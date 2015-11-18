works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/namedCapture@a31a38be12d4dad4aec6257a47c0e2307679589c")

pair.df <- read.table("pairs.txt", header=TRUE, as.is=TRUE)

pattern <- paste0(
  "(?<prefix>.*?)",
  "/profile/",
  "(?<profile>.*?)",
  "/",
  "(?<chr>[^/]+)")

match.mat <- str_match_named(pair.df$chrom.url, pattern)

profile.mat <- unique(match.mat[, c("prefix", "profile")])
data.by.profile <- list()
for(profile.i in 1:nrow(profile.mat)){
  match.vec <- profile.mat[profile.i,]
  profile.id <- match.vec[["profile"]]
  cat(sprintf("%4d / %4d profiles %s\n",
              profile.i, nrow(profile.mat), profile.id))
  url.vec <- structure(c(
    sprintf("%s/secret/%s.bedGraph.gz",
      match.vec[["prefix"]],
      profile.id),
    sprintf("%s/export/tdhock5@gmail.com/%s/%s/csv/",
      match.vec[["prefix"]],
      profile.id,
            c("regions", "copies"))),
                       names=c("probes.bedGraph.gz",
                         "regions.csv", "segments.csv"))
  profile.dir <- file.path("profiles", profile.id)
  dir.create(profile.dir, showWarnings=FALSE, recursive=TRUE)
  profile.data <- list()
  for(local.base in names(url.vec)){
    local.path <- file.path(profile.dir, local.base)
    if(!file.exists(local.path)){
      u <- url.vec[[local.base]]
      download.file(u, local.path)
    }
    data.type <- sub("[.].*", "", local.base)
    new.names <- paste0(sub("s$", "", data.type), c("Start", "End"))
    dt <- if(local.base == "probes.bedGraph.gz"){
      df <- read.table(local.path, skip=1,
                 col.names=c("chrom", new.names, "logratio"))
      data.table(df)
    }else{
      tmp <- fread(local.path)
      setnames(tmp, c("min", "max"), new.names)
      tmp[, chrom := paste0("chr", chromosome)]
      tmp
    }
    setkeyv(dt, c("chrom", new.names))
    profile.data[[data.type]] <- dt
  }
  data.by.profile[[profile.id]] <- profile.data
}

data.by.type.lists <- list()
for(label.i in 1:nrow(match.mat)){
  match.vec <- match.mat[label.i,]
  label.row <- pair.df[label.i,]
  label.chrom <- paste0("chr", match.vec[["chr"]])
  profile.id <- match.vec[["profile"]]
  profile.data <- data.by.profile[[profile.id]]
  chrom.data <- list()
  for(data.type in names(profile.data)){
    dt <- profile.data[[data.type]]
    data.by.type.lists[[data.type]][[paste(label.i)]] <-
      chrom.data[[data.type]] <- 
      data.table(label.i, dt[label.chrom])
  }
  breakpoint <- chrom.data$segments[-1, segmentStart]
  break.dt <- data.table(label.i, breakpoint, breakpoint0=breakpoint)
  setkey(break.dt, breakpoint, breakpoint0)
  setkey(chrom.data$regions, regionStart, regionEnd)
  chrom.data$regions[, region.i := seq_along(regionStart)]
  over.dt <- foverlaps(break.dt, chrom.data$regions, nomatch=0L)
  setkey(over.dt, region.i)
  larger <- over.dt[J(region.i=label.row$larger), breakpoint]
  smaller <- over.dt[J(region.i=label.row$smaller), breakpoint]
  data.by.type.lists$breakpoints[[paste(label.i)]] <- break.dt
  data.by.type.lists$labels[[paste(label.i)]] <-
    data.table(label.i, smaller, larger)
}

data.by.type <- list()
for(data.type in names(data.by.type.lists)){
  data.by.type[[data.type]] <- do.call(rbind, data.by.type.lists[[data.type]])
}

save(data.by.type, file="data.by.type.RData")
