works_with_R("3.2.2",
             data.table="1.9.6",
             bcp="4.0.0",
             EBS="3.0",
             postCP="1.7.2") # Archived.

load("data.by.type.RData")

key.dt <- data.table(label.i=5)
setkey(data.by.type$probes, label.i)
setkey(data.by.type$segments, label.i)
probes <- data.by.type$probes[key.dt]
setkey(probes, probeStart, probeEnd)
probes[, probe.i := seq_along(probeStart) ]
segs <- data.by.type$segments[key.dt]
setkey(segs, segmentStart, segmentEnd)
over.dt <- foverlaps(probes, segs, nomatch=0L)
stopifnot(nrow(over.dt) == nrow(probes))
last.dt <- over.dt[, list(
  last.i=max(probe.i)
  ), by=segmentStart]
last.index.vec <- last.dt[-.N, last.i]
data.vec <- probes$logratio

## Takes 8GB of RAM and is slow...
set.seed(1)
data.vec <- rnorm(17147)
fit <- EBSegmentation(data.vec, model=2, Kmax=5)

set.seed(1)
data.vec <- c(rpois(125,1),rpois(100,5),rpois(50,1),rpois(75,5),rpois(50,1))
fit <- EBSegmentation(data.vec,model=1,Kmax=20)
## Li: a matrix of size Kmax*(length+1). Element [i,j] is the
## log-probability of interval [1,j[ being segmented in i segments

## Col: a matrix of size (length+1)*Kmax. Element [i,j] is the
## log-probability of interval [i,n] being segmented in j segments
model.segments <- 5
break.i <- 2
t.vec <- (break.i + 1):(length(data.vec) - (model.segments - break.i))
Aux <- fit@Li[break.i, t.vec] + fit@Col[t.vec, model.segments - break.i]
ma <- max(Aux)
Aux2 <- exp(Aux - ma)
prob.vec <- c(
  rep(0, break.i),
  Aux2/sum(Aux2),
  rep(0, model.segments - break.i)
  )

## bcp has linear time complexity but only returns one probability
## vector, so can not be used to compute an error band for each
## breakpoint.
fit <- bcp(data.vec)
plot(fit$posterior.prob, type="l")
plot(fit)

## I tried versions of postCP in the CRAN archive, from v1.0 to
## v1.7.2, but none of them work. They all give a matrix of zeros for
## fit$cp.est. I emailed Nuel and Luong and asked them if they can
## help.
fit <- postCP(data.vec, last.index.vec, model=2, ci=0.9)
print(fit)

