works_with_R("3.2.2",
             data.table="1.9.6",
             "tdhock/animint@4257e8cf76eb5021a98010b6629b779a4f383b24",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2")

load("data.by.type.RData")

ann.colors <- c(
    "1breakpoint"="#ff7d7d",
    "0breakpoints"='#f6f4bf'
  )

gg <- 
ggplot()+
  geom_tallrect(aes(xmin=regionStart/1e6, xmax=regionEnd/1e6,
                    fill=annotation),
                color="grey",
                alpha=0.5,
                data=data.by.type$regions)+
  scale_fill_manual(values=ann.colors)+
  geom_point(aes(probeStart/1e6, logratio),
             shape=1,
             data=data.by.type$probes)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(label.i ~ ., scales="free", labeller=function(var, val){
    paste("label", val)
  })+
  xlab("position on chromosome (mega bases)")+
  geom_segment(aes(smaller/1e6, Inf,
                   xend=larger/1e6, yend=Inf),
               color="blue",
               size=3,
               data=data.by.type$labels)+
  geom_point(aes(larger/1e6, Inf),
               color="blue",
               shape=1,
               size=3,
               data=data.by.type$labels)+
  coord_cartesian(xlim=range(data.by.type$probes$probeStart/1e6))+
  ylab("logratio (approximate DNA copy number)")+
  ggtitle("Blue circle on larger confidence band")

png("figure-labels.png", height=7, width=100, res=200, units="in")
print(gg)
dev.off()

setkey(data.by.type$probes, label.i)
setkey(data.by.type$regions, label.i)
setkey(data.by.type$labels, label.i)
## TODO: figure zoomed just to arrow.
for(label.i in unique(data.by.type$probes$label.i)){

  key.dt <- data.table(label.i)
  regions <- data.by.type$regions[key.dt]
  lim.vec <-
    c(min(regions$regionStart),
      max(regions$regionEnd))/1e6

  gg <-
    
    ggplot()+
      geom_tallrect(aes(xmin=regionStart/1e6, xmax=regionEnd/1e6,
                        fill=annotation),
                    color="grey",
                    alpha=0.5,
                    data=regions)+
      scale_fill_manual(values=ann.colors)+
      geom_point(aes(probeStart/1e6, logratio),
                 shape=1,
                 data=data.by.type$probes[key.dt])+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      xlab("position on chromosome (mega bases)")+
      geom_segment(aes(smaller/1e6, Inf,
                       xend=larger/1e6, yend=Inf),
                   color="blue",
                   size=3,
                   data=data.by.type$labels[key.dt])+
      geom_point(aes(larger/1e6, Inf),
                 color="blue",
                 shape=1,
                 size=3,
                 data=data.by.type$labels[key.dt])+
      coord_cartesian(xlim=lim.vec)+
      ylab("logratio (approximate DNA copy number)")+
      ggtitle(paste("Blue circle on larger confidence band,",
                     "label", label.i))

  png(sprintf("figure-labels-%d.png", label.i),
      height=7, width=20, res=200, units="in")
  print(gg)
  dev.off()
}  
