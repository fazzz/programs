felwFillConMapwpath <- function (data.name.pmf,data.name.path,label.x="",label.y="",title="",level=c(0,2,4,6,8,10,12,14,16,18,20),kt="") {
  data <- read.table(data.name.pmf)
  n.XtimesY <- nrow(data)
  n.X <- 1
  temp.X <- data[1,1]
  for ( i in 1:n.XtimesY) {
    if(data[i,1] != temp.X ) {
      n.X <- n.X + 1
      temp.X <- data[i,1]
    }
  }
  n.Y <- n.XtimesY / n.X
  n.bin <- sqrt(nrow(data))
  new.data1 <- matrix(data$V3,nrow=n.Y,ncol=n.X)
  x <- data$V1[seq(1,n.XtimesY,n.Y)]
  y <- data$V2[1:n.Y]
  new.data <- t(new.data1)
  y <- y
  FillConMap(x,y,new.data,color=topo.colors,level=level, #nlevels=10,
             plot.title=title(main=title,xlab=label.x,ylab=label.y,line=c(3.0,3.0),cex=2.0)
#             ,col=gray((8:0)/8)
             )

  a <- read.table(data.name.path)

#  lines(a[,1],a[,2],lwd=4,col="red")
  points(a[,1],a[,2],lwd=0.5,col="red")
}
