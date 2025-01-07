fel_persp <- function (data.name,label.x="",label.y="",title="",level=c(0,2,4,6,8,10,12,14,16,18,20),kt="",xrange=xrange,yrange=yrange,mar=c(5,5,2,2),xaflag="T",yaflag="F",norm="F",plot.axis) {
  data <- read.table(data.name)
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
  x <- x*fact.x

  y <- y*fact.y

  if (norm=="T") {
    Nx <- length(x)
    Ny <- length(y)

    x <- seq(-1.0,1.0,length=Nx)
    y <- seq(-1.0,1.0,length=Ny)    
  }
  
  new.data <- t(new.data1)

  y <- y
  par(mar=mar)

  persp(x,y,new.data,
        axes="FALSE",
#        axes="TRUE",
        box="TRUE",
#        box="FALSE",
#        theta=60,
#        phi=60,
        theta=80,
        phi=10,
#        border=NA,
#        shade=0.1,
#        col=rainbow(50),
        col=heat.colors(10),
#        col=grey(0:10/11),
        expand=0.3,
        xlab="",ylab="",zlab=""
#        plot.title=title(main=title,xlab=label.x,ylab=label.y,line=c(3.0,3.0),cex=2.0),
#        xrange=xrange,yrange=yrange
        )
}
