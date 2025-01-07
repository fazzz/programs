
felwFillConMapwrangeworwoaxiswpoints <- function (data.name,label.x="",label.y="",title="",level=c(0,2,4,6,8,10,12,14,16,18,20),kt="",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(1.0,1.0,0.2,0.2),plot.axis,xaflag="T",yaflag="T",norm="F",labelflag="T",pdata,iro) {
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
#  par(mar=mar)
  par(mai=mai)
  FillConMapwrangeworwoaxis(x,y,new.data,color=topo.colors,level=level, #nlevels=10,
                            plot.title=title(main=title,xlab=label.x,ylab=label.y,line=c(3.0,3.0),cex=2.0),
                            xrange=xrange,yrange=yrange,plot.axis=plot.axis,wxaxis=xaflag,wyaxis=yaflag,
                            xrange.axis=xrange.axis,yrange.axis=yrange.axis,labelflag=labelflag,
                            )

  a <- read.table(pdata)
  a <- a*fact.p
  n <- length(a)+1
  
  for ( i in 1:n ) {                        
    points(a[i,1],a[i,2],lwd=1.0,col=iro[i])
  }                                         
  
}
