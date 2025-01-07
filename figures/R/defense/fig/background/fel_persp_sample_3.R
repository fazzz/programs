require(MASS)
set.seed(125)
x <- rnorm(150,mean=3*rbinom(150,prob=.5,size=1),sd=1)
y <- rnorm(150,mean=4*rbinom(150,prob=.5,size=2),sd=1)
d <- kde2d(x,y,n=50)

name.out <- "/home/yamamori/defense/eps/fel_sample_3"
file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.0,height=4.0,horizontal=FALSE,onefile=FALSE,paper="special")

kde2dplot <- function(d,                # a 2d density computed by kde2D
                      ncol=50,          # the number of colors to use
                      zlim=c(0,max(z)), # limits in z coordinates
                      nlevels=20,       # see option nlevels in contour
                      theta=30,         # see option theta in persp
#                      phi=30)           # see option phi in persp
#                       phi=160)           # see option phi in persp
#                       phi=90)           # see option phi in persp
#                       phi=120)           # see option phi in persp
#                       phi=60)           # see option phi in persp
#                       phi=240)           # see option phi in persp
#                       phi=200)           # see option phi in persp
                      phi=190)           # see option phi in persp
                      {
z   <- d$z
nrz <- nrow(z)
ncz <- ncol(z)

couleurs  <- tail(topo.colors(trunc(1.4 * ncol)),ncol)
#couleurs  <- tail(heat.colors(trunc(1.4 * ncol)),ncol)
fcol      <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol      <- fcol[-nrz,-ncz]

persp(d,col=fcol,zlim=zlim,theta=theta,phi=phi,box="FALSE",axes="TRUE",xlab="",ylab="",zlab="")

}

