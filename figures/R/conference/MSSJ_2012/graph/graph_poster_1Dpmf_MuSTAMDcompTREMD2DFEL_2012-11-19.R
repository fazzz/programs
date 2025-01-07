#source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange.R")

dir01 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase1 <- paste(dir01,"/e_CG-FG_NH_2012-07-27",sep="")

dir02 <- "~/calspa/refcalc/REMD/AD"
dirbase2 <- paste(dir02,"/s_REVAC_2012-07-23_",ff,sep="")

dir03 <- "~/calspa/refcalc/UmbSam/AD"
dirbase3 <- paste(dir03,"/s_UmbSam_vac_2012-11-12_",ff2,sep="")

nTZs<-1

title=name.title

if ( phsiflag == "psi" )
  label.x<-expression(paste(phi))
if ( phsiflag == "phi" )
  label.x<-expression(paste(psi))

label.y="pmf"

xrange <- c(-3.14,3.14,4)
xrange.axis <- c(-2,3,5)
yrange.axis <- c(0,height,10)
height.ext <- height * 1.2
yrange <- c(0,height.ext,10)

fact.x <- 1
fact.y <- 1
ave.x <- 1

file.name <- paste(name.out,".tiff",sep='')
cat(file.name,'\n')
#tiff(file.name,width=400,height=250)
#tiff(file.name,width=750,height=350)
tiff(file.name,width=700,height=500)
#tiff(file.name,width=500,height=500)
#tiff(file.name,width=500,height=400)
#tiff(file.name,width=400,height=350)

par(mfrow=c(4,5))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

dx <- 6.28 / num
cat(dx,'\n')

s <- 1
for (t in 1:20) {
  name <- NULL
  name[1] <- paste(dirbase1,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname1,"/nEX=",numEX1,"/fq=",TLbase1,"/pmf/1D/pmf_pymbar_1_TAA=",TAA,"_TCG=",TCG,"_TZ=",TZs[1],"_KZAAo=",KZAAo,"_KZCGo=",KZCGo,"_",AACG,phsiflag,"@",t,sep="")
  name[2] <- paste(dirbase2,"/",pname2,"/nEX=",numEX2,"/freq=",TLbase2,"ps","/pmf/1D/pmf_p1",phsiflag,"@",t,sep="")
  if ( phsiflag=="phi" ) {
    xy <- "x"
  }
  if ( phsiflag=="psi" ) {
    xy <- "y"
  }
  name[3] <- paste(dirbase3,"/",pname3,"/pmf/1D/pmf_1D_UmbMD_vac_",TLns,"ns_",xy,"=",t,".txt",sep="")

  cat(name[1],'\n')
  cat(name[2],'\n')

  id.xs <- c(1,1,1)
  id.ys <- c(2,2,2)
  ids.ys <- c(3,3,3)
  iro <- c(2,3,1)
  senshu <- c(1,1,1)
#  tenshu <- c(1,1)
  tenshu <- c(3,3,19)
  hutosa <- c(1,1,1)
  is.leg <- 0

  plmGeneralwsrange(data.names=name,
                    sd.names=name,
                    id.ys=id.ys,
                    ids.ys=ids.ys,
                    is.sen=rep(0,20),
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    warrow="T")
  val = -3.0 + 0.3 * (t -1)
  txt <- paste(phsiflag,"=",val,sep='')
  text(0,19.0,txt,cex=2.0)
  box(lwd=2.0)
  
  k <- t %% numy
  if ( k == 1 ) {
    axis(2,yaxp=yrange.axis,lwd=2.0,cex.axis=1.5)
    mtext(outer=T,label.y,side=2,line=4.0,cex=1.5)
  }
  
  if ( t > numy*(numx-1) ) {
    axis(1,xaxp=xrange.axis,lwd=2.0,cex.axis=1.5)
    mtext(outer=T,label.x,side=1,line=4.0,cex=1.5)
#    s <- s+1
  }
}

dev.off()
