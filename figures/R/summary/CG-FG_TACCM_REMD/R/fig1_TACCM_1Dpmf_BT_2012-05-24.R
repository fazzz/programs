
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(-180,180,12)
xrange.axis1 <- c(-180,180,6)
#xrange.axis2 <- c(-150,180,11)

yrange <- c(-6.0,1.0,2)    
yrange.axis1 <- c(-6.0,1.0,7)
#yrange.axis2 <- c(-6.0,1.0,7)

label.x <- expression(paste("dihedral angle "))
label.y <- expression(paste("pmf [kcal/mol]"))

name.prot <- c("BT")
num.prot <- length(name.prot)

tau <- 0.1

T    <-  c( "300"   ) 
T2    <- c(  300    ) 
nT   <- length(T)

TB    <-  c( "600" )
nTB   <- length(TB)

#KZ <- c( "60.00", "600.00")
KZ <- c( "600.00" )
nKZ <- length(KZ)

mZ <- c( "100000.00", "500000.00" )
nmZ <- length(mZ)

#par(mfrow=c(nKZ,nmZ))
par(mfrow=c(1,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

filex.names<-NULL
filey.names<-NULL
  
for (n in 1:nT) {
  for (p in 1:nKZ) {
    for (q in 1:nmZ) {

      filex.names<-NULL
      filey.names<-NULL
      label.names<-NULL
      sen <- NULL
      senshu <- NULL
      id.ys <- NULL
      id.xs <- NULL
      iro <-NULL

      for (o in 1:nTB) {

        sen[o]<-c(1)
        senshu[o] <- 1
        id.ys[o]<-c(2)
        iro[o]<-c(1+o)

        dirbase.name <- paste("/home/yamamori/calspa/TACCM/BT/equ_vT_vpara_Amber_TACCM_MD_NH_2012-02-27/tau=",tau,"//TB=",TB[o],"/KZ=",KZ[p],"/mZ=",mZ[q],"/anl/",sep="")

        filex.names[o]<-paste(dirbase.name,name.prot,"_equ_T=",T[n],".Dhist",sep='')
        filey.names[o]<-paste(dirbase.name,name.prot,"_equ_T=",T[n],".Dhist",sep='')
        label.names[o]<-paste("TB=",TB[o],sep='')
        
        ave.x[o] <- 1
      }

      o <- nTB +1
      sen[o]<-c(1)
      senshu[o] <- 2
      id.ys[o]<-c(2)
      iro[o]<-"black"

      
      filex.names[o]<-paste("~/calspa/TACCM/BT/pmf_ns=20_tl=100000000_fc=10.0_temp=300.txt")
      filey.names[o]<-paste("~/calspa/TACCM/BT/pmf_ns=20_tl=100000000_fc=10.0_temp=300.txt")
      label.names[o]<-paste("ref",sep='')                                                                   
      ave.x[o] <- 1                                                                                         

      plmGeneralwsrangewdrangewdiffx(datax.names=filex.names,
                                     datay.names=filey.names,
                                     label.names=label.names,
                                     id.ys=id.ys,
                                     label.size=0.5,axis.size=2.0,
                                     iro=iro,axis.ft="F",is.header="T",
                                     sdiro=iro,
                                     xrange=xrange,yrange=yrange,
                                     sdyrange=yrange,
                                     is.sen=sen,width=10.0,
                                     ,warrow="F",ave.x=ave.x)
               
      box(lwd=2.0)
      if ( p == 1 && q==1 ) {
#        text(150,1.0,"(a)")
      }
      if ( p == 1 && q==2 ) {
 #       text(150,1.0,"(b)")
      }
      if ( p == 2 && q==1 ) {
  #      text(150,1.0,"(c)")
      }
      if ( p == 2 && q==2 ) {
   #     text(150,1.0,"(d)")
      }
      
      if ( p == nKZ && q==1 ) {
        axis(1,xaxp=xrange.axis2,lwd=2.0)
        mtext(outer=T,label.x,side=1,line=5.0,cex=1.0)
      }
      else if ( p == nKZ && q==2 ) {
        axis(1,xaxp=xrange.axis1,lwd=2.0)
        mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
      }
      
#      if ( q == 1 && p == 1 ) {
      axis(2,yaxp=yrange.axis1,lwd=2.0)
      mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
#     }
      ######################################################
      ## else if ( q == 1 && p == 2 ) {                   ##
      ##   axis(2,yaxp=yrange.axis2,lwd=2.0)              ##
      ##   mtext(outer=T,label.y,side=2,line=3.0,cex=1.0) ##
      ## }                                                ##
      ######################################################
    }
  }
}

name.out=paste("~/summary/CG-FG_TACCM_REMD/R/fig1_TACCM_1Dpmf_BT_2012-05-24",sep='')
#OutTiff(name.out,w_val=700,h_val=500)
OutTiff(name.out,w_val=400,h_val=300)
