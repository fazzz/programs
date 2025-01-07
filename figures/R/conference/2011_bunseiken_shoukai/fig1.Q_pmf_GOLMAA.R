source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

dirbase.name <- paste("/home/yamamori/calspa/GOLMAA/SH3/equ_1nlo_vT_vpara_ABAMD_GOLMAA_PROTEINS2008_AAMD_NH_2012-01-05/tau=1.0/ep=",ep[1],"/coff=",cutoff[1],"/",sep="")

par(mfrow=c(1,2))
par(mar=c(0.0,6.0,6.0,0.0))
par(oma=c(7.5,0.0,0.0,2.0))


leg.pos="top"
is.leg <- 0

xrange <- c(300,400,10)
xrange.axis <- c(300,400,10)

label.x <- expression(paste("T (K)"))

name.prot <- c("SH3_1nlo")
num.prot <- length(name.prot)

ep <- c(  "0.3" )
cutoff <- c( "4.7" )

filex.names<-NULL
filey.names<-NULL
yrange <- c(0,5.0,5)
yrange.axis <- yrange
label.y <- expression(paste("Cv"))
  

leg.pos="topleft"
is.leg <- 0

xrange <- c(0,10000000,10)
xrange.axis <- c(0,10000000,10)

label.x <- expression(paste("time step "))

name.prot <- c("SH3")
num.prot <- length(name.prot)

T    <-  c( "355" )
nT   <- length(T)

ep <- c( "0.3" )
nep <- length(ep)

nb <- c( "3" )
nnb <- length(nb)

cutoff <- c(  "4.7" )
ncutoff <- length(cutoff)

filex.names<-NULL
filey.names<-NULL
yrange <- c(-0.05,1.05,10)
yrange.axis <- c(0,1,10)
label.y <- expression(paste("Q"))
  
for (k in 1:nnb) {
  for (i in 1:nep) {
    for (l in 1:ncutoff) {
      for (j in 1:nT) {
        sen[j]<-c(1)
        senshu[j] <- 1
        iro[j]<-c(1+j)
        id.ys[j]<-c(1)

        dirbase.name <- paste("/home/yamamori/calspa/GOLMAA/SH3/equ_1nlo_vT_vpara_ABAMD_GOLMAA_PROTEINS2008_AAMD_NH_2012-01-05/tau=1.0/ep=",ep[i],"/coff=",cutoff[l],"/",sep="")
        
        filex.names[j]<-paste(dirbase.name,name.prot,"_equ_T=",T[j],".out",sep='')
        filey.names[j]<-paste(dirbase.name,"anl/",name.prot,"_equ_T=",T[j],".qca",sep='')
        label.names[j]<-paste("T=",T[j],sep='')
        ave.x[j] <- 1
      }
    
      plmGeneralwsrangewdrangewdiffx(datax.names=filex.names,
                                     datay.names=filey.names,
                                     label.names=label.names,
                                     sd.names=file.names,
                                     id.ys=id.ys,
                                     ids.ys=ids.ys,
                                     label.size=0.5,axis.size=2.0,
                                     iro=iro,axis.ft="F",is.header="F",
                                     sdiro=iro,
                                     xrange=xrange,yrange=yrange,
                                     sdyrange=yrange,
                                     is.sen=sen,width=10.0,
                                     ,warrow="F",ave.x=ave.x)
    
      box(lwd=3.5)
  
      if ( k == nnb ) {
        axis(1,xaxp=xrange.axis,lwd=3.5)
        mtext(label.x,side=1,line=4.0,cex=1.5)
      }
      
      if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=3.5)
        mtext(label.y,side=2,line=3.0,cex=1.5)
      }
    }
  }
}

leg.pos="top"
is.leg <- 0

xrange <- c(-0.05,1.1,10)
xrange.axis <- c(0,1,10)

label.x <- expression(paste("Q "))

name.prot <- c("SH3")
num.prot <- length(name.prot)

T    <-  c(  "352", "355", "360" ) 
T2    <- c(   352 ,  355 ,  360  ) 
nT   <- length(T)

ep <- c(  "0.3" )
cutoff <- c( "4.7" )
width <- "0.05"

dirbase.name <- paste("/home/yamamori/calspa/GOLMAA/SH3/equ_1nlo_vT_vpara_ABAMD_GOLMAA_PROTEINS2008_AAMD_NH_2012-01-05/tau=1.0/ep=",ep[1],"/coff=",cutoff[1],"/",sep="")

filex.names<-NULL
filey.names<-NULL
yrange <- c(0,7,7)
yrange.axis <- yrange
label.y <- expression(paste("pmf"))
  
for (j in 1:nT) {
  sen[j]<-c(1)
  senshu[j] <- 1
  iro[j]<-c(1+j)
  id.ys[j]<-c(2)
  
  filex.names[j]<-paste(dirbase.name,"pmf/pmf_T=",T[j],"_",width,sep='')
  filey.names[j]<-paste(dirbase.name,"pmf/pmf_T=",T[j],"_",width,sep='')
  label.names[j]<-paste("T=",T[j],sep='')
  ave.x[j] <- 1
}
    
plmGeneralwsrangewdrangewdiffx(datax.names=filex.names,           
                               datay.names=filey.names,           
                               label.names=label.names,           
                               sd.names=file.names,               
                               id.ys=id.ys,                       
                               ids.ys=ids.ys,                     
                               label.size=0.5,axis.size=2.0,      
                               iro=iro,axis.ft="F",is.header="F", 
                               sdiro=iro,                         
                               xrange=xrange,yrange=yrange,       
                               sdyrange=yrange,                   
                               is.sen=sen,width=10.0,             
                               ,warrow="F",ave.x=ave.x)           
    
box(lwd=3.5)


axis(1,xaxp=xrange.axis,lwd=3.5)       
mtext(label.x,side=1,line=4.0,cex=1.5) 
                                       
axis(2,yaxp=yrange.axis,lwd=3.5)       
mtext(label.y,side=2,line=3.0,cex=1.5) 

name.out=paste("~/gakkai/2011_bunseiken_shoukai/fig1.Q_pmf_GOLMAA",sep='')
OutTiff(name.out,w_val=600,h_val=300)
