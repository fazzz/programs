source("~/Rspa/ini.R")
source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,3,3)
yrange <- c(-3,-1,2)

name.prot <- c("G20")
name.prot2 <- c("GLY_20")
numres.prot <- c(20,35)
num.prot <- length(name.prot)
dirbase.name <- "/home/yamamori/calspa/refcalc/"

for ( k in 1:4  ) {
  sen[k]<-c(3)
  n<-2+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  tenshu[k]<-c(20)
}

par(mfrow=c(1,num.prot))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

type <- c("rmsd","aved")

file.names<-NULL
for (t in 1:2) {
  cat("numres ave std","\n",file = outfilname)
  for (i in 1:num.prot) {
    file.names[1]=paste(dirbase.name,name.prot[i],"/econ_wdi_NVE_100ps/",name.prot2[i],"v_econ_300_100ps_1-10_1ntc",".econ.",type[t],".av2",sep='')
    file.names[2]=paste(dirbase.name,name.prot[i],"/econ_wdi_NVE_100ps/",name.prot2[i],"v_econ_300_100ps_1-10_2ntc",".econ.",type[t],".av2",sep='')

    plmGeneralwsrange(data.names=file.names,
                      sd.names=file.names,
                      id.ys=id.ys,
                      ids.ys=ids.ys,
                      label.size=0.5,axis.size=2.0,
                      iro=iro,axis.ft="F",is.header="T",
                      sdiro=iro,
                      xrange=xrange,yrange=yrange,
                      sdyrange=yrange,
                      is.sen=sen,width=10.0,
                      ,warrow="T")
    text(5,-3,name.prot[i])
    box(lwd=2.0)

    ###############################################################
    ## legend(0.1,-2.9,legend=c("without constraint"),bty="n",   ##
    ##        lty=c(1),cex=1.5,y.intersp=2,                      ##
    ##        col=iro[1],ncol=leg.col)                           ##
    ## legend(0.1,-3.4,legend=c("with bond constraint"),bty="n", ##
    ##        lty=c(1),cex=1.5,y.intersp=2,                      ##
    ##        col=iro[2],ncol=leg.col)                           ##
    ###############################################################
    
    axis(1,xaxp=xrange,lwd=2.0,cex.axis=1.5)
    mtext(outer=T,label.x,side=1,line=4.0,cex=2.0)
  
    if (i==1 ) {
      axis(2,yaxp=yrange,lwd=2.0,cex.axis=1.5)
#      mtext(outer=T,label.y,side=2,line=4.0,cex=2.0)
    }
  }
  name.out=paste("~/papers/TAMD/fig.tyuukan_1_econ_ref_",type[t],sep='')
  OutTiff(name.out,w_val=500,h_val=400)
}


