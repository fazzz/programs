source("~/Rspa/ini.R")
source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,7,7)
yrange <- c(-7,0,7)

name.res <- c("ALA","ASP","ASN","ARG","CYS","GLN", "GLY","GLU","HIE","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TYR","TRP","VAL")

num.res <- length(name.res)

dirbase <- "/home/yamamori/TAMD/dipeptides/"

for ( k in 1:4  ) {
  sen[k]<-c(3)
  n<-2+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  tenshu[k]<-c(20)
}

par(mfrow=c(4,5))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.res) {
  file.names[1]=paste(dirbase.name,name.res[i],"/econ_wdi_NVE_AMBER_100ps_termon/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10.econ.rmsd.av",sep='')
  file.names[2]=paste(dirbase.name,name.res[i],"/econ_wdi_NVE_AMBER_ai_debug_100ps_termon/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10.econ.aved.av",sep='')  
  file.names[3]=paste("/home/yamamori/calspa/refcalc/dipep/",name.res[i],"/econ_wdi_NVE_100ps/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10_2ntc",".econ.aved.av",sep='')

  plmGeneralwsrange(data.names=file.names,
                      sd.names=file.names,
#                      label.names=label.names,
                      id.ys=id.ys,
                      ids.ys=ids.ys,
                      label.size=0.5,axis.size=2.0,
                      iro=iro,axis.ft="F",is.header="T",
                      sdiro=iro,
                      xrange=xrange,yrange=yrange,
                      sdyrange=yrange,
                      is.sen=sen,width=10.0,
                      ,warrow="T")
  text(5,-6,name.res[i])
  box(lwd=2.0)
  if (i==1) {
    ###############################################
    ## legend(3,-5.0,legend=c("vV-H'"),bty="n",  ##
    ##        lty=c(1),cex=1.0,y.intersp=2,      ##
    ##        col=iro[3],ncol=leg.col)           ##
    ## legend(3,-4.0,legend=c("vV-H"),bty="n",   ##
    ##        lty=c(1),cex=1.0,y.intersp=2,      ##
    ##        col=iro[4],ncol=leg.col)           ##
    ## legend(3,-3.0,legend=c("pc6-H"),bty="n",  ##
    ##        lty=c(1),cex=1.0,y.intersp=2,      ##
    ##        col=iro[1],ncol=leg.col)           ##
    ## legend(3,-2.0,legend=c("pc6-H'"),bty="n", ##
    ##        lty=c(1),cex=1.0,y.intersp=2,      ##
    ##        col=iro[2],ncol=leg.col)           ##
    ###############################################
    legend(0.2,-2.0,legend=c("pc6'"),bty="n",  
           lty=c(1),cex=1.0,y.intersp=2,      
           col=iro[1],ncol=leg.col)           
    legend(0.2,-1.0,legend=c("vV"),bty="n",   
           lty=c(1),cex=1.0,y.intersp=2,      
           col=iro[2],ncol=leg.col)           
    legend(0.2,-0.0,legend=c("CMD"),bty="n",  
           lty=c(1),cex=1.0,y.intersp=2,      
           col=iro[3],ncol=leg.col)           
  }

  a <- read.table(file.names[3],header=TRUE)
  cat(a[2,2])
  lines(c(1,6),c(a[2,2],a[2,2]),lwd=2,lty="dashed",col=iro[4])

  
  if (i==16 || i==17 || i==18 || i==19 || i==20) {
    axis(1,xaxp=xrange,lwd=2.0)
    mtext(outer=T,label.x,side=1,line=5.0,cex=0.8)
  }
  
  if (i==1 || i==6 || i==11 || i==16 ) {
    axis(2,yaxp=yrange,lwd=2.0)
    mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
  }

}
name.out=paste("~/papers/TAMD/fig2.a.econ_NVE_20_aminoacids",sep='')
OutTiff(name.out,w_val=880,h_val=660)
