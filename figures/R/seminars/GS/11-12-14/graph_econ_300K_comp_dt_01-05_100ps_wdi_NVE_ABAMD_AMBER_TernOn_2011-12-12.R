source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,8,8)
xrange.axis <- c(0,7,7)
yrange <- c(-6,1,7)
yrange.axis <- c(-6,0,6)

prot.name <- NULL
prot.name2 <- NULL

prot.name <- c("HP35","UBIQ")         
prot.name2 <- c("HP35","UBIQ")
num.prot <- length(prot.name)

numini=10

for ( i in 1:num.prot  ) {
  dirbase.name[i] <- paste("/home/yamamori/calspa/TAMD/",prot.name2[i],sep="")
}

for ( k in 1:4  ) {
  sen[k]<-c(3)
  n<-1+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

par(mfrow=c(1,num.prot))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.prot) {
  file.names[1]=paste(dirbase.name[i],"/econ_vdt_wdi_NVE_AMBER_TermOn_11-11-20_pc6/",prot.name[i],"v_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(dirbase.name[i],"/econ_vdt_wdi_NVE_AMBER_TermOn_11-11-20_mvV/",prot.name[i],"v_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')  
  file.names[3]=paste("/home/yamamori/calspa/refcalc/",prot.name2[i],"/econ_wdi_NVE_100ps/",prot.name[i],"v_econ_",temp[1],"_100ps_1-10_2ntc",".econ.rmsd.av",sep='')

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
  text(5,0.5,prot.name[i])
  box(lwd=2.0)
  ###################################################
  ## if (i==1) {                                   ##
  ##   legend(2.5,-5.5,legend=c("pc6'"),bty="n",   ##
  ##          lty=c(1),cex=1.0,y.intersp=2,        ##
  ##          col=iro[1],ncol=leg.col)             ##
  ##   legend(2.5,-5.0,legend=c("vV"),bty="n",     ##
  ##          lty=c(1),cex=1.0,y.intersp=2,        ##
  ##          col=iro[2],ncol=leg.col)             ##
  ##   legend(2.5,-4.5,legend=c("CMD"),bty="n",    ##
  ##          lty=c(1),cex=1.0,y.intersp=2,        ##
  ##          col=iro[3],ncol=leg.col)             ##
  ## }                                             ##
  ###################################################

  if (file.access(file.names[3]) == 0){
    a <- read.table(file.names[3],header=TRUE)
    cat(a[2,2])
    lines(c(1,6),c(a[2,2],a[2,2]),lwd=2,lty="dashed",col="gray")
  }
  
  axis(1,xaxp=xrange.axis,lwd=2.0)
  mtext(outer=T,label.x,side=1,line=5.0,cex=0.8)
  
  if (i==1) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
    mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
  }
}

name.out=paste("~/seminars/GS/11-12-14/graph_econ_300K_comp_dt_01-05_100ps_wdi_NVE_ABAMD_AMBER_TernOn_2011-12-14",sep='')
OutTiff(name.out,w_val=600,h_val=300)
