source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

xrange <- c(0,8,8)
xrange.axis <- c(0,7,7)
yrange <- c(-7,0,7)
yrange.axis <- c(-7,-1,6)

name.res <- c( "ASN", "GLN", "LEU","SER", "LYS","ILE","ALA", "CYS", "MET", "THR", "VAL", "ARG", "ASP", "GLY", "HIE","PHE","PRO", "THR","TYR","TRP" )

numini=10

num.res <- length(name.res)

dirbase.name <- "/home/yamamori/calspa/TAMD-s-ivy3/TAMD/dipeptides/"
dirbase21.name <-paste("~/calspa/TAMD-s-ivy3/TAMD/mod_parm/dipep",sep='')
dirbase22.name <- paste("econ_wdi_NVE_AMBER_TermOn_wtune_2012-02-29_pc6",sep='')

name.out=paste("~/thesis/TAMD/eps/fig_Chapt5_econdependence_dt=1-5fs_T=300K_comp_TAMDvsCMD_20dipeptides.eps",sep="")
postscript(name.out,width=6.0,height=6.0,horizontal=FALSE,onefile=FALSE,paper="special")

hutosa <- 1

id.ys <- NULL
ids.ys <- NULL

for ( k in 1:4  ) {
  sen[k]<-c(1)
  n<-1+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

par(mfrow=c(4,5))
par(oma=c(7.5,4.0,2.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

file.names<-NULL
for (i in 1:num.res) {
  file.names[1]=paste(dirbase.name,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_11-11-17_pc6/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(dirbase.name,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_11-11-17_mvV/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')  
  file.names[3]=paste("/home/yamamori/calspa/refcalc/dipep/",name.res[i],"/econ_wdi_NVE_100ps/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10_2ntc",".econ.rmsd.av",sep='')

  senshu[1]<-1
  senshu[2]<-1
  
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
  box(lwd=2.0)
  text(5,-6,name.res[i])

  a <- read.table(file.names[3],header=TRUE)
  cat(a[2,2])
  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")

  if (i==16 || i==17 || i==18 || i==19 || i==20) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
#    mtext(outer=T,label.x,side=1,line=5.0,cex=0.8)
  }
  
  if (i==1 || i==6 || i==11 || i==16 ) {
    axis(2,yaxp=yrange,lwd=2.0)
#    mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
  }
}

label.x <- expression(paste("time step (fs)"))
mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)        
label.y <- expression(paste("energy (kcal/mol)"))
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
