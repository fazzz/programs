source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

#xrange <- c(0,8,8)
xrange <- c(0,10,10)
xrange.axis <- c(0,7,7)
yrange <- c(-7,0,7)
yrange <- c(-7,0,7)

numini=10

name.res <- c( "ALA", "ARG", "ASN", "ASP", "CYS",  "GLN", "GLU", "GLY", "HIE", "ILE", "LEU", "LYS", "MET", "SER", "TRP", "THR", "TYR", "PHE", "PRO", "VAL" )
num.res <- length(name.res)

dirbase.name <- "/home/yamamori/calspa/TAMD_econ/dipep/"
tune <- c(0.1)
ntune <- length(tune)

for ( k in 1:4  ) {
  sen[k]<-c(1)
  n<-1+k
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

iro[1]<-"red"
iro[2]<-"blue"

par(mfrow=c(4,5))
par(oma=c(7.5,4.0,2.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

file.names<-NULL
k2<-1
for (i in 1:num.res) {
  for (j in 1:ntune) {
    file.names[1]=paste(dirbase.name,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_2012-04-17_pc6/dtune=",tune[j],"/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')
    file.names[2]=paste(dirbase.name,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_afclustwoH_2012-04-17_pc6/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')

                                        #    file.names[2]=paste("/home/yamamori/calspa/refcalc/dipep/",name.res[i],"/econ_wdi_NVE_100ps_afclustwoH_2012-04-19/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10_2ntc",".econ.rmsd.av",sep='')                                        ##    file.names[2]=paste("/home/yamamori/calspa/refcalc/dipep/",name.res[i],"/econ_wdi_NVE_100ps/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10_2ntc",".econ.rmsd.av",sep='')

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
    text(5,-6,name.res[i])
    box(lwd=2.0)

    lines(c(1,9),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")

#    a <- read.table(file.names[3],header=TRUE)
#    cat(a[2,2])

    if (i>=16 ) {
      axis(1,xaxp=xrange.axis,lwd=2.0)
      mtext(outer=T,label.x,side=1,line=5.0,cex=0.8)
    }
  
    if (i==1 || i==6  || i==11 || i==16 ) {
      axis(2,yaxp=yrange,lwd=2.0)
      mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
    }

  }
}

name.out=paste("~/calspa/TAMD_econ/dipep/R/graph_econ_vdt_100ps_wdi_NVE_ABAMD_Ambr_TermOn_dihed_tune_afclustwoH_2012-04-17_2",sep='')
OutTiff(name.out,w_val=850,h_val=500)

