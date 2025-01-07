source("~/Rspa/ini.R")
source("~/Rspa/set_spe.R")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,600,6)
yrange <- c(0,2.0,2)

name.res <- c("ALA","ASP","ASN","ARG","CYS","GLN", "GLY","GLU","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TYR","TRP","VAL")

num.res <- length(name.res)

name.res.ref <- c()

tp <- c("30")
num.tp <- length(tp)

dirbase.nma <-"/home/yamamori/calspa/TAMD/spectrum/nmadas/"
dirbase.tamd <- "/misc/ivy3/yamamori/ajisai/TAMD/spectrum/tamd_ECEPP/testspace_debug_af_optq/"

sen=c(1,0)
iro<-c(3,2)
id.ys<-c(2,2)
tenshu<-c(1,1)

par(mfrow=c(4,5))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,4.0,2.0,2.0))

for (i in 1:num.res) {
  file.names.nma=paste(dirbase.nma,name.res[i],'_3_13_af_sim_2_4/',name.res[i],'.freq_hist',sep='')
  file.names.tamd=paste(dirbase.tamd,name.res[i],'_spe_step_equ_2_4/',tp[1],'_0.1/spectrum/spe_',name.res[i],'Dv_c_a.txt',sep='')

  file.names=c(file.names.tamd,file.names.nma)

  plmGeneralwsrange(data.names=file.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",
                    xrange=xrange,yrange=yrange,
                    is.sen=sen,width=10.0,is.header="F")
  text(500,1.2,name.res[i])
  if (i==1) {
    legend(300,1.9,legend=c("tamd"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[1],ncol=leg.col)
    legend(300,1.7,legend=c("nmadas"),bty="n",
           cex=1.0,y.intersp=2,
           col=iro[2],ncol=leg.col,pch=c(1))
  }
  box(lwd=2.0)

  if (i==1 || i==6 || i==11 || i==16 ) {
    axis(2,yaxp=yrange,lwd=2.0)
  }

  if (i==16 || i==17 || i==18 || i==19 || i==20) {
    axis(1,xaxp=xrange,lwd=2.0)
  }
  
  if (i==16) {
    mtext(outer=T,label.x,side=1,line=2.5,cex=0.8)
  }
  if (i==1) {
    mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
  }
}
name.out=paste("~/papers/TAMD/fig1.comp_spec_tamd_nmadas_20_aminoacids",sep='')
OutTiff(name.out,w_val=880,h_val=660)
