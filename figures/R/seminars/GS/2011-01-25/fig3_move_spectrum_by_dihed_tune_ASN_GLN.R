source("ini.R")
source("set_spe.R")
source("set_res.R")
source("plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 1

name.res <- c("ASN", "GLN")
num.res <- length(name.res)

#id.ys.gen <- c(8,9,7,14,9,11,7,10,9,12,10,12,10,9,5,9,10,10,9,6)

id.ys.gen <- c(7)
temp <- "300"
id.ys.wogen <- c(3,4,5,6,8,9,10,11)
n.id.ys.wogen <- length(id.ys.wogen)
temp <- "300"

label.x="fraction of V_n"

xrange <- c(0,1.1,11)
xrange.axis <- c(0.1,1.0,9)

yrange <- c(0,700,7)
yrange.axis <- c(0,700,7)

par(mfrow=c(1,1))
par(oma=c(7.5,4.0,2.0,2.0))
par(mar=c(3.0,3.0,3.0,3.0))

dirbase.name <-paste("/home/yamamori/calspa/TAMD/dipeptides/",sep='')

label.names<-NULL
file.names<-NULL

iro[1] <- "red"
iro[2] <- "blue"

for ( i in 1:num.res  ) {
  j<-1
  cat(id.ys)
  sen[i]<-c(3)
  id.ys[i] <- 2
  file.names[i] <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-01-16_pc6/spe_dtrj_",name.res[i],"Dv_T=300_tau=1_1fs_100ps_max_period_sum.txt",sep='')
  label.names[i] <- name.res[i]
}  

plmGeneralwsrange(data.names=file.names,
                  label.names=label.names,
                  id.ys=id.ys,
                  label.size=0.5,axis.size=2.0,
                  iro=iro,axis.ft="F",is.header="F",
                  sdiro=iro,
                  xrange=xrange,yrange=yrange,
                  sdyrange=yrange,
                  is.sen=sen,width=10.0,
                  ,warrow="F")
text(300,0.4,name.res[i])
box(lwd=2.0)

axis(1,xaxp=xrange.axis,lwd=2.0)
mtext(label.x,side=1,line=5.0,cex=0.8)                                                                    
axis(2,yaxp=yrange.axis,lwd=2.0)
mtext(label.y,side=2,line=2.5,cex=0.8)                                                                 

name.out=paste("~/seminars/GS/2011-01-25/fig3_move_spectrum_by_dihed_tune_ASN_GLN",sep="")
OutTiff(name.out,w_val=600,h_val=500)

