source("ini.R")
source("set_spe.R")
source("set_res.R")
source("plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

name.res <- c("ASN")
num.res <- length(name.res)

id.ys.gen <- c(7)
temp <- "300"
id.ys.wogen <- c(7)
n.id.ys.wogen <- length(id.ys.wogen)
temp <- "300"

para.dihed.V <- c( "0.2", "0.4", "0.6", "0.8"  )
num.para.dihed.V <- length(para.dihed.V)

xrange <- c(0,1000,10)
xrange.axis <- c(0,900,9)

yrange <- c(0.0,0.2,6)
yrange.axis <- yrange

num.para.dihed.V <- length(para.dihed.V)

par(mfrow=c(1,num.para.dihed.V))
par(oma=c(7.5,4.0,6.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

dirbase.name <-paste("/home/yamamori/calspa/TAMD/dipeptides/",sep='')
dirbase21.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase22.name <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')

label.names<-NULL
file.names<-NULL
title<-NULL

for ( k in 1:num.para.dihed.V  ) {
  for ( i in 1:num.res  ) {
    j<-1
    cat(id.ys)
    sen[2]<-c(1)
    n<-1+j
    iro[1]<- "red"
    id.ys[1] <- 7
    file.names[j] <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-01-16_pc6/n=",para.dihed.V[k],"/anl/spe_dtrj_ASNDv_T=300_tau=1_1fs_100ps.txt",sep='')

    file.names[2] <- paste(dirbase21.name,"/",name.res[i],"/",dirbase22.name,"/n=",para.dihed.V[k],"/anl/spe_dtrj_",name.res[i],"Dv_T=",temp[i],"_tau=1_1fs_100ps.txt",sep='')

    id.ys[2] <- 7
    iro[2]<- "blue"

    senshu[1] <- 1
    senshu[2] <- 2
    
    plmGeneralwsrange(data.names=file.names,
                      label.names=para.dihed.V,
                      id.ys=id.ys,
                      label.size=0.5,axis.size=2.0,
                      iro=iro,axis.ft="F",is.header="F",
                      sdiro=iro,
                      xrange=xrange,yrange=yrange,
                      sdyrange=yrange,
                      is.sen=sen,width=10.0,
                      warrow="F")
    box(lwd=2.0)

    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
}
name.out=paste("~/seminars/GS/2012-02-29/figS1.spectra_chang_dihed.tune_comp_programs_2012-02-29",sep="")
OutTiff(name.out,w_val=600,h_val=350)
