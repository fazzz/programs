source("ini.R")
source("set_spe.R")
source("set_res.R")
source("plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

xrange <- c(0,600,10)
xrange.axis <- c(0,500,5)

yrange <- c(0,0.5,1)
yrange.axis <- c(-1,1,1)

name.res <- c("ASN","GLN","ASP","ALA")
num.res <- length(name.res)
is.ys.gen <- c(7,8,4,3)

temp <- "300"
num.temp <- length(temp)
num.dt <- length(dt)
tau <- "1"
dt <- "1"

par(mfrow=c(2,1))
par(mar=c(5.0,0.0,0.0,0.0))
par(oma=c(7.5,4.0,2.0,2.0))

dirbase.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase2.name <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')

label.names<-NULL
file.names<-NULL

for ( k in 1:2  ) {
  sen[k]<-c(1)
  n<-1+k
}
iro[1]<-"red"
iro[2]<-"blue"

file.names <- NULL
file.name.max.period <- NULL
id.ys <- NULL

k<-1
for ( i in 1:2  ) {
  file.names <- NULL
  file.name.max.period <- NULL
  id.ys <- NULL
  j <- i
  for ( j in 1:2  ) {
    file.names[j] <- paste(dirbase.name,"/",name.res[k],"/",dirbase2.name,"/n=1.0/anl/spe_dtrj_",name.res[k],"Dv_T=",temp,"_tau=1_1fs_100ps.txt",sep='')
    file.name.max.period[j] <- paste(dirbase.name,"/",name.res[k],"/",dirbase2.name,"/n=1.0/anl/spe_dtrj_",name.res[k],"Dv_T=",temp,"_tau=1_1fs_100ps_max_period.txt",sep='')

    b <- read.table(file.name.max.period[j])
    id.ys[j] <- b[,1]
    k <- k+1
  }

#  cat(file.names,"\n")
  cat(id.ys,"\n")
  plmGeneralwsrange(data.names=file.names,             ##
                    id.ys=id.ys,                       ##
                    label.size=0.5,axis.size=2.0,      ##
                    iro=iro,axis.ft="F",is.header="T", ##
                    sdiro=iro,                         ##
                    xrange=xrange,yrange=yrange,       ##
                    sdyrange=yrange,                   ##
                    is.sen=sen,width=10.0,             ##
                    ,warrow="F")                       ##
  box(lwd=2.0)                                         ##

  axis(1,xaxp=xrange.axis,lwd=2.0)
  
  axis(2,yaxp=yrange.axis,lwd=2.0)
}

name.out=paste("~/seminars/GS/2012-02-29/fig2.spectra_faster_dihed_T=300K_4dipep_2012-02-29",sep="")
OutTiff(name.out,w_val=200,h_val=500)
