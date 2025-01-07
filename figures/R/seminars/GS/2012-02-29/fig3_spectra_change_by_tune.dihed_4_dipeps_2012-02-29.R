
source("ini.R")
source("set_spe.R")
source("set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

name.res <- c( "ASN", "GLN", "LYS", "SER" )
num.res <- length(name.res)

id.terminal.dihed.num <- c(1,1,1,1)
id.terminal.dihed <- c(5,6,7,4)

temp <- c("30","300","30","30" )

para.dihed.V <- "0.3"
num.para.dihed.V <- length(para.dihed.V)

xrange <- c(0,1000,10)
xrange.axis <- c(0,1000,10)

yrange <- c(-0.6,0.6,5)
yrange.axis <- yrange

num.para.dihed.V <- length(para.dihed.V)

par(mfrow=c(2,2))
par(oma=c(7.5,4.0,6.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

dirbase.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase2.name <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')

label.names<-NULL
file.names<-NULL
title<-NULL

k <- 1
for ( i in 1:num.res  ) {
  j<-1

  cat(i,name.res[i],"\n")

  file.name <- paste(dirbase.name,"/",name.res[i],"/",dirbase2.name,"/n=",para.dihed.V,"/anl/spe_diff_dtrj_",name.res[i],"Dv_T=",temp[i],"_tau=1_1fs_100ps.txt",sep='')

  file.name.max.period <- paste(dirbase.name,"/",name.res[i],"/",dirbase2.name,"/n=1.0/anl/spe_dtrj_",name.res[i],"Dv_T=",temp[i],"_tau=1_1fs_100ps_max_period.txt",sep='')

  a <- read.table(file.name)
  n.id.ys.wogen <- ncol(a)

  b <- read.table(file.name.max.period)
    

  id.xs<-NULL
  id.ys<-NULL
  ave.x<-NULL
  sen<-NULL
  file.names <- NULL

  for ( j in 1:4 ) {
    sen[j]<-c(1)                  
    n<-1+j                        
    iro[j] <- "gray"

    id.xs[j] <- 1                 
    id.ys[j] <- j+2                 
    ave.x[j] <- 1                 
    file.names[j] <- file.name    
  }

  k <- k+1
  cat("k=",k,"\n")
  cat("id.xs=",id.xs,"\n")
  cat("id.ys=",id.ys,"\n")

  iro[4] <- "red"
  id.ys[4] <- id.terminal.dihed[i]+2
  
  plmGeneralwsrange(data.names=file.names,
                    id.xs=id.xs,
                    id.ys=id.ys,
                    ave.x=ave.x,
                    ids.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="F",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,
                    width=10.0,
                    warrow="F")
  label <- paste(name.res[i],"_n=0.1")
  box(lwd=2.0)
  
  if ( i == 3 || i == 4 ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
}

name.out=paste("~/seminars/GS/2012-02-29/fig3_spectra_change_by_tune.dihed_4_dipeps_2012-02-29",sep="")
OutTiff(name.out,w_val=600,h_val=500)
