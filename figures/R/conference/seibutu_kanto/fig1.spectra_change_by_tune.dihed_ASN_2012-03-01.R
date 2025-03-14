
source("ini.R")
source("set_spe.R")
source("set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

name.res <- c( "ASN" )
num.res <- length(name.res)

temp <- "30"

id.ys.gen <- c(7)
id.ys.wogen <- c(3,4,5,6,8,9,10,11,7)
n.id.ys.wogen <- length(id.ys.wogen)

para.dihed.V <- "0.3"
num.para.dihed.V <- length(para.dihed.V)

xrange <- c(0,1000,10)
xrange.axis <- c(0,1000,10)

yrange <- c(-0.8,0.8,5)
yrange.axis <- c(-10,10,5)

num.para.dihed.V <- length(para.dihed.V)

par(mfrow=c(1,2))
#par(oma=c(7.5,4.0,6.0,2.0))
par(oma=c(3.0,1.0,2.0,0.5))
par(mar=c(0.0,0.0,0.0,0.0))

dirbase.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase2.name <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')

label.names<-NULL
file.names<-NULL
title<-NULL

k <- 1
j<-1

cat(i,name.res[1],"\n")

file.name <- paste(dirbase.name,"/",name.res[1],"/",dirbase2.name,"/n=",para.dihed.V,"/anl/spe_diff_dtrj_",name.res[1],"Dv_T=",temp,"_tau=1_1fs_100ps.txt",sep='')

id.xs<-NULL
id.ys<-NULL
ave.x<-NULL
sen<-NULL
file.names <- NULL

for ( j in 1:n.id.ys.wogen ) {
  sen[j]<-c(1)                  
  n<-1+j                        
  iro[j] <- "gray"
  if ( j==n.id.ys.wogen )
    iro[j] <- "red"
  
  id.xs[j] <- 1                 
  id.ys[j] <- id.ys.wogen[j]
  ave.x[j] <- 1                 
  file.names[j] <- file.name    
}

cat(iro,"\n")

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
label <- paste(name.res[1],"_n=0.1")
box(lwd=2.0)

axis(1,xaxp=xrange.axis,lwd=2.0)

name.res <- c( "GLN" )
num.res <- length(name.res)

temp <- "300"

id.ys.gen <- c(8)
id.ys.wogen <- c(3,4,5,6,7,9,10,11,12,8)
n.id.ys.wogen <- length(id.ys.wogen)

para.dihed.V <- "0.3"
num.para.dihed.V <- length(para.dihed.V)

xrange <- c(0,1000,10)
xrange.axis <- c(0,1000,10)

yrange <- c(-0.8,0.8,5)
yrange.axis <- c(-10,10,5)

num.para.dihed.V <- length(para.dihed.V)

dirbase.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase2.name <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')

label.names<-NULL
file.names<-NULL
title<-NULL

k <- 1
j<-1

cat(i,name.res[1],"\n")

file.name <- paste(dirbase.name,"/",name.res[1],"/",dirbase2.name,"/n=",para.dihed.V,"/anl/spe_diff_dtrj_",name.res[1],"Dv_T=",temp,"_tau=1_1fs_100ps.txt",sep='')

id.xs<-NULL
id.ys<-NULL
ave.x<-NULL
sen<-NULL
file.names <- NULL

for ( j in 1:n.id.ys.wogen ) {
  sen[j]<-c(1)                  
  n<-1+j                        
  iro[j] <- "gray"
  if ( j==n.id.ys.wogen )
    iro[j] <- "red"
  
  id.xs[j] <- 1                 
  id.ys[j] <- id.ys.wogen[j]
  ave.x[j] <- 1                 
  file.names[j] <- file.name    
}

cat(iro,"\n")

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
label <- paste(name.res[1],"_n=0.1")
box(lwd=2.0)

axis(1,xaxp=xrange.axis,lwd=2.0)

name.out=paste("~/gakkai/seibutu_kanto/fig1.spectra_change_by_tune.dihed_ASN_2012-03-01",sep="")
OutTiff(name.out,w_val=450,h_val=250)
