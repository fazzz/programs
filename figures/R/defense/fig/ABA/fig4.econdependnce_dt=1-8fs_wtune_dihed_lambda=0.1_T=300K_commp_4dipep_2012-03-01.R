source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"
temp2 <- "300"

xrange <- c(0,9,9)
yrange <- c(-7,0,7)

lambda <- "0.1"

name.res <- c( "ASN","GLN","LEU","SER")

numini=10

num.res <- length(name.res)

dirbase.name <- "/home/yamamori/calspa/TAMD/dipeptides/"
dirbase21.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase22.name <- paste("econ_wdi_NVE_AMBER_TermOn_wtune_2012-02-29_pc6",sep='')

sen<-NULL
id.ys<-NULL
ids.ys<-NULL
tenshu<-NULL

for ( k in 1:4  ) {
  sen[k]<-c(1)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

#iro[1] <- "magenta"
iro[1] <- "red"
iro[2] <- "red"
iro[3] <- "blue"

par(mfrow=c(2,3))
par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
k2<-1
for (i in 1:num.res) {
  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]<- paste(dirbase21.name,"/",name.res[i],"/",dirbase22.name,"/tune=",lambda,"/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]<- paste(dirbase21.name,"/",name.res[i],"/",dirbase22.name,"/tune=1.0/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')

  file.names[3]=paste("/home/yamamori/calspa/refcalc/dipep/",name.res[i],"/econ_wdi_NVE_100ps/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10_2ntc",".econ.rmsd.av",sep='')

  senshu[1]<-3
  senshu[2]<-1
  senshu[3]<-1
  
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
  lines(c(1,7),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")

  if ( i==3 )   {
    xrange.axis <- c(0,8,8)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  if ( i==4 )   {
    xrange.axis <- c(0,9,9)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  
  if ( i==1 ) {
    yrange.axis <- c(-7,0,7)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  if ( i==3 ) {
    yrange.axis <- c(-7,-1,6)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }

  if (i==2 || i==4) {                                                                                               
#    par(mar=c(7.0,10.0,0.0,0.0))
    par(mar=c(0.0,10.0,0.0,0.0))
    dirbase.name2 <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
    dirbase2.name2 <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')
    xrange2 <- c(0,600,6)
    xrange.axis2 <- c(0,600,6)
    yrange2 <- c(-0.15,0.15,1)
    yrange.axis2 <- c(-1,1,1)

    label.names2<-NULL
    file.names2<-NULL

    sen2 <- NULL
    iro2 <- NULL
    sen2[1] <- 1
    sen2[2] <- 1
    senshu[1]<-1
    senshu[2]<-2
    iro2[1]<-"black"
    iro2[2]<-"black"

    file.names2 <- NULL
    file.name.max.period2 <- NULL
    id.ys2 <- NULL

    i2 <- i/2

    file.names2 <- NULL
    file.name.max.period2 <- NULL
    id.ys2 <- NULL
    j2 <- i
    for ( j2 in 1:2  ) {
      file.names2[j2] <- paste(dirbase.name2,"/",name.res[k2],"/",dirbase2.name2,"/n=0.1/anl/spe_diff_dtrj_",name.res[k2],"Dv_T=",temp2,"_tau=1_1fs_100ps.txt",sep='')
      
      file.name.max.period2[j2] <- paste(dirbase.name2,"/",name.res[k2],"/",dirbase2.name2,"/n=1.0/anl/spe_dtrj_",name.res[k2],"Dv_T=",temp2,"_tau=1_1fs_100ps_max_period.txt",sep='')
#      b2 <- read.table(file.name.max.period2[j2])
#      id.ys2[j2] <- b2[,1]

      k2 <- k2+1
    }

    if (i==4 )
      file.names2[1] <- paste(dirbase.name2,"/",name.res[3],"/",dirbase2.name2,"/n=0.1/anl/spe_diff_dtrj_",name.res[3],"Dv_T=30_tau=1_1fs_100ps.txt",sep='')
    if (i==4 )
      file.names2[2] <- paste(dirbase.name2,"/",name.res[4],"/",dirbase2.name2,"/n=0.1/anl/spe_diff_dtrj_",name.res[4],"Dv_T=30_tau=1_1fs_100ps.txt",sep='')

    if (i==2) {
      id.ys2 <- c(7,8)
    }
    if (i==4) {
      id.ys2 <- c(7,6)
    }

                                        #  cat(file.names,"\n")
    cat(id.ys2,"\n")
    plmGeneralwsrange(data.names=file.names2,
                      id.ys=id.ys2,
                      label.size=0.5,axis.size=2.0,
                      iro=iro2,axis.ft="F",is.header="T",
                      sdiro=iro2,
                      xrange=xrange2,yrange=yrange2,
                      sdyrange=yrange2,
                      is.sen=sen,width=10.0,
                      ,warrow="F")
    box(lwd=2.0)
    if ( i==4 ) {
      axis(1,xaxp=xrange.axis2,lwd=2.0)
    }
    axis(2,yaxp=yrange.axis2,lwd=2.0)
  }
}

name.out=paste("~/gakkai/seibutu_kanto/fig4.econdependnce_dt=1-8fs_wtune_dihed_lambda=0.1_T=300K_commp_4dipep_2012-03-01",sep='')
#OutTiff(name.out,w_val=600,h_val=480)
OutTiff(name.out,w_val=850,h_val=480)
