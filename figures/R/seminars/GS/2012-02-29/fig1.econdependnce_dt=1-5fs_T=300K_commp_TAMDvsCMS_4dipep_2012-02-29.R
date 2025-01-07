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

yrange <- c(-7,0,7)

name.res <- c( "ASN","GLN","ASP","GLY")

numini=10

num.res <- length(name.res)

dirbase.name <- "/home/yamamori/calspa/TAMD/dipeptides/"

for ( k in 1:4  ) {
  sen[k]<-c(1)
  n<-1+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

par(mfrow=c(2,3))

par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
k2<-1
for (i in 1:num.res) {
  par(mar=c(7.0,0.0,0.0,0.0))

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

  a <- read.table(file.names[3],header=TRUE)
  cat(a[2,2])
  lines(c(1,6),c(a[2,2],a[2,2]),lwd=2,lty="dashed",col="gray")

  axis(1,xaxp=xrange.axis,lwd=2.0)
  
  if (i==1 || i==3 ) {
    axis(2,yaxp=yrange,lwd=2.0)
  }

  if (i==2 || i==4) {
    par(mar=c(7.0,10.0,0.0,0.0))

    dirbase.name2 <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
    dirbase2.name2 <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')

    xrange2 <- c(0,600,10)
    xrange.axis2 <- c(0,500,5)

    yrange2 <- c(0,0.5,1)
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
      file.names2[j2] <- paste(dirbase.name2,"/",name.res[k2],"/",dirbase2.name2,"/n=1.0/anl/spe_dtrj_",name.res[k2],"Dv_T=",temp,"_tau=1_1fs_100ps.txt",sep='')
      file.name.max.period2[j2] <- paste(dirbase.name2,"/",name.res[k2],"/",dirbase2.name2,"/n=1.0/anl/spe_dtrj_",name.res[k2],"Dv_T=",temp,"_tau=1_1fs_100ps_max_period.txt",sep='')

      b2 <- read.table(file.name.max.period2[j2])
      id.ys2[j2] <- b2[,1]
      k2 <- k2+1
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
    
    axis(1,xaxp=xrange.axis2,lwd=2.0)
    axis(2,yaxp=yrange.axis2,lwd=2.0)

  }
}


name.out=paste("~/seminars/GS/2012-02-29/fig1.econdependnce_dt=1-5fs_T=300K_commp_TAMDvsCMS_4dipep_2012-02-29",sep='')
OutTiff(name.out,w_val=850,h_val=500)
