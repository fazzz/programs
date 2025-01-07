source("ini.R")
source("set_spe.R")
source("set_res.R")
source("plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

name.res <- c("ASN")
num.res <- length(name.res)

#id.ys.gen <- c(8,9,7,14,9,11,7,10,9,12,10,12,10,9,5,9,10,10,9,6)

id.ys.gen <- c(7)
temp <- "300"
id.ys.wogen <- c(3,4,5,6,8,9,10,11,7)
n.id.ys.wogen <- length(id.ys.wogen)
temp <- "300"

para.dihed.V <- c( "0.2", "0.6", "0.8"  )
num.para.dihed.V <- length(para.dihed.V)
para.dihed.n <- c( "0.2", "0.6", "0.8"  )
num.para.dihed.n <- length(para.dihed.n)
para.14es <- c( "0.0001", "0.001", "0.01" ) 
num.para.14es <- length(para.14es)
para.14LJ <- c( "0.0001", "0.001", "0.01" )
num.para.14LJ <- length(para.14LJ)
para.14esLJ <- c( "0.0001", "0.001", "0.01" ) 
num.para.14esLJ <- length(para.14esLJ)

xrange <- c(0,1000,10)
xrange.axis <- c(0,1000,10)

yrange <- c(-0.3,0.3,6)
yrange.axis <- yrange

num.para.dihed.V <- length(para.dihed.V)

par(mfrow=c(num.para.dihed.V,5))
par(oma=c(7.5,4.0,6.0,2.0))
par(mar=c(0.0,1.0,0.0,0.0))

dirbase.name <-paste("/home/yamamori/calspa/TAMD/dipeptides/",sep='')

label.names<-NULL
file.names<-NULL
title<-NULL

for ( k in 1:num.para.dihed.V  ) {
  for ( m in 1:5  ) {
    for ( i in 1:num.res  ) {
      j<-1
      cat(id.ys)
      for ( j in 1:n.id.ys.wogen  ) {
        sen[j]<-c(1)
        n<-1+j
        if ( j==9 )
          iro[j] <- "red"
        else
          iro[j]<- "gray"
        id.ys[j] <- id.ys.wogen[j]
        if (m==1) {
          file.names[j] <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-01-16_pc6/n=",para.dihed.V[k],"/anl/spe_diff_dtrj_ASNDv_T=300_tau=1_1fs_100ps.txt",sep='')
          lam <- expression(paste(lambda))
          label=paste(para.dihed.V[k],sep="")
          if ( k==1 )
            title <- "V"
          else 
            title <- NULL
          }
        else if (m==2) {
          file.names[j] <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_wdtune_dihed_term_2012-01-16_pc6/n=",para.dihed.n[k],"/anl/spe_diff_dtrj_ASNDv_T=300_tau=1_1fs_100ps.txt",sep='')
          lam <- expression(paste(lambda))
          label=paste(para.dihed.n[k],sep="")
          if ( k==1 )
            title <- "n"
          else 
            title <- NULL
        }
        else if (m==3) {
          file.names[j] <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_wdtune_14es_term_2012-01-16_pc6/n=",para.14es[k],"/anl/spe_diff_dtrj_ASNDv_T=300_tau=1_1fs_100ps.txt",sep='')
          lam <- expression(paste(lambda))
          label=paste(para.14es[k],sep="")
          if ( k==1 )
            title <- "1-4 es"
          else 
            title <- NULL
        }
        else if (m==4) {
          file.names[j] <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_wdtune_14LJ_term_2012-01-16_pc6/n=",para.14LJ[k],"/anl/spe_diff_dtrj_ASNDv_T=300_tau=1_1fs_100ps.txt",sep='')
          lam <- expression(paste(psi))
          label=paste(para.14LJ[k],sep="")
          if ( k==1 )
            title <- "1-4 LJ"
          else 
            title <- NULL
        }
        else if (m==5) {
          file.names[j] <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_wdtune_14esLJ_term_2012-01-16_pc6/n=",para.14esLJ[k],"/anl/spe_diff_dtrj_ASNDv_T=300_tau=1_1fs_100ps.txt",sep='')
          lam <- expression(paste(psi))
          label=paste(para.14esLJ[k],sep="")
          if ( k==1 )
            title <- "1-4 es & LJ"
          else 
            title <- NULL
        }
      }

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
      text(300,0.2,label)
      box(lwd=2.0)

      if ( k == 1 ) {
        mtext(title,side=3,line=1.0,cex=0.8)
      }
      
      if ( k == num.para.dihed.V ) {
        axis(1,xaxp=xrange.axis,lwd=2.0)
        mtext(outer=T,label.x,side=1,line=5.0,cex=0.8)
      }
      ########################################
      ## if ( m == 1 ) {                    ##
      ##   axis(2,yaxp=yrange.axis,lwd=2.0) ##
      ## }                                  ##
      ########################################
      mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
    }
  }
}
name.out=paste("~/seminars/GS/2011-01-25/fig2_move_spectrum_by_tune_ASN",sep="")
OutTiff(name.out,w_val=900,h_val=500)
