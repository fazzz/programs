 source("~/Rspa/ini.R")
 source("~/Rspa/set_enecon.sh")
 source("~/Rspa/set_res.R")
 source("~/Rspa/plmGeneral_wsrange.R")

 leg.pos="topright"
 is.leg <- 0

# xrange <- c(-pi,pi,6)
 xrange <- c(-4,4,8)
 xrange.axis <- c(-4,3,7)
 yrange <- c(-16,160,17)
 yrange.axis <- yrange

 num.K <- 6
 num.h <- 4

# pmf <- array( c( "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=10-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=10-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=10-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=10-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=20-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=20-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=20-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=20-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=30-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=30-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=30-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=30-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=40-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=40-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=40-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=40-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=50-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=50-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=50-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=50-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=60-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=60-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=60-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=60-h=0.05.txt" ) , dim=c(num.K,num.h)  )
 pmf <- array( c( "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=10-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=10-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=10-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=10-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=20-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=20-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=20-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=20-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=30-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=30-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=30-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=30-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=40-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=40-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=40-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=40-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=50-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=50-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=50-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=50-h=0.05.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=60-h=0.001.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=60-h=0.005.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=60-h=0.01.txt" , "/home/yamamori/calspa/GBRW/Ala-dipeptide_oneD/s_UmbSam_vac_2014-10-24_ff99SB/Umb36_alPhi_K=40/pmf_US_vac_K=60-h=0.05.txt" ) , dim=c(num.h,num.K)  )

 file.name=paste("/home/yamamori/Report/2015-02/eps/fig_GBRW_US_Ala-dipeptide_ondD_2015-02-02.eps",sep='')
 postscript(file.name,width=15.30,height=3.75,horizontal=FALSE,onefile=FALSE,paper="special")

 label.x <- "dihedral angle (rad.)"
 label.y <- "pmf"

 par(mfrow=c(1,6))
 par(omi=c(1.50,1.20,0.3,0.3))

 for (i in 1:num.K) {
  par(mar=c(0.0,0.0,0.0,0.0))

  file.names<-NULL
  hutosa <- NULL
  sen <- NULL
  id.ys  <- NULL
  tenshu <- NULL
  senshu <- NULL

  for (j in 1:num.h) {
#    file.names[j]=paste(pmf[i,j],sep='')
    file.names[j]=paste(pmf[j,i],sep='')
    hutosa[j] <- 1
    sen[j] <- 1
    id.ys[j]  <- 3
    tenshu[j] <- 1
    senshu[j] <- 1
  }

  plmGeneralwsrange(data.names=file.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    xrange=xrange,yrange=yrange,
                    is.sen=sen)

  txt <- paste("# of base =",i,0,spe="")
  text(0,145,txt,cex=1.00)
  box(lwd=2.0)

  if ( i > -1 ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  ######################################
  # if ( i %% 6 == 1 || i == 1) {  #
  #   axis(2,yaxp=yrange.axis,lwd=2.0) #
  # }				       #
  ######################################
    
}

mtext(outer=T,label.x,side=1,line=6.0,ce=1.0)
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
