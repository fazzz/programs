#!~/bin/sh

nx=6
ny=1

dir=~/calspa/GBRW/Ala-dipeptide_oneD
dirOUT=~/Report/2015-02

parafilename=${dir}/para/para_UMB_vac_ff99SB_36_2014-08-04_2.sh

source ${parafilename}

basename=fig_GBRW_US_Ala-dipeptide_ondD_2015-02-02

paraRfil=${dirOUT}/graph/${basename}.R
figfilename=${dirOUT}/eps/${basename}.eps

dirpmfbase=${dir}/s_UmbSam_vac_2014-10-24_${ff}/
dirpmf=${dirpmfbase}${pname}

pmf=()

k=1
for i in `seq 1 ${numK}`; do
    for j in `seq 1 ${numh}`; do
	pmf[$k]=\"${dirpmf}/pmf_US_vac_K=${K[$i]}-h=${h[$j]}.txt\"
	k=`expr ${k} + 1`
	if [ $i -ne ${numK} -o $j -ne ${numh} ]; then
	    pmf[$k]=,
	    k=`expr ${k} + 1`
	fi
    done
done

#omixl=0.75
##############
# omixl=1.00 #
# omiyl=0.80 #
# omixs=0.2  #
# omiys=0.2  #
# hen=1.5    #
##############

omixl=1.50
omiyl=1.20
omixs=0.3
omiys=0.3
hen=2.25

width=$( echo ${nx} \* ${hen} + ${omixl} + ${omixs} | bc  )
echo ${width}

height=$( echo ${ny} \* ${hen} + ${omiyl} + ${omiys} | bc  )
echo ${height}

nymo=`expr ${ny} - 1`
nxxnymo=`expr ${nx} \* ${nymo}`
nxxnymo=`expr ${nxxnymo} - 1`

cat <<EOF > ${paraRfil}
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

 num.K <- ${numK}
 num.h <- ${numh}

# pmf <- array( c( ${pmf[*]} ) , dim=c(num.K,num.h)  )
 pmf <- array( c( ${pmf[*]} ) , dim=c(num.h,num.K)  )

 file.name=paste("${figfilename}",sep='')
 postscript(file.name,width=${width},height=${height},horizontal=FALSE,onefile=FALSE,paper="special")

 label.x <- "dihedral angle (rad.)"
 label.y <- "pmf"

 par(mfrow=c(${ny},${nx}))
 par(omi=c(${omixl},${omiyl},${omixs},${omixs}))

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

  if ( i > ${nxxnymo} ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  ######################################
  # if ( i %% ${nx} == 1 || i == 1) {  #
  #   axis(2,yaxp=yrange.axis,lwd=2.0) #
  # }				       #
  ######################################
    
}

mtext(outer=T,label.x,side=1,line=6.0,ce=1.0)
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
EOF

Rscript ${paraRfil}
