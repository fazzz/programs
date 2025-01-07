#!/bin/sh

dirOUT=~/seminars/GS/2015-03-18

parafilename=${dirOUT}/fig/para_fig_econ_vdt_100ps_wdi_NVE_ABAMD_AMBER_Termon_comp_PC6-mLF_2015-03-18.sh

nx=9
ny=3

source ${parafilename}

dir=~/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged

direconbase=econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning

basename=fig_econ_vdt_100ps_wdi_NVE_ABAMD_AMBER_Termon_comp_PC6-mLF_2015-03-18

paraRfil=${dirOUT}/graph/${basename}.R
figfilename=${dirOUT}/eps/${basename}.eps

RESname=( )
proname=( )
direconPC6=( )
direconmLF=( )

NRESM=`expr ${NRES} - 1`
k=1
for  ires in `seq 1 ${NRESM}`; do
    i=`expr ${ires} - 1`
    RESname[$i]=\"${RES[$ires]}\", 
    proname[$i]=\"${RES[$ires]}Dv\", 
    direconPC6[$i]=\"${dir}/${RES[$ires]}/${direconbase}${integtypePC6}_${clustname[$k]}\", 
    k=`expr ${k} + 1 `
    direconmLF[$i]=\"${dir}/${RES[$ires]}/${direconbase}${integtypemLF}_${clustname[$k]}\", 
    k=`expr ${k} + 1 `
done

i=`expr ${i} + 1`
ires=${NRES}

RESname[$i]=\"${RES[$ires]}\"
proname[$i]=\"${RES[$ires]}Dv\"
direconPC6[$i]=\"${dir}/${RES[$ires]}/${direconbase}${integtypePC6}_${clustname[$k]}\"
k=`expr ${k} + 1 `
direconmLF[$i]=\"${dir}/${RES[$ires]}/${direconbase}${integtypemLF}_${clustname[$k]}\"
k=`expr ${k} + 1 `

echo ${proname[*]}

omixl=0.75
omiyl=0.80
omixs=0.2
omiys=0.2
hen=1.5
#omixl=1.50
#omiyl=1.60
#omixs=0.4
#omiys=0.4
#hen=3.0

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
 source("~/Rspa/plmGeneral_wsrange_logu.R")

 leg.pos="topright"
 is.leg <- 0

 T <- "600"
 temp <- "300"

 xrange <- c(0,16,16)
# xrange.axis <- c(0,16,16)
# yrange <-  c(-8,0,8) #c(-14,0,14)
 yrange <-  c(-8,2,10) #c(-14,0,14)
# yrange.axis <- yrange

 num.res=${NRES}
 numini=${numini}
 numsim=${numsim}

 direconPC6 <- c(	${direconPC6[*]}  )
 direconmLF <- c(	${direconmLF[*]}  )
 RESname <- c(	${RESname[*]}  )
 proname <- c(	${proname[*]}  )

 file.name=paste("${figfilename}",sep='')
 postscript(file.name,width=${width},height=${height},horizontal=FALSE,onefile=FALSE,paper="special")

 label.x <- "dt (fs)"
 label.y <- "log(dE)"

 par(mfrow=c(${ny},${nx}))
 par(omi=c(${omixl},${omiyl},${omixs},${omixs}))

 hutosa <- c(2, 2, 2)
 sen <- c(3, 3, 3)
 id.ys  <- c(2, 2, 2)
 ids.ys <- c(3, 3, 3)
 tenshu <- c(20, 20, 20)
 senshu <- c(1, 1, 1)
 iro    <-c( "red", "blue", "green" )
 file.names<-NULL

 for (i in 1:num.res) {
  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=paste(direconPC6[i],"/",proname[i],"_econ_${temp0}_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(direconmLF[i],"/",proname[i],"_econ_${temp0}_100ps_1-",numini,".econ.rmsd.av",sep='')

  plmGeneralwsrangelogu(data.names=file.names,
                       sd.names=file.names,
                       id.ys=id.ys,
                       ids.ys=ids.ys,
                       label.size=0.5,axis.size=2.0,
                       iro=iro,axis.ft="F",is.header="T",
                       sdiro=iro,
                       xrange=xrange,yrange=yrange,
                       sdyrange=yrange,
                       is.sen=sen,width=10.0,
                       warrow="T")
  text(8,-7,RESname[i],cex=1.75)
  box(lwd=2.0)

#  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")
#  lines(c(7,7),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(8,8),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(10,10),c(-10,0),lwd=2,lty="dashed",col="gray")

  if ( i == 19 ) {
    xrange.axis <- c(0,15,15)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  else if ( i > 18 && i < 27 ) {
    xrange.axis <- c(1,16,15)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  else if ( i > 18 ) {
    xrange.axis <- c(1,16,15)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if ( i == 1) {
    yrange.axis <- c(-8,2,10)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  else if ( i == 10) {
    yrange.axis <- c(-8,1,9)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  else if ( i == 19) {
    yrange.axis <- c(-8,1,9)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}

#mtext(outer=T,label.x,side=1,line=4.0,cex=1.0)
#mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)

EOF

Rscript ${paraRfil}
