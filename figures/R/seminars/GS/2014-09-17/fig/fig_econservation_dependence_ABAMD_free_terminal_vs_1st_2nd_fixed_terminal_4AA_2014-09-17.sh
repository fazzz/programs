#!~/bin/sh

parafilename=~/seminars/GS/2014-09-17/para/para_sander_26_dipeptides_2014-09-14.sh

source ${parafilename}

dir=~/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged

dirOUT=~/seminars/GS/2014-09-17
direconbase=econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning

basename=fig_econservation_dependence_ABAMD_free_terminal_vs_1st_2nd_fixed_terminal_4AA_2014-09-17

paraRfil=${dirOUT}/graph/${basename}.R
figfilename=${dirOUT}/eps/${basename}.eps

proname=( )
direcon_free=( )
direcon_fixed1=( )
direcon_fixed2=( )

NRESM=`expr ${NRES} - 1`
k=1
for  ires in `seq 1 ${NRESM}`; do
    i=`expr ${ires} - 1`
    proname[$i]=\"${RES[$ires]}Dv\", 
    direcon_free[$i]=\"${dir}/${RES[$ires]}/${direconbase}_${clustname[$k]}\", 
    k=`expr ${k} + 1 `
    direcon_fixed1[$i]=\"${dir}/${RES[$ires]}/${direconbase}_${clustname[$k]}\", 
    k=`expr ${k} + 1 `
    if [ ${clustnum[$ires]} -eq 3 ]; then
	direcon_fixed2[$i]=\"${dir}/${RES[$ires]}/${direconbase}_${clustname[$k]}\", 
	k=`expr ${k} + 1 `
    else
	direcon_fixed2[$i]=\"XX\", 
    fi
done

i=`expr ${i} + 1`
ires=${NRES}

proname[$i]=\"${RES[$ires]}Dv\"
direcon_free[$i]=\"${dir}/${RES[$ires]}/${direconbase}_${clustname[$k]}\"
k=`expr ${k} + 1 `
direcon_fixed1[$i]=\"${dir}/${RES[$ires]}/${direconbase}_${clustname[$k]}\"
k=`expr ${k} + 1 `
if [ ${clustnum[$ires]} -eq 3 ]; then
    direcon_fixed2[$i]=\"${dir}/${RES[$ires]}/${direconbase}_${clustname[$k]}\"
    k=`expr ${k} + 1 `
else
    direcon_fixed2[$i]=\"XX\"
fi

echo ${proname[*]}

omixl=0.75
omiyl=0.80
omixs=0.2
omiys=0.2
hen=1.5

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

 xrange <- c(0,14,14)
 xrange.axis <- c(0,14,14)
 yrange <- c(-14,0,14) #c(-8,0,8)
 yrange.axis <- yrange

 num.res=${NRES}
 numini=${numini}
 numsim=${numsim}

 direconfree <- c(	${direcon_free[*]}  )
 direconfixed <- c(	${direcon_fixed1[*]}  )
 direconfixed2 <- c(	${direcon_fixed2[*]}  )
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
 iro    <-c( "red", "black", "green" )
 file.names<-NULL

 for (i in 1:num.res) {
  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=paste(direconfree[i],"/",proname[i],"_econ_${temp0}_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(direconfixed[i],"/",proname[i],"_econ_${temp0}_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[3]=paste(direconfixed2[i],"/",proname[i],"_econ_${temp0}_100ps_1-",numini,".econ.rmsd.av",sep='')

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
  text(8,-11,proname[i],cex=1.75)
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")
  lines(c(7,7),c(-10,0),lwd=2,lty="dashed",col="gray")
  lines(c(8,8),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(10,10),c(-10,0),lwd=2,lty="dashed",col="gray")

  if ( i > ${nxxnymo} ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if ( i %% ${nx} == 1 || i == 1) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}

mtext(outer=T,label.x,side=1,line=4.0,cex=1.0)
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)

EOF

Rscript ${paraRfil}
