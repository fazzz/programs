#!~/bin/sh

opt=(dummy name )
nopt=${#opt[*]}
if [ $# -le `expr ${nopt} - 2`  ]; then
    echo "USAGE " $0  ${opt[*]:1:${nopt}}
    echo $*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

nx=5
ny=4
mindt=1 
maxdt=15 
inc=1 
temp0=300 
numini=10

dt=( dummy ${mindt}  )
i=1
while [ ${dt[$i]} -le ${maxdt}  ]; do
    i=`expr ${i} + 1`
    j=`expr ${i} - 1`
    dt[$i]=`expr ${dt[$j]} + ${inc}`
done
numsim=$i

parafilename=~/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/para/para_sander_fixed_end_ALL_2014-03-12.sh

source ${parafilename}

dir=~/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize
dirfixed=~/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end
dirOUT=~/seminars/GS/2014-04-08/
direconbase=econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning

basename=fig_econdependence_fixed_Term_vs_free_Term_20-dipeptides_2014-04-08
    
paraRfil=${dirOUT}/graph/${basename}.R
figfilename=${dirOUT}/eps/${basename}.eps

proname=( )
direcon=( )
direconfixed=( )

NRESM=`expr ${NRES} - 1`
echo "NRES="${NRES}
echo "NRESM="${NRESM}
for  ires in `seq 1 ${NRESM}`; do
    i=`expr ${ires} - 1`
    proname[$i]=\"${RES[$ires]}Dv\", 
    direcon[$i]=\"${dir}/${RES[$ires]}/${direconbase}${clustname0[$ires]}\", 
    direconfixed[$i]=\"${dirfixed}/${RES[$ires]}/${direconbase}_${clustname[$ires]}\", 
    direconfixed2[$i]=\"${dirfixed}/${RES[$ires]}/${direconbase}_${clustname2[$ires]}\", 
done

i=`expr ${i} + 1`
#ires=`expr ${ires} + 1`
ires=${NRES}

proname[$i]=\"${RES[$ires]}Dv\" 
direcon[$i]=\"${dir}/${RES[$ires]}/${direconbase}${clustname0[$ires]}\" 
direconfixed[$i]=\"${dirfixed}/${RES[$ires]}/${direconbase}_${clustname[$ires]}\" 
direconfixed2[$i]=\"${dirfixed}/${RES[$ires]}/${direconbase}_${clustname2[$ires]}\" 

echo ${proname[*]}

#omixl=0.80
omixl=0.75
#omiyl=0.5
omiyl=0.80
omixs=0.2
#omixs=0.4
omiys=0.2
#hen=2.0
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
 yrange <- c(-8,0,8)
 yrange.axis <- yrange

 num.res=${NRES}
 numini=${numini}
 numsim=${numsim}

 direcon <- c(	${direcon[*]}  )
 direconfixed <- c(	${direconfixed[*]}  )
 direconfixed2 <- c(	${direconfixed2[*]}  )
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

  file.names[1]=paste(direcon[i],"/",proname[i],"_econ_${temp0}_100ps_1-",numini,".econ.rmsd.av",sep='')
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
  text(8,-6,proname[i],cex=1.0)
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")
  lines(c(7,7),c(-10,0),lwd=2,lty="dashed",col="gray")
  lines(c(8,8),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(10,10),c(-10,0),lwd=2,lty="dashed",col="gray")

  if ( i > 15 ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  yrange.axis <- c(-7,0,7)
  if ( i == 1 ) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  yrange.axis <- c(-8,-1,7)
  if (i == 6 || i == 11) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  yrange.axis <- c(-8,-1,7)
  if ( i == 16) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}

mtext(outer=T,label.x,side=1,line=2.0,cex=1.0)
mtext(outer=T,label.y,side=2,line=2.0,cex=1.0)

EOF

Rscript ${paraRfil}
