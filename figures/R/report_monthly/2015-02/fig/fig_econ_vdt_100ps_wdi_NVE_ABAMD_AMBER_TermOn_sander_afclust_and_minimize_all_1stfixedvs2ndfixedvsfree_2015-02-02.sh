#!~/bin/sh

nx=9
ny=3

dir=~/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged

parafilename=~/Report/2015-02/fig/para_sander_27dipeptides_2015-02-02_fig.sh


source ${parafilename}

direconbase1st=econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning
direconbase2nd=econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning
direconbasefree=econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning

dirOUT=~/Report/2015-02

basename=fig_econ_vdt_100ps_wdi_NVE_ABAMD_AMBER_TermOn_sander_afclust_and_minimize_all_1stfixedvs2ndfixedvsfree_2015-02-02

paraRfil=${dirOUT}/graph/${basename}.R
figfilename=${dirOUT}/eps/${basename}.eps

RESIDUE=( )
econ1st=( )
econ2nd=( )
econfree=( )

k=1
for  ires in `seq 1 ${NRES}`; do
    RESIDUE[$k]=\"${RES[$ires]}\"
    econ1st[$k]=\"${dir}/${RES[$ires]}/${direconbase1st}_${clustname1st[$ires]}/${RES[$ires]}Dv_econ_300_100ps_1-10.econ.rmsd.av\"
    econ2nd[$k]=\"${dir}/${RES[$ires]}/${direconbase2nd}_${clustname2nd[$ires]}/${RES[$ires]}Dv_econ_300_100ps_1-10.econ.rmsd.av\"
    econfree[$k]=\"${dir}/${RES[$ires]}/${direconbasefree}_${clustnamefree[$ires]}/${RES[$ires]}Dv_econ_300_100ps_1-10.econ.rmsd.av\"

    k=`expr ${k} + 1 `
    if [ $ires -ne ${NRES} ]; then
	RESIDUE[$k]=,
	econ1st[$k]=,
	econ2nd[$k]=,
	econfree[$k]=,
	k=`expr ${k} + 1`
    fi
done

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

 xrange <- c(0,16,16)
 yrange <- c(-8,2,10)

 num.res=${NRES}

 RES <- c( ${RESIDUE[*]}  )
 econ1st <- c(	${econ1st[*]}  )
 econ2nd <- c(	${econ2nd[*]}  )
 econfree <- c(	${econfree[*]}  )

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

  file.names[1]=econfree[i]
  file.names[2]=econ1st[i]
  file.names[3]=econ2nd[i]

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
  text(8,-7,RES[i],cex=1.75)
  box(lwd=2.0)

  if ( i >= 19 && i < 27 ) {
    xrange.axis <- c(0,15,15)
    axis(1,xaxp=xrange.axis,lwd=2.0,cex=0.75)
  }
  else if ( i == 27 ) {
    xrange.axis <- c(0,16,16)
    axis(1,xaxp=xrange.axis,lwd=2.0,cex=0.75)
  }

  if ( i == 1) {
    yrange.axis=c( -7,2,9 )
    axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.75)
  }
  else if ( i == 10) {
    yrange.axis=c( -7,2,9 )
    axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.75)
  }
  else if ( i == 19) {
    yrange.axis=c( -8,2,10 )
    axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.75)
  }
    
}

mtext(outer=T,label.x,side=1,line=4.0,cex=1.0)
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)

EOF

Rscript ${paraRfil}
