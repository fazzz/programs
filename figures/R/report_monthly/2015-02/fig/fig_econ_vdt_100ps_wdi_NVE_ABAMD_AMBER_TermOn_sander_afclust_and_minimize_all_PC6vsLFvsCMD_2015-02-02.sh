#!~/bin/sh

nx=9
ny=3

dir=~/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged
dirCMD=~/calspa/refcalc-s-ivy3/refcalc/dipep

parafilename=${dir}/para/para_sander_27_dipeptide_comp_PC6-mLP_2015-01-22_fig.sh

source ${parafilename}

direconbasePC6=econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning${integtypePC6}
direconbasemLF=econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning${integtypemLF}
direconbaseCMD=econ_wdi_NVE_100ps

dirOUT=~/Report/2015-02

basename=fig_econ_vdt_100ps_wdi_NVE_ABAMD_AMBER_TermOn_sander_afclust_and_minimize_all_PC6vsLFvsCMD_2015-02-02

paraRfil=${dirOUT}/graph/${basename}.R
figfilename=${dirOUT}/eps/${basename}.eps

RESIDUE=( )
econPC6=( )
econmLF=( )
econCMD=( )

k=1
for  ires in `seq 1 ${NRES}`; do
    RESIDUE[$k]=\"${RES[$ires]}\"
    econPC6[$k]=\"${dir}/${RES[$ires]}/${direconbasePC6}_${clustname[$k]}/${RES[$ires]}Dv_econ_300_100ps_1-10.econ.rmsd.av\"
    econmLF[$k]=\"${dir}/${RES[$ires]}/${direconbasemLF}_${clustname[$k]}/${RES[$ires]}Dv_econ_300_100ps_1-10.econ.rmsd.av\"
    econCMD[$k]=\"${dirCMD}/${RES[$ires]}/${direconbaseCMD}/${RES[$ires]}Dv_econ_300_100ps_1-10_2ntc.econ.rmsd.av\"

    k=`expr ${k} + 1 `
    if [ $ires -ne ${NRES} ]; then
	RESIDUE[$k]=,
	econPC6[$k]=,
	econmLF[$k]=,
	econCMD[$k]=,
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
 econPC6 <- c(	${econPC6[*]}  )
 econmLF <- c(	${econmLF[*]}  )
 econCMD <- c(	${econCMD[*]}  )

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
 iro    <-c( "red", "blue", "black" )
 file.names<-NULL

 for (i in 1:num.res) {
  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=econPC6[i]
  file.names[2]=econmLF[i]
  file.names[3]=econCMD[i]

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
