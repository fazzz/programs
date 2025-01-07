#!~/bin/sh

opt=(dummy height width1 tiffheight tiffwidth nb )
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

proname=AD
dir=~/calspa/refcalc/UmbSam/AD
dirout=~/Report/2012-11

mode=2
width=0.3
filenamepmf=( dummy )

#parafile=( dummy ~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=5_2012-11-12.sh ~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh ~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=20_2012-11-12.sh  )

parafile=( dummy ~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh )

np=`expr ${#parafile[*]} - 1`

n=${np}

n1=`expr 1 + ${n}`

for i in `seq 1 ${np}`; do
    source ${parafile[$i]}

    TLns=`expr ${TLbase} / 1000`
    direqubase=${dir}/s_UmbSam_vac_2012-11-12_${ff}/
    dirpmf=${direqubase}${pname}/pmf
    pmf=${dirpmf}/pmf_UmbMD_vac_${TLns}ns_nbx=${nb}_nby=${nb}.txt

    filenamepmf[$i]=${pmf}
done

paraRfile=${dirout}/graph/para_s_UmbSam_vac_2012-11-14.R

cat <<EOF > ${paraRfile}
dir <- "${dir}"
dirout <- "${dirout}"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dirout,"/tiff/","fig_s_UmbSam_vac_2012-11-14",sep='')
   
level <- seq(0.0,${height},${width1})

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=${tiffwidth},height=${tiffheight})

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c($(echo -n 1)$(for i in `seq 2 ${n1}`; do echo -n ,${i} ;done)),1,${n1},byrow=TRUE),c($(for i in `seq 1 ${n}`; do echo -n "3," ;done)1),c(1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

EOF
for i in `seq 1 ${n}`; do 
    echo "name.pmf <- paste(\"${filenamepmf[$i]}\",sep='')" >> ${paraRfile}
    echo "felwFillConMap(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt=\"F/kBT\")" >> ${paraRfile}
done 

cat <<EOF >> ${paraRfile}
felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
EOF

Rscript ${paraRfile}; echo ${paraRfile}


