#!~/bin/sh

opt=(dummy name parabasefile height width1 tiffheight tiffwidth )
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

dir=~/calspa/TACCM_CGAAREMD/AD

proname=AD

source ${parabasefile}

dir=~/calspa/MFEP/AD/

n1=`expr 1 + ${n}`

parafile=~/calspa/TACCM_CGAAREMD/AD/fig/fig_e_CG-FG_NH_2012-10-29_${name}.R
cat <<EOF > ${parafile}
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "${dir}"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("${filenamepath}",sep='')
name.out <- paste(dir,"/tiff/","fig_e_CG-FG_NH_2012-10-29",sep='')
   
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
    echo "name.pmf <- paste(\"${filenamepmf[$i]}\",sep='')" >> ${parafile}
    echo "felwFillConMap(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt=\"F/kBT\")" >> ${parafile}
done 

cat <<EOF >> ${parafile}
felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
EOF

Rscript ${parafile}

display ${dir}//tiff/fig_e_CG-FG_NH_2012-10-29.tiff &
