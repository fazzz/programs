#!~/bin/sh

opt=(dummy nameout name height widthFEL )
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

dir=~/defense/

paraRfil=${dir}/graph/fel_pers.R

cat <<EOF > ${paraRfil}

level <- seq(0,${height},${widthFEL})

fact.x <- 1

fact.y <- 1

name.title <- " "

name.inp <- "${name}"
name.out <- "${dir}/tiff/graph_fel_pers"

par(oma = c(0,0,0,0) )
source("~/defense/fig/fel_persp.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1),c(1))

title=name.title

label.x <- expression(paste("CV1"))
label.y <- expression(paste("CV2"))

xrange <- c(0.0,1.0,10)
yrange <- c(0.0,1.0,10)

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=780,height=640)

fel_persp(name.inp,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,5,1,1),plot.axis="yes")

EOF

Rscript ${paraRfil}; echo ${paraRfil}

display ${dir}/tiff/graph_fel_pers.tiff &
