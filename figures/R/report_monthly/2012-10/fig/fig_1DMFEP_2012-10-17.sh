#!~/bin/sh

K=17
mode=2
xi=-1.4
yi=1.1
xf=1.1
yf=-0.8
height=6
width=6
tiffheight=200
tiffwidth=300

proname=AD

dir=~/calspa/MFEP/AD

source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=750_fq=10ps_99SB.sh

filename=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/path2_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}_@${xi},${yi}-@${yi},${yf}.txt

paraRfil=~/Report/2012-10/graph/para_fig_1DMFEP_2012-10-17.R

cat <<EOF > ${paraRfil}
source("~/Rspa/plmGeneral_wsrange.R")

name.title <- NULL

name.out <- "~/Report/2012-10//tiff/fig_1DMFEP_2012-10-17"

label.x=" "
label.y="pmf"

xrange <- c(0,1.0,10)
xrange.axis <- xrange
yrange <- c(0,${height},${width})
yrange.axis <- yrange

fact.x <- 1
fact.y <- 1
ave.x <- 1

file.name <- paste(name.out,".tiff",sep='')
cat(file.name,'\n')
tiff(file.name,width=${tiffwidth},height=${tiffheight})

par(mfrow=c(1,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(3.0,6.0,2.0,2.0))

name <- "${filename}"

id.xs <- c(1,1)
id.ys <- c(2,2)
ids.ys <- c(3,3)
iro <- c(1,2)
senshu <- c(1,1)
tenshu <- c(3,3)
hutosa <- c(1,1)
is.leg <- 0

plmGeneralwsrange(data.names=name,
                  sd.names=name,
                  id.ys=id.ys,
                  ids.ys=ids.ys,
                  is.sen=rep(0,20),
                  label.size=0.5,axis.size=2.0,
                  iro=iro,axis.ft="F",is.header="T",
                  sdiro=iro,
                  xrange=xrange,yrange=yrange,
                  sdyrange=yrange,
                  warrow="T")
box(lwd=2.0)

axis(2,yaxp=yrange.axis,lwd=2.0,cex.axis=1.5)
mtext(outer=T,label.y,side=2,line=4.0,cex=1.5)

axis(1,xaxp=xrange.axis,lwd=2.0,cex.axis=1.5)
mtext(outer=T,label.x,side=1,line=4.0,cex=1.5)

EOF

Rscript ${paraRfil}; echo ${paraRfil}


