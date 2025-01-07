#!~/bin/sh

K=( dummy 3 4 5 6 )
nK=5
mode=2
tiffwidth=420
tiffheight=200
height=20
width=4

proname=AD

#source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=350_fq=10ps_99SB_KZmax=1000.sh
source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=750_fq=10ps_99SB.sh

dir=~/calspa/MFEP/AD

filename=( dummy )

filename[1]=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_0.3_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21

source ~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh

filename[2]=~/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_${ff}/${pname}/nEX=${numEX}/freq=${TLbase}ps/pmf/pmf_ADv

parafile=~/gakkai/MSSJ_2012/graph/para_fig_abst_2DFEL_2012-10-22.R
cat <<EOF > ${parafile}
level <- seq(0,10,1)

name.title <- NULL

MyColor <- function(n,alpha=1)
{			
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%3		
        k <- n%/%3		
        i <- n - j - k		
        c(if (i > 0) hsv(h = seq.int(from = 40/60, to = 25/60,
            length.out = i), alpha = alpha), if (j > 0) hsv(h = seq.int(from = 23/60,
            to = 11/60, length.out = j), alpha = alpha), if (k > 
            0) hsv(h = seq.int(from = 8/60, to = 0/60, length.out = k-1),
            alpha = alpha, s = seq.int(from = 1, to = 0.9, length.out = k-1),
            v = 1),hsv(0,0,1))
    }			
    else character(0L)	
}			

dir <- "~/gakkai/MSSJ_2012"

title <- NULL

name <- paste("${filename}",sep='')

level <- seq(0.0,${height},${width})

name.out <- paste(dir,"/tiff/","fig_abst_2DFEL_2012-10-22",sep='')
file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=${tiffwidth},height=${tiffheight})

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/gakkai/MSSJ_2012/graph/fel_wFillConMap.R")
source("~/gakkai/MSSJ_2012/graph/fel_wFillConBar.R")
						
nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(3,3,1),c(1))

#par(mar =c(5,5,2,2) )
par(mar =c(4,4,2,2) )
#par(mar =c(3,3,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)
#par(cex.lab=1.0)
label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

felwFillConMap("${filename[2]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConMap("${filename[1]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar("${filename[1]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
EOF

Rscript ${parafile}
