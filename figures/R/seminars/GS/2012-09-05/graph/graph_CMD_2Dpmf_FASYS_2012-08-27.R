 
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "/home/yamamori/calspa/refcalc/CMD/FiveAtomSys"
dirbase <- paste(dir0,"/e_CMD_NH_2012-08-20",sep="")

title=name.title

label.x <- expression(paste(theta,1,sep=""))
label.y <- expression(paste(theta,2,sep=""))

name <- paste(dirbase,"/ff=",ffname,"/tau=",tau,"/TL=",TL,"/anl/pmf_FASYSv_T",T,"_",width,sep="")
    
cat(name,'\n')

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=tiffwidth,height=tiffheight)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()

