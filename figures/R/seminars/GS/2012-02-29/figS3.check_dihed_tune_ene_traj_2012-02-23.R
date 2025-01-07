
source("~/Rspa/ini.R")
source("set_res.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

leg.pos="top"
is.leg <- 1

xrange <- c(0,10000,10)
xrange.axis <- c(0,8000,8)

label.x <- expression(paste("time step "))

T    <- c( "30" )

name.res <- c("ALA","ASP","ASN","ARG","CYS","GLN", "GLY","GLU","HIE","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","VAL")
num.res <- length(name.res)

dirbase.name <- "/home/yamamori/calspa/TAMD/mod_parm/dipep/"
dircheck.name <- "/check_dihed_tune_10ps_ABAMD_AMBER_Termon_2012-02-23_pc6/"

par(mfrow=c(4,5))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

file.names<-NULL
file.names<-NULL

yrange <- c(-180,40,9)
yrange.axis <- c(-160,0,8)
label.y <- expression(paste("E"))

for ( i in 1:num.res  ) {
  j<-1
  id.ys[1] <- 2
  id.ys[2] <- 2
  file.names[1] <- paste(dirbase.name,name.res[i],dircheck.name,name.res[i],"Dv_ref_check_T=",T,"_1fs_10ps.out",sep='')
  label.names[1] <- "ref"
  file.names[2] <- paste(dirbase.name,name.res[i],dircheck.name,name.res[i],"Dv_mod_check_T=",T,"_1fs_10ps.out",sep='')
  label.names[2] <- "mod"
  
  plmGeneralwsrange(data.names=file.names,
                    label.names=label.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,width=10.0,
                    ,warrow="F")
  text(5000,0.4,name.res[i])
  box(lwd=2.0)

  if (i==16 || i==17 || i==18 || i==19 || i==20) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
    mtext(outer=T,label.x,side=1,line=5.0,cex=0.8)
  }
  
  if (i==1 || i==6 || i==11 || i==16 ) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
    mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
  }
}

name.out=paste("~/seminars/GS/2012-02-29/figS3.check_dihed_tune_ene_traj_2012-02-23",sep='')
OutTiff(name.out,w_val=1000,h_val=700)




