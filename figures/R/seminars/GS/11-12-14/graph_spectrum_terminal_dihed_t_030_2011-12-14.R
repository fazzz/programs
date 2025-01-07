source("ini.R")
source("set_spe.R")
source("set_res.R")
source("plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

xrange <- c(0,600,10)
xrange.axis <- c(0,500,5)

yrange <- c(0,0.5,1)
yrange.axis <- c(-1,1,1)

name.res <- c("SER","ASN","GLN","ILE","LYS","THR", "LEU","GLU","PRO","VAL","GLY","CYS","PHE","TRP", "MET","ALA","HIE","TYR","ASP", "ARG")

num.res <- length(name.res)
#id.ys.gen <- c(5, 14, 6, 9, 6, 11, 10, 7, 9, 5, 7, 9, 10, 6, 9, 6, 10, 9, 5, 6)
id.ys.gen <- c(6,7,11,5,9,6,7,10,5,6,7,6,9,9,10,5,9,10,9,14)
temp <- "30"
num.temp <- length(temp)
num.dt <- length(dt)
tau <- "1"
dt <- "1"

par(mfrow=c(4,5))
par(oma=c(7.5,4.0,2.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

dirbase.name <-paste("/home/yamamori/calspa/TAMD/dipeptides/",sep='')

label.names<-NULL
file.names<-NULL

for ( k in 1:4  ) {
  sen[k]<-c(1)
  n<-1+k
  iro[k]<-c(n)
}

for ( i in 1:num.res  ) {
  j<-1
  id.ys <- id.ys.gen[i]
  cat(id.ys)
  file.names <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_11-11-20_pc6/anl/","spe_dtrj_",name.res[i],"Dv_T=",temp,"_tau=",tau,"_",dt,"fs_100ps.txt",sep='')
  
  plmGeneralwsrange(data.names=file.names,
#                    sd.names=file.names,
                    id.ys=id.ys,
#                    ids.ys=ids.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,width=10.0,
                    ,warrow="F")
  text(300,0.4,name.res[i])
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
name.out=paste("~/seminars/GS/11-12-14/graph_spectrum_terminal_dihed_t_",temp,"_2011-12-14",sep="")
OutTiff(name.out,w_val=600,h_val=500)
