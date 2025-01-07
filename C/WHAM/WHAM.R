
WHAM.multi.T<- function(data.names,is.header="T",id.xs=rep(2,100),Ts=c(300,400),Tobj=300) {

  KB=1.98723e-3
  
  n.data <- length(data.names)
  N <- 1:n.data
  U <- 1:n.data
  f <- 1:n.data
  N.total <- 0
  for ( i in 1:n.data) {
    if (file.access(data.names[i]) != 0){
      print(paste(data.names[i],"does not exist."))
    }
    if (is.header=="T") {
      a <- read.table(data.names[i],header=TRUE)
    }
    else {
      a <- read.table(data.names[i],header=FALSE)
    }
    N[i] <- length(a[,id.xs[i]])
    N.total <- N.total+N[i]
    if (i==1) U <- a[,id.xs[i]]
    else U <- rbind(U,a[,id.xs[i]])
  }    

  cat("N=",N,"\n")
  cat("N.total=",N.total,"\n")

  temp <- N.total*n.data
  C <- NULL
  N.now <- 0
  for ( i in 1:n.data) {
    beta <- 1.0/KB*(1.0/Ts[i]-1.0/Tobj)
    for ( j in 1:n.data) {
      for ( k in 1:N[j]) {
        N.now <- N.now+1
        C[N.now] <- exp(-1.0*beta*U[j,k])
      }
    }
  }
  C <- matrix(C,nrow=N.total)
  
  ntotal <- N.total
  nsim <- n.data

  A <- function(gs) {                   
    n1 <- 0.0                               
    for ( i in 2:nsim ) {                   
      n1 <- n1 + N[i]*gs[i-1]               
    }                                       
    n2 <- 0.0                          
    for ( l in 1:ntotal ) {            
      din <- N[1]*C[l,1]                        
      for ( j in 2:nsim  )             
        din <- din + N[j]*C[l,j]*exp(gs[j-1])
      browser()
      n2 <- n2 + log(1.0/din)          
    }                                  
    browser()
    -n1-n2                             
  }                                    

  grA <- function(gs) {
    dA <- NULL
    for ( i in 2:nsim ) {
      n <- 0.0
      for ( l in 1:ntotal ) {
        din <-N[1]*C[l,1]
        for ( k in 2:nsim  )
          din <- din + N[k]*C[l,k]*exp(gs[k-1])
        n <- n + C[l,i]/din
      }
      dA[i-1] <- N[i]*(exp(gs[i-1])*n-1.0)
    }
    dA
  }
  
  A0 <-  A(rep(0.0,n.data-1))
  cat("A0=",A0,"\n")
  dA0 <- grA(rep(0.0,n.data-1))
  cat("grA0=",dA0,"\n")
  
  result<-optim(rep(1.0,n.data-1),A,grA,method="BFGS")
  result
}
