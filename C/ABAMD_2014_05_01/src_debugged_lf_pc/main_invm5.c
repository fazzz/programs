#include <stdio.h>
#include <stdlib.h>

void invm5(double *a, double *inv_a, int n);

int main(int argc, char *argv[])
{
  int i,j;
  double *a;
  double *inv_a; 

  a=(double *)malloc(sizeof(double)*4*4);
  inv_a=(double *)malloc(sizeof(double)*4*4);

  a[0*4+0]=1;   a[0*4+1]=2;   a[0*4+2]=0; a[0*4+3]=-1;

  a[1*4+0]=-1;  a[1*4+1]=1;  a[1*4+2]=2;  a[1*4+3]=0;

  a[2*4+0]=2;   a[2*4+1]=0;  a[2*4+2]=1;  a[2*4+3]=1;

  a[3*4+0]=1;   a[3*4+1]=-2; a[3*4+2]=-1;  a[3*4+3]=1;

  invm5(a, inv_a, 4);

  for(i=0;i<4;i++){             
    for(j=0;j<4;j++){           
      printf(" %f",inv_a[i*4+j]);
    }
    printf("\n");
  }                          

  return 0;
}

void invm5(double *a, double *inv_a, int n) {
  //  double a[4][4]={{1,2,0,-1},{-1,1,2,0},{2,0,1,1},{1,-2,-1,1}};
  //  double inv_a[4][4]; 
  double buf; 
  int i,j,k; 
  //  int n=4;  

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      inv_a[i*n+j]=(i==j)?1.0:0.0;
    }
  }

  for(i=0;i<n;i++){
    buf=1/a[i*n+i];
    for(j=0;j<n;j++){
      a[i*n+j]*=buf;
      inv_a[i*n+j]*=buf;
    }
    for(j=0;j<n;j++){
      if(i!=j){
	buf=a[j*n+i];
	for(k=0;k<n;k++){
	  a[j*n+k]-=a[i*n+k]*buf;
	  inv_a[j*n+k]-=inv_a[i*n+k]*buf;
	}
      }
    }
  }

/**********************************/
/* for(i=0;i<n;i++){		  */
/*   for(j=0;j<n;j++){		  */
/*     printf(" %f",inv_a[i][j]); */
/*   }				  */
/*   printf("\n");		  */
/**********************************/
}

/* 
 2.000000 2.000000 -1.000000 3.000000
 -4.000000 -5.000000 3.000000 -7.000000
 3.000000 4.000000 -2.000000 5.000000
 -7.000000 -8.000000 5.000000 -11.000000
*/
