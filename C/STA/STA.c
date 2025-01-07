
double STA_ave(double *data, int numdata) {
  int i,j;

  for (k=0;k<numcolum;++k)
      sum[j*numcolum+k]+=data[k];
  }
  for (i=0;i<numave;++i)
    for (j=0;j<numcolum;++j)
      ave[i*numcolum+j]=sum[i*numcolum+j]/numinterval;


}

double STA_var(double *data, int numdata) {
  for (k=0;k<numcolum;++k)
      sum[j*numcolum+k]+=data[k];
  }
  for (i=0;i<numave;++i)
    for (j=0;j<numcolum;++j)
      ave[i*numcolum+j]=sum[i*numcolum+j]/numinterval;

}
