#ifndef INCLUDE_DCA_M
#define INCLUDE_DCA_M

double DCAm_cW(double *W,double *VS,double *V);
double DCAm_cGamma(double *gamma, double *beta, double *W,double *V, double Q, double *VS);
double DCAm_cV(double *V, double *PA, double *PB);
double DCAm_cBeta(double *beta, double *b2A, double *b1B, double *Coacc);
double DCAm_cABI(double *PC1, double *PC2,
		 double *PC12,double *PC21,
		 double *bC1, double *bC2,
		 double *PA1, double *PA2,
		 double *PA12,double *PA21,
		 double *bA1, double *bA2,
		 double *PB1, double *PB2, 
		 double *PB12,double *PB21,
		 double *bB1, double *bB2,
		 double *W, double *gamma);
double DCAm_mainpass(double *PC1, double *PC2, 
		     double *PC12,double *PC21,
		     double *bC1, double *bC2,
		     double *PA1, double *PA2,
		     double *PA12,double *PA21,
		     double *bA1, double *bA2,
		     double *PB1, double *PB2,
		     double *PB12,double *PB21,
		     double *bB1, double *bB2,
		     double *Coacc, double Q);

#endif
