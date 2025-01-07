#ifndef INCLUDE_ABAb_s_rst
#define INCLUDE_ABAb_s_rst

void ABAbs_restat_read(char *inputvelofilename,int numclut,
		      double *correct,double **correctt_Term,double **correctt_Term2,double correct_s[6]
		      ,int MODE,int TERMMODE);

void ABAbs_restat_write_vel(char *restatvelofilename,int numclut,
			   double *correct,double **correctt_Term,double **correctt_Term2,double correct_s[6]
			   ,int MODE,int TERMMODE);

void ABAbs_restat_read_new(char *inputvelofilename,int numclut,
			  double *correct,double **correctt_Term,double **correctt_Term2,double correct_gzi[5],
			  int MODE,int TERMMODE);

void ABAbs_restat_write_vel_new(char *restatvelofilename,int numclut,
			       double *correct,double **correctt_Term,double **correctt_Term2,double correct_gzi[5],
			       int MODE,int TERMMODE);

void ABAbs_restat_read_new_mvV(char *inputvelofilename,int numclut,
			      double *qvel_b1,double *qvel_b2,
			      double *vel_Term_b1,double *vel_Term_b2,
			      double zeta,
			      int MODE,int TERMMODE);

void ABAbs_restat_write_vel_new_mvV(char *restatvelofilename,int numclut,
				   double *qvel_b1,double *qvel_b2,
				   double *vel_Term_b1,double *vel_Term_b2,
				   double zeta,
				   int MODE,int TERMMODE);

void ABAbs_restat_read_chain(char *inputvelofilename,int numclut,
			    double *correct,double **correctt_Term,double **correctt_Term2,double **correct_gzi,
			    int M,
			    int MODE,int TERMMODE);

void ABAbs_restat_write_vel_chain(char *restatvelofilename,int numclut,
				 double *correct,double **correctt_Term,double **correctt_Term2,double **correct_gzi,
				 int M,
				 int MODE,int TERMMODE);

#endif
