
struct clust {
  int numclt;
  int nainiclt[100];
  int termflag[100];
  int natom[100];
  int nbranch[100];
  int type[100];
  int nafinclt[100][4];
  int parent[100];
  int child[100][4];
  int index[100];
};

int io_dismisstrj_Amberform(FILE *inputfile,int numstep,int numatom);
int io_dismissdata(FILE *inputfile,int numstep,int numdata);
int io_inputtrj_Amberform(FILE *inputfile,double *data);
int io_scanconf(FILE *inputfile,int numatom,double *data,int flag);
int io_scanconfwj(FILE *inputfile,int numatom,double *data,int flag);
int io_scanconf_Amber_rst(FILE *inputfile,double *data);
int io_scanconf_Amber_ini(FILE *inputfile,int numatom, double *data);
int io_scan_aw(FILE *inputfile,int numstep,double *data);
int io_scan_atomtype_traj_aw(FILE *inputfile,int numstep,double *data,int typeflag,char atomtype[20][4], int numtype,int numatomtype);
int io_scan_data(FILE *inputfile,int numdimension,double *data,int flag);
int io_scantimeseries(FILE *inputfile,int numstep,int numdimension,double *data,int flag);
int io_scanodtimeseries(FILE *inputfile,int numcolum,int numspecificcolum, int numintialraw, int numfinalrow, double *data);
int io_scanslidetimeseries(FILE *inputfile,int numstep, int numslide,int numdimension,double *data,int flag);
int io_scanmatrix(FILE *inputfile,int numrow,int numcolum,double *data);
int io_scantraj_aw(FILE *inputfile,int numstep,int numatom,double *traj);
int io_scandtraj(FILE *inputfile,int numstep,int numdihed,double *traj);
int io_scanspdtraj(FILE *inputfile,int numstep,int numdihed,int numspdihed,double *traj);
int io_scanspcalphatraj_aw(FILE *inputfile,int numstep,int numatom,int numres,int numspres,double *traj);
int io_scansptraj_aw(FILE *inputfile,int numstep,int numatom,int numspatom,double *traj);
int io_scnadcoldata(FILE *inputfile,int numf,int numi,int numcol,int xcol,int ycol, double *data);
int io_scantrj(FILE *inputfile,int numatom,int numstep,double **data);
int io_scanconf_atomtype(FILE *inputfile,char *atomname,int numname, int numatom,double *crd);
int io_scancatrj(FILE *inputfile,int numatom,int numstep,double **data);
//int io_scancaconf(FILE *inputfile,int numatom,double *crd);
int io_scantrj2(FILE *inputfile,int numatom,int numstep,double *data);
int io_scancaconf(FILE *inputfile,double *data);
int io_scancoldata(FILE *inputfile,int numf,int numi,int numcol,int col, double *data);
int io_scnascoldata(FILE *inputfile,int numf,int numi,int numcol,int scol, double *data);
int io_scnasdcoldata(FILE *inputfile,int numf,int numi,int numcol,int scol,int scol2, double *data);
double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data);
int io_scancamass(FILE *inputfile,int numatom,double *mass);
int io_chop_timeseries_data(FILE *inputfile,FILE *outputfile,int numlength, int numdata);
int io_chop_trj_data(FILE *inputfile,FILE *outputfile,int numlength, int numatom);

int io_outputconf(FILE *outputfile,int numatom,double *data,int flag);
int io_outputmatrix(FILE *inputfile,int numrow,int numcolum,double *data);
int io_outputdata(FILE *outputfile, int numdimension, double *data);
int io_outputdata_f(FILE *outputfile, int numdimension, double *data);
int io_outputtimeseries(FILE *outputfile,int numstep,int numdimension,double *data);
int io_outputtimeseries_f(FILE *outputfile,int numstep,int numdimension,double *data);
int io_outputt_spdim(FILE *outputfile,int numstep,int numdimension,int spdim,double *data);
int io_outputconf_Amberform(FILE *inputfile,int numatom,double *data);

int io_thin_down_data(FILE *inputfile,FILE *outputfile,int numinterval);

int io_outtrj(FILE *inputfile,int numatom,int numstep,double **data);

void *io_inputtrj_for_gather_Amberform(FILE *inputfile,double *trj,int numstep, int numatom);
int io_addtrj_Amberform(FILE *outputfile,double *trj,int numstep,int numatom);
int io_outtrj2(FILE *outputfile,int numatom,int numstep,double *data);

int io_scanmatrix_ifv(FILE *inputfile,int numcolum,int numrow,int numinirow, int numfinrow,double **data);

int io_scanclust(FILE *clustfile, struct clust *clt);
int io_scanconfwc(FILE *inputfile,int numatom,double *data);
