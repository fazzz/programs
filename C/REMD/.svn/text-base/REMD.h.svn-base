
#define MAXREP 40

struct rep_HREMD_vac_SDFF {
  char *parmtop;
  char rstname[100];
  char dir[100];
  char dirmin[100];
  char dirrex[100];
  char dirsam[100];
};

int HreplicaExchange(int numjudge,int numreplica,char *crd,char *vel, char *top,char *clt,char *pn,char *inpmd[MAXREP],char *rstin[MAXREP],char *pepca[MAXREP],char *dirbase,double temp,int waittime );

int  runMD(char *inpmd,char *crd, char *vel,char *top, char *clt,char *rst, char *out, char *pn, char *dir, char *pepca,char *qsubfilename);

int  runMD_Amber(char *inpmd,char *crd,char *ref, char *vel,char *top, char *rst, char *trj, char *out, char *pn, char *dir,char *qsubfilename, char *simtype);

double calcdelta(char *crdnamei,char *crdnamej,char *eignamei,char *eignamej,char *parmname,double temp);

int HREMD_vac_SDFF(int numexchange,int numreplica,char *minin,char *rexin, char *samin, char **crd,char **top,char *pn,char *dirbase,char *indexfilename,int numnb, int num14,double temp,int waittime,int numstep,int numatom);

double CD_SDFF(char *crdnamei,char *crdnamej,char *parmname1,char *parmname2,char *indexfilename,int numnb, int num14,double temp);

void gatherdata(int numexchange,int numreplica,int numstep,int numatom,char *dirbase,char *pn);








