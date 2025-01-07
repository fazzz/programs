
struct data{
  int numstep;
  double dt;
  double *Temp;
  double *enet;
  double *enek;
  double *enep;
  double *enekv;
  double *enepv;
};



int pickMDdata(FILE *mddata, FILE *out);
