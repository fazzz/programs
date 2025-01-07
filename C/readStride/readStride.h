
typedef struct Stridedata SD;

struct Stridedata {
  int resnum;
  char resname[3];
  char OL_type_2nd_st;
  char type_2nd_st[10];

  double phi;
  double psi;
  double area;
};

int readStride(FILE *stridefile,SD *SDdata,int numres);



