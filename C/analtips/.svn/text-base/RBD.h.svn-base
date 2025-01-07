
typedef struct RigidBodyData {
  int num_origin;
  int num_terminal[10];
  int flag_terminal;
  int numatom_body;
  int num_branch;
  int hingmat;
  int num_parent_body;
  int num_child_body[10];
}RigidBody;

typedef struct RigidBodyCollection {
  int numbody;
  int *indexofABAcyc;

  RigidBody *RB;
}RigidBodyColl;

int RBD_setInertia(double *crd,double *crd_org, double *mass,int numatom,double *Inertia);
double RBD_setcom(double *crd,double *mas,int numatom,double *com);
double RBD_removecom(double *crd,double *mas,int numatom);
int RBD_setmassweightedcrd(double *crd,double *mas,int numatom,double *mwcrd);

void RBD_trans_body_coordinate(int natom1,int natom2,int natom3,int numatom,double *coord,double *Rotmat, double *coord_body);
int RBD_raedRBdata(FILE *inputfile2,RigidBodyColl RBC);
