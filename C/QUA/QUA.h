
void qua_rot(double q[4],double r[4],double roted[4]);
void qua_product(double r[4],double q[4], double p[4]);
void qua_conjugate(double q[4],double qs[4]);
void qua_trans_omgtodqua(double omg[3],double q[4],double dq[4]);
void trans_omgtodqua(double omg[3],double q[4],double dq[4]);
void qua_setrotmat(double q[4],double rotmat[3][3]);
void qua_euler(double q[4],double euler[3]);
