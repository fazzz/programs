/////////////////////////////////////////
//                                     //
//           /MD_NA/src/BD.h           //
//                                     //
/////////////////////////////////////////

// �萔
int DYNMODE;

// �֐�
// ���C�W���̐ݒ���s���֐��S
// ���̖̂��C�W���̐ݒ���s���֐�
void set_friction_tensor(void);
// ���C�W���̐ݒ�̕⏕���s���֐�
void sub_set_friction_tensor(int nNumClut);
// ���̂̋��ߎ����s���֐�
double pick_clust_radius(int nNumClut);
// ���̂̋��ߎ��̕⏕���s���֐�_1
double sub_pick_clust_radius_1(int nNumClut);
// ���̂̋��ߎ��̕⏕���s���֐�_2
double sub_pick_clust_radius_2(int nNumClut);
// ���̂̋��ߎ��̕⏕���s���֐�_3
double sub_pick_clust_radius_3(int nNumClut);
// ���S�̒T�����s���֐�
void mc_move(double origin[3], double step_limit);
// �����ȕψʂ̔������s���֐�
double randum_delta_q(double step_limit);

// �h���̌v�Z���s���֐��S
// ���̗̂h���̌v�Z���s���֐�
void Calc_Brownian(void);
// �g�U�W���̐ݒ���s���֐�
void set_diffusion_tensor(int nNumClut);
void Calc_Brownian_force(int nNumClut);

void calc_d_theta_cycle(void);
void calc_sp_velo_cycle(int nNumClut);
void sub_calc_sp_velo_cycle(int nNumClut);

