/////////////////////////////////////////
//                                     //
//           /MD_NA/src/BD.h           //
//                                     //
/////////////////////////////////////////

// ’è”
int DYNMODE;

// ŠÖ”
// –€CŒW”‚Ìİ’è‚ğs‚¤ŠÖ”ŒS
// „‘Ì‚Ì–€CŒW”‚Ìİ’è‚ğs‚¤ŠÖ”
void set_friction_tensor(void);
// –€CŒW”‚Ìİ’è‚Ì•â•‚ğs‚¤ŠÖ”
void sub_set_friction_tensor(int nNumClut);
// „‘Ì‚Ì‹…‹ß—‚ğs‚¤ŠÖ”
double pick_clust_radius(int nNumClut);
// „‘Ì‚Ì‹…‹ß—‚Ì•â•‚ğs‚¤ŠÖ”_1
double sub_pick_clust_radius_1(int nNumClut);
// „‘Ì‚Ì‹…‹ß—‚Ì•â•‚ğs‚¤ŠÖ”_2
double sub_pick_clust_radius_2(int nNumClut);
// „‘Ì‚Ì‹…‹ß—‚Ì•â•‚ğs‚¤ŠÖ”_3
double sub_pick_clust_radius_3(int nNumClut);
// ’†S‚Ì’Tõ‚ğs‚¤ŠÖ”
void mc_move(double origin[3], double step_limit);
// ”÷¬‚È•ÏˆÊ‚Ì”­¶‚ğs‚¤ŠÖ”
double randum_delta_q(double step_limit);

// —h“®‚ÌŒvZ‚ğs‚¤ŠÖ”ŒS
// „‘Ì‚Ì—h“®‚ÌŒvZ‚ğs‚¤ŠÖ”
void Calc_Brownian(void);
// ŠgUŒW”‚Ìİ’è‚ğs‚¤ŠÖ”
void set_diffusion_tensor(int nNumClut);
void Calc_Brownian_force(int nNumClut);

void calc_d_theta_cycle(void);
void calc_sp_velo_cycle(int nNumClut);
void sub_calc_sp_velo_cycle(int nNumClut);

