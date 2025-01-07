////////////////////////////////
//                            //
// mas_MD/src/gen_rand_num.h  //
//                            //
////////////////////////////////

// �萔
#define M32(x) (((1UL << 32) - 1) & (x))

// �ϐ�
static unsigned long seed = 1;
static int jrnd;
static unsigned long x[521];

// �֐�
// ��l�����̐��������̕��@�ōs���֐��S
// ���`�����@��p���Ĉ�l�����̐������s���֐�
void init_rnd_lcm(unsigned long x);
unsigned long irnd_lcm(void);
double rnd_lcm(void);

// M �n�@��p���Ĉ�l�����̐������s���֐�
void rnd521_Ms(void);
void init_rnd_Ms(unsigned long seed);
unsigned long irnd_Ms(void);
double rnd_Ms(void);

// rand() ��p����
// 0 �ȏ� 1 �ȉ��̎�����
// ��l�����̐������s���֐�
double rnd(void);

// ���K���z�̐������s���֐�
double nrmd_Box_Muller_method(int num, 
                              double ave, double vari);

double rnd_n(int n);

