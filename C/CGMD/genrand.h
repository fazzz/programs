////////////////////////////////
//                            //
// mas_MD/src/gen_rand_num.h  //
//                            //
////////////////////////////////

// 定数
#define M32(x) (((1UL << 32) - 1) & (x))

// 変数
static unsigned long seed = 1;
static int jrnd;
static unsigned long x[521];

// 関数
// 一様乱数の生成を種種の方法で行う関数郡
// 線形合同法を用いて一様乱数の生成を行う関数
void init_rnd_lcm(unsigned long x);
unsigned long irnd_lcm(void);
double rnd_lcm(void);

// M 系法を用いて一様乱数の生成を行う関数
void rnd521_Ms(void);
void init_rnd_Ms(unsigned long seed);
unsigned long irnd_Ms(void);
double rnd_Ms(void);

// rand() を用いて
// 0 以上 1 以下の実数の
// 一様乱数の生成を行う関数
double rnd(void);

// 正規分布の生成を行う関数
double nrmd_Box_Muller_method(int num, 
                              double ave, double vari);

double rnd_n(int n);

