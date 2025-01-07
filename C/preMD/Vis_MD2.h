/////////////////////////////////////////
//                                     //
//         /Vis_MD/src/DataBase.h      //
//                                     //
/////////////////////////////////////////


// 定数

#define MAXNATOM 100   /* 原子数の上限*/
#define MAXNCLUST 30   /* 自由度の上限*/
#define MAXNDH 30      /* 二面角の上限*/
#define MAXNTD 10       /* 二面種の上限*/
#define MAXNRESIDUE 10  /* 残基数の上限*/
#define MAXNBH 4        /* 分岐数の上限*/
#define NRESIDUE 30     /* 残基種の数 */
#define NHASH 26
#define MMAXNATOM  100;


#define PI 3.141592653589793238

// 変数

///////////////////////////////////////////////////////////////////////////////
// 原子タイプのデータ
struct atomtype
{
	// 原子のタイプ番号
	int atomnum;
	// VDW データ_1
	double gamma;
	// VDW データ_2
	double eata;
	// 電荷
	double e_elect;
};


// 原子のデータ
struct atomdata
{
	// 原子種 ID
	struct atomtype atomID;
	// 座標
	double coord[3];

	// 結合している原子の番号
	int bond_interacting[8];
	// 結合している原子の数
	int num_bond_interacting;

	// 非結合性相互作用しない組
	int not_interacting[20];
	int not_interacting_kminusone[20];
	int not_interacting_kplusone[20];

};



// 二面角のデータ
struct dihedang{
	// 二面角を構成する原子の番号
	// とデータへのインデックス
	int atomnum[6];
};

// 二面角のデータベース
struct dihedang_parm{
	// 相互作用の大きさ
	double V_2;
	// 位相
	double theta;
	// 多重度
	int num;
};

// 剛体のデータ
struct clustdata{
	// 剛体の番号
	int clustID;
	// 原点
	int origin;
	// 終点
	int term[MAXNBH];

	// 終端
	int termflag;
	// 剛体中の原子数
	int nNumAtom;
	// 剛体中の枝の数
	int nNumBranch;
	// 自由度
	int nDegOfFoo;
	// "子"の剛体の番号
	int nNumChildClust[MAXNBH];
	// "親"側の剛体の番号
	int nNumParentClust;

	// 剛体中の二面角の数
	int nNumDihed;
	// 剛体中の二面角の数
	int nNumDiehd_f[MAXNBH];
	// 剛体中の二面角の数
	int nNumDihed_b;
	// 剛体中の二面角の数
	int nNumAtom_f;
	// 剛体中の二面角の数
	int nNumAtom_b;

	double xoord_clust[MAXNBH][MAXNATOM][3];
};
///////////////////////////////////////////////////////////////////////////////

// 残基のデータベース
typedef struct residuedata residuedata;

struct residuedata{
	int nNumAtom;
	int nNumClut;
	int nNumDihed;

	// 原子のデータ
	struct atomdata atm[MAXNATOM];
	// 剛体のデータ
	struct clustdata clt[MAXNCLUST];
	// 二面角のデータ
	struct dihedang dihedang[MAXNDH];

	// 他の残基と結合の原子
	int HeadAtom;
	int TailAtom;
	// 他の残基と結合の剛体
	int HeadClut;
	int TailClut;
};

struct residuedata dataofresidueontable[26];

///////////////////////////////////////////////////////////////////////////////

int NumResInStdResData;

typedef struct elementoftable elementoftable;

struct elementoftable {
	char *nameofthisresidue;

	residuedata dataofthisresdiue;
	elementoftable *next;
};

elementoftable *table[NHASH];

///////////////////////////////////////////////////////////////////////////////

double xoord[MAXNATOM][3];

double MatTransAbsoToLocal[3][3];

///////////////////////////////////////////////////////////////////////////////

// タンパク質のデータ
struct protein{
	// 残基配列のデータ
	char *Sequence[MAXNRESIDUE];
	// タンパク質名
	char *NameProt;

	// 原子種数
	int totalatomtype;
	// 原子数
	int nNumAtom;
	// 原子種数
	int nNumAtomType;
	// 剛体数
	int nNumClut;
	// 残基数
	int nNumResidue;
	// 二面角数
	int nNumDihedAng;
	// 二面角種数
	int nNumDihedType;
	// 非結合性相互作用へのインデックス
	double atmtypeIndex[40];
	// LJ パーム A
	double parmA[40];
	// LJ パーム B
	double parmB[40];

	// タンパク質中の原子のデータ
	struct atomdata atm[MAXNATOM];
	// タンパク質中の原子のタイプのデータ
	struct atomtype atmtype[MAXNATOM];
	// タンパク質中の剛体のデータ
	struct clustdata clt[MAXNCLUST];
	// タンパク質中の二面角のデータ
	struct dihedang dihedang[MAXNDH];
	// タンパク質中の二面角のデータベース
	struct dihedang_parm dihedang_parm[MAXNTD];
}prot;

///////////////////////////////////////////////////////////////////////////////

// 関数
///////////////////////////////////////////////////////////////////////////////
// データの取得を行う関数
int PickData(void);

// インプットデータの作成を行う関数
void CreateData(void);

void TransAbsoToLocal(int nNumAtom1, int nNumAtom2, int nSumNumAtom);
void TransLocalToAbso(int nNumAtom1, int nNumAtom2, int nSumNumAtom, 
                      double delta_dihed);

// インプットの作成を行う関数
void CreateInPut(void);
void CreateCoordFile(FILE *input);
void CreateSeqFile(FILE *input);
void CreateClustFile(FILE *input);
void CreateTopFile(FILE *input);

// ハッシュ中の残基データをの取得を行う関数
elementoftable *LURData(char *nameofthisresidue,
									 int create,
									 residuedata *dataofthisresdiue);
//void LURData(char *nameofthisresidue,
//									 int create,
//									 residuedata *dataofthisresdiue,
//	                                                                 elementoftable *element);

// ハッシュ関数
unsigned int hash(char * str);

// 標準残基データ中にインプットした
// 残基名はあるかの判定を行う関数
int MatchResName(char *nameofthisresidue);
///////////////////////////////////////////////////////////////////////////////


