/////////////////////////////////////////
//                                     //
//         /Vis_MD/src/DataBase.h      //
//                                     //
/////////////////////////////////////////


// ���

#define MAXNATOM 100   /* ���ҿ��ξ��*/
#define MAXNCLUST 30   /* ��ͳ�٤ξ��*/
#define MAXNDH 30      /* ���̳Ѥξ��*/
#define MAXNTD 10       /* ���̼�ξ��*/
#define MAXNRESIDUE 10  /* �Ĵ���ξ��*/
#define MAXNBH 4        /* ʬ�����ξ��*/
#define NRESIDUE 30     /* �Ĵ��ο� */
#define NHASH 26
#define MMAXNATOM  100;


#define PI 3.141592653589793238

// �ѿ�

///////////////////////////////////////////////////////////////////////////////
// ���ҥ����פΥǡ���
struct atomtype
{
	// ���ҤΥ������ֹ�
	int atomnum;
	// VDW �ǡ���_1
	double gamma;
	// VDW �ǡ���_2
	double eata;
	// �Ų�
	double e_elect;
};


// ���ҤΥǡ���
struct atomdata
{
	// ���Ҽ� ID
	struct atomtype atomID;
	// ��ɸ
	double coord[3];

	// ��礷�Ƥ��븶�Ҥ��ֹ�
	int bond_interacting[8];
	// ��礷�Ƥ��븶�Ҥο�
	int num_bond_interacting;

	// ��������ߺ��Ѥ��ʤ���
	int not_interacting[20];
	int not_interacting_kminusone[20];
	int not_interacting_kplusone[20];

};



// ���̳ѤΥǡ���
struct dihedang{
	// ���̳Ѥ������븶�Ҥ��ֹ�
	// �ȥǡ����ؤΥ���ǥå���
	int atomnum[6];
};

// ���̳ѤΥǡ����١���
struct dihedang_parm{
	// ��ߺ��Ѥ��礭��
	double V_2;
	// ����
	double theta;
	// ¿����
	int num;
};

// ���ΤΥǡ���
struct clustdata{
	// ���Τ��ֹ�
	int clustID;
	// ����
	int origin;
	// ����
	int term[MAXNBH];

	// ��ü
	int termflag;
	// ������θ��ҿ�
	int nNumAtom;
	// ������λޤο�
	int nNumBranch;
	// ��ͳ��
	int nDegOfFoo;
	// "��"�ι��Τ��ֹ�
	int nNumChildClust[MAXNBH];
	// "��"¦�ι��Τ��ֹ�
	int nNumParentClust;

	// ����������̳Ѥο�
	int nNumDihed;
	// ����������̳Ѥο�
	int nNumDiehd_f[MAXNBH];
	// ����������̳Ѥο�
	int nNumDihed_b;
	// ����������̳Ѥο�
	int nNumAtom_f;
	// ����������̳Ѥο�
	int nNumAtom_b;

	double xoord_clust[MAXNBH][MAXNATOM][3];
};
///////////////////////////////////////////////////////////////////////////////

// �Ĵ�Υǡ����١���
typedef struct residuedata residuedata;

struct residuedata{
	int nNumAtom;
	int nNumClut;
	int nNumDihed;

	// ���ҤΥǡ���
	struct atomdata atm[MAXNATOM];
	// ���ΤΥǡ���
	struct clustdata clt[MAXNCLUST];
	// ���̳ѤΥǡ���
	struct dihedang dihedang[MAXNDH];

	// ¾�λĴ�ȷ��θ���
	int HeadAtom;
	int TailAtom;
	// ¾�λĴ�ȷ��ι���
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

// ����ѥ����Υǡ���
struct protein{
	// �Ĵ�����Υǡ���
	char *Sequence[MAXNRESIDUE];
	// ����ѥ���̾
	char *NameProt;

	// ���Ҽ��
	int totalatomtype;
	// ���ҿ�
	int nNumAtom;
	// ���Ҽ��
	int nNumAtomType;
	// ���ο�
	int nNumClut;
	// �Ĵ��
	int nNumResidue;
	// ���̳ѿ�
	int nNumDihedAng;
	// ���̳Ѽ��
	int nNumDihedType;
	// ��������ߺ��ѤؤΥ���ǥå���
	double atmtypeIndex[40];
	// LJ �ѡ��� A
	double parmA[40];
	// LJ �ѡ��� B
	double parmB[40];

	// ����ѥ�����θ��ҤΥǡ���
	struct atomdata atm[MAXNATOM];
	// ����ѥ�����θ��ҤΥ����פΥǡ���
	struct atomtype atmtype[MAXNATOM];
	// ����ѥ�����ι��ΤΥǡ���
	struct clustdata clt[MAXNCLUST];
	// ����ѥ���������̳ѤΥǡ���
	struct dihedang dihedang[MAXNDH];
	// ����ѥ���������̳ѤΥǡ����١���
	struct dihedang_parm dihedang_parm[MAXNTD];
}prot;

///////////////////////////////////////////////////////////////////////////////

// �ؿ�
///////////////////////////////////////////////////////////////////////////////
// �ǡ����μ�����Ԥ��ؿ�
int PickData(void);

// ����ץåȥǡ����κ�����Ԥ��ؿ�
void CreateData(void);

void TransAbsoToLocal(int nNumAtom1, int nNumAtom2, int nSumNumAtom);
void TransLocalToAbso(int nNumAtom1, int nNumAtom2, int nSumNumAtom, 
                      double delta_dihed);

// ����ץåȤκ�����Ԥ��ؿ�
void CreateInPut(void);
void CreateCoordFile(FILE *input);
void CreateSeqFile(FILE *input);
void CreateClustFile(FILE *input);
void CreateTopFile(FILE *input);

// �ϥå�����λĴ�ǡ�����μ�����Ԥ��ؿ�
elementoftable *LURData(char *nameofthisresidue,
									 int create,
									 residuedata *dataofthisresdiue);
//void LURData(char *nameofthisresidue,
//									 int create,
//									 residuedata *dataofthisresdiue,
//	                                                                 elementoftable *element);

// �ϥå���ؿ�
unsigned int hash(char * str);

// ɸ��Ĵ�ǡ�����˥���ץåȤ���
// �Ĵ�̾�Ϥ��뤫��Ƚ���Ԥ��ؿ�
int MatchResName(char *nameofthisresidue);
///////////////////////////////////////////////////////////////////////////////


