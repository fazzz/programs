

#ifndef INCLUDE_mkTA
#define INCLUDE_mkTA

#define S1 0
#define S2 1
#define S3 2

#define PHI 0
#define PSI 1
#define OMEGA 2
#define KAI 3
#define ALL 10

#define NTERM 1
#define NOTERM 2
#define CTERM 3

int mkTACCMinput_set_atomnumi_atomnuml(int atomnumj,int atomnumk,
				       char *atomnamej, char *atomnamek,
				       char *resnamej,  char *resnamek,
				       int *atomnumi,  int *atomnuml, int *mesg, int termflag);

int judge_dihedtype(char *atomnamej, char *atomnamek/*, char *treej, char *treek*/, char *resnamej,  char *resnamek,
		    char *atomnamei, char *atomnamel/*, char *treei, char *treel*/,int termflag);

#endif
