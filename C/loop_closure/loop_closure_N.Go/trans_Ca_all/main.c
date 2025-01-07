#include <stdio.h>
#include <stdlib.h>
#include "loop_closure.h"


struct protein
{
  int atomnum;
  double coord[100][3];
  double temfac[100];
  double num[100];
  char *atomname[100];
  char *atomtype[100];
  char *resname[100];
  char *secname[100];
  int resnum[100];
};

void pick_data(char *inputfilename,double s[3],double u[3],double v[3]);
void build_backbone(char *ouputfilename,double psi[3],double phi[3]);
void write_pdb(struct protein prot, char *outputfilename);

int main(int argc, char *argv[])
{
	char *inputfilename;
	char *outputfilename;
	int i,num_root;
	double s[3], u[3], v[3];
	double dihedang[4][6];
	double psi[3],phi[3];
	double psi1[4],psi2[4],psi3[4],phi1[4],phi2[4],phi3[4];
	double coord_calpha[100][3];


	if (argc>1)
	{
	  inputfilename=argv[1];
	}
	else
	{
	  printf("no input file name !! \n");
	  exit(1);
	}
	sprintf(outputfilename,"%s_out",inputfilename);
    
	pick_data(inputfilename,s,u,v);

	num_root=loop_closure(s,u,v,psi1,psi2,psi3,phi1,phi2,phi3);

	for (i=0;i<num_root;++i)
	{
	  psi[0]=psi1[i]*pi/180.0;
	  psi[1]=psi2[i]*pi/180.0;
	  psi[2]=psi3[i]*pi/180.0;
	  phi[0]=phi1[i]*pi/180.0;
	  phi[1]=phi2[i]*pi/180.0;
	  phi[2]=phi3[i]*pi/180.0;
	  build_backbone(outputfilename,psi,phi);
	}

	//	return 0;
}

void pick_data(char *inputfilename,double s[3],double u[3],double v[3])
{
    int i;
    double x,y,z;
    FILE *input1;


    if ((input1=fopen(inputfilename,"r")) == NULL)
    {
        printf("can't open %s",inputfilename);
        exit(1);
    }
    
    
    fscanf(input1,"%lf %lf %lf",&x,&y,&z);
    s[0]=x;s[1]=y;s[2]=z;    
    fscanf(input1,"%lf %lf %lf",&x,&y,&z);
    u[0]=x;u[1]=y;u[2]=z;    
    fscanf(input1,"%lf %lf %lf",&x,&y,&z);
    v[0]=x;v[1]=y;v[2]=z;    
       
    fclose(input1);
    
}

void build_backbone(char *outputfilename,double psi[3],double phi[3])
{
    FILE *out;
    int i, j, k;
    double vectCatoCp[4][3],vectCatoNn[4][3],vectCatoCan[4][3];
    double Talpha[3][3],Tbeta[3][3];
    double Rpsi1[3][3],Rphi1[3][3],Rpsi2[3][3],Rphi2[3][3],Rpsi3[3][3],Rphi3[3][3];
    double TaRs1[3][3],TaRs1Tb[3][3],TaRs1TbRh1[3][3],TaRs1TbRh1Ta[3][3],TaRs1TbRh1TaRs2[3][3],TaRs1TbRh1TaRs2Tb[3][3],TaRs1TbRh1TaRs2TbRh2[3][3],TaRs1TbRh1TaRs2TbRh2Ta[3][3],TaRs1TbRh1TaRs2TbRh2TaRs3[3][3],TaRs1TbRh1TaRs2TbRh2TaRs3Tb[3][3],TaRs1TbRh1TaRs2TbRh2TaRs3TbRh3[3][3];
 
    struct protein prot;

    vectCatoCp[0][0]=1.530;
    vectCatoCp[0][1]=0.000;
    vectCatoCp[0][2]=0.000;
    
    vectCatoCan[0][0]=3.519;
    vectCatoCan[0][1]=1.436;
    vectCatoCan[0][2]=0.000;

   vectCatoNn[0][0]=2.067;
   vectCatoNn[0][1]=1.206;
   vectCatoNn[0][2]=0.000;
        
   fomTmat(cos(alpha), sin(alpha), Talpha);
   fomTmat(cos(beta), sin(beta), Tbeta);
   fomRmat(cos(psi[0]), sin(psi[0]), Rpsi1);
   fomRmat(cos(psi[1]), sin(psi[1]), Rpsi2);
   fomRmat(cos(psi[2]), sin(psi[2]), Rpsi3);
   fomRmat(cos(phi[0]), sin(phi[0]), Rphi1);
   fomRmat(cos(phi[1]), sin(phi[1]), Rphi2);
   fomRmat(cos(phi[2]), sin(phi[2]), Rphi3);

   mmult(Talpha,Rpsi1,TaRs1);
   mmult(TaRs1,Tbeta,TaRs1Tb);
   mmult(TaRs1Tb,Rphi1,TaRs1TbRh1);
   mmult(TaRs1TbRh1,Talpha,TaRs1TbRh1Ta);
   mmult(TaRs1TbRh1Ta,Rpsi2,TaRs1TbRh1TaRs2);
   mmult(TaRs1TbRh1TaRs2,Tbeta,TaRs1TbRh1TaRs2Tb);
   mmult(TaRs1TbRh1TaRs2Tb,Rphi2,TaRs1TbRh1TaRs2TbRh2);
   mmult(TaRs1TbRh1TaRs2TbRh2,Talpha,TaRs1TbRh1TaRs2TbRh2Ta);
   mmult(TaRs1TbRh1TaRs2TbRh2Ta,Rpsi3,TaRs1TbRh1TaRs2TbRh2TaRs3);
   mmult(TaRs1TbRh1TaRs2TbRh2TaRs3,Tbeta,TaRs1TbRh1TaRs2TbRh2TaRs3Tb);
   mmult(TaRs1TbRh1TaRs2TbRh2TaRs3Tb,Rphi3,TaRs1TbRh1TaRs2TbRh2TaRs3TbRh3);

   mvult(TaRs1TbRh1,vectCatoCp[0],vectCatoCp[1]);
   mvult(TaRs1TbRh1,vectCatoNn[0],vectCatoNn[1]);
   mvult(TaRs1TbRh1,vectCatoCan[0],vectCatoCan[1]);
   mvult(TaRs1TbRh1TaRs2TbRh2,vectCatoCp[0],vectCatoCp[2]);
   mvult(TaRs1TbRh1TaRs2TbRh2,vectCatoNn[0],vectCatoNn[2]);
   mvult(TaRs1TbRh1TaRs2TbRh2,vectCatoCan[0],vectCatoCan[2]);
   mvult(TaRs1TbRh1TaRs2TbRh2TaRs3TbRh3,vectCatoCp[0],vectCatoCp[3]);
   mvult(TaRs1TbRh1TaRs2TbRh2TaRs3TbRh3,vectCatoNn[0],vectCatoNn[3]);
   mvult(TaRs1TbRh1TaRs2TbRh2TaRs3TbRh3,vectCatoCan[0],vectCatoCan[3]);

   for (j=0;j<3;++j)
   {
     prot.coord[0][j]= 0.0;
    }
   for (j=0;j<3;++j)
   {
     prot.coord[1][j]= vectCatoCp[0][j];
    }
   for (j=0;j<3;++j)
   {
     prot.coord[2][j]= vectCatoNn[0][j];
   }
   for (j=0;j<3;++j)
   {
     prot.coord[3][j]= vectCatoCan[0][j];
   }
   for (j=0;j<3;++j)
   {
     prot.coord[4][j]= vectCatoCan[0][j]+vectCatoCp[1][j];
   }
   for (j=0;j<3;++j)
   {
     prot.coord[5][j]= vectCatoCan[0][j]+vectCatoNn[1][j];
   }
   for (j=0;j<3;++j)
   {
     prot.coord[6][j]= vectCatoCan[0][j]+vectCatoCan[1][j];
   }
   for (j=0;j<3;++j)
   {
      prot.coord[7][j]= vectCatoCan[0][j]+vectCatoCan[1][j]+vectCatoCp[2][j];
   }
   for (j=0;j<3;++j)
   {
      prot.coord[8][j]= vectCatoCan[0][j]+vectCatoCan[1][j]+vectCatoNn[2][j];
   }
   for (j=0;j<3;++j)
   {
      prot.coord[9][j]= vectCatoCan[0][j]+vectCatoCan[1][j]+vectCatoCan[2][j];
   }
   for (j=0;j<3;++j)
   {
     prot.coord[10][j]= vectCatoCan[0][j]+vectCatoCan[1][j]+vectCatoCan[2][j]+vectCatoCp[3][j];
   }
   for (j=0;j<3;++j)
   {
      prot.coord[11][j]= vectCatoCan[0][j]+vectCatoCan[1][j]+vectCatoCan[2][j]+vectCatoNn[3][j];
   }
   for (j=0;j<3;++j)
   {
      prot.coord[12][j]= vectCatoCan[0][j]+vectCatoCan[1][j]+vectCatoCan[2][j]+vectCatoCan[3][j];
   }

   prot.atomname[0]="C";
   prot.resname[0]="ALA";
   prot.secname[0]="A";
   prot.resnum[0]=1;
   prot.temfac[0]=1.00;
   prot.num[0]=6.00;
   prot.atomtype[0]="C";

   for (i=0;i<4;++i)
   {
     prot.atomname[3*i+1]="C";
     prot.resname[3*i+1]="ALA";
     prot.secname[3*i+1]="A";
     prot.resnum[3*i+1]=1;
     prot.temfac[3*i+1]=1.00;
     prot.num[3*i+1]=6.00;
     prot.atomtype[3*i+1]="C";

     prot.atomname[3*i+2]="C";
     prot.resname[3*i+2]="ALA";
     prot.secname[3*i+2]="A";
     prot.resnum[3*i+2]=1;
     prot.temfac[3*i+2]=1.00;
     prot.num[3*i+2]=6.00;
     prot.atomtype[3*i+2]="C";

     prot.atomname[3*i+3]="C";
     prot.resname[3*i+3]="ALA";
     prot.secname[3*i+3]="A";
     prot.resnum[3*i+3]=1;
     prot.temfac[3*i+3]=1.00;
     prot.num[3*i+3]=6.00;
     prot.atomtype[3*i+3]="C";
   }

   prot.atomnum=13;
   for(i=0;i<prot.atomnum;++i)
   {
     printf("%13.3lf%8.3lf%8.3lf\n",prot.coord[i][0],prot.coord[i][1],prot.coord[i][2]);
    }
   write_pdb(prot,outputfilename);
       
}

void write_pdb(struct protein prot, char *outputfilename)
{
  int i;
  FILE *out;

    if ((out=fopen(outputfilename,"a")) == NULL)
    {
        printf("can't open %s",outputfilename);
        exit(1);
    }

    printf("********************************************** \n");
    for(i=0;i<prot.atomnum;++i)
    {
      fprintf(out,"ATOM%7d  %-4s%-5s%-1s%4d%11.3lf%8.3lf%8.3lf%6.2lf%6.2lf           %s\n",
	      i,prot.atomname[i],prot.resname[i],prot.secname[i],prot.resnum[i],prot.coord[i][0],prot.coord[i][1],prot.coord[i][2],prot.temfac[i],prot.num[i],prot.atomtype[i]);

      printf("ATOM%7d  %-4s%-5s%-1s%4d%11.3lf%8.3lf%8.3lf%6.2lf%6.2lf           %s\n",
	      i,prot.atomname[i],prot.resname[i],prot.secname[i],prot.resnum[i],prot.coord[i][0],prot.coord[i][1],prot.coord[i][2],prot.temfac[i],prot.num[i],prot.atomtype[i]);
    }

    fclose(out);
}
