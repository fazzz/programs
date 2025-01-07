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
    double coord_backbone[10][3];
    double p1dummy[3],p2dummy[3],p3dummy[3];
    double Talpha[3][3],Tbeta[3][3],Rphi[3][3],Rpsi[3][3],TransMat[3][3],TransMatdummy[3][3],TransMatdummy2[3][3],TransMatdummy3[3][3],TransMatdummy4[3][3];
    double thispsi,thisphi;
    double coordummy[3];
    struct protein prot;

    p1[0]=1.470;
    p1[1]=0.000;
    p1[2]=0.000;
    
    //    p2[0]=2.067;
    //    p2[1]=1.206;
    //    p2[2]=0.000;

    p2[0]=3.519;
    p2[1]=1.436;
    p2[2]=0.000;

    p3[0]=1.980;
    p3[1]=1.443;
    p3[2]=0.000;
    
   for (i=0;i<10;++i)
   {
        for (j=0;j<3;++j)
        {
             coord_backbone[i][j]=0.0;
        }
    }
   
   for (i=0;i<10;++i)
   {
        for (j=0;j<3;++j)
        {
            if (i==j)
                TransMat[i][j]=1.0;
            else
                TransMat[i][j]=0.0;
        }
    }

   fomTmat(cos(alpha),sin(alpha),Talpha);
   fomTmat(cos(beta),sin(beta),Tbeta);
    
   printf("********************************************** \n");
   printf("%d  %8.3lf  %8.3lf  %8.3lf \n",0,coord_backbone[0][0],coord_backbone[0][1],coord_backbone[0][2]);
   for (i=0;i<3;++i)
   {


	thispsi=psi[i];
        thisphi=phi[i];
        
        fomRmat(cos(thisphi),sin(thisphi),Rphi);
        fomRmat(cos(thispsi),sin(thispsi),Rpsi);
        
	//        // Calpha
	//         mvult(TransMat,p1,p1dummy);
	//        for (j=0;j<3;++j)
	//        {
	//            coord_backbone[3*i+1][j]=coord_backbone[3*i][j]+p1dummy[j];
	//        }
	//       printf("%d  %8.3lf  %8.3lf  %8.3lf \n",3*i+1,coord_backbone[3*i+1][0],coord_backbone[3*i+1][1],coord_backbone[3*i+1][2]);
	//       // Cprime
	//      mvult(TransMat,p3,p3dummy);
	//       for (j=0;j<3;++j)
	//        {
	//           coord_backbone[3*i+2][j]=coord_backbone[3*i][j]+p3dummy[j];
	//        }
	//       printf("%d  %8.3lf  %8.3lf  %8.3lf \n",3*i+2,coord_backbone[3*i+2][0],coord_backbone[3*i+2][1],coord_backbone[3*i+2][2]);
	//        // N
	//        mmult(Talpha,Rpsi,TransMatdummy);
	//        mmult(Tbeta,Rphi,TransMatdummy2);
	//        mmult(TransMatdummy,TransMatdummy2,TransMatdummy3);	
	//	mmult(TransMat,TransMatdummy3,TransMatdummy4);
	//	mvult(TransMatdummy4,p2,p2dummy);
	//        for (j=0;j<3;++j)
	//      {
	//          coord_backbone[3*i+3][j]=coord_backbone[3*i+1][j]+p2dummy[j];
	//     }
	//       printf("%d  %8.3lf  %8.3lf  %8.3lf \n",3*i+3,coord_backbone[3*i+3][0],coord_backbone[3*i+3][1],coord_backbone[3*i+3][2]);

	//        for (j=0;j<3;++j)
	//        {
	//           for (k=0;k<3;++k)
	//           {
	//	     TransMat[j][k]=TransMatdummy4[j][k];
	//	   }        
	//	}
	//
	//
	//
	//        // Calpha
	//        mmult(Talpha,Rpsi,TransMatdummy);
	//        mmult(Tbeta,Rphi,TransMatdummy2);
	//        mmult(TransMatdummy,TransMatdummy2,TransMatdummy3);	
	//        mmult(TransMat,TransMatdummy3,TransMatdummy4);
	//        mvult(TransMat,p2,p2dummy);
	//        for (j=0;j<3;++j)
	//        {
	//            coord_backbone[3*i+4][j]=coord_backbone[3*i+1][j]+p2dummy[j];
	//        }
	//        printf("%d  %8.3lf  %8.3lf  %8.3lf \n",3*i+3,coord_backbone[3*i+3][0],coord_backbone[3*i+3][1],coord_backbone[3*i+3][2]);
	//
	//        for (j=0;j<3;++j)
	//        {
	//            for (k=0;k<3;++k)
	//            {
	//                TransMat[j][k]=TransMatdummy4[j][k];
	//            }        
	//        }
	//
	//     for (j=0;j<3;++j)
	//     {
	//       prot.coord[3*i][j]=coord_backbone[3*i][j];
	//       prot.coord[3*i+1][j]=coord_backbone[3*i+1][j];
	//       prot.coord[3*i+2][j]=coord_backbone[3*i+2][j];
	//     }
	double Rpsi1[3][3],Rphi1[3][3],Rpsi2[3][3],Rphi2[3][3];
	double TaRs1[3][3],TaRs1Tb[3][3],TaRs1TbRh1[3][3],TaRs1TbRh1Ta[3][3],TaRs1TbRh1TaRs2[3][3],TaRs1TbRh1TaRs2Tb[3][3],TaRs1TbRh1TaRs2TbRh2[3][3],TaRs1TbRh1TaRs2TbRh2Ta[3][3];
	double s1[3],s2[3];

	fomTmat(cos(alpha), sin(alpha), Talpha);
	fomTmat(cos(beta), sin(beta), Tbeta);
	fomRmat(cos(psi[0]), sin(psi[0]), Rpsi1);
	fomRmat(cos(psi[1]), sin(psi[1]), Rpsi2);
	fomRmat(cos(phi[0]), sin(phi[0]), Rphi1);
	fomRmat(cos(phi[1]), sin(phi[1]), Rphi2);

	mmult(Talpha,Rpsi1,TaRs1);
	mmult(TaRs1,Tbeta,TaRs1Tb);
	mmult(TaRs1Tb,Rphi1,TaRs1TbRh1);
	mmult(TaRs1TbRh1,Talpha,TaRs1TbRh1Ta);
	mmult(TaRs1TbRh1Ta,Rpsi2,TaRs1TbRh1TaRs2);
	mmult(TaRs1TbRh1TaRs2,Tbeta,TaRs1TbRh1TaRs2Tb);
	mmult(TaRs1TbRh1TaRs2Tb,Rphi2,TaRs1TbRh1TaRs2TbRh2);
	mmult(TaRs1TbRh1TaRs2TbRh2,Talpha,TaRs1TbRh1TaRs2TbRh2Ta);

	mvult(TaRs1TbRh1,p2,s1);
	mvult(TaRs1TbRh1TaRs2TbRh2,p2,s2);
	
	printf("%8.3lf %8.3lf %8.3lf\n",p2[0],p2[1],p2[2]);
	printf("%8.3lf %8.3lf %8.3lf\n",p2[0]+s1[0],p2[1]+s1[1],p2[2]+s1[2]);
	printf("%8.3lf %8.3lf %8.3lf\n",p2[0]+s1[0]+s2[0],p2[1]+s1[1]+s2[1],p2[2]+s1[2]+s2[2]);

     prot.atomname[3*i]="N";
     prot.resname[3*i]="ALA";
     prot.secname[3*i]="A";
     prot.resnum[i]=1;
     prot.temfac[3*i]=1.00;
     prot.num[3*i]=6.00;
     prot.atomtype[3*i]="N";
     prot.atomname[3*i+1]="C";
     prot.resname[3*i+1]="ALA";
     prot.secname[3*i+1]="A";
     prot.resnum[i]=1;
     prot.temfac[3*i+1]=1.00;
     prot.num[3*i+1]=6.00;
     prot.atomtype[3*i+1]="C";
     prot.atomname[3*i+2]="C";
     prot.resname[3*i+2]="ALA";
     prot.secname[3*i+2]="A";
     prot.resnum[i]=1;
     prot.temfac[3*i+2]=1.00;
     prot.num[3*i+2]=6.00;
     prot.atomtype[3*i+2]="C";
   }

    prot.atomnum=9;
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

    for(i=0;i<prot.atomnum;++i)
    {
      fprintf(out,"ATOM%7d  %-4s%-5s%-1s%4d%12.3lf%8.3lf%8.3lf%6.2lf%6.2lf           %s\n",
	      i,prot.atomname[i],prot.resname[i],prot.secname[i],prot.resnum[i],prot.coord[i][0],prot.coord[i][1],prot.coord[i][2],prot.temfac[i],prot.num[i],prot.atomtype[i]);

      printf("ATOM%7d  %-4s%-5s%-1s%4d%12.3lf%8.3lf%8.3lf%6.2lf%6.2lf           %s\n",
	      i,prot.atomname[i],prot.resname[i],prot.secname[i],prot.resnum[i],prot.coord[i][0],prot.coord[i][1],prot.coord[i][2],prot.temfac[i],prot.num[i],prot.atomtype[i]);
    }

    fclose(out);
}
