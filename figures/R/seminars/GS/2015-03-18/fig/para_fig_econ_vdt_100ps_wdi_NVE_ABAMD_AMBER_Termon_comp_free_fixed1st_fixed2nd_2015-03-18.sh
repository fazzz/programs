#!/bin/sh

RES=( dummy ) clustnum=( dummy ) clustname=( dummy )
RES[1]=GLY  clustnum[1]=2 clustname[1]=Term_free  clustname[2]=CH3_Term_fixed
RES[2]=ALA  clustnum[2]=2 clustname[3]=Term_free  clustname[4]=CH3_Term_fixed
RES[3]=VAL  clustnum[3]=2 clustname[5]=Term_fixed clustname[6]=CH3_Term_fixed
RES[4]=LEU  clustnum[4]=2 clustname[7]=Term_free  clustname[8]=CH3_Term_fixed
RES[5]=ILE  clustnum[5]=2 clustname[9]=Term_free  clustname[10]=CH3_Term_fixed
RES[6]=PRO  clustnum[6]=2 clustname[11]=Term_free clustname[12]=CH3_Term_fixed

RES[7]=MET  clustnum[7]=2  clustname[13]=Term_S_free clustname[14]=CH3_Term_fixed
RES[8]=CYS  clustnum[8]=3  clustname[15]=Term_free  clustname[16]=CH3_Term_fixed clustname[17]=CH3_SH_Term_fixed
RES[9]=CYX  clustnum[9]=2  clustname[18]=Term_free  clustname[19]=CH3_OH_Term_fixed
RES[10]=CYM clustnum[10]=2 clustname[20]=Term_free  clustname[21]=CH3_Term_fixed 
RES[11]=SER clustnum[11]=3 clustname[22]=Term_free  clustname[23]=CH3_Term_fixed clustname[24]=CH3_OH_Term_fixed
RES[12]=THR clustnum[12]=3 clustname[25]=Term_free  clustname[26]=CH3_Term_fixed clustname[27]=CH3_OH_Term_fixed

RES[13]=ASN clustnum[13]=3 clustname[28]=Term_free  clustname[29]=CH3_Term_fixed clustname[30]=CH3_NH2_Term_fixed
RES[14]=ASH clustnum[14]=3 clustname[31]=Term_free  clustname[32]=CH3_Term_fixed clustname[33]=CH3_OH_Term_fixed
RES[15]=GLN clustnum[15]=3 clustname[34]=Term_free  clustname[35]=CH3_Term_fixed clustname[36]=CH3_NH2_Term_fixed
RES[16]=GLH clustnum[16]=3 clustname[37]=Term_free  clustname[38]=CH3_Term_fixed clustname[39]=CH3_OH_Term_fixed
RES[17]=ASP clustnum[17]=2 clustname[40]=Term_free  clustname[41]=CH3_Term_fixed
RES[18]=GLU clustnum[18]=2 clustname[42]=Term_free  clustname[43]=CH3_Term_fixed

RES[19]=LYS clustnum[19]=2 clustname[44]=Term_free  clustname[45]=CH3_NH3_Term_fixed
RES[20]=LYN clustnum[20]=3 clustname[46]=Term_free  clustname[47]=CH3_Term_fixed clustname[48]=CH3_NH2_Term_fixed
RES[21]=ARG clustnum[21]=3 clustname[49]=Term_C_free clustname[50]=Term_C_free_CH3_fixd clustname[51]=CH3_NH2_Term_fixed
RES[22]=HIE clustnum[22]=2 clustname[52]=Term_free  clustname[53]=CH3_Term_fixed
RES[23]=HIP clustnum[23]=2 clustname[54]=Term_free  clustname[55]=CH3_Term_fixed
RES[24]=HID clustnum[24]=2 clustname[56]=Term_free  clustname[57]=CH3_Term_fixed

RES[25]=PHE clustnum[25]=2 clustname[58]=Term_free  clustname[59]=CH3_Term_fixed
RES[26]=TY3 clustnum[26]=3 clustname[60]=Term_free  clustname[61]=CH3_Term_fixed clustname[62]=CH3_OH_Term_fixed
RES[27]=TRP clustnum[27]=2 clustname[63]=Term_free  clustname[64]=CH3_Term_fixed

NRES=$(expr ${#RES[*]} - 1)
#NRES=26

#nx=5
#ny=6

mindt=1
maxdt=15

inc=1

temp0=300
numini=10
numsim=15
