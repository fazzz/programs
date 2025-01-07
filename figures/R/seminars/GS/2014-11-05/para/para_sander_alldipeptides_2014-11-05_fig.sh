#!/bin/sh

RES=( dummy ) clustnum=( dummy ) clustname=( dummy )
RES[1]=GLY    clustnum[1]=2  clustname[1]=Term_free   clustname[2]=CH3_Term_fixed
RES[2]=ALA    clustnum[2]=2  clustname[3]=Term_free   clustname[4]=CH3_Term_fixed
RES[3]=ASP    clustnum[3]=2  clustname[5]=Term_free   clustname[6]=CH3_Term_fixed
RES[4]=GLU    clustnum[4]=2  clustname[7]=Term_free   clustname[8]=CH3_Term_fixed
RES[5]=LEU    clustnum[5]=2  clustname[9]=Term_free   clustname[10]=CH3_Term_fixed
RES[6]=ILE    clustnum[6]=2  clustname[11]=Term_free  clustname[12]=CH3_Term_fixed
RES[7]=ASN    clustnum[7]=3  clustname[13]=Term_free  clustname[14]=CH3_Term_fixed clustname[15]=CH3_NH2_Term_fixed
RES[8]=GLN    clustnum[8]=3  clustname[16]=Term_free  clustname[17]=CH3_Term_fixed clustname[18]=CH3_NH2_Term_fixed
RES[9]=VAL    clustnum[9]=2  clustname[19]=Term_fixed clustname[20]=CH3_Term_fixed
RES[10]=HIE   clustnum[10]=2 clustname[21]=Term_free  clustname[22]=CH3_Term_fixed
RES[11]=SER   clustnum[11]=3 clustname[23]=Term_free  clustname[24]=CH3_Term_fixed clustname[25]=CH3_OH_Term_fixed
RES[12]=THR   clustnum[12]=3 clustname[26]=Term_free  clustname[27]=CH3_Term_fixed clustname[28]=CH3_OH_Term_fixed
RES[13]=MET   clustnum[13]=2 clustname[29]=XX         clustname[30]=CH3_Term_fixed
RES[14]=CYS   clustnum[14]=3 clustname[31]=Term_free  clustname[32]=CH3_Term_fixed clustname[33]=CH3_SH_Term_fixed
RES[15]=PRO   clustnum[15]=2 clustname[34]=Term_free  clustname[35]=CH3_Term_fixed
RES[16]=ARG   clustnum[16]=3 clustname[36]=Term_free  clustname[37]=CH3_Term_fixed clustname[38]=CH3_NH2_Term_fixed
RES[17]=LYS   clustnum[17]=2 clustname[39]=Term_free  clustname[40]=CH3_NH3_Term_fixed
RES[18]=PHE   clustnum[18]=2 clustname[41]=Term_free  clustname[42]=CH3_Term_fixed
RES[19]=TY3   clustnum[19]=3 clustname[43]=Term_free  clustname[44]=CH3_Term_fixed clustname[45]=CH3_OH_Term_fixed
RES[20]=TRP   clustnum[20]=2 clustname[46]=Term_free  clustname[47]=CH3_Term_fixed
RES[21]=ASH   clustnum[21]=3 clustname[48]=Term_free  clustname[49]=CH3_Term_fixed clustname[50]=CH3_OH_Term_fixed
RES[22]=GLH   clustnum[22]=3 clustname[51]=Term_free  clustname[52]=CH3_Term_fixed clustname[53]=CH3_OH_Term_fixed
RES[23]=HIP   clustnum[23]=2 clustname[54]=Term_free  clustname[55]=CH3_Term_fixed
RES[24]=HID   clustnum[24]=2 clustname[56]=Term_free  clustname[57]=CH3_Term_fixed
RES[25]=CYX   clustnum[25]=2 clustname[58]=Term_free  clustname[59]=CH3_OH_Term_fixed
RES[26]=CYM   clustnum[26]=1 clustname[60]=Term_free  clustname[61]=CH3_Term_fixed 
RES[27]=LYN   clustnum[27]=3 clustname[62]=Term_free  clustname[63]=CH3_Term_fixed clustname[64]=CH3_NH2_Term_fixed

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
