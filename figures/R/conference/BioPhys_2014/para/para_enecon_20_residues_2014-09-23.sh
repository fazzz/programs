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
RES[13]=MET   clustnum[13]=2 clustname[29]=Term_free  clustname[30]=CH3_Term_fixed
RES[14]=CYS   clustnum[14]=3 clustname[31]=Term_free  clustname[32]=CH3_Term_fixed clustname[33]=CH3_SH_Term_fixed
RES[15]=PRO   clustnum[15]=2 clustname[34]=Term_free  clustname[35]=CH3_Term_fixed
RES[16]=ARG   clustnum[16]=2 clustname[36]=Term_free clustname[37]=CH3_NH2_Term_fixed
RES[17]=LYS   clustnum[17]=2 clustname[38]=Term_free  clustname[39]=CH3_NH3_Term_fixed
RES[18]=PHE   clustnum[18]=2 clustname[40]=Term_free  clustname[41]=CH3_NH3_Term_fixed
RES[19]=TY3   clustnum[19]=3 clustname[42]=Term_free  clustname[43]=CH3_Term_fixed clustname[44]=CH3_OH_Term_fixed
RES[20]=TRP   clustnum[20]=2 clustname[45]=Term_free  clustname[46]=CH3_Term_fixed

NRES=$(expr ${#RES[*]} - 1)
#NRES=11

nx=5
ny=4

mindt=1
maxdt=15

inc=1

temp0=300
numini=10
numsim=15
