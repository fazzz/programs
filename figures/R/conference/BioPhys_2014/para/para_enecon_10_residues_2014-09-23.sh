#!/bin/sh

GLY ALA ASP LYS ILE VAL LEU ASN THR CYS 

RES=( dummy ) clustnum=( dummy ) clustname=( dummy )
RES[1]=GLY    clustnum[1]=2  clustname[1]=Term_free   clustname[2]=CH3_Term_fixed
RES[2]=ALA    clustnum[2]=2  clustname[3]=Term_free   clustname[4]=CH3_Term_fixed
RES[3]=ASP    clustnum[3]=2  clustname[5]=Term_free   clustname[6]=CH3_Term_fixed
RES[4]=LYS    clustnum[4]=2  clustname[7]=Term_free   clustname[8]=CH3_NH3_Term_fixed
RES[5]=ILE    clustnum[5]=2  clustname[9]=Term_free   clustname[10]=CH3_Term_fixed
RES[6]=VAL    clustnum[6]=2  clustname[11]=Term_fixed clustname[12]=CH3_Term_fixed
RES[7]=LEU    clustnum[7]=2  clustname[13]=Term_free  clustname[14]=CH3_Term_fixed
RES[8]=ASN    clustnum[8]=3  clustname[15]=Term_free  clustname[16]=CH3_Term_fixed clustname[17]=CH3_NH2_Term_fixed
RES[9]=THR    clustnum[9]=3  clustname[18]=Term_free  clustname[19]=CH3_Term_fixed clustname[20]=CH3_OH_Term_fixed
RES[10]=CYS   clustnum[10]=3 clustname[21]=Term_free  clustname[22]=CH3_Term_fixed clustname[23]=CH3_SH_Term_fixed

NRES=$(expr ${#RES[*]} - 1)
#NRES=11

nx=5
ny=2

mindt=1
maxdt=15

inc=1

temp0=300
numini=10
numsim=15
