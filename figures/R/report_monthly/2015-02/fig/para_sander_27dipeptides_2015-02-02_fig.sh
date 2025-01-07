#!/bin/sh

RES=( dummy ) clustnamefree=( dummy ) clustname1st=( dummy ) clustname2nd=( dummy )
RES[1]=GLY clustnamefree[1]=Term_free  clustname1st[1]=CH3_Term_fixed clustname2nd[1]=""
RES[2]=ALA clustnamefree[2]=Term_free  clustname1st[2]=CH3_Term_fixed clustname2nd[2]=""
RES[3]=VAL clustnamefree[3]=Term_fixed clustname1st[3]=CH3_Term_fixed clustname2nd[3]=""
RES[4]=LEU clustnamefree[4]=Term_free  clustname1st[4]=CH3_Term_fixed clustname2nd[4]=""
RES[5]=ILE clustnamefree[5]=Term_free  clustname1st[5]=CH3_Term_fixed clustname2nd[5]=""
RES[6]=PRO clustnamefree[6]=Term_free  clustname1st[6]=CH3_Term_fixed clustname2nd[6]=""

RES[7]=MET  clustnamefree[7]=Term_S_free clustname1st[7]=CH3_Term_fixed    clustname2nd[7]=""
RES[8]=CYS  clustnamefree[8]=Term_free   clustname1st[8]=CH3_Term_fixed    clustname2nd[8]=CH3_SH_Term_fixed
RES[9]=CYX  clustnamefree[9]=Term_free   clustname1st[9]=CH3_OH_Term_fixed clustname2nd[9]=""
RES[10]=CYM clustnamefree[10]=Term_free  clustname1st[10]=CH3_Term_fixed   clustname2nd[10]=""
RES[11]=SER clustnamefree[11]=Term_free  clustname1st[11]=CH3_Term_fixed   clustname2nd[11]=CH3_OH_Term_fixed
RES[12]=THR clustnamefree[12]=Term_free  clustname1st[12]=CH3_Term_fixed   clustname2nd[12]=CH3_OH_Term_fixed

RES[13]=ASN clustnamefree[13]=Term_free  clustname1st[13]=CH3_Term_fixed clustname2nd[13]=CH3_NH2_Term_fixed
RES[14]=ASH clustnamefree[14]=Term_free  clustname1st[14]=CH3_Term_fixed clustname2nd[14]=CH3_OH_Term_fixed
RES[15]=GLN clustnamefree[15]=Term_free  clustname1st[15]=CH3_Term_fixed clustname2nd[15]=CH3_NH2_Term_fixed
RES[16]=GLH clustnamefree[16]=Term_free  clustname1st[16]=CH3_Term_fixed clustname2nd[16]=CH3_OH_Term_fixed
RES[17]=ASP clustnamefree[17]=Term_free  clustname1st[17]=CH3_Term_fixed clustname2nd[17]=""
RES[18]=GLU clustnamefree[18]=Term_free  clustname1st[18]=CH3_Term_fixed clustname2nd[18]=""

RES[19]=LYS clustnamefree[19]=Term_free   clustname1st[19]=CH3_NH3_Term_fixed   clustname2nd[19]=""
RES[20]=LYN clustnamefree[20]=Term_free   clustname1st[20]=CH3_Term_fixed       clustname[20]=CH3_NH2_Term_fixed
RES[21]=ARG clustnamefree[21]=Term_C_free clustname1st[21]=Term_C_free_CH3_fixd clustname[21]=CH3_NH2_Term_fixed
RES[22]=HIE clustnamefree[22]=Term_free   clustname1st[22]=CH3_Term_fixed       clustname2nd[22]=""
RES[23]=HIP clustnamefree[23]=Term_free   clustname1st[23]=CH3_Term_fixed	clustname2nd[23]=""
RES[24]=HID clustnamefree[24]=Term_free   clustname1st[24]=CH3_Term_fixed	clustname2nd[24]=""

RES[25]=PHE clustnamefree[25]=Term_free  clustname1st[25]=CH3_Term_fixed clustname2nd[25]=""
RES[26]=TY3 clustnamefree[26]=Term_free  clustname1st[26]=CH3_Term_fixed clustname2nd[26]=CH3_OH_Term_fixed
RES[27]=TRP clustnamefree[27]=Term_free  clustname1st[27]=CH3_Term_fixed clustname2nd[27]=""

clustname=( dummy )
clustnamemod[1]=tune_s_0.1    # GLY
clustnamemod[2]=tune_s_0.1    # ALA
clustnamemod[3]=tune_s_0.1    # VAL
clustnamemod[4]=tune_s_0.1    # LEU
clustnamemod[5]=tune_s_0.1    # ILE
clustnamemod[6]=tune_s_0.1    # PRO

clustnamemod[7]=tune_s_0.1    # MET
clustnamemod[8]=tune_s_0.1    # CYS
clustnamemod[9]=tune_s_0.1    # CYX
clustnamemod[10]=tune_s_0.1   # CYM
clustnamemod[11]=tune_s_0.1   # SER
clustnamemod[12]=tune_s_0.1   # THR

clustnamemod[13]=tune_s_0.1   # ASN
clustnamemod[14]=tune_s_0.1   # ASH
clustnamemod[15]=tune_s_0.1   # GLN
clustnamemod[16]=tune_s_0.1   # GLH
clustnamemod[17]=tune_s_0.1   # ASP
clustnamemod[18]=tune_s_0.1   # GLU

clustnamemod[19]=tune_s_0.1   # LYS
clustnamemod[20]=tune_s_0.1   # LYN
clustnamemod[21]=tune_s_0.1   # ARG
clustnamemod[22]=tune_s_0.1   # HIE
clustnamemod[23]=tune_s_0.1   # HIP
clustnamemod[24]=tune_s_0.1   # HID

clustnamemod[25]=tune_s_0.1   # PHE
#clustnamemod[26]=tune_s_0.1  # TY3
clustnamemod[27]=tune_s_0.1   # TRP

NRES=$(expr ${#RES[*]} - 1)

mindt=1
maxdt=15

inc=1

temp0=300
numini=10
numsim=15
