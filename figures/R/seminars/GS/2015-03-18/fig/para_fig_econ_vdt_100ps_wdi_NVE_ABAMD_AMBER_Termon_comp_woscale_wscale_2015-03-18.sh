#!/bin/sh

integtype=PC6
integflag=--predictcorrect

RES=( dummy )
RES[1]=GLY
RES[2]=ALA
RES[3]=VAL
RES[4]=LEU
RES[5]=ILE
RES[6]=PRO

RES[7]=MET
RES[8]=CYS
RES[9]=CYX
RES[10]=CYM
RES[11]=SER
RES[12]=THR

RES[13]=ASN
RES[14]=ASH
RES[15]=GLN
RES[16]=GLH
RES[17]=ASP
RES[18]=GLU

RES[19]=LYS
RES[20]=LYN
RES[21]=ARG
RES[22]=HIE
RES[23]=HIP
RES[24]=HID

RES[25]=PHE
RES[26]=TY3
RES[27]=TRP


clustname=( dummy )
clustname[1]=tune_s_0.1          # GLY
clustname[2]=tune_s_0.1          # ALA
clustname[3]=tune_s_0.1          # VAL
clustname[4]=tune_s_0.1          # LEU
clustname[5]=tune_s_0.1          # ILE
clustname[6]=tune_s_0.1          # PRO

clustname[7]=tune_s_0.1          # MET
clustname[8]=tune_s_0.1          # CYS
clustname[9]=tune_s_0.1          # CYX
clustname[10]=tune_s_0.1         # CYM
clustname[11]=tune_s_0.1         # SER
clustname[12]=tune_s_0.1         # THR

clustname[13]=tune_s_0.1         # ASN
clustname[14]=tune_s_0.1         # ASH
clustname[15]=tune_s_0.1         # GLN
clustname[16]=tune_s_0.1         # GLH
clustname[17]=tune_s_0.1         # ASP
clustname[18]=tune_s_0.1         # GLU

clustname[19]=tune_s_0.1         # LYS
clustname[20]=tune_s_0.1         # LYN
clustname[21]=tune_s_0.1         # ARG
clustname[22]=tune_s_0.1         # HIE
clustname[23]=tune_s_0.1         # HIP
clustname[24]=tune_s_0.1         # HID

clustname[25]=tune_s_0.1         # PHE
clustname[26]=tune_s_0.1         # TY3
clustname[27]=tune_s_0.1         # TRP


NRES=$(expr ${#RES[*]} - 1)

mindt=1
maxdt=15

inc=1

temp0=300
numini=10
numsim=15
