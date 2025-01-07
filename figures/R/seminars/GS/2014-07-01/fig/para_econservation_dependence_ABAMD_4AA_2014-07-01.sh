#!/bin/sh

RES=( dummy ALA LEU ILE GLN  )
NRES=$(expr ${#RES[*]} - 1)

T=600

numini=10
nCyc=200

clustname_free=(   dummy )
clustname_fixed1=( dummy )
clustname_fixed2=( dummy )

clustname_free[1]=Term_free         # ALA
clustname_fixed1[1]=CH3_Term_fixed  # ALA
clustname_fixed2[1]=""              # ALA

clustname_free[2]=Term_free         # LEU
clustname_fixed1[2]=CH3_Term_fixed  # LEU
clustname_fixed2[2]=""              # LEU

clustname_free[3]=Term_free         # ILE
clustname_fixed1[3]=CH3_Term_fixed  # ILE
clustname_fixed2[3]=""              # ILE

clustname_free[4]=Term_free         # GLN
clustname_fixed1[4]=CH3_Term_fixed  # GLN
clustname_fixed2[4]=NH2_Term_fixed  # GLN

phi=( dummy -60 -135 )
psi=( dummy -45  150 )

Ncrds=$(expr ${#phi[*]} - 1)

TCMD=100
       
NUM=5
