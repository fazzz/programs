#!/bin/sh

RES=( dummy ALA )
NRES=$(expr ${#RES[*]} - 1)

T=600

numini=10
nCyc=200

clustname=( dummy )

clustname_free[1]=Term_free         # ALA
clustname_fixed1[1]=CH3_Term_fixed  # ALA
clustname_fixed2[1]=""              # ALA

phi=( dummy -60 -135 )
psi=( dummy -45  150 )

Ncrds=$(expr ${#phi[*]} - 1)

TCMD=100
       
NUM=5
