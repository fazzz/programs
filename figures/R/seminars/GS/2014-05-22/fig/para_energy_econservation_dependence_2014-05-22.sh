#!/bin/sh

RES=( dummy ALA VAL LEU ILE THR \
            SER  \
            ASN GLN \
            LYS \
            CYS \
            ASP ARG GLU GLY HIE MET TRP TY3 PHE PRO \
            ASH CYM CYX GLH HIP LYN )
    
NRES=$(expr ${#RES[*]} - 1)

T=600

numini=10
nCyc=200

clustname0=( dummy "" "" "" "" "" \
                   "" "" \
                   "" "" \
                   "" \
 		   "" \
		   "" "" "" "" "" "" "" "" "" \
                   "" "" "" "" "" "" "" "" )  

clustname=( dummy CH3_Term_fixed CH3_Term_fixed CH3_Term_fixed CH3_Term_fixed \
                  CH3_Term_fixed \
                  Term_fixed  \
                  NH2_Term_fixed NH2_Term_fixed \
                  NH3_Term_fixed \
		  SH_Term_fixed \
    		  Term_fixed Term_fixed Term_fixed Term_fixed Term_fixed Term_fixed Term_fixed Term_fixed Term_fixed Term_fixed \
                  Term_fixed Term_fixed Term_fixed Term_fixed Term_fixed Term_fixed )

clustname2=( dummy XX XX XX XX \
                   CH3_OH_Term_fixed \
                   OH_Term_fixed \
                   XX XX \
                   XX \
 		   XX \
		   XX XX XX XX XX XX XX OH_Term_fixed XX 2014-03-17 \
                   OH_Term_fixed XX XX OH_Term_fixed XX XX  )

psi=( dummy -45  150 )

Ncrds=$(expr ${#phi[*]} - 1)

TCMD=100
       
NUM=5
