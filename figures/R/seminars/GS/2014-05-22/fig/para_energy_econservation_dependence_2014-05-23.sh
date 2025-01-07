#!/bin/sh

RES=( dummy \
    GLY \
    ALA \
    ASP \
    GLU \
    LEU \
    ILE \
    ASN \
    GLN \
    VAL \
    HIE \
    SER \
    THR \
    MET \
    CYS \
    PRO \
    ARG \
    LYS \
    PHE \
    TY3 \
    TRP \
    ASH \
    GLH \
    HIP \
    CYX \
    CYM \
    LYN )
    
NRES=$(expr ${#RES[*]} - 1)

clustname0=( dummy "" \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    XX 		      \
    "" 		      \
    "" 		      \
    XX 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    "" 		      \
    ""                \
    "" ) 

clustname=( dummy Term_fixed \
              CH3_Term_fixed \
                  Term_fixed \
                  Term_fixed \
              CH3_Term_fixed \
              CH3_Term_fixed \
              CH3_Term_fixed \
              CH3_Term_fixed \
              CH3_Term_fixed \
                  Term_fixed \
                  Term_fixed \
              CH3_Term_fixed \
              CH3_Term_fixed \
              CH3_Term_fixed \
                  Term_fixed \
                  XX         \
              NH3_Term_fixed \
                  Term_fixed \
                  Term_fixed \
                  Term_fixed \
                  Term_fixed \
                  Term_fixed \
                  Term_fixed \
                  Term_fixed \
                  Term_fixed \
                  Term_fixed )

clustname2=( dummy XX \
    XX 	              \
    XX 		      \
    XX		      \
    XX                \
    XX                \
    NH2_Term_fixed    \
    NH2_Term_fixed    \
    XX 		      \
    XX 		      \
    OH_Term_fixed     \
    CH3_OH_Term_fixed \
    Term_fixed        \
    SH_Term_fixed     \
    XX                \
    Term_fixed        \
    XX 		      \
    XX                \
    OH_Term_fixed     \
    XX                \
    OH_Term_fixed     \
    OH_Term_fixed     \
    XX                \
    XX 		      \
    SH_Term_fixed     \
    NH2_Term_fixed )  \
