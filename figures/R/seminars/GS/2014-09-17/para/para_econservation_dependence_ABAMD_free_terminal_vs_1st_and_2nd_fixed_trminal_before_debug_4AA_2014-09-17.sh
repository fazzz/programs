#!/bin/sh

RES=( dummy  GLY  ALA  ASP  GLU )
    
NRES=$(expr ${#RES[*]} - 1)

clustname0=( dummy ""  ""  ""  "" ) 

clustname=( dummy Term_fixed CH3_Term_fixed Term_fixed Term_fixed  )

clustname2=( dummy XX  XX  XX  XX  )
