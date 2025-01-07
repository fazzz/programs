#!/usr/bin/gawk -f

BEGIN{
	usage = "modify the shape of this file to PDNformat"
}

$1 ~ /(^/.*\..*/)/ {
	for(i=2; i<=NF; ++i){
			if($i ~ /^[ H  HD  Hd  LI  Li  BE  Be  B  C  N  O  F  NA  Na  MG  Mg  AL  Al  SI  Si  P  S  CL  Cl  K  CA  Ca  SC  Sc  TI  Ti  V  CR  Cr  MN  Mn  FE  Fe  CO  Co  NI  Ni  CU  Cu  ZN  Zn  GA  Ga  GE  Ge  AS  As  SE  Se  BR  Br  RB  Rb  SR  Sr  Y  ZR  Zr  NB Nb  MO  Mo  TC  Tc  RU  Ru  RH  Rh  Pd  PD  AG  Ag  CD  Cd  IN  In  SN  Sn  SB  Sb  TE  Te  I  CS  Cs  BA  Ba  HF  Hf  TA  Ta  W  RE  Re  OS  Os  IR  Ir  PT  Pt  AU  Au  HG  Hg  TL  Tl  PB  Pb  BI  Bi  PO  Po  AT  At]/ ){
				PDBatomname[num_atom] = $i
				flag = 1
				break
			if (flag == 1){
				break
			}
		}
	++num_atom
	}

	for(j=i+1; j<=NF; ++j){
			if($j ~ /^[ALA GRY TRY CYS VAL PHE PRO MET ILE LEU ASP GLU LYS ARG SER THR TYR HIS ASN GLN TRP]/ ){
				PDBresiduename[num_atom] = $j
				if ($(j+1) ~ /[A B]/)
				{
					PDBresiduenum0[num_atom]=$(j+1)
					j=j+1
				}
				break
		}
	}

	for (p=j+1;p <= NF;++p){
		if($p ~ /[0-9]*/ ){
			PDBresiduenum[num_atom]=$p
			break
		}
	}


	for (k=2;k <= NF;++k){
		if($k ~ /.*\..*/ && $k+1 ~ /.*\..*/ && $k+2 ~ /.*\..*/ ){
			PDBcoordx[num_atom] = $k
			PDBcoordy[num_atom] = $(k+1)
			PDBcoordz[num_atom] = $(k+2)
			break
		}
	}

	for(l=k+3; l<=NF; ++l){
		if($k ~ /.*\..*/ ){
			PDBsenyuu[num_atom] = $l
			break
		}
	}

	for(m=l+1; m<=NF; ++m){
		if($k ~ /.*\..*/ ){
			PDBondoinsi[num_atom] = $m
			break
		}
	}

}
