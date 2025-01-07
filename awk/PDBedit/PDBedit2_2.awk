#!/usr/bin/gawk -f

BEGIN{
	usage = "usage S to show,D to delete,F to flag,H to delete specific atom,W to save, Q to quit"
}

$1 ~ /(^(AT)|^(HET))/ {
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

END{
	#print NF
	if (num_atom == 0){
		print "PDBfile cannot read"
	}
	else{
	for(;;){
		print "PDBediter > "
		getline option < "/dev/stdin"
		option = toupper(option)
		if( "" == option){
		}
		else if (option == "H"){
			PDBhdelete(PDBdata)
		}
		else if (option == "N"){
			PDBnum(PDBdata)
		}
		else if (option == "D"){
			PDBdelete(PDBdata)
		}
		else if (option == "L"){
			PDBdelete2(PDBdata)
		}
		else if (option == "A"){
			PDBadd(PDBdata)
		}
		else if (option == "S"){
			PDBshow(PDBdata)
		}
		else if (option == "F"){
			PDBflag(PDBdata)
		}
		else if (option == "W"){
			PDBwrite(PDBdata)
			writeflag = 1
		}
		else if (option == "Q"){
			if (writeflag == 0){
				print "you don't save any change !!"
				print "Will you really quit ?"
				print "yes, no"
				getline userswill < "/dev/stdin"
				userswill = toupper(userswill)
				if (userswill== "YES"){
					exit
				}
				else if (userswill== "NO"){
				}
			}
			else {
				exit
			}
		}
		else{
			print usage
		}
	}
	}
}

func PDBhdelete(PDBdata)
{
	print "to delete specific atom \n set atomname: "
	getline nameofatom  < "/dev/stdin"
	nameofatom = toupper(nameofatom)

	num_h = 0
	j=0
	for (i=0; i<= num_atom; ++i){
		if (PDBatomname[i] !~ nameofatom){
			PDBatomname[j] = PDBatomname[i]
			PDBresiduename[j] = PDBresiduename[i]
			PDBresiduenum[j] = PDBresiduenum[i]
			PDBresiduenum0[j] = PDBresiduenum0[i]
			PDBcoordx[j] = PDBcoordx[i]
			PDBcoordy[j] = PDBcoordy[i]
			PDBcoordz[j] = PDBcoordz[i]
			PDBsenyuu[j] = PDBsenyuu[i]
			PDBondoinsi[j] = PDBondoinsi[i]
			++j
			++newnum_atom
		}
		else {
			++num_h
		}
	}

	for (i=0; i<= num_atom-num_h; ++i){
		PDBdata[i] = "ATOM" "      " i "  " PDBatomname[i] "   " PDBresiduename[i] " " PDBresiduenum0[i] " " PDBresiduenum[i] "      " PDBcoordx[i] " " PDBcoordy[i] " " PDBcoordz[i] " " PDBsenyuu[i] " " PDBondoinsi[i]
	}

	if(num_h == 0) {
		print "any such atom is not contained!!"
	}
}

func PDBdelete(PDBdata)
{
	for (i=0; i<= num_atom; ++i){
		PDBresiduename[i] =  "\b"
		PDBresiduenum[i] =  "\b"
		PDBresiduenum0[i] = "\b"
		PDBsenyuu[i] = "\b"
		PDBondoinsi[i]  = "\b"
	}

	for (i=0; i<= num_atom; ++i){
		PDBdata[i] = "ATOM" "      " i "  " PDBatomname[i] "   " PDBresiduename[i] " " PDBresiduenum0[i] " " PDBresiduenum[i] "      " PDBcoordx[i] " " PDBcoordy[i] " " PDBcoordz[i] " " PDBsenyuu[i] " " PDBondoinsi[i]
	}
	print "look sample !!"
	print PDBdata[1]
}

func PDBdelete2(PDBdata)
{
	for (i=0; i<= num_atom; ++i){
		PDBsenyuu[i] = "\b"
		PDBondoinsi[i]  = "\b"
	}

	for (i=0; i<= num_atom; ++i){
		PDBdata[i] = "ATOM" "      " i "  " PDBatomname[i] "   " PDBresiduename[i] " " PDBresiduenum0[i] " " PDBresiduenum[i] "      " PDBcoordx[i] " " PDBcoordy[i] " " PDBcoordz[i] " " PDBsenyuu[i] " " PDBondoinsi[i]
	}
	print "look sample !!"
	print PDBdata[1]
}


func PDBadd(PDBdata)
{
	for (i=0; i<= num_atom; ++i){
		PDBsenyuu[i] = "1.00 "
		PDBondoinsi[i]  = "0.00 "
	}

	for (i=0; i<= num_atom; ++i){
		PDBdata[i] = "ATOM" "      " i "  " PDBatomname[i] "   " PDBresiduename[i] " " PDBresiduenum0[i] " " PDBresiduenum[i] "      " PDBcoordx[i] " " PDBcoordy[i] " " PDBcoordz[i] " " PDBsenyuu[i] " " PDBondoinsi[i]
	}
	print "look sample !!"
	print PDBdata[1]
}


func PDBshow(PDBdata)
{
	print "specify number of the atom to look"
	getline num < "/dev/stdin"
	for (i=0; i<= num_atom; ++i){
		PDBdata[i] = "ATOM" "      " i "  " PDBatomname[i] "   " PDBresiduename[i] " " PDBresiduenum0[i] " " PDBresiduenum[i] "      " PDBcoordx[i] " " PDBcoordy[i] " " PDBcoordz[i] " " PDBsenyuu[i] " " PDBondoinsi[i]
	}
	print PDBdata[num]
}

func PDBflag(PDBdata)
{
	for(;;){
		print "set atom_num or finish: "
		getline atom_num < "/dev/stdin"
		num_num[++p]=atom_num
		if (atom_num ~ /^f/) {
			ffflag = 1
		}
		else if(atom_num <= num_atom){
			for (i=0; i< num_atom; ++i){
				if (i == atom_num){
					PDBatomname[i] = PDBatomname[i] "#"
					PDBdata[i] = "ATOM" "      " i "  " PDBatomname[i] "   " PDBresiduename[i] " " PDBresiduenum0[i] " " PDBresiduenum[i] "      " PDBcoordx[i] " " PDBcoordy[i] " " PDBcoordz[i] " " PDBsenyuu[i] " " PDBondoinsi[i]
				}
				else {
					PDBdata[i] = PDBdata[i]
				}
			}
		}
		if (ffflag == 1){
			break
		}
	}
	for (i=1;i<p;++i){
		print PDBdata[num_num[i]]
	}
}

func PDBwrite(PDBdata)
{
	print "set outputfilneme: "
	getline filename  < "/dev/stdin"
	for (i=1; i<=num_atom; ++i){
			printf "ATOM %7d  %4-s %-3s%-2s %-4d%12.3f%8.3f%8.3f%8.3f%8.3f \n",i,PDBatomname[i],PDBresiduename[i],PDBresiduenum0[i],PDBresiduenum[i],PDBcoordx[i],PDBcoordy[i],PDBcoordz[i],PDBsenyuu[i],PDBondoinsi[i] >> filename".pdb"
	}
}

func PDBnum(PDBdata)
{
	j=0
	print "to number specific atom \n set atomname: "
	getline nameofatom  < "/dev/stdin"
	nameofatom = toupper(nameofatom)

	for (i=0; i<= num_atom; ++i){
		if (PDBatomname[i] ~ nameofatom){
			++j
			PDBatomname[i] = PDBatomname[i] " " j
			PDBdata[i] = "ATOM" " " i " " PDBatomname[i] " " PDBresiduename[i] " "  PDBcoordx[i] " " PDBcoordy[i] " " PDBcoordz[i]
		}
		else{
			PDBdata[i] = "ATOM" " " i " " PDBatomname[i] " " PDBresiduename[i] " "  PDBcoordx[i] " " PDBcoordy[i] " " PDBcoordz[i]
		}
	}




}

