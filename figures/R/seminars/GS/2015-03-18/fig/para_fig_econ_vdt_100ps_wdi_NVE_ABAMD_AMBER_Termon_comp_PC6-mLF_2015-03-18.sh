#!/bin/sh

integtypePC6=""
integtypemLF=_LeapFrog

RES=( dummy ) clustname=( dummy )
RES[1]=ALA    clustname[1]=Term_free clustname[2]=Term_free

RES[1]=GLY    clustname[1]=Term_free        clustname[2]=Term_free        
RES[2]=ALA    clustname[3]=Term_free        clustname[4]=Term_free        
RES[3]=VAL    clustname[5]=Term_fixed       clustname[6]=Term_fixed       
RES[4]=LEU    clustname[7]=Term_free        clustname[8]=Term_free        
RES[5]=ILE    clustname[9]=Term_free        clustname[10]=Term_free        
RES[6]=PRO    clustname[11]=Term_free       clustname[12]=Term_free        
	                                                                  
RES[7]=MET    clustname[13]=Term_S_free     clustname[14]=Term_S_free      
RES[8]=CYS    clustname[15]=Term_free       clustname[16]=Term_free        
RES[9]=CYX    clustname[17]=Term_free        clustname[18]=Term_free        
RES[10]=CYM   clustname[19]=Term_free       clustname[20]=Term_free       
RES[11]=SER   clustname[21]=Term_free       clustname[22]=Term_free       
RES[12]=THR   clustname[23]=Term_free       clustname[24]=Term_free       
	                                                                  
RES[13]=ASN   clustname[25]=Term_free       clustname[26]=Term_free       
RES[14]=ASH   clustname[27]=Term_free       clustname[28]=Term_free       
RES[15]=GLN   clustname[29]=Term_free       clustname[30]=Term_free       
RES[16]=GLH   clustname[31]=Term_free       clustname[32]=Term_free       
RES[17]=ASP   clustname[33]=Term_free       clustname[34]=Term_free       
RES[18]=GLU   clustname[35]=Term_free       clustname[36]=Term_free       
	                                                                  
RES[19]=LYS   clustname[37]=Term_free       clustname[38]=Term_free       
RES[20]=LYN   clustname[39]=Term_free       clustname[40]=Term_free       
RES[21]=ARG   clustname[41]=Term_C_free     clustname[42]=Term_C_free     
RES[22]=HIE   clustname[43]=Term_free       clustname[44]=Term_free       
RES[23]=HIP   clustname[45]=Term_free       clustname[46]=Term_free       
RES[24]=HID   clustname[47]=Term_free       clustname[48]=Term_free       
	                                                                  
RES[25]=PHE   clustname[49]=Term_free       clustname[50]=Term_free       
RES[26]=TY3   clustname[51]=Term_free       clustname[52]=Term_free       
RES[27]=TRP   clustname[53]=Term_free       clustname[54]=Term_free       

NRES=$(expr ${#RES[*]} - 1)
#NRES=26

#nx=5
#ny=6

mindt=1
maxdt=15

inc=1

temp0=300
numini=10
numsim=15
