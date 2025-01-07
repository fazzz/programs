#!/bin/sh

cp ../DCA/DCA.c                ABA.c               
cp ../DCA/DCA_pick_data.c      ABA_pick_data.c
cp ../DCA/DCA_prepass.c	       ABA_prepass.c
cp ../DCA/DCA_set_imat.c       ABA_set_imat.c
cp ../DCA/DCA_set_lref.c       ABA_set_lref.c
cp ../DCA/DCA_set_tmat.c       ABA_set_tmat.c
cp ../DCA/DCAs_setfrc.c	       ABA_set_frc.c
cp ../DCA/DCAs_trans.c	       ABA_set_trans.c
cp ../DCA/mainDCA.c	       mainABA.c

cp ../DCA/DCA.h                ABA.h
cp ../DCA/DCA_pick_data.h      ABA_pick_data.h
cp ../DCA/DCA_prepass.h	       ABA_prepass.h
cp ../DCA/DCA_set_imat.h       ABA_set_imat.h
cp ../DCA/DCA_set_lref.h       ABA_set_lref.h
cp ../DCA/DCA_set_tmat.h       ABA_set_tmat.h
cp ../DCA/DCAs_setfrc.h	       ABA_set_frc.h
cp ../DCA/DCAs_trans.h	       ABA_set_trans.h

sed -i "s/DCA/ABA/g" ABA.c          
sed -i "s/DCA/ABA/g" ABA_pick_data.c
sed -i "s/DCA/ABA/g" ABA_prepass.c  
sed -i "s/DCA/ABA/g" ABA_set_imat.c 
sed -i "s/DCA/ABA/g" ABA_set_lref.c 
sed -i "s/DCA/ABA/g" ABA_set_tmat.c 
sed -i "s/DCA/ABA/g" ABA_set_frc.c  
sed -i "s/DCA/ABA/g" ABA_set_trans.c
sed -i "s/DCA/ABA/g" mainABA.c      
               
sed -i "s/DCA/ABA/g" ABA.h          
sed -i "s/DCA/ABA/g" ABA_pick_data.h
sed -i "s/DCA/ABA/g" ABA_prepass.h  
sed -i "s/DCA/ABA/g" ABA_set_imat.h 
sed -i "s/DCA/ABA/g" ABA_set_lref.h 
sed -i "s/DCA/ABA/g" ABA_set_tmat.h 
sed -i "s/DCA/ABA/g" ABA_set_frc.h  
sed -i "s/DCA/ABA/g" ABA_set_trans.h
