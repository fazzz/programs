#
# Makefile
#

# definitions
CC      = gcc -g
OBJS	= mainABAMD_NH_GOLMAA_dbasin.o ABA_set_lref.o ABA_set_imat.o ABA_set_tmat.o ABA_pick_data.o ABA_set_trans.o ABA_set_frc.o ABA_mainpass.o ABA_backpass.o ABA_update.o ABA_prepass.o ABA_integ.o ABA_calcattfrc.o ABA_calcKineE.o ABA_Nose-Hoover.o ABA.o
LIBS	= ABA.h
TARGET	= ABAMD_NH_GOLMAA

# rules of generations

$(TARGET):  $(OBJS) $(LIBS)
	$(CC) -o $@ $(OBJS) -L${HOME}/lib -I${HOME}/include  -L${HOME}/mylib -I${HOME}/myinclude -Xlinker -rpath -Xlinker ${HOME}/mylib -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0 -I${HOME}/lib/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -lstringOTF -lMD -lABA_Nose-Hoover -lABA_gtree -lABA_Inverse_mainpass -lABA_Inverse -lABA_Inverse_backpass -lnetcdf_mineL -lFFL -lquaternion -lABA -lABA_calcKineE -lABA_calcattfrc -lABA_integ -lABA_update -lABA_prepass -lABA_backpass -lABA_mainpass -lABA_pick_data -lABA_set_frc -lABA_set_trans -lABA_set_imat -lABA_set_tmat -lABA_set_lref -lDCA -lDCAs_trans -lDCAs_assembletree -lDCA_backpass  -lDCA_mainpass -lDCA_pick_data  -lDCA_set_imat -lDCA_set_tmat -lDCA_set_lref -lGOLM_check -lGOLMAA_dbasin  -lGOLMAA_set -lGOLMAA -lGOLM_check -lGOLM_set -lGOLM -lPTL -lRMSF -lIOV2  -lLEASQDV -lreadAOUT  -lENECON  -lPROTOPO -lKE  -lCHCOLLA -lCA -lCB -lNC -lPDB -lSBFF -lSTRING  -lnetcdf_mine -lECEPE -lMBAR  -lFASYS -lCSI -lbestfit -lTOPO -lBOXMULL  -lPMF -lREMD -lMB -lPT -lPDBf -lSSL -lFF -lEF -lINTG -lQUA -lLA -lHIST  -lmymath -lIO -ldPCA -lPCA  -lSPE -lMC  -lrmsd -lSHELL -lRAND  -lRBD -llapack_LINUX_kusu -lblas_LINUX_kusu  -lf2c_kusu -lfftw3  -lm -lgc -lnetcdf  -lglib-2.0 ; \
	cp $@ ~/mybin/; \

.c.o:
	$(CC) -c -L${HOME}/lib -I${HOME}/include  -L${HOME}/mylib -I${HOME}/myinclude -Xlinker -rpath -Xlinker ${HOME}/mylib -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0 -I${HOME}/lib/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -L/lib64 $<;

ABA_backpass.o:ABA.h
ABA.o:ABA.h
ABA_mainpass.o:ABA.h
ABA_pick_data.o:ABA.h
ABA_prepass.o:ABA.h
ABA_set_frc.o:ABA.h
ABA_set_imat.o:ABA.h
ABA_set_lref.o:ABA.h
ABA_set_tmat.o:ABA.h
ABA_set_trans.o:ABA.h

clean:
	rm *.o;
