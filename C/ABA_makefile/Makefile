#
# Makefile
#

# definitions
CC      = gcc -g
OBJS	= mainABAMD_NH_new.o ABA_set_lref.o ABA_set_imat.o ABA_set_tmat.o ABA_pick_data.o ABA_set_trans.o ABA_set_frc.o ABA_mainpass.o ABA_backpass.o ABA_update.o ABA_prepass.o ABA_integ.o ABA_calcattfrc.o ABA_calcKineE.o ABA_Nose-Hoover_new.o ABA_Nose-Hoover.o ABA.o ABA_set_rst.o LA.o EF.o MB.o PTL.o FFL.o TOPO.o quaternion.o netcdf_mineL.o

LIBS	= ABA.h
TARGET	= ABAMD_NH

# rules of generations

$(TARGET):  $(OBJS) $(LIBS)
	$(CC) -o $@ $(OBJS) -L${HOME}/lib -I${HOME}/include  -L${HOME}/mylib -I${HOME}/myinclude -Xlinker -rpath -Xlinker ${HOME}/mylib -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0 -I${HOME}/lib/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -lstringOTF -llapack_LINUX_kusu -lblas_LINUX_kusu  -lf2c_kusu -lfftw3  -lm -lgc -lnetcdf  -lglib-2.0 ; \
#	cp $@ ~/mybin/; \

.c.o:
	$(CC) -c -L${HOME}/lib -I${HOME}/include  -L${HOME}/mylib -I${HOME}/myinclude -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0 -I${HOME}/lib/glib-2.0/include  $<;

EF.o:EF.h
LA.o:LA.h
MB.o:MB.h
FFL.o:FFL.h
PTL.o:PTL.h
TOPO.o:TOPO.h
mymath.o:mymath.h
netcdf_mineL.o:netcdf_mineL.h
quaternion.o:quaternion.h

ABA_mainpass.o:ABA.h
ABA_backpass.o:ABA.h
ABA_pick_data.o:ABA.h
ABA_prepass.o:ABA.h
ABA_set_frc.o:ABA.h
ABA_set_imat.o:ABA.h
ABA_set_lref.o:ABA.h
ABA_set_tmat.o:ABA.h
ABA_set_trans.o:ABA.h
ABA_set_rst.o:ABA.h
ABA.o:ABA.h
ABA_Nose-Hoover.o:ABA.h
ABA_Nose-Hoover_new.o:ABA.h
mainABAMD_NH.o:ABA.h
mainABAMD_NH_new.o:ABA.h

clean:
	rm *.o;
