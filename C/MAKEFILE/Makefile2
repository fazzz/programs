#
# Makefile
#

# definitions
#CC	= gcc
#CC	 = gcc -lm
CC      = gcc -g
OBJS	= VisualMD.o LookUpData.o PickData.o Hash.o Usage.o MatchResName.o CreateData.o CreateInPut.o emalloc.c
LIBS	= Vis_MD.h
TARGET	= Vis_mass_MD.exe

# rules of generations

$(TARGET):  $(OBJS) $(LIBS)
	$(CC) -lm -o $@ $(OBJS); \
#	cp $@ ../../../bin/; \

.c.o:
	$(CC) -c $<;

CreateData.o: Vis_MD.h
CreateInPut.o: Vis_MD.h
Hash.o: Vis_MD.h
LookUpData.o: Vis_MD.h
MatchResName.o: Vis_MD.h
PickData.o: Vis_MD.h
Usage.o: Vis_MD.h
VisualMD.o: Vis_MD.h
emalloc.o: Vis_MD.h

clean:
	rm *.o; \
	cd $bin/; \
	rm Vis_mass_MD.exe;
