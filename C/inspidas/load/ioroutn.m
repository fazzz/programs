CGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOG
C
C  *** IOROUTN  ***
C      THIS FILE INCLUDES I/O ROUTINES.
C
C      RDANGL INSTALLED BY AKIO KITAO, APR. 27, 1990
C      WTANGL INSTALLED BY AKIO KITAO, APR. 27, 1990
C
CGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOGOG
C
C     THIS SUBROUTINE READS DIHEDRAL ANGLE DATA FROM ANG.DATA FILE
C     ( ECEPP FORMAT)
C
      SUBROUTINE RDANGL(IREAD)
C
C     *** DECLARE ***
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (PI=3.141 592 653 589 793D0,RAD=PI/180.D0)
C
CMSP  INCLUDE (MAXSIZE)
      PARAMETER(MAXATM=2050,MAXVAR=800,MAXRES=150,MAXSS=10,
     >          MXRS1=MAXRES+1,MXRS10=10*MAXRES,MAXLEN=6000,
     >          MAXDBL=1000,MAXDP=700,MAXDPL=MAXATM+MAXDBL/2,
     >          MAXDP2=MAXDP*2,MXPAIR=MAXVAR*(MAXVAR+1)/2,
     >          MXPAR1=(MAXVAR+1)*(MAXVAR+2)/2,MAXPER=16000,
     >          MAXCNS=900,MAXJ=2*MAXVAR,MAXSS2=MAXSS*2,
     >          MAXINT=40000,MAXI1=MAXINT/2,MAXI2=MAXINT-MAXI1,
     >          MXATM1=MAXATM+1,MVAR1=MAXVAR+1,
     >          MWORK=100000,MIWORK=100000,MCHWOR=100000)
C
C  MAXAT:  MAXIMUM NUMBER OF ATOMS IN A PROTEIN
C  MAXVAR: MAXMUM NUMBER OF DIHEDRAL ANGLE VARIABLES
C  MAXRS:  MAXMUM NUMBER OF RESIDUES
C  MAXSS:  MAXMUM NUMBER OF S-S CROSSLINK
C  MXRS1:  MAXRS+1
C  MXRS10: MAXRS*10
C  MAXLEN: ARRAY SIZE OF INTER IN PREP IN PRECEP ?
C  MAXDBL: MAXMUM NUMBER OF DIPOLES (?)
C  MWORK:  MAXIMUM SIZE OF REAL WORK ARRAY
C  MIWORK: MAXIMUM SIZE OF INTEGER WORK ARRAY
C  MCWORK: MAXIMUM SIZE OF CHARACTER WORK ARRAY
CMSP  INCLUDE (SIZE)
C...Translated by FPP 6.0 (3.06G3) 02/26/96  10:07:07   -dc
C SIZE
C  THIS COMMON BLOCK STORES THE NUMBERS RELATED TO PROTEIN SIZE
C
      COMMON/SIZE/NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      COMMON/SIZE1/NATM1,NVAR1,NPAIR1
C
C  NUMATM:   NUMBER OF ATOMS IN A PROTEIN
C  NUMVAR:   NUMBER OF (DIHEDRAL) ANGLE VARIABLES
C  NN2:      NUMVAR*(NUMVAR+1)/2
C  NUMINT:     ?
C  NSS:      NUMBER OF S-S CROSSLINK (?)
C  NUMRES:   NUMBER OF RESIDUES
C  NATM1=NUMATM+1
C  NVAR1=NUMVAR+1
C  NPAIR1=(NUMVAR+1)*(NUMVAR+2)/2
CMSP  INCLUDE (COMMAND)
C COMMAND.CMN
C
C PURPOSE : COMMON BLOCKS FOR PARSER
C
      COMMON/COMMAN /COMLYN(2000)
      COMMON/COMMSI  /COMLEN,MXCMSZ
      CHARACTER*1 COMLYN
      INTEGER COMLEN
      EQUIVALENCE (COMLYN(1),COMLY2 )
      CHARACTER*2000 COMLY2
C
C     COMLYN  - STORES CHARACTERS IN INPUT (COMMAND) LINE, UP TO 2000
C               CHARACTERS
C     COMLEN  - ACTUAL NUMBER OF NON-BLANK CHARACTERS
C     MXCMSZ  - MAX NUMBER OF CHARACTERS = 2000
C     COMLYN2 - A CHARACTER*2000 REPRESENTATION OF THE COMMAND LINE
C
CMSP  INCLUDE (UNITNUMS)
      COMMON/UNITNU  /IN1,IN2,IN3,IN4,IFN4,
     *               IOUT1,IOUT2,IOUT3,
     * IN01,IOUT02,IOUT06,
C
C PREIN
C
     * IN11,IN13,IOUT12,IOUT14,IOUT16,IN15,IN17,IOUT18,IN19,IOUT20,
C
C PRECEP
C
     * IN21,IN23,IN25,IN27,IN29,
     * IOUT22,IOUT24,IOUT26,IOUT28,IOUT30,
C
C MINIMIZATION
C
     * IN31,IN33,IN35,IN37,IN39,
     * IOUT32,IOUT34,IOUT36,IOUT38,IOUT39,IOUT40,
C
C NORMAL MODE
C
     * IN41,IN43,IN45,IN47,IN49,
     * IOUT42,IOUT44,IOUT46,IOUT48,IOUT50,IOUT52,IOUT54,IOUT56,IOUT58,
     * IOUT60,IOUT62,IOUT64,IOUT66,
C
C MONTE CARLO
C
     * IN71,IN73,IN75,IN77,IN79,
     * IOUT72,IOUT74,IOUT76,IOUT78
C
C PRIMARY INPUT FILE ----- IN01
C PRIMARY OUTPUT FILE ---- IOUT02
C PRIMARY LOG FILE ------- IOUT06
C
C ODD NUMBERS ARE RESERVED FOR INPUT,
C EVEN NUMBERS ARE RESERVED FOR OUTPUT
C FOR EXPLICIT DEFINITIONS OF INPUT AND OUTPUT FILES, PLEASE SEE MAIN
CMSP  INCLUDE (VARANG)
      COMMON/VARANG/ANGRAD(10,MAXRES),VAR(MAXVAR),INDXV(MAXVAR),
     >              ANGLES(10,MAXRES)
C
C  ANGRAD :  DIHEDRAL ANGLES IN ECCEP ORDER (RADIAN)
C  VAR:    DIHEDRAL ANGLES IN WAKO ORDER (RADIAN)
C  INDXV : INDVX=NNN*100+MMM (MMM:1-10)
C            NNN: RESIDUE NUMBER
C            MMM: DIHEDRAL ANGLE NUMBER IN THE RESIDUE(NNN) IN ECEPP
C                 ORDER
C  ANGLES: DIHEDRAL ANGLES IN ECEPP ORDER (DEGREE)
C          USED AS THE WORK]]
C
C     *** READ DIHEDRAL ANGLE DATA FROM ANG.DATA FILE ***
C
      DO 10 I=1,NUMRES
        READ(IREAD,1100)(ANGLES(J,I),J=1,10)
 1100   FORMAT(10F8.3)
   10 CONTINUE
C
CDIR@ IVDEP
      DO I = 1, 10*NUMRES
         ANGRAD(I,1) = ANGLES(I,1)*RAD
      END DO
C
      RETURN
      END
C
C     THIS SUBROUTINE WRITE DIHEDRAL ANGLE DATA
C     ( ECEPP FORMAT)
C
      SUBROUTINE WTANGL(IWRITE)
C
C     *** DECLARE ***
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
CMSP  INCLUDE (MAXSIZE)
      PARAMETER(MAXATM=2050,MAXVAR=800,MAXRES=150,MAXSS=10,
     >          MXRS1=MAXRES+1,MXRS10=10*MAXRES,MAXLEN=6000,
     >          MAXDBL=1000,MAXDP=700,MAXDPL=MAXATM+MAXDBL/2,
     >          MAXDP2=MAXDP*2,MXPAIR=MAXVAR*(MAXVAR+1)/2,
     >          MXPAR1=(MAXVAR+1)*(MAXVAR+2)/2,MAXPER=16000,
     >          MAXCNS=900,MAXJ=2*MAXVAR,MAXSS2=MAXSS*2,
     >          MAXINT=40000,MAXI1=MAXINT/2,MAXI2=MAXINT-MAXI1,
     >          MXATM1=MAXATM+1,MVAR1=MAXVAR+1,
     >          MWORK=100000,MIWORK=100000,MCHWOR=100000)
C
C  MAXAT:  MAXIMUM NUMBER OF ATOMS IN A PROTEIN
C  MAXVAR: MAXMUM NUMBER OF DIHEDRAL ANGLE VARIABLES
C  MAXRS:  MAXMUM NUMBER OF RESIDUES
C  MAXSS:  MAXMUM NUMBER OF S-S CROSSLINK
C  MXRS1:  MAXRS+1
C  MXRS10: MAXRS*10
C  MAXLEN: ARRAY SIZE OF INTER IN PREP IN PRECEP ?
C  MAXDBL: MAXMUM NUMBER OF DIPOLES (?)
C  MWORK:  MAXIMUM SIZE OF REAL WORK ARRAY
C  MIWORK: MAXIMUM SIZE OF INTEGER WORK ARRAY
C  MCWORK: MAXIMUM SIZE OF CHARACTER WORK ARRAY
CMSP  INCLUDE (SIZE)
C...Translated by FPP 6.0 (3.06G3) 02/26/96  10:07:07   -dc
C SIZE
C  THIS COMMON BLOCK STORES THE NUMBERS RELATED TO PROTEIN SIZE
C
      COMMON/SIZE/NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      COMMON/SIZE1/NATM1,NVAR1,NPAIR1
C
C  NUMATM:   NUMBER OF ATOMS IN A PROTEIN
C  NUMVAR:   NUMBER OF (DIHEDRAL) ANGLE VARIABLES
C  NN2:      NUMVAR*(NUMVAR+1)/2
C  NUMINT:     ?
C  NSS:      NUMBER OF S-S CROSSLINK (?)
C  NUMRES:   NUMBER OF RESIDUES
C  NATM1=NUMATM+1
C  NVAR1=NUMVAR+1
C  NPAIR1=(NUMVAR+1)*(NUMVAR+2)/2
CMSP  INCLUDE (COMMAND)
C COMMAND.CMN
C
C PURPOSE : COMMON BLOCKS FOR PARSER
C
      COMMON/COMMAN /COMLYN(2000)
      COMMON/COMMSI  /COMLEN,MXCMSZ
      CHARACTER*1 COMLYN
      INTEGER COMLEN
      EQUIVALENCE (COMLYN(1),COMLY2 )
      CHARACTER*2000 COMLY2
C
C     COMLYN  - STORES CHARACTERS IN INPUT (COMMAND) LINE, UP TO 2000
C               CHARACTERS
C     COMLEN  - ACTUAL NUMBER OF NON-BLANK CHARACTERS
C     MXCMSZ  - MAX NUMBER OF CHARACTERS = 2000
C     COMLYN2 - A CHARACTER*2000 REPRESENTATION OF THE COMMAND LINE
C
CMSP  INCLUDE (UNITNUMS)
      COMMON/UNITNU  /IN1,IN2,IN3,IN4,IFN4,
     *               IOUT1,IOUT2,IOUT3,
     * IN01,IOUT02,IOUT06,
C
C PREIN
C
     * IN11,IN13,IOUT12,IOUT14,IOUT16,IN15,IN17,IOUT18,IN19,IOUT20,
C
C PRECEP
C
     * IN21,IN23,IN25,IN27,IN29,
     * IOUT22,IOUT24,IOUT26,IOUT28,IOUT30,
C
C MINIMIZATION
C
     * IN31,IN33,IN35,IN37,IN39,
     * IOUT32,IOUT34,IOUT36,IOUT38,IOUT39,IOUT40,
C
C NORMAL MODE
C
     * IN41,IN43,IN45,IN47,IN49,
     * IOUT42,IOUT44,IOUT46,IOUT48,IOUT50,IOUT52,IOUT54,IOUT56,IOUT58,
     * IOUT60,IOUT62,IOUT64,IOUT66,
C
C MONTE CARLO
C
     * IN71,IN73,IN75,IN77,IN79,
     * IOUT72,IOUT74,IOUT76,IOUT78
C
C PRIMARY INPUT FILE ----- IN01
C PRIMARY OUTPUT FILE ---- IOUT02
C PRIMARY LOG FILE ------- IOUT06
C
C ODD NUMBERS ARE RESERVED FOR INPUT,
C EVEN NUMBERS ARE RESERVED FOR OUTPUT
C FOR EXPLICIT DEFINITIONS OF INPUT AND OUTPUT FILES, PLEASE SEE MAIN
CMSP  INCLUDE (VARANG)
      COMMON/VARANG/ANGRAD(10,MAXRES),VAR(MAXVAR),INDXV(MAXVAR),
     >              ANGLES(10,MAXRES)
C
C  ANGRAD :  DIHEDRAL ANGLES IN ECCEP ORDER (RADIAN)
C  VAR:    DIHEDRAL ANGLES IN WAKO ORDER (RADIAN)
C  INDXV : INDVX=NNN*100+MMM (MMM:1-10)
C            NNN: RESIDUE NUMBER
C            MMM: DIHEDRAL ANGLE NUMBER IN THE RESIDUE(NNN) IN ECEPP
C                 ORDER
C  ANGLES: DIHEDRAL ANGLES IN ECEPP ORDER (DEGREE)
C          USED AS THE WORK]]
C
      PARAMETER (PI=3.141 592 653 589 793D0,RAD=PI/180.D0)
C
C     *** WRITE DIHEDRAL ANGLE DATA
C
CDIR@ IVDEP
      DO I = 1, 10*NUMRES
         ANGLES(I,1) = ANGRAD(I,1)/RAD
      END DO
C
      DO 20 I=1,NUMRES
        WRITE(IWRITE,1100)(ANGLES(J,I),J=1,10)
 1100   FORMAT(10F8.3)
   20 CONTINUE
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  GTANGL (IREAD,II)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
CMSP  INCLUDE (MAXSIZE)
C MAXSIZE
C     MAXSIZE DEFINES MAXIMUM ARRAY SIZE
C
      PARAMETER(MAXATM=2050,MAXVAR=800,MAXRES=150,MAXSS=10,
     >          MXRS1=MAXRES+1,MXRS10=10*MAXRES,MAXLEN=6000,
     >          MAXDBL=1000,MAXDP=700,MAXDPL=MAXATM+MAXDBL/2,
     >          MAXDP2=MAXDP*2,MXPAIR=MAXVAR*(MAXVAR+1)/2,
     >          MXPAR1=(MAXVAR+1)*(MAXVAR+2)/2,MAXPER=16000,
     >          MAXCNS=900,MAXJ=2*MAXVAR,MAXSS2=MAXSS*2,
     >          MAXINT=40000,MAXI1=MAXINT/2,MAXI2=MAXINT-MAXI1,
     >          MXATM1=MAXATM+1,MVAR1=MAXVAR+1,
     >          MWORK=100000,MIWORK=100000,MCHWOR=100000)
C
C  MAXAT:  MAXIMUM NUMBER OF ATOMS IN A PROTEIN
C  MAXVAR: MAXMUM NUMBER OF DIHEDRAL ANGLE VARIABLES
C  MAXRS:  MAXMUM NUMBER OF RESIDUES
C  MAXSS:  MAXMUM NUMBER OF S-S CROSSLINK
C  MXRS1:  MAXRS+1
C  MXRS10: MAXRS*10
C  MAXLEN: ARRAY SIZE OF INTER IN PREP IN PRECEP ?
C  MAXDBL: MAXMUM NUMBER OF DIPOLES (?)
C  MWORK:  MAXIMUM SIZE OF REAL WORK ARRAY
C  MIWORK: MAXIMUM SIZE OF INTEGER WORK ARRAY
C  MCWORK: MAXIMUM SIZE OF CHARACTER WORK ARRAY
CMSP  INCLUDE (SIZE)
C SIZE
C  THIS COMMON BLOCK STORES THE NUMBERS RELATED TO PROTEIN SIZE
C
      COMMON/SIZE/NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      COMMON/SIZE1/NATM1,NVAR1,NPAIR1
C
C  NUMATM:   NUMBER OF ATOMS IN A PROTEIN
C  NUMVAR:   NUMBER OF (DIHEDRAL) ANGLE VARIABLES
C  NN2:      NUMVAR*(NUMVAR+1)/2
C  NUMINT:     ?
C  NSS:      NUMBER OF S-S CROSSLINK (?)
C  NUMRES:   NUMBER OF RESIDUES
C  NATM1=NUMATM+1
C  NVAR1=NUMVAR+1
C  NPAIR1=(NUMVAR+1)*(NUMVAR+2)/2
CMSP  INCLUDE (VARANG)
C VARANG
C
      COMMON/VARANG/ANGRAD(10,MAXRES),VAR(MAXVAR),INDXV(MAXVAR),
     >              ANGLES(10,MAXRES)
C
C  ANGRAD :  DIHEDRAL ANGLES IN ECCEP ORDER (RADIAN)
C  VAR:    DIHEDRAL ANGLES IN WAKO ORDER (RADIAN)
C  INDXV : INDVX=NNN*100+MMM (MMM:1-10)
C            NNN: RESIDUE NUMBER
C            MMM: DIHEDRAL ANGLE NUMBER IN THE RESIDUE(NNN) IN ECEPP
C                 ORDER
C  ANGLES: DIHEDRAL ANGLES IN ECEPP ORDER (DEGREE)
C          USED AS THE WORK]]
CMSP  INCLUDE (UNITNUMS)
C  UNITNUMS.CMN
C
C  THIS FILE CONTAINS THE  COMMON BLOCK NECESSARY
C  TO ASSIGN ALL I/0 UNIT NUMBERS. A SUMMARY OF ALL I/O FILES USED BY
C  IMPACT APPEARS BELOW, ORGANIZED BY TASK.
C
C
      COMMON/UNITNU  /IN1,IN2,IN3,IN4,IFN4,
     *               IOUT1,IOUT2,IOUT3,
     * IN01,IOUT02,IOUT06,
C
C PREIN
C
     * IN11,IN13,IOUT12,IOUT14,IOUT16,IN15,IN17,IOUT18,IN19,IOUT20,
C
C PRECEP
C
     * IN21,IN23,IN25,IN27,IN29,
     * IOUT22,IOUT24,IOUT26,IOUT28,IOUT30,
C
C MINIMIZATION
C
     * IN31,IN33,IN35,IN37,IN39,
     * IOUT32,IOUT34,IOUT36,IOUT38,IOUT39,IOUT40,
C
C NORMAL MODE
C
     * IN41,IN43,IN45,IN47,IN49,
     * IOUT42,IOUT44,IOUT46,IOUT48,IOUT50,IOUT52,IOUT54,IOUT56,IOUT58,
     * IOUT60,IOUT62,IOUT64,IOUT66,
C
C MONTE CARLO
C
     * IN71,IN73,IN75,IN77,IN79,
     * IOUT72,IOUT74,IOUT76,IOUT78
C
C PRIMARY INPUT FILE ----- IN01
C PRIMARY OUTPUT FILE ---- IOUT02
C PRIMARY LOG FILE ------- IOUT06
C
C ODD NUMBERS ARE RESERVED FOR INPUT,
C EVEN NUMBERS ARE RESERVED FOR OUTPUT
C FOR EXPLICIT DEFINITIONS OF INPUT AND OUTPUT FILES, PLEASE SEE MAIN
C     COMMON/NUMBER/  NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
C     DIMENSION  VAR(NUMVAR)
C  FOR FRESH START THE INPUT FILE CONTAINS DIHEDRAL ANGLES ONLY.
C  FOR RESTART THE INPUT FILE CONTAINS DIHEDRAL ANGLES AND THE
C  ITERATION NUMBER.
C
      II=1
      READ(IREAD)  (VAR(I),I=1,NUMVAR)
      READ(IREAD,END=10)  II
        WRITE(IOUT06,*) 'READ FROM A RESTART FILE'
        WRITE(IOUT02,*) 'READ FROM A RESTART FILE'
        II = II + 1
  10  CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  PTANGL (IWRITE,VAR,II)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
CMSP  INCLUDE (MAXSIZE)
C MAXSIZE
C     MAXSIZE DEFINES MAXIMUM ARRAY SIZE
C
      PARAMETER(MAXATM=2050,MAXVAR=800,MAXRES=150,MAXSS=10,
     >          MXRS1=MAXRES+1,MXRS10=10*MAXRES,MAXLEN=6000,
     >          MAXDBL=1000,MAXDP=700,MAXDPL=MAXATM+MAXDBL/2,
     >          MAXDP2=MAXDP*2,MXPAIR=MAXVAR*(MAXVAR+1)/2,
     >          MXPAR1=(MAXVAR+1)*(MAXVAR+2)/2,MAXPER=16000,
     >          MAXCNS=900,MAXJ=2*MAXVAR,MAXSS2=MAXSS*2,
     >          MAXINT=40000,MAXI1=MAXINT/2,MAXI2=MAXINT-MAXI1,
     >          MXATM1=MAXATM+1,MVAR1=MAXVAR+1,
     >          MWORK=100000,MIWORK=100000,MCHWOR=100000)
C
C  MAXAT:  MAXIMUM NUMBER OF ATOMS IN A PROTEIN
C  MAXVAR: MAXMUM NUMBER OF DIHEDRAL ANGLE VARIABLES
C  MAXRS:  MAXMUM NUMBER OF RESIDUES
C  MAXSS:  MAXMUM NUMBER OF S-S CROSSLINK
C  MXRS1:  MAXRS+1
C  MXRS10: MAXRS*10
C  MAXLEN: ARRAY SIZE OF INTER IN PREP IN PRECEP ?
C  MAXDBL: MAXMUM NUMBER OF DIPOLES (?)
C  MWORK:  MAXIMUM SIZE OF REAL WORK ARRAY
C  MIWORK: MAXIMUM SIZE OF INTEGER WORK ARRAY
C  MCWORK: MAXIMUM SIZE OF CHARACTER WORK ARRAY
CMSP  INCLUDE (SIZE)
C SIZE
C  THIS COMMON BLOCK STORES THE NUMBERS RELATED TO PROTEIN SIZE
C
      COMMON/SIZE/NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      COMMON/SIZE1/NATM1,NVAR1,NPAIR1
C
C  NUMATM:   NUMBER OF ATOMS IN A PROTEIN
C  NUMVAR:   NUMBER OF (DIHEDRAL) ANGLE VARIABLES
C  NN2:      NUMVAR*(NUMVAR+1)/2
C  NUMINT:     ?
C  NSS:      NUMBER OF S-S CROSSLINK (?)
C  NUMRES:   NUMBER OF RESIDUES
C  NATM1=NUMATM+1
C  NVAR1=NUMVAR+1
C  NPAIR1=(NUMVAR+1)*(NUMVAR+2)/2
C     COMMON/NUMBER/  NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      DIMENSION  VAR(MAXVAR)
      WRITE(IWRITE) (VAR(I),I=1,NUMVAR)
      WRITE(IWRITE) II
      REWIND  IWRITE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  PTCOOR (IWRITE)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
CMSP  INCLUDE (MAXSIZE)
C MAXSIZE
C     MAXSIZE DEFINES MAXIMUM ARRAY SIZE
C
      PARAMETER(MAXATM=2050,MAXVAR=800,MAXRES=150,MAXSS=10,
     >          MXRS1=MAXRES+1,MXRS10=10*MAXRES,MAXLEN=6000,
     >          MAXDBL=1000,MAXDP=700,MAXDPL=MAXATM+MAXDBL/2,
     >          MAXDP2=MAXDP*2,MXPAIR=MAXVAR*(MAXVAR+1)/2,
     >          MXPAR1=(MAXVAR+1)*(MAXVAR+2)/2,MAXPER=16000,
     >          MAXCNS=900,MAXJ=2*MAXVAR,MAXSS2=MAXSS*2,
     >          MAXINT=40000,MAXI1=MAXINT/2,MAXI2=MAXINT-MAXI1,
     >          MXATM1=MAXATM+1,MVAR1=MAXVAR+1,
     >          MWORK=100000,MIWORK=100000,MCHWOR=100000)
C
C  MAXAT:  MAXIMUM NUMBER OF ATOMS IN A PROTEIN
C  MAXVAR: MAXMUM NUMBER OF DIHEDRAL ANGLE VARIABLES
C  MAXRS:  MAXMUM NUMBER OF RESIDUES
C  MAXSS:  MAXMUM NUMBER OF S-S CROSSLINK
C  MXRS1:  MAXRS+1
C  MXRS10: MAXRS*10
C  MAXLEN: ARRAY SIZE OF INTER IN PREP IN PRECEP ?
C  MAXDBL: MAXMUM NUMBER OF DIPOLES (?)
C  MWORK:  MAXIMUM SIZE OF REAL WORK ARRAY
C  MIWORK: MAXIMUM SIZE OF INTEGER WORK ARRAY
C  MCWORK: MAXIMUM SIZE OF CHARACTER WORK ARRAY
CMSP  INCLUDE (SIZE)
C SIZE
C  THIS COMMON BLOCK STORES THE NUMBERS RELATED TO PROTEIN SIZE
C
      COMMON/SIZE/NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      COMMON/SIZE1/NATM1,NVAR1,NPAIR1
C
C  NUMATM:   NUMBER OF ATOMS IN A PROTEIN
C  NUMVAR:   NUMBER OF (DIHEDRAL) ANGLE VARIABLES
C  NN2:      NUMVAR*(NUMVAR+1)/2
C  NUMINT:     ?
C  NSS:      NUMBER OF S-S CROSSLINK (?)
C  NUMRES:   NUMBER OF RESIDUES
C  NATM1=NUMATM+1
C  NVAR1=NUMVAR+1
C  NPAIR1=(NUMVAR+1)*(NUMVAR+2)/2
CMSP  INCLUDE (COORD)
C COORD
C
      COMMON/COORD/CO(3,MAXATM)
C
C  CO: VARIABLE COORDINATES (ECEPP ORDER)
C
      WRITE(IWRITE) ((CO(J,I),J=1,3),I=1,NUMATM)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RDSTD(FL,LEVEL,IDENT,ATTRIB,NR,RECS)
C-----------------------------------------------------------------------
      INTEGER FL, LEVEL, NR
      CHARACTER IDENT*8, ATTRIB*8, RECS*(*)
      COMMON /BUFFER/BUF
      CHARACTER*80 BUF
      CALL RDSTDM(FL,BUF,LEVEL,IDENT,ATTRIB,NR,RECS)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RDSTDM(FL,BUF,LEVEL,IDENT,ATTRIB,NR,RECS)
C-----------------------------------------------------------------------
      INTEGER FL, LEVEL, NR
      CHARACTER IDENT*8, ATTRIB*8, RECS*(*), BUF*80
      CHARACTER*8 LEV0,LEV1,LEV2,LEV3,LEV4,LEV5,ENDL
      DATA LEV0/'0       '/, LEV1/'1       '/, LEV2/'2       '/, LEV3/
     &'3       '/, LEV4/'4       '/, LEV5/'5       '/, ENDL/'END     '/
      IF(.NOT.(NR .EQ. 0))GOTO  3000
      READ(FL,'(A80)',END=99) BUF
 3000 CONTINUE
      IF(.NOT.(BUF(1:8) .EQ. ENDL))GOTO  3002
      LEVEL = -1
      RETURN
 3002 CONTINUE
      I = 1
 3004 CONTINUE
      IF(.NOT.(I+79 .LE. LEN(RECS)))GOTO  3007
      RECS(I:I+79) = BUF
      I = I + 80
 3007 CONTINUE
      READ (FL,'(A80)',END=99) BUF
 3005 IF(.NOT.(BUF (1:8) .EQ. LEV0 .OR.BUF (1:8) .EQ. LEV1 .OR.BUF (1:8)
     & .EQ. LEV2 .OR.BUF (1:8) .EQ. LEV3 .OR.BUF (1:8) .EQ. LEV4 .OR.
     &BUF (1:8) .EQ. LEV5 .OR.BUF (1:8) .EQ. ENDL ))GOTO  3004
      IF(.NOT.(RECS(1:8) .EQ. LEV0))GOTO  3009
      LEVEL = 0
 3009 CONTINUE
      IF(.NOT.(RECS(1:8) .EQ. LEV1))GOTO  3011
      LEVEL = 1
 3011 CONTINUE
      IF(.NOT.(RECS(1:8) .EQ. LEV2))GOTO  3013
      LEVEL = 2
 3013 CONTINUE
      IF(.NOT.(RECS(1:8) .EQ. LEV3))GOTO  3015
      LEVEL = 3
 3015 CONTINUE
      IF(.NOT.(RECS(1:8) .EQ. LEV4))GOTO  3017
      LEVEL = 4
 3017 CONTINUE
      IF(.NOT.(RECS(1:8) .EQ. LEV5))GOTO  3019
      LEVEL = 5
 3019 CONTINUE
      IDENT = RECS(9:16)
      ATTRIB = RECS(17:24)
      NR = (I-1)/80
      RETURN
 99    LEVEL = -1
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTINFO(KEY,NI,INFORM,NR,RECS)
C-----------------------------------------------------------------------
      INTEGER NI, NR
      CHARACTER KEY*8, INFORM*(*), RECS*(*)
      CHARACTER*8 CONT,SPCE,LEV0,LEV1,LEV2,LEV3,LEV4,LEV5
      DATA CONT/'-       '/, SPCE/'        '/, LEV0/'0       '/, LEV1/
     &'1       '/, LEV2/'2       '/, LEV3/'3       '/, LEV4/'4       '/,
     & LEV5/'5       '/
      NI = 0
      I=1
 3000 IF(.NOT.(I.LE.NR))GOTO  3002
      IB = 80*(I-1)
      IF(.NOT.(RECS(IB+1:IB+8 ) .EQ. SPCE .OR.RECS(IB+1:IB+8 ) .EQ.
     &LEV0 .OR.RECS(IB+1:IB+8 ) .EQ. LEV1 .OR.RECS(IB+1:IB+8 ) .EQ.
     &LEV2 .OR.RECS(IB+1:IB+8 ) .EQ. LEV3 .OR.RECS(IB+1:IB+8 ) .EQ.
     &LEV4 .OR.RECS(IB+1:IB+8 ) .EQ. LEV5 ))GOTO  3003
      IF(.NOT.(RECS(IB+ 9:IB+16) .EQ. KEY))GOTO  3005
      GO TO 1
 3005 CONTINUE
      IF(.NOT.(RECS(IB+41:IB+48) .EQ. KEY))GOTO  3007
      GO TO 2
 3007 CONTINUE
 3003 CONTINUE
      I=I+1
      GOTO  3000
 3002 CONTINUE
      RETURN
 3009 IF(.NOT.(I.LE.NR))GOTO  3011
      IB = 80*(I-1)
      IF(.NOT.(RECS(IB+9:IB+16) .NE. CONT))GOTO  3012
      GOTO  3011
 3012 CONTINUE
 1     INFORM = INFORM(1:NI) // RECS(IB+17:IB+40)
      NI = NI + 24
      IF(.NOT.(RECS(IB+41:IB+48) .NE. CONT))GOTO  3014
      GOTO  3011
 3014 CONTINUE
 2     INFORM = INFORM(1:NI) // RECS(IB+49:IB+72)
      NI = NI + 24
      I=I+1
      GOTO  3009
 3011 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WTEIG(IOUT)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
CMSP  INCLUDE (MAXSIZE)
C MAXSIZE
C     MAXSIZE DEFINES MAXIMUM ARRAY SIZE
C
      PARAMETER(MAXATM=2050,MAXVAR=800,MAXRES=150,MAXSS=10,
     >          MXRS1=MAXRES+1,MXRS10=10*MAXRES,MAXLEN=6000,
     >          MAXDBL=1000,MAXDP=700,MAXDPL=MAXATM+MAXDBL/2,
     >          MAXDP2=MAXDP*2,MXPAIR=MAXVAR*(MAXVAR+1)/2,
     >          MXPAR1=(MAXVAR+1)*(MAXVAR+2)/2,MAXPER=16000,
     >          MAXCNS=900,MAXJ=2*MAXVAR,MAXSS2=MAXSS*2,
     >          MAXINT=40000,MAXI1=MAXINT/2,MAXI2=MAXINT-MAXI1,
     >          MXATM1=MAXATM+1,MVAR1=MAXVAR+1,
     >          MWORK=100000,MIWORK=100000,MCHWOR=100000)
C
C  MAXAT:  MAXIMUM NUMBER OF ATOMS IN A PROTEIN
C  MAXVAR: MAXMUM NUMBER OF DIHEDRAL ANGLE VARIABLES
C  MAXRS:  MAXMUM NUMBER OF RESIDUES
C  MAXSS:  MAXMUM NUMBER OF S-S CROSSLINK
C  MXRS1:  MAXRS+1
C  MXRS10: MAXRS*10
C  MAXLEN: ARRAY SIZE OF INTER IN PREP IN PRECEP ?
C  MAXDBL: MAXMUM NUMBER OF DIPOLES (?)
C  MWORK:  MAXIMUM SIZE OF REAL WORK ARRAY
C  MIWORK: MAXIMUM SIZE OF INTEGER WORK ARRAY
C  MCWORK: MAXIMUM SIZE OF CHARACTER WORK ARRAY
CMSP  INCLUDE (SIZE)
C SIZE
C  THIS COMMON BLOCK STORES THE NUMBERS RELATED TO PROTEIN SIZE
C
      COMMON/SIZE/NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      COMMON/SIZE1/NATM1,NVAR1,NPAIR1
C
C  NUMATM:   NUMBER OF ATOMS IN A PROTEIN
C  NUMVAR:   NUMBER OF (DIHEDRAL) ANGLE VARIABLES
C  NN2:      NUMVAR*(NUMVAR+1)/2
C  NUMINT:     ?
C  NSS:      NUMBER OF S-S CROSSLINK (?)
C  NUMRES:   NUMBER OF RESIDUES
C  NATM1=NUMATM+1
C  NVAR1=NUMVAR+1
C  NPAIR1=(NUMVAR+1)*(NUMVAR+2)/2
CMSP  INCLUDE (EIGEN)
C EIGEN
        COMMON/EIGEN/ EIGVAL(MAXVAR),EIGVEC(MAXVAR,MAXVAR)
C EIGVAL: EIGEN VALUES
C EIGVEC : EIGEN VECTORS
CMSP  INCLUDE (UNITNUMS)
C  UNITNUMS.CMN
C
C  THIS FILE CONTAINS THE  COMMON BLOCK NECESSARY
C  TO ASSIGN ALL I/0 UNIT NUMBERS. A SUMMARY OF ALL I/O FILES USED BY
C  IMPACT APPEARS BELOW, ORGANIZED BY TASK.
C
C
      COMMON/UNITNU  /IN1,IN2,IN3,IN4,IFN4,
     *               IOUT1,IOUT2,IOUT3,
     * IN01,IOUT02,IOUT06,
C
C PREIN
C
     * IN11,IN13,IOUT12,IOUT14,IOUT16,IN15,IN17,IOUT18,IN19,IOUT20,
C
C PRECEP
C
     * IN21,IN23,IN25,IN27,IN29,
     * IOUT22,IOUT24,IOUT26,IOUT28,IOUT30,
C
C MINIMIZATION
C
     * IN31,IN33,IN35,IN37,IN39,
     * IOUT32,IOUT34,IOUT36,IOUT38,IOUT39,IOUT40,
C
C NORMAL MODE
C
     * IN41,IN43,IN45,IN47,IN49,
     * IOUT42,IOUT44,IOUT46,IOUT48,IOUT50,IOUT52,IOUT54,IOUT56,IOUT58,
     * IOUT60,IOUT62,IOUT64,IOUT66,
C
C MONTE CARLO
C
     * IN71,IN73,IN75,IN77,IN79,
     * IOUT72,IOUT74,IOUT76,IOUT78
C
C PRIMARY INPUT FILE ----- IN01
C PRIMARY OUTPUT FILE ---- IOUT02
C PRIMARY LOG FILE ------- IOUT06
C
C ODD NUMBERS ARE RESERVED FOR INPUT,
C EVEN NUMBERS ARE RESERVED FOR OUTPUT
C FOR EXPLICIT DEFINITIONS OF INPUT AND OUTPUT FILES, PLEASE SEE MAIN
C
      WRITE(IOUT,*) NUMVAR
      WRITE(IOUT,1000) (EIGVAL(I),I=1,NUMVAR)
      WRITE(IOUT,1000) ((EIGVEC(I,J),I=1,NUMVAR),J=1,NUMVAR)
 1000 FORMAT(1X,10E12.4)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE PTEIG(IOUT)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
CMSP  INCLUDE (MAXSIZE)
C MAXSIZE
C     MAXSIZE DEFINES MAXIMUM ARRAY SIZE
C
      PARAMETER(MAXATM=2050,MAXVAR=800,MAXRES=150,MAXSS=10,
     >          MXRS1=MAXRES+1,MXRS10=10*MAXRES,MAXLEN=6000,
     >          MAXDBL=1000,MAXDP=700,MAXDPL=MAXATM+MAXDBL/2,
     >          MAXDP2=MAXDP*2,MXPAIR=MAXVAR*(MAXVAR+1)/2,
     >          MXPAR1=(MAXVAR+1)*(MAXVAR+2)/2,MAXPER=16000,
     >          MAXCNS=900,MAXJ=2*MAXVAR,MAXSS2=MAXSS*2,
     >          MAXINT=40000,MAXI1=MAXINT/2,MAXI2=MAXINT-MAXI1,
     >          MXATM1=MAXATM+1,MVAR1=MAXVAR+1,
     >          MWORK=100000,MIWORK=100000,MCHWOR=100000)
C
C  MAXAT:  MAXIMUM NUMBER OF ATOMS IN A PROTEIN
C  MAXVAR: MAXMUM NUMBER OF DIHEDRAL ANGLE VARIABLES
C  MAXRS:  MAXMUM NUMBER OF RESIDUES
C  MAXSS:  MAXMUM NUMBER OF S-S CROSSLINK
C  MXRS1:  MAXRS+1
C  MXRS10: MAXRS*10
C  MAXLEN: ARRAY SIZE OF INTER IN PREP IN PRECEP ?
C  MAXDBL: MAXMUM NUMBER OF DIPOLES (?)
C  MWORK:  MAXIMUM SIZE OF REAL WORK ARRAY
C  MIWORK: MAXIMUM SIZE OF INTEGER WORK ARRAY
C  MCWORK: MAXIMUM SIZE OF CHARACTER WORK ARRAY
CMSP  INCLUDE (SIZE)
C SIZE
C  THIS COMMON BLOCK STORES THE NUMBERS RELATED TO PROTEIN SIZE
C
      COMMON/SIZE/NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      COMMON/SIZE1/NATM1,NVAR1,NPAIR1
C
C  NUMATM:   NUMBER OF ATOMS IN A PROTEIN
C  NUMVAR:   NUMBER OF (DIHEDRAL) ANGLE VARIABLES
C  NN2:      NUMVAR*(NUMVAR+1)/2
C  NUMINT:     ?
C  NSS:      NUMBER OF S-S CROSSLINK (?)
C  NUMRES:   NUMBER OF RESIDUES
C  NATM1=NUMATM+1
C  NVAR1=NUMVAR+1
C  NPAIR1=(NUMVAR+1)*(NUMVAR+2)/2
CMSP  INCLUDE (EIGEN)
C EIGEN
        COMMON/EIGEN/ EIGVAL(MAXVAR),EIGVEC(MAXVAR,MAXVAR)
C EIGVAL: EIGEN VALUES
C EIGVEC : EIGEN VECTORS
CMSP  INCLUDE (UNITNUMS)
C  UNITNUMS.CMN
C
C  THIS FILE CONTAINS THE  COMMON BLOCK NECESSARY
C  TO ASSIGN ALL I/0 UNIT NUMBERS. A SUMMARY OF ALL I/O FILES USED BY
C  IMPACT APPEARS BELOW, ORGANIZED BY TASK.
C
C
      COMMON/UNITNU  /IN1,IN2,IN3,IN4,IFN4,
     *               IOUT1,IOUT2,IOUT3,
     * IN01,IOUT02,IOUT06,
C
C PREIN
C
     * IN11,IN13,IOUT12,IOUT14,IOUT16,IN15,IN17,IOUT18,IN19,IOUT20,
C
C PRECEP
C
     * IN21,IN23,IN25,IN27,IN29,
     * IOUT22,IOUT24,IOUT26,IOUT28,IOUT30,
C
C MINIMIZATION
C
     * IN31,IN33,IN35,IN37,IN39,
     * IOUT32,IOUT34,IOUT36,IOUT38,IOUT39,IOUT40,
C
C NORMAL MODE
C
     * IN41,IN43,IN45,IN47,IN49,
     * IOUT42,IOUT44,IOUT46,IOUT48,IOUT50,IOUT52,IOUT54,IOUT56,IOUT58,
     * IOUT60,IOUT62,IOUT64,IOUT66,
C
C MONTE CARLO
C
     * IN71,IN73,IN75,IN77,IN79,
     * IOUT72,IOUT74,IOUT76,IOUT78
C
C PRIMARY INPUT FILE ----- IN01
C PRIMARY OUTPUT FILE ---- IOUT02
C PRIMARY LOG FILE ------- IOUT06
C
C ODD NUMBERS ARE RESERVED FOR INPUT,
C EVEN NUMBERS ARE RESERVED FOR OUTPUT
C FOR EXPLICIT DEFINITIONS OF INPUT AND OUTPUT FILES, PLEASE SEE MAIN
C
      WRITE(IOUT) NUMVAR,(EIGVAL(I),I=1,NUMVAR),
     &            ((EIGVEC(I,J),I=1,NUMVAR),J=1,NUMVAR)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTEIG(IIN)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
CMSP  INCLUDE (MAXSIZE)
C MAXSIZE
C     MAXSIZE DEFINES MAXIMUM ARRAY SIZE
C
      PARAMETER(MAXATM=2050,MAXVAR=800,MAXRES=150,MAXSS=10,
     >          MXRS1=MAXRES+1,MXRS10=10*MAXRES,MAXLEN=6000,
     >          MAXDBL=1000,MAXDP=700,MAXDPL=MAXATM+MAXDBL/2,
     >          MAXDP2=MAXDP*2,MXPAIR=MAXVAR*(MAXVAR+1)/2,
     >          MXPAR1=(MAXVAR+1)*(MAXVAR+2)/2,MAXPER=16000,
     >          MAXCNS=900,MAXJ=2*MAXVAR,MAXSS2=MAXSS*2,
     >          MAXINT=40000,MAXI1=MAXINT/2,MAXI2=MAXINT-MAXI1,
     >          MXATM1=MAXATM+1,MVAR1=MAXVAR+1,
     >          MWORK=100000,MIWORK=100000,MCHWOR=100000)
C
C  MAXAT:  MAXIMUM NUMBER OF ATOMS IN A PROTEIN
C  MAXVAR: MAXMUM NUMBER OF DIHEDRAL ANGLE VARIABLES
C  MAXRS:  MAXMUM NUMBER OF RESIDUES
C  MAXSS:  MAXMUM NUMBER OF S-S CROSSLINK
C  MXRS1:  MAXRS+1
C  MXRS10: MAXRS*10
C  MAXLEN: ARRAY SIZE OF INTER IN PREP IN PRECEP ?
C  MAXDBL: MAXMUM NUMBER OF DIPOLES (?)
C  MWORK:  MAXIMUM SIZE OF REAL WORK ARRAY
C  MIWORK: MAXIMUM SIZE OF INTEGER WORK ARRAY
C  MCWORK: MAXIMUM SIZE OF CHARACTER WORK ARRAY
CMSP  INCLUDE (SIZE)
C SIZE
C  THIS COMMON BLOCK STORES THE NUMBERS RELATED TO PROTEIN SIZE
C
      COMMON/SIZE/NUMATM,NUMVAR,NN2,NUMINT,NSS,NUMRES
      COMMON/SIZE1/NATM1,NVAR1,NPAIR1
C
C  NUMATM:   NUMBER OF ATOMS IN A PROTEIN
C  NUMVAR:   NUMBER OF (DIHEDRAL) ANGLE VARIABLES
C  NN2:      NUMVAR*(NUMVAR+1)/2
C  NUMINT:     ?
C  NSS:      NUMBER OF S-S CROSSLINK (?)
C  NUMRES:   NUMBER OF RESIDUES
C  NATM1=NUMATM+1
C  NVAR1=NUMVAR+1
C  NPAIR1=(NUMVAR+1)*(NUMVAR+2)/2
CMSP  INCLUDE (EIGEN)
C EIGEN
        COMMON/EIGEN/ EIGVAL(MAXVAR),EIGVEC(MAXVAR,MAXVAR)
C EIGVAL: EIGEN VALUES
C EIGVEC : EIGEN VECTORS
CMSP  INCLUDE (UNITNUMS)
C  UNITNUMS.CMN
C
C  THIS FILE CONTAINS THE  COMMON BLOCK NECESSARY
C  TO ASSIGN ALL I/0 UNIT NUMBERS. A SUMMARY OF ALL I/O FILES USED BY
C  IMPACT APPEARS BELOW, ORGANIZED BY TASK.
C
C
      COMMON/UNITNU  /IN1,IN2,IN3,IN4,IFN4,
     *               IOUT1,IOUT2,IOUT3,
     * IN01,IOUT02,IOUT06,
C
C PREIN
C
     * IN11,IN13,IOUT12,IOUT14,IOUT16,IN15,IN17,IOUT18,IN19,IOUT20,
C
C PRECEP
C
     * IN21,IN23,IN25,IN27,IN29,
     * IOUT22,IOUT24,IOUT26,IOUT28,IOUT30,
C
C MINIMIZATION
C
     * IN31,IN33,IN35,IN37,IN39,
     * IOUT32,IOUT34,IOUT36,IOUT38,IOUT39,IOUT40,
C
C NORMAL MODE
C
     * IN41,IN43,IN45,IN47,IN49,
     * IOUT42,IOUT44,IOUT46,IOUT48,IOUT50,IOUT52,IOUT54,IOUT56,IOUT58,
     * IOUT60,IOUT62,IOUT64,IOUT66,
C
C MONTE CARLO
C
     * IN71,IN73,IN75,IN77,IN79,
     * IOUT72,IOUT74,IOUT76,IOUT78
C
C PRIMARY INPUT FILE ----- IN01
C PRIMARY OUTPUT FILE ---- IOUT02
C PRIMARY LOG FILE ------- IOUT06
C
C ODD NUMBERS ARE RESERVED FOR INPUT,
C EVEN NUMBERS ARE RESERVED FOR OUTPUT
C FOR EXPLICIT DEFINITIONS OF INPUT AND OUTPUT FILES, PLEASE SEE MAIN
C
      READ(IIN) NUMVAR,(EIGVAL(I),I=1,NUMVAR),
     &          ((EIGVEC(I,J),I=1,NUMVAR),J=1,NUMVAR)
      RETURN
      END
