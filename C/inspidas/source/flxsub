      SUBROUTINE WRSTD1(FL,LEVEL,IDENT,ATTRIB,NR,RECS)
      CHARACTER IDENT*8, ATTRIB*8, RECS*(*)
      INTEGER FL, LEVEL, NR
      CHARACTER*8 LEV0,LEV1,LEV2,LEV3,LEV4,LEV5,ENDL
      DATA LEV0 /'0       '/, LEV1 /'1       '/, LEV2 /'2       '/,
     &LEV3 /'3       '/, LEV4 /'4       '/, LEV5 /'5       '/, ENDL /
     &'END     '/
      IF(.NOT.(NR .EQ. 0))GOTO  3000
      NR = 1
 3000 CONTINUE
      IF(.NOT.(LEVEL .LT. 0))GOTO  3002
      RECS = ENDL
      NR = 1
      GOTO  3003
 3002 CONTINUE
      IF(.NOT.(LEVEL .EQ. 0))GOTO  3004
      RECS(1:8) = LEV0
 3004 CONTINUE
      IF(.NOT.(LEVEL .EQ. 1))GOTO  3006
      RECS(1:8) = LEV1
 3006 CONTINUE
      IF(.NOT.(LEVEL .EQ. 2))GOTO  3008
      RECS(1:8) = LEV2
 3008 CONTINUE
      IF(.NOT.(LEVEL .EQ. 3))GOTO  3010
      RECS(1:8) = LEV3
 3010 CONTINUE
      IF(.NOT.(LEVEL .EQ. 4))GOTO  3012
      RECS(1:8) = LEV4
 3012 CONTINUE
      IF(.NOT.(LEVEL .EQ. 5))GOTO  3014
      RECS(1:8) = LEV5
 3014 CONTINUE
      RECS( 9:16) = IDENT
      RECS(17:24) = ATTRIB
      RECS(25:40) = ' '
 3003 CONTINUE
      I=1
 3016 IF(.NOT.(I.LE.NR))GOTO  3018
      WRITE (FL,'(A80)')RECS(80*(I-1)+1:80*I)
      I=I+1
      GOTO  3016
 3018 CONTINUE
      RETURN
      END
      SUBROUTINE ADINF1(KEY,NI,INFORM,NR,RECS)
      CHARACTER KEY*8, INFORM*(*), RECS*(*)
      INTEGER NI, NR
      CHARACTER*8 CONT,SPCE
      DATA CONT/'-       '/, SPCE/'        '/
      IF(.NOT.(NR .EQ. 0))GOTO  3000
      I=1
      GO TO 2
 3000 CONTINUE
      IF(.NOT.(NR .EQ. 1))GOTO  3002
      IF(.NOT.(RECS(41:48) .EQ. SPCE))GOTO  3004
      I=1
      GO TO 2
 3004 CONTINUE
      I=2
      GO TO 1
 3002 CONTINUE
      IF(.NOT.(RECS(80*(NR-1)+1:80*(NR-1)+8) .NE. SPCE))GOTO  3006
      I=NR+1
      GO TO 1
 3006 CONTINUE
      IF(.NOT.(RECS(80*(NR-1)+41:80*(NR-1)+48) .EQ. SPCE))GOTO  3008
      I=NR
      GO TO 2
 3008 CONTINUE
      I=NR+1
      GO TO 1
 1    RECS(80*(I-1)+ 1:80*(I-1)+80) = ' '
      RECS(80*(I-1)+ 9:80*(I-1)+16) = KEY
      RECS(80*(I-1)+17:80*(I-1)+40) = INFORM
      NNI = 25
      GO TO 4
 2     RECS(80*(I-1)+41:80*(I-1)+48) = KEY
      RECS(80*(I-1)+49:80*(I-1)+72) = INFORM
      NNI = 25
      I = I+1
      GO TO 3
 3010 IF(.NOT.(I.LE.LEN(RECS)/80))GOTO  3012
 3     IF(.NOT.(NNI .GT. NI))GOTO  3013
      GOTO  3012
 3013 CONTINUE
      RECS(80*(I-1)+ 1:80*(I-1)+80) = ' '
      RECS(80*(I-1)+ 9:80*(I-1)+16) = CONT
      RECS(80*(I-1)+17:80*(I-1)+40) = INFORM(NNI:NI)
      NNI = NNI + 24
 4     IF(.NOT.(NNI .GT. NI))GOTO  3015
      I=I+1
      GOTO  3012
 3015 CONTINUE
      RECS(80*(I-1)+41:80*(I-1)+48) = CONT
      RECS(80*(I-1)+49:80*(I-1)+72) = INFORM(NNI:NI)
      NNI = NNI + 24
      I=I+1
      GOTO  3010
 3012 CONTINUE
      NR = I-1
      RETURN
      END
