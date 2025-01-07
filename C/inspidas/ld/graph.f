C
C =============================================================
C       Graph making subroutine
C                                 Programmed by Koji Oda 1990.5
C =============================================================
C
      SUBROUTINE GRAPH(IOUT,TYPE,TITLE,COLMS,LINES,NDATA,
     &                 X,DX,XNAME,
     &                 Y,DY,YNAME)
C ==========
C parameters
C ==========
      implicit real*8 (A-H,O-Z)
      INTEGER IOUT
      CHARACTER TYPE*4
      CHARACTER*(*) TITLE
      INTEGER COLMS, LINES
      INTEGER NDATA
      REAL*8 X(128),Y(128)
      INTEGER DX,DY
      CHARACTER*(*) XNAME, YNAME
C
      CHARACTER SCR(64)*128
      CHARACTER XFMT*8, YFMT*8
      INTEGER XFMTL, YFMTL
      INTEGER XFLOR, YFLOR, XFRAC, YFRAC
      INTEGER LFTMAR, RGTMAR, TOPMAR, BTMMAR
      REAL*8 SCALE
      REAL*8 XMIN,XMAX, YMIN,YMAX
      REAL*8 DMY
      REAL*8 XX,YY
      INTEGER CX,CY
      REAL*8 XRANGE, YRANGE, CXR, CYR
      REAL*8 PX,PY
      GETX(PX) = NINT(LFTMAR             + CXR * (PX - XMIN) / XRANGE)
      GETY(PY) = NINT(LINES - BTMMAR + 1 - CYR * (PY - YMIN) / YRANGE)
C
C
C ====================
C initialize variables
C ====================
C
C
C initialize SCR()
C ----------------
      DO 301 I=1,LINES
        DO 302 J=1,COLMS
          SCR(I)(J:J)=' '
  302   CONTINUE
  301 CONTINUE
C
C get maximum & minimum value of each data
C ----------------------------------------
      XMIN = X(1)
      XMAX = X(1)
      YMIN = Y(1)
      YMAX = Y(1)
      DO 201 I = 2, NDATA
        IF (X(I).GT.XMAX) THEN
          XMAX = X(I)
        ELSE IF (X(I).LT.XMIN) THEN
            XMIN = X(I)
        END IF
        IF (Y(I).GT.YMAX) THEN
          YMAX = Y(I)
        ELSE IF (Y(I).LT.YMIN) THEN
            YMIN = Y(I)
        END IF
  201 CONTINUE
C
C     XRANGE = XMAX - XMIN
      YRANGE = YMAX - YMIN
C
C
C rearrange maximum & minimum
C ---------------------------
C     XMAX = XMAX + XRANGE * 0.1
C     XMIN = XMIN - XRANGE * 0.1
      YMAX = YMAX + YRANGE * 0.1
      YMIN = YMIN - YRANGE * 0.1
      IF (TYPE.EQ.'HIST') THEN
        XMIN = 0.0
        YMIN = 0.0
      END IF
      XRANGE = XMAX - XMIN
      YRANGE = YMAX - YMIN
C
C
C set format
C ----------
      IF (ABS(XMAX).GT.ABS(XMIN)) THEN
        DMY=XMAX
      ELSE
        DMY=XMIN
      END IF
C
      XFLOR =  INT(LOG10(ABS(DMY)))+1
      IF (XFRAC.LE.0) THEN
        XFRAC = 0
      END IF
C
      IF (ABS(YMAX).GT.ABS(YMIN)) THEN
        DMY=YMAX
      ELSE
        DMY=YMIN
      END IF
      YFLOR =  INT(LOG10(ABS(DMY)))+1
      YFRAC = -INT(LOG10(ABS(YRANGE))) + 1
      IF (YFRAC.LE.0) THEN
        YFRAC = 0
      END IF
      XFMTL = XFLOR + XFRAC + 2
      YFMTL = YFLOR + YFRAC + 2
C
      XFMT = '(F  .  )'
      YFMT = '(F  .  )'
      WRITE(XFMT(3:4),'(I2)') XFMTL
      WRITE(XFMT(6:7),'(I2)') XFRAC
      WRITE(YFMT(3:4),'(I2)') YFMTL
      WRITE(YFMT(6:7),'(I2)') YFRAC
C
C set margin
C ----------
C              YNAME  Y-VARIABLE'S FORMAT  Y-AXIS('I')
      LFTMAR = 2    + YFMTL              + 2
      RGTMAR = XFMTL / 2
C              TITILE
      TOPMAR = 2
C              XNAME  X-AXIS('-')  X-VARIABLE
      BTMMAR = 1    + 1          + 1
C
      CXR = DBLE(COLMS - LFTMAR - RGTMAR)
      CYR = DBLE(LINES - TOPMAR - BTMMAR)
C
C
C ============
C main routine
C ============
C
C
C   put Origin
C   ----------
      SCR(LINES - BTMMAR + 1)(LFTMAR:LFTMAR) = '+'
C
C   draw X-axis
C   -----------
      DO 101 I = LFTMAR + 1, COLMS
        SCR(LINES - BTMMAR + 1)(I:I) = '-'
  101 CONTINUE
C
C   draw Y-axis
C   -----------
      DO 102 I = TOPMAR + 1, LINES - BTMMAR
        SCR(I)(LFTMAR:LFTMAR) = '!'
  102 CONTINUE
C
C   put title
C   ---------
      I0 = INT(COLMS / 2) - INT(LEN(TITLE) / 2)
      SCR(1)(I0:I0 + LEN(TITLE) - 1) = TITLE
C
C   put X-name
C   ----------
      I0 = LFTMAR + INT((COLMS - LFTMAR) / 2)
     &     - INT(LEN(XNAME) / 2)
      SCR(LINES)(I0:I0 + LEN(XNAME) - 1) = XNAME
C
C   put Y-name
C   ----------
      I0 = TOPMAR + INT((LINES - (TOPMAR + BTMMAR)) / 2)
     &     - INT(LEN(YNAME) / 2)
      DO 103 I = 1, LEN(YNAME)
        SCR(I0 + I)(1:1) = YNAME(I:I)
  103 CONTINUE
C
C   put X-variables
C   ---------------
      SCALE = 10.0**INT(LOG10(ABS(XRANGE)/DX))
      DO 106 I = 0,DX
        XX = XMIN + I * XRANGE / DX
        XX = INT(XX / SCALE) * SCALE
        CX = GETX(XX)
        IF (CX.GE.0.AND.CX.LE.COLMS) THEN
          SCR(LINES-BTMMAR+1)(CX:CX) = '+'
          WRITE(SCR(LINES-BTMMAR+2)
     &         (CX-XFMTL/2:CX-XFMTL/2+XFMTL-1),XFMT) XX
       END IF
  106 CONTINUE
C
C   put Y-variables
C   ---------------
      SCALE = 10.0**INT(LOG10(ABS(YRANGE)/DY))
      DO 105 I = 0,DY
        YY = YMIN + I * YRANGE / DY
        YY = INT(YY / SCALE) * SCALE
        CY = GETY(YY)
        IF (CY.GE.0.AND.CY.LE.(LINES-BTMMAR)) THEN
          SCR(CY)(LFTMAR:LFTMAR) = '+'
          WRITE(SCR(CY)(3:LFTMAR - 1),YFMT) YY
        END IF
  105 CONTINUE
C
C   put data
C   --------
      DO 104 I = 1, NDATA
        CX = GETX(X(I))
        CY = GETY(Y(I))
        SCR(CY)(CX:CX) = '*'
        IF (TYPE.EQ.'HIST') THEN
          DO 107 J = CY+1, LINES - BTMMAR
            SCR(J)(CX:CX) = '*'
  107     CONTINUE
        END IF
  104 CONTINUE
C
C =========
C print out
C =========
C
      WRITE(IOUT,601)
      WRITE(IOUT,602) (SCR(I)(1:COLMS), I = 1, LINES)
  601 FORMAT(1H1)
  602 FORMAT(' ',A)
      RETURN
      END
