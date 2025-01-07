C             NICER               LEVEL=1        DATE=84.12.08
      SUBROUTINE NSHOUS(A,NMAX,NN,NNE,NNV,EPS,IORD,E,V,ILL,
     *           W1,W2,W3,W4,W5,W6,W7)
CI--------------------------------------------------------------------I
CI     SUBPROGRAM FOR STANDARD EIGEN-PROBLEM , A*V=V*E ,   BY         I
CI (HOUSEHOLDER)-(BISECTION & NO-ROOT-QR)-(INVERSE-ITERATION) METHOD  I
CI-----------------------------------------------(VERSION-2,LEVEL-1)--I
      DIMENSION A(NMAX,NN),E(NN),V(NMAX,NNV)      ,LS(20),KS(20)
      DIMENSION W1(NN),W2(NN),W3(NN),W4(NN),W5(NN),W6(NN),W7(1)
C     W7(I) IS NEVER USED IN VERSION-2.
      DATA ZERO,HALF,FIBONA,ONE /0.E 0,0.5E 0,0.6180339E 0,1.E 0/
      DATA EXPM30,EXPM20,EXPM6  /1.E-30,1.E-20,1.E-6/
      DATA EXPP6,EXPP12,EXPP18  /1.E 6,1.E12,1.E18/
C     IAP=0 FOR SCALAR MACHINE   &   IAP=1 FOR VECTOR MACHINE
      DATA IAP /0/
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.0) RETURN
      IF(NE.GT.N.OR.NV.GT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
      NM2=N-2
      NE8=NE*8
      IF(N.EQ.2)                  GO TO 130
C  TRI-DIAGONALIZATION                      ( BY A.HOUSEHOLDER )
      DO 120 K=1,NM2
      KP1=K+1
      E(K)=A(K,K)
      SUM=ZERO
      DO 10 J=KP1,N
      E(J)=A(K,J)
   10 SUM=E(J)*E(J)+SUM
      S=SIGN(SQRT(SUM),E(KP1))
      W1(K)=-S
      E(KP1)=E(KP1)+S
      A(K,KP1)=E(KP1)
      H=E(KP1)*S
      IF(H.EQ.ZERO)     GO TO 110
      SUMM=ZERO
      DO 40 I=KP1,NM1
      SUM=ZERO
      DO 20 J=KP1,I
   20 SUM=A(J,I)*E(J)+SUM
      IP1=I+1
      DO 30 J=IP1,N
   30 SUM=A(I,J)*E(J)+SUM
      W1(I)=SUM/H
      SUMM=W1(I)*E(I)+SUMM
   40 CONTINUE
      SUM=ZERO
      DO 50 J=KP1,N
   50 SUM=A(J,N)*E(J)+SUM
      W1(N)=SUM/H
      SUMM=W1(N)*E(N)+SUMM
      U=SUMM*HALF/H
      IF(IAP.EQ.1)  GO TO 70
      DO 60 J=KP1,N
      W1(J)=E(J)*U-W1(J)
      DO 60 I=KP1,J
   60 A(I,J)=W1(I)*E(J)+W1(J)*E(I)+A(I,J)
      GO TO 110
   70 CONTINUE
      DO 100 J=KP1,N
      W1(J)=E(J)*U-W1(J)
      EJ=E(J)
      DO 80 I=KP1,J
   80 A(I,J)=W1(I)*EJ+A(I,J)
      W1J=W1(J)
      DO 90 I=KP1,J
   90 A(I,J)=E(I)*W1J+A(I,J)
  100 CONTINUE
  110 A(K,K)=H
  120 CONTINUE
  130 E(NM1)=A(NM1,NM1)
      E(N)=A(N,N)
      W1(NM1)=A(NM1,N)
      W1(N)=ZERO
      GERSCH=ABS(E(1))+ABS(W1(1))
      DO 140 I=1,NM1
      SUM=ABS(E(I+1))+ABS(W1(I))+ABS(W1(I+1))
      IF(SUM.GT.GERSCH) GERSCH=SUM
      IF(NV.NE.0)  V(I,NV)=E(I)
  140 CONTINUE
      IF(NV.NE.0)  V(N,NV)=E(N)
      DEL=EPS*GERSCH
      IF(DEL.EQ.ZERO)             RETURN
      IF(NE8.LT.N)                GO TO 360
C  NO-ROOT-QR METHOD FOR EIGENVALUES (BY PAL-WALKER-KAHAN & M.SHIMASAKI)
      DD=DEL*DEL
      DO 150 I=1,NM1
  150 W3(I+1)=W1(I)*W1(I)
      K=N
  160 KM1=K-1
      IF(W3(K).LT.DD) GO TO 190
      EPE=(E(KM1)+E(K))*HALF
      EME=E(K)-EPE
      QRSHIF=EPE-SIGN(SQRT(W3(K)+EME*EME),EPE)
      CC=ONE
      SS=ZERO
      G=E(1)-QRSHIF
      PP=G*G
      DO 180 I=1,KM1
      BB=W3(I+1)
      TT=PP+BB
      W3(I)=SS*TT
      OLDCC=CC
      SS=BB/TT
      CC=PP/TT
      OLDG=G
      IF(CC.EQ.ZERO) GO TO 170
      G=(E(I+1)-QRSHIF)*CC-OLDG*SS
      E(I)=E(I+1)+OLDG-G
      PP=(G*G)/CC
      GO TO 180
  170 G=-OLDG
      E(I)=E(I+1)+OLDG+OLDG
      PP=BB*OLDCC
  180 CONTINUE
      E(K)=G+QRSHIF
      W3(K)=SS*PP
      GO TO 160
  190 K=K-1
      IF(K.GT.1) GO TO 160
C  QUICK SORT OF EIGENVALUES                ( BY C.HOARE )
      ISP=0
      L=1
      K=N
  200 IF(K-L.LT.16)          GO TO 290
      M=(K+L)/2
      MAX=K
      IF(E(M).GT.E(K))   MAX=M
      IF(E(L).GT.E(MAX)) MAX=L
      IF(MAX.EQ.K)           GO TO 210
      BK=E(MAX)
      E(MAX)=E(K)
      E(K)=BK
  210 IF(E(L).GE.E(M))       GO TO 220
      BK=E(L)
      E(L)=E(M)
      E(M)=BK
  220 BK=E(L)
      I=L
      J=K
      GO TO 250
  230 E(J)=E(I)
      E(I)=BK
  240 J=J-1
  250 IF(BK.LT.E(J))         GO TO 240
      IF(J.LE.I)             GO TO 270
      E(I)=E(J)
      E(J)=BK
  260 I=I+1
      IF(E(I).LT.BK)         GO TO 260
      IF(J.GT.I)             GO TO 230
  270 ISP=ISP+1
      IF(K-I.GE.I-L)         GO TO 280
      LS(ISP)=L
      KS(ISP)=I-1
      L=I+1
      GO TO 200
  280 LS(ISP)=I+1
      KS(ISP)=K
      K=I-1
      GO TO 200
  290 IF(K-L.LT.1)           GO TO 330
      J=K
  300 BK=E(J-1)
      I=J
  310 IF(E(I).GE.BK)         GO TO 320
      E(I-1)=E(I)
      I=I+1
      IF(I.LE.K)             GO TO 310
  320 E(I-1)=BK
      J=J-1
      IF(J.GT.L)             GO TO 300
  330 IF(ISP.EQ.0)           GO TO 340
      L=LS(ISP)
      K=KS(ISP)
      ISP=ISP-1
      GO TO 200
  340 IF(IORD.LT.0)               GO TO 500
      NP1=N+1
      NE2=N/2
      DO 350 K=1,NE2
      TEMP=E(K)
      E(K)=E(NP1-K)
  350 E(NP1-K)=TEMP
      GO TO 500
C  BISECTION METHOD FOR EIGENVALUES         ( BY J.GIVENS)
  360 CONTINUE
      DO 370 I=1,NE
      W2(I)=E(I)
      E(I)=-GERSCH
      W4(I)=GERSCH
  370 W3(I)=-W1(I)*W1(I)
      NEP1=NE+1
      DO 380 I=NEP1,N
      W2(I)=E(I)
  380 W3(I)=-W1(I)*W1(I)
      IF(IORD.GT.0)               GO TO 400
      DO 390 I=1,N
  390 W2(I)=-W2(I)
  400 CONTINUE
      DO 490 K=1,NE
  410 X=(E(K)+W4(K))*HALF
      IF(ABS(W4(K)-X).LE.DEL)     GO TO 490
      NAG=0
      I=1
  420 S=W2(I)-X
  430 IF(S.GE.ZERO) NAG=NAG+1
      IF(ABS(S).LT.EXPM30)        GO TO 440
      I=I+1
      IF(I.GT.N)                  GO TO 450
      S=W3(I-1)/S+W2(I)-X
      GO TO 430
  440 I=I+2
      IF(I.LE.N)                  GO TO 420
  450 IF(NAG.GE.K)                GO TO 470
      DO 460 J=K,NE
      IF(X.LT.W4(J)) W4(J)=X
  460 CONTINUE
      GO TO 410
  470 MG=NAG
      IF(NE.LT.MG) MG=NE
      DO 480 J=K,MG
  480 E(J)=X
      GO TO 410
  490 E(K)=X
  500 IF(NV.EQ.0)                 GO TO 810
      IF(NE8.GE.N .OR. IORD.GT.0) GO TO 520
      DO 510 I=1,N
      W1(I)=-W1(I)
  510 V(I,NV)=-V(I,NV)
C   INVERSE ITERATION FOR EIGENVECTORS      ( BY H.WIELANDT )
  520 FN=FLOAT(N)
      SN=SQRT(FN)
      SEPS=SQRT(EPS)
      EPS1=(GERSCH*EXPM6)/(FN*SEPS)
      EPS2=EXPP12*FN/SEPS
      TM6N=EXPM6*SN
      TP6N=EXPP6*SN
      TP18N=EXPP18*SN
      RN=ZERO
      RA=EPS*FIBONA
      DO 530 J=1,N
      W3(J)=ZERO
      W4(J)=W1(J)
      W5(J)=V(J,NV)-E(1)
      RN=RN+RA
      IF(RN.GE.EPS) RN=RN-EPS
  530 W6(J)=RN
      IG=1
      DO 750 I=1,NV
      IM1=I-1
      IF(I.EQ.1)                       GO TO 550
      DO 540 J=1,NM1
      W3(J)=ZERO
      W4(J)=W1(J)
      W5(J)=V(J,NV)-E(I)
  540 W6(J)=V(J+1,IM1)
      W3(N)=ZERO
      W4(N)=W1(N)
      W5(N)=V(N,NV)-E(I)
      W6(N)=V(1,IM1)
  550 CONTINUE
      DO 580 J=1,NM1
      IF(ABS(W5(J)).GE.ABS(W1(J)))     GO TO 560
      W2(J)=-W5(J)/W1(J)
      W5(J)=W1(J)
      T=W5(J+1)
      W5(J+1)=W4(J)
      W4(J)=T
      W3(J)=W4(J+1)
      IF(W3(J).EQ.ZERO) W3(J)=DEL
      W4(J+1)=ZERO
      GO TO 570
  560 IF(W5(J).EQ.ZERO) W5(J)=DEL
      W2(J)=-W1(J)/W5(J)
  570 W4(J+1)=W3(J)*W2(J)+W4(J+1)
  580 W5(J+1)=W4(J)*W2(J)+W5(J+1)
      IF(W5(N).EQ.ZERO) W5(N)=DEL
      ITELIM=2
      IF(EPS.LT.EXPM20)                          ITELIM=3
      IF(I.GT.1.AND.ABS(E(I)-E(IM1)).LT.EPS1)    ITELIM=3
      DO 640 IT=1,ITELIM
      IF(IT.EQ.1)       GO TO 600
      DO 590 J=1,NM1
      IF(W3(J).EQ.ZERO) GO TO 590
      T=W6(J)
      W6(J)=W6(J+1)
      W6(J+1)=T
  590 W6(J+1)=W6(J)*W2(J)+W6(J+1)
  600 W6(N)=W6(N)/W5(N)
      W6(NM1)=(W6(NM1)-W6(N)*W4(NM1))/W5(NM1)
      VN=ABS(W6(N))+ABS(W6(NM1))
      IF(N.EQ.2)                       GO TO 620
      DO 610 KK=2,NM1
      K=N-KK
      W6(K)=(W6(K)-W6(K+1)*W4(K)-W6(K+2)*W3(K))/W5(K)
  610 VN=ABS(W6(K))+VN
  620 IF(VN.GT.TM6N.AND.VN.LT.TP6N)    GO TO 640
      IF(IT.EQ.ITELIM)                 GO TO 650
      IF(VN.GT.EPS2)                   GO TO 650
      DUMP=FN/VN
      DO 630 J=1,N
  630 W6(J)=W6(J)*DUMP
  640 CONTINUE
C   RE-ORTHOGONALIZATION                    ( BY GRAM & SCHMIDT )
  650 CONTINUE
      DO 660 J=IG,I
      IF(ABS(E(J)-E(I)).LT.EPS1)       GO TO 670
  660 CONTINUE
      J=I
  670 IG=J
      IF(IG.EQ.I.AND.VN.LT.TP18N)      GO TO 720
      DUMP=FN/VN
      DO 680 J=1,N
  680 W6(J)=W6(J)*DUMP
      IF(IG.EQ.I)                      GO TO 720
      DO 710 K=IG,IM1
      SUM=ZERO
      DO 690 J=1,N
  690 SUM=V(J,K)*W6(J)+SUM
      S=-SUM
      DO 700 J=1,N
  700 W6(J)=V(J,K)*S+W6(J)
  710 CONTINUE
C   NORMALIZATION
  720 SUM=ZERO
      DO 730 J=1,N
  730 SUM=W6(J)*W6(J)+SUM
      SINV=ONE/SQRT(SUM)
      DO 740 J=1,N
  740 V(J,I)=W6(J)*SINV
  750 CONTINUE
C   BACK-TRANSFORMATION OF EIGEN-VECTORS
      IF(N.EQ.2)                       GO TO 810
      DO 800 J=1,NM2
      K=N-J-1
      IF(A(K,K).EQ.ZERO)               GO TO 800
      DO 760 KK=K,N
  760 W1(KK)=A(K,KK)
      KP1=K+1
      DO 790 I=1,NV
      SUM=ZERO
      DO 770 KK=KP1,N
  770 SUM=V(KK,I)*W1(KK)+SUM
      S=-SUM/W1(K)
      DO 780 KK=KP1,N
  780 V(KK,I)=W1(KK)*S+V(KK,I)
  790 CONTINUE
  800 CONTINUE
  810 ILL=0
      IF(NE8.GE.N .OR. IORD.GT.0) RETURN
      DO 820 I=1,NE
  820 E(I)=-E(I)
      RETURN
      END
      SUBROUTINE NSJENS(A,NMAX,NN,NNE,NNV,EPS,ITER,ESHIFT,
     *                  E,V,ILL,W1,W2,U)
CI--------------------------------------------------------------------I
CI       SUBPROGRAM FOR STANDARD EIGEN-PROBLEM , A*V=V*E , BY         I
CI      JENNINGS METHOD                                               I
CI-----------------------------------------------(VERSION-2,LEVEL-1)--I
      DIMENSION  A(NMAX,NN),E(NNE),V(NMAX,NNV)
      DIMENSION  W1(NN),W2(NN),U(NMAX,NNV)
      DATA ZERO,HALF,V0985,ONE /0.E 0,0.5E 0,0.985E 0,1.E 0/
      DATA TWO,TEN,HUND        /2.E 0,1.E 1,1.E 2/
      DATA INDEX/0/
      IF(INDEX.EQ.0)  WRITE(6,1000)
 1000 FORMAT(1H //5X,70HPACKAGE-NAME : NICER(NAGOYA ITERATIVE COMPUTATIO
     *N EIGENVALUE ROUTINES)     /20X,
     * 50H(VERSION-2,LEVEL-1) MODIFIED ON JAN. 1983         / 5X,
     * 90HREFERENCE    : Y.BEPPU AND I.NINOMIYA;COMPUTER PHYSICS COMMUNI
     *CATIONS, V.23,Y.1981,P.123       /)
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.1) RETURN
      IF(NE.GT.N.OR.NV.LT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
      NEP1=NE+1
      NVM1=NV-1
      IF(INDEX.NE.0)                   GO TO 10
      SEPS=SQRT(EPS*HUND)
      SEPS10=SEPS*TEN
      SSEPS=SQRT(SEPS10)
      INDEX=1
   10 SEPS2=SEPS*( ABS(A(1,1)) + ABS(A(N,N)) )
      ITELIM=ITER
      ITERE=0
      IF(ABS(ESHIFT).LT.SEPS2)         GO TO 40
      DO 20 I=1,NE
      A(I,I)=A(I,I)-ESHIFT
   20 E(I)=E(I)-ESHIFT
      DO 30 I=NEP1,N
   30 A(I,I)=A(I,I)-ESHIFT
CC   TRANSFORMATION OF LARGE SQUARE-MATRIX-A TO SMALL SQUARE-MATRIX-V
CC  BY MATRIX-MULTIPLICATION ; FIRST U=A*V , THEN V=TRANS(V)*U .
   40 CONTINUE
      DO 80 I=1,NM1
      DO 50 K=1,I
   50 W1(K)=A(K,I)
      IP1=I+1
      DO 60 K=IP1,N
   60 W1(K)=A(I,K)
      DO 80 J=1,NV
      SUM=ZERO
      DO 70 K=1,N
   70 SUM=V(K,J)*W1(K)+SUM
      U(I,J)=SUM
   80 CONTINUE
      DO 100 J=1,NV
      SUM=ZERO
      DO 90 K=1,N
   90 SUM=A(K,N)*V(K,J)+SUM
      U(N,J)=SUM
  100 CONTINUE
      DO 130 J=1,NV
      DO 110 K=1,N
  110 W1(K)=V(K,J)
      DO 130 I=J,NV
      SUM=ZERO
      DO 120 K=1,N
  120 SUM=U(K,I)*W1(K)+SUM
      V(I,J)=SUM
  130 CONTINUE
CC  SOLVING EIGEN-PROBLEM IN SUB-SPACE
      ERROR=ABS(E(1)-V(1,1))
      E(1)=V(1,1)
      IF(NV.EQ.1)                      GO TO 630
      DO 140 K=2,NV
      ERROR=ABS(E(K)-V(K,K))+ERROR
  140 E(K)=V(K,K)
      E1ENE=ABS(E(1)) + ABS(E(NE))
      ERROR=ERROR/( E1ENE*FLOAT(NV) )
C
      IF(ERROR.GT.SSEPS)               GO TO 260
C          BY APPROXIMATE FORMULA
      DEGENE=SEPS*E1ENE
      K=1
  150 CONTINUE
      ND=1
  160 KND=K+ND
      IF(KND.GT.NV)                    GO TO 170
      IF(ABS(E(KND)-E(K)).GT.DEGENE)   GO TO 170
      ND=ND+1
      GO TO 160
  170 IF(ND.EQ.1)                      GO TO 190
      DD=FLOAT(ND)*SEPS10
      ND2=ND/2
      LI=K
      LF=K+ND2-1
      LM=ND-1
      DW=FLOAT(LM)
      DO 180 L=LI,LF
      DDDW=DD*DW
      V(L,L)=   DDDW*V(L,L)
      LL=LI+LM
      V(LL,LL)=-DDDW*V(LL,LL)
      DW=DW-TWO
      LM=LM-1
  180 CONTINUE
  190 K=K+ND
      IF(K.LT.NV)                      GO TO 150
C
      DO 200 I=1,NVM1
      IP1=I+1
      DO 200 J=IP1,NV
      V(I,J)=V(J,I) / ( V(J,J)-V(I,I) )
  200 V(J,I)=-V(I,J)
C
      SUM=ZERO
      DO 210 J=2,NV
  210 SUM=V(J,1)*V(J,1)+SUM
      IF(SUM.GT.ONE)    GO TO 230
      S1=SQRT(ONE-SUM)
      SUM=ZERO
      DO 220 J=1,NVM1
  220 SUM=V(J,NV)*V(J,NV)+SUM
      IF(SUM.GT.ONE)    GO TO 230
      SNV=SQRT(ONE-SUM)
      VDIAGO=HALF*(S1+SNV)
      ILLIND=0
      GO TO 240
  230 VDIAGO=V0985
      ILLIND=200
  240 CONTINUE
      DO 250 K=1,NV
  250 V(K,K)=VDIAGO
      GO TO 510
C          BY HOUSEHOLDER-QR FORMULA
  260 NVM2=NV-2
      DO 270 J=2,NV
      J1=J-1
      DO 270 I=1,J1
  270 V(I,J)=V(J,I)
      IF(NV.LE.2) GO TO 360
      DO 350 K=1,NVM2
      K1=K+1
      E(K)=V(K,K)
      SUM=ZERO
      DO 280 J=K1,NV
      E(J)=V(K,J)
  280 SUM=E(J)*E(J)+SUM
      SUM=SIGN(SQRT(SUM),E(K1))
      W1(K)=-SUM
      E(K1)=E(K1)+SUM
      V(K,K1)=E(K1)
      C=E(K1)*SUM
      IF(C.LE.ZERO) GO TO 340
      SUMM=ZERO
      DO 320 I=K1,NV
      SUM=ZERO
      DO 290 J=K1,I
  290 SUM=V(J,I)*E(J)+SUM
      IF(I.GE.NV) GO TO 310
      IP1=I+1
      DO 300 J=IP1,NV
  300 SUM=V(I,J)*E(J)+SUM
  310 W1(I)=SUM/C
  320 SUMM=E(I)*W1(I)+SUMM
      W=SUMM*HALF/C
      DO 330 J=K1,NV
      W1(J)=E(J)*W-W1(J)
      DO 330 I=K1,J
  330 V(I,J)=E(J)*W1(I)+E(I)*W1(J)+V(I,J)
  340 V(K,K)=C
  350 CONTINUE
  360 E(NV)=V(NV,NV)
      V(NV,NV)=ONE
      IF(NV.EQ.1) GO TO 510
      E(NVM1)=V(NVM1,NVM1)
      W1(NVM1)=V(NVM1,NV)
      W1(NV)=ZERO
      V(NVM1,NVM1)=ONE
      V(NVM1,NV)=ZERO
      V(NV,NVM1)=ZERO
      IF(NV.EQ.2) GO TO 420
      DO 410 L=1,NVM2
      K=NVM1-L
      K1=K+1
      C=-V(K,K)
      V(K,K)=ONE
      IF(C.GE.ZERO) GO TO 390
      DO 380 J=K1,NV
      SUM=ZERO
      DO 370 I=K1,NV
  370 SUM=V(I,J)*V(K,I)+SUM
      S=SUM/C
      DO 380 I=K1,NV
  380 V(I,J)=S*V(K,I)+V(I,J)
  390 CONTINUE
      DO 400 I=K1,NV
      V(K,I)=ZERO
  400 V(I,K)=ZERO
  410 CONTINUE
  420 D=ABS(E(1))
      DO 430 J=2,NV
      IF(( ABS(W1(J-1))+ABS(E(J)) ).GT.D) D=ABS(W1(J-1))+ABS(E(J))
  430 CONTINUE
      IF(D.EQ.ZERO) GO TO 510
      D=D*EPS
      K=NV
  440 L=K
  450 IF(ABS(W1(L-1)).LT.D) GO TO 460
      L=L-1
      IF(L.GT.1) GO TO 450
  460 IF(L.EQ.K) GO TO 500
      WWW=(E(K-1)+E(K))*HALF
      R=E(K)-WWW
      Z=WWW-SIGN(SQRT(W1(K-1)*W1(K-1)+R*R),WWW)
      EE=E(L)-Z
      E(L)=EE
      FF=W1(L)
      R=SQRT(EE*EE+FF*FF)
      J=L
      GO TO 480
  470 R=SQRT(E(J)*E(J)+W1(J)*W1(J))
      W1(J-1)=S*R
      EE=C*E(J)
      FF=C*W1(J)
  480 C=E(J)/R
      S=W1(J)/R
      WWW=E(J+1)-Z
      E(J)=(FF*C+WWW*S)*S+EE+Z
      E(J+1)=C*WWW-S*FF
      DO 490 I=1,NV
      R=V(I,J+1)
      V(I,J+1)=R*C-V(I,J)*S
  490 V(I,J)=V(I,J)*C+R*S
      J=J+1
      IF(J.LT.K) GO TO 470
      W1(K-1)=E(K)*S
      E(K)=E(K)*C+Z
      GO TO 440
  500 K=K-1
      IF(K.GT.1) GO TO 440
C        STRAIGHT INSERTION SORT OF EIGENVALUES
  510 J=NV
  520 CONTINUE
      DO 530 L=1,NV
  530 W1(L)=V(L,J-1)
      EK=E(J-1)
      DO 550 I=J,NV
      IF(ABS(E(I)).LE.ABS(EK)) GO TO 570
      DO 540 L=1,NV
  540 V(L,I-1)=V(L,I)
  550 E(I-1)=E(I)
      DO 560 L=1,NV
  560 V(L,NV)=W1(L)
      E(NV)=EK
      GO TO 590
  570 CONTINUE
      DO 580 L=1,NV
  580 V(L,I-1)=W1(L)
      E(I-1)=EK
  590 J=J-1
      IF(J.GT.1) GO TO 520
C
CC  NORM-ORTHOGONALIZATION OF EIGEN-VECTORS
      DO 620 I=1,N
      DO 600 K=1,NV
  600 W1(K)=U(I,K)
      DO 620 J=1,NV
      SUM=ZERO
      DO 610 K=1,NV
  610 SUM=V(K,J)*W1(K)+SUM
      U(I,J)=SUM
  620 CONTINUE
C
  630 CONTINUE
      DO 660 J=1,NV
      DO 640 K=1,N
  640 W1(K)=U(K,J)
      DO 660 I=1,J
      SUM=ZERO
      DO 650 K=1,N
  650 SUM=U(K,I)*W1(K)+SUM
      V(I,J)=SUM
  660 CONTINUE
C
      IF(V(1,1).LT.ZERO)      ILLIND=200
      IF(V(NV,NV).LT.ZERO)    ILLIND=200
      IF(V(1,1).LT.V(NV,NV))  ILLIND=200
      W2(1)=ONE/SQRT(V(1,1))
      IF(NV.EQ.1)             GO TO 710
      V(1,2)=V(1,2)*W2(1)
      W2(2)=ONE/SQRT(  V(2,2)-V(1,2)**2  )
      IF(NV.EQ.2)             GO TO 710
      DO 700 J=3,NV
      V(1,J)=V(1,J)*W2(1)
      JM1=J-1
      DO 680 I=2,JM1
      SUM=ZERO
      IM1=I-1
      DO 670 L=1,IM1
  670 SUM=V(L,I)*V(L,J)+SUM
      V(I,J)=(  V(I,J)-SUM  )*W2(I)
  680 CONTINUE
      SUM=ZERO
      DO 690 L=1,JM1
  690 SUM=V(L,J)*V(L,J)+SUM
      W2(J)=ONE/SQRT(  V(J,J)-SUM  )
  700 CONTINUE
C
  710 CONTINUE
      C11=W2(1)
      DO 720 J=1,N
  720 V(J,1)=U(J,1)*C11
      IF(NV.EQ.1)  GO TO 760
      DO 750 I=2,NV
      CII=W2(I)
      IM1=I-1
      DO 730 K=1,IM1
  730 W1(K)=V(K,I)
      DO 750 J=1,N
      SUM=ZERO
      DO 740 K=1,IM1
  740 SUM=V(J,K)*W1(K)+SUM
      V(J,I)=(U(J,I)-SUM)*CII
  750 CONTINUE
  760 ITERE=ITERE+1
      IF(ITERE.LT.ITELIM .AND. ERROR.GT.SEPS10)  GO TO 40
      IF(ITERE.GE.ITELIM .AND. ERROR.GT.SEPS10)  GO TO 770
C
      ILL=0
      IF(ILLIND.EQ.200)                          GO TO 770
      GO TO 780
  770 ILL=200
  780 CONTINUE
      ITER=ITERE
      IF(ABS(ESHIFT).LT.SEPS2) RETURN
      DO 790 I=1,NE
      A(I,I)=A(I,I)+ESHIFT
  790 E(I)=E(I)+ESHIFT
      DO 800 I=NEP1,N
  800 A(I,I)=A(I,I)+ESHIFT
      RETURN
      END
      SUBROUTINE NGHOUS(AB,NNMAX,NN,NNE,NNV,EPS,IORD,ICHO,BD,E,V,ILL,
     *                  W1,W2,W3,W4,W5,W6,W7)
CI--------------------------------------------------------------------I
CI       SUBPROGRAM FOR GENERALIZED EIGEN-PROBLEM , A*V=B*V*E , BY    I
CI  ( SIMULTANEOUSLY TRIANGULAR DECOMPOSITION )-                      I
CI  (HOUSEHOLDER)-(BISECTION & NO-ROOT-QR)-(INVERSE-ITERATION) METHOD I
CI-----------------------------------------------(VERSION-2,LEVEL-1)--I
      DIMENSION AB(NNMAX,NN),BD(NN),V(NNMAX,NNV),E(NNE)
      DIMENSION W1(NN),W2(NN),W3(NN),W4(NN),W5(NN),W6(NN),W7(1)
C     W7(I) IS NEVER USED IN VERSION-2.
      DATA ZERO,ONE /0.E 0,1.E 0/
      NMAX=NNMAX
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.0) RETURN
      IF(NE.GT.N.OR.NV.GT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
      DO 10 I=1,NM1
      IP1=I+1
      DO 10 J=IP1,N
      T=AB(J,I)
      AB(J,I)=AB(I,J)
   10 AB(I,J)=T
CC   CHOLESKY TRANSFORMATION OF A BY DECOMPOSED B
      IF(ICHO.EQ.1)          GO TO 60
      ILL=100
      IF(BD(1).LT.EPS)       RETURN
      BD(1)=ONE/SQRT(BD(1))
      AB(1,2)=AB(1,2)*BD(1)
      PPP=BD(2)-AB(1,2)*AB(1,2)
      IF(PPP.LT.EPS)         RETURN
      BD(2)=ONE/SQRT(PPP)
      IF(N.EQ.2)             GO TO 60
      DO 50 J=3,N
      AB(1,J)=AB(1,J)*BD(1)
      JM1=J-1
      DO 30 I=2,JM1
      SUM=ZERO
      IM1=I-1
      DO 20 L=1,IM1
   20 SUM=AB(L,I)*AB(L,J)+SUM
      AB(I,J)=(  AB(I,J)-SUM  )*BD(I)
   30 CONTINUE
      SUM=ZERO
      DO 40 L=1,JM1
   40 SUM=AB(L,J)*AB(L,J)+SUM
      PPP=BD(J)-SUM
      IF(PPP.LT.EPS)         RETURN
      BD(J)=ONE/SQRT(PPP)
   50 CONTINUE
C
   60 CONTINUE
      DO 70 J=1,N
   70 AB(J,1)=AB(J,1)*BD(1)
      DO 90 I=2,N
      BDI=BD(I)
      IM1=I-1
      DO 90 J=I,N
      SUM=ZERO
      DO 80 K=1,IM1
   80 SUM=AB(K,I)*AB(J,K)+SUM
      AB(J,I)=(AB(J,I)-SUM)*BDI
   90 CONTINUE
C
      AB(1,1)=AB(1,1)*BD(1)
      DO 110 I=2,N
      SUM=ZERO
      IM1=I-1
      DO 100 K=1,IM1
  100 SUM=AB(K,1)*AB(K,I)+SUM
      AB(I,1)=(AB(I,1)-SUM)*BD(I)
  110 CONTINUE
C
      DO 150 J=2,NM1
      JM1=J-1
      JP1=J+1
      SUM=ZERO
      DO 120 K=1,JM1
  120 SUM=AB(K,J)*AB(J,K)+SUM
      AB(J,J)=(AB(J,J)-SUM)*BD(J)
      DO 150 I=JP1,N
      SUM=ZERO
      DO 130 K=1,JM1
  130 SUM=AB(K,I)*AB(J,K)+SUM
      IM1=I-1
      DO 140 K=J,IM1
  140 SUM=AB(K,J)*AB(K,I)+SUM
      AB(I,J)=(AB(I,J)-SUM)*BD(I)
  150 CONTINUE
      SUM=ZERO
      DO 160 K=1,NM1
  160 SUM=AB(K,N)*AB(N,K)+SUM
      AB(N,N)=(AB(N,N)-SUM)*BD(N)
C#  SOLVING STANDARD EIGEN-PROBLEM FOR TRANSFORMED A
      DO 170 I=1,NM1
      IP1=I+1
      DO 170 J=IP1,N
      T=AB(J,I)
      AB(J,I)=AB(I,J)
  170 AB(I,J)=T
      CALL NSHOUS(AB,NMAX,N,NE,NV,EPS,IORD,E,V,ILL,W1,W2,W3,W4,W5,W6,W7)
      IF(NV.EQ.0.OR.ILL.NE.0) RETURN
C#  CONVERSION OF EIGEN-VECTORS( STANDARD FORM -> GENERALIZED FORM )
      DO 190 I=1,NV
      V(N,I)=V(N,I)*BD(N)
      DO 190 JJ=1,NM1
      J=N-JJ
      SUM=ZERO
      JP1=J+1
      DO 180 K=JP1,N
  180 SUM=AB(K,J)*V(K,I)+SUM
      V(J,I)=(V(J,I)-SUM)*BD(J)
  190 CONTINUE
      RETURN
      END
      SUBROUTINE NGJENS(AB,NNMAX,NN,NNE,NNV,EPS,BD,IUV,ITER,ESHIFT,
     *                  E,V,U,ILL,W1,W2)
CI--------------------------------------------------------------------I
CI       SUBPROGRAM FOR GENERALIZED EIGEN-PROBLEM , A*V=B*V*E , BY    I
CI      (SIMULTANEOUSLY TRIANGULAR DECOMPOSITION)-(JENNINGS) METHOD   I
CI-----------------------------------------------(VERSION-2,LEVEL-1)--I
      DIMENSION AB(NNMAX,NN),BD(NN),E(NNE),V(NNMAX,NNV),U(NNMAX,NNV)
      DIMENSION W1(NN),W2(NN)
      DATA ZERO,ONE /0.E 0,1.E 0/
      NMAX=NNMAX
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.1) RETURN
      IF(NE.GT.N.OR.NV.LT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
C   CHOLESKY TRASFORMATION OF A BY DECOMPOSED B
      DO 10 J=1,N
   10 AB(1,J)=AB(1,J)*BD(1)
      DO 30 I=2,N
      BDI=BD(I)
      IM1=I-1
      DO 30 J=I,N
      SUM=ZERO
      DO 20 K=1,IM1
   20 SUM=AB(I,K)*AB(K,J)+SUM
      AB(I,J)=(AB(I,J)-SUM)*BDI
   30 CONTINUE
C
      AB(1,1)=AB(1,1)*BD(1)
      DO 50 I=2,N
      SUM=ZERO
      IM1=I-1
      DO 40 K=1,IM1
   40 SUM=AB(1,K)*AB(I,K)+SUM
      AB(1,I)=(AB(1,I)-SUM)*BD(I)
   50 CONTINUE
C
      DO 90 J=2,NM1
      JM1=J-1
      JP1=J+1
      SUM=ZERO
      DO 60 K=1,JM1
   60 SUM=AB(K,J)*AB(J,K)+SUM
      AB(J,J)=(AB(J,J)-SUM)*BD(J)
      DO 90 I=JP1,N
      SUM=ZERO
      DO 70 K=1,JM1
   70 SUM=AB(K,J)*AB(I,K)+SUM
      IM1=I-1
      DO 80 K=J,IM1
   80 SUM=AB(J,K)*AB(I,K)+SUM
      AB(J,I)=(AB(J,I)-SUM)*BD(I)
   90 CONTINUE
      SUM=ZERO
      DO 100 K=1,NM1
  100 SUM=AB(K,N)*AB(N,K)+SUM
      AB(N,N)=(AB(N,N)-SUM)*BD(N)
C   TRANSFORMATION OF INITIAL EIGEN-VECTORS
      IF(IUV.EQ.1)           GO TO 150
      DO 130 I=1,NM1
      W2(I)=ONE/BD(I)
      IP1=I+1
      DO 110 K=IP1,N
  110 W2(K)=AB(K,I)
      DO 130 J=1,NV
      SUM=ZERO
      DO 120 K=I,N
  120 SUM=V(K,J)*W2(K)+SUM
      U(I,J)=SUM
  130 CONTINUE
      BNN=ONE/BD(N)
      DO 140 J=1,NV
  140 U(N,J)=V(N,J)*BNN
C#  SOLVING STANDARD EIGEN-PROBLEM FOR TRANSFORMED A
  150 CALL NSJENS(AB,NMAX,N,NE,NV,EPS,ITER,ESHIFT,E,U,ILLE,W1,W2,V)
C#  CONVERSION OF EIGEN-VECTORS( STANDARD FORM -> GENERALIZED FORM )
      DO 170 I=1,NV
      V(N,I)=U(N,I)*BD(N)
      DO 170 JJ=1,NM1
      J=N-JJ
      SUM=ZERO
      JP1=J+1
      DO 160 K=JP1,N
  160 SUM=AB(K,J)*V(K,I)+SUM
      V(J,I)=(U(J,I)-SUM)*BD(J)
  170 CONTINUE
      ILL=ILLE
      RETURN
      END
      SUBROUTINE NSHOUD(A,NMAX,NN,NNE,NNV,EPS,IORD,E,V,ILL,
     *           W1,W2,W3,W4,W5,W6,W7)
CI--------------------------------------------------------------------I
CI     SUBPROGRAM FOR STANDARD EIGEN-PROBLEM , A*V=V*E ,   BY         I
CI (HOUSEHOLDER)-(BISECTION & NO-ROOT-QR)-(INVERSE-ITERATION) METHOD  I
CI-----------------------------------------------(VERSION-2,LEVEL-1)--I
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NMAX,NN),E(NN),V(NMAX,NNV)      ,LS(20),KS(20)
      DIMENSION W1(NN),W2(NN),W3(NN),W4(NN),W5(NN),W6(NN),W7(1)
C     W7(I) IS NEVER USED IN VERSION-2.
      DATA ZERO,HALF,FIBONA,ONE /0.D 0,0.5D 0,0.6180339D 0,1.D 0/
      DATA EXPM30,EXPM20,EXPM6  /1.D-30,1.D-20,1.D-6/
      DATA EXPP6,EXPP12,EXPP18  /1.D 6,1.D12,1.D18/
C     IAP=0 FOR SCALAR MACHINE   &   IAP=1 FOR VECTOR MACHINE
      DATA IAP /0/
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.0) RETURN
      IF(NE.GT.N.OR.NV.GT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
      NM2=N-2
      NE8=NE*8
      IF(N.EQ.2)                  GO TO 130
C  TRI-DIAGONALIZATION                      ( BY A.HOUSEHOLDER )
      DO 120 K=1,NM2
      KP1=K+1
      E(K)=A(K,K)
      SUM=ZERO
      DO 10 J=KP1,N
      E(J)=A(K,J)
   10 SUM=E(J)*E(J)+SUM
      S=DSIGN(DSQRT(SUM),E(KP1))
      W1(K)=-S
      E(KP1)=E(KP1)+S
      A(K,KP1)=E(KP1)
      H=E(KP1)*S
      IF(H.EQ.ZERO)     GO TO 110
      SUMM=ZERO
      DO 40 I=KP1,NM1
      SUM=ZERO
      DO 20 J=KP1,I
   20 SUM=A(J,I)*E(J)+SUM
      IP1=I+1
      DO 30 J=IP1,N
   30 SUM=A(I,J)*E(J)+SUM
      W1(I)=SUM/H
      SUMM=W1(I)*E(I)+SUMM
   40 CONTINUE
      SUM=ZERO
      DO 50 J=KP1,N
   50 SUM=A(J,N)*E(J)+SUM
      W1(N)=SUM/H
      SUMM=W1(N)*E(N)+SUMM
      U=SUMM*HALF/H
      IF(IAP.EQ.1)  GO TO 70
      DO 60 J=KP1,N
      W1(J)=E(J)*U-W1(J)
      DO 60 I=KP1,J
   60 A(I,J)=W1(I)*E(J)+W1(J)*E(I)+A(I,J)
      GO TO 110
   70 CONTINUE
      DO 100 J=KP1,N
      W1(J)=E(J)*U-W1(J)
      EJ=E(J)
      DO 80 I=KP1,J
   80 A(I,J)=W1(I)*EJ+A(I,J)
      W1J=W1(J)
      DO 90 I=KP1,J
   90 A(I,J)=E(I)*W1J+A(I,J)
  100 CONTINUE
  110 A(K,K)=H
  120 CONTINUE
  130 E(NM1)=A(NM1,NM1)
      E(N)=A(N,N)
      W1(NM1)=A(NM1,N)
      W1(N)=ZERO
      GERSCH=DABS(E(1))+DABS(W1(1))
      DO 140 I=1,NM1
      SUM=DABS(E(I+1))+DABS(W1(I))+DABS(W1(I+1))
      IF(SUM.GT.GERSCH) GERSCH=SUM
      IF(NV.NE.0)  V(I,NV)=E(I)
  140 CONTINUE
      IF(NV.NE.0)  V(N,NV)=E(N)
      DEL=EPS*GERSCH
      IF(DEL.EQ.ZERO)             RETURN
      IF(NE8.LT.N)                GO TO 360
C  NO-ROOT-QR METHOD FOR EIGENVALUES (BY PAL-WALKER-KAHAN & M.SHIMASAKI)
      DD=DEL*DEL
      DO 150 I=1,NM1
  150 W3(I+1)=W1(I)*W1(I)
      K=N
  160 KM1=K-1
      IF(W3(K).LT.DD) GO TO 190
      EPE=(E(KM1)+E(K))*HALF
      EME=E(K)-EPE
      QRSHIF=EPE-DSIGN(DSQRT(W3(K)+EME*EME),EPE)
      CC=ONE
      SS=ZERO
      G=E(1)-QRSHIF
      PP=G*G
      DO 180 I=1,KM1
      BB=W3(I+1)
      TT=PP+BB
      W3(I)=SS*TT
      OLDCC=CC
      SS=BB/TT
      CC=PP/TT
      OLDG=G
      IF(CC.EQ.ZERO) GO TO 170
      G=(E(I+1)-QRSHIF)*CC-OLDG*SS
      E(I)=E(I+1)+OLDG-G
      PP=(G*G)/CC
      GO TO 180
  170 G=-OLDG
      E(I)=E(I+1)+OLDG+OLDG
      PP=BB*OLDCC
  180 CONTINUE
      E(K)=G+QRSHIF
      W3(K)=SS*PP
      GO TO 160
  190 K=K-1
      IF(K.GT.1) GO TO 160
C  QUICK SORT OF EIGENVALUES                ( BY C.HOARE )
      ISP=0
      L=1
      K=N
  200 IF(K-L.LT.16)          GO TO 290
      M=(K+L)/2
      MAX=K
      IF(E(M).GT.E(K))   MAX=M
      IF(E(L).GT.E(MAX)) MAX=L
      IF(MAX.EQ.K)           GO TO 210
      BK=E(MAX)
      E(MAX)=E(K)
      E(K)=BK
  210 IF(E(L).GE.E(M))       GO TO 220
      BK=E(L)
      E(L)=E(M)
      E(M)=BK
  220 BK=E(L)
      I=L
      J=K
      GO TO 250
  230 E(J)=E(I)
      E(I)=BK
  240 J=J-1
  250 IF(BK.LT.E(J))         GO TO 240
      IF(J.LE.I)             GO TO 270
      E(I)=E(J)
      E(J)=BK
  260 I=I+1
      IF(E(I).LT.BK)         GO TO 260
      IF(J.GT.I)             GO TO 230
  270 ISP=ISP+1
      IF(K-I.GE.I-L)         GO TO 280
      LS(ISP)=L
      KS(ISP)=I-1
      L=I+1
      GO TO 200
  280 LS(ISP)=I+1
      KS(ISP)=K
      K=I-1
      GO TO 200
  290 IF(K-L.LT.1)           GO TO 330
      J=K
  300 BK=E(J-1)
      I=J
  310 IF(E(I).GE.BK)         GO TO 320
      E(I-1)=E(I)
      I=I+1
      IF(I.LE.K)             GO TO 310
  320 E(I-1)=BK
      J=J-1
      IF(J.GT.L)             GO TO 300
  330 IF(ISP.EQ.0)           GO TO 340
      L=LS(ISP)
      K=KS(ISP)
      ISP=ISP-1
      GO TO 200
  340 IF(IORD.LT.0)               GO TO 500
      NP1=N+1
      NE2=N/2
      DO 350 K=1,NE2
      TEMP=E(K)
      E(K)=E(NP1-K)
  350 E(NP1-K)=TEMP
      GO TO 500
C  BISECTION METHOD FOR EIGENVALUES         ( BY J.GIVENS)
  360 CONTINUE
      DO 370 I=1,NE
      W2(I)=E(I)
      E(I)=-GERSCH
      W4(I)=GERSCH
  370 W3(I)=-W1(I)*W1(I)
      NEP1=NE+1
      DO 380 I=NEP1,N
      W2(I)=E(I)
  380 W3(I)=-W1(I)*W1(I)
      IF(IORD.GT.0)               GO TO 400
      DO 390 I=1,N
  390 W2(I)=-W2(I)
  400 CONTINUE
      DO 490 K=1,NE
  410 X=(E(K)+W4(K))*HALF
      IF(DABS(W4(K)-X).LE.DEL)    GO TO 490
      NAG=0
      I=1
  420 S=W2(I)-X
  430 IF(S.GE.ZERO) NAG=NAG+1
      IF(DABS(S).LT.EXPM30)       GO TO 440
      I=I+1
      IF(I.GT.N)                  GO TO 450
      S=W3(I-1)/S+W2(I)-X
      GO TO 430
  440 I=I+2
      IF(I.LE.N)                  GO TO 420
  450 IF(NAG.GE.K)                GO TO 470
      DO 460 J=K,NE
      IF(X.LT.W4(J)) W4(J)=X
  460 CONTINUE
      GO TO 410
  470 MG=NAG
      IF(NE.LT.MG) MG=NE
      DO 480 J=K,MG
  480 E(J)=X
      GO TO 410
  490 E(K)=X
  500 IF(NV.EQ.0)                 GO TO 810
      IF(NE8.GE.N .OR. IORD.GT.0) GO TO 520
      DO 510 I=1,N
      W1(I)=-W1(I)
  510 V(I,NV)=-V(I,NV)
C   INVERSE ITERATION FOR EIGENVECTORS      ( BY H.WIELANDT )
  520 FN=DFLOAT(N)
      SN=DSQRT(FN)
      SEPS=DSQRT(EPS)
      EPS1=(GERSCH*EXPM6)/(FN*SEPS)
      EPS2=EXPP12*FN/SEPS
      TM6N=EXPM6*SN
      TP6N=EXPP6*SN
      TP18N=EXPP18*SN
      RN=ZERO
      RA=EPS*FIBONA
      DO 530 J=1,N
      W3(J)=ZERO
      W4(J)=W1(J)
      W5(J)=V(J,NV)-E(1)
      RN=RN+RA
      IF(RN.GE.EPS) RN=RN-EPS
  530 W6(J)=RN
      IG=1
      DO 750 I=1,NV
      IM1=I-1
      IF(I.EQ.1)                       GO TO 550
      DO 540 J=1,NM1
      W3(J)=ZERO
      W4(J)=W1(J)
      W5(J)=V(J,NV)-E(I)
  540 W6(J)=V(J+1,IM1)
      W3(N)=ZERO
      W4(N)=W1(N)
      W5(N)=V(N,NV)-E(I)
      W6(N)=V(1,IM1)
  550 CONTINUE
      DO 580 J=1,NM1
      IF(DABS(W5(J)).GE.DABS(W1(J)))   GO TO 560
      W2(J)=-W5(J)/W1(J)
      W5(J)=W1(J)
      T=W5(J+1)
      W5(J+1)=W4(J)
      W4(J)=T
      W3(J)=W4(J+1)
      IF(W3(J).EQ.ZERO) W3(J)=DEL
      W4(J+1)=ZERO
      GO TO 570
  560 IF(W5(J).EQ.ZERO) W5(J)=DEL
      W2(J)=-W1(J)/W5(J)
  570 W4(J+1)=W3(J)*W2(J)+W4(J+1)
  580 W5(J+1)=W4(J)*W2(J)+W5(J+1)
      IF(W5(N).EQ.ZERO) W5(N)=DEL
      ITELIM=2
      IF(EPS.LT.EXPM20)                          ITELIM=3
      IF(I.GT.1.AND.DABS(E(I)-E(IM1)).LT.EPS1)   ITELIM=3
      DO 640 IT=1,ITELIM
      IF(IT.EQ.1)       GO TO 600
      DO 590 J=1,NM1
      IF(W3(J).EQ.ZERO) GO TO 590
      T=W6(J)
      W6(J)=W6(J+1)
      W6(J+1)=T
  590 W6(J+1)=W6(J)*W2(J)+W6(J+1)
  600 W6(N)=W6(N)/W5(N)
      W6(NM1)=(W6(NM1)-W6(N)*W4(NM1))/W5(NM1)
      VN=DABS(W6(N))+DABS(W6(NM1))
      IF(N.EQ.2)                       GO TO 620
      DO 610 KK=2,NM1
      K=N-KK
      W6(K)=(W6(K)-W6(K+1)*W4(K)-W6(K+2)*W3(K))/W5(K)
  610 VN=DABS(W6(K))+VN
  620 IF(VN.GT.TM6N.AND.VN.LT.TP6N)    GO TO 640
      IF(IT.EQ.ITELIM)                 GO TO 650
      IF(VN.GT.EPS2)                   GO TO 650
      DUMP=FN/VN
      DO 630 J=1,N
  630 W6(J)=W6(J)*DUMP
  640 CONTINUE
C   RE-ORTHOGONALIZATION                    ( BY GRAM & SCHMIDT )
  650 CONTINUE
      DO 660 J=IG,I
      IF(DABS(E(J)-E(I)).LT.EPS1)      GO TO 670
  660 CONTINUE
      J=I
  670 IG=J
      IF(IG.EQ.I.AND.VN.LT.TP18N)      GO TO 720
      DUMP=FN/VN
      DO 680 J=1,N
  680 W6(J)=W6(J)*DUMP
      IF(IG.EQ.I)                      GO TO 720
      DO 710 K=IG,IM1
      SUM=ZERO
      DO 690 J=1,N
  690 SUM=V(J,K)*W6(J)+SUM
      S=-SUM
      DO 700 J=1,N
  700 W6(J)=V(J,K)*S+W6(J)
  710 CONTINUE
C   NORMALIZATION
  720 SUM=ZERO
      DO 730 J=1,N
  730 SUM=W6(J)*W6(J)+SUM
      SINV=ONE/DSQRT(SUM)
      DO 740 J=1,N
  740 V(J,I)=W6(J)*SINV
  750 CONTINUE
C   BACK-TRANSFORMATION OF EIGEN-VECTORS
      IF(N.EQ.2)                       GO TO 810
      DO 800 J=1,NM2
      K=N-J-1
      IF(A(K,K).EQ.ZERO)               GO TO 800
      DO 760 KK=K,N
  760 W1(KK)=A(K,KK)
      KP1=K+1
      DO 790 I=1,NV
      SUM=ZERO
      DO 770 KK=KP1,N
  770 SUM=V(KK,I)*W1(KK)+SUM
      S=-SUM/W1(K)
      DO 780 KK=KP1,N
  780 V(KK,I)=W1(KK)*S+V(KK,I)
  790 CONTINUE
  800 CONTINUE
  810 ILL=0
      IF(NE8.GE.N .OR. IORD.GT.0) RETURN
      DO 820 I=1,NE
  820 E(I)=-E(I)
      RETURN
      END
      SUBROUTINE NSJEND(A,NMAX,NN,NNE,NNV,EPS,ITER,ESHIFT,
     *                  E,V,ILL,W1,W2,U)
CI--------------------------------------------------------------------I
CI       SUBPROGRAM FOR STANDARD EIGEN-PROBLEM , A*V=V*E , BY         I
CI      JENNINGS METHOD                                               I
CI-----------------------------------------------(VERSION-2,LEVEL-1)--I
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  A(NMAX,NN),E(NNE),V(NMAX,NNV)
      DIMENSION  W1(NN),W2(NN),U(NMAX,NNV)
      DATA ZERO,HALF,V0985,ONE /0.D 0,0.5D 0,0.985D 0,1.D 0/
      DATA TWO,TEN,HUND        /2.D 0,1.D 1,1.D 2/
      DATA INDEX/0/
      IF(INDEX.EQ.0)  WRITE(6,1000)
 1000 FORMAT(1H //5X,70HPACKAGE-NAME : NICER(NAGOYA ITERATIVE COMPUTATIO
     *N EIGENVALUE ROUTINES)     /20X,
     * 50H(VERSION-2,LEVEL-1) MODIFIED ON JAN. 1983         / 5X,
     * 90HREFERENCE    : Y.BEPPU AND I.NINOMIYA;COMPUTER PHYSICS COMMUNI
     *CATIONS, V.23,Y.1981,P.123       /)
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.1) RETURN
      IF(NE.GT.N.OR.NV.LT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
      NEP1=NE+1
      NVM1=NV-1
      IF(INDEX.NE.0)                   GO TO 10
      SEPS=DSQRT(EPS*HUND)
      SEPS10=SEPS*TEN
      SSEPS=DSQRT(SEPS10)
      INDEX=1
   10 SEPS2=SEPS*( DABS(A(1,1)) + DABS(A(N,N)) )
      ITELIM=ITER
      ITERE=0
      IF(DABS(ESHIFT).LT.SEPS2)        GO TO 40
      DO 20 I=1,NE
      A(I,I)=A(I,I)-ESHIFT
   20 E(I)=E(I)-ESHIFT
      DO 30 I=NEP1,N
   30 A(I,I)=A(I,I)-ESHIFT
CC   TRANSFORMATION OF LARGE SQUARE-MATRIX-A TO SMALL SQUARE-MATRIX-V
CC  BY MATRIX-MULTIPLICATION ; FIRST U=A*V , THEN V=TRANS(V)*U .
   40 CONTINUE
      DO 80 I=1,NM1
      DO 50 K=1,I
   50 W1(K)=A(K,I)
      IP1=I+1
      DO 60 K=IP1,N
   60 W1(K)=A(I,K)
      DO 80 J=1,NV
      SUM=ZERO
      DO 70 K=1,N
   70 SUM=V(K,J)*W1(K)+SUM
      U(I,J)=SUM
   80 CONTINUE
      DO 100 J=1,NV
      SUM=ZERO
      DO 90 K=1,N
   90 SUM=A(K,N)*V(K,J)+SUM
      U(N,J)=SUM
  100 CONTINUE
      DO 130 J=1,NV
      DO 110 K=1,N
  110 W1(K)=V(K,J)
      DO 130 I=J,NV
      SUM=ZERO
      DO 120 K=1,N
  120 SUM=U(K,I)*W1(K)+SUM
      V(I,J)=SUM
  130 CONTINUE
CC  SOLVING EIGEN-PROBLEM IN SUB-SPACE
      ERROR=DABS(E(1)-V(1,1))
      E(1)=V(1,1)
      IF(NV.EQ.1)                      GO TO 630
      DO 140 K=2,NV
      ERROR=DABS(E(K)-V(K,K))+ERROR
  140 E(K)=V(K,K)
      E1ENE=DABS(E(1)) + DABS(E(NE))
      ERROR=ERROR/( E1ENE*DFLOAT(NV) )
C
      IF(ERROR.GT.SSEPS)               GO TO 260
C          BY APPROXIMATE FORMULA
      DEGENE=SEPS*E1ENE
      K=1
  150 CONTINUE
      ND=1
  160 KND=K+ND
      IF(KND.GT.NV)                    GO TO 170
      IF(DABS(E(KND)-E(K)).GT.DEGENE)  GO TO 170
      ND=ND+1
      GO TO 160
  170 IF(ND.EQ.1)                      GO TO 190
      DD=DFLOAT(ND)*SEPS10
      ND2=ND/2
      LI=K
      LF=K+ND2-1
      LM=ND-1
      DW=DFLOAT(LM)
      DO 180 L=LI,LF
      DDDW=DD*DW
      V(L,L)=   DDDW*V(L,L)
      LL=LI+LM
      V(LL,LL)=-DDDW*V(LL,LL)
      DW=DW-TWO
      LM=LM-1
  180 CONTINUE
  190 K=K+ND
      IF(K.LT.NV)                      GO TO 150
C
      DO 200 I=1,NVM1
      IP1=I+1
      DO 200 J=IP1,NV
      V(I,J)=V(J,I) / ( V(J,J)-V(I,I) )
  200 V(J,I)=-V(I,J)
C
      SUM=ZERO
      DO 210 J=2,NV
  210 SUM=V(J,1)*V(J,1)+SUM
      IF(SUM.GT.ONE)    GO TO 230
      S1=DSQRT(ONE-SUM)
      SUM=ZERO
      DO 220 J=1,NVM1
  220 SUM=V(J,NV)*V(J,NV)+SUM
      IF(SUM.GT.ONE)    GO TO 230
      SNV=DSQRT(ONE-SUM)
      VDIAGO=HALF*(S1+SNV)
      ILLIND=0
      GO TO 240
  230 VDIAGO=V0985
      ILLIND=200
  240 CONTINUE
      DO 250 K=1,NV
  250 V(K,K)=VDIAGO
      GO TO 510
C          BY HOUSEHOLDER-QR FORMULA
  260 NVM2=NV-2
      DO 270 J=2,NV
      J1=J-1
      DO 270 I=1,J1
  270 V(I,J)=V(J,I)
      IF(NV.LE.2) GO TO 360
      DO 350 K=1,NVM2
      K1=K+1
      E(K)=V(K,K)
      SUM=ZERO
      DO 280 J=K1,NV
      E(J)=V(K,J)
  280 SUM=E(J)*E(J)+SUM
      SUM=DSIGN(DSQRT(SUM),E(K1))
      W1(K)=-SUM
      E(K1)=E(K1)+SUM
      V(K,K1)=E(K1)
      C=E(K1)*SUM
      IF(C.LE.ZERO) GO TO 340
      SUMM=ZERO
      DO 320 I=K1,NV
      SUM=ZERO
      DO 290 J=K1,I
  290 SUM=V(J,I)*E(J)+SUM
      IF(I.GE.NV) GO TO 310
      IP1=I+1
      DO 300 J=IP1,NV
  300 SUM=V(I,J)*E(J)+SUM
  310 W1(I)=SUM/C
  320 SUMM=E(I)*W1(I)+SUMM
      W=SUMM*HALF/C
      DO 330 J=K1,NV
      W1(J)=E(J)*W-W1(J)
      DO 330 I=K1,J
  330 V(I,J)=E(J)*W1(I)+E(I)*W1(J)+V(I,J)
  340 V(K,K)=C
  350 CONTINUE
  360 E(NV)=V(NV,NV)
      V(NV,NV)=ONE
      IF(NV.EQ.1) GO TO 510
      E(NVM1)=V(NVM1,NVM1)
      W1(NVM1)=V(NVM1,NV)
      W1(NV)=ZERO
      V(NVM1,NVM1)=ONE
      V(NVM1,NV)=ZERO
      V(NV,NVM1)=ZERO
      IF(NV.EQ.2) GO TO 420
      DO 410 L=1,NVM2
      K=NVM1-L
      K1=K+1
      C=-V(K,K)
      V(K,K)=ONE
      IF(C.GE.ZERO) GO TO 390
      DO 380 J=K1,NV
      SUM=ZERO
      DO 370 I=K1,NV
  370 SUM=V(I,J)*V(K,I)+SUM
      S=SUM/C
      DO 380 I=K1,NV
  380 V(I,J)=S*V(K,I)+V(I,J)
  390 CONTINUE
      DO 400 I=K1,NV
      V(K,I)=ZERO
  400 V(I,K)=ZERO
  410 CONTINUE
  420 D=DABS(E(1))
      DO 430 J=2,NV
      IF(( DABS(W1(J-1))+DABS(E(J)) ).GT.D) D=DABS(W1(J-1))+DABS(E(J))
  430 CONTINUE
      IF(D.EQ.ZERO) GO TO 510
      D=D*EPS
      K=NV
  440 L=K
  450 IF(DABS(W1(L-1)).LT.D) GO TO 460
      L=L-1
      IF(L.GT.1) GO TO 450
  460 IF(L.EQ.K) GO TO 500
      WWW=(E(K-1)+E(K))*HALF
      R=E(K)-WWW
      Z=WWW-DSIGN(DSQRT(W1(K-1)*W1(K-1)+R*R),WWW)
      EE=E(L)-Z
      E(L)=EE
      FF=W1(L)
      R=DSQRT(EE*EE+FF*FF)
      J=L
      GO TO 480
  470 R=DSQRT(E(J)*E(J)+W1(J)*W1(J))
      W1(J-1)=S*R
      EE=C*E(J)
      FF=C*W1(J)
  480 C=E(J)/R
      S=W1(J)/R
      WWW=E(J+1)-Z
      E(J)=(FF*C+WWW*S)*S+EE+Z
      E(J+1)=C*WWW-S*FF
      DO 490 I=1,NV
      R=V(I,J+1)
      V(I,J+1)=R*C-V(I,J)*S
  490 V(I,J)=V(I,J)*C+R*S
      J=J+1
      IF(J.LT.K) GO TO 470
      W1(K-1)=E(K)*S
      E(K)=E(K)*C+Z
      GO TO 440
  500 K=K-1
      IF(K.GT.1) GO TO 440
C        STRAIGHT INSERTION SORT OF EIGENVALUES
  510 J=NV
  520 CONTINUE
      DO 530 L=1,NV
  530 W1(L)=V(L,J-1)
      EK=E(J-1)
      DO 550 I=J,NV
      IF(DABS(E(I)).LE.DABS(EK)) GO TO 570
      DO 540 L=1,NV
  540 V(L,I-1)=V(L,I)
  550 E(I-1)=E(I)
      DO 560 L=1,NV
  560 V(L,NV)=W1(L)
      E(NV)=EK
      GO TO 590
  570 CONTINUE
      DO 580 L=1,NV
  580 V(L,I-1)=W1(L)
      E(I-1)=EK
  590 J=J-1
      IF(J.GT.1) GO TO 520
C
CC  NORM-ORTHOGONALIZATION OF EIGEN-VECTORS
      DO 620 I=1,N
      DO 600 K=1,NV
  600 W1(K)=U(I,K)
      DO 620 J=1,NV
      SUM=ZERO
      DO 610 K=1,NV
  610 SUM=V(K,J)*W1(K)+SUM
      U(I,J)=SUM
  620 CONTINUE
C
  630 CONTINUE
      DO 660 J=1,NV
      DO 640 K=1,N
  640 W1(K)=U(K,J)
      DO 660 I=1,J
      SUM=ZERO
      DO 650 K=1,N
  650 SUM=U(K,I)*W1(K)+SUM
      V(I,J)=SUM
  660 CONTINUE
C
      IF(V(1,1).LT.ZERO)      ILLIND=200
      IF(V(NV,NV).LT.ZERO)    ILLIND=200
      IF(V(1,1).LT.V(NV,NV))  ILLIND=200
      W2(1)=ONE/DSQRT(V(1,1))
      IF(NV.EQ.1)             GO TO 710
      V(1,2)=V(1,2)*W2(1)
      W2(2)=ONE/DSQRT(  V(2,2)-V(1,2)**2  )
      IF(NV.EQ.2)             GO TO 710
      DO 700 J=3,NV
      V(1,J)=V(1,J)*W2(1)
      JM1=J-1
      DO 680 I=2,JM1
      SUM=ZERO
      IM1=I-1
      DO 670 L=1,IM1
  670 SUM=V(L,I)*V(L,J)+SUM
      V(I,J)=(  V(I,J)-SUM  )*W2(I)
  680 CONTINUE
      SUM=ZERO
      DO 690 L=1,JM1
  690 SUM=V(L,J)*V(L,J)+SUM
      W2(J)=ONE/DSQRT(  V(J,J)-SUM  )
  700 CONTINUE
C
  710 CONTINUE
      C11=W2(1)
      DO 720 J=1,N
  720 V(J,1)=U(J,1)*C11
      IF(NV.EQ.1)  GO TO 760
      DO 750 I=2,NV
      CII=W2(I)
      IM1=I-1
      DO 730 K=1,IM1
  730 W1(K)=V(K,I)
      DO 750 J=1,N
      SUM=ZERO
      DO 740 K=1,IM1
  740 SUM=V(J,K)*W1(K)+SUM
      V(J,I)=(U(J,I)-SUM)*CII
  750 CONTINUE
  760 ITERE=ITERE+1
      IF(ITERE.LT.ITELIM .AND. ERROR.GT.SEPS10)  GO TO 40
      IF(ITERE.GE.ITELIM .AND. ERROR.GT.SEPS10)  GO TO 770
C
      ILL=0
      IF(ILLIND.EQ.200)                          GO TO 770
      GO TO 780
  770 ILL=200
  780 CONTINUE
      ITER=ITERE
      IF(DABS(ESHIFT).LT.SEPS2) RETURN
      DO 790 I=1,NE
      A(I,I)=A(I,I)+ESHIFT
  790 E(I)=E(I)+ESHIFT
      DO 800 I=NEP1,N
  800 A(I,I)=A(I,I)+ESHIFT
      RETURN
      END
      SUBROUTINE NGHOUD(AB,NNMAX,NN,NNE,NNV,EPS,IORD,ICHO,BD,E,V,ILL,
     *                  W1,W2,W3,W4,W5,W6,W7)
CI--------------------------------------------------------------------I
CI       SUBPROGRAM FOR GENERALIZED EIGEN-PROBLEM , A*V=B*V*E , BY    I
CI  ( SIMULTANEOUSLY TRIANGULAR DECOMPOSITION )-                      I
CI  (HOUSEHOLDER)-(BISECTION & NO-ROOT-QR)-(INVERSE-ITERATION) METHOD I
CI-----------------------------------------------(VERSION-2,LEVEL-1)--I
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AB(NNMAX,NN),BD(NN),V(NNMAX,NNV),E(NNE)
      DIMENSION W1(NN),W2(NN),W3(NN),W4(NN),W5(NN),W6(NN),W7(1)
C     W7(I) IS NEVER USED IN VERSION-2.
      DATA ZERO,ONE /0.D 0,1.D 0/
      NMAX=NNMAX
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.0) RETURN
      IF(NE.GT.N.OR.NV.GT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
      DO 10 I=1,NM1
      IP1=I+1
      DO 10 J=IP1,N
      T=AB(J,I)
      AB(J,I)=AB(I,J)
   10 AB(I,J)=T
CC   CHOLESKY TRANSFORMATION OF A BY DECOMPOSED B
      IF(ICHO.EQ.1)          GO TO 60
      ILL=100
      IF(BD(1).LT.EPS)       RETURN
      BD(1)=ONE/DSQRT(BD(1))
      AB(1,2)=AB(1,2)*BD(1)
      PPP=BD(2)-AB(1,2)*AB(1,2)
      IF(PPP.LT.EPS)         RETURN
      BD(2)=ONE/DSQRT(PPP)
      IF(N.EQ.2)             GO TO 60
      DO 50 J=3,N
      AB(1,J)=AB(1,J)*BD(1)
      JM1=J-1
      DO 30 I=2,JM1
      SUM=ZERO
      IM1=I-1
      DO 20 L=1,IM1
   20 SUM=AB(L,I)*AB(L,J)+SUM
      AB(I,J)=(  AB(I,J)-SUM  )*BD(I)
   30 CONTINUE
      SUM=ZERO
      DO 40 L=1,JM1
   40 SUM=AB(L,J)*AB(L,J)+SUM
      PPP=BD(J)-SUM
      IF(PPP.LT.EPS)         RETURN
      BD(J)=ONE/DSQRT(PPP)
   50 CONTINUE
C
   60 CONTINUE
      DO 70 J=1,N
   70 AB(J,1)=AB(J,1)*BD(1)
      DO 90 I=2,N
      BDI=BD(I)
      IM1=I-1
      DO 90 J=I,N
      SUM=ZERO
      DO 80 K=1,IM1
   80 SUM=AB(K,I)*AB(J,K)+SUM
      AB(J,I)=(AB(J,I)-SUM)*BDI
   90 CONTINUE
C
      AB(1,1)=AB(1,1)*BD(1)
      DO 110 I=2,N
      SUM=ZERO
      IM1=I-1
      DO 100 K=1,IM1
  100 SUM=AB(K,1)*AB(K,I)+SUM
      AB(I,1)=(AB(I,1)-SUM)*BD(I)
  110 CONTINUE
C
      DO 150 J=2,NM1
      JM1=J-1
      JP1=J+1
      SUM=ZERO
      DO 120 K=1,JM1
  120 SUM=AB(K,J)*AB(J,K)+SUM
      AB(J,J)=(AB(J,J)-SUM)*BD(J)
      DO 150 I=JP1,N
      SUM=ZERO
      DO 130 K=1,JM1
  130 SUM=AB(K,I)*AB(J,K)+SUM
      IM1=I-1
      DO 140 K=J,IM1
  140 SUM=AB(K,J)*AB(K,I)+SUM
      AB(I,J)=(AB(I,J)-SUM)*BD(I)
  150 CONTINUE
      SUM=ZERO
      DO 160 K=1,NM1
  160 SUM=AB(K,N)*AB(N,K)+SUM
      AB(N,N)=(AB(N,N)-SUM)*BD(N)
C#  SOLVING STANDARD EIGEN-PROBLEM FOR TRANSFORMED A
      DO 170 I=1,NM1
      IP1=I+1
      DO 170 J=IP1,N
      T=AB(J,I)
      AB(J,I)=AB(I,J)
  170 AB(I,J)=T
      CALL NSHOUD(AB,NMAX,N,NE,NV,EPS,IORD,E,V,ILL,W1,W2,W3,W4,W5,W6,W7)
      IF(NV.EQ.0.OR.ILL.NE.0) RETURN
C#  CONVERSION OF EIGEN-VECTORS( STANDARD FORM -> GENERALIZED FORM )
      DO 190 I=1,NV
      V(N,I)=V(N,I)*BD(N)
      DO 190 JJ=1,NM1
      J=N-JJ
      SUM=ZERO
      JP1=J+1
      DO 180 K=JP1,N
  180 SUM=AB(K,J)*V(K,I)+SUM
      V(J,I)=(V(J,I)-SUM)*BD(J)
  190 CONTINUE
      RETURN
      END
      SUBROUTINE NGJEND(AB,NNMAX,NN,NNE,NNV,EPS,BD,IUV,ITER,ESHIFT,
     *                  E,V,U,ILL,W1,W2)
CI--------------------------------------------------------------------I
CI       SUBPROGRAM FOR GENERALIZED EIGEN-PROBLEM , A*V=B*V*E , BY    I
CI      (SIMULTANEOUSLY TRIANGULAR DECOMPOSITION)-(JENNINGS) METHOD   I
CI-----------------------------------------------(VERSION-2,LEVEL-1)--I
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AB(NNMAX,NN),BD(NN),E(NNE),V(NNMAX,NNV),U(NNMAX,NNV)
      DIMENSION W1(NN),W2(NN)
      DATA ZERO,ONE /0.D 0,1.D 0/
      NMAX=NNMAX
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.1) RETURN
      IF(NE.GT.N.OR.NV.LT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
C   CHOLESKY TRASFORMATION OF A BY DECOMPOSED B
      DO 10 J=1,N
   10 AB(1,J)=AB(1,J)*BD(1)
      DO 30 I=2,N
      BDI=BD(I)
      IM1=I-1
      DO 30 J=I,N
      SUM=ZERO
      DO 20 K=1,IM1
   20 SUM=AB(I,K)*AB(K,J)+SUM
      AB(I,J)=(AB(I,J)-SUM)*BDI
   30 CONTINUE
C
      AB(1,1)=AB(1,1)*BD(1)
      DO 50 I=2,N
      SUM=ZERO
      IM1=I-1
      DO 40 K=1,IM1
   40 SUM=AB(1,K)*AB(I,K)+SUM
      AB(1,I)=(AB(1,I)-SUM)*BD(I)
   50 CONTINUE
C
      DO 90 J=2,NM1
      JM1=J-1
      JP1=J+1
      SUM=ZERO
      DO 60 K=1,JM1
   60 SUM=AB(K,J)*AB(J,K)+SUM
      AB(J,J)=(AB(J,J)-SUM)*BD(J)
      DO 90 I=JP1,N
      SUM=ZERO
      DO 70 K=1,JM1
   70 SUM=AB(K,J)*AB(I,K)+SUM
      IM1=I-1
      DO 80 K=J,IM1
   80 SUM=AB(J,K)*AB(I,K)+SUM
      AB(J,I)=(AB(J,I)-SUM)*BD(I)
   90 CONTINUE
      SUM=ZERO
      DO 100 K=1,NM1
  100 SUM=AB(K,N)*AB(N,K)+SUM
      AB(N,N)=(AB(N,N)-SUM)*BD(N)
C   TRANSFORMATION OF INITIAL EIGEN-VECTORS
      IF(IUV.EQ.1)           GO TO 150
      DO 130 I=1,NM1
      W2(I)=ONE/BD(I)
      IP1=I+1
      DO 110 K=IP1,N
  110 W2(K)=AB(K,I)
      DO 130 J=1,NV
      SUM=ZERO
      DO 120 K=I,N
  120 SUM=V(K,J)*W2(K)+SUM
      U(I,J)=SUM
  130 CONTINUE
      BNN=ONE/BD(N)
      DO 140 J=1,NV
  140 U(N,J)=V(N,J)*BNN
C#  SOLVING STANDARD EIGEN-PROBLEM FOR TRANSFORMED A
  150 CALL NSJEND(AB,NMAX,N,NE,NV,EPS,ITER,ESHIFT,E,U,ILLE,W1,W2,V)
C#  CONVERSION OF EIGEN-VECTORS( STANDARD FORM -> GENERALIZED FORM )
      DO 170 I=1,NV
      V(N,I)=U(N,I)*BD(N)
      DO 170 JJ=1,NM1
      J=N-JJ
      SUM=ZERO
      JP1=J+1
      DO 160 K=JP1,N
  160 SUM=AB(K,J)*V(K,I)+SUM
      V(J,I)=(U(J,I)-SUM)*BD(J)
  170 CONTINUE
      ILL=ILLE
      RETURN
      END
