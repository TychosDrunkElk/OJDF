
C
C       *********************************************************
C       *****   OPPOSED-JET DIFFUSION FLAME (AXISYMMETRIC)  *****
C       *********************************************************
C
        PROGRAM OJDFS2
C         (VERSION -S2)
C
C       ----------------------------------------------------------
C       |                                                        |
C       |  PROGRAM OJDF WAS DEVELOPED BY JAMES S. T'IEN TO MODEL |
C       |  THE DIFFUSION FLAME IN THE STAGNATION POINT REGION OF |
C       |  OPPOSING OXIDIZER AND FUEL FLOWS.  THIS VERSION, -S2, |
C       |  IS FOR A SOLID FUEL FACING A GASEOUS OXIDIZER STREAM. |
C       |  THIS VERSION IS MODIFIED FROM A PREVIOUS VERSION -S1  |
C       |  BY INCLUDING ARREHNIUS PYROLYSIS LAW AT THE SOLID     |
C       |  SURFACE AND SPECIFICALLY DESIGNED TO INVESTIGATE THE  |
C       |  EFFECT OF RADIATIVE LOSS AT SMALL STRETCH RATE ON     |
C       |  FLAME EXTINCTION.                                     |
C       |                                                        |
C       ----------------------------------------------------------
C
        DIMENSION T(25),F(25),U(25),YF(25),YO(25),W(25)
        DIMENSION TN(25),FN(25),UN(25),YFN(25),YON(25),WN(25)
C
        DATA CSCP,  XL,  XNO  /1.32,4.32,-1.92/
        DATA TW,TCE/2.5,1.0/
        DATA M/25/
        DATA SC,PR/0.7,0.7/
        DATA RO,XMU/1.176E-03,1.85E-04/
        DATA EW,SB/50.3,2.32E+06/
        DATA B/5.27E+07/
C	B=frequency factor in reaction in gas reaction ????
C	F? what is F??????????????
C       RO=AMBIENT DENSITY (G/CM3)
C       XMU=AMBIENT VISCOSITY (G/CM-SEC)
C       SB=PRE-EXPONENTIAL FACTOR, PYROLYSIS LAW (G/CM2-SEC)
C       EW=ACTIVATION ENERGY ,     PYROLYSIS LAW (NONDIM.)
C
        DATA T/2.5,3.3,4.1,4.9,5.8,6.6,7.1,6.7,5.4,4.2,3.2,2.5,2.0,
     1         1.6,1.4,1.3,1.2,1.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0/
        DATA F/-1.1,-.98,-.74,-.36,.13,.69,1.3,1.9,2.5,3.0,3.5,3.9,
     1         4.3,4.6,4.9,5.3,5.6,5.9,6.2,6.5,6.8,7.1,7.4,7.7,8.0/
        DATA YF/.5125,.4552,.3920,.3240,.2536,.1850,.1243,.0807,
     1          .0565,.0402,.0279,.0190,.0126,.0082,.0052,.0033,
     1          .0020,.0012,.0007,.0004,.0002,.0001,.0001,.0,.0/
        DATA YO/.0014,.0016,.0017,.0020,.0026,.0050,.0136,.0400,
     1          .0860,.1270,.1590,.1827,.1994,.2110,.2188,.2240,
     1          .2272,.2293,.2306,.2314,.2314,.2314,.2314,.2314,.2314/
C
  002   FORMAT (I4)
  003   FORMAT ('  INPUT NO. OF TIME STEPS, I4')
        WRITE (*,003)
        READ (*,002) N
  005   FORMAT ('  INPUT dt, dy')
        WRITE (*, 005)
        READ (*,*) DT, DY
  006   FORMAT ('  INPUT HEAT OF COMB. Q, OXYGEN MASS FRAC, YOE')
        WRITE (*,006)
        READ(*,*) Q, YOE
  007   FORMAT ('  INPUT ACTI.ENERGY E,  VELOSITY GRADIENT (A)')
        WRITE (*,007)
        READ (*,*) E, A
  008   FORMAT ('  INPUT SURF. RADIA. EMMISIVITY,RAD. ABSORP. QA')
        WRITE (*,008)
        READ (*,*) EPSONS, QA
C
        M1=M-1
        M2=M-2
        N1=N-1
C
C       ---------- INITIAL CONDITIONS ---------
C
        U(1) =0.
        U(M) =1.0
        W(1) =D*YF(1)*YO(1)/EXP(E/T(1))
C
        DO 55 I=2,M1
   55   U(I) =(F(I+1)-F(I-1))/4./DY
        DA=DT/DY
        DB=DA/DY
C
        KOUNT=1
  101   CONTINUE
  
        TN(M) =1.
        YON(M) = YOE
        YFN(M) =0.
        UN(M) =1.0
        D =B/A
        S =0.2304*EPSONS/SQRT(A)
C       THE ABOVE NUMERICAL VALUE IS FOR TE=300K
        BETA =SB/SQRT(RO*XMU*A)
C
C       ----------- DIFFERENTIAL EQUATIONS ----------
        DO 32 I=2,M1
        XI=I-1
        ACT =E/T(I)
        ARH=EXP(-1*ACT)
   36   W(I) =D*YF(I)*YO(I)*ARH
C
C       UPWIND SCHEME FOR T,YO,YF,U IN CONVECTIVE TERM
        IF (F(I).GE.0.0) GO TO 61
        TL =T(I-1)
        TR =T(I)
        YOL =YO(I-1)
        YOR =YO(I)
        YFL =YF(I-1)
        YFR =YF(I)
        UL =U(I-1)
        UR =U(I)
        GO TO 62
   61   CONTINUE
        TL =T(I)
        TR =T(I+1)
        YOL=YO(I)
        YOR=YO(I+1)
        YFL=YF(I)
        YFR=YF(I+1)
        UL =U(I)
        UR =U(I+1)
   62   CONTINUE
        TN(I) = T(I) +DB*(T(I+1)-2.*T(I)+T(I-1))/PR+DA*F(I)*(TR-TL)
     1    +DT*Q*W(I)
        YFN(I)=(YF(I)+DB*(YF(I+1)-2.*YF(I)+YF(I-1))/SC +DA*F(I)*
     1    (YFR-YFL)) /(1.0+ DT*D*YO(I)*ARH)
        YON(I)=(YO(I)+DB*(YO(I+1)-2.*YO(I)+YO(I-1))/SC +DA*F(I)*
     1    (YOR-YOL)) /(1.0-DT*XNO*D*YF(I)*ARH)
        UN(I) =U(I) +DA*F(I)*(UR-UL) +DB* (U(I+1)-2.*U(I)+U(I-1))
     1    +    DT*(T(I)-U(I)*U(I))
   32   CONTINUE
C
C
C       -------- BOUNDARY CONDITION AT SURFACE ----------
        TN(1)=-EW/ALOG(-F(1)/BETA)
        
        ACT1=E/TN(1)
        ARH1 =1./EXP(ACT1)
        YF0J=F(1)*2.0*DY*(YF(1)-1.0)*SC +YF(2)
        YFN(1)=(YF(1)+DB*(YF(2)-2.*YF(1)+YF0J)/SC+.5*DA*(YF(2)
     1     -YF0J)*F(1)) /(1.0 +DT*D*YO(1)*ARH1)
        YO0J=F(1)*2.0*DY*YO(1)*SC +YO(2)
        YON(1)=(YO(1)+DB*(YO(2)-2.*YO(1)+YO0J)/SC+.5*DA*(YO(2)
     1     -YO0J)*F(1))/(1.0-XNO*DT*D*YF(1)*ARH1)
        WN(1) = D*YFN(1)*YON(1)*ARH1
C
C       -------- TO OBTAIN FW ---------
        TG = TN(1) -CSCP*TCE +XL
        TIMAG = (-TN(2)+2.*TN(1) -.5*F(1)*TN(2)*PR*DY-Q*WN(1)*DY*DY*PR) / (1.-.5*F(1)*PR*DY)
        FN(1) = ((TN(2)-TIMAG)/2./DY - S*TN(1)**4. +QA)/(-TG)/PR
C
C       ------- INTEGRATION OF VELOCITY TO OBTAIN F -------
        FN(2) =FN(1)+DY*(UN(2)+UN(1))
        DO 51 I =1,M2
        FN(I+2) =FN(I) +2.*DY*(UN(I+2) +4.*UN(I+1) +UN(I))/3.0
  51    CONTINUE
C
        IF(TN(2).LT.TW) GO TO 99
C
        DO 88 I=1,M
        T(I)=TN(I)
        YF(I)=YFN(I)
        YO(I)=YON(I)
        U (I)=UN (I)
        F (I)=FN (I)
  88    CONTINUE
C
C       ------- PRINT OUT -------
        IF (MOD(KOUNT,50).NE.0) GO TO 31
C so every 50 kounts, if Kount lt N, repeat calculations
        XJ=KOUNT
        TAU=XJ*DT
        WRITE(*,11) TAU,A,F(1)
  11    FORMAT('  TAU =',F10.5,5X,'A =',F12.5,5X,'FW =',F10.5)
        WRITE(*,12) D,S,BETA
  12    FORMAT('  D =',E12.7,5X,'S =',F10.5,5X,'BETA =',E12.7)
        WRITE (*,21) (T (I), I=1,M)
        WRITE (*,23) (YF(I), I=1,M)
        WRITE (*,24) (YO(I), I=1,M)
        WRITE (*,25) (W (I), I=1,M)
        WRITE (*,26) (U (I), I=1,M)
        WRITE (*,20)
  20    FORMAT(' ')
  21    FORMAT(' T= ',25f12.5)
  23	FORMAT(' YF= ',25f12.5)
  24	FORMAT(' YO= ',25f12.5)
  25	FORMAT(' W= ',25f12.5)
  26	FORMAT(' U= ',25f12.5)

  22    FORMAT(///)
        print *,"  number   "
C
  31    CONTINUE
        KOUNT =KOUNT +1
        IF (KOUNT.GT.N) GO TO 99
        GO TO 101
  99    CONTINUE
        STOP
        END
