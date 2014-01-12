C                                                                       
C       *********************************************************       
C       *****   OPPOSED-JET DIFFUSION FLAME (2-DIMENSIONAL) *****       
C       *********************************************************       
C                                                                       
C      PROGRAM OJDFS2                                                   
C           (VERSION -S2)                                               
C                                                                       
C       ----------------------------------------------------------      
C       !                                                        !      
C       !  PROGRAM OJDF WAS DEVELOPED BY JAMES S. T'IEN TO MODEL !      
C       !  THE DIFFUSION FLAME IN THE STAGNATION POINT REGION OF !      
C       !  OPPOSING OXIDIZER AND FUEL FLOWS.  THE CONFIGURATION  !      
C       !  IS FOR A  GAS  FUEL FACING A GASEOUS OXIDIZER STREAM. !      
C       !  THIS VERSION INCLUDE CONTROLLED FUEL EJECTION VELOCITY!      
C       !  SPECIFICALLY  IT  IS     DESIGNED TO INVESTIGATE THE  !      
C       !  EFFECT OF RADIATIVE LOSS AT SMALL STRETCH RATE ON     !      
C       !  FLAME EXTINCTION.  A PRINTOUT OF FLAME STRUCTURE AT   !      
C       !  EVERY SO MANY ITERATIONS.                             !      
C       !               NOVEMBER, 1985                           !      
C       ----------------------------------------------------------      
C********************************************************************   
C                                                                       
C			THE MODIFICATION OF THE SIMULATION              
C        FOR THE GAS FUEL DIFFUSION FLAME WITH MIXED                    
C	CONVECTION FOR 2-D AND AXISYMMETRIC TSUJI-TYPE BURNER           
C                                                                       
C      1. To modify the original version in order to                    
C         run it in the UNIX machine;                                   
C	 2. To make the program work for both 2-D and axisymmetric      
C		geometries, the ONLY cotrolling parameter is CONFI;     
C	 3. To add the gravity force in the momentum equation;          
C	 4. To introduce the mixed convection and mixed stretch rate    
C		with the Froude number (TSAI);                          
C	 5. To change the boundary condition to include the buoyancy    
C		effects;                                                
C                                                                       
C	                          Bai Han (02/22/02)                    
C********************************************************************   
      PROGRAM TIEN                                                      
C*****precision > double                                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)                
C*****END precision > double                                            
C*****precision > single                                                
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)                           
C*****END precision > single                                            
	CHARACTER G                                                     
C                                                                       
      DIMENSION T(40),F(40),U(40),YF(40),YO(40),W(40)                   
      DIMENSION TN(40),FN(40),UN(40),YFN(40),YON(40),WN(40)             
C                                                                       
      DATA  XNO  /-4.00/                                                
      DATA M,N/40,99999/                                                
      DATA SCF,SCO,PR/0.56,0.56,0.70/                                   
      DATA RO,XMU/1.177E-03,1.846E-04/                                  
      DATA B,E/1.801E+09,54.19/                                         
      DATA XNUM/1.168/                                                  
      DATA Q,YOE/122.0,0.2325/                                          
      DATA EPSONS/0.25/                                                 
	DATA TREF/4.5/                                                  
      DATA DT1,DT2,DY/0.0005,0.0015,0.2000/                             
C ****************************************************************      
C  DT1 is the initial time step (For first 100 iterations)              
C  DT2 is the final time step                                           
C  DY is the spatial grid size of Y, it is usually about 0.15           
C ****************************************************************      
C       Q=HEAT OF COMBUSTION (DIMENSIONLESS)                            
C       YOE=MASS FRACTION OF O2 IN ATMOSPHERE (DIMENSIONLESS)           
C       EPSONS=SURFACE EMMISIVITY (DIMENSIONLESS)                       
C       RO=AMBIENT DENSITY (G/CM3)                                      
C       XMU=AMBIENT VISCOSITY (G/CM-SEC)                                
C       XNUM=KINEMATIC VISCOSITY EVALUATED AT MEAN TEMPERATURE(CM2/S)   
C       E=ACTIVATION ENERGY IN GAS REACTION (DIMENSIONLESS)             
C       B=FREQUENCY FACTOR  IN GAS REACTION (DIMENSIONLESS)             
C       XMDOT=MASS FLOW RATE OF FUEL (G/CM-SEC)                         
C                                                                       
      DATA T/2.5,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.3,5.7,6.3,6.3,       
     1   7.2,6.4,4.2,3.2,2.5,1.9,1.6,1.4,1.3,1.2,1.1,1.,1.,1.,1.,1.,1., 
     2   1.,1.,1.,1.,1.,1.,1.,1.,1.,1./                                 
      DATA F/-2.00,-1.62,-1.5,-1.38,-1.10,-.93,-.76,-.54,-.32,-.15,     
     1  0.03,.2,.45,.65,1.1,1.8,2.5,3.,3.5,4.3,4.7,5.3,5.9,6.6,7.1,8.,  
     2  8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8.,8./                      
      DATA YF/.9982,.9955,.984,.9433,.8947,.8327,.7711,.7167,.6050,     
     1           .5058,.4127,.3063,.2010,.1209,.0926,.0645,             
     2       .035,.011,.0015,.0005,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 
     3       0.,0.,0.,0.,0.,0.,0.,0./                                   
      DATA YO/.0002,.0004,.0005,.0008,.0016,.0023,.0034,.0048,.0058,    
     1          .0069,.0085,.01,.012,.016,.032,.04,.05,.063,.08,        
     2   .12,.1333,.178,.2,.2222,.2325,.2325,.2325,.2325,.2325,.2325,   
     3     .2325,.2325,.2325,.2325,.2325,                               
     4     .2325,.2325,.2325,.2325,.2325/                               
C                                                                       
      PARAMETER (LOUT1=4, LOUT=5)                                       
C *******************************************************************   
C      PARAMETER FOR THE CONFIGAURATION:                                
C	   2-D			: CONFI=0.                              
C        axisymmetric	: CONFI=1.                                      
C *******************************************************************	
C	PARAMETER (CONFI=0.)                                            
C                                                                       
C            open the user's file                                       
                                                                        
C*****unix                                                              
C      OPEN (LIN, FORM='FORMATTED', FILE='input')                       
      OPEN (LOUT, FORM='FORMATTED', FILE='output')                      
      OPEN (LOUT1, FORM='FORMATTED', FILE='out1')                       
C                                                                       
C*****unix                                                              
      WRITE (LOUT,601)                                                  
      WRITE (LOUT,*) '####PLEASE CHOOSE THE CONFIGAURATION:'            
	WRITE (LOUT,602)                                                
C*****PC                                                                
	WRITE (*,601)                                                   
      WRITE (*,*) '####PLEASE CHOOSE THE CONFIGAURATION:'               
	WRITE (*,602)                                                   
                                                                        
      READ(*,'(I3)') NCONFI                                             
	IF (NCONFI.EQ.1) CONFI=1.                                       
	IF (NCONFI.EQ.0) CONFI=0.                                                                                                              
C*****unix                                                              
      WRITE (LOUT,*) '####DO YOU WANT TO INCLUDE THE GRAVITY FORCE? '   
	WRITE (LOUT,603)                                                
	WRITE (LOUT,604)                                                
                                                                        
C*****PC                                                                
      WRITE (*,*) '####DO YOU WANT TO INCLUDE THE GRAVITY FORCE? '      
	WRITE (*,603)                                                   
                                                                        
      READ(*,'(A1)') 	G                                               
                                                                        
	WRITE (*,604)                                                   
	READ (*,*) R,UOXY,FW                                            
                                                                        
C                                                                       
C007   FORMAT (' INPUT VEL. GRAD. (A), FUEL FEED RATE (FW):2F6.3')      
C      WRITE (*,007)                                                    
C      READ (*,*) A,FW                                                  
C       AF = 1.98                                                       
C ********************************************************************  
                                                                        
C Note: A is assigned directly in the original version!                 
C       Here the Af is calculated by the potential flow theory.         
C  The velocity of the forced air flow                                  
C	UOXY = 3.   ! (m/s)                                             
C  The radius of the cylinder (2D) or sphere (Axisym) is                
C	R = 2.5     ! (m)                                               
C  The forced flow velocity gradient is given by                        
   	 AF = (4.-CONFI)/2.*UOXY/R                                      
C The stretch rate induced by the natural convection                    
	AN =SQRT((1.176-.261)/1.176*(9.8/R))                            
C      AN=0.0                                                           
C      TSAI=0.0                                                         
C For purely forced flow. We can treat g=0. instead of 9.8!             
                                                                        
C ********************************************************************  
C  The nondimesional fuel injection rate:                               
C       FW =-10.                                                        
C ********************************************************************* 
C  To introduce the mixed convection and mixed stretch rate             
C		with the Froude number (TSAI);                          
C      Forced flow: TSAI=0. and AF=A                                    
C	 Natural convection: TSAI>10000, and AF=0.                      
C ********************************************************************* 
C	TSAI= 0.                                                        
C	A   =AF*SQRT(1.+TSAI)                                           
C	IF (AF.EQ.0..AND.TSAI.GE.10000.) THEN                           
C	A  =SQRT((1.176-.261)/1.176*(9.8/R))                            
C      ENDIF                                                            
C ********************************************************************  
C Note:                                                                 
C 1. The ambient density is assumed as the same as ambient air;         
C        1.176 kg/m3                                                    
C 2. The characteristic density of the boundary layer is assumed the    
C    same as the air at 1400K (0.261 kg/m3)                             
C 3. TSAI is the square of the ratio between AN and AF                  
C                                                   Bai Han (03/01/02)  
C ********************************************************************  
      IF (AF.EQ.0.)	THEN                                            
	A  = AN                                                         
	TSAI=50000.                                                     
	ELSE IF (G.EQ.'n'.OR.G.EQ.'N') THEN                             
	A  =AF                                                          
	TSAI=0.                                                         
	ELSE                                                            
	TSAI=(AN/AF)*(AN/AF)                                            
      A  =AF*SQRT(1.+TSAI)                                              
	END IF                                                          
                                                                        
013   FORMAT(F6.2,1X,F7.4)                                              
C                                                                       
      M1=M-1                                                            
      M2=M-2                                                            
      N1=N-1                                                            
C                                                                       
559   FORMAT(' -----------------------------------------------------')  
      WRITE(4,559)                                                      
C560   FORMAT(' A, FW, Q, EPSONS, N =',4(3X,F8.4),3X,I5)                
C      WRITE (4,560) A,FW,Q,EPSONS,N                                    
C      WRITE(5,560) A,FW,Q,EPSONS,N                                     
560   FORMAT(' Af, Tsai, A, FW, Q, EPSONS, N =',6(3X,F8.4),3X,I5)       
      WRITE (4,560) AF,TSAI,A,FW,Q,EPSONS,N                             
      WRITE(5,560) AF,TSAI,A,FW,Q,EPSONS,N                              
      WRITE(4,23) DT1,DT2,DY                                            
561   FORMAT(' SCF,SCO,PR =',3(5X,F8.4))                                
      WRITE(4,561) SCF,SCO,PR                                           
      WRITE(4,559)                                                      
      WRITE(4,558)                                                      
C                                                                       
C       ---------- INITIAL CONDITIONS ---------                         
C                                                                       
      U(1) =0.                                                          
C      U(M) =1.0                                                        
C To change the boundary condition to include the buoyancy effects;     
	U(M) =1.0/SQRT(1.+TSAI)                                         
                                                                        
      T(1)= 1. +0.5/SQRT(-FW)                                           
C      XMDOT = -FW*SQRT(RO*XMU*A)                                       
C *******************************************************************   
C     The Mass Flow Rate of Fuel (XMDOT):                               
C	   2-D			: CONFI=0.                              
C        axisymmetric	: CONFI=1.                                      
C *******************************************************************	
	XMDOT = -FW*SQRT(RO*XMU*A)*SQRT((1.+CONFI)/2.)                  
      TMAX1 = 1.0                                                       
C                                                                       
      DO 55 I=2,MAX(2,M1)                                               
55    U(I) =(F(I+1)-F(I-1))/4./DY                                       
C                                                                       
      KOUNT=1                                                           
101   CONTINUE                                                          
      IF(KOUNT.LT.100) DT=DT1                                           
      IF(KOUNT.GE.100) DT=DT2                                           
      DA=DT/DY                                                          
      DB=DA/DY                                                          
      TN(M) =1.                                                         
      YON(M) = YOE                                                      
      YFN(M) =0.                                                        
C      UN(M) =1.0                                                       
C To change the boundary condition to include the buoyancy effects;     
	UN(M) =1.0/SQRT(1.+TSAI)                                        
                                                                        
C      D =B/A                                                           
C *******************************************************************   
C     The Damkohler Number (D):                                         
C	   2-D			: CONFI=0.                              
C        axisymmetric	: CONFI=1.                                      
C *******************************************************************	
	D =2.*B/((1.+CONFI)*A)                                          
      S =0.2304*EPSONS/SQRT(A)                                          
C       THE ABOVE NUMERICAL VALUE IS FOR TE=300K                        
C                                                                       
C       ----------- DIFFERENTIAL EQUATIONS ----------                   
      DO 32 I=2,MAX(2,M1)                                               
      XI=I-1                                                            
      ACT =E/T(I)                                                       
      IF(ACT.GT.35.) GO TO 35                                           
      ARH=1./EXP(ACT)                                                   
      GO TO 36                                                          
35    CONTINUE                                                          
      ARH=0.                                                            
36    W(I) =D*YF(I)*YO(I)*ARH/T(I)                                      
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
61    CONTINUE                                                          
      TL =T(I)                                                          
      TR =T(I+1)                                                        
      YOL=YO(I)                                                         
      YOR=YO(I+1)                                                       
      YFL=YF(I)                                                         
      YFR=YF(I+1)                                                       
      UL =U(I)                                                          
      UR =U(I+1)                                                        
62    CONTINUE                                                          
      TN(I) = T(I) +DB*(T(I+1)-2.*T(I)+T(I-1))/PR+DA*F(I)*(TR-TL)       
     1    +DT*Q*W(I)                                                    
      YFN(I)=(YF(I)+DB*(YF(I+1)-2.*YF(I)+YF(I-1))/SCF +DA*F(I)*         
     1    (YFR-YFL)) /(1.0+ DT*D*YO(I)*ARH/T(I))                        
      YON(I)=(YO(I)+DB*(YO(I+1)-2.*YO(I)+YO(I-1))/SCO +DA*F(I)*         
     1    (YOR-YOL)) /(1.0-DT*XNO*D*YF(I)*ARH/T(I))                     
                                                                        
C *******************************************************************   
C      CONFIGAURATION PARAMETER IN MOMENTUM EQUATION:                   
C	   2-D			: CONFI=0.                              
C        axisymmetric	: CONFI=1.                                      
C                                                                       
C	 To add the BUOYANT force in the momentum equation:             
C	  IF TSAI=0. IT IS THE PURELY FORCED FLOW CASE!                 
C *******************************************************************	
                                                                        
C      UN(I) =U(I) +DA*F(I)*(UR-UL) +DB* (U(I+1)-2.*U(I)+U(I-1))        
C     1    + 2.*DT*(T(I)-U(I)*U(I))                                     
      UN(I) =U(I) +DA*F(I)*(UR-UL) +DB* (U(I+1)-2.*U(I)+U(I-1))         
     1    + 2./(1.+CONFI)*DT*((T(I)+(T(I)-1.)/(TREF-1.)*TREF            
     2					    *TSAI)/(1.+TSAI)            
     3    -U(I)*U(I))                                                   
                                                                        
32    CONTINUE                                                          
C                                                                       
C                                                                       
C       -------- BOUNDARY CONDITION AT SURFACE ----------               
      UN(1)=0.0                                                         
      FN(1)=FW                                                          
      ACT1=E/T(1)                                                       
      IF(ACT1.GT.35.) GO TO 45                                          
      ARH1 =1./EXP(ACT1)                                                
      GO TO 46                                                          
45    CONTINUE                                                          
      ARH1=0.0                                                          
46    CONTINUE                                                          
      YF0J=FN(1)*2.0*DY*(YF(1)-1.0)*SCF +YF(2)                          
      YFN(1)=(YF(1)+DB*(YF(2)-2.*YF(1)+YF0J)/SCF+.5*DA*(YF(2)           
     1     -YF0J)*FN(1)) /(1.0 +DT*D*YO(1)*ARH1/T(1))                   
      YO0J=FN(1)*2.0*DY*YO(1)*SCO +YO(2)                                
      YON(1)=(YO(1)+DB*(YO(2)-2.*YO(1)+YO0J)/SCO+.5*DA*(YO(2)           
     1     -YO0J)*FN(1))/(1.0-XNO*DT*D*YF(1)*ARH1/T(1))                 
      W(1) = D*YFN(1)*YON(1)*ARH1/T(1)                                  
C                                                                       
C       -------- TO OBTAIN TW ---------                                 
      TG = T(1) -TN(M)                                                  
      TIMAG =T(2)-2.*DY*(-FN(1))*PR*TG-2.*DY*S*(T(1)**4.-1.)            
      TN(1) =(T(2)+TIMAG)/2.                                            
      TN(1) = (TN(1)+T(1))/2.                                           
C                                                                       
C       ------- INTEGRATION OF VELOCITY TO OBTAIN F -------             
      FN(2) =FN(1)+DY*(UN(2)+UN(1))                                     
      DO 51 I=1,MAX(1,M2)                                               
      FN(I+2) =FN(I) +2.*DY*(UN(I+2) +4.*UN(I+1) +UN(I))/3.0            
51    CONTINUE                                                          
C                                                                       
      TMAX =TN(I)                                                       
      DO 555 I=1,MAX(1,M1)                                              
      TMAX = DMAX1(TMAX, TN(I+1))                                       
555   CONTINUE                                                          
C                                                                       
      DO 88 I=1,MAX(1,M)                                                
      T(I)=TN(I)                                                        
      YF(I)=YFN(I)                                                      
      YO(I)=YON(I)                                                      
      U (I)=UN (I)                                                      
      F (I)=FN (I)                                                      
88    CONTINUE                                                          
C                                                                       
C       ------- PRINT OUT -------                                       
      IF (MOD(KOUNT,1500).NE.0) GO TO 31                                
558   FORMAT('       TMAX        TWALL  ')                              
      WRITE(4,70) TMAX, T(1)                                            
70    FORMAT(2(5X,E8.3))                                                
      XJ = KOUNT                                                        
      TAU = XJ*DT                                                       
      WRITE(5,25)                                                       
25    FORMAT( '---------------------------------------------------------
     *----------------------------------------')                        
      WRITE(5,11) TAU,A,F(1)                                            
      WRITE(5,24) D,S,XMDOT                                             
11    FORMAT( ' TAU = ',F10.5,5X,' A = ',F6.2,5X,' FW = ',F7.3)         
24    FORMAT(' D = ',E12.7,5X,' S = ',F10.5,5X,' XMDOT = ',E12.7)       
      WRITE(5,23) DT1,DT2,DY                                            
23    FORMAT(' DT1 = ',F8.4,5X,' DT2 = ',F8.4,5X,' DY = ',F8.4)         
      WRITE(5,561) SCF,SCO,PR                                           
      WRITE(5,20)                                                       
      WRITE(5,60)                                                       
      WRITE (UNIT=5,FMT=21) (T(I),I=1,MAX(1,M))                         
      WRITE(5,20)                                                       
      WRITE(5,65)                                                       
      WRITE (UNIT=5,FMT=21) (YF(I),I=1,MAX(1,M))                        
      WRITE(5,20)                                                       
      WRITE(5,66)                                                       
      WRITE (UNIT=5,FMT=21) (YO(I),I=1,MAX(1,M))                        
      WRITE(5,20)                                                       
      WRITE(5,63)                                                       
      WRITE (UNIT=5,FMT=21) (W(I),I=1,MAX(1,M))                         
      WRITE(5,20)                                                       
      WRITE(5,64)                                                       
      WRITE (UNIT=5,FMT=21) (U(I),I=1,MAX(1,M))                         
      WRITE(5,67)                                                       
      WRITE (UNIT=5,FMT=68) (F(I),I=1,MAX(1,M))                         
      WRITE(5,22)                                                       
60    FORMAT('  TEMPERATURE PROFILE')                                   
65    FORMAT('  FUEL MASS FRACTION PROFILE')                            
66    FORMAT('  OXYGEN MASS FRACTION PROFILE')                          
63    FORMAT('  REACTION RATE DISTRIBUTION')                            
64    FORMAT('  VELOCITY PROFILE')                                      
67    FORMAT('  STREAM FUNCTION')                                       
68    FORMAT(13(1X,F7.3))                                               
20    FORMAT(/)                                                         
21    FORMAT(13(1X,E8.3))                                               
22    FORMAT(///)                                                       
      IF(ABS(TMAX-TMAX1).LT.0.001) GO TO 91                             
      IF(TMAX.LT.2.5) GO TO 98                                          
      TMAX1 = TMAX                                                      
C                                                                       
31    CONTINUE                                                          
      KOUNT =KOUNT +1                                                   
      IF (KOUNT.GT.N) GO TO 99                                          
      GO TO 101                                                         
98    CONTINUE                                                          
      WRITE(4,559)                                                      
      WRITE(4,97)                                                       
      WRITE(5,97)                                                       
97    FORMAT('  TMAX IS LESS THAN 2.5.  COMPUTATION STOPPED.')          
      GO TO 99                                                          
91    CONTINUE                                                          
      WRITE(4,559)                                                      
      WRITE(4,96)                                                       
      WRITE(5,96)                                                       
96    FORMAT( ' TMAX HAS NOT CHANGED SINCE LAST PRINTOUT. STOPPED. ')   
      GO TO 99                                                          
99    CONTINUE                                                          
      WRITE(5,25)                                                       
      WRITE(5,560) A,FW,Q,EPSONS,N                                      
      WRITE(5,25)                                                       
      WRITE(5,11) TAU,A,F(1)                                            
      WRITE(5,24) D,S,XMDOT                                             
      WRITE(5,23) DT1,DT2,DY                                            
      WRITE(5,561) SCF,SCO,PR                                           
601   FORMAT(//5X,'The Simulation for Gas Fuel Diffusion Flame '/,      
     1  'with Mixed Convection for 2-D and Axisymmetric ',              
     2   'Tsuji-type Burner'/,/6X, 'Prof.J.S. Tien'/,/10X,              
     3    'Modified by Bai Han in March, 2002'//)                       
602	FORMAT(' Please type in "0" for 2-D case OR ',                  
     1  'type in "1" for axisymmetric case, ',                          
     2   'Then press the ENTER')                                        
603	FORMAT(' Please type in "n" for purely forced flow OR ',        
     1  'type in "y" for the case with normal gravity: ',               
     2   'Then press the ENTER')                                        
604	FORMAT (//' INPUT Radius of Burner(m), Ambient Velosity(m/s)',  
     1        ', and FUEL FEED RATE (FW,negative):3F6.3')               
      STOP                                                              
      END                                                               
                                                                        
