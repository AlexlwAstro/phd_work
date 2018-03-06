C*****************************************************************************
C  QUESTA ROUTINE CONTIENE: EOS - NEUTR - NKAPPA - NKAPPANODIFF - NKAPPADIFF - 
C                           NSUPERA - SIGMANEW 
C*****************************************************************************
C
C
C===================================================================
C====              ROUTINE FOR THE EOS TABLES                   ====
C===================================================================
C====   version 2.0, 24-02-2015                                 ====
C===================================================================
      SUBROUTINE EOS(MESH,NP,NE,ND,P,T,XX,ROUT,RO,GRAD,CP,PMU,
     #G1,DEL,MUE,NVAR)
C==== variables: mesh number, ***, number of elements, number of ***,
C==== pressure, temperature, elemental abundances, output label,
C==== density, adiabatic
C==== gradient, heat capacity at const. P, ***, ***, ***, mean molecular
C==== mass per electron, number of ***
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
C===================================================================
      CHARACTER*1 TESTO(82),ROUT(8),TFIN(14)
      INTEGER*4 IS,IT,IP(4),NT,NPG,NT0,MESH,NP,NE,ND,INI,IFIN,     
     #IDT(4),IDP(4),ERRORE,NVAR
C Santi - linea originale di sotto
cc      REAL*8 P(1),T(1),XX(NE,NP),RO(NP,3),GRAD(NP,3),CP(NP,3),
C===============================================================
      REAL*8 P(1),T(1),XX(NE,LIM),RO(NP,3),GRAD(NP,3),CP(NP,3),
     #PMU(NP,3),ROGAR,ADGAR,PMUGAR,CPGAR,TEMP,PGAS,TTT,PPP
      REAL*8 G1GAR,DELGAR,MUEGAR,G1(NP,3),DEL(NP,3),MUE(NP,3)
      REAL*8 YTAB(6),DP(4),DELTAT,DELTAP,XTAB(6)    
      REAL*8 TERMO(6,90,230,10), THERMO(6,90,230,10),TT(230),PP(230)
      REAL*8 OUTP(6,4,10),OUTT(6,10),ADGAR1,ADGAR2
      REAL*8 HYD,HE3,HE4,CAR,NIT,OXY,X,Y,XA(4),YA(4),ZINI,DIFF
      LOGICAL OK,IDEBUG   
C===================================================================
      DATA IS/0/,NPG/90/,NT/230/,NQT/6/      
      DATA PASSOLT/0.02D0/,PASSOHT/0.04D0/,PASSOP/0.2D0/DPMAX/17.8D0/
C HT=high temp., LT=low temp., P=pressure, DP=
      DATA YTAB/0.0D0,0.2D0,0.3D0,0.4D0,0.7D0,1.0D0/
      DATA DELTAT/1.D-6/,DELTAP/1.D-6/
c      DATA XTAB/1.0D0,0.8D0,0.7D0,0.6D0,0.3D0,0.0D0/

      data ippo/0/,ippa/0/,ico/0/
            
      IDEBUG=.TRUE.
C=================================================================== 
  101 FORMAT(1X,82A1)
  151 FORMAT('EoS:',82A1)
  102 FORMAT(12X,F8.5,38X,F10.5)  
  103 FORMAT(11X,3F12.5,3D14.5,14A1)  
C===================================================================      
C============== READING EOS TABLES =================================
C===================================================================
      IF(IS.EQ.0) THEN
         OPEN(3,FILE='eos1ext.dat',FORM='FORMATTED',STATUS='OLD') 
        OPEN(20,FILE='provaeos.deb',FORM='FORMATTED')
         REWIND(3)

c     NOTE that only 1 EOS filename is accepted for input - need
c     to select an EOS file and copy into eos1ext.dat

C===================================================================
C============== READING THE HEADERS ================================
C===================================================================
         READ(3,101) TESTO
         WRITE(2,101) TESTO ! to stampe (used in main.f)
         WRITE(*,101) TESTO
            WRITE(10,151) TESTO ! to grafi2p1 (used in main.f)
         DO J=1,17
          READ(3,101) TESTO     
         END DO
C===================================================================         
         DO K=1,10
          READ(3,101) TESTO
          WRITE(*,101) TESTO         
C         
         NT0=1
            IF(K.GE.3.AND.K.LE.6) NT0=151
C            
          DO J=NT0,NT
           READ(3,102) TT(J),PP(J) ! get log(T) and log(Pmin)
                   IF(IDEBUG) WRITE(20,102) TT(J),PP(J)
           DO N=1,NPG
            READ(3,103) (TERMO(I,N,J,K),I=1,NQT),TFIN
                IF(IDEBUG) WRITE(20,103) (TERMO(I,N,J,K),I=1,NQT)   
              END DO
          END DO
         END DO         
C         
         CLOSE(3)
         CLOSE(20)
C===================================================================
C=== RIORGANIZATION OF THE ABUNDANCES IN THE TABLES              ===
C===================================================================
      DO J=1,NT
       DO I=1,NQT
        DO N=1,NPG
        THERMO(I,N,J,1)=TERMO(I,N,J,1) 
        END DO
       END DO
      END DO
C===
      DO L=2,5
       K=L+5
        DO J=1,NT
        DO I=1,NQT
          DO N=1,NPG
          THERMO(I,N,J,L)=TERMO(I,N,J,K) 
          END DO
        END DO
        END DO
       END DO
C===       
      DO J=1,NT
       DO I=1,NQT
        DO N=1,NPG
        THERMO(I,N,J,6)=TERMO(I,N,J,2) 
        END DO
       END DO
      END DO
C===      
      DO L=7,10
       K=L-4
        DO J=151,NT
        DO I=1,NQT
          DO N=1,NPG
          THERMO(I,N,J,L)=TERMO(I,N,J,K) 
          END DO
        END DO
        END DO
       END DO
       
       do i=1,10
         do j=1,4
          do k=1,6            
            OUTP(k,j,i)=0.d0
          enddo
         enddo
       enddo
       
C========================================================================
C=== NOW: Y=0.0, 0.2, 0.3, 0.4, 0.7, 1-Z, PURE HE, PURE C, PURE N, PURE O
C========================================================================
      YTAB(6)=1.0D0-ZINI
      STEPLT=1.0D0/PASSOLT  ! IT IS USED FOR THE TABLE POSITIONING 
      STEPHT=1.0D0/PASSOHT  !          ,,  ,,  ,,
      STEPP=1.0D0/PASSOP    !          ,,  ,,  ,,
      IS=1
      ENDIF
C========================================================================
C====  MAIN LOOP
C========================================================================
      DO KK=1,NP
         IND=KK+MESH-1   !CHIMICA SHIFTATA SE MESH>1
         HYD = XX(1,IND)
         HE3 = XX(2,IND)
         HE4 = XX(3,IND)
         CAR = XX(4,IND)
         NIT = XX(5,IND)
         OXY = XX(6,IND)

C====
         DIFF = 1.D0-(HYD+HE3+HE4+NIT+CAR+OXY)
         OXY  = OXY+DIFF ! = 1 - sum(all other elements)

C====
       IF(HE4.GT.YTAB(6)) THEN
        HE4=YTAB(6)
       ENDIF
       
C====      
       IF(HYD.GE.1.0D-30) THEN
          OK=.TRUE.
          K = 1
          DO WHILE (OK)
           K = K+1
            IF(YTAB(K).GE.HE4) THEN
               L  = K
               LL = L-1
               OK = .FALSE.
             ENDIF
             IF(K.EQ.6) OK=.FALSE.
          END DO
C====
          AY2=YTAB(K-1)
          AY1=YTAB(K)
C==== 
          BETA1=(AY2-HE4)/(AY2-AY1)        ! LL ---> K-1
          BETA2=(HE4-AY1)/(AY2-AY1)        ! L  ---> K
C====
          INI=LL
         IFIN=L
        GO TO 33
       ENDIF
      IF(HE4.GE.1.0D-30) THEN      
       IF(G(6,1).LT.0.D0.AND.G(4,1).LT.150.D0.
     #AND.XXX(4,1).LT.0.2D0) THEN
        INI=6
        IFIN=6
         BETA1=0.D0
         BETA2=1.0D0
         HYD=0.D0
        if(ippo.eq.0) then
         ippo=1
         write(*,*) mesh, he4, ini,ifin, beta1, beta2
        endif
       ELSE
         INI=7
        IFIN=10
        HYD=0.D0
       if(ippa.eq.0) then
       ippa=1
       write(*,*) mesh,he4,car,oxy,nit,ini,ifin
       write(*,*) he4+car+oxy+nit
       write(*,*) g(6,1),g(4,1)
c       pause
       endif
       ENDIF         
      ELSE
         INI=7
        IFIN=10
        HYD=0.D0
        HE4=0.D0
        NIT=0.D0
      ENDIF    
C===================================================================        
 33   CONTINUE  
C===================================================================
C==== LOOP ON T AND P
C===================================================================
         DO KD=1,ND
            TTT = T(KK) ! constant assignment here w.r.t. KD
            PPP = P(KK)
            IF (KD.EQ.2) THEN
               PPP = PPP*(1.0D0+DELTAP)
            ELSE IF (KD.EQ.3) THEN
               TTT = TTT*(1.0D0+DELTAT)
            END IF
            PGAS = PPP-CN3*TTT**4 ! Note that gas pressure used
            IF (PGAS.LE.0.D0) THEN
               ERRORE=3
               CALL STOPERR(PPP,TTT,ICHIM,MESH,ERRORE,ROUT)
            END IF
            TEMP = CONV*DLOG(TTT)
            PGAS = CONV*DLOG(PGAS)

            IF(TEMP.LT.3.D0.OR.TEMP.GT.9.12D0) THEN
               ERRORE=2
               CALL STOPERR(PGAS,TEMP,ICHIM,MESH,ERRORE,ROUT)
            END IF
C===================================================================
C==== POSITIONING IN TEMPERATURE                                ====
C===================================================================
           IF(T(KK).LT.1.D+06) THEN       ! case log(T)<6
              DTEMP = TEMP-3.0D0
              IT = INT(DTEMP*STEPLT)+1 
           ELSE
              DTEMP = TEMP-6.0D0
             IT=INT(DTEMP*STEPHT)+151   ! 151 index of logT=6 
            ENDIF ! significance of 151?
C===================================================================
C==== POSITIONING IN PRESSURE                                   ====
C===================================================================
            IF(IT.LT.151) THEN
             DO J=1,4
             DP(J)=PGAS-PP(IT-2+J)
             END DO
           ELSEIF(IT.EQ.151) THEN
             DP(1)=PGAS-PP(IT-2)
            DP(2)=PGAS-PP(IT)
            DP(3)=PGAS-PP(IT+2)
            DP(4)=PGAS-PP(IT+3)
           ELSEIF(IT.EQ.152) THEN
             DP(1)=PGAS-PP(IT-1)
            DP(2)=PGAS-PP(IT+1)
            DP(3)=PGAS-PP(IT+2)
            DP(4)=PGAS-PP(IT+3)
            ELSE
              DO J=1,4
              DP(J)=PGAS-PP(IT-1+J)
              END DO
            ENDIF        
C====
             DO J=1,4      
              IF(DP(J).LT.0.D0.OR.DP(J).GT.DPMAX) THEN
              ERRORE=2
              CALL STOPERR(PGAS,TEMP,ICHIM,MESH,ERRORE,ROUT)
               ENDIF
              
             IP(J)=INT(DP(J)*STEPP)+1
             END DO

C===================================================================
C==== START INTERPOLATING                                       ==== 
C===================================================================
C==== INTERPOLATION IN PRESSURE                                 ====
C===================================================================
      DO L=INI,IFIN
       DO J=1,6
        OUTT(J,L)=1.0D0
       END DO
      END DO
C====     
      DO IJ=INI,IFIN   

       IF(INI.EQ.8.AND.IJ.EQ.9) GO TO 34

       DO M=1,NVAR
            
       X=PGAS
            
            DO I=1,4
         IF(IT.LT.151) THEN
          IDT(I)=IT-2+I
          ELSEIF(IT.GE.229)THEN
          IDT(I)=226+I
         ELSE         
           IDT(I)=IT-1+I
         ENDIF  


         IF(IT.EQ.151) THEN 
          IF(I.EQ.1) IDT(1)=IT-2
          IF(I.EQ.2) IDT(2)=IT
          IF(I.EQ.3) IDT(3)=IT+2
          IF(I.EQ.4) IDT(4)=IT+3
         ELSEIF(IT.EQ.152) THEN
          IF(I.EQ.1) IDT(1)=IT-1
          IF(I.EQ.2) IDT(2)=IT+1
          IF(I.EQ.3) IDT(3)=IT+2
          IF(I.EQ.4) IDT(4)=IT+3
         ENDIF  
       
          DO N=1,4   
           
           IF(IP(I).EQ.1) THEN
            IDP(N)=N
           ELSEIF(IP(I).GE.76)THEN
            IDP(N)=73+N
           ELSE
            IDP(N)=IP(I)-2+N
           ENDIF
          
           XA(N)=PP(IDT(I))+PASSOP*DFLOAT(IDP(N)-1) !FLOAT(IDP(N)-1)
           YA(N)=THERMO(M,IDP(N),IDT(I),IJ) 
           END DO
            CALL POLINT(XA,YA,4,X,Y,DY)
            OUTP(M,I,IJ)=Y
        END DO
              
       END DO
 34   CONTINUE             
      END DO  
C===================================================================
C==== INTERPOLATION IN TEMPERATURE                              ====
C===================================================================
      DO IJ=INI,IFIN

       IF(INI.EQ.8.AND.IJ.EQ.9) GO TO 35

       DO M=1,NVAR
            
       X=TEMP
            
           DO I=1,4
        XA(I)=TT(IDT(I))
        YA(I)=OUTP(M,I,IJ)
       END DO    
         CALL POLINT(XA,YA,4,X,Y,DY)
         OUTT(M,IJ)=Y
       END DO
 35   CONTINUE              
      END DO
      
      DO JJ=INI,IFIN
        OUTT(1,JJ)=DEXP(COND*OUTT(1,JJ))
      IF(NVAR.EQ.1) GO TO 36       
        OUTT(3,JJ)=DEXP(COND*OUTT(3,JJ))
 36   CONTINUE       
      END DO
C===================================================================
      IF(INI.LE.6) THEN
C===================================================================
C==== ADDITIVE VOLUME LAW                                       ====
C===================================================================       

       ROGAR= 1.0D0/((BETA1/OUTT(1,IFIN))+(BETA2/OUTT(1,INI)))
       IF(NVAR.EQ.1) GO TO 37
      
       CPGAR= BETA1*OUTT(3,IFIN)+BETA2*OUTT(3,INI)

       ADGAR1=OUTT(3,IFIN)*BETA1
       ADGAR1=ADGAR1*OUTT(2,IFIN)
       
       ADGAR2=OUTT(3,INI)*BETA2
       ADGAR2=ADGAR2*OUTT(2,INI)

       ADGAR=ADGAR1+ADGAR2
      
       ADGAR =ADGAR/CPGAR
       IF(NVAR.GT.1.AND.NVAR.LE.3) GO TO 38
       
       G1GAR = BETA1*OUTT(4,IFIN)+BETA2*OUTT(4,INI)
       DELGAR = BETA1*OUTT(5,IFIN)+BETA2*OUTT(5,INI)
       MUEGAR = BETA1*OUTT(6,IFIN)+BETA2*OUTT(6,INI)
              
C       IF(NVAR.GT.3.AND.NVAR.LE.6) GO TO 39
           
      ELSE
      
       ROGAR = 1.0D0/(HE4/OUTT(1,7)+CAR/OUTT(1,8)+
     #NIT/OUTT(1,9)+OXY/OUTT(1,10))
     
       IF(NVAR.EQ.1) GO TO 37
       
       CPGAR= HE4*OUTT(3,7)+CAR*OUTT(3,8)+       
     #NIT*OUTT(3,9)+OXY*OUTT(3,10)
     
       ADGAR= (OUTT(3,7)*HE4*OUTT(2,7)+
     #OUTT(3,8)*CAR*OUTT(2,8)+OUTT(3,9)*NIT*OUTT(2,9)+
     #OUTT(3,10)*OXY*OUTT(2,10))/CPGAR

       IF(NVAR.GT.1.AND.NVAR.LE.3) GO TO 38

       G1GAR= HE4*OUTT(4,7)+CAR*OUTT(4,8)+       
     #NIT*OUTT(4,9)+OXY*OUTT(4,10)
       DELGAR= HE4*OUTT(5,7)+CAR*OUTT(5,8)+       
     #NIT*OUTT(5,9)+OXY*OUTT(5,10)
       MUEGAR= HE4*OUTT(6,7)+CAR*OUTT(6,8)+       
     #NIT*OUTT(6,9)+OXY*OUTT(6,10)

            
       ENDIF
       
       G1(KK,KD)   = G1GAR
       DEL(KK,KD)  = DELGAR  
       MUE(KK,KD)  = MUEGAR
       
 38   CONTINUE
       CP(KK,KD)   = CPGAR
       GRAD(KK,KD) = ADGAR
       PMU(KK,KD)  = ((ROGAR*TTT)/PPP)*(CPGAR*ADGAR)   
C    (rho*T/P)*(C_P)*(dlnT/dlnP)   
 37    CONTINUE
       
       RO(KK,KD) = ROGAR      

C===================================================================              
C===================================================================              
         END DO

         OXY = OXY - DIFF

      END DO
      
      RETURN
      END
C===================================================================              
C===================================================================              
C===================================================================              

      SUBROUTINE NEUTR(RO,T,XX,QT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
C************* PRODUZIONE DI NEUTRINI **********************************
      DIMENSION XX(MAXNE)
      DIMENSION A0(3),A1(3),A2(3),B1(3),B2(3),B3(3),C1(3),F1(3)
      DATA A0/6.002D+19,4.886D+10,2.320D-07/
      DATA A1/2.084D+20,7.580D+10,8.449D-08/
      DATA A2/1.872D+21,6.023D+10,1.787D-08/
      DATA B1/9.383D-01,6.290D-03,2.581D-02/
      DATA B2/-4.141D-1,7.483D-03,1.734D-02/
      DATA B3/5.829D-02,3.061D-04,6.990D-04/
      DATA C1/5.5924D0,1.5654D0,0.56457D0/
C******  NN NUMERO DI NEUTRINI OLTRE QUELLI ELETTRONICI ****************
C *****  NN = number of neutrino types other than electron neutrinos
      DATA NN/2/,IALL/0/
      QT=0.D0
      IF(T.LT.1.D+08.AND.XX(1).GT.1.D-7) RETURN
      CVPA=1.122356D0+NN*0.254356D0
      CVMA=0.622356D0-NN*0.245644D0
      SUM=0.D0
      DO IC=1,MAXNE-1
         SUM=SUM+XX(IC)*ZET(IC)/ATW(IC)
      END DO
      VME=1.D0/SUM
C***********************************************************************
C PRODUZIONE NEUTRINI DA:    PLASMA (Q1) - PHOTO (Q2) - COPPIE (Q3)
C        MUNAKATA,H., KOHYAMA,Y., ITOH,N. AP.J. 296,197 - 1985
C                + ERRATA DEL 1986
C***********************************************************************
      VLA=T/5.93097D+09
      VLA2=VLA*VLA
      VLA3=VLA*VLA2
      VLA4=VLA2*VLA2
      VLOLA=DLOG(VLA)
      ROME=RO/VME
      CSI=((ROME/1.D+09)**(1.D0/3.D0))/VLA
      DO K=1,3
         F1(K)=0.D0
         CHK=-C1(K)*CSI
         IF(DABS(CHK).LE.40.D0) THEN
            F1(K)=(A0(K)+A1(K)*CSI+A2(K)*CSI*CSI)*DEXP(CHK)
            F1(K)=F1(K)/(CSI*CSI*CSI+B1(K)/VLA+B2(K)/VLA2+B3(K)/VLA3)
         ENDIF
      END DO
C***************************  COPPIE  **********************************
      Q3=0.D0
      CHK2=DABS(2.D0/VLA)
      IF(CHK2.LE.40.D0) THEN
         G1=1.D0-13.04D0*VLA2+133.5D0*VLA4+1534.D0*VLA4*VLA2+918.6D0
     &      *VLA4*VLA4
         VLA05=DEXP(0.5D0*VLOLA)
         QPA=1.D0/(10.748D0*VLA2+.3967D0*VLA05+1.005D0)*(1.D0+ROME/
     &       (7.692D7*VLA3+9.715D6*VLA05))**(-.3)
         Q3=.5D0*CVPA*(1.D0+CVMA/CVPA*QPA)*G1*DEXP(-CHK2)*F1(1)
      ENDIF
C*************************  FOTONEUTRINI  ******************************
      Q2=0.D0
      VCSI=DLOG(CSI)
      CSI448=DEXP(4.48D0*VCSI)
      CSI348=DEXP(3.48D0*VCSI)
      CHK3=(.556D0*CSI448/(150.D0+CSI348))
      IF(DABS(CHK3).LE.40.D0) THEN
         VLA355=DEXP(3.555D0*VLOLA)
         QPHO=.666D0*(1.D0+2.045D0*VLA)**(-2.066)/(1.D0+ROME/(1.875D8
     &        *VLA+1.653D8*VLA2+8.499D8*VLA3-1.604D8*VLA4))
         Q2=.5D0*CVPA*(1.D0-CVMA/CVPA*QPHO)*VLA4*VLA4*1.D9*CSI**3*F1(2)
     #   *.893D0*(1.D0+143.8D0*VLA355)**.3516*DEXP(CHK3)
      ENDIF
C***************************  PLASMA  **********************************
C      Q1=(0.872356D0+NN*4.356D-03)*ROME**3*F1(3)
C***********************************************************************
C   FORMULA PER PLASMA NEUTRINI DA HAFT, RAFFELT & WEISS 94 ApJ 425,222
C***********************************************************************
C***************************  PLASMA  **********************************
      Q1=0.D0
      XGAM2=(1.1095D11*ROME)/(T**2.D0*(1.D0+
     #        (1.019D-6*ROME)**(2.D0/3.D0))**0.5D0)
      XGAM=XGAM2**0.5D0
      FT=2.4D0+0.6D0*XGAM**0.5D0+0.51D0*XGAM+1.25D0*XGAM**(3.D0/2.D0)
      FL=(8.6D0*XGAM2+1.35D0*XGAM**(7.D0/2.D0))/
     #        (225.0D0-17.0D0*XGAM+XGAM2)
      XICS=(1.D0/6.D0)*(17.5D0 + DLOG10(2.D0*ROME) - 3.D0*DLOG10(T))
      XY=(1.D0/6.D0)*(-24.5D0+DLOG10(2.D0*ROME)+3.D0*DLOG10(T))
      DXX=DABS(XICS)
      VAR=XY-1.6D0+1.25D0*XICS
      IF(DXX.GT.0.7D0.OR.XY.LT.0.D0) THEN
      FXY=1.D0
      ELSE
      ZERO=0.D0   
      FXY= 1.05D0+(0.39D0-1.25D0*XICS-0.35D0*DSIN(4.5D0*XICS)-
     #        0.3D0*DEXP(-(4.5D0*XICS+0.9D0)**2.D0))*
     #        DEXP(-(DMIN1(ZERO,VAR)/(0.57D0-0.25D0*XICS))**2.D0)
      ENDIF
      QAPPR=3.00D+21*(VLA**9.D0)*(XGAM**6.D0)*DEXP(-XGAM)*(FT+FL)*FXY
      Q1=0.9325D0*QAPPR
C*************************** RICOMBINAZIONE ****************************
      ZM=0.D0
      AM=0.D0
      DO IC=1,MAXNE-1
         ZM=ZM+XX(IC)*ZET(IC)
         AM=AM+XX(IC)*ATW(IC)
      END DO
      S1=1.85D-04*(ZM**6)*RO*RO*VLA2/(AM*AM)
      S2=-1.57D+05*ZM*ZM/T-2.428D-05*(RO**(.666666))/VLA
      Q4=S1*DEXP(S2)
C******* NEUTRINI DI BREMSSTRAHLUNG - DICUS ET AL. (1976) AP.J. 210,481
      T8=T/1.D+08
      EFER= (1.018D-04 * (ROME**(2.D0/3.D0)) + 1.D0 )**.5
      BY=1.D0/DSQRT(1.D0-1.D0/(EFER*EFER))
      BA2=BY*BY
      VLB=DLOG(ABS((BY+1.D0)/(BY-1.D0)))
      D2=BY/215.D0
      VLD2=DLOG((2.D0+D2)/D2)
      BB1=-2.D0/3.D0+BA2+.5D0*BY*(1.D0-BA2)*VLB
      BB2=-4.D0+2.D0*(1.D0+D2)*VLD2
      BB3=VLD2-2.D0/(2.D0+D2)
      FB=0.14D0*(BB1*BB2-BB3*(BA2-1.D0)*(2.D0/3.D0+.5D0*BA2-BY/4.D0*
     &   (BA2+1.D0)*VLB))
      GB=0.14D0*BB3*(BA2-1.D0)*(2.D0/3.D0-5.D0/2.D0*BA2-BY/4.D0*
     &   (3.D0-5.D0*BA2)*VLB)
      Q5=2.D0*ZM*ZM/AM*T8**6*(.5D0*CVPA*FB-.5D0*CVMA*GB)*RO
      QT=-(Q1+Q2+Q3+Q4+Q5)/RO
      IF(IALL.EQ.0)RETURN
      WRITE(2,100)RO,T,Q1,Q2,Q3,Q4
      WRITE(2,100)EFER,BY,BA2,VLB,D2,VLD2,BB1,BB2,BB3,FB,GB
  100 FORMAT(1P,11E10.3)
      RETURN
      END
C ===== NKAPPA JUST CHOOSES BETWEEN NKAPPADIFF AND NKAPPANODIFF
      SUBROUTINE NKAPPA(MESH,NP,NE,ND,RHO,T,XX,CAP)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INTEGER*4    NP,NE,ND,MESH
      REAL*8       RHO(NP,3),T(1),XX(NE,LIM),CAP(NP,3)
      COMMON/OPACITY/INDIFFINI

      IF(INDIFF.NE.0) THEN
        CALL NKAPPADIFF(MESH,NP,NE,ND,RHO,T,XX,CAP)
c        CALL NKAPPANODIFF(MESH,NP,NE,ND,RHO,T,XX,CAP)
      ELSE
        CALL NKAPPANODIFF(MESH,NP,NE,ND,RHO,T,XX,CAP)
      ENDIF
      RETURN
      END
C      
      SUBROUTINE NKAPPANODIFF(MESH,NP,NE,ND,RHO,T,XX,CAP)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'consts.2p1'
      INCLUDE 'maincom.2p1'
C ***********  DEFINIZIONE DELLE VARIABILI  *****************************
      CHARACTER*1  TESTOL(53),TESTOH(53),TESTO1(200),COMPO(33)
      INTEGER*4    NP,NE,ND,J,K,I,L,JY(3),VF,ICHIM,KK,KT,IK,
     ^             JK,JY1V,JY2V,IIV,JYF,M,JYY,N,JY1,JY2,
     ^             ICHI,ITEM,INDI,MESH,IND,JYYV(3),INDR0
      REAL*8       RHO(NP,3),T(1),XX(NE,LIM),CAP(NP,3)
      
      REAL*8       TK(87),OP(31,87,13),R,TEMP,T63,CATD(3),
     ^             HYD,HE3,HE4,CAR,OXY,RHP(4),TT(4),TTV,DELTAT,
     ^             A(4),B(4),TDEL,AX(3),AY(3),T63DEL,
     ^             AAA(12),BBB(12),CCC(12),DDD(12),CAPR(4),CAPJ(3),
     ^             CAPS(3),COEFF(4),YTAB(10),SUMMA,CAPDP(3),
     ^             CATDT(3),KPDP(4),KTDT(4),RP,RT,HE
      real*8       atest(4)
      LOGICAL      ESCI,IZ,FAST,HBURN,DEBUG
C ************  COMMON RELATIVO ALLA CUB  ******************************
      COMMON/CUBICA/A,B,COEFF
      COMMON/OPACITY/INDIFFINI
C **********************************************************************
      DATA IZ/.TRUE./
      DATA DEBUG/.FALSE./
      DATA YTAB/1.D0,0.9D0,0.8D0,0.65D0,0.5D0,0.3D0,0.2D0,0.1D0,
     #0.05D0,0.0D0/
C
      DATA JYYV/0.d0,0.d0,0.d0/
      DATA JY1V/0/,JY1/0/,JY2V/0/,JY2/0/,IIV/0/,TTV/0./
      DATA DELTAT/1.D-6/
C **********************************************************************
  100 FORMAT(1X,53A1)
  151 FORMAT('Opacity: ',19A1,' ',53A1,10A1)
  200 FORMAT(200A1)
  299 FORMAT(36A1)
  300 FORMAT(3X,'READING TABLE LOW T FOR MIX:',2X,36A1)
  301 FORMAT(F5.3,29F8.3)
  400 FORMAT(3X,'READING TABLE HIGH T FOR MIX:',2X,36A1)
  401 FORMAT(F5.3,31F8.3)
C ***************  LETTURA TABELLE DI OPACITA' *************************
      IF(IZ) THEN
c         WRITE(*,'(d9.2)')zini
        DO i=1,9
          YTAB(I)=YTAB(I)-ZINI
        ENDDO
         OPEN(3,FILE='opacity.dat',STATUS='OLD')
        IF(DEBUG)open(unit=18,file='lello')
C   **********************  T <= 4.10  ************************
         READ(3,100) TESTOL
        WRITE(2,100) TESTOL
         WRITE(*,100) TESTOL
        IF(DEBUG)WRITE(18,100) TESTOL
         READ(3,200)  TESTO1
        IF(DEBUG)WRITE(18,200)  TESTO1
         DO J=1,10
            READ(3,299)  COMPO
            WRITE(*,300) COMPO
           IF(DEBUG)WRITE(18,300) COMPO
            DO K=1,23
               READ(3,301) TK(K),(OP(L,K,J),L=3,31)
              IF(DEBUG)WRITE(18,301) TK(K),(OP(L,K,J),L=3,31)
            END DO
         END DO
C   **********************  T >= 4.15  ************************
         READ(3,100) TESTOH
         WRITE(2,100) TESTOH
         WRITE(*,100) TESTOH
        IF(DEBUG)WRITE(18,100) TESTOH
         READ(3,200)  TESTO1
         READ(3,200)  TESTO1
        IF(DEBUG)WRITE(18,200)  TESTO1
         DO J=1,10
            READ(3,299)  COMPO
            WRITE(*,400) COMPO
           IF(DEBUG)WRITE(18,400) COMPO
            DO K=24,87    
               READ(3,401) TK(K),(OP(L,K,J),L=1,31)
              IF(DEBUG)WRITE(18,401) TK(K),(OP(L,K,J),L=1,31)
            END DO
         END DO
        WRITE(10,151)
     #(TESTOH(i),i=1,19),TESTOL,(COMPO(I),I=23,32)
C   **** LEGGE TABELLE DI PURO CARBONIO E PURO OSSIGENO ****
         DO I=12,13
            READ(3,299)  COMPO
            WRITE(*,400) COMPO
           IF(DEBUG)WRITE(18,400) COMPO
            DO K=24,87 
               READ(3,401) TEMP,(OP(L,K,I),L=1,31)
              IF(DEBUG)WRITE(18,401) TEMP,(OP(L,K,I),L=1,31)
            END DO
         END DO
C   ******************* CREA TABELLA DI PURO ELIO + Z ******
         DO J=1,87  
            DO K=1,31
               OP(K,J,11)=OP(K,J,1)
            END DO
           IF(DEBUG)write(18,'(31f8.3)')(OP(K,J,11),k=1,31)
         END DO
         IZ=.FALSE.
         CLOSE(3)
      ENDIF
      IF(DEBUG)close(18)
C************************* CALCOLO OPACITA'*****************************
      DO VF=1,NP
         IND = VF+MESH-1  !CHIMICA SHIFTATA SE MESH>1
         HYD = XX(1,IND)
         HE3 = XX(2,IND)
         HE4 = XX(3,IND)
         CAR = XX(4,IND)
         OXY = XX(6,IND)
 
         HE=HE3+HE4

C     ******* SE STA BRUCIANDO H SCEGLIE LE DUE TABELLE IN Y ***

         IF(HYD.GT.1.D-30) THEN
   
            DO J=1,9
             IF(HE.LE.YTAB(J).AND.HE.GE.YTAB(J+1)) THEN
              IF(J.EQ.9) THEN
               JY(1)=8
                JY(2)=9
                JY(3)=10
                JYF=3
             ELSE
                JY(1)=J
                JY(2)=J+1
                JY(3)=J+2
                JYF=3
              ENDIF
             HBURN=.TRUE.
              GO TO 1
             ENDIF
           ENDDO
         ELSE
            JY(1)   = 11
            JY(2)   = 12
            JY(3)   = 13
            JYF   =  3
            HBURN = .FALSE.
         ENDIF
  1      CONTINUE

        T63   = (T(VF)/1.0D6)**3.D0 
        T63DEL = ((T(VF)*(1.D0+DELTAT))/1.0D6)**3.D0 
         R     = CONV*DLOG(RHO(VF,1)/T63)
        TEMP  = CONV*DLOG(T(VF))
        IF(TEMP.LT.3.D0.OR.TEMP.GT.9.065D0) THEN
          WRITE(*,*)'FUORI RANGE TEMPERATURA IN OPACITA`'
          STOP
        ENDIF
        IF(TEMP.LT.4.15D0) THEN
          IF(R.LT.-7.D0.OR.R.GT.7.D0) THEN
            WRITE(*,*)'FUORI RANGE R --> LogT<4.15'
            STOP
          ENDIF
        ELSE
          IF(R.LT.-8.D0.OR.R.GT.7.D0) THEN
            WRITE(*,*)'FUORI RANGE R --> LogT>4.15'
            STOP
          ENDIF
        ENDIF
         IF(ND.GT.1) THEN
             RT=CONV*DLOG(RHO(VF,3)/T63DEL) 
             RP=CONV*DLOG(RHO(VF,2)/T63)
         ENDIF
C ******************* CERCA LE QUATTRO TEMPERATURE *********************
c ******************* = SEARCH FOR THE 4 TEMPERATURES
         IF(TEMP.LT.6.0D0) THEN
            K=1
            ESCI=.TRUE.
            DO WHILE (ESCI)
               K=K+1
              IF(K.EQ.61) ESCI=.FALSE.
              IF(K.LT.61.AND.TK(K).GT.TEMP) ESCI=.FALSE.
            END DO
            KT=K-2
            IF(KT.LT.1) KT=1
            KK=KT
            DO K=1,4
               TT(K)=TK(KK)
               KK=KK+1
            END DO
         ELSEIF(TEMP.GE.6.0D0.AND.TEMP.LT.8.1D0) THEN
            JK=60+INT((TEMP-6.0D0)/0.1D0)+1
            KT=JK-1
            IF(JK.GT.80) KT=79
            KK=KT
            DO IK=1,4
               TT(IK)=TK(KK)
               KK=KK+1
            END DO
         ELSEIF(TEMP.GE.8.1D0) THEN
            JK=82+INT((TEMP-8.1D0)/0.2D0)+1
           KT=JK-1
           IF(JK.GT.85) KT=84
            KK=KT
            DO IK=1,4
               TT(IK)=TK(KK)
               KK=KK+1
            END DO
        ENDIF
C ******************* CERCA LE QUATTRO DENSITA' ************************
c ******************* = SEARCH FOR THE 4 DENSITIES
        INDR0=INT((R+8.D0)/0.5D0)+1
        IF(TEMP.LE.4.15D0) THEN
          IF(INDR0.LT.4) INDR0=4
        ELSE
          IF(INDR0.LT.2) INDR0=2
        ENDIF
        IF(INDR0.GT.29) INDR0=29
C **********************************************************************
C * SE GRIGLIATO NON E' CAMBIATO DALLA CHIAMATA PRECEDENTE METTE
C * FAST = TRUE
C * = IF GRAVITY IS NOT CHANGED FROM THE PREVIOUSLY CALLED SET
C * FAST = TRUE
         IF(JY(1).EQ.JYYV(1).AND.JY(2).EQ.JYYV(2).AND.
     #       JY(3).EQ.JYYV(3).AND.TT(1).EQ.TTV.AND.INDR0.EQ.
     #   IIV) THEN
            FAST=.TRUE.
         ELSE
            FAST=.FALSE.
         ENDIF
C **************** INTERPOLA CUBICAMENTE SU R *************************
C **************** = INTERPOLATE CUBICALLY ON R
         DO JYY=1,JYF
         
            
            IF(FAST.EQV..FALSE.) THEN
               DO M=1,4
                  DO N=1,4
                   A(N)=OP((INDR0-2+N),(KT-1+M),JY(JYY))
                     B(N)=-8.d0+(INDR0-3+N)*0.5D0
                  END DO
                 CALL NCUB
                 CAPR(M)=COEFF(1)*R*R*R+COEFF(2)*R*R+
     ^                   COEFF(3)*R+COEFF(4)
                 IF(ND.GT.1) THEN
                    KPDP(M)=COEFF(1)*RP*RP*RP+COEFF(2)*RP*RP+
     ^                      COEFF(3)*RP+COEFF(4)
                    KTDT(M)=COEFF(1)*RT*RT*RT+COEFF(2)*RT*RT+
     ^                      COEFF(3)*RT+COEFF(4)
                 ENDIF
                 INDI = (JYY-1)*4 + M
                 AAA(INDI)=COEFF(1)
                 BBB(INDI)=COEFF(2)
                 CCC(INDI)=COEFF(3)
                 DDD(INDI)=COEFF(4)
               END DO
            ELSE
               DO M=1,4
                  INDI = (JYY-1)*4 + M
                  CAPR(M)=AAA(INDI)*R*R*R+BBB(INDI)*R*R+
     ^            CCC(INDI)*R+DDD(INDI)
                  IF(ND.GT.1) THEN
                     KPDP(M)=AAA(INDI)*RP*RP*RP+BBB(INDI)*RP*RP+
     ^                       CCC(INDI)*RP+DDD(INDI)
                     KTDT(M)=AAA(INDI)*RT*RT*RT+BBB(INDI)*RT*RT+
     ^                       CCC(INDI)*RT+DDD(INDI)
                  ENDIF
               END DO
            ENDIF
C********  INTERPOLA CUBICAMENTE SU T E CALCOLA DK/DR E DK/DT **********
C ******* = INTERPOLATE CUBICALLY ON T & CALCULATE DK/DR & DK/DT
            DO J=1,ND
               DO K=1,4
                  IF(J.EQ.1) THEN
                     A(K)=CAPR(K)
                  ELSE IF(J.EQ.2) THEN
                     A(K)=KPDP(K)
                  ELSE IF(J.EQ.3) THEN
                     A(K)=KTDT(K)
                  ENDIF
                  B(K)=TT(K)
               END DO
               CALL NCUB
               CAPJ(JYY)=COEFF(1)*TEMP*TEMP*TEMP+COEFF(2)*TEMP*TEMP+
     ^                   COEFF(3)*TEMP+COEFF(4)
               IF(J.EQ.1) THEN
                  CAPS(JYY)=DEXP(COND*(CAPJ(JYY)))
               ELSE IF(J.EQ.2) THEN
                      CAPDP(JYY)=DEXP(COND*(CAPJ(JYY)))
               ELSE IF(J.EQ.3) THEN
                      TDEL=CONV*DLOG(T(VF)*(1.D0+DELTAT))
                      CATDT(JYY)=COEFF(1)*TDEL*TDEL*TDEL+COEFF(2)*
     ^                           TDEL*TDEL+COEFF(3)*TDEL+COEFF(4)
                      CATDT(JYY)=DEXP(COND*(CATDT(JYY)))
               ENDIF
            END DO
         END DO
C ********************** INTERPOLAZIONE SU CHIMICA *********************
c ********************** = INTERPOLATION ON CHEMISTRY
         IF (HBURN) THEN
              DO K=1,3
               AY(K)=CAPS(K)
               AX(K)=YTAB(JY(K))
              END DO
           CAP(VF,1)=AY(1)+(HE-AX(1)) * ((AY(2)-AY(1))/(AX(2)-AX(1))+
     #    (HE-AX(2))/(AX(3)-AX(2)) * ((AY(3)-AY(1))/(AX(3)-AX(1)) - 
     #    (AY(2)-AY(1))/(AX(2)-AX(1))) )
         


            IF(ND.GT.1) THEN

              DO K=1,3
               AY(K)=CATDT(K)
               AX(K)=YTAB(JY(K))
              END DO
          CAP(VF,3)=AY(1)+(HE-AX(1)) * ((AY(2)-AY(1))/(AX(2)-AX(1)) +
     #    (HE-AX(2))/(AX(3)-AX(2)) * ((AY(3)-AY(1))/(AX(3)-AX(1)) - 
     #    (AY(2)-AY(1))/(AX(2)-AX(1))) )

              DO K=1,3
               AY(K)=CAPDP(K)
               AX(K)=YTAB(JY(K))
              END DO
          CAP(VF,2)=AY(1)+(HE-AX(1)) * ((AY(2)-AY(1))/(AX(2)-AX(1)) +
     #    (HE-AX(2))/(AX(3)-AX(2)) * ((AY(3)-AY(1))/(AX(3)-AX(1)) - 
     #    (AY(2)-AY(1))/(AX(2)-AX(1))) )

            ENDIF
C ****************** QUADRATIC INTERPOLATION ON Y **********************
         ELSE
C ************ SE H=0 COMBINAZIONE LINEARE SU Y - C - O ****************
C **** IF NO H BURNING, LINEAR COMBINATION OF Y, C, O
            SUMMA   = (HE4+CAR+OXY)
            CAP(VF,1) = (CAPS(1)*HE4+CAPS(2)*CAR+CAPS(3)*OXY)/SUMMA
            IF(ND.GT.1) THEN
              CAP(VF,3) = (CATDT(1)*HE4+CATDT(2)*CAR+CATDT(3)*OXY)/SUMMA
              CAP(VF,2) = (CAPDP(1)*HE4+CAPDP(2)*CAR+CAPDP(3)*OXY)/SUMMA
            END IF
         ENDIF
        DO I=1,3
           JYYV(I)=JY(I)
        ENDDO
         TTV =TT(1)
         IIV =INDR0
      END DO
      RETURN
      END
      SUBROUTINE NKAPPADIFF(MESH,NP,NE,ND,RHO,T,XX,CAP)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(MX=10)
      INCLUDE 'consts.2p1'
      INCLUDE 'maincom.2p1'
C ***********  DEFINIZIONE DELLE VARIABILI  *****************************
      CHARACTER*1  TESTOL(53),TESTOH(53)
      CHARACTER*1  TESTO(72),TESTO1(200),COMPO(33),testata(20)
      INTEGER*4    NP,NE,ND,J,K,I,L,JY(3),VF,ICHIM,KK,KT,IK,
     ^             JK,JY1V,JY2V,IIV,JYF,JOTA,M,JYY,N,JY1,JY2,
     ^             ICHI,ITEM,INDI,MESH,IND,JYYV(3),INDR0,
     #             jj,lij,ki,jfinale
      INTEGER*4    MMM,IOP,NUM,IZi
      REAL*8       RHO(NP,3),T(1),XX(NE,LIM),CAP(NP,3)
      
      REAL*8       TK(87),OP(31,87,13),R,TEMP,T63,CATD(3),
     ^             HYD,HE3,HE4,CAR,OXY,RHP(4),TT(4),TTV,DELTAT,
     ^             A(4),B(4),TDEL,AX(3),AY(3),T63DEL,
     ^             AAA(12),BBB(12),CCC(12),DDD(12),CAPR(4),CAPJ(3),
     ^             CAPS(3),COEFF(4),YTAB(10),SUMMA,CAPDP(3),
     ^             CATDT(3),KPDP(4),KTDT(4),RP,RT,HE,coef,feo,
     #             OPU(31,87,13,10),ZETAVAL(10),ZMA,zmarn,zillav(10),
     #             absolferro
      REAL*4       X,Y,DERIV,ZIN,FVALUE,
     #             FDERIV
      LOGICAL      ESCI,IZ,FAST,HBURN,DEBUG
C ************  COMMON RELATIVO ALLA CUB  ******************************
      COMMON/CUBICA/A,B,COEFF
C **********************************************************************
C ************  COMMON E DATA RELATIVI ALLA SPLINk    ******************
      COMMON/SP/X(MX),Y(MX),DERIV(MX,2),ZIN(MX),FVALUE(MX),
     #             FDERIV(MX,2)
      DATA MMM,MMC,IOP,NUM/1,1,-1,10/
      DATA IZ/.TRUE./
      DATA DEBUG/.TRUE./
      DATA YTAB/1.D0,0.9D0,0.8D0,0.65D0,0.5D0,0.3D0,0.2D0,0.1D0,
     #0.05D0,0.0D0/
C
      DATA JYYV/0.d0,0.d0,0.d0/
      DATA JY1V/0/,JY1/0/,JY2V/0/,JY2/0/,IIV/0/,TTV/0./
      DATA DELTAT/1.D-6/
C **********************************************************************
  100 FORMAT(1X,72A1)
  151 FORMAT('Opacity: ',19A1,' ',53A1,10A1)
  200 FORMAT(200A1)
  299 FORMAT(36A1)
  300 FORMAT(3X,'READING TABLE LOW T FOR MIX:',2X,36A1)
  301 FORMAT(F5.3,29F8.3)
  333 FORMAT('Z OUT OF TABLE: Z=',D9.2,' Z(1)=',D9.2,' Z(10)=',D9.2,
     #' mesh=',I7)
  400 FORMAT(3X,'READING TABLE HIGH T FOR MIX:',2X,36A1)
  401 FORMAT(F5.3,31F8.3)
  402 FORMAT(14X,D8.1)
  403 FORMAT(20A1,14X,f10.7)
  404 FORMAT('LETTA TABELLA PER METALLICITA` Z=',D8.1)
C ***************  LETTURA TABELLE DI OPACITA' *************************
      IF(IZ) THEN
c       WRITE(*,'(d9.2)')zini
       DO i=1,9
        YTAB(I)=YTAB(I)-ZINI
       ENDDO
       OPEN(3,FILE='opacyall.dat',STATUS='OLD')
       IF(DEBUG)open(unit=18,file='lello')
C   **********************  T <= 4.10  ************************
        read(3,403)testata
       write(*,403)testata
       DO N=1,10
        READ(3,402)ZETAVAL(N)
        WRITE(*,404)ZETAVAL(N)
        READ(3,100) TESTOL
        IF(N.EQ.10)WRITE(2,100) TESTOL
         IF(N.EQ.10)WRITE(*,100) TESTOL
        IF(DEBUG)WRITE(18,100) TESTOL
         READ(3,200)  TESTO1
        IF(DEBUG)WRITE(18,200)  TESTO1
         DO J=1,10
            READ(3,299)  COMPO
            IF(N.EQ.10)WRITE(*,300) COMPO
           IF(DEBUG)WRITE(18,300) COMPO
            DO K=1,23
               READ(3,301) TK(K),(OPU(L,K,J,N),L=3,31)
              IF(DEBUG)WRITE(18,301) TK(K),(OPU(L,K,J,N),L=3,31)
            END DO
         END DO
C   **********************  T >= 4.15  ************************
         READ(3,100) TESTOH
         IF(N.EQ.10)WRITE(2,100) TESTOH
         IF(N.EQ.10)WRITE(*,100) TESTOH
        IF(DEBUG)WRITE(18,100) TESTOH
         READ(3,200)  TESTO1
         READ(3,200)  TESTO1
        IF(DEBUG)WRITE(18,200)  TESTO1
         DO J=1,10
            READ(3,299)  COMPO
            IF(N.EQ.10)WRITE(*,400) COMPO
           IF(DEBUG)WRITE(18,400) COMPO
            DO K=24,87    
               READ(3,401) TK(K),(OPU(L,K,J,N),L=1,31)
              IF(DEBUG)WRITE(18,401) TK(K),(OPU(L,K,J,N),L=1,31)
            END DO
         END DO
        IF(N.EQ.10) then
C   **** LEGGE TABELLE DI PURO CARBONIO E PURO OSSIGENO ****
           DO I=12,13
             READ(3,299)  COMPO
             WRITE(*,400) COMPO
            IF(DEBUG)WRITE(18,400) COMPO
             DO K=24,87 
                READ(3,401) TEMP,(OP(L,K,I),L=1,31)
               IF(DEBUG)WRITE(18,401) TEMP,(OP(L,K,I),L=1,31)
             END DO
           END DO
        ENDIF
        ENDDO
C   ******************* CREA TABELLA DI PURO ELIO + Z ******
        DO N=1,10
        DO J=1,87  
            DO K=1,31
               OPU(K,J,11,N)=OPU(K,J,1,N)
            END DO
           IF(DEBUG)write(18,'(31f8.3)')(OPU(K,J,11,N),k=1,31)
         END DO
         IZ=.FALSE.
         CLOSE(3)
        ENDDO
      ENDIF
      IF(DEBUG)close(18)
C   ******************* FINE FASE LETTURA *********************      
      DO VF=1,NP
         IND = VF+MESH-1  !CHIMICA SHIFTATA SE MESH>1
         HYD = XX(1,IND)
         HE3 = XX(2,IND)
         HE4 = XX(3,IND)
         CAR = XX(4,IND)
         OXY = XX(6,IND)
         FEO = XX(26,IND)  !VEDERE DIFFMAU PER XX PASSATI
c************************************************************************
c****************uso il ferro come indicatore di metallicita'************
         ZMA1=FEO/ABSOL(26) ! iron abundance
         ZMA2=1.d0-HYD-HE3-HE4 ! general metallicity
         IF(ZMA2.GT.0.D0)THEN
           ZMA=DMIN1(ZMA1,ZMA2) ! min of zma1, zma2
         ELSE
           ZMA=ZMA1
         ENDIF
        if(zma.le.0.d0)write(*,'(5f15.10)')ZMA,HYD,HE3,HE4
c        pause
        ZIN(1)=dlog10(ZMA)
         DO i=2,10
           x(i)=dlog10(ZETAVAL(i))
         ENDDO
         x(1)=-30.D0 
c************************************************************************
         HE=HE3+HE4

C     ******* SE STA BRUCIANDO H SCEGLIE LE TRE TABELLE IN Y ***

         IF(HYD.GT.1.D-30) THEN
   
            DO J=1,9
             IF(HE.LE.YTAB(J).AND.HE.GE.YTAB(J+1)) THEN
              IF(J.EQ.9) THEN
               JY(1)=8
                JY(2)=9
                JY(3)=10
                JYF=3
             ELSE
                JY(1)=J
                JY(2)=J+1
                JY(3)=J+2
                JYF=3
              ENDIF
             HBURN=.TRUE.
              GO TO 1
             ENDIF
           ENDDO
         ELSE
            JY(1)   = 11
            JY(2)   = 12
            JY(3)   = 13
            JYF   =  3
            HBURN = .FALSE.
         ENDIF
  1      CONTINUE

        T63   = (T(VF)/1.0D6)**3.D0 
        T63DEL = ((T(VF)*(1.D0+DELTAT))/1.0D6)**3.D0 
         R     = CONV*DLOG(RHO(VF,1)/T63)
        TEMP  = CONV*DLOG(T(VF))
        IF(TEMP.LT.3.D0.OR.TEMP.GT.8.7D0) THEN
          WRITE(*,*)'FUORI RANGE TEMPERATURA IN OPACITA`'
          STOP
        ENDIF
        IF(TEMP.LT.4.15D0) THEN
          IF(R.LT.-7.D0.OR.R.GT.7.D0) THEN
            WRITE(*,*)'FUORI RANGE R --> LogT<4.15'
            STOP
          ENDIF
        ELSE
          IF(R.LT.-8.D0.OR.R.GT.7.D0) THEN
            WRITE(*,*)'FUORI RANGE R --> LogT>4.15'
            STOP
          ENDIF
        ENDIF
         IF(ND.GT.1) THEN
             RT=CONV*DLOG(RHO(VF,3)/T63DEL) 
             RP=CONV*DLOG(RHO(VF,2)/T63)
         ENDIF
C ******************* CERCA LE QUATTRO TEMPERATURE *********************
         IF(TEMP.LT.6.0D0) THEN
            K=1
            ESCI=.TRUE.
            DO WHILE (ESCI)
               K=K+1
              IF(K.EQ.61) ESCI=.FALSE.
              IF(K.LT.61.AND.TK(K).GT.TEMP) ESCI=.FALSE.
            END DO
            KT=K-2
            IF(KT.LT.1) KT=1
            KK=KT
            DO K=1,4
               TT(K)=TK(KK)
               KK=KK+1
            END DO
         ELSEIF(TEMP.GE.6.0D0.AND.TEMP.LT.8.1D0) THEN
            JK=60+INT((TEMP-6.0D0)/0.1D0)+1
            KT=JK-1
            IF(JK.GT.80) KT=79
            KK=KT
            DO IK=1,4
               TT(IK)=TK(KK)
               KK=KK+1
            END DO
         ELSEIF(TEMP.GE.8.1D0) THEN
            JK=82+INT((TEMP-8.1D0)/0.2D0)+1
           KT=JK-1
           IF(JK.GT.85) KT=84
            KK=KT
            DO IK=1,4
               TT(IK)=TK(KK)
               KK=KK+1
            END DO
        ENDIF
C ******************* CERCA LE QUATTRO DENSITA' ************************
        INDR0=INT((R+8.D0)/0.5D0)+1
        IF(TEMP.LE.4.15D0) THEN
          IF(INDR0.LT.4) INDR0=4
        ELSE
          IF(INDR0.LT.2) INDR0=2
        ENDIF
        IF(INDR0.GT.29) INDR0=29
c*****************interpolazione in metallicita'************************
        IF(ZIN(1).LT.X(1).OR.ZIN(1).GT.X(10)) THEN
             WRITE(*,333)ZMA,X(1),X(10),MESH
             WRITE(*,'(5D18.7)')HYD,HE3,HE4,1.D0-HYD-HE3-HE4,ZMA2
             WRITE(*,'(5D18.7)')zma1
             WRITE(*,*)'pippppppa'
             STOP
c          IF(ZIN(1).GT.X(10))ZIN(1)=X(10)
        ENDIF
c*****************controlla se punto griglia************************
           DO ki=1,10 
             if(ZIN(1).eq.X(ki)) then
               if(JY(1).eq.11)then
                jfinale=11
              else
                jfinale=JY(3)
              endif
               do lij=JY(1),jfinale
                 do k=kt,kt+3
                    do jj=INDR0-1,INDR0+2           
                       OP(jj,k,lij)=OPU(jj,k,lij,ki)
                    enddo    
                 enddo       
               enddo        
               goto 339
             end if
           ENDDO
c****** interpolazione in metallicita' lineare se Z<0.0001 **************
        IF(ZIN(1).LT.X(2)) THEN
c          write(*,*)'interpolazione lineare in opacita'
	  DX=(ZIN(1)-X(1))/(X(2)-X(1))
          DO LIJ=JY(1),JY(3) 
             DO K=KT,KT+3
                DO JJ=INDR0-1,INDR0+2 
                   OP(JJ,K,LIJ)=OPU(JJ,K,LIJ,1)+
     #                         (OPU(JJ,K,LIJ,2)-OPU(JJ,K,LIJ,1))*DX
                ENDDO    
             ENDDO      
          ENDDO
       ELSE
c************************************************************************
              if(JY(1).eq.11)then
               jfinale=11
             else
               jfinale=JY(3)
             endif
              do lij=JY(1),JY(3) 
                 do k=kt,kt+3
                    do jj=INDR0-1,INDR0+2 
                      DO IZI=1,10
                       Y(IZI)=OPU(JJ,K,LIJ,IZI)
                     ENDDO
                     CALL SPLINK(NUM,NUM,MMM,MMC,IOP)
                     OP(JJ,K,LIJ)=FVALUE(1)
                  enddo    
                 enddo      
              enddo
        ENDIF
C **********************************************************************          
 339     CONTINUE       
C **********************************************************************
C * SE GRIGLIATO NON E' CAMBIATO DALLA CHIAMATA PRECEDENTE METTE
C * FAST = TRUE
c        IF(JY1.EQ.JY1V.AND.JY2.EQ.JY2V.AND.TT(1).EQ.TTV.AND.II(1).EQ.
c    #   IIV) THEN
c           FAST=.TRUE.
c        ELSE
            FAST=.FALSE.
c        ENDIF
C **************** INTERPOLA CUBICAMENTE SU RO *************************
         DO JYY=1,JYF
         
            JOTA=JY(JYY)
            
            IF(FAST.EQV..FALSE.) THEN
               ICHI  = (JOTA-1)*31*87
               DO M=1,4
                  ITEM = ICHI + ((KT-1)-1+M)*31 + INDR0-2
                  DO N=1,4
                     B(N)=-8.d0+(INDR0-3+N)*0.5D0
                   A(N)=OP((INDR0-2+N),(KT-1+M),JY(JYY))
                  END DO
                 CALL NCUB
                 CAPR(M)=COEFF(1)*R*R*R+COEFF(2)*R*R+
     ^                   COEFF(3)*R+COEFF(4)
                 IF(ND.GT.1) THEN
                    KPDP(M)=COEFF(1)*RP*RP*RP+COEFF(2)*RP*RP+
     ^                      COEFF(3)*RP+COEFF(4)
                    KTDT(M)=COEFF(1)*RT*RT*RT+COEFF(2)*RT*RT+
     ^                      COEFF(3)*RT+COEFF(4)
                 ENDIF
                 INDI = (JYY-1)*4 + M
                 AAA(INDI)=COEFF(1)
                 BBB(INDI)=COEFF(2)
                 CCC(INDI)=COEFF(3)
                 DDD(INDI)=COEFF(4)
               END DO
            ELSE
               DO M=1,4
                  INDI = (JYY-1)*4 + M
                  CAPR(M)=AAA(INDI)*R*R*R+BBB(INDI)*R*R+
     ^            CCC(INDI)*R+DDD(INDI)
                  IF(ND.GT.1) THEN
                     KPDP(M)=AAA(INDI)*RP*RP*RP+BBB(INDI)*RP*RP+
     ^                       CCC(INDI)*RP+DDD(INDI)
                     KTDT(M)=AAA(INDI)*RT*RT*RT+BBB(INDI)*RT*RT+
     ^                       CCC(INDI)*RT+DDD(INDI)
                  ENDIF
               END DO
            ENDIF
C********  INTERPOLA CUBICAMENTE SU T E CALCOLA DK/DR E DK/DT **********
            DO J=1,ND
               DO K=1,4
                  IF(J.EQ.1) THEN
                     A(K)=CAPR(K)
                  ELSE IF(J.EQ.2) THEN
                     A(K)=KPDP(K)
                  ELSE IF(J.EQ.3) THEN
                     A(K)=KTDT(K)
                  ENDIF
                  B(K)=TT(K)
               END DO
               CALL NCUB
               CAPJ(JYY)=COEFF(1)*TEMP*TEMP*TEMP+COEFF(2)*TEMP*TEMP+
     ^                   COEFF(3)*TEMP+COEFF(4)
               IF(J.EQ.1) THEN
                  CAPS(JYY)=DEXP(COND*(CAPJ(JYY)))
               ELSE IF(J.EQ.2) THEN
                      CAPDP(JYY)=DEXP(COND*(CAPJ(JYY)))
               ELSE IF(J.EQ.3) THEN
                      TDEL=CONV*DLOG(T(VF)*(1.D0+DELTAT))
                      CATDT(JYY)=COEFF(1)*TDEL*TDEL*TDEL+COEFF(2)*
     ^                           TDEL*TDEL+COEFF(3)*TDEL+COEFF(4)
                      CATDT(JYY)=DEXP(COND*(CATDT(JYY)))
               ENDIF
            END DO
         END DO
C ********************** INTERPOLAZIONE SU CHIMICA *********************
         IF (HBURN) THEN
              DO K=1,3
               AY(K)=CAPS(K)
               AX(K)=YTAB(JY(K))
              END DO
           CAP(VF,1)=AY(1)+(HE-AX(1)) * ((AY(2)-AY(1))/(AX(2)-AX(1))+
     #    (HE-AX(2))/(AX(3)-AX(2)) * ((AY(3)-AY(1))/(AX(3)-AX(1)) - 
     #    (AY(2)-AY(1))/(AX(2)-AX(1))) )
         


            IF(ND.GT.1) THEN

              DO K=1,3
               AY(K)=CATDT(K)
               AX(K)=YTAB(JY(K))
              END DO
          CAP(VF,3)=AY(1)+(HE-AX(1)) * ((AY(2)-AY(1))/(AX(2)-AX(1)) +
     #    (HE-AX(2))/(AX(3)-AX(2)) * ((AY(3)-AY(1))/(AX(3)-AX(1)) - 
     #    (AY(2)-AY(1))/(AX(2)-AX(1))) )

              DO K=1,3
               AY(K)=CAPDP(K)
               AX(K)=YTAB(JY(K))
              END DO
          CAP(VF,2)=AY(1)+(HE-AX(1)) * ((AY(2)-AY(1))/(AX(2)-AX(1)) +
     #    (HE-AX(2))/(AX(3)-AX(2)) * ((AY(3)-AY(1))/(AX(3)-AX(1)) - 
     #    (AY(2)-AY(1))/(AX(2)-AX(1))) )

            ENDIF
C ****************** QUADRATIC INTERPOLATION ON Y **********************
         ELSE
C ************ SE H=0 COMBINAZIONE LINEARE SU Y - C - O ****************
            SUMMA   = (HE4+CAR+OXY)
            CAP(VF,1) = (CAPS(1)*HE4+CAPS(2)*CAR+CAPS(3)*OXY)/SUMMA
            IF(ND.GT.1) THEN
              CAP(VF,3) = (CATDT(1)*HE4+CATDT(2)*CAR+CATDT(3)*OXY)/SUMMA
              CAP(VF,2) = (CAPDP(1)*HE4+CAPDP(2)*CAR+CAPDP(3)*OXY)/SUMMA
            END IF
         ENDIF
        DO I=1,3
           JYYV(I)=JY(I)
        ENDDO
         TTV =TT(1)
         IIV =INDR0
      END DO
      RETURN
      END 
      
      SUBROUTINE NSUPERA(RO,P,T,CAP,Q,CP,ADIAB,GI,ALPHA,GRSAD,ACCO,RAD) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'consts.2p1'
C***********************************************************************
C ===== CHECK IF PHYSICAL QUANTITIES ARE REALISTIC
      IF(P.LE.0.D0.OR.RO.LE.0.D0.OR.GI.EQ.0.D0.OR.T.LE.0.D0) 
     # STOP'MUOIO IN SUPERA'
      HP=P/(RO*GI) 
C* HP DEVE ESSERE < O = DELLA DISTANZA (IN PRESS.) DAL BORDO EST. CONVET 
C     HPMAX=(P/PTOP-1.)/1.72 
C     HP=DMIN1(HP,HPMAX) 
C*********************************************************************** 
      AML=ALPHA*HP 
      STAB=RAD-ADIAB 
      E=(DSQRT(Q)*CP*CAP*GI*RO**2.5*AML**2)/(48.D0*DSQRT(2.D0) 
     %*STBOLZ*P**.5*T**3) 
      A0=9.D0/4.D0 
      B=(E**2*STAB/A0)**(1.D0/3.D0) 
      IF(B.EQ.0) THEN 
         WRITE(*,*)'B=0'
         WRITE(2,*)'B=0'
         STOP 
      ENDIF 
      DELTA=1.D-12 
      SDELTA=.1D0*DELTA ! STEP SIZE = 10% OF DELTA
      IF(B.LT.1.D0) THEN 
        A0B2=A0*B*B 
        Z1=(A0B2**3)*(1.D0-3.D0*B*A0B2-(3.D0-(9.D0/A0))*A0B2**3) 
        Z=Z1 
        Y0=Z**(1.D0/3.D0) 
        L=0 
   10   Y1=Y0 ! ITERATION LOOP = '10'
        GUESS1=Y1+B*Y1*Y1+A0B2*Y1*Y1*Y1-A0B2 
        VAR=DABS(GUESS1) 
        IF(VAR.LE.DELTA) GO TO 20 
        Y2=Y1*(1.D0+SDELTA) 
        GUESS2=Y2+B*Y2*Y2+A0B2*Y2*Y2*Y2-A0B2 
        IF((GUESS2-GUESS1).EQ.0.D0) THEN
           WRITE(2,*)'NSUPERA B<1'
           WRITE(*,*)'NSUPERA B<1'
           SDELTA=SDELTA*2.D0 ! CHANGE STEP SIZE TO 1% OF DELTA
           GO TO 10 
        ENDIF 
        Y0=Y1-GUESS1*(Y2-Y1)/(GUESS2-GUESS1) 
        IF(L.GT.100) GO TO 30 
        L=L+1 
        GO TO 10 
      ELSE 
        FAC1=(1.D0+B)/A0/B**2 
        FAC2=(1.D0+2.D0*B)/3.D0/(1.D0+B) 
        Z2=1.D0-FAC1*(1.D0-FAC1*FAC2+(1.D0/9.D0)*(9.D0*FAC2**2-1.D0)* 
     #     FAC1**2) 
        Z=Z2 
        Y0=Z**(1.D0/3.D0) 
        L=0 
   12   Y1=Y0 
        GUESS1=1.D0-Y1*Y1*Y1-(Y1+B*Y1*Y1)/A0/B**2 
        VAR=DABS(GUESS1) 
        IF(VAR.LE.DELTA) GO TO 20 
        Y2=Y1*(1.D0+SDELTA) 
        GUESS2=1.D0-Y2*Y2*Y2-(Y2+B*Y2*Y2)/A0/B**2 
        IF((GUESS2-GUESS1).EQ.0.D0) THEN 
           WRITE(2,*)'NSUPERA B>1'
           WRITE(*,*)'NSUPERA B>1'
           SDELTA=SDELTA*2.D0 
           GO TO 12 
        ENDIF 
        Y0=Y1-GUESS1*(Y2-Y1)/(GUESS2-GUESS1) 
C        WRITE(*,*)'B>1 ',L,RO,T,GUESS1 
        IF(L.GT.100) GO TO 30 
        L=L+1 
        GO TO 12 
      ENDIF 
   20 CONTINUE 
      Z=Y0**3 
      GRSAD=(1.D0-Z)*RAD+Z*ADIAB 
      CV=B*Y0 
      GRAD1=GRSAD-(CV/E)**2 
      ACCO=GI*Q**.5*(RO/P)**.5*AML*(GRSAD-GRAD1)**.5/(2.D0*DSQRT(2.D0)) 
      RETURN 
   30 CONTINUE 
      IF(VAR.GT.1.D-10) WRITE(*,100) 
  100 FORMAT(1X,'ATTENZIONE MIX-LEN NON CONVERGE') 
      END           

      SUBROUTINE SIGMANEW (IFASE,NPRO,T,RO,YF,CROS)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      DIMENSION YF(MAXNE),CROS(MAXPR,2)
      T9=T*1.D-9
      WLOT9=DLOG(T9)
      T912=DEXP(0.5D0*WLOT9) ! = exp(log(T9^0.5))
      T913=DEXP(1.D0/3.D0*WLOT9) ! = exp(log(T9^(1/3)))
      T9=DEXP(WLOT9)
      T92 =T9*T9
      T93 =T92*T9
      T94 =T93*T9
      T95 =T94*T9
      T923=T913*T913
      T932=T9*T912
      T943=T9*T913
      T953=T9*T923
c****************************************************************
C  16-10-02 
C  INSERITE LE SEZIONI D'URTO NACRE (Angulo et al. 1999) PER
C    LE COMBUSTIONI DI H e He SALVO ALTRA INDICAZIONE. 
C    inserito rate Kunz et al. 2002 per C12(a,g)O16 - 
C  I RATE COMMENTATI SONO QUELLI "VECCHI".
c****************************************************************
      DO I=1,NPRO
         CROSSO=0.D0
c        
C######################################################################
C***********************  CATTURE PROTONICHE  *************************
C######################################################################
         IF(ICROS(IFASE,I).EQ.1) THEN
C************************(1)* P(P,E+NU)D ******************************
          CROSSO=4.08D-15/T923*DEXP(-3.381D0/T913)*
     #       (1.D0+3.82D0*T9+1.51D0*T92+0.144D0*T93-1.14D-02*T94)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.2) THEN
C********************(2)*  D + P => HE3 + G  **************************  
           IF(T9.LE.0.11D0)THEN
             CROSSO=(1.81D03/T923)*DEXP(-3.721D0/T913)*
     #           (1.D0+14.3D0*T9-90.5D0*T92+395.D0*T93)
           ELSE
             CROSSO=(2.58D03/T923)*DEXP(-3.721D0/T913)*
     #           (1.D0+3.96D0*T9+0.116D0*T92)
           ENDIF
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.3) THEN
******************(3)*  LI7 + P => HE4 + HE4  *************************
            T9L=T9**0.576D0
            CROSSO=(7.20D08/T923)*DEXP(-8.473D0/T913-T92/42.25D0)*
     #      (1.D0+1.05D0*T9-0.653D0*T92+0.185D0*T93-2.12D-02*T94+
     #       9.30D-04*T95)+9.85D06*T9L*DEXP(-10.415D0/T9)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.4) THEN
**************(4)*  BE7 + P => BE8 + E + NU  => HE4 + HE4  ************
            CROSSO=(2.61D05/T923)*DEXP(-10.264D0/T913)*
     #      (1.D0-5.11D-02*T9+4.68D-02*T92-6.60D-03*T93+3.12D-04*T94)+
     #       (2.05D03/T932)*DEXP(-7.345D0/T9)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.5) THEN
C*********************(5)* C12(P,G)N13=>C13 ***************************
            CROSSO=(2.D07/T923)*DEXP(-13.692D0/T913-T92/0.2116D0)*
     #      (1.D0+9.89D0*T9-59.8D0*T92+266.D0*T93)+(1.D05/T932)*
     #      DEXP(-4.913D0/T9)+(4.24D05/T932)*DEXP(-21.62D0/T9)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.6) THEN
*********************(6)* C13(P,G)N14 *********************************
            T9C13=T9**(-0.864D0)
            CROSSO=(9.57D07/T923*(1.D0+3.56D0*T9)*
     #             DEXP(-13.72D0/T913-T92)+
     #             (1.5D06/T932)*DEXP(-5.930D0/T9)+
     #             6.83D05*T9C13*DEXP(-12.057D0/T9))*
     #            (1.D0-2.07D0*DEXP(-37.938D0/T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.7) THEN
C*********************(7)* N14(P,G)O15=>N15 ***************************
c            T9N=T9**0.380D0
c            CROSSO=(4.83D07/T923)*DEXP(-15.231D0/T913-T92/0.64D0)*
c     #      (1.D0-2.D0*T9+3.41D0*T92-2.43D0*T93)+
c     #       (2.36D03/T932)*DEXP(-3.01D0/T9)+
c     #       6.72D03*T9N*DEXP(-9.53D0/T9)
CC********************FORMICOLA et AL. 2003  ****************************
          T9N=T9**0.0682D0
          CROSSO=(3.12D07/T923)*DEXP(-15.193D0/T913-T92/(0.486D0**2))*
     #      (0.782D0-1.5D0*T9+17.97D0*T92-3.32D0*T93)+
     #      (2.11D03/T932)*DEXP(-2.998D0/T9)+
     #       8.42D02*T9N*DEXP(-4.891D0/T9)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.8) THEN
C********************(8)** N15(P,G)O16 ********************************
            T9N152=T9**0.095D0
            IF(T9.LE.3.5D0)THEN
               CROSSO=(1.08D09/T923)*DEXP(-15.254D0/T913-T92/0.1156D0)*
     #         (1.D0+6.15D0*T9+16.4D0*T92)+
     #          (9.23D03/T932)*DEXP(-3.597D0/T9)+
     #          (3.27D06/T932)*DEXP(-11.024D0/T9)
            ELSE
              CROSSO=3.54D04*T9N152*DEXP(-2.306D0/T9)
            ENDIF
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.9) THEN
C*********************(9)* N15(P,A)C12 ********************************
            T9N15=T9**0.917
            IF(T9.LE.2.5D0) THEN
               CROSSO=1.12D12/T923*DEXP(-15.253D0/T913-T92/0.0784D0)*
     #         (1.D0+4.95D0*T9+143D0*T92)+
     #          (1.01D08/T932)*DEXP(-3.643/T9)+
     #          (1.19D09/T932)*DEXP(-7.406/T9)
            ELSE
              CORSSO=4.17D07*T9N15*DEXP(-3.292D0/T9)
            END IF
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.10) THEN
C********************(10)** O16(P,G)F17 (==>O17)***********************
            T9O16=T9**(-0.82D0)
            CROSSO=(7.37D07*T9O16*DEXP(-16.696D0/T913))*
     #             (1.D0+202.D0*DEXP(-70.348D0/T9-0.161D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.11) THEN
C*********************(11)* O17(P,A)N14 *******************************
            T9O17=T9**1.591D0
            T9O172=T9**0.95D0
            IF(T9.LE.6.D0)THEN
              CROSSOG=9.20D08/T923*DEXP(-16.715/T913-T92/0.0036D0)*
     #              (1.D0-80.31D0*T9+2211.D0*T92)+
     #              (9.13D-04/T932)*DEXP(-0.7667D0/T9)+
     #              (9.68D0/T932)*DEXP(-2.083D0/T9)+
     #              (8.13D06/T932)*DEXP(-5.685D0/T9)+
     #               1.85D06*T9O17*DEXP(-4.848D0/T9)
             ELSE
              CROSSOG=8.73D06*T9O172*DEXP(-7.508D0/T9)
              write(*,*)'t9>6'
             ENDIF
             CROSSO=CROSSOG*
     #             (1.D0+1.033D0*DEXP(-10.034D0/T9-0.165D0*T9))
C********************LANDRE` et AL. 1990  ****************************
c      PIPPO1=0.2D0
c      PIPPO2=1.0D0
c      CROSSO=1.53D7/T923*DEXP(-16.712D0/T913-(T9/0.565D0)**2.D0)*
c     #(1.D0+.025D0*T913+5.39D0*T923+0.94D0*T9+13.5D0*T943+5.98D0*T953)
c     #+2.92D6*T9*DEXP(-4.247D0/T9)
c     #+1.78D5/T923*DEXP(-16.67D0/T913)/(0.479D0*T923+0.00312D0)**2
c     #+PIPPO1*
c     #(2.8D11*T9*DEXP(-16.67D0/T913-(T9/0.04D0)**2)+2.94D-3/T932*
c     #DEXP(-0.767D0/T9))
c     #+PIPPO2*98.D0/T932*DEXP(-2.077D0/T9)      
cC**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.12) THEN
C ***********(12)**   NE20 + P => NA21 + G   **************************
            T9Ne203=T9**(-1.84D0)
            T9Ne204=T9**(-0.641D0)
            CROSSO=(2.35D07*T9Ne203*DEXP(-19.451D0/T913)*
     #            (1.D0+10.8D0*T9)+18.D0/T932*DEXP(-4.247D0/T9)+
     #            (9.83D0/T932)*DEXP(-4.619D0/T9)+
     #             6.67D04*T9Ne204*DEXP(-11.922D0/T9))*
     #            (1.D0-7.929D0*DEXP(-20.108D0/T9-0.327D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.13) THEN
C ************(13)*   NE21 + P => NA22 + G   **************************
            T9Ne21=T9**(-0.128D0)
            T9Ne211=T9**0.42D0
            IF(T9.LE.2.D0)THEN
              CROSSOG=(4.68D08/T923)*DEXP(-19.465/T913-T92/0.04D0)+
     #               (8.18D-04/T932)*DEXP(-1.085D0/T9)+
     #               (6.11D0/T932)*DEXP(-1.399D0/T9)+
     #               (1.34D04/T932)*DEXP(-3.009D0/T9)+
     #               1.26D05*T9Ne21*DEXP(-4.962D0/T9)
            ELSE
              CROSSOG=3.04D04*T9Ne211*DEXP(-2.65D0/T9)
            END IF
            CROSSO=CROSSOG*(1.D0-0.708D0*DEXP(-3.851D0/T9-0.156D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.14) THEN
C ************(14)*   NE22 + P => NA23 + G   **************************
            T9Ne22=T9**0.725D0
            T9Ne222=T9**0.816D0
            IF(T9.LE.2.D0)THEN
              CROSSOG=(1.11D-09/T932)*DEXP(-0.422D0/T9)+
     #               (6.83D-05/T932)*DEXP(-0.81D0/T9)+
     #               (9.76D-03/T932)*DEXP(-1.187D0/T9)+
     #               (1.06D-01/T932)*DEXP(-1.775D0/T9)+
     #               8.51D04*T9Ne22*DEXP(-4.315/T9)
            ELSE
              CROSSOG=6.3D04*T9Ne222*DEXP(-3.91D0/T9)
            END IF
            CROSSO=CROSSOG*
     #            (1.D0-1.41D0*DEXP(-14.651D0/T9-0.02D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.15) THEN
C **************(15)*     NA23 + P => NE20 + A  ***********************
            T9Na23=T9**(-1.48D0)
            T9Na232=T9**1.456D0
            T9Na233=T9**1.291D0
            IF(T9.LE.5.D0)THEN
               CROSSO=(8.39D09/T923)*DEXP(-20.77D0/T913-T92/0.01)*
     #             (1.D0+45.2D0*T9)+3.09D-13/T932*DEXP(-0.42D0/T9)+
     #             (8.12D-03/T932)*DEXP(-1.601/T9)+
     #             (4.37D0/T932)*DEXP(-1.934/T9)+
     #             7.5D03*T9Na23*DEXP(-3.15D0/T9)+
     #             1.05D06*T9Na232*DEXP(-4.482/T9)
            ELSE
               CROSSO=3.93D06*T9Na233*DEXP(-9.277D0/T9)
            END IF
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.16) THEN
C **************(16)*     NA23 + P => MG24 + G  ***********************
            T9Na234=T9**1.112D0
            IF(T9.LE.5.D0)THEN
               CROSSOG=9.55D07/T923*DEXP(-20.77D0/T913-T92/0.09D0)*
     #                (1.D0-10.8D0*T9+61.08D0*T92)+
     #                 (8.2D-02/T932)*DEXP(-1.601D0/T9)+
     #                 (85.2D0/T932)*DEXP(-2.808D0/T9)+
     #                 (1.7D04/T932)*DEXP(-3.458D0/T9)+
     #                 5.94D04*DEXP(-5.734D0/T9)
            ELSE
              CROSSOG=5.6D03*T9Na234*DEXP(-2.337/T9)
            END IF
            CROSSO=CROSSOG*(1.D0-0.56D0*DEXP(-5.119D0/T9-0.05D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.17) THEN
C **************(17)* MG24 + P => MG25 + G (CFZ 88) *******************
            CROSSO=(5.6D08/T923*DEXP(-22.019D0/T913)*
     #       (1.D0+0.019D0*T913-0.173D0*T923-0.023D0*T9)+
     #       (1.48D03/T932)*DEXP(-2.484/T9)+4.D03*DEXP(-4.18D0/T9))/
     #       (1.D0+5.D0*DEXP(-15.882D0/T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.18) THEN
C*************(18)* MG25 + P => AL26 + G (AL26=>MG26) *****************
            T9Mg25=T9**0.647D0
            T9Mg252=T9**1.362D0
            T9Mg253=T9**1.262D0
            IF(T9.LE.2.D0)THEN
               CROSSO1=(3.07D-16/T932)*DEXP(-0.435D0/T9)+
     #           (3.7D-08/T932)*DEXP(-0.673D0/T9)+
     #           (1.6D-05/T932)*DEXP(-1.074D0/T9)+
     #           1.27D04*T9Mg25*DEXP(-3.055D0/T9)
               CROSSO2=(8.15D-17/T932)*DEXP(-0.435D0/T9)+
     #           (8.68D-09/T932)*DEXP(-0.673D0/T9)+
     #           (2.82D-06/T932)*DEXP(-1.074D0/T9)+
     #           3.48D03*T9Mg252*DEXP(-2.906/T9)
            ELSE
              CROSSO1=8.75D03*T9*DEXP(-2.997D0/T9)
              CROSSO2=3.91D03*T9Mg253*DEXP(-3.229D0/T9)
            END IF
            CROSSO=(CROSSO1+CROSSO2)*
     #             (1.D0-0.352D0*DEXP(-7.221D0/T9+0.068D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.19) THEN
C*************(19)* MG26 + P => AL27 + G ******************************
            T9Mg26=T9**(-1.565D0)
            T9Mg262=T9**0.215D0
            T9Mg263=T9**1.068D0
            IF(T9.LE.3.5D0)THEN
               CROSSOG=(8.54D-12/T932)*DEXP(-0.605D0/T9)+
     #            (1.93D-06/T932)*DEXP(-1.044D0/T9)+
     #            (9.67D-03/T932)*DEXP(-1.726D0/T9)+
     #            (9.50D04/T932)*DEXP(-3.781D0/T9)+
     #            10.2D0*T9Mg26*DEXP(-2.521D0/T9)+
     #            7.07D04*T9Mg262*DEXP(-3.947D0/T9)
            ELSE
               CROSSOG=3.95D04*T9Mg263*DEXP(-4.990D0/T9)
            END IF
            CROSSO=CROSSOG*
     #            (1.D0-1.259D0*DEXP(-20.076D0/T9+0.069D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.20) THEN
C************(20)** AL27 + P => SI28 + G ******************************
            T9Al27=T9**(-0.2D0)
            T9Al272=T9**1.12D0
            T9Al273=T9**0.251D0
            T9Al274=T9**0.549D0
            IF(T9.LE.6.D0)THEN
               CROSSOG=(2.51D-11/T932)*DEXP(-0.839D0/T9)+
     #            48.2D0*T9Al27*DEXP(-2.223D0/T9)+
     #            1.76D03*T9Al272*DEXP(-3.196/T9)+
     #            3.25D04*T9Al273*DEXP(-5.805D0/T9)
            ELSE
               CROSSOG=1.62D05*T9Al274*DEXP(-17.222/T9)
            END IF
            CROSSO=CROSSOG*
     #            (1.D0-0.669D0*DEXP(-10.426D0/T9+0.008D0*T9))
C**********************************************************************
C######################################################################
C************************  CATTURA ALFA  ******************************
C######################################################################
         ELSEIF(ICROS(IFASE,I).EQ.21) THEN
C*******************(21)*  HE3 + HE4 => BE7 + G  **********************
c            CROSSO=(5.46D06/T923)*DEXP(-12.827D0/T913)*
c     #         (1.D0-0.307D0*T9+8.81D-02*T92-1.06D-02*T93+4.46D-04*T94)
C******************Cyburt e Devis (2008) Phys. rev. C *****************
C**********************************************************************
            CROSSO=DEXP(15.609867D0-12.827077d0/T913+2.d0/3.d0*dlogT9)*
     #             (1.d0-0.020478d0*t923+0.211995*T943)/
     #             (1.d0-0.255059d0*t923+0.338573*T943)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.22) THEN
C ***************(22)*  C12 + A => O16 + G  (NACRE) *******************
c            CROSSO1=6.66D07/T92*DEXP(-32.123/T913-T92/21.16D0)*
c     #             (1.D0+2.54D0*T9+1.04D0*T92-0.226D0*T93)+
c     #              (1.39D03/T932)*DEXP(-28.930/T9)
c            CROSSO2=(6.56D07/T92)*DEXP(-32.123/T913-T92/1.69D0)*
c     #             (1.D0+9.23D0*T9-13.7D0*T92+7.4D0*T93)
c            CROSSO3=19.2D0*T92*DEXP(-26.9D0/T9)
c            CROSSO=CROSSO1+CROSSO2+CROSSO3
c************************************************************************
c************** KUNZ et al. 2002*****************************************
c************************************************************************
c      CROSSO=(1.21D8/(T9*T9*(1.d0+6.06D-2/T923)**2))*
c     #          DEXP(-(32.12D0/T913)-(T9/1.7D0)**2)+
c     #          7.4D8/(T9*T9*(1.d0+0.47D0/T923)**2)*DEXP(-32.12D0/T913)+
c     #          (1.53D4/T923)*(1.d0+2.0D6*T913)*DEXP(-38.534D0/T913)
c************** Hammer et al. 2005 **************************************
cc************************************************************************
      CROSSO=(1.51D8/(T9*T9*(1.d0+6.66D-2/T923)**2))*
     #          DEXP(-(32.12D0/T913)-(T9/1.03D0)**2)+
     #          1.11D9/(T9*T9*(1.d0+0.735D0/T923)**2)*
     #          DEXP(-32.12D0/T913)+
     #          (1.62D4/T923)*(1.d0+2.19D6*T913)*DEXP(-38.814D0/T913)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.23) THEN
C ************(23)*   O16 + A => NE20 + G   ***************************
            T9O162=T9**2.966D0
            CROSSO=(2.68D10/T923)*DEXP(-39.76D0/T913-T92/2.56D0)+
     #             (51.1D0/T932)*DEXP(-10.32D0/T9)+
     #             (616.1D0/T932)*DEXP(-12.2D0/T9)+
     #             0.41D0*T9O162*DEXP(-11.9D0/T9)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.24) THEN
C ***********(24)*    N14 + A => F18 + G  (F18 ==> 018) ***************
            T9N14=0.344D0
            T9N142=T9**1.567D0
            IF(T9.LE.2.D0) THEN
               CROSSOG=(7.93D11/T923)*DEXP(-36.035D0/T913-T92/0.0049)+
     #               (1.85D-10/T932)*DEXP(-2.75D0/T9)+
     #               (2.62D0/T932)*DEXP(-5.045D0/T9)+
     #               (2.93D03*T9N14)*DEXP(-10.561D0/T9)
            ELSE
               CROSSOG=1.52D02*T9N142*DEXP(-6.315D0/T9)
            END IF
            CROSSO=CROSSOG*(1.D0-0.34D0*DEXP(-26.885D0/T9-0.012D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.25) THEN
C ***********(25)*    O18 + A => NE22 + G     *************************
            T9O18=T9**(-0.221D0)
            IF(T9.LE.6.D0) THEN
              CROSSOG=(1.95D-13/T932)*DEXP(-2.069D0/T9)+
     #              (1.56D-02/T932)*DEXP(-4.462D0/T9)+
     #              (10.1D0/T932)*DEXP(-6.391D0/T9)+
     #              (44.1D0/T932)*DEXP(-7.389D0/T9)+
     #              (3.44D05*(T9**(-0.5D0)))*DEXP(-22.103D0/T9)
            ELSE
              CROSSOG=3.31D05*T9O18*DEXP(-24.99D0/T9)
            END IF
            CROSSO=CROSSOG*
     #            (1.D0-1.411D0*DEXP(-20.533D0/T9-0.038D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.26) THEN
C ************(26)*   NE20 + A => MG24 + G   **************************
            T9Ne20=T9**(-0.532D0)
            T9Ne202=T9**2.229D0
            IF(T9.LE.1.D0) THEN
              CROSSO=8.72D0*T9Ne20*DEXP(-8.995D0/T9)
            ELSE
              CROSSO=3.74D02*T9Ne202*DEXP(-12.681D0/T9)
            END IF
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.27) THEN
C ************(27)*   NE22 + A => MG25 + N    *************************
            T9Ne226=T9**0.83D0
            T9Ne227=T9**(2.78D0)
            T9Ne228=T9**0.892D0
            T9Ne229=T9**2.879D0
            IF(T9.LE.2.D0)THEN
              CROSSOG=7.4D0*DEXP(-7.79D0/T9)+
     #               1.30-04*T9Ne226*DEXP(-5.52D0/T9)+
     #               9.41D03*T9Ne227*DEXP(-11.7D0/T9)+
     #               8.59D06*T9Ne228*DEXP(-24.4D0/T9)
            ELSE
              CROSSOG=1.51D05*T9Ne229*DEXP(-16.717D0/T9)
            END IF
            CROSSO=CROSSOG*
     #            (1.D0+2.674D0*DEXP(-15.025D0/T9-0.321D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.28) THEN
C ************(28)*   NE22 + A => MG26 + G    *************************
            T9Ne223=T9**(-1.064D0)
            T9Ne224=T9**(-2.556D0)
            T9Ne225=T9**3.322D0
            IF(T9.LE.1.25D0)THEN
              CROSSOG=3.55D-09/T932*DEXP(3.927D0/T9)+
     #               7.07D-01*T9Ne223*DEXP(-7.759D0/T9)+
     #               1.27D-03*T9Ne224*DEXP(-6.555D0/T9)
            ELSE
              CROSSOG=1.76D0*T9Ne225*DEXP(-12.412D0/T9)
            END IF
            CROSSO=CROSSOG*(1.D0-0.005D0*DEXP(-5.109D0/T9+0.373D0*T9))
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.29) THEN
C **************(29)* MG24 + A => SI28 + G (CFZ88)*********************
            GT9=1.D0+5.D0*DEXP(-15.882D0/T9)
            CROSSO=((4.78D+01/T932*DEXP(-13.506/T9)+
     #            2.38D3/T932*DEXP(-15.218/T9)+2.47D+02*
     #            T932*DEXP(-15.147/T9)+0.01D0*1.72D-09/
     #            T932*DEXP(-5.028/T9)+1.25D-03/T932*
     #            DEXP(-7.929/T9)+2.43D+01/T9*DEXP(-11.523/T9))/GT9)
C**********************************************************************
C######################################################################
C***************************  CORPI IDENTICI  *************************
C######################################################################
         ELSEIF(ICROS(IFASE,I).EQ.30) THEN
C******************(30)*  HE3 + HE3 => HE4 + 2P  **********************
            CROSSO=5.59D10/T923*DEXP(-12.277/T913)*
     #            (1.D0-0.135D0*T9+2.54D-02*T92-1.29D-03*T93)
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.31) THEN
C **************(31)*  C12 + C12 ----- NE20 + A  (CFZ88)   ************
            T9A=T9/(1.0D0+0.067*T9)
            T9A13=T9A**0.3333D0
            T9A23=T9A**0.6667D0
            T9A56=T9A**0.8333D0
            IF(T9. LT. 3.D0) THEN
               BRA=0.5D0
            ELSE
               BRA=0.65D0
            ENDIF
            SCC=(1.26D27*T9A56/T932*DEXP(-84.165/T9A13)
     #         /(DEXP(-0.01*T9A**4)+5.56D-03*DEXP(1.685*T9A23)))
            CROSSO=BRA*SCC
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.32) THEN
C **************(32)*  C12 + C12 ----- NA23 + P  (CFZ88)   ************
            T9A=T9/(1.0D0+0.067*T9)
            T9A13=T9A**0.3333D0
            T9A23=T9A**0.6667D0
            T9A56=T9A**0.8333D0
            IF(T9. LT. 3.D0) THEN
               BRA=0.5D0
            ELSE
               BRA=0.65D0
            ENDIF
            SCC=(1.26D27*T9A56/T932*DEXP(-84.165/T9A13)
     #         /(DEXP(-0.01*T9A**4)+5.56D-03*DEXP(1.685*T9A23)))
            CROSSO=(1.D0-BRA)*SCC
C**********************************************************************
         ELSEIF(ICROS(IFASE,I).EQ.33) THEN
C *************(33)* O16 + O16 => S32 + G  (CFZ 88)  ******************
            CROSSO=(7.10D+36/T923*DEXP(-135.93D0/T913-0.629D0*
     #             T923-0.445D0*T943+0.0103D0*T9**2))
C**********************************************************************
C######################################################################
C***************************   3 ALFA  ********************************
C######################################################################
         ELSEIF(ICROS(IFASE,I).EQ.34) THEN
C ***********(34)*    HE4 + 2A => C12 + G *****************************
            T93A=T9**(-0.65D0)
            CROSSO1=(2.43D09/T923)*DEXP(-13.49D0/T913-T92/0.0225D0)*
     #             (1.D0+74.5D0*T9)+(6.09D05/T932)*DEXP(-1.054D0/T9)
            CROSSO2=(2.76D07/T923)*DEXP(-23.57D0/T913-T92/0.16D0)*
     #             (1.D0+5.47D0*T9+326.D0*T92)+
     #              (130.7D0/T932)*DEXP(-3.338D0/T9)+
     #              (2.51D4/T932)*DEXP(-20.307D0/T9)
            IF(T9.LE.0.03D0) THEN
              CROSSO=CROSSO1*CROSSO2*3.07D-16*(1.D0-29.1D0*T9+
     #               1308.D0*T92)
            ELSE
              CROSSO=CROSSO1*CROSSO2*3.44D-16*(1.D0+0.0158D0*T93A)
            END IF
C**********************************************************************
C######################################################################
C***********************  DECADIMENTO BETA -  *************************
C######################################################################
         ELSEIF(ICROS(IFASE,I).EQ.35) THEN
C**************(35)*  BE7 => LI7 + E + NU  **  (CFZ 88)  **************
            CROSSO=1.34D-10/T912*(1.D0-0.537D0*T913+3.86D0*T923+
     #             0.0027D0/T9*DEXP(2.515D-03/T9))
C**********************************************************************
C######################################################################
C***********************  FOTODISINTEGRAZIONE  ************************
C######################################################################
         ELSEIF(ICROS(IFASE,I).EQ.36) THEN
C ************(36)*   NE20 + G => O16 + A  (CFZ 88)  ******************
            CROSSO=5.65D+10*T932*DEXP(-54.937D0/T9)
C**********************************************************************
         END IF
         CROS(I,1)=CROSSO
      END DO
      RETURN
      END

            
            
