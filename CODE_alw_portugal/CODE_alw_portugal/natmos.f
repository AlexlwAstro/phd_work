C*****************************************************************************
C  QUESTA ROUTINE CONTIENE: CONDCON - LEGGI - NATMOS
C**********************************************************************
      SUBROUTINE CONDCON(VLUM,TEP,VMAS,ZZZ,PT,TGA)
      PARAMETER(NUMT=60,NUMG=7)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'consts.2p1'
      COMMON/DATIBOUND/T(NUMT,NUMG),G(NUMT,NUMG),AAX(NUMT,NUMG,4),
     #NT,NG 
      COMMON/CUBICA/A,B,COEFF
      DIMENSION A(4),B(4),COEFF(4),AFI(2),AX(NUMG,2),ARA(2,6)
      DATA INI/0/
C***********************************************************
      IF(INI.EQ.0) THEN
c       CALL LEGGI(ZZZ)
       INI=1
      ENDIF
C***********************************************************
      IPIP=0
C***********************************************************
      VVL=10.D0**VLUM
      AT=10.D0**TEP
      AG=DLOG10((NGATCOST*VMAS*(AT**4))/VVL)
C**********SE ANDIAMO FUORI GRIGLIATO SI FERMA**************
C if temperature/gravity value lies outside current grid limits
C***********************************************************
      IF(AT.GT.T(NT,1)) THEN
      WRITE(*,*) 'FUORI RANGE TEMP. HIGH-->',AT
      STOP
      ENDIF
      IF(AT.LT.T(1,1)) THEN
      WRITE(*,*) 'FUORI RANGE TEMP. LOW -->',AT
      STOP
      ENDIF
C===========================================================
C      IF(AG.GT.5.8D0) THEN
      IF(AG.GT.G(1,NG)) THEN
        WRITE(*,*) 'FUORI RANGE GRAV. HIGH-->',AG
        STOP
      ENDIF
C      IF(AG.LT.3.5D0) THEN
      IF(AG.LT.G(1,1)) THEN
        WRITE(*,*) 'FUORI RANGE GRAV. LOW-->',AG
        STOP
      ENDIF
C*****************FINE CONTROLLO****************************
      DO 2 K=1,NT
      IF(AT.EQ.T(K,1)) THEN
      DO 1 J=1,NG
       DO 3 I=1,2
        AX(J,I)=AAX(K,J,I)
   3  CONTINUE
   1  CONTINUE
      IPIP=1
      GO TO 11 
      ELSE
      GO TO 2
      ENDIF
   2  CONTINUE
  11  CONTINUE 
      IF(IPIP.EQ.1) GO TO 33
C*********************************************************************
      DO 7 I=1,2
      DO 4 J=1,NG
      DO 5 K=1,NT-1
      IF(AT.GT.T(K,J).AND.AT.LT.T(K+1,J)) THEN
      INI=K
      IF(K.EQ.1) INI=2
      IF(K.EQ.NT-1) INI=NT-2
      GO TO 30
      ELSE
      GO TO 5
      ENDIF
   5  CONTINUE
  30  CONTINUE
      DO 6 L=1,4
      B(L)=T(INI+L-2,J)
      A(L)=AAX(INI+L-2,J,I)
   6  CONTINUE
      CALL NCUB
      AX(J,I)=COEFF(1)*AT*AT*AT+COEFF(2)*AT*AT+COEFF(3)*AT+COEFF(4)
C==================================================================
   4  CONTINUE
   7  CONTINUE
C**********************************************************
  33  CONTINUE
C**********************************************************
      DO 55 K=1,NG
      IF(AG.EQ.G(1,K)) THEN      
       TGA=AX(K,1)
       PT=AX(K,2)
       GO TO 333
      ELSE
       GO TO 55
      ENDIF
  55  CONTINUE           
C==========================================================
      DO 44 J=1,NG-1
      IF(AG.GT.G(1,J).AND.AG.LT.G(1,J+1)) THEN
      IGI=J
      IF(J.EQ.1) IGI=2
      IF(J.EQ.NG-1) IGI=NG-2
      GO TO 65
      ELSE
      GO TO 44
      ENDIF
  44  CONTINUE    
  65  CONTINUE
C==========================================================
      DO 77 I=1,2  
      DO 88 L=1,4
      B(L)=G(1,IGI+L-2)
      A(L)=AX(IGI+L-2,I)
  88  CONTINUE
      CALL NCUB
      AFI(I)=COEFF(1)*AG*AG*AG+COEFF(2)*AG*AG+COEFF(3)*AG+COEFF(4)
  77  CONTINUE
C**********************************************************
      PT=AFI(2)
      TGA=AFI(1)
C==========================================================      
 333  CONTINUE    
      RETURN
      END

      SUBROUTINE LEGGI(ZZZ)
      PARAMETER(NUMT=60,NUMG=7)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80 TESTO
      COMMON/DATIBOUND/T(NUMT,NUMG),G(NUMT,NUMG),AAX(NUMT,NUMG,4),
     #NT,NG 
C====================================================================
      OPEN(7,FILE='BOUNDALLAH.DAT',STATUS='OLD')
      READ(7,80) TESTO
      WRITE(*,80) TESTO
      READ(7,80) TESTO
      WRITE(*,80) TESTO
      READ(7,80) TESTO
      WRITE(*,80) TESTO
      READ(7,90) NT,NG
      WRITE(*,95) nt,ng
      READ(7,80) TESTO
      READ(7,80) TESTO
      WRITE(*,*)NT,NG      
      DO 60 J=1,NT
        DO 70 K=1,NG
       READ(7,107) T(J,K),G(J,K),(AAX(J,K,I),I=1,4)
  70   CONTINUE
  60  CONTINUE
       CLOSE(7)
  80  FORMAT(A80)
  90  FORMAT(I3,I2)
  95  FORMAT('numero T --->',I3,' numero gravita --->',I2)  
 107  FORMAT(2X,F8.1,3X,F4.1,2X,4E12.5)
C===================================================================
      RETURN
      STOP
      END
      SUBROUTINE NATMOS(LU,TE,MAS,X,Z,FRAZ,ALFA,RATM,PATM,TATM,MATM,NATM
     $,IPR,IDEBUG,IDER,ISING,MAXNE,nmd)
C***********************************************************************
C    INPUT:
C       LU=>LOG(L/LSUN) - TE=>LOG(TEFF) -  MAS=>MASSA/MSUN
C        X=>CHIMICA DELL'ATMOSFERA
C     FRAZ=>FRAZIONE DI MASSA INCLUSA NELL'ATMOSFERA
C     ALFA=>ALTEZZA DI SCALA DI PRESSSIONE L/HP
C     IPR=>SE=1 STAMPA DETTAGLI DELL'ATMOSFERA
C   OUTPUT:
C        P=>PRESSIONE ALLA BASE DELL'ATMOSFERA
C        T=>TEMPERATURA ALLA BASE DELL'ATMOSFERA
C        R=>RAGGIO ALLA BASE DELL'ATMOSFERA
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'consts.2p1'
      REAL *8 LU,TE,MAS,X,Z,FRAZ,ALFA,PTOT,PGAS,PRAD,RAGGIO,VSOUND
      REAL *8 MTOT,RTOT,LTOT,TEFF,FACT,FACR,DGRA,T,DP,SERV,RAPPORTO
      REAL *8 PASSOP,PASSOT,PASSOR,GRAD,GRADV,BORDOMASSA,PASSO,FACTS
      REAL *8 TAU,TAUV,GIINI,T4,GI,O1,O2,O3,VEELCO,atau(3000)
      REAL *8 RO,RADIAT,CAP,PATM,TATM,RATM,MATM,DTAUP,DTAU,AD,DEP
      REAL *8 RSER,PSER,TSER,MSER,CP,MU,MUT,DR,DM,DT,ERMAX,PGASV,PT
      REAL *8 HM,HR,HP,HT,HRO,HMU,HRAD,HAD,HGRA,HKAP,HGI,T1,Q,RIP,HBET
      REAL *8 C1(3),C2(3),C3(3),FATTO,grt,xhhhh,radiat1,G1,DEL,EMUE
      real *8 RAG,ahlim,ahelim,rlim,APP(26,3000),radr,bpp(26,3000)
      REAL *8 DDR,CON,CON1,CONF,confv,RTOT1,DDR1,emcoatm,difgra,difgrav
      REAL *8 Psave,Tsave,Rsave,Msave,deltaR,deltaP,deltaM,deltaT
      LOGICAL*1 CONTINUA,ULTIMO,PIPPO
      INTEGER *4 IPRESSSTART,IPRESSSTOP,IPRESSSWITCH,nmd1
      INTEGER *4 IPR,NUM,K,IDEBUG,STAMPA,NATM,HNU,IDER,ISING,MAXNE
      INTEGER *4 indatm,nmd,nmlimatm,nmlim,lll,iii,jjj,npa,icatm,ipippo
      integer *4 KK,II
      DIMENSION X(MAXNE),RATM(1),PATM(1),TATM(1),MATM(1)
      common/NGATM/xhhhh,nmlim,indatm,nmlimatm
      common /difatm/ emcoatm,icatm
      
      common/NOSUB/IPRESSSTART,IPRESSSTOP,IPRESSSWITCH
      common/astros/NPA,BPP

      DATA FACT/3.D-02/,FACR/4.D-02/,PASSOP/5.D-02/,DTAU/1.D-01/
      DATA FACTS/3.D-02/
      DATA ERMAX/1.D-10/
C***********************************************************************
      
      RAG  = STBSUN*DSQRT(10.D0**LU)/((10.D0**TE)**2)
      LLL=0
      IF(ISING.EQ.1.AND.IDER.EQ.0) THEN
         IF(TE.GT.3.8D0) THEN
            FACT  = 3.D-02
            FACTS = 2.D-02
            FACR  = 3.D-02
            PASSOP= 3.D-02
            DTAU  = 3.D-02
         ELSE
           FACT  = 3.D-03
            FACTS = 2.D-03
            FACR  = 3.D-03
            PASSOP= 3.D-03
            DTAU  = 2.D-02
         ENDIF
      ENDIF
      MTOT=MAS*SUNMg33                                          !M(1.E33)
      LTOT=(10.D0**LU)*SUNLe*1.D-32                             !L(1.E32)
      TEFF=10.D0**TE                                            !T(KELVI)
      RTOT=DSQRT(1.D12/(4.D0*PI*STBOLZ))*(DSQRT(LTOT))/(TEFF**2)

      RAGGIO=RTOT/(SUNRcm*1.D-10)                              !R/(R_SUN)
      RTOT1=RTOT*1.D10
      IF(IPR.EQ.1.OR.IDEBUG.EQ.1) THEN
         WRITE(2,100)MAS,LU,TE,RAGGIO,X(3),Z,ALFA,FRAZ
      ENDIF
c      if(nmlimatm.eq.0)then
c         grt=dlog10(GRAVC0*1.D+13*MTOT/(RTOT*RTOT))
c         if(grt.gt.4.d0)then
c            write(*,*)'grt=',grt
c            nmlimatm=nmd
c            write(*,*)'nmlimatm=',nmlimatm
c            indatm=0
c         endif
c      endif
      indatm=1
C***********************************************************************
      IF(indatm.GT.0.5D0) THEN
c***********************************************************************                   
C*************CALCOLO CONDIZIONE AL CONTORNO: P0 DI TAU0 ***************
        TAU   = 1.D-10
        TAUV  = TAU
        IF(IDER.EQ.0) THEN
           IF(TE.LT.4.6D0) THEN
              FATTO=0.1D0
           ELSE
              FATTO=0.03D0
           ENDIF
           PGAS  = FATTO*CN3*(0.75D0*TEFF)**4        ! CGS
           PGASV = PGAS
        ELSE
           PGAS  = PGASV
        ENDIF
       GIINI = GRAVC0*1.D+13*MTOT/(RTOT*RTOT)         ! CGS
        GI    = GIINI
        NUM   = 0


        DO WHILE (DABS((TAUV-TAU)/TAU).GT.1.D-08.AND.NUM.LT.200.OR.
     $   NUM.EQ.0)
           TAUV=TAU
           NUM=NUM+1
C***********************KRISHNA SWAMY 1966***************************
c           T4=0.75D0*(TEFF**4)*(TAU+1.390D0-0.815D0*DEXP(-2.54D0*TAU)-
c     $      0.025D0*DEXP(-30.D0*TAU))
C*********************** VERRAZZA *********************************
           T4=0.75D0*(TEFF**4)*(TAU+1.017D0-0.3D0*DEXP(-2.54D0*TAU)-
     $      0.291D0*DEXP(-30.D0*TAU))
c**********HOLWEGER-MUELLER
c
c       T4=TEFF**4*(0.2477d0+(1.106700d0/((0.844-dlog10(TAU))**1.52d0)))
c
C*************EDDINGTON******************************************   
C
c           T4=0.75D0*(TEFF**4)*(TAU+2.D0/3.D0)
C***********************************************************************
           T=T4**0.25D0
           PRAD=CN3*T4
           PTOT=PGAS+PRAD
           CALL EOS(1,1,MAXNE,1,PTOT,T,X,'ATMOS   ',RO,O1,O2,O3,G1,
     #DEL,EMUE,1)
           CALL NKAPPA(1,1,MAXNE,1,RO,T,X,CAP)
           RADIAT=(3.D0*1.D+12)/(64.D0*PI*STBOLZ)*CAP*LTOT*PTOT/
     $      (RTOT*RTOT*GIINI*T4)
           GI=GIINI*(1.D0-4.D0*RADIAT*PRAD/PTOT)
           TAU=PGAS*CAP/GI
           IF(IDEBUG.EQ.1) THEN
              IF(NUM.EQ.1) THEN
                 WRITE(2,101)
              ENDIF
              WRITE(2,102)NUM,TAU,PGAS,PRAD,T,GI,RADIAT,RO,CAP
           ENDIF
        END DO


        IF(NUM.GE.200) THEN
           WRITE(2,103)
           WRITE(*,103)
        ENDIF
        IF(TAU.LE.1.D-08)TAU=1.D-08
        IF(IPR.EQ.1.OR.IDEBUG.EQ.1) THEN
           SERV=PRAD/PTOT
           WRITE(2,104)TAU,PTOT,T,SERV,NUM
        ENDIF
C******* INTEGRAZIONE IN TAU FINO A TAU=2/3 O FINO A CONVETTIVO ********
        NUM=0
       grad=0.d0
        CONTINUA=.TRUE.
        DO WHILE (CONTINUA)
           NUM=NUM+1
           IF(NUM.GT.10000) THEN
              WRITE(*,*)' TROPPI PASSI IN TAU'
              WRITE(2,*)' TROPPI PASSI IN TAU'
              STOP
           ENDIF
           DTAUP=DTAU*TAU
           IF((TAU+DTAUP).GE.(2.D0/3.D0)) DTAUP=(2.D0/3.D0-TAU)+1.D-12
           TAU=TAU+DTAUP
           DEP=GI*DTAUP/CAP
           PGAS=PGAS+DEP
C***********************KRISHNA SWAMY 1966***************************
c           T4=0.75D0*(TEFF**4)*(TAU+1.390D0-0.815D0*DEXP(-2.54D0*TAU)-
c     $      0.025D0*DEXP(-30.D0*TAU))
C*********************** VERRAZZA *********************************
           T4=0.75D0*(TEFF**4)*(TAU+1.017D0-0.3D0*DEXP(-2.54D0*TAU)-
     $      0.291D0*DEXP(-30.D0*TAU))
c**********HOLWEGER-MUELLER
c
c       T4=TEFF**4*(0.2477d0+(1.106700d0/((0.844-dlog10(TAU))**1.52d0)))
c
C*************EDDINGTON******************************************   
C
c           T4=0.75D0*(TEFF**4)*(TAU+2.D0/3.D0)
C***********************************************************************
           T=T4**0.25D0
           PRAD=CN3*T4
           PTOT=PGAS+PRAD
           CALL EOS(1,1,MAXNE,1,PTOT,T,X,'ATMOS   ',RO,AD,CP,MU,G1,
     #DEL,EMUE,6)
           CALL NKAPPA(1,1,MAXNE,1,RO,T,X,CAP)
           RADIAT=(3.D0*1.D+12)/(64.D0*PI*STBOLZ)*CAP*LTOT*PTOT/  
     $      (RTOT*RTOT*GIINI*T4)
          RADIAT1=(3.D0*1.D-01)/(64.D0*PI*STBOLZ*GRAVC0)*CAP*LTOT*PTOT/
     $             (Mtot*T**4)
           GI=GIINI*(1.D0-4.D0*RADIAT*PRAD/PTOT)
           IF(RADIAT.GE.AD) then
            CONTINUA=.FALSE.
c             write(*,*)'il convettivo blocca int in atmosfera '
c            pause
          ENDIF
          IF(IDEBUG.EQ.1) THEN
              IF(NUM.EQ.1) THEN
                 WRITE(2,109)
              ENDIF
              WRITE(2,105)NUM,TAU,PGAS,PRAD,T,RO,CAP,GI,RADIAT,AD
           ENDIF
c      IF(RAG.GT.(rlim-1.d-2).AND.RAG.LT.(rlim+1.d-2))then
       IF(ipr.eq.2)then
               BPP(1,num)=rtot*1.d10
              BPP(2,num)=MTOT/(MAS*SUNMg33)
              BPP(3,num)=T
              BPP(4,num)=PGAS+PRAD
              BPP(5,num)=RO
              BPP(6,num)=x(1)
              BPP(7,num)=(10**LU)*SUNLE
              BPP(8,num)=CAP
              BPP(9,num)=0.D0
              BPP(10,num)=AD
              BPP(11,num)=1.d0-x(1)-X(2)-X(3)
              BPP(12,num)=0.D0
              BPP(13,num)=0.D0
              BPP(14,num)=0.D0
              BPP(15,num)=X(2)
              BPP(16,num)=X(3)
              BPP(17,num)=X(4)
              BPP(18,num)=X(10)
              BPP(19,num)=X(5)
              BPP(20,num)=X(6)
              BPP(21,num)=grad
              BPP(22,num)=radiat
              BPP(23,num)=G1
              BPP(24,num)=CP
              BPP(25,num)=DEL
              BPP(26,num)=EMUE
              atau(num)=tau
              
       ENDIF
     
               IF(IPR.EQ.1) THEN
                  STAMPA=NUM-(NUM/1)*1
               ELSE
                  STAMPA=1
               ENDIF
               IF(IPR.EQ.1.AND.(NUM.EQ.1.OR.ULTIMO.EQV..TRUE.)) STAMPA=0
               IF(IDEBUG.EQ.1) STAMPA=0
               
          IF(TAU.GE.(2.D0/3.D0)) CONTINUA=.FALSE.
        END DO
               npa=num
       do lll=npa,2,-1
         BPP(1,lll-1)=BPP(1,lll)+(atau(lll)-atau(lll-1))/
     #   (((BPP(8,lll)+BPP(8,lll-1))/2.d0)*
     #   (BPP(5,lll)+BPP(5,lll-1))/2.d0)
c         write(29,'(4e16.9)')BPP(23,lll),BPP(1,lll)
c         write(*,'(4e16.9)')BPP(23,lll),BPP(1,lll)
       enddo
       lll=0

        IF(IPR.EQ.1.OR.IDEBUG.EQ.1) THEN
           WRITE(2,106)TAU,T,PTOT,PGAS,PRAD,CAP,NUM
        ENDIF
c***********************************************************************
c***************MODIFICA PER LETTURA BOUNDALLAH**************************        
      ELSE
       GIINI=GRAVC0*1.D+13*MTOT/(RTOT*RTOT)
        CALL CONDCON(LU,TE,MAS,Z,PT,T)
        T4=T**4
        PRAD=CN3*T4
        PGAS=PT
        PTOT=PGAS+PRAD
        CALL EOS(1,1,MAXNE,1,PTOT,T,X,'ATMOS   ',RO,AD,O2,O3,G1,
     #DEL,EMUE,3)
        CALL NKAPPA(1,1,MAXNE,1,RO,T,X,CAP)
        RADIAT=(3.D0*1.D+12)/(64.D0*PI*STBOLZ)*CAP*LTOT*PTOT/
     $      (RTOT*RTOT*GIINI*T4)
        GI=GIINI*(1.D0-4.D0*RADIAT*PRAD/PTOT)
      ENDIF
c===== modifica per eliminare subatmosfera
      Rsave = RTOT
      Psave = PTOT            ! CGS
      Tsave = T               ! CGS
      Msave = MTOT
      ipippo=0
c      if(nmd.gt.40)then
c       IPRESSSWITCH=1
c       goto 666
c      endif
      difgra=radiat-ad     !inizializza diff gradienti come vengono da tau=2/3
      icatm=0
      LLL=0  !SERVE PER ASTROSISMOLOGIA SE ABBIAMO LA SUBATMOSFERA
      if (nmd.gt.IPRESSSTOP.and.IPRESSSWITCH.eq.0) goto 666
c===== FINE modifica per eliminare subatmosfera
C*************** INTEGRAZIONE IN PRES FINO A MASSA=1-FRAZ **************
      RATM(1) = RTOT
      PATM(1) = PTOT            ! CGS
      TATM(1) = T               ! CGS
      MATM(1) = MTOT
      BORDOMASSA=MTOT*(1.D0-FRAZ)
      GRADV=RADIAT
      NUM=1
      HNU=0
      CONTINUA=.TRUE.
      ULTIMO=.FALSE.
c      PASSO=PASSOP/10.D0
      PASSO=PASSOP/5000.D0
      DP=PASSO*PTOT
      IF(IPR.EQ.1.OR.IDEBUG.EQ.1) WRITE(2,107)
      DO WHILE (CONTINUA)
         IF(NUM.GT.1995) THEN
            WRITE(*,*)' TROPPI PASSI IN ATMOSFERA'
            WRITE(2,*)' TROPPI PASSI IN ATMOSFERA'
            STOP
         ENDIF
         RSER=RATM(NUM)
         PSER=PATM(NUM)
         TSER=TATM(NUM)
         MSER=MATM(NUM)
         DO K=1,3
            VEELCO=0.D0
            CALL EOS(1,1,MAXNE,1,PSER,TSER,X,'ATMOS   ',RO,AD,CP,MU,G1,
     #DEL,EMUE,3)
            CALL NKAPPA(1,1,MAXNE,1,RO,TSER,X,CAP)
            RADIAT=(3.D0*1.D-01)/(64.D0*PI*STBOLZ*GRAVC0)*CAP*LTOT*PSER/
     $             (MSER*TSER**4)
           IF(RADIAT.GE.AD) THEN
              GI=GRAVC0*1.D+13*MSER/(RSER*RSER)
               Q=MU
               IF(Q.LT.0.D0) Q=1.D-20
              CALL NSUPERA(RO,PSER,TSER,CAP,Q,CP,AD,GI,ALFA,GRAD,VEELCO,
     $                     RADIAT)
            ELSE
               GRAD=RADIAT
            ENDIF
           C1(K)=-(1.D-23/GRAVC0)*RSER*RSER/(MSER*RO)*DP
           C2(K)=-(4.D0*PI*1.D-26/GRAVC0)*(RSER**4)/MSER*DP
            C3(K)=TSER/PSER*GRAD*DP
            IF(K.EQ.1) THEN
               RIP=0.5D0
            ELSE
               RIP=1.0D0
            ENDIF
            PSER=PATM(NUM)+DP*RIP
            RSER=RATM(NUM)+C1(K)*RIP
            MSER=MATM(NUM)+C2(K)*RIP
            TSER=TATM(NUM)+C3(K)*RIP
            IF(K.EQ.1.AND.NUM.NE.HNU) THEN
               PRAD=CN3*TATM(NUM)**4
               RAPPORTO=PRAD/PATM(NUM)
               VSOUND=DSQRT(PATM(NUM)/RO)
              HM  = MATM(NUM)/(MAS*SUNMg33)
              HR  = RATM(NUM)/(RAGGIO*SUNRcm*1.d-10)
               HP  = DLOG10(PATM(NUM))
               HT  = DLOG10(TATM(NUM))
               HRO = DLOG10(RO)
               HGI = DLOG10(GI)
               HKAP= DLOG10(CAP)
               HMU = MU
C               HBET= RAPPORTO
               HBET= DLOG10((VEELCO+1.D-10)/VSOUND)
               HRAD= RADIAT
               HAD = AD
               HGRA= GRAD
               HNU = NUM
c************ Questo blocco serve per la diffusione *******************
               difgrav=difgra
               difgra=radiat-ad
               if(difgrav.ge.0.d0.and.difgra.lt.0.d0) then
                 emcoatm=matm(num)/1.989d0
                 icatm=1
               endif
C***********************************************************************
              
              
C***********************************************************************
               IF(IPR.EQ.1) THEN
                  STAMPA=NUM-(NUM/1)*1
               ELSE
                  STAMPA=1
               ENDIF
               IF(IPR.EQ.1.AND.(NUM.EQ.1.OR.ULTIMO.EQV..TRUE.)) STAMPA=0
               IF(IDEBUG.EQ.1) STAMPA=0
               
              IF(STAMPA.EQ.0) THEN
                  WRITE(2,108)HM,HR,HP,HT,HRO,HMU,HKAP,HRAD,HAD,HGRA,
     $                  HBET,HGI,HNU
               ENDIF
              IF(IPR.EQ.2)THEN
               LLL=LLL+1
               APP(1,LLL)=RATM(NUM)*1.d10
               APP(2,LLL)=MATM(NUM)/(MAS*SUNMg33)
               APP(3,LLL)=TATM(NUM)
               APP(4,LLL)=PATM(NUM)
               APP(5,LLL)=RO
               APP(6,LLL)=x(1)
               APP(7,LLL)=(10**LU)*SUNLE
               APP(8,LLL)=CAP
               APP(9,LLL)=0.d0
               APP(10,LLL)=AD
               APP(11,LLL)=1.d0-x(1)-X(2)-X(3)
               APP(12,LLL)=RAG*SUNRCM-RATM(NUM)*1.d10
               APP(13,LLL)=0.D0
               APP(14,LLL)=0.D0
               APP(15,LLL)=X(2)
               APP(16,LLL)=X(3)
               APP(17,LLL)=X(4)
               APP(18,LLL)=X(10)
               APP(19,LLL)=X(5)
               APP(20,LLL)=X(6)
               APP(21,LLL)=grad
               APP(22,LLL)=radiat
               APP(23,LLL)=G1
               APP(24,LLL)=CP
               APP(25,LLL)=DEL
               APP(26,LLL)=EMUE
              ENDIF
C********************** DETERMINA NUOVO PASSO IN PRESSIONE *************
            ENDIF
         END DO
         IF(ULTIMO.EQV..FALSE.) THEN
            DR=(C1(1)+2.*C1(2)+C1(3))/4.D0
            DM=(C2(1)+2.*C2(2)+C2(3))/4.D0
            DT=(C3(1)+2.*C3(2)+C3(3))/4.D0
            MSER=((MATM(NUM)+DM)-BORDOMASSA)/BORDOMASSA
            IF(MSER.GT.(-1.D0*ERMAX)) THEN
               RATM(NUM+1)=RATM(NUM)+DR
               PATM(NUM+1)=PATM(NUM)+DP
               TATM(NUM+1)=TATM(NUM)+DT
               MATM(NUM+1)=MATM(NUM)+DM
               NUM=NUM+1
               DGRA=DABS(GRADV-GRAD)
               IF(DGRA.GT.0.D0) THEN
                  IF(GRAD.GT.0.40) THEN
                     PASSOT=FACTS*PASSO/DGRA
                  ELSE
                     PASSOT=FACT*PASSO/DGRA
                  ENDIF
               ELSE
                  PASSOT=PASSOP
               ENDIF
               PASSOR=DABS(FACR*PASSO*RATM(NUM-1)/DR)
               PASSO=DMIN1(PASSOP,PASSOT,PASSOR)
c                  passo=passo/10.d0
               IF(PASSO.LT.1.D-03)PASSO=1.D-03
               GRADV=GRAD
               DP=PASSO*PATM(NUM)
               IF(MSER.LT.ERMAX) ULTIMO=.TRUE.
               NATM=NUM
C***********************************************************************
            ELSE
               PASSO=PASSO*DABS((MATM(NUM)-BORDOMASSA)/DM)
               DP=PASSO*PATM(NUM)
            ENDIF
         ELSE
            CONTINUA=.FALSE.
         ENDIF
      END DO
c      if(ipr.eq.2)
c        WRITE(39,'(a9,2x,i4)')'ATMOSFERA',lll
c       DO JJJ=LLL,1,-1
c        write(39,'(1P,25D17.9)')(APP(III,JJJ),III=1,22)
c       ENDDO
c      endif
c       WRITE(39,'(a9,2x,i4)')'photosfe',npa
c       DO JJJ=npa,1,-1
cc       BPP(12,jjj)=rtot*1.d10-BPP(1,JJJ)
c        write(39,'(1P,25D17.9)')(BPP(III,JJJ),III=1,22)
c       ENDDO
c===== modifica per eliminare subatmosfera
666   if (nmd.ge.IPRESSSTART.and.nmd.le.IPRESSSTOP.and.
     $    IPRESSSWITCH.eq.0.and.ipippo.lt.1) then
        deltaR=dlog10(RATM(NATM))-dlog10(Rsave)
        deltaP=dlog10(PATM(NATM))-dlog10(Psave)
        deltaT=dlog10(TATM(NATM))-dlog10(Tsave)
        deltaM=dlog10(MATM(NATM))-dlog10(Msave)
       nmd1=nmd
c       if(nmd.gt.10)nmd1=10
c        RATM(1) = 10**(dlog10(Rsave)+deltaR/
c     $         DFLOAT(IPRESSSTOP-IPRESSSTART+1)*
c     $         DFLOAT(IPRESSSTOP+1-nmd1))
c        PATM(1) = 10**(dlog10(Psave)+deltaP/
c     $         DFLOAT(IPRESSSTOP-IPRESSSTART+1)*
c     $         DFLOAT(IPRESSSTOP+1-nmd1))
c           
c        TATM(1) = 10**(dlog10(Tsave)+deltaT/
c     $         DFLOAT(IPRESSSTOP-IPRESSSTART+1)*
c     $         DFLOAT(IPRESSSTOP+1-nmd1))
c        
c        MATM(1) = 10**(dlog10(Msave)+deltaM/
c     $         DFLOAT(IPRESSSTOP-IPRESSSTART+1)*
c     $         DFLOAT(IPRESSSTOP+1-nmd1))
        RATM(1) = 10**(dlog10(Rsave)+deltaR/
     $         49.d0*
     $         DFLOAT(49-nmd1))
        PATM(1) = 10**(dlog10(Psave)+deltaP/
     $         49.d0*
     $         DFLOAT(49-nmd1))           
        TATM(1) = 10**(dlog10(Tsave)+deltaT/
     $         49.d0*
     $         DFLOAT(49-nmd1))                   
        MATM(1) = 10**(dlog10(Msave)+deltaM/
     $         49.d0*
     $         DFLOAT(49-nmd1))                   
        NATM=1
      elseif (nmd.gt.IPRESSSTART.and.IPRESSSWITCH.eq.0) then
        RATM(1) = Rsave
        PATM(1) = Psave           ! CGS
        TATM(1) = Tsave               ! CGS
        MATM(1) = Msave
        NATM=1
      end if
c===== FINE modifica per eliminare subatmosfera
      IF(ipr.eq.2)then
        IF(LLL.GT.0)THEN
         DO II=2,LLL
           DO KK=1,26
             BPP(KK,II+NPA-1)=APP(KK,II)
           ENDDO
         ENDDO
         NPA=NPA+LLL-1
       ENDIF
c       WRITE(40,'(a9,2x,i4)')'photosfe',npa
       DO JJJ=npa,1,-1
       BPP(12,jjj)=rtot*1.d10-BPP(1,JJJ)
c        write(40,'(1P,25D17.9)')(BPP(III,JJJ),III=1,22)
       ENDDO
      ENDIF

      RETURN
C***********************************************************************
  100 FORMAT(/,1X,'M=',F7.4,' logL=',F7.4,' logTe=',F6.4,' R/Ro=',
     $1P,E9.2,' Y=',0P,F6.3,' Z=',1P,E10.3,0P,' alfa_ml=',F7.4,
     $' Fraz=',1P,E10.3,/)
  101 FORMAT(' ITER  TAU       PGAS      PRAD      TEMP      GRAVITY   G
     $RA-RAD   RO        KAP')
  102 FORMAT(1X,I3,1P,8E10.3)
  103 FORMAT(1X,'******** CALCOLO P0 DI TAU0 NON CONVERGE ********')
  104 FORMAT(' P(TAU0=',1P,E10.3,')=',E10.3,' T=',E10.3,' PRAD/PTOT=',
     $E10.3,' N.ITER=',I3,/)
  105 FORMAT(I5,1P,9E9.2)
  106 FORMAT(1X,'FINE ATMOSFERA TAU=',1P,E10.3,1X,'T=',E10.3,1X,'PTOT=',
     $E10.3,1X,'PGAS=',E10.3,1X,'PRAD=',E10.3,1X,'OPAC=',E10.3,1X,
     $'NITER=',I3,/)
  107 FORMAT(5X,'MASSA',5X,'RAGGIO',6X,'LOG P',3X'LOG T',2X,'LOG RO',
     $4X,'MU',5X,'LOG K',3X,'GRA-RAD',5X'GRA-AD',4X,'GRA-EFF',
     $1X,'LOG(Vc/Vsound)',1X,'LOG G',1X,'NPAS')
  108 FORMAT(1P,2E12.5,0P,5F8.4,1P,3E11.3,0P,F10.3,F10.3,I4)
  109 FORMAT(/,1X,'  N   TAU      PGAS     PRAD     TEMP     RO       KA
     $P      GRAVITY  GRA-RAD  GRA-AD')
      END
