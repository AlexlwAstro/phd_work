C*****************************************************************************
C  QUESTA ROUTINE CONTIENE: LEVMES - NEWMES - OPTIM  - PASTEM - 
C                           QUATM  
C*****************************************************************************
C
      SUBROUTINE LEVMES(N) 
C 
C    QUESTA SUBROUTINE SERVE PER TOGLIERE UN MESH 
C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      DIMENSION PX(MAXNE) 
      DO J=1,MAXNE-1 
         PX(J)=XXX(J,N-1)*(G(5,N)-G(5,N-1))+XXX(J,N)*(G(5,N+1)-G(5,N)) 
      END DO 
      DM=G(5,N+1)-G(5,N-1) 
      DO J=1,MAXNE-1 
         XXX(J,N-1)=PX(J)/DM 
      END DO 
      KK=N+1 
      DO K=KK,MAXME 
         DO J=1,7 
            G(J,K-1)=G(J,K) 
         END DO 
         DO J=1,MAXNE-1 
            XXX(J,K-1)=XXX(J,K) 
         END DO 
         CARXY(1,K-1)=CARXY(1,K) 
         CARXY(2,K-1)=CARXY(2,K) 
      END DO 
      MAXME=MAXME-1 
C     WRITE(*,*)N,G(5,N) 
      RETURN 
      END 

      SUBROUTINE NEWMES(N) 
C 
C    QUESTA SUBROUTINE AGGIUNGE UN MESH 
C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      DIMENSION PX(MAXNE),PG(7) 
      DO J=1,7 
         PG(J)=(G(J,N-1)+G(J,N))/2.D0 
      END DO 
c      PG(1)=DEXP((DLOG(G(1,N-1))+DLOG(G(1,N)))/2.D0)
      DO J=1,MAXNE-1 
c         PX(J)=XXX(J,N-1) 
         PX(J)=(XXX(J,N-1)+XXX(J,N))/2.d0 
      END DO 
      XID=CARXY(1,N-1) 
      XEL=CARXY(2,N-1) 
      DO K=N,MAXME 
         L=MAXME-K+N 
         DO I=1,7 
            G(I,L+1)=G(I,L) 
         END DO 
         DO J=1,MAXNE-1 
            XXX(J,L+1)=XXX(J,L) 
         END DO 
         CARXY(1,L+1)=CARXY(1,L) 
         CARXY(2,L+1)=CARXY(2,L) 
      END DO 
      DO I=1,7 
         G(I,N) = PG(I) 
      END DO 
      DO I=1,MAXNE-1 
         XXX(I,N) = PX(I) 
      END DO 
      CARXY(1,N)=XID 
      CARXY(2,N)=XEL 
      MAXME=MAXME+1 
c      WRITE(*,100)N,G(1,N),G(2,N),G(3,N),G(4,N),G(5,N),G(6,N),G(7,N) 
c      WRITE(*,100)N,XXX(1,N),XXX(2,N),XXX(3,N)
  100 FORMAT(I4,1P,7D10.3) 
      RETURN 
      END
      SUBROUTINE OPTIM
C
C    QUESTA SUBROUTINE SERVE PER ALTERARE IL NUMERO DEI MESH
C
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL*1 SALTA
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      DIMENSION VMM(5),URA(4),ULA(4),UPA(4),UTA(4),UMA(4)
      DIMENSION UGRA(4)
C DEFINIZIONE DIVERSI INTERVALLI DI MASSA
C     DATA VMM/ 0.0D+00 , 1.400D-01 , 1.000D0 , 1.000D0 , 1.000D0 /
      DATA VMM/ 0.0D+00 , 0.995D 00 , 1.000D0 , 1.000D0 , 1.000D0 /
c      DATA VMM/ 0.0D+00 , 1.000D00 , 1.000D0 , 1.000D0 , 1.000D0 /
C INTERVALLO:      #1        #2       #3       #4
C PRE FITTARE IL SOLE .10 .01 .10 .03 .01
C PRE PICCOLE MASSE   .20 .02 .10 .05 .01
      DATA URA /  0.30D0 ,  0.10D0 ,  0.30D0,  0.10D0/
      DATA ULA /  0.02D0 ,  0.006D0 , 0.02D0,  0.03D0/
      DATA UPA /  0.15D0 ,  0.050 ,   0.15D0,  0.15D0/
      DATA UTA /  0.05D0 ,  0.0150 ,  0.02D0,  0.10D0/
      DATA UMA /  0.01D0 ,  0.003D0 , 0.01D0,  0.01D0/
      DATA UGRA/  0.20D0 ,  0.06D0 ,  0.05D0,  0.05D0/
      DATA IS/0/
      EMSO=EMTOT/SUNMg33
      YM1=0.D0
      YM2=0.D0
      SCALO=0.D0
      if(is.eq.0)then
       nmd0=nmd
       is=1
      endif
c===introdotto per presequenza=============
      IF(XXX(7,1).LT.1.D-15)THEN
c===introdotto per presequenza=============
       IF(BB(1).GT.1.D-03.AND.XXX(1,1).GT.0.D0.AND.EMSO.GE.1.D0) THEN
         YM1=GICO/EMTOT-.02D0
         IF(YM1.LT.0.D0)YM1=0.D0
         YM2=YM1+0.02D0
       ENDIF
      ENDIF
csanti
      IF(BB(1).GT.1.D-02.AND.XXX(1,1).LE.0.D0.AND.XXX(3,1).GT.0.D0) THEN
        YM1=GICO/EMTOT-.01D0
        IF(YM1.LT.0.D0)YM1=0.D0
       DYM2=(20.D0-(EMTOT/SUNMg33))*.005D0
        IF( DYM2.LE.0.D0 )  DYM2 = 0.D0
        IF(EMSO.GT.2.D0) THEN
          YM2=YM1+DYM2
        ELSE
          YM2=YM1+0.05D0+DYM2
        ENDIF
      ENDIF
csanti
C********************** INIZIA PARTE INFITTIMENTO **********************
C      IF(G(4,1).LT.100.D0) THEN
C         N=8
C      ELSE
C         N=3
C      ENDIF
CSANTI
      N=2
c      write(*,'(3i8)')nmd0,nmd-nmd0,nmd
      if(nmd-nmd0.ge.1)n=4
      if(nmd-nmd0.ge.2)n= 5
      if(nmd-nmd0.ge.3)n=  8
C migliore configurazione per calcolo 1Mo asterosismologica
c      N=2
c      if(nmd.ge.4)n=4
c      if(nmd.ge.5)n=5
c      if(nmd.ge.6)n= 8
CSANTI
    1 N=N+1
      IF(MAXME.GT.(LIM-5)) GO TO 10
      GIG=G(5,N)/EMTOT
      DO I=1,4
         IF(GIG.GE.VMM(I).AND.GIG.LT.VMM(I+1)) GO TO 3
      END DO
      IF(I.GT.4)I=4
    3 CONTINUE
      UL=ULA(I)
      UP=UPA(I)
      UT=UTA(I)
      UR=URA(I)
      UM=UMA(I)
      UGR=UGRA(I)
      IF(XXX(1,1).GT.0.D0.AND.NSHUT.GE.1)THEN
        UL=ULA(1)/2.D0
        UP=UPA(1)/2.D0
        UT=UTA(1)/2.D0
        UR=URA(1)/2.D0
        UM=UMA(1)/2.D0
        UGR=UGRA(1)/2.D0
      ENDIF
      IF(ISTART.EQ.5)THEN
        UL=ULA(I)*6.D0
        UP=UPA(I)*6.D0
        UT=UTA(I)*6.D0
        UR=URA(I)*6.D0
        UM=UMA(I)*6.D0
        UGR=UGRA(I)*6.D0
      ENDIF
c      IF(ISTART.EQ.3.and.eta.gt.0.d0)THEN
c        UL=ULA(I)/6.D0
c        UP=UPA(I)/6.D0
c        UT=UTA(I)/6.D0
c        UR=URA(I)/6.D0
c        UM=UMA(I)/6.D0
c        UGR=UGRA(I)/6.D0
c      ENDIF

CSANTI
C      IF((G(2,N)/G(2,MAXME)).GT.1.5D0.AND.FL3A.GT.1.5D0) THEN
C        DENOM=G(2,N-1)*8.D0
C      ELSE
C        DENOM=G(2,MAXME)
C      ENDIF
CSANTI
      DENOM=G(2,MAXME)
CSANTI      
      P1=DABS(G(1,N)-G(1,N-1))/G(1,N-1)
      P2=DABS(G(2,N)-G(2,N-1))/DENOM
      P3=DABS(G(3,N)-G(3,N-1))/G(3,N-1)
      P4=DABS(G(4,N)-G(4,N-1))/G(4,N-1)
      P5=DABS(G(5,N)-G(5,N-1))/EMTOT
      IF(FLCAR.GT.1.D-1.AND.XXX(1,N).LE.0.D0.AND.XXX(3,N).LE.0.D0)
     #   P2=0.D0
      IF(EMSO.GT.30.D0.AND.G(6,N).GT.0.D0.AND.G(4,N).LT.10.D0) THEN
         P6=DABS((G(7,N)-G(7,N-1))/G(7,N-1))
      ELSE
         P6=0.D0
      ENDIF
      IF(N.EQ.MAXME) P6=0.D0
CSANTI
C      IF(BB(1).GT.1.D-3.AND.EMSO.GT.1.5D0.AND.GIG.GT.YM1.AND.GIG.LT.YM2.
C     $   AND.XXX(1,1).GT.0.D0) THEN
C         UR=URA(I)/5.D0
C         UL=ULA(I)/5.D0
C         UP=UPA(I)/3.D0
C         UT=UTA(I)/5.D0
C         UM=UMA(I)/10.D0
C      ENDIF
CSANTI

      IF(BB(1).GT.0.D0.AND.GIG.GT.YM1.AND.GIG.LT.YM2) THEN
c      write(*,*)'bb1----1'
       IF(XXX(1,1).GT.0.D0) THEN
         UR=URA(I)/5.D0
         UL=ULA(I)/5.D0
         UP=UPA(I)/3.D0
         UT=UTA(I)/5.D0
         UM=UMA(I)/10.D0
       ELSE
         UR=URA(I)/10.D0
         UL=ULA(I)/10.D0
         UP=UPA(I)/5.D0
         UT=UTA(I)/10.D0
         UM=UMA(I)/10.D0
       ENDIF
      ENDIF 
CSANTI

CSANTI
      IF(XXX(1,N).LE.0.D0.AND.XXX(3,N-1).GT.1.D-8.AND.BB(1).GT.0.03D0
     &   .AND.DABS(G(6,N)).LT.0.04) THEN
c       write(*,*)'bb1----2'
Cc         UP=UPA(I)/5.D0  ! originale
Cc         UT=UTA(I)/5.D0  ! originale
Cc         UM=UMA(I)/20.D0 ! originale
CC************* CONTROLLO SU ABBONDANZA HE CENTRALE (Santi) **************
       IF(XXX(3,1).LT.0.15D0) THEN
         UP=UPA(I)/10.D0
         UT=UTA(I)/10.D0
         UM=UMA(I)/30.D0
       ELSE
         UP=UPA(I)/5.D0
         UT=UTA(I)/5.D0
         UM=UMA(I)/20.D0
       ENDIF  
C************************************************************************
      ENDIF
CSANTI

C SANTI X DIFF
      IF(INDIFF.EQ.1.AND.GIG.LT.0.99D0.and.xxx(1,1).gt.0.d0) THEN
       SCALO=0.5D0
      ELSE
       SCALO=1.0D0
      ENDIF     

      UL=UL*SCALO
      UP=UP*SCALO
      UT=UT*SCALO
      UR=UR*SCALO
      UM=UM*SCALO
      UGR=UGR*SCALO

C SANTI X DIFF
      pgr=0.d0
      if(G(6,N).ge.0.d0.and.G(6,N-1).lt.0.d0 )then
       pgr=dabs(G(6,N-1)-G(6,N))
c       write(*,*)'env ',G(6,N-1),G(6,N),n-1,n,pgr
c       if(pgr.lt.1.d-3)then
c        write(*,*)'scatta per env'
c        pause
c       endif
      endif

      SALTA=.FALSE.
      IF(FLCAR.LE.0.D0) THEN
         IF(XXX(3,N).LE.0.D0.AND.HT1.LT.1.D+00) SALTA=.TRUE.
      ELSE
         IF(XXX(3,N).GT.0.D0.AND.HT1.LT.1.D+00) SALTA=.TRUE.
      ENDIF
      IF(SALTA.EQV..TRUE.) GO TO 6
       IF(P1.LT.UR.AND.P2.LT.UL.AND.P3.LT.UP.AND.P4.LT.UT.AND.P5.LT.UM.
     $ AND.P6.LT.UGR.and.pgr.lt.1.d-4) GO TO 6
      CALL NEWMES(N)
      n=n+1
    6 CONTINUE
      IF(N.LT.MAXME) GO TO 1
C****************** FINE FASE INFITTIMENTO MESHPOINTS ******************
   10 CONTINUE
C****************** INIZIA FASE SFOLTIMENTO MESHPOINTS *****************
      N= 1
   11 N=N+1
      GIG=G(5,N)/EMTOT
      DO I=1,4
         IF(GIG.GE.VMM(I).AND.GIG.LT.VMM(I+1)) GO TO 13
      END DO
      IF(I.GT.4)I=4
   13 CONTINUE
c      UL=ULA(I)*0.8D0
c      UP=UPA(I)*0.8D0
c      UT=UTA(I)*0.8D0
c      UR=URA(I)*0.8D0
c      UM=UMA(I)*0.8D0
c      UGR=UGRA(I)*0.8D0
CSANTI      
      IF( XXX(1,1).LE.0.D0 .AND. XXX(3,1).GT.1.D-08 .AND.
     #DABS(G(6,N)).LE.0.04D0) GO TO 16
      IF( XXX(1,N-1).LE. 0.D0 . AND .XXX(1,N)  .GT.0.D0 ) GO TO 16
      IF( XXX(1,N)  .LE. 0.D0 . AND .XXX(1,N+1).GT.0.D0 ) GO TO 16
      GMX=G(5,N)/EMTOT
      IF(BB(1).GT.0.D0.AND.GIG.GE.YM1.AND.GIG.LE.YM2.AND.
     #XXX(1,1).LE.0.D0)GO TO 16

      IF(BB(1).GT.0.D0.AND.GIG.GT.YM1.AND.GIG.LT.YM2.
     $   AND.XXX(1,1).GT.0.D0) GO TO 16
CSANTI          
C      IF(BB(1).GT.1.D-3.AND.EMSO.GT.1.5D0.AND.GIG.GT.YM1.AND.GIG.LT.YM2.
C     $   AND.XXX(1,1).GT.0.D0) THEN
C         UR=URA(I)/5.D0*0.8D0
C         UL=ULA(I)/5.D0*0.8D0
C         UP=UPA(I)/3.D0*0.8D0
C         UT=UTA(I)/5.D0*0.8D0
C         UM=UMA(I)/10.D0*0.8D0
C      ENDIF
      IF(XXX(1,N).LE.0.D0.AND.XXX(3,N-1).GT.1.D-8.AND.BB(1).GT.0.03D0
     $   .AND.DABS(G(6,N)).LT.0.04D0) THEN
c       write(*,*)'bb1----3'
       IF(XXX(3,1).LT.0.15D0) THEN
         UP=UPA(I)/10.D0*0.8D0
         UT=UTA(I)/10.D0*0.8D0
         UM=UMA(I)/30.D0*0.8D0
       ELSE
         UP=UPA(I)/5.D0*0.8D0
         UT=UTA(I)/5.D0*0.8D0
         UM=UMA(I)/20.D0*0.8D0
       ENDIF
      ENDIF
C      UH=0.05D0
C      UHE=0.05D0
C      WH=0.D0
C      WHE=0.D0
CSANTI
      IF(BB(1).GT.0.D0.AND.GIG.GT.YM1.AND.GIG.LT.YM2) THEN
c       write(*,*)'bb1----4'
       IF(XXX(1,1).GT.0.D0) THEN
         UR=(URA(I)/5.D0)*0.8D0
         UL=(ULA(I)/5.D0)*0.8D0
         UP=(UPA(I)/3.D0)*0.8D0
         UT=(UTA(I)/5.D0)*0.8D0
         UM=(UMA(I)/10.D0)*0.8D0
       ELSE
         UR=(URA(I)/10.D0)*0.8D0
         UL=(ULA(I)/10.D0)*0.8D0
         UP=(UPA(I)/5.D0)*0.8D0
         UT=(UTA(I)/10.D0)*0.8D0
         UM=(UMA(I)/10.D0)*0.8D0
       ENDIF
      ENDIF 




CSANTI
      UH=0.001D0
      UHE=0.001D0
      WH=0.D0
      WHE=0.D0
      
csanti      
c      IF(XXX(1,N-1).GT.1.D-6) WH=DABS((XXX(1,N-1)-XXX(1,N+1))/XXX(1,N))
c      IF(XXX(3,N-1).GT.1.D-6) WHE=DABS((XXX(3,N-1)-XXX(3,N+1))/XXX(3,N))
c      IF(BB(1).GT.1.D-3.AND.XXX(1,N).LE.0.D0.AND.XXX(3,N).GT.1.D-08.AND.
c     #   DABS(G(6,N)).LE.2.D-2.AND.WHE.GT.UHE) GO TO 16
c      IF(N.LT.(MAXME-6)) THEN
c         IF(XXX(1,N).LE.0.D0.AND.XXX(1,N+5).GT.0.D0) GO TO 16
c         IF(XXX(1,N).GT.0.D0.AND.XXX(1,N-5).LE.0.D0) GO TO 16
c         IF(XXX(3,N).LE.0.D0.AND.XXX(3,N+1).GT.0.D0) GO TO 16
c         IF(XXX(3,N).GT.0.D0.AND.XXX(3,N-1).LE.0.D0) GO TO 16
c      ENDIF
c      IF((G(2,N)/G(2,MAXME)).GT.1.5D0.AND.FL3A.GT.1.5D0) THEN
c         DENOM=G(2,N)*8.D0
c      ELSE
c         DENOM=G(2,MAXME)
c      ENDIF
c      P1=DABS(G(1,N-1)-G(1,N+1))/G(1,N)
c      P2=DABS(G(2,N-1)-G(2,N+1))/DENOM
c      P3=DABS(G(3,N-1)-G(3,N+1))/G(3,N)
c      P4=DABS(G(4,N-1)-G(4,N+1))/G(4,N)
c      P5=DABS(G(5,N-1)-G(5,N+1))/EMTOT
c      IF(EMSO.GT.30.D0.AND.G(6,N).GT.0.D0.AND.G(4,N).LT.10.D0) THEN
c         P6=DABS((G(7,N-1)-G(7,N+1))/G(7,N))
c      ELSE
c         P6=0.D0
c      ENDIF
c      SALTA=.FALSE.
c      IF(FLCAR.LE.0.D0) THEN
c         IF(XXX(3,N).LE.0.D0.AND.HT1.LT.1.D+00) SALTA=.TRUE.
c      ELSE
c         IF(XXX(3,N).GT.0.D0.AND.HT1.LT.1.D+00) SALTA=.TRUE.
c      ENDIF
c      IF(SALTA.EQV..TRUE.) GO TO 16
c      IF(P1.GT.UR.OR.P2.GT.UL.OR.P3.GT.UP.OR.P4.GT.UT.OR.P5.GT.UM.OR.
c     $ P6.GT.UGR) GO TO 16
c      IF(WH.GT.UH.OR.WHE.GT.UHE) GO TO 16
c      
csanti      
C SANTI X DIFF
      IF(INDIFF.EQ.1.AND.GIG.LT.0.99D0.and.xxx(1,1).gt.0.d0) THEN
       SCALO=0.5D0
      ELSE
       SCALO=1.0D0
      ENDIF     
      SCALO=1.0D0
      UL=UL*SCALO
      UP=UP*SCALO
      UT=UT*SCALO
      UR=UR*SCALO
      UM=UM*SCALO
      UGR=UGR*SCALO

C SANTI X DIFF
      
csanti 
      IF(SALTA.EQV..TRUE.) GO TO 6
      UH=0.0001D0
      UHE=0.0001D0
      WH=0.D0
      WHE=0.D0
      WH=DABS(XXX(1,N-1)-XXX(1,N))/XXX(1,N)
      WHE=DABS(XXX(3,N-1)-XXX(3,N))/XXX(3,N)

      IF(WH.GT.UH.OR.WHE.GT.UHE) GO TO 16
c      pgr=1.d0
c      if(G(6,N).ge.0.d0.and.G(6,N-1).lt.0.d0 )then
c       pgr=dabs(G(6,N-1)-G(6,N))
c       write(*,*)'l env ',G(6,N-1),G(6,N),n-1,n,pgr
c       if(prg.gt.1.d-4)goto 16
c      elseif(G(6,N).lt.0.d0.and.G(6,N-1).ge.0.d0 )then
c       pgr=dabs(G(6,N-1)-G(6,N))
c       write(*,*)'l core ',G(6,N-1),G(6,N),n-1,n,pgr
c       if(prg.gt.1.d-4)goto 16
c      endif
      
      DENOM=G(2,MAXME)
      P1=DABS(G(1,N-1)-G(1,N+1))/G(1,N)
      P2=DABS(G(2,N-1)-G(2,N+1))/DENOM
      P3=DABS(G(3,N-1)-G(3,N+1))/G(3,N)
      P4=DABS(G(4,N-1)-G(4,N+1))/G(4,N)
      P5=DABS(G(5,N-1)-G(5,N+1))/EMTOT

      IF(EMSO.GT.30.D0.AND.G(6,N).GT.0.D0.AND.G(4,N).LT.10.D0) THEN
         P6=DABS((G(7,N-1)-G(7,N+1))/G(7,N))
      ELSE
         P6=0.D0
      ENDIF
      IF(P1.GT.UR.OR.P2.GT.UL.OR.P3.GT.UP.OR.P4.GT.UT.OR.P5.GT.UM.OR.
     $ P6.GT.UGR) GO TO 16

             
      SALTA=.FALSE.       
csanti    

      CALL LEVMES(N)
   16 CONTINUE
      IF(N.LT.(MAXME-1)) GO TO 11
C********************* FINE FASE SFOLTIMENTO ***************************
C
CC********** LEVA    MESH CHE DISTANO MENO DI 1.D-012 IN MASSA ***********
   17 CONTINUE
      MAM=MAXME-1
      DO 18 K=1,MAM
c      IF(((G(5,K+1)-G(5,K))/(G(5,K+1)+G(5,K))).GT.1.D-12) GO TO 18
      IF((G(5,K+1)-G(5,K)).GT.1.D-13) GO TO 18
      DO LL=K,MAM
         DO NN=1,7
            G(NN,LL)=G(NN,LL+1)
         END DO
         DO NN=1,MAXNE-1
            XXX(NN,LL)=XXX(NN,LL+1)
         END DO
         CARXY(1,LL)=CARXY(1,LL+1)
         CARXY(2,LL)=CARXY(2,LL+1)
      END DO
      MAXME=MAXME-1
      GO TO 17
   18 CONTINUE
       G5MIN=G(5,MAXME-1)-G(5,1)
       DO K=1,MAXME-1
        DELG5=G(5,K+1)-G(5,K)
        IF(DELG5.LE.G5MIN)THEN
          G5MIN=DELG5
          KMINIMO=K
        ENDIF
       ENDDO	
        WRITE(22,'(2I5,1P,4D20.13)')NMD,K,G5MIN,G(5,KMINIMO),
     #                             G(5,KMINIMO+1)
CC****** IN CASO DI DEBUG STAMPA FISICA E CHIMICA MODIFICATA ************
      IF(IPRALL.EQ.0)RETURN
      WRITE(2,100)
      DO K=1,MAXME
         WRITE(2,101)K,(G(J,K),J=1,6),(XXX(I,K),I=1,6)
      END DO
  100 FORMAT(1X,'FISICA DOPO REZONING: K - R - L - P - T - M - DGRAD - H
     #- HE3 - HE4 - C - N - O',/)
  101 FORMAT(1X,I4,1P,12E10.3)
      RETURN
      END

      SUBROUTINE PASTEM(IRSTA)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      
      common/NOSUB/IPRESSSTART,IPRESSSTOP,IPRESSSWITCH


C **********************************************************************
  101 FORMAT(1X,1P,9D8.1)
  110 FORMAT(2X,'RAGGIO',2X,'X-CEN',3X,'PRESS',3X,'TEMP',4X,
     *'ELIO',4X,'CNO',5X,'LUM',5X,'DMMAX',3X,'MPUNTO')
  120 FORMAT(2X,'RAGGIO',2X,'X-CEN',3X,'PRESS',3X,'TEMP',4X,
     *'ELIO',4X,'CAR',5X,'LUM',5X,'DMMAX',3X,'MPUNTO')
C  *********************************************************************
      HT1V=HT1
C     TIM=GRAVC0*DSQRT(4.D0*PI*STBOLZ)*(EMTOT*EMTOT)*(10.**(2*TEFF))/
C    #    (G(2,MAXME)**(3.D0/2.D0))*3.170979D+10  ! TEMPO SCALA KH
      TIM=1.D+10  ! VALORE GRANDE A KAZZO PER FARE IL MINIMO ALLA FINE
      TIMR=TIM
      TIML=TIM
      TIMMA=TIM
      TIMX=TIM
      TIMP=TIM
      TIMT=TIM
      TMPMA=TIM
      TIMHEL=TIM
      TIMCNO=TIM
      SMASS=EMTOT/SUNMg33
C******* CONTROLLA VARIAZIONI MASSIME DI R - L - P - T *****************
      DIPV=0.D0
      DITV=0.D0
      DO 7 K=2,MAXME-1
         DO 2 L=2,MAXMV
            IF(G(5,K)-GG(5,L))3,3,2
    2    CONTINUE
         GO TO 7
    3    K1=L-1
         IF((G(6,K-1)*G(6,K)).LE.0.D0) GO TO 7
         IF((G(6,K+1)*G(6,K)).LE.0.D0) GO TO 7
         IF((G(6,K )*GG(6,L)).LE.0.D0) GO TO 7
         RAT=(G(5,K)-GG(5,K1))/(GG(5,L)-GG(5,K1))
         PP=GG(3,K1)+RAT*(GG(3,L)-GG(3,K1))
         TT=GG(4,K1)+RAT*(GG(4,L)-GG(4,K1))
         DIP=DABS((PP-G(3,K))/PP)*100.D0
         DIT=DABS((TT-G(4,K))/TT)*100.D0
         IF(DIP.GT.DIPV)THEN
            DIPV=DIP
         ENDIF
         IF(DIT.GT.DITV)THEN
            DITV=DIT
         ENDIF
    7 CONTINUE
C******* CONTROLLI IN COMBUSTIONE CENTRALE DI IDROGENO  ****************
      IF(XXX(1,1).GT.0.D0) THEN
         VALR=3.D-4
        IF(FRAZ.LT.0.9991D0) THEN
            VALP=20.D0
            VALUE=10.D0
         ELSE
            VALP=150.D0
            VALUE=150.D0
         ENDIF
         IF(G(6,1).GT.0.D0) THEN
           VALX=.4D0-(1.D0-XXX(1,1))/5.D0
         ELSE
           VALX=.15D0-(1.D0-XXX(1,1))/10.D0
         ENDIF
        IF(ETA.GT.0.d0.AND.XXX(1,1).LT.XINI*0.97D0)THEN
          valp= 48.d0
          value=19.d0
          VALL=0.0006D0
          VALX=.7D0-(1.D0-XXX(1,1))/5.D0
         IF(XXX(1,1).LT.1.D-1.and.SMASS.GE.1.0D0.and.SMASS.LT.1.2D0)THEN
          valp= 43.d0
          value=20.d0
          VALL=0.001D0
          VALX=.7D0-(1.D0-XXX(1,1))/25.D0
         ENDIF
         IF(XXX(1,1).LT.1.D-1.and.SMASS.GE.1.2D0)THEN
          valp= 55.d0
          value=35.d0
          VALL=0.001D0
          VALX=.7D0-(1.D0-XXX(1,1))/25.D0
         ENDIF
         IF(XXX(1,1).LT.1.D-3.and.SMASS.GE.1.2D0)THEN
          valp= 65.d0
          value=45.d0
          VALL=0.001D0
          VALX=.7D0-(1.D0-XXX(1,1))/25.D0
         ENDIF
c        IF(ZINI.GT.1.d-2.AND.XXX(1,1).LT.1.D-3)THEN
c         write(*,*)'pipppppo'
c         valp= 22.d0
c         value=10.d0
c         VALL=0.001D0
c         VALX=.7D0-(1.D0-XXX(1,1))/15.D0
c        ENDIF
        ENDIF
         VALL=0.002D0
         IF(G(4,1).LT.8.D0) VALL=.02D0

!         DXDTM=(XCE(1,2)-XCE(1,1))/HT1V
c############################  per limitare passo temporale (vd. main) #####
         DXDTM=(XNM(1)-XCE(1,1))/HT1V
c###########################################################################
         IF(EMTOT.GT.12.D0) THEN
            TAGX=1.D-15
         ELSE
            TAGX=1.D-05
         ENDIF
         IF(DABS(DXDTM).GT.0.D0) THEN
            IF(XXX(1,1).GT.TAGX) THEN
               TIMX=DABS(VALX*XNM(1)/DXDTM)
            ELSE
               TIMX=DABS(VALX*TAGX/DXDTM)
            ENDIF
         ENDIF
         goto 456
C###################################################################
C#######riduce il passo temporale dopo pres.######################## 
         dhe3=(xxx(2,1)-xxv(2,1))/ht1v
         xl=0.12d-13
         dhe3lim=1.d-14
             if(nmlim.eq.0)then
                xxxl=3.0d-04
             else
                xxxl=0.1d-03
             endif
         if(smass.le.0.75d0.and.smass.gt.0.58d0)then
             if(nmlim.eq.0)then
                xxxl=3.3d-04
             else
                xxxl=0.1d-03
             endif
             xl=0.24d-13
             dhe3lim=1.d-15
         elseif(smass.le.0.58d0.and.smass.gt.0.53d0)then
             if(nmlim.eq.0)then
                xxxl=0.55d-03
             else
                xxxl=0.3d-03
             endif
             xl=0.12d-12
             dhe3lim=1.d-16
             if(xxx(2,1).gt.0.8d-03)xl=1.d-13
         elseif(smass.le.0.53d0.and.smass.gt.0.36d0)then
             xxxl=0.7d-03
             xl=0.12d-12
             dhe3lim=1.d-16
         elseif(smass.le.0.36d0.and.smass.ge.0.30d0)then
             if(nmlim.eq.0)then
                xxxl=1.6d-03
             else
                xxxl=3.0d-03
             endif
             xl=0.12d-12
             dhe3lim=1.d-15
         elseif(smass.lt.0.30d0)then
             xxxl=0.2d-03
             xl=0.12d-13
             dhe3lim=1.d-15
         endif
         if(xxx(2,1).gt.xxxl)then
            if(nmlim.eq.0)then
               dhe3f=dhe3
               if(dhe3.le.dhe3lim)then
                  nmlim=nmd
                  dhe3f=xl
               endif
            endif
            IF(SMASS.Lt.1.D0.AND.SMASS.GT.0.56D0)THEN
               vall=vall*(1.d13*dabs(dhe3f))
               valp=valp*(3.d11*dabs(dhe3f))
               value=value*(2.d11*dabs(dhe3f))
            ELSEIF(SMASS.Le.0.56D0.AND.SMASS.GT.0.32D0)THEN
               vall=vall*(1.d12*dabs(dhe3f)*smass/0.5d0)
               valp=valp*(1.d10*dabs(dhe3f)*smass/0.5d0)
               value=value*(1.d10*dabs(dhe3f)*smass/0.5d0)
            ELSEIF(SMASS.Le.0.32D0)THEN
               vall=vall*(1.d12*dabs(dhe3f)*smass/0.5d0)
               valp=valp*(1.d10*dabs(dhe3f)*smass/0.5d0)
               value=value*(1.d11*dabs(dhe3f)*smass/0.5d0)                     
            ENDIF
         endif
 456    continue
C#######################################################################   
C*******   CONTROLLI IN COMBUSTIONE CENTRALE DI ELIO *******************
      ELSEIF(XXX(3,1).GT.0.D0) THEN
         IF(TEFF.LT.4.0D0) THEN
            VALP=7.D0
         ELSE
            VALP=3.D0
         ENDIF
         VALUE=5.D0
         VALR=3.0D-4
c=========================== prova elimina ht1/2
         if(g(6,1).lt.0.d0)then
          VALP=3.5D0
          VALR=1.5D-4
        endif
c========================end prova elimina ht1/2
         IF(XXX(1,MAXME-1).LT.0.4D0) THEN
            VALP=80.D0
            VALUE=50.D0
         ENDIF
         IF(XXX(1,MAXME-1).LE.0.D0) THEN
            VALP=80.D0
            VALUE=50.D0
         ENDIF
         VALX=.3D0-(1.D0-XXX(3,1))/10.D0
         VALL=0.005D0
c        IF(SMASS.GT.5.0D0.and.ZINI.LE.4.d-5)THEN
c              vall=vall*0.05d0
c              valp=valp*0.02d0
c              value=value*0.05d0
c             VALX=VALX*0.02d0
c        ENDIF
        IF(ETA.GT.0.d0)THEN
          valp= 35.d0
          value=15.d0
          VALL=0.0006D0
c          IF(STARTMASS.GE.1.8D0)THEN
c            valp= 75.d0
c            value=35.d0
c            VALL=0.006D0
c          ENDIF
        ENDIF
         DXDTM=(XNM(3)-XCE(3,1))/HT1V
         IF(DABS(DXDTM).GT.0.D0) THEN
            IF(XXX(3,1).GT.1.D-5) THEN
               TIMX=DABS((VALX*XNM(3))/DXDTM)  !/3.D0
            ELSE
               TIMX=DABS(VALX*1.D-5/DXDTM)
            ENDIF
         ENDIF

C******* CONTROLLI IN COMBUSTIONE CENTRALE DI CARBONIO  ****************
      ELSEIF(XXX(4,1).GT.0.D0) THEN
         VALR=3.D-4
        VALP=2.D0
         VALUE=1.D0
         IF(EMTOT.GT.24.D0) THEN
            VALP=2.D0
            VALUE=1.D0
         ENDIF
        IF(ETA.GT.0.d0)THEN
          valp= 45.d0
          value=28.d0
        ENDIF
         IF(G(6,1).GE.0.D0) THEN
           VALX=.4D0-(1.D0-XXX(4,1))/5.D0
         ELSE
           VALX=.15D0-(1.D0-XXX(4,1))/10.D0
         ENDIF
         VALL=0.001D0
         DXDTM=(XCE(4,2)-XCE(4,1))/HT1V
         IF(DABS(DXDTM).GT.0.D0) THEN
            IF(XXX(4,1).GT.1.0D-03) THEN
               TIMX=DABS(VALX*XCE(4,2)/DXDTM)
            ELSE
               TIMX=DABS(VALX*1.D-3/DXDTM)
            ENDIF
         ENDIF
      ELSE
         VALR=3.D-4
        VALP=2.D0
         VALUE=1.D0
         VALL=0.001D0
      ENDIF
C***********************************************************************

c      DMDT=(G(5,MAXME)-GG(5,MAXMV))/HT1V
      DMDT=(4.D-13)*REIMAS*ETA*(10.D0**(ELLOG*1.5D0))/
     #                           (EMTOV*(10.D0**(2.D0*TEFFV)))      

      DLDTM=(DLOG10(G(2,MAXME))-DLOG10(GG(2,MAXMV)))/HT1V
      DRDTM=(TEFF-TEFFV)/HT1V
      DPDTM=DIPV/HT1V
      DTDTM=DITV/HT1V
C*******L=.03-R=.002-X=.01-P=.1-T=.03**********************************
      DMP=DMDT
      LOLO=LOBOR
      lillo=maxme-int(dfloat(maxme)/6.5d0)
c      IF((MAXME-LOBOR).LT.30) LOLO=MAXME-30

      DMT=G(5,MAXME)-G(5,LOLO)
c      write(*,*)lolo,lillo
c      write(*,'(i6,3d20.12)')NMD,DMP,dmt
      F3ALFA=dlog10(FL3A*(10.d0**ellog))
      AMAS=2.d0
      if(nshut.gt.0)AMAS=1.7d0
      IF(ETA.GT.0.D0.AND.DMP.GT.0.d0.and.xxx(1,1).lt.1.d-29.
     #   and.DMT.gt.2.d-8.and.STARTMASS.LT.AMAS
     #   .and.xxx(3,1).gt.0.d0)THEN
         IF(LOLO.LT.lillo) THEN
            TMPMA=0.30*DMT/(DMP*startmass)
         ELSE
c           write(*,*)'di sotto'
            TMPMA=0.25*DMT/(DMP*startmass)
         ENDIF
      ENDIF
      IF(DABS(DMDT).GT.0.D0) THEN
         IF(XXX(1,MAXME-1).GT.0.D0) THEN
            IF(XXX(1,1).GT.0.D0) THEN
             TIMMA=DABS(.0005*G(5,MAXME)/(DMDT*startmass))
           ELSE
c             write(*,*)ellog
             TIMMA=DABS(.000005*G(5,MAXME)/(DMDT*startmass))
           ENDIF
           IF(ZINI.GT.0.6D-2.AND.FL3A.GT.7.D-6.AND.
     #         STARTMASS.LT.AMAS)TIMMA=TIMMA*2.D0
         ELSE
            TIMMA=DABS(.002*G(5,MAXME)/(DMDT*startmass))
         ENDIF
      ENDIF
      IF(DABS(DLDTM).GT.0.D0) TIML=DABS(VALL/DLDTM)
      IF(DABS(DRDTM).GT.0.D0) TIMR=DABS(VALR*TEFF/DRDTM)
      IF(DABS(DPDTM).GT.0.D0) TIMP=DABS(VALP/DPDTM)
      IF(DABS(DTDTM).GT.0.D0) TIMT=DABS(VALUE/DTDTM)
C***********************************************************************
      IF(FL3A.GT.1.D-03) THEN
         PPRO=(GENHE+GENHEV)/2.D0
         DEHEL=DABS((GENHE-GENHEV)/(HT1V*PPRO))
         COHEL=0.10D0
         IF(XXX(3,1).LE.0.D0) COHEL=.03D0
         IF(DEHEL.GT.0.D0)TIMHEL=COHEL/DEHEL
      ENDIF
CCC      IF(FLCAR.GE.1.D-3) THEN
      IF(FLCAR.GE.1.D-4) THEN
         VALMED=(GENCARV+GENCAR)/2.D0
         DECAR=DABS((GENCARV-GENCAR)/(HT1V*VALMED))
         COCAR=.008D0
         IF(XXX(4,1).LE.0.D0)COCAR=.01D0
         IF(DECAR.GT.0.D0)TIMCAR=COCAR/DECAR
      ENDIF
      IF(FLCNO.GT.0.03D0) THEN
         ECN=(GENCNO+GENCNOV)/2.D0
         DECNO=DABS((GENCNO-GENCNOV)/(HT1V*ECN))
         COCNO=0.15D0
         IF(DECNO.GT.0.D0)TIMCNO=COCNO/DECNO
      END IF
CCC      IF(FLCAR. LT. 1.D-03) THEN
      IF(FLCAR. LT. 1.D-4) THEN
         WRITE(2,110)
         WRITE(*,110)
         HT1=DMIN1(TIM,TIMR,TIMX,TIMP,TIMT,TIMHEL,TIMCNO,TIML,TIMMA,
     $             TMPMA)
         WRITE(*,101)TIMR,TIMX,TIMP,TIMT,TIMHEL,TIMCNO,TIML,TMPMA,TIMMA
         WRITE(2,101)TIMR,TIMX,TIMP,TIMT,TIMHEL,TIMCNO,TIML,TMPMA,TIMMA
      ELSE
         WRITE(2,120)
         WRITE(*,120)
         HT1=DMIN1(TIM,TIMR,TIMX,TIMP,TIMT,TIMHEL,TIMCAR,TIML,TIMMA,
     $             TMPMA)
         WRITE(*,101)TIMR,TIMX,TIMP,TIMT,TIMHEL,TIMCAR,TIML,TMPMA,TIMMA
         WRITE(2,101)TIMR,TIMX,TIMP,TIMT,TIMHEL,TIMCAR,TIML,TMPMA,TIMMA         
      ENDIF
C*********** CONTROLLA VARIAZIONI BRUSCHE NEL PASSO TEMPORALE **********
         IF(IGONG.NE.0.AND.TEMPO.GT.1.0D8.AND.ETA.EQ.0)then
          WRITE(*,*)'LIMITING STEP TO HT/10 - FGONG FILES IN A GIVEN R R
     #ANGE'
       HT1=HT1/10.0d0
         ENDIF
C***********************************************************************
C---- per calibrazione solare
c      if(tempo.gt.4.45d09)ht1=ht1/10.d0
c      if(tempo.gt.4.5d09)ht1=ht1/2.3d0
c      if(tempo+ht1.gt.4.57d9)then
c        ht1=(4.57d9-tempo)
c	return
c      endif
c      if(tempo.gt.4.57d09)stop
C---- per calibrazione solare
c      IF(XXX(3,1).LE.0.D0.and.ETA.GT.0.d0)HT1=2.0d0*HT1

c      if(xxx(1,1).lt.1.d-4.and.xxx(1,1).gt.0.d0.and
c     #   .nshut.ne.0)ht1=ht1/2.d0
c      if(xxx(1,1).gt.1.d-4.and.nshut.ne.0.and.ht1.gt.1.d7)ht1=ht1/2.d0
c      if(xxx(3,1).lt.1.d-10.and.xxx(3,1).gt.0.d0)ht1=ht1/2.0d0
c      if(indiff.gt.0.and.xxx(1,1).lt.4.d-1.and.ht1.gt.5.d6)then
c         ht1=5.d6
c      endif
c      IF(XXX(3,1).LE.0.D0)HT1=HT1/2.d0
c      IF(XXX(3,1).GT.0.8D0.AND.FL3A.GT.2.D-5)
c     #HT1=HT1*2.5D0
      IF(XXX(3,1).LE.0.0D0.AND.FL3A.LT.1.D0.AND.FLCNOV.LT.1.D-3)
     #HT1=HT1*1.8D0
      IF(XXX(3,1).LE.0.0D0.AND.FL3A.LT.1.D0.AND.FLCNOV.GT.1.D-3)
     #HT1=HT1/4.8D0
c      IF(XXX(3,1).GT.0.8D0.AND.FL3A.GT.1.D-4.AND.ZINI.LT.4.d-5)
c     #HT1=HT1*0.5D0
C---- modifica per D-burning in presequenza
c      IF(startmass.LT.0.7d0)THEN
c       IF(XXX(1,1).GE.XINI*(1.D0-1.D-4)) ht1=ht1/(40.04D0-114.85D0*smass
c    # +98.669D0*(smass*smass))
c      ELSE
c       IF(XXX(1,1).GE.XINI*(1.D0-1.D-4))ht1=ht1/7.d0
c      ENDIF
C------------------------------------------------
      IF(startmass.LT.0.7d0)THEN
       IF(XXX(1,1).GE.XINI*(1.D0-1.D-4)) ht1=ht1/(40.04D0-114.85D0*smass
     # +98.669D0*(smass*smass))
      ELSE
       IF(XXX(7,1).GE.1.D-15)THEN
        ht1=ht1/7.d0
       ELSE
        IF(XXX(1,1).GE.XINI*(1.D0-1.D-4)) ht1=ht1/3.D0
       ENDIF 
      ENDIF
c       IF(XXX(1,1).LE.0.7d0.and.ht1.gt.1.0d7) ht1=ht1/4.d0
C------------------------------------------------
      IF(XXX(3,1).LE.0.D0)HT1=HT1*4.D0

      HT2=HT1V*1.5D0
      IF(startmass.gt.12.d0)ht1=ht1/2.d0
      IF(XXX(3,1).LE.0.95D0.AND.XXX(3,1).GE.0.7D0.
     #and.g(6,1).gt.0.d0) HT2=HT1V*1.2D0
      IF(XXX(3,1).LE.0.D0)HT2=HT1V*1.2D0
      IF(XXX(3,1).LE.0.D0.AND.FL3A.GT.2.D0)HT2=HT1V*1.05D0
      IF(HT1.GT.HT2)HT1=HT2
      HT3=HT1V/3.0D0
      IF(XXX(3,1).LE.0.D0.AND.FL3A.GT.2.D0)HT3=HT1V/1.05D0
      IF(HT1.LT.HT3)HT1=HT3
c      if(dabs(dm).ge.2.d-5)ht1=ht1*dabs(1.d-5/dm)
c      if(xxx(1,1).le.0.d0.and.xxx(3,1).gt.9.d-1.and.g(6,1).lt.0.d0.
c     #and.fl3a.lt.1.d-1.and.ht1.lt.5.d2.and.eta.gt.0.d0)ht1=5.d2
c===== modifica per eliminare subatmosfera
      if (nmd.lt.IPRESSSTOP+1.and.IPRESSSWITCH.eq.0) then
        HT1=1.D+1
        print *,'*****************TIME***************'  
      end if       
c===== FINE modifica per eliminare subatmosfera
c       if(nmd.eq.25)fraz=(fraz+1.d0)/2.d0
c       if(ht1.gt.1.d7)ht1=1.d7

      RETURN
      END

      SUBROUTINE QUATM(NABLA) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'atmosfer.2p1'
      INCLUDE 'consts.2p1'
      INTEGER IPRESSSTART,IPRESSSTOP,IPRESSSWITCH
      LOGICAL*1 CALCOLA,VARIATO 
      COMMON/CUBICA/A(4),B(4),ALFA,BETA,GAMMA,DELTA 

      common/NOSUB/IPRESSSTART,IPRESSSTOP,IPRESSSWITCH

      DIMENSION ELLA(4),ETTA(4),SERV(4),TT1(2),EL1(2), 
     $          C1(4),C2(4),AT(4,2,2),PTR(3,4,4),X(MAXNE) 
      DIMENSION RATM(2000),PATM(2000),TATM(2000),VMATM(2000) 
      DIMENSION VBP(3),VBT(3),VBR(3) 
      DATA RAPPOL/0.20/,RAPLMA/0.10/,RAPLMI/0.01/ 
      DATA RAPPOT/0.10/,RAPTMA/0.02/,RAPTMI/0.001/ 
      DATA DRATT/0.0005/,DRAPP/0.005/ 
      DATA DELLE/0.00001/,DETTE/0.000005/ 
      DATA LALLA/1/,INI/0/,XIDROV/0.0000D0/ 
      DO K=1,MAXNE-1 
         X(K)=XXX(K,MAXME-1) 
      END DO 
c      write(*,*)fraz
      BOTT=1.D0-FRAZ 
      VARIATO=.FALSE. 
c===== modifica per eliminare subatmosfera

      IF (nmd.le.IPRESSSTOP.and.IPRESSSWITCH.eq.0) then
        INI=0
      END IF
c===== FINE  modifica per eliminare subatmosfera
      DM=0.d0
      IF(ETA.GT.0.D0.AND.NMD.GE.IPRESSSTOP)THEN
       DM=-(4.D-13)*REIMAS*HT1*ETA*(10.D0**(ELLOG*1.5D0))/(EMTOT*
     #                           (10.D0**(2.D0*TEFF)))      
       DM=DABS(DM)
      ENDIF
C*********************************************************************** 
      IF(ISUB.EQ.1.OR.NABLA.EQ.3) THEN 
        EMSOLA=EMTOT/SUNMg33
         CALL NATMOS(ELLOG,TEFF,EMSOLA,X,ZINI,BOTT,ALPHA, 
     $               RATM,PATM,TATM,VMATM,NATM,1,1,0,ISINGOLA,MAXNE,nmd) 
         RETURN 
      ENDIF 
C*********************************************************************** 
C*********************************************************************** 
      IF(ISINGOLA.EQ.1) THEN 
C        IF(X(1).LT.XIDROV) THEN 
C           XGAR=X(1)+X(3) 
C           XVAF=X(1) 
C           X(1)=XIDROV*0.997D0 
C           IF(X(1).LT.XVAF)X(1)=XVAF 
C           X(3)=XGAR-X(1) 
C        ENDIF 
        EMSOLA=EMTOT/SUNMg33
         DO K=1,3 
            IF(K.EQ.1) THEN 
               ELLE=ELLOG 
               TEMP=TEFF 
            ELSEIF(K.EQ.2) THEN 
               ELLE=ELLOG+DELLE 
               TEMP=TEFF 
            ELSEIF(K.EQ.3) THEN 
               ELLE=ELLOG 
               TEMP=TEFF+DETTE 
            ENDIF 
            KK=K-1 
            CALL NATMOS(ELLE,TEMP,EMSOLA,X,ZINI,BOTT,ALPHA, 
     $            RATM,PATM,TATM,VMATM,NATM,0,0,KK,ISINGOLA,MAXNE,nmd) 
            VBP(K)=DLOG10(PATM(NATM)) 
            VBT(K)=DLOG10(TATM(NATM)) 
            VBR(K)=DLOG10(RATM(NATM)) 
         END DO 
C        IF(ITERAZ.EQ.1) XIDROV=X(1) 
         PB=10.**(VBP(1)-17.D0) 
         TB=10.**(VBT(1)-6.D0) 
         RB=10.**VBR(1) 
         DPL=(VBP(2)-VBP(1))/DELLE 
         DTL=(VBT(2)-VBT(1))/DELLE 
         DRL=(VBR(2)-VBR(1))/DELLE 
         DPT=(VBP(3)-VBP(1))/DETTE 
         DTT=(VBT(3)-VBT(1))/DETTE 
         DRT=(VBR(3)-VBR(1))/DETTE 
C        WRITE(2,233)ELLOG,TEFF,PB,TB,RB,DPL,DPT,DRL,DRT,DTL,DTT 
C        WRITE(*,233)ELLOG,TEFF,PB,TB,RB,DPL,DPT,DRL,DRT,DTL,DTT 
C 233    FORMAT(2F8.5,1P,9E10.3) 
         RETURN 
      ENDIF 
C*********************************************************************** 
      CALCOLA=.FALSE. 
      IF(INI.EQ.0) THEN 
         IF(LALLA.EQ.0) THEN 
            EL(1)=4.51 
            EL(2)=4.53 
            EL(3)=4.55 
            EL(4)=4.57 
            TE(1)=3.560 
            TE(2)=3.580 
            TE(3)=3.600 
            TE(4)=3.620 
           PM1=EMTOT/SUNMg33
           IF(ISTART.EQ.5) THEN
            PM2=PM1*1.1D0
           ELSE
             PM2=PM1*0.95 
            ENDIF
           CALCOLA=.TRUE. 
         ELSEIF(ISTART.NE.1) THEN 
            EL(1)=ELLOG-0.090  
            EL(2)=EL(1)+0.06 
            EL(3)=EL(2)+0.06 
            EL(4)=EL(3)+0.06 
            TE(1)=TEFF-0.030 
            TE(2)=TE(1)+0.02 
            TE(3)=TE(2)+0.02 
            TE(4)=TE(3)+0.02 
           PM1=EMTOT/SUNMg33
           IF(ISTART.EQ.5) THEN
            PM2=PM1*1.1D0
           ELSE
             PM2=PM1*0.95D0 
            ENDIF
            CALCOLA=.TRUE. 
cc===== modifica per eliminare subatmosfera
         ELSEIF(ISTART.EQ.1.and.IPRESSSWITCH.EQ.0) THEN             
            EL(1)=ELLOG-0.090 
            EL(2)=EL(1)+0.06 
            EL(3)=EL(2)+0.06 
            EL(4)=EL(3)+0.06 
            TE(1)=TEFF-0.030 
            TE(2)=TE(1)+0.02 
            TE(3)=TE(2)+0.02 
            TE(4)=TE(3)+0.02 
           PM1=EMTOT/SUNMg33
            PM2=PM1*0.95D0 
            CALCOLA=.TRUE. 
c===== FINE  modifica per eliminare subatmosfera
         ENDIF 
         INI=1 
      ENDIF 
      
      IF(ISTART.EQ.5.AND.EMTOT/SUNMg33.GT.PM2) THEN
         PM1=EMTOT/SUNMg33
         PM2=PM1*1.1D0
         CALCOLA=.TRUE.
      ENDIF 

      IF(ISTART.NE.5.AND.EMTOT/SUNMg33.LE.PM2) THEN
	 PM1=EMTOT/SUNMg33
          PM2=PM1*0.95D0 
          CALCOLA=.TRUE.
      ENDIF

      IF(CALCOLA.EQV..TRUE.) THEN 
         DO JL=1,4 
            DO JT=1,4 
               CALL NATMOS(EL(JL),TE(JT),PM1,X,ZINI,BOTT,ALPHA, 
     $              RATM,PATM,TATM,VMATM,NATM,0,0,0,ISINGOLA,MAXNE,nmd) 
               U(1,JL,JT)=DLOG10(PATM(NATM)) 
               U(2,JL,JT)=DLOG10(TATM(NATM)) 
               U(3,JL,JT)=DLOG10(RATM(NATM)) 
               WRITE(2,333)EL(JL),TE(JT),(U(IT,JL,JT),IT=1,3),PM1 
c               WRITE(*,333)EL(JL),TE(JT),(U(IT,JL,JT),IT=1,3),PM1 
               IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0) THEN 
                  CALL NATMOS(EL(JL),TE(JT),PM2,X,ZINI,BOTT,ALPHA, 
     $                  RATM,PATM,TATM,VMATM,NATM,0,0,0,ISINGOLA,
     $                  MAXNE,nmd) 
                  V(1,JL,JT)=DLOG10(PATM(NATM)) 
                  V(2,JL,JT)=DLOG10(TATM(NATM)) 
                  V(3,JL,JT)=DLOG10(RATM(NATM)) 
                  WRITE(2,333)EL(JL),TE(JT),(V(IT,JL,JT),IT=1,3),PM2 
c                  WRITE(*,333)EL(JL),TE(JT),(V(IT,JL,JT),IT=1,3),PM2 
               ENDIF 
            END DO 
         END DO 
      ENDIF 
C**********************CONTROLLO SULLA LUMINOSITA'********************** 
      DO WHILE (ELLOG.GT.EL(3).OR.ELLOG.LT.EL(2)) 
         VARIATO=.TRUE. 
         IF(ELLOG.LT.EL(2)) THEN 
            L=1 
         ELSE 
            L=0 
         END IF 
         DO I=1,3 
            JJ=I+1+L*(3-2*I) 
            JL=I*(1-L)+L*(JJ+1) 
            DO JT=1,4 
               DO JB=1,3 
                  U(JB,JL,JT)=U(JB,JJ,JT) 
                  IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0)THEN
                   V(JB,JL,JT)=V(JB,JJ,JT) 
                  ENDIF
	       END DO 
            END DO 
            EL(JL)=EL(JJ) 
         END DO 
         IB=3-L 
         IA=IB+1-2*L 
         DO J=1,4 
            ELLA(J)=EL(IB)+(1-2*L)*RAPPOL*DABS((EL(3)-EL(2))/ 
     $           (U(2,2,J)-U(2,3,J))) 
            IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0)THEN
              ELLA(J)=EL(IB)+(1-2*L)*RAPPOL* 
     $                    DABS((EL(3)-EL(2))/(V(2,2,J)-V(2,3,J))) 
            ENDIF
         END DO 
         ORMAX=EL(IB)+RAPLMA*(1-2*L) 
         ORMIN=EL(IB)+RAPLMI*(1-2*L) 
         IF(L.EQ.0) THEN 
            EL(4)=DMIN1(ELLA(1),ELLA(2),ELLA(3),ELLA(4),ORMAX) 
            EL(4)=DMAX1(EL(4),ORMIN) 
         ELSE IF(L.EQ.1) THEN 
            EL(1)=DMAX1(ELLA(1),ELLA(2),ELLA(3),ELLA(4),ORMAX) 
            EL(1)=DMIN1(EL(1),ORMIN) 
         ENDIF 
         IELIA=EL(IA)*10000.D0 
         EL(IA)=(IELIA)/10000.D0 
         WRITE(2,931)EL(IA),ELLOG,TEFF 
         WRITE(*,931)EL(IA),ELLOG,TEFF 
         DO J=1,4 
            CALL NATMOS(EL(IA),TE(J),PM1,X,ZINI,BOTT,ALPHA, 
     $            RATM,PATM,TATM,VMATM,NATM,0,0,0,ISINGOLA,MAXNE,nmd) 
            U(1,IA,J)=DLOG10(PATM(NATM)) 
            U(2,IA,J)=DLOG10(TATM(NATM)) 
            U(3,IA,J)=DLOG10(RATM(NATM)) 
            WRITE(2,333)EL(IA),TE(J),(U(IT,IA,J),IT=1,3),PM1 
c            WRITE(*,333)EL(IA),TE(J),(U(IT,IA,J),IT=1,3),PM1 
            IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0) THEN 
               CALL NATMOS(EL(IA),TE(J),PM2,X,ZINI,BOTT,ALPHA, 
     $              RATM,PATM,TATM,VMATM,NATM,0,0,0,ISINGOLA,MAXNE,nmd) 
               V(1,IA,J)=DLOG10(PATM(NATM)) 
               V(2,IA,J)=DLOG10(TATM(NATM)) 
               V(3,IA,J)=DLOG10(RATM(NATM)) 
               WRITE(2,333)EL(IA),TE(J),(V(IT,IA,J),IT=1,3),PM2 
c               WRITE(*,333)EL(IA),TE(J),(V(IT,IA,J),IT=1,3),PM2 
            ENDIF 
         END DO 
      END DO 
C******************CONTROLLO SULLA TEMPERATURA************************** 
      DO WHILE (TEFF.GT.TE(3).OR.TEFF.LT.TE(2)) 
         VARIATO=.TRUE. 
         IF(TEFF.LT.TE(2)) THEN 
            L=1 
         ELSE 
            L=0 
         ENDIF 
         DO I=1,3 
            JJ=I+1+L*(3-2*I) 
            JT=I*(1-L)+L*(JJ+1) 
            DO JL=1,4 
               DO JB=1,3 
                  U(JB,JL,JT)=U(JB,JL,JJ) 
                  IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0)THEN
                   V(JB,JL,JT)=V(JB,JL,JJ) 
                  ENDIF
               END DO 
            END DO 
            TE(JT)=TE(JJ) 
         END DO 
         IB=3-L 
         IA=IB+1-2*L 
         DO J=1,4 
            ETTA(J)=TE(IB)+(1-2*L)*RAPPOT*DABS((TE(3)-TE(2))/ 
     $              (U(2,J,2)-U(2,J,3))) 
            IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0)THEN
               ETTA(J)=TE(IB)+(1-2*L)*RAPPOT* 
     $                      DABS((TE(3)-TE(2))/(V(2,J,2)-V(2,J,3))) 
            ENDIF
         END DO 
         ORMAX=TE(IB)+RAPTMA*(1-2*L) 
         ORMIN=TE(IB)+RAPTMI*(1-2*L) 
         IF(L.EQ.0) THEN 
            TE(4)=DMIN1(ETTA(1),ETTA(2),ETTA(3),ETTA(4),ORMAX) 
            TE(4)=DMAX1(TE(4),ORMIN) 
         ENDIF 
         IF(L.EQ.1) THEN 
            TE(1)=DMAX1(ETTA(1),ETTA(2),ETTA(3),ETTA(4),ORMAX) 
            TE(1)=DMIN1(TE(1),ORMIN) 
         ENDIF 
         IELIA=TE(IA)*10000.D0 
         TE(IA)=(IELIA)/10000.D0 
         WRITE(2,901) TE(IA),ELLOG,TEFF 
         WRITE(*,901) TE(IA),ELLOG,TEFF 
         DO J=1,4 
            CALL NATMOS(EL(J),TE(IA),PM1,X,ZINI,BOTT,ALPHA, 
     $           RATM,PATM,TATM,VMATM,NATM,0,0,0,ISINGOLA,MAXNE,nmd) 
            U(1,J,IA)=DLOG10(PATM(NATM)) 
            U(2,J,IA)=DLOG10(TATM(NATM)) 
            U(3,J,IA)=DLOG10(RATM(NATM)) 
            WRITE(2,333)EL(J),TE(IA),(U(IT,J,IA),IT=1,3),PM1 
c            WRITE(*,333)EL(J),TE(IA),(U(IT,J,IA),IT=1,3),PM1 
            IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0) THEN 
               CALL NATMOS(EL(J),TE(IA),PM2,X,ZINI,BOTT,ALPHA, 
     $               RATM,PATM,TATM,VMATM,NATM,0,0,0,ISINGOLA,
     $               MAXNE,nmd) 
               V(1,J,IA)=DLOG10(PATM(NATM)) 
               V(2,J,IA)=DLOG10(TATM(NATM)) 
               V(3,J,IA)=DLOG10(RATM(NATM)) 
               WRITE(2,333)EL(J),TE(IA),(V(IT,J,IA),IT=1,3),PM2 
c               WRITE(*,333)EL(J),TE(IA),(V(IT,J,IA),IT=1,3),PM2 
            ENDIF 
         END DO 
      END DO 
C****************STAMPA NUOVO QUATM SE CAMBIATO************************* 
      IF(VARIATO.EQV..TRUE.) THEN 
         WRITE(2,334) 
         DO JL=1,4 
            DO JT=1,4 
               WRITE(2,333)EL(JL),TE(JT),(U(IT,JL,JT),IT=1,3),PM1 
               IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0) THEN 
                  WRITE(2,333)EL(JL),TE(JT),(V(IT,JL,JT),IT=1,3),PM2 
               ENDIF 
            END DO 
         END DO 
         WRITE(2,334) 
      ENDIF 
C**************************CALCOLO GRIGLIATO**************************** 
      IF(DM.GT.1.D-10.AND.IPRESSSWITCH.EQ.0) THEN 
        RAT=(EMTOT/SUNMg33-PM2)/(PM1-PM2)
         DO J=1,3 
            DO K=1,4 
               DO L=1,4 
                  PTR(J,K,L)=V(J,K,L)+RAT*(U(J,K,L)-V(J,K,L)) 
               END DO 
            END DO 
         END DO 
      ELSE 
         DO J=1,3 
            DO K=1,4 
               DO L=1,4 
                  PTR(J,K,L)=U(J,K,L) 
               END DO 
            END DO 
         END DO 
      ENDIF 
C****************CALCOLO GRIGLIATINO************************************ 
      DDRATT=2.D0*DRATT 
      DDRAPP=2.D0*DRAPP 
      TT1(1)=TEFF+DRATT 
      TT1(2)=TEFF-DRATT 
      EL1(1)=ELLOG+DRAPP 
      EL1(2)=ELLOG-DRAPP 
      DO ITE=1,2 
         DO ILU=1,2 
            DO I=1,3 
               DO J=1,4 
                  DO K=1,4 
                     A(K)=PTR(I,J,K) 
                     B(K)=TE(K) 
                  END DO 
                  CALL NCUB 
                  SERV(J)=ALFA*TT1(ITE)*TT1(ITE)*TT1(ITE)+ 
     $                    BETA*TT1(ITE)*TT1(ITE)+GAMMA*TT1(ITE)+DELTA 
               END DO 
               DO J=1,4 
                  A(J)=SERV(J) 
                  B(J)=EL(J) 
               END DO 
               CALL NCUB 
               AT(I,ITE,ILU)=ALFA*EL1(ILU)*EL1(ILU)*EL1(ILU)+ 
     $                  BETA*EL1(ILU)*EL1(ILU)+GAMMA*EL1(ILU)+DELTA 
            END DO 
         END DO 
      END DO 
      RATT=(TEFF-TT1(2))/DDRATT 
      RATL=(ELLOG-EL1(2))/DDRAPP 
      DO K=1,4 
         C1(K)=AT(K,2,1)+RATT*(AT(K,1,1)-AT(K,2,1)) 
         C2(K)=AT(K,2,2)+RATT*(AT(K,1,2)-AT(K,2,2)) 
      END DO 
      PB=10.**(C2(1)+RATL*(C1(1)-C2(1))-17.D0) 
      TB=10.**(C2(2)+RATL*(C1(2)-C2(2))-6.D0) 
      RB=10.**(C2(3)+RATL*(C1(3)-C2(3))) 
      DPL=(C1(1)-C2(1))/DDRAPP 
      DTL=(C1(2)-C2(2))/DDRAPP 
      DRL=(C1(3)-C2(3))/DDRAPP 
      DO K=1,3 
         C1(K)=AT(K,1,2)+RATL*(AT(K,1,1)-AT(K,1,2)) 
         C2(K)=AT(K,2,2)+RATL*(AT(K,2,1)-AT(K,2,2)) 
      END DO 
      DPT=(C1(1)-C2(1))/DDRATT 
      DTT=(C1(2)-C2(2))/DDRATT 
      DRT=(C1(3)-C2(3))/DDRATT 
      IF(IPRALL.EQ.1)WRITE(2,222)PB,TB,RB,DPL,DPT,DTL,DTT,DRL,DRT 
       IF (nmd.ge.IPRESSSTART.and.nmd.le.IPRESSSTOP.
     #and.ipressswitch.eq.0.and.fraz.lt.(1-1.d-6))THEN
        IF(nmd.ge.5)THEN
          g(5,maxme)=vmatm(1)
          gg(5,maxme)=vmatm(1)
          fraz=g(5,maxme)/emtot
c         IPRESSSTOP=nmd+10
c        write(*,*)IPRESSSTOP,'--------------'
        ENDIF
       ENDIF
      RETURN 
  222 FORMAT(1X,'P',1P,E10.3,'  T',E10.3,'  R',E10.3,'  DP',2E10.3,'  DR 
     $',2E10.3,'  DT',2E10.3) 
  333 FORMAT(6F8.4,1P,3E10.3) 
  334 FORMAT(///) 
  931 FORMAT(3X,'SQUATM IN LUMINOSITA',3F9.4) 
  901 FORMAT(3X,'SQUATM IN T EFFETTIVA',3F9.4) 
      END 
