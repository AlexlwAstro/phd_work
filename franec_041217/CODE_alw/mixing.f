C*****************************************************************************
C  QUESTA ROUTINE CONTIENE: BCONV - MIXINGALL - MIXINGHE - OVERSH
C*****************************************************************************
C
      SUBROUTINE BCONV(NCONV,JBCONV)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      
      common/NOSUB/IPRESSSTART,IPRESSSTOP,IPRESSSWITCH
      COMMON/OVER/FACTOR

C***********************************************************************
      DIMENSION JBCONV(2,50)
      LOGICAL*1 CONTINUA
C***********************************************************************
      CONTINUA=.TRUE.
C********************  NSHUT=1 ATTIVA OVERSHOOTING  ********************
      DO NCONV=1,50
         JBCONV(1,NCONV)=0
         JBCONV(2,NCONV)=0
      END DO
      NCONV=0
      KBORD=0
      GICO=0.D0
C****************************  CORE CONVETTIVO  ************************
      IF(G(6,1).GT.0.D0) THEN
c         IF(XXX(1,1).GT.0.D0) THEN                    ! CORE DI IDROGENO
            K=1
            DO WHILE (G(6,K).GE.0.D0.AND.K.LT.MAXME)
               K=K+1
            END DO
            IF(K.GT.1) THEN
C***********************************************************************
C KBORD=K PER CONSISTENZA CON FRA3P0 - ATTENZIONE IL MESH K E' RADIATIVO
c KBORD=K FOR CONSISTENCY WITH FRA3PO - ATTENTION - MESH K IS RADIATIVE
c***********************************************************************
               KBORD=K
            ELSE
               KBORD=0
            ENDIF
            IF(KBORD.GE.MAXME-1) KBORD=MAXME-1
C*******************  EVOLUZIONE CON OVERSHOOTING  *********************
c           IF(NSHUT.GT.0.AND.XXX(1,1).GT.1.D-08.AND.BB(1).GT.1.D-2)THEN
           IF(XXX(1,1).GT.1.D-08.AND.BB(1).GT.1.D-2.
     #AND.FACTOR.GT.0.D0)THEN
               CONVCLAS=G(5,KBORD+1)/EMTOT
               RCONCLAS=G(1,KBORD+1)
               CALL OVERSH(K-1,KBORD,AMBDA,VMOVER)
               RCONOVER=G(1,KBORD+1)
               CONVOVER=G(5,KBORD+1)/EMTOT
               WRITE(2,323)RCONCLAS,CONVCLAS,RCONOVER,CONVOVER,KBORD
               WRITE(*,323)RCONCLAS,CONVCLAS,RCONOVER,CONVOVER,KBORD
            ENDIF
         NCONV=1
         JBCONV(1,1)=1
         JBCONV(2,1)=KBORD
         GICO=G(5,KBORD+1)
         K=KBORD
      ELSE
         K=1
         NCONV=0
      ENDIF
C**************************  SHELL CONVETTIVE  *************************
      DO WHILE(CONTINUA)
         MINT=0
         MEST=0
         DO WHILE (G(6,K).LT.0.D0.AND.K.LE.(MAXME-1))
            K=K+1
         END DO
         IF(K.LT.MAXME-1) MINT=K-1  ! BORDO INT. ZONA NCONV-ESIMA CONV.
         DO WHILE (G(6,K).GE.0.D0.AND.K.LE.(MAXME-1))
            K=K+1
         END DO
         MEST=K-1                     ! BORDO EXT. ZONA NCONV-ESIMA CONV.
c===== modifica per eliminare subatmosfera
        IF (nmd.gt.4.and.IPRESSSWITCH.eq.0) then
          solmas=fraz-fraz/5.d6
c         write(*,'(2f15.10)')G(5,MEST)/emtot,solmas
          if(G(5,MEST)/emtot.ge.solmas)MEST=MAXME-1
        ENDIF
c===== fine modifica per eliminare subatmosfera

C===========================================================================
         IF(MINT.GT.0) THEN
            IF((MEST-MINT).GT.3) THEN         ! (PER SHELL DI 3 MESH...)
C===========================================================================         
           if(nconv.gt.0) then
               IF(MINT.GT.(JBCONV(2,NCONV)+1)) THEN
                  NCONV=NCONV+1
                  JBCONV(1,NCONV)=MINT
               ENDIF
            
           else
                 NCONV=NCONV+1
                 JBCONV(1,NCONV)=MINT
            endif                 
C===========================================================================

               JBCONV(2,NCONV)=MEST
            ENDIF
         ENDIF
         IF(K.GE.MAXME-1) CONTINUA=.FALSE.
         IF(NCONV.GE.50) CONTINUA=.FALSE.
      END DO
C***********************************************************************
      IF(NCONV.GT.0) THEN
         WRITE(*,100)NCONV,((JBCONV(LL,II),LL=1,2),II=1,NCONV)
      END IF
C***********************************************************************
cC* EVOLUZIONE CON "UNDERSHOOTING"  *************************************
c      KENV1=MAXME-1
c      IF(NCONV.GT.1)KENV1=JBCONV(1,NCONV)
c      IF(NCONV.EQ.1.AND.G(6,1).LT.0.D0)KENV1=JBCONV(1,NCONV)
c      write(*,*)KENV1
c      IF(KENV1.LT.MAXME-1)THEN
c          CONVCLAS=G(5,KENV1)/EMTOT
c          RCONCLAS=G(1,KENV1)
c          CALL UNDERSH(KENV1,KENV2,AMBDA,VMOVER)
c         write(*,*)KENV2
c          RCONOVER=G(1,KENV2)
c          CONVOVER=G(5,KENV2)/EMTOT
c          WRITE(2,323)RCONCLAS,CONVCLAS,RCONOVER,CONVOVER,KENV2
c          WRITE(*,323)RCONCLAS,CONVCLAS,RCONOVER,CONVOVER,KENV2
c         JBCONV(1,NCONV)=KENV2
c      ENDIF
c***********************************************************************
C**********************  FORMATI E SCRITTURE  **************************
  100 FORMAT(1X,'# CONV.:',I5,' BORDI:',101I5)
  323 FORMAT('CLASSIC BORDER R:',1P,D10.2,' M:',D10.2,
     #  ' OVER BORDER R:',D10.2,' M:',D10.2,' MESH:',0P,I4)
C***********************************************************************
      RETURN
      END

      SUBROUTINE MIXINGALL(NCONV,JBCONV)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      DIMENSION JBCONV(2,50),XTOT(MAXNE)
C
C============================================================     
      DO JJ=1,NCONV
         JINI=JBCONV(1,JJ)
         JFIN=JBCONV(2,JJ)
        IF(JINI.EQ.1.AND.XXX(1,1).EQ.0.D0) GOTO 123
         DO J=1,MAXNE-1
            XTOT(J)=0.D0
         END DO
        DMC=0.D0
         DO L=JINI,JFIN
             DO I=1,MAXNE-1
                XTOT(I)=XTOT(I)+XXX(I,L)*(G(5,L+1)-G(5,L))
             END DO
             DMC=DMC+G(5,L+1)-G(5,L)
         END DO
         DO N=JINI,JFIN
            DO I=1,MAXNE-1
               XXX(I,N)=XTOT(I)/DMC
            END DO
         END DO
 123  CONTINUE
      END DO        
      RETURN
      END

      SUBROUTINE MIXINGHE
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL*1 CONTINUA
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      COMMON/FREMEM/ XTOT(MAXNE),PG(6),XSERV(MAXNE,LIM),
     &               O((MAXNE+1)*LIM-MAXNE-MAXNE*LIM-6)
C*************IBLK=1 BLOCCA I PULSI: IBLK=0 NO *************************
C* SE IBLK=1 ALLORA ISMOTH=1 AMMORBIDISCE I BP, ISMOTH = 0 BLOCCA I BP
      DATA LBCON/2/,IBLK/1/,ISMOTH/0/
      EMSOL=EMTOT/SUNMg33
      K=0
      IGRO=0
      IBOR=0
      KCORMIX=0
      GIKO=0.D0
      GSEMI=0.D0
      L1=MAXME
      LBIT=1
C********* TRATTAMENTO CONVEZIONE IN COMBUSTIONE CENTRALE DI ELIO ******
         IF( G(6,1).LT.0.D0 ) GO TO 25
         DO KK=1,MAXME
            DO JJ=1,MAXNE-1
               XSERV(JJ,KK)=XXX(JJ,KK)
            END DO
         END DO
         IF(EMSOL.LT.10.D0) THEN
           LBCON=64
           LBIT=64
         ELSE
            LBCON=16
            LBIT=16
         ENDIF
         L1=1
         L2=LBCON
         IF(L2.LT.2) L2=2
   10    CONTINUE
         IF( LBIT.LT.1) LBIT=1
         IF( LBCON.LT.2) LBCON=2
C************************** MIXING DA L1 A L2 **************************
         DO J=1,MAXNE-1
            XTOT(J)=0.D0
         END DO
         DMC=0.D0
         DO K=L1,L2
            DO J=1,MAXNE-1
               XTOT(J)=XTOT(J)+XXX(J,K)*(G(5,K+1)-G(5,K))
           END DO
            DMC=DMC+G(5,K+1)-G(5,K)
         END DO
         DO K=L1,L2
            DO J=1,MAXNE-1
               XXX(J,K)=XTOT(J)/DMC
            END DO
         END DO
C***********************************************************************
         CONTINUA=.TRUE.
         DO WHILE (CONTINUA)
C******************** TEST: RAD > AD DA L1 A L2 ************************
            K=L1-1
            ICONTI=0
            DO WHILE(ICONTI.EQ.0)
               K=K+1
               ZL=(G(2,K)+G(2,K+1))/2.D0
               PP=(G(3,K)+G(3,K+1))/2.D0
               TT=(G(4,K)+G(4,K+1))/2.D0
               EM=(G(5,K)+G(5,K+1))/2.D0
               TM=TT*1.D+06
               PR=PP*1.D+17
               CALL EOS(K,1,MAXNE,1,PR,TM,XXX,'MIXING  ',
     &                     RHO,DAD,CSP,O1,G1,DEL,EMUE,3)
               CALL NKAPPA(K,1,MAXNE,1,RHO,TM,XXX,CAP)
              DRAD=CRAD*CAP*PP*ZL/(EM*(TT**4))
C              DADO=DRAD-DAD
C              WRITE(*,340)K,L1,L2,XXX(3,K),DADO,XXX(5,K)
C              WRITE(2,340)K,L1,L2,XXX(3,K),DADO,XXX(5,K)
C 340          FORMAT(1X,3I4,1P,3D12.4)
               IF(DRAD.LT.DAD.OR.K.GE.L2) ICONTI=1
            END DO
            IF(K.GE.L2) THEN
C***************** AGGIUNGE LBIT MESH AL CONVETTIVO ********************
               IF(L2.GE.(MAXME-1)) RETURN
               L2=L2+LBIT
C************************** MIXING DA L1 A L2 **************************
               IF(LBIT.EQ.1) THEN
                  IF(XXX(1,L2).GT.0.D0) THEN
                     WRITE(*,*)'HYD ',L2
                     WRITE(2,*)'HYD ',L2
                  ENDIF
                  DO J=1,MAXNE-1
                     XTOT(J)=XTOT(J)+XXX(J,L2)*(G(5,L2+1)-G(5,L2))
                  END DO
                  DMC=DMC+G(5,L2+1)-G(5,L2)
                  DO K=L1,L2
                     DO J=1,MAXNE-1
                        XXX(J,K)=XTOT(J)/DMC
                     END DO
                  END DO
               ELSE
                  DO J=1,MAXNE-1
                     XTOT(J)=0.D0
                  END DO
                  DMC=0.D0
                  DO K=L1,L2
                     DO J=1,MAXNE-1
                        XTOT(J)=XTOT(J)+XXX(J,K)*(G(5,K+1)-G(5,K))
                     END DO
                     DMC=DMC+G(5,K+1)-G(5,K)
                  END DO
                  DO K=L1,L2
                     DO J=1,MAXNE-1
                        XXX(J,K)=XTOT(J)/DMC
                     END DO
                  END DO
               ENDIF
            ELSE
               CONTINUA=.FALSE.
            ENDIF
         END DO
C**************** MESH K RADIATIVO *************************************
         IF(L1.EQ.1.AND.LBIT.GT.1) THEN
            L2=L2-LBIT/2
            LBIT=LBIT/2
            DO KK=1,MAXME
               DO JJ=1,MAXNE-1
                  XXX(JJ,KK)=XSERV(JJ,KK)
               END DO
            END DO
            GO TO 10
         END IF
c        WRITE(2,234)L1,L2,K
c        WRITE(*,234)L1,L2,K
c 234    FORMAT(1X,'*** MIX DA:',I3,'  A:',I3,'  BORDO A:',I3,' ***')
         IF(L1.EQ.1) THEN
            KCORMIX = K+1
            GIKO = G(5,K+1)
         ENDIF
         XXMAX=XCE(3,1)
C***********************************************************************
         IF(L1.EQ.1.AND.IBLK.EQ.1.AND.XXX(3,1).GT.XXMAX.AND.
     #      XXMAX.GT.0.D0.AND.XXMAX.LT.0.15D0) THEN
            IBOR=1
         ELSE
            IBOR=0
         ENDIF
c        ibor=0  !  breathing pulses yes
         IF(IBOR.EQ.1) THEN
            DO J=1,MAXNE-1
               XTOT(J)=0.D0
            END DO
            DMC=0.D0
            DO KK=1,MAXME
               DO JJ=1,MAXNE-1
                  XXX(JJ,KK)=XSERV(JJ,KK)
               END DO
            END DO
            DO KL=1,L2
               DO J=1,MAXNE-1
                  XTOT(J)=XTOT(J)+XXX(J,KL)*(G(5,KL+1)-G(5,KL))
               END DO
               DMC=DMC+G(5,KL+1)-G(5,KL)
               GROHE=XTOT(3)/DMC
C              WRITE(*,222)KL,L2,GROHE,XCE(3,1)
C 222          FORMAT(1X,2I4,1P,2D12.4)
               IF(GROHE.GT.XCE(3,1)) THEN
                  IF(ISMOTH.EQ.1) THEN
C********** QUESTO BLOCCHETTO SERVE PER SMOOTTARE IL CONVETTIVO ********
                     IGRO=KL
                     GICA=GIKO
                     GIKO=G(5,KL+1)
                  ELSE
C********** QUESTO BLOCCHETTO SERVE PER BLOCCARE IL CONVETTIVO *********
                     IGRO=KL-1
                     GICA=GIKO
                     GIKO=G(5,KL)
C***********************************************************************
                  ENDIF
                  WRITE(*,200)L2,K,GICA,KL,GIKO
                  WRITE(2,200)L2,K,GICA,KL,GIKO
                  GO TO 32
               ENDIF
            END DO
            STOP
   32       CONTINUE
C********** QUESTO BLOCCHETTO SERVE PER BLOCCARE IL CONVETTIVO *********
            IF(IBLK.EQ.1.AND.ISMOTH.EQ.0) THEN
               DO J=1,MAXNE-1
                 XTOT(J)=XTOT(J)-XXX(J,IGRO+1)*(G(5,IGRO+2)-G(5,IGRO+1))
               END DO
               DMC=DMC-(G(5,IGRO+2)-G(5,IGRO+1))
            ENDIF
C***********************************************************************
            DO KL=1,IGRO
               DO J=1,MAXNE-1
                  XXX(J,KL)=XTOT(J)/DMC
               END DO
            END DO
         ENDIF
C******** CONTROLLA PRESENZA DI SHELL CONVETTIVE ***********************
   25    CONTINUE
         LBIT=1
         IF(IBOR.EQ.1) GO TO 43
         LA=K+1
         I=0
         MAX=MAXME -1
         DO 11 K=LA,MAX
          IF(G(6,1).LT.0.D0.AND.G(6,K).LT.0.D0) GO TO 11
          IF(G(6,1).GE.0.D0.AND.G(6,K).LT.0.D0.AND.K.GT.L2) GO TO 11
          IF(G(6,1).LT.0.D0.AND.G(6,K).GE.0.D0.AND.G(4,K).GE.80.d0)THEN
            L1=K
            L2=L1+1
            GO TO 10
          END IF
          IF( G(6,K).GE.0.D0 .AND. XXX(1,K).GT.0.D0) GO TO 11
          IF(G(6,1).LT.0.D0.AND.G(6,K).GE.0.D0.
     #           AND.G(4,K).LT.80.d0)GO TO 11
          ZL = ( G( 2 , K ) + G( 2 , K+1 ) ) / 2.D0
          PP = ( G( 3 , K ) + G( 3 , K+1 ) ) / 2.D0
          TT = ( G( 4 , K ) + G( 4 , K+1 ) ) / 2.D0
          EM = ( G( 5 , K ) + G( 5 , K+1 ) ) / 2.D0
          TM = TT * 1.D+06
          PR = PP * 1.D+17 

          CALL EOS(K,1,MAXNE,1,PR,TM,XXX,'MIXING  ',
     &                       RHO,DAD,CSP,O1,G1,DEL,EMUE,3)
          CALL NKAPPA(K,1,MAXNE,1,RHO,TM,XXX,CAP)
          DRAD=CRAD*CAP*PP*ZL/(EM*(TT**4))

          IF( DRAD . LT . DAD ) GO TO 11
          IF(XXX(1,K).GT.0.D0) GO TO 11
          L1 = K
          IF(G(6,1).LT.0.d0) THEN
            L2=L1+1
          ELSE
            L2=L2+1
            L2=MAX0(L2,L1+1)
          END IF
          GO TO 10
   11    CONTINUE
   43   CONTINUE
      RETURN
  200 FORMAT(1X,'************* ELIO CENTRALE SUPERATO *************',/,
     #' MIX. FINO A:',I4,' BORDO CC:',I4,'M-CC:',1P,E12.5,' STOP A:',I4
     #,' M-CC:',1P,E12.4)
  323   FORMAT('CLASSIC BORDER R:',1p,D10.2,' M:',D10.2,' MESH:',0p,i4,
     #  1p,' OVER BORDER R:',D10.2,' M:',D10.2,' MESH:',0p,i4)
      END

      SUBROUTINE OVERSH(KCLAS,KOVER,AMBDA,VMOVER) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      COMMON/OVER/FACTOR
c      FACTOR IS DEFINED IN MAIN.F
      
      write(*,123)FACTOR
      KOVER=MAXME-1 
      P1  = G(3,KCLAS+1) 
      T1  = G(4,KCLAS+1) 
      TT1 = T1 * 1.D+06 
      PG1 = P1 * 1.D+17 
      CALL EOS(KCLAS,1,MAXNE,1,PG1,TT1,XXX,'OVERSH  ', 
     &            RHO1,DAD,CSP,O1,G1,DEL,EMUE,1) 
      GRAV=G(5,KCLAS+1)/(G(1,KCLAS+1)*G(1,KCLAS+1))
      AMBDA=OVERCOST*P1/(RHO1*GRAV)
      BOEST=G(1,KCLAS+1)+AMBDA*FACTOR 
      DO K=(KCLAS+1),(MAXME-1) 
        IF(G(1,K).GT.BOEST) THEN 
          KOVER=K-1 
          GO TO 1 
        ENDIF 
      END DO 
   1  CONTINUE 
      VMOVER=G(5,KOVER) 
      RETURN 
  123   FORMAT('LAMBDAover=',f10.5)
      END       

      SUBROUTINE UNDERSH(KCLAS,KOVER,AMBDA,VMOVER) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      FACTOR=0.5D+00
      write(*,123)FACTOR
      KOVER=MAXME-1 
      P1  = G(3,KCLAS) 
      T1  = G(4,KCLAS) 
      TT1 = T1 * 1.D+06 
      PG1 = P1 * 1.D+17 
      CALL EOS(KCLAS,1,MAXNE,1,PG1,TT1,XXX,'OVERSH  ', 
     &            RHO1,DAD,CSP,O1,G1,DEL,EMUE,1) 
      GRAV=G(5,KCLAS)/(G(1,KCLAS)*G(1,KCLAS))
      AMBDA=OVERCOST*P1/(RHO1*GRAV)
      BOEST=G(1,KCLAS)-AMBDA*FACTOR 
      write(*,'(3f15.7)')G(1,KCLAS),AMBDA,BOEST
      DO K=(KCLAS),1,-1 
        IF(G(1,K).LT.BOEST) THEN 
          KOVER=K-1 
          GO TO 1 
        ENDIF 
      END DO 
   1  CONTINUE 
      VMOVER=G(5,KOVER) 
      RETURN 
  123   FORMAT('LAMBDAunder=',f10.5)
      END              
