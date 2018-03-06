C*****************************************************************************
C  QUESTA ROUTINE CONTIENE: ENERGY - EPSI - EPSIG - EVOLUT - NCROSS -
C   INDICIZZA - RAP - READNET - PLASMA - SC - SMALLRAP.
C*****************************************************************************
C  
C   30-09-2004 ULTIMA MODIFICA "RAP" (INTRODOTTO DYDYM(MAXNE)
C   20-01-2006 cambiato burning YM=(YF+YV)/2
C
      SUBROUTINE ENERGY (T,RO,YF,IFASE,CROS,NNP)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      DIMENSION CROS(MAXPR,2),YF(MAXNE)
      DO J=1,NNP
         CROS(J,2)=QVAL(IFASE,J)
      END DO
      RETURN
      END

      SUBROUTINE EPSI (RO,T,JF,EPSN,EPSPP,EPSCNO,EPSA,ECAR,IDER)
c ********* INPUT / OUTPUT ************
c * T RO temperatura & densita'       *
c * YF abbondanze in numero           *
c * EPSN Epsilon nucleare             *
c * IMPORTANTE!!!! UNITA' CGS         *
c *************************************
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      INCLUDE 'consts.2p1'
      DIMENSION CROS(MAXPR,2),YF(MAXNE)
      COMMON/FLUX/FLUSSO(20)
      DO J=1,MAXNE
         YF(J)=XXX(J,JF)/ATW(J)
      END DO
      EPSPP=0.D0
      EPSCNO=0.D0
      EPSA=0.D0
      ECAR=0.D0
      EPSN=0.D0
c######### SKIPPA FASI INUTILI E RICERCA FASE EVOLUTIVA ###############
      CALL SKIPFASE (T,RO,YF,JF,IFASE)
      IF(IFASE.EQ.0) RETURN
C #####################################################################
      NNP = NP(IFASE)
      CALL NCROSS (RO,T,YF,IFASE,NNP,CROS,IDER)
      CALL ENERGY (T,RO,YF,IFASE,CROS,NNP)
      DO J=1,NNP
         EPSN=EPSN+YF(ICA(IFASE,NET(IFASE,J,1)))*YF(ICA(IFASE,
     $   NET(IFASE,J,2)))*CROS(J,1)*CROS(J,2)/DK(IFASE,J)
      FLUSSO(J)=FLUSSO(J)+YF(ICA(IFASE,NET(IFASE,J,1)))*YF(ICA(IFASE,
     $   NET(IFASE,J,2)))*CROS(J,1)*CROS(J,2)/DK(IFASE,J)*ERGEV*NAvo*RO
      IF(IFASE.EQ.1.AND.J.EQ.7) EPSPP=EPSN*ERGEV*NAvo*RO
      FLUSSO(8)=FLUSSO(8)+EPSPP
      END DO
      EPSN=EPSN*ERGEV*NAvo*RO
      IF(IFASE.EQ.1) THEN
         EPSCNO=EPSN-EPSPP
      ELSE IF(IFASE.EQ.2) THEN
         EPSA=EPSN
      ELSE IF(IFASE.EQ.3) THEN
         ECAR=EPSN
      ENDIF
      RETURN
      END   
   
      SUBROUTINE EPSIG(P,T,CP,GRAVI,DAD,PV,TV,JF) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
C======================================================================
C==== CALCOLO MODELLI DI ZAHB =========================================
      IF(ISTART.EQ.5) THEN
        GRAVI=0.D0
        RETURN
      ENDIF      
C======================================================================    
      IF(NMD.GT.80.OR.ISTART.EQ.2) THEN 
         DP=P-PV 
         DT=T-TV 
         PPUNTO=DP*T*DAD/P 
         TPUNTO=-DT 
        GRAVI=((PPUNTO+TPUNTO)*CP)/(SECANNI*HT1*1.D-6)
      ENDIF
c      if(xxx(3,1).lt.0.15d0.and.xxx(3,1).gt.0.d0)then
c        jlim=1
c       do jj=2,maxme-1
c         IF( XXX(1,JJ-1).LE.0.D0.AND.XXX(1,JJ).GT.0.D0 )jlim=jj
c       enddo
c        if(jf.le.jlim)then
c          GRAVI=0.d0
c       endif
c      endif
      RETURN 
      END 

      SUBROUTINE NCROSS (RO,T,YF,IFA,NNP,CROS,IDER)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      REAL    *8 SC,CROS(MAXPR,2),YF(MAXNE),CROSER(MAXPR,2),
     %           Z1,Z2,RO,T
      INTEGER *4 J,IFA,NNP
C*********************************************************************      
      CALL PLASMA (RO,T,YF)             !   CARATTERISTICHE DEL PLASMA
      IF(IDER.NE.2)THEN
         CALL SIGMANEW (IFA,NNP,T,RO,YF,CROSER)  !<NA*SIGMA*V>
      END IF
      DO J=1,NNP
         IF(IFA.EQ.2.AND.J.EQ.1)THEN
            CROS(J,1)=SC(RO,T,YF,2.D0,2.D0,2,100)*
     %                SC(RO,T,YF,2.D0,4.D0,2,100)*CROSER(J,1)
         ELSE
C#####################################################################
            IF(ICA(IFA,NET(IFA,J,2)).EQ.MAXNE) THEN
               CROS(J,1)=CROSER(J,1)   !cattura elettronica
            ELSE
C#####################################################################            
               Z1=ZET(ICA(IFA,NET(IFA,J,1)))
               Z2=ZET(ICA(IFA,NET(IFA,J,2)))
               CROS(J,1)=SC(RO,T,YF,Z1,Z2,IFA,J)*CROSER(J,1)
            ENDIF   
         END IF
         if(cros(j,1).lt.0.d0) then
            write(*,*)'ATTENZIONE! SEZIONE D`URTO NEGATIVA'
            write(*,*)cros(j,1),SC(RO,T,YF,Z1,Z2,IFA,J),croser(j,1),j
            stop
         end if
      END DO
      RETURN
      END

      SUBROUTINE EVOLUT
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      INCLUDE 'consts.2p1'
      dimension xymax(maxme),ZMETAT(maxme)
C###################################################################      
      LOGICAL*1 CONVETTA
C###################################################################
      DIMENSION YF(MAXNE)
c      DATA SECANNI/3.15587D+07/
      PAS=HT1*SECANNI
C ************ MAIN LOOP *************
      DO 1 K=1,MAXME-1
         if(k.eq.1) then
           do ll=1,maxne-1
              xnm(ll)=xxx(ll,1)
           end do
        end if
        xymax(k)=0.d0
        ZMETAT(k)=0.0d0
         do jj=4,maxne-1
           ZMETAT(k)=ZMETAT(k)+xxx(jj,k)
         end do
         xymax(k)=1.0d0-zmetat(k)
         KEQUIL(K)=0
         DO J=1,MAXNE
            YF(J)=XXX(J,K)/ATW(J)
         END DO
C ************** CALCOLA RO **********
         T=1.0D+06*(G(4,K)+G(4,K+1))/2.D0
         P=1.0D+17*(G(3,K)+G(3,K+1))/2.D0
c         T=1.0D+06*G(4,K)
c         P=1.0D+17*G(3,K)
         CALL EOS(K,1,MAXNE,1,P,T,XXX,'EVOLUT  ',RO,O1,O2,O3,
     #G1,DEL,EMUE,1)
C#####RICERCA FASE EVOLUTIVA E SKIPPA FASI INUTILI###################         
         CALL SKIPFASE (T,RO,YF,K,IFASE)
         IF(IFASE.EQ.0) GO TO 1
C####################################################################         
         NNE=NE(IFASE)
         NNP=NP(IFASE)
C ******************** CONTROLLO INTEGRAZIONE: **********************
C ** CHIMICA MIXED -> SERIE DI TAYLOR   CHIMICA LOCALE -> RAPNEW ****
C **************  SOLO SE CENTRALE ****************************
C####################################################################
         CONVETTA=.FALSE.
         IF(KBORD.GT.0)THEN
            IF(XXX(1,1).GT.0.D0.OR.XXX(3,1).GT.0.D0)THEN
               IF(K.LE.KBORD)CONVETTA=.TRUE.  !VERIFICARE SE K=KBORD CONVETTIVO
            ENDIF
         ENDIF
C####################################################################
         CALL RAP(RO,T,PAS,IFASE,NNE,NNP,YF,K,CONVETTA)
         DO J=1,MAXNE
            XXX(J,K)=YF(J)*ATW(J)
         END DO
         IF(IFASE. EQ. 3) THEN
            CARXY(1,K)=XXX(1,K)
            CARXY(2,K)=XXX(3,K)
            XXX(1,K)=0.D0
            XXX(3,K)=0.D0
         ENDIF
         if(k.eq.1) then
           do ll=1,maxne-1
              xnm(ll)=xxx(ll,1)
           end do
        end if
        IF(IFASE.EQ.1)THEN
            XTEST=XXX(1,K)+XXX(3,K)
            IF(XTEST.GT.XYMAX(K))XXX(3,K)=XYMAX(K)-XXX(1,K)
         ENDIF
    1 CONTINUE
      RETURN
      END

      SUBROUTINE INDICIZZA(PROC)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      CHARACTER*1 TESTATA(72)
      CHARACTER*4 namel(maxne),nomi(3,maxne),proc(3,maxpr,4)
      CHARACTER*16 REAZIONI(36)
      character*16 STRINGA
      LOGICAL *1 DEBUG
      DATA DEBUG/.false./
      DATA REAZIONI/
!*******************  CATTURE PROTONICHE  ************************
     #   'H   H   g   D   ','D   H   g   He3 ','Li7 H   He4 He4 ',
     #   'Be7 H   He4 He4 ','C12 H   g   C13 ','C13 H   g   N14 ',
     #   'N14 H   g   N15 ','N15 H   g   O16 ','N15 H   He4 C12 ',
     #   'O16 H   g   O17 ','O17 H   He4 N14 ','Ne20H   g   Ne21',
     #   'Ne21H   g   Ne22','Ne22H   g   Na23','Na23H   He4 Ne20',
     #   'Na23H   g   Mg24','Mg24H   g   Mg25','Mg25H   g   Mg26',
     #   'Mg26H   g   Al27','Al27H   g   Si28',
!**********************  CATTURE ALPHA  **************************
     #   'He3 He4 g   Be7 ','C12 He4 g   O16 ','O16 He4 g   Ne20',
     #   'N14 He4 g   O18 ','O18 He4 g   Ne22','Ne20He4 g   Mg24',
     #   'Ne22He4 n   Mg25','Ne22He4 g   Mg26','Mg24He4 g   Si28',
!***********************  CORPI IDENTICI  ************************
     #   'He3 He3 He4 H   ','C12 C12 He4 Ne20','C12 C12 H   Na23',
     #   'O16 O16 He4 Si28',
!*************************  3ALPHA  ******************************
     #   'He4 Dumyg   C12 ',
!**********************    DECADIMENTI  **************************
     #   'Be7 bet-g   Li7 ',
!*******************   FOTODISINTEGRAZIONI   *********************
     #   'Ne20DumyHe4 O16 '/
!*****************************************************************
  100 FORMAT(2X,'***** ATTENZIONE!  PROCESSO ',I4,2X,
     #A16,' IN FASE ',I2,' NON TROVATO!  ********')
!*****************************************************************
      DO I=1,3
         DO J=1,NP(I)
            STRINGA=PROC(I,J,1)//PROC(I,J,2)
     #           //PROC(I,J,3)//PROC(I,J,4)
            KINI=0
            KFIN=0
            IF(PROC(I,J,2).EQ.'bet-') THEN
               KINI=35
               KFIN=35
            ELSEIF(PROC(I,J,2).EQ.'Dumy') THEN
              IF(I.EQ.2) THEN
                 KINI=34
                 KFIN=34
              ELSEIF(I.EQ.3)THEN
                 KINI=36
                 KFIN=36
              ELSE
                 WRITE(*,100)J,STRINGA,I
                 STOP
              END IF
            ELSEIF(PROC(I,J,2).EQ.'H   ') THEN
               KINI=1
               KFIN=20
            ELSEIF(PROC(I,J,2).EQ.'He4 ') THEN
               KINI=21
               KFIN=29
            ELSEIF(PROC(I,J,1).EQ.PROC(I,J,2)) THEN
               KINI=30
               KFIN=33
            ENDIF
            KK=KINI
            DO WHILE(STRINGA.NE.REAZIONI(KK).AND.KK.LE.KFIN)
               KK=KK+1
            END DO
            IF(KK.GT.KFIN) THEN
               WRITE(*,100)J,STRINGA,I
               STOP
            END IF
            ICROS(I,J)=KK
            IF(DEBUG)WRITE(*,'(I4,1X,A16,2X,A16)')ICROS(I,J),
     #             REAZIONI(KK),STRINGA
          END DO
       END DO
       RETURN
       END            

      SUBROUTINE RAP(RO,T,PAS,IFASE,NNE,NNP,YF,JF,CONVETTA)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      LOGICAL*1 CONVETTA
      DIMENSION  YF(MAXNE),YM(MAXNE),YV(MAXNE),CROS(MAXPR,2)
      DIMENSION  B(MAXNE),A(MAXNE*MAXNE),DXX(MAXNE,MAXNE),DX(MAXNE),
     #           DYDYM(MAXNE)
      NI=0
      BETA=1.D0
      GAMMA=1.D0
      ESCI=1.D0
      NEL=NNE
      CALL NCROSS(RO,T,YF,IFASE,NNP,CROS,1)
      DO I=1,NNE
         YV(ICA(IFASE,I))=YF(ICA(IFASE,I))
        DYDYM(ICA(IFASE,I))=1.D0

c=== test problem R
         IF(IFASE.EQ.1) THEN
          if(ica(ifase,I).eq.1) DYDYM(ICA(IFASE,I))=0.5D0
         ELSEIF(IFASE.EQ.2) THEN
          if(ica(ifase,I).eq.3) DYDYM(ICA(IFASE,I))=0.5D0
         ENDIF
          
c===

      END DO
      YM(MAXNE)=YF(MAXNE)
C######################################################################      
      IF(IFASE.EQ.1)THEN
         I=1
         DO WHILE (ICA(IFASE,I).NE.3.AND.I.LE.NNE)
           I=I+1 
         ENDDO
         ICONSER=I
      ELSEIF(IFASE.EQ.2)THEN
         I=1
         DO WHILE (ICA(IFASE,I).NE.6.AND.I.LE.NNE)
           I=I+1 
         ENDDO
         ICONSER=I   
      ENDIF
      SUMOLD=0.D0
      DO I=1,MAXNE-1
         SUMOLD=SUMOLD+YF(I)*ATW(I)
      ENDDO
C#####################################################################                 
c ***************** MAIN LOOP *************
      desci=1.D-8
      DO WHILE (ESCI. GE. DESCI)
        if(ni.ge.20)then
          desci=1.D-7
c          write(*,*)'ni maggiore di 20'
c         pause
        endif
        if(ni.ge.100)then
          desci=1.D-6
c          write(*,*)'ni maggiore di 100'
c         pause
        endif
        if(ni.ge.50)then
          desci=1.D-5
c          write(*,*)'ni maggiore di 500'
c         pause
        endif
        DO I=1,NNE
           DX(I)=0.D0
           DO K=1,NNE
             DXX(I,K)=0.D0
           END DO
C#########################################################           
           YM(ICA(IFASE,I))=YF(ICA(IFASE,I))
          if(ifase.eq.1) then  
             if(ica(ifase,I).eq.1) then
          YM(ICA(IFASE,I))=(YF(ICA(IFASE,I))+YV(ICA(IFASE,I)))/2.d0
              endif
          elseif(ifase.eq.2) then
             if(ica(ifase,I).eq.3) then
          YM(ICA(IFASE,I))=(YF(ICA(IFASE,I))+YV(ICA(IFASE,I)))/2.d0
             endif          
          endif   
        END DO
        IF(CONVETTA)THEN
           IF(IFASE.EQ.1.AND.XXV(1,1).gT.1.D-6)THEN
              YM(1)=YV(1)
             DYDYM(1)=0.D0
           ELSEIF(IFASE.EQ.2.AND.XXV(3,1).GT.1.D-6)THEN
              YM(3)=YV(3)
             DYDYM(3)=0.D0
           ENDIF
        ENDIF
        IF(IFASE. EQ. 2) YM(MAXNE)=YM(3)**2*RO/6.D0
        DO J=1,NNP
           NET1=NET(IFASE,J,1)
           NET2=NET(IFASE,J,2)
           NET3=NET(IFASE,J,3)
           NET4=NET(IFASE,J,4)
c ******* CALCOLO DERIVATE DY/DT *********
           CROSSO=RO*CROS(J,1)/DK(IFASE,J)
           YM1=YM(ICA(IFASE,NET1))*YM(ICA(IFASE,NET2))*CROSSO
           DX(NET1)=DX(NET1)-YM1*CC(IFASE,J,1)
           IF(NET2.LE.NNE) DX(NET2)=DX(NET2)-YM1*CC(IFASE,J,2)
           DX(NET3)=DX(NET3)+YM1*CC(IFASE,J,3)
           IF(NET4.NE.0) DX(NET4)=DX(NET4)+YM1*CC(IFASE,J,4)
        END DO
        DO J=1,NNP
           NET1=NET(IFASE,J,1)
           NET2=NET(IFASE,J,2)
           NET3=NET(IFASE,J,3)
           NET4=NET(IFASE,J,4)
c ******* CALCOLO DERIVATE DF/DY *********
           CROSSO=RO*CROS(J,1)/DK(IFASE,J)
           YM1=YM(ICA(IFASE,NET1))*CROSSO
           YM2=YM(ICA(IFASE,NET2))*CROSSO
           DXX(NET1,NET1)=DXX(NET1,NET1)+YM2*CC(IFASE,J,1)*DYDYM(NET1)
           DXX(NET3,NET1)=DXX(NET3,NET1)-YM2*CC(IFASE,J,3)*DYDYM(NET1)
           IF(NET2.LE.NEL) THEN
             DXX(NET1,NET2)=DXX(NET1,NET2)+YM1*CC(IFASE,J,1)*DYDYM(NET2)
             DXX(NET2,NET2)=DXX(NET2,NET2)+YM1*CC(IFASE,J,2)*DYDYM(NET2)
             DXX(NET3,NET2)=DXX(NET3,NET2)-YM1*CC(IFASE,J,3)*DYDYM(NET2)
             DXX(NET2,NET1)=DXX(NET2,NET1)+YM2*CC(IFASE,J,2)*DYDYM(NET1)
             IF(NET4.NE.0) DXX(NET4,NET2)=DXX(NET4,NET2)-
     #            YM1*CC(IFASE,J,4)*DYDYM(NET2)
           END IF
           IF(NET4.NE.0) DXX(NET4,NET1)=DXX(NET4,NET1)-
     #            YM2*CC(IFASE,J,4)*DYDYM(NET1)
        END DO
        DO I=1,NEL
           DXX(I,I)=1.D0/PAS+DXX(I,I)
        END DO  
c **** RISOL. SET EQUAZIONI LINEARI DF=-B ****
        DO I=1,NEL
          B(I)=-(YF(ICA(IFASE,I))-YV(ICA(IFASE,I)))/PAS+DX(I)
        END DO
C##################################################################        
        DO I=1,NNE
           DXX(ICONSER,I)=ATW(ICA(IFASE,I))
        ENDDO
        SUMNEW=0.D0
        DO I=1,MAXNE-1
           SUMNEW=SUMNEW+YF(I)*ATW(I)
        ENDDO   
        B(ICONSER)=-(SUMNEW-SUMOLD)
C##################################################################        
        DO IE=1,NEL
          DO II=1,NEL
            IA=(IE-1)*NEL+II
            A(IA)=DXX(II,IE)
          END DO
        END DO
        CALL KERNEL (A,B,NEL,KS)
        IF(KS. EQ. 1) THEN
          WRITE(*,130)
          STOP
        ENDIF
        DO I=1,NEL
           YF(ICA(IFASE,I))=YF(ICA(IFASE,I))+B(I)
        END DO
        NI=NI+1
c ** VARIAZIONI PERCENTUALI E VERIFICA EQUAZ. **
        IF(NI. NE. 1) THEN
          BETA=0.D0
          GAMMA=0.D0
          DO I=1,NEL
            IF(YV(ICA(IFASE,I)).GE.1.D-15) THEN
C######################################################            
               IF(I.EQ.ICONSER)THEN
                  SUMNEW=0.D0
                  DO II=1,MAXNE-1
                     SUMNEW=SUMNEW+YF(II)*ATW(II)
                  ENDDO
                  ERR1=DABS((SUMNEW-SUMOLD)/SUMOLD)
               ELSE
                  DELTA=DABS((YF(ICA(IFASE,I))-YV(ICA(IFASE,I)))/
     #             YV(ICA(IFASE,I)))
                  IF(DX(I).NE.0.D0.AND.DELTA.GT.1.D-20) then !08)THEN
                  ERR1=DABS(1.D0-(YF(ICA(IFASE,I))-YV(ICA(IFASE,I)))/
     #             (PAS*DX(I)))
                  ELSE 
                  ERR1=0.D0
                  ENDIF
               ENDIF   
C######################################################                  
              ERR2=DABS(B(I)/YF(ICA(IFASE,I)))
              IF(ERR1. GT. BETA) THEN
                 IB=I
                 BETA=ERR1
              ENDIF
              IF(ERR2. GT. GAMMA) THEN
                 IG=I
                 GAMMA=ERR2
              ENDIF
            ENDIF
          END DO
          ESCI=DMIN1(BETA,GAMMA)
          IF(NI. GE. 1000) THEN
             WRITE(*,444)IB,BETA,IG,GAMMA,NI,JF
             WRITE(2,444)IB,BETA,IG,GAMMA,NI,JF
             WRITE(*,250)
             WRITE(*,113)T,RO,JF
             WRITE(*,111)(YF(ica(IFASE,K)),K=1,6)
             WRITE(*,111)(Yv(ica(IFASE,K)),K=1,6)
  111        FORMAT(1P,6E9.2)
             WRITE(2,250)
             WRITE(2,113)T,RO,JF
             WRITE(2,111)(YF(K),K=1,6)
             STOP
          ENDIF
        ENDIF
      END DO
      RETURN
  113 FORMAT(1X,'T =',1P,E10.3,2X,'RO =',1P,E10.3,2X,'MESH =',I4)
  130 FORMAT(1X,'ATTENZIONE QUALCOSA NON VA NELLA EPSI. IL DET = 0 NEL
     #           RAPHSON/NEWTON')
  250 FORMAT(1X,'MORTO: ITERA TROPPO NEL RAPHSON/NEWTON ')
  444 FORMAT(1X,'ITERAZIONI EPSI: DF =',I3,E10.3,' DX =',I3,E10.3,
     #          ' N. ITER.=',I3,' N. MESH=',I4)
      END




      SUBROUTINE READNET
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      DIMENSION LCC(3,MAXPR,4)
      CHARACTER*1 TESTATA(72),TMIST(40)
      CHARACTER*4 namel(maxne),nomi(3,maxne),proc(3,maxpr,4)
      LOGICAL *1 DEBUG,SPANNA
      DATA DEBUG/.false./
C************ LETTURA PROCESSI E NETWORK **************
      OPEN (19,FILE='network',STATUS='OLD')
      READ(19,101)TESTATA
      WRITE(*,101)TESTATA
      READ(19,101)TESTATA
      WRITE(*,101)TESTATA
      WRITE(10,153)TESTATA
  101 FORMAT(72A1)
  151 FORMAT(25X,40A1)
  152 FORMAT('Mixture: ',40A1)
  153 FORMAT('Nuclear Cross Sections: ',72A1)
  102 FORMAT(4I4,4F4.0,F8.3)
  103 FORMAT(10F5.1)
  110 FORMAT(3X,'********  EL:',I4,' IN FASE:',I2,1X,A4,
     #         ' NON TROVATO  ********')
  120 FORMAT(3X,'********  EL:',I4,1X,A4,' IN FASE:',I2,
     #         ' NON EVOLUTO  ********')
  130 FORMAT(3X,'********  EL:',I4,1X,A4,' IN PROC.:',I2,
     #         ' IN FASE:',I4,' NON TROVATO  ********')
      DO I=1,3
         READ(19,'(I4,66A1)')NE(I),(TESTATA(LL),LL=1,66)
         IF(DEBUG)WRITE(*,'(I4,66A1)')NE(I),(TESTATA(LL),LL=1,66)
         DO J=1,NE(I)+1
            READ(19,'(5X,A4)')NOMI(I,J)
            IF(DEBUG)WRITE(*,'(I4,1X,A4)')J,NOMI(I,J)
         END DO
         READ(19,'(I4,66A1)')NP(I),(TESTATA(LL),LL=1,66)
         IF(DEBUG)WRITE(*,'(I4,66A1)')NP(I),(TESTATA(LL),LL=1,66)
         DO J=1,NP(I)
            READ(19,'(6X,I1,A4,1X,I1,A4,1X,I1,A4,1X,I1,
     #         A4,2X,F6.3)')(LCC(I,J,LL),PROC(I,J,LL),LL=1,4),QVAL(I,J)
            DO LL=1,4
               CC(I,J,LL)=DFLOAT(LCC(I,J,LL))
            END DO
            IF(DEBUG)
     #      WRITE(*,'(I4,2X,I1,A4,1X,I1,A4,1X,I1,A4,1X,
     #            I1,A4,2X,F6.3)')J,
     #            (LCC(I,J,LL),PROC(I,J,LL),LL=1,4),QVAL(I,J)
         END DO
      END DO
      WRITE(16,170)(NOMI(1,J),J=1,NE(1))
      WRITE(66,170)(NOMI(1,J),J=1,NE(1))
      WRITE(56,156)(NOMI(1,J),J=1,NE(1))
  170 format('#NMD',8X,'LOG(T)',9X,'LOG(L/LSUN)',4X,'LOG(TE)',5X,
     # 'LOG(R/RSUN)',4X,'LOG(GI)',8X,'MCC',7X,'LOG(RHOc)',5X,
     # 'LOG(Tc)',8X,'Mec',7X,'LOG(ROec)',5X,'LOG(Tec)',2X,
     #  12(5X,A4,4X))
  156 format('#NMD',8X,'LOG(T)',12X,'dT',12X,'tHe3',12X,'t_smr',
     #  6X,'LOG(L/LSUN)',4X,'LOG(TE)',4X,
     # 'LOG(R/RSUN)',4X,'LOG(GI)',8X,'MCC',7X,'LOG(RHOc)',5X,
     # 'LOG(Tc)',8X,'Mec',7X,'LOG(ROec)',5X,'LOG(Tec)',2X,
     #  12(5X,A4,4X))
C#####################################################################
C il 12 qui sopra e` il # di elementi coinvolti nella reazione,
C cambiando il network va modificato. 
C##################################################################### 
      READ(19,151)TMIST
      WRITE(*,152)TMIST
      WRITE(10,152)TMIST
      READ(19,101)TESTATA
      DO JV=1,MAXNE
         ABSOL(JV)=0.D0
      END DO
      DO I=1,MAXNE
         READ(19,'(6X,A4,2X,F7.4,D8.1,2X,D12.5)')
     #         NAMEL(I),ATW(I),ZET(I),ABSOL(I)
      END DO
      CLOSE(19)
      DO I=1,3
         DO J=1,NE(I)+1
            KK=1
            DO WHILE(NAMEL(KK).NE.NOMI(I,J).AND.KK.LE.MAXNE)
               KK=KK+1
            END DO
            IF(KK.GT.MAXNE) THEN
               WRITE(*,110)J,I,NOMI(I,J)
               STOP
            END IF
            SPANNA=.TRUE.
            JJ=1
            DO WHILE (SPANNA)
               LL=1
               DO WHILE (PROC(I,JJ,LL).NE.NOMI(I,J).AND.
     #                 LL.LE.4)
               LL=LL+1
              if(LL.eq.5) then
               ll=4
              go to 131
              endif
               END DO
 131           continue              
C=================================================

               IF(LL.LE.4)SPANNA=.FALSE.               
               JJ=JJ+1
               IF(JJ.GT.NP(I))SPANNA=.FALSE.
            END DO
            IF(JJ.GT.NP(I).AND.LL.GT.4.AND.J.LE.NE(I)) THEN
               WRITE(*,120)J,NOMI(I,J),I 
               STOP
            END IF  
            ICA(I,J)=KK
         END DO
         DO J=1,NP(I)
            DO LL=1,4
               IF(PROC(I,J,LL).EQ.'g   ')THEN
               ELSEIF(PROC(I,J,LL).EQ.'bet-')THEN
               ELSEIF(PROC(I,J,LL).EQ.'bet+')THEN
               ELSEIF(PROC(I,J,LL).EQ.'Dumy')THEN
               ELSE   
                  SPANNA=.TRUE.
                  JJ=1
                  DO WHILE(SPANNA)
                     IF(PROC(I,J,LL).EQ.NOMI(I,JJ))SPANNA=.FALSE.
                     JJ=JJ+1
                     IF(JJ.GT.NE(I)+1)SPANNA=.FALSE.
                  END DO
                  IF(JJ.GT.NE(I)+1) THEN
                     WRITE(*,130)LL,PROC(I,J,LL),J,I
                     STOP
                  END IF
               END IF
            END DO
         END DO
      END DO
      CALL INDICIZZA(PROC)
      DO I=1,3
         DO J=1,NP(I)
            DO LL=1,4
               IF(PROC(I,J,LL).EQ.'bet+ ')PROC(I,J,LL)='Dumy'
               IF(PROC(I,J,LL).EQ.'bet- ')PROC(I,J,LL)='Dumy'
               IF(PROC(I,J,LL).EQ.'g    ') THEN
                  NET(I,J,LL)=0
               ELSE
                  KK=1
                  DO WHILE(NAMEL(KK).NE.PROC(I,J,LL).AND.KK.LE.MAXNE)
                     KK=KK+1
                  END DO
                  IF(KK.GT.MAXNE) THEN
                     WRITE(*,110)J,I,PROC(I,J,LL)
                     STOP
                  END IF
                  KK=1
                  DO WHILE(PROC(I,J,LL).NE.NOMI(I,KK).AND.KK.LE.NE(I))
                     KK=KK+1
                  END DO
                  NET(I,J,LL)=KK
               END IF
            END DO           
            IF(DEBUG) THEN
            IF(NET(I,J,3).EQ.0) THEN
            WRITE(*,'(I4,2X,I1,A4,1X,I1,A4,1X,I1,A4,1X,I1,A4,2X,F6.3)')
     #            J,(LCC(I,J,LL),NAMEL(ICA(I,NET(I,J,LL))),LL=1,2),
     #            LCC(I,J,3),'g   ',LCC(I,J,4),
     #             NAMEL(ICA(I,NET(I,J,4))),QVAL(I,J)
            ELSE
            WRITE(*,'(I4,2X,I1,A4,1X,I1,A4,1X,I1,A4,1X,I1,A4,2X,F6.3)')
     #      J,(LCC(I,J,LL),NAMEL(ICA(I,NET(I,J,LL))),LL=1,4),QVAL(I,J)
            ENDIF
            ENDIF
            IF(NET(I,J,3).EQ.0) THEN
               NET(I,J,3)=NET(I,J,4)
               NET(I,J,4)=0
            END IF
            IF(NET(I,J,1). EQ. NET(I,J,2)) THEN
                DK(I,J)=2
            ELSE
                DK(I,J)=1
            ENDIF
         END DO
      END DO
      DO IFA=1,3          !FATTORI UTILIZZATI NELLA SCREENING
         DO J=1,NP(IFA)   !ATTENZIONE AL CASO DELLE 3 ALFA - VEDI NCROSS
            Z1=ZET(ICA(IFA,NET(IFA,J,1)))
            Z2=ZET(ICA(IFA,NET(IFA,J,2)))
            PASSIO(IFA,J) = (Z1+Z2)**(1.86)-Z1**(1.86)-Z2**(1.86)
            PASSI(1,IFA,J)=(Z1+Z2)**(5.D0/3.D0)-
     #                      Z1**(5.D0/3.D0)-Z2**(5.D0/3.D0)
            PASSI(2,IFA,J)=(Z1+Z2)**(4.D0/3.D0)-
     #                      Z1**(4.D0/3.D0)-Z2**(4.D0/3.D0)
            PASSI(3,IFA,J)=(Z1+Z2)**(2.D0/3.D0)-
     #                      Z1**(2.D0/3.D0)-Z2**(2.D0/3.D0)
         END DO
      END DO
      RETURN
      END       

      SUBROUTINE PLASMA(RHO,T,YF) 
      IMPLICIT REAL*8(A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      DIMENSION YF(MAXNE),F(MAXNE) 
!  *****************   CALCOLO PESO MOLECOLARE   *************** 
      SUM=0.D0 
      DO K=1,MAXNE-1 
         SUM=SUM+YF(K) 
      END DO 
      PEM=1.D0/SUM 
!  *****************   CALCOLO DI TETA   *********************** 
      FER=5.4885D+07*RHO*(1.D0+YF(1)*ATW(1))/(T*DSQRT(T)) 
      IF(FER.LE.3.D0) THE=1.D0/(1.D0+.39716D0*FER-.00929D0*FER*FER) 
      IF(FER.GT.3.D0) 
     %   THE=1.1447D0/((FER**(2.D0/3.D0))*(1.D0+.28077D0/FER)) !Teta 
!  ***********   CALCOLO ZETA MEDIO E ZETA TILDE  ************** 
      ZM=0.D0 
      ZT=0.D0 
      DO K=1,MAXNE-1 
         F(K)=YF(K)*PEM                 ! ni/nI 
         ZM=ZM+ZET(K)*F(K)              ! z Medio 
         ZT=ZT+ZET(K)*F(K)*(ZET(K)+THE) 
      END DO 
      ZT=DSQRT(ZT)                      ! z Tilde 
      ZT058=ZT**(0.58)                  ! QUANTITA' CHE USA NELLA SC 
      ZM028=ZM**(0.28) 
      ZM13=ZM**(1.D0/3.D0) 
!  ***********   CALCOLO  < z**(3b-1) >  ********************** 
      ZAV=0.D0 
      DO K=1,MAXNE-1 
         ZAV=ZAV+F(K)*(ZET(K)**(1.58))         ! < z**(3b-1) > 
      END DO 
      RETURN 
      END 
      FUNCTION SC(RHO,T,YF,Z1,Z2,IFASE,IPRO)
      IMPLICIT REAL *8(A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      INCLUDE 'consts.2p1'
      REAL    *8 RHO,T,YF(MAXNE),SC,SUM,PEM,ZM,ZT,
     %           FL,ZAV,FL12,SCR1,SCR2,OVER
      REAL    *8 HTAG,A1,A2,G12,MU12,T12,CHK,Z1,Z2,INFSTR
      DATA OVER/200.D0/
      DATA INFSTR/0.2D0/
C **********************************************************************
C  WEAK - INTERMEDIATE AND INTERMEDIATE-STRONG SCREENING:
C     GRABOSKE, DE WITT, GROSSMAN AND COOPER, AP. J., 181, 457-474, 1973
C     DEWITT, GRABOSKE AND COOPER, AP. J., 181, 439-456, 1973
C  STRONG SCREENING:
C     ITOH,TOTSUJI AND ICHIMARU, AP. J., 218, 477-483, 1977
C     ITOH, TOTSUJI, ICHIMARU AND DE WITT, AP. J., 234, 1079-1084 , 1979
C  ********************   TOPPA PER IL BRANCHING DEL BE7  *************
      HTAG=HPLANCK/PI
      
      IF(IPRO.EQ.100) THEN
                PASS0 = (Z1+Z2)**(1.86)-Z1**(1.86)-Z2**(1.86)
                PASS1 = (Z1+Z2)**(5.D0/3.D0)
     %                  -Z1**(5.D0/3.D0)-Z2**(5.D0/3.D0)
                PASS2 = (Z1+Z2)**(4.D0/3.D0)
     %                  -Z1**(4.D0/3.D0)-Z2**(4.D0/3.D0)
                PASS3 = ( (Z1+Z2)**(2.D0/3.D0)
     %                  -Z1**(2.D0/3.D0)-Z2**(2.D0/3.D0) )
      ELSE
                PASS0 = PASSIO(IFASE,IPRO)
                PASS1 = PASSI(1,IFASE,IPRO)
                PASS2 = PASSI(2,IFASE,IPRO)
                PASS3 = PASSI(3,IFASE,IPRO)
      ENDIF
C  *******************************************************************
!  *************    CALCOLO LAMBDA ZERO    ********************
      FL=1.88D+08*DSQRT(RHO/(PEM*T**3))
!  ************   CALCOLO LAMBDA 12  *************************
      FL12=FL*Z1*Z2*ZT           ! Lambda 12
!  **********     CALCOLO SCREENING  *************************
      IF(FL12.LE..1D0) THEN
          SC=DEXP(ZT*FL*Z1*Z2)              ! WEAK
      ELSE
          IF(FL12.LE.5.D0) THEN
             FLV=DLOG(FL)
             FL086=DEXP(0.86*FLV)
             SCR1 = 0.38D0*ZAV/(ZT058*ZM028)*PASS0*FL086
             IF(FL12.LE.2.D0) THEN
                SC=DEXP(SCR1)               ! INTERMEDIATE
             ELSE
                FL23=DEXP((2.D0/3.D0)*FLV)
                SCR2 = 0.624D0*ZM13*FL23*( PASS1+
     %                 0.316D0*ZM13*PASS2+
     %                 0.737D0/ZM*PASS3/FL23 )
                IF(SCR2.GT.OVER)SCR2=OVER
                SC = DMIN1(SCR1,SCR2)      ! INTERMEDIATE-STRONG
                SC = DEXP(SC)
             ENDIF
          ELSE
             SUM=0.D0
             DO K=1,MAXNE-1
                SUM = SUM + ZET(K)*YF(K)
             END DO
            A1   = ((3.D0*Z1)/(4.D0*PI*RHO/Umass*SUM))**(1.D0/3.D0)
            A2   = ((3.D0*Z2)/(4.D0*PI*RHO/Umass*SUM))**(1.D0/3.D0)
            G12  = (Z1*Z2*ECharg**2)/(0.5D0*(A1+A2)*KB*T)
            MU12 = (2.D0*Z1*Z2)/(Z1+Z2)*Umass 
             T12 =((27.D0*PI**2/4.D0)*(2.D0*MU12*Z1**2*Z2**2*ECharg**4)/
     %              (KB*T*HTAG**2))**(1.D0/3.D0)
             CHK  = 3.D0*G12/T12
             IF(CHK.GE.INFSTR) THEN            ! STRONG (ITHO)
c                 write(*,*)'strong screening itoh'
                 SC = DEXP(1.25D0*G12-0.095D0*T12*(3.D0*G12/T12)**2)
             ELSE
c                 write(*,*)'strong screening GRABOSKE'
                 FL23=FL**(2.D0/3.D0)
                 SCR2 = 0.624D0*ZM13*FL23*( PASS1+
     %                  0.316D0*ZM13*PASS2+
     %                  0.737D0/ZM*PASS3/FL23 )
                 IF(SCR2.GT.OVER)then
                    write(*,*)scr2,'SCR2 MAGGIORE DI OVER'
                    SCR2=OVER
                 endif
                 SC = DEXP(SCR2)               ! STRONG (GRABOSKE)
             END IF
         END IF
      END IF
      RETURN
      END

      SUBROUTINE SMALLRAP(NCONV,JBCONV)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      INCLUDE 'consts.2p1'
C###################################################################      
      LOGICAL*1 CONVETTA
      LOGICAL*1 EQUIHE3
C###################################################################
      DIMENSION YF(MAXNE),JBCONV(2,50),CROS(MAXPR,2)
C      DATA SECANNI/3.15587D+07/
      PAS=ht1*SECANNI/100.D0
c       PASSO=100.d0*SECANNI 
      EQUIHE3=.FALSE.
C ************ MAIN LOOP *************
c      write(*,*) 'entro in smallrap'
      pasanni=0.d0
      the3=0.d0
      isomma=0

      if(JBCONV(1,1).eq.1)then
        DO J=1,MAXNE
         YF(J)=XXX(J,1)/ATW(J)
        END DO
        coeff1=0.d0
       coeff2=0.d0
       coeff3=0.d0
       coeff4=0.d0
       coeff5=0.d0
       do K=1,JBCONV(2,1)-1

        T=1.0D+06*(G(4,K)+G(4,K+1))/2.D0
         P=1.0D+17*(G(3,K)+G(3,K+1))/2.D0

         CALL EOS(K,1,MAXNE,1,P,T,XXX,'SMALLRAP',RHO,O1,O2,O3,G1,DEL,
     #EMUE,1)
         CALL SKIPFASE (T,RHO,YF,JF,IFASE)
        IF(IFASE.EQ.1) then
           NNP = NP(IFASE)
            CALL NCROSS (RHO,T,YF,IFASE,NNP,CROS,1)
             coeff1=coeff1+(rho*cros(3,1)*(G(5,K+1)-G(5,K)))
            coeff2=coeff2+(rho*cros(4,1)*(G(5,K+1)-G(5,K)))
c        ELSEIF(IFASE.EQ.2) then
c          NNP = NP(IFASE)
c            CALL NCROSS (RHO,T,YF,IFASE,NNP,CROS,1)
c        if(k.eq.1)write(*,*)cros(1,1),cros(2,1),cros(3,1),cros(5,1)
c        if(k.eq.JBCONV(2,1)-1)
c     #write(*,*)cros(1,1),cros(2,1),cros(3,1),cros(5,1)
c             coeff3=coeff3+(rho*cros(2,1)*(G(5,K+1)-G(5,K)))
c            coeff4=coeff4+(rho*cros(3,1)*(G(5,K+1)-G(5,K)))
c            coeff5=coeff5+(rho*cros(5,1)*(G(5,K+1)-G(5,K)))
        ENDIF
        enddo
       if(coeff1.le.0.d0.and.coeff2.le.0.d0)then
        THE3=1.d99*secanni
       else
        THE3=G(5,JBCONV(2,1))/(coeff1*YF(2)+coeff2*YF(3))
        endif
         if(ht1*secanni.gt.the3)then
             isomma=1
          if(pas.lt.5.d0*the3)then
           xlim=xini-xini*0.2d0
           write(*,'(d14.7)')xlim
           if(xxx(1,1).lt.xlim)then
             write(*,*)'----------------------------------------'
             write(*,*)'ht1/100 < 5the3 quando l`H centrale e` >'
             write(*,'(d14.7)')xlim
             write(*,*)'----------------------------------------'
           else
             yhini=yf(1)
              EQUIHE3=.true.
             pas=5.d0*the3
             write(*,*)'ht1 > the3 con pas < 5the3'
           endif
          else
             pas=the3*5.d0
             write(*,*)'ht1 > the3 con pas > 5the3'
          endif
        endif
c         write(*,*)'the3      pasanni'
c         write(*,*)the3/secanni,pas/secanni
c        if(coeff3.le.0.d0.and.coeff4.le.0.d0.and.coeff5.le.0.d0)then
c        Tc12=1.d99*secanni
c        To16=1.d99*secanni
c        To18=1.d99*secanni
c       else
c        Tc12=G(5,JBCONV(2,1))/(coeff3*YF(3))
c        To16=G(5,JBCONV(2,1))/(coeff4*YF(3))
c        To18=G(5,JBCONV(2,1))/(coeff5*YF(3))
c        endif
c      write(*,*)pas/secanni,Tc12/secanni,To16/secanni,To18/secanni
      endif

c      if(JBCONV(1,1).eq.1)then
c        DO J=1,MAXNE
c         YF(J)=XXX(J,1)/ATW(J)
c        END DO
c       K=JBCONV(2,1)-1
c
c        T=1.0D+06*(G(4,K)+G(4,K+1))/2.D0
c         P=1.0D+17*(G(3,K)+G(3,K+1))/2.D0
c
c         CALL EOS(K,1,MAXNE,1,P,T,XXX,'SMALLRAP',RHO,O1,O2,O3,G1,DEL,
c     #EMUE,1)
c         CALL SKIPFASE (T,RHO,YF,JF,IFASE)
c        IF(IFASE.EQ.1) then
c           NNP = NP(IFASE)
c            CALL NCROSS (RHO,T,YF,IFASE,NNP,CROS,1)
c           the3=1/((YF(2)*cros(3,1)+YF(3)*cros(4,1))*rho)
c        ENDIF
c         if(ht1*secanni.gt.the3)then
c          if(pas.lt.the3)then
c            yhini=yf(1)
c             EQUIHE3=.true.
c            pas=the3*1.5d0
c          else
c            if(pas.gt.the3*1.5d0)pas=the3*1.5d0
c          endif
c        endif
c         write(*,*)the3/secanni,pas/secanni
c      endif
           
           
           
           
      DO JJ=1,NCONV
         JINI=JBCONV(1,JJ)
         JFIN=JBCONV(2,JJ)
       DO 1 K=JINI,JFIN
         KEQUIL(K)=0
         DO J=1,MAXNE
            YF(J)=XXX(J,K)/ATW(J)
         END DO
C ************** CALCOLA RO **********
         T=1.0D+06*(G(4,K)+G(4,K+1))/2.D0
         P=1.0D+17*(G(3,K)+G(3,K+1))/2.D0
c        T=1.0D+06*G(4,K)
c         P=1.0D+17*G(3,K)
         CALL EOS(K,1,MAXNE,1,P,T,XXX,'SMALLRAP',RO,O1,O2,O3,
     #G1,DEL,EMUE,1)
C ####### SKIPPA FASI INUTILI E RICERCA FASE EVOLUTIVA###############
         CALL SKIPFASE (T,RO,YF,K,IFASE)
         IF(IFASE.NE.1)goto 1
C ###################################################################
         NNE=NE(IFASE)
         NNP=NP(IFASE)
C ******************** CONTROLLO INTEGRAZIONE: **********************
C ** CHIMICA MIXED -> SERIE DI TAYLOR   CHIMICA LOCALE -> RAPNEW ****
C **************  SOLO SE CENTRALE ****************************
C####################################################################

         CONVETTA=.FALSE.
         CALL RAP(RO,T,PAS,IFASE,NNE,NNP,YF,K,CONVETTA)
        
        if(EQUIHE3.and.jini.eq.1)then
          yxh=atw(1)*(yf(1)-yhini)
          yf(1)=yhini
          XXX(3,K)=Yf(3)*atw(3)+yxh
          yf(3)=XXX(3,K)/atw(3)
          if(k.eq.1)then
            write(*,*)yxh
          endif
        endif        
        
C####################################################################               
         DO J=1,MAXNE
            XXX(J,K)=YF(J)*ATW(J)
         END DO
         IF(IFASE. EQ. 3) THEN
            CARXY(1,K)=XXX(1,K)
            CARXY(2,K)=XXX(3,K)
            XXX(1,K)=0.D0
            XXX(3,K)=0.D0
         ENDIF
         do ml=1,maxne-1
            if(xxx(ml,k).lt.0.d0) then
               if(dabs(xxx(ml,k)).gt.1.d-10) then
                   write(*,*)'DT in SMALLRAP TROPPO LUNGO!!!!!!'
                   write(*,'(2i5,3x,e13.5)')k,ml,xxx(ml,k)
                 stop
              end if
            end if
         end do
        
        
    1 CONTINUE
      ENDDO
      
c      do jj=1,maxme-1
c         do ml=1,maxne-1
c            if(xxx(ml,jj).lt.0.d0) then
c               if(dabs(xxx(ml,jj)).gt.1.d-10) then
c                   write(*,*)'DT in SMALLRAP TROPPO LUNGO!!!!!!'
c                   write(*,'(2i5,3x,e13.5)')jj,ml,xxx(ml,jj)
c                 pause
c              end if
cc                 stop
c            end if
c         end do
c      end do
c################################################
      IF(ISOMMA.GT.0)THEN
c        write(*,*)ht1,PAS
        PASANNI=PAS/SECANNI
        HT1=HT1+PASANNI
       if(equihe3)then
         HT1=HT1-PASANNI
       endif
c       write(*,*)ht1,PAS
      ENDIF
c################################################      
      RETURN
      END
