      PROGRAM BIGTAB2P1
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 AGE,FNAME
      integer igio,imes,ianno
      CHARACTER*51 vER
      CHARACTER*100 CROSSEC,MIXT,EOS,OPAC
      CHARACTER*150 SPA
      CHARACTER*6 EVFASE
      CHARACTER*6 DIF
      DIMENSION A(58)
C=========================================================================
C  Questa versione deve essere utilizzata solo ed esclusivamente con i 
C  grafi generati con i codici FRANEC la cui versione e` almeno 2.1.0
C=========================================================================
C
C   A(1)= PASSO TEMPORALE   
C   A(2)= LUMINOSITA` SUPERF  
C   A(3)= TEMP EFF.
C   A(4)= RAGGIO TOT  
C   A(5)= MASSA TOT
C   A(6)= PRESS. CENTR  
C   A(7)= TEMP CENTR
C   A(8)= DENSITA CENTR
C   A(9)= H CENTR
C   A(10)= He CENTR
C   A(11)= C CENTR
C   A(12)= N CENTR
C   A(13)= O CENTR
C   A(14)= LUMINOSITA` P-P
C   A(15)= LUMINOSITA` CNO
C   A(16)= LUMINOSITA` 3-ALFA
C   A(17)= LUMINOSITA` GRAVITA
C   A(18)= PERIOD SPACING
C
C   A(19)= MASSA per TEMP MAX FUORI CENTRO
C   A(20)= PRESSIONE per TEMP MAX FUORI CENTRO
C   A(21)= TEMP MAX FUORI CENTRO
C   A(22)= DENSITA per TEMP MAX FUORI CENTRO
c
C   A(23)= MASSA per MAX ENERGIA H-BURNING
C   A(24)= MAX ENERGIA H-BURNING
C   A(25)= RAGGIO/RTOT per MAX ENERGIA H-BURNING
C   A(26)= LUMINOSITA  per MAX ENERGIA H-BURNING
C   A(27)= PRESSIONE per MAX ENERGIA H-BURNING
C   A(28)= TEMPERATURA  per MAX ENERGIA H-BURNING
C   A(29)= DENSITA  per MAX ENERGIA H-BURNING
C   A(30)= MASSA SHELL H
C
C   A(31)= MASSA per MAX ENERGIA He-BURNING
C   A(32)= MAX ENERGIA He-BURNING
C   A(33)= RAGGIO/RTOT per MAX ENERGIA He-BURNING
C   A(34)= LUMINOSITA  per MAX ENERGIA He-BURNING
C   A(35)= PRESSIONE per MAX ENERGIA He-BURNING
C   A(36)= TEMPERATURA  per MAX ENERGIA He-BURNING
C   A(37)= DENSITA  per MAX ENERGIA He-BURNING
C   A(38)= MASSA SHELL He
C        CORE CONVETTIVO 
C   A(39)= MASSA CORE CONVETTIVO / MTOT
C   A(40)= RAGGIO CORE CONVETTIVO / RTOT
C   A(41)= LUMINOSITA CORE CONVETTIVO / LUMINOSITA SUPERFICE
C   A(42)= PRESSIONE CORE CONVETTIVO
C   A(43)= TEMPERATURA CORE CONVETTIVO
C   A(44)= H  CORE CONVETTIVO
C   A(45)= He  CORE CONVETTIVO
C   A(46)= C  CORE CONVETTIVO
C   A(47)= N  CORE CONVETTIVO
C   A(48)= O  CORE CONVETTIVO
C        BASE INVILUPPO CONVETTIVO
C   A(49)= MASSA  BASE INVILUPPO CONVETTIVO / MTOT
C   A(50)= RAGGIO  BASE INVILUPPO CONVETTIVO / RTOT
C   A(51)= LUMINOSITA  BASE INVILUPPO CONVETTIVO
C   A(52)= PRESSIONE  BASE INVILUPPO CONVETTIVO
C   A(53)= TEMPERATURA  BASE INVILUPPO CONVETTIVO
C   A(54)= H  BASE INVILUPPO CONVETTIVO
C   A(55)= He  BASE INVILUPPO CONVETTIVO
C   A(56)= N  BASE INVILUPPO CONVETTIVO
C   A(57)= C  BASE INVILUPPO CONVETTIVO
C   A(58)= O  BASE INVILUPPO CONVETTIVO
C
C

C========================================================================
      
      OPEN(UNIT=1,FILE='grafi2p1',STATUS='OLD')
      OPEN(UNIT=7,FILE='asttab2p1',STATUS='NEW')
      OPEN(UNIT=2,FILE='bigtab2p1',STATUS='NEW')
      READ( 1,'(18X,2I3,I5,10x,a12)')igio,imes,ianno,VER
      WRITE(*,'(18X,2I3,I5,10x,a12)')igio,imes,ianno,VER
      READ( 1, '(A150)' )SPA
C      WRITE(*, '(A150)')SPA
      READ( 1 ,'(A100)')CROSSEC
      WRITE(*,'(A100)')CROSSEC
      READ( 1 ,'(A100)')MIXT
      WRITE(*,'(A100)')MIXT
      READ( 1 ,'(A100)')EOS
      WRITE(*,'(A100)')EOS
      READ( 1 ,'(A100)')OPAC
      WRITE(*,'(A100)')OPAC
      READ( 1, '(A150)')SPA
C      WRITE(*, '(A150)' )SPA
      READ( 1, 770)AM,Y,Z,AML,ETA,OV,DIF,EVFASE,PROGM
      WRITE(*,772)AM,Y,Z,AML,ETA,OV,DIF,EVFASE,PROGM
      READ( 1, '(A150)')SPA
C      WRITE(*, '(A150)' )SPA
      WRITE(2,132) igio,imes,ianno,ver,AM,Y,Z,AML,ETA,OV,DIF,EVFASE,
     #PROGM
      WRITE(2,103)
      WRITE(2,133)      
      WRITE(7,132) igio,imes,ianno,ver,AM,Y,Z,AML,ETA,OV,DIF,EVFASE,
     #PROGM
      WRITE(7,107)
      WRITE(7,133)      
    1 CONTINUE
C indented stuff: 1 int, 5 doubles
      READ( 1 , 100 , END = 3 ) NMOD , AGE , ( A(L) , L =  1 , 3 )
C non-indented stuff: 11 lines x 5 numbers per line = 55 objects
      READ( 1 , 101 , END = 3 )              ( A(L) , L =  4 , 58 )
      A(49)=A(49)*A(5)
      A(39)=A(39)*A(5)
      A(19)=A(19)*A(5)
      BB=1.-A(14)-A(15)-A(16)-A(17)
      IF(A(9).GE.1.E-30) THEN
        WRITE(2,104)NMOD,AGE,A(9),A(2),A(3),A(7),A(8),A(39),A(23)
     $  ,A(31),A(49),A(14),A(15),A(16),A(17),A(55),A(5),A(21)
         WRITE(7,108)NMOD,AGE,A(9),A(2),A(3),A(4),A(39)/A(5),A(40),
     #               A(49)/A(5),A(50),A(5),A(18)
      ELSE
        WRITE(2,104)NMOD,AGE,A(10),A(2),A(3),A(7),A(8),A(39),A(23)
     $  ,A(31),A(49),A(14),A(15),A(16),A(17),A(55),A(5),A(21)
         WRITE(7,108)NMOD,AGE,A(10),A(2),A(3),A(4),A(39)/A(5),A(40),
     #               A(49)/A(5),A(50),A(5),A(18)
      ENDIF
      GO TO 1
    3 STOP
  100 FORMAT(I5,1PD22.15,3D16.8)
  101 FORMAT(1P5D16.8)
  102 FORMAT(A8)
  132 FORMAT(2I3,I5,' == ',a12,' | Mass=',F8.4,1X,'|',1x,'Y=',F6.3,1X,'|
     #',1x,'Z=',1P,E9.2,1X,'|',1x,'ML=',0P,F5.2,1X,'|',1x,'Mass loss=',
     #F4.2,1X,'|',1x,'Oversh= ',f5.3,1X,'|',1x,'Diff: ',A3,1X,'|',1x,
     #'from: ',a6,1X,'|',1x,'Mass of prog.=',F8.4)

  133 FORMAT('==========================================================
     #==================================================================
     #================================================')
  107 FORMAT(' #mod',5X,'log(t)',5x,'H/HE',5x,'log(L/Lo)',3x,'logTe',7x,
     #'log(R/Ro)',4x,'Mcc/M',7x,'Rcc/R',7x,'Mce/M',7x,'Rce/R',6x,
     #'Mtot/Mo',6x,'P-sp')
  103 FORMAT(' #mod',4X,'log(t)',5x,'H/HE',6x,'logL',6x,'logTe',3x,'logT
     #c',2x,'logRc',4x,'Mcc',4x,'M_cHe',3x,'M_cCO',4x,'Mce',5x,
     #'Lpp/Ls',5x,'Lcno/Ls',6x,'L3a/Ls',6x,'Lgr/Ls',5x,
     #'He_sup',5x,'Mtot',5x,'log(Tmax)')
  104 FORMAT(I6,F12.8,1P,E9.2,0P,2F10.6,2F7.3,4F8.4,1P,4(1X,E11.4),
     # 1P,E10.2,0P,2F10.6)
  108 FORMAT(I6,F12.8,1P,E9.2,0P,8F12.8,1P,E11.3)

  770 FORMAT (5X,F8.4,5X,F6.3,5X,1P,E9.2,
     #6X,0P,F5.2,13X,F4.2,11X,f5.3,9X,A3,9X,a6,17X,F8.4)

  772 FORMAT ('Mass=',F8.4,1X,'|',1x,'Y=',F6.3,1X,'|',1x,'Z=',1P,E9.2,
     #1X,'|',1x,'ML=',0P,F5.2,1X,'|',1x,'Mass loss=',F4.2,1X,'|',1x,
     #'Oversh= ',f5.3,1X,'|',1x,'Diff: ',A3,1X,'|',1x,'from: ',a6,
     #1X,'|',1x,'Mass of prog.=',F8.4)
      END
