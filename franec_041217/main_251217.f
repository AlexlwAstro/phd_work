C*****************************************************************************
C  QUESTA ROUTINE CONTIENE: MAIN - DEJA - FATO - FITTA - GTIME - HBREZO - 
C     HENYEY - INNES - MASLOS - NEWBVAL - PLOTTA -RESNUC - RWMOD - SKIPFASE - 
C     STAMPA - STOPERR - VEIOVE, also RSTAM
C*****************************************************************************
C
C
C****************************************************************************
C NELLA ROUTINE "ENERGY" CI SONO:
C         ENERGY - EPSI - EPSIG - EVOLUT - NCROSS -
C         INDICIZZA - RAP - READNET - PLASMA - SC - SMALLRAP.
C****************************************************************************
C****************************************************************************
C NELLA ROUTINE "FISICA" CI SONO:
C         EOS - NEUTR - NKAPPA - NSUPERA - SIGMANEW 
C****************************************************************************
C****************************************************************************
C NELLA ROUTINE "MATH" CI SONO:
C          KERNEL - NCUB - PARAB - POLINT 
C****************************************************************************
C****************************************************************************
C NELLA ROUTINE "MIXING" CI SONO:
C         BCONV - MIXINGALL - MIXINGHE - OVERSH
C****************************************************************************
C****************************************************************************
C NELLA ROUTINE "NATMOS" CI SONO:
C         CONDCON - LEGGI - NATMOS
C****************************************************************************
C****************************************************************************
C NELLA ROUTINE "OPTIMA" CI SONO:
C             LEVMES - NEWMES - OPTIM - PASTEM - QUATM 
C****************************************************************************
C
      PROGRAM FRANEC    
C
C
C  VERSIONE 2.1.0 - last change - 7-05-2015
C
C 1) eliminate alcune chiamate inutili alla state per la nuova def. di Q
C 
C 2) corretto errore posizionamento EOS nel caso puro idrogeno 
C 
C 3) 30-09-2004  MODIFICA "RAP" (INTRODOTTO DYDYM(MAXNE))
C
C 4) 20-01-2006 modificati:
C   schema di burning in rap  YM=(YF+YV)/2 
C   dimezzato passo temporale in RGB e cambiato 1.d-20 in 1.d-5 nel 
C        controllo passo temporale in He burning
C   optim per la parte che riguarda il # di mesh per il convettivo
C        centrale alla fine dell'He burning
C
C 5) 30-10-2006 introdotta nuova routine eos.f al posto della NSTATE
C
C 6) 12-4-2007 modifiche come da quaderno Santi; eliminata NSTATE e INDCONV
C
C 7) 16-4-2007 eliminata divisione passo temporale in RGB; introdotto VALR e
C              ridefiniti VALP, VALUE in combustione di He nella PASTEM
C
C 8) 17-4-2007 eliminati commenti e chiamate superflue
C              
C 9) 18-4-2007 eliminato INDI = ITEM + N in NKAPPANODIFF
C10) 18-4-2007 eliminato EQUIVALENCE in NKAPPADIFF
C              
C11) 18-4-2007 imposto HT1=5000.D0/(EMSOL**3*G(4,1)**2/225.D0) invece di
C                HT1=50000.D0/(... per VLM da ZAMS
C12) 23-4-2007 inclusa eos.f in fisica.f
C
C13) 02-7-2007 
C     modificata nella subroutine eos la condizione
C        IF(G(6,1).LT.0.D0.AND.G(4,1).LT.150.D0)THEN
C     con
C        IF(G(6,1).LT.0.D0.AND.G(4,1).LT.150.D0.
C     #AND.XXX(4,1).LT.0.2D0) THEN
C14) 29-1-2008 
C     modificata la subroutine SMALL RAP per
C     equilibrio He3 e genera nuovo file tempi.dat 
C15)  1-2-2008 
C     modificata la subroutine SMALL RAP per
C     controllare abbondanze negative elementi solo nei convettivi 
C16)  4-2-2008 
C     limite H per diffusione a1.d-10 invece che 0.d0
C17) 27-1-2009 
C     Corretto posizionamento in pressione al bordo superiore nella EOS
C     Corretto posizionamento in temperatura al bordo superiore nella EOS
C     Corretto posizionamento in chimica nel caso nucleo CO INI=8 --> INI=7 
C18) 04-4-2011 
C     Corretta linea 947 nella energy.f per eliminare errata abbondanza di elio
C     nell'inviluppo convettivo in fase di H burninc in core convettivo:
C     if(EQUIHE3)then --> if(EQUIHE3.and.jini.eq.1)then
C19) 07-06-2013 
C     modificata OPTIM per migliorare la ripartizione in mesh
C20) 07-06-2013
c     aggiunto RSTAM per scrittura file per astrosismologia, adesso legge
C     il file rlim.txt in cui é indicato il raggio limite per cui scrive in
c     fort.39 le grandezze per astrosismologia da dare come input al programma 
c     che da Gamma1 e altre grandezze nel formato fgong
c21) 07-06-2013
C     inserita nuova routine per diffusione (DIFFUS) che tratta le zone 
C     centrali con regressione polinomiale delle velocitá di diffusione
c22) 07-06-2013
C     modificati i formati per mixing lenght e elio iniziale nel modstart
C     per avere maggiore accuratezza nel calcolo del MSS
c23) 07-06-2013 
c     modificato const.2p0 adottando le formule per il calcolo di STB10,
c     STBSUN e REIMAS.
c24) 07-06-2013 
c     ripristinato lo sviluppo in serie nella henyey (questo consente di
c     ottenere risultati migliori per il dln(M)/dln(R) per le regioni centrali
c25) 10-10-2014 
c     modificata OPTIM considerando parametri diversi per le diverse zone
c     ripartite in massa per migliorare la ripartizione in mesh nelle zone
c     esterne
c26) 04-11-2014 
c     adottate le nuove sezioni d´urto per:
c       He3+He4 -> Be7+G      Cyburt e Devis (2008) Phys. rev. C 
c       N14+p -> O15+g        Formicola et al. (2003)
c       C12+ALPFA -> O16+G    Hammer et al. 2005esterne
c27) 05-11-2014 
c     Modificata la INNES per includere le abbondanze primordiali di
c      D=3.90000D-05
c      He3=2.30000D-05
c      Li7=2.60000D-09
c     solo partendo dalla presequenza
c28) 22-03-2015 
c     Modificata la EOS per includere le grandezze termodinamiche finalizzate
c     al calcolo delle strutture per l´astrosismologia; le modifiche tengono
c     conto anche del differente formato in cui viene scritta la eos, che 
c     é stata estesa ad alte pressioni per il calcolo delle VLM
c     solo partendo dalla presequenza
c29) 12-04-2015 
c     Introdotta la routine RSTAM per la scrittura dei files fgong 
c     in cui calcoliamo tutte le grandezze necessarie per scrivere il file fgong
c     
c     A tal fine e' stata modificata anche la natmos introducendo la matrice
c     BPP che viene poi passata alla RSTAM attraverso un common.
c30) 23-04-2015 
c     Modificata la OPTIM per migliorare la definizione dei mesh nel 
c     passaggio radiativo-convettivo nell envelope convettivo; questo elimina 
c     il problema del BUMP dell RGB che per strutture prossime alla scomparsa, 
c     presentavano un andamento non monotono decrescente con la massa.
c     E' stato inoltre ripristinato il controllo in massa sui mesh, per 
c     evitare che la differenza in massa tra due mesh possa essere inferiore
c     a 1.D-12 limite di risoluzione nella scrittura dei files di ripartenza.
c     (senza questo controllo alla ripartenze si puo avere durante le iterazioni
c     per la convergenza che E(2,mesh)=infinity
c
C KSA = IPRALL
C fix number of models stored in grafi?
C KSB=N PRINTA OGNI N MODELLI
C number of models in stampe?
C KSC=0 TEMPO INFINITO - = N NORE DI ESECUZIONE
C 0 = no stop (9999)
C KSD=K SCRIVE OGNI K MESH
C stampe - store at each k mesh number
C KSE=1 NON USA LA OPTIM
C at 1, does not use optimisation of mesh number
C KSF=1 USA IL REZONING DELL'ATMOSFERA
C 1 = modify mesh point distribution in atmosphere
C KSG=N GRAFICA TEMPORALE OGNI N MODELLI
C number of models stored in grafi
C KSH=M GRAFICA STRUTTURA OGNI M MODELLI
C same but in stampe (i.e. structure)
C KSW=1 PER LAVORI IN BATCH
C 1 = write to screen (never used)
C ISINGOLA=1 ATMOSFERA SINGOLA
C 1 = force the computation of new atmosphere (don't use)
C NSHUT=1 LAVORA CON OVER
C 1 = account for overshooting, 0 = don't
C INDIFF=1 INCLUDE DIFFUSIONE
C 1 = consider atomic diffusion
C KFGONG=1 STORE FGONG FILE - 0 NO FGONG FILE
C fgong = asteroseismology file - store/don't store
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MZ=100)
      LOGICAL*1   IFLAG,MACINA
      CHARACTER*11 NOME(MZ),MODELLO
      CHARACTER*10 STRIFGONG
      CHARACTER*16 NOMEFGONG
      CHARACTER*16 NOMEFGONGD
      CHARACTER*18 NOMEFGONGETA
      CHARACTER*18 NOMEFGONGETAD
      CHARACTER*12 VERSIONE
      CHARACTER*8 FMT 
      CHARACTER*5 X1
      CHARACTER*4 X2
      CHARACTER*2 Z2
      CHARACTER*3 ZET,ELIO,Y1
      CHARACTER*25 FILENAME ! format complete name fgong file name+nmd
      
      DIMENSION ADX(3),FM(3)
      DIMENSION ahlim(2,10),ahelim(2,10),rlim(2,10)
      integer IPRESSSTART,IPRESSSTOP,IPRESSSWITCH
      INCLUDE 'maincom.2p1'
      INCLUDE 'atmosfer.2p1'
      INCLUDE 'consts.2p1'
      COMMON/FREMEM/GGV(5,LIM),O(MAXNE+1-5,LIM)
      COMMON/OVER/FACTOR
      COMMON/OPACITY/INDIFFINI
      COMMON/REWRMOD/NMOD,KSA,KSB,KSC,KSD,KSE,KSG,KSH,KSW,NITER,NPER,
     #KFGONG
      
      common/NOSUB/IPRESSSTART,IPRESSSTOP,IPRESSSWITCH

      
      DIMENSION TRAS(5),JBCONV(2,50)
      DIMENSION AMZAHB(MZ),ROR(3),AMOMIN(20),AMOMAX(20),ZZZ(20)
      DATA TRAS/1.D-10,1.D-32,1.D-17,1.D-06,1.D-33/    
      DATA NABLA/-1/,NUM1/0/
      DATA AMOMIN/1.30D0,1.30D0,1.30D0,1.30D0,1.30D0,1.30D0,1.32D0,
     #1.32D0,1.24D0,1.21D0,1.17D0,1.13D0,1.10D0,1.09D0,1.08D0,1.08D0,
     #1.09D0,1.11D0,1.10D0,1.09D0/
      DATA AMOMAX/2.09D0,1.78D0,1.68D0,1.59D0,1.54D0,1.50D0,1.47D0,
     #1.45D0,1.44D0,1.43D0,1.42D0,1.42D0,1.42D0,1.42D0,1.42D0,1.42D0,
     #1.43D0,1.47D0,1.42D0,1.40D0/
      DATA ZZZ/1.D-5,5.D-5,1.D-4,1.98D-4,3.14D-4,4.43D-4,6.24D-4,
     #7.87D-4,9.9D-4,1.397D-3,1.97D-3,3.11D-3,3.9D-3,6.14D-3,7.7D-3,
     #9.64D-3,1.721D-2,2.081D-2,2.865D-2,3.905D-2/
C***********************************************************************
      OPEN(UNIT=2,FILE='stampe')
      OPEN(UNIT=4,FILE='modstart')
C     OPEN(UNIT=9,FILE='CNO')
C     REWIND(9)
      OPEN(UNIT=10,FILE='grafi2p1')
C     OPEN(UNIT=14,FILE='grchi-1')
      OPEN(UNIT=16,FILE='chimsup2p1')
      OPEN(UNIT=66,FILE='chimcen2p1')
      OPEN(UNIT=46,FILE='flussi2p1.dat')
      OPEN(UNIT=56,FILE='tempi2p1.dat')
      WRITE(46,*)'  AGE            L/Lsun         Fpp            Fdp            
     #   Fhe3HE3         Fhe3he4          Fbe           Fli7p            
     #Fbe7p         Ftotpp'
C***********************************************************************
      IFLAG=.TRUE.
C     CALL UNDER0(IFLAG)
      WRITE(*,100)
      WRITE(2,100)
      VERSIONE='FRANEC 2.1.2'
      MODELLO='modant2'
      IPRALL=0
      DMAS=0.D0
      IPP=0
      NMOD=0
      TEMPO=1.D-09
      DAM=0.D0
      ITERAZ=0
      ERR=0.D0
      IIST=1
      ISTA=1
C######### PER INDIVIDUARE  inizio He CENTRALE ###########
      NMDHEI=0
C######### PER INDIVIDUARE  FINE He CENTRALE ###########
      NMDHE=0
C######### PER INDIVIDUARE 100 MODELLI DOPO FINE He CENTRALE ###########
      NMDCO=1
C#######################################################################
C#######################################################################
      NNMD=0
C#########per il passo temporale########################################
      nmlim=0
C#######################################################################
C############per andare con condizioni al contorno######################
      indatm=1
      nmlimatm=0
C#######################################################################
C*******   ISTART => 1(MODANT) - 2(PRESEQUENZA) - 3(MS) - 4(HB)   ******
C*******          => 5 (ZAHB)                                     ******
      READ(4,111)ISTART,NITER,NPER,ISUB
      READ(4,951)KSA,KSB,KSC,KSD,KSE,KSF,KSG,KSH,KSW,ISINGOLA,NSHUT
     #           ,INDIFFINI,KFGONG

      INDIFF=0
C#######################################################################    
      CALL gtime(igio,imes,ianno)
c===== writes in date and version
      WRITE(10,770)IGIO,IMES,IANNO,VERSIONE
      WRITE(10,773)
 770  FORMAT ('Model computed on ',2I3,I5,' by using ',A12)
 773  FORMAT ('=========================================================
     #==================================================================
     #============')

      CALL READNET
c===== modifica per eliminare subatmosfera
      
        IPRESSSWITCH=1   
        IPRESSSTART=0
        IPRESSSTOP=0
      IF(ISTART.EQ.2)THEN
        IPRESSSWITCH=0   ! SE DIVERSO DA 0 NON ELIMINA SUBATMOSFERA
        IPRESSSTART=5
        IPRESSSTOP=50+IPRESSSTART-1
      ENDIF
c===== fine  modifica per eliminare subatmosfera



      IF(ISTART.EQ.2) THEN
C****************************** FITTING PRESEQUENZA ********************
C**************************** pre-main sequence fitting
         ITCC=0
         CALL INNES(ITCC,IPP)
         EN1=0.D0
         EN2=0.D0
         ITCC=1
         ELLOGV=ELLOG
         TEFFV=TEFF
         G(5,MAXME+1)=G(5,MAXME)
         DO K=1,MAXME
            DO J=1,5
               G(J,K)=G(J,K)*TRAS(J)
               GG(J,K)=G(J,K)
               GGV(J,K)=G(J,K)
            END DO
         END DO
       
         DO K=1,MAXME-1
            EN1=EN1+G(4,K)*(G(5,K+1)-G(5,K))
         END DO
         CALL INNES(ITCC,IPP)
         DO K=1,MAXME-1
            EN2=EN2+1.D-06*G(4,K)*(GGV(5,K+1)-GGV(5,K))
         END DO
         HT1=DABS(EN1-EN2)/(120.55*(10.D0**ELLOG+10.D0**ELLOGV))*2.5D+8
         IPP=1
         CALL INNES(ITCC,IPP)
         RAZ=(10.D0**ELLOGV)/(10.D0**ELLOG)
         DO K=1,MAXME
            DO J=1,5
               G(J,K)=G(J,K)*TRAS(J)
            END DO
            GG(2,K)=G(2,K)*RAZ
            GGV(2,K)=G(2,K)*RAZ
         END DO
       
         HT1V=HT1
         TEMPO=HT1
         XCE(3,1)=XXX(3,1)
         MAXMV=MAXME
         EMTOV=EMTOT
         LUCA=1
         NMD=1
C***********************************************************************
      ELSE IF(ISTART.EQ.3.OR.ISTART.EQ.5) THEN
C***************************** FITTING MS ED HB & ZAHB *****************
C**************** MS and HB & ZAHB fitting
         ITCC=0
         CALL INNES(ITCC,IPP)
       EMSOL=EMTOT/(SUNMg33)
         DO K=1,MAXME
            DO J=1,5
               G(J,K)=G(J,K)*TRAS(J)
               GG(J,K)=G(J,K)
               GGV(J,K)=G(J,K)
            END DO
         END DO
         TEFFV=TEFF
       HT1=5000.D0/(EMSOL**3*G(4,1)**2/225.D0)
         IF(ISTART.EQ.4.OR.ISTART.EQ.5)HT1=HT1/500.D0
         IF(ISTART.EQ.3.AND.XXX(1,1).LE.1.D-5)HT1=HT1/1000.D0
         HT1V=HT1
         MAXMV=MAXME
         EMTOV=EMTOT
         LUCA=1
         NMD=1
C=======================================================================
C==== LETTURA FILE DEFINIZIONE MASSE DI ZAHB E NOMI FILES ==============
         IF(ISTART.EQ.5) THEN
        OPEN(UNIT=18, FILE='initial.zahb')
        READ(18,'(2I6,F12.8,F7.4)') NTRK, NRIP, TFLASH, PROGM ! N. MODELLI ZAHB, T_flash, massa prog   
         DO K = 1, NTRK
              READ(18,434) AMZAHB(K),NOME(K)  ! MASSE ZAHB E NOME FILE
         END DO 
         TFLASH=(10.D0**TFLASH)-1.0D6 ! ETA AL FLASH - 1Myr
       IIST=0
       ENDIF  
C***********************************************************************
      ELSE IF(ISTART.EQ.1.OR.ISTART.EQ.4) THEN
C***************************** LETTURA MODANT **************************
         CALL RWMOD(0,'NNNN',INDIFFINI)
         LUCA=1
         NMD=NMOD+LUCA
         IF(IPRESSSWITCH.EQ.0)THEN
          IPRESSSTART=5
          IPRESSSTOP=50+IPRESSSTART-1
         ENDIF
c         IF(XXX(1,1).LE.0.d0)INDIFF=0
C***********************************************************************
      END IF
C**************************   FINE  LETTURA  ***************************
c===== modifica per eliminare subatmosfera in caso di modelli di HB senza mass loss
c      IF(ISTART.EQ.4.AND.ETA.EQ.0.D0)THEN
c        IPRESSSWITCH=0   ! SE DIVERSO DA 0 NON ELIMINA SUBATMOSFERA
c        IPRESSSTART=5
c        IPRESSSTOP=50+IPRESSSTART-1
c      ENDIF
c===== fine modifica per eliminare subatmosfera in caso di modelli di HB senza mass loss

      EMSOL=EMTOT/(SUNMg33)
      YINI=1.D0-XINI-ZINI
      WRITE(2,444)EMSOL,YINI,ZINI,ALPHA
      WRITE(*,444)EMSOL,YINI,ZINI,ALPHA
      NABLA=3
      
C***************** here compute lambda overshooting (FACTOR) **********************
       FACTOR=0.D0
       IF(nshut.gt.0)THEN
       SMIN=0.D0  !MASSA Max with CORE-CONV (SOLAR UNIT) lower than 0.04Mo
       SMAX=0.D0  !MASSA Min with CORE-CONV (SOLAR UNIT) always gtreater than 0.04Mo
c        WRITE(*,'(1P,D10.3)')ZINI
	DO I=1,20
         IF(ZINI.EQ.ZZZ(I))THEN
	   SMIN=AMOMIN(I)
	   SMAX=AMOMAX(I)
	   GOTO 80
	 ELSEIF(ZINI.GT.ZZZ(I).AND.ZINI.LT.ZZZ(I+1))THEN
	   DMZ=(ZINI-ZZZ(I))/(ZZZ(I+1)-ZZZ(I))
	   SMIN=AMOMIN(I)+(AMOMIN(I+1)-AMOMIN(I))*DMZ
	   SMAX=AMOMAX(I)+(AMOMAX(I+1)-AMOMAX(I))*DMZ
c            WRITE(*,'(1P,2D10.3)')SMIN,SMAX
	   GOTO 80
	 ELSEIF(ZINI.GT.ZZZ(20).or.ZINI.LT.ZZZ(1))THEN
	   WRITE(*,*)'WARNING: OVER BUT NO OVER'
	 ENDIF
	ENDDO 
 80    CONTINUE
        IF(STARTMASS.GT.SMIN)THEN
          IF(STARTMASS.LE.SMAX)THEN
           FACTOR=0.2d0/(SMAX-SMIN)*(STARTMASS-SMIN)
          ELSE
           FACTOR=0.2d0
          ENDIF
        ELSE
	  FACTOR=0.d0
	ENDIF
       ENDIF
C***************** end computing lambda overshooting (FACTOR) **********************


C#######################################################################
C####### KFGONG = 0 NO STORE FGONG FILE - KFGONG = 1 YES          ######
C####### IMODFG = 0 store fgong files each KFG models             ######
C####### IMODFG = 1 store fgong files when radius is in KFG interval ###
C#######################################################################
C########### NO STORE FGONG FILE FOR MASS LT 0.7mO OR GT 4.0Mo  ########
      IF(STARTMASS.LT.0.7D0.OR.STARTMASS.GT.4.D0)KFGONG=0
C#######################################################################
 47   FORMAT(3E13.4)
      IF(KFGONG.EQ.1) THEN     
      OPEN(UNIT=49,FILE='initial.fgong')
        READ(49,'(2I3)') IMODFG, KFG
         READ(49,'(A10)') STRIFGONG
         if(zini.lt.0.9d-4)then
           kzz=nint(zini*100000.d0)
	   kzz=kzz*10
           FMT = '(I2.2)' ! an integer of width 2 
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
           ZET=TRIM(Z2)//'5'
         elseif(zini.ge.0.9d-4.and.zini.lt.3.5d-4)then
           kzz=nint(zini*10000.d0)
	   kzz=kzz*10
           FMT = '(I2.2)' ! an integer of width 2 
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
	   ZET=TRIM(Z2)//'4'
         elseif(zini.ge.3.5d-4.and.zini.lt.5.5d-4)then
           kzz=nint(zini*100000.d0)
           FMT = '(I2.2)' ! an integer of width 2 
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
	   ZET=TRIM(Z2)//'4'
         elseif(zini.ge.5.5d-4.and.zini.lt.0.9d-3)then
           kzz=nint(zini*10000.d0)
	   kzz=kzz*10
           FMT = '(I2.2)' ! an integer of width 2 
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
	   ZET=TRIM(Z2)//'4'
         elseif(zini.ge.0.9d-3.and.zini.lt.1.3d-3)then
           kzz=nint(zini*1000.d0)
	   kzz=kzz*10
           FMT = '(I2.2)' ! an integer of width 2
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
           ZET=TRIM(Z2)//'3'
         elseif(zini.ge.1.3d-3.and.zini.lt.1.8d-3)then
           kzz=nint(zini*10000.d0)
c	   kzz=kzz*10
           FMT = '(I2.2)' ! an integer of width 2
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
           ZET=TRIM(Z2)//'3'
         elseif(zini.ge.1.8d-3.and.zini.lt.9.d-3)then
           kzz=nint(zini*1000.d0)
	   kzz=kzz*10
           FMT = '(I2.2)' ! an integer of width 2
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
           ZET=TRIM(Z2)//'3'
         elseif(zini.ge.9.d-3.and.zini.lt.1.0d-2)then
           kzz=nint(zini*100.d0)
	   kzz=kzz*10
           FMT = '(I2.2)' ! an integer of width 2
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
           ZET=TRIM(Z2)//'2'
         elseif(zini.ge.1.0d-2)then
           kzz=nint(zini*1000.d0)
c	   kzz=kzz*10
           FMT = '(I2.2)' ! an integer of width 2
           WRITE (Z2,FMT) kzz ! converting integer to string using a 'internal file'
           ZET=TRIM(Z2)//'2'
         endif
	 
	 kel=nint(yini*1000.d0)
         FMT = '(I3.3)' ! an integer of width 3
         WRITE (Y1,FMT) kel ! converting integer to string using a 'internal file'
	 ELIO=TRIM(Y1)
	 if(yini.ge.2.69d-1.and.yini.le.2.71d-1)ELIO='sun'
	 if(ISTART.eq.4.and.
     #      zini.ge.1.72d-2.and.zini.lt.1.722d-2)ELIO='sun'
	 IF(FACTOR.GT.0.D0)THEN
	  STRIFGONG='z'//ZET//'y'//ELIO//'ao'
	 ELSE
          STRIFGONG='z'//ZET//'y'//ELIO//'ae'
         ENDIF
	 MMIN=NINT(STARTMASS*1.d2)
         FMT = '(I4.4)' ! an integer of width 4
         WRITE (X2,FMT) MMIN ! converting integer to string using a 'internal file'


	 IF(ETA.GT.0)then
           IF(INDIFFINI.GT.0)THEN
             NOMEFGONGETAD=TRIM(X2)//STRIFGONG//'de1.'
           ELSE
             NOMEFGONGETA=TRIM(X2)//STRIFGONG//'ce1.'
c          write(*,*)NOMEFGONGETA
           ENDIF
         ELSE 
           IF(INDIFFINI.GT.0)THEN
	    NOMEFGONGD=TRIM(X2)//STRIFGONG//'d.'
           ELSE
            NOMEFGONG=TRIM(X2)//STRIFGONG//'c.'
c          write(*,*)NOMEFGONG
           ENDIF
         ENDIF
	IF(IMODFG.EQ.0) THEN
         CLOSE(49)
         GO TO 1
        ELSE 
         READ(49,*)
         DO I=1,KFG
          READ(49,47)ahlim(1,I),ahelim(1,I),rlim(1,I)
          WRITE(*,47)ahlim(1,I),ahelim(1,I),rlim(1,I)
          READ(49,47)ahlim(2,I),ahelim(2,I),rlim(2,I)
          WRITE(*,47)ahlim(2,I),ahelim(2,I),rlim(2,I)
          READ(49,*)
         ENDDO
         CLOSE(49) 
        END IF 
       END IF
 1     CONTINUE
C#######################################################################
      IRSTA=0    ! added because at this stage it does not enter in RSTAM
      
      CALL STAMPA(DMAS,NABLA,KSD,KSG,KSH,ISTA,IRSTA,INDIFFINI)
      ITEMPOV=0
      DO WHILE (LUCA.LE.NITER)
C**********  TRASFERIMENTI DELLE VARIE QUANTITA NEI MODELLI VECCHI  ****
      write(*,*)fraz
         FLCNOV  = FLCNO
         FL3AV   = FL3A
         FLCARV  = FLCAR
         GENCNOV = GENCNO
         GENHEV  = GENHE
         GENCARV = GENCAR
         DO K=1,MAXMV
            DO J=1,5
               GGV(J,K)=GG(J,K)
            END DO
         END DO
         DO K=1,MAXME
            DO J=1,7
               GG(J,K)=G(J,K)
            END DO
         END DO
         do k=1,maxme
            do j=1,maxne-1
              xxv(j,k)=xxx(j,k)
            end do
         end do     
         MAXMV=MAXME
         TEFFV=TEFF
C##########################################################
         IF(INDIFFINI.GT.0.D0.AND.NMDHEI.EQ.0)THEN 
c          IF(XXX(1,1).LE.0.d0.or.XXX(1,1).GT.XINI*0.99d0)THEN
          IF(XXX(1,1).LE.XINI*0.995d0)THEN
            INDIFF=INDIFFINI
c             INDIFF=indiffini
          ENDIF
          IF(XXX(1,1).LE.0.D0.AND.XXX(3,1).GT.0.d0.AND.G(6,1).GT.0)THEN
            INDIFF=0
            NMDHEI=1
          ENDIF

c          IF(XXX(1,1).LE.0.d0)INDIFF=0
         ENDIF
c         WRITE(*,*)indiff
C##########################################################
         IF(NMD.GE.3.AND.KSE.EQ.0.AND.LUCA.GE.2) CALL OPTIM
          CALL BCONV(NCONV,JBCONV)
C**************  EVOLUZIONE CHIMICA  ***********************************
C********* DIFFUSIONE SE H>0 - INDIFF=1 ********************************
         write(*,*)'entro'
         IF(NMD.GT.4.AND.XXX(1,1).GT.1.D-15.and.tempo.gt.1.d7.
     #      and.tempo.lt.2.d10.AND.INDIFF.EQ.1)CALL DIFFUS
c         IF(NMD.GT.4.AND.XXX(1,1).GT.1.D-3.and.tempo.gt.1.d7.
c     #      and.tempo.lt.2.d10.AND.INDIFF.EQ.1)CALL DIFFUS
         write(*,*)'uscito'
C***********************************************************************       
       DO I=1,MAXNE-1
            XCE(I,1)=XXX(I,1)
         END DO
       if(nmd.gt.2)then
          CALL EVOLUT
          IF(NCONV.GT.0) THEN
         CALL MIXINGALL(NCONV,JBCONV)
          IF(G(6,1).GT.0.D0.AND.XCE(1,1).LE.0.D0) THEN
           CALL MIXINGHE
          ENDIF 
          ENDIF
          IF(NCONV.GT.0)CALL SMALLRAP(NCONV,JBCONV)
         end if
         DO K=1,MAXME
            DO J=1,MAXNE-1
             IF(XXX(J,K).LT.0.D0.AND.DABS(XXX(J,K)).GT.1.D-10)THEN
                  WRITE(*,*)'MUOIO AL MESH',K,' ELEM',J,XXX(J,K)
c                  STOP
               ENDIF 
C#######################################################################
c      AZZERAMENTO ABBONDANZE CHIMICHE SE INFERIORI A...
C#######################################################################               
               IF(XXX(J,K).LT.1.D-30) XXX(J,K)=0.D0
c               IF(XXX(J,K).LT.1.D-15) XXX(J,K)=0.D0
c            if(j.eq.3)then
c              IF(XXX(J,K).LT.1.D-25) XXX(J,K)=0.D0
c            else
c              IF(XXX(J,K).LT.1.D-30) XXX(J,K)=0.D0
c            endif
C#######################################################################             
            END DO
         END DO
         DO I=1,MAXNE-1
            XCE(I,2)=XXX(I,1)
         END DO
         TEMPO=TEMPO+HT1
         EMTOV=EMTOT
c         write(*,*)MAXME,MAXMV
         CALL MASLOS(DMAS)
         EMSOL=EMTOT/(SUNMg33)
         NABLA=-1
         ITERAZ=0
         IPRALL=KSA
         NUMAX= 100
         MACINA=.TRUE.
         COR=0.D0
         VER=0.D0
C***************************** ITERAZIONI ******************************
         DO WHILE(MACINA)
            ITERAZ=ITERAZ+1
            IF(ITERAZ.GT.NUMAX) THEN
C******** SE LE ITERAZIONI SONO TROPPE SI FERMA*************************
               WRITE(2,4444)
               WRITE(*,4444)
               STOP
            ELSE
               CALL QUATM(NABLA)
               CALL HENYEY(ERR,NUM0,IERR)
               IF(ITERAZ.LE.1) THEN
                  WRITE(2,9111)IERR,NUM0,ERR,ITERAZ
                  WRITE(*,9111)IERR,NUM0,ERR,ITERAZ
               ELSE
                  WRITE(2,3456)NUM1,NUM2,DAM,IERR,NUM0,ERR,ITERAZ
                  WRITE(*,3456)NUM1,NUM2,DAM,IERR,NUM0,ERR,ITERAZ
               ENDIF
            END IF
C****************** CONTROLLO CONVERGENZA ******************************
            IF(ITERAZ.GT.1) THEN
               VER=DABS(ERR)
               COR=DABS(DAM)
               IF(VER.LT.5.D-4.AND.COR.LT.5.D-4)MACINA=.FALSE.
               IF(VER.LT.1.D-4.OR. COR.LT.1.D-4)MACINA=.FALSE.
               IF(ITERAZ.GT.6) THEN
                  IF(COR.LE.1.D-3.OR.VER.LE.1.D-3)MACINA=.FALSE.
               END IF
               IF(ITERAZ.GT.9) THEN
                  IF(COR.LE.5.D-2.OR.VER.LE.5.D-2)MACINA=.FALSE.
               END IF
             IF(ISTART.EQ.5.AND.ITERAZ.GT.10) MACINA=.FALSE.
            ENDIF
C*************************************************************************
            IF(MACINA.EQV..TRUE.) THEN
               FRAT=1.D0
               IF(XXX(1,1).LE.0..AND.XXX(3,1).LE.0..AND.COR.GT.1.D+1)
     #         FRAT=0.5D0
               IF(NMD.LE.3.AND.ITERAZ.LT.3)FRAT=.1D00
               CALL FITTA(MAXME,DAM,NUM1,NUM2,FRAT)
            ENDIF
         END DO

C******** CONTROLLA SE RAGGIO E/O PRESSIONE PRESENTANO INVERSIONI********
         MUX=MAXME-1
         K=1
         DO WHILE (G(1,K).LT.G(1,K+1).AND.K.LE.MUX)
            K=K+1
         END DO
         IF(K.LT.MAXME) THEN
            WRITE(2,4414)K
            WRITE(*,4414)K
            STOP
         ENDIF
         K=1
         DO WHILE (G(3,K).GT.G(3,K+1).AND.K.LE.MUX)
            K=K+1
         END DO
         IF(K.LT.MAXME) THEN
            WRITE(2,4415)K
            WRITE(*,4415)K
            STOP
         ENDIF
C***********   MODELLO  CONVERGIUTO  OK  *******************************
         LUCA=LUCA+1                        !INDICE MODELLO INTERNO
         NMD=NMOD+LUCA                      !INDICE MODELLO TOTALE
         IF(ISTART.EQ.4.AND.NNMD.EQ.0)THEN
c           IF(ETA.GT.0.D0)THEN
            IF(NMD.EQ.3)THEN
	     NNMD=NRIP/10000
	     NMD=NRIP-NNMD*10000
	     nmod=nmd-luca
            ENDIF
c           ELSE
c            IF(NMD.EQ.IPRESSSTOP)THEN
c	     NNMD=NRIP/10000
c	     NMD=NRIP-NNMD*10000
c	     nmod=nmd-luca
c            ENDIF
c           ENDIF
	 ENDIF
         IF(NMD.GT.9998) THEN
	    NMOD=100
            LUCA=19
         ENDIF
C***************************** STAMPA RISULTATI ************************
C***********************************************************************
C*************** CONTROLLI PER DECIDERE SE SALVARE I FILE FGONG ********
         IF(XXX(3,1).LE.0.D0.AND.NMDHE.EQ.0)NMDHE=NMD+100
         IF(NMDHE.NE.0.AND.NMD.GE.NMDHE.AND.NMDCO.NE.0)NMDCO=0
	 NMDL=200
	 IF(ISTART.EQ.4)NMDL=55
         IF(KFGONG.EQ.1.AND.NMD.GT.NMDL.AND.NMDCO.NE.0.AND.
     #ISTART.NE.5) THEN
          IGONG=0
C***********************************************************************
C***** Put in the following line IRSTA=2 in case the asympt period 
C***** spacing is needed but not fgong file
C***********************************************************************
C         IRSTA=0
          IRSTA=2   
c          write(*,*)'kfg -->  ',kfg, nmd
          IF(IMODFG.EQ.0) THEN
           IF((NMD-(NMD/KFG)*KFG).EQ.0) THEN ! SAVE FGONG @EACH KFG MODEL
            IRSTA=1
            GO TO 55
           END IF 
          ELSE         
C****** SAVE FGONG IF THE CONDITION ON THE STELLAR RADIUS IS FULFILLED               
           RAG  = STBSUN*DSQRT(10.D0**ELLOG)/((10.D0**TEFF)**2)
           DO I=1,KFG
            IF(AHLIM(1,I).GT.0.D0)THEN
             IF(XXX(1,1).LE.AHLIM(1,I).AND.XXX(1,1).GE.AHLIM(2,I))THEN
              IF(RAG.LE.rlim(1,I).AND.RAG.GE.rlim(2,I)) IRSTA=1
             ENDIF
            ELSE
             IF(XXX(1,1).LE.0.D0.AND.XXX(3,1).LE.AHELIM(1,I).
     #                               AND.XXX(3,1).GE.AHELIM(2,I))THEN
              IF(RAG.LE.rlim(1,I).AND.RAG.GE.rlim(2,I))IRSTA=1
             ENDIF
            ENDIF
            IF(IRSTA.EQ.1) IGONG=1 ! CONTROL FOR TIME STEP REDUCTION ACTIVATED
           END DO
          END IF 
         ELSE
          IRSTA=2
         END IF
 55    CONTINUE
C          IRSTA=2   
       IF(IRSTA.GT.0)THEN
C***********************************************************************
C***     DEFINING THE NAME OF THE OUTPUT FGONG FILE - STRING+NMD   *****
C***********************************************************************
         FMT = '(I5.5)' ! an integer of width 5 with zeros at the left
         I1 = NMD+10000*NNMD
         WRITE (X1,FMT) I1 ! converting integer to string using a 'internal file'
         IF(ETA.GT.0.d0)THEN
           IF(INDIFFINI.GT.0)THEN
            FILENAME=NOMEFGONGETAD//X1
           ELSE
            FILENAME=NOMEFGONGETA//X1
           ENDIF
         ELSE
           IF(INDIFFINI.GT.0)THEN
            FILENAME=NOMEFGONGD//X1
           ELSE
            FILENAME=NOMEFGONG//X1
           ENDIF
         ENDIF
C***********************************************************************
          RTOT = RAG*SUNRcm
          
          rph=(G(1,maxme-1)+G(1,maxme))*1.d10/2.d0   ! SERVE???
C***********************************************************************
C**** COMPUTING THE SECOND DERIVATIVE OF P & RHO AT THE CENTER wrt RADIUS          
C***********************************************************************     
      DO JF=1,3
        TM = G( 4 , JF ) * 1.D+06
        PR = G( 3 , JF ) * 1.D+17
        CALL EOS(JF,1,MAXNE,1,PR,TM,XXX,'MAIN    ',RHO,DAD,CSPA,PMAU,G1,
     #DEL,EMUE,3)
        ROR(JF)=RHO
      ENDDO
      DP1=(G(3,2)-G(3,1))*1.D17/((G(1,2)-G(1,1))*1.D10)
      DP2=(G(3,3)-G(3,2))*1.D17/((G(1,3)-G(1,2))*1.D10)
      DPR=(DP2-DP1)/((G(1,2)-G(1,1))*1.D10)
      DPR=DPR*RTOT**2/(G(3,1)*1.D17)
      
      DRO1=(ROR(2)-ROR(1))/((G(1,2)-G(1,1))*1.D10)
      DRO2=(ROR(3)-ROR(2))/((G(1,3)-G(1,2))*1.D10)
      DRORO=(DRO2-DRO1)/((G(1,2)-G(1,1))*1.D10)
      DRORO=DRORO*RTOT**2/ROR(1)     
C***********************************************************************
C               
Cc       write(*,'(4e16.9)')ELLOG,(10**ELLOG)*sunle,sunle
C        ddf=1.d-2
Cc       IF(RAG.GT.(rlim-ddf).AND.RAG.LT.(rlim+ddf))then
C         WRITE(39,'(i7)')NMD
C         WRITE(39,'(100a)')'Nmesh     R/Rsun          Rph            
C     # Mtot           Ls                Z                Xin        alph
C     #a_ml         d2Pc/dr2         d2Roc/dr2           Age'
C         WRITE(39,'(i5,1p20e16.9)')MAXME-1,RAG,RTOT,
C     #EMTOT*1.d33,10**ELLOG*sunle,zini,xini,alpha,DPR,DRORO,tempo
C         WRITE(41,'(i7)')NMD
C         WRITE(41,'(100a)')'Nmesh     R/Rsun          Rph            
C     # Mtot           Ls                Z                Xin        alph
C     #a_ml         d2Pc/dr2         d2Roc/dr2           Age'
C         WRITE(41,'(i5,1p20e16.9)')MAXME-1,RAG,RTOT,
C     #EMTOT*1.d33,10**ELLOG*sunle,zini,xini,alpha,DPR,DRORO,tempo
C     
C***********************************************************************     
         
c          write(*,*)'IRSTA', IRSTA
          CALL RSTAM(FILENAME,IRSTA,INDIFFINI)
         END IF
C***********************************************************************
C***** END CHANGES FOR STORING DATA IN FGONG FILES                  ****
C***********************************************************************
C***********************************************************************
       IF((NMD-(NMD/KSB)*KSB).EQ.0) NABLA=3
         IF(KSB.EQ.98.AND.(NMD-(NMD/KSB)*KSB).EQ.0) REWIND 02
         CALL STAMPA(DMAS,NABLA,KSD,KSG,KSH,ISTA,IRSTA,INDIFFINI)
C***********************************************************************
C     CONTROLLO PER STOPPARE IL CALCOLO DOPO UN CERTO NUMERO DI 
C     PULSI TERMICI - L'ELIO AL CENTRO DEVE ESSERE ESAURITO
C
C         WRITE(*,'(D9.4,f10.4)')FL3A,ellog
cc     F3ALFA=dlog10(FL3A*(10.d0**ellog))
C      WRITE(*,'(D9.4)')f3alfa
c         fermati=dlog10(10.d0**5*(emsol/2.8d0))
cc         fermati=6.d0
cc     IF(XXX(3,1).LE.0.D0.AND.F3ALFA.GT.fermati)THEN
cc       write(*,*)'luminosita` 3 alfa superiore a 10.d6 Lsun'
cc       write(*,*)'controllo per i pulsi termici'
cc       stop
cc     endif
C
C         WRITE(*,'(D9.4,f10.4)')FL3A,ellog
c      F3ALFA=dlog10(FL3A*(10.d0**ellog))
C      WRITE(*,'(D9.4)')f3alfa
c         fermati=dlog10(10.d0**5*(emsol/2.8d0))
         fermati=70.d0
       IF(NMD.GT.9998) THEN
          NNMD=NNMD+1
       ENDIF
       IF(ISTART.EQ.4.AND.NNMD.GT.3)THEN
         write(*,*)'oltre terzo giro...'
         stop
       ENDIF
       IF(XXX(3,1).LE.0.D0.AND.FL3A.GT.fermati)THEN
         write(*,*)'luminosita` FLUSSO 3 alfa superiore a 70'
         write(*,*)'controllo per i pulsi termici'
         stop
       endif
      CONTROLLO PER STOPPARE IL CALCOLO a Hc=0.3
c       if(xxx(1,1).lt.0.3d0)stop
C****************    PASSO TEMPORALE E PUNCH  **************************
         IF(ISTART.EQ.5) THEN
        HT1=1.0D-4
       ELSE 
          CALL PASTEM(IGONG)
       ENDIF 
c         WRITE(*,*)'========================'
c      WRITE(*,'(2i10,e25.15)')nmd,maxme,G(5,maxme)
c      WRITE(34,'(2i10,e25.15)')nmd,maxme,G(5,maxme)
c      WRITE(*,*)'========================'

c      pause
         IOK=LUCA-(LUCA/NPER)*NPER
         IF(IOK.EQ.0) THEN
            IF (MODELLO.EQ.'modant1') THEN
               MODELLO='modant2'
            ELSE
               MODELLO='modant1'
            END IF
            CALL RWMOD(1,MODELLO,INDIFFINI)
          IF(ISTART.EQ.5) THEN
           ETA=1.0D0
           HT1=1.0D-4
           HT1V=HT1
           TEMPO=1.0D-9  
          ENDIF
         END IF
C=======================================================================
C=============== PUNCH MODELLI DI ZAHB =================================
         IF(ISTART.EQ.5) THEN
       DIFMAS=0.D0
        DO J=1,NTRK
        DIFMAS=DABS(EMSOL-AMZAHB(J))
         IF(DIFMAS.LT.1.0D-10) THEN
         write(*,435) AMZAHB(J)
           KFGONG=1
	   MODELLO=NOME(J)
           TEMPO=TFLASH
           NMDP=NMD
	   NMD=2
             CALL RWMOD(1,MODELLO,INDIFFINI)
           NMD=NMDP
           DIFMAS=0.D0
           ETA=1.0D0
           HT1=1.0D-4
           HT1V=HT1
           TEMPO=1.0D-9
           IF(J.EQ.NTRK) STOP
           GO TO 757
         ENDIF  
          END DO
       ENDIF 
  757  CONTINUE
C***********************************************************************
C        CALL TIMER(ITEMPO)
C        IDELTA=ITEMPO-ITEMPOV
C        ITEMPOV=ITEMPO
C        WRITE(*,*)IDELTA
      END DO             !LOOP PRINCIPALE SU LUCA
      WRITE(2,101)
      WRITE(*,101)
      STOP
C***********************************************************************
  100 FORMAT(10X,'OR VA`, CH`UN SOL VOLERE E` D`AMBEDUE.',/,
     $       10X,'TU DUCA, TU SEGNORE E TU MAESTRO;',/,
     $       10X,'COSI` LI DISSI E POI CHE MOSSO FUE',/,
     $       10X,'ENTRAI PER LO CAMMINO ALTO E SILVESTRO.',/)
  101 FORMAT(10X,'PERCHE GLI OCCHI DELL UOM CERCAN MORENDO IL SOLE')
  111 FORMAT(4I5)
  199 FORMAT(1X,1P,3D12.4)
  434 FORMAT(F7.4,3X,A11)
  435 FORMAT('######## PUNCHO MODELLO ZAHB - MASSA= ',F7.4,' ########')
  444 FORMAT(/,1X,'MASSA=',F7.4,2X,'Y=',F6.3,2X,'Z=',1P,E9.2,2X,'ML=',
     *0P,F5.2,/)
  951 FORMAT(14I5)
 4444 FORMAT(1X,'SUPIN RICADDE E PIU` NON PARVE FORA')
 4414 FORMAT(1X,'INVERSIONE DEL RAGGIO AL MESH',I5)
 4415 FORMAT(1X,'INVERSIONE DELLA PRESSIONE AL MESH',I5)
 4810 FORMAT(80I1)
 9111 FORMAT(1X,'ITERIAMO...',11X,'E(',I1,',',I4,')=',1P,E9.2,I4)
 3456 FORMAT(1X,'DX(',I1,',',I4,')=',1P,E9.2,'  E(',I1,',',I4,')=',
     *1P,E9.2,I4)
C***********************************************************************
      END

      SUBROUTINE DEJA(LUMI,TEFF,HT1,MPUN) 
      IMPLICIT NONE 
      INCLUDE 'consts.2p1'
      REAL*8 LUMI,TEFF,MPUN,HT1 
      REAL*8 A(0:5,0:4),X,Y,LOMP,B,C 
      INTEGER*4 N,I,J 
      DATA  A/6.34916, 3.41678,-1.08683, 0.13095, 0.22427, 0.11968, 
     $       -5.04240, 0.15629, 0.41952,-0.09825, 0.46591, 0.00000, 
     $       -0.83426, 2.96244,-1.37272, 0.13025, 0.00000, 0.00000, 
     $       -1.13925, 0.33659,-1.07493, 0.00000, 0.00000, 0.00000, 
     $       -0.12202, 0.57576, 0.00000, 0.00000, 0.00000, 0.00000/ 
      X=(TEFF-4.05D0)/0.75D0  
      Y=(LUMI-4.60D0)/2.10D0  
      IF(X.GT.+1.D0) X=+1.D0 
      IF(X.LT.-1.D0) X=-1.D0 
      IF(Y.GT.+1.D0) Y=+1.D0 
      IF(Y.LT.-1.D0) Y=-1.D0 
      LOMP=0.D0 
      DO N=0,5 
         DO I=0,N 
            J=N-I 
            IF(J.LE.4) THEN 
               B=DCOS(I*DACOS(X)) 
               C=DCOS(J*DACOS(Y)) 
               LOMP=LOMP+A(I,J)*B*C 
C              WRITE(2,100)N,I,J,A(I,J),B,C,LOMP 
            ENDIF 
         END DO 
      END DO 
      MPUN=-(SUNMg33)*HT1*10.D0**(-LOMP)
      RETURN 
C 100 FORMAT(3I4,4F9.5) 
      END 

      SUBROUTINE FATO(ELL,TEF,PCEN,TCEN,AS,AT,AU,AV,IPP,INDU) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
C ********************************************************************** 
      COMMON/DE/EL0,EL1,EL2,P0,P1,P2,R0,R1,R2,T0,T1,T2,DPC,DTC,DELL,DTEF 
     *,EL0P,EL1P,EL2P,P0P,P1P,P2P,R0P,R1P,R2P,T0P,T1P,T2P,CORZ 
C ********************************************************************** 
  547 FORMAT(7F8.4,2I4,2I2) 
  789 FORMAT(///) 
  555 FORMAT(15X,'DL',7X,'DTE',7X,'DPC',7X,'DTC') 
  333 FORMAT(11X,10E10.3) 
      INDU=0 
      WRITE(2,555) 
      WRITE(*,555) 
      AA=(EL2-EL0)/DPC 
      AB=(R2-R0)/DPC 
      AC=(P2-P0)/DPC 
      AD=(T2-T0)/DPC 
      IF(ISTART.EQ.2.AND.IPP.EQ.0)GO TO 50 
      AE=(EL1-EL0)/DTC 
      AF=(R1-R0)/DTC 
      AG=(P1-P0)/DTC 
      AH=(T1-T0)/DTC 
   50 CONTINUE 
      AI=(EL2P-EL0P)/DELL 
      AL=(R2P-R0P)/DELL 
      AM=(P2P-P0P)/DELL 
      AN=(T2P-T0P)/DELL 
      AO=(EL1P-EL0P)/DTEF 
      AP=(R1P-R0P)/DTEF 
      AQ=(P1P-P0P)/DTEF 
      AR=(T1P-T0P)/DTEF 
      IF(ISTART.EQ.2.AND.IPP.EQ.0)GO TO 51 
      ALF1=AL-AF*AI/AE 
      ALF2=AM-AG*AL/AF 
      ALF3=AN-AH*AM/AG 
      BET1=AP-AF*AO/AE 
      BET2=AQ-AG*AP/AF 
      BET3=AR-AH*AQ/AG 
      GAM1=AF*AA/AE-AB 
      GAM2=AG*AB/AF-AC 
      GAM3=AH*AC/AG-AD 
      ETA1=AT-AF*AS/AE 
      ETA2=AU-AG*AT/AF 
      ETA3=AV-AH*AU/AG 
      ZOT1=ALF2-GAM2*ALF1/GAM1 
      ZOT2=BET2-GAM2*BET1/GAM1 
      ZOT3=ETA2-GAM2*ETA1/GAM1 
      ZOT4=ALF3-GAM3*ALF2/GAM2 
      ZOT5=BET3-GAM3*BET2/GAM2 
      ZOT6=ETA3-GAM3*ETA2/GAM2 
      DISC1=ZOT1*ZOT5-ZOT2*ZOT4 
      ERREL=(ZOT5*ZOT3-ZOT2*ZOT6)/DISC1 
      ERRTE=(ZOT1*ZOT6-ZOT3*ZOT4)/DISC1 
      ERRPC=ETA1/GAM1-ALF1*ERREL/GAM1-BET1*ERRTE/GAM1 
      ERRTC=AI*ERREL/AE+AO*ERRTE/AE-AA*ERRPC/AE-AS/AE 
      WRITE(2,333)ERREL,ERRTE,ERRPC,ERRTC 
      WRITE(*,333)ERREL,ERRTE,ERRPC,ERRTC 
      ZUM1=ERREL/ELL 
      ZUM2=ERRTE/TEF 
      ZUM3=ERRPC/PCEN 
      ZUM4=ERRTC/TCEN 
      WRITE(2,333)ZUM1,ZUM2,ZUM3,ZUM4 
      WRITE(*,333)ZUM1,ZUM2,ZUM3,ZUM4 
      WRITE(2,789) 
      IF(DABS(ZUM1).GT.2.)GO TO 48 
      IF(DABS(ZUM2).GT.2.)GO TO 48 
      IF(DABS(ZUM3).GT.2.)GO TO 48 
      IF(DABS(ZUM4).GT.2.)GO TO 48 
      IF(DABS(ERREL)-.10*ELL)40,40,41 
   41 ERREL=.10*ELL*(ERREL/DABS(ERREL)) 
   40 ELL=ELL+ERREL*CORZ 
      IF(DABS(ERRTE)-.10*TEF)42,42,43 
   43 ERRTE=.10*TEF*(ERRTE/DABS(ERRTE)) 
   42 TEF=TEF+ERRTE*CORZ 
      IF(DABS(ERRPC)-.10*PCEN)44,44,45 
   45 ERRPC=.10*PCEN*(ERRPC/DABS(ERRPC)) 
   44 PCEN=PCEN+ERRPC*CORZ 
      IF(DABS(ERRTC)-.10*TCEN)46,46,47 
   47 ERRTC=.10*TCEN*(ERRTC/DABS(ERRTC)) 
   46 TCEN=TCEN+ERRTC*CORZ 
      ELLOG=DLOG10(ELL/SUNLe)
      PPL=DLOG10(PCEN) 
      PTL=DLOG10(TCEN) 
      TEFF=DLOG10(TEF) 
      TEFF=TEFF 
      EMMU=EMTOT/SUNMg33
      WRITE(2,547)EMMU,ELLOG,PPL,PTL,TEFF 
      WRITE(*,547)EMMU,ELLOG,PPL,PTL,TEFF 
   49 RETURN 
   48 INDU=1 
      GO TO 49 
   51 CONTINUE 
      FAT1=AC*AL/AB-AM 
      FAT2=AC*AP/AB-AQ 
      FAT3=AC*AT/AB-AU 
      FAT4=AD*AL/AB-AN 
      FAT5=AD*AP/AB-AR 
      FAT6=AD*AT/AB-AV 
      DISC2=FAT1*FAT5-FAT2*FAT4 
      ERREL=(FAT3*FAT5-FAT2*FAT6)/DISC2 
      ERRTE=(FAT1*FAT6-FAT3*FAT4)/DISC2 
      ERRPC=AL*ERREL/AB+AP*ERRTE/AB-AT/AB 
      ZUM1=ERREL/ELL 
      ZUM2=ERRTE/TEF 
      ZUM3=ERRPC/PCEN 
      WRITE(2,333)ZUM1,ZUM2,ZUM3 
      WRITE(*,333)ZUM1,ZUM2,ZUM3 
      WRITE(2,789) 
      IF(DABS(ZUM1).GT.2)GO TO 67 
      IF(DABS(ZUM2).GT.2)GO TO 67 
      IF(DABS(ZUM3).GT.2)GO TO 67 
      IF(DABS(ERREL)-.10*ELL)60,60,61 
   61 ERREL=.10*ELL*(ERREL/DABS(ERREL)) 
   60 ELL=ELL+ERREL*CORZ 
      IF(DABS(ERRTE)-.10*TEF)62,62,63 
   63 ERRTE=.10*TEF*(ERRTE/DABS(ERRTE)) 
   62 TEF=TEF+ERRTE*CORZ 
      IF(DABS(ERRPC)-.10*PCEN)64,64,65 
   65 ERRPC=.10*PCEN*(ERRPC/DABS(ERRPC)) 
   64 PCEN=PCEN+ERRPC*CORZ 
      ELLOG=DLOG10(ELL/SUNLe)
      PPL=DLOG10(PCEN) 
      TEFF=DLOG10(TEF) 
      PTL=DLOG10(TCEN) 
      TEFF=TEFF 
      EMMU=EMTOT/SUNMg33
      WRITE(2,547)EMMU,ELLOG,PPL,PTL,TEFF 
      WRITE(*,547)EMMU,ELLOG,PPL,PTL,TEFF 
   66 RETURN 
   67 INDU=1 
      GO TO 66 
      END 
c---------- Versione per Santi
c      subroutine gtime(igio,imes,ianno)
cc***********************************************************************
cc      Fornisce giorno, mese, anno
cc***********************************************************************
cc
c      character*3 cmes,mesi(12)
c      data mesi/'Jan','Feb','Mar','Apr','May','Jun','Jul',
c     *      'Aug','Sep','Oct','Nov','Dec'/
cc
c      imes=0
c      call system('date > ccc')
c      open (99,file='ccc',status='unknown')
c      read (99,10)cmes,igio,ianno
c        write(*,10)cmes,igio,ianno
c      close(99)
c      call system('rm ccc')
c 10   format(4x,a3,1x,i2,15x,i4)
cc======================================================================
c      do i=1,12
c      imes=imes+1
c      if(cmes.eq.mesi(i)) go to 30
c      end do
cc
c30    return
c      end
c---------------------------------------------------------------------
c---------versione per Adriano
	subroutine gtime(igio,imes,ianno)
c***********************************************************************
c      Fornisce giorno, mese, anno
c***********************************************************************
c
	character*3 cmes,mesi(12)
	data mesi/'jan','feb','mar','apr','may','jun','jul',
     *	'aug','sep','oct','nov','dec'/
c
        imes=0
	call system('date > ccc')
	open (99,file='ccc',status='unknown')
	read (99,10)igio,cmes,ianno
c        write(*,10)igio,cmes,ianno
	close(99)
	call system('rm ccc')
 10	format(4x,i2,1x,a3,14x,i4)
C***** example of current format output: "    12 dec              2017"
c======================================================================
	do i=1,12
	imes=imes+1
	if(cmes.eq.mesi(i)) go to 30
	end do
c
30	return
	end
              
      SUBROUTINE HBREZO 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
C ********************************************************************** 
      COMMON/HORIZ/AMCO,YCO,CCO 
C ********************************************************************** 
      DIMENSION NMESH(6),CC(MAXNE) 
      DATA NMESH/20,10,5,5,3,2/ 
      DO K=1,MAXNE-1 
         CC(K)=XXX(K,MAXME-1) 
      END DO 
      CORE=AMCO*SUNMg
      DO 1 K=1,MAXME 
         IF(G(5,K).GE.CORE)GO TO 2 
    1 CONTINUE 
    2 IF(K.GT.(MAXME-3))GO TO 15 
      RAT=CORE/G(5,K) 
      KCOR=K 
      DO J=1,K 
         G(5,J)=G(5,J)*RAT 
      END DO 
      IND=0 
      DO 12 ITER=1,2 
      DO 11 I=1,3 
      IND=IND+1 
      DM=G(5,K+1)-G(5,K) 
      IF(I-2)4,5,6 
    4 NME=NMESH(IND) 
      GO TO 7 
    5 NME=NMESH(IND) 
      GO TO 7 
    6 NME=NMESH(IND) 
    7 DM=DM/ FLOAT(NME+1) 
      KK=K+1 
      NITER=MAXME-K 
      DO N=KK,MAXME 
         K1=N+NME 
         GG(5,K1)=G(5,N) 
      END DO 
      DO N=1,NME 
         K1=K+N 
         G(5,K1)=G(5,K)+DM* FLOAT(N) 
      END DO 
      MAXME=MAXME+NME 
      DO N=1,NITER 
         K1=K+NME+N 
         G(5,K1)=GG(5,K1) 
      END DO 
   11 K=K+NME+1 
      K=KCOR 
   12 CONTINUE 
      DO 13 L=1,K 
      XXX(1,L)=0.D0 
      XXX(2,L)=0.D0 
      XXX(3,L)=YCO 
   13 XXX(4,L)=CCO 
      K2=K+1 
      DO 14 L=K2,MAXME 
      DO 14 I=1,MAXNE-1 
   14    XXX(I,L)=CC(I) 
      IF(IPRALL.NE.1)RETURN 
      WRITE(2,222)(G(5,NUM),NUM=1,MAXME) 
      WRITE(2,222)(XXX(1,NUM),NUM=1,MAXME)
C      WRITE(2,222)(PMAU(1,NUM),NUM=1,MAXME)
  222 FORMAT(1P,13E10.3) 
      RETURN 
   15 WRITE(2,111) 
  111 FORMAT(12X,'MASSA DI CORE TROPPO GRANDE, MUOIO DISPERATO',//) 
      STOP 
      END

      SUBROUTINE HENYEY( ERR , INUM , IERR ) ! Henyey method
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      COMMON/FREMEM/ALF(4,LIM),BET(4,LIM),GAM(4,LIM),DG(5,LIM),
     &              O((MAXNE+1-20)*LIM),PV(LIM),TV(LIM),ENGRAV(LIM)
      COMMON/SYSTEMA/ALF1(4),BET1(4),GAM1(4),E(4),B(4,4),C(4,4)
      DIMENSION ESP(3) , CAP(3) , DRAD(3) , DAD(3)
      DIMENSION EPSN(3),EPNEU(3),GRAVI(3),RHO(3),PMAU(3),CSPA(3),GRSAD
     #(3)
      DIMENSION CICCIO(2)
      ERR  = 0.D0
      IERR = 0
      INUM = 0
      DO K = 1 , 4
         ALF(K,1) = 0.D0
         BET(K,1) = 0.D0
         GAM(K,1) = 0.D0
      END DO
      BET(3,1) = 1.D0
      GAM(4,1) = 1.D0
      DO K=1,MAXME
         KK=K+1
         DO J=1,5
            DG(J,K)=G(J,KK)-G(J,K)
         END DO
      END DO
      DO JF=1,MAXME-1
         JFF=JF+1
         EMAS=(G(5,JF)+G(5,JFF))/2.D0
         L=2
         DO WHILE(EMAS.GT.GG(5,L))
            L=L+1
         END DO
         RAT=(EMAS-GG(5,L-1))/(GG(5,L)-GG(5,L-1))
         PV(JF)=GG(3,L-1)+RAT*(GG(3,L)-GG(3,L-1))
         TV(JF)=GG(4,L-1)+RAT*(GG(4,L)-GG(4,L-1))
      END DO
C ************************* LOOP PRINCIPALE ****************************
      MAX=MAXME-1
      DO JF=1,MAX
         JFF=JF+1
C ****** QUANTITA' FISICHE DEFINITI NEL PUNTO JF + 1/2 *****************
C ****** midpoint mean of zone (JF, JF+1)
         RR = ( G( 1 , JF ) + G( 1 , JFF ) ) / 2.D0
         ZL = ( G( 2 , JF ) + G( 2 , JFF ) ) / 2.D0
         PP = ( G( 3 , JF ) + G( 3 , JFF ) ) / 2.D0
         TT = ( G( 4 , JF ) + G( 4 , JFF ) ) / 2.D0
         EM = ( G( 5 , JF ) + G( 5 , JFF ) ) / 2.D0
C * CHIMICA DEL MESH JF E' QUELLA DELLA REGIONE COMPRESA TRA JF E JF+1 *
C **** Chemistry of mesh JF is that of the region between JF and JF+1
C ******************** DEFINIZIONE DELLE DERIVATE **********************
         DELP = 1.D-06*PP
         DELT = 1.D-06*TT
         DERL = 1.D-06*ZL
         DERR = 1.D-06*RR
C ************ CALCOLO DERIVATE DELLE VARIE QUANTITA' FISICHE **********
C *** calculate derivative of various physical quantities
         TM = TT*1.D+06
         PR = PP*1.D+17
c      write(*,*)'hen',tm,pr
        CALL EOS(JF,1,MAXNE,3,PR,TM,XXX,'HENYEY  ',RHO,DAD,CSPA,PMAU,G1,
     #DEL,EMUE,3)
         CALL NKAPPA(JF,1,MAXNE,3,RHO,TM,XXX,CAP)
       GI=GRAVC0*1.D+13*EM/(RR*RR)
C ***** GI = magnitude of grav. accn. = GM/(R*R)
         DO I=1,3
            IF(I.EQ.2) THEN
               PP=PP+DELP
            ELSEIF(I.EQ.3) THEN
               TT=TT+DELT
               PP=PP-DELP
            ENDIF
            EPSN(I) = 0.D0
            EPNEU(I)= 0.D0
            GRAVI(I)= 0.D0
            TM=TT*1.D+06
            PR=PP*1.D+17
            CALL EPSI(RHO(I),TM,JF,EPSN(I),V1,V2,V3,V4,I)
            CALL NEUTR(RHO(I),TM,XXX(1,JF),EPNEU(I))
            CALL EPSIG(PP,TT,CSPA(I),GRAVI(I),DAD(I),PV(JF),TV(JF),JF)
          ESP(I)=EPSN(I)+EPNEU(I)+GRAVI(I)
          DRAD(I)=CRAD*CAP(I)*PP*ZL/(EM*(TT**4))
            GRSAD(I)=DAD(I)
            IF(DRAD(1).GE.DAD(1)) THEN
               IF(DRAD(I).GE.DAD(I)) THEN
               TVA=TM
                     CALL EOS(JF,1,MAXNE,1,PR,TVA,XXX,'HEN-SUPE',
     &                           O1,O2,O3,PKM,G1,DEL,EMUE,3)
                     Q=PKM
                     IF(Q.LT.0.D0) Q=0.D0
                     IF(I.EQ.1) QQQ=Q
                  CALL NSUPERA(RHO(I),PR,TM,CAP(I),Q,CSPA(I),DAD(I),GI,
     $                      ALPHA,GRSAD(I),O1,DRAD(I))
               ELSE
                  GRSAD(I)=DAD(I)
                  IF(I.EQ.1) QQQ=1.D0
               ENDIF
            ENDIF
         END DO
         TT=TT-DELT
         CICCIO(1)=0.D0
         CICCIO(2)=0.D0
         TM=TT*1.D+06
         PR=PP*1.D+17
         IF(DRAD(1).GE.DAD(1)) THEN
            Q=QQQ
            DO JJ=1,2
               IF(JJ.EQ.1) THEN
                  DRADA=DRAD(1)*(1.D0+DERL/ZL)
               ELSE
                  DRADA=DRAD(1)
                  RR=RR+DERR
               ENDIF
             GI=GRAVC0*1.D+13*EM/(RR*RR)
               CALL NSUPERA(RHO(1),PR,TM,CAP(1),Q,CSPA(1),DAD(1),GI,
     $                   ALPHA,CICCIO(JJ),O1,DRADA)
            END DO
            RR=RR-DERR
         END IF
C ********************** FINE CALCOLO DERIVATE *************************
         DDELT = 2.D0*DELT
         DDELP = 2.D0*DELP
         DDELL = 2.D0*DERL
         DDERR = 2.D0*DERR
         DEPT  = (ESP(3)-ESP(1))/DDELT
         DEPP  = (ESP(2)-ESP(1))/DDELP
         DROT  = (RHO(3)-RHO(1))/DDELT
         DROP  = (RHO(2)-RHO(1))/DDELP
         DGRRT = (DRAD(3)-DRAD(1))/DDELT
         DGRRP = (DRAD(2)-DRAD(1))/DDELP
       DGRRL = (CRAD/2.D0)*CAP(1)*PP/(EM*(TT**4))
C ** SCEGLIE IL GRADIENTE IN BASE AL CRITERIO DI SCHWARTZCHILD *********
c ** CHOOSE THE GRADIENT AT THE BASE ACCORDING TO THE SCHWARTZSCHILD
c ** CRITERION
         G(6,JF)=DRAD(1)-DAD(1)
         VAGRA=G(6,JF)
         IF(VAGRA.LT.0.D0) THEN
            GRAD = DRAD(1)
            DGRT = DGRRT
            DGRP = DGRRP
            DGRL = DGRRL
            DGRR = 0.D0
         ELSE
            GRAD = GRSAD(1)
            DGRT = (GRSAD(3)-GRSAD(1))/DDELT
            DGRP = (GRSAD(2)-GRSAD(1))/DDELP
            DGRL = (CICCIO(1)-GRSAD(1))/DDELL
            DGRR = (CICCIO(2)-GRSAD(1))/DDERR
         ENDIF
         G(7,JF)=GRAD
         ENGRAV(JF)=GRAVI(1)
C ********************** CALCOLA LE DERIVATE ***************************
         DDROP =  DROP / RHO( 1 )
         DDROT =  DROT / RHO( 1 )
         RR1   =  1.D0 / RR
         DG1   =  1.D0 / DG( 1 , JF )
         DG3   =  1.D0 / DG( 3 , JF )
         DG10P = -10.D0 * DG( 5 , JF ) * DEPP
         DG10T = -10.D0 * DG( 5 , JF ) * DEPT
         IF(JF.EQ.1) THEN    ! SVILUPPI CENTRALI
            B( 1 , 1 ) =  0.D0
            C( 1 , 1 ) =  2.D0
            B( 2 , 1 ) =  0.D0
            C( 2 , 1 ) =  0.D0
            B( 3 , 1 ) =  G( 3 , JF  ) * (  DG3 + 2.d0*DDROP )
            C( 3 , 1 ) =  G( 3 , JFF ) * ( -DG3 + 2.d0*DDROP )
            B( 4 , 1 ) =  G( 4 , JF  ) * (  2.D0*DDROT )
            C( 4 , 1 ) =  G( 4 , JFF ) * (  2.D0*DDROT )
            B( 1 , 2 ) =  0.D0
            C( 1 , 2 ) =  3.D0
            B( 2 , 2 ) =  0.D0
            C( 2 , 2 ) =  0.D0
            B( 3 , 2 ) =  G( 3 , JF  ) * DDROP
            C( 3 , 2 ) =  G( 3 , JFF ) * DDROP
            B( 4 , 2 ) =  G( 4 , JF  ) * DDROT
            C( 4 , 2 ) =  G( 4 , JFF ) * DDROT
         ELSE
            B( 1 , 1 ) =  G( 1 , JF  ) * (  RR1 + DG1   )
            C( 1 , 1 ) =  G( 1 , JFF ) * (  RR1 - DG1   )
            B( 2 , 1 ) =  0.D0
            C( 2 , 1 ) =  0.D0
            B( 3 , 1 ) =  G( 3 , JF  ) * ( -DG3 - DDROP )
            C( 3 , 1 ) =  G( 3 , JFF ) * (  DG3 - DDROP )
            B( 4 , 1 ) =  G( 4 , JF  ) * ( -DDROT )
            C( 4 , 1 ) =  G( 4 , JFF ) * ( -DDROT )
            B( 1 , 2 ) =  G( 1 , JF  ) * (  RR1 - DG1   )
            C( 1 , 2 ) =  G( 1 , JFF ) * (  RR1 + DG1   )
            B( 2 , 2 ) =  0.D0
            C( 2 , 2 ) =  0.D0
            B( 3 , 2 ) =  G( 3 , JF  ) * DDROP
            C( 3 , 2 ) =  G( 3 , JFF ) * DDROP
            B( 4 , 2 ) =  G( 4 , JF  ) * DDROT
            C( 4 , 2 ) =  G( 4 , JFF ) * DDROT
         END IF
c           B( 1 , 1 ) =  G( 1 , JF  ) * (  RR1 + DG1   )
c            C( 1 , 1 ) =  G( 1 , JFF ) * (  RR1 - DG1   )
c            B( 2 , 1 ) =  0.D0
c            C( 2 , 1 ) =  0.D0
c            B( 3 , 1 ) =  G( 3 , JF  ) * ( -DG3 - DDROP )
c            C( 3 , 1 ) =  G( 3 , JFF ) * (  DG3 - DDROP )
c            B( 4 , 1 ) =  G( 4 , JF  ) * ( -DDROT )
c            C( 4 , 1 ) =  G( 4 , JFF ) * ( -DDROT )
c            B( 1 , 2 ) =  G( 1 , JF  ) * (  RR1 - DG1   )
c            C( 1 , 2 ) =  G( 1 , JFF ) * (  RR1 + DG1   )
c            B( 2 , 2 ) =  0.D0
c            C( 2 , 2 ) =  0.D0
c            B( 3 , 2 ) =  G( 3 , JF  ) * DDROP
c            C( 3 , 2 ) =  G( 3 , JFF ) * DDROP
c            B( 4 , 2 ) =  G( 4 , JF  ) * DDROT
c            C( 4 , 2 ) =  G( 4 , JFF ) * DDROT
         B( 1 , 3 ) =  0.D0
         C( 1 , 3 ) =  0.D0
         B( 2 , 3 ) = -G( 2 , JF  )
         C( 2 , 3 ) =  G( 2 , JFF )
         B( 3 , 3 ) =  G( 3 , JF  ) * DG10P
         C( 3 , 3 ) =  G( 3 , JFF ) * DG10P
         B( 4 , 3 ) =  G( 4 , JF  ) * DG10T
         C( 4 , 3 ) =  G( 4 , JFF ) * DG10T
         B( 1 , 4 ) = -G( 1 , JF  ) * DGRR
         C( 1 , 4 ) = -G( 1 , JFF ) * DGRR
         B( 2 , 4 ) = -G( 2 , JF  ) * DGRL
         C( 2 , 4 ) = -G( 2 , JFF ) * DGRL
         B(3,4)=G(3,JF )*(DG(4,JF)*(.5D0*DG(3,JF)+PP)/(TT*DG(3,JF)**2)-
     $       DGRP)
         C(3,4)=G(3,JFF)*(DG(4,JF)*(.5D0*DG(3,JF)-PP)/(TT*DG(3,JF)**2)-
     $       DGRP)
         B(4,4)=G(4,JF )*(PP*(-.5D0*DG(4,JF)-TT)/(DG(3,JF)*TT*TT)-DGRT)
         C(4,4)=G(4,JFF)*(PP*(-.5D0*DG(4,JF)+TT)/(DG(3,JF)*TT*TT)-DGRT)
C ******************* VERIFICA EQUAZIONI *******************************
         IF(JF.LE.1) THEN       ! SVILUPPI CENTRALI
           E(1)=DLOG(DABS(2.D0*4.D0*PI*GRAVC0*1.D05/3.d0*RHO(1)**2*
     #          G(1,2)**2/DG(3,1)))   
c ** E(1) = LOG(|2*(4*PI/3)*10^5*GRAVCO*((RHO(1)*G(1,2))^2)/DG(3,1)|)
         E(2)=DLOG(DABS(4.D0*PI*1.D-03/3.d0*G(1,2)**3*RHO(1)/G(5,2)))
         ELSE
        E(1)=DLOG(DABS(DG(3,JF)*RR*RR/(DG(1,JF)*
     #         GRAVC0*1.D06*EM*RHO(1))))
        E(2)=DLOG(DABS(4.D0*PI*1.D-03*RR*RR*RHO(1)*DG(1,JF)/DG(5,JF)))
         ENDIF
C=======================================================================
c          E(1)=DLOG(DABS(DG(3,JF)*RR*RR/(DG(1,JF)*
c     #         GRAVC0*1.D06*EM*RHO(1))))
c         E(2)=DLOG(DABS(4.D0*PI*1.D-03*RR*RR*RHO(1)*DG(1,JF)/DG(5,JF)))
         E(3)=DG( 2 , JF ) - ESP(1) * 10.D0 * DG( 5 , JF )
         E(4)=DG( 4 , JF ) * PP / ( DG( 3 , JF ) * TT ) - GRAD
         DO N = 1 , 4
            ALF1(N) = -E(N) - B(1,N) * ALF(1,JF) - B(2,N) * ALF(2,JF)
            BET1(N) = -B(1,N) * BET(1,JF) - B(2,N) * BET(2,JF) - B(3,N)
            GAM1(N) = -B(1,N) * GAM(1,JF) - B(2,N) * GAM(2,JF) - B(4,N)
         END DO
         IF(E(3).EQ.0.D0.AND.DG(2,JF).EQ.0.D0.AND.ESP(1).EQ.0.D0) THEN
           E(3)=0.D0
         ELSE
           E(3)=E(3)/DSQRT(DG(2,JF)**2+100.D0*ESP(1)*ESP(1)*DG(5,JF)**2)
         ENDIF
        E(4)=E(4)/DSQRT(PP*PP*DG(4,JF)**2/(TT*TT*DG(3,JF)**2)+GRAD*GRAD)
C *********** CALCOLA QUANTO SONO VERIFICATE LE EQUAZIONI **************
         DO N = 1 , 4
            IF( DABS( E(N) ) . GT . DABS( ERR ) ) THEN
               ERR  = E( N )
               IERR = N
               INUM = JF
            ENDIF
         END DO
C ************** IN CASO DI DEBUG STAMPA VARIE QUANTITA' ***************
         IF(IPRALL.NE.0) THEN
            WRITE(2,110) JF
            WRITE(2,90)RR,ZL,PP,TT,EM,DEPP,DEPT,DROP,DROT,DGRL,DGRP,DGRT
      WRITE(2,90)(EPSN(NN),NN=1,3),(EPNEU(NN),NN=1,3),(GRAVI(NN),NN=1,3)
         WRITE(2,90)(RHO(NN),NN=1,3),(CAP(NN),NN=1,3),(DRAD(NN),NN=1,3),
     &           (ESP(NN),NN=1,3)
         WRITE(2,90)G(1,JF),G(1,JFF),DG(1,JF),G(2,JF),G(2,JFF),DG(2,JF),
     &            G(3,JF),G(3,JFF),DG(3,JF),G(4,JF),G(4,JFF),DG(4,JF)
C           WRITE(2,90)(XX(IP),IP=1,6),G(5,JF),G(5,JFF),DG(5,JF)
            WRITE(2,90)E(1),E(2),E(3),E(4)
C           WRITE(2,100)C(1,1),C(3,1),C(4,1),B(1,1),B(3,1),B(4,1),E(1)
C           WRITE(2,100)C(1,2),C(3,2),C(4,2),B(1,2),B(3,2),B(4,2),E(2)
C           WRITE(2,100)C(2,3),C(3,3),C(4,3),B(2,3),B(3,3),B(4,3),E(3)
C           WRITE(2,100)C(2,4),C(3,4),C(4,4),B(2,4),B(3,4),B(4,4),E(4)
         ENDIF
         CALL RESNUC(JFF)
      END DO
      G(6,MAXME)=G(6,MAXME-1)
      G(7,MAXME)=G(7,MAXME-1)
C ********************** FINE LOOP PRINCIPALE **************************
C*********************** END OF MAIN LOOP
      RETURN
   90 FORMAT(1X,1P,12E10.3)
  100 FORMAT(1X,1P,E10.3,'DX',E10.3,'DP',E10.3,'DT=', E10.3,'DXV',
     *E10.3,'DPV',E10.3,'DTV',E10.3)
  110 FORMAT(/,1X,'MESH N.',I4)
      END

      SUBROUTINE FITTA(M,DAM,NUM1,NUM2,FRAT) 
      IMPLICIT REAL*8(A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'atmosfer.2p1'
      INCLUDE 'consts.2p1'
      COMMON/FREMEM/ALF(4,LIM),BET(4,LIM),GAM(4,LIM),DG(5,LIM),
     #       O((MAXNE+1-17)*LIM)
      DIMENSION DA(4),DAV(4) 
      DDR=(DLOG(RB)-DLOG(G(1,M))) 
      DDP=(DLOG(PB)-DLOG(G(3,M))) 
      DDT=(DLOG(TB)-DLOG(G(4,M))) 
      DEN=1.D0-BET(2,M)*DPL-GAM(2,M)*DTL 
      TAO=(BET(2,M)*DPT+GAM(2,M)*DTT)/DEN 
      DAO=(ALF(2,M)+BET(2,M)*DDP+GAM(2,M)*DDT)/DEN 
      CAO=BET(1,M)*(DPL*TAO+DPT)+GAM(1,M)*(DTL*TAO+DTT)-DRL*TAO-DRT 
      DTEF=-ALF(1,M)-DAO*(BET(1,M)*DPL+GAM(1,M)*DTL-DRL)+DDR-BET(1,M)* 
     #DDP-GAM(1,M)*DDT 
      DTEF=DTEF/(CAO) 
      DEL1=TAO*DTEF+DAO 
      TEFF=TEFF+DTEF*FRAT/COND
      DA(3)=DDP+DPL*DEL1+DPT*DTEF 
      DA(4)=DDT+DTL*DEL1+DTT*DTEF 
      DAM=0.D0 
      DO K=1,M 
         L=M-K+1 
         DA(1)=ALF(1,L)+BET(1,L)*DA(3)+GAM(1,L)*DA(4) 
         DA(2)=ALF(2,L)+BET(2,L)*DA(3)+GAM(2,L)*DA(4) 
         DAV3=DA(3) 
         DAV4=DA(4) 
         DO I=1,4 
            IF(DABS(DA(I)).GT.DABS(DAM)) THEN 
               DAM=DA(I)*1.00000000001D0 
               NUM1=I 
               NUM2=L 
            ENDIF 
         END DO 
         DO J=1,4 
            DA(J)=DA(J)*G(J,L) 
            G(J,L)=G(J,L)+DA(J)*FRAT 
            IF(L.NE.M) DG(J,L)=DG(J,L)+DAV(J)*FRAT-DA(J)*FRAT 
         END DO 
         IF(IPRALL.EQ.1) THEN 
            WRITE(2,112) L,(G(IT,L),IT=1,4),(DA(IT),IT=1,4) 
            WRITE(2,112) L,(DG(IT,L),IT=1,4) 
         ENDIF 
         DO I=1,4 
            DAV(I)=DA(I) 
         END DO 
         DA(3)=ALF(3,L)+BET(3,L)*DAV3+GAM(3,L)*DAV4 
         DA(4)=ALF(4,L)+BET(4,L)*DAV3+GAM(4,L)*DAV4 
      END DO 
      ELLOG=DLOG10(G(2,M)/(SUNLe*1.D-32))
      RETURN 
  112 FORMAT(1X,I5,1P,8E13.6) 
      END 

      SUBROUTINE INNES(ITCC,IPP)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'atmosfer.2p1'
      INCLUDE 'consts.2p1'
      DIMENSION ZN(280),ZNV(280),PV(LIM),TV(LIM)
c      DIMENSION ZN(180),ZNV(180),PV(LIM),TV(LIM)
C **********************************************************************
      COMMON/HORIZ/AMCO,YCO,CCO
      COMMON/DE/EL0,EL1,EL2,P0,P1,P2,R0,R1,R2,T0,T1,T2,DPC,DTC,DELL,DTEF
     *,EL0P,EL1P,EL2P,P0P,P1P,P2P,R0P,R1P,R2P,T0P,T1P,T2P,CORZ
C **********************************************************************
      DATA ZN/0.D+0,1.28D-7,5.10D-7,1.540D-6,
     #4.90D-6,1.50D-5,4.6D-5,.000133,.00035,.0008
     *,.0013,.002,.0025,.003,.004,.005,.006,.007,.008,.009,.01,.011,.012
     *,.013,.014,.015,.016,.017,.018,.019,.02,.021,.022,.023,.024,.025,
     *.026,.027,.028,.029,.03,.032,.034,.036,.038,.04,.042,.044,.046,.
     *048,.05,.052,.054,.056,.058,.06,.062,.064,.066,.068,.072,.076,.08,
     *.084,.088,.092,.096,.1,.104,.108,.112,.116,.12,.124,.128,.132,.136
     *,.14,.144,.148,.152,.156,.16,.164,.168,.172,.176,.18,.184,.188,.19
     *2,.196,.2,.205,.21,.215,.22,.225,.23,.234,.242,.25,.26,.27,.28,.29
     *,.30,.3075,.315,.33,.345,.36,.375,.39,.405,.42,.435,.45,.465,.48,.
     *49,.5,.515,.53,.54,.55,.56,.57,.58,.59,.6,.61,.62,.63,.64,.65,.66,
     *.67,.68,.69,.7,.71,.72,.73,.74,.75,.76,.77,.78,.79,.8,.81,.82,.83,
     *.84,.85,.855,.86,.865,.87,.875,.88,.885,.89,.895,.9,.905,.91,.915,
     *.92,.9225,.925,.9275,.93,.932,.934,.936,.938,.939,.94,
     *.941,.942,.943,.944,.945,.946,.947,.948,.949,.950,
     *.951,.952,.953,.954,.955,.956,.957,.958,.959,.960,
     *.961,.962,.963,.964,.965,.966,.967,.968,.969,.970,
     *.971,.972,.973,.974,.975,.976,.977,.978,.979,.980,
     *.981,.982,.983,.984,.985,.986,.987,.988,.989,.990,
     *.9902,.9904,.9906,.9908,.9910,.9912,.9914,.9916,.9918,.9920,
     *.9922,.9924,.9926,.9928,.9930,.9932,.9934,.9936,.9938,.9940,
     *.9942,.9944,.9946,.9948,.9950,.9952,.9954,.9956,.9958,.9960,
     *.9962,.9964,.9966,.9968,.9970,.9972,.9974,.9976,.9978,.9980,
     *.9981,.9982,.9983,.9984,.9985,.9986,.9987,.9988,.9989,.9990D0/
      DATA ELLPC/0.D0/,ELLTC/0.D0/,EMMU/0.D0/,ZNV/280*0.D0/
c      DATA VER/1.D-10/ ORIGINALE
      DATA VER/1.D-9/
C **********************************************************************
c  555 FORMAT(5D10.3)
  555 FORMAT(1D10.3,f10.7,3D10.3) ! for solar calibration
  888 FORMAT(/,1X,'CONVERGENZA RAGGIUNTA',/) ! convergence reached
  665 FORMAT(5F8.4)
c  666 FORMAT(2D10.3,5F8.4,I4,F5.2)
  666 FORMAT(f11.9,D10.3,5F8.4,I4,F5.2)    ! for solar calibration
  333 FORMAT(///)
  111 FORMAT(1X,'INTERNO',6X,'R',E12.5,3X,'L',E12.5,3X,'P',E12.5,3X,'T'
     *,E12.5)
  222 FORMAT(1X,'ESTERNO',6X,'R',E12.5,3X,'L',E12.5,3X,'P',E12.5,3X,'T'
     *,E12.5)
  444 FORMAT(5X,'M',7X,'L',6X,'TE',5X,'Y',7X,'Z',8X,'P',7X,'T',
     *6X,'R',5X,'RO',/)
C *********************************************************************
      MAXME=280
      MAXMV=280
C ********MAXME E MAXMV PER INNES MODIFICATA****************************

      IF(ISUB.EQ.1) WRITE(*,444)
      IF(ITCC.NE.1) THEN ! read in from modstart
         READ(4,555)FARZ,ALPHA,AMCO,ETA,CCO
         FRAZ=1.D0-FARZ/1.D+02
         READ(4,666)YINI,ZINI,EMMU,ELLOG,TEFF,ELLPC,ELLTC,LAST,CORZ
         STARTMASS=EMMU
         IF(CORZ.LE..09)CORZ=1.D0
         DO K=1,MAXME
            ZNV(K)=ZN(K)*FRAZ/.999D0  
         END DO
         XINI=1.D0-YINI-ZINI-ABSOL(2)
         IF(ISTART.EQ.4.OR.ISTART.EQ.5)YCO=1.D0-CCO-ZINI
       IF(ISTART.EQ.3.OR.ISTART.EQ.5)THEN
         ABSOL(2)=0.D0
         ABSOL(7)=0.D0
         ABSOL(8)=0.D0
       ELSEif(ISTART.EQ.2)THEN
         ABSOL(7)=ABSOL(7)/ZINI
         ABSOL(8)=ABSOL(8)/ZINI
       ENDIF
      END IF
      IF (IPP.EQ.1) THEN
         DO L=1,MAXME-1
            PV(L)=(GG(3,L)+GG(3,L+1))/2.D0
            TV(L)=(GG(4,L)+GG(4,L+1))/2.D0
         END DO
      ENDIF
      DO KM=1,MAXME
         XXX(1,KM)=XINI     !-ABSOL(7)*ZINI
         XXX(2,KM)=ABSOL(2)
         XXX(3,KM)=YINI
      END DO
      SUMABSOL=0.D0-ABSOL(7)-absol(8)
      DO J=4,MAXNE-2
        SUMABSOL=SUMABSOL+ABSOL(J)
c     write(*,*)j,SUMABSOL
      ENDDO
      ABSOL(MAXNE-1)=1.D0-SUMABSOL  !SCARICO SUL FERRO PER NORMALIZZARE
c      write(*,*)SUMABSOL,absol(26)
c      pause
      DO KC=4,MAXNE-1
         DO KM=1,MAXME
            XXX(KC,KM)=ZINI*ABSOL(KC)
         END DO
      END DO
c      write(*,*)XXX(26,1)
c      pause
      TEF=10.D0**TEFF
      ELL=(10.D0**ELLOG)*SUNLe
      EMTOT=EMMU*(SUNMg33)
      EMTAT=EMMU*SUNMg
      PCEN=10.D0**ELLPC
      TCEN=10.D0**ELLTC
      DO K=1,MAXME
         G(5,K)=EMTAT*ZNV(K)
      END DO 
      IF(ISUB.EQ.1) THEN
         CALL QUATM(-1)
         WRITE(*,*)'CALCOLO SINGOLA ATMOSFERA TERMINATO'
         WRITE(2,*)'CALCOLO SINGOLA ATMOSFERA TERMINATO'
         STOP
      ENDIF
      IF(ISTART.EQ.4.OR.ISTART.EQ.5)CALL HBREZO
      IF(ITCC.NE.0) THEN
         IF(IPP.EQ.0)TCEN=TCEN+CCO
         IF(IPP.EQ.1)PCEN=G(3,1)
         ELLPC=DLOG10(PCEN)
         ELLTC=DLOG10(TCEN)
      ENDIF
      WRITE(2,665)EMMU,ELLOG,ELLPC,ELLTC,TEFF
      WRITE(*,665)EMMU,ELLOG,ELLPC,ELLTC,TEFF
      DELL= .0005D0 * ELL
      DTEF= .0005D0 * TEF
      DPC = .0005D0 * PCEN
      DTC = .0005D0 * TCEN
  1   CONTINUE
      DO KI=1,3
         IF(KI.EQ.1) THEN
            TTC=TCEN
            PPC=PCEN
            VLS=ELL
            TES=TEF
         ELSE IF(KI.EQ.2) THEN
            TTC=TCEN+DTC
            TES=TEF+DTEF
         ELSE IF(KI.EQ.3) THEN
            TTC=TCEN
            PPC=PCEN+DPC
            TES=TEF
            VLS=ELL+DELL
         END IF
         CALL VEIOVE(VLS,PPC,TTC,EM,R,ELLE,P,T,PV,TV,1,IPP,KI,LAST)
         G(6,1)=G(6,2)
         GG(6,1)=GG(6,2)
         IF(P.LE.0.D0.OR.T.LE.0.D0) CALL NEWBVAL
         IF(KI.EQ.1) THEN
            EL0=ELLE
            R0=R
            P0=P
            T0=T
            WRITE(2,111)R,ELLE,P,T
            WRITE(*,111)R,ELLE,P,T
         ELSE IF(KI.EQ.2) THEN
            EL1=ELLE
            R1=R
            P1=P
            T1=T
         ELSE IF(KI.EQ.3) THEN
            EL2=ELLE
            R2=R
            P2=P
            T2=T
         ENDIF
       ELLOG=DLOG10(VLS/SUNLe)
         TEFF=DLOG10(TES)
         CALL QUATM(-1)
         R=RB*1.D+10
         P=PB*1.D+17
         T=TB*1.D+06
         ELLE=VLS
         EM=EMTAT*FRAZ
         IF(KI.EQ.1) THEN
            G(1,1)=0.D0
            G(2,1)=0.D0
            G(3,1)=PCEN
            G(4,1)=TCEN
            G(1,MAXME)=R
            G(2,MAXME)=VLS
            G(3,MAXME)=P
            G(4,MAXME)=T
         ENDIF
         CALL VEIOVE(VLS,PPC,TTC,EM,R,ELLE,P,T,PV,TV,2,IPP,KI,LAST)
         G(6,MAXME)=G(6,MAXME-1)
         GG(6,MAXME)=GG(6,MAXME-1)
         IF(R.LE.0.D0.OR.ELLE.LE.0.D0) CALL NEWBVAL
         IF(KI.EQ.1) THEN
            EL0P=ELLE
            R0P=R
            P0P=P
            T0P=T
            WRITE(2,222)R,ELLE,P,T
            WRITE(*,222)R,ELLE,P,T
            WRITE(2,333)
         ELSEIF (KI.EQ.2) THEN
            EL1P=ELLE
            R1P=R
            P1P=P
            T1P=T
         ELSEIF (KI.EQ.3) THEN
            EL2P=ELLE
            R2P=R
            P2P=P
            T2P=T
         ENDIF
         IF(KI.EQ.1) THEN
            AS=(EL0-EL0P)
            AT=(R0-R0P)
            AU=(P0-P0P)
            AV=(T0-T0P)
            ASS=DABS(AS/EL0)
            ATS=DABS(AT/R0)
            AUS=DABS(AU/P0)
            AVS=DABS(AV/T0)
            IF(ISTART.EQ.2.AND.IPP.EQ.0) ASS=0.D0
            IF(ASS.LT.VER.AND.ATS.LT.VER.AND.AUS.LT.VER.AND.AVS.LT.VER)
     $      THEN
               WRITE(2,888)
               WRITE(*,888)
             ELLOG=DLOG10(G(2,MAXME)/SUNLe)
               TEFF=DLOG10(TEF)
               RETURN
            ENDIF
         ENDIF
      END DO
      CALL FATO(ELL,TEF,PCEN,TCEN,AS,AT,AU,AV,IPP,INDU)
      IF(INDU.EQ.1) CALL NEWBVAL
      GO TO 1
      END
      SUBROUTINE MASLOS(DM) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'

      DIMENSION RATM(2000),PATM(2000),TATM(2000),VMATM(2000),X(MAXNE) 
      DIMENSION GGA(5,lim)
      common/NOSUB/IPRESSSTART,IPRESSSTOP,IPRESSSWITCH

      IPIPPA=1 
      LOBOR=-1
      IF(ETA.LE.0.D0) RETURN 
      IF(xxx(1,1).gt.XINI*0.9850d0.and.startmass.gt.2.d0) RETURN 
      IF(xxx(1,1).le.0.d0.and.xxx(3,1).gt.0.d0.and.
     #   xxx(3,1).lt.0.9d0.and.startmass.gt.1.8d0) RETURN
      IF(xxx(1,1).le.0.d0.and.xxx(3,1).gt.0.9d0.and.
     #     FL3A.gt.1.d-3.and.g(6,1).le.0.d0) RETURN
 
cc===== modifica per eliminare subatmosfera
c      IF(nmd.LE.IPRESSSTOP+55.and.IPRESSSWITCH.eq.0)then
      IF(nmd.LE.100.and.IPRESSSWITCH.eq.0)then
       RETURN 
      ENDIF
       IF(HT1.LT.2.5d1.AND.ISTART.NE.5)RETURN  
       IF(HT1.LT.2.5d2.AND.XXX(3,1).LT.0.5D0.AND.
     #XXX(3,1).GT.0.D0)RETURN 
       IF(XXX(3,1).GT.0.D0.AND.XXX(3,1).LT.5.D-4.AND.
     #    STARTMASS.GT.1.8d0)RETURN 
       dmaxme=2.d0*dfloat(maxme-maxmv)/dfloat(maxme)
c       write(*,*)'dmaxme',dmaxme
       if(dmaxme.le.-0.08d0.AND.ISTART.NE.5)then
c        pause
        return
       endif
c       if(maxmv-maxme.gt.300)return
c      if (nmd.lt.120.and.istart.eq.3) then
c       return
c      endif
c      i=1
c      DO WHILE ((G(5,i).lt.(0.99*EMTOT)).and.(i.lt.MAXME-30))
c       i=i+1
c      enddo
c      ISTEP=MAXME-i-1
c
C=======================================================================
C=== ACCRESCIMENTO PER CREAZIONE MODELLI DI ZAHB =======================
      IF(ISTART.EQ.5) THEN
        IF (IPRESSSWITCH.ne.0) then  
         DM=0.0010D0*SUNMg33
         if(nmd.le.10)dm=0.0001D0*SUNMg33
         EMTOT=EMTOT+DM
         FINMAS=EMTOT*FRAZ
         RAT=(FINMAS-G(5,MAXME-30))/(G(5,MAXME)-G(5,MAXME-30))
         LA=MAXME-29
          DO L=LA,MAXME
            G(5,L)=G(5,MAXME-30)+RAT*(G(5,L)-G(5,MAXME-30))
          END DO
        ELSE

          DM=5.D-4*SUNMg33
          EMTOT=EMTOT+DM
          FINMAS=EMTOT
          RAT=(FINMAS-G(5,MAXME-ISTEP))/(G(5,MAXME)-G(5,MAXME-ISTEP))
          LA=MAXME-ISTEP+1
          DO L=LA,MAXME
            G(5,L)=G(5,MAXME-ISTEP)+RAT*(G(5,L)-G(5,MAXME-ISTEP))
          END DO
          if (G(5,LA).lt.(0.99*EMTOT).and.ISTEP.gt.30) ISTEP=ISTEP-10 
        END IF

        RETURN
      ENDIF
c===== fine modifica per eliminare subatmosfera
C=======================================================================     
c      write(*,'(e13.5)')dm
C**********************************************************************
C         MASS-LOSS BY REIMERS (1975) IN Mo/anno
C
      DM=-(4.D-13)*REIMAS*HT1*ETA*(10.D0**(ELLOG*1.5D0))/(EMTOT*
     #                           (10.D0**(2.D0*TEFF)))
      DM=DM*SUNMg33  !*SUNMg33 PERCHE IN USCITA DIVIDO PER SUNMg33 !!!!
C**********************************************************************      
      EMTOT=EMTOT+DM 
c===== modifica per eliminare subatmosfera
      if (IPRESSSWITCH.eq.0) goto 660
c===== fine modifica per eliminare subatmosfera

      FINMAS=EMTOT*FRAZ
      IPIPPA=0 
      IF(IPIPPA.EQ.0) THEN 
C***************************************************************** 
         DO K=1,MAXNE-1 
            X(K)=XXX(K,MAXME-1) 
         END DO 
         BOTT=1.D0-FRAZ 
       EMSOLA=EMTOT/(SUNMg33)
         CALL NATMOS(ELLOG,TEFF,EMSOLA,X,ZINI,BOTT,ALPHA, 
     $               RATM,PATM,TATM,VMATM,NATM,0,0,1,1,MAXNE,nmd) 
         PB=PATM(NATM)*1.D-17 
         TB=TATM(NATM)*1.D-6 
         RB=RATM(NATM) 
C***************************************************************** 
         K=1 
         DO WHILE (G(5,K).LT.FINMAS.AND.K.LT.MAXME) 
            K=K+1 
         END DO 
         KMAX=K 
         IF(DABS((FINMAS-G(5,K-1))/G(5,K-1)).LT.1.D-8)KMAX=K-1 
         G(1,KMAX)=RB 
         G(3,KMAX)=PB 
         G(4,KMAX)=TB 
         G(5,KMAX)=FINMAS 
         MAXME=KMAX 
         LA=MAXME-11 
C= Santi = added the following definition of lobor to avoid problem in pastem =
         LOBOR=MAXME-11 
C============================================================================== 
         LU=LA+1 
         DO LL=1,4 
            IF(LL.NE.2) THEN 
               RAT=(G(LL,MAXME)-G(LL,LA))/(G(5,MAXME)-G(5,LA)) 
               DO L=LU,MAXME 
                  G(LL,L)=G(LL,LA)+RAT*(G(5,L)-G(5,LA)) 
               END DO 
            ENDIF 
         END DO 
      ELSE 
         LA=LOBOR 
         LU=LA+1 
         RAT=(FINMAS-G(5,LA))/(G(5,MAXME)-G(5,LA)) 
         DO L=LU,MAXME 
            G(5,L)=G(5,LA)+RAT*(G(5,L)-G(5,LA)) 
         END DO 
      ENDIF 
        RETURN
c===== modifica per eliminare subatmosfera
 660   CONTINUE 
      FINMAS=EMTOT*FRAZ
      
      K=1 
      DMLIM=5.d0*DM
      IF(startmass.gt.1.8d0.AND.xxx(1,1).LE.0.d0)DMLIM=5.d0*DM
C originale      DO WHILE (G(5,K).LT.FINMAS+4.d0*DM.AND.K.LT.MAXME)
      DO WHILE (G(5,K).LT.FINMAS+DMLIM.AND.K.LT.MAXME)
        K=K+1 
      END DO 
      KMAX=K
      
      IF(FINMAS-G(5,KMAX).LT.1.D-8)THEN
        write(*,*)'finmas-g5 < 1.d8'
        K=1 
        DO WHILE (G(5,K).LT.FINMAS.AND.K.LT.MAXME) 
           K=K+1 
        END DO 
        KMAX=K 
        DO K=1,MAXNE-1 
          X(K)=XXX(K,MAXME-1) 
        END DO 
        BOTT=1.D0-FRAZ 
        EMSOLA=EMTOT/(SUNMg33)
        CALL NATMOS(ELLOG,TEFF,EMSOLA,X,ZINI,BOTT,ALPHA, 
     $               RATM,PATM,TATM,VMATM,NATM,0,0,1,1,MAXNE,nmd) 
        PB=PATM(NATM)*1.D-17 
        TB=TATM(NATM)*1.D-6 
        RB=RATM(NATM) 
        MAXME=KMAX 
	lobor=kmax
	RETURN
      ENDIF
      LOBOR=KMAX
      RAT=(EMTOT*FRAZ-G(5,KMAX))/(G(5,MAXME)-G(5,KMAX))
      DO K=KMAX+1,MAXME
        G(5,K)=G(5,KMAX)+RAT*(G(5,K)-G(5,KMAX))
        GG(5,K)=G(5,K)    ! um Lverlust zu ermoeglichen
      END DO
c===== fine modifica per eliminare subatmosfera

      RETURN 
      END
      SUBROUTINE NEWBVAL 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      INCLUDE 'consts.2p1'
      WRITE(*,*)'QUALCOSA E MINORE DI ZERO NELLA INNES'
      STOP 
      END 
      SUBROUTINE PLOTTA(RTOT) 
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'consts.2p1' 
      REAL*8 VA 
      CHARACTER*1 A(51,120),SIMB(10),PNUM(11),BLANK,STAR,SII,POINT 
      INCLUDE 'maincom.2p1'
      DIMENSION VAR(10),VARM(10),IVAV(10) 
      DATA BLANK/' '/,STAR/'*'/,SII/'I'/,POINT/'.'/ 
      DATA SIMB/'R','L','P','T','H','3','4','C','N','O'/ 
      DATA PNUM/'0','9','8','7','6','5','4','3','2','1','0'/ 
  111 FORMAT(120A1) 
  222 FORMAT(11X,'I---------I---------I---------I---------I---------I--- 
     *------I---------I---------I---------I---------I') 
  333 FORMAT(10X,'0.0',7X,'0.1',7X,'0.2',7X,'0.3',7X,'0.4',7X,'0.5',7X,'
     *0.6',7X,'0.7',7X,'0.8',7X,'0.9',7X,'1.0',3X,'M/MTOT') 
  444 FORMAT(//////) 
  555 FORMAT(1H1) 
  666 FORMAT(//,1X,'R=',1PE9.2,2X,'L=',E9.2,2X,'P=',E9.2,2X,'T=',E9.2,2X 
     *,'H=',E9.2,2X,'He3=',E9.2,2X,'He4=',E9.2,2X,'C=',E9.2,2X,'N=',E9.2,2X, 
     *'O=',E9.2) 
      WRITE(2,555) 
      DO 20 K=1,51 
      DO 22 L=1,120 
   22 A(K,L)=BLANK 
   20 A(K,12)=SII 
      N=1 
      DO 21 K=1,11 
      A(N,7)=PNUM(1) 
      A(N,8)=POINT 
      A(N,9)=PNUM(K) 
   21 N=N+5 
      A(1,7)=PNUM(10) 
      VARM(1)=RTOT 
      DO 7 J=2,4 
      VARM(J)=G(J,1) 
      DO 8 JJ=2,MAXME 
      PAP=G(J,JJ) 
      VARM(J)=DMAX1(PAP,VARM(J)) 
    8 CONTINUE 
    7 CONTINUE 
      DO 9 J=1,6 
      JJJ=J+4 
      VARM(JJJ)=XXX(J,1) 
      DO 10 JJ=2,MAXME 
      VARM(JJJ)=DMAX1(VARM(JJJ),XXX(J,JJ)) 
   10 CONTINUE 
    9 CONTINUE 
      N=0 
      DO 1 K=1,101 
      EM=.01* FLOAT(K-1) 
      EM=EM*EMTOT 
      IF(EM.GT.G(5,MAXME))GO TO 1 
      N=N+1 
      DO 2 L=2,MAXME 
      IF(EM-G(5,L))3,3,2 
    2 CONTINUE 
    3 LL=L-1 
      IF(LL.GT.(MAXME-1))LL=MAXME-1 
      RAT=(EM-G(5,LL))/(G(5,LL+1)-G(5,LL)) 
      DO 4 M=1,4 
    4 VAR(M)=G(M,LL)+RAT*(G(M,LL+1)-G(M,LL)) 
      DO 5 M=1,6 
      MM=M+4 
    5 VAR(MM)=XXX(M,LL)+RAT*(XXX(M,LL+1)-XXX(M,LL)) 
      MAC=K+11 
      DO 15 M=1,10 
      IF(VARM(M).EQ.0) THEN 
      VA=50.26 
      ELSE 
      VA=(VAR(M)/VARM(M))*50.+.26 
      END IF 
      LIVA=VA 
      IVA=51-(LIVA) 
      IF(IVA.LT.1)IVA=1 
      IF(M.EQ.1)GO TO 12 
      MU=M-1 
      DO 11 MM=1,MU 
      IF(IVA.EQ.IVAV(MM))GO TO 14 
   11 CONTINUE 
   12 A(IVA,MAC)=SIMB(M) 
      GO TO 13 
   14 A(IVA,MAC)=STAR 
   13 IVAV(M)=IVA 
   15 CONTINUE 
    1 CONTINUE 
      N=N+11 
      WRITE(2,333) 
      WRITE(2,222) 
      DO 16 M=1,51 
   16 WRITE(2,111)(A(M,K),K=1,N) 
      WRITE(2,222) 
      WRITE(2,333) 
      VARM(1)=VARM(1)*1.E 10   
      VARM(2)=VARM(2)*1.E 32 
      VARM(3)=VARM(3)*1.E 17 
      VARM(4)=VARM(4)*1.E 06   
      WRITE(2,666)(VARM(K),K=1,10) 
      WRITE(2,444) 
      RETURN 
      END
      SUBROUTINE RESNUC(JFF) 
      IMPLICIT REAL*8(A-H,O-Z) 
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      COMMON/FREMEM/ALF(4,LIM),BET(4,LIM),GAM(4,LIM),DG(5,LIM),
     #              O((MAXNE+1-17)*LIM)
      COMMON/SYSTEMA/ALF1(4),BET1(4),GAM1(4),E(4),B(4,4),C(4,4) 
      DATA  IPRQUI/0/ 
      AA=C(1,1)/C(1,2) 
      AB=C(1,4)/C(1,2) 
      DD=C(2,4)/C(2,3) 
      FF=((C(3,4)-AB*C(3,2))-DD*C(3,3))/(C(3,1)-AA*C(3,2)) 
      H=(C(4,4)-AB*C(4,2)-DD*C(4,3))-FF*(C(4,1)-AA*C(4,2)) 
      V1=((ALF1(4)-AB*ALF1(2)-DD*ALF1(3))-FF*(ALF1(1)-AA*ALF1(2)))/H 
      V2=((BET1(4)-AB*BET1(2)-DD*BET1(3))-FF*(BET1(1)-AA*BET1(2)))/H 
      V3=((GAM1(4)-AB*GAM1(2)-DD*GAM1(3))-FF*(GAM1(1)-AA*GAM1(2)))/H 
      H=C(3,1)-AA*C(3,2) 
      V4=((ALF1(1)-AA*ALF1(2))-V1*(C(4,1)-AA*C(4,2)))/H 
      V5=((BET1(1)-AA*BET1(2))-V2*(C(4,1)-AA*C(4,2)))/H 
      V6=((GAM1(1)-AA*GAM1(2))-V3*(C(4,1)-AA*C(4,2)))/H 
      V7=(ALF1(3)-C(3,3)*V4-C(4,3)*V1)/C(2,3) 
      V8=(BET1(3)-C(3,3)*V5-C(4,3)*V2)/C(2,3) 
      V9=(GAM1(3)-C(3,3)*V6-C(4,3)*V3)/C(2,3) 
      V10=(ALF1(2)-C(3,2)*V4-C(4,2)*V1)/C(1,2) 
      V11=(BET1(2)-C(3,2)*V5-C(4,2)*V2)/C(1,2) 
      V12=(GAM1(2)-C(3,2)*V6-C(4,2)*V3)/C(1,2) 
      H=V3*V5-V6*V2 
      ALF(3,JFF)=-(V3*V4-V6*V1)/H 
      BET(3,JFF)=V3/H 
      GAM(3,JFF)=-V6/H 
      ALF(4,JFF)=-(V1+V2*ALF(3,JFF))/V3 
      BET(4,JFF)=-V2*BET(3,JFF)/V3 
      GAM(4,JFF)=(1.D0-V2*GAM(3,JFF))/V3 
      ALF(1,JFF)=V10+V11*ALF(3,JFF)+V12*ALF(4,JFF) 
      BET(1,JFF)=V11*BET(3,JFF)+V12*BET(4,JFF) 
      GAM(1,JFF)=V11*GAM(3,JFF)+V12*GAM(4,JFF) 
      ALF(2,JFF)=V7+V8*ALF(3,JFF)+V9*ALF(4,JFF) 
      BET(2,JFF)=V8*BET(3,JFF)+V9*BET(4,JFF) 
      GAM(2,JFF)=V8*GAM(3,JFF)+V9*GAM(4,JFF) 
      IF(IPRALL.EQ.0.AND.IPRQUI.EQ.0)RETURN 
      WRITE(2,111)ALF(1,JFF),BET(1,JFF),GAM(1,JFF),ALF(2,JFF),BET(2,JFF) 
     &,GAM(2,JFF),ALF(3,JFF),BET(3,JFF),GAM(3,JFF),ALF(4,JFF),BET(4,JFF) 
     &,GAM(4,JFF) 
  111 FORMAT(1X,1P,12E10.3) 
      RETURN 
      END           
      SUBROUTINE RWMOD(IRW,MODELLO,INDIFFINI) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER*11 MODELLO 
      INCLUDE 'atmosfer.2p1'
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      COMMON/REWRMOD/NMOD,KSA,KSB,KSC,KSD,KSE,KSG,KSH,KSW,NITER,NPER,
     #KFGONG
      common/NOSUB/IPRESSSTART,IPRESSSTOP,IPRESSSWITCH

C*********************************************************************** 
  100 FORMAT(4I5) 
  104 FORMAT(4I5,f8.4,2I3,I6,I3) 
  101 FORMAT(13I5) 
  102 FORMAT(1P,4D22.15)  ! FORMATO ESTESO 
  103 FORMAT(/,1X,'MODELLO PUNCHATO',/) 
C*********************************************************************** 
      IF(IRW.EQ.0) THEN 
C***************************** LETTURA MODANT ************************** 
         READ(4,104)MAXME,MAXMV,NMOD,KBORD,PROGM,IIST,NNMD,NRIP,
     #               IPRESSSWITCH
         READ(4,102)(((U(J,K,L),L=1,4),K=1,4),J=1,3), 
     $   (((V(J,K,L),L=1,4),K=1,4),J=1,3),(EL(L),L=1,4),(TE(L),L=1,4), 
     $   ELLOG,ETA,EMTOT,EMTOV,TEMPO,XINI,ZINI,ALPHA,HT1,HT1V,TEFF,TEFFV 
     $   ,PM1,PM2,FRAZ,FL3A,FLCNO,VEMCO,ABSOL(MAXNE-1) 
         READ(4,102)((G(J,K),J=1,7),K=1,MAXME) 
         READ(4,102)((XXX(J,K),J=1,MAXNE-1),K=1,MAXME) 
         READ(4,102)((GG(J,K),J=1,7),K=1,MAXMV) 
	 STARTMASS=PROGM
C*********************************************************************** 
      ELSE 
C***************************** PERFORA ********************************* 
C*********************************************************************** 
C************* INIZIALIZZAZIONE VARIABILI DA CALCOLO ZAHB ************** 
         IHJK=1
         IF(ISTART.EQ.5) THEN
          ETA=0.D0
c         ETA=0.2D0
          HT1=1.0D1
          HT1V=HT1
          IHJK=4
         ENDIF
	 IF(ISTART.EQ.4)IHJK=4
C*********************************************************************** 
         OPEN(UNIT=7,FILE=MODELLO) 
         REWIND(7) 
         IMOD=NMD-1 
         WRITE(7,100) IHJK,NITER,NPER,ISUB
         WRITE(7,101) KSA,KSB,KSC,KSD,KSE,KSF,KSG,KSH,KSW,ISINGOLA,NSHUT
     #                ,INDIFFINI,KFGONG
         WRITE(7,104)MAXME,MAXMV,IMOD,KBORD,PROGM,IIST,NNMD,NRIP,
     #               IPRESSSWITCH
         WRITE(7,102)(((U(J,K,L),L=1,4),K=1,4),J=1,3), 
     $   (((V(J,K,L),L=1,4),K=1,4),J=1,3),(EL(J),J=1,4),(TE(J),J=1,4), 
     $   ELLOG,ETA,EMTOT,EMTOV,TEMPO,XINI,ZINI,ALPHA,HT1,HT1V,TEFF,TEFFV 
     $   ,PM1,PM2,FRAZ,FL3A,FLCNO,VEMCO,ABSOL(MAXNE-1) 
         WRITE(7,102)((G(J,K),J=1,7),K=1,MAXME) 
         WRITE(7,102)((XXX(J,K),J=1,MAXNE-1),K=1,MAXME) 
         WRITE(7,102)((GG(J,K),J=1,7),K=1,MAXMV) 
         WRITE(2,103) 
         WRITE(*,103) 
         CLOSE(7) 
C*********************************************************************** 
      ENDIF 
C***************************** ESCE ************************************ 
      RETURN 
      END
      SUBROUTINE SKIPFASE (T,RO,YF,KK,IFASE)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE'nuclear.2p1'
      INCLUDE 'consts.2p1'
      DIMENSION YF(MAXME)

      IFASE=0
c ********* SKIPPA FASI INUTILI *********
      IF(T. LT. 1.D+05) RETURN
c      IF(T. LT. 8.D+06) RETURN
      IF(XXX(1,1).LT.XINI*(1.D0-1.D-4).and.T.LT.1.D+06.
     #AND.NMD.GE.40) RETURN
c      IF(T. LT. 1.D+06. AND. NMD. GE. 40) RETURN
      IF(YF(1). LE. 0.D0. AND. YF(3). GT. 0.D0. AND. T. LE. 4.D7) RETURN
      IF(YF(1). LE. 0.D0. AND. YF(3). LE. 0.D0. AND. T. LE. 2.D9) RETURN
      IF(T. GT. 2.D09) THEN
         WRITE(*,150)
         WRITE(2,150)
         STOP
      ENDIF
C ******************* RICERCA FASE EVOLUTIVA: **************************
C ** IFASE=1 -> COMB. IDR. IFASE=2 -> COMB. ELIO  IFASE=3 -> AVANZATE **
      IF(T.LE.1.D08) THEN
         IF(YF(1).GT.0.D0) THEN
            IFASE=1
            YF(MAXNE)=1.D0/RO  
         ELSE
            IFASE=2
            YF(MAXNE)=YF(3)**2*RO/6.D0    ! 3-alfa
         ENDIF
      ELSE IF(T.LE.4.D08) THEN
         IF(YF(3).GT.0.D0) THEN
            IFASE=2
            YF(MAXNE)=YF(3)**2*RO/6.D0    ! 3-alfa
         ELSE
            IFASE=3
            YF(MAXNE)=1.D0/RO         ! fotodis. del neon
            YF(1)=CARXY(1,KK)
            YF(3)=CARXY(2,KK)/ATW(3)
         ENDIF
      ELSE IF(T.GT.4.D08) THEN
         IFASE=3
         YF(MAXNE)=1.D0/RO         ! fotodis. del neon
      ENDIF
C***********************************************************************
  150 FORMAT(10X,'NON PIU` ANDRAI FARFALLONE AMOROSO',/,
     $       10X,'NOTTE E GIORNO D`INTORNO GIRANDO',/,
     $       10X,'DELLE BELLE TURBANDO IL RIPOSO',/,
     $       10X,'NARCISETTO ADONCINO D`AMOR!')
      RETURN
      END

      SUBROUTINE STAMPA(DMAS,NABLA,KSD,IFVA,IFMO,ISTA,IRSTA,INDIFFINI)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      INCLUDE 'consts.2p1'
      DIMENSION A(8),AA(4),BE(12),CW(11,50),STAR(11,LIM),
C-FABIO     &DD(6,10),FF(8),HH(8),CV(9),RPQ(9),ES(20)
     &DD(6,50),FF(8),HH(8),CV(9),RPQ(9),ES(20)
      DIMENSION EE(14),DFIS(4),DFIX(4),IVME(4),IMINT(50),IMEST(50)
      DIMENSION DRAD(LIM),XX(MAXNE),xmassa(lim)
      DIMENSION CROS(MAXPR,2),YF(MAXNE)
      CHARACTER*6  EVFASE
      CHARACTER*3  DIF
      COMMON/FREMEM/DG(5,LIM),O((MAXNE+1-5-1)*LIM),ENGRAV(LIM)
      COMMON/FLUX/FLUSSO(20)
      COMMON/OVER/FACTOR
      DATA INEU/0/,ICHI/0/,IPICCI/1/
      DATA MEMHEB/0/
      
      DO K=1,50
         IMINT(K)=0
         IMEST(K)=0
         DO J=1,11
            CW(J,K)=0.D0
            IF(J.LE.6) DD(J,K)=0.D0
         END DO
      END DO
      DO 51 K=1,12
         BE( K ) = 0.D0
         BB( K ) = 0.D0
         IF( K . GT . 9 ) GO TO 51
         CV( K ) = 0.D0
         RPQ(K ) = 0.D0
         IF( K . GT . 8 ) GO TO 51
         FF( K ) = 0.D0
         HH( K ) = 0.D0
         IF( K . GT . 4 ) GO TO 51
         AA( K ) = 0.D0
   51 CONTINUE
C************************* QUANTITA' VARIE ****************************
      IF(IRSTA.GT.0) THEN
       ASPS1=ASYPS1
      ELSE
       ASPS1=0.D0
      END IF
C***********************************************************************
      AGE  = DLOG10(TEMPO)
      A(1) = DLOG10(HT1)
      A(2) = ELLOG
      A(3) = TEFF
      RAG  = STBSUN*DSQRT(10.D0**ELLOG)/((10.D0**TEFF)**2)
      RTOT = RAG*SUNRcm*1.D-10
      EMSOL= EMTOT/(SUNMg33)
      A(4) = DLOG10(RAG)
      A(6) = EMSOL
      A(7) = ETA
      A(8) = DMAS/(SUNMg33*HT1)
      WRITE(2,100) AGE,A,ITERAZ,MAXME,NMD
C *************** CERCA ZONE CONVETTIVE ********************************
C  ***************   CERCA CORE CONVETTIVO  *************************
      KCORE= 0
      ICOV = 0
      K=1
      DO WHILE (G(6,K).GE.0.D0.AND.K.LE.(MAXME-1))
         K=K+1
      END DO
      IF(K.GT.1.AND.K.LE.(MAXME)) THEN
         ICOV=1
         KCORE=K   !PRIMO MESH RADIATIVO
      ENDIF
C  ***************   CERCA BORDO ENVELOPE CONVETTIVO  ****************
      IENV = 0
      KENV = 0
      K=MAXME-1
      DO WHILE (G(6,K).GE.0.D0.AND.K.GT.(KCORE+1))
         K=K-1
      END DO
      IF(K.GT.1.AND.K.LT.(MAXME-1)) THEN
         IENV=1
         KENV=K   !MESH BORDO INFERIORE ENVELOPE CONVETTIVO (RADIAT.)
      ENDIF
C  ***************   CERCA SHELL CONVETTIVE  ************************
      IF(KCORE.NE.0) THEN
         K2=KCORE+1
      ELSE
         K2=2
      ENDIF
      IF(KENV.NE.0) THEN
         K3=KENV
      ELSE
         K3=MAXME-1
      ENDIF
      L1 = 0
      N1 = 0
      DO K=K2,K3
         IF(G(6,K).GE.0.D0.AND.G(6,K-1).LT.0.D0) THEN
            L1=L1+1
            IMINT(L1)=K-1  !BORDI INTERNI SHELL CONVETTIVE
         ELSE IF(G(6,K).LT.0.D0.AND.G(6,K-1).GE.0.D0) THEN
            N1=N1+1
            IMEST(N1)=K    !BORDI ESTERNI SHELL CONVETTIVE
         ENDIF
      END DO
C  ***************   VERIFICA BORDO ENVELOPE CONVETTIVO  ****************
      IF(KENV.EQ.0.AND.L1.NE.0)THEN
       app=(G(5,MAXME)-G(5,IMEST(N1)))/G(5,MAXME)
       
c       WRITE(*,*)'---------------'
c       WRITE(*,'(4I7)')L1,IMINT(L1),IMEST(N1),MAXME
c       WRITE(*,'(F10.5)')APP
c       WRITE(*,*)'---------------'
       IF(APP.LT.5.D-3)THEN
         IENV=1
       KENV=IMINT(L1)
       ENDIF
      ENDIF
C ************* STAMPA ANDAMENTI FISICA ALL'INTERNO ********************
      KSHHE  = 0
      KSHHY  = 0
      NTMAX  = 1
      NHMAX  = 1
      NHEMAX = 1
      L      = 1
      L3     = 1
      N      = 1
      N3     = 1
      TMAX   = G( 4 , 1 )
      AHMAX  = 0.D0
      AHEMAX = 0.D0
      AGRMAX = 0.D0
      ANUMAX = 0.D0
      GENPP  = 0.D0
      GENCAR = 0.D0
      GENCNO = 0.D0
      GENHE  = 0.D0
      GENGRA = 0.D0
      GENNU  = 0.D0
      CORCO  = 0.D0
      CORHE  = 0.D0
      CORCONV= 0.D0
      ENVCONV= 0.D0
      DO J=1,6
         G(J,MAXME+1)=G(J,MAXME)
      END DO
      DO K=1,MAXME
         KK=K+1
         DO J=1,5
            DG(J,K)=G(J,KK)-G(J,K)
         END DO
      END DO
      ELSU = G(2,MAXME)
      TEMP = DEXP(TEFF*COND)
      RTOT = STB10*DSQRT(ELSU)/(TEMP*TEMP)
      
C **************** LOOP PRINCIPALE *************************************
c         if(xxx(3,1).gt.0.d0)then
c          RAG  = STBSUN*DSQRT(10.D0**ELLOG)/((10.D0**TEFF)**2)
c          RAG1  = RAG 
c       ddf=1.d-2
c       IF(RAG.GT.(rlim-ddf).AND.RAG.LT.(rlim+ddf))then 
cc       WRITE(*,'(2e16.9)')RAG1
cc       pause
c        DO K=1,MAXNE-1
c            XX(K)=XXX(K,MAXME-1)
c         END DO
c         BOT=1.D0-FRAZ
c         CALL NATMOS(ELLOG,TEFF,EMSOL,XX,ZINI,BOT,ALPHA,
c     $       O(1),O(2001),O(4001),O(6001),NATM,1,0,0,ISINGOLA,
c     $        MAXNE,nmd)
c       ENDIF
c         endif
      IF(NABLA.EQ.3) THEN
         DO K=1,MAXNE-1
            XX(K)=XXX(K,MAXME-1)
         END DO
         BOT=1.D0-FRAZ
         CALL NATMOS(ELLOG,TEFF,EMSOL,XX,ZINI,BOT,ALPHA,
     $       O(1),O(2001),O(4001),O(6001),NATM,1,0,0,ISINGOLA,
     $        MAXNE,nmd)
         WRITE(2,108)
      ENDIF
      do kk=1,10
         flusso(kk)=0.d0
      enddo
      DO 1 JF = 1 , MAXME
      JFF = JF + 1
      RR = ( G( 1 , JF ) + G( 1 , JFF ) ) / 2.D0
      ZL = ( G( 2 , JF ) + G( 2 , JFF ) ) / 2.D0
      PP = ( G( 3 , JF ) + G( 3 , JFF ) ) / 2.D0
      TT = ( G( 4 , JF ) + G( 4 , JFF ) ) / 2.D0
      EM = ( G( 5 , JF ) + G( 5 , JFF ) ) / 2.D0
      TM = TT * 1.D+06
      PR = PP * 1.D+17
      EPSN  = 0.D0
      GRAVI = 0.D0
      VEELCO= 0.D0
C      PMAU = 0.D0
      CALL EOS(JF,1,MAXNE,1,PR,TM,XXX,'STAMPA  ',RHO,DAD,CSP,PMAU,G1,DEL
     #,EMUE,6)
      WRITE(*,*) PMAU
      if(emue.lt.0.d0)then
       write(*,*)jf,pr,tm,g1,emue
       pause
      endif
      if(jf.eq.1.and.icov.eq.0)then
        pasanni=0.d0
      DO J=1,MAXNE
         YF(J)=XXV(J,1)/ATW(J)
        END DO
        CALL SKIPFASE (TM,RHO,YF,JF,IFASE)
        IF(IFASE.EQ.1) then
         NNP = NP(IFASE)
           CALL NCROSS (RHO,TM,YF,IFASE,NNP,CROS,1)
c         write(*,*)'icov stampe=',icov
c         write(*,*)icov,cros(3,1),cros(4,1),rho
            THE3=1.d0/((yf(2)*cros(3,1)+yf(3)*cros(4,1))*rho)
c         pause
      ENDIF
      endif


       if(jf.eq.1.)then
            THE3=THE3/secanni
c         write(*,*)'telio      passetto'
c         write(*,*)the3,pasanni
       endif

      CALL EPSI( RHO , TM , JF , EPSN , EPSPP , EPSCNO , EPSA , ECAR,1)
      GRAVI=ENGRAV(JF)
      CALL NEUTR( RHO , TM , XXX(1,JF) , EPNEU )
      CALL NKAPPA(JF,1,MAXNE,1,RHO,TM,XXX,CAP)
      DRAD(JF)= CRAD * CAP * PP * ZL / ( EM * ( TT**4 ) )
      GRSAD=0.D0
      IF(DRAD(JF).GE.DAD) THEN
       GI=GRAVC0*1.D+13*EM/(RR*RR)
         Q=PMAU
         IF(Q.LT.0.D0) Q=0.D0
         CALL NSUPERA(RHO,PR,TM,CAP,Q,CSP,DAD,GI,ALPHA,GRSAD,VEELCO,
     $                DRAD(JF))
      ENDIF
      IF(EPSPP.EQ.0.D0.AND.EPSCNO.EQ.0.D0) THEN
         DELTA=0.D0
      ELSE
         DELTA=EPSPP/(EPSPP+EPSCNO)
      ENDIF
      EPSH   = EPSPP + EPSCNO
      GENPP  = GENPP  + DELTA * EPSH * 10.D0 * DG( 5 , JF )
      GENCNO = GENCNO + (1.D0 - DELTA ) * EPSH * 10.D0 * DG( 5 , JF )
      GENHE  = GENHE  + EPSA  * 10.D0 * DG( 5 , JF )
      GENCAR = GENCAR + ECAR  * 10.D0 * DG( 5 , JF )
      GENGRA = GENGRA + GRAVI * 10.D0 * DG( 5 , JF )
      GENNU  = GENNU  + EPNEU * 10.D0 * DG( 5 , JF )
      IF( JF . EQ . 1 ) ROC   = RHO
      IF( JF . EQ . 1 ) DELCE = DELTA
C *****************  CERCA MASSIMO DELLE 3A   **************************
      IF( EPSA . GT . AHEMAX  ) THEN
         AHEMAX = EPSA
         NHEMAX = JF
         ROHEMX = RHO
      ENDIF
C *****************  CERCA MASSIMO DELL' H    **************************
      IF( EPSH . GT .  AHMAX  ) THEN
       AHMAX = EPSH
         NHMAX = JF
         ROHMX = RHO
      ENDIF

C *****************  CERCA MASSIMO DEI NEUTRINI  ***********************
      IF( EPNEU . GT . ANUMAX ) THEN
         ANUMAX = EPNEU
      ENDIF
C *****************  CERCA MASSIMO GRAVITAZIONALE  *********************
      IF( DABS( GRAVI ) . LE . DABS( AGRMAX ) ) GO TO 6
      AGRMAX = GRAVI
    6 CONTINUE
C *****************  CERCA TEMPERATURA MASSIMA  ************************
      IF( TT . LE . TMAX ) GO TO 7
      TMAX    = TT
      NTMAX   = JF
      AA( 1 ) = EM / EMTOT
      AA( 2 ) = DLOG10( PP ) + 17.D0
      AA( 3 ) = DLOG10( TT ) +  6.D0
      AA( 4 ) = DLOG10( RHO )
    7 CONTINUE
C *********  CALCOLA QUANTITA' AL BORDO DEL CORE CONVETTIVO  ***********
      IF( JF.EQ.KCORE) THEN
         BB( 1 ) = EM / EMTOT
         BB( 2 ) = RR / RTOT
         BB( 3 ) = ZL / G( 2 , MAXME )
         BB( 4 ) = DLOG10( PP ) + 17.D0
         BB( 5 ) = DLOG10( TT ) +  6.D0
         BB( 6 ) = DLOG10( RHO )
        IF(KCORE.LT.MAXME)THEN
         BB( 7 ) = XXX( 1 , JF + 1 )
         BB( 8 ) = XXX( 2 , JF + 1 )
         BB( 9 ) = XXX( 3 , JF + 1 )
         BB(10 ) = XXX( 4 , JF + 1 )
         BB(11 ) = XXX( 5 , JF + 1 )
         BB(12 ) = XXX( 6 , JF + 1 )
        ELSE
         BB( 7 ) = XXX( 1 , JF  )
         BB( 8 ) = XXX( 2 , JF  )
         BB( 9 ) = XXX( 3 , JF  )
         BB(10 ) = XXX( 4 , JF  )
         BB(11 ) = XXX( 5 , JF  )
         BB(12 ) = XXX( 6 , JF  )
        ENDIF
      ENDIF
C *** CALCOLA GRANDEZZE RELATIVE A BORDO INTERNO SHELL CONVETTIVE ******
      IF(L1.GT.0.AND.JF.EQ.IMINT(L)) THEN
         CW( 1 , L ) = EM / EMTOT
         CW( 2 , L ) = RR / RTOT
         CW( 3 , L ) = ZL / G( 2 , MAXME )
         CW( 4 , L ) = DLOG10( PP ) + 17.D0
         CW( 5 , L ) = DLOG10( TT ) +  6.D0
         CW( 6 , L ) = DLOG10( RHO )
         CW( 7 , L ) = XXX( 1 , JF+1)   !CHIMICA ZONA CONVETTIVA
         CW( 8 , L ) = XXX( 3 , JF+1)
         CW( 9 , L ) = XXX( 4 , JF+1)
         CW(10 , L ) = XXX( 5 , JF+1)
         CW(11 , L ) = XXX( 6 , JF+1)
         L = L + 1
      ENDIF
C *** CALCOLA GRANDEZZE RELATIVE A BORDO ESTERNO SHELL CONVETTIVE ******
      IF(N1.GT.0.AND.JF.EQ.IMEST(N)) THEN
         DD( 1 , N ) = EM / EMTOT
         DD( 2 , N ) = RR / RTOT
         DD( 3 , N ) = ZL / G( 2 , MAXME )
         DD( 4 , N ) = DLOG10( PP ) + 17.D0
         DD( 5 , N ) = DLOG10( TT ) +  6.D0
         DD( 6 , N ) = DLOG10( RHO )
         N = N + 1
      ENDIF
C *** CALCOLA GRANDEZZE RELATIVE ALL'ENVELOPE CONVETTIVO  **************
      IF( JF.EQ.KENV ) THEN
         BE( 1 ) = EM / EMTOT
         BE( 2 ) = RR / RTOT
         BE( 3 ) = ZL / G( 2 , MAXME )
         BE( 4 ) = DLOG10( PP ) + 17.D0
         BE( 5 ) = DLOG10( TT ) +  6.D0
         BE( 6 ) = DLOG10( RHO )
         BE( 7 ) = XXX( 1 , MAXME-1 )
         BE( 8 ) = XXX( 2 , MAXME-1 )
         BE( 9 ) = XXX( 3 , MAXME-1 )
         BE(10 ) = XXX( 4 , MAXME-1 )
         BE(11 ) = XXX( 5 , MAXME-1 )
         BE(12 ) = XXX( 6 , MAXME-1 )
      ENDIF
C ******** STAMPA VARIE QUANTITA' SE NABLA=3 OGNI KSD MESHPOINTS *******
      IF(MEMHEB.GT.0.AND.NABLA.NE.3) GO TO 58
      IF( NABLA . NE . 3 ) GO TO 1
      IF( ICOV.EQ.1.AND.JF.EQ.1 )             WRITE(2,101)
      IF( ICOV.EQ.1.AND.JF.EQ.KCORE )         WRITE(2,102)
      IF( IENV.EQ.1.AND.JF.EQ.(KENV+1) )      WRITE(2,107)
      IF(IMINT(L3).GT.0) THEN
         IF( JF.EQ.(IMINT(L3)+1) ) THEN
            WRITE(2,103)  !BORDO INTERNO SHELL CONVETTIVA
            L3=L3+1
         END IF
      ENDIF
      IF(IMEST(N3).GT.0) THEN
         IF( JF.EQ.(IMEST(N3)) ) THEN
            WRITE(2,104)  !BORDO ESTERNO SHELL CONVETTIVA
            N3=N3+1
         END IF
      END IF
      IF( JF.GT.1 ) THEN  !CORE DI ELIO O DI CO
         IF( XXX(3,JF-1).LE.0.D0.AND.XXX(3,JF).GT.0.D0 ) WRITE(2,105)
         IF( XXX(1,JF-1).LE.0.D0.AND.XXX(1,JF).GT.0.D0 ) WRITE(2,106)
      END IF
      JRONA = 0
      JRANA = JF - ( JF / KSD ) * KSD
      IF( JF.EQ.1.OR.JF.EQ.MAXME ) JRONA = 1
      IF( JRANA.NE.0.AND.JRONA.EQ.0 ) GO TO 1
   58 CONTINUE
      BETA=(CN3*TM**4)/PR 
      VSOUND=DSQRT(PP*1.D+17/RHO)
      EE( 1 )  = EM / EMTOT
      EE( 2 )  = RR / RTOT
      EE( 3 )  = DLOG10( PP ) + 17.D0
      EE( 4 )  = DLOG10( TT ) +  6.D0
      EE( 5 )  = DLOG10( RHO )
      EE( 6 )  = ZL / G( 2 , MAXME )
      EE( 7 )  = EPSN
      EE( 8 )  = GRAVI
      EE( 9 )  = EPNEU
      EE( 10 ) = GRSAD
      EE( 11 ) = DRAD(JF)
      EE( 12 ) = DAD
C     EE( 13 ) = CAP
      EE( 13 ) = VEELCO/VSOUND
      EE( 14 ) = PMAU
C********* STAMPA VALORI PER S-PROCESS CODE SE MEMHEB > 0 **************
      IF(XXX(1,1).GT.0.0D0) GO TO 59
      IF(DRAD(JF).LT.DAD.OR.EE(4).LT.7.8D0.OR.MEMHEB.LE.0) GO TO 59
      IF(XXX(3,1).GT.0.7D0.OR.XXX(3,1).LE.0.D0) GO TO 59
   59 CONTINUE
  124 FORMAT(I5)
C***********************************************************************
      IF(NABLA.EQ.3)WRITE(2,109) EE , JF
    1 CONTINUE         !FINE LOOP PRINCIPALE SU JF
C*********8************ End of main loop over JF
      FLPP  = GENPP  / G( 2 , MAXME )
      FLCNO = GENCNO / G( 2 , MAXME )
      FL3A  = GENHE  / G( 2 , MAXME )
      FLGRA = GENGRA / G( 2 , MAXME )
      FLNEU = GENNU  / G( 2 , MAXME )
      FLCAR = GENCAR / G( 2 , MAXME )
C **************** CALCOLA QUANTITA' AL CENTRO   ***********************
      CV( 1 ) = DLOG10( G( 3 , 1 ) ) + 17.D0
      CV( 2 ) = DLOG10( G( 4 , 1 ) ) +  6.D0
      CV( 3 ) = DLOG10( ROC )
      CV( 4 ) = XXX( 1 , 1 )
      CV( 5 ) = XXX( 2 , 1 )
      CV( 6 ) = XXX( 3 , 1 )
      CV( 7 ) = XXX( 4 , 1 )
      CV( 8 ) = XXX( 5 , 1 )
      CV( 9 ) = XXX( 6 , 1 )
C ********* CALCOLA GRANDEZZE RELATIVE ALLA SHELL DI IDROGENO **********
      IF( XXX( 1 , 1 ) . GT . 0.D0 ) GO TO 22
      XHA = 0.01D0 * XXX( 1 , MAXME )
      XHB = 0.90D0 * XXX( 1 , MAXME )
      DO 23 K = 2 , MAXME
      IF( XXX( 1 , K ) . GT . 0.D0 ) GO TO 24
   23 CONTINUE
      IF( K . GE . MAXME ) GO TO 88
   24 CONTINUE
      CORHE = G( 5 , K ) / (SUNMg33)
      DO 25 LMA = 2 , MAXME
      IF( XXX( 1 , LMA ) . GE . XHA ) GO TO 26
   25 CONTINUE
      IF( LMA . GT . MAXME ) LMA = MAXME
   26 CONTINUE
      DO 27 LMB = LMA , MAXME
      IF( XXX( 1 , LMB ) . GE . XHB ) GO TO 28
   27 CONTINUE
      IF( LMB . GT . MAXME ) LMB = MAXME
   28 CONTINUE
c      FF(  1 ) = CORHE
      FF(  1 ) = G( 5 , NHMAX ) / EMTOT * EMSOL
      FF(  2 ) = AHMAX
      FF(  3 ) = G( 1 , NHMAX ) / RTOT
      FF(  4 ) = G( 2 , NHMAX ) / G( 2 , MAXME )
      FF(  5 ) = DLOG10( G( 3 , NHMAX ) ) + 17.D0
      FF(  6 ) = DLOG10( G( 4 , NHMAX ) ) +  6.D0
      FF(  7 ) = DLOG10( ROHMX )
      FF(  8 ) = ( G( 5 , LMB ) - G( 5 , LMA ) ) / EMTOT * EMSOL
      KSHHY = 1
   88 CONTINUE
C ********* CALCOLA GRANDEZZE RELATIVE ALLA SHELL DI ELIO **************
      IF( XXX( 3 , 1 ) . GT . 0.D0 ) GO TO 22
      XEA = 0.01D0
      XEB = 0.90D0
      DO 30 K = 2 , MAXME
      IF( XXX( 3 , K ) . GT . 0.D0 ) GO TO 31
   30 CONTINUE
      IF( K . GT . MAXME ) K = MAXME
   31 CONTINUE
      CORCO = G( 5 , K ) / (SUNMg33)
      DO 32 LEA = 2 , MAXME
      IF( XXX( 3 , LEA ) . GE . XEA ) GO TO 33
   32 CONTINUE
      IF( LEA . GT . MAXME ) LEA = MAXME
   33 CONTINUE
      DO 34 LEB = LEA , MAXME
      IF( XXX( 3 , LEB ) . GE . XEB ) GO TO 35
   34 CONTINUE
      IF( LEB . GT . MAXME ) LEB = MAXME
   35 CONTINUE
      HH( 1 ) = G( 5 , NHEMAX ) / EMTOT * EMSOL
      HH( 2 ) = AHEMAX
      HH( 3 ) = G( 1 , NHEMAX ) / RTOT
      HH( 4 ) = G( 2 , NHEMAX ) / G( 2 , MAXME )
      HH( 5 ) = DLOG10( G( 3 , NHEMAX ) ) + 17.D0
      HH( 6 ) = DLOG10( G( 4 , NHEMAX ) ) +  6.D0
      HH( 7 ) = DLOG10( ROHEMX )
      HH( 8 ) = ( G( 5 , LEB ) - G( 5 , LEA ) ) / EMTOT * EMSOL
      KSHHE = 1
   22 CONTINUE
C **************** CHIMICA E VARIAZIONI TEMPORALI FISICA ***************
      IF( NABLA . EQ . 3 ) WRITE(2,111)
      DO 36 K = 2 , MAXME
      DO 37 J = 1 , MAXMV
      IF( GG( 5 , J ) . GE . G( 5 , K ) ) GO TO 38
   37 CONTINUE
   38 IF( J . GT . MAXMV ) J = MAXMV
      RAT = (G( 5 , K ) - GG( 5 , J-1 )) / (GG( 5 , J ) - GG( 5 , J-1 ))
      DO 39 I = 1 , 4
      FISV = GG( I , J-1 ) + RAT * ( GG( I , J ) - GG( I , J-1 ) )
      IF( DABS( FISV ) . LT . 1.D-14 ) GO TO 39
      DFIS( I ) = ( G( I , K ) - FISV ) / FISV * 1.D+02
   39 CONTINUE
      DO 40 I = 1 , 4
      IF( K . GT . 2 ) GO TO 41
      IVME( I ) = K
      DFIX( I ) = DFIS( I )
   41 CONTINUE
      IF( DABS( DFIS( I ) ) . LE . DABS( DFIX( I ) ) ) GO TO 40
      DFIX( I ) = DFIS( I )
      IVME( I ) = K
   40 CONTINUE
      IF( NABLA . NE . 3 ) GO TO 36
      JRONA = 0
      JRANA = K - ( K / KSD ) * KSD
      IF( K . EQ . 2 . OR . K . EQ . MAXME ) JRONA = 1
      IF( JRANA . NE . 0 . AND . JRONA . EQ . 0 ) GO TO 36
      KM = K - 1
      x1mas=(G(5,k)+G(5,K-1))/2.d0
      x1mas=x1mas/emtot
      XI = XXX( 1 , KM )
      XHE3 = XXX( 2 , KM )
      XHE4 = XXX( 3 , KM )
      XC = XXX( 4 , KM )
      XN = XXX( 5 , KM )
      XO = XXX( 6 , KM )
      XD = XXX( 7 , KM )
      XLi7 = XXX( 8 , KM )
      XBe7 = XXX( 9 , KM )
      XC13 = XXX( 10 , KM )
      XN15 = XXX( 11 , KM )
      XO17 = XXX( 12 , KM )
      XO18 = XXX( 13 , KM )
      XFe56=XXX(26,KM)
      IF( XI . LE . 0. . AND . XHE4 . LE . 0. ) THEN
      XI   = CARXY( 1 , KM )
      XHE4 = CARXY( 2 , KM )
C      XPMAU = PMAU(KM)
      ENDIF
      WRITE(2,112) x1mas,XI,XD,XHE3,XHE4,XLi7,XBe7,XC,XC13,XN,
     # XN15,XO,XO17,XO18,XFe56,KM
   36 CONTINUE
C ********************* STAMPA FISICA VARIA ****************************
      IF( NABLA.EQ.3 ) CALL PLOTTA( RTOT )
      WRITE(2,113) CV
      IF( KCORE.GT.0 ) WRITE(2, 114 ) BB
      IF( NTMAX . GT . 1 ) WRITE(2, 115 ) AA
      IF( KSHHE . EQ . 1 ) WRITE(2, 116 ) HH
      IF( N1.GT.0 ) THEN
         DO  LL = 1 , N1
            WRITE(2, 117 ) ( CW( I , LL ) , I = 1 , 11 )
            WRITE(2, 118 ) ( DD( I , LL ) , I = 1 ,  6 )
         END DO
      ENDIF
      IF( KSHHY.EQ.1 ) WRITE(2, 119 ) FF
      IF( KENV.GT.0 ) THEN
         WRITE(2,120)(BE(I),I=1,12)
      ELSE
         WRITE(2,126)(XXX(J,MAXME-1),J=1,13)
      ENDIF
      WRITE(2,121)GENPP,FLPP,GENCNO,FLCNO,GENHE,FL3A,GENGRA,FLGRA,GENNU,
     &FLNEU,GENCAR,FLCAR
      WRITE(2,122) ( IVME(LL) , DFIX(LL) , LL = 1 , 4 )
C***********************************************************************
C stampe schermo
C***********************************************************************
      WRITE(*,199)
      WRITE(*,200)AGE,A(1),A(2),A(3),A(4),CV(1),CV(2),CV(3),EMSOL,NMD
      WRITE(*,200)AGE,(A(L),L=1,4),CV(1),CV(2),CV(3),EMSOL,NMD
C------------- formato per calibrazione solare
c     WRITE(*,210)AGE,A(1),A(2),A(3),A(4),CV(1),CV(2),CV(3),EMSOL,NMD
C-------------------------------------------------------

      IF(XXX(1,1).GT.0.D0) THEN
         XCEN = XXX(1,1)
         WRITE(*,201)XCEN,FLPP,FLCNO,FL3A,FLGRA,FLNEU,MAXME
      ELSEIF(XXX(3,1).GT.0.D0) THEN
         XCEN = XXX(3,1)
         WRITE(*,202)XCEN,FLPP,FLCNO,FL3A,FLGRA,FLNEU,MAXME
      ELSEIF(XXX(4,1).GT.0.D0) THEN
         XCEN = XXX(4,1)
         WRITE(*,203)XCEN,FLCAR,FLCNO,FL3A,FLGRA,FLNEU,MAXME
      ENDIF
      CORCONV=BB(1)*EMSOL
      AMMAX=AA(1)*EMSOL
      ZX=(XXX(26,MAXME-1)/ABSOL(26))/XXX(1,MAXME-1)

      WRITE(*,204) CORCONV,CORHE,CORCO,BE(1),AMMAX,AA(3)
      WRITE(*,125) (XXX(J,MAXME-1),J=1,6),XXX(26,MAXME-1),ZX
      WRITE(*,130) ASPS1
      WRITE(*,122) ( IVME(LL) , DFIX(LL) , LL = 1 , 4 )
C ************ MEMORIZZA VARIE QUANTITA' PER GRAFICI *******************
c ******* STORE MISCELLANEOUS QUANTITIES FOR GRAPHICS
      JRANA = NMD - ( NMD / IFVA ) * IFVA
      IF( JRANA . NE . 0 ) GO TO 52
      IF( L1 . LT . 1 ) L1 = 1
      RPQ( 1 ) = XXX( 1 , 1 )
      RPQ( 2 ) = XXX( 3 , 1 )
      RPQ( 3 ) = XXX( 4 , 1 )
      RPQ( 4 ) = XXX( 5 , 1 )
      RPQ( 5 ) = XXX( 6 , 1 )
      RPQ( 6 ) = FLPP
      RPQ( 7 ) = FLCNO
      RPQ( 8 ) = FL3A
      RPQ( 9 ) = FLGRA
      SUMMA = 0.D0
      DO JJ=1,20
         ES(JJ) = XXX(JJ,1)
         SUMMA = SUMMA + XXX(JJ,1)
      END DO
      IF( ES(1).LE.0.D0.AND.ES(3).LE.0.D0) THEN
      ES(1) = CARXY( 1 , 1 )
      ES(3) = CARXY( 2 , 1 )
      SUMMA = SUMMA + ES(1) + ES(2)
      ENDIF
C####################################################################
      WRITE(16,'(I4,D22.15,10E13.5,100E13.5)')NMD,AGE,
     # (A(J),J=2,5),BB(1)*EMSOL,CV(3),CV(2),BE(1)*EMSOL,BE(6),BE(5),
     # (XXX(ICA(1,J),MAXME-1),J=1,NE(1))
C#################################################################### 
C####################################################################
      WRITE(66,'(I4,D22.15,10E13.5,100E13.5)')NMD,AGE,
     # (A(J),J=2,5),BB(1)*EMSOL,CV(3),CV(2),BE(1)*EMSOL,BE(6),BE(5),
     # (XXX(ICA(1,J),1),J=1,NE(1))
      WRITE(56,'(I4,4D16.7,10E13.5,100E13.5)')NMD,AGE,ht1,the3,
     #pasanni,(A(J),J=2,5),BB(1)*EMSOL,CV(3),CV(2),BE(1)*EMSOL,BE(6),
     #BE(5),(XXX(ICA(1,J),1),J=1,NE(1))
C####################################################################
      WRITE(46,'(f15.10,9e15.7)')age,a(2),(FLUSSO(l),L=1,8)
C####################################################################    
      IF(ISTA.EQ.1)THEN
       IF(ISTART.EQ.1)THEN
         IF(IIST.EQ.0)THEN
         EVFASE=' ZAHB '
       ELSE
         EVFASE='MODANT'
       ENDIF
       ELSEIF(ISTART.EQ.2)THEN
         EVFASE='PRE-MS'
       ELSEIF(ISTART.EQ.3)THEN
         EVFASE=' ZAMS '
       ELSEIF(ISTART.EQ.4)THEN
         EVFASE=' ZAHB '
       ELSEIF(ISTART.EQ.5)THEN
         EVFASE='for HB'
       ENDIF
       OVE=0.d0
       DIF=' no'
       if(INDIFFINI.GT.0)DIF='yes'
       write(10,773)
       if(istart.ge.4.or.istart.eq.1)then
        write(10,772)EMSOL,YINI,ZINI,ALPHA,ETA,FACTOR,DIF,EVFASE,PROGM
       else
        PROGM=EMSOL
        write(10,772)EMSOL,YINI,ZINI,ALPHA,ETA,FACTOR,DIF,EVFASE,PROGM
       endif
       ista=0
       write(10,773)
      ENDIF

 771  FORMAT ('Mass=',F8.4,1X,'|',1x,'Y=',F6.3,1X,'|',1x,'Z=',1P,E9.2,
     #1X,'|',1x,'ML=',0P,F5.2,1X,'|',1x,'Mass loss=',F4.2,1X,'|',1x,
     #'Oversh= ',f5.3,1X,'|',1x,'Diff: ',A3,1X,'|',1x,'from: ',a6)
 772  FORMAT ('Mass=',F8.4,1X,'|',1x,'Y=',F6.3,1X,'|',1x,'Z=',1P,E9.2,
     #1X,'|',1x,'ML=',0P,F5.2,1X,'|',1x,'Mass loss=',F4.2,1X,'|',1x,
     #'Oversh= ',f5.3,1X,'|',1x,'Diff: ',A3,1X,'|',1x,'from: ',a6,
     #1X,'|',1x,'Mass of prog.=',F8.4)
 773  FORMAT ('=========================================================
     #==================================================================
     #============')

      WRITE( 10,300 ) NMD+10000*NNMD,AGE,A(1),A(2),A(3),A(4),A(6),CV(1)
     &,CV(2),CV(3),RPQ,ASPS1
      WRITE( 10,301 ) AA,FF,HH,BB(1),BB(2),BB(3),BB(4),BB(5),BB(7),
     #BB(9),BB(10),BB(11),BB(12),BE(1),BE(2),BE(3),BE(4),BE(5),BE(7),
     #BE(9),BE(10),BE(11),BE(12)
c      WRITE( 10,301 ) AA,FF,HH,(BB(I),I=1,5),BB(7),(BB(K),K=9,12),
c     &                (BE(I),I=1,5),BE(7),(BE(I),I=9,12)
      IF (ICHI.EQ.1 ) THEN
         WRITE( 14,300 ) NMD,AGE,SUMMA,(ES(K),K=1,16)
         WRITE( 14,301 ) (ES(K),K=17,20)
      ENDIF
C     IF (INEU.EQ.1) THEN
C        WRITE(15,303 ) AGE,(PPNEU(K),K=1,6),TCRIT
C        WRITE( *,304 ) AGE,(PPNEU(K),K=1,6),TCRIT
C     ENDIF
   52 CONTINUE
      IF(IPICCI.EQ.1) RETURN
      IF( ( NMD - ( NMD / IFMO ) * IFMO ) . NE . 0 ) RETURN
      WRITE(11,302) NMD,MAXME,AGE,A(2),A(3),XINI,ZINI,EMTOT
      DO 48 K = 1 , MAXME
      DO 48 I = 1 , 11
      IF( I . GT . 5 ) GO TO 49
      STAR( I , K ) =   G(  I , K )
      GO TO 48
   49 CONTINUE
      II = I - 5
      STAR( I , K ) = XXX( II , K )
   48 CONTINUE
      WRITE( 11 , 301 ) (( STAR(I,K) , I = 1 , 11 ) , K = 1 , MAXME )
      RETURN
  100 FORMAT(/,' AGE:',F15.11,' DT:',F8.4,' L:',F8.4,' TE:',F8.4,' R:',F
     &8.4,' GR:',1P,D10.2,' M:',0P,F9.5,' E:',F6.3,' DM:',1P,D10.3,I2
     &,2I5)
  101 FORMAT(10X,'********** CORE CONVETTIVO **********')
  102 FORMAT(10X,'********** BORDO CONVETTIVO CENTRALE **********')
  103 FORMAT(10X,'********** BORDO INTERNO SHELL CONVETTIVA **********')
  104 FORMAT(10X,'********** BORDO ESTERNO SHELL CONVETTIVA **********')
  105 FORMAT(10X,'********** CORE DI CARBONIO - OSSIGENO **********')
  106 FORMAT(10X,'********** CORE DI ELIO **********')
  107 FORMAT(10X,'********** BORDO CONVEZIONE SUPERFICIALE **********')
  108 FORMAT(///,'   M/MTOT      R/RTOT  LOG(P)  LOG(T)  LOG(RO)   L/LSU
     &P     ENUC      EGRAV     ENEU      DELTA     RADIAT    ADIAB  V-S
     &OUND  PMAU  MESH',/)
  109 FORMAT(F12.9,F10.7,F8.4,F7.4,F8.4,
     &1P,D13.5,1P,3D10.3,1P,D9.2,1P,D11.3,
     &0P,F6.3,1P,D9.2,1x,F12.9,I5)
C 111 FORMAT(1H1,/,'   M/MTOT     IDROGENO    ELIO3       ELIO4       CA
C    &RBONIO  AZOTO     OSSIGENO      DR/DT     DL/DT   DP/DT   DTEMP/DT
C    &MESH ')
  111 FORMAT(1H1,/,'      Mr         H        D        HE3     HE4            
     #Li7      Be7     C12      C13       N14      N15      O16      O17        
     &     O18      Fe56    MESH')
C 112 FORMAT(F12.9,1P,6D12.5,4D10.2,I5)
  112 FORMAT(1X,f12.9,14D9.2,I5)
  113 FORMAT(1X,'VALORI CENTRALI:',3F8.4,1P,6D12.5)
  114 FORMAT(1X,'CORE CONVETTIVO:',F10.7,F9.6,1P,D10.2,0P,3F7.3,1P,6D10.
     #3)
  115 FORMAT(1X,'TEMPERATURA MASSIMA NON CENTRALE:',F12.9,3F8.4)
  116 FORMAT(1X,'HE-SHELL:',F10.7,1P,3D10.2,0P,3F7.3,1P,D10.3)
  117 FORMAT(1X,'BORDO INT. SH. CONVET.:',
     &2F10.7,1P,D10.2,0P,3F7.3,1P,5D10.3)
  118 FORMAT(1X,'BORDO ESTERNO SHELL CONVETTIVA:',
     &2F10.7,1P,D10.2,0P,3F7.3)
  119 FORMAT(1X,' H-SHELL:',F10.7,1P,3D10.2,0P,3F7.3,1P,D10.3)
  120 FORMAT(1X,'BORDO ENV. CONV.:',
     &2F10.7,1P,D10.2,0P,3F7.3,1P,6D10.3)
  121 FORMAT(1X,1P,12D10.2)
  122 FORMAT(1X,' DR(',I4,') =',F8.4,
     &          ' DL(',I4,') =',F8.4,
     &          ' DP(',I4,') =',F8.4,
     &          ' DT(',I4,') =',F8.4)
  126 FORMAT(' CHI-SUP:',1P,13E9.2)
  125 FORMAT(' XSUP ','H:',1P,E9.2,' HE3:',1P,E9.2,' HE4:',1P,E9.2,
     &' C:',1P,E9.2,' N:',1P,E9.2,' O:',1P,E9.2,' Fe:',1P,E9.2,
     #' Z/X:',0P,F6.4)
  128 FORMAT(3f11.7,1P,2E10.2,E11.4,3E10.2,E11.4)
  129 FORMAT(' log(L/Lo)    logTe  ',3x,' log(Tc)',7x,'H',7x,'HE3',7x,'H
     #E4',8x,' C',9x,'N',9x,'O',8x,'Z')
  130 FORMAT(' Asymptotic Period Spacing l=1 (sec) ---> ',F9.2)   
  199 FORMAT('************************************************************
     &***********************************************')

  200 FORMAT(' t:',F8.5,' Dt:',F7.4,' L:',F7.4,' Te:',F7.4,' R:',F7.4,
     &' Pc:',F7.3,' Tc:',F6.3,' ROc:',F6.3,' M:',F9.5,' #:',I5)
  210 FORMAT(' t:',F8.5,' Dt:',F7.4,' L:',F11.8,' Te:',F11.8,' R:'
     $,F11.8,' Pc:',F7.3,' Tc:',F6.3,' ROc:',F6.3,' M:',F9.5,' #:',I5)
     
  201 FORMAT(' H :',1P,D9.2,' - Lpp:',D9.2,' - Lcno:',D9.2,' - L3a:',
     &D9.2,' - Lgr:',D9.2,' - Lneu:',D9.2,0P,' mesh:',I6)
  202 FORMAT(' HE:',1P,D9.2,' - Lpp:',D9.2,' - Lcno:',D9.2,' - L3a:',
     &D9.2,' - Lgr:',D9.2,' - Lneu:',D9.2,0P,' mesh:',I6)
  203 FORMAT(' C:',1P,D9.2,' - Lpp:',D9.2,' - Lcno:',D9.2,' - L3a:',
     &D9.2,' - Lgr:',D9.2,' - Lneu:',D9.2,0P,' mesh:',I6)
  204 FORMAT(' Mcc/Mo:',F8.4,' - Mhe/Mo:',F8.4,' - Mco/Mo:',F8.4,
     &' - Mr_envc:',F8.4,' - MTmax/Mo:',F8.4,' - Tmax:',F6.4)
C*********************** IBM *******************************************
C 300 FORMAT(I4,A8,17A4)
C 301 FORMAT(20A4)
C 302 FORMAT(2I5,A8,5A4)
C***********************************************************************
C
C************************ FPS ******************************************
  300 FORMAT(I5,1P,D22.15,3E16.8,/,5E16.8,/,5E16.8,/,5E16.8)
  301 FORMAT(1P,5E16.8)
  302 FORMAT(2I5,1P,D22.15,3E16.8,/,2E16.8)
  303 FORMAT(F12.8,1P,7D10.3)
  304 FORMAT(F12.8,1P,7D9.2)
C***********************************************************************
      END ! end of subroutine

      SUBROUTINE STOPERR(PRESS,TE,ELEMENTO,MESH,ERROR,ROUT)
      IMPLICIT NONE
      CHARACTER*1 ROUT(8)
      CHARACTER*70 etadifover
      INTEGER *4 ERROR,ELEMENTO,MESH
      REAL    *8 PRESS,TE
  100 FORMAT(' USCITO FUORI TABELLE FISICHE',/,1X,
     #'P=',1PD12.4,1X,'T=',1PD12.4,1X,'MESH=',I4,/,1X,'STOP IN: ',8A1)
  101 FORMAT(' ATTENZIONE! MORTO NELLA STATE. CONTROLLARE LA CHIMICA',
     #/,1X,'ELEMENTO=',I3,1X,'MESH=',I4,/,1X,'STOP IN: ',8A1)
  102 FORMAT(' ATTENZIONE! PRESSIONE DEL GAS NEGATIVA',
     #/,1X,'P=',1PD12.4,1X,'T=',1PD12.4,1X,'MESH=',I4,/,1X,'STOP IN: '
     #,8A1)
      IF (ERROR.EQ.1) THEN
         WRITE(*,101)ELEMENTO,MESH,ROUT
         WRITE(2,101)ELEMENTO,MESH,ROUT
         STOP
      ELSE IF (ERROR.EQ.2) THEN
              WRITE(*,100)PRESS,TE,MESH,ROUT
              WRITE(2,100)PRESS,TE,MESH,ROUT
              STOP
      ELSE IF (ERROR.EQ.3) THEN
              WRITE(*,102)PRESS,TE,MESH,ROUT
              WRITE(2,102)PRESS,TE,MESH,ROUT
              STOP
      END IF
      END

      SUBROUTINE VEIOVE(ELL,PCEN,TCEN,EM,R,EL,P,T,PV,TV,
     &                  INTEST,IPP,KI,MFI)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      INCLUDE 'consts.2p1'
      DIMENSION PV(LIM),TV(LIM)
      DATA EPS/0.D0/    
c      IF(INTEST.EQ.1) THEN
      IF(INTEST.le.1) THEN
         EM=0.D0
         EL=0.D0
         R=0.D0
         P=PCEN
         T=TCEN
         CALL EOS(1,1,MAXNE,1,P,T,XXX,'VEIOVE  ',ROCEN,O1,CPP,
     #                 O2,G1,DEL,EMUE,3)
         NUM=MFI-1
      ELSE
         EL=ELL
         NUM=MAXME-MFI
      ENDIF
      DO K=1,NUM
         K2=K
         K1=K2+1
         IF(INTEST.EQ.2) THEN
            K2=MAXME-K+1
            K1=K2-1
         ENDIF
         DM=G(5,K1)-G(5,K2)
         EM=EM+DM
         CALL EOS(K2,1,MAXNE,1,P,T,XXX,'VEIOVE  ',RO,DAD,
     #                  CPP,O1,G1,DEL,EMUE,3)
c     WRITE(*,*)'DEBUG ', RO
c     PAUSE
         IF(K.GT.1.OR.INTEST.GT.1) THEN
            DR=DM/(4.D0*PI*RO*(R**2))
         ELSE
            DR=(3.D0*DM/(4.D0*PI*RO))**(1.D0/3.D0)
         ENDIF
         R=R+DR
         IF(R.LE.0.D0) then
            write(*,*)'raggio minore di zero al mesh:',k1,maxme
            RETURN
         end if
         CALL NKAPPA(K2,1,MAXNE,1,RO,T,XXX,CAPPA)
         GRAVI=0.D0
         PPP=P*1.D-17
         TTT=T*1.D-06
         IF(IPP.EQ.1)CALL EPSIG(PPP,TTT,CPP,GRAVI,DAD,PV(K1),TV(K1),K1)
         CALL EPSI(RO,T,K1,EPS,V1,V2,V3,V4,1)
         CALL NEUTR(RO,T,XXX(1,K2),QT)
         EPS=(EPS+QT+GRAVI)
         DL=DM*EPS
         EL=EL+DL
         IF(EL.LE.0.D0.AND.INTEST.EQ.2) then
            write(*,*)'luminosita minore di zero al mesh:',k1,maxme
            write(*,*)p,g(2,k1),dl
            write(*,*) dm, eps,p,t,ro
            write(*,*)g(5,k1),g(5,k2),k1,k2
            write(*,*)intest
            RETURN
         end if
       DP=-GRAVC0*EM*RO*DR/(R**2)
         P=P+DP
         IF(P.LE.0.D0) then
            write(*,*)'pression minore di zero al mesh:',k1,maxme
            RETURN
         end if
       DRAD=CRAD*1.D8*CAPPA*(P/(T**4))*(EL/EM)
         IF(ISTART.EQ.2.AND.IPP.EQ.0) THEN
             GRAD=DAD
         ELSE
            IF(DRAD.GE.DAD) THEN
               GRAD=DAD
            ELSE
               GRAD=DRAD
            ENDIF
         ENDIF
         DT=GRAD*T*DP/P
         T=T+DT
         IF(T.LE.0.D0) then
            write(*,*)'temperature minore di zero al mesh:',k1,maxme
            RETURN
         end if
         IF(KI.EQ.1) THEN
            G(1,K1)=R
            G(2,K1)=EL
            G(3,K1)=P
            G(4,K1)=T
            G(6,K1)=DRAD-DAD
            GG(6,K1)=DRAD-DAD
c            WRITE(2,222)K1,(G(JL,K1),JL=1,6),(XX(KL),KL=1,6),
c     &                  EPS,QT,GRAVI
c 222        FORMAT(1X,I4,1P,12D10.3,/,1X,3D12.5)
         ENDIF
      END DO
      RETURN
      END
C========================================================================            
C correz apportate a Aarhus 
C 1) la 15 quantita' che tabulavamo era la gravita'
C superficiale! loro vogliono la costante di gravitazione!
C
C========================================================================              
      SUBROUTINE RSTAM(FILENAME,IRSTA,INDIFFINI)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*25 FILENAME
      CHARACTER*17 STDIF,STML,STTAU
      CHARACTER*19 STOVER
      CHARACTER*12 STML1
      CHARACTER*5 EML,APP
            
      INCLUDE 'maincom.2p1'
      INCLUDE 'nuclear.2p1'
      INCLUDE 'consts.2p1'
      PARAMETER (MAX=15000)
      DIMENSION A(8),AA(4),CV(9)
      DIMENSION DRAD(LIM),XX(MAXNE),xmassa(lim)
      DIMENSION CROS(MAXPR,2),YF(MAXNE),BPP(26,3000),AST(26,LIM)
      DIMENSION FP(3),FRO(3),FRDO(3),AX(3),ADX(3),FM(3)
      DIMENSION DP(LIM),DPP(LIM),A5(LIM),A2(LIM),A3(LIM)
      DIMENSION GLOC(LIM),ENNEBV(LIM)
      
      COMMON/FREMEM/DG(5,LIM),O((MAXNE+1-5-1)*LIM),ENGRAV(LIM)
      common/astros/NPA,BPP


      ELSU = G(2,MAXME)
      TEMP = DEXP(TEFF*COND)
          RAG  = STBSUN*DSQRT(10.D0**ELLOG)/((10.D0**TEFF)**2)
          RTOT = RAG*SUNRcm
c      write(*,'(4e16.9)')RTOT,RTOT,RTOT,RTOT
C **************** LOOP PRINCIPALE *************************************
c       WRITE(39,773)
c       WRITE(39,774)
c       WRITE(40,775)
c       WRITE(40,776)
      DO 1 JF = 1 , MAXME-1
      JFF = JF + 1
      RR = ( G( 1 , JF ) ) ! + G( 1 , JFF ) ) / 2.D0
      ZL = ( G( 2 , JF ) ) ! + G( 2 , JFF ) ) / 2.D0
      PP = ( G( 3 , JF ) ) ! + G( 3 , JFF ) ) / 2.D0
      TT = ( G( 4 , JF ) ) ! + G( 4 , JFF ) ) / 2.D0
      EM = ( G( 5 , JF ) ) ! + G( 5 , JFF ) ) / 2.D0
      TM = TT * 1.D+06
      PR = PP * 1.D+17
      EPSN  = 0.D0
      GRAVI = 0.D0
      VEELCO= 0.D0
      CALL EOS(JF,1,MAXNE,1,PR,TM,XXX,'STAMPA  ',RHO,DAD,CSP,PMAU,G1,DEL
     #,EMUE,6)

      CALL EPSI( RHO , TM , JF , EPSN , EPSPP , EPSCNO , EPSA , ECAR,1)
      GRAVI=ENGRAV(JF)
      CALL NEUTR( RHO , TM , XXX(1,JF) , EPNEU )
      CALL NKAPPA(JF,1,MAXNE,1,RHO,TM,XXX,CAP)
      DRAD(JF)= CRAD * CAP * PP * ZL / ( EM * ( TT**4 ) )
      GRSAD=0.D0
      IF(DRAD(JF).GE.DAD) THEN
       GI=GRAVC0*1.D+13*EM/(RR*RR)
         Q=PMAU
         IF(Q.LT.0.D0) Q=0.D0
         CALL NSUPERA(RHO,PR,TM,CAP,Q,CSP,DAD,GI,ALPHA,GRSAD,VEELCO,
     $                DRAD(JF))
      ENDIF
      IF(EPSPP.EQ.0.D0.AND.EPSCNO.EQ.0.D0) THEN
         DELTA=0.D0
      ELSE
         DELTA=EPSPP/(EPSPP+EPSCNO)
      ENDIF
      EPSH   = EPSPP + EPSCNO
      GENPP  = GENPP  + DELTA * EPSH * 10.D0 * DG( 5 , JF )
      GENCNO = GENCNO + (1.D0 - DELTA ) * EPSH * 10.D0 * DG( 5 , JF )
      GENHE  = GENHE  + EPSA  * 10.D0 * DG( 5 , JF )
      GENCAR = GENCAR + ECAR  * 10.D0 * DG( 5 , JF )
      GENGRA = GENGRA + GRAVI * 10.D0 * DG( 5 , JF )
      GENNU  = GENNU  + EPNEU * 10.D0 * DG( 5 , JF )
      IF( JF . EQ . 1 ) ROC   = RHO
      IF( JF . EQ . 1 ) DELCE = DELTA
    
C********** STAMPA VARIE QUANTITA' OGNI MESHPOINTS **********************

      BETA=(CN3*TM**4)/PR 
      VSOUND=DSQRT(PP*1.D+17/RHO)
      
      AST( 1 , JF) = RR *1.d10
      AST( 2 , JF) = EM / EMTOT
      AST( 3 , JF) = TM
      AST( 4 , JF) = PR
      AST( 5 , JF) = RHO 
      AST( 6 , JF) = XXX(1,JF)
      AST( 7 , JF) = ZL*1.d32
      AST( 8 , JF) = CAP
      AST( 9 , JF) = EPNEU+EPSN+GRAVI 
      AST( 10, JF) = DAD
      AST( 11, JF) = 1.d0-XXX(1,JF)-XXX(2,JF)-XXX(3,JF)
      AST( 12, JF) = RTOT-RR*1.d10
      AST( 13, JF) = GRAVI
      AST( 14, JF) = GENGRA
      AST( 15, JF) = XXX(2,JF)
      AST( 16, JF) = XXX(3,JF)
      AST( 17, JF) = XXX(4,JF)
      AST( 18, JF) = XXX(10,JF)
      AST( 19, JF) = XXX(5,JF)
      AST( 20, JF) = XXX(6,JF)
      AST( 21, JF) = GRSAD      ! grad_rad - grad_ad
      AST( 22, JF) = DRAD(JF)
      AST( 23, JF) = G1
      AST( 24, JF) = CSP
      AST( 25, JF) = DEL
      AST( 26, JF) = EMUE
C***********************************************************************
c       WRITE(39,109) (AST(L,JF),L=1,22),XXX(8,JF)
c       WRITE(40,109) (AST(L,JF),L=1,6),XXX(3,JF),(AST(J,JF),J=23,26)

    1 CONTINUE         !FINE LOOP PRINCIPALE SU JF

       DO K=1,MAXNE-1
            XX(K)=XXX(K,MAXME-1)
         END DO
         BOT=1.D0-FRAZ
       EMSOL= EMTOT/(SUNMg33)

         CALL NATMOS(ELLOG,TEFF,EMSOL,XX,ZINI,BOT,ALPHA,
     $       O(1),O(2001),O(4001),O(6001),NATM,2,0,0,ISINGOLA,
     $        MAXNE,nmd)
c       write(*,200)ELLOG,TEFF,DLOG10(RAG)

      DO JJJ=npa-1,1,-1
       DO II=1,26
         AST(II,JF+NPA-JJJ-1)=BPP(II,JJJ)
       ENDDO
      ENDDO
      NTOT=JF+NPA-2

C***********************************************************************
C  a questo punto tutte le grandezze necessarie sono negli ast(*,*); 
C  NTOT é il numero di mesh totali, comprensivi anche della parte 
C  della fotosfera
C***********************************************************************
       LINT=0
c       write(42,'(i7,f12.8)')nmd,rtot/SUNRcm
       DO M=2,NTOT
      IF(M.EQ.NRI)THEN
        DO I=1,3
          FP(I)  = dlog10(AST(4,M+I-3))
          FRO(I) = dlog10(AST(5,M+I-3))
          FRDO(I) = dlog(AST(5,M+I-3))
            AX(I)  = dlog10(AST(1,M+I-3))
          ADX(I)  = dlog(AST(1,M+I-3))
          FM(I) = dlog(AST(2,M+I-3)*EMTOT)
c         write(*,'(i4,1p5e16.9)')m,FP(I),FM(I),ADX(I)
c         pause
        ENDDO
      ELSEIF(M.GE.10.AND.M.LE.NRI-10)THEN
        DO I=1,3
         ydat=0.d0
         do k1=-5,5
          ydat=ydat+dlog10(AST(5,M+k1+I-2))
         enddo
          FRO(I) = ydat/11.d0
          FP(I)  = dlog10(AST(4,M+I-2))
          FRDO(I) = dlog(AST(5,M+I-2))
            AX(I)  = dlog10(AST(1,M+I-2))
          ADX(I)  = dlog(AST(1,M+I-2))
          FM(I) = dlog(AST(2,M+I-2)*EMTOT)
        ENDDO
      ELSE
        DO I=1,3
          FP(I)  = dlog10(AST(4,M+I-2))
          FRO(I) = dlog10(AST(5,M+I-2))
          FRDO(I) = dlog(AST(5,M+I-2))
            AX(I)  = dlog10(AST(1,M+I-2))
          ADX(I)  = dlog(AST(1,M+I-2))
          FM(I) = dlog(AST(2,M+I-2)*EMTOT)
C         write(*,'(i4,1p5e16.9)')m,FP(I),FRO(I),AX(I)
C         pause
        ENDDO
        ENDIF

      RRR=DLOG10(AST(1,M))
            A2(M) =  (2.D0*RRR-(AX(2)+AX(3)))*FP(1)/
     #        ((AX(1)-AX(2))*(AX(1)-AX(3)))
        A2(M) = A2(M)+(2.D0*RRR-(AX(1)+AX(3)))*FP(2)/
     #          ((AX(2)-AX(1))*(AX(2)-AX(3)))
        A2(M) = A2(M)+(2.D0*RRR-(AX(1)+AX(2)))*FP(3)/
     #                ((AX(3)-AX(1))*(AX(3)-AX(2)))  

      A3(M) =  (2.D0*RRR-(AX(2)+AX(3)))*FRO(1)/
     #           ((AX(1)-AX(2))*(AX(1)-AX(3)))
        A3(M) = A3(M)+(2.D0*RRR-(AX(1)+AX(3)))*FRO(2)/
     #                ((AX(2)-AX(1))*(AX(2)-AX(3)))
        A3(M) = A3(M)+(2.D0*RRR-(AX(1)+AX(2)))*FRO(3)/
     #                ((AX(3)-AX(1))*(AX(3)-AX(2)))
        DP(M)=A2(M)/AST(23,M)-A3(M)
      DPP(M)=DP(M)

        IF(AST(21,M).GT.0.D0)THEN
          A5(M)=-GRAVC0*(AST(21,M)-AST(10,M))*AST(25,M)*AST(2,M)*EMTOT*
     #                  1.d33*AST(5,M)/(AST(4,M)*AST(1,M))
        ELSE
        A5(M)=-GRAVC0*(AST(22,M)-AST(10,M))*AST(25,M)*AST(2,M)*EMTOT*
     #                  1.d33*AST(5,M)/(AST(4,M)*AST(1,M))
        ENDIF
      b3=(dlog(AST(2,M+1)*EMTOT)-dlog(AST(2,M)*EMTOT))/
     #        (dlog(AST(1,M+1))-dlog(AST(1,M)))
      UU=4.d0*pi*AST(5,M)*AST(1,M)**3/(EMTOT*1.d33*AST(2,M))
       IF(b3.gt.2.998d0)LINT=M+1
c      write(42,109)DP(M),A5(M),DPP(M),a2(m),a3(m),
c     #        AST(1,M)/RTOT,AST(2,m),b3,UU
       ENDDO

       NHEST=0
       NMCC=0
       DO M=NTOT-10,2,-1
         IF(AST(6,M+6)-AST(6,M).GT.1.D-10)THEN
           NHEST=M
        GOTO 10
       ENDIF
       ENDDO
 10    CONTINUE
       DO M=1,NTOT
         IF(AST(6,M).gt.1.d-20)THEN
        NHINT=M
        goto 11
       ENDIF
       ENDDO
 11    CONTINUE
       AHEL=1.d0-(1.05d0*zini)
       DO M=NHINT,1,-1
         IF(AST(16,M).lt.AHEL)THEN
        NHEINT=M-1
        goto 12
       ENDIF
       ENDDO
 12    CONTINUE
       IF(AST(6,1).LE.1.D-30)THEN
         DO M=2,NHEINT+10
         IF(AST(21,M).GT.0)THEN
           NMCC=M
         ENDIF
       ENDDO
       ENDIF
c       IF(NMCC.NE.0)THEN
c         DRA2=(A2(NHEINT+6)-A2(NMCC-2))/
c     #          (AST(1,NHEINT+6)-AST(1,NMCC-2))
c         DRA3=(A3(NHEINT+6)-A3(NMCC-2))/
c     #          (AST(1,NHEINT+6)-AST(1,NMCC-2))
c         DO M=NMCC+1,NHEINT+5
c           A2(M)=A2(NMCC-2)+DRA2*(AST(1,M)-AST(1,NMCC))
c           A3(M)=A3(NMCC-2)+DRA3*(AST(1,M)-AST(1,NMCC))
c           DP(M)=A2(M)/AST(23,M)-A3(M)
c           DPP(M)=DP(M)
c         ENDDO
c       ENDIF
       
      IF(LINT.LT.3)LINT=3
       NRIGHE=NTOT-LINT
       DO M=LINT,NTOT
        IF(AST(6,1).GT.0.D0)THEN
         IF(M.GT.NHINT.AND.M.LT.NHINT+3)DP(M)=DPP(M)
         IF(M.GE.NHEST+10)DP(M)=A5(M)
        ELSE
         IF(NMCC.NE.0)THEN
           IF(M.LT.NMCC)DP(M)=A5(M)
         ELSE
         IF(M.LT.NHEINT)DP(M)=A5(M)
        ENDIF
         IF(M.GT.NHINT.AND.M.LT.NHINT+3)DP(M)=DPP(M)
         IF(M.GE.NHEST+10)DP(M)=A5(M)
        ENDIF
       ENDDO
c       WRITE(41,777)
c       WRITE(41,773)
c       DO M=LINT,NTOT
c         WRITE(41,109)DP(M),A5(M),DPP(M),AST(1,M)/RTOT,
c     #(AST(I,M),I=2,26)
c       ENDDO
C=======================================================================
C=== Computing the Brunt-Vaisala frequency N ===========================
C=======================================================================
      ASYPS1=0.D0 ! ASYMPTOTIC PERIOD SPACING FOR L=1
      write(*,'(d15.7,2i8)')ASYPS1,LINT,NTOT
      DO J=LINT,NTOT
       GLOC(J)=(GRAVC0*AST(2,J)*(EMTOT*1.d33))/(AST(1,J)**2) ! GRAV LOCALE
c       write(*,*) j,gloc(j),AST(2,J),AST(1,J)
C=== N^2=g*A5*(M/Mtot)/(r/Rtot)^3 - from Dennis           
cc       ENNEBV(J)=(GLOC(J)*DP(J)*AST(2,J))/((AST(1,J)/RTOT)**3)
       ENNEBV(J)=(GLOC(J)*DP(J))/AST(1,J)
       IF(ENNEBV(J).LT.0.D0) ENNEBV(J)=0.D0
       ENNEBV(J)=DSQRT(ENNEBV(J))    
      END DO
      
       DO J=LINT,NTOT-1 
        DELTARR= 2.0D0*(AST(1,J+1)-AST(1,J))/(AST(1,J+1)+AST(1,J))
        IF(ENNEBV(J).LT.1.D-20) GO TO 3 
        ASYPS1=ASYPS1 + (ENNEBV(J)*DELTARR)
  3    CONTINUE        
       END DO
      write(*,'(d15.7)')ASYPS1
       ASYPS1=((2.D0*PI**2)/DSQRT(2.D0))*1.0D0/ASYPS1
       
       IF(IRSTA.EQ.2) RETURN
       
C=======================================================================     
c   SCRITTURA FILE FORMATO FGONG:
C   SCRIVE LE PRIME 4 RIGHE: LA PRIMA IDENTIFICA IL NOME DEL MODELLO, DATA
C   DI CREAZIONE, CODICE USATO E LUOGO;
C   DALLA 2 ALLA 4 SONO DESCRITTIVE DEL MODELLO
C=======================================================================     
      CALL gtime(igio,imes,ianno)
      OPEN(UNIT=33,FILE=filename)
      STTAU='T(tau):  K-S 66  '  
      IF(INDIFFINI.EQ.0)THEN
        STDIF='-  Diffusion NO  '
      ELSE
        STDIF='-  Diffusion YES '        
      ENDIF
      IF(NSHUT.EQ.0)THEN
        STOVER='- Overshooting NO  '
      ELSE
        STOVER='- Overshooting YES '       
      ENDIF
      IF(ETA.LE.0.D0)THEN
        STML='- Mass Loss NO   '
      ELSE
        STML1='- Mass Loss '
        WRITE(APP,'(F5.3)') ETA
        READ(APP,'(A5)') EML
        STML=STML1//EML
      ENDIF
      write(33,110)FILENAME,ast(6,1),ast(16,1),rtot/SUNRcm,
     #igio,imes,ianno
      write(33,'(100a)')'EoS: EOS1 by IRWIN  - Opacity: OPAL (logT>4)  F
     #erguson et al. 2005  - Mixture: Caffau et al. 2011'
      write(33,'(3A17,A19)')STTAU,STML,STDIF,STOVER
      write(33,'(100a)')'Note: var(12,*)= -(partial ln rho(T,P)/ln T; 9.
     #999999999E+99 indicates not computed value'

 110  FORMAT(A23,' Model with Hc=',f8.4,1x,'Hec=',f8.4,1x,'R/Ro=',f8.4
     #,1x,' computed on',2I3,I5,'  by using FRANEC CODE 2.1.2 at INAF-OA
     #Te')
 112  FORMAT (1P5E16.9)
    
      ICONST=15
      IVERS=210
      IVAR=25
      FI=9.D0/4.D0
      EP=1.D0/162.D0
      BET=1.D0
      ALAM=1.D0
      VN=9.999999999D99
      
      DP1=(AST(4,2)-AST(4,1))/((AST(1,2)-AST(1,1)))
      DP2=(AST(4,3)-AST(4,2))/((AST(1,3)-AST(1,2)))
      DPRE=(DP2-DP1)/((AST(1,2)-AST(1,1)))
      DPRE=DPRE*RTOT**2/AST(4,1)
      
      DRO1=(AST(5,2)-AST(5,1))/((AST(1,2)-AST(1,1)))
      DRO2=(AST(5,3)-AST(5,2))/((AST(1,3)-AST(1,2)))
      DRORO=(DRO2-DRO1)/((AST(1,2)-AST(1,1)))
      DRORO=DRORO*RTOT**2/AST(5,1)

C santi aarhus - 050515      
c      ag=GRAVC0*(EMTOT*1.d33)/(RTOT**2)
      ag=GRAVC0
      
      
       write(33,'(4I10)')NRIGHE,ICONST,IVAR,IVERS    
       write(33,112)EMTOT*1.d33,RTOT,10**ELLOG*SUNLE,ZINI,XINI+absol(2),
     #ALPHA,FI,EP,BET,ALAM,DPRE,DRORO,TEMPO,10**TEFF,ag
      
       DO M=NTOT,LINT,-1
         ZZZ=1.D0-AST(6,M)-AST(16,M)
         write(33,112)AST(1,M),DLOG(AST(2,M)),(AST(L,M),L=3,9),
     #               AST(23,M),AST(10,M),AST(25,M),AST(24,M),AST(26,M),
     #               DP(M),VN,ZZZ,(AST(L,M),L=12,15),(AST(L,M),L=17,20)
       ENDDO
      
      
      
C ********************* STAMPA FISICA VARIA ****************************
 200  FORMAT(' L:',F10.7,' Te:',F10.7,' R:',F10.7)
 773  FORMAT ('=========================================================
     #==================================================================
     #==================================================================
     #==================================================================
     #==================================================================
     #==============================================')

 774  FORMAT ('      R (cm)           M/Mtot           Temp             
     #Pre              Ro                  H              L             
     #    k                E              Grad_ad            Z          
     #   R-r              E_g                L_g                He3     
     #     He4              C12               C13                 N14   
     #         O16          GRSAD            G_RAD       Li7')
 775  FORMAT ('=========================================================
     #==================================================================
     #==================================================')
 776  FORMAT ('      R (cm)           M/Mtot           Temp             
     #Pre              Ro                  H               He           
     # Gamma1               CP             Der           1/MU_e')
 777  FORMAT ('        DP               A5               DPP            
     #  R (cm)        M/Mtot           Temp             Pre             
     # Ro                    H              L                 k         
     #     E              Grad_ad            Z             R-r          
     #    E_g                L_g                He3          He4        
     #      C12               C13                 N14            O16    
     #        GRSAD            G_RAD        GAM_1              CP       
     #   1/Mu_e')
 109  FORMAT (1P,30D17.9)
      RETURN
      END
