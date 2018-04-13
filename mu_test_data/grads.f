C*****************************************************************************
C ROUTINE FOR PLOTTING GRADIENTS AND GRAPHS OF STAMPE DATA
C
C Author: Alex Lisboa-Wright
C N.B.: Directory for input & output is FeHp006_alw/mu_test_data
C
C*****************************************************************************
C SUBROUTINE TO AVOID GRADIENT NANS DUE TO 0/0
      SUBROUTINE AVOIDNAN(RES,NUM,DEN)
       Real*8 NUM,DEN,RES ! NUMERATOR, DENOMINATOR, RESULT OF DIVISION
       IF (NUM .EQ. 0.0 .AND. DEN .EQ. 0.0) THEN
        RES = 0.0
       ELSE
        RES = NUM/DEN
       END IF
      END



C **********************************************************************
C          MAIN PROGRAM
C **********************************************************************
      PROGRAM GRADS

      IMPLICIT REAL*8(A-H,O-Z)
C      include 'CODE_alw/maincom.2p1'
      include 'consts.2p1'
      character*30 abund_in,model_in,head_in,file_out,eq_out,time_out
      integer*4 kcount,MAXME,NMODEL,KMESH,MESH
      PARAMETER (MAXME=1424)
      PARAMETER (NSTEPS = 100)
      real*8 Kthm
      Real*8 XI(MAXME),XHE3(MAXME),XHE4(MAXME),XC(MAXME),XN(MAXME),
     # XO(MAXME),XD(MAXME),XLi7(MAXME),XBe7(MAXME),XC13(MAXME),
     # XN15(MAXME),XO17(MAXME),XO18(MAXME),XFe56(MAXME),x1mas(MAXME),
     # MFRAC(MAXME),rfrac(MAXME),logP(MAXME),logT(MAXME),
     # logDens(MAXME),lum(MAXME),enuc(MAXME),egrav(MAXME),
     # eneu(MAXME),DELTA(MAXME),RADIAT(MAXME),ADIAB(MAXME),
     # VSOUND(MAXME),opac(maxme),cp(maxme),DXN_DR_K(MAXME)
     # ,DEL_MU_K(MAXME), XN_KT_OUT(NSTEPS),XN_KT_OUT_VAR(NSTEPS)

C NOTE ABOUT READING:
C For '_model' files, copy from NoMachine FIRST, then delete all
C strings from the file BEFORE running grads2p1p4.go

c    CHARACTER(*),PARAMETER :: fileplace = "FeHp006_alw/mu_test_data/"
      abund_in = 'diff_mu_logL=2_1231_abund'
      model_in = 'diff_mu_logL=2_1231_model'
C      head_in = 'diff_mu_logL=2_1231_header'
      file_out = 'diff_mu_logL=2_1231_grads'
      eq_out = 'diff_mu_logL=2_1231_eq'
      time_out = 'diff_mu_logL=2_1231_time'
C fileplace//      
c      OPEN(UNIT=21,FILE=head_in,status='OLD')
      OPEN(UNIT=28,FILE=abund_in,status='OLD')
      OPEN(UNIT=29,FILE=model_in,status='OLD')
      OPEN(UNIT=38,FILE=file_out,status='REPLACE')
      OPEN(UNIT=44,FILE=eq_out,status='REPLACE')
      OPEN(UNIT=45,FILE=time_out,status='REPLACE')
      
      WRITE(*,*) file_out
      write(38,110)
      write(38,120)
 110  FORMAT(1H1,/,'      Mr            MU           d(mu)/dr       dln(
     &mu)/dr   dln(mu)/dln(p)  d(X_Li7)/dr   d(X_C)/dr      d(X_N)/dr  
     &   d(X_He3)/dr     d(X_He4)/dr    Dthl     MESH')

Cx1mas,MU,DMU_DR ,DLNMU_DR ,DEL_MU ,
C     # DXLi7_DR ,DXC_DR ,DXN_DR ,KM
 120  FORMAT('==========================================================
     #==================================================================
     #===========================================================')

      kcount = 0
C READING DO-LOOP
 
      DO K = 1,MAXME
       READ(28,100,end=40) x1mas(k),XI(k),XD(k),XHE3(k),XHE4(k),
     # XLi7(k),XBe7(k),XC(k),XC13(k),XN(k),XN15(k),XO(k),
     # XO17(k),XO18(k),XFe56(k),Kmesh
C       write(*,*) Kmesh
       
C MU,DMU_DR,DLNMU_DR,DEL_MU,DXLi7_DR,# DXC_DR,DXN_DR,KM
       READ(29,101,end=40) mfrac(k),rfrac(k),logP(k),logT(k),
     # logDens(k),lum(k),enuc(k),egrav(k),eneu(k),DELTA(k),
     # RADIAT(k),ADIAB(k),opac(k),cp(k),VSOUND(k),MESH
       kcount = kcount + 1
      END DO
   40 WRITE(*,*)'Reading finished: number of lines: ', kcount
C      WRITE(*,*) x1mas(1)

C CALCULATION DO-LOOP
      DO K = 1,MAXME

       KM = K - 1
       KP = K + 1

C SINGLE-MESH QUANTITIES

C    CALCULATE FINVMU = INVERSE OF FINAL_MU = 1/MU = 1/MU_I + 1/MU_E
C SO: 1/MU = SUM_I ((Z_I + 1)*X_I/A_I), ASSUMING FULL IONIZATION
C N.B.: FINVMU_B DEFINED AT KM = K - 1
      FINVMU_B = (XI(km)*2/1)+(XHE3(km)*3/3)+(XHE4(km)*3/4)
     $ +(XC(km)*7/12)+(XN(km)*8/14)+(XO(km)*9/16)+(XD(km)*2/2)
     $ +(XLi7(km)*4/7)+(XBe7(km)*5/7)+(XC13(km)*7/13)+(XN15(km)*8/16)
     $ +(XO17(km)*9/17)+(XO18(km)*9/18)+(XFe56(km)*26/56)

C CALCULATE MEAN MOLECULAR WEIGHT, FINAL_MU_B, BY INVERTING FINVMU
      FINAL_MU_B = 1/FINVMU_B

C DEFINE FINVMU_K, THE EQUIVALENT OF FINVMU, BUT AT K = KM + 1
C errorbait
      FINVMU_K = (XI(k)*2/1)+(XHE3(k)*3/3)+(XHE4(k)*3/4)
     $ +(XC(k)*7/12)+(XN(k)*8/14)+(XO(k)*9/16)+(XD(k)*2/2)
     S +(XLi7(k)*4/7)+(XBe7(k)*5/7)+(XC13(k)*7/13)+(XN15(k)*8/16)
     $ +(XO17(k)*9/17)+(XO18(k)*9/18)+(XFe56(k)*26/56)

C IN THE SAME WAY, DEFINE FINAL_MU_K
      FINAL_MU_K = 1/FINVMU_K

C DEFINE FINVMU_F, THE EQUIVALENT OF FINVMU, BUT AT KP = K + 1
      FINVMU_F = (XI(kp)*2/1)+(XHE3(kp)*3/3)+(XHE4(kp)*3/4)
     $ +(XC(kp)*7/12)+(XN(kp)*8/14)+(XO(kp)*9/16)+(XD(kp)*2/2)
     S +(XLi7(kp)*4/7)+(XBe7(kp)*5/7)+(XC13(kp)*7/13)+(XN15(kp)*8/16)
     $ +(XO17(kp)*9/17)+(XO18(kp)*9/18)+(XFe56(kp)*26/56)

C IN THE SAME WAY, DEFINE FINAL_MU_F
      FINAL_MU_F = 1/FINVMU_F

C MULTIPLE-MESH QUANTITIES - NEED TO TREAT K = 1 & MAXME AS SPECIAL
C CASES
      if (k .eq. 1) then
C forward derivatives for k = 1
C TAKE DIFFERENCE OF THE MU'S
      DMU = FINAL_MU_F - FINAL_MU_K
C Focus on individual isotopes
      DXC = XC(KP) - XC(K)
      DXLi7 = XLi7(KP) - XLi7(K)
      DXN = XN(KP) - XN(K)
      DXHE3 = XHE3(KP) - XHE3(K)
      DXHE4 = XHE4(KP) - XHE4(K)
C CALCULATE THE (LINEAR) MOLECULAR GRADIENT
C BY GETTING THE KTH-LAYER WIDTHS FOR DIFFERENT PARAMETERS
      DR = (2.21E+01)*(RFRAC(KP)-RFRAC(K))
      DLNP = DLOG(10**logP(KP)) - DLOG(10**logP(K))
      DLOGP = logP(kp)-logP(k)
      DLNT = DLOG(10**logT(KP)) - DLOG(10**logT(K))
      DLOGT = logT(kp)-logT(k)
C      DMASS = ()

      elseif (k .eq. MAXME) then
C backward derivatives for k = MAXME
C TAKE DIFFERENCE OF THE MU'S
      DMU = FINAL_MU_K - FINAL_MU_B
C Focus on individual isotopes
      DXC = XC(K) - XC(KM)
      DXLi7 = XLi7(K) - XLi7(KM)
      DXN = XN(K) - XN(KM)
      DXHE3 = XHE3(K) - XHE3(KM)
      DXHE4 = XHE4(K) - XHE4(KM)
C CALCULATE THE (LINEAR) MOLECULAR GRADIENT
C BY GETTING THE KTH-LAYER WIDTHS FOR DIFFERENT PARAMETERS
      DR = (2.21E+01)*(RFRAC(K)-RFRAC(KM))
      DLNP = DLOG(10**logP(K)) - DLOG(10**logP(KM))
      DLOGP = logP(k)-logP(km)
      DLNT = DLOG(10**logT(K)) - DLOG(10**logT(KM))
      DLOGT = logT(k)-logT(km)

      else
C CENTRAL DERIVATIVE METHOD FOR REST
C TAKE DIFFERENCE OF THE MU'S
      DMU = FINAL_MU_F - FINAL_MU_B
C Focus on individual isotopes
      DXC = XC(KP) - XC(KM)
      DXLi7 = XLi7(KP) - XLi7(KM)
      DXN = XN(KP) - XN(KM)
      DXHE3 = XHE3(KP) - XHE3(KM)
      DXHE4 = XHE4(KP) - XHE4(KM)
C CALCULATE THE (LINEAR) MOLECULAR GRADIENT
C BY GETTING THE KTH-LAYER WIDTHS FOR DIFFERENT PARAMETERS
      DR = (2.21E+01)*(RFRAC(KP)-RFRAC(KM))
      DLNP = DLOG(10**logP(KP)) - DLOG(10**logP(KM))
      DLOGP = logP(kp)-logP(km)
      DLNT = DLOG(10**logT(KP)) - DLOG(10**logT(KM))
      DLOGT = logT(kp)-logT(km)
C      DMASS = ()
      endif

C TAKE DIFFERENCE OF NATURAL LOG OF THE MU'S
C BACKWARD DERIVATIVE METHOD
C      DLNMU_B = DLOG(FINAL_MU_K) - DLOG(FINAL_MU_B)
C      DLOGMU_B = DLOG10(FINAL_MU_K) - DLOG10(FINAL_MU_B)
C FORWARD DERIVATIVE METHOD
C      DLNMU_F = DLOG(FINAL_MU_F) - DLOG(FINAL_MU_K)
C      DLOGMU_F = DLOG10(FINAL_MU_F) - DLOG10(FINAL_MU_K)
C CENTRAL DERIVATIVE METHOD
      DLNMU = DLOG(FINAL_MU_F) - DLOG(FINAL_MU_B)
      DLOGMU = DLOG10(FINAL_MU_F) - DLOG10(FINAL_MU_B)

C DEFINE RELEVANT LINEAR GRADIENTS
C CENTRAL DERIVATIVE METHOD
      CALL AVOIDNAN(DMU_DR,DMU,DR) ! d(mu)/dr
C      CALL AVOIDNAN(DMU_DM,DMU,DMASS) ! d(mu)/dm
      CALL AVOIDNAN(DEL_MU,DLNMU,DLNP) ! d(ln(mu))/d(ln(P))
      CALL AVOIDNAN(DLNMU_DR,DLNMU,DR) ! d(ln(mu))/dr
      CALL AVOIDNAN(DXC_DR,DXC,DR) ! C12
      CALL AVOIDNAN(DXN_DR,DXN,DR) ! N14
      CALL AVOIDNAN(DXLi7_DR,DXLi7,DR) ! Li7
      CALL AVOIDNAN(DXHE3_DR,DXHE3,DR) ! He3
      CALL AVOIDNAN(DXHE4_DR,DXHE4,DR) ! He4
c assign thermohaline-relevant numbers to arrays
      DXN_DR_K(K) = (DXN_DR/SUNRcm)
      DEL_MU_K(K) = DEL_MU
C write out numerical results to files

      write(38,102) x1mas(K),FINAL_MU_K,DMU_DR,DLNMU_DR,DEL_MU,
     # DXLi7_DR,DXC_DR,DXN_DR,DXHE3_DR,DXHE4_DR,Dthl,K
      
 100  FORMAT(1X,F12.9,14D14.7,105X,I5)
 101  FORMAT(F12.9,F10.7,F10.6,F8.5,F8.4,D13.5,3D11.3,D9.2,
     &D11.3,F6.3,D12.4,D12.4,D9.2,I5)
 102  FORMAT(1X,F12.9,9D15.7,D12.4,I5)

c ,XI(kP),XI(kM),KP ,2D15.7,I5 ,DR ,DLNP_2
      enddo
   
C Thermohaline mixing loop
      do k = 1,MAXME
c       write(*,*) 'XN(1),XN(2) = ',XN(1),',',XN(2)
       km = k - 1
       kp = k + 1
C THERMOHALINE MIXING IMPLEMENTATION

C get extant parameters into loop 
c physical radius in cgs
      rphys_m = (2.21E+01*SUNRcm)*(RFRAC(KM))
      rphys_k = (2.21E+01*SUNRcm)*(RFRAC(K))
      rphys_p = (2.21E+01*SUNRcm)*(RFRAC(KP))

c density
      dens_m = 10**logDens(KM)
      dens_k = 10**logDens(K)
      dens_p = 10**logDens(KP)

      th_dxn_dr_m = DXN_DR_K(KM)
      th_dxn_dr_k = DXN_DR_K(K)
      th_dxn_dr_p = DXN_DR_K(KP)


C Special cases for boundary conditions: dx/dr = 0 at rfrac = 0,1

      if (k .eq. 1) then
       th_dxn_dr_m = 0
       th_dxn_dr_k = 0
      elseif(k .eq. 2) then
       th_dxn_dr_m = 0
      elseif(k .eq. (MAXME-1)) then
       th_dxn_dr_p = 0
      elseif(k .eq. MAXME) then
       th_dxn_dr_p = 0
       th_dxn_dr_k = 0
      endif



C CN3 = (radiation constant)/3 in cgs units (erg/(cm**3 K**4))
C Clight = light speed in cgs units (cm/s**2)
C Temperature conversion
      temp = 10**(logT(k))
C Heat capacity cp(k) comes from stampe 'model' file

c First, get thermal diffusivity (cgs units)
      Kthm = (4*CN3*Clight*(temp**3))/(OPAC(k)*(dens_k**2)*cp(k))

C Final result (Cantiello & Langer (2010) & Maeder et al. (2013))
C with (phi/delta) = 1 (must always be +ve from first principles)
      Cthl = 1000 ! Maeder et al. (2013), from Ulrich (1972)


C LEDOUX CRITERION (CANTIELLO, LANGER, 2010):
C  0 >= DEL-DEL_AD >= (PHI/DELTA)*DEL_MU
C WITH PHI = PART(DLN(RHO)/DLN(MU))_F,T
c AND DELTA = PART(DLN(RHO)/DLN(T))_F,MU
C If loop comes from Lattanzio et al. (2015): D = 0 unless del_MU < 0
      if ((DEL_MU_K(K) .LT. 0) .AND. (DEL_MU_K(K) .GT.
     # (RADIAT(K)-ADIAB(K)))) then
        Dthl = Cthl*Kthm*(DEL_MU_K(K))/(RADIAT(K)-ADIAB(K))
c        write(*,*) 
      else
        Dthl = 0.0
      end if

      if (DEL_MU_K(K) .LT. 0) then
        Dthl_0 = Cthl*Kthm*(DEL_MU_K(K))/(RADIAT(K)-ADIAB(K))
      else
        Dthl_0 = 0.0
      end if

      if (Dthl_0 .ne. Dthl) then
        write(*,*) 'Ledoux-Schwartzschild criterion mismatch of: '
        write(*,*) Dthl_0, Dthl,Kthm
c        write(*,*) temp,opac(k),dens_k,cp(k), k
        write(*,*) DEL_MU_K(k),(RADIAT(K)-ADIAB(K)),k
      end if

C DIFFUSION EQUATION CALCULATION - construction (linear)
C using concept & equations from:
C https://me.ucsb.edu/~moehlis/APC591/tutorials/tutorial5/node3.html
C and:
C www.math.utep.edu/Faculty/sewell/tamu/664/prob4b.pdf

C CORRECT WAY!!!
C dX/dt = (1/rho*r*r)*d/dr(D*rho*r*r*dX/dr)
C Test case: Nitrogen-14

C get a (fixed) test value of D_thm for the equation
      Dthl_fixed = 0.1457D+06

C N.B.: IF CHANGING ELEMENT(S), DO IT HERE!
c collect terms in rightmost bracket for k,k+1,k-1 into 1 variable
      BR_PROD_M = ((rphys_m**2)*dens_m*Dthl_fixed)*th_dxn_dr_m
c      BR_PROD_K = DXN_DR*(rphys_k*rphys_k*dens_k*Dthl_fixed)
      BR_PROD_P = ((rphys_p**2)*dens_p*Dthl_fixed)*th_dxn_dr_p
      
c variable coefficient
      BR_PROD_M_var = ((rphys_m**2)*dens_m*Dthl)*th_dxn_dr_m
c      BR_PROD_K_var = DXN_DR*(rphys_k*rphys_k*dens_k*Dthl_fixed)
      BR_PROD_P_var = ((rphys_p**2)*dens_p*Dthl)*th_dxn_dr_p

c central difference for second radius-differential
      call avoidnan(dbrprod_dr,(BR_PROD_P-BR_PROD_M),
     # (rphys_p-rphys_m))
c same for variable coefficient
      call avoidnan(dbrprod_dr_var,(BR_PROD_P_var-BR_PROD_M_var),
     # (rphys_p-rphys_m))

c      IF (k .eq. 1 .or. k .eq. MAXME) then
c       write(*,*) 'k = ',k
c       write(*,*) mfrac(kp),mfrac(k),mfrac(km)
c       write(*,*) x1mas(kp),x1mas(k),x1mas(km)
c       write(*,*) rfrac(kp),rfrac(k),rfrac(km),dr
c       write(*,*) XN(kp),XN(k),XN(km),DXN,DXN_DR
c      ELSEIF (k .eq. (MAXME)) then
c       write(*,*) rfrac(kp),rfrac(k),rfrac(km),dr
c      ELSE
c      END IF

C timestep size in years, then seconds
      dt_yr = 500
      dt_s = dt_yr*secanni

C Diffusion equation - gives change in mass fraction
c due to time only, in cgs
      dXN_dt = dbrprod_dr/((rphys_k**2)*dens_k)
      dXN_dt_var = dbrprod_dr_var/((rphys_k**2)*dens_k)
C density of N in layer at time dt later
      XN_K_tpdt = XN(K) + dt_s*dXN_dt
C use non-fixed coefficients
      XN_K_tpdt_var = XN(K) + dt_s*dXN_dt_var

      if ((XN_K_tpdt_var .ne. XN(K)) .and. Dthl .eq. 0) then
       write(*,*) 'equation output errors',k
      endif

      if (k .eq. MAXME) then
	write(*,*) 'k = MAXME = ', k
      end if

      if (k .eq. 1) then
       write(44,103) Dthl_fixed,nsteps,dt_s,dt_yr
       write(44,104)
       write(44,105)

       write(45,130)
       write(45,131)
      end if
C write '_eq' data to file
      write(44,106) x1mas(K),rfrac(k),Dthl,BR_PROD_M,BR_PROD_M_var,
     # BR_PROD_P,BR_PROD_P_var,XN(k),XN_K_tpdt,XN_K_tpdt_var,K


C Write header to '_time' file
C time variation: test with 'nsteps' time-steps
C first step already done above: introduce iterating (array) variable
      XN_KT_OUT(1) = XN_K_tpdt
      XN_KT_OUT_var(1) = XN_K_tpdt_var

      do n = 2,nsteps
       XN_KT_OUT(N) = XN_KT_OUT(N-1) + dt_s*dXN_dt
       XN_KT_OUT_var(N) = XN_KT_OUT_var(N-1) + dt_s*dXN_dt_var
       if (XN_KT_OUT_var(N) .le. 0) then
        XN_KT_OUT_var(N) = 0
       elseif (XN_KT_OUT(N) .le. 0) then
        XN_KT_OUT(N) = 0
       endif
      enddo
C write
C FORMAT FROM BGT2P1.F: NMOD , AGE , ( A(L) , L =  1 , 3 )
      WRITE(45,132) x1mas(K),XN(K),XN_KT_OUT(1),(XN_KT_OUT(L), L =  10,
     #100,10 ),K

C '_eq' file format statements
 103  FORMAT('fixed Dthl =',D11.4,' yr, number of steps = ',I4,', step s
     #ize =',D11.4,' sec =',D11.4,' yr' )
 104  FORMAT('      Mr         R/R*         Dthl        BR_PROD_M    BR_
     #PROD_M_var    BR_PROD_P     BR_PROD_P_var     XN(t)        XN(t+dt
     #)     XN(t+dt)_var  MESH')
 105  FORMAT('==========================================================
     #==================================================================
     #=========================')
 106  FORMAT(1X,F12.9,F10.7,8D15.7,I5)

C '_time' file format statements
 130  FORMAT('      Mr            XN(t)        XN(t+dt)      XN(t+10dt) 
     #    XN(t+20dt)     XN(t+30dt)     XN(t+40dt)     XN(t+50dt)     XN
     #(t+60dt)    XN(t+70dt)    XN(t+80dt)     XN(t+90dt)    XN(t+100dt)
     #   MESH')
 131  FORMAT('==========================================================
     #==================================================================
     #==================================================================
     #================')
 132  FORMAT(1X,F12.9,12D15.7,I5)

      enddo
      

      WRITE(*,*) 'END OF WRITING'
c      close(38)

      END

