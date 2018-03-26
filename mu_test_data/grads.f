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

      PROGRAM GRADS

      IMPLICIT REAL*8(A-H,O-Z)
C      include 'CODE_alw/maincom.2p1'
      include 'consts.2p1'
      character*30 abund_in,model_in,head_in,file_out,eq_out
      integer*4 kcount,MAXME,NMODEL,KMESH,MESH
      PARAMETER (MAXME=1424)
      real*8 Kthm
      Real*8 XI(MAXME),XHE3(MAXME),XHE4(MAXME),XC(MAXME),XN(MAXME),
     # XO(MAXME),XD(MAXME),XLi7(MAXME),XBe7(MAXME),XC13(MAXME),
     # XN15(MAXME),XO17(MAXME),XO18(MAXME),XFe56(MAXME),x1mas(MAXME),
     # MFRAC(MAXME),rfrac(MAXME),logP(MAXME),logT(MAXME),
     # logDens(MAXME),lum(MAXME),enuc(MAXME),egrav(MAXME),
     # eneu(MAXME),DELTA(MAXME),RADIAT(MAXME),ADIAB(MAXME),
     # VSOUND(MAXME),opac(maxme),cp(maxme),RHO_PHYS(MAXME),
     # RHO_N(MAXME)
c      WRITE(*,*) 'HELLO WORLD'


C NOTE ABOUT READING:
C For '_model' files, copy from NoMachine FIRST, then delete all
C strings from the file BEFORE running grads2p1p4.go

c    CHARACTER(*),PARAMETER :: fileplace = "FeHp006_alw/mu_test_data/"
      abund_in = 'diff_mu_logL=2_1231_abund'
      model_in = 'diff_mu_logL=2_1231_model'
C      head_in = 'diff_mu_logL=2_1231_header'
      file_out = 'diff_mu_logL=2_1231_grads'
      eq_out = 'diff_mu_logL=2_1231_eq'
C fileplace//      
c      OPEN(UNIT=21,FILE=head_in,status='OLD')
      OPEN(UNIT=28,FILE=abund_in,status='OLD')
      OPEN(UNIT=29,FILE=model_in,status='OLD')
      OPEN(UNIT=38,FILE=file_out,status='REPLACE')
      OPEN(UNIT=44,FILE=eq_out,status='REPLACE')
      
      WRITE(*,*) file_out
      write(38,110)
      write(38,120)
 110  FORMAT(1H1,/,'      Mr            MU           d(mu)/dr       dln(
     &mu)/dr   dln(mu)/dln(p)  d(X_Li7)/dr   d(X_C)/dr      d(X_N)/dr  
     &   d(X_He3)/dr     d(X_He4)/dr    Dthm_0   MESH')

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
C       if (km .eq. 1) then
C	 write(*,*) 'km = 1' 
C       else 
C       write(*,*) x1mas(k),',',XI(k),',',mfrac(k),',',rfrac(k)
       kcount = kcount + 1
      END DO
   40 WRITE(*,*)'Reading finished: number of lines: ', kcount
C      WRITE(*,*) x1mas(1)

C CALCULATION DO-LOOP
      DO K = 1,MAXME

       KM = K - 1
       KP = K + 1

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

C TAKE DIFFERENCE OF THE MU'S
C CENTRAL DERIVATIVE METHOD - CHOSEN METHOD FOR SINGLE DERIVATIVE
      DMU = FINAL_MU_F - FINAL_MU_B

C Focus on individual isotopes
C CENTRAL DERIVATIVE METHOD
      DXC = XC(KP) - XC(KM)
      DXLi7 = XLi7(KP) - XLi7(KM)
      DXN = XN(KP) - XN(KM)
      DXHE3 = XHE3(KP) - XHE3(KM)
      DXHE4 = XHE4(KP) - XHE4(KM)

C CALCULATE THE (LINEAR) MOLECULAR GRADIENT
C BY GETTING THE KTH-LAYER WIDTHS FOR DIFFERENT PARAMETERS

C CENTRAL DERIVATIVE METHOD
      DR = (2.21E+01)*(RFRAC(KP)-RFRAC(KM))
      DLNP = DLOG(10**logP(KP)) - DLOG(10**logP(KM))
      DLOGP = logP(kp)-logP(km)
      DLNT = DLOG(10**logT(KP)) - DLOG(10**logT(KM))
      DLOGT = logT(kp)-logT(km)
C      DMASS = ()


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

C THERMOHALINE MIXING IMPLEMENTATION
C LEDOUX CRITERION (CANTIELLO, LANGER, 2010):
C  0 >= DEL-DEL_AD >= (PHI/DELTA)*DEL_MU
C WITH PHI = PART(DLN(RHO)/DLN(MU))_F,T
c AND DELTA = PART(DLN(RHO)/DLN(T))_F,MU

c Calculate thermohaline diffusion coefficient
C conversions to get physical parameters in cgs units
C Density conversion
      density = 10**(logDens(k))
C CN3 = (radiation constant)/3 in cgs units (erg/(cm**3 K**4))
C Clight = light speed in cgs units (cm/s**2)
C Temperature conversion
      temp = 10**(logT(k))
C Heat capacity cp(k) comes from stampe 'model' file

c First, get thermal diffusivity (cgs units)
      Kthm = (4*CN3*Clight*(temp**3))/(OPAC(k)*(density**2)*cp(k))

C Final result (Cantiello & Langer (2010) & Maeder et al. (2013))
C with (phi/delta) = 1 (must always be +ve from first principles)
      Cthl = 1000 ! Maeder et al. (2013), from Ulrich (1972)
C If loop comes from Lattanzio et al. (2015): D = 0 unless del_MU < 0
C      if ((DEL_MU .LT. 0) .AND. (DEL_MU .GT. (RADIAT(K)-ADIAB(K))))
C     # then
C        Dthl = Cthl*Kthm*(DEL_MU)/(RADIAT(K)-ADIAB(K))
C      else
C        Dthl = 0.0
C      end if

      if (DEL_MU .LT. 0) then
        Dthl_0 = Cthl*Kthm*(DEL_MU)/(RADIAT(K)-ADIAB(K))
      else
        Dthl_0 = 0.0
      end if

C DIFFUSION EQUATION CALCULATION - d2(X_i)/dr2 construction (linear)
C using concept & equations from:
C https://me.ucsb.edu/~moehlis/APC591/tutorials/tutorial5/node3.html
C and:
C www.math.utep.edu/Faculty/sewell/tamu/664/prob4b.pdf

C diffusion equation uses number density as the variable being
C differetiated - comes from continuity equation. For 1 species, with
C no fusion modelled, can represent the number density with (X_i)*rho,
C since each particle of the species has a constant mass

C First, get the physical density from the log(rho) data from stampe
C and convert to element-specific masses: test case = N14
      RHO_N(KM) = XN(KM)*(10**logDens(KM))
      RHO_N(K) = XN(K)*(10**logDens(K))
      RHO_N(KP) = XN(KP)*(10**logDens(KP))

C Second, calulate the numerator
C Numerator = difference between forward and backward differentials
C i.e., (dT/dr)_(k+1/2) - (dT/dr)_(k-1/2)
C Define backward, forward difference in T
      DRHO_N_B = RHO_N(K) - RHO_N(KM)
      DRHO_N_F = RHO_N(KP) - RHO_N(K)

      DXN_F = XN(KP) - XN(K)
      DXN_B = XN(K) - XN(KM)

C Define backward, forward differences in physical distance, in this
C (1D) case the radius, r
      DR_B = (2.21E+01*SUNRcm)*(rfrac(k) - rfrac(km))
      DR_F = (2.21E+01*SUNRcm)*(rfrac(kp) - rfrac(k))

C Divide to get the numerator components, then take their difference
C to get the final numerator
      CALL AVOIDNAN(DRHO_N_DR_B,DRHO_N_B,DR_B)
      CALL AVOIDNAN(DRHO_N_DR_F,DRHO_N_F,DR_F)
      D2RHO_N_DR = DRHO_N_DR_F - DRHO_N_DR_B

      CALL AVOIDNAN(DXN_DR_B,DXN_B,DR_B)
      CALL AVOIDNAN(DXN_DR_F,DDXN_F,DR_F)
      D2XN_DR = DXN_DR_F - DXN_DR_B

C Denominator = central difference for r
      DR_2 = (DR*SUNRcm)/2

C Combine the numerator and denominator to obtain the
C second-order differential
      CALL AVOIDNAN(D2RHO_N_DR2,D2RHO_N_DR,DR_2)

      CALL AVOIDNAN(D2XN_DR2,D2XN_DR,DR_2)

C timestep calculation, using logL = 2.1231 and ages at (k-1),(k+1),
C then divide by 2 to approximate to (k-1/2),(k+1/2), then divide by 100
c      dtime_yr = ((10**10.10616736) - (10**10.10616654))/(2*100)
c      dtime_s = dtime_yr*secanni

C get a (fixed) test value of D_therm for the equation
      Dthm_eq = 0.1457D+06

C Diffusion equation - gives change in physical
C density rho_N(r,t) due to time only, in cgs
      dRHO_N_dt = Dthm_eq*D2RHO_N_DR2
C get 'normalising' time, tau = Rtot^2/D in years
      t_norm = ((2.21E+01*SUNRcm)**2)/(Dthm_eq*secanni)

C write out numerical results
      if (k .eq. MAXME) then
	write(*,*) 'km = (MAXME-1) = ', km
        write(*,*) RHO_N(km),RHO_N(k)
      end if

      write(38,102) x1mas(K),FINAL_MU_K,DMU_DR,DLNMU_DR,DEL_MU,
     # DXLi7_DR,DXC_DR,DXN_DR,DXHE3_DR,DXHE4_DR,Dthl_0,K
      
      if (k .eq. 1) then
       write(44,103) Dthm_eq,t_norm
       write(44,104)
       write(44,105)
      end if

      if (dr_2 .eq. 0.0) then
       write(*,*) rfrac(kp),rfrac(km)
      end if

      write(44,106) x1mas(K),rfrac(k),XN(k),D2XN_DR2,DRHO_N_F,DRHO_N_B,
     # DR_2,D2RHO_N_DR,D2RHO_N_DR2,dRHO_N_dt,K

 100  FORMAT(1X,F12.9,14D14.7,105X,I5)
 101  FORMAT(F12.9,F10.7,F10.6,F8.5,F8.4,D13.5,3D11.3,D9.2,
     &D11.3,F6.3,D12.4,D12.4,D9.2,I5)
 102  FORMAT(1X,F12.9,9D15.7,D12.4,I5)

 103  FORMAT('Dthm = ',D11.4,', normalisation time = ',D11.4)
 104  FORMAT('      Mr         R/R*         X_N          d2X_N/dr2      
     #drho_N_rf       drho_N_rb        dr          d2rho_N/dr    d2rho_N
     #/dr2    drho_N/dt   MESH')
 105  FORMAT('==========================================================
     #==================================================================
     #=========================')
 106  FORMAT(1X,F12.9,F10.7,8D15.7,I5)

c ,XI(kP),XI(kM),KP ,2D15.7,I5 ,DR ,DLNP_2
      
      enddo
   

      WRITE(*,*) 'END OF WRITING'
      close(38)

      END

