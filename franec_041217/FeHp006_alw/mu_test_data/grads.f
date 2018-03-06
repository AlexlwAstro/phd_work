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
c      subroutine gradient

      IMPLICIT REAL*8(A-H,O-Z)
C      include 'CODE_alw/maincom.2p1'
      include 'consts.2p1'
      character*30 abund_in,model_in,head_in,file_out
      integer*4 kcount,MAXME,NMODEL,KMESH,MESH,x,xm,xp
      PARAMETER (MAXME=1424)
      real*8 Kthm
      Real*8 XI(MAXME),XHE3(MAXME),XHE4(MAXME),XC(MAXME),XN(MAXME),
     # XO(MAXME),XD(MAXME),XLi7(MAXME),XBe7(MAXME),XC13(MAXME),
     # XN15(MAXME),XO17(MAXME),XO18(MAXME),XFe56(MAXME),x1mas(MAXME),
     # MFRAC(MAXME),rfrac(MAXME),logP(MAXME),logT(MAXME),
     # logDens(MAXME),lum(MAXME),enuc(MAXME),egrav(MAXME),
     # eneu(MAXME),DELTA(MAXME),RADIAT(MAXME),ADIAB(MAXME),
     # VSOUND(MAXME),opac(maxme),cp(maxme)
c      WRITE(*,*) 'HELLO WORLD'
      
c    CHARACTER(*),PARAMETER :: fileplace = "FeHp006_alw/mu_test_data/"
      abund_in = 'diff_mu_logL=2_1231_abund'
      model_in = 'diff_mu_logL=2_1231_model'
C      head_in = 'diff_mu_logL=2_1231_header'
      file_out = 'diff_mu_logL=2_1231_grads'
C fileplace//      
c      OPEN(UNIT=21,FILE=head_in,status='OLD')
      OPEN(UNIT=28,FILE=abund_in,status='OLD')
      OPEN(UNIT=29,FILE=model_in,status='OLD')
      OPEN(UNIT=38,FILE=file_out,status='REPLACE')
      
C      read(21,180) MAXME,NMODEL
C       180  FORMAT(117X, 2I5)
C      WRITE(*,*) MAXME,',',NMODEL
      write(38,110)
      write(38,120)
 110  FORMAT(1H1,/,'      Mr            MU        d(mu)/dr_2    dln(mu)/   
     &dr_2 dln(mu)/dln(p)_2   d(X_Li7)/dr_2   d(X_C)/dr_2  d(X_N)/dr_2   
     &  d(X_He3)/dr_2  d(X_He4)/dr_2    Dthm      Dthm_0     MESH')

Cx1mas,MU,DMU_DR_2,DLNMU_DR_2,DEL_MU_2,
C     # DXLi7_DR_2,DXC_DR_2,DXN_DR_2,KM
 120  FORMAT('==========================================================
     #==================================================================
     #=====================================================')
c      IF (K .EQ. 1) THEN

c      END IF
      kcount = 0
C READING DO-LOOP
C 
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
      WRITE(*,*) x1mas(1)

C CALCULATION DO-LOOP
      DO K = 1,MAXME

       KM = K - 1
       KP = K + 1

C    CALCULATE FINVMU = INVERSE OF FINAL_MU = 1/MU = 1/MU_I + 1/MU_E
C SO: 1/MU = SUM_I ((Z_I + 1)*X_I/A_I), ASSUMING FULL IONIZATION
C N.B.: FINVMU DEFINED AT KM = K - 1
      FINVMU = (XI(km)*2/1)+(XHE3(km)*3/3)+(XHE4(km)*3/4)+(XC(km)*7/12)
     $ +(XN(km)*8/14)+(XO(km)*9/16)+(XD(km)*2/2)+(XLi7(km)*4/7)
     $ +(XBe7(km)*5/7)+(XC13(km)*7/13)+(XN15(km)*8/16)+(XO17(km)*9/17)
     $ +(XO18(km)*9/18)+(XFe56(km)*26/56)

C CALCULATE MEAN MOLECULAR WEIGHT, FINAL_MU, BY INVERTING FINVMU
      FINAL_MU = 1/FINVMU

C DEFINE FINVMU_K, THE EQUIVALENT OF FINVMU, BUT AT K = KM + 1
C errorbait
      FINVMU_K = (XI(k)*2/1)+(XHE3(k)*3/3)+(XHE4(k)*3/4)
     $ +(XC(k)*7/12)+(XN(k)*8/14)+(XO(k)*9/16)+(XD(k)*2/2)
     S +(XLi7(k)*4/7)+(XBe7(k)*5/7)+(XC13(k)*7/13)+(XN15(k)*8/16)
     $ +(XO17(k)*9/17)+(XO18(k)*9/18)+(XFe56(k)*26/56)

C IN THE SAME WAY, DEFINE FINAL_MU_K
      FINAL_MU_K = 1/FINVMU_K

C DEFINE FINVMU_P, THE EQUIVALENT OF FINVMU, BUT AT KP = K + 1
      FINVMU_P = (XI(kp)*2/1)+(XHE3(kp)*3/3)+(XHE4(kp)*3/4)
     $ +(XC(kp)*7/12)+(XN(kp)*8/14)+(XO(kp)*9/16)+(XD(kp)*2/2)
     S +(XLi7(kp)*4/7)+(XBe7(kp)*5/7)+(XC13(kp)*7/13)+(XN15(kp)*8/16)
     $ +(XO17(kp)*9/17)+(XO18(kp)*9/18)+(XFe56(kp)*26/56)

C IN THE SAME WAY, DEFINE FINAL_MU_P
      FINAL_MU_P = 1/FINVMU_P

C TAKE DIFFERENCE OF THE MU'S
C 1-MESH POINT WIDTH DERIVATIVE METHOD
C      DMU = FINAL_MU_K - FINAL_MU
C 1-MESH POINT WIDTH DERIVATIVE METHOD
C      DMU_P = FINAL_MU_P - FINAL_MU_K
C 2-POINT WIDTH METHOD
      DMU_2 = FINAL_MU_P - FINAL_MU

C Focus on individual isotopes
C 1-POINT METHOD
C      DXC = XC(K) - XC(KM)
C      DXLi7 = XLi7(K) - XLi7(KM)
C      DXN = XN(K) - XN(KM)
C 1-POINT METHOD
C      DXC_P = XC(KP) - XC(K)
C      DXLi7_P = XLi7(KP) - XLi7(K)
C      DXN_P = XN(KP) - XN(K)
C 2-POINT METHOD
      DXC_2 = XC(KP) - XC(KM)
      DXLi7_2 = XLi7(KP) - XLi7(KM)
      DXN_2 = XN(KP) - XN(KM)
      DXHE3_2 = XHE3(KP) - XHE3(KM)
      DXHE4_2 = XHE4(KP) - XHE4(KM)

C CALCULATE THE (LINEAR) MOLECULAR GRADIENT
C GET THE KTH-LAYER WIDTHS FOR DIFFERENT PARAMETERS
C 1-MESH POINT WIDTH DERIVATIVE METHOD -
C      DR = (2.21E+01)*(RFRAC(K)-RFRAC(KM))
C      DLNP = DLOG(10**logP(K)) - DLOG(10**logP(KM))
C      DLOGP = logP(K)-logP(KM)
C      DLNT = DLOG(10**logT(K)) - DLOG(10**logT(KM))
C      DLOGT = logT(K)-logT(KM)
C      DMASS = ()/emtot
C 1-MESH POINT WIDTH DERIVATIVE METHOD +
C      DR_P = (2.21E+01)*(RFRAC(KP)-RFRAC(K))
C      DLNP_P = DLOG(10**logP(KP)) - DLOG(10**logP(K))
C      DLOGP_P = logP(kp)-logP(k)
C      DLNT_P = DLOG(10**logT(KP)) - DLOG(10**logT(K))
C      DLOGT_P = logT(kp)-logT(k)
C      DMASS_P = ()/emtot
C 2-POINT WIDTH METHOD
      DR_2 = (2.21E+01)*(RFRAC(KP)-RFRAC(KM))
      DLNP_2 = DLOG(10**logP(KP)) - DLOG(10**logP(KM))
      DLOGP_2 = logP(kp)-logP(km)
      DLNT_2 = DLOG(10**logT(KP)) - DLOG(10**logT(KM))
      DLOGT_2 = logT(kp)-logT(km)
c      DMASS_2 = ()


C TAKE DIFFERENCE OF NATURAL LOG OF THE MU'S
C 1-POINT METHOD -
C      DLNMU = DLOG(FINAL_MU_K) - DLOG(FINAL_MU)
C      DLOGMU = DLOG10(FINAL_MU_K) - DLOG10(FINAL_MU)
C 1-POINT METHOD -
C      DLNMU_P = DLOG(FINAL_MU_P) - DLOG(FINAL_MU_K)
C      DLOGMU_P = DLOG10(FINAL_MU_P) - DLOG10(FINAL_MU_K)
C 2-POINT METHOD
      DLNMU_2 = DLOG(FINAL_MU_P) - DLOG(FINAL_MU)
      DLOGMU_2 = DLOG10(FINAL_MU_P) - DLOG10(FINAL_MU)

C DEFINE RELEVANT LINEAR GRADIENTS
C 1-POINT METHOD (K-1 to K)
C      DMU_DR = DMU/DR ! d(mu)/dr
C      DMU_DM = DMU/DMASS ! d(mu)/dm
C      DEL_MU = DLNMU/DLNP ! d(ln(mu))/d(ln(P))
C      DLNMU_DR = DLNMU/DR ! d(ln(mu))/dr
C      DXC_DR = DXC/DR
C      DXN_DR = DXN/DR
C      DXLi7_DR = DXLi7/DR

C DEFINE RELEVANT LINEAR GRADIENTS
C 1-POINT METHOD (K to K+1)
C      DMU_DR_P = DMU_P/DR_P ! d(mu)/dr
C      DMU_DM_P = DMU/DMASS ! d(mu)/dm
C      DEL_MU_P = DLNMU_P/DLNP_P ! d(ln(mu))/d(ln(P))
C      DLNMU_DR_P = DLNMU_P/DR_P ! d(ln(mu))/dr
C      DXC_DR_P = DXC_P/DR_P
C      DXN_DR_P = DXN_P/DR_P
C      DXLi7_DR_P = DXLi7_P/DR_P


C 2-POINT METHOD
      CALL AVOIDNAN(DMU_DR_2,DMU_2,DR_2) ! d(mu)/dr
C      CALL AVOIDNAN(DMU_DM_2,DMU_2,DMASS_2) ! d(mu)/dm
      CALL AVOIDNAN(DEL_MU_2,DLNMU_2,DLNP_2) ! d(ln(mu))/d(ln(P))
      CALL AVOIDNAN(DLNMU_DR_2,DLNMU_2,DR_2) ! d(ln(mu))/dr
      CALL AVOIDNAN(DXC_DR_2,DXC_2,DR_2) ! C12
      CALL AVOIDNAN(DXN_DR_2,DXN_2,DR_2) ! N14
      CALL AVOIDNAN(DXLi7_DR_2,DXLi7_2,DR_2) ! Li7
      CALL AVOIDNAN(DXHE3_DR_2,DXHE3_2,DR_2) ! He3
      CALL AVOIDNAN(DXHE4_DR_2,DXHE4_2,DR_2) ! He4

C THERMOHALINE MIXING IMPLEMENTATION
C LEDOUX CRITERION (CANTIELLO, LANGER, 2010):
C  0 >= DEL-DEL_AD >= (PHI/DELTA)*DEL_MU
C WITH PHI = PART(DLN(RHO)/DLN(MU))_P,T
c AND DELTA = PART(DLN(RHO)/DLN(T))_P,MU

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
C If loop comes from Lattanzio et al. (2015): D = 0 unless del_mu < 0
      if ((DEL_MU_2 .LT. 0) .AND. (DEL_MU_2 .GT. (RADIAT(K)-ADIAB(K))))
     # then
        Dthl = Cthl*Kthm*(DEL_MU_2)/(RADIAT(K)-ADIAB(K))
      else
        Dthl = 0.0
      end if

      if (DEL_MU_2 .LT. 0) then
        Dthl_0 = Cthl*Kthm*(DEL_MU_2)/(RADIAT(K)-ADIAB(K))
      else
        Dthl_0 = 0.0
      end if

      if (k .eq. MAXME) then
	write(*,*) 'km = (MAXME-1) = ', km
C       else
c	write(*,*) 'km = ',km
        
      end if
c       WRITE(*,*) 'WRITING K'
C       write(*,*) XXX( 1 , K )
      write(38,102) x1mas(K),FINAL_MU_K,DMU_DR_2,DLNMU_DR_2,DEL_MU_2,
     # DXLi7_DR_2,DXC_DR_2,DXN_DR_2,DXHE3_DR_2,DXHE4_DR_2,Dthl,Dthl_0,K
 100  FORMAT(1X,F12.9,14D14.7,105X,I5)
 101  FORMAT(F12.9,F10.7,F10.6,F8.5,F8.4,D13.5,3D11.3,D9.2,
     &D11.3,F6.3,D12.4,D12.4,D9.2,I5)
C 101  FORMAT(F12.9,F10.7,F8.4,F7.4,F8.4,1P,D13.5,1P,3D11.3,1P,D9.2,1P,
C     &D11.3,0P,F6.3,1P,D12.4,D12.4,D9.2,I5)
 102  FORMAT(1X,F12.9,9D15.7,2D12.4,I5)

c ,XI(kP),XI(kM),KP ,2D15.7,I5 ,DR_2,DLNP_2
      
      enddo
   

      WRITE(*,*) 'END OF WRITING'
      close(38)

      END

