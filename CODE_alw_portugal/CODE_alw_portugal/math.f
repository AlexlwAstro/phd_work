C*****************************************************************************
C  QUESTA ROUTINE CONTIENE: KERNEL - NCUB - PARAB - POLINT -SPLINK
C*****************************************************************************
C
      SUBROUTINE KERNEL (A,B,N,KS) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION B(1),A(1) 
      TOL=0.D0 
      KS=0 
      JJ=-N 
      DO 65 J=1,N 
      JY=J+1 
      JJ=JJ+N+1 
      BIGA=0.D0 
      IT=JJ-J 
      DO 30 I=J,N 
      IJ=IT+I 
      IF(DABS(BIGA)-DABS(A(IJ))) 20,30,30 
   20 BIGA=A(IJ) 
      IMAX=I 
   30 CONTINUE 
      IF(DABS(BIGA)-TOL) 35,35,40 
   35 KS=1 
      RETURN 
   40 I1=J+N*(J-2) 
      IT=IMAX-J 
      DO 50 K=J,N 
      I1=I1+N 
      I2=I1+IT 
      SAVE=A(I1) 
      A(I1)=A(I2) 
      A(I2)=SAVE 
   50 A(I1)=A(I1)/BIGA 
      SAVE=B(IMAX) 
      B(IMAX)=B(J) 
      B(J)=SAVE/BIGA 
      IF(J-N) 55,70,55 
   55 IQS=N*(J-1) 
      DO 65 IX=JY,N 
      IXJ=IQS+IX 
      IT=J-IX 
      DO 60 JX=JY,N 
      IXJX=N*(JX-1)+IX 
      JJX=IXJX+IT 
   60 A(IXJX)=A(IXJX)-(A(IXJ)*A(JJX)) 
   65 B(IX)=B(IX)-(B(J)*A(IXJ)) 
   70 NY=N-1 
      IT=N*N 
      DO 80 J=1,NY 
      IA=IT-J 
      IB=N-J 
      IC=N 
      DO 80 K=1,J 
      B(IB)=B(IB)-A(IA)*B(IC) 
      IA=IA-N 
   80 IC=IC-1 
      RETURN 
      END 

      SUBROUTINE NCUB 
      IMPLICIT REAL*8 (A-H,O-Z) 
      COMMON/CUBICA/A1,A2,A3,A4,B1,B2,B3,B4,ALFA,BETA,GAMMA,DELTA 
      D1=B2*B2+B2*B1 
      C1=B3*B3+B3*B1-D1 
      C2=B3-B2 
      D2=(A2-A1)/(B2-B1) 
      C3=(A3-A1)/(B3-B1)-D2 
      C4=B4*B4+B4*B1-D1 
      C5=B4-B2 
      C6=(A4-A1)/(B4-B1)-D2 
      D3=C1*C5-C2*C4 
      ALFA=(C3*C5-C2*C6)/D3 
      BETA=(C1*C6-C4*C3)/D3 
      GAMMA=D2-ALFA*(D1+B1*B1)-BETA*(B2+B1) 
      DELTA=A1-ALFA*B1*B1*B1-BETA*B1*B1-GAMMA*B1 
      RETURN 
      END 

      SUBROUTINE PARAB(X1,X2,X3,Y1,Y2,Y3,AA,BB,CC) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DX1=X2-X1 
      DX2=X3-X1 
      DQ1=X2**2-X1**2 
      DQ2=X3**2-X1**2 
      DY1=Y2-Y1 
      DY2=Y3-Y1 
      RAP=DX1/DX2 
      AA=(DY1-RAP*DY2)/(DQ1-RAP*DQ2) 
      BB=(DY1-AA*DQ1)/DX1 
      CC=Y1-AA*X1**2-BB*X1 
      RETURN 
      END 

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
c      write(*,*) 'nn  ',n
c      pause
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
          write(*,*) 'son in polint'
          write(*,*) x,(xa(k),k=1,4)
          pause 'failure in polint'
          endif 
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      SUBROUTINE SPLINK(N,NC,M,MC,IOP)
      PARAMETER(MX=10)
      COMMON/SP/X(MX),Y(MX),DERIV(MX,2),ZIN(MX),FVALUE(MX),
     $FDERIV(MX,2)
      COMMON/SPAPPR/SECD1,SECDN,VOFINT,IERR,NXY
      DATA ZERO,HALF,ONE,THREE/0.,.5,1.,3./
      DATA THIRD,SIXTH/0.33333333,0.16666667/
C      DO IM=1,N
C            WRITE(*,*)Y(IM),X(IM),ZIN(1)
C      ENDDO
C      PAUSE
 1000 IF (IOP.GT.0) GO TO 1110
      IERR=0
      IF (N.GE.4) GO TO 1010
      IERR=1
      GO TO 2000
 1010 IF (IOP.NE.-1) GO TO 1015
      SECD1=ZERO
      SECDN = ZERO
      BET1=ONE/(ONE+HALF*(X(2)-X(1))/(X(3)-X(2)))
      ALF1=BET1*(ONE- ((X(2)-X(1))/(X(3)-X(2)))**2)
      BETN=ONE/(ONE+HALF*(X(N)-X(N-1))/(X(N-1)-X(N-2)))
      ALFN=BETN*(ONE- ((X(N)-X(N-1))/(X(N-1)-X(N-2)))**2)
 1015 DERIV(1,2)=SECD1
      DERIV(N,2)=SECDN
      DERIV(1,1)=ZERO
      DXPLUS=X(2)-X(1)
      IF ( DXPLUS.GT.ZERO) GO TO 1020
      IN=1
      IERR=2
      GO TO 2000
 1020 DYPLUS=(Y(2)-Y(1))/DXPLUS
      IU=N-1
      DO 1040 I=2,IU
      DXMIN =DXPLUS
      DYMIN =DYPLUS
      DXPLUS=X(I+1)-X(I)
      IF (DXPLUS.GT.ZERO) GO TO 1030
      IN=I
      IERR=2
      GO TO 2000
 1030 DXINV =ONE/(DXPLUS+DXMIN)
      DYPLUS=(Y(I+1)-Y(I))/DXPLUS
      DIVDIF=DXINV*(DYPLUS-DYMIN)
      ALF   =HALF*DXINV*DXMIN
      BET   =HALF-ALF
      IF (I.EQ.2)  DIVDIF=DIVDIF-THIRD*ALF*DERIV(1,2)
      IF (I.EQ.IU) DIVDIF=DIVDIF-THIRD*BET*DERIV(N,2)
      IF (I.EQ.2) ALF=ZERO
      IF (IOP.NE.-1) GO TO 1035
      IF (I.NE.2) GO TO 1032
      BET=BET*ALF1
      DIVDIF=DIVDIF*BET1
      GO TO 1035
 1032 IF (I.NE.IU) GO TO 1035
      ALF=ALF*ALFN
      DIVDIF=DIVDIF*BETN
 1035 DXINV =ONE/(ONE+ALF*DERIV(I-1,1))
      DERIV(I,1)=-DXINV*BET
      DERIV(I,2)= DXINV*(THREE*DIVDIF-ALF*DERIV(I-1,2))
 1040 CONTINUE
 1050 DO 1060 I=2,IU
      J=N-I
      DERIV(J,2)=DERIV(J,1)*DERIV(J+1,2)+DERIV(J,2)
 1060 CONTINUE
      IF (IOP.NE.-1) GO TO 1070
      DERIV(1,2)=((X(3)-X(1))/(X(3)-X(2)))*DERIV(2,2)-((X(2)-X(1))/(X(3)
     1-X(2)))*DERIV(3,2)
      DERIV(N,2)=-((X(N)-X(N-1))/(X(N-1)-X(N-2)))*DERIV(N-2,2)+((X(N)-X(
     1N-2))/(X(N-1)-X(N-2)))*DERIV(N-1,2)
 1070 VOFINT=ZERO
      DO 1080 I=1,IU
      DXPLUS=X(I+1)-X(I)
      DYPLUS=Y(I+1)-Y(I)
      DIVDIF=DYPLUS/DXPLUS
      DERIV(I,1)=DIVDIF-DXPLUS*(THIRD*DERIV(I,2)+SIXTH*DERIV(I+1,2))
      DXPLUS=HALF*DXPLUS
      VOFINT=VOFINT+DXPLUS*(Y(I+1)+Y(I)-THIRD*(DERIV(I+1,2)+DERIV(I,2))*
     1DXPLUS**2)
 1080 CONTINUE
      DXPLUS=X(N)-X(N-1)
      DYPLUS=Y(N)-Y(N-1)
      DIVDIF=DYPLUS/DXPLUS
      DERIV(N,1)=DIVDIF+DXPLUS*(SIXTH*DERIV(N-1,2)+THIRD*DERIV(N,2))
      NXY=N
 1110 IF (M.LT.1) RETURN
      XL=X(1)
      XU=X(2)
      IP=3
      IL=0
 1120 DO 1160 J=1,M
      ARG=ZIN(J)
      IF (ARG.GT.XU) THEN
      IF(IP.GT.NXY) GO TO 1185
      IPP=IP
      DO 1180 I=IPP,NXY
      IF (ARG.GT.X(I)) GO TO 1180
      XL=X(I-1)
      XU=X(I)
      IP=I+1
      IL=0
      GO TO 1140
 1180 CONTINUE
 1185 IERR=3
      IP=NXY+1
      GO TO 2010
      ENDIF
      IF (ARG.LT.XL) THEN
      IPP=IP
      DO 1200 I=1,IPP
      II=IP-I-2
      IF (II.EQ.0) GO TO 1210
      IF (ARG.LT.X(II)) GO TO 1200
      XL=X(II)
      XU=X(II+1)
      IP=II+2
      IL=0
      GO TO 1140
 1200 CONTINUE
 1210 IERR=4
      IP=3
      GO TO 2010
      ENDIF
 1130 IF (IL.GT.0) GO TO 1150
 1140 II=IP-2
      A0=Y(II)
      A1=DERIV(II,1)
      A4=DERIV(II,2)
      A6=(DERIV(II+1,2)-A4)/(XU-XL)
      A2=HALF*A4
      A3=SIXTH*A6
      A5=HALF*A6
      IL=1
 1150 ARG=ARG-XL
      FVALUE(J)=((A3*ARG+A2)*ARG+A1)*ARG+A0
      FDERIV(J,1)=(A5*ARG+A4)*ARG+A1
      FDERIV(J,2)=A6*ARG+A4
 1155 CONTINUE
 1160 CONTINUE
      RETURN
 2000 CONTINUE
      IF(IERR.EQ.1) WRITE(*,3001)
      IF(IERR.EQ.2) WRITE(*,3002)IN,X(IN),X(IN+1),Y(IN),Y(IN+1)
      RETURN
 2010 CONTINUE
      IF(IERR.EQ.3) WRITE(*,3003)ARG
      IF(IERR.EQ.4) WRITE(*,3004)ARG
      FVALUE(J)=ZERO
      FDERIV(J,1)=ZERO
      FDERIV(J,2)=ZERO
      II=IP-2
      XL=X(II)
      XU=X(II+1)
      IL=0
      RETURN
 3001 FORMAT('*** NUMBER OF POINTS LESS THEN 4 ***')
 3002 FORMAT('*** X VALUES NOT INCREASING ***',I5,1P,4E11.4)
 3003 FORMAT('*** Z VALUE ABOVE RANGE ***',1P,E15.8)
 3004 FORMAT('*** Z VALUE BELOW RANGE ***',1P,E15.8)
      END
                
