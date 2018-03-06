      SUBROUTINE DIFFUS  !diffusione con bordo esterno all'inviluppo conv. 9-2-98
C     DIFFUSION WITH EXTERNAL BOUNDARY TO THE CONVECTIVE ENVELOPE
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'maincom.2p1'
      parameter (nele=13) ! 13 elemental species considered
      dimension em(lim),r(lim),p(lim),t(lim),ro(lim),app(nele,lim)
      dimension zatom(nele),atom(nele),y(nele,lim),zot(lim)
c     statements of species' atomic mass (atom) and number (zatom)
      data atom/1.d0,4.d0,12.d0,14.d0,16.d0,7.d0,13.d0,15.d0,17.d0,
     #21.d0,24.d0,27.d0,56.d0/
c     #18.d0,56.d0/
      data zatom/1.d0,2.d0,6.d0,7.d0,8.d0,3.d0,6.d0,7.d0,8.d0,
     #10.d0,12.d0,13.0D0,26.d0/
c     #8.d0,26.d0/
      solmas=fraz-fraz/5.d6
      do k=maxme-1,1,-1
       if(G(5,k)/emtot.lt.solmas)then
c       write(*,'(i6,3f15.8)')k,G(5,k)/emtot,fraz
c       pause
         if(G(6,k).lt.0.d0) go to 1
       endif
      end do
   1  maxo=k
c      write(*,*)'maxo, maxme:  ',maxo,maxme
C ************** CALCOLA RO **********
       do l=1,nele
        xxx(l,maxme+1)=xxx(l,maxme)
       enddo
       do l=1,5
         G(l,maxme+1)=G(l,maxme)
       enddo
      DO K=1,maxme
       ! takes simple average of kth and (k+1)th terms as the defintion
       ! for layer at mesh point k
       T(K)=1.0D+06*(G(4,K)+G(4,K+1))/2.D0
        P(K)=1.0D+17*(G(3,K)+G(3,K+1))/2.D0
        PS=P(K)
        TS=T(K)
c        CALL NSTATE(K,1,maxne,1,PS,TS,XXX,'EVOLUT  ',ROS,O1,O2,O3)
       CALL EOS(K,1,MAXNE,1,PS,TS,XXX,'DIFFUSIO',ROS,DAD,CSP,PMAU,G1,
     #DEL,EMUE,1)
        em(k)=1.0D33*(G(5,K)+G(5,K+1))/2.D0
        r(k)=1.0D10*(G(1,K)+G(1,K+1))/2.D0
        RO(K)=ROS
        y(1,k)=xxx(1,k)    ! Idrogeno
        y(2,k)=xxx(3,k)    ! Elio
        y(3,k)=xxx(4,k)    ! Carbonio
        y(4,k)=xxx(5,k)    ! azoto
        y(5,k)=xxx(6,k)    ! Ossigeno
       y(6,k)=xxx(8,k)    ! Li7
       y(7,k)=xxx(10,k)   ! C13
       y(8,k)=xxx(11,k)   ! N15  
       y(9,k)=xxx(12,k)   ! O17 
       y(10,k)=xxx(14,k)  ! Ne20
       y(11,k)=xxx(18,k)  ! Mg24  
       y(12,k)=xxx(21,k)  ! Al27  
       y(13,k)=xxx(26,k)  ! Ferro
        zot(k)=1.d0
        do j=1,nele
           if(y(j,k).le.1.d-30)y(j,k)=1.d-30
          zot(k)=zot(k)-y(j,k) ! avanzo al mesh k
        end do
C      if(k.le.15)then
C       write(34,'(i5,6e16.8)')k,T(K),P(K),r(k),RO(K),em(k)
C      endif
       do j=1,nele
        app(j,k)=y(j,k)
       enddo
      end do
c      write(*,*)'Fe prima'
c      write(*,'(7e18.10)')(y(4,i)/absol(26),i=1,maxme)
c      write(*,*)'Fe prima'
      call tbl(y,atom,zatom,em,r,p,t,ro,maxo)
c      write(*,*)'Fe dopo'
c      write(*,'(7e18.10)')(y(4,i)/absol(26),i=1,maxme)
c      write(*,*)'Fe dopo'
      do k=1,maxme
       do j=1,nele
         if(y(j,k).gt.app(j,k)*1.10d0.or
     #      .y(j,k).lt.app(j,k)*0.90d0)then
           do l=1,nele
             y(l,k)=app(l,k)
           enddo
         endif
       enddo
      end do



      do k=1,maxme
         sumz=0.d0
         do j=1,nele
           if(j.ne.2) then
             if(y(j,k).ge.0.95d0)then
c              do l=1,nele
               y(l,k)=app(l,k)
c              enddo
             endif
c  originale             if(y(j,k).le.1.d-10)y(j,k)=1.d-10
             if(y(j,k).le.1.d-10)y(j,k)=1.d-10
             sumz=sumz+y(j,k)
           end if
         end do
         y(2,k)=1.d0-sumz-zot(k)
         if(y(2,k).le.1.d-10) then
           write(*,*)'elio negativo: ',k,y(2,k)
           write(*,*)k,(y(l,k),l=1,nele)
           write(2,*)'elio negativo: ',k,y(2,k)
           y(2,k)=1.d-10
         end if
        xxx(1,k)=y(1,k)
        xxx(3,k)=y(2,k)
        xxx(4,k)=y(3,k)
        xxx(5,k)=y(4,k)
        xxx(6,k)=y(5,k)
        xxx(8,k)=y(6,k)
        xxx(10,k)=y(7,k)
        xxx(11,k)=y(8,k)
        xxx(12,k)=y(9,k)
        xxx(14,k)=y(10,k)
        xxx(18,k)=y(11,k)
        xxx(21,k)=y(12,k)
        xxx(26,k)=y(13,k)
      end do
      sumz=0.d0
      do j=4,maxne-1
         sumz=sumz+xxx(j,maxme-1)
      end do
      zsx=sumz/xxx(1,maxme-1)
      write(*,111)sumz,zsx
      write(2,111)sumz,zsx
      return
 111  format(1x,'<--- DIFFUSION --->  Z e Z/X : ',2f10.5)
      END
      
      subroutine tbl(x,a,z,em,ri,pri,ti,roi,maxo)
      implicit real*8(a-h,o-z)
      include 'maincom.2p1'
      include 'consts.2p1'
      parameter (nele=13,nr=12)
      common /difatm/emcoatm,icatm
      dimension x(nele,lim),u(nele,lim),dxo(nele)
      dimension f(nele,lim),gx(nele)
      dimension ri(lim),pri(lim),roi(lim),ti(lim),em(lim)
      dimension z(nele),a(nele),cl(nele+1,nele+1)
      dimension zserv(nele+1),aserv(nele+1),xserv(nele+1)
      dimension ap(nele+1),at(nele+1),ax(nele+1,nele+1)
      dimension cc(nele,lim),beta(nr+1),ydat(lim),xdat(lim)
      dimension dri(lim),xsi(nele,lim),dxco(nele),fatm(nele)
      dimension emtr(lim,nr+1),emtrt(nr+1,lim),alfa(nr+1,nr+1)
      dimension segno(13)
      data idebug/1/
      data segno/1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1./
c** IMPORTANTE !!! fisica e chimica definiti al punto intermedio **
c** definizione variabili adimensionali nel punto centrale ********
      e=ECharg
      if(idebug.eq.1) then
        open(7,file='fff.dat')
        write(7,*)'modello', nmd,'    massa, raggio, H, C, N, Fe'
      end if
c      write(*,'(7e18.10)')(x(1,i),i=1,6)
      m=nele+1
      dti=ht1/6.d13              ! passo in t in unita' interne
      min=0
      mxx=0
      do k=1,maxme
         ti(k)=ti(k)/1.d7
         roi(k)=roi(k)/1.d2
         ri(k)=ri(k)/SUNRcm
         if(ri(k).lt.0.5d0)then
          min=k
        endif
         em(k)=em(k)/SUNMg
        if(G(5,k)/emtot.lt.0.999d0)then
          mxx=k
        endif
      end do
c** definisce u e gli intervalli in raggio *******
C   define u and the ranges in radius
      do k=1,maxme
        if(k.eq.1) then
           dri(1)=ri(1)
        else
           dri(k)=ri(k)-ri(k-1)
        end if
        sumz=0.d0
        do j=1,nele
           sumz=sumz+x(j,k)*z(j)/a(j)
        end do
        do j=1,nele
           u(j,k)=x(j,k)*roi(k)*ri(k)**2
           cc(j,k)=x(j,k)/a(j)/sumz
        end do
      end do
      if(maxo.lt.(maxme-1)) then
         maximo=maxo          ! primo radiativo sotto inviluppo convettivo
      else
         maximo=maxme-2       ! se il bordo del convettivo e' oltre o a maxme-1
      end if
c******** calcolo f(k) *********
      do k=2,maxme-1
         tcgs=ti(k)*1.d7
         rocgs=roi(k)*1.d2
         sum=0.
         sumz=0.
         sumz2=0.
         do j=1,nele
            xserv(j)=x(j,k)
            aserv(j)=a(j)
            zserv(j)=z(j)
            sum=sum+xserv(j)/aserv(j)
            sumz=sumz+xserv(j)*zserv(j)/aserv(j)
            sumz2=sumz2+xserv(j)*(zserv(j)**2)/aserv(j)
         end do
         aserv(m)=5.446d-4
         zserv(m)=-1.d0
         sum=sum*rocgs*NAvo
         sumz=sumz*rocgs*NAvo
         sumz2=sumz2*rocgs*NAvo
         a0=(3.d0/(4.d0*PI*sum))**0.333333d0
         dl=dsqrt(kb*tcgs/(4.d0*PI*(e**2)*(sumz2+sumz)))
         amda=max(a0,dl)
         do i=1,m
           do j=1,m
              z1=zserv(i)
              if(z1.le.0.)z1=-z1
              z2=zserv(j)
              if(z2.le.0.)z2=-z2
              cl(i,j)=0.81245d0*dlog(1.d0+0.18769d0*
     %                   (4.*kb*tcgs*amda/(Z1*Z2*e**2))**1.2d0)
           end do
         end do
         call coedif(m,aserv,zserv,xserv,cl,ap,at,ax)
         gp=(dlog(pri(k+1))-dlog(pri(k-1)))/(dri(k)+dri(k+1))
         gt=(dlog(ti(k+1))-dlog(ti(k-1)))/(dri(k)+dri(k+1))
         do i=1,nele
           if(i.ne.2) then
             gx(i)=(dlog(cc(i,k+1))-dlog(cc(i,k-1)))/(dri(k)+dri(k+1))
           end if
         end do
         do j=1,nele
           if(j.ne.2) then
               gchim=0.d0
               do i=1,nele
                 if(i.ne.2) then
                    gchim=gchim+ax(j,i)*gx(i)
                 end if
               end do
               xsi(j,k)=at(j)*gt+ap(j)*gp+gchim
               f(j,k)=ti(k)**2.5d0*xsi(j,k)*u(j,k)/roi(k)
           end if
         end do
      end do
c ****************** smoot derivata prima ************************
      if(idebug.eq.1) then
         do k=1,maximo+1
           write(7,444)em(k),ri(k),f(1,k)
         end do
      end if
 444  format(1p,6E12.5)
 449  format(i5,1p,6E12.5)

      do j=1,nele
         if(j.ne.2) then
            f(j,1)=f(j,3)+(dri(1)+dri(2))/dri(2)*(f(j,2)-f(j,3))
             do k=1,maximo+1
                ctrl=f(j,k)/segno(j)
                if(ctrl.lt.0.d0)f(j,k)=0.d0    ! elimina spike
             end do
            do k=6,maximo-4
               ydat(k)=0.d0
               do k1=-5,5
                  ydat(k)=ydat(k)+f(j,k+k1)  ! moving average - window=5
               end do
               ydat(k)=ydat(k)/11.d0
               xdat(k)=em(k)
              
                if(f(j,k)/segno(j).le.0.d0)f(j,k)=ydat(k)
c         if (idebug.eq.1.and.j.eq.1)write(7,449)k,xdat(k-5),ydat(k-5)
            end do
            ntot=maximo-10
            call polino(ntot,nr,xdat,ydat,beta,emtr,emtrt,alfa)
            do k=1,min
               f(j,k)=0.d0
               do iex=1,nr
                  f(j,k)=f(j,k)+beta(iex+1)*em(k)**iex
               end do
            end do
         end if
      end do

      if(idebug.eq.1) then
         write(7,*)'bordi....  ',maxme, maximo
         do k=1,maximo+1
          write(7,444)f(1,k)
         end do
      end if
c ******* correzioni abbondanze **************
c **** 213=R(sun)**3*rho0*4*PI/M(sun)  *******
c      write(*,*)'calcolo correzioni: inizio   '

      if(mxx.lt.maximo)then
       do k=2,maxme
        do j=1,nele           ! zona radiativa
           if(j.ne.2) then
             if(k.le.mxx)then
             dx=dti/(em(k+1)-em(k-1))*(f(j,k+1)-f(j,k-1))*213.d0
              x(j,k)=x(j,k)-dx
            else
             dx=dti/(em(mxx+1)-em(mxx-1))*
     #               (f(j,mxx+1)-f(j,mxx-1))*213.d0
c             dx=dti/(emtot/1.989d0-em(mxx))*(f(j,mxx))*213.d0
              x(j,k)=x(j,k)-dx
            endif
           end if
        end do
       end do

      else
       do k=2,maximo
        do j=1,nele           ! zona radiativa
           if(j.ne.2) then
             dx=dti/(em(k+1)-em(k-1))*(f(j,k+1)-f(j,k-1))*213.d0
             x(j,k)=x(j,k)-dx
           end if
        end do
       end do

       if(icatm.eq.1) then
c       write(*,*)'emcoatm',emco,emcoatm
        emco=emtot/SUNMg33-emcoatm        ! massa env convettivo se in atm
       else
        emco=emtot/SUNMg33-em(maximo+1)   ! massa env convettivo se in henyey
c        write(*,*)'emcoatm convett',em(maximo+1),em(maxme),emtot,
c     #emtot/1.989d0
       if(emco.le.0.d0)then
         write(*,*)emco
         pause
       endif 
       end if
       do j=1,nele              ! mesh centrale (sviluppo in serie)
        if(j.ne.2) then
           x(j,1)=x(j,3)+(em(1)-em(3))/(em(2)-em(3))*(x(j,2)-x(j,3))
           dxco(j)=f(j,maximo+1)*213.d0*dti/emco
        end if
       end do
 10    format(1p,5e13.5)
       do k=maximo+1,maxme
         do j=1,nele           ! envelope convettivo
          if(j.ne.2) then
            x(j,k)=x(j,k)+dxco(j)
c             dx=dti/emco*f(j,k)*213.d0
c             if(k.le.maxme-5)then
c             x(j,k)=x(j,k)+dx
c            else
c             x(j,k)=x(j,k)+dxco(j)
c            endif
          end if
         end do
       end do
      endif
c      write(2,445)fatm(1)
c      write(2,446)fatm(5),fatm(6)
c      write(*,445)fatm(1)
c      write(*,446)fatm(5),fatm(6)
 445  format('fH,fC,fN: ',1p,3e13.6)
c 446  format('fO,fFe:   ',1p,2e13.6)
      return
      end

      subroutine polino(maxme,nr,xdat,ydat,bb,x,xt,aa)
c****** regressione polinomiale *********
      implicit real*8(a-h,o-z)
      dimension x(maxme,nr+1),xt(nr+1,maxme),aa(nr+1,nr+1),bb(nr+1,1)
      dimension xdat(maxme),ydat(maxme)
      nd=nr+1
c***** costruzione della matrice X *****
      do i=1,maxme
         x(i,1)=1.d0 ! set order 0 polynomial term in 1st column
         do J=2,nd ! for order > 0, use x(dat) value as root of term
            x(i,j)=xdat(i)**(j-1)
         end do
      end do
c***** costruzione della matrice X'***** = transpose of X
      do i=1,nd
        do j=1,maxme
          xt(i,j)=x(j,i)
        end do
      end do
c***** calcolo matrice AA=X'*X ********
      call xmatr(nd,nd,maxme,xt,x,aa)
c***** calcolo matrice BB=X'*Y ********
      call xmatr(nd,1,maxme,xt,ydat,bb)
c***** risoluzione equazione AA*C=BB **
      call kernel (aa,bb,nd,ks)
      if(ks.eq.1) then
         write(*,*)'qualcosa non va nella kernel'
         stop
      end if
c      write(2,*)'coefficienti'
c      do i=1,nd
c        write(2,*)bb(i,1)
c      end do
      return
      end

      subroutine xmatr(nr,nc,maxme,x,y,aa)
c*********** prodotto matriciale  X*Y=AA ****************
c*  nr    = numero righe matrice X                      *
c*  nc    = numero colonne matrice Y                    *
c*  maxme = numero colonne matrice X o righe matrice Y  *
c********************************************************
      implicit real*8(a-h,o-z)
      dimension x(nr,maxme),y(maxme,nc),aa(nr,nc)
      do i=1,nr
         do j=1,nc
            aa(i,j)=0.d0
         end do
      end do
      do i=1,nr
         do k=1,nc
            do j=1,maxme
               aa(i,k)=aa(i,k)+x(i,j)*y(j,k)
            end do
         end do
      end do
      return
      end
C*************************************************************
C This routine was written by Anne A. Thoul, at the Institute
C for Advanced Study, Princeton, NJ 08540.
C See Thoul et al., Ap.J. 421, p. 828 (1994)
C The subroutines LUBKSB and LUDCMP are from Numerical Recipes.
C*************************************************************
C This routine inverses the burgers equations.
C
C The system contains N equations with N unknowns.
C The equations are: the M momentum equations,
C                    the M energy equations,
C                    two constraints: the current neutrality
C                                     the zero fluid velocity.
C The unknowns are: the M diffusion velocities,
C                   the M heat fluxes,
C                   the electric field E
C                   the gravitational force g.
C
C**************************************************
      SUBROUTINE coedif(M,A,Z,X,CL,AP,AT,AX)

C The parameter M is the number of species considered.
C
C Fluid 1 is the hydrogen
C Fluid 2 is the helium
C Fluids 3 to M-1 are the heavy elements
C Fluid M is the electrons
C
C The vectors A,Z and X contain the atomic mass numbers,
C the charges (ionization), and the mass fractions, of the elements.
C NOTE: Since M is the electron fluid, its mass and charge must be
C      A(M)=m_e/m_u
C      Z(M)=-1.
C
C The array CL contains the values of the Coulomb Logarithms.
C The vector AP, AT, and array AX contains the results for the diffusion
C coefficients.

      IMPLICIT NONE

      INTEGER M,N,I,J,L,MMAX,NMAX
      PARAMETER (MMAX=20,NMAX=42)
      INTEGER INDX(NMAX)
      REAL*8 A(M),Z(M),X(M),AP(M),AT(M),AX(M,M),CL(M,M)
      REAL*8 C(MMAX),CC,AC,XX(MMAX,MMAX),Y(MMAX,MMAX),YY(MMAX,MMAX),
     $     K(MMAX,MMAX)
      REAL*8 ALPHA(NMAX),NU(NMAX),GAMMA(NMAX,NMAX),DELTA(NMAX,NMAX),
     $     GA(NMAX)
      REAL*8 TEMP,KO,D

C The vector C contains the concentrations
C CC is the total concentration: CC=sum(C_s)
C AC is proportional to the mass density: AC=sum(A_s C_s)
C The arrays XX,Y,YY and K are various parameters which appear in
C Burgers equations.
C The vectors and arrays ALPHA, NU, GAMMA, DELTA, and GA represent
C the "right- and left-hand-sides" of Burgers equations, and later
C the diffusion coefficients.


C Initialize parameters:

      KO=2.
      N=2*M+2
      DO I=1,M
         C(I)=0.
      ENDDO
      CC=0.
      AC=0.

C Calculate concentrations from mass fractions:

      TEMP=0.
      DO I=1,M-1
         TEMP=TEMP+Z(I)*X(I)/A(I)
      ENDDO
      DO I=1,M-1
         C(I)=X(I)/A(I)/TEMP
      ENDDO
      C(M)=1.

C Calculate CC and AC:

      DO I=1,M
         CC=CC+C(I)
         AC=AC+A(I)*C(I)
      ENDDO

C Calculate the mass fraction of electrons:

      X(M)=A(M)/AC

C Calculate the coefficients of the burgers equations

      DO I=1,M
         DO J=1,M
            XX(I,J)=A(J)/(A(I)+A(J))
            Y(I,J)=A(I)/(A(I)+A(J))
            YY(I,J)=3.0*Y(I,J)+1.3*XX(I,J)*A(J)/A(I)
            K(I,J)=1.*CL(I,J)*
     $           dSQRT(A(I)*A(J)/(A(I)+A(J)))*C(I)*C(J)*
     $           Z(I)**2*Z(J)**2
         ENDDO
      ENDDO

C Write the burgers equations and the two constraints as
C alpha_s dp + nu_s dT + sum_t(not 2 or M) gamma_st dC_t
C                     = sum_t delta_st w_t

      DO I=1,M
         ALPHA(I)=C(I)/CC
         NU(I)=0.
         DO J=1,M
            GAMMA(I,J)=0.
         ENDDO
         DO J=1,M
            IF ((J.NE.2).AND.(J.NE.M)) THEN
               GAMMA(I,J)=-C(J)/CC+C(2)/CC*Z(J)*C(J)/Z(2)/C(2)
               IF (J.EQ.I) THEN
                  GAMMA(I,J)=GAMMA(I,J)+1.
               ENDIF
               IF (I.EQ.2) THEN
                  GAMMA(I,J)=GAMMA(I,J)-Z(J)*C(J)/Z(2)/C(2)
               ENDIF
               GAMMA(I,J)=GAMMA(I,J)*C(I)/CC
            ENDIF
         ENDDO

         DO J=M+1,N
            GAMMA(I,J)=0.
         ENDDO
      ENDDO

      DO I=M+1,N-2
         ALPHA(I)=0.
         NU(I)=2.5*C(I-M)/CC
         DO J=1,N
            GAMMA(I,J)=0.
         ENDDO
      ENDDO

      ALPHA(N-1)=0.
      NU(N-1)=0.
      DO J=1,N
         GAMMA(N-1,J)=0.
      ENDDO

      ALPHA(N)=0.
      NU(N)=0.
      DO J=1,N
         GAMMA(N,J)=0.
      ENDDO

      DO I=1,N
         DO J=1,N
            DELTA(I,J)=0.
         ENDDO
      ENDDO

      DO I=1,M
         DO J=1,M
            IF (J.EQ.I) THEN
               DO L=1,M
                  IF(L.NE.I) THEN
                     DELTA(I,J)=DELTA(I,J)-K(I,L)
                  ENDIF
               ENDDO
            ELSE
               DELTA(I,J)=K(I,J)
            ENDIF
         ENDDO

         DO J=M+1,N-2
            IF(J-M.EQ.I) THEN
               DO L=1,M
                  IF (L.NE.I) THEN
                     DELTA(I,J)=DELTA(I,J)+0.6*XX(I,L)*K(I,L)
                  ENDIF
               ENDDO
            ELSE
               DELTA(I,J)=-0.6*Y(I,J-M)*K(I,J-M)
            ENDIF
         ENDDO

         DELTA(I,N-1)=C(I)*Z(I)

         DELTA(I,N)=-C(I)*A(I)
      ENDDO

      DO I=M+1,N-2
         DO J=1,M
            IF (J.EQ.I-M) THEN
               DO L=1,M
                  IF (L.NE.I-M) THEN
                     DELTA(I,J)=DELTA(I,J)+1.5*XX(I-M,L)*K(I-M,L)
                  ENDIF
               ENDDO
            ELSE
               DELTA(I,J)=-1.5*XX(I-M,J)*K(I-M,J)
            ENDIF
         ENDDO

         DO J=M+1,N-2
            IF (J-M.EQ.I-M) THEN
               DO L=1,M
                  IF (L.NE.I-M) THEN
                     DELTA(I,J)=DELTA(I,J)-Y(I-M,L)*K(I-M,L)*
     $                    (1.6*XX(I-M,L)+YY(I-M,L))
                  ENDIF
               ENDDO
               DELTA(I,J)=DELTA(I,J)-0.8*K(I-M,I-M)
            ELSE
               DELTA(I,J)=2.7*K(I-M,J-M)*XX(I-M,J-M)*Y(I-M,J-M)
            ENDIF
         ENDDO

         DELTA(I,N-1)=0.

         DELTA(I,N)=0.
      ENDDO

      DO J=1,M
         DELTA(N-1,J)=C(J)*Z(J)
      ENDDO
      DO J=M+1,N
         DELTA(N-1,J)=0.
      ENDDO

      DO J=1,M
         DELTA(N,J)=C(J)*A(J)
      ENDDO
      DO J=M+1,N
         DELTA(N,J)=0.
      ENDDO




C Inverse the system for each possible right-hand-side, i.e.,
C if alpha is the r.h.s., we obtain the coefficient A_p
C if nu    ---------------------------------------- A_T
C if gamma(i,j) ----------------------------------- A_Cj
C
C If I=1, we obtain the hydrogen diffusion velocity
C If I=2, ------------- helium   ------------------
C If I=3,M-1, --------- heavy element -------------
C If I=M, ------------- electrons -----------------
C For I=M,2M, we get the heat fluxes
C For I=N-1, we get the electric field
C For I=N, we get the gravitational force g

      CALL LUDCMP(DELTA,N,NMAX,INDX,D)

      CALL LUBKSB(DELTA,N,NMAX,INDX,ALPHA)
      CALL LUBKSB(DELTA,N,NMAX,INDX,NU)
      DO J=1,N
         DO I=1,N
            GA(I)=GAMMA(I,J)
         ENDDO
         CALL LUBKSB(DELTA,N,NMAX,INDX,GA)
         DO I=1,N
            GAMMA(I,J)=GA(I)
         ENDDO
      ENDDO

C The results for the coefficients must be multiplied by p/K_0:

      DO I=1,M
         ALPHA(I)=ALPHA(I)*KO*AC*CC
         NU(I)=NU(I)*KO*AC*CC
         DO J=1,M
            GAMMA(I,J)=GAMMA(I,J)*KO*AC*CC
         ENDDO
      ENDDO

      DO I=1,M
         AP(I)=ALPHA(I)
         AT(I)=NU(I)
         DO J=1,M
            AX(I,J)=GAMMA(I,J)
         ENDDO
      ENDDO


      RETURN

      END

*********************************************************************
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER N,NP
C     ..
C     .. Array Arguments ..
      REAL*8 A(NP,NP),B(N)
      INTEGER INDX(N)
C     ..
C     .. Local Scalars ..
      REAL*8 SUM
      INTEGER I,II,J,LL
C     ..
      II = 0
      DO 12 I = 1,N
          LL = INDX(I)
          SUM = B(LL)
          B(LL) = B(I)
          IF (II.NE.0) THEN
              DO 11 J = II,I - 1
                  SUM = SUM - A(I,J)*B(J)
   11         CONTINUE

          ELSE IF (SUM.NE.0.) THEN
              II = I
          END IF

          B(I) = SUM
   12 CONTINUE
      DO 14 I = N,1,-1
          SUM = B(I)
          IF (I.LT.N) THEN
              DO 13 J = I + 1,N
                  SUM = SUM - A(I,J)*B(J)
   13         CONTINUE
          END IF

          B(I) = SUM/A(I,I)
   14 CONTINUE
      RETURN

      END

*********************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

      IMPLICIT NONE

C     .. Parameters ..
      INTEGER NMAX
      REAL*8 TINY
      PARAMETER (NMAX=100,TINY=1.0D-20)
C     ..
C     .. Scalar Arguments ..
      REAL*8 D
      INTEGER N,NP
C     ..
C     .. Array Arguments ..
      REAL*8 A(NP,NP)
      INTEGER INDX(N)
C     ..
C     .. Local Scalars ..
      REAL*8 AAMAX,DUM,SUM
      INTEGER I,IMAX,J,K
C     ..
C     .. Local Arrays ..
      REAL*8 VV(NMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      D = 1.
      DO 12 I = 1,N
          AAMAX = 0.
          DO 11 J = 1,N
              IF (ABS(A(I,J)).GT.AAMAX) AAMAX = ABS(A(I,J))
   11     CONTINUE
          IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
          VV(I) = 1./AAMAX
   12 CONTINUE
      DO 19 J = 1,N
          IF (J.GT.1) THEN
              DO 14 I = 1,J - 1
                  SUM = A(I,J)
                  IF (I.GT.1) THEN
                      DO 13 K = 1,I - 1
                          SUM = SUM - A(I,K)*A(K,J)
   13                 CONTINUE
                      A(I,J) = SUM
                  END IF

   14         CONTINUE
          END IF

          AAMAX = 0.
          DO 16 I = J,N
              SUM = A(I,J)
              IF (J.GT.1) THEN
                  DO 15 K = 1,J - 1
                      SUM = SUM - A(I,K)*A(K,J)
   15             CONTINUE
                  A(I,J) = SUM
              END IF

              DUM = VV(I)*ABS(SUM)
              IF (DUM.GE.AAMAX) THEN
                  IMAX = I
                  AAMAX = DUM
              END IF

   16     CONTINUE
          IF (J.NE.IMAX) THEN
              DO 17 K = 1,N
                  DUM = A(IMAX,K)
                  A(IMAX,K) = A(J,K)
                  A(J,K) = DUM
   17         CONTINUE
              D = -D
              VV(IMAX) = VV(J)
          END IF

          INDX(J) = IMAX
          IF (J.NE.N) THEN
              IF (A(J,J).EQ.0.) A(J,J) = TINY
              DUM = 1./A(J,J)
              DO 18 I = J + 1,N
                  A(I,J) = A(I,J)*DUM
   18         CONTINUE
          END IF

   19 CONTINUE
      IF (A(N,N).EQ.0.) A(N,N) = TINY
      RETURN
      END

C      SUBROUTINE THERMOHALINE(MU)
C      INCLUDE 'maincom.2p1'
C      DIMENSION MU(MAXME)
C      DO K = 1, MAXME-1
C         DR = G(1,K+1) - G(1,K)
C         LOGMU(K) = DLOG(MU(K))
C         DLOGMU = LOGMU(K+1) - LOGMU(K)
C         GRADIENT = DLOGMU/DR
C      END DO
C      
C      END
