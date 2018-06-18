      program integra
c****************** declarations ******************
      dimension ff(21,19000), fs(2000), wavef(21,19000), waves(2000)
      dimension ff2(21,19000), fvega(2000), bcV(21)
      dimension wavemicr(2000), ratio(2000), yw(2000)
      dimension ndat(21)
      character namein(49)*100
c
c************** units *******************************
c      open(unit=1, file='WFPC2widefilters.dat')
c      open(unit=2, file='INPUT')
      open(unit=2, file='fm10ak2odfnew.pck')
      open(unit=3, file='check')
      open(unit=4, file='spettroVega')
      open(unit=9, file='OUTPUT')
      open (unit=11, file='INPUTfilters', status='old' )
c****************** insert extinction A(V)*************
      write(*,*) 'Insert A(V) and Rv=A(V)/E(B-V) (standard value 3.1)'
      read(*,*) avl, rv
c***************************************************
c*********** read number and names filters 
       
       read (11,*) nfilt ! number filters

       do 35 irt=1,nfilt ! number filters
         read(11,'(a)') namein(irt)       
c         write(*,*) irt, 'leggo nome',  namein(irt)
 35    continue

       close(11)
c************** read filter response curves*******************

      do 2222 ik=1,nfilt

      open(unit=1, file=namein(ik), status='old' ) 

      do 1 i=1, 5
        read(1,*)
 1    continue  

       ndat(ik)=0

      do 2 i=1, 19000
c       read(1,*, end=21) wavef(ik, i),ff(ik,i)
        read(1,*, end=21) indp, wavef(ik, i),ff(ik,i)
c        write(*,*) 'arrivo qui', indp
        wavef(ik,i)=wavef(ik,i)/10.e0        ! lambda in nm
c       write(3,*) i, wavef(ik,i), ff(ik,i)
 2    continue       

 21   ndat(ik)=i-1

      write(*,*) ik, namein(ik), ndat(ik)
 2222 continue
      close(1)

c********** read Vega spectrum **************************
      do 13 i=1, 2
        read(4,*)
 13    continue      
      do 14 i=1, 1221
        read(4,101) waves(i), fvega(i)    ! lambda in nm
c        write(3,*) i, waves(i), fvega(i)
c******** convert to Flambda for integration
        fvega(i)=4.0e0*2.99792458e17*fvega(i)/(waves(i)**2.e0)
c********* dilution factor ***********
        fvega(i)=fvega(i)*6.247e-17
c*********************************** 
 14   continue    

 101  format(9x,f9.2,20x,e13.4)
      close(4)
c103  format (4x,f8.0,9x,f8.5)
 103  format (4x,f8.0,9x,f8.5)
 104  FORMAT(8e10.4)
 105  FORMAT(5e10.4)
c
c*******
c******interpolo S(filtro) sugli stessi lambda dello spettro
c*******
      do 2223 ilj=1,nfilt

      do 5 i=1, 1221
        do 6 j=1, ndat(ilj)
          if(waves(i).lt.wavef(ilj,1)) then 
             ff2(ilj,i)=0.0e0
             write(3,*) waves(i), ff2(ilj,i)
          goto 5
          end if
          if(waves(i).gt.wavef(ilj,ndat(ilj))) then 
             ff2(ilj,i)=0.0e0
             write(3,*) waves(i), ff2(ilj,i)
          goto 5
          end if
          if(waves(i).eq.wavef(ilj,j)) then
           ff2(ilj,i)=ff(ilj,j)     
           write(3,*) waves(i), ff2(ilj,i)
           goto 5
          end if
        if(waves(i).lt.wavef(ilj,j).and.j.gt.1.and.j.lt.ndat(ilj)) then
           coe=(ff(ilj,j)-ff(ilj,j-1))/(wavef(ilj,j)-wavef(ilj,j-1))
           ff2(ilj,i)=coe*(waves(i)-wavef(ilj,j-1))+ff(ilj,j-1)     
           write(3,*) waves(i), ff2(ilj,i)
           goto 5
          end if
 6      continue
 5    continue  

 2223 continue
c      
      close(3)
c
c
c
c********** read input spectra and start loop ***********
c
c waves(i)    wavelength in nm
c Flambda=4*fs*c/wavelength^2  c=speed of light 
c
c
c      rewind(2)
c
c
c
c     do 699 i=1,153
c      write(*,*) 'leggo'
c      read(2,*) 
c 699  continue
c
      do 698 i=1,175
      read(2,*) 
 698  continue

  
c************ loop over spectra ******************
c
c
      do 700 iz=1,500
      write(*,*) 'arrivo qui'

c******* read Teff and log(g)

      read(2,103) teff, zlogg      
      write(*,*) teff, zlogg
 
c********** read fluxes *****************
 
c       write(*,*) 'arrivo qui', iz
c       write(*,*), iz
      
       do 701 j=1,152

c      write(*,*) 'leggo'

c       read(2,*)
 
       read (2,104) (fs(8*(j-1)+il), il=1,8)

c       write(*,*) (fs(8*(j-1)+il), il=1,8)
c       pause

 701  continue


      read(2,105) (fs(1216+il), il=1,5)

c      write(*,*) 'arrivo qui'
c********** skip continuum fluxes ************
      do 702 ij=1,153
       read(2,*)
 702  continue 

c******** convert to Flambda for integration
      do 703 ip=1,1221 
       fs(ip)=4.0e0*2.99792458e17*fs(ip)/(waves(ip)**2.e0)
 703  continue 
c***********************************************
c
c******************* Extinction law*************
c  Cardelli et al. (1989) with Rv=3.1 ***********
c   
      do 3211 i=1, 1221
      wavemicr(i)=1.e3/waves(i) ! (1/lambda in micron)
c
      if (wavemicr(i).lt.0.3e0) ratio(i)=0.0e0 
c
      if(wavemicr(i).ge.0.3e0.and.wavemicr(i).le.1.1e0) then 
      ratio(i)=0.574e0*((wavemicr(i))**1.61e0)-(0.527e0/3.1e0)*
     #((wavemicr(i))**1.61e0)   
      end if
c
      if(wavemicr(i).ge.1.1e0.and.wavemicr(i).le.3.3e0) then
      yw(i)=wavemicr(i)-1.82e0      
      app=1.e0+0.17699e0*yw(i)-0.50447e0*yw(i)**2.e0-0.02427e0*yw(i)**
     #3.e0+0.72085e0*yw(i)**4.e0+0.01979e0*yw(i)**5.e0
     #-0.7753e0*yw(i)**6.e0+0.32999e0*yw(i)**7.e0
      bpp=1.41338e0*yw(i)+2.28305e0*yw(i)**2.e0+1.07233e0*yw(i)**
     #3.e0-5.38434e0*yw(i)**4.e0-0.62251e0*yw(i)**5.e0
     #+5.3026*yw(i)**6.e0-2.09002e0*yw(i)**7.e0
      ratio(i)=app+(bpp/rv)
      end if
c
      if(wavemicr(i).ge.3.3e0.and.wavemicr(i).le.8.0e0) then

      if(wavemicr(i).ge.5.9e0.and.wavemicr(i).le.8.0) then
      zfa=-0.04473e0*((wavemicr(i)-5.9e0)**2.e0)-0.009779e0*
     #((wavemicr(i)-5.9e0)**3.e0)
      zfb=0.2130e0*((wavemicr(i)-5.9e0)**2.e0)+0.1207e0*
     #((wavemicr(i)-5.9e0)**3.e0)
      end if 
      if(wavemicr(i).lt.5.9e0) then
      zfa=0.0e0
      zfb=0.0e0
      end if
      app=zfa+1.752e0-0.316e0*wavemicr(i)-0.104e0/
     #(((wavemicr(i)-4.67e0)**2.e0)+0.341e0)
      bpp=zfb-3.090e0+1.825e0*wavemicr(i)+1.206e0/
     #(((wavemicr(i)-4.62e0)**2.e0)+0.263e0)
      ratio(i)=app+(bpp/rv)
      end if
c
      if(wavemicr(i).ge.8.0e0.and.wavemicr(i).le.10.0e0) then 
      app=-1.073e0-0.628e0*(wavemicr(i)-8.0e0)
     #+0.137e0*((wavemicr(i)-8.0e0)**2.e0)
     #-0.070e0*((wavemicr(i)-8.0e0)**3.e0)
      bpp=13.67e0+4.257e0*(wavemicr(i)-8.0e0)
     #-0.42e0*((wavemicr(i)-8.0e0)**2.e0)
     #+0.374e0*((wavemicr(i)-8.0e0)**3.e0)
      ratio(i)=app+(bpp/rv)
      end if
      ratio(i)=ratio(i)*avl
 3211 continue

c**********************************************************
c
c***************** start Girardi et al. (2000) formulas ********
      aa1=4.77e0 
      aa2a=alog10(4.e0*3.1415e0*100.e0)
      aa2c=alog10(5.67051e-5*(teff**4.e0)/(3.844e33))
      aa2b=2.e0*(alog10(3.0857e0)+18.e0)
      aa3=-2.5e0*(aa2a+aa2b+aa2c)
      aa4=0.0e0  ! all Vega magnitudes=0 in Vegamag HST
c      write(*,*) aa2a,aa2b,aa2c
c************* calcolo flusso totale (not necessary for stellar BCs
c
c      flussotot=0.0e0
c      do 7 i=2, 1221
c      flussomedio=(fs(i-1)+fs(i))/2.e0 
c      wavemedia=(waves(i)+waves(i-1))/2.0e0
c      flussotot=flussotot+waves(i)*
c     #flussomedio*(waves(i)-waves(i-1))    
c 7    continue
c      write(*,*) -2.5*alog10(flussotot)
c      
c************* calcolo flusso nel filtro
c

      do 2224 ikj=1,nfilt

      flussotot2=0.0e0
      flussotot2vega=0.0e0


      do 8 i=2, 1221
      flussomedio=(fs(i-1)+fs(i))/2.e0 
      rispostamedia=(ff2(ikj,i-1)+ff2(ikj,i))/2.e0
      wavemedia=(waves(i)+waves(i-1))/2.0e0
      flussotot2=flussotot2+waves(i)*flussomedio
     #*(10.e0**(-0.4e0*ratio(i)))*
     #rispostamedia*(waves(i)-waves(i-1))    
      flussomedio2=(fvega(i-1)+fvega(i))/2.e0
      flussotot2vega=flussotot2vega+waves(i)*flussomedio2*
     #rispostamedia*(waves(i)-waves(i-1))    
 8    continue
c**********
      bcV(ikj)=aa1+aa3+aa4+2.5e0*alog10(flussotot2/flussotot2vega)
c      write(5,*) teff, zlogg, bcV
c**********
c     write(*,*) -2.5e0*alog10(flussotot2)
c     write(*,*) -2.5e0*alog10(flussotot2vega) 
c      write(*,*) iz, 'correzione bolometrica=',bcV(ikj)

 2224 continue

      write(9,444) teff, zlogg, (bcV(ikj), ikj=1,nfilt)

 700  continue

 444  format(f8.0, f6.2, 51f8.3)
      close(2)
      close(9)
      stop
      end
