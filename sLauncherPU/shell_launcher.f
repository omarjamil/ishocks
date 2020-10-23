      program simul8
c  program to simulate a lightcurve, using the method of Timmer & Konig
c  assuming various types of power-spectral model
c  Phil Uttley, Feb 2001
      implicit none
      character*30 slfl,lcfl,inpowfl
      character*1 simquery,gamq(10),massq(10),cmass, cwidth
      double precision t(10000000),dflux(10000000),bsize,mean
      double precision gflux(10000000),mflux(10000000)
      double precision a(10),f0(10),p(10,7),c,var,slrate
      double precision f,f1,f2,a_mod(10),mult,df,intpow,powgen
      double precision ran2,mtotintpow,gtotintpow,dp(7),rms(10)
      double precision mrms,grms,k_m,k_g,m_ind,g_ind
      double precision mfluxav,gfluxav,mfluxsqav,gfluxsqav
      double precision mvar,gvar,intmass,avgam,width, width_factor, wid
      integer npoints,idum,mod(10),i,ncomp,j

      open(unit=2,file="launch.in",status='old')

      read (2,*), npoints,bsize,k_m,m_ind,k_g,g_ind,idum,cmass,slrate
      read (2,*), slfl,lcfl
      read (2,*), width_factor
      read (2,*), cwidth,wid
      do i=1,10      
        read (2,*,end=10), a(i),f0(i),p(i,1),p(i,2),p(i,3),p(i,4),
     *p(i,5),p(i,6),p(i,7),mod(i),rms(i),gamq(i),massq(i)
      enddo
      
      close (unit=2)
 10   ncomp=i-1

      simquery='n'

c  c is just a constant can be used to add noise, not used in this code
      c=0.0
c  this bit works out the integrated powers for each PSD component and 
c  compares them with the required rms (if it is set) to rescale
c  PSD component normalisations to values that produce required rms.  
c  Also works out total powers, necessary for rescaling to take
c  account of effect of exponentiation, which increases the rms
      f2=1./(2.*bsize)
      f1=1./(real(npoints)*bsize)
      mult=10**(log10(f2/f1)/100000.0)
      mtotintpow=0.0
      gtotintpow=0.0
      do i=1,ncomp
        intpow=0.0
        do j=1,7
          dp(j)=p(i,j)
        enddo
        do j=1,100000
          f=f1*mult**(real(j))
          df=f-(f1*mult**(real(j-1)))
          intpow=intpow+(powgen(f,a(i),f0(i),dp,mod(i),c)*df)
        enddo
        if (rms(i).eq.0.0) then
          a_mod(i)=a(i)
        else if (rms(i).ne.0.0) then
          a_mod(i)=a(i)*((rms(i)**2)/intpow)
        endif
        if (massq(i).eq.'y') then
          if (rms(i).eq.0.0) then
            mtotintpow=mtotintpow+intpow
          else if (rms(i).ne.0.0) then
            mtotintpow=mtotintpow+rms(i)**2
          endif
        endif
        if (gamq(i).eq.'y') then
          if (rms(i).eq.0.0) then
            gtotintpow=gtotintpow+intpow
          else if (rms(i).ne.0.0) then
            gtotintpow=gtotintpow+rms(i)**2
          endif
        endif
      enddo
c  now, since we are exponentiating, modify a_mod again so that simulated
c  lightcurve will have desired variance/PS normalisation, and
c  modify mean also (must use fractional rms^2 units in PS or quote
c  fractional variance)


c  Now set time,mass and gamma time series to zero so we can sum components
      do i=1,npoints
        t(i)=real(i)*bsize
        mflux(i)=0.0
        gflux(i)=0.0
      enddo
c  Now sum the time series over all the contributing components
      do i=1,ncomp
        do j=1,7
          dp(j)=p(i,j)
        enddo
c  This subroutine makes the lightcurves using the Timmer & Koenig
c  (1995) method, they have mean 0 and are linear
        call simlc(t,npoints,dflux,bsize,0.0d0,a_mod(i),f0(i),dp,
     *idum,mod,c,1,var,simquery,inpowfl)
        print *, "Simulation variance = ",var
c  Sum and include corrections to variance to compensate for 
c  exponentiation later        
        do j=1,npoints
          if (massq(i).eq.'y') then
          mflux(j)=mflux(j)+(dflux(j)*(log(var+1.0)/var))
          endif
          if (gamq(i).eq.'y') then
          gflux(j)=mflux(j)+(dflux(j)*(log(var+1.0)/var))
          endif
        enddo
      enddo
c  Now exponentiate to get the final time series and record
c  mean and flux^2 to get variance for rescaling
      mfluxav=0.0
      gfluxav=0.0
      do i=1,npoints
c  Incorporate index for mass at this stage:
        mflux(i)=(exp(mflux(i)))**m_ind
        gflux(i)=(exp(gflux(i)))**g_ind
        mfluxav=mfluxav+mflux(i)
        gfluxav=gfluxav+gflux(i)
        mfluxsqav=mfluxsqav+mflux(i)**2
        gfluxsqav=gfluxsqav+gflux(i)**2
      enddo
      mfluxav=mfluxav/real(npoints)
      gfluxav=gfluxav/real(npoints)
      mfluxsqav=mfluxsqav/real(npoints)
      gfluxsqav=gfluxsqav/real(npoints)
      mvar=mfluxsqav-mfluxav**2
      gvar=gfluxsqav-gfluxav**2

      open (unit=2,file=lcfl,status='new')

c  Now rescale to mean=1
      do i=1,npoints
        mflux(i)=mflux(i)*(1./mfluxav)
        gflux(i)=1.0+((k_g-1.0)*(gflux(i)/gfluxav))
        write (2,*), t(i),mflux(i),gflux(i)
      enddo
      close (unit=2)
c  Now determine shell launch times and parameters, shells are launched
c  with either constant mass, or at constant intervals.  In either case
c  the key parameter is the shell launch rate, slrate, which in the case of
c  constant intervals just defines the time interval between launches,
c  (with mass determined by the flux integrated between those times)
c  whereas for a constant mass, the rate is the *average*, and is defined
c  so that for constant mflux(i)=1 there are slrate launches per unit time,
c  thus the mass per shell is set to be 1/slrate, and shells are launched
c  whenever the integrated mass since the last launch exceeds this value.
c  The final shell mass launched is then scaled by the normalisation
c  factor k_m.

      open (unit=2,file=slfl,status='new')

      intmass=0.0
      avgam=0.0
      j=0
      do i=1,npoints
c  integrate the mass
        j=j+1
        intmass=intmass+(mflux(i)*bsize)      
c  determine average gamma weighted by the present mass flux
        avgam=avgam+(mflux(i)*bsize*gflux(i))
c        avgam=avgam+(bsize*gflux(i))
c  first, for constant mass case:
        if (cmass.eq.'y'.and.intmass.ge.(1./slrate)) then
        avgam=avgam/intmass
c  use average gamma to determine shell width (in m)
        if (cwidth.eq.'n') then
        width=width_factor*sqrt(1-(1/avgam)**2)*real(j)*bsize*3*10**6
        else if(cwidth.eq.'y') then
        width=wid
        end if
        write (2,*), real(i)*bsize,k_m*intmass,avgam,width
        intmass=0.0
        avgam=0.0
        j=0
c  now for constant launch interval:
        else if (cmass.eq.'n'.and.(real(j)*bsize).ge.(1./slrate)) then
        avgam=avgam/intmass
        if (cwidth.eq.'n') then
        width=width_factor*sqrt(1-(1/avgam)**2)*real(j)*bsize*3*10**8
        else if(cwidth.eq.'y') then
        width=wid
        end if
        write (2,*), real(i)*bsize,k_m*intmass,avgam,width
        intmass=0.0
        avgam=0.0
        j=0
        endif
      enddo
      
      end
c
c
c
      subroutine simlc(t,npoints,flux,bsize,utbin,a,f0,p,idum,
     *      mod,c,tmult,var,simquery,inpowfl)
c  subroutine which simulates an underlying lightcurve of known power 
c  spectrum and then samples it with the same sampling of the observed 
c  lightcurve.
c  Uses powgen,four1,gasdev
      implicit none
      character*30 inpowfl,simquery*1
      double precision tmax,nyqf,a,f0,p(7),f,newt,c,var,utbin,bsize
      double precision minf,per(10000000),powgen,pow,t(10000000)
      double precision flux(10000000),dflux(10000000),avflux,tbin
      double precision infreq(0:100000),inpow(0:100000),merr,perr
      double precision delfreq,delpow
      double precision gasdev
      integer nf,i,ii,idum,j,npoints,k,mod,tmult,nut,ninfreq
      integer jj

      
c      if (utbin.eq.0.0) then
c      tbin=bsize
c      else
c      tbin=utbin
c      endif
      tbin=bsize
      nyqf=tbin*2.0
      tmax=nyqf*2.0**(real
     *   (int(log10((t(npoints)-t(1))/nyqf)/log10(2.0))+tmult))
      print *, nyqf,tmult
      minf=1/tmax
      nf=nint(tmax/nyqf)
      per(1)=0.0
      per(2)=0.0
        
      if (simquery.eq.'y') then
      open (unit=4,file=inpowfl,status='old')
      do i=1,1000000
        read (4,*,end=15) infreq(i),inpow(i),merr,perr
         inpow(i)=inpow(i)/infreq(i)
         inpow(i)=a*inpow(i)
c        print *, infreq(i),inpow(i),merr,perr       
      enddo
 15   ninfreq=i-1
      close (unit=4)
      endif
      infreq(0)=0.0
      inpow(0)=inpow(1)
      print *, simquery,nf,minf
      jj=0
      do i=3,((2*nf)+1),2
        f=(real(i-1)/2.0)*minf
        if (simquery.eq.'y') then
        if (f.gt.infreq(jj)) then
        jj=jj+1
        endif
        if (jj.le.ninfreq) then
        delfreq=infreq(jj)-infreq(jj-1)
        delpow=inpow(jj)-inpow(jj-1)
        pow=inpow(jj-1)+(((f-infreq(jj-1))/delfreq)*delpow)
c        print *, delfreq,delpow,pow
        else
        pow=inpow(ninfreq)
        endif
c        pow=pow/f
        else
        pow=powgen(f,a,f0,p,mod,c)
        endif
        pow=pow*(2*nf/nyqf)
        per(i)=gasdev(idum)*sqrt(0.5*pow)
        per(i+1)=gasdev(idum)*sqrt(0.5*pow)
        ii=(4*nf)-(i-2)
        per(ii)=per(i)     
        per(ii+1)=-per(i+1)
      enddo
 
      per((2*nf)+1)=0.0  
        
      call four1(per,2*nf,-1)
        
      avflux=0.0
      k=0
      j=1
      do i=1,4*nf-1,2 
        k=k+1
        per(i)=per(i)/real(2*nf)
        newt=((i-1)/2)*(nyqf/2.0)
        dflux(k)=per(i)
c        if (utbin.eq.0.) then
c        if (newt.ge.(t(j)).and.((t(j)-t(j-1)).
c        if (newt.ge.(t(j)-t(1)).and.((t(j)-t(j-1)).
c     *    ge.(nyqf/2.0).or.j.eq.1)) then
        
        flux(j)=per(i)
        avflux=avflux+flux(j)
        j=j+1
c        else if ((t(j)-t(j-1)).lt.(nyqf/2.0).and.j.ne.1) then
c        flux(j)=flux(j-1)
c        avflux=avflux+flux(j)
c        j=j+1
c        endif
c        else
c        if (newt.ge.(t(j)-t(1)).and.newt.lt.(t(j)-t(1)+bsize)) then
c        nut=nut+1        
c        flux(j)=flux(j)+per(i)
c        else if (newt.ge.(t(j)-t(1)+bsize)) then
c        flux(j)=flux(j)/real(nut)
c        j=j+1
c        flux(j)=0.0
c        nut=0
c        endif        
c        endif

        if (j.gt.npoints) goto 444
      enddo
     
 444  continue
        
      var=0.0
      avflux=avflux/real(npoints)
      do i=1,npoints
        var=var+(flux(i)-avflux)**2
      enddo
      var=var/(npoints-1)

      return
        
      end


      function powgen(f,a,f0,p,mod,c)
      implicit none
      double precision f,a,f0,p(7),c,qpo,sig
      double precision powgen,mmed,mp1,mp2,mf0
      double precision fcent1,fcent2,fcent3,q1,q2,q3
      integer mod
 
      if (mod.eq.1) then
      if (f.lt.p(1).and.f.ge.p(2)) then
      powgen=(a*(f/f0)**(-p(4)))+c
      else if (f.lt.p(2)) then
      powgen=a*((p(2)/f0)**(p(5)-p(4)))*(f/f0)**(-p(5))+c
      else
      powgen=a*((p(1)/f0)**(p(3)-p(4)))*(f/f0)**(-p(3))+c
      endif
 
      else if (mod.eq.2) then
      mmed=(1.+(f/p(2))**p(4))**(-(p(3)/p(4)))
      powgen=(a*mmed)+c 

      else if (mod.eq.3) then
      
      mf0=f0/p(6)
      mp1=p(1)/p(6)
      mp2=p(2)/p(6)
      
      if (f.lt.mp1.and.f.ge.mp2) then
      powgen=(a*(f/mf0)**(-p(4)))+c
      else if (f.lt.mp2) then
      powgen=a*((mp2/mf0)**(p(5)-p(4)))*(f/mf0)**(-p(5))+c
      else
      powgen=a*((mp1/mf0)**(p(3)-p(4)))*(f/mf0)**(-p(3))+c
      endif

      else if (mod.eq.4) then

c     q=p(7) so work out FWHM sigma

      sig=p(6)/p(7)
      qpo=(p(5)/3.14159)*((sig/2.)/((f-p(6))**2.+(sig/2.)**2.))

      powgen=qpo
      
c      if (f.lt.p(1).and.f.ge.p(2)) then
c      powgen=(a*(f/f0)**(-p(4)))+c+qpo
c      else if (f.lt.p(2)) then
c      powgen=a*((p(2)/f0)**(-p(4)))*(f/f0)**(0.0)+c+qpo
c      else
c      powgen=a*((p(1)/f0)**(p(3)-p(4)))*(f/f0)**(-p(3))+c+qpo
c      endif
      
c Double Lorentzian peak
      else if (mod.eq.9) then
      q1=0.01
      q2=0.125
      q3=5

      fcent1=p(1)/sqrt(1+1/(4*q1**2))
      fcent2=p(2)/sqrt(1+1/(4*q2**2))
      fcent3=p(3)/sqrt(1+1/(4*q3**2))
      powgen=a*q1*fcent1/(fcent1**2 + 4*q1**2*(f-fcent1)**2)
     *+a*p(4)*q2*fcent2/(fcent2**2+4*q2**2*(f-fcent2)**2)+c
     *+a*p(5)*q3*fcent3/(fcent3**2+4*q3**2*(f-fcent3)**2)+c
      endif
      
      return
      end

      subroutine four1(data,nn,isign)
      integer isign,nn
      double precision data(2*nn)
      integer i,istep,j,m,mmax,n
      double precision tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp
      
      n=2*nn
      j=1
      do 11 i=1,n,2
        if (j.gt.i) then
        tempr=data(j)
        tempi=data(j+1)
        data(j)=data(i)
        data(j+1)=data(i+1)
        data(i)=tempr
        data(i+1)=tempi
        endif
        m=n/2
 1      if ((m.ge.2).and.(j.gt.m)) then
        j=j-m
        m=m/2
        goto 1
        endif
        j=j+m
 11   enddo
      mmax=2
 2    if (n.gt.mmax) then
      istep=2*mmax
      theta=6.28318530717959d0/(isign*mmax)
      wpr=-2.d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.d0
      wi=0.d0
      do 13 m=1,mmax,2
        do 12 i=m,n,istep
          j=i+mmax
          tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
          tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
          data(j)=data(i)-tempr
          data(j+1)=data(i+1)-tempi
          data(i)=data(i)+tempr
          data(i+1)=data(i+1)+tempi
 12     enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
 13   enddo
      mmax=istep
      goto 2
      endif
      return
      end


      function gasdev(idum)
      integer idum
      double precision gasdev
      integer iset
      double precision fac,gset,rsq,v1,v2,ran2
      save iset,gset
      data iset/0/
      if (iset.eq.0) then
 1    v1=2.*ran2(idum)-1.
      v2=2.*ran2(idum)-1.
      print *, "Random 1 = ",v1
      print *, "Random 2 = ",v2
      rsq=v1**2+v2**2
      if (rsq.ge.1..or.rsq.eq.0.) goto 1
      fac=sqrt(-2.*log(rsq)/rsq)
      gset=v1*fac
      gasdev=v2*fac
      print *, "gas, gset = ",gasdev,gset
      iset=1
      else
      gasdev=gset
      print *, "gasdev= ",gasdev
      iset=0
      endif
      return
      end


      function ran2(idum)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      double precision ran2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *   ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     *   ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-6,rnmx=1.-eps)

      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if (idum.lt.0) idum=idum+im1
          if (j.le.ntab) iv(j)=idum
 11     continue
        iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if (iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)
      return
      end

      function poidev(xm,idum)
      integer idum
      double precision poidev,xm,pi
      parameter (pi=3.141592654)
      double precision alxm,em,g,oldm,sq,t,y,gammln,ran2
      save alxm,g,oldm,sq
      data oldm /-1./

      if (xm.lt.12.) then
      if (xm.ne.oldm) then
      oldm=xm
      g=exp(-xm)
      endif
      em=-1
      t=1.
 2    em=em+1.
      t=t*ran2(idum)
      if (t.gt.g) goto 2
      else
      if (xm.ne.oldm) then
      oldm=xm
      sq=sqrt(2.*xm)
      alxm=log(xm)
      g=xm*alxm-gammln(xm+1.)
      endif
 1    y=tan(pi*ran2(idum))  
      em=sq*y+xm
      if (em.lt.0.) goto 1
      em=int(em)
      t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)   
      if (ran2(idum).gt.t) goto 1
      endif
      poidev=em
      return
      end

**************************************************************************



      FUNCTION GAMMLN(XX)
**************************************************************************



C  This function returns the value log gamma(XX) for XX>0.
C  Actually this function works only for XX>1.0. (It needs more for
C  XX<1).

      IMPLICIT NONE

C  Variables for the argument of the function and the result.
      double precision XX,GAMMLN
      
C  Other variables.
      double precision COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      INTEGER J 

      DATA COF,STP/76.1800,-86.5053,24.0141,-1.2317,
     +             0.1209E-2,-0.5364E-5,2.5066/
      DATA HALF,ONE,FPF/0.5,1.0,5.5/

**********************************************************************
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
         X=X+ONE
         SER=SER+COF(J)/X
   11 CONTINUE
      GAMMLN=TMP+LOG(STP*SER)

      RETURN
      END



