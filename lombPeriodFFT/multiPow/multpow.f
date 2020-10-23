c    program multpow
c  program to make power spectrum from multiple power spectra
c  obtained from individual segments of a lightcurve, allows selection
c  of segemnts according to their flux.  RXTE specific, needs background
c  lightcurve and observed lightcurve
      implicit none
      character*50 obsfl,bgdfl,psfl,limquery*1
      real*8 t(100000),flux(100000),ff(0:100000),err(100000)
      real*8 oflux(100000),bflux(100000),dum,logmean
      real*8 freq(100000),dff(100000),dflux(100000)
      real*8 bfreq(100000),bpow(100000),bperr(100000),bklev
      real*8 dt(100000),tlim,bsize,df,dav,minflux,maxflux
      real*8 bfac
      real*8 avoflux,avofluxtot,avflux,avfluxtot,tseg
      integer i,nhi,npoints,np,nf,j,minp,minb,nptot,nbfreq
      integer nseg

c..imh addition to read params from file

      open(unit=2,file='powparams',status='old')
      read (2,*), obsfl,bgdfl,bklev,bsize,bfac,minb,
     + tlim,minp,limquery,minflux,maxflux,psfl
      close (unit=2)

c      goto 888 

c      print *, "What input observed lightcurve shall I use?"
c      read *, obsfl
c      print *, "What is the name of the background lightcurve?"
c      read *, bgdfl
c      if (bgdfl.eq."none") then
c      print *, "What is the background level?"
c      read *, bklev 
c      endif
c      print *, "What is the bin time, and what PS bin factor and"
c      print *, "minimum number of bins shall I use?"
cc..bins means number of individual psd measurements which will be 
cc..averaged to make a binned-up bin, typically at least 10
c      read *, bsize,bfac,minb
c      print *, "What maximum gap size and minimum no. of points"
c      print *, "shall I use?"
cc..min no of datapoints in lightcurve from which to make powerspectrum
cc..at least a few tens
c      read *, tlim,minp
c      print *, "Shall I impose any flux limits on the parts of the"
c      print *, "lightcurve I use to make a power spectrum? (y/n)"
c      read *, limquery
c      if (limquery.eq."y") then
c      print *, "What is the lower and upper flux limit?"
c      read *, minflux,maxflux
c      else
c      minflux=-1E6
c      maxflux=1E6
c      endif
c      print *, "What file shall I write the PSD to?"
c      read *, psfl

c 888  continue

      print *, "What maximum segment size shall I use?"
      read *, tseg


c  first read in observed and background lightcurves and calculate
c  background subtracted lightcurve
      open (unit=2,file=obsfl,status='old')
      if (bgdfl.ne."none") then
      open (unit=3,file=bgdfl,status='old')
      endif
      do i=1,100000
        read (2,*,end=5), t(i),oflux(i),err(i)
        if (bgdfl.ne."none") then
        read (3,*), dum,bflux(i),dum
        else
        bflux(i)=bklev
        endif
        flux(i)=oflux(i)-bflux(i)       
      enddo
 5    npoints=i-1
      close (unit=2)
      if (bgdfl.ne."none") then
      close (unit=3)
      endif
c  set up initial values
      dav=0.0
      np=0
      nf=0
      nptot=0
      avflux=0.0
      avfluxtot=0.0
      avoflux=0.0
      avofluxtot=0.0
      nseg=0
c  now go through lightcurve and make power spectra
      do i=1,npoints
        np=np+1
        dt(np)=t(i)-t(i-np+1)
        dflux(np)=flux(i)
        dav=dav+dflux(np)
        avflux=avflux+flux(i)
        avoflux=avoflux+oflux(i)
        if (t(i+1)-t(i).ge.tlim.or.i.eq.npoints.or.dt(np).
     *      gt.tseg) then
        if (np.lt.minp) goto 100
        dav=dav/real(np)
        if (dav.lt.minflux.or.dav.ge.maxflux) goto 100
        avofluxtot=avofluxtot+avoflux
        avfluxtot=avfluxtot+avflux
        nptot=nptot+np
        df=1./dt(np)
c        print *,df,bsize,npoints,nf
c	      print *,"bfac,minb",bfac,minb
        call powcal(bsize,df,dt,dflux,np,ff,nhi,dav)
c        print *,"NHI outside powcal",nhi
        do j=1,nhi
          dff(j+nf)=ff(j)/dav**2
c          dff(j+nf)=ff(j)
          freq(j+nf)=df*j
        enddo
c	print *,nf
        nseg=nseg+1
        print *, nseg," segments recorded"
        nf=nf+nhi
 100    np=0 
c        dav=0.0
        avoflux=0.0
        avflux=0.0
        endif
      enddo

c  now sort the resulting array of power spectra 
      call sort2(nf,freq,dff)
c  now bin up the sorted array
      call binps(nf,freq,dff,minb,bfac,bfreq,bpow,bperr,nbfreq)
c  now write to the output file
      logmean=2*log10(avfluxtot/real(nptot))
      open (unit=2,file=psfl,status='new')      
c      write (2,*), "read serr 2"
      do i=1,nbfreq
c        write (2,*), bfreq(i),bpow(i)-logmean,bperr(i)
        write (2,*), bfreq(i),bpow(i),bperr(i)
      enddo
      close (unit=2)
c  finally, print out average count rate  of lightcurve used to make PS
      print *, "Average total count rate = ",avofluxtot/real(nptot)
      print *, "Average source count rate = ",avfluxtot/real(nptot)
      print *, "AVERAGE FLUX", dav
      end

      subroutine binps(nf,freq,dff,minb,bfac,bfreq,bff,bperr,nbfreq)
       implicit none
      real*8 freq(100000),dff(100000),bfac,flast,fprev,ffav,logftot
      real*8 bfreq(100000),bff(100000),bperr(100000)
      integer jj,j,nbfreq,nf,minb,i,k

      jj=0
      j=0
      logftot=0.0
      ffav=0.0
      flast=bfac*freq(1)
      fprev=freq(1)  
       
c       print *,"flast,fprev",flast,fprev
        
      do i=1,nf
         ffav=ffav+log10(dff(i))+0.253
         j=j+1
         logftot=logftot+log10(freq(i))
         if ((freq(i+1).gt.flast.or.i.eq.nf).and.j.ge.minb) then
            jj=jj+1
            bff(jj)=ffav/real(j)
            bfreq(jj)=logftot/real(j)
            bperr(jj)=0.0
            do k=i-(j-1),i
               bperr(jj)=bperr(jj)+(log10(dff(k))+0.253-bff(jj))**2
            enddo
c     The  following is the error on the mean.
            bperr(jj)=sqrt(bperr(jj)/(j*real(j-1)))
c     The following is simply the spread of the distribution.
c            bperr(jj)=sqrt(bperr(jj)/(j-1))
            fprev=freq(i)
            flast=bfac*freq(i)
            logftot=0.0
            ffav=0.0
            j=0  
         endif
      enddo
      nbfreq=jj
      
      return        
      end

   

      subroutine powcal(bsize,df,t,flux,npoints,ff,nhi,dav)
      implicit none
      real*8 t(1000000),flux(1000000),ff(0:1000000),a,pi
      real*8 fr(0:1000000),fi(0:1000000)
      real*8 c,s,df,dav,bsize
      integer i,j,npoints,nhi,counter
      parameter (pi=3.1415926536)

      nhi=int(0.5/(bsize*df))
       counter =0 
      do i=0,nhi
      counter=counter+1
        fr(i)=0.
        fi(i)=0.
        a=2*pi*i*df
	

        do j=1,npoints
          c=cos(a*t(j))
          s=sin(a*t(j))
          fr(i)=fr(i)+((flux(j)-dav)*c)
          fi(i)=fi(i)+((flux(j)-dav)*s)
        enddo

        ff(i)=((nhi*df)**(-1))*(fr(i)**2+fi(i)**2)/(real(npoints))
        
      enddo

	do j=0,20
c	print *,j,ff(j)
	enddo
c	print *,"counter",counter
      return
      end

      subroutine sort2(n,arr,brr)
      integer n,m,nstack
      real*8 arr(n),brr(n)
      parameter (m=7,nstack=50)
      integer i,ir,j,jstack,k,l,istack(nstack)
      real a,b,temp
      jstack=0
      l=1
      ir=n
 1    if (ir-l.lt.m) then
      do 12 j=l+1,ir
        a=arr(j)
        b=brr(j)
        do 11 i=j-1,1,-1
          if (arr(i).le.a) goto 2
          arr(i+1)=arr(i)
          brr(i+1)=brr(i)
 11     enddo
        i=0
 2      arr(i+1)=a
        brr(i+1)=b
 12   enddo
      if (jstack.eq.0) return
      ir=istack(jstack)
      l=istack(jstack-1)
      jstack=jstack-2

      else

      k=(l+ir)/2
      temp=arr(k)     
      arr(k)=arr(l+1)
      arr(l+1)=temp
      temp=brr(k)
      brr(k)=brr(l+1)
      brr(l+1)=temp

      if (arr(l+1).gt.arr(ir)) then
      temp=arr(l+1)
      arr(l+1)=arr(ir)
      arr(ir)=temp
      temp=brr(l+1)
      brr(l+1)=brr(ir)
      brr(ir)=temp
      endif

      if (arr(l).gt.arr(ir)) then
      temp=arr(l)
      arr(l)=arr(ir)
      arr(ir)=temp
      temp=brr(l)
      brr(l)=brr(ir)
      brr(ir)=temp
      endif

      if (arr(l+1).gt.arr(l)) then
      temp=arr(l+1)
      arr(l+1)=arr(l)
      arr(l)=temp
      temp=brr(l+1)
      brr(l+1)=brr(l)
      brr(l)=temp
      endif

      i=l+1
      j=ir
      a=arr(l)
      b=brr(l)

 3    continue
      i=i+1
      if (arr(i).lt.a) goto 3
 4    continue
      j=j-1
      if (arr(j).gt.a) goto 4
      if (j.lt.i) goto 5
      temp=arr(i)
      arr(i)=arr(j)
      arr(j)=temp
      temp=brr(i)
      brr(i)=brr(j)
      brr(j)=temp
      goto 3
 5    arr(l)=arr(j)
      arr(j)=a
      brr(l)=brr(j)
      brr(j)=b
      jstack=jstack+2

      if (jstack.gt.nstack) then 
         write (*,*), 'nstack too small in sort2'
         read (*,'()')
      endif
      if (ir-i+1.ge.j-l) then
      istack(jstack)=ir
      istack(jstack-1)=i
      ir=j-1
      else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
      endif
      endif
      goto 1

      end         

