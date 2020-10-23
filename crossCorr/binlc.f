      program binlc
c
      implicit none
      character*50, inlc,outlc
      double precision t(10000000),flux(10000000),err(10000000)
      double precision derr(1000000),data(1000000),tstart,tbin
      double precision bt(4000000),bflux(4000000),mu,averr
      double precision berr(4000000),dum,oldbin
      integer i,j,k,npoints,ndata,nbins

C$$$      print *, "What lightcurve shall I bin up?"
C$$$      read *, inlc
C$$$      print *, "What is the bin time and new bin time?"
C$$$      read *, oldbin,tbin
C$$$      print *, "What file shall I write to?"
C$$$      read *, outlc
      open (unit=1, file="rebin.par", status='old')
      read (1, *), inlc
      read (1, *), oldbin,tbin
      read (1, *), outlc

      open (unit=2,file=inlc,status='old')
      open (unit=3,file=outlc,status='new')

      do i=1,10000000
c          print *, i
          read (2,*,end=10),t(i),flux(i),err(i)
      enddo
 10   npoints=i-1
 
      nbins=0
      tstart=t(1)
      ndata=0
      j=0
      do i=1,npoints
        j=j+1
        data(j)=flux(i)
        k=int((t(i)-tstart)/tbin)+1
        derr(j)=err(i)
        if ((int((t(i+1)-tstart)/tbin)+1).gt.k.or.i.eq.npoints) then
c        if ((t(i+1)-t(i-j+1)).lt.(tbin/2.0)) goto 50
c        if ((real(j)*oldbin).lt.(tbin/2.0)) goto 50
        ndata=j
        call wterrmu(data,derr,ndata,mu,averr)
        nbins=nbins+1
        bflux(nbins)=mu
        bt(nbins)=(k*tbin)+tstart
        berr(nbins)=averr
 50   j=0      
        endif   
      enddo

C$$$  do i=1,nbins
C$$$        write (3,*),(bt(i)-(tbin/2.0)),bflux(i),berr(i)
C$$$      enddo

C$$$  Not writing the error column - OJ
      do i=1,nbins
        write (3,*),(bt(i)-(tbin/2.0)),bflux(i)
      enddo
      print *, npoints," points read"
      print *, nbins," points written"
      end


      subroutine wterrmu(data,err,ndata,mu,averr)
c  given input values and errors, produces a single weighted output mean
c  and error
      implicit none
      double precision mutot,errtot,data(100000),err(100000),mu,averr
      integer i, ndata
        
      
      mutot=0.0
      errtot=0.0
      do i=1,ndata
        mutot=mutot+data(i)
        errtot=errtot+err(i)**2
      enddo 
      
      mu=mutot/real(ndata)
      averr=sqrt(errtot)/real(ndata)
        
      end
