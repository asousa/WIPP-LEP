       program newray
c      new version of main program-ray tracing     wjbmar73
c-------
c	brought onto star-vax by beth pridgen, mar 2 1983
c	revisited by ngo d. hoc, april 9, 1985.
c	last modified by D. Lauben, June 26, 1995; Aug 18, 1997
c-------
c	N.D. Hoc modifications:
c		input and output units use files
c		funct modified to prevent math errors, exponents too large for
c		vax precision.
c	D. Lauben modifications:
c		write to unit 8 now includes phase time (tt) and L-value (ele)
c		added stop when ele reaches -tgfina when tgfina < 0
c           write an and ani(i),i=1,4 to unit 8 in ADAMS()
c		eliminated plot calls as necessary to compile on Nova5 (unix)
c		also altered calls to date(),time() and secnds() to become 
c		   fdate(), ctime(time()), and dtime(), as req'd.
c		   resized dates,times,xtime0,xtime1 arrays as req'd.
c		upgraded to implicit real*8 with double precision math calls
c               altered auto-step algorithm to upper bound h to its initial 
c                  value (thru hmax variable)
c               WILL clobber existing output files (alternate is frustrating)
c               provided explicit unit 8 name(newray.dat) and format statement
c                 to avoid annoying carriage-returns on non-DEC unix.
c-------
c previous comments from Ngo Hoc:
c
c	input data file renamed to newray.in, instead of raytrace.in
c	output data file rename to newray.out
c       add "numres" to reduce step size when resonance occurs.
c	print date and time on output.   
c       use free-format read except character read.
c       check geocentric distance p(1,j) at the beginning of funct 
c       to prevent potential overflow in dens at lines 120, 130 and 140.      
c       delete variable "kn" in input data.        
c       add "nsuppr" in input data to suppress (if nsuppr=1) all output points
c       whose latitudes lie outside satellite trajectory.
c	add interactive option and put "interac, numres, nsuppr, spelat" 
c	on the first data card.
c	read only "distre" and "l", no "dummmm". 
c	handle all raytracing plots in one call to rayplt, thus eliminate
c	the extra variable "kplot" in read statement: "num kskip mode ..."
c-------
	implicit real*8 (a-h)
	implicit real*8 (o-z)

c      dimension glab1(50),glab2(50)
c      character dates*9, times*8 
      character*24 dates, times, fdate, ctime
      integer*4 time

      real*8 xdra(2000),ydra(2000)
      real*8 l,l0,lk,latitu,grarad,gf,latmax,latmin
      real xtime0,xtime1

      common /resona/ nreson, ndecre, numres
      common /ryplt/ freq, altinj
      common /x/ latitu,latmax,latmin, xtime0(2), xtime1(2), nsuppr
      common a,absb,al,alpha(4),alpha0(4),am,an,an1,an2,ane0,ani(4),
     1 ap,apl,apr,aps,ar,arl,as,b,c,cosz2,cosz22,cpsi,cpsi2,dd(10),ddk,
     2 def(10),delt,dlgndr(4),dlgndt(4),egfeq,expk,f0e,fkc,gf,h,hmin,
     3 l,l0(10),lk,pi,print,psi,psires,r0,radgra,rbase,rconsn,rdiv,relb,
     4 rstop,rzero,scbot,scr,sidedu(10),sinz2,sinz22,spsi,spsi2,tgfina,
     5 therm,tiniti,x(4),x0(5),y(4),y01,y02,z(5),clight,consf0,
     6 hl2n(10),hl2s(10),hu2n(10),hu2s(10),rducln(10),rducls(10),
     7 rducun(10),rducus(10),durbe,refalt,hmax,
     8 init,kducts,kinit,kj,kount,kr,kskip,kt,ktape,lines,mode,num,
     9 numray
	
	dimension hducls(10), hducln(10), hducus(10), hducun(10)
	
c      data pi/3.141592653589793/,clight/2.99792501 e05/,
c     1      consf0/80.61638449/,r0 /6370./

      pi = 3.141592653589793
      clight = 2.99792501 e05
      consf0 = 80.61638449
      r0 = 6370.

c
c     constants       (consf0 = q**2 / 4*pi*pi*e0*m)
c-----
      open (unit=4, file='newray.in', status='old')
      open (unit=7, file='newray.out', status='unknown')
      open (unit=8, file='newray.dat', status='unknown')
      open (unit=3, file='Lprofile.dat', status='unknown')
      open (unit=9, file='latprofile.dat', status='unknown')

      print *, 'newray v1.8 Aug 97'

      zero = 0.0
    
      radgra = 180./pi
      grarad = 1./radgra
      am = 1.67252 e-27 / 9.10955854 e-31
      read (4,*) intera, numres, nsuppr, spelat
      if (intera .eq. 0) goto 50
         write (*,12)
   12    format (' please enter maximum number of tracing ',
     1   'continuations at resonant cone')
	 read (*,*) numres
	 write (*,14)
   14	 format (' enter 1/0 if suppressing output outside ',
     1   'satellite range is or is not desired')
         read (*,*) nsuppr
   50 kinit=2
      init=-1
c
c      call date (dates)
c      call time (times)
c      xtime0=secnds(0.0)

	dates = fdate()
	times = ctime( time() )
	dummmm = dtime(xtime0)
	
c     read in the satellite coordinates
c
      count=0.0
      write(*,90)
   90 format (' read input')
      latmax=-100.
      latmin=100.
  105 continue
c  	  print *, 'past 105'
  	  read (4,*,end=9501) distre, latitu
      dummmm = 0.0
      if (distre .le. -1.) dummmm = 99999.0
      distre = distre/r0
      if (latitu .le. latmin) latmin=latitu
      if (latitu .ge. latmax) latmax=latitu
  130    continue
c  	  print *, 'past 130'
            if (count .gt. 0.0) read(8,*,end=9501) dummmm,distre,latitu,
     1         (zz,i=1,11),gf,ap,ar,al
            if (count .lt. 1.0) then
c234567890c234567890c234567890c234567890c234567890c234567890c23456789012
               write(8,88) real(dummmm),real(distre),real(latitu),
     1           real(pi),real(zero),real(zero),real(zero),
     2           real(zero),real(zero),real(an),(real(ani(i)),i=1,4),
     3           real(gf),real(ap),real(ar),real(al)
c               write(8,*) char(10)
 88            format(18(1x,e15.8))
            endif
            if (dummmm .gt. 99998.0) go to 140
               if (count .gt .0.0) go to 130
                  if (count .lt. 1.0) go to 105
  140 continue
c     if latmax=latmin, then no satellite trajectory is given, disable
c     variable "nsuppr", latmax and latmin play no role.
      if (abs(latmax-latmin) .le. 1.e-5) nsuppr=0
      latmax=latmax+3.
      latmin=latmin-3.

c     input of model data
      read(4,*,end=9501) num,kskip,mode,kount,kducts,ktape,
     1    refalt,dsrrng,dsrlat,dsdens
  890     format (5i5, i2, f11.2, 3f10.3, i5)
          if (num .eq. 0) go to 9501
          write(7,410) dates, times
  410	  format (1h1,20x,'===== newray.out ===== date : ',a9,
     1        ',   time: ',a8,' ====='//)
          write(7,31)
   31     format(1x,23hdata cards for this run/)
          write(7,895) intera, numres, nsuppr, spelat
  895     format (3i5,f8.2)
          write(7,890)num,kskip,mode,kount,kducts,ktape,
     1         refalt,dsrrng,dsrlat,dsdens

          kfile=ktape
  910     read (4,*,end=9501) egfeq,therm,hm,absb,relb
          write(7,900) egfeq,therm,hm,absb,relb
  900     format(4f15.3,f15.4)

          read (4,*,end=9501) rbase,ane0,(alpha0(i),i=2,4)
          write(7,900)rbase,ane0,(alpha0(i),i=2,4)

          read (4,*,end=9501) rzero,scbot,rstop,rdiv,hmin
          write(7,900) rzero,scbot,rstop,rdiv,hmin

          if (kducts.eq.0) go to 500
             read (4,*,end=9501) lk,expk,ddk,rconsn,scr
             write(7,900) lk,expk,ddk,rconsn,scr

             if (kducts.eq.1) go to 500
                do 935 k=2,kducts
                   read (4,*,end=9501) l0(k),def(k),dd(k),
     1             rducln(k),hducln(k),rducun(k),hducun(k),
     2             rducls(k),hducls(k),rducus(k),hducus(k),sidedu(k)
c 940              format (3f5.3,9f7.1)    
                   write(7,940) l0(k),def(k),dd(k),rducln(k),hducln(k),
     1             rducun(k),hducun(k),rducls(k),hducls(k),rducus(k),
     2             hducus(k),sidedu(k)
  940              format (4f15.3,2f15.4,/,4f15.3,2f15.4)
 	           hl2n(k)=hducln(k)**2
	           hl2s(k)=hducls(k)**2
     	           hu2n(k)=hducun(k)**2
  935           hu2s(k)=hducus(k)**2
  500 read (4,*,end=9501) pstalt,palt1,palt2,platit,pstlat,plat1,plat2,
     1    paltit
      write(7,*) pstalt,palt1,palt2,platit,pstlat,plat1,plat2, paltit
c
c    setting ane0 to desired value at dsrrng,dsrlat
      z(2)=(90.00-dsrlat)*grarad
      z(1)=dsrrng*r0
      call dens
      ane0=ane0*dsdens/ani(1)
c
      if (pstalt .eq. 0.0) go to 520
c        density profiles
         z(2) = (90.00-platit)*grarad
         z(1) = palt1*r0
         pstep=pstalt*r0
         pfinal=palt2*r0
         write(7,505) platit
  505    format('1',14x,'alt(re)',15x,'l',15x,
     1   'el/cm3   density profile at latitude=',f4.1,' deg')
         i=0
  510    call dens
            palt=z(1)/r0
c           write(7,515) palt,l,ani(1)
  515       format(3f20.3)
            z(1)=z(1)+pstep
c
c           plotting of density profile on ncar cdc7600
            i=i+1
            xdra(i)=palt
            ydra(i)=ani(1)
            write(3,*) xdra(i),ydra(i)     ! altitude(Re),Ne
         if (z(1) .le. pfinal) go to 510

c            call agseti(6HY/LOG.,1)
c            call agsetf(11HLABEL/NAME.,1HL)
c            call agseti(12HLINE/NUMBER.,100)
c            call agsetp(10HLINE/TEXT.,25HELECTRON DENSITY (EL/CC)$,1)
c            call agsetf(11HLABEL/NAME.,1HB)
c            call agseti(12HLINE/NUMBER.,-100)
c            call agsetp(10HLINE/TEXT.,21HRADIAL DISTANCE (RE)$,1)
c            encode(42,32,glab1)platit
   32       format(' density profile at ',f4.1,' deg latitude     ')
c            call ezxy(xdra,ydra,i,glab1)

  520 if (pstlat .eq. 0.0) go to 60
         z(1)=paltit+r0
         z(2)=(90.00+plat1)*grarad
         pstep=pstlat*grarad
         pfinal=(90.00+plat2)*grarad
         write(7,525) paltit
  525    format('1',16x,'lat',17x,'l',15x,
     1      'el/cm3   density profile at altitude=',f7.1,' km')
         i=0
  530    call dens
            plat=z(2)*radgra-90.00
c           write(7,515) plat,l,ani(1)
            z(2)=z(2)+pstep
c
c           plotting of density profile on ncar cdc7600
            i=i+1
            xdra(i)=plat
            ydra(i)=ani(1)
            write(9,*) xdra(i),ydra(i)     ! latitude,Ne
         if ( z(2) .le. pfinal) go to 530

c            call agseti(6HY/LOG.,0)
c            call agsetf(11HLABEL/NAME.,1HL)
c            call agseti(12HLINE/NUMBER.,100)
c            call agsetp(10HLINE/TEXT.,25HELECTRON DENSITY (EL/CC)$,1)
c            call agsetf(11HLABEL/NAME.,1HB)
c            call agseti(12HLINE/NUMBER.,-100)
c            call agsetp(10HLINE/TEXT.,27HGEOMAGNETIC LATITUDE (DEG)$,1)
c            encode(40,33,glab2)paltit
   33       format(' density profile at ',f7.1,' km altitude ')
c            call ezxy(xdra,ydra,i,glab2)
c...  numray = number of rays to be traced.  
      numray=0     

      print *, 'execution'
   60 continue
      read (4,*,end=9501) fkc,x0(1),latitu,delt,psi,tiniti,x0(5),tgfina
c
c     input of ray data
c     frequency in khz, and initial values of geocentric distande in km,
c     latitude, delta -- angle from upward vertical to wave normal,
c     psi -- angle from geomagnetic field to wave normal (pos. inward),
c     initial phase time, initial group time, + maximum group time.
c     all angles in degrees, times in seconds. mag field up at pos lats.
c     the value of delta is used if it is different from 360 deg,
c     psi is used if delta is set equal to 360 deg.
c
      specul=0.0
      numray=numray+1
   63 durbe=0.0
      nreson=0  
      no584=0
c
c     use nreson as a counter to count number of times correction on
c     the psi angle is made (when /psi/ > /psiresonance/) for one ray.
c
      if (fkc .ge. 0.0) goto 61
         tgfina=tgfold
         go to 9000
   61 if (fkc .eq. 0.0) go to 50
      freq=fkc
      altinj=x0(1)-r0
      tgfold=tgfina
      init=-1
      x0(2)=90.00-latitu

c     output of heading heading heading
c      call date (dates) 
c      call time (times)
	dates = fdate()
	times = ctime(time())

      write(7,961) fkc,x0(1),latitu,delt,psi,tiniti,x0(5),tgfina,
     1   dates,times	
  961 format ('1',8f11.3,5x,'===== ',a9,',   ',a8)
c      xtime1=secnds(0.0)
      dummmm = dtime(xtime1)
      write(7,970) egfeq,therm,hm,absb,relb,kskip,mode,num
  970 format ('0','el gf eq=',f6.1,'    temp=',f7.1,'    hm=',f7.3,
     1'    absb=',f6.3,'    relb=',e9.2,'    kskip=',i2,'  mode=',i2,
     2'   no. comp.=',i2)
      write(7,980) rbase,ane0,(alpha0(i),i=2,4 )
  980 format (' ','rbase=',f10.4,'    density=',f15.4,'    h+=',f15.8,
     1'    he+=',f15.8,'    o+=',f15.8)
      write(7,981) rzero,scbot,rstop,rdiv,kducts
  981 format (' ','rzero=',f8.2,'    scbot=',f8.2,'    rstop=',
     1   f8.2,'    rdiv=' ,f8.2,'    kducts=',i3)
      lines=5

      if (kducts .eq. 0)  go to 998
c        plasmapause
         lines=lines+3
         s2t=x0(1)/(r0*lk)
         if (s2t .lt. 1.0) goto 600
            theta=0.
            thet=0.
            go to 991
  600    tant=sqrt((1.0-s2t)/s2t)
         theta=datan(tant)
         thet=radgra*theta
         if (x0(2).lt.90.0) go to 991
            thet=180.00-thet
  991    s2t=x0(1)/(r0*(lk+ddk))
         if (s2t.lt.1.0) go to 620
            width=0.
            go to 630
  620    tant=sqrt((1.-s2t)/s2t)
         width=x0(1)*(datan(tant)-theta)
  630    write(7,992) lk,thet,expk,ddk,width,rconsn,scr
  992    format('0','knee',f10.5,f15.5,6x,'(r**-',f5.3,')',f24.6,f11.4,
     1'     r of const. dens. ',f7.1,'    scr=',f7.1)
         if  (kducts .eq.1) go to 700
            if (l0(2) .gt. 0.0) go to 999
c              sinusoidal density perturbation
  800          al02=-l0(2)
               kinit=3
               write(7,810) def(2),dd(2),al02,rducln(2),
     1             rducls(2),rducun(2),rducus(2),hducln(2),
     2             hducls(2),hducun(2),hducus(2),sidedu(2)
  810          format(' ','sinusoidal dens pert   ampl=',f7.4,
     1         '   period=',f7.4,' l, with rel max at l=',f7.4/9f10.2)
               lines=lines+1
               if (kducts.eq.2) goto 700
c                 duct(s)
  999       write(7,990) x0(1)						
  990       format (' ','duct',4x,'l',4x,'thet(',f7.1,') enhcmt ',
     1         'smwth eq  smwth r   rducln  rducls  rducun  rducus  ',
     2         'hducln  hducls  hducun  hducus  sidedu')		
            lines=lines+1

            do 993 kduc=kinit,kducts
               s2t=x0(1)/(r0*l0(kduc))
               if (s2t.lt.1.0) go to 640
                  theta=0.
                  thet=0.
                  go to 997
  640          tant=sqrt((1.-s2t)/s2t)
               theta=datan(tant)
               thet=radgra*theta
               if (x0(2) .lt. 90.0) go to 997
                  thet=180.00-thet
  997          s2t=x0(1)/(r0*(l0(kduc)+dd(kduc)))
               if (s2t.lt.1.0) go to 650
                  width=0.
                  go to 993
  650          tant=sqrt((1.-s2t)/s2t)
               width=x0(1)*(datan(tant)-theta)
	       lines=lines+1
  993          write(7,994) kduc,l0(kduc),thet,def(kduc),dd(kduc),width,
     1         rducln(kduc),rducls(kduc),rducun(kduc),rducus(kduc),
     2         hducln(kduc),hducls(kduc),hducun(kduc),hducus(kduc),
     3         sidedu(kduc)
  994          format (' ',i3,f9.3,f11.3,2f9.3,f10.3,1x,8f8.1,f5.1)
               lines=lines+1
  700          write(7,701)
  701          format(' the density model is not good in the ',
     1         'region which is both below the r of const.',
     2         '  dens.  and  outside the knee')

c              recalculate delta using psi if delta equals 360.0
  998 if (delt.ne.360.0) go to 995
         alat=latitu*grarad
         beta=90.00-radgra*datan(2.0*dsin(alat)/dcos(alat))
         delt=beta+psi

  995 write(7,996) fkc,x0(1),latitu,delt
  996 format('0','frequency,khz=',f12.4,'   initial     r=',f10.4,
     1'     latitude=',f10.5,'   delta=',f10.5)
      lines=lines+2
      delt=delt*grarad
   74 x0(2)=x0(2)*grarad
  100 kr=1
c     transfer to next statement if kj=2 ,radial dist less than rdiv
c     - see 50 and 130 in funct
  110 h=3.12500*hm/sqrt(fkc)
      hmax = h
      print=16*iabs(kount)*h
      printt=print/clight
  150 if (kr .ne. 2) go to 210
  160    kr=1
  170    go to 290
  210 z(1)=x0(1)
  220 z(2)=x0(2)

  230 call mu

  260 if (x0(1) .lt. rdiv) go to 290
  270    kr=2
  280    go to 360
  290 if(lines.lt.54)go to 292
         write(7,322)
  322    format(1h1)
         write(7,321)
         lines=lines-57
  292 continue

      write(7,291) hm,h,hmin,kount,print,ktape,kfile,nsuppr
  291 format(' hm=',f8.3,'   h=',f8.2,'   hmin=',f8.2,'   kount=',i3,
     1 '   print=',f8.2,'   ktape=',i3,'   kfile=',i3,
     2 '   nsuppr=',i2)
      if (kount .eq. 0 .and. nsuppr .eq. 1) goto 315
      if (kount) 305,300,310
  300    write(7,301)
  301    format(' complete output   all computed points printed'//)
         go to 320
  305 write(7,306)
  306 format(' abbreviated output   only critical points printed'//)
      go to 320
  310 write(7,311) printt
  311 format(' bounded output   printed points spaced by dt=',f7.4,//)
      goto 320
  315 write(7,316)
  316 format (' suppressed output:  only points in satellite trajectory'
     1   ,' range printed'//)
  320 write(7,321)
  321 format('    t     tg   dist(re)  alt(km)   lat     l     f/gfeq'
     1 /'f/gf     dens     lhr    mu   mucos  delta   psi   psiray '
     2 /'psires    e//(kev)')
      lines=lines+6
c
	write(6,330) numray,fkc,x0(1),latitu,radgra*delt,radgra*psi
  330	format(1x,i4,' fkc=',f7.2,' R=',f8.1,' lat=',f7.2,
     1    ' delt=',f7.2,' psi=',f7.2)
  
  340 call adams

      if (durbe .lt. 0.0) go to 63
      if (durbe.gt.0.00 .and. numres .eq. 0) go to 1060
c        if nreson > numres, goto 1060 to process another ray input
         if (nreson .gt. numres) goto 1060
c        if kj=1,2,or 3 go to 360,110,or 1060 respectively upon 
c        return fm ads
  341    if (kj-2) 360,110,1060
  360       h=50.00*hm/sqrt(fkc)
            hmax = h
  400       go to 290
 1060 if (ktape.eq.0) go to 40
         f99=99999.0
c234567890c234567890c234567890c234567890c234567890c234567890c23456789012
         write(8,88) real(f99),real(l),real(lo),real(lk),real(zero),
     +     real(zero),real(zero),real(zero),real(zero),real(an),
     +	   (real(ani(i)),i=1,4),real(gf),real(ap),real(ar),real(al)
c         write(8,*) char(10)
         kfile=kfile+1
   40 continue
      go to 60
 9000 continue
c      call rayplt
 9500 continue
      count=1.0
      go to 130
 9501 continue
c      xtimes=secnds(xtime0)
      xtimes = dtime(xtime0)
      write (*,9560) xtimes
      write(7,9560) xtimes
 9560 format (///,' total execution time for this set of rays is ',
     1   f10.3,' sec')
      stop
c    end of main program main main main main main main main main
      end
c
c ----------------------------------------------------------------------
      subroutine adams
c     adams predictor-corrector ne.  adams adams adams
	implicit real*8 (a-h)
	implicit real*8 (o-z)
      real*8 l,l0,lk,latitu,mucos,latmax,latmin
      real*8 hl2n,hl2s,hu2n,hu2s
      real xtime0,xtime1

      common /resona/ nreson, ndecre, numres
      common /x/ latitu,latmax,latmin, xtime0(2), xtime1(2), nsuppr
      common a,absb,al,alpha(4),alpha0(4),am,an,an1,an2,ane0,ani(4),
     1 ap,apl,apr,aps,ar,arl,as,b,c,cosz2,cosz22,cpsi,cpsi2,dd(10),ddk,
     2 def(10),delt,dlgndr(4),dlgndt(4),egfeq,expk,f0e,fkc,gf,h,hmin,
     3 l,l0(10),lk,pi,print,psi,psires,r0,radgra,rbase,rconsn,rdiv,relb,
     4 rstop,rzero,scbot,scr,sidedu(10),sinz2,sinz22,spsi,spsi2,tgfina,
     5 therm,tiniti,x(4),x0(5),y(4),y01,y02,z(5),clight,consf0,
     6 hl2n(10),hl2s(10),hu2n(10),hu2s(10),rducln(10),rducls(10),
     7 rducun(10),rducus(10),durbe,refalt,hmax,
     8 init,kducts,kinit,kj,kount,kr,kskip,kt,ktape,lines,mode,num,
     9 numray

      dimension ae(5),ak(5,4),e(5),f(5,5),p(5,5),xp(5,5),
     1 walpha(3),wan(3),wbeta(3),wf0e(3),wgama(3),wgf(3),
     2 wpsi(3),wpsire(3),wray(3),wsin22(3)
c     set up initial values for each ray. h is step size of integration
c
      n=5
      istop=0
      kt=0
      bound=-1.e6
      nspeci=0
      nreson=0
      oldlat=0.
      olddis=0.
   70 if (kt .eq. 1 .or. ndecre .eq. 1) h=h/4.0
      if (h .ge. hmin) goto 71
	 write(7,601) hmin
	 durbe=1.0
	 goto 1350
   71 if (numres .eq. 0 .or. kt .eq. 1) goto 78
      if (ndecre .eq. 0) goto 78
         oldh=h*4.
	 htt=h/clight         
         no584=1
         if (nreson .gt. numres) goto 1325
         ndecre=0
   72    if (lines .ge. 60) call title (lines)
         if (kount .eq. 0 .and. nsuppr .eq. 1) goto 78
         write(7,74) oldh, h, htt
   74    format (' action: current h =', f13.6,' is reduced ',
     1   'to h =',f13.6,',   dt=h/c=',f13.10,' sec'/)
         lines=lines+2
   78 do 90 i=1,n
   80     p(i,1)=x0(i)
   90 continue
  130 t=tiniti*clight
      j=1
      call funct(t,p,f,j,irtrn)
      if (nspeci .ne. 1) goto 131
c	nspeci = -1 : latitude adjustment has just been done
	nspeci = -1
	h = oldh
        if (h .gt. hmax) h = hmax
	htt = h/clight
	write(7,586) h, htt
	lines=lines+2
  131 if (ndecre .eq. 1) goto 70
      if(irtrn-1) 132,70,1400
  132 z2=p(2,j)*radgra
      latitu=90.00-z2
      delta=datan2(z(4),z(3)) *radgra
      alt=p(1,j)-r0
      ele=p(1,j)/(r0*sinz22)
      fgf=fkc/gf
      psirg=psires*radgra
      psig=psi*radgra
      sumq=(alpha(2)+alpha(3)/4.00+alpha(4)/16.00 )/am
      flhr=sqrt(sumq/(1.0/(f0e*f0e)+1.0/(gf*gf)))
      rdth=z(1)*f(2,j)
      delray=radgra*datan2(rdth,f(1,j))
      psiray=delray-delta+psig
      mucos=an*abs(cpsi)
      vloc=(1.0-fgf)/(fgf*mucos)
c     enerpb = kinetic energy (kevolt) of motion // b-field
      enerpb = 2.8428464e-6*(vloc*clight)**2
      distre=p(1,j)/r0
      tt=t/clight
      fgfeq = fkc/(egfeq/(ele*ele*ele))
      intan=(psiray+360.)/90.
      iquad=1+mod(intan,4)
c     psiray lies in quadrant iquad=1 if 0 .lt. psiray .lt. 90
      ihemi=dsign(1.d0,latitu)
      if (psiray .lt. -180.) psiray = psiray + 360.
      if (kount .ne. 0) goto 136
         if (nsuppr .eq. 0) goto 137
         if (latitu .lt. latmin .or. latitu .gt. latmax) goto 139
         call mskip (latitu,oldlat,olddis,lines)
         goto 137
  136 if (t .lt. bound) go to 140
  137       if (lines .ge. 60) call title(lines)
	    write(7,*)'t=',tt
  138       write(7,1187) p(5,j),distre,alt,latitu,ele,fgfeq,
     1        fgf,ani(1), flhr,an,mucos,delta,psig,psiray,psirg,enerpb
            lines=lines+1
  139       bound=t+print-1.e-6
            if (ktape .eq. 0) go to 140
c234567890c234567890c234567890c234567890c234567890c234567890c23456789012
            write(8,88) real(p(5,j)),real(distre),real(latitu),
     +        real(delta),real(tt),real(ele),real(psig),
     +        real(psiray),real(psirg),real(an),(real(ani(i)),i=1,4),
     +        real(gf),real(ap),real(ar),real(al)
 88            format(18(1x,e15.8))
c            write(8,*) char(10)
c
c this writes TG, DRE, LAT, DELTA, TP, ELE, PSIG, PSIRAY, PSIRG, AN, ANI(I)
c tg = group time
c dre = distance in earth radii
c lat = latitude
c delta = accuracy check
c tp = phase time
c ele = L shell
c psig = angle of k-vector (of group?)
c psiray = angle of ray
c psirg = angle of resonance cone for this freq
c an = refractive index
c ani(1) = e-  density /cc
c ani(2) = H+  density
c ani(3) = He+ density
c ani(4) = O+  density

  140 if (kskip .eq. 0) go to 180
  150    ia=2
  160    ib=4
  170    go to 300
  180 reltes =14.200*relb
  190 abstes =14.200*absb
  200 factor=relb/absb
  210 db=reltes /200.00
  220 h=2.0*h
c.....runge=kutta starting method.
  270    ia=2
  280   ib=2
  300 do 510 j=ia,ib

  310    call funct(t,p,f,j-1,irtrn)

	 if (numres .eq. 0 .and. durbe.gt.0.0) go to 1400
         if (nreson .gt. numres) goto 1325
         if (ndecre .eq. 1) goto 70
         if (irtrn-1) 312,70,1400
  312    if (j .lt. 3) go to 320
  311       walpha(j-1)=alpha(2)
            wbeta(j-1)=alpha(3)
            wgama(j-1)=alpha(4)
            wan(j-1)=an
            wpsi(j-1)=psi
            wpsire(j-1)=psires
            wf0e(j-1)=f0e
            wgf(j-1)=gf
            wsin22(j-1)= sinz22
            rdth=z(1)*f(2,j-1)
            wray(j-1)=radgra*datan2(rdth,f(1,j-1))
  320    do 350 i=1,n
  330       ak(i,1)=h*f(i,j-1)
  340       p(i,j)=p(i,j-1)+ak(i,1)/2.
  350    continue
  360    ttemp=t+h/2.

  370    call funct(ttemp,p,f,j,irtrn)

	 if (numres .eq. 0 .and. durbe .gt. 0.0) go to 1400
         if (nreson .gt. numres) goto 1325
         if (ndecre .eq. 1) goto 70
         if (irtrn-1) 380,70,1400
  380    do 410 i=1,n
  390       ak(i,2)=h*f(i,j)
  400       p(i,j)=p(i,j-1)+ak(i,2)/2.
  410    continue

  411    call funct(ttemp,p,f,j,irtrn)

     	 if (numres .eq. 0 .and. durbe.gt.0.00) go to 1400
         if (nreson .gt. numres) goto 1325
         if (ndecre .eq. 1) goto 70
         if (irtrn-1) 420,70,1400
  420    do 450 i=1,n
  430      ak(i,3)=h*f(i,j)
  440      p(i,j)=p(i,j-1)+ak(i,3)
  450    continue
  460    t=t+h

  470    call funct(t,p,f,j,irtrn)

	 if (numres .eq. 0 .and. durbe .gt. 0.) goto 1400  
   	 if (nreson .gt. numres) goto 1325
         if (ndecre .eq. 1) goto 70
         if (irtrn-1) 480,70,1400
  480    do 505 i=1,n
  490       ak(i,4)=h*f(i,j)
  500       p(i,j)=p(i,j-1)+(ak(i,1)+2.*ak(i,2)+2.*ak(i,3)+ak(i,4))/6.
  505    continue
  510 continue 
c
c....... end runge-kutta method.
c
  520 if (ib .ne. 2) go to 640
  530    do 550 i=1,n
  540       xp(i,5)=p(i,2)
  550    continue
c        xp(i) are the values with double interval for error analysis
  560    t=t-h
  570    h=h/2.
  580    if (kount .lt. 0)go to 590
  582       if (no584 .eq. 1) goto 590
               no584=0
               htt=h/clight
               if (kount .eq. 0 .and. nsuppr .eq. 1) goto 590
                  if (lines .ge. 60) call title(lines)
  584             write(7,586) h,htt
  586             format(  ' in the following calculation h=',f13.6,
     1            ',   dt=h/c=',f13.10,' sec'/)
                  lines=lines+2
  590    if (h.ge.hmin) go to 620
  600       write(7,601) hmin
  601       format(' h less than hmin =',f13.6/, ' to continue, '
     1           ,'change hmin or relax error conditions absb, relb')
     		print *, 'h smaller than hmin'
     		print *, h,hmin
            istop=1
  620    ib=3
  630    go to 300
  640 if (ib .ne. 3) go to 810
c        is accuracy criterion met
  650    j=3
  660    do 661 i=1,5
            ae(i)=xp(i,5)-p(i,j)
  661    continue
         do 760 i=1,5
  670       e(i)=abs(ae(i))
  680       if (e(i) .ge. abs(p(i,j))*reltes) go to 690
  681          e(i)=e(i)/abs(p(i,j))
  682          go to 760
  690       if (e(i) .ge. abstes) go to 700
  691          e(i)=e(i)*factor
  692          go to 760
  700       t=t-h
            if (kount .eq. 0 .and. nsuppr .eq. 1) goto 710
               if (kount .lt. 0) go to 710
                  if (lines .ge. 60) call title(lines)
                  write(7,701) i,(ae(ij),ij=1,5)
  701             format(   ' i=',i2,'  absolute errors '/,5f20.10)
                  lines=lines+2
  710       if (j .ne. 5) go to 530
  720       t=t-h
            do 740 iii=1,n
  730          p(iii,1)=p(iii,3)
  740       continue
  750       go to 270
  760    continue
  770    if (j .eq. 5) go to 1080
  780    ia=4
  790    ib=4
  800    go to 300
  810 t=t-3.00*h
  820 do 921 j=2,3
  830    t=t+h
c        print-out for runge-kutta
  862    z2=p(2,j)*radgra
         temlat = 90.00-z2
	 if (nspeci .eq. -1) goto 863
	    ratlat = (spelat-oldlat)/(temlat-oldlat)
	    if (ratlat .lt. .004 .or. ratlat .gt. .996) goto 1172
	       oldh = h
	       h = h*ratlat
               if (h .gt. hmax) h = hmax
	       htt = h/clight
	       nspeci = 1
	       goto 72
  863    latitu=temlat
         delta=datan2(p(4,j),p(3,j))*radgra
  871    alt=p(1,j)-r0
         dense=wf0e(j)*wf0e(j)/consf0
  872    ele=p(1,j)/(r0*wsin22(j))
         fgf=fkc/wgf(j)
  875    psirg=wpsire(j)*radgra
         psig=wpsi(j)*radgra
         sumq=(walpha(j)+wbeta(j)/4.00+wgama(j)/16.00 )/am
         flhr=sqrt(sumq/(1./(wf0e(j)*wf0e(j)) + 1./(wgf(j)*wgf(j))))
         gf=wgf(j)
         an=wan(j)
         delray=wray(j)
         psiray=delray-delta+psig
         cospsi=abs(dcos(wpsi(j)))
         mucos=an*cospsi
         vloc=(1.0-fgf)/(fgf*mucos)
         enerpb = 2.8428464e-6*(vloc*clight)**2
         distre=p(1,j)/r0
         tt=t/clight
         fgfeq = fkc/(egfeq/(ele*ele*ele))
         iquada=iquad
         iquad=1+mod(int((psiray+360.00)/90.00),4)
         ihemia=ihemi
         ihemi=dsign(1.d0,latitu)
	 if (psiray .lt. -180.) psiray=psiray+360.
	     if (tgfina .lt. 0.0) then
	     	if (ele .ge. (-tgfina)) then
                   write(6,801) -tgfina,ele,latitu
 801               format(1x,'Lstop=',f5.2,' L=',f5.2,' lat=',f7.2)
	     		istop = 1
	     	endif
	     else
            if (p(1,j).lt.rstop .or. p(5,j).gt.tgfina) then
               write(6,802) tgfina,p(5,j),rstop,p(1,j)
 802           format(6x,'Tgstop',f5.2,' Tg=',f5.2,'  Rlostop=',f8.1,
     1              ' R=',f8.1) 
            	istop = 1
            endif
         endif
         if (kount .ne. 0) goto 878
            if (nsuppr .eq. 0) goto 880
            if (latitu .lt. latmin .or. latitu .gt. latmax) goto 885
            call mskip (latitu,oldlat,olddis,lines)
            goto 880
  878    if (iquad .ne. iquada .or. ihemi.ne.ihemia
     1         .or. istop.eq.1)goto 880
               if (t .lt. bound)  go to 921
                  if (kount .lt. 0) go to 885
  880    if (lines .ge. 60) call title(lines)
	    write(7,*)'t=',tt
  882       write(7,1187) p(5,j),distre,alt,latitu,ele,fgfeq,
     1        fgf,dense,flhr,an,mucos,delta,psig,psiray,psirg,enerpb
         lines=lines+1
  885    bound=t+print-1.e-6
         ndecre=0
         if (ktape .eq. 0) go to 890
c234567890c234567890c234567890c234567890c234567890c234567890c23456789012
         write(8,88) real(p(5,j)),real(distre),real(latitu),
     +     real(delta),real(tt),real(ele),real(psig),real(psiray),
     +     real(psirg),real(an),(real(ani(i)),i=1,4),real(gf),
     +     real(ap),real(ar),real(al)
c         write(8,*) char(10)
  890    if (istop .eq. 1) go to 1350
c
c......  the ray specularly reflects when alt.lt.refalt
         if (specul.gt.0.0 .and. alt.gt.refalt) specul=0.0
         if (specul.gt.0.0) go to 53
         if (alt .gt. refalt) go to 51
             write(7,52) alt
   52        format(/23h specular reflection at, f9.2,12h km altitude)
             tiniti=tt
             x0(5)=p(5,j)
             x0(1)=distre*r0
             delt=180.0-delta
             if (delt .gt. 180.0) delt=delt-360.0
             durbe=-1.0
             kj=3
             specul=1.0
             go to 1400
c
   51    continue
  921 continue
c
      t=t+h
  922 konly=0
c     begin adams method.

  930 call funct(t,p,f,4,irtrn)

      if (numres .eq. 0 .and. durbe.gt.0.0) go to 1400
         if (nreson .gt. numres) goto 1325
         if (ndecre .eq. 1) goto 70
         if (irtrn-1) 1170,70,1400
c           print-out for adams
 1170    j=4
         z2=p(2,j)*radgra
         temlat = 90.00-z2
	 nspeci = 0
	 ratlat = (spelat-oldlat)/(temlat-oldlat)
	 if (ratlat .lt. .004 .or. ratlat .gt. .996) goto 1172
	    oldh = h
	    h = h*ratlat
            if (h .gt. hmax) hmax = h
            htt = h/clight
	    nspeci = 1
	    goto 72
 1172    latitu=temlat
         delta=datan2(z(4),z(3)) *radgra
         alt=p(1,j)-r0
         ele=p(1,j)/(r0*sinz22)
         fgf =fkc/gf
 1175    psig=psi*radgra
 1176    psirg=psires*radgra
 1177    sumq=(alpha(2)+alpha(3)/4.00+alpha(4)/16.00 )/am
 1178    flhr=sqrt(sumq/(1.0/(f0e*f0e)+1.0/(gf*gf)))
         rdth=z(1)*f(2,j)
         delray=radgra*datan2(rdth,f(1,j))
         psiray=delray-delta+psig
         mucos=an*abs(cpsi)
         vloc=(1.0-fgf)/(fgf*mucos)
         enerpb = 2.8428464e-12*(vloc*clight)**2
         distre=p(1,j)/r0
         tt=t/clight
         fgfeq = fkc/(egfeq/(ele*ele*ele))
         iquada=iquad
         iquad=1+mod(int((psiray+360.00)/90.00),4)
         ihemia=ihemi
         ihemi=dsign(1.d0,latitu)
	     if (psiray .lt. -180.) psiray = psiray + 360.
	     if (tgfina .lt. 0.0) then
	     	if (ele .ge. (-tgfina)) then
                   write(6,801) -tgfina,ele,latitu
                   istop = 1
	     	endif
	     else
            if (p(1,j).lt.rstop .or. p(5,j).gt.tgfina) then
               write(6,802) tgfina,p(5,j),rstop,p(1,j)
            	istop = 1
            endif
         endif
         if (kount .ne. 0) goto 3500
            if (nsuppr .eq. 0) goto 1179
            if (latitu .lt. latmin .or. latitu .gt. latmax) goto 1188
            call mskip (latitu,oldlat,olddis,lines)
            goto 1179
 3500    if(iquad.ne.iquada .or. ihemi.ne.ihemia 
     1      .or. istop.eq.1)goto 1179
            if (t .lt. bound)  go to 1190
            if (kount .lt. 0)go to 1188
 1179    if (lines .ge. 60) call title(lines)
	    write(7,*)'t=',tt
 1186        write(7,1187) p(5,j),distre,alt,latitu,ele,fgfeq,
     1      fgf,ani(1),flhr,an,mucos,delta,psig,psiray,psirg,enerpb
 1187    format(1x,f7.3,7g10.3)
         lines=lines+1
 1188    bound=t+print-1.e-6
         if (ktape .eq. 0) go to 1189
c234567890c234567890c234567890c234567890c234567890c234567890c23456789012
         write(8,88) real(p(5,j)),real(distre),real(latitu),
     +     real(delta),real(tt),real(ele),real(psig),real(psiray),
     +     real(psirg),real(an),(real(ani(i)),i=1,4),real(gf),
     +     real(ap),real(ar),real(al)
c         write(8,*) char(10)

 1181       format(4f10.4)
 1189    if (istop .eq. 1) go to 1350
c
c        the ray specularly reflects when alt.lt.refalt
         if (specul.gt.0.0 .and. alt.gt.refalt) specul=0.0
         if (specul.gt.0.0) go to 53
         if (alt .gt. refalt) go to 53
         write(7,54)alt
   54    format(/23h specular reflection at,f9.2,12h km altitude)
         tiniti=tt
         x0(5)=p(5,j)
         x0(1)=distre*r0
         delt=180.0-delta
         if (delt .gt. 180.0) delt=delt-360.0
         durbe=-1.0
         kj=3
         specul=1.0
         go to 1400
c
   53    continue
 1190    do 1180 i=1,5
 1180    x0(i)=p(i,3)
         tiniti=(t-h)/clight
         if (konly .eq. 1) go to 1250
c        xp(i,5) is predictor value, p(i,5) is corrector value
  940    do 960 i=1,n
  950       xp(i,5)=p(i,4)+h*(55.00*f(i,4)-59.00*f(i,3)
     1              +37.00*f(i,2)-9.00 *f(i,1))/24.00
  960    continue
  970    t=t+h
  980    call funct(t,xp,f,5,irtrn)
         if (numres .eq. 0 .and. durbe.gt.0.0) go to 1400
         if (nreson .gt. numres) goto 1325
         if (ndecre .eq. 1) goto 70
         if (irtrn-1) 990,70,1400
  990    do 1010 i=1,n
 1000       p(i,5)=p(i,4)+h*(9.00*f(i,5)+19.00*f(i,4)
     1             -5.00*f(i,3)+f(i,2))/24.00
 1010    continue
 1020    if (kskip .ne. 0) go to 1080
 1030       j=5
 1040       go to 660
 1080    do 1120 i=1,n
            p(i,3)=p(i,4)
 1090       p(i,4)=p(i,5)
 1100       do 1120 j=2,5
 1110          f(i,j-1)=f(i,j)
 1120    continue
 1210    if (kskip .gt. 0) go to 922
c           test whether the interval can be doubled.
 1220    do 1240 i=1,5
 1230       if (e(i) .gt. db) go to 922
 1240    continue
         konly=1
      go to 930
 1250 konly=0
      do 1270 i=1,n
 1260    p(i,1)=p(i,4)
 1270 continue
 1280 h=4.00*h
      if (h .gt. hmax) h = hmax
      if (kount .eq. 0 .and. nsuppr .eq. 1) goto 1300
         if (kount .lt. 0) go to 1300
            if (lines .ge. 60) call title(lines)
            write(7,1291) h
 1291       format(' double h (current h = ',f13.6,' )')
            lines=lines+1
 1300 go to 270
c
 1325 write(7,1330) numres
 1330 format (' number of step size reductions due to resonance ',
     1  'exceeds allowed maximum =',i3,/)
 1350 write(7,1351) p(1,j),latitu,delta,tt,p(5,j)
 1351 format('0end of raytracing.'/' r,lat,delt,t,tg ',5f15.8,/)
c      xtime=secnds(xtime1)
c      xtime2=secnds(xtime0)
      xtime = dtime(xtime1)
      xtime2 = dtime(xtime0)
      write(7,1360) xtime,xtime2
 1360 format (' execution time for this ray    = ',f10.3,' sec',/,
     1        ' total execution time up to now = ',f10.3,' sec')
      kj=3
 1400 return
c    end adams end adams end adams end adams end adams
      end
c
      subroutine title (lines)
	implicit real*8 (a-h)
	implicit real*8 (o-z)
      write(7,20)
 20   format('1   t     tg   dist(re)  alt(km)   lat     l     f/gfeq
     1 f/gf     dens     lhr    mu   mucos  delta   psi   psiray ',
     2 'psires    e//(kev)'/)
      lines=lines-58
      return
      end
c
      subroutine mskip (xlatit,oldlat,olddis,lines)
	implicit real*8 (a-h)
	implicit real*8 (o-z)
      distan=xlatit-oldlat
      if (distan*olddis .gt. 1.e-6) goto 20 
      write(7,10)
  10  format (/)
      lines=lines+1
  20  oldlat=latitu
      olddis=distan
      return
      end
c
      subroutine funct(t,p,dzdt,j,irtrn)
c     differential equations to be solved by adams.   funct funct funct
c     dz(i)/dt  z(1)= geocentric distance, z(2)= colatitude,
c               z(3)= n cos delta, z(4)= n sin delta, z(5)= group time
	implicit real*8 (a-h)
	implicit real*8 (o-z)
      real*8 l,l0,lk
      real*8 hl2n,hl2s,hu2n,hu2s
      real xtime0,xtime1

      common /resona/ nreson, ndecre, numres
      common /x/ xlatit, xlatma, xlatmi, xtime0(2), xtime1(2), nsuppr
      common a,absb,al,alpha(4),alpha0(4),am,an,an1,an2,ane0,ani(4),
     1 ap,apl,apr,aps,ar,arl,as,b,c,cosz2,cosz22,cpsi,cpsi2,dd(10),ddk,
     2 def(10),delt,dlgndr(4),dlgndt(4),egfeq,expk,f0e,fkc,gf,h,hmin,
     3 l,l0(10),lk,pi,print,psi,psires,r0,radgra,rbase,rconsn,rdiv,relb,
     4 rstop,rzero,scbot,scr,sidedu(10),sinz2,sinz22,spsi,spsi2,tgfina,
     5 therm,tiniti,x(4),x0(5),y(4),y01,y02,z(5),clight,consf0,
     6 hl2n(10),hl2s(10),hu2n(10),hu2s(10),rducln(10),rducls(10),
     7 rducun(10),rducus(10),durbe,refalt,hmax,
     8 init,kducts,kinit,kj,kount,kr,kskip,kt,ktape,lines,mode,num,
     9 numray

      dimension dzdt(5,5),p(5,5)
 
      kt=0
      irtrn=0
c.... if p(1,j) < rdiv, go back to adams to reduce step size by setting
c     irtrn = 1 (to goto line 70, adams), kt = 1 (to reduce h by 1/4).
      if (p(1,j) .gt. rdiv) goto 33
	 kt=1
	 irtrn=1
	 return
   33 do 35 i=1,5
   34 z(i)=p(i,j)
   35 continue
   40 call mu
      if (ndecre .eq. 1) return
      if (numres .eq. 0 .and. durbe .gt. 0.0) go to 1240
c     upon returning from mu, go back to calling program if
c     ndecre=1 to decrease step size
      if (kt .ne. 1) goto 50
      irtrn=1
      return
   50 if (z(1) .lt. rdiv) go to 140
   60    if (kr .ne. 1) go to 130
   70       kr=2
   80       tiniti  =t/clight
   90       do 110 i=1,5
  100          x0(i)=z(i)
  110       continue
  120       kj=1
  121       irtrn=2
            return
  130    if (z(1) .gt. rdiv) go to 200
  140 if (kr .ne. 2) go to 200
  150    tiniti  =t/clight
  160    do 180 i=1,5
  170       x0(i)=z(i)
  180    continue
  190    kj=2
  191    irtrn=2
         return
  200 ratio=an/an1
      if (fkc.gt.gf/am) go to 210
      if (mode.eq.+1) go to 221
  210    arg=-ap/as
  220    if (arg .ge. 0.0) go to 230
  221       psires=pi/2.00  * dsign(1.d0,psi)
  222       go to 400
c
  230 psires=datan(sqrt(arg))
  240 apsi=abs(psi)
  250 snpsi=psi/apsi
  260 if (apsi .ge. pi/2.00  ) go to 290
  270    psires=snpsi*psires
  280    go to 300
  290 psires=snpsi*(pi-psires)
  300 if (abs(psires-psi) .le. pi/36.00) go to 320
  310 if (abs(apsi-pi/2.0) .gt. pi/36.0) go to 400
  320    an12=an1*an1
  330    an14=an12*an12
  340    tpsi12=-(ap*an14-2.0*aps*an12+c)/(as*an14-(arl+aps)*an12+c)
  350    if (tpsi12 .ge. 0.0) go to 440
  360       if (an2 .gt. 0.0) go to 400
  370          write(7,371)
  371          format(' no solution for n ',/)
  380          kj=3
  390          irtrn=2
               return
c              correction of refractive index vector.
  400 do 420 i=3,4
  410    z(i)=z(i)*ratio
  411    p(i,j)=z(i)
  420 continue
  430 go to 750
c
  440 if (abs(ratio-1.0) .le. 1.e-07) go to 750
  450    psi1=datan(sqrt(tpsi12))
  460    if (apsi .ge. pi/2.00  ) go to 490
  470       psi1=snpsi*psi1
  480       go to 500
  490    psi1=snpsi*(pi-psi1)
  500    delpsi=psi1-psi
  510    if (an2 .gt. 0.0) go to 540
  520       psi=psi1
  530       go to 700
  540    delan=an-an1
  550    x1=an1*delpsi
  560    tana=x1/delan
         tana2=tana*tana
  570    cosa2=1.0/(1.0+tana2)
  580    sina2= tana2 *cosa2
  590    x2=x1*cosa2
  600    x4=delan*sina2
  610    dpsi3=x2/an1
  620    dan2=x4*delan
  630    if (dan2 .le. an1*an1*1.e-05) go to 680
            if (kount .eq. 0 .and. nsuppr .eq. 1) goto 650
               if (kount .lt. 0) go to 650
                  if (lines .ge. 60) call title(lines)
  640             write(7,641) h
  641             format(' correction vector is big. step h = ',f13.6)
                  lines=lines+1
  650       if (kskip .ne. 0) go to 680
               kt=1
  670          irtrn=1
               return
  680    an1=an1+x4
  690    psi=psi+dpsi3
  700    cpsi=dcos(psi)
         cpsi2=cpsi*cpsi
  710    spsi=dsin(psi)
         spsi2=spsi*spsi
  720    an=an1
  721    an2=an1*an1
  730    z(3)=an1*(y01*cpsi-y02*spsi)
  740    z(4)=an1*(y02*cpsi+y01*spsi)
         p(3,j)=z(3)
         p(4,j)=z(4)
  750 s2psi=2.0*spsi*cpsi
  760 dadpsi=s2psi*(as-ap)
  770 dbdpsi=s2psi*(arl-aps)
c     dcdpsi=0.0
  771 an1=an
  780 an12=an*an
  790 an13=an1*an12
  800 an14=an12*an12
  810 denom1=1.0/(an12*spsi)
  820 dpsidg  =(z(3)*y02-z(4)*y01)*denom1*an1
      denoms=4.0*an13*a-2.0*an1*b
      if (abs(denoms) .ge. 1.e-35) goto 830
         irtrn=1
         kt=1
         return
  830 denom2=1.0/denoms
  840 dndpsi=-(an14*dadpsi-an12*dbdpsi)*denom2
  850 dpsidz=(cpsi*z(3)-an1*y01)*denom1
  860 dandz3=dndpsi*dpsidz
  870 dpsidz=(cpsi*z(4)-an1*y02)*denom1
  880 dandz4=dndpsi*dpsidz
c     dgdz1=0.0, dpsdz1=0.0
  890 dgdz2  =1.0/(2.00*cosz22+sinz22/2.00)
  910 dpsdz2=dpsidg*dgdz2
  920 dandz1=0.0
  930 dandz2=dndpsi*dpsdz2
  931 dardf=0.0
  932 daldf=0.0
  933 dapdf=0.0
  970 dapdx=-1.0
  940 do 1120 i=1,num
  950    dardx=-1.0/(1.0+y(i))
  960    daldx=-1.0/(1.0-y(i))
  980    dadx=(dardx+daldx)*spsi2/2.00+dapdx*cpsi2
  990    dbdx=(ar*daldx+al*dardx)*spsi2+(ap*(dardx+daldx)/2.00+
     1        as*dapdx)*(1.0+cpsi2)
 1000    dcdx=apr*daldx+apl*dardx+arl*dapdx
 1010    dandx=(-dadx*an14+dbdx*an12-dcdx)*denom2
         ymns=1.0-y(i)
         ypls=1.0+y(i)
 1020    dardy= x(i)/(ypls*ypls)
 1030    daldy=-x(i)/(ymns*ymns)
c        dapdy(i)=0.0
 1040    dady=(dardy+daldy)*spsi2/2.
 1050    dbdy=(ar*daldy+al*dardy)*spsi2+ap*(dardy+daldy)*(1.0+cpsi2)/2.
 1060    dcdy=apl*dardy+apr*daldy
 1061    dxdz1=dlgndr(i) *x(i)
c        to prevent underflow of dxdz2
         check=1.e18*dlgndt(i)*x(i)
         if (abs(check).gt.1.e-20) goto 1062
            dxdz2=1.e-38
            goto 1070
 1062    dxdz2=check*1.e-18
 1070    dandy=(-an14*dady+an12*dbdy-dcdy)*denom2
 1080    dydz1=-3.e0*y(i)/z(1)
 1090    dydz2=-3.e0*y(i)*sinz2*cosz2/(1.0+3.00*cosz22)
 1100    dandz1=dandz1+dandx*dxdz1+dandy*dydz1
 1110    dandz2=dandz2+dandx*dxdz2+dandy*dydz2
 1111    dardf=(2.0+y(i))*dardy+dardf
 1112    daldf=-(2.0-y(i))*daldy+daldf
 1113    dapdf=x(i)+dapdf
 1120 continue
 1124 dapdf=2.0*dapdf
 1125 dadf=(dardf+daldf)*spsi2/2.00+dapdf*cpsi2
 1126 dbdf=(al*dardf+ar*daldf)*spsi2+((ar+al)*dapdf/2.00+
     1      ap*(dardf+daldf)/2.00)*(1.0+cpsi2)
 1127 dcdf=arl*dapdf+apl*dardf+apr*daldf
 1128 dandf=(-an14*dadf+an12*dbdf-dcdf)*denom2
 1129 gn=an+dandf
c
 1130 dzdt(1,j)=(z(3)-an1*dandz3  )/an12
 1140 dzdt(2,j)=(z(4)-an1*dandz4  )/(an12*z(1))
 1150 dzdt(3,j)=dandz1  /an1+z(4)*dzdt(2,j)
 1160 dzdt(4,j)=(dandz2  /an1-z(4)*dzdt(1,j))/z(1)
 1170 dzdt(5,j)=gn /(an*clight)
 1240 return
c     end of funct subroutine.
      end
c
      subroutine mu
c     refractive index mu mu mu mu mu mu mu mu mu  mu mu
	implicit real*8 (a-h)
	implicit real*8 (o-z)
      real*8 l,l0,lk
      real*8 hl2n,hl2s,hu2n,hu2s

      common /resona/ nreson, ndecre, numres
      common a,absb,al,alpha(4),alpha0(4),am,an,an1,an2,ane0,ani(4),
     1 ap,apl,apr,aps,ar,arl,as,b,c,cosz2,cosz22,cpsi,cpsi2,dd(10),ddk,
     2 def(10),delt,dlgndr(4),dlgndt(4),egfeq,expk,f0e,fkc,gf,h,hmin,
     3 l,l0(10),lk,pi,print,psi,psires,r0,radgra,rbase,rconsn,rdiv,relb,
     4 rstop,rzero,scbot,scr,sidedu(10),sinz2,sinz22,spsi,spsi2,tgfina,
     5 therm,tiniti,x(4),x0(5),y(4),y01,y02,z(5),clight,consf0,
     6 hl2n(10),hl2s(10),hu2n(10),hu2s(10),rducln(10),rducls(10),
     7 rducun(10),rducus(10),durbe,refalt,hmax,
     8 init,kducts,kinit,kj,kount,kr,kskip,kt,ktape,lines,mode,num,
     9 numray
c
      init=init+1
      ndecre=0
   40 call dens
   50 f0e=sqrt(consf0 *ani(1))
   60 f0efkc=f0e/fkc
      x(1)=f0efkc*f0efkc
   70 cons=x(1)/am
   80 x(2)=cons*alpha(2)
   90 x(3)=cons*alpha(3)/4.
  100 x(4)=cons*alpha(4)/16.
      do 109 i=2,4
         if (x(i).eq. 0.0) go to 109
            if (x(i) .gt. 1.0e-35) go to 109
               x(i)=1.0e-35
  109 continue
  150 r0z1=r0/z(1)
      gf=egfeq*r0z1*r0z1*r0z1*sqrt(1.0+3.00*cosz22)
  160 y(1)=-gf/fkc
  170 y(2)=-y(1)/am
  180 y(3)=y(2)/4.
  190 y(4)=y(3)/4.
  200 ar=1.0
  210 al=1.0
  220 ap=1.0
  230 do 270 i=1,num
  240    ar=ar-x(i)/(1.0+y(i))
  250    al=al-x(i)/(1.0-y(i))
  260    ap=ap-x(i)
  270 continue
  280 as=(ar+al)/2.
      ad=(ar-al)/2.
      if (ad.lt.0.0) mode=1
      if (ad.ge.0.0) mode=-1
  290 taz2=sinz2/cosz2
  300 y01=dsign(1.d0,taz2)/sqrt(1.0+taz2*taz2/4.00)
  310 y02=sqrt(1.0-y01*y01)
  320 if (init .gt. 1) go to 370
  330    cpsi=y01*dcos(delt)+y02*dsin(delt)
  340    spsi=y01*dsin(delt)-y02*dcos(delt)
  360    go to 400
  370       an1=sqrt(z(3)*z(3)+z(4)*z(4))
  380       cpsi=(y01*z(3)+y02*z(4))/an1
  390       spsi=(y01*z(4)-y02*z(3))/an1
  400 cpsi2=cpsi*cpsi
  410 spsi2=spsi*spsi
  420 arl=ar*al
  430 aps=ap*as
  431 apl=ap*al
  432 apr=ap*ar
  440 a=as*spsi2+ap*cpsi2
  450 b=arl*spsi2+aps*(1.0+cpsi2)
  460 c=ap*arl
  470 if (fkc .ge. gf/am ) go to 500
  471    if (mode .lt. 0) go to 500
  472       if (c .ge. 0.0) go to 480
               atemp=(arl-aps)*spsi2
  473          f2=atemp*atemp+4.00*(ap*ad)*(ap*ad)*cpsi2
  474          an2=(2.0*c)/(b-sqrt(f2))
  475          go to 550
  480       kt=1
            if (kount .lt. 0) go to 490
  481          write(7,482)
  482          format('   an2  was  negative  here'  )
               lines=lines+1
  490       return
  500 f2=b*b -4.0*a*c
  510 if (b .gt. 0.0) go to 540
  520    an2=(b-sqrt(f2))/(2.0*a)
  530    go to 550
  540 an2=(2.0*c)/(b+sqrt(f2))
  550 if (an2 .ge. 0.0) go to 570
         if (kount .lt. 0) go to 480
c
  560       if (lines .ge. 60) call title(lines)
            write(7,562)
  562       format(  ' psi outside resonance cone')
	    lines=lines+1
c           set ndecre=1 to go back to statement 70 in adams to
c           decrease step size by a factor 1/4
            durbe=1.0
            if (numres .eq. 0) goto 630
               ndecre=1
               nreson=nreson+1
               goto 630
  570 an=sqrt(abs(an2))
  580 if (init .gt. 1) go to 620
  590    x0(3)=an*dcos(delt)
  600    x0(4)=an*dsin(delt)
  601    an1=an
  620 psi=datan2(spsi,cpsi)     
  630 return
c    end mu end mu end mu end mu end mu end mu end mu end mu
      end
c
      subroutine dens
c    density density density density density density
	implicit real*8 (a-h)
	implicit real*8 (o-z)
      real*8 l,l0,lk,latitu,latmax,latmin
      real*8 hl2n,hl2s,hu2n,hu2s
      real xtime0,xtime1

      common /x/ latitu,latmax,latmin, xtime0(2), xtime1(2), nsuppr
      common a,absb,al,alpha(4),alpha0(4),am,an,an1,an2,ane0,ani(4),
     1 ap,apl,apr,aps,ar,arl,as,b,c,cosz2,cosz22,cpsi,cpsi2,dd(10),ddk,
     2 def(10),delt,dlgndr(4),dlgndt(4),egfeq,expk,f0e,fkc,gf,h,hmin,
     3 l,l0(10),lk,pi,print,psi,psires,r0,radgra,rbase,rconsn,rdiv,relb,
     4 rstop,rzero,scbot,scr,sidedu(10),sinz2,sinz22,spsi,spsi2,tgfina,
     5 therm,tiniti,x(4),x0(5),y(4),y01,y02,z(5),clight,consf0,
     6 hl2n(10),hl2s(10),hu2n(10),hu2s(10),rducln(10),rducls(10),
     7 rducun(10),rducus(10),durbe,refalt,hmax,
     8 init,kducts,kinit,kj,kount,kr,kskip,kt,ktape,lines,mode,num,
     9 numray

      dimension exnor(4),qi(4),sh(4)
c
      cosz2 = dcos(z(2))
      sinz2 = dsin(z(2))
      cosz22 = cosz2*cosz2
      sinz22 = sinz2*sinz2
   40 if (init .gt. 0) go to 110
   50    rb7370 = rbase/7370.
         sh(2) = 1.150600 * therm * rb7370*rb7370
   90    sh(3) = sh(2)/4.
  100    sh(4) = sh(3)/4.
         do 105 i = 2,4
  105    alpha(i) = 0.0
c           electron and ion densities.
  110 gph = rbase*(1.0-rbase/z(1))
  120 exnor(2) = dexp(-gph/sh(2))
  130 exnor(3) = exnor(2)*exnor(2)*exnor(2)*exnor(2)
  140 exnor(4) = exnor(3)*exnor(3)*exnor(3)*exnor(3)
  150 q = 0.0
  160 sumi = 0.0
  170 do 210 i = 2,num
  180    qi(i) = alpha0(i)*exnor(i)
  190    q = q+qi(i)
c        sum1 = -dqdz
  200    sumi = sumi+qi(i)/sh(i)
  210 continue
      do 211 i = 2,num
  211 alpha(i) = qi(i)/q
c     anr = ne / ne0      radial dependence only - de model
  220 anr = sqrt(q)
c     anli takes into account the lower ionosphere
  221 arg = (z(1)- rzero)/scbot
      if (arg.lt.13.0) go to 222
         arg = 13.0
  222 exarg = dexp(-arg*arg)
  223 anli = 1.0-exarg
  224 dlnlid  = arg*exarg*2./(scbot*anli)
  231 l = z(1)/(r0*sinz22)
      dlnldr = 0
      dlnldt = 0
      ani(1) = ane0*anr
      if (init .ne. 0) go to 233
         write(7,232)z(1),ani(1),(alpha(i),i = 2,4)
  232    format(' de model at r = ',f10.4,'    dens = ',e15.8,'   h+ = '
     1            ,e15.8,'   he+ = ',e15.8,'   o+ = ',e15.8)
         lines = lines+1
  233 ani(1) = ani(1)*anli
      if (kducts.eq.0)go to 250
c        plasmapause
         cotz2 = cosz2/sinz2
         deltal = l-lk
         if (deltal.lt.0.0) go to 410
            d2 = ddk*ddk
            argl = deltal*deltal/(d2*2.0)
            if (argl.lt.80.00) go to 400
               argl = 80.00
  400       f = dexp(-argl)
            trm = (rconsn/z(1))**expk
            argr = (z(1)-rconsn)/scr
            if (argr.lt.12.50) go to 405
               argr = 12.5
  405       fr = dexp(-argr*argr)
            trmodl = trm+(1.-trm)*fr
            dtrmdr = -expk*trm*(1.-fr)/z(1)-(1.-trm)*fr*2.*argr/scr
            anlk = f+trmodl*(1.0-f)
            factor = deltal*f*l*(1.0-trmodl)/(d2*anlk)
            dlnldr = dlnldr-factor/z(1) + (1.-f)*dtrmdr/anlk
            dlnldt = dlnldt+2.0*factor*cotz2
            ani(1) = ani(1)*anlk
  410    if (kducts.eq.1) go to 250
            if (l0(2) .gt. 0.0) go to 1000
c              sinusoidal density perturbation
  900          kinit = 3
               deltal = l+l0(2)
               if (deltal*sidedu(2) .ge. 0.0) go to 920
  910             deltal = 0
  920          delk = -l0(2)-(lk+ddk)+dd(2)/2
               critl = (lk+ddk)+mod(delk,dd(2))
c              *critl* is the first zero in the sinu pert beyond
c              the knee --  sinusoidal perturbation begins here
               if (l .le. critl) go to 990
                  argl = 2.0*pi*deltal/dd(2)
                  delnl = (def(2)/2.)*(1.+dcos(argl))
  930             if (latitu.le.0 .and. z(1).le.rducus(2)) go to 950
                  if (latitu.ge.0 .and. z(1).le.rducun(2)) go to 950
c                    continue only if point lies above rducun and rducus
                     if (latitu.ge.0) delr = z(1)-rducun(2)
                     if (latitu.le.0) delr = z(1)-rducus(2)
                     if (latitu.le.0) arglr = delr*delr/hu2s(2)
                     if (latitu.ge.0) arglr = delr*delr/hu2n(2)
                     if (arglr.ge.75.0) go to  990
                        frduct = dexp(-arglr)
                        if (latitu.le.0) delroh = 2.0*delr/hu2s(2)
                        if (latitu.ge.0) delroh = 2.0*delr/hu2n(2)
                        go to 960
  950             if (latitu.le.0 .and. z(1).ge.rducls(2)) go to 970
                  if (latitu.ge.0 .and. z(1).ge.rducln(2)) go to 970
c                     continue here only if point lies below 
c                     rducln and rducls
                      if (latitu.ge.0) delr = z(1)-rducln(2)
                      if (latitu.le.0) delr = z(1)-rducls(2)
                      if (latitu.le.0) arglr = delr*delr/hl2s(2)  
                      if (latitu.ge.0) arglr = delr*delr/hl2n(2)
                      if (arglr.ge.75.0) go to  990
                         frduct = dexp(-arglr)
                         if (latitu.ge.0) delroh = 2.0*delr/hl2n(2)
                         if (latitu.le.0) delroh = 2.0*delr/hl2s(2)
  960                    delnl = delnl*frduct
                         anl = 1.0+delnl
                         fac = pi*def(2)*l*dsin(argl)/(dd(2)*anl)
                         onedut = fac*2.0*cotz2
                         onedur = -fac/z(1) - delnl*delroh/anl
                         go to 980
  970              anl = 1.0+delnl
                   fac = pi*def(2)*l*dsin(argl)/(dd(2)*anl)
                   onedut = fac*2.0*cotz2
                   onedur = -fac/z(1)
  980              dlnldr = dlnldr+onedur
                   dlnldt = dlnldt+onedut
                   ani(1) = ani(1)*anl
  990              if (kducts.eq.2) goto 250
c              duct(s)
 1000       do 240 kduc = kinit,kducts
               deltal = l-l0(kduc)
               if (deltal*sidedu(kduc) .ge. 0.0) go to 1200
 1100             deltal = 0
 1200          d2 = dd(kduc)*dd(kduc)
               argl = deltal*deltal/(d2*2.0)
               if (argl.gt. 80.0 ) go to 240
  234             delnl = def(kduc)*dexp(-argl)
                  if (latitu.ge.0 .and. z(1).le.rducun(kduc)) goto 1500
                  if (latitu.le.0 .and. z(1).le.rducus(kduc)) goto 1500
c                    continue only if point lies above rducs and rducun
                     if (latitu.ge.0) delr = z(1)-rducun(kduc)
                     if (latitu.le.0) delr = z(1)-rducus(kduc)
                     if (latitu.ge.0) arglr = delr*delr/hu2n(kduc)
                     if (latitu.le.0) arglr = delr*delr/hu2s(kduc)
                     if (arglr.ge.75.0) go to 240
                        frduct = dexp(-arglr)
                        if (latitu.ge.0) delroh = 2.0*delr/hu2n(kduc)
                        if (latitu.le.0) delroh = 2.0*delr/hu2s(kduc)
                        go to 1600
 1500             if (latitu.ge.0 .and. z(1).ge.rducln(kduc)) go to 235
                  if (latitu.le.0 .and. z(1).ge.rducls(kduc)) go to 235
c                    continue here only if point lies
c                    below rducln and rducls
                     if (latitu.ge.0) delr = z(1)-rducln(kduc)
                     if (latitu.le.0) delr = z(1)-rducls(kduc)
                     if (latitu.ge.0) arglr = delr*delr/hl2n(kduc)
                     if (latitu.le.0) arglr = delr*delr/hl2s(kduc)
                     if (arglr.ge.75.0) go to 240
                        frduct = dexp(-arglr)
                        if (latitu.ge.0) delroh = 2.0*delr/hl2n(kduc)
                        if (latitu.le.0) delroh = 2.0*delr/hl2s(kduc)
 1600                   delnl = delnl*frduct
                        anl = 1.0+delnl
                        fac = delnl*deltal*l/(anl*d2)
                        onedut = fac*2.0*cotz2
                        onedur = -fac/z(1) - delnl*delroh/anl
                        go to 1700
  235             anl = 1.0+delnl
  236             onedur = -delnl*deltal*l/(anl*d2*z(1))
  237             onedut = -onedur*2.0*z(1)*cotz2
 1700             dlnldr = dlnldr+onedur
                  dlnldt = dlnldt+onedut
                  ani(1) = ani(1)*anl
  240       continue
  250 do 270 i = 2,num
  260    ani(i) = ani(1)*alpha(i)
  270 continue
c     vzs = dzdr
  290 vzs = (rbase/z(1))*(rbase/z(1))
c     dlnrdr  = d nr(de) dr / nr
  300 dlnrdr = -sumi*vzs/(2.0*q)
  310 do 330 i = 2,num
  320    dlgndr(i) = -dlnrdr-vzs/sh(i)+dlnldr +dlnlid
  321    dlgndt(i) = dlnldt
  330 continue
      dlgndr(1) = dlnrdr+dlnldr  +dlnlid
  341 dlgndt(1) = dlnldt
  350 return
c    end density end density end density end density
      end
c
c   ----- Ancient plotting routine ------------  
c   commenting out, since there's no need for it
c   --------------------------------------------
C       subroutine rayplt
C c       dd80 plotting routine for ncar 7700
C c this code retained since it would be realitively easy to activate the 
C c NCAR plot calls.  But we have Matlab, and who has the time [dsl 6/95].
C c
C 	implicit real*8 (a-h)
C 	implicit real*8 (o-z)
C       real*8 l0,lp(10),lat,latitu,l,ll,lk,gf,grarad
C       real*8 xdra(150),ydra(150), xnew,ynew,xc,xd,yc,yd,dt,xa,ya,y1
C       character  endplt*51, noplot*7
C c      dimension lab1(30),lab2(30),ifield(2000)
C       logical ymult,plus, grlab(51)
C       real*8 hl2n,hl2s,hu2n,hu2s

C       common /ryplt/ freq, altinj
C       common a,absb,al,alpha(4),alpha0(4),am,an,an1,an2,ane0,ani(4),
C      1 ap,apl,apr,aps,ar,arl,as,b,c,cosz2,cosz22,cpsi,cpsi2,dd(10),ddk,
C      2 def(10),delt,dlgndr(4),dlgndt(4),egfeq,expk,f0e,fkc,gf,h,hmin,
C      3 l,l0(10),lk,pi,print,psi,psires,r0,radgra,rbase,rconsn,rdiv,relb,
C      4 rstop,rzero,scbot,scr,sidedu(10),sinz2,sinz22,spsi,spsi2,tgfina,
C      5 therm,tiniti,x(4),x0(5),y(4),y01,y02,z(5),clight,consf0,
C      6 hl2n(10),hl2s(10),hu2n(10),hu2s(10),rducln(10),rducls(10),
C      7 rducun(10),rducus(10),durbe,refalt,hmax,
C      8 init,kducts,kinit,kj,kount,kr,kskip,kt,ktape,lines,mode,num,
C      9 numray

C       data ymult/1hs/,plus/1h+/
C c
C 	  print *, 'rayplt'
C       grarad = 3.14159/180.
C       write  (6,20)
C    20 format (1h1,26h data read by plot routine/)
C    30 format (51a1)
C    40 format (2x, 51a1)
C    50 format (a51)
C    55 format (a7)
C       ktold = ktape
C c     looping from 60 to 310 to process each plot. to end plotting, 
C c     the data card for grlab must contain the string end
C    60 rewind 8
C 	 read (4,30) grlab
C          backspace 4
C          read (4,50) endplt
C          if (endplt .eq. 'end') return
C          backspace 4
C          read (4,55) noplot
C 	    write(7,40) grlab
C 	    ktape  = ktold
C       kfile  = 1
C       read (4,*)   kplot,xscale,yscale,flbot,fltop,(lp(j),j = 1,5)
C       write(7,126) kplot,xscale,yscale,flbot,fltop,(lp(j),j = 1,5)
C   126 format (i5, 4f5.1, 5f5.2, i5)
C       if (noplot .eq. 'noplot') goto 60
C       if (kplot .le. 0) return
C       if (kplot .eq. 2) dis = .7      
C       if (kplot .le. 10) goto 128
C c	if kplot .ge. 11, then the wave normal will be plotted at every
C c	distance equal to dis = .1*(kplot-10).
C 	dis = .1*(kplot-10)
C 	kplot = 2
C  128  if (kplot .ge. 3) go to 200
C c
C c     meridional plane plot
C       xc = -5./yscale
C       xd = -xc
C       yc = 0.
C       yd = 10.0/yscale
C       if (kfile.gt.1) go to 163
C c_         call flash1 (ifield,10000)
C c         call set(.1,.9,.1,.9,xc,xd,yc,yd,1)
C          ltop = fltop
C          dt = .05/yscale
C c         call frstpt(-dt,1.)
C c         call vector(dt,1.)
C          lmt = ltop-1
C          do 111 i = 1,lmt
C             ynew = i
C             y1 = i+1
C c            call frstpt(0.,ynew)
C c            call vector(0.,y1)
C c            call frstpt(-dt,y1)
C c            call vector(dt,y1)
C   111    continue
C c         call frstpt(1.,0.)
C c         call vector(-1.,0.)
C          latitu = -90.
C          do 130 i = 1,19
C             lat = latitu*grarad
C             xnew = dsin(lat)
C             ynew = dcos(lat)
C c            call vector(xnew,ynew)
C   130    latitu = latitu+10.
C          j = 0
C   140    continue
C c  	   call frstpt(0.,1.)
C          j = j+1
C          if (lp(j) .eq. 0.0) go to 160
C             latitu = -75.
C             ii = -1
C             do 150 i = 1,151
C                lat = latitu*grarad
C                c = dcos(lat)
C                r = lp(j)*c*c
C                if (r.lt.1.0) go to 150
C                   xnew = r*dsin(lat)
C                   ynew = r*c
C c                  if (ii.lt.0) call frstpt(xnew,ynew)
C c                  if (ii.gt.0) call vector(xnew,ynew)
C                   ii = -ii
C   150       latitu = latitu+1.
C             go to 140
C c_  160    call flash2(1,len)
C  160            continue
C          write(7,162) len
C   162    format(' ',i10)
C c_         call flash3(1)
C c         encode(51,3332,lab1)grlab
C  3332    format(51a1)
C c         call pwrx(512,974,lab1,51,.8,0.,1)
C c         encode(66,3333,lab1)freq,altinj
C  3333 format(14h"pru"f"pru"  = ,f5.1,32h"pru" khz"pru"   ,   injected at
C      1,f7.1,8h"pru" km)
C c      call pwrx(512,50,lab1,66,.8,0.,1)
C c      encode(128,3334,lab2)alpha0(2),alpha0(3),alpha0(4),expk,lk
C  3334 format(14h"pru"de-model ,f4.1,10h h"s"+"n" ,f3.1,21h h"prl"e"pru""
C      1s"+"n" ,f3.1,39h o"s"+"n"   "prl"n"pgu","prl"r"s""pru"-,f2.0,30h"n
C      2""prl" for "pru"l"prl"9"pru",f2.0)
C c      call pwrx(512,10,lab2,128,.8,0.,1)
C c
C   163 continue
C c      call frstpt(0.,0.)
C       xp = 0.
C       yp = 0.
C       ii = 1
C       tgbf = 0.0
C       i = 0
C       count = 1.0
C   170 read(8) dummmm,distre,latitu,delta,zero,zero,zero,zero,zero,
C      +        zero,zero,zero,zero,zero
C   171    format(10x,4f10.4) 
C          if (xdra(1).gt.99998.0) go to 300
C c            if (dummmm.gt.99998.0 .and. kfile.eq.1)
C c     1      call points (xdra,ydra,i,ymult,1)
C             if (dummmm-99998.0) 201,201,300
C c
C   201          continue
C                if ( dummmm.gt.tgfina) go to 170
C                   if (distre.lt.1.0) go to 170
C                      if (dummmm .lt. tgbf) go to 170
C                         tgbf = dummmm
C                         lat = -latitu*grarad
C                         xnew = distre*dsin(lat)
C                         ynew = distre*dcos(lat)
C                         if (kfile.ne.1)go to 175
C                            i = i+1
C                            xdra(i) = xnew
C                            ydra(i) = ynew
C                            if (kfile.eq.1) go to 170
C   175                   continue
C                         if (ynew.gt.fltop .or. 
C      1                  abs(xnew).gt.abs(xc)) go to 163
C c                           if (ii.eq.1) call frstpt(xnew,ynew)
C                            ii = 0
C c                           call vector (xnew,ynew)
C                            if (kfile .eq. 1) go to 170
C                               if (kplot.ne.1) go to 172
C                                  time = count*0.5
C                                  if (dummmm.lt.time) go to 172
C c                                    call points (xnew,ynew,1,plus,1)
C                                     count = count+1.0
C   172                         continue
C                               if (kplot.eq.1) go to 170
C c                                wave normal angle
C                                  dx = xnew-xp
C        			         dy = ynew-yp
C 			         d = sqrt(dx*dx+dy*dy)
C c				 plot wave normal every distance dis:
C 		 	         if (d.lt.dis) go to 170
C 		                    wn = (-latitu+delta)*grarad
C 			            xa = xnew+0.2*dsin(wn)
C 			            ya = ynew+0.2*dcos(wn)
C c			            call frstpt (xa,ya)
C c			            call vector (xnew,ynew)
C 			            xp = xnew
C 			            yp = ynew
C 			            go to 170
C c.... l vs lat plot
C   200 xc = 80.0/xscale
C       xd = -xc
C       yc = flbot/yscale
C       yd = fltop/yscale
C c      call set (.1,.9,.1,.9,xc,xd,yc,yd,1)
C       idl = fltop-flbot
C c      encode(51,3332,lab1)grlab
C c      call pwrit (512,974,lab1,51,2,0,0)
C c      encode (29,3335,lab1)
C  3335 format(29hgeomagnetic latitude (deg)   )
C c      call pwrit(512,10,lab1,29,2,0,0)
C c      call labmod(6h(f4.0),6h(f4.0),4,4,1,1,0,0,0)
C c      call halfax(8,1,idl,1,0.,flbot,1,1)
C       j = 0
C c FRF - go to shouldn't go to the end of a do loop, not allowed by 
C c modern fortran compilers.  Changed to "go to 240".
C c      if (kplot.eq.4) go to 230
C       if (kplot.eq.4) go to 240
C          if (kplot.eq.3) go to 220
C             latitu = -75.
C             ii = -1
C             do 210 i = 1,151
C                c = dcos(latitu*grarad)
C                ll = 1./(c*c)
C                if (ll.gt.fltop .or. ll.lt.flbot) go to 205
C                   xnew = latitu
C                   ynew = ll-flbot+1
C c                  if (ii.lt.0) call frstpt(xnew,ynew)
C c                  if (ii.gt.0) call vector(xnew,ynew)
C                   ii = -ii
C                   ii = 0
C c                  call vector(xnew,ynew)
C                   go to 210
C   205          continue
C c               call frstpt(0.,0.)
C   210       latitu = latitu+1.
C   220       continue
C c            call frstpt(0.,0.)
C             j = j+1
C             if (lp(j) .eq. 0.0) go to 240
C                invlat = datan(sqrt(lp(j)-1.))/grarad
C                latitu = -invlat
C                ii = -1
C                do 230 i = 1,151
C                   xnew = latitu
C                   ynew = lp(j)-flbot+1
C c                  if (ii.lt.0) call frstpt(xnew,ynew)
C c                  if (ii.gt.0) call vector(xnew,ynew)
C                   ii = -ii
C                   latitu = latitu+1.
C                   if(latitu.ge.invlat) go to 220
C   230 continue
C c
C   240 continue
C c      call frstpt(0.,0.)
C       ii = 1
C       tgbf = 0.0
C   250 read(8) dummmm,distre,latitu,delta,zero,zero,zero,zero,zero,
C      +        zero,zero,zero,zero,zero
C   251    format(10x,4f10.4)
C          if (dummmm-99998.0) 202,202,300
C   202       continue
C             if (dummmm.gt.tgfina) go to 250
C                if (distre.lt.1.0) go to 250
C                   if (kfile.eq.1) go to 250
C                      if (dummmm.lt.tgbf) go to 250
C                         tgbf = dummmm
C                         c = dcos(latitu*grarad)
C                         l = distre/(c*c)
C                         xnew = latitu
C                         ynew = l
C                         if (l.gt.fltop .or. l.lt.flbot) go to 240
C c                           if (ii.eq.1) call frstpt(xnew,ynew)
C c                           call vector(xnew,ynew)
C                            ii = 0
C                            go to 250
C   300 continue
C c      call frstpt(xnew,ynew)
C       kfile = kfile+1
C       ktst = mod(kfile,numray)
C       nspace = 0
C       if (ktst.eq.1) nspace  = 1
C       write(7,500) kfile,ktst
C   500 format(' kfile  = ',i5,'    ktst  = ',i5)
C       kks = kfile-ktape-1
C c
C       if (kks) 400,330,400
C   330     ktape = kfile
C           if (nspace .ne. 0) go to 310
C              if (kplot.ge.3) go to 240
C                 go to 163
C   310     continue
C c          call frame
C 	  goto 60

C   400 return
C       end
