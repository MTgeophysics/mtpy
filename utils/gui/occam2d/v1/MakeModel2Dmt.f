C taken from
C http://mahi.ucsd.edu/Steve/Occam/code.html


	program MakeModel2Dmt
C OCCAM2DMT 2.1: Steven Constable IGPP/SIO La Jolla CA 92093-0225
c Program Revision 2.0, 14 Jan 1993
c
c constructs the standard 2D MT model, mesh and startup files from a minimum of 
c information.
c 
c
        implicit none

c      include 'imp.inc'
	integer maxdim
	parameter (maxdim = 500)
	character*80 modfil, datfil, itform , datetm, descr
	real sides(maxdim), bot(maxdim), dlz(maxdim), dly(maxdim)
	real width(maxdim), thick(maxdim), boffs
	integer nfe(maxdim), nfev(maxdim)
	real tolrq, pmu, rlast, tobt, aldp, t1, wmax, trigger, spcing
	real sitloc(maxdim), sitlok(maxdim),spcng1,spcng2, t, d1, d2
c sitlok() = dummy array to be used like sitloc()
	integer nrk, maxitr, iruf, nit, idebug, ifftol, nlay 
	integer n, nodey, nodez, ncol0, nbot, nside, ncol, nrc, np, nr
	integer j,i,k, ndata, nfre, mcol
c nrk = number of dummy site locations
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c set the number of elements in the side blocks (this controls the relaxation
c of the 2D model into the 1D boundary conditions).
	nside = 7
c set the number of elements in the bottom layer
	nbot = 4
c
c read in the data file
	write(*,*) ' enter data file name'
	read (*,'(a)') datfil
	nr = maxdim
      call jnputd(datfil, nr, sitloc, nrc, nfre, ndata)
      write(*,*) nrc,' site locations read'
c
	itform = 'OCCAMITER_1.0'
	modfil = 'INMODEL'
	descr = 'Startup file from MakeModel2Dmt_Linux32bit'
c	
	write(*,*) ' enter maxiumum number of iterations'
	read (*,*) maxitr
	write(*,*) ' enter required RMS misfit'
	read (*,*) tolrq
	iruf = 1
	idebug = 1
	nit = 0
	pmu = 5.0
	rlast = 1.0e+07
	tobt = 100.
	ifftol = 0
	
	write(*,*) ' enter number of layers'
	read (*,*) nlay
	write(*,*) ' enter number of layers per decade'
	read (*,*) aldp
	write(*,*) ' enter thickness of top layer (m)'
	read (*,*) t1
	write(*,*) ' enter maximum width of surface bricks'
	read (*,*) wmax
	write(*,*) ' enter trigger level for brick doubling (0.75 sugg.)'
	read (*,*) trigger
c
c We assume the site locations are ordered.  Sweep through them computing
c the intersite spacing and construct the mesh around these.
c At the moment there is one brick per station, two elements per station,
c with the node placed at the station.  
c Assumes at least 3 sites exist.
c
c I want to restrict the largest width the bricks can get.  It turned out 
c that one simple method for doing this was to add dummy stations.
c The first station is copied to the dummy array:
	j = 1
	sitlok(j) = sitloc(1)
c sweep through the other stations
	do 20 i = 2,nrc
c find the station spacing
	  spcing = sitloc(i) - sitloc(i-1)
c find out by how many extra site locations we want (plus the real one):
	  n = int(spcing/wmax) + 1
c add them
 	  do 15 k = 1, n
	    j = j + 1
15	    sitlok(j) = sitloc(i-1) + float(k)/float(n)*spcing
20	continue
	nrk = j
	write(*,*) j,' dummy station locations defined'
c 
c j tracks the number of finite elephants
c k tracks the number of regularization bricks
c boffs tracks the distance from the right edge of the left layer to the first
c station.
c First the left side elements 
	spcng1 = (sitlok(2) - sitlok(1))/2.
	sides(1) = 3.*spcng1
	do 30 i = 2,nside
	  sides(i) = sides(i-1)*3.
30      if (sides(i) .gt. 1000000.) sides(i) = 1000000.
	j = 0
	k = 1
	width(k) = 0.0
	do 31 i = nside,1,-1
	  j = j+1
	  width(k) = width(k) + sides(i)
31	  dly(j) = sides(i)
	nfe(k) = nside
c put one column of bricks on the left side for good measure
      dly(j+1) = spcng1
      dly(j+2) = spcng1
	j = j + 2
      nfe(k+1) = 2
	width(k+1) = 2.*spcng1
	boffs = width(k+1)
	k = k + 1
c and set the left hand station
      dly(j+1) = spcng1
      dly(j+2) = spcng1
	j = j + 2
      nfe(k+1) = 2
	width(k+1) = 2.*spcng1
	boffs = boffs + spcng1
	k = k + 1
c now sweep through the rest of the stations up to the right one
	do 100 i = 2, nrk-1
	  spcng2 = (sitlok(i+1) - sitlok(i))/2.
	  dly(j+1) = spcng1
	  dly(j+2) = spcng2
	  j = j + 2
	  nfe(k+1) = 2
	  width(k+1) = spcng1 + spcng2
 	  k = k + 1
	  spcng1 = spcng2
100	continue
c now fix the right hand station
      dly(j+1) = spcng2
      dly(j+2) = spcng2
	j = j + 2
      nfe(k+1) = 2
	width(k+1) = 2.*spcng2
	k = k + 1
c and add an extra brick
      dly(j+1) = spcng2
      dly(j+2) = spcng2
	j = j + 2
      nfe(k+1) = 2
	width(k+1) = 2.*spcng2
	k = k + 1
c and set the right side elements
	width(k) = 0.0
	sides(1) = 3.*spcng2
	do 130 i = 2,nside
	  sides(i) = sides(i-1)*3.
130      if (sides(i) .gt. 1000000.) sides(i) = 1000000.
	do 133 i=1,nside
	  j = j+1
	  width(k) = width(k) + sides(i)
133	  dly(j) = sides(i) 
 	k = k + 1
      nfe(k) = nside
	nodey = j+1
	ncol0 = k
c set the binding offset. 
	boffs = sitlok(1) - boffs
c
c now attend to the vertical coordinate and set up the array of 
c layer thicknesses
c Get a depth multiplier from the layers per decade
	t = 10.**(1./aldp)
c set the thickness of the first layer and initialize the depth counter
	thick(1) = t1
	d1 = t1
	do 25 i = 2, nlay-1
c increment the depth and get a thickness from the last two depths
	  d2 = d1*t
	  thick(i) = d2-d1
c don't allow any layer to get thicker than the top layer
	  if (thick(i) .lt. t1) thick(i) = t1
25	  d1 = d1 + thick(i)
c set the thickness of the elements in the bottom layer
	bot(1) = 3.*thick(nlay-1)
	do 35 i = 2,nbot
35	  bot(i) = bot(i-1)*3.
c now set up the finite element partioning
	k = 0
	dlz(1) = thick(1)/2.
	dlz(2) = thick(1)/2.
	nfev(1) = 2
	k = k + 2
	dlz(3) = thick(2)/2.
	dlz(4) = thick(2)/2.
	nfev(2) = 2
	k = k + 2
	do 36 i = 3,nlay-1
	  k = k+1
	  nfev(i) = 1
36	  dlz(k) = thick(i)
	do 37 i = 1, nbot
	  k = k+1
37	  dlz(k) = bot(i)
	nfev(nlay) = nbot
	nodez = k+1
c
c write the mesh filr
	open (12, file='MESH') 
	write(12,'(a)') ' MESH FILE FROM MAKE2DMODEL'
	write(12,'(6i6)') 0, nodey, nodez, 0, 0, 2
	write(12,'(8f10.1)') (dly(i), i=1,nodey-1)
	write(12,'(8f10.1)') (dlz(i), i=1,nodez-1)
	write(12,'(i6)') 0
	do 50 i = 1,nodez-1
	  do 50 k = 1,4
50	    write(12,'(170a1)') ('?', j = 1,nodey-1)
	close (12)
c
c write the model file
	open (12, file=modfil)
	write(12,'(a)') 'FORMAT:           '//'OCCAM2MTMOD_1.0'
	write(12,'(a)') 'MODEL NAME:       '//'MODELFILE FROM MAKE2DMODEL'
	write(12,'(a)') 'DESCRIPTION:      '//'SIMPLE INVERSION'
	write(12,'(a)') 'MESH FILE:        '//'MESH'
	write(12,'(a)') 'MESH TYPE:        '//'PW2D'
	write(12,'(a)') 'STATICS FILE:     none'
	write(12,'(a)') 'PREJUDICE FILE:   none'
	write(12,'(a18,f10.1)') 'BINDING OFFSET: ', boffs
	write(12,'(a18,i4)') 'NUM LAYERS:    ', nlay
c
c initialize the parameter counter
	np = 0
	ncol = ncol0
	do 260 i = 1, nlay
c sweep through all the bricks in this row to see we can double any up
c (note that this applies to the surface and bottom layers too; I think  
c this is reasonable)
	  k = 1
205	  k = k+1
210	  if (k .eq. ncol-2) goto 250
c (end of row (we don't want to mess with the right side block);
c go write this stuff out)
c
c test to see if the combined width of the next two bricks is greater than 
c the layer thickness     
	  if (thick(i) .gt. trigger*(width(k) + width(k+1)) ) then
c it is, double up the bricks
	    width(k) = width(k) + width(k+1)
	    nfe(k) = nfe(k) + nfe(k+1)
c condense the rest of the width() and nfe() arrays
	    do 220 j = k+2, ncol
	      width(j-1) = width(j)
220		nfe(j-1) = nfe(j)
c subtract one from the column count
	    ncol = ncol-1
c go test k again to see if this condensation bought us to the end of the row
	    goto 210
	  end if
c increment k and test again
	  goto 205
c this is where we dump the row out to file
250	  write(12,'(2i6)') nfev(i), ncol
	  write(12,'(50i3)') (nfe(j), j=1,ncol)
	  np = np + ncol
	  if (i.eq.1) mcol = ncol
c end of main loop
260	continue
c	
	write(12,'(a)') 'NO. EXCEPTIONS:   0'
c	
	close (12) 
c
c write the startup file
	open (12, file='startup')
c
      call datime(datetm)
c write stuff (it would be neater to call itrout(), part of occam, but
c the occam package is not linked in and itrout() is designed open a file
c called 'ITERxx', while here we want 'startup'.
	write(12,'(a)') 'FORMAT:           '//itform
	write(12,'(a)') 'DESCRIPTION:      '//descr
	write(12,'(a)') 'MODEL FILE:       '//modfil
	write(12,'(a)') 'DATA FILE:        '//datfil
	write(12,'(a)') 'DATE/TIME:        '//datetm
	write(12,'(a,i3)')    'MAX ITER:       ',maxitr
	write(12,'(a,f6.3)')  'REQ TOL:         ',tolrq
	write(12,'(a,i3)')    'IRUF:           ',iruf
	write(12,'(a,i3)')    'DEBUG LEVEL:    ',idebug
	write(12,'(a,i3)')    'ITERATION:      ',nit
	write(12,'(a,g15.7)') 'PMU:           ',pmu
	write(12,'(a,g15.7)') 'RLAST:          ',rlast
	write(12,'(a,g15.7)') 'TLAST:         ',tobt
	write(12,'(a,i3)')    'IFFTOL:         ',ifftol
	write(12,'(a,i6)')    'NO. PARMS:       ',np
	write(12,'(4g15.7)') (2.0, i=1,np)
c
	close(12)
c	
c now indicate what the dimensions should be:
      write(*,*) ' Dimensions:'
      write(*,*) '     npp, ndd (occamdim.inc) = ', np,',', ndata
      write(*,*) '     ipy, ipz (mt2Ddim.inc) = ',nodey,',', nodez+10
      write(*,*) '     ipi (mt2Ddim.inc) = ', np
      write(*,*) '     nl (mt2Ddim.inc) = ', nlay
      write(*,*) '     nf, nr (mt2Ddim.inc) = ', nfre,',', nrc
      write(*,*) '     mcol (mt2Ddim.inc) = ', mcol
c
	stop
	end
c-----------------------------------------------------------------------
      subroutine jnputd(datfil, nr, sitloc, nrc, nfre, ndata)
c-----------------------------------------------------------------------
c
c OCCAM2DMT 2.0: Steven Constable IGPP/SIO La Jolla CA 92093-0225
c Subroutine Revision 1.00, 13 Jan 1993
c
c jnputd opens the data file datfil and reads in just the site locations
c
c on input:
c   datfil = data file name 
c   nr = dimension of sitloc
c
c on output:
c   sitloc() = site locations
c   nrc = number of sites
c   ndata = number of data
c
c calls:  
c   filerr, idcode
c
c includes:
c      include 'imp.inc'
c
c parameters:
      integer iof 
      parameter (iof = 15)
c   iof = io unit number for file operations
c arguments:
      character*(*) datfil
	integer nrc, nr, ndata, nfre
	real sitloc(nr)
c local variables
      integer lerr, i
      character*80  string
	real trash
      integer idcode
c   idcode = function to perform the equivalent of a free-format internal
c
c-----------------------------------------------------------------------
c open data file 
      open (iof, file=datfil, status='old', iostat=lerr)
      if (lerr .ne. 0) then
        write(*,*) ' Error opening data file'
        stop
      end if
c read stuff
      read(iof,'(18x,a)', end=198, err=199) string
      if (string(1:16) .ne. 'OCCAM2MTDATA_1.0') call filerr(
     *  ' Data file not supported',iof)
      read(iof,'(18x,a)', end=198, err=199) string
      read(iof,'(18x,a)', end=198, err=199) string
      nrc = idcode(string, iof)
      if (nrc .gt. nr)    call filerr (' Too many sites', iof)
c skip past the site names
      do 5 i = 1, nrc
5        read(iof,'(a)', end=198, err=199) string
      read(iof,'(a)', end=198, err=199) string
c get the site locations
      read(iof,*, end=198, err=199) (sitloc(i), i=1,nrc)
c read past the frequencies
      read(iof,'(18x,a)', end=198, err=199) string
      nfre = idcode(string, iof)
      read(iof,*, end=198, err=199) (trash, i=1,nfre)
      read(iof,'(18x,a)', end=198, err=199) string
      ndata = idcode(string, iof)
c      
      return
c
198      call filerr (' Data file ended prematurely', iof)
199      call filerr (' Error reading data file', iof)
c
      end
c
C-----------------------------------------------------------------------
      subroutine filerr(mssg, io1)
C-----------------------------------------------------------------------
C
C OCCAM 2.0: Steven Constable IGPP/SIO La Jolla CA 92093-0225
c Subroutine Revision 2.00, 13 May 1992
C
c filerr prints an error message, closes file io1, and stops
c
c on input:
c    mssg = character string containing error message
c    i01 = unit number of file to close (0 = no file open)
c on output:
c    outputs nothing
c calls:
c    no other routines
c
c
c includes:
c      include 'imp.inc'
c input arguments:
      character*(*) mssg
      integer io1
c
C-----------------------------------------------------------------------
      write(*,*) mssg
      if (io1 .gt. 0) close (io1)
      stop
      end
C-----------------------------------------------------------------------
      integer function idcode(string, iof2)
C-----------------------------------------------------------------------
C
C OCCAM 2.0: Steven Constable IGPP/SIO La Jolla CA 92093-0225
c Subroutine Revision 2.00, 13 Jan 1993
C
c idcode provides a patch to get around the lack of an internal
c free-format (list-directed) read on some compilers, notably the
c cray.  This version for integers.  Butchered from decode by Bob Parker.
c 
c on input:
c    string = character string containing number to be read
c    iof2 = unit number of open file, to be closed on error 
c on output:
c    idcode contains an integer number from string
c calls:
c    filerr
c
c includes:
c      include 'imp.inc'
c input arguments:
        implicit none
        character*(*) string
	integer iof2
c local variables
      character*40 local
	integer k, k1, kn, l
	real realno
c
c  Terminate line at % sign
      kn=index(string//'%', '%') - 1
      k1=1
      do 1100 k=k1, kn
        if (string(k:k) .ne. ' ') goto 1200
 1100 continue
      call filerr (' No number found by idcode', iof2)
 1200 do 1300 l=k, kn
        if (string(l:l).eq. ',' .or. string(l:l) .eq. ' ') goto 1500
 1300 continue
 1500 local=string(k:l-1)
      read (local, '(f40.0)', err=1900) realno
	idcode = int(realno + 0.001)
	return
 1900 call filerr (' Error reading string in idcode', iof2)
      end
C-----------------------------------------------------------------------
      subroutine datime(datetm)
c-----------------------------------------------------------------------
c OCCAM 2.0: Steven Constable IGPP/SIO La Jolla CA 92093-0225
c Subroutine Revision 1.00, 18 May 1992
c
c Dummy date and time utility for OCCAM 2.0.  
c
c on input:
c   nothing
c on output:
c   datetm = 80 character string containing message saying that date and time
c           are not available
c      
      character*80 datetm
c      
c-----------------------------------------------------------------------
      datetm = 'No date and time available on this system'      
      return
      end
c-----------------------------------------------------------------------
      
