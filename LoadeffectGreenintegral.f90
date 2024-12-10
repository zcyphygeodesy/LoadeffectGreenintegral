!  LoadeffectGreenintegral.f90 
!
!  FUNCTIONS:
!  LoadeffectGreenintegral - Entry point of console application.
!
!****************************************************************************

      program LoadeffectGreenintegral
      implicit none
	character*800::loadgrfl
	character*80000::line
	integer i,j,k,n,m,nn,kk,sn,kln,astat(5)
      real*8 tdn(24),rec(8000),GF(8000,9),BLH(3),direct(10),indrct(14)
	real*8 GRS(6),dcin(40),dout(40),hd(6),ae,dr,L0,RAD
	integer nlon,nlat,maxn,kndld,nk,hgt,kpar(12)!
	integer::status=0
	real*8,allocatable::ewh(:,:)
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0
      dr=300.d3!Green's integral radius (m)
      GRS(6)=1.d0!=1.d0 land water or sea level variation load, =-1.d0 surface atmosphere load
      !read the load EWH grid (cm)
      open(unit=8,file="swscSEP2018041112.dat",status="old",iostat=status)
      if(status/=0) goto 902
      read(8,'(a)') line
      call PickReclong(line,kln,rec,sn)
      if(sn<6)then
         close(8);goto 902
      endif
      hd(1:6)=rec(1:6)
	nlat=nint((hd(4)-hd(3))/hd(6))
	nlon=nint((hd(2)-hd(1))/hd(5))
	hd(5)=(hd(2)-hd(1))/dble(nlon)
	hd(6)=(hd(4)-hd(3))/dble(nlat)
 	allocate(ewh(nlat,nlon), stat=astat(1))
 	do i=1,nlat
	   read(8,*,end=905)(ewh(i,j),j=1,nlon)
          do j=1,nlon
            if(ewh(i,j)>9.d3)ewh(i,j)=0.d0
          enddo
      enddo
905   close(8)
      RAD=datan(1.d0)/45.d0;ewh=ewh*1.d-2!cm->m
      !read load green functions for indirect load effect 负荷格林函数(间接影响)
	write(loadgrfl,*) "LoadGreen.txt"
      call LGrnFunc(loadgrfl,GF)
      open(unit=8,file="calcpnt.txt",status="old",iostat=status)
      if(status/=0)goto 904
      open(unit=10,file="reslt.txt",status="replace")!Output file 输出文件
      read(8,'(a)') line
      write(10,'(a)')trim(line)
      kk=0
      do while(.not.eof(8))  
         read(8,'(a)') line
         call PickReclong(line,kln,rec,sn)
         if(sn<4)goto 906
         BLH(1)=rec(3);BLH(2)=rec(2);BLH(3)=rec(4)
         kk=kk+1;tdn(1:14)=0.d0;direct=0.d0;indrct=0.d0
         call rntGreenintegral(BLH,ewh,hd,nlat,nlon,GF,direct,indrct,GRS,dr)
         tdn(1:7)=direct(1:7)+indrct(1:7);tdn(8:10)=indrct(8:10)
         tdn(11)=tdn(10)-tdn(1);tdn(12:14)=direct(8:10)+indrct(12:14)
         write(10,'(a,40F14.4)')trim(line),(tdn(i),i=1,14)
         if(kk/500*500==kk)write(*, '(a,i9)'), '    Calculated point number: ',kk
906      continue
	enddo
 	close(8)
      close(10)
904   deallocate(ewh)
902	continue
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      pause
      end
