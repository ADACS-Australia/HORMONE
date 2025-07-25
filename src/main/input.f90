module input_mod
 implicit none

 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE READ_MESA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To read MESA file

subroutine read_mesa(mesafile,r,m,rho,pres,mu,comp,comp_list)

 use constants,only:msun,rsun
 use composition_mod,only:get_imu

 implicit none

 character(len=100),intent(in):: mesafile
 character(len=30),allocatable::header(:),dum(:)
 real(8),allocatable,dimension(:),intent(out):: r, m, rho, pres
 real(8),allocatable,dimension(:,:),intent(out),optional:: mu, comp
 character(len=10),allocatable,intent(out),optional::comp_list(:)
 real(8),allocatable,dimension(:,:):: dat,Xfrac
 character(len=10000):: dumc
 integer:: nn,ui, i,j, lines, columns, nrel, nel,istat
 character(len=10),allocatable:: element_list(:)


!-----------------------------------------------------------------------------

! list of relevant elements ! ------------------------------------------------
 nrel = 20 ! number of elements in the list
 allocate(element_list(nrel))
 element_list( 1) = 'h1'
 element_list( 2) = 'he3'
 element_list( 3) = 'he4'
 element_list( 4) = 'c12'
 element_list( 5) = 'n14'
 element_list( 6) = 'o16'
 element_list( 7) = 'ne20'
 element_list( 8) = 'mg24'
 element_list( 9) = 'si28'
 element_list(10) = 's32'
 element_list(11) = 'ar36'
 element_list(12) = 'ca40'
 element_list(13) = 'ti44'
 element_list(14) = 'cr48'
 element_list(15) = 'cr60'
 element_list(16) = 'fe52'
 element_list(17) = 'fe54'
 element_list(18) = 'fe56'
 element_list(19) = 'co56'
 element_list(20) = 'ni56'

! reading data from datafile ! -----------------------------------------------
 open(newunit=ui,file=mesafile,status='old',iostat=istat)
 if(istat/=0)then
  print*,'Error: Input MESA file cannot be found.'
  print*,'Make sure to specify the correct file name.'
  print'(3a)',"mesafile='",trim(mesafile),"'"
  stop
 end if
 read(ui,'()');read(ui,'()')
 read(ui,*) lines, lines
 read(ui,'()');read(ui,'()')
 read(ui,'(a)') dumc

! counting columns
 allocate(dum(500)) ; dum = 'aaa'
 read(dumc,*,end=101) dum
101 do i = 1, 500
  if(dum(i)=='aaa')then
   columns = i - 1
   exit
  end if
 end do

 allocate(header(1:columns),dat(1:lines,1:columns))
 header(1:columns) = dum(1:columns)
 deallocate(dum)

! find relevant chemical elements in the headers
 if(present(comp))then
  nel = 0
  do j = 1, nrel
   do i = 1, columns
    if(trim(header(i))==trim(element_list(j)))then
     nel = nel + 1 ! first count how many elements
     exit
    end if
   end do
  end do
  allocate(comp_list(nel),comp(1:nel,0:lines))
  nn = 1
  element_loop: do j = 1, nrel
   do i = 1, columns
    if(trim(header(i))==trim(element_list(j)))then
     comp_list(nn) = trim(header(i))
     nn = nn + 1
     if(nn>nel)exit element_loop
     exit
    end if
   end do
  end do element_loop
 end if

 do i = 1, lines
  read(ui,*) dat(lines-i+1,1:columns)
 end do

 allocate(m(0:lines))
 allocate(r,rho,pres,mold=m)
 if(present(mu))then
  allocate(mu(0:lines,0:1))
  allocate(Xfrac,mold=mu)
  mu(1,0) = -1d0
 end if

 do i = 1, columns
  select case(trim(header(i)))
  case('mass')
   m(1:lines) = dat(1:lines,i) * msun
  case('logM','log_mass')
   m(1:lines) = 10d0**dat(1:lines,i) * msun
  case('density','rho')
   rho(1:lines) = dat(1:lines,i)
  case('logRho')
   rho(1:lines) = 10d0**dat(1:lines,i)
  case('radius_cm')
   r(1:lines) = dat(1:lines,i)
  case('radius_km')
   r(1:lines) = dat(1:lines,i) * 1d5
  case('radius')
   r(1:lines) = dat(1:lines,i) * rsun
  case('logR')
   r(1:lines) = 10d0**dat(1:lines,i) * rsun
  case('logR_cm')
   r(1:lines) = 10d0**dat(1:lines,i)
  case('pressure')
   pres(1:lines) = dat(1:lines,i)
  case('logP')
   pres(1:lines) = 10d0**dat(1:lines,i)
  case('mu')
   if(present(mu)) mu(1:lines,1) = dat(1:lines,i)
  case('h1','x_mass_fraction_H')
   if(present(mu)) Xfrac(1:lines,0) = dat(1:lines,i)
  case('he4','y_mass_fraction_He')
   if(present(mu)) Xfrac(1:lines,1) = dat(1:lines,i)
  end select
  if(present(comp))then
   do j = 1, nel
    if(trim(header(i))==trim(comp_list(j))) comp(j,1:lines) = dat(1:lines,i)
   end do
  end if
 end do

 close(ui)

! Set central values
 r(0) = 0d0
 m(0) = 0d0
 rho(0) = rho(1)
 pres(0) = pres(1)
 if(present(mu))then
  if(mu(1,0)<0d0)then ! compute mu from X and Y if not outputted
   do i = 1, lines
    mu(i,1) = 1d0/get_imu(Xfrac(i,0:1))
   end do
  end if
  mu(:,0) = m(:)
  mu(0,1) = mu(1,1)
 end if
 if(present(comp))comp(:,0) = comp(:,1)

return
end subroutine read_mesa


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE ERROR_EXTRAS
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output error message and stop simulation when failing to find extras file.

subroutine error_extras(simutype,extrasfile)

 character(len=*),intent(in):: simutype,extrasfile

!-----------------------------------------------------------------------------

 print*,'Error: Model parameter file cannot be found.'
 print'(5a)','Copy over "../para/extras_',simutype,'" to "',trim(extrasfile),&
             '" and specify model parameters'
 stop

end subroutine error_extras

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                             SUBROUTINE ERROR_NML
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output error message and stop simulation when relevant namelist cannot be found in extras file

subroutine error_nml(simutype,extrasfile)

 character(len=*),intent(in)::simutype,extrasfile

!-----------------------------------------------------------------------------

 print*,'Error: extras file does not contain relevant namelist'
 print'(5a)','Copy over contents of "../para/extras_',simutype,'" to "',&
             trim(extrasfile),'" and specify model parameters'
 stop

return
end subroutine error_nml

end module input_mod
