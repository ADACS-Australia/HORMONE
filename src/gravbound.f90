!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE GRAVBOUND
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary conditions for miccg with companion star

 subroutine gravbound

  use settings,only:bc3is
  use grid
  use constants
  use physval
  use gravmod

  implicit none

  integer error1
  real*8 dphiii, dphiio, dphi1o, dphi3i, dphi3o

!------------------------------------------------------------------------------

! cylindrical (axial symmetry) ################################################
  if(crdnt==1.and.je==1.and.bc3is/=1)then

   call multipole

   error1=0
   do i = gie+1,gie+2
    do k=gks,gke
     phi1o(i,k) = 0.0d0
     do ll=0,llmax
      dphi1o = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**dble(ll)/rdis(i,k)
      if(ll/=0.and.abs(dphi1o/phi1o(i,k)) < grverr) exit
      phi1o(i,k) = phi1o(i,k) + dphi1o
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi1o(i,k) = G * phi1o(i,k)! + pripot1o(i,k)
    end do !k-loop

    if(error1==1)then
     write(6,*)"Error from gravbound 1: Number of terms in multipole expansion &
         &is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do k = gks, gke
    do i = gie+1, gie+2
     grvphi(i,js,k) = phi1o(i,k)
    end do
   end do

   error1=0
   do k = gks-2, gks-1
    do i=gis,gie
     phi3i(i,k) = 0.0d0
     do ll=0,llmax
      dphi3i = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**dble(ll)/rdis(i,k)
      if(ll/=0.and.abs(dphi3i/phi3i(i,k)) < grverr) exit
      phi3i(i,k) = phi3i(i,k) + dphi3i
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi3i(i,k) = G*phi3i(i,k)! + pripot3i(i)
    end do !i-loop

    if(error1==1)then
     write(6,*)"Error from gravbound3i: Number of terms in multipole expansion &
        &  is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   error1=0
   do k = gke+1,gke+2
    do i=gis,gie
     phi3o(i,k) = 0.0d0
     do ll=0,llmax
      dphi3o = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**dble(ll)/rdis(i,k)
      if(ll/=0.and.abs(dphi3o/phi3o(i,k)) < grverr) exit
      phi3o(i,k) = phi3o(i,k) + dphi3o
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi3o(i,k) = G*phi3o(i,k)! + pripot3o(i)
    end do !i-loop

    if(error1==1)then
     write(6,*)"Error from gravbound3o: Number of terms in multipole expansion &
         & is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do i = gis, gie
    do k = gke+1, gke+2
     grvphi(i,js,k) = phi3o(i,k)     
    end do
    do k = gks-2, gks-1
     grvphi(i,js,k) = phi3i(i,k)
    end do
   end do

! cylindrical (equatorial+axial symmetry) #####################################
  elseif(crdnt==1.and.je==1.and.bc3is==1)then

   call multipole

   error1=0
   do i = gie+1, gie+2
    do k=gks,gke
     phi1o(i,k) = 0.0d0
     do ll=0,llmax,2
      dphi1o = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**dble(ll)/rdis(i,k)
      if(ll/=0.and.abs(dphi1o/phi1o(i,k)) < grverr) exit
      phi1o(i,k) = phi1o(i,k) + dphi1o
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi1o(i,k) = 2d0*G* phi1o(i,k)! + pripot1o(i,k)
    end do !k-loop

    if(error1==1)then
     write(6,*)"Error from gravbound 1: Number of terms in multipole expansion &
         & is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do i = ie+1, ie+2
    do k = ks, ke
     grvphi(i,js,k) = phi1o(i,k)
    end do
   end do

   error1=0
   do k = gke+1, gke+2
    do i=gis,gie
     phi3o(i,k) = 0.0d0
     do ll=0,llmax,2
      dphi3o = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**dble(ll)/rdis(i,k)
      if(ll/=0.and.abs(dphi3o/phi3o(i,k)) < grverr) exit
      phi3o(i,k) = phi3o(i,k) + dphi3o
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi3o(i,k) = 2d0*G*phi3o(i,k)! + pripot3o(i,k)
    end do !i-loop

    if(error1==1)then
     write(6,*)"Error from gravbound3o: Number of terms in multipole expansion &
         & is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do i = gis, gie
    do k = gke+1, gke+2
     grvphi(i,js,k) = phi3o(i,k)
    end do
   end do
   grvphi(is:ie,js,ks-1) = grvphi(is:ie,js,ks)
   grvphi(is:ie,js,ks-2) = grvphi(is:ie,js,ks+1)
   phi3i = 0d0
   grvphi(is-1,js,ks:ke) = grvphi(is,js,ks:ke)
   grvphi(is-2,js,ks:ke) = grvphi(is+1,js,ks:ke)

! spherical coordinates (axial symmetry) #######################################
  elseif(crdnt==2.and.ke==1)then

   call multipole

   error1=0
   do i = gie+1, gie+2
    do j = gjs,gje
     phiio(i,j) = 0.0d0
     do ll=0,llmax
      dphiio = -G*Pl(ll,j)*ml(ll)/x1(i)
      if(ll/=0.and.abs(dphiio/phiio(i,j)) < grverr) exit
      phiio(i,j) = phiio(i,j) + dphiio
     end do !ll-loop
     phiio(i,j) = phiio(i,j) - G*mc(is-1)/x1(i)
     if(ll>=llmax) error1=1
    end do

    if(error1==1)then
     write(6,*)"Error from gravbound i: Number of terms in multipole expansion &
        & is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do i = gie+1, gie+2
    do j = gjs, gje
     grvphi(i,j,ks) = phiio(i,j)
    end do
   end do

   if(xi1s>0d0)then
    error1=0
    do i = gis-2, gis-1
     do j = gjs, gje
      phiii(i,j) = 0.0d0
      do ll=0,llmax
       call multipoleinner
       dphiii = -G*Pl(ll,j)*ml(ll)/x1(is-1)
       if(ll/=0.and.abs(dphiii/phiii(i,j)) < grverr) exit
       phiii(i,j) = phiii(i,j) + dphiii
      end do !ll-loop
      phiii(i,j) = phiii(i,j) - G*mc(is-1)/x1(is-1)
      if(ll>=llmax) error1=1
     end do

     if(error1==1)then
      write(6,*)"Error from gravbound i: Number of terms in inner multipole &
           & expansion is not enough. Computation stopped. tn=",tn,"i=",i
      stop
     end if
    end do
   end if


   if(mc(is-1)==0.d0.and.xi1s/=0.d0)error1=1
   if(error1==1)then
    write(6,*)"Error from gravbound i: inner boundary is wrong",xi1s,mc(is-1)
    stop
   end if

  end if

return
end subroutine gravbound


!------------------------------------------------------------------------------
 subroutine multipole

   use grid
   use physval
   use gravmod

   implicit none

   integer ii, jj, kk

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(crdnt==1.and.je==1)then
   jj = js
   ml = 0d0
!$omp parallel do private (ll,kk,ii)
   do ll = 0, llmax
    do kk = ks,ke
     do ii = is,ie
      ml(ll) = ml(ll) + d(ii,jj,kk)*(rdis(ii,kk)/rdis(ie,ke))**dble(ll) &
             * Plc(ll,ii,kk)*dvol(ii,jj,kk)
     end do
    end do
   end do
!$omp end parallel do
   
! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elseif(crdnt==2.and.ke==1)then
   kk = ks
   do ll = 0, llmax
    ml(ll) = 0.0d0
    do jj = js,je
     do ii = is,ie
      ml(ll) = ml(ll) + d(ii,jj,kk)*(x1(ii)/x1(ie+1))**dble(ll) &
             * Pl(ll,jj) * dvol(ii,jj,kk)
     end do
    end do
   end do

  end if

  return
  end subroutine multipole



!------------------------------------------------------------------------------
 subroutine multipoleinner

   use grid
   use physval
   use gravmod

   implicit none

   integer ii, jj, kk

!------------------------------------------------------------------------------

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(crdnt==2.and.ke==1)then
   kk = ks
   ml(ll) = 0.0d0
   do jj = js,je
    do ii = is,ie
     ml(ll) = ml(ll) + d(ii,jj,kk)*(x1(is-1)/x1(ii))**(ll+1) &
              * Pl(ll,jj) * dvol(ii,jj,kk)
    end do
   end do

  end if

  return
  end subroutine multipoleinner