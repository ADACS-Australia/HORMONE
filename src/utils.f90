module utils
 implicit none

contains

! convert polar to cartesian coordinates
 pure function polcar(xp) result(x)
  implicit none
  real*8,intent(in):: xp(1:3)
  real*8:: x(1:3)
  x(1) = xp(1)*sin(xp(2))*cos(xp(3))
  x(2) = xp(1)*sin(xp(2))*sin(xp(3))
  x(3) = xp(1)*cos(xp(2))
 end function polcar

! convert cartesian to polar coordinates
 pure function carpol(x) result(xp)
  implicit none
  real*8,intent(in):: x(1:3)
  real*8:: xp(1:3)
  xp(1) = norm2(x)
  xp(2) = acos(x(3)/xp(1))
  xp(3) = atan2(x(2),x(1))
 end function carpol


 pure function softened_potential(r,h) result(phi)
! From Price & Monaghan 2007
 implicit none
 real*8,intent(in):: r,h
 real*8:: phi,q
 q=r/h
 if(q>=2d0)then
  phi = -1d0/r
 elseif(q>=1d0)then
  phi = (4d0/3d0*q**2-q**3+0.3d0*q**4-q**5/30d0-1.6d0+1/(15d0*q))/h
 elseif(q>=0d0)then
  phi = (2d0/3d0*q**2-0.3d0*q**4+0.1d0*q**5-1.4d0)/h
 end if
 
end function softened_potential


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE GET_VCAR
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate the Cartesian vector components from polar coordinates

subroutine get_vcar(xcar,x3,v1,v2,v3,vcar)

 use constants,only:pi
 
 implicit none

 real*8,intent(in):: xcar(1:3),x3,v1,v2,v3
 real*8,intent(out)::vcar(1:3)
 real*8,dimension(1:3)::uvec1,uvec2,uvec3

!-----------------------------------------------------------------------------

 uvec1 = xcar/norm2(xcar)
 uvec2 = rotz(roty(rotz(uvec1,-x3),0.5*pi),x3)
 uvec3 = cross(uvec1,uvec2)
 
 vcar = v1*uvec1 + v2*uvec2 + v3*uvec3 

return
end subroutine get_vcar


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE GET_VPOL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate the Polar vector components from Cartesian coordinates

subroutine get_vpol(xcar,x3,vcar,v1,v2,v3)

 use constants,only:pi
 
 implicit none

 real*8,intent(in):: xcar(1:3),x3,vcar(1:3)
 real*8,intent(out)::v1,v2,v3
 real*8,dimension(1:3)::uvec1,uvec2,uvec3

!-----------------------------------------------------------------------------

 uvec1 = xcar/norm2(xcar)
 uvec2 = rotz(roty(rotz(uvec1,-x3),0.5*pi),x3)
 uvec3 = cross(uvec1,uvec2)
 
 v1 = dot_product(vcar,uvec1)
 v2 = dot_product(vcar,uvec2)
 v3 = dot_product(vcar,uvec3)

return
end subroutine get_vpol

! rotate about the x axis
 pure function rotx(x,theta) result(xp)
  implicit none
  real*8,intent(in):: x(1:3), theta
  real*8:: xp(1:3)
  xp(1) = x(1)
  xp(2) = cos(theta)*x(2) - sin(theta)*x(3)
  xp(3) = sin(theta)*x(2) + cos(theta)*x(3)
 end function rotx

! rotate about the y axis
 pure function roty(x,theta) result(xp)
  implicit none
  real*8,intent(in):: x(1:3), theta
  real*8:: xp(1:3)
  xp(1) = cos(theta)*x(1) + sin(theta)*x(3)
  xp(2) = x(2)
  xp(3) =-sin(theta)*x(1) + cos(theta)*x(3)
 end function roty

! rotate about the z axis
 pure function rotz(x,theta) result(xp)
  implicit none
  real*8,intent(in):: x(1:3), theta
  real*8:: xp(1:3)
  xp(1) = cos(theta)*x(1) - sin(theta)*x(2)
  xp(2) = sin(theta)*x(1) + cos(theta)*x(2)
  xp(3) = x(3)
 end function rotz

! cross product
 pure function cross(x,y) result(z)
  implicit none
  real*8,intent(in):: x(1:3),y(1:3)
  real*8:: z(1:3)
  z(1) = x(2)*y(3)-x(3)*y(2)
  z(2) = x(3)*y(1)-x(1)*y(3)
  z(3) = x(1)*y(2)-x(2)*y(1)
 end function cross

 ! linear interpolation
 pure function intpol(x,y,z) result(val)
  implicit none
  real*8,intent(in):: x(1:2), y(1:2)
  real*8,intent(in):: z
  real*8:: val
  val = ((x(2)-z)*y(1)+(z-x(1))*y(2))/(x(2)-x(1))
 end function intpol

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                      SUBROUTINE GEOMETRICAL_SERIES
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate dxi's in a geometrical series.

subroutine geometrical_series(dxi,xmin,is,ie,xis,xie)

 implicit none

 integer,intent(in):: is,ie
 real*8,intent(in):: xis,xie,xmin
 real*8,intent(inout),allocatable:: dxi(:)
 integer i
 real*8 xrng, irng, xr, xrnew, xrmax, err, maxerr, fx, dfx

!-----------------------------------------------------------------------------

 xrmax = 1.15d0
 maxerr = 1d-10

 xr = 1.01d0
 xrng = xie - xis ; irng = dble(ie - is + 1)

 if(xrng/irng<xmin)then
  print *,"Error from geometrical_series ;"
  print *,"xmin should be smaller or uniform mesh should be chosen",xmin
  stop
 end if

 do i = 1, 10000000
  fx = (xr-1d0)*xrng - xmin * (xr**irng-1d0)
  dfx = xrng - irng * xmin * xr**(irng-1d0)

  xrnew = xr - fx/dfx

  if(abs((xrnew-xr)/xr)<maxerr)then
   xr = xrnew ; exit
  end if
  if(xrnew<1d0)xrnew = 2d0

  xr = xrnew
 end do

 if(xr>xrmax)then
  print *,"xmin too small", xmin, xr
  stop
 end if

 dxi(is) = xmin
 do i = is+1, ie
  dxi(i) = dxi(i-1) * xr
 end do
 dxi(is-1) = dxi(is) ; dxi(is-2) = dxi(is+1)
 dxi(ie+1) = dxi(ie)*xr ; dxi(ie+2) = dxi(ie)*xr*xr

 if(xr-1d0<maxerr) dxi = (xie-xis) / irng

return
end subroutine geometrical_series


end module utils