module ejtfilemod

  implicit none

  character*30 ejtfile,ejtfilebin
  real*8 t_ejt_out, dt_ejt_out, remmass, inimass, belmass, totmass, totangmom, totangmoma, totangmomb, totenergy, totenergya, totenergyb, ToW, ToWb, ToWt, totintene, totinteneb, totintenet, totgrvene,totgrveneb,totgrvenet,totkinene,totkineneb,totkinenet, centre_of_mass,centre_of_massb,vel,velb

end module ejtfilemod
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE OUTPUT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output data

subroutine output

  use settings,only:outstyle,gravswitch,compswitch,spn,eq_sym,mag_on
  use grid
  use physval
  use constants
  use gravmod
  use ejtfilemod
  use particle_mod

  implicit none

  logical extflag
  character*35 pltfile, binfile, ptcfile, bptfile
  real*8 phih, vsq, pr, pt
  real*8,allocatable:: lapphi(:,:,:)

!----------------------------------------------------------------------------


!gridfile---------------------------------------------------------------
  if(tn==0)then
   open(unit=40,file='data/gridfile.bin',status='replace',form='unformatted')

   write(40) x1  (gis-2:gie+2), xi1(gis-2:gie+2), &
             dxi1(gis-2:gie+2), dx1(gis-2:gie+2), &
             x2  (gjs-2:gje+2), xi2(gjs-2:gje+2), &
             dxi2(gjs-2:gje+2), dx2(gjs-2:gje+2), &
             x3  (gks-2:gke+2), xi3(gks-2:gke+2), &
             dxi3(gks-2:gke+2), dx3(gks-2:gke+2)
   close(40)

   open(unit=50,file='data/gridfile.dat',status='replace')
   write(50,'()')
   write(50,'(2(a5),3(a15))')'#  i','k','x1(3)','x3(4)','dvol(5)'

   if(crdnt==2)then
    k=ks;j=js-1
    do i = is, ie
     write(50,'(2(i5),3(1x,1PE14.6e2))')i,j,x1(i),0d0,dvol(i,j,k)
    end do
    write(50,*)
   end if

   do k = ks,ke,2
    do j = js,je
     do i = is,ie,2
      write(50,'(2(i5),3(1x,1PE14.6e2))')i,k,x1(i),x3(k),dvol(i,j,k)
     end do
     if(ie/=1)write(50,*)
    end do
    !   if(ie==1)write(30,*)
   end do

   if(crdnt==2)then
    k=ks;j=je+1
    do i = is, ie
     if(eq_sym)then
      write(50,'(2(i5),3(1x,1PE14.6e2))')i,j,x1(i),0.5d0*pi,dvol(i,j,k)
     else
      write(50,'(2(i5),3(1x,1PE14.6e2))')i,j,x1(i),pi,dvol(i,j,k)
     end if
    end do
    write(50,*)
   end if

  end if


  close(50)


!pltfile---------------------------------------------------------------------
  if(outstyle==1)then
   write(pltfile,'(a8,i11.11,a5)')'data/plt',int(time),'s.dat'
   if(time>2147483647d0)then
    write(pltfile,'(a8,i9.9,a7)')'data/plt',int(time*0.01d0),'00s.dat'
   end if
  elseif(outstyle==2)then
   write(pltfile,'(a8,i8.8,a4)')'data/plt',tn,'.dat'
  end if

  open(unit=30,file = pltfile, status='replace')

  write(30,'(a,i7,2(a,1PE12.4e2))') '#tn =',tn,'  time= ',time
  write(30,'(19(a15))')'d(6)','e(7)','p(8)','v1(9)','v2(10)','phi(11)','T(12)','mu(13)','X_h1(14)','X_He4(15)','X_He3(16)','X_C12(17)','X_N14(18)','X_O16(19)','X_Ne20(20)','X_Mg(21)','shock(22)'

  if(gravswitch==1)then
   j = js ; k = ks
   do i = is, ie-1
    grvphi(i,js:je,k) = 0d0
    do n = i+1, ie
     grvphi(i,js:je,k) = grvphi(i,js:je,k) - sum(d(n,js:je,k)*dvol(n,js:je,k))/sum(dvol(n,js:je,k))*x1(n)*dxi1(n)
    end do
    grvphi(i,js:je,k) = G*(-mc(i)/x1(i)+4d0*pi*grvphi(i,js:je,k))
   end do
   grvphi(ie,js:je,k) = -G*mc(ie)/x1(ie)
  end if

  call shockfind
  
  if(crdnt==2)then
   j =js;k=ks
   do i = is, ie
    write(30,'(16(1x,1PE14.6e2),i2)')&
         d(i,j,k),e(i,j,k),p(i,j,k),v1(i,j,k),v2(i,j,k), &
         grvphi(i,j,k)+extgrv(i,j,k),T(i,j,k),1d0/imu(i,j,k),spc(1:spn,i,j,k),shock(i,j,k)
   end do
   write(30,'()')
  end if

  do k = ks,ke,2
   do j = js,je
    do i = is,ie,2
     write(30,'(16(1x,1PE14.6e2),i9)')&
          d(i,j,k),e(i,j,k),p(i,j,k),v1(i,j,k),v3(i,j,k), &
          grvphi(i,j,k)+extgrv(i,j,k),T(i,j,k),1d0/imu(i,j,k),spc(1:spn,i,j,k),shock(i,j,k)!,lapphi(i,j,k)
    end do
    if(ie/=1)write(30,*)
   end do
!   if(ie==1)write(30,*)
  end do

  if(crdnt==2)then
   j =je;k=ks
   do i = is, ie
    write(30,'(16(1x,1PE14.6e2),i9)')&
         d(i,j,k),e(i,j,k),p(i,j,k),v1(i,j,k),v2(i,j,k), &
         grvphi(i,j,k)+extgrv(i,j,k),T(i,j,k),1d0/imu(i,j,k),spc(1:spn,i,j,k),shock(i,j,k)
   end do
   write(30,'()')
  end if


  close(30)


!binfile----------------------------------------------------------------
  if(outstyle==1)then
   write(binfile,'(a8,i11.11,a5)')'data/bin',int(time),'s.dat'
   if(time>2147483647d0)then
    write(binfile,'(a8,i9.9,a7)')'data/bin',int(time*0.01d0),'00s.dat'
   end if
  elseif(outstyle==2)then
   write(binfile,'(a8,i8.8,a4)')'data/bin',tn,'.dat'   
  end if

  open(unit=10,file=binfile,status='replace',form='unformatted')

  write(10)tn,time,mc(is-1)!,de_dt,domega_dt
  write(10) d (is:ie,js:je,ks:ke), &
            v1(is:ie,js:je,ks:ke), &
            v2(is:ie,js:je,ks:ke), &
            v3(is:ie,js:je,ks:ke), &
            e (is:ie,js:je,ks:ke)
  if(gravswitch>=2)write(10)grvphi(gis:gie,gjs:gje,gks:gke)
  if(gravswitch==3)then
   write(10)grvphiold(gis:gie,gjs:gje,gks:gke), &
            dt_old
  end if
  if(compswitch>=2)write(10)spc(1:spn,is:ie,js:je,ks:ke)
  if(mag_on)then
   write(10) b1(is:ie,js:je,ks:ke), &
             b2(is:ie,js:je,ks:ke), &
             b3(is:ie,js:je,ks:ke), &
             phi(is:ie,js:je,ks:ke)
  end if


  close(10)

  if(include_particles)then
!ptcfile----------------------------------------------------------------
  if(outstyle==1)then
   write(ptcfile,'(a8,i11.11,a5)')'data/ptc',int(time),'s.dat'
   if(time>2147483647d0)then
    write(ptcfile,'(a8,i9.9,a7)')'data/ptc',int(time*0.01d0),'00s.dat'
   end if
  elseif(outstyle==2)then
   write(ptcfile,'(a8,i8.8,a4)')'data/ptc',tn,'.dat'
  end if

  if(crdnt==1.and.je==1)then
   do n = 1, np
    extflag = .false.
    do k = ks, ke
     if(ptcx(2,n)>=xi3(k-1))then;if(ptcx(2,n)<xi3(k))then
      do i = is, ie
       if(ptcx(1,n)>=xi1(i-1))then;if(ptcx(1,n)<xi1(i))then
        if(e(i,js,k)+grvphi(i,js,k)*d(i,js,k)>0d0)then
         ptci(2,n) = 1
        else
         ptci(2,n) = 0
        end if
        extflag=.true.
        exit
       end if;end if
      end do
      if(extflag)exit
     end if;end if
    end do
   end do
  elseif(crdnt==2.and.ke==1)then
   do n = 1, np
    extflag = .false.
    pr = sqrt(ptcx(1,n)**2d0+ptcx(2,n)**2d0)
    pt = acos(ptcx(2,n)/pr)
    do j = js, je
     if(pt>=xi2(j-1))then;if(pt<xi2(j))then
      do i = is, ie
       if(pr>=xi1(i-1))then;if(pr<xi1(i))then
        if(e(i,j,ks)+grvphi(i,j,ks)*d(i,j,ks)>0d0)then
         ptci(2,n) = 1
        else
         ptci(2,n) = 0
        end if
        extflag=.true.
        exit
       end if;end if
      end do
      if(extflag)exit
     end if;end if
    end do
   end do
  end if

  open(unit=70,file = ptcfile, status='replace')

  write(70,'(a,i7,a,1PE12.4e2,2(a5,i9))')&
       '#tn =',tn,'  time= ',time,'np= ',np,'npl=',npl
  write(70,'(a7,2a3,3a15)')'label','ej','ub','mass','x1','x3'

  do i = 1, np
   write(70,'(i7,2i3,3(1PE15.7e2))')ptci(0:2,i),ptcx(0:2,i)
  end do

  close(70)

!bptfile----------------------------------------------------------------
  if(outstyle==1)then
   write(bptfile,'(a8,i11.11,a5)')'data/bpt',int(time),'s.dat'
   if(time>2147483647d0)then
    write(bptfile,'(a8,i9.9,a7)')'data/bpt',int(time*0.01d0),'00s.dat'
   end if
  elseif(outstyle==2)then
   write(bptfile,'(a8,i8.8,a4)')'data/bpt',tn,'.dat'
  end if

  open(unit=80,file = bptfile, status='replace',form='unformatted')

  write(80)np,npl
  write(80)ptci(0:2,1:np),ptcx(0:2,1:np),ptc_in(1:jmax)

  close(80)

  end if

!remmassfile---------------------------------------------------------------
  if(tn==0)then
   open(unit=60,file='data/remmass.dat',status='replace')
   inimass = 0d0
   do k = ks, ke ; do j = js, je ; do i = is, ie
    if(e(i,j,k)+d(i,j,k)*grvphi(i,j,k)<=0d0)then
     inimass = inimass + dvol(i,j,k)*d(i,j,k) + coremass
    end if
   end do ; end do ; end do
   write(60,'(a,1PE16.8e2,a15,1PE16.8e2)') &
    'initial_env_mass = ',inimass/msun,'core_mass = ',coremass/msun
   write(60,'(7a16)')'time','remmass','belmass','v_kick','v_kickb'
  else
   open(unit=60,file='data/remmass.dat',status='old',position='append')
  end if

  remmass = coremass ; belmass = coremass
  centre_of_mass = 0d0 ; centre_of_massb = 0d0
  vel = 0d0 ; velb = 0d0
  do k = ks, ke
   do j = js, je
    do i = is, ie
     if(e(i,j,k)+d(i,j,k)*(grvphi(i,j,k)+extgrv(i,j,k))<=0d0)then
      remmass = remmass + d(i,j,k)*dvol(i,j,k)
      vel = vel + d(i,j,k)*dvol(i,j,k)*v3(i,j,k)
     end if
     if(v1(i,j,k)**2d0+v3(i,j,k)**2d0+(grvphi(i,j,k)+extgrv(i,j,k))<=0d0)then
      belmass = belmass + d(i,j,k)*dvol(i,j,k)
      velb = velb + d(i,j,k)*dvol(i,j,k)*v3(i,j,k)
     end if
    end do
   end do
  end do
  write(60,'(7(1PE16.8e2))') time, remmass/msun, belmass/msun, vel/remmass, velb/belmass
  call flush(60)

return
end subroutine output
