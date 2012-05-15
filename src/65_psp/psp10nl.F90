!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp10nl
!! NAME
!! psp10nl
!!
!! FUNCTION
!! Hartwigsen-Goedecker-Hutter nonlocal pseudopotential (from preprint of 1998).
!! Uses Gaussians for fully nonlocal form, analytic expressions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, FRD, XG, GMR, PT, SC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  hij(0:lmax,3,3)=factor defining strength of (max 3) projectors for each
!!   angular momentum channel l among 0, 1, ..., lmax
!!  lmax=maximum angular momentum
!!  mproj=maximum number of projectors in any channel
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=number of grid points for qgrid
!!  nproj(1:lmax+1)=number of projectors in any channel
!!  qgrid(mqgrid)=array of |G| values
!!  rr(0:lmax)=core radius for each 0<l<lmax channel (bohr)
!!
!! OUTPUT
!!  ekb(mpsang,mproj)=Kleinman-Bylander energies
!!  ffspl(mqgrid,2,mpssang,mproj)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projectors
!!
!! PARENTS
!!      psp10in
!!
!! CHILDREN
!!      leave_new,spline,wrtout,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp10nl(ekb,ffspl,hij,lmax,mproj,mpsang,mqgrid,nproj,qgrid,rr)

 use m_profiling

 use defs_basis
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp10nl'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,mproj,mpsang,mqgrid
!arrays
 integer,intent(in) :: nproj(mpsang)
 real(dp),intent(in) :: hij(0:lmax,3,3),qgrid(mqgrid),rr(0:lmax)
 real(dp),intent(out) :: ekb(mpsang,mproj),ffspl(mqgrid,2,mpsang,mproj)

!Local variables-------------------------------
!scalars
 integer :: info,ipack,iproj,iqgrid,jproj,ll,numproj
 real(dp) :: qmax,rrl
 character(len=500) :: message
 character :: jobz,uplo
!arrays
 real(dp) :: ap(2,9),rwork1(9),work1(2,9),ww(3),yp1j(3),ypnj(3)
 real(dp),allocatable :: ppspl(:,:,:,:),uu(:,:),work(:),zz(:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' psp10nl : enter '
!stop
!ENDDEBUG

 ABI_ALLOCATE(ppspl,(mqgrid,2,mpsang,mproj))
 ABI_ALLOCATE(work,(mqgrid))

 qmax=qgrid(mqgrid)
 jobz='v'
 uplo='u'
 ekb(:,:)=zero

 lloop: do ll=0,lmax
   ap(:,:)=zero
   numproj=nproj(ll+1)

!  Fill up the matrix in packed storage
   prjloop: do jproj=1,numproj
     priloop: do iproj=1,jproj
       ipack=iproj+(jproj-1)*jproj/2; if(mod((jproj-1)*jproj,2)/=0)stop "odd"
       ap(1,ipack)=hij(ll,iproj,jproj)
     end do priloop
   end do prjloop

   if(numproj/=0)then

     ABI_ALLOCATE(uu,(numproj,numproj))
     ABI_ALLOCATE(zz,(2,numproj,numproj))

     if (numproj > 1) then
       call ZHPEV(jobz,uplo,numproj,ap,ww,zz,numproj,work1,rwork1,info)
       uu(:,:)=zz(1,:,:)
     else
       ww(1)=hij(ll,1,1)
       uu(1,1)=one
     end if

!    Initialization of ekb, and spline fitting

     if (ll==0) then ! s channel

       rrl=rr(0)
       do iproj=1,numproj
         ekb(1,iproj)=ww(iproj)*32.d0*(rrl**3)*(pi**2.5d0)/(4.d0*pi)**2
         if(iproj==1)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,1,1)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2)
           end do
           yp1j(1)=zero
           ypnj(1)=-(two_pi*rrl)**2*qmax*exp(-0.5d0*(two_pi*qmax*rrl)**2)
         else if(iproj==2)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,1,2)=2.0d0/sqrt(15.0d0)     &
&             *exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) &
&             *( 3.d0-(two_pi*qgrid(iqgrid)*rrl)**2 )
           end do
           yp1j(2)=zero
           ypnj(2)=2.0d0/sqrt(15.0d0)*(two_pi*rrl)**2*qmax &
&           *exp(-0.5d0*(two_pi*qmax*rrl)**2) * (-5.d0+(two_pi*qmax*rrl)**2)
         else if(iproj==3)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,1,3)=(4.0d0/3.0d0)/sqrt(105.0d0)*&
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * &
&             (15.0d0-10.0d0*(two_pi*qgrid(iqgrid)*rrl)**2 + &
&             (two_pi*qgrid(iqgrid)*rrl)**4)
           end do
           yp1j(3)=zero
           ypnj(3)=(4.0d0/3.0d0)/sqrt(105.0d0)*exp(-0.5d0*(two_pi*qmax*rrl)**2) * &
&           (two_pi*rrl)**2*qmax*(-35.0d0+14d0*(two_pi*qmax*rrl)**2-(two_pi*qmax*rrl)**4)
         end if
         call spline(qgrid,ppspl(:,1,1,iproj),mqgrid,&
&         yp1j(iproj),ypnj(iproj),ppspl(:,2,1,iproj))
       end do

     else if (ll==1) then ! p channel

       rrl=rr(1)
       do iproj=1,numproj
         ekb(2,iproj)=ww(iproj)*64.d0*(rrl**5)*(pi**2.5d0)/(4.d0*pi)**2
         if(iproj==1)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,2,1)=(1.0d0/sqrt(3.0d0))* &
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * (two_pi*qgrid(iqgrid))
           end do
           yp1j(1)=two_pi*(1.0d0/sqrt(3.0d0))
           ypnj(1)=-two_pi*((two_pi*qmax*rrl)**2-1.d0)*exp(-0.5d0*(two_pi*qmax*rrl)**2)*&
&           (1.0d0/sqrt(3.0d0))
         else if(iproj==2)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,2,2)=(2.0d0/sqrt(105.0d0))* &
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * &
&             (two_pi*qgrid(iqgrid))*(5.0d0-(two_pi*qgrid(iqgrid)*rrl)**2)
           end do
           yp1j(2)=(5.0d0*two_pi)*(2.0d0/sqrt(105.0d0))
           ypnj(2)=(2.0d0/sqrt(105.0d0))*two_pi*exp(-0.5d0*(two_pi*qmax*rrl)**2)* &
&           (-8*(two_pi*qmax*rrl)**2 + (two_pi*qmax*rrl)**4 + 5.0d0)
         else if(iproj==3)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,2,3)=(4.0d0/3.0d0)/sqrt(1155d0)*&
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * &
&             (two_pi*qgrid(iqgrid))*&
&             (35.0d0-14.0d0*(two_pi*qgrid(iqgrid)*rrl)**2+(two_pi*qgrid(iqgrid)*rrl)**4)
           end do
           yp1j(3)=(35.0d0*two_pi)*(4.0d0/3.0d0)/sqrt(1155.0d0)
           ypnj(3)=(4.0d0/3.0d0)/sqrt(1155.0d0)*two_pi*exp(-0.5d0*(two_pi*qmax*rrl)**2)* &
&           (35.0d0-77.0d0*(two_pi*qmax*rrl)**2+19.0d0*(two_pi*qmax*rrl)**4 - &
&           (two_pi*qmax*rrl)**6)
         end if
         call spline(qgrid,ppspl(:,1,2,iproj),mqgrid,&
&         yp1j(iproj),ypnj(iproj),ppspl(:,2,2,iproj))
       end do

     else if (ll==2) then ! d channel

!      If there is a third projector. Warning : only two projectors are allowed.
       if ( numproj>2 ) then
         write(message, '(a,a,a,a,a,a)' ) ch10,&
&         ' psp10nl : ERROR -',ch10,&
&         '  only two d-projectors are allowed ',ch10,&
&         '  Action : check your pseudopotential file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if

       rrl=rr(2)
       do iproj=1,numproj
         ekb(3,iproj)=ww(iproj)*128.d0*(rrl**7)*(pi**2.5d0)/(4.d0*pi)**2
         if(iproj==1)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,3,1)=(1.0d0/sqrt(15.0d0))* &
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * (two_pi*qgrid(iqgrid))**2
           end do
           yp1j(1)=zero
           ypnj(1)=(1.0d0/sqrt(15.0d0))*(two_pi**2)*&
&           exp(-0.5d0*(two_pi*qmax*rrl)**2)*qmax*(2d0-(two_pi*qmax*rrl)**2)
         else if(iproj==2)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,3,2)=(2.0d0/3.0d0)/sqrt(105.0d0)* &
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * &
&             ((two_pi*qgrid(iqgrid))**2)*(7.0d0-(two_pi*qgrid(iqgrid)*rrl)**2)
           end do
           yp1j(2)=zero
           ypnj(2)=(2.0d0/3.0d0)/sqrt(105.0d0)*exp(-0.5d0*(two_pi*qmax*rrl)**2)* &
&           qmax*(two_pi**2)*( (two_pi*qmax*rrl)**4 - 11.0d0*(two_pi*qmax*rrl)**2 + 14.0d0)
         end if
         call spline(qgrid,ppspl(:,1,3,iproj),mqgrid,&
&         yp1j(iproj),ypnj(iproj),ppspl(:,2,3,iproj))
       end do

     else if (ll==3) then ! f channel

!      If there is a second projector. Warning : only one projector is allowed.
       if ( numproj>1 ) then
         write(message, '(a,a,a,a,a,a)' ) ch10,&
&         ' psp10nl : ERROR -',ch10,&
&         '  only one f-projector is allowed ',ch10,&
&         '  Action : check your pseudopotential file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if

       rrl=rr(3)
       ekb(4,1)=ww(1)*(256.0d0/105.0d0)*(rrl**9)*(pi**2.5d0)/(4.d0*pi)**2
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,4,1)=(two_pi*qgrid(iqgrid))**3* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2)
       end do
!      Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
       yp1j(1)=zero
       ypnj(1)=(two_pi**3)*qmax**2*exp(-0.5d0*(two_pi*qmax*rrl)**2)*&
&       (3.0d0-(two_pi*qmax*rrl)**2)
!      Fit spline to get second derivatives by spline fit
       call spline(qgrid,ppspl(:,1,4,1),mqgrid,&
&       yp1j(1),ypnj(1),ppspl(:,2,4,1))

     else; stop "lmax>3?"
     end if

!    Linear combination using the eigenvectors
     ffspl(:,:,ll+1,:)=zero
     do jproj=1,numproj
       do iproj=1,numproj
         do iqgrid=1,mqgrid
           ffspl(iqgrid,1:2,ll+1,jproj)=ffspl(iqgrid,1:2,ll+1,jproj) &
&           +uu(iproj,jproj)*ppspl(iqgrid,1:2,ll+1,iproj)
         end do
       end do
     end do

     ABI_DEALLOCATE(uu)
     ABI_DEALLOCATE(zz)

!    End condition on numproj(/=0)
   end if

 end do lloop

!DEBUG
!write(std_out,*)' psp10nl : after lloop '
!stop
!ENDDEBUG

 ABI_DEALLOCATE(ppspl)
 ABI_DEALLOCATE(work)

!DEBUG
!write(std_out,*)' psp10nl : exit '
!stop
!ENDDEBUG

end subroutine psp10nl
!!***
