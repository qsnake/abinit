!{\src2tex{textfont=tt}}
!!****f* ABINIT/ramansus
!!
!! NAME
!! ramansus
!!
!! FUNCTION
!! Compute the raman susceptibilities of zone-center phonons
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  d2cart = second order derivatives of the energy wrt all perturbations
!!  dchide(3,3,3) = non-linear optical coefficients from dtchi
!!  dchidt(natom,3,3,3) = first-order change of the electronic dielectric
!!     tensor induced by an individual atomic displacement
!!  displ = phonon mode atomic displacements
!!  mpert = maximum number of perturbations
!!  natom = number of atoms
!!  phfrq = phonon frequencies
!!  qphnrm=(described below)
!!  ucvol = unit cell volume
!!
!! OUTPUT
!!  qphon(3)= to be divided by qphnrm, give the phonon wavevector;
!!     if qphnrm==0.0_dp, then the wavevector is zero (Gamma point)
!!     and qphon gives the direction of
!!     the induced electric field; in the latter case, if qphon is
!!     zero, no non-analytical contribution is included.
!!  rsus
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ramansus(d2cart,dchide,dchidt,displ,mpert,&
& natom,phfrq,qphon,qphnrm,rsus,ucvol)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ramansus'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: mpert,natom
 real(dp),intent(in) :: qphnrm,ucvol
!arrays
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert),dchide(3,3,3)
 real(dp),intent(in) :: dchidt(natom,3,3,3),displ(2,3*natom,3*natom)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(inout) :: qphon(3)
 real(dp),intent(out) :: rsus(3*natom,3,3)

!Local variables-------------------------------
!scalars
 integer :: analyt,i1,i1dir,i1pert,i2dir,iatom,idir,imode
 real(dp) :: epsq,fac,qphon2
 character(len=500) :: message
!arrays
 real(dp) :: dijk_q(3,3)
 real(dp),allocatable :: zeff(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)'ramansus : enter' ; stop
!ENDDEBUG

 ABI_ALLOCATE(zeff,(3,natom))

 rsus(:,:,:) = 0._dp
 epsq        = 0._dp
 zeff(:,:)   = 0._dp
 dijk_q(:,:) = 0._dp

!Determine the analyticity of the matrix.
 analyt=1
 if(abs(qphnrm)<tol8)analyt=0
 if(abs(qphon(1))<tol8.and.abs(qphon(2))<tol8.and.abs(qphon(3))<tol8)&
& analyt=1



!In the case the non-analyticity is required :
 if(analyt == 0) then

!  Normalize the limiting direction
   qphon2=qphon(1)**2+qphon(2)**2+qphon(3)**2
   qphon(1)=qphon(1)/sqrt(qphon2)
   qphon(2)=qphon(2)/sqrt(qphon2)
   qphon(3)=qphon(3)/sqrt(qphon2)

!  Get the dielectric constant for the limiting direction
   epsq= 0._dp
   do i1dir=1,3
     do i2dir=1,3
       epsq=epsq+qphon(i1dir)*qphon(i2dir)*&
&       d2cart(1,i1dir,natom+2,i2dir,natom+2)
     end do
   end do

!  Check if epsq > 0
   if (epsq < tol8) then
     write(message,'(a,a,a,a,e13.6)')ch10,&
&     ' ramansus : BUG -',ch10,&
&     '  The value of epsq must be > 0 while it is found to be',epsq
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if


!  Get the effective charges for the limiting direction
   do i1dir=1,3
     do i1pert=1,natom
       zeff(i1dir,i1pert)=0.0_dp
       do i2dir=1,3
         zeff(i1dir,i1pert)=zeff(i1dir,i1pert)+qphon(i2dir)*&
&         d2cart(1,i1dir,i1pert,i2dir,natom+2)
       end do
     end do
   end do

!  Get the NLO tensor for the limiting direction
!  $\sum_{k} d_{ijk} \cdot q_k$

   dijk_q(:,:) = 0._dp
   do i1dir = 1, 3
     do i2dir = 1, 3
       do idir = 1, 3
         dijk_q(i1dir,i2dir) = dijk_q(i1dir,i2dir) + &
&         dchide(i1dir,i2dir,idir)*qphon(idir)
       end do
     end do
   end do

   fac = 16._dp*pi/(ucvol*epsq)
   do imode = 1, 3*natom
     do iatom = 1, natom
       do idir = 1, 3
         i1=idir + (iatom - 1)*3
         rsus(imode,:,:) = rsus(imode,:,:) + &
&         (dchidt(iatom,idir,:,:) - fac*zeff(idir,iatom)*dijk_q(:,:))* &
&         displ(1,i1,imode)
       end do  ! disp
     end do  ! iatom
   end do  ! imode
   rsus(:,:,:) = rsus(:,:,:)*sqrt(ucvol)

 else

   do imode = 1, 3*natom
     do iatom = 1, natom
       do idir = 1, 3
         i1=idir + (iatom - 1)*3
         rsus(imode,:,:) = rsus(imode,:,:) + &
&         dchidt(iatom,idir,:,:)*displ(1,i1,imode)
       end do  ! disp
     end do  ! iatom
   end do  ! imode
   rsus(:,:,:) = rsus(:,:,:)*sqrt(ucvol)

 end if      ! analyt == 0

 if (analyt == 0) then

   write(ab_out,*) ch10
   write(ab_out, '(a,/,a,3f9.5)' )&
&   ' Raman susceptibility of zone-center phonons, with non-analyticity in the',&
&   '  direction (cartesian coordinates)',qphon(1:3)+tol10
   write(ab_out,'(a)')&
&   ' -----------------------------------------------------------------------'
   write(ab_out,*) ch10

 else

   write(ab_out,*) ch10
   write(ab_out,*)' Raman susceptibilities of transverse zone-center phonon modes'
   write(ab_out,*)' -------------------------------------------------------------'
   write(ab_out,*) ch10

 end if

 do imode = 1, 3*natom
   write(ab_out,'(a4,i3,2x,a2,f7.2,a6)')' Mode',imode,&
&   ' (',phfrq(imode)*Ha_cmm1,' cm-1)'
   do idir = 1,3
     write(ab_out,'(a,4x,3(f16.9,2x))')';',rsus(imode,idir,:)
   end do
   write(ab_out,*)
 end do

 ABI_DEALLOCATE(zeff)

end subroutine ramansus
!!***
