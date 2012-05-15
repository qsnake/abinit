!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtchi
!!
!! NAME
!! dtchi
!!
!! FUNCTION
!! Reads the non-linear optical susceptibility tensor and the
!! first-order change in the linear dielectric susceptibility
!! induced by an atomic displacement in the
!! Gamma Block coming from the Derivative Data Base
!! (third-order derivatives).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert)= matrix of third-order energies
!! natom= number of atoms in unit cell
!! mpert =maximum number of ipert
!! ramansr= if /= 0, impose sum rule on first-order derivatives
!!                   of the electronic susceptibility with respect
!!                   to atomic displacements
!!
!!
!! OUTPUT
!! dchide(3,3,3) = non-linear optical coefficients
!! dchidt(natom,3,3,3) = first-order change of the electronic dielectric
!!   tensor induced by an individual atomic displacement
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dtchi(blkval,dchide,dchidt,mpert,natom,ramansr)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtchi'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ramansr
!arrays
 real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(out) :: dchide(3,3,3),dchidt(natom,3,3,3)

!Local variables -------------------------
!scalars
 integer :: depl,elfd1,elfd2,elfd3,iatom,ivoigt
 real(dp) :: wttot
!arrays
 integer :: voigtindex(6,2)
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert),dvoigt(3,6),sumrule(3,3,3)
 real(dp) :: wghtat(natom)

! *********************************************************************

 d3cart(1,:,:,:,:,:,:) = reshape(blkval(1,:),&
& shape = (/3,mpert,3,mpert,3,mpert/))
 d3cart(2,:,:,:,:,:,:) = reshape(blkval(2,:),&
& shape = (/3,mpert,3,mpert,3,mpert/))


!Extraction of non-linear optical coefficients

 do elfd1 = 1,3
   do elfd2 = 1,3
     do elfd3 = 1,3
       dchide(elfd1,elfd2,elfd3) = &
&       d3cart(1,elfd1,natom+2,elfd2,natom+2,elfd3,natom+2)
     end do
   end do
 end do

!Transform to Voigt notations

 voigtindex(:,1) = (/1,2,3,2,1,1/)
 voigtindex(:,2) = (/1,2,3,3,3,2/)
 do ivoigt = 1, 6
   elfd2 = voigtindex(ivoigt,1)
   elfd3 = voigtindex(ivoigt,2)
   do elfd1 = 1, 3
     dvoigt(elfd1,ivoigt) = &
&     0.5_dp*(dchide(elfd1,elfd2,elfd3) + dchide(elfd1,elfd3,elfd2))
   end do
 end do

!Transform to pm/V

 dvoigt(:,:) = dvoigt(:,:)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb


!Extraction of $\frac{d \chi}{d \tau}$

 do iatom = 1, natom
   do depl = 1,3
     do elfd1 = 1,3
       do elfd2 = 1,3
         dchidt(iatom,depl,elfd1,elfd2) = &
&         d3cart(1,depl,iatom,elfd1,natom+2,elfd2,natom+2)
       end do
     end do
   end do
 end do

 wghtat(:) = 0._dp
 if (ramansr == 1) then

   wghtat(:) = 1._dp/dble(natom)

 else if (ramansr == 2) then

   wttot = 0._dp
   do iatom = 1, natom

     do depl = 1,3
       do elfd1 = 1,3
         do elfd2 = 1,3
           wghtat(iatom) = wghtat(iatom) + abs(dchidt(iatom,depl,elfd1,elfd2))
         end do
       end do
     end do
     wttot = wttot + wghtat(iatom)

   end do

   wghtat(:) = wghtat(:)/wttot

 end if


 write(ab_out,*)ch10
 write(ab_out,*)'Non-linear optical coefficients d (pm/V)'
 write(ab_out,'(6f12.6)')dvoigt(1,:)
 write(ab_out,'(6f12.6)')dvoigt(2,:)
 write(ab_out,'(6f12.6)')dvoigt(3,:)

 if (ramansr /= 0) then

   write(ab_out,*)ch10
   write(ab_out,*)'The violation of the Raman sum rule'
   write(ab_out,*)'by the first-order electronic dielectric tensors ',&
&   'is as follows'
   write(ab_out,*)'    atom'
   write(ab_out,*)' displacement'

   sumrule(:,:,:) = 0._dp
   do elfd2 = 1,3
     do elfd1 = 1,3
       do depl = 1,3
         do iatom = 1, natom
           sumrule(depl,elfd1,elfd2) = sumrule(depl,elfd1,elfd2) + &
&           dchidt(iatom,depl,elfd1,elfd2)
         end do
         do iatom = 1, natom
           dchidt(iatom,depl,elfd1,elfd2) = dchidt(iatom,depl,elfd1,elfd2) - &
&           wghtat(iatom)*sumrule(depl,elfd1,elfd2)
         end do
       end do
     end do
   end do

   do depl = 1,3
     write(ab_out,'(6x,i2,3(3x,f16.9))') depl,sumrule(depl,1,1:3)
     write(ab_out,'(8x,3(3x,f16.9))') sumrule(depl,2,1:3)
     write(ab_out,'(8x,3(3x,f16.9))') sumrule(depl,3,1:3)
     write(ab_out,*)
   end do

 end if    ! ramansr

 write(ab_out,*)ch10
 write(ab_out,*)' First-order change in the electronic dielectric '
 write(ab_out,*)' susceptibility tensor (Bohr^-1)'
 write(ab_out,*)' induced by an atomic displacement'
 if (ramansr /= 0) then
   write(ab_out,*)' (after imposing the sum over all atoms to vanish)'
 end if
 write(ab_out,*)'  atom  displacement'
 do iatom = 1,natom
   do depl = 1,3

     write(ab_out,'(1x,i4,9x,i2,3(3x,f16.9))')iatom,depl,&
&     dchidt(iatom,depl,1,:)
     write(ab_out,'(16x,3(3x,f16.9))')&
&     dchidt(iatom,depl,2,:)
     write(ab_out,'(16x,3(3x,f16.9))')&
&     dchidt(iatom,depl,3,:)

   end do

   write(ab_out,*)

 end do

!DEBUG
!sumrule(:,:,:) = 0._dp
!do elfd2 = 1,3
!do elfd1 = 1,3
!do depl = 1,3
!do iatom = 1, natom
!sumrule(depl,elfd1,elfd2) = sumrule(depl,elfd1,elfd2) + &
!&     dchidt(iatom,depl,elfd1,elfd2)
!end do
!end do
!end do
!end do
!do depl = 1,3
!write(ab_out,'(6x,i2,3(3x,f16.9))') depl,sumrule(depl,1,1:3)
!write(ab_out,'(8x,3(3x,f16.9))') sumrule(depl,2,1:3)
!write(ab_out,'(8x,3(3x,f16.9))') sumrule(depl,3,1:3)
!write(ab_out,*)
!end do
!ENDEBUG




end subroutine dtchi
!!***
