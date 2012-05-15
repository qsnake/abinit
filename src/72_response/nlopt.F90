!{\src2tex{textfont=tt}}
!!****f* ABINIT/nlopt
!! NAME
!! nlopt
!!
!!
!! FUNCTION
!! Output of all quantities related to third-order derivatives
!! of the energy.
!! Compute the permutations of the three perturbations, then
!! write out the whole matrix of third order derivatives
!! in reduced coordinates. Finally, compute the non-linear optical
!! susceptibility d and the first-order change in the dielectric
!! susceptibility tensor induced by an atomic displacement.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blkflg(3,mpert,3,mpert,3,mpert)= ( 1 if the element of the 3dte
!!   has been calculated ; 0 otherwise )
!!  d3(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!  gprimd(3,3)=dimensional primitive translations for
!!              reciprocal space(bohr^-1)
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  ucvol=unit cell volume (bohr^3)
!!
!! OUTPUT
!! carflg(3,mpert,3,mpert,3,mpert)=1 if the element of d3cart has been
!!   calculated, 0 otherwise
!! d3cart(2,3,mpert,3,mpert,3,mpert)=matrix of third-order energy
!!   derivatives in cartesian coordinates
!!
!!
!! PARENTS
!!      nonlinear,rdddb9
!!
!! CHILDREN
!!      cart39
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine nlopt(blkflg,carflg,d3,d3cart,gprimd,mpert,natom,rprimd,ucvol)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nlopt'
 use interfaces_72_response, except_this_one => nlopt
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
 integer,intent(out) :: carflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert),gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: d3cart(2,3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)

! *******************************************************************

!DEBUG
!write(std_out,*)'nlopt : enter'
!write(std_out,*)'ucvol = ',ucvol ; stop
!ENDDEBUG

!Compute the permutations of the perturbations

 d3cart(:,:,:,:,:,:,:) = 0._dp

 do i1pert = 1,mpert
   do i2pert = 1,mpert
     do i3pert = 1,mpert
       do i1dir=1,3
         do i2dir=1,3
           do i3dir=1,3

!            Check if all elements are available

             if ((blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=0).and. &
&             (blkflg(i1dir,i1pert,i3dir,i3pert,i2dir,i2pert)/=0).and. &
&             (blkflg(i2dir,i2pert,i1dir,i1pert,i3dir,i3pert)/=0).and. &
&             (blkflg(i2dir,i2pert,i3dir,i3pert,i1dir,i1pert)/=0).and. &
&             (blkflg(i3dir,i3pert,i1dir,i1pert,i2dir,i2pert)/=0).and. &
&             (blkflg(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert)/=0)) then

               d3cart(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = &
&               (  d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) + &
&               d3(:,i1dir,i1pert,i3dir,i3pert,i2dir,i2pert) + &
&               d3(:,i2dir,i2pert,i1dir,i1pert,i3dir,i3pert) + &
&               d3(:,i2dir,i2pert,i3dir,i3pert,i1dir,i1pert) + &
&               d3(:,i3dir,i3pert,i1dir,i1pert,i2dir,i2pert) + &
&               d3(:,i3dir,i3pert,i2dir,i2pert,i1dir,i1pert))*sixth

             end if
           end do
         end do
       end do
     end do
   end do
 end do


!Transform to cartesian coordinates

 carflg(:,:,:,:,:,:) = 0

 do i1pert = 1, mpert
   do i2pert = 1, mpert
     do i3pert = 1, mpert

       do i2dir = 1, 3
         do i3dir = 1, 3

           vec1(:) = d3cart(1,:,i1pert,i2dir,i2pert,i3dir,i3pert)
           flg1(:) = blkflg(:,i1pert,i2dir,i2pert,i3dir,i3pert)
           call cart39(flg1,flg2,gprimd,i1pert,natom,rprimd,vec1,vec2)
           d3cart(1,:,i1pert,i2dir,i2pert,i3dir,i3pert) = vec2(:)
           carflg(:,i1pert,i2dir,i2pert,i3dir,i3pert) = flg2(:)

         end do
       end do

       do i1dir = 1, 3
         do i3dir = 1, 3

           vec1(:) = d3cart(1,i1dir,i1pert,:,i2pert,i3dir,i3pert)
           flg1(:) = blkflg(i1dir,i1pert,:,i2pert,i3dir,i3pert)
           call cart39(flg1,flg2,gprimd,i2pert,natom,rprimd,vec1,vec2)
           d3cart(1,i1dir,i1pert,:,i2pert,i3dir,i3pert) = vec2(:)
           carflg(i1dir,i1pert,:,i2pert,i3dir,i3pert) = flg2(:)

         end do
       end do

       do i1dir = 1, 3
         do i2dir = 1, 3

           vec1(:) = d3cart(1,i1dir,i1pert,i2dir,i2pert,:,i3pert)
           flg1(:) = blkflg(i1dir,i1pert,i2dir,i2pert,:,i3pert)
           call cart39(flg1,flg2,gprimd,i3pert,natom,rprimd,vec1,vec2)
           d3cart(1,i1dir,i1pert,i2dir,i2pert,:,i3pert) = vec2(:)
           carflg(i1dir,i1pert,i2dir,i2pert,:,i3pert) = flg2(:)

         end do
       end do


     end do
   end do
 end do


!Compute non linear-optical coefficients d_ijk (atomic units)

 i1pert = natom+2
 d3cart(:,:,i1pert,:,i1pert,:,i1pert) = &
& -3._dp*d3cart(:,:,i1pert,:,i1pert,:,i1pert)/(ucvol*2._dp)


!Compute first-order change in the electronic dielectric
!susceptibility (Bohr^-1) induced by an atomic displacement

 d3cart(1:2,1:3,1:natom,1:3,natom + 2,1:3,natom + 2) = &
& -6._dp*d3cart(1:2,1:3,1:natom,1:3,natom + 2,1:3,natom + 2)/ucvol



end subroutine nlopt
!!***
