!{\src2tex{textfont=tt}}
!!****f* ABINIT/d3sym
!! NAME
!! d3sym
!!
!!
!! FUNCTION
!! Given a set of calculated elements of the 3DTE matrix,
!! build (nearly) all the other matrix elements that can be build using
!! symmetries.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  indsym(4,nsym,natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  nsym=number of space group symmetries
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  Input/Output
!!  blkflg(3,mpert,3,mpert,3,mpert)= matrix that indicates if an
!!   element of d3 is available (1 if available, 0 otherwise)
!!  d3(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!
!! PARENTS
!!      nonlinear,rdddb9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine d3sym(blkflg,d3,indsym,mpert,natom,nsym,&
& symrec,symrel)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd3sym'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3(2,3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: found,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,idisy1,idisy2,idisy3
 integer :: ipesy1,ipesy2,ipesy3,isym,ithree
 real(dp) :: sumi,sumr
!arrays
 integer :: sym1(3,3),sym2(3,3),sym3(3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)'d3sym : enter'
!do i1dir = 1, 3
!do i2dir = 1, 3
!do i3dir = 1, 3
!write(std_out,*)i1dir,i2dir,i3dir,&
!&   blkflg(i1dir,natom+2,i2dir,natom+2,i3dir,natom+2)
!end do
!end do
!end do
!stop
!ENDDEBUG

!First, take into account the permuations symmetry of
!(i1pert,i1dir) and (i3pert,i3dir)

 do i1pert = 1, mpert
   do i2pert = 1, mpert

     do i3pert = 1, mpert

       do i1dir = 1, 3
         do i2dir = 1, 3
           do i3dir = 1, 3

             if ((blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1).and. &
&             (blkflg(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert)/=1)) then

               d3(:,i3dir,i3pert,i2dir,i2pert,i1dir,i1pert) = &
               d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)

               blkflg(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert) = 1

             end if

           end do
         end do
       end do

     end do
   end do
 end do



!Big Big Loop : symmetrize three times, because
!of some cases in which one element is not yet available
!at the first pass, and even at the second one !

 do ithree=1,3

!  Loop over perturbations

   do i1pert = 1, mpert
     do i2pert = 1, mpert
       do i3pert = 1, mpert

         do i1dir = 1, 3
           do i2dir = 1, 3
             do i3dir = 1, 3

!              Will get element (idir1,ipert1,idir2,ipert2)
!              so this element should not yet be present ...

               if(blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=1)then


                 d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 0_dp

                 do isym = 1, nsym

                   found = 1

                   if (i1pert <= natom) then
                     ipesy1 = indsym(4,isym,i1pert)
                     sym1(:,:) = symrec(:,:,isym)
                   else if (i1pert == natom + 2) then
                     ipesy1 = i1pert
                     sym1(:,:) = symrel(:,:,isym)
                   else
                     found = 0
                   end if

                   if (i2pert <= natom) then
                     ipesy2 = indsym(4,isym,i2pert)
                     sym2(:,:) = symrec(:,:,isym)
                   else if (i2pert == natom + 2) then
                     ipesy2 = i2pert
                     sym2(:,:) = symrel(:,:,isym)
                   else
                     found = 0
                   end if

                   if (i3pert <= natom) then
                     ipesy3 = indsym(4,isym,i3pert)
                     sym3(:,:) = symrec(:,:,isym)
                   else if (i3pert == natom + 2) then
                     ipesy3 = i3pert
                     sym3(:,:) = symrel(:,:,isym)
                   else
                     found = 0
                   end if

                   sumr = 0_dp ; sumi = 0_dp;
                   do idisy1 = 1, 3
                     do idisy2 = 1, 3
                       do idisy3 = 1, 3

                         if ((sym1(i1dir,idisy1) /=0).and.(sym2(i2dir,idisy2) /=0) &
&                         .and.(sym3(i3dir,idisy3) /=0)) then

                           if (blkflg(idisy1,ipesy1,idisy2,ipesy2,idisy3,ipesy3) == 1) then

                             sumr = sumr + sym1(i1dir,idisy1)*sym2(i2dir,idisy2)*&
&                             sym3(i3dir,idisy3)*d3(1,idisy1,ipesy1,idisy2,ipesy2,idisy3,ipesy3)
                             sumi = sumi + sym1(i1dir,idisy1)*sym2(i2dir,idisy2)*&
&                             sym3(i3dir,idisy3)*d3(2,idisy1,ipesy1,idisy2,ipesy2,idisy3,ipesy3)

                           else

                             found = 0

                           end if

                         end if

                       end do
                     end do
                   end do

                   if (found == 1) then
                     d3(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumr
                     d3(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumi
                     blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
                   end if

                 end do  ! isym

               end if  ! blkflg


!              Close loop over perturbations

             end do
           end do
         end do
       end do
     end do
   end do

 end do  ! close loop over ithree



end subroutine d3sym
!!***
