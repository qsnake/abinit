!{\src2tex{textfont=tt}}
!!****f* ABINIT/d3output
!! NAME
!! d3output
!!
!!
!! FUNCTION
!! Write the matrix of third-order derivatives to the output file
!! and the DDB
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
!!  mpert =maximum number of ipert
!!  unddb = unit number for DDB output
!!
!! OUTPUT
!!
!! NOTES
!!  d3 holds the third-order derivatives before computing
!!  the permutations of the perturbations.
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      create_ddb_blk,destroy_ddb_blk,nullify_ddb_blk,write_blok8,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine d3output(blkflg,d3,mband,mpert,nkpt,unddb)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_ddb_blk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd3output'
 use interfaces_14_hidewrite
 use interfaces_72_response, except_this_one => d3output
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mband,mpert,nkpt,unddb
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: choice,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,index,msize
 character(len=500) :: message
 type(ddb_blk_type), pointer :: ddb_blk
!arrays

!*************************************************************************




 ABI_ALLOCATE(ddb_blk,)
 call nullify_ddb_blk(ddb_blk)

 msize = 27*mpert*mpert*mpert
 call create_ddb_blk(msize, 1, ddb_blk)

 choice = 2

 ddb_blk%typ = 3
 ddb_blk%nrm = one
 ddb_blk%qpt = zero   ! this has to be changed in case anharmonic
!force constants have been computed


!Write blok of third-order derivatives to ouput file

 write(message,'(a,a,a,a,a)')ch10,&
& ' Matrix of third-order derivatives (reduced coordinates)',ch10,&
& ' before computing the permutations of the perturbations',ch10
 call wrtout(ab_out,message,'COLL')

 write(ab_out,*)'    j1       j2       j3              matrix element'
 write(ab_out,*)' dir pert dir pert dir pert           real part           imaginary part'

 do i1pert=1,mpert
   do i1dir=1,3

     do i2pert=1,mpert
       do i2dir=1,3

         do i3pert=1,mpert
           do i3dir=1,3

             index = i1dir + &
&             3*((i1pert-1)+mpert*((i2dir-1) + &
&             3*((i2pert-1)+mpert*((i3dir-1) + 3*(i3pert-1)))))
             ddb_blk%flg(index,1) = blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
             ddb_blk%val(:,index,1)= d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)

             if (blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=0) then

               write(ab_out,'(3(i4,i5),2f22.10)')&
&               i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,&
&               d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)

             end if

           end do
         end do

       end do
     end do

   end do
 end do


!Write blok of third-order derivatives to DDB
 call write_blok8(ddb_blk,1,choice,mband,mpert,&
& msize,nkpt,unddb)

 call destroy_ddb_blk(ddb_blk)
 ABI_DEALLOCATE(ddb_blk)

end subroutine d3output
!!***
