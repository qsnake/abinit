!{\src2tex{textfont=tt}}
!!****f* ABINIT/carttransf
!! NAME
!! carttransf
!!
!! FUNCTION
!! Transform a second-derivative matrix (EIG2D) from reduced
!! coordinates to cartesian coordinates. 
!! 
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (PB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  blkflg(msize,nblok)=
!!   ( 1 if the element of the dynamical matrix has been calculated ;
!!     0 otherwise )
!!  gprimd(3,3)=basis vector in the reciprocal space
!!  iqpt  = number of the Q-point currently used
!!  mband = maximal number of bands
!!  mpert = maximum number of ipert
!!  msize = size of the EIG2D arrays (3*mpert*3*mpert)
!!  natom = number of atom
!!  nblok = number of bloks in blkflg
!!  nkpt  = number of K-points
!!  rprimd(3,3) = basis vector in the real space
!!
!! OUTPUT
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  EIG2D matrix has been calculated correctly ; 0 otherwise )
!!
!! SIDE EFFECT
!! blkval2(2,msize,mband,nkpt)=Second order eigenvalues (EIG2D)
!! is transformed from reduced coordinates to cartesian coordinates
!!
!! PARENTS
!!      thmeig
!!
!! CHILDREN
!!      carteig2d
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine carttransf(blkflg,blkval2,carflg,gprimd,iqpt,mband,&
& mpert,msize,natom,nblok,nkpt,rprimd)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'carttransf'
 use interfaces_77_ddb, except_this_one => carttransf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,msize
 integer,intent(in) :: iqpt
 integer,intent(in) :: mpert,nblok
 integer,intent(inout) :: natom,nkpt
!arrays
 integer,intent(in) :: blkflg(msize,nblok)
 integer,intent(out) :: carflg(3,mpert,3,mpert)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(inout) :: blkval2(2,msize,mband,nkpt)
!Local variables-------------------------------
!scalars
integer :: iatom1,iatom2,iband,idir1,idir2,ikpt
integer :: index
!arrays
real(dp),allocatable :: blkflgtmp(:,:,:,:,:)
real(dp),allocatable :: blkval2tmp(:,:,:,:,:,:)
real(dp),allocatable :: d2cart(:,:,:,:,:)

! *********************************************************************

!Start by allocating local arrays
 ABI_ALLOCATE(blkflgtmp,(3,mpert,3,mpert,1))
 ABI_ALLOCATE(blkval2tmp,(2,3,mpert,3,mpert,1))
 ABI_ALLOCATE(d2cart,(2,3,mpert,3,mpert))

!Begin by formating the arrays to be compatible with cart29
!Then call cart29 to transform the arrays in cartesian coordinates
!Finally reformat the cartesian arrays in old format 
 do ikpt=1,nkpt
   do iband=1,mband

     do idir1=1,3
       do iatom1=1,mpert
         do idir2=1,3
           do iatom2=1,mpert
             index = idir1 + 3*((iatom1 - 1) + natom * ((idir2-1)+3*(iatom2-1)))
             blkflgtmp(idir1,iatom1,idir2,iatom2,1) = blkflg(index,iqpt)
             blkval2tmp(:,idir1,iatom1,idir2,iatom2,1) = blkval2(:,index,iband,ikpt)
           end do
         end do
       end do
     end do

!    The 1sin the argument of cart29 are respectively iblok and nblok. We are doing only one blok.
     call carteig2d(blkflg(:,iqpt),blkval2tmp,carflg,d2cart,gprimd,1,mpert,natom,1,rprimd)

     do idir1=1,3
       do iatom1=1,mpert
         do idir2=1,3
           do iatom2=1,mpert
             index = idir1 + 3*((iatom1 - 1) + natom * ((idir2-1)+3*(iatom2-1)))
             blkval2(:,index,iband,ikpt) = d2cart(:,idir1,iatom1,idir2,iatom2)
           end do
         end do
       end do
     end do

   end do
 end do

!Deallocating local arrays
 ABI_DEALLOCATE(blkflgtmp)
 ABI_DEALLOCATE(blkval2tmp)
 ABI_DEALLOCATE(d2cart)

end subroutine carttransf
!!***
