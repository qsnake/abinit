!{\src2tex{textfont=tt}}
!!****f* ABINIT/cart29
!! NAME
!! cart29
!!
!!
!! FUNCTION
!! Transform a second-derivative matrix from reduced
!! coordinates to cartesian coordinates, and also
!! 1) add the ionic part of the effective charges,
!! 2) normalize the electronic dielectric tensor, and
!!     add the vacuum polarisation
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blkflg(3,mpert,3,mpert,nblok)=
!!   ( 1 if the element of the dynamical matrix has been calculated ;
!!     0 otherwise )
!!  blkval(2,3,mpert,3,mpert,nblok)=DDB values
!!  gprimd(3,3)=basis vector in the reciprocal space
!!  iblok=number of the blok that will be transformed
!!  mpert =maximum number of ipert
!!  natom=number of atom
!!  nblok=number of blocks (dimension of blkflg and blkval)
!!  ntypat=number of atom types
!!  rprimd(3,3)=basis vector in the real space
!!  typat(natom)=integer label of each type of atom (1,2,...)
!!  ucvol=unit cell volume
!!  zion(ntypat)=charge corresponding to the atom type
!!
!! OUTPUT
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  d2cart(2,3,mpert,3,mpert)=
!!    dynamical matrix, effective charges, dielectric tensor,....
!!    all in cartesian coordinates
!!
!! SIDE EFFECTS
!
!!
!! NOTES
!
!!
!! PARENTS
!!      gath3,rdddb9
!!
!! CHILDREN
!!      cart39
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine cart29(blkflg,blkval,carflg,d2cart,&
& gprimd,iblok,mpert,natom,nblok,ntypat,rprimd,typat,ucvol,zion)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cart29'
 use interfaces_72_response, except_this_one => cart29
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ntypat
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok),typat(natom)
 integer,intent(out) :: carflg(3,mpert,3,mpert)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok),gprimd(3,3),rprimd(3,3)
 real(dp),intent(in) :: zion(ntypat)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)

! *********************************************************************

!First, copy the data blok in place.
 d2cart(:,:,:,:,:)=blkval(:,:,:,:,:,iblok)

!Cartesian coordinates transformation (in two steps)
!First step
 do ipert1=1,mpert
   do ipert2=1,mpert
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           vec1(idir2)=d2cart(ii,idir1,ipert1,idir2,ipert2)
!          Note here blkflg
           flg1(idir2)=blkflg(idir1,ipert1,idir2,ipert2,iblok)
         end do
         call cart39(flg1,flg2,&
&         gprimd,ipert2,natom,rprimd,vec1,vec2)
         do idir2=1,3
           d2cart(ii,idir1,ipert1,idir2,ipert2)=vec2(idir2)
!          And here carflg
           carflg(idir1,ipert1,idir2,ipert2)=flg2(idir2)
         end do
       end do
     end do
   end do
 end do

!Second step
 do ipert1=1,mpert
   do ipert2=1,mpert
     do ii=1,2
       do idir2=1,3
         do idir1=1,3
           vec1(idir1)=d2cart(ii,idir1,ipert1,idir2,ipert2)
!          Note here carflg
           flg1(idir1)=carflg(idir1,ipert1,idir2,ipert2)
         end do
         call cart39(flg1,flg2,&
&         gprimd,ipert1,natom,rprimd,vec1,vec2)
         do idir1=1,3
           d2cart(ii,idir1,ipert1,idir2,ipert2)=vec2(idir1)
!          And here carflg again
           carflg(idir1,ipert1,idir2,ipert2)=flg2(idir1)
         end do
       end do
     end do
   end do
 end do

!For the dielectric tensor, takes into account the volume
!of the unit cell, and add the unit matrix (polarization of
!the vacuum)
 do idir1=1,3
   do idir2=1,3
     do ii=1,2
       d2cart(ii,idir1,natom+2,idir2,natom+2)=&
&       -four_pi/ucvol*d2cart(ii,idir1,natom+2,idir2,natom+2)
     end do
   end do
 end do

 do idir1=1,3
   d2cart(1,idir1,natom+2,idir1,natom+2)=&
&   1.0_dp+d2cart(1,idir1,natom+2,idir1,natom+2)
 end do

!Add the ionic charges to delta z to get the effective charges
 do ipert1=1,natom
   do idir1=1,3
     d2cart(1,idir1,ipert1,idir1,natom+2)=&
&     zion(typat(ipert1))+d2cart(1,idir1,ipert1,idir1,natom+2)
   end do
 end do
 do ipert2=1,natom
   do idir2=1,3
     d2cart(1,idir2,natom+2,idir2,ipert2)=&
&     zion(typat(ipert2))+d2cart(1,idir2,natom+2,idir2,ipert2)
   end do
 end do

!For the piezoelectric tensor, takes into account the volume
!of the unit cell
 do ipert2=natom+3,natom+4
   do idir1=1,3
     do idir2=1,3
       do ii=1,2
         d2cart(ii,idir1,natom+2,idir2,ipert2)=&
&         (1.0_dp/ucvol)*d2cart(ii,idir1,natom+2,idir2,ipert2)
       end do
     end do
   end do
 end do

end subroutine cart29
!!***
