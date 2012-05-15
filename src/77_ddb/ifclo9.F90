!{\src2tex{textfont=tt}}
!!****f* ABINIT/ifclo9
!!
!! NAME
!! ifclo9
!!
!! FUNCTION
!! Convert from cartesian coordinates to local coordinates
!! the 3*3 interatomic force constant matrix
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ifccar(3,3)= matrix of interatomic force constants in cartesian
!!  coordinates
!! vect1(3)= cartesian coordinates of the first local vector
!! vect2(3)= cartesian coordinates of the second local vector
!! vect3(3)= cartesian coordinates of the third local vector
!!
!! OUTPUT
!! ifcloc(3,3)= matrix of interatomic force constants in local
!!  coordinates
!!
!! PARENTS
!!      rsiaf9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ifclo9(ifccar,ifcloc,vect1,vect2,vect3)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifclo9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!arrays
 real(dp),intent(in) :: ifccar(3,3),vect1(3),vect2(3),vect3(3)
 real(dp),intent(out) :: ifcloc(3,3)

!Local variables -------------------------
!scalars
 integer :: ii,jj
!arrays
 real(dp) :: work(3,3)

! *********************************************************************

 do jj=1,3
   do ii=1,3
     work(jj,ii)=0.0_dp
   end do
   do ii=1,3
     work(jj,1)=work(jj,1)+ifccar(jj,ii)*vect1(ii)
     work(jj,2)=work(jj,2)+ifccar(jj,ii)*vect2(ii)
     work(jj,3)=work(jj,3)+ifccar(jj,ii)*vect3(ii)
   end do
 end do

 do jj=1,3
   do ii=1,3
     ifcloc(ii,jj)=0.0_dp
   end do
   do ii=1,3
     ifcloc(1,jj)=ifcloc(1,jj)+vect1(ii)*work(ii,jj)
     ifcloc(2,jj)=ifcloc(2,jj)+vect2(ii)*work(ii,jj)
     ifcloc(3,jj)=ifcloc(3,jj)+vect3(ii)*work(ii,jj)
   end do
 end do

end subroutine ifclo9
!!***
