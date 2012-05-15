!{\src2tex{textfont=tt}}
!!****f* ABINIT/recip
!! NAME
!! recip
!!
!! FUNCTION
!! Computes the reciprocal unit cell
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (XG,RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! hkl(3)=Miller indices of the plane
!! rprimd(3,3)=orientation of the unit cell in 3D
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! x1(3)=point coordinates
!!
!! PARENTS
!!      planeint,volumeint
!!
!! CHILDREN
!!      matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine recip(x1,hkl,rprimd)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recip'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments--------------------------------------------------------------
!arrays
 integer,intent(in) :: hkl(3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: x1(3)

!Local variables--------------------------------------------------------
!scalars
 integer :: ii
!arrays
 real(dp) :: mminv(3,3)

! *************************************************************************

 call matr3inv(rprimd,mminv)

!write(std_out,*) 'rec.lattice (cartesian):',mminv(1,1),mminv(1,2),mminv(1,3)
!write(std_out,*) 'rec.lattice (cartesian):',mminv(2,1),mminv(2,2),mminv(2,3)
!write(std_out,*) 'rec.lattice (cartesian):',mminv(3,1),mminv(3,2),mminv(3,3)

 do ii=1,3
   x1(ii)=mminv(ii,1)*hkl(1) + mminv(ii,2)*hkl(2) + mminv(ii,3)*hkl(3)
 end do

end subroutine recip

!!***
