!{\src2tex{textfont=tt}}
!!****f* ABINIT/reduce
!! NAME
!! reduce
!!
!! FUNCTION
!! Transforms coordinates of an input point
!! from cartesian to crystallographic
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (RC,XG,LSI)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! rcart(3)=position vector in crystallographic coordinates
!! rprimd(3,3)=orientation of the unit cell in 3D
!!
!! OUTPUT
!! r(3)=position vector in cartesian coordinates
!!
!! PARENTS
!!      lineint,planeint,pointint,volumeint
!!
!! CHILDREN
!!      matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine reduce(r,rcart,rprimd)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'reduce'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments-------------------------------------------------------------
!arrays
 real(dp),intent(in) :: rcart(3),rprimd(3,3)
 real(dp),intent(out) :: r(3)

!Local variables--------------------------------------------------------
!scalars
!arrays
 real(dp) :: mminv(3,3)

! *************************************************************************

 call matr3inv(rprimd,mminv)
 r(1)=rcart(1)*mminv(1,1)+rcart(2)*mminv(2,1)+rcart(3)*mminv(3,1)
 r(2)=rcart(1)*mminv(1,2)+rcart(2)*mminv(2,2)+rcart(3)*mminv(3,2)
 r(3)=rcart(1)*mminv(1,3)+rcart(2)*mminv(2,3)+rcart(3)*mminv(3,3)

end subroutine reduce
!!***
