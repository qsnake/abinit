!{\src2tex{textfont=tt}}
!!****f* ABINIT/vdot
!! NAME
!! vdot
!!
!! FUNCTION
!! Computes the cross product of two vectors
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! x1(3)=first vector
!! x2(3)=second vector
!!
!! OUTPUT
!! x3(3)=cross product of x1 * x2
!!
!! PARENTS
!!      planeint,volumeint
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vdot(x1,x2,x3)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vdot'
!End of the abilint section

 implicit none

!Arguments-------------------------------------------------------------
!arrays
 real(dp),intent(in) :: x1(3),x2(3)
 real(dp),intent(out) :: x3(3)

!Local variables-------------------------------

! *************************************************************************

 x3(1)=x1(2)*x2(3)-x2(2)*x1(3)
 x3(2)=x1(3)*x2(1)-x2(3)*x1(1)
 x3(3)=x1(1)*x2(2)-x2(1)*x1(2)

 return
end subroutine vdot
!!***
