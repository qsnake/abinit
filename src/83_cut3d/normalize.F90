!{\src2tex{textfont=tt}}
!!****f* ABINIT/normalize
!! NAME
!! normalize
!!
!! FUNCTION
!! Normalizes the value of v
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  v = on entry, vector to be normalized
!!
!! OUTPUT
!!  v = on exit, vector normalized

!! SIDE EFFECTS
!!   v=value to be normalized
!!
!! PARENTS
!!      lineint,planeint,volumeint
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine normalize(v)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'normalize'
!End of the abilint section

 implicit none

!Arguments-------------------------------------------------------------
!arrays
 real(dp),intent(inout) :: v(3)

!Local variables--------------------------------------------------------
!scalars
 integer :: idir
 real(dp) :: norm

! *************************************************************************

 norm=0.0
 do idir=1,3
   norm=norm+v(idir)**2
 end do
 norm=sqrt(norm)

 do idir=1,3
   v(idir)=v(idir)/norm
 end do

 return
end subroutine normalize
!!***
