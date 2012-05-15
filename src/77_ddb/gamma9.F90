!{\src2tex{textfont=tt}}
!!****f* ABINIT/gamma9
!!
!! NAME
!! gamma9
!!
!! FUNCTION
!! This small routine checks if the wavevector qphon and the
!! corresponding normalisation factor represent a phonon at Gamma.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! qphon(3)=wavevector
!! qphnrm=normalisation factor
!! qtol=tolerance
!!
!! OUTPUT
!! gamma= if 1, means that the wavevector is indeed at Gamma
!!  otherwise 0.
!!
!! PARENTS
!!      gtblk9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gamma9(gamma,qphon,qphnrm,qtol)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gamma9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(out) :: gamma
 real(dp),intent(in) :: qphnrm,qtol
!arrays
 real(dp),intent(in) :: qphon(3)

!Local variables-------------------------------

! *********************************************************************

 if( (  abs(qphon(1))<qtol .and.  &
& abs(qphon(2))<qtol .and.        &
& abs(qphon(3))<qtol      ) .or.  &
& abs(qphnrm)<qtol )then
   gamma=1
 else
   gamma=0
 end if

end subroutine gamma9
!!***
