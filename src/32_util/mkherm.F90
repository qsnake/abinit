!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkherm
!! NAME
!! mkherm
!!
!! FUNCTION
!! Make the complex array(ndim,ndim) hermitian,
!!     by adding half of it to its hermitian conjugate.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! ndim=dimension of the matrix
!! array= complex matrix
!!
!! SIDE EFFECTS
!! array= hermitian matrix made by adding half of array to its hermitian conjugate
!!
!! PARENTS
!!      anaddb,phfrq3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkherm(array,ndim)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkherm'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ndim
!arrays
 real(dp),intent(inout) :: array(2,ndim,ndim)

!Local variables -------------------------
!scalars
 integer :: i1,i2

! *********************************************************************

 do i1=1,ndim
   do i2=1,i1
     array(1,i1,i2)=(array(1,i1,i2)+array(1,i2,i1))*half
     array(2,i1,i2)=(array(2,i1,i2)-array(2,i2,i1))*half
     array(1,i2,i1)=array(1,i1,i2)
     array(2,i2,i1)=-array(2,i1,i2)
   end do
 end do

end subroutine mkherm
!!***
