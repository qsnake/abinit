!{\src2tex{textfont=tt}}
!!****f* ABINIT/overlap_g
!! NAME
!! overlap_g
!!
!! FUNCTION
!! Compute the scalar product between WF at two different k-points
!! < u_{n,k1} | u_{n,k2} >
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! mpw = maximum dimensioned size of npw
!! npw_k1 = number of plane waves at k1
!! pwind_k = array required to compute the scalar product (see initberry.f)
!! vect1 = wavefunction at k1: | u_{n,k1} >
!! vect2 = wavefunction at k1: | u_{n,k2} >
!!
!! OUTPUT
!! doti = imaginary part of the scalarproduct
!! dotr = real part of the scalarproduct
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!! In case a G-vector of the basis sphere of plane waves at k1
!! does not belong to the basis sphere of plane waves at k2,
!! pwind = 0. Therefore, the dimensions of vect1 &
!! vect2 are (1:2,0:mpw) and the element (1:2,0) MUST be
!! set to zero.
!!
!! The current implementation if not compatible with TR-symmetry (i.e. istwfk/=1) !
!!
!! PARENTS
!!      bec3,die3,ebp3,edie3,gbefd3,gradberry3,qmatrix,smatrix,smatrix_paw
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine overlap_g(doti,dotr,mpw,npw_k1,pwind_k,vect1,vect2)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'overlap_g'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpw,npw_k1
 real(dp),intent(out) :: doti,dotr
!arrays
 integer,intent(in) :: pwind_k(mpw)
 real(dp),intent(in) :: vect1(1:2,0:mpw),vect2(1:2,0:mpw)

!Local variables-------------------------------
!scalars
 integer :: ipw,jpw
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)'overlap_g : enter'
!ENDDEBUG

 dotr = zero ; doti = zero

!Check if vect1(:,0) = 0 and vect2(:,0) = 0

 if ((abs(vect1(1,0)) > tol12).or.(abs(vect1(2,0)) > tol12).or. &
& (abs(vect2(1,0)) > tol12).or.(abs(vect2(2,0)) > tol12)) then
   write(message,'(a,a,a)')' overlap_g : BUG -',ch10,&
&   '   vect1(:,0) and/or vect2(:,0) are not equal to zero'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Compute the scalar product

!$OMP PARALLEL DO ORDERED PRIVATE(ipw,jpw) REDUCTION(+:doti,dotr) &
!$OMP SHARED(pwind_k,vect1,vect2,npw_k1)
 do ipw = 1, npw_k1
   jpw = pwind_k(ipw)
   dotr = dotr + vect1(1,ipw)*vect2(1,jpw) + vect1(2,ipw)*vect2(2,jpw)
   doti = doti + vect1(1,ipw)*vect2(2,jpw) - vect1(2,ipw)*vect2(1,jpw)
 end do
!$OMP END PARALLEL DO

end subroutine overlap_g
!!***
