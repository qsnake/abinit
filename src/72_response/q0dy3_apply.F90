!{\src2tex{textfont=tt}}
!!****f* ABINIT/q0dy3_apply
!! NAME
!! q0dy3_apply
!!
!! FUNCTION
!! Takes care of the inclusion of the ewald q=0 term in the dynamical
!! matrix - corrects the dyew matrix provided as input
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG, MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dyewq0(3,3,natom) = part needed to correct
!!    the dynamical matrix for atom self-interaction.
!!  natom= number of atom in the unit cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dyew(2,3,natom,3,natom)= dynamical matrix corrected on output
!!
!! NOTES
!! Should be used just after each call to ewald3, for both
!! q==0 and the real wavelength.
!!
!! The q0dy3_apply should be used in conjunction with the subroutine
!! ewald3 (or ewald9):
!! First, the call of ewald3 with q==0 should be done,
!!   then the call to q0dy3_calc will produce
!!   the dyewq0 matrix from the (q=0) dyew matrix
!! Second, the call of ewald3 with the real q (either =0 or diff 0)
!!   should be done, then the call to q0dy3_apply
!!   will produce the correct dynamical matrix dyew starting from
!!   the previously calculated dyewq0 and the bare(non-corrected)
!!   dyew matrix
!!
!! PARENTS
!!      gtdyn9,mkifc9,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine q0dy3_apply(natom,dyewq0,dyew)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'q0dy3_apply'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: dyewq0(3,3,natom)
 real(dp),intent(inout) :: dyew(2,3,natom,3,natom)

!Local variables -------------------------
!scalars
 integer :: ia,mu,nu

! *********************************************************************

 do mu=1,3
   do nu=1,3
     do ia=1,natom
       dyew(1,mu,ia,nu,ia)=dyew(1,mu,ia,nu,ia)-dyewq0(mu,nu,ia)
     end do
   end do
 end do

end subroutine q0dy3_apply
!!***
