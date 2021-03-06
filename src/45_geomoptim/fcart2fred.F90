!{\src2tex{textfont=tt}}
!!****f* ABINIT/fcart2fred
!!
!! NAME
!! fcart2fred
!!
!! FUNCTION
!! Convert cartesian forces into reduced forces
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!  natom=Number of atoms in the unitary cell 
!!  rprimd(3,3)=dimensional primitive
!!
!! OUTPUT
!!  fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!    Unlike fred, fcart has been corrected by enforcing
!!    the translational symmetry, namely that the sum of force
!!    on all atoms is zero.
!!
!! PARENTS
!!      mover,prec_simple
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fcart2fred(fcart,fred,rprimd,natom)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fcart2fred'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: fcart(3,natom)
 real(dp),intent(out) :: fred(3,natom)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu
!arrays

! *************************************************************************

 do iatom=1,natom
   do mu=1,3
     fred(mu,iatom)= - (rprimd(mu,1)*fcart(1,iatom)+&
&     rprimd(mu,2)*fcart(2,iatom)+&
&     rprimd(mu,3)*fcart(3,iatom))
   end do
 end do

end subroutine fcart2fred
!!***
