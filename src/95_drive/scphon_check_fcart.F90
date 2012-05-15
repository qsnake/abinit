!!****f* ABINIT/scphon_check_fcart
!! NAME
!! scphon_check_fcart
!!
!! FUNCTION
!! Check that the cartesian displacements presently being imposed on the atoms
!! in the supercell are consistent with the calculated forces and the hypothesis
!! that the phonons are harmonic, or at least that the equilibrium structure is
!! a minimum of energy: takes the scalar product of the displacement by the
!! force, for each atom, and this should be negative
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cartesian_displacements= displacement of all atoms from their initial
!!   equilibrium positions
!! fcart= forces on all atoms, in cartesian coordinates
!! natom= number of atoms in the full supercell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_check_fcart(cartesian_displacements,fcart,natom)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_check_fcart'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: cartesian_displacements(3,natom),fcart(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom
 real(dp) :: scprod

! *************************************************************************

 write(std_out,*) 'natom=',natom

 write(std_out,'(a)') 'cartesian_displacements(3,natom) = '
 do iatom=1,natom
   write(std_out,'(3E20.10)')  cartesian_displacements(:,iatom)
 end do

 write(std_out,'(a)') 'fcart(3,natom) = '
 do iatom=1,natom
   write(std_out,'(3E20.10)')  fcart(:,iatom)
 end do

 write(std_out,*) ' Calculate scalar product of force times displacement.'
 write(std_out,*) '  Should be negative'
 do iatom=1,natom
   scprod=sum(fcart(:,iatom)*cartesian_displacements(:,iatom))
   write(std_out,*) 'atom ', iatom, ' F.dR = ', scprod
 end do


end subroutine scphon_check_fcart
!!***

