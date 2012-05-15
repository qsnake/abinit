!!****f* ABINIT/scphon_update_xcart
!! NAME
!! scphon_update_xcart
!!
!! FUNCTION
!! From normal mode displacements, calculate the cartesian displacements
!! for all atoms in the supercell, and update xcart
!!
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! sqrt_amass_pcell= masses of the atoms in the primitive unit cell, in atomic units
!! natom= number of atoms in the full supercell
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! normal_mode_displacements= calculated displacements of canonical coordinates
!!   (normal modes) of phonons at desired temperature.
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! pcell_atom_in_supercell= mapping of atoms to an atom index in the primitive
!!   unit cell
!! phonon_eigvec_ref=reference phonon eigenvectors, from the anaddb equil run
!! phononq= phonon q vectors used in anaddb run (reduced coordinates)
!! supercell_vectors= vector for each atom in the supercell, which points
!!   to the unit cell it is contained in, in integer units of the primitive cell
!!   lattice vectors
!! xcart0= initial positions of atoms, before perturbations
!!
!! OUTPUT
!! cartesian_displacements= displacements of all atoms in the supercell, in
!!   cartesian coordinates
!! xcart= old, then new positions of all atoms, in cartesian coordinates
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!
!! SOURCE
! update xcart with difference between old and new normal mode displacements

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_update_xcart (sqrt_amass_pcell,cartesian_displacements,natom,&
&   natom_primitive_cell,normal_mode_displacements,&
&   nphononq,pcell_atom_in_supercell,phonon_eigvec_ref,phononq,supercell_vectors,xcart,xcart0)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_update_xcart'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,natom_primitive_cell,nphononq
!arrays
 integer,intent(in) :: pcell_atom_in_supercell(natom)
 real(dp),intent(in) :: normal_mode_displacements(3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: phononq(3,nphononq)
 real(dp),intent(in) :: sqrt_amass_pcell(natom_primitive_cell)
 real(dp),intent(in) :: supercell_vectors(3,natom),xcart0(3,natom)
 real(dp),intent(inout) :: cartesian_displacements(3,natom),xcart(3,natom)

!Local variables-------------------------------
! the instantaneous change in the cartesian displacements
!scalars
 integer :: iatom,iatom_in_pcell,imode_primitive_cell,iq
 real(dp) :: argument,cosarg,sinarg
!arrays
 real(dp),allocatable :: delta_cartesian_displacements(:,:)

! *************************************************************************

 ABI_ALLOCATE(delta_cartesian_displacements,(3,natom))

!
!NOTE: only real part of cartesian displacement is taken here.
!what happens if there is an imaginary part?
!
 delta_cartesian_displacements=zero
 do iatom=1,natom
   iatom_in_pcell=pcell_atom_in_supercell(iatom)
   do iq=1,nphononq
     argument=-two_pi * sum(phononq(:,iq)*supercell_vectors(:,iatom))
     cosarg=cos(argument)
     sinarg=sin(argument)

     do imode_primitive_cell=1,3*natom_primitive_cell
       delta_cartesian_displacements(:,iatom)=delta_cartesian_displacements(:,iatom) + &
&       normal_mode_displacements(imode_primitive_cell,iq)*&
&       (phonon_eigvec_ref(1,(iatom_in_pcell-1)*3+1:iatom_in_pcell*3,imode_primitive_cell,iq)*cosarg  &
&       -phonon_eigvec_ref(2,(iatom_in_pcell-1)*3+1:iatom_in_pcell*3,imode_primitive_cell,iq)*sinarg) &
&       / sqrt_amass_pcell(iatom_in_pcell)
     end do
   end do
 end do
!Normalization chosen in PRL and confirmed by Souvatzis
 delta_cartesian_displacements = delta_cartesian_displacements &
& /sqrt(dble(nphononq))
!TODO: this normalization should be by number of ATOMS
!delta_cartesian_displacements = delta_cartesian_displacements &
!& /dble(nphononq)

!!this is needed for normalization of cartesian coordinates
!! NO - the phonon displacement vectors are already normalized correctly
!delta_cartesian_displacements = delta_cartesian_displacements &
!& /sqrt(dble(natom_primitive_cell))

 xcart = xcart0 + delta_cartesian_displacements
 cartesian_displacements=delta_cartesian_displacements

 write(std_out,*)
 write(std_out,*) '  new xcart = '
 write(std_out,'(3E20.10)')  xcart
 write(std_out,*)

end subroutine scphon_update_xcart
!!***



