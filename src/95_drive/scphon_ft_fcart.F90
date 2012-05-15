!!****f* ABINIT/scphon_ft_fcart
!! NAME
!! scphon_ft_fcart
!!
!! FUNCTION
!! Fourier Transform cartesian forces on all supercell atoms, with respect to
!! the supercell lattice vectors (ie multiples of the primitive unit cell which
!! are contained in the supercell). This returns a force on each atom in the
!! primitive unit cell, for each q-vector in the dual grid of the supercell.
!!
!! The force is divided by the square root of the mass of the appropriate atom, in prevision of the
!! calculation of new frequencies.
!!
!! The dual grid should be the same as the input q-point grid on which the
!! equilibrium phonons were calculated.
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
!! fcart= forces on all supercell atoms, in cartesian coordinates
!! natom= number of atoms in the full supercell
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! phononq= phonon q vectors used in anaddb run (reduced coordinates)
!! pcell_atom_in_supercell= mapping of atoms to an atom index in the primitive
!!   unit cell
!! supercell_vectors= vector for each atom in the supercell, which points
!!   to the unit cell it is contained in, in integer units of the primitive cell
!!
!! OUTPUT
!! forces_on_atoms_ft= FT of cartesian forces on atoms, wrt the superlattice
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_ft_fcart(sqrt_amass_pcell,fcart,natom,natom_primitive_cell,nphononq,phononq,&
&   forces_on_atoms_ft,pcell_atom_in_supercell,supercell_vectors)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_ft_fcart'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalar
!scalars
 integer,intent(in) :: natom,natom_primitive_cell,nphononq
!arrays
 integer,intent(in) :: pcell_atom_in_supercell(natom)
 real(dp),intent(in) :: fcart(3,natom),phononq(3,nphononq)
 real(dp),intent(in) :: sqrt_amass_pcell(natom_primitive_cell)
 real(dp),intent(in) :: supercell_vectors(3,natom)
 real(dp),intent(out) :: forces_on_atoms_ft(2,3*natom_primitive_cell,nphononq)

!Local variables-------------------------------
!scalars
 integer :: iatom,iatom_in_pcell,idir,indx_pcell,iq
 real(dp) :: argument
 character(len=500) :: message

! *************************************************************************

 if (natom_primitive_cell*nphononq /= natom) then
   write (message,'(a,a)') 'Error: number of phonon q times number of atoms ',&
&   ' in unit cell should be equal to number of atoms in supercell'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

 forces_on_atoms_ft = zero
 do iq=1,nphononq
   do iatom=1,natom
     iatom_in_pcell=pcell_atom_in_supercell(iatom)

     argument=sum(supercell_vectors(1:3,iatom)*phononq(1:3,iq))
!    abinit does not add 2 pi factor into k vectors
     argument=argument*two_pi

     do idir=1,3
       indx_pcell=idir+(iatom_in_pcell-1)*3
!      Presumes fcart is real (which is certain), so its FT has some inversion
!      symmetry
       forces_on_atoms_ft(1,indx_pcell,iq)=forces_on_atoms_ft(1,indx_pcell,iq)+&
&       fcart(idir,iatom)*cos(argument)/sqrt_amass_pcell(iatom_in_pcell)
       forces_on_atoms_ft(2,indx_pcell,iq)=forces_on_atoms_ft(2,indx_pcell,iq)+&
&       fcart(idir,iatom)*sin(argument)/sqrt_amass_pcell(iatom_in_pcell)
     end do
   end do
 end do

!looks like this is the choice of normalization in the PRL
 forces_on_atoms_ft = forces_on_atoms_ft/sqrt(dble(nphononq))

 write(std_out,*) ' Re(FT of fcart) / sqrt(M) : '
 write(std_out,'(3(E20.10,2x))') forces_on_atoms_ft(1,:,:)
 write(std_out,*)

end subroutine scphon_ft_fcart
!!***

