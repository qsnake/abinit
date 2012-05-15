!!****f* ABINIT/scphon_supercell_vectors_init
!! NAME
!! scphon_supercell_vectors_init
!!
!! FUNCTION
!! Calculate the integer vectors, for each atom in the supercell, which point to
!! the primitive unit cell the atom is contained in. Also output an array which
!! gives the equivalent atom in the primitive cell basis, if there are several
!! atoms in the primitive cell.
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
!! natom= number of atoms in the full supercell
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! pcell=container type with ancillary variables and dimensions from anaddb run
!! supercell_multiplicity=number of times the primitive unit cell is repeated
!!   along each axis, in the supercell
!! xred= reduced coordinates of all atoms in the supercell
!!
!! OUTPUT
!! pcell_atom_in_supercell= mapping of atoms to an atom index in the primitive
!!   unit cell
!! supercell_vectors= vector for each atom in the supercell, which points
!!   to the unit cell it is contained in, in integer units of the primitive cell
!!   lattice vectors
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_supercell_vectors_init(natom,natom_primitive_cell,&
&   pcell,pcell_atom_in_supercell,&
&   supercell_multiplicity,supercell_vectors,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_primcell_ddb_info

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_supercell_vectors_init'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,natom_primitive_cell
 type(primcell_ddb_info),intent(inout) :: pcell
!arrays
 integer,intent(in) :: supercell_multiplicity(3)
 integer,intent(out) :: pcell_atom_in_supercell(natom)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: supercell_vectors(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,iatom_primcell,idir,ii
 character(len=500) :: message
!arrays
 real(dp) :: relative_position(3)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 write (message,'(a)') ' supercell vectors :  '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 write (message,*) ' supercell multiplicity:',supercell_multiplicity(:)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 write(ab_out,*) 'xred:'
 do ii=1,natom
   write(ab_out,*) xred(:,ii)
 end do

 do iatom=1,natom
!  for each atom find unit cell which it should belong to.
   do idir=1,3
     supercell_vectors(idir,iatom)=dble(floor(&
&     xred(idir,iatom)*supercell_multiplicity(idir) + 0.01))
     if (supercell_vectors(idir,iatom) >= supercell_multiplicity(idir) .or. &
&     supercell_vectors(idir,iatom) < 0) then
       write (message,'(a,2I6,2x,I6)') 'error in atomic disposition ',&
&       iatom,idir,supercell_vectors(idir,iatom)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       stop
     end if
   end do
!  find which atom in the primitive cell corresponds to the present one in the
!  supercell
   pcell_atom_in_supercell(iatom)=-1
!  1. if the ordering is preserved, just modulo natom_primitive_cell
   pcell_atom_in_supercell(iatom) = mod(iatom-1,natom_primitive_cell) + 1

!  2. otherwise really need to seek primcell atom
   relative_position=xred(:,iatom)-dble(supercell_vectors(:,iatom))/dble(supercell_multiplicity(:))
   do iatom_primcell=1,natom_primitive_cell
     if (sum((relative_position(:)-pcell%xred(:,iatom_primcell))**2) < tol8) then
       pcell_atom_in_supercell(iatom) = iatom_primcell
       exit
     end if
   end do
   if (pcell_atom_in_supercell(iatom) < 1  .or. &
&   pcell_atom_in_supercell(iatom) > natom_primitive_cell) then
     write(std_out,*) 'Error: pcell atom index is out of bounds '
     stop
   end if

   write (message,'(3E20.8)') supercell_vectors (:,iatom)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,I6)') ' equiv primcell atom ', pcell_atom_in_supercell (iatom)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end do

 write(message,'(a)') 'supercell_vectors:'
 call wrtout(ab_out,message,'COLL')
 do ii=1,natom
   write(ab_out,'(3E20.10)') supercell_vectors(:,ii)
   call wrtout(ab_out,message,'COLL')
 end do

end subroutine scphon_supercell_vectors_init
!!***


