!!****f* ABINIT/scphon_phonon_init
!! NAME
!! scphon_phonon_init
!!
!! FUNCTION
!! Return phonon frequencies and eigenvectors initialized in
!! phonon_eigvec_ref and phonon_eigval_ref, read in from files
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtfil= datatype for file names and units
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!!
!! OUTPUT
!! phonon_eigval_ref=phonon eigenfrequencies, from the anaddb equil run
!! phonon_eigvec_ref=reference phonon eigenvectors, from the anaddb equil run
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


subroutine scphon_phonon_init (fnameabi_phfrq,fnameabi_phvec,&
&  natom_primitive_cell,nphononq,phonon_eigvec_ref,&
&  phonon_eigval_ref)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_phonon_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom_primitive_cell,nphononq
 character(len=fnlen),intent(in) :: fnameabi_phfrq
 character(len=fnlen),intent(in) :: fnameabi_phvec

!arrays
 real(dp),intent(out) :: phonon_eigval_ref(3*natom_primitive_cell,nphononq)
 real(dp),intent(out) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)

!Local variables-------------------------------
!scalars
 integer :: imode_primitive_cell,indx,iq,jmode_primitive_cell
 integer :: phfrq_unit,phvec_unit
!character(len=500) :: message

! ************************************************************************

!read in _PHFRQ and _PHVEC files from anaddb output
!

!first eigenvalues
 phfrq_unit=300
 open (unit=phfrq_unit,file=fnameabi_phfrq)

!read in header of 3 lines
 read(phfrq_unit,*)
 read(phfrq_unit,*)
 read(phfrq_unit,*)

!for each qpoint, read in index frq_1 frq_2 ... frq_3natom
 do iq=1,nphononq
   read (phfrq_unit,'(i6)',ADVANCE='NO') indx
   do imode_primitive_cell=1,3*natom_primitive_cell
     read (phfrq_unit,'(e16.8,2x)',ADVANCE='NO') phonon_eigval_ref(imode_primitive_cell,iq)
   end do
   read (phfrq_unit,*)
 end do

 close (phfrq_unit)

!now eigenvectors
 phvec_unit=300
 open (unit=phvec_unit,file=fnameabi_phvec)

!read in header of 3 lines
 read(phvec_unit,*)
 read(phvec_unit,*)
 read(phvec_unit,*)

!for each qpoint, read in index frq_1 frq_2 ... frq_3natom
 do iq=1,nphononq
   read (phvec_unit,'(i6)',ADVANCE='NO') indx
   do imode_primitive_cell=1,3*natom_primitive_cell
     do jmode_primitive_cell=1,3*natom_primitive_cell
       read (phvec_unit,'(2e25.15,2x)',ADVANCE='NO') phonon_eigvec_ref(:,jmode_primitive_cell,imode_primitive_cell,iq)
     end do
   end do
   read (phvec_unit,*)
 end do

 close (phvec_unit)

end subroutine scphon_phonon_init
!!***

