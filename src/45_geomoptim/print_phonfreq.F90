!!****f* ABINIT/print_phonfreq
!! NAME
!! print_phonfreq
!!
!! FUNCTION
!!   print phonon frequencies to standardized file SCphon_TIMx_PHFRQ
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! istep=iteration number in SC phonon run -1 gives TIM0, as for Broyden algorithm
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! phonon_eigval=phonon eigenfrequencies, updated inside SC phonon run
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  writes to file
!!
!! PARENTS
!!      scphon,scphon_new_frequencies
!!
!! CHILDREN
!!      fappnd
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine print_phonfreq(istep,natom_primitive_cell,nphononq,phonon_eigval)

 use m_profiling

 use defs_basis
 use m_io_tools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_phonfreq'
 use interfaces_45_geomoptim, except_this_one => print_phonfreq
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,natom_primitive_cell,nphononq
!arrays
 real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)

!Local variables-------------------------------
! scfrq_unit= file unit for output of SC phonon frequencies as a function of iteration
! FIXME: replace by call to getunit
! should become input from files file out radix
!scalars
 integer :: imode_primitive_cell,iq,scfrq_unit
 character(len=fnlen) :: outfilename_radix
 character(len=fnlen) :: phonon_freq_filename

! *************************************************************************

 outfilename_radix="SCphon"
 call fappnd(phonon_freq_filename,outfilename_radix,istep)
 phonon_freq_filename = trim(phonon_freq_filename)//"_PHFRQ"
 scfrq_unit = get_unit()
 open (unit=scfrq_unit,file=phonon_freq_filename)
 write (scfrq_unit,*) '#'
 write (scfrq_unit,*) '# phonon frequencies (in Ha) on qph1l list of qpoints'
 write (scfrq_unit,*) '#'

 do iq=1,nphononq
   write (scfrq_unit,'(i6)',ADVANCE='NO') iq
   do imode_primitive_cell=1,3*natom_primitive_cell
     write (scfrq_unit,'(E16.8,2x)',ADVANCE='NO') phonon_eigval(imode_primitive_cell,iq)
   end do
   write (scfrq_unit,*)
 end do

 close (scfrq_unit)

end subroutine print_phonfreq
!!***
