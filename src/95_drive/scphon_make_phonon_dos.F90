!!****f* ABINIT/scphon_make_phonon_dos
!! NAME
!! scphon_make_phonon_dos
!!
!! FUNCTION
!! Simple Gaussian summation of the density of states for the phonon frequencies
!! given in input. This is robust, but superceded by
!! scphon_interpolate_phonon_and_dos which does much better by interpolating the
!! phonon frequencies explicitly before summing the DOS.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dos_smearing= smearing width used to smooth the phonon DOS
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! nfreq_int= number of frequencies in the grid on which the phonon DOS is
!!   calculated
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! pcell=container type with ancillary variables and dimensions from anaddb run
!! phonon_eigval=phonon eigenfrequencies, updated inside SC phonon run
!!
!! OUTPUT
!! maxfreq=maximum frequency for which the DOS is calculated
!! minfreq=minimum frequency for which the DOS is calculated
!! phonon_dos=array containing the phonon DOS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!      simpson_int
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_make_phonon_dos (dos_smearing,natom_primitive_cell,&
&    nfreq_int,nphononq,maxfreq,minfreq,phonon_dos,phonon_eigval)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_make_phonon_dos'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom_primitive_cell,nfreq_int,nphononq
 real(dp),intent(in) :: dos_smearing
 real(dp),intent(out) :: maxfreq,minfreq
!arrays
 real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
 real(dp),intent(out) :: phonon_dos(nfreq_int)

!Local variables-------------------------------
!scalars
 integer :: ifreq,imode,iq,unit_phondos
 real(dp) :: domega,freq,gauss_prefactor,inv_dos_smearing
!arrays
 real(dp) :: phonon_dos_int(nfreq_int)

! *************************************************************************

 inv_dos_smearing = one / dos_smearing
 gauss_prefactor = one/dos_smearing/sqrt(pi)

 maxfreq=1.1_dp*maxval(phonon_eigval)
 minfreq=minval(phonon_eigval)

!add a few negative frequencies to account for smearing of acoustic modes
 minfreq=minfreq-dos_smearing*five

 domega=(maxfreq-minfreq)/dble(nfreq_int-1)

!calculate phonon DOS from (small) number of frequencies we have
 phonon_dos= zero
 do iq=1,nphononq
   do imode=1,3*natom_primitive_cell
     do ifreq=1,nfreq_int
       freq=(minfreq+dble(ifreq-1)*(maxfreq-minfreq)/dble(nfreq_int-1) &
&       - phonon_eigval(imode,iq))*inv_dos_smearing
       if (abs(freq) > seven) cycle
       phonon_dos(ifreq) = phonon_dos(ifreq) + gauss_prefactor*exp(-freq*freq)
     end do
   end do
 end do

!normalize for number of qpoints
 phonon_dos=phonon_dos/dble(nphononq)

 unit_phondos=401
 open (unit=unit_phondos,file='phonondos.dat')
 write (unit_phondos,*) '#'
 write (unit_phondos,*) '#  phonon dos calculated from self-consistent phonon spectrum'
 write (unit_phondos,*) '#'
 do ifreq=1,nfreq_int
   freq=minfreq+dble(ifreq-1)*(maxfreq-minfreq)/dble(nfreq_int-1)
   write (unit_phondos,*) freq, phonon_dos(ifreq)
 end do

 phonon_dos_int=zero
 call simpson_int(nfreq_int,domega,phonon_dos,phonon_dos_int)

 write (unit_phondos,'(a,F10.3)') '# integral = ', phonon_dos_int(nfreq_int)
 close(unit_phondos)

end subroutine scphon_make_phonon_dos
!!***


