!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_prtscphon
!! NAME
!!  m_prtscphon
!!
!! FUNCTION
!!  This module contains routines to write out self-consistent phonon auxiliary files
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt
!!
!! PARENTS 
!! 
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_prtscphon

 use m_profiling

 implicit none 

 private

 integer :: phfrq_unit, phvec_unit, iphl1, outscphon

 public :: beginprtscphon
 public :: prtscphon
 public :: endprtscphon
!!***

contains

!****f* ABINIT/beginprtscphon
!!
!! NAME
!! beginprtscphon
!!
!! FUNCTION
!! Open the file units for scphon outputs and print header
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  phonon_freq_filename = file name for frequencies
!!  phonon_vec_filename = file name for eigenvectors
!!  outscphon = flag to output the vectors (or not)
!!
!! OUTPUTS
!!
!! NOTES
!!  opens 2 file units contained in present module
!! 
!! PARENTS
!!      mkphbs,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine beginprtscphon(phonon_freq_filename, phonon_vec_filename, outscphon_in)

 use defs_basis
 use m_io_tools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'beginprtscphon'
!End of the abilint section

 implicit none

 integer, intent(in) :: outscphon_in
 character(len=fnlen), intent(in) :: phonon_freq_filename, phonon_vec_filename

!  open and write header for phonon frequency file
 phfrq_unit=get_unit()
 open (unit=phfrq_unit,file=phonon_freq_filename)
 write (phfrq_unit,*) '#'
 write (phfrq_unit,*) '# phonon frequencies (in Ha) on fineqpath list of qpoints'
 write (phfrq_unit,*) '#'

 if (outscphon > 0) then
!  open and write header for phonon frequency file
   phvec_unit=get_unit()
   open (unit=phvec_unit,file=phonon_vec_filename)
   write (phvec_unit,*) '#'
   write (phvec_unit,*) '# phonon eigenvectors (dimensionless) on fineqpath list of qpoints'
   write (phvec_unit,*) '#'
 end if

 iphl1=1

 outscphon = outscphon_in

end subroutine beginprtscphon
!!***

!****f* ABINIT/prtscphon
!!
!! NAME
!! prtscphon
!!
!! FUNCTION
!! Print out one line of phonon frequencies resp. eigenvectors to appropriate files, for 1 q-point
!! Used for output in context of self consistent phonons
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! eigvec(2,3,natom,3,natom) = phonon eigenvectors for 1 q-point
!! natom = number of atoms
!! phfrq = phonon frequencies for 1 q-point
!!
!! OUTPUTS
!!
!! NOTES
!!  writes to files open on phfrq_unit and phvec_unit
!! 
!! PARENTS
!!      mkphbs,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine prtscphon(eigvec, natom, phfrq, qphon)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtscphon'
!End of the abilint section

 implicit none

!scalars
 integer, intent(in) :: natom
!arrays
 real(dp), intent(in) :: phfrq(3*natom), eigvec(2,3,natom,3,natom), qphon(3)

!local
 integer :: ii, idir, jj, jdir
 
 write (phfrq_unit,'(i6)',ADVANCE='NO') iphl1
 do ii=1,3*natom
   write (phfrq_unit,'(E16.8,2x)',ADVANCE='NO') phfrq(ii)
 end do
 write (phfrq_unit,'(a,3E20.10)') '# qpt ', qphon

 if (outscphon > 0) then
! write phonon eigenvectors to file
   write (phvec_unit,'(i6)',ADVANCE='NO') iphl1
     do ii=1,natom
     do idir=1,3
       do jj=1,natom
         do jdir=1,3
           write (phvec_unit,'(2E25.15,2x)',ADVANCE='NO') eigvec(:,jdir,jj,idir,ii)
         end do
       end do
     end do
   end do
   write (phvec_unit,*)
 end if ! outscphon

 iphl1 = iphl1 + 1

end subroutine prtscphon
!!***

!****f* ABINIT/endprtscphon
!!
!! NAME
!! endprtscphon
!!
!! FUNCTION
!! Close the file units for scphon outputs
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUTS
!!
!! NOTES
!!  closes 2 file units contained in present module
!! 
!! PARENTS
!!      mkphbs,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine endprtscphon()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'endprtscphon'
!End of the abilint section

 implicit none

 close (phfrq_unit)
 if (outscphon > 0) close (phvec_unit)

end subroutine endprtscphon
!!***


end module m_prtscphon
!!***
