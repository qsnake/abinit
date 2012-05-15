!{\src2tex{textfont=tt}}
!!****f* ABINIT/init9
!!
!! NAME
!! init9
!!
!! FUNCTION
!! Initialize the code ppddb9  :
!! write heading and make the first i/os
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!! character(len=fnlen) filnam(7)=character strings giving file names
!!
!! NOTES
!! 1. Should be executed by one processor only.
!! 2. File names refer to following files, in order:
!!     (1) Formatted input file
!!     (2) Formatted output file
!!     (3) Input Derivative Database
!!     (4) Output Molecular Dynamics
!!     (5) Input electron-phonon matrix elements
!!     (6) Root name for electron-phonon file names
!!     (7) Name of file containing the 3 ddk filenames and the GS wf file name
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine init9(filnam)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!arrays
 character(len=fnlen),intent(out) :: filnam(7)

!Local variables-------------------------------

! *********************************************************************

!Read the file names
 write(std_out,*)' Give name for      formatted input file : '
 read(05, '(a)' ) filnam(1)
 write(std_out,'(a,a)' )'-   ',trim(filnam(1))
 write(std_out,*)' Give name for     formatted output file : '
 read(05, '(a)' ) filnam(2)
 write(std_out,'(a,a)' )'-   ',trim(filnam(2))
 write(std_out,*)' Give name for input derivative database : '
 read(05, '(a)' ) filnam(3)
 write(std_out,'(a,a)' )'-   ',trim(filnam(3))
 write(std_out,*)' Give name for output molecular dynamics : '
 read(05, '(a)' ) filnam(4)
 write(std_out,'(a,a)' )'-   ',trim(filnam(4))
 write(std_out,*)' Give name for input elphon matrix elements  (GKK file) : '
 read(05, '(a)' ) filnam(5)
 write(std_out,'(a,a)' )'-   ',trim(filnam(5))
 write(std_out,*)' Give root name for elphon output files : '
 read(05, '(a)' ) filnam(6)
 write(std_out,'(a,a)' )'-   ',trim(filnam(6))
 write(std_out,*)' Give name for file containing ddk filenames for elphon/transport : '
 read(05, '(a)' ) filnam(7)
 write(std_out,'(a,a)' )'-   ',trim(filnam(7))

end subroutine init9
!!***
