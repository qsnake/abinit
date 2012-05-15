!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffWriteDataRec
!! NAME
!! WffWriteDataRec
!!
!! FUNCTION
!! Generic subroutines to write data in one record of a wavefunction file
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (XG,MB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
!!  We use several procedures with the same generic name
!!    WffWriteDataRec  contains
!!               WffWriteDataRec_int2d: write integer 2D array
!!               WffWriteDataRec_dp1d : write double precision 1D array
!!               WffWriteDataRec_dp2d : write double precision 2D array
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!***


!!****f* ABINIT/WffWriteDataRec_int2d
!! NAME
!! WffWriteDataRec_int2d
!!
!! FUNCTION
!! Subroutine to write data in one record of a wavefunction file
!! Handles integer 2D arrays
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (XG,MB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! intarray=array of integer numbers
!! n1,n2=sizes of the integer array to be written
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffWriteDataRec_int2d(intarray,ierr,n1,n2,wff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_wffile
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffWriteDataRec_int2d'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) ::  n1,n2
 integer,intent(out) :: ierr
 integer,intent(in) :: intarray(n1,n2)

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   write(wff%unwff,iostat=ierr) intarray(1:n1,1:n2)

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then
   call xderiveWRecInit(wff,ierr)
   call xderiveWrite(wff,intarray,n1,n2,MPI_COMM_SELF,ierr)
   call xderiveWRecEnd(wff,ierr)
#endif
 end if

end subroutine WffWriteDataRec_int2d
!!***


!!****f* ABINIT/WffWriteDataRec_dp1d
!! NAME
!! WffWriteDataRec_dp1d
!!
!! FUNCTION
!! Subroutine to write data in one record of a wavefunction file
!! Handles double precision 1D arrays
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (XG,MB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! dparray=array of double precision numbers
!! ndp=size of the double precision array to be written
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffWriteDataRec_dp1d(dparray,ierr,ndp,wff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_wffile
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffWriteDataRec_dp1d'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) ::  ndp
 integer,intent(out) :: ierr
 real(dp),intent(in) :: dparray(ndp)

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   write(wff%unwff,iostat=ierr) dparray(1:ndp)

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then
   call xderiveWRecInit(wff,ierr)
   call xderiveWrite(wff,dparray,ndp,MPI_COMM_SELF,ierr)
   call xderiveWRecEnd(wff,ierr)
#endif
 end if

end subroutine WffWriteDataRec_dp1d
!!***


!!****f* ABINIT/WffWriteDataRec_dp2d
!! NAME
!! WffWriteDataRec_dp2d
!!
!! FUNCTION
!! Subroutine to write data in one record of a wavefunction file
!! Handles double precision 2D arrays
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (XG,MB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! dparray=array of double precision numbers
!! n1,n2=sizes of the double precision array to be written
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine WffWriteDataRec_dp2d(dparray,ierr,n1,n2,wff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_wffile
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffWriteDataRec_dp2d'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) ::  n1,n2
 integer,intent(out) :: ierr
 real(dp),intent(in) :: dparray(n1,n2)

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   write(wff%unwff,iostat=ierr) dparray(1:n1,1:n2)

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then
   call xderiveWRecInit(wff,ierr)
   call xderiveWrite(wff,dparray,n1,n2,MPI_COMM_SELF,ierr)
   call xderiveWRecEnd(wff,ierr)
#endif
 end if

end subroutine WffWriteDataRec_dp2d
!!***
