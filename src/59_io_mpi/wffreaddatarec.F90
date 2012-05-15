!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffReadDataRec
!! NAME
!! WffReadDataRec
!!
!! FUNCTION
!! Generic subroutines to read data in one record of a wavefunction file
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (XG,MB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!!
!! NOTES
!!  We use several procedures with the same generic name
!!    WffReadDataRec  contains
!!               WffReadDataRec_dp   : read double precision 1D array
!!               WffReadDataRec_dp2d : read double precision 2D array
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!***


!!****f* ABINIT/WffReadDataRec_dp1d
!! NAME
!! WffReadDataRec_dp1d
!!
!! FUNCTION
!! Subroutine to read data in one record of a wavefunction file
!! Handles double precision 1D arrays
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (XG,MB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.

!! INPUTS
!! ndp=size of the double precision array to be read
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! dparray=array of double precision numbers
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      xderiveread,xderiverrecend,xderiverrecinit
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffReadDataRec_dp1d(dparray,ierr,ndp,wff)

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
#define ABI_FUNC 'WffReadDataRec_dp1d'
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
 real(dp),intent(out) :: dparray(ndp)

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   read (wff%unwff,iostat=ierr) dparray(1:ndp)

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then
   call xderiveRRecInit(wff,ierr)
   call xderiveRead(wff,dparray,ndp,MPI_COMM_SELF,ierr)
   call xderiveRRecEnd(wff,ierr)
#endif
 end if

end subroutine WffReadDataRec_dp1d
!!***


!!****f* ABINIT/WffReadDataRec_dp2d
!! NAME
!! WffReadDataRec_dp2d
!!
!! FUNCTION
!! Subroutine to read data in one record of a wavefunction file
!! Handles double precision 2D arrays
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (XG,MB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.

!! INPUTS
!! n1,n2=sizes of the double precision array to be read
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! dparray=array of double precision numbers
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      xderiveread,xderiverrecend,xderiverrecinit
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine WffReadDataRec_dp2d(dparray,ierr,n1,n2,wff)

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
#define ABI_FUNC 'WffReadDataRec_dp2d'
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
 real(dp),intent(out) :: dparray(n1,n2)

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
   read (wff%unwff,iostat=ierr) dparray(1:n1,1:n2)

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then
   call xderiveRRecInit(wff,ierr)
   call xderiveRead(wff,dparray,n1,n2,MPI_COMM_SELF,ierr)
   call xderiveRRecEnd(wff,ierr)
#endif
 end if

end subroutine WffReadDataRec_dp2d
!!***
