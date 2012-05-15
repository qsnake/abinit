!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffDelete
!! NAME
!! WffDelete
!!
!! FUNCTION
!! This subroutine closes a Wf file, and delete it.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! PARENTS
!!      gstate,loper3,outwf,respfn
!!
!! CHILDREN
!!      mpi_file_close
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine WffDelete(wff,ier)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile

#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffDelete'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer, intent(out) :: ier

!Local variables-------------------------------

! *************************************************************************

 ier=0
 if (wff%accesswff==IO_MODE_FORTRAN) then !  All processors see a local file
   close(unit=wff%unwff,status='delete')

 else if (wff%accesswff==IO_MODE_FORTRAN_MASTER)then !  Only the master processor see a local file
   if (wff%master==wff%me) close (unit=wff%unwff,status='delete')

#if defined HAVE_MPI_IO
 else if (wff%accesswff==IO_MODE_MPI)then
   if ( wff%fhwff /= -1 )then
     call MPI_FILE_CLOSE(wff%fhwff,ier)
   end if
   if (wff%master==wff%me ) then
     close(unit=wff%unwff,status='delete')
     wff%fhwff = -1
   end if
   wff%offwff=0;wff%off_recs=0;wff%lght_recs=0
   wff%nbOct_recMarker=-1
   wff%kgwff=-1
#endif

 end if

end subroutine WffDelete
!!***
