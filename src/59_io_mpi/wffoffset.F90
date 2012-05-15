!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffOffset
!! NAME
!! WffOffset
!!
!! FUNCTION
!! Tool to manage WF file in the MPI/IO case : broadcast the offset of
!! the first k-point data block
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  wff <type(wffile_type)> = structured info about the wavefunction file
!!  sender = id of the sender
!!  spaceComm = id of the space communicator handler
!!
!! OUTPUT
!!  ier = error code returned by the MPI call
!!
!! PARENTS
!!      outwf
!!
!! CHILDREN
!!      mpi_bcast,xmax_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine WffOffset(wff,sender,spaceComm,ier)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile
 use m_errors
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffOffset'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer          ,intent(inout) :: sender
 integer          ,intent(in)    :: spaceComm
 integer          ,intent(out)   :: ier

!Local variables ------------------------------
#if defined HAVE_MPI_IO
 integer :: icom
 integer(kind=MPI_OFFSET_KIND) :: ima
#endif

! *********************************************************************

#if defined HAVE_MPI_IO
 if (wff%accesswff == IO_MODE_MPI) then
   call xmax_mpi(sender,icom,spaceComm,ier)
   if (icom>=0)then
     ima=wff%offwff
     call MPI_BCAST(ima,1,wff%offset_mpi_type,icom,spaceComm,ier)
     wff%offwff=ima
   end if
 end if ! accesswff
#else
 ier = 0
 ABI_UNUSED((/wff%accesswff,sender,spaceComm/))
#endif

end subroutine WffOffset

!!***
