!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_atom
!! NAME
!!  clnmpi_atom
!!
!! FUNCTION
!!  Cleans-up the mpi informations for the parallelism over atoms (PAW).
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2012 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      gstate,respfn
!!
!! CHILDREN
!!      mpi_comm_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine clnmpi_atom(mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_atom'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
#if defined HAVE_MPI
 integer :: ierr
#endif

! ***********************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_MPI
 if (mpi_enreg%comm_atom/=MPI_COMM_NULL.and.mpi_enreg%comm_atom/=MPI_COMM_SELF) then
   call MPI_COMM_FREE(mpi_enreg%comm_atom,ierr)
 end if
#endif

 if (associated(mpi_enreg%atom_indx))  then
   ABI_DEALLOCATE(mpi_enreg%atom_indx)
 end if
 mpi_enreg%nproc_atom=1

 DBG_EXIT("COLL")

end subroutine clnmpi_atom
!!***
