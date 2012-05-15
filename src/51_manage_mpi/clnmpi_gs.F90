!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_gs
!! NAME
!!  clnmpi_gs
!!
!! FUNCTION
!!  Cleans-up the mpi informations for the ground-state datasets
!!  (mostly deallocate parts of mpi_enreg).
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2012 ABINIT group (AR, XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!    If the cpp option MPI is activated, deallocate:
!!      - mpi_enreg%proc_distrb
!!      - mpi_enreg%kpt_comm
!!
!! TODO
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      mpi_comm_free,mpi_group_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine clnmpi_gs(mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_gs'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: group,ierr,nkptspol
#endif

! ***********************************************************************

!DEBUG
!write(std_out,*)' clnmpi_gs : enter'
!stop
!ENDDEBUG

#if defined HAVE_MPI
 if (associated(mpi_enreg%proc_distrb))  then
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 end if
 if (associated(mpi_enreg%kpt_comm)) then
   if (mpi_enreg%paralbd >= 1) then
     nkptspol=size(mpi_enreg%kpt_comm)
     do group=1,nkptspol
       if (mpi_enreg%kpt_comm(group)/=MPI_COMM_NULL.and.mpi_enreg%kpt_comm(group)/=MPI_COMM_SELF) then
         call MPI_COMM_FREE(mpi_enreg%kpt_comm(group),ierr)
       end if
     end do
     call MPI_GROUP_FREE(mpi_enreg%world_group,ierr)
   end if
   ABI_DEALLOCATE(mpi_enreg%kpt_comm)
 end if
#endif

end subroutine clnmpi_gs
!!***
