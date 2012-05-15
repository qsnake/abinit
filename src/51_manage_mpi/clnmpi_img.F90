!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_img
!! NAME
!!  clnmpi_img
!!
!! FUNCTION
!!  Cleans-up the mpi informations for parallelism over images of the cell (npimage>1).
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MT,GG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      driver,invars1,invars2m,m_results_out
!!
!! CHILDREN
!!      mpi_comm_free,mpi_comm_size,mpi_group_free,mpi_group_size
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine clnmpi_img(mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi

#if defined HAVE_MPI && defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_img'
!End of the abilint section

 implicit none
#if defined HAVE_MPI && defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
#if defined HAVE_MPI
         integer :: ierr,csize
#endif

! ***********************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_MPI


 if (mpi_enreg%group_one_img/=mpi_enreg%world_group.and.mpi_enreg%group_one_img/=xmpi_group_null) then
   call MPI_GROUP_SIZE(mpi_enreg%group_one_img,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_GROUP_FREE(mpi_enreg%group_one_img,ierr)
   end if
 end if
 mpi_enreg%group_one_img=mpi_enreg%world_group

 if (mpi_enreg%comm_one_img/=xmpi_world.and.mpi_enreg%comm_one_img/=xmpi_comm_null.and. &
& mpi_enreg%comm_one_img/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%comm_one_img,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%comm_one_img,ierr)
   end if
 end if
 mpi_enreg%comm_one_img=xmpi_world

 if (mpi_enreg%comm_img/=xmpi_world.and.mpi_enreg%comm_img/=xmpi_comm_null.and. &
& mpi_enreg%comm_img/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%comm_img,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%comm_img,ierr)
   end if
 end if
 mpi_enreg%comm_img=xmpi_comm_null
#endif

 if (associated(mpi_enreg%index_img))  then
   ABI_DEALLOCATE(mpi_enreg%index_img)
 end if
 if (associated(mpi_enreg%distrb_img))  then
   ABI_DEALLOCATE(mpi_enreg%distrb_img)
 end if

 mpi_enreg%paral_img=0
 mpi_enreg%nimage=1
 mpi_enreg%me_img=0
 mpi_enreg%me_one_img=0

 DBG_EXIT("COLL")

end subroutine clnmpi_img
!!***
