!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_world
!! NAME
!!  initmpi_world
!!
!! FUNCTION
!!  Initializes the mpi information for world.
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2012 ABINIT group (FJ, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! OUTPUT
!!
!!  mpi_enreg=informations about MPI parallelization
!! 
!!
!! SIDE EFFECTS
!! xmpi_world is redifined for the number of processors on which ABINIT is launched
!!
!! TODO
!!
!! PARENTS
!!      finddistrproc
!!
!! CHILDREN
!!      abi_io_redirect,mpi_comm_create,mpi_comm_group,mpi_comm_rank
!!      mpi_comm_size,mpi_group_free,mpi_group_incl
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_world(mpi_enreg,nproc)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_world'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in)::nproc
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
#if defined HAVE_MPI
!scalars
 integer :: new_world_group, world_group,ierr,ii

!arrays
 integer,allocatable :: ranks(:)
#endif
! ***********************************************************************

!DBG_ENTER("COLL")

#if defined HAVE_MPI
 if(nproc==mpi_enreg%nproc) return
!Creation of groups of communicators

 ABI_ALLOCATE(ranks,(0:nproc-1))
 do ii=0,nproc-1
   ranks(ii)=ii
 end do
 call MPI_COMM_GROUP(MPI_COMM_WORLD,world_group,ierr)
 call MPI_GROUP_INCL(world_group,nproc,ranks,new_world_group,ierr)
 call MPI_COMM_CREATE(MPI_COMM_WORLD,new_world_group,xmpi_world,ierr)
 call MPI_GROUP_FREE(world_group,ierr)
 call MPI_GROUP_FREE(new_world_group,ierr)
 ABI_DEALLOCATE(ranks)
 if(mpi_enreg%me<nproc)  then
   mpi_enreg%world_comm=xmpi_world
!  mpi_enreg%world_group=MPI_GROUP_NULL
   call MPI_COMM_GROUP(mpi_enreg%world_comm,mpi_enreg%world_group,ierr)
   call MPI_COMM_RANK(xmpi_world,mpi_enreg%me,ierr)
   call MPI_COMM_SIZE(xmpi_world,mpi_enreg%nproc,ierr)
   
   mpi_enreg%comm_respfn=mpi_enreg%world_comm
   call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)
 else
   mpi_enreg%me=-1
 end if
#endif
 
!DBG_EXIT("COLL")

end subroutine initmpi_world
!!***
