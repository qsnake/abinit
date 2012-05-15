!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_grid
!! NAME
!!  clnmpi_grid
!!
!! FUNCTION
!!  Cleans-up the mpi informations for parallelism over grid (kpt/band/fft).
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
!!      clnmpi_fft,invars1
!!
!! CHILDREN
!!      mpi_comm_free,mpi_comm_size
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine clnmpi_grid(mpi_enreg)

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
#define ABI_FUNC 'clnmpi_grid'
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

 if (associated(mpi_enreg%sizecart))  then
   ABI_DEALLOCATE(mpi_enreg%sizecart)
 end if
 if (associated(mpi_enreg%coords))    then
   ABI_DEALLOCATE(mpi_enreg%coords)
 end if

 if (mpi_enreg%commcart_4d/=xmpi_world.and.mpi_enreg%commcart_4d/=xmpi_comm_null.and. &
& mpi_enreg%commcart_4d/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%commcart_4d,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%commcart_4d,ierr)
   end if
 end if

 if (mpi_enreg%commcart_3d/=xmpi_world.and.mpi_enreg%commcart_3d/=xmpi_comm_null.and. &
& mpi_enreg%commcart_3d/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%commcart_3d,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%commcart_3d,ierr)
   end if
 end if

 if (mpi_enreg%commcart/=xmpi_world.and.mpi_enreg%commcart/=xmpi_comm_null.and. &
& mpi_enreg%commcart/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%commcart,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%commcart,ierr)
   end if
 end if

 if (mpi_enreg%comm_spinfft/=xmpi_world.and.mpi_enreg%comm_spinfft/=xmpi_comm_null.and. &
& mpi_enreg%comm_spinfft/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%comm_spinfft,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%comm_spinfft,ierr)
   end if
 end if

 if (mpi_enreg%comm_bandspin/=xmpi_world.and.mpi_enreg%comm_bandspin/=xmpi_comm_null.and. &
& mpi_enreg%comm_bandspin/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%comm_bandspin,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%comm_bandspin,ierr)
   end if
 end if


 if (mpi_enreg%comm_fft/=xmpi_world.and.mpi_enreg%comm_fft/=xmpi_comm_null.and. &
& mpi_enreg%comm_fft/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%comm_fft,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%comm_fft,ierr)
   end if
 end if

 if (mpi_enreg%comm_band/=xmpi_world.and.mpi_enreg%comm_band/=xmpi_comm_null.and. &
& mpi_enreg%comm_band/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%comm_band,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%comm_band,ierr)
   end if
 end if

 if (mpi_enreg%comm_spin/=xmpi_world.and.mpi_enreg%comm_spin/=xmpi_comm_null.and. &
& mpi_enreg%comm_spin/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%comm_spin,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%comm_spin,ierr)
   end if
 end if

 if (mpi_enreg%comm_kpt/=xmpi_world.and.mpi_enreg%comm_kpt/=xmpi_comm_null.and. &
& mpi_enreg%comm_kpt/=xmpi_self) then
   call MPI_COMM_SIZE(mpi_enreg%comm_kpt,csize,ierr)
   if (ierr==0.and.csize>0) then
     call MPI_COMM_FREE(mpi_enreg%comm_kpt,ierr)
   end if
 end if

#endif

 DBG_EXIT("COLL")

end subroutine clnmpi_grid
!!***
