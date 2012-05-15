!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_fft
!! NAME
!!  clnmpi_fft
!!
!! FUNCTION
!!  Cleans-up the mpi informations for FFT or BAND-FFT parallelism
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
!!   If the cpp option MPI is activated, deallocate
!!   mpi_enreg%fft_comm
!!
!! TODO
!!
!! PARENTS
!!      gstate,invars2m,respfn
!!
!! CHILDREN
!!      clnmpi_grid,mpi_comm_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine clnmpi_fft(mpi_enreg)

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
#define ABI_FUNC 'clnmpi_fft'
 use interfaces_51_manage_mpi, except_this_one => clnmpi_fft
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
 integer :: ierr,group
#endif

! ***********************************************************************

!DEBUG
!write(std_out,*)' clnmpi_fft : enter'
!ENDDEBUG

#if defined HAVE_MPI
 if (mpi_enreg%fft_master_comm/=MPI_COMM_NULL.and.mpi_enreg%fft_master_comm/=MPI_COMM_SELF) then
   call MPI_COMM_FREE(mpi_enreg%fft_master_comm,ierr)
 end if

 if (associated(mpi_enreg%fft_comm)) then
   do group=1,size(mpi_enreg%fft_comm)
     if (mpi_enreg%fft_comm(group)/=MPI_COMM_NULL.and.mpi_enreg%fft_comm(group)/=MPI_COMM_SELF) then
       call MPI_COMM_FREE(mpi_enreg%fft_comm(group),ierr)
     end if
   end do
   ABI_DEALLOCATE(mpi_enreg%fft_comm)
 end if
#endif

 if (associated(mpi_enreg%nplanes_fft))  then
   ABI_DEALLOCATE(mpi_enreg%nplanes_fft)
 end if
 if (associated(mpi_enreg%ind_fft_planes))  then
   ABI_DEALLOCATE(mpi_enreg%ind_fft_planes)
 end if

 call clnmpi_grid(mpi_enreg)

!DEBUG
!write(std_out,*)' clnmpi_fft : end'
!ENDDEBUG

end subroutine clnmpi_fft
!!***
