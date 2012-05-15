!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_band
!! NAME
!!  clnmpi_band
!!
!! FUNCTION
!!  Cleans-up the mpi informations for the BAND parallelism (paralbd=1).
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!    mpi_enreg%band_comm(nkpt*nsppol)=comm array of BAND set
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      mpi_comm_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine clnmpi_band(mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_band'
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
         integer :: group,ierr,nkptspol
#endif

! ***********************************************************************

!DEBUG
!write(std_out,*)' clnmpi_band : enter'
!ENDDEBUG

#if defined HAVE_MPI
 if (mpi_enreg%has_band_comm==1) then
   nkptspol=size(mpi_enreg%band_comm)
   do group=1,nkptspol
     if (mpi_enreg%band_comm(group)/=MPI_COMM_NULL.and.mpi_enreg%band_comm(group)/=MPI_COMM_SELF) then
       call MPI_COMM_FREE(mpi_enreg%band_comm(group),ierr)
     end if
   end do
 end if
#endif
 if (associated(mpi_enreg%band_comm))  then
   ABI_DEALLOCATE(mpi_enreg%band_comm)
 end if

!DEBUG
!write(std_out,*)' clnmpi_band : end'
!ENDDEBUG

end subroutine clnmpi_band
!!***
