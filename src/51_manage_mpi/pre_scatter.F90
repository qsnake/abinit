!{\src2tex{textfont=tt}}
!!****f* ABINIT/pre_scatter
!!
!! NAME
!!  pre_scatter
!!
!! FUNCTION
!!  Scatters data among FFT processors.
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2012 ABINIT group (?,YP)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  n1,n2,n3= FFT grid dimensions
!!  n4= n3/mpi_enreg%nproc_fft
!!  array_allgather= data to scatter among FFT procs
!!
!! OUTPUT
!!  array= scattered data
!!
!! SIDE EFFECTS
!!  mpi_enreg= information about MPI parallelization
!!
!! PARENTS
!!      fresid
!!
!! CHILDREN
!!      xcomm_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine pre_scatter(array,array_allgather,n1,n2,n3,n4,mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pre_scatter'
 use interfaces_51_manage_mpi, except_this_one => pre_scatter
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,n3,n4
 real(dp),intent(out) :: array(n1,n2,n4,1)
 real(dp),intent(in) :: array_allgather(n1,n2,n3,1)
 type(mpi_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: spaceComm
 integer ::  old_paral_level

! *********************************************************************

 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
!Perform the reverse operation
 array(:,:,:,:)=array_allgather(:,:,&
& n3/mpi_enreg%nproc_fft*mpi_enreg%me_fft+1:n3/mpi_enreg%nproc_fft*(mpi_enreg%me_fft+1),:)
 mpi_enreg%paral_level=old_paral_level
end subroutine pre_scatter
!!***
