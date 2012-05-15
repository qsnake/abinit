!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_gs
!! NAME
!!  initmpi_gs
!!
!! FUNCTION
!!  Initializes the MPI information for the ground-state datasets.
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2012 ABINIT group (AR, XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!    mpi_enreg%paralbd=option for parallelisation over the bands
!!    If the cpp option MPI is activated, also initialize
!!      mpi_enreg%proc_distrb(nkpt,mband,nsppol)  array that
!!        describes the work of each processor (allocated here)
!!      mpi_enreg%nproc=number of processors
!!      mpi_enreg%me=my processor number
!!
!! TODO
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      distrb2,mpi_comm_create,mpi_comm_rank,mpi_comm_size,mpi_group_free
!!      mpi_group_incl
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine initmpi_gs(dtset,mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_gs'
 use interfaces_51_manage_mpi, except_this_one => initmpi_gs
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
#if defined HAVE_MPI
 integer :: kpt_group,mband,nbdblock,nkpt,nsppol
!Variables introduced for MPI version
 integer :: group,ierr,iikpt,iisppol,irank,max_proc,me,min_proc
 integer :: nband_k,nproc,spacecomm,spacegroup
 integer,allocatable :: ranks(:)
#endif

! ***********************************************************************

!DEBUG
!write(std_out,*)' initmpi_gs : enter'
!stop
!ENDDEBUG

 mpi_enreg%paralbd=0
 mpi_enreg%paral_level=2
 
!Set up information for FFT parallelism
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1
 if(dtset%paral_kgb==1) mpi_enreg%nproc_kpt=dtset%npkpt
 if(dtset%paral_kgb==0) mpi_enreg%paral_fft=0
!MPIWF Should be modified for real FFT parallelism ...

#if defined HAVE_MPI

!Take into account possible parallelization over images
 if (mpi_enreg%paral_img==0) then
   spacecomm=mpi_enreg%comm_one_img
   spacegroup=mpi_enreg%group_one_img
   me=mpi_enreg%me_one_img
   nproc=mpi_enreg%nproc_one_img
 else
   spacecomm=mpi_enreg%world_comm
   spacegroup=mpi_enreg%world_group
   me=mpi_enreg%me
   nproc=mpi_enreg%nproc
 end if

 if(dtset%paral_kgb==1) nproc=mpi_enreg%nproc_kpt

 mband=dtset%mband
 nbdblock=dtset%nbdblock
 nkpt=dtset%nkpt
 nsppol=dtset%nsppol
 ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt,mband,nsppol))
 if (nkpt >= nproc) then
   mpi_enreg%paralbd=0
 else
   if (nbdblock == 1) then
     mpi_enreg%paralbd=0
   else
     mpi_enreg%paralbd=nbdblock
   end if
 end if

 if(dtset%paral_kgb ==1) mpi_enreg%proc_distrb(:,:,:)=0 ! temporary

!This routine sets up the array mpi_enreg%proc_distrb, and also mpi_enreg%me
 call distrb2(mband, dtset%nband, nkpt, nsppol, mpi_enreg)

 if (mpi_enreg%paralbd >= 1) then
!  Creation of groups of communicators
   ABI_ALLOCATE(mpi_enreg%kpt_comm,(nkpt*nsppol))
   ABI_ALLOCATE(ranks,(mpi_enreg%nproc_per_kpt))
   do iisppol=1,nsppol
     do iikpt=1,nkpt
       group=iikpt+(iisppol-1)*nkpt
       nband_k=dtset%nband(iikpt+(iisppol-1)*nkpt)
       min_proc=minval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
       max_proc=maxval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
       do irank=1,mpi_enreg%nproc_per_kpt
         ranks(irank)=min_proc+irank-1
         if (ranks(irank)==mpi_enreg%me) then
           mpi_enreg%num_group=group
         end if
       end do
       call MPI_GROUP_INCL(spacegroup,mpi_enreg%nproc_per_kpt,ranks,kpt_group,ierr)
       call MPI_COMM_CREATE(spacecomm,kpt_group,mpi_enreg%kpt_comm(group),ierr)
       call MPI_GROUP_FREE(kpt_group,ierr)
     end do
   end do
   if ((nproc > mpi_enreg%nproc_per_kpt*nkpt*nsppol) .and. &
&   (me >= mpi_enreg%nproc_per_kpt*nkpt*nsppol)) then
     mpi_enreg%num_group=0
     mpi_enreg%me_group=-1
     mpi_enreg%nproc_group=-1
   else
     call MPI_COMM_RANK(mpi_enreg%kpt_comm(mpi_enreg%num_group), &
&     mpi_enreg%me_group,ierr)
     call MPI_COMM_SIZE(mpi_enreg%kpt_comm(mpi_enreg%num_group), &
&     mpi_enreg%nproc_group,ierr)
   end if
   ABI_DEALLOCATE(ranks)
 end if

#endif

!DEBUG
!write(std_out,*)' initmpi_gs : exit '
!stop
!ENDDEBUG

end subroutine initmpi_gs
!!***
