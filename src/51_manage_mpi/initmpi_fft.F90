!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_fft
!! NAME
!!  initmpi_fft
!!
!! FUNCTION
!!  Initializes the mpi information for FFT or BAND-FFT parallelism.
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2012 ABINIT group (AR, XG, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=informations about MPI parallelization
!!
!! OUTPUT
!!  Not up to date ; to be updated !
!!  mpi_enreg=informations about MPI parallelization
!!    mpi_enreg%fft_comm(nkpt)=comm array of FFT set
!!    mpi_enreg%me_fft=index of the processor in the FFT set
!!    mpi_enreg%nproc_fft=number of processors int the FFT set
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! PARENTS
!!      gstate,invars2m,respfn
!!
!! CHILDREN
!!      initmpi_grid,mpi_comm_create,mpi_comm_rank,mpi_comm_size,mpi_group_free
!!      mpi_group_incl
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_fft(dtset,mpi_enreg)

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
#define ABI_FUNC 'initmpi_fft'
 use interfaces_51_manage_mpi, except_this_one => initmpi_fft
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
#if defined HAVE_MPI
!scalars
 integer :: nkpt,nsppol
 character(len=500) :: msg
 integer :: fft_group,ierr,iikpt,iikpt_modulo,iproc,irank,isppol,me,nproc
 integer :: SpaceComm,SpaceGroup
 integer ::    nproc_tmp
!arrays
 integer,allocatable :: ranks(:)
#endif
! ***********************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_MPI
!This is to take into account possible parallelization over images of the cell
 if (mpi_enreg%paral_img==0) then
   me=mpi_enreg%me
   nproc=mpi_enreg%nproc
   SpaceComm=mpi_enreg%world_comm
   SpaceGroup=mpi_enreg%world_group
 else
   me=mpi_enreg%me_one_img
   nproc=mpi_enreg%nproc_one_img
   SpaceComm=mpi_enreg%comm_one_img
   SpaceGroup=mpi_enreg%group_one_img
 end if

!Set up information for FFT parallelism
 if(dtset%paral_kgb == 1) then
   mpi_enreg%mode_para='b'
   mpi_enreg%paralbd=0
   mpi_enreg%paral_fft=1
   mpi_enreg%fft_option_lob=1
   if (dtset%fft_opt_lob /= 0) mpi_enreg%fft_option_lob=dtset%fft_opt_lob
   mpi_enreg%nproc_fft=nproc
 else
   mpi_enreg%nproc_fft = 1
   mpi_enreg%num_group_fft=0
 end if

 nkpt=dtset%nkpt
 nsppol=dtset%nsppol

 if (nproc>0) then

!  Creation of groups of communicators
   if(dtset%paral_kgb == 1) then

!    ====== MT june 2011 : is this part useful ? Abinit works without it !
     ABI_ALLOCATE(mpi_enreg%fft_comm,(nkpt*nsppol))
     ABI_ALLOCATE(ranks,(mpi_enreg%nproc_fft))
     do isppol=1,nsppol
       do iikpt=1,nkpt
         do irank=1,mpi_enreg%nproc_fft
           ranks(irank)=irank-1
           if (ranks(irank)==me) mpi_enreg%num_group_fft=iikpt
         end do
         call MPI_GROUP_INCL(SpaceGroup,mpi_enreg%nproc_fft,ranks,fft_group,ierr)
         call MPI_COMM_CREATE(SpaceComm,fft_group,mpi_enreg%fft_comm(iikpt+(isppol-1)*nkpt),ierr)
         call MPI_GROUP_FREE(fft_group,ierr)
       end do
     end do
     call MPI_COMM_RANK(mpi_enreg%fft_comm(mpi_enreg%num_group_fft),mpi_enreg%me_fft,ierr)
     if (mpi_enreg%me_fft==0) then
       mpi_enreg%master_fft=me
     else
       mpi_enreg%master_fft=-1
     end if
     call MPI_COMM_SIZE(mpi_enreg%fft_comm(mpi_enreg%num_group_fft),mpi_enreg%nproc_fft,ierr)
     ABI_DEALLOCATE(ranks)
!    ===============================================

     if(mpi_enreg%mode_para=='b') then
       mpi_enreg%nproc_fft  = dtset%npfft
       mpi_enreg%nproc_band = dtset%npband
       mpi_enreg%nproc_kpt  = dtset%npkpt
       mpi_enreg%bandpp     = dtset%bandpp

       if((dtset%use_gpu_cuda==1).and.(mpi_enreg%nproc_fft/=1))then
         write(msg,'(3a,i5)') &
&         '  When use_gpu_cuda is on, the number of FFT processors, npfft, must be 1',ch10,&
&         '  However, npfft=',mpi_enreg%nproc_fft
         MSG_ERROR(msg)
       end if

       if(modulo(dtset%ngfft(2),mpi_enreg%nproc_fft)/=0)then
         write(msg,'(5a,i5,a,i5)') &
&         '  The number of FFT processors, npfft, should be',ch10,&
&         '  a multiple of the number of ngfft(2).',ch10,&
&         '  However, npfft=',mpi_enreg%nproc_fft,' and ngfft(2)=',dtset%ngfft(2)
         MSG_BUG(msg)
       end if

       do iikpt=1,nkpt*nsppol
         iikpt_modulo = modulo(iikpt,nkpt)+1
         if ((dtset%istwfk(iikpt_modulo)==2).and.(dtset%ngfft(7)==401)) then
           if ((mpi_enreg%bandpp==0).or. &
           ((mpi_enreg%bandpp/=1).and.(modulo(mpi_enreg%bandpp,2)/=0))) then
             write(msg,'(3a,i5)') &
&             '  The number bandpp should be 1 or a multiple of 2',ch10,&
&             '  However, bandpp=',mpi_enreg%bandpp
             MSG_BUG(msg)
           end if
           if(modulo(dtset%nband(iikpt),mpi_enreg%nproc_band*mpi_enreg%bandpp)/=0)then
             write(msg,'(5a,i5,a,i5)') &
&             '  The number of band for the k-point, nband_k, should be',ch10,&
&             '  a multiple of the number nproc_band*bandpp.',ch10,&
&             '  However, nband_k=',dtset%nband(iikpt),' and nproc_band*bandpp=', &
&             mpi_enreg%nproc_band* mpi_enreg%bandpp
             MSG_BUG(msg)
           end if
         elseif ((dtset%istwfk(iikpt_modulo)==2) .and. (dtset%ngfft(7)==400)) then
           msg='  The fftalg=400 with istwfk=2 is not valid'
           MSG_BUG(msg)
         else
           if(modulo(dtset%nband(iikpt),mpi_enreg%nproc_band*mpi_enreg%bandpp)/=0)then
             write(msg,'(5a,i5,a,i5)') &
&             '  The number of band for the k-point, nband_k, should be',ch10,&
&             '  a multiple of the number nproc_band*bandpp.',ch10,&
&             '  However, nband_k=',dtset%nband(iikpt),' and nproc_band*bandpp=', &
&             mpi_enreg%nproc_band* mpi_enreg%bandpp
             MSG_BUG(msg)
           end if
           if ((mpi_enreg%bandpp==0)) then
             write(msg,'(a,i5,2a,i5,2a,i5)')&
&             '  The number bandpp should not be 0 with fftalg=',dtset%ngfft(7),ch10,&
&             ' and istwfk=',dtset%istwfk(iikpt_modulo),ch10,&
&             '  However, bandpp=',mpi_enreg%bandpp
             MSG_BUG(msg)
           end if
         end if
       end do

       if (mpi_enreg%paral_compil_kpt==1) then
         if(modulo(nkpt*nsppol,mpi_enreg%nproc_kpt)/=0)then
           write(msg,'(5a,i5,a,i5)') &
&           '  The number of KPT processors, npkpt, should be',ch10,&
&           '  a multiple of the number of nkpt*nsppol.',ch10,&
&           '  However, npkpt=',mpi_enreg%nproc_kpt,' and nkpt*nsppol=',nkpt*nsppol
           MSG_WARNING(msg)
         end if
       end if

!      Call to main routine to create cartesian grid of processors
       call initmpi_grid(mpi_enreg)

     else  ! mpi_enreg%mode_para/='b'
       mpi_enreg%nproc_fft   = 1
       mpi_enreg%nproc_band  = 1
       mpi_enreg%comm_fft    = MPI_COMM_SELF
       mpi_enreg%comm_band   = MPI_COMM_SELF
       mpi_enreg%comm_kpt    = SpaceComm
       mpi_enreg%me_fft      = 0
       mpi_enreg%me_band     = 0
       mpi_enreg%me_kpt      = me
       mpi_enreg%commcart    = MPI_COMM_SELF
       mpi_enreg%commcart_3d = SpaceComm
       mpi_enreg%me_cart_2d  = 0
       mpi_enreg%nproc_kpt   = nproc
       mpi_enreg%bandpp      = 1
     end if
   end if

 else ! nproc==0
   mpi_enreg%nproc_fft   = 0
   mpi_enreg%nproc_band  = 0
   mpi_enreg%nproc_kpt   = 0
   mpi_enreg%me_fft      = 0
   mpi_enreg%me_band     = 0
   mpi_enreg%me_kpt      = 0
   mpi_enreg%me_cart_2d  = 0
   mpi_enreg%comm_fft    = MPI_COMM_NULL
   mpi_enreg%comm_band   = MPI_COMM_NULL
   mpi_enreg%comm_kpt    = MPI_COMM_NULL
   mpi_enreg%commcart    = MPI_COMM_NULL
   mpi_enreg%commcart_3d = MPI_COMM_NULL
   mpi_enreg%bandpp      = 1
 end if

!Creation of communicator for master
 if (mpi_enreg%paral_fft==0.and.nproc>0) then
   ABI_ALLOCATE(ranks,(nproc))
   do iproc=0,nproc-1
     ranks(iproc+1)=iproc
   end do
!  This did not work with the g95 compiler XG081106: nproc was replaced by 0 . Very strange
!  call MPI_GROUP_INCL(SpaceGroup,nproc,ranks,fft_group,ierr)
   nproc_tmp=nproc
   call MPI_GROUP_INCL(SpaceGroup,nproc_tmp,ranks,fft_group,ierr)
   call MPI_COMM_CREATE(SpaceComm,fft_group,mpi_enreg%fft_master_comm,ierr)
   call MPI_GROUP_FREE(fft_group,ierr)
   ABI_DEALLOCATE(ranks)
 else if (nproc>0) then
!  only one proc per group_fft
   ABI_ALLOCATE(ranks,(1))
   ranks(1)=0
   call MPI_GROUP_INCL(SpaceGroup,1,ranks,fft_group,ierr)
   call MPI_GROUP_FREE(fft_group,ierr)
   ABI_DEALLOCATE(ranks)
   mpi_enreg%fft_master_comm=MPI_COMM_SELF
 else
   mpi_enreg%fft_master_comm=MPI_COMM_NULL
 end if

!MT june 2011 : these pieces of code are never used
!They seem to be there for a future? parallelization over FFT
!if (nproc>0.and.mpi_enreg%paral_compil_fft==1) then
!if(dtset%paral_kgb ==0) then
!allocate(mpi_enreg%fft_comm(nkpt*nsppol))
!mpi_enreg%master_fft=-1
!do isppol=1,nsppol
!do iikpt=1,nkpt
!iproc_min=minval(mpi_enreg%proc_distrb(iikpt,:,isppol))
!iproc_max=maxval(mpi_enreg%proc_distrb(iikpt,:,isppol))
!if (me==iproc_min) mpi_enreg%master_fft=me
!allocate(ranks(iproc_max-iproc_min+1))
!iiproc=1
!do iproc=iproc_min,iproc_max
!ranks(iiproc)=iproc
!iiproc=iiproc+1
!end do
!!       With MPI on SGI machine "Spinoza", there is a limitation
!!       of the number of groups that can be defined. When FFT // is off
!!       paral_kgb=0, there is actually no reason to define
!!       the following groups (in the present implementation). So,
!!       these lines can be skipped.
!#if !defined FC_MIPSPRO
!call MPI_GROUP_INCL(SpaceGroup,iproc_max-iproc_min+1,ranks,fft_group,ierr)
!call MPI_COMM_CREATE(SpaceComm,fft_group,mpi_enreg%fft_comm(iikpt+(isppol-1)*nkpt),ierr)
!call MPI_GROUP_FREE(fft_group,ierr)
!#endif
!deallocate(ranks)
!end do
!end do
!end if
!if (dtset%mgfft /=0) then
!allocate(mpi_enreg%nplanes_fft(dtset%nkpt))
!allocate(mpi_enreg%ind_fft_planes(dtset%nkpt,dtset%ngfft(2)))
!end if
!end if

#endif

 DBG_EXIT("COLL")

end subroutine initmpi_fft
!!***
