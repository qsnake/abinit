!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_respfn
!! NAME
!!  initmpi_respfn
!!
!! FUNCTION
!!  Creates groups for Parallelization over Perturbations.
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2012 ABINIT group (franm, nimi)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!  spaceComm=Communicator out of which the Groups should be created
!!
!! OUTPUT
!!  spaceComm_noparalrf=Communicator without the parallelization over perturbations
!!                      (this is the value of the comm at entry of the routine)
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! TODO
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      abi_io_redirect,leave_new,mpi_comm_create,mpi_comm_group,mpi_comm_rank
!!      mpi_comm_size,mpi_group_incl,mpi_group_range_incl
!!      mpi_group_translate_ranks,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine initmpi_respfn(mpi_enreg,spaceComm,spaceComm_noparalrf)

 use m_profiling

 use defs_basis
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_respfn'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(MPI_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: spaceComm
 integer,intent(out) :: spaceComm_noparalrf

!Local variables-------------------------------
!no_abirules
#if defined HAVE_MPI
        integer:: igroup_respfn, nproc_spaceComm, me_spaceComm, nrest, spaceGroup
        integer:: range,ierr,master
        integer:: spaceRange(3,1)
        integer:: masterRank(1),masterRank2(1)
        integer,allocatable    :: masterGroup(:)
        character*500 :: message
#endif

! ***********************************************************************

!DEBUG
!write(std_out,*)' initmpi_respfn : enter'
!stop
!ENDDEBUG

!Save communicator
 spaceComm_noparalrf=mpi_enreg%comm_respfn

#if defined HAVE_MPI

!set master-id
 master=0

!get MPI envriroment
 call MPI_COMM_RANK(spaceComm,me_spaceComm,ierr)
 if(ierr /= MPI_SUCCESS) then
   write(message,'(2a,i4)') ch10, 'initmpi_respfn: Error on MPI_COMM_RANK: ', ierr
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 call MPI_COMM_SIZE(spaceComm,nproc_spaceComm,ierr)
 if(ierr /= MPI_SUCCESS) then
   write(message,'(2a,i4)') ch10, 'initmpi_respfn: Error on MPI_COMM_SIZE: ', ierr
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
!write(std_out,*) me_spaceComm, '/', nproc_spaceComm, ':', 'MPI_Comm_group...'
 call MPI_COMM_GROUP(spaceComm,spaceGroup,ierr)
 if(ierr /= MPI_SUCCESS) then
   write(message,'(2a,i4)') ch10, 'initmpi_respfn: Error on MPI_COMM_GROUP: ', ierr
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!write(std_out,*) me_spaceComm, '/', nproc_spaceComm, ':', 'Allocate Memory...'
 ABI_ALLOCATE(masterGroup,(mpi_enreg%ngroup_respfn))
 ABI_ALLOCATE(mpi_enreg%respfn_group,(mpi_enreg%ngroup_respfn))

!write(std_out,*) me_spaceComm, '/', nproc_spaceComm, ':', 'respfn_group allocated with size: ', sizeof(mpi_enreg%respfn_group)
 ABI_ALLOCATE(mpi_enreg%respfn_comm,(mpi_enreg%ngroup_respfn))

!write(std_out,*) me_spaceComm, '/', nproc_spaceComm, ':', 'respfn_comm allocated with size: ', sizeof(mpi_enreg%respfn_comm)

!calculate how many process per group and rest
 nrest=modulo(nproc_spaceComm,mpi_enreg%ngroup_respfn)
 range=nproc_SpaceComm / mpi_enreg%ngroup_respfn
!init spaceRange(2,1) with -1 so first spaceRange(1,1)=0
 spaceRange(2,1)=-1
!spaceRange(3,1)=1 so stepping is 1 -> select each process in range
 spaceRange(3,1)=1

!loop over all groups
 do igroup_respfn=1,mpi_enreg%ngroup_respfn
!  set first process of group = last process of last group + 1
   spaceRange(1,1)=spaceRange(2,1)+1
!  set last process of group
   spaceRange(2,1)=spaceRange(1,1)+range-1
!  add one process to group if nrest > 0
   if(nrest>0) then
     spaceRange(2,1)=spaceRange(2,1)+1
     nrest=nrest-1
   end if

!  create MPI group
   call MPI_GROUP_RANGE_INCL(spaceGroup,1,spaceRange,mpi_enreg%respfn_group(igroup_respfn),ierr)
   if(ierr /= MPI_SUCCESS) then
     write(message,'(2a,i4)') ch10, 'initmpi_respfn: Error on creating group nr. ',igroup_respfn
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
!  create MPI communicator
   call MPI_COMM_CREATE(spaceComm,mpi_enreg%respfn_group(igroup_respfn),mpi_enreg%respfn_comm(igroup_respfn),ierr)
   if(ierr /= MPI_SUCCESS) then
     write(message,'(2a,i4)') ch10, 'initmpi_respfn: Error while creating communicator from group ',igroup_respfn
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  if me is in actual range set my_respfn_group
   if((spaceRange(1,1)<=me_spaceComm).AND.(spaceRange(2,1)>=me_spaceComm)) then

     mpi_enreg%my_respfn_group=mpi_enreg%respfn_group(igroup_respfn)

     mpi_enreg%my_respfn_comm=mpi_enreg%respfn_comm(igroup_respfn)

!    number of processors in my group
     mpi_enreg%nproc_respfn = spaceRange(2,1) - spaceRange(1,1) + 1
!    number of my processor in my group
     call MPI_COMM_RANK(mpi_enreg%respfn_comm(igroup_respfn),mpi_enreg%me_respfn,ierr)
     if(ierr /= MPI_SUCCESS) then
       write(message,'(2a,i4)') ch10, 'initmpi_respfn: Error on MPI_COMM_RANK: ', ierr
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     mpi_enreg%comm_respfn = mpi_enreg%my_respfn_comm
   end if
 end do

!create MPI master group
 do igroup_respfn=1,mpi_enreg%ngroup_respfn
   masterRank(1)=0
   call MPI_GROUP_TRANSLATE_RANKS(mpi_enreg%respfn_group(igroup_respfn),1,masterRank,spaceGroup,masterRank2,ierr)
   masterGroup(igroup_respfn)=masterRank2(1)
 end do

 call MPI_GROUP_INCL(spaceGroup,mpi_enreg%ngroup_respfn,masterGroup,mpi_enreg%respfn_master_group,ierr)
 if(ierr /= MPI_SUCCESS) then
   write(message,'(2a)') ch10, 'initmpi_respfn: Error while creating master group '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 call MPI_COMM_CREATE(spaceComm,mpi_enreg%respfn_master_group,mpi_enreg%respfn_master_comm,ierr)

 if(ierr /= MPI_SUCCESS .or. mpi_enreg%respfn_master_group == MPI_COMM_NULL) then
   write(message,'(2a)') ch10, 'initmpi_respfn: Error while creating communicator from master group '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 ABI_DEALLOCATE(masterGroup)
#endif

!Has also to change the communicator used in leave_test
 call abi_io_redirect(new_leave_comm=mpi_enreg%comm_respfn)

!DEBUG
!write(std_out,*)' initmpi_respfn : exit '
!stop
!ENDDEBUG

end subroutine initmpi_respfn
!!***
