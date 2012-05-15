!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_respfn
!! NAME
!!  clnmpi_respfn
!!
!! FUNCTION
!!  Cleans-up the mpi informations for parallelization over perturbations.
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2012 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  spaceComm_noparalrf=Communicator without the parallelization over perturbations
!!                      (to be put in mpi_enreg%comm_respfn in this routine)
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      abi_io_redirect,leave_test,mpi_comm_free,mpi_group_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine clnmpi_respfn(mpi_enreg,spaceComm_noparalrf)

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
#define ABI_FUNC 'clnmpi_respfn'
 use interfaces_51_manage_mpi, except_this_one => clnmpi_respfn
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 integer,intent(in) :: spaceComm_noparalrf
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: ierr,igroup_cnt,ngroups,spaceComm
 character(len=500) :: msg
#endif

! ***********************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_MPI
 if(mpi_enreg%paral_compil_respfn == 1) then

!  Reset communicators
   spaceComm=mpi_enreg%comm_respfn
   mpi_enreg%comm_respfn=spaceComm_noparalrf
   call abi_io_redirect(new_leave_comm=spaceComm_noparalrf)

!  Test whether _all_ groups are working properly
   call leave_test()

!  Free groups
   if (spaceComm/=MPI_COMM_NULL.and.spaceComm/=MPI_COMM_SELF) then
     ngroups=size(mpi_enreg%respfn_comm)
     if (ngroups>0) then
       do igroup_cnt=1,ngroups
         if (mpi_enreg%respfn_comm(igroup_cnt)==spaceComm) then
           call MPI_COMM_FREE(mpi_enreg%respfn_comm(igroup_cnt),ierr)
           if (ierr/=MPI_SUCCESS) then
             write(unit=msg,fmt='(a,i3)') '  Error on releasing Communicator for group nr ',igroup_cnt
             MSG_BUG(msg)
           end if
           call MPI_GROUP_FREE(mpi_enreg%respfn_group(igroup_cnt),ierr)
           if (ierr/=MPI_SUCCESS) then
             write(unit=msg,fmt='(a,i3)') '  Error on releasing group nr ',igroup_cnt
             MSG_BUG(msg)
           end if
         end if
       end do
     end if
   end if

   ABI_DEALLOCATE(mpi_enreg%respfn_group)
   ABI_DEALLOCATE(mpi_enreg%respfn_comm)
   mpi_enreg%me_respfn=-1
 end if
#endif

 DBG_EXIT("COLL")

end subroutine clnmpi_respfn
!!***
