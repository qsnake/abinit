!{\src2tex{textfont=tt}}
!!****f* ABINIT/leave_test
!! NAME
!!  leave_test
!!
!! FUNCTION
!!  Routine that tests whether exit must be done,
!!  because of eventual problems encountered by another processor.
!!  In this case, will make a clean exit.
!!  In the sequential case, return.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2012 ABINIT group (GMR, XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abinit,alloc_hamilt_gpu,clnmpi_respfn,dyfnl3,eltfrkin3,eltfrnl3,energy
!!      forstrnps,gstateimg,initylmg,inwffil3,iofn1,ladielmt,lavnl,memana,mkrho
!!      mkrho3,mlwfovlp,mlwfovlp_pw,newkpt,nselt3,nstdy3,nstpaw3,outwf
!!      partial_dos_fractions_paw,pawmkrhoij,prctfvw1,prctfvw2,pspheads_comm
!!      rhofermi3,scfcv,scfcv3,suscep_dyn,suscep_kxc_dyn,suscep_stat,tddft
!!      vtorho,vtorho3,wfsinp
!!
!! CHILDREN
!!      leave_myproc,mpi_allreduce,mpi_barrier,timab,wrtout,xmpi_abort
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine leave_test()

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi, only : xmpi_abort

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'leave_test'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: gl_check_bit,ierr,my_check_bit
 real(dp) :: tsec(2)
 character(len=500) :: message
#endif

! **********************************************************************

#if defined HAVE_MPI
 call timab(48,1,tsec)
!Synchronize
 call MPI_BARRIER(abinit_comm_leave,ierr)
 call timab(48,2,tsec)
 write(message, '(a)' ) ' leave_test : synchronization done...'
 call wrtout(std_out,message,'PERS')

!Everything is allright for me
 my_check_bit=0
 call timab(48,1,tsec)
!See what about the others
 call MPI_ALLREDUCE(my_check_bit,gl_check_bit,1,MPI_INTEGER,MPI_SUM,abinit_comm_leave,ierr)
 call timab(48,2,tsec)

!Check for exit
!MT 2011-11-28: when is this test activated ? (gl_check_bit>0) should be never true...
!(at least for static processes)
 if(gl_check_bit>0) then
   write(message, '(a,I6,a)' ) ' leave_test : error - ', gl_check_bit,&
&   ' processors are not answering. Exiting...'
   call wrtout(std_out,message,'PERS')
   call leave_myproc(option=1)
   call xmpi_abort()
!  call MPI_FINALIZE(ierr)
 end if
#endif

end subroutine leave_test
!!***
