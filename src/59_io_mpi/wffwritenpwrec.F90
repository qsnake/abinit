!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffWriteNpwRec
!! NAME
!! WffWriteNpwRec
!!
!! FUNCTION
!! This subroutine writes the npw record of a wavefunction file
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (XG,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!! nband_disk=number of bands
!! npw=number of plane waves
!! nspinor=number of spinorial components of the wavefunctions
!! opt_paral=(optional argument, default=1, only used for MPI-IO)
!!           1: all procs in the communicator write the data
!!           2: only master in the communicator writes the data
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      rwwf,vtowfk3
!!
!! CHILDREN
!!      mpi_barrier,mpi_bcast,xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine WffWriteNpwRec(ierr,nband_disk,npw,nspinor,wff,&
&                         opt_paral) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_wffile
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffWriteNpwRec'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: nband_disk,npw,nspinor
 integer,intent(in),optional :: opt_paral
 integer,intent(out) :: ierr

!Local variables-------------------------------
 integer :: opt_paral_
#if defined HAVE_MPI_IO
 integer :: me
#endif
  
! *************************************************************************

 ierr=0
 opt_paral_=1;if (present(opt_paral)) opt_paral_=opt_paral

 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff ==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then
!  $if (wff_ireadf90(wff)) then
   write(wff%unwff,iostat=ierr) npw,nspinor,nband_disk

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then
   me=-1;if (opt_paral_==2) me=wff%me_mpiio
   if ((me==-1.and.opt_paral_==1).or.(me==0.and.opt_paral_==2)) then
     call xderiveWRecInit(wff,ierr)
     call xderiveWrite(wff,npw,ierr)
     call xderiveWrite(wff,nspinor,ierr)
     call xderiveWrite(wff,nband_disk,ierr)
     call xderiveWRecEnd(wff,ierr)
   end if
   if (opt_paral_==2.and.wff%spaceComm_mpiio/=MPI_COMM_SELF) then
     call MPI_BARRIER(wff%spaceComm_mpiio,ierr)
     call MPI_BCAST(wff%offwff,1,wff%offset_mpi_type,0,wff%spaceComm_mpiio,ierr)
   end if
#endif
 end if

end subroutine WffWriteNpwRec
!!***

