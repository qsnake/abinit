!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffClose
!! NAME
!! WffClose
!!
!! FUNCTION
!! This subroutine closes a Wf file.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! PARENTS
!!      conducti_nc,conducti_paw,conducti_paw_core,emispec_paw,gstate,ioarr
!!      kss2wfk,linear_optics_paw,loop3dte,loper3,m_wfs,nonlinear,nstdy3
!!      nstpaw3,optic,optics_paw,optics_paw_core,optics_vloc,outwf,respfn,rrho
!!      scfcv3,suscep,uderiv,wffile,wfk_read_ene
!!
!! CHILDREN
!!      etsf_io_low_close,etsf_io_low_error_to_str,mpi_file_close
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine WffClose(wff,ier)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile
 use m_errors
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffClose'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type), intent(inout) :: wff
 integer, intent(out) :: ier

!Local ------------------------------------
#if defined HAVE_TRIO_ETSF_IO
  type(etsf_io_low_error) :: error
  logical                 :: lstat
  character(len = etsf_io_low_error_len)   :: errmess
#endif

! *************************************************************************

 ier=0
 if(wff%accesswff==IO_MODE_FORTRAN) then ! All processors see a local file
   close(unit=wff%unwff)

#if defined HAVE_TRIO_ETSF_IO
 else if(wff%accesswff == IO_MODE_ETSF)then

   call etsf_io_low_close(wff%unwff, lstat, error_data = error)
   if (.not. lstat) then
     call etsf_io_low_error_to_str(errmess, error)
     MSG_ERROR(errmess)
   end if

#endif
 else if(wff%accesswff==IO_MODE_FORTRAN_MASTER)then !  Only the master processor see a local file
   if(wff%master==wff%me) close (unit=wff%unwff)    ! VALGRIND complains buf points to uninitialized bytes

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then
   call MPI_FILE_CLOSE(wff%fhwff,ier)
   if (wff%master==wff%me ) close(unit=wff%unwff)
   wff%offwff=0;wff%off_recs=0;wff%lght_recs=0
   wff%nbOct_recMarker=-1
   wff%kgwff=-1
#endif

 end if

end subroutine WffClose
!!***
