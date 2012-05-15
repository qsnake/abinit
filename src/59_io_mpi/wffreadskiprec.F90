!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffReadSkipRec
!! NAME
!! WffReadSkipRec
!!
!! FUNCTION
!! This subroutine move forward or backward in a Wf file by nrec records.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (ZL,DCA,XG,GMR,MB,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nrec=number of records
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! TODO
!! For the future : one should treat the possible errors of backspace
!!
!! PARENTS
!!      gstate,nstdy3,nstpaw3,nstwf3,nstwf4,randac,rwwf,vtowfk3,wfkfermi3
!!
!! CHILDREN
!!      mvrecord,rwrecordmarker
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine WffReadSkipRec(ierr,nrec,wff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_wffile
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools,   only : mvrecord

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffReadSkipRec'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in)  :: nrec
 integer,intent(out) :: ierr
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: irec
 integer(kind=MPI_OFFSET_KIND) :: delim_record,offset
#endif

! *************************************************************************

 ierr=0
 if( wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) then

   call mvrecord(ierr,nrec,wff%unwff)

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then

   if (nrec>0) then ! Move forward nrec records
     do irec=1,nrec
       wff%off_recs = wff%offwff
       call rwRecordMarker(1,wff%offwff,delim_record,wff,ierr)
       wff%lght_recs = delim_record
     end do
   else             ! Move backward -nrec records
     do irec=1,-nrec
       offset = wff%offwff-wff%nbOct_recMarker
       call rwRecordMarker(1,offset,delim_record,wff,ierr)
       wff%lght_recs = delim_record
       wff%offwff = wff%offwff - delim_record - 2*wff%nbOct_recMarker
       wff%off_recs = wff%offwff
     end do
   end if
#endif
 end if ! wff%accesswff==0,1 or -1

end subroutine WffReadSkipRec
!!***
