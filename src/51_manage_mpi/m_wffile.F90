!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_wffile
!! NAME
!!  m_wffile
!!
!! FUNCTION
!!  This module provides the definition of the wffile_type used to WF file data.
!!  As the type contains MPI-dependent fields, it has to be declared in a MPI-managed directory.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!! wffile_type : a handler for dealing with the IO of a wavefunction file
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_wffile

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

 private

!public procedures.
 public :: getRecordMarkerLength_wffile
 public :: xnullifyOff
 public :: xderiveWRecEnd
 public :: xderiveWRecInit
 public :: xderiveRRecEnd
 public :: xderiveRRecInit
#if defined HAVE_MPI_IO
 public :: rwRecordMarker
#endif
 public :: clsopn
 public :: wff_usef90
!!***

!!****t* m_wffile/wffile_type
!! NAME
!! wffile_type
!!
!! FUNCTION
!! This structure datatype is a handler for dealing with the IO of a
!! wavefunction file.
!! It contains, among other things, the method of access to the file
!! (standard F90 read/write, or NetCDF call, or MPI IO), the unit number
!! if applicable, the filename, the information on the
!! parallelism, etc ...
!!
!! NOTES
!!
!! SOURCE

 type, public :: wffile_type

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar
  integer :: unwff
   ! unwff  unit number of unformatted wavefunction disk file

  integer :: accesswff
   ! Method to access the wavefunction file
   !   IO_MODE_FORTRAN for usual Fortran IO routines
   !   IO_MODE_FORTRAN_MASTER if usual Fortran IO routines, but only the master node in the parallel case
   !   IO_MODE_MPI if MPI/IO routines (this access method is only available in parallel)
   !   IO_MODE_NETCDF if NetCDF routines (obsolete, do not use)
   !   IO_MODE_ETSF, NetCDF format read via etsf-io.

  integer :: formwff
   ! formwff=format of the eigenvalues
   !   -1 => not used
   !    0 => vector of eigenvalues
   !    1 => hermitian matrix of eigenvalues

  integer :: headform
   ! headform=format of the header

  integer ::  kgwff
   ! kgwff  if 1 , read or write kg_k ; if 0, do not care about kg_k

! Character
  character(len=fnlen) :: fname
   ! filename (if available)

! In case of MPI parallel use
  integer :: master
   ! index of the processor master of the IO procedure when the WffOpen call is issued

  integer :: me
   ! index of my processor in the spaceComm communicator

  integer :: me_mpiio
   ! index of my processor in the spaceComm_mpiio communicator

  integer :: nproc
   ! number of processors that will have access to the file

  integer :: spaceComm
   ! space communicator for the standard FORTRAN access to the file

  integer :: spaceComm_mpiio
   ! space communicator for the MPI/IO access to the file

! In case of MPI/IO : additional information
#if defined HAVE_MPI_IO
  integer :: fhwff
   ! file handle used to access the file with MPI/IO.

  integer(kind=MPI_OFFSET_KIND) :: nbOct_int,nbOct_dp,nbOct_ch
   ! nbOct_int byte number of int value
   ! nbOct_dp  byte number of dp value
   ! nbOct_ch  byte number of character value

  integer(kind=MPI_OFFSET_KIND) :: nbOct_recMarker
   ! byte number of Fortran file record markers

  integer(kind=MPI_OFFSET_KIND) :: lght_recs
   ! length of record

  integer :: marker_mpi_type
   ! MPI Datatype for Fortran record markers

  integer(kind=MPI_OFFSET_KIND)  :: offwff,off_recs
   ! offwff   offset position of unformatted wavefunction disk file
   ! off_recs offset position of start record
   ! (used in parallel MPI-IO)

  integer :: offset_mpi_type
   ! MPI Datatype for INTEGER(kind=MPI_OFFSET_KIND)
#endif

 end type wffile_type


CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/getRecordMarkerLength_wffile
!! NAME
!!  getRecordMarkerLength_wffile
!!
!! FUNCTION
!!  Get the record marker length of the FORTRAN header of a file to access it in MPI/IO.
!!  This routine assumes that the header has been written (and flushed) in the file.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! SIDE EFFECTS
!!  wff=<type(wffile_type)>=structured info for reading/writing the wavefunctions
!!      only%nbOct_recMarker is changed
!!
!! PARENTS
!!      hdr_skip
!!
!! CHILDREN
!!      xnullifyoff
!!
!! SOURCE

subroutine getRecordMarkerLength_wffile(wff)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getRecordMarkerLength_wffile'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: headform,ierr,ii,iimax
 integer(kind=MPI_OFFSET_KIND)  :: posit,rml
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)
#endif

!************************************************************************

#if defined HAVE_MPI_IO

 if (wff%nbOct_recMarker>0) return

!wff%nbOct_recMarker=4;return
!call flush(wff%unwff)
!call MPI_FILE_SYNC(wff%fhwff,ierr)

!Only master do that
 ierr=0
 if (wff%master==wff%me) then

! Define number of INTEGER types to be tested
#if defined HAVE_FC_INT_QUAD
   iimax=4
#else
   iimax=3
#endif

! Try to read headform
   rml=-1;ii=0
   do while (wff%nbOct_recMarker<=0.and.ii<iimax)
     ii=ii+1
     if (ii==1) rml=4
     if (ii==2) rml=8
     if (ii==3) rml=2
     if (ii==4) rml=16
     posit=rml+6*wff%nbOct_ch
     call MPI_FILE_READ_AT(wff%fhwff,posit,headform,1,MPI_INTEGER,statux,ierr)
     if (ierr==MPI_SUCCESS) then
       if (headform==wff%headform) wff%nbOct_recMarker=rml
     end if
    end do

    if (ierr/=MPI_SUCCESS) then
     MSG_BUG("Header problem")
    end if

   if (ii==iimax.and.wff%nbOct_recMarker<=0) then
     if (iimax>=4) then
       write(msg,'(3a)') &
&        ' Your architecture is not able to handle 16, 8, 4 or 2-bytes FORTRAN file record markers !',ch10,&
&        ' You cannot use ABINIT and MPI/IO.'
     else
       write(msg,'(3a)') &
&        '  Your architecture is not able to handle 8, 4 or 2-bytes FORTRAN file record markers !',ch10,&
&        '  You cannot use ABINIT and MPI/IO.'
     end if
     MSG_ERROR(msg)
   else
     write(msg,'(a,i0)') &
&     '  MPI/IO accessing FORTRAN file header: detected record mark length=',wff%nbOct_recMarker
     MSG_COMMENT(msg)
   end if

 end if  ! me=master

!Broadcast record marker length
 if (wff%spaceComm/=MPI_COMM_SELF) then
   call MPI_BCAST(wff%nbOct_recMarker,1,wff%offset_mpi_type,wff%master,wff%spaceComm,ierr)
 end if

!Select MPI datatype for markers
 if (wff%nbOct_recMarker==4) then
   wff%marker_mpi_type=MPI_INTEGER4
 else if (wff%nbOct_recMarker==8) then
   wff%marker_mpi_type=MPI_INTEGER8
#if defined HAVE_FC_INT_QUAD
 else if (wff%nbOct_recMarker==16) then
   wff%marker_mpi_type=MPI_INTEGER16
#endif
 else if (wff%nbOct_recMarker==2) then
   wff%marker_mpi_type=MPI_INTEGER2
 end if

#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine getRecordMarkerLength_wffile
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/rwRecordMarker
!! NAME
!!  rwRecordMarker
!!
!! FUNCTION
!!  Read/Write a record marker in a FORTRAN file at a given file pointer position.
!!  This is needed to access data in a FORTRAN file with MPI/IO.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  option=1 for reading by current proc
!!         2 for writing by current proc
!!         3 for reading by all procs
!!         4 for writing by all procs
!!  posit= position of the MPI/IO file pointer
!!  wff=<type(wffile_type)>=structured info for reading/writing
!!     Use here only:
!!       wff%fhwff= handle of the MPI/IO file
!!       wff%nbOct_recMarker= length of Fortran record markers
!!
!! OUTPUT
!!  ierr= error code
!!
!! SIDE EFFECTS
!!  posit= position of the MPI/IO file pointer
!!         updated after the reading (with the length of the record)
!!  recordmarker= content of the record marker
!!
!! PARENTS
!!      hdr_skip,m_wffile,wffreadskiprec
!!
!! CHILDREN
!!      xnullifyoff
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine rwRecordMarker(option,posit,recordmarker,wff,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rwRecordMarker'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option
 integer(kind=MPI_OFFSET_KIND),intent(inout) :: posit,recordmarker
 integer,intent(out) :: ierr
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
!scalars
 integer*2  :: delim_record2
 integer*4  :: delim_record4
 integer*8  :: delim_record8
#if defined HAVE_FC_INT_QUAD
 integer*16 :: delim_record16
#endif
!character(len=500) :: msg
!arrays
 integer  :: statux(MPI_STATUS_SIZE)

!************************************************************************

 ierr=0

 if (option==1) then
   if (wff%nbOct_recMarker==4) then
     call MPI_FILE_READ_AT(wff%fhwff,posit,delim_record4 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record4
   else if (wff%nbOct_recMarker==8) then
     call MPI_FILE_READ_AT(wff%fhwff,posit,delim_record8 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record8
#if defined HAVE_FC_INT_QUAD
   else if (wff%nbOct_recMarker==16) then
     call MPI_FILE_READ_AT(wff%fhwff,posit,delim_record16,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record16
#endif
   else if (wff%nbOct_recMarker==2) then
     call MPI_FILE_READ_AT(wff%fhwff,posit,delim_record2 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record2
   else
     MSG_BUG('Wrong record marker length!')
   end if

 else if (option==2) then
   if (wff%nbOct_recMarker==4) then
     delim_record4 = recordmarker
     call MPI_FILE_WRITE_AT(wff%fhwff,posit,delim_record4 ,1,wff%marker_mpi_type,statux,ierr)
   else if (wff%nbOct_recMarker==8) then
     delim_record8 = recordmarker
     call MPI_FILE_WRITE_AT(wff%fhwff,posit,delim_record8 ,1,wff%marker_mpi_type,statux,ierr)
#if defined HAVE_FC_INT_QUAD
   else if (wff%nbOct_recMarker==16) then
     delim_record16 = recordmarker
     call MPI_FILE_WRITE_AT(wff%fhwff,posit,delim_record16,1,wff%marker_mpi_type,statux,ierr)
#endif
   else if (wff%nbOct_recMarker==2) then
     delim_record2 = recordmarker
     call MPI_FILE_WRITE_AT(wff%fhwff,posit,delim_record2 ,1,wff%marker_mpi_type,statux,ierr)
   else
     MSG_BUG('Wrong record marker length!')
   end if

 else if (option==3) then
   if (wff%nbOct_recMarker==4) then
     call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,delim_record4 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record4
   else if (wff%nbOct_recMarker==8) then
     call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,delim_record8 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record8
#if defined HAVE_FC_INT_QUAD
   else if (wff%nbOct_recMarker==16) then
     call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,delim_record16,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record16
#endif
   else if (wff%nbOct_recMarker==2) then
     call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,delim_record2 ,1,wff%marker_mpi_type,statux,ierr)
     recordmarker = delim_record2
   else
     MSG_BUG('Wrong record marker length !')
   end if

 else if (option==4) then
   if (wff%nbOct_recMarker==4) then
     delim_record4 = recordmarker
     call MPI_FILE_WRITE_AT_ALL(wff%fhwff,posit,delim_record4 ,1,wff%marker_mpi_type,statux,ierr)
   else if (wff%nbOct_recMarker==8) then
     delim_record8 = recordmarker
     call MPI_FILE_WRITE_AT_ALL(wff%fhwff,posit,delim_record8 ,1,wff%marker_mpi_type,statux,ierr)
#if defined HAVE_FC_INT_QUAD
   else if (wff%nbOct_recMarker==16) then
     delim_record16 = recordmarker
     call MPI_FILE_WRITE_AT_ALL(wff%fhwff,posit,delim_record16,1,wff%marker_mpi_type,statux,ierr)
#endif
   else if (wff%nbOct_recMarker==2) then
     delim_record2 = recordmarker
     call MPI_FILE_WRITE_AT_ALL(wff%fhwff,posit,delim_record2 ,1,wff%marker_mpi_type,statux,ierr)
   else
     MSG_BUG('Wrong record marker length!')
   end if

 else
   MSG_BUG('Wrong value for option!')
 end if

 posit = posit + recordmarker + 2*wff%nbOct_recMarker

end subroutine rwRecordMarker
#endif
!!***

!------------------------------------------------------------------------------------

!!****f* m_wffile/xnullifyOff
!! NAME
!!  xnullifyOff
!!
!! FUNCTION
!!  In case of MPI I/O, nullify the offset of a WF file
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wff=<type(wffile_type)>=structured info for reading/writing
!!
!! PARENTS
!!      m_wffile
!!
!! CHILDREN
!!      xnullifyoff
!!
!! SOURCE

subroutine xnullifyOff(wff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xnullifyOff'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff

! *************************************************************************
 
#if defined HAVE_MPI_IO
 wff%offwff    = 0
 wff%off_recs  = 0
 wff%lght_recs = 0
#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine xnullifyOff
!!***

!------------------------------------------------------------------------------------

!!****f* m_wffile/xderiveWRecEnd
!! NAME
!!  xderiveWRecEnd
!!
!! FUNCTION
!!  Writes the first and last wavefunction block marker using MPI/IO
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  me_proc= (optional argument) index of current proc
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      ioarr,outxfhist,rwwf,wffwritedatarec,wffwritenpwrec
!!
!! CHILDREN
!!      xnullifyoff
!!
!! NOTES
!!  We assume that:
!!    wff%offwff contains the position of the end of the record
!!    wff%off_recs contains the position of the beginning of the record
!!
!! SOURCE

subroutine xderiveWRecEnd(wff,ierr,me_proc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveWRecEnd'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in),optional :: me_proc
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: me
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit
!arrays
#endif

! *************************************************************************

 ierr=0

#if defined HAVE_MPI_IO
 me=-1;if (present(me_proc)) me=me_proc
 if (me==-1.or.me==0) then

   delim_record=wff%offwff-wff%off_recs-wff%nbOct_recMarker

!  Write the first word of the record
   posit=wff%off_recs
   call rwRecordMarker(2,posit,delim_record,wff,ierr)

!  Write the last word of the record
   posit=wff%offwff
   call rwRecordMarker(2,posit,delim_record,wff,ierr)

 end if

 wff%offwff = wff%offwff + wff%nbOct_recMarker
#endif

 RETURN
 ABI_UNUSED((/wff%me,me_proc/))

end subroutine xderiveWRecEnd
!!***

!------------------------------------------------------------------------------

!!****f* m_wffile/xderiveWRecInit
!! NAME
!!  xderiveWRecInit
!!
!! FUNCTION
!!  Writes the first wavefunction block marker using MPI/IO.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  me_proc= (optional argument) index of current proc
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      ioarr,outxfhist,rwwf,wffwritedatarec,wffwritenpwrec
!!
!! CHILDREN
!!      xnullifyoff
!!
!! NOTES
!!  We assume that:
!!    wff%offwff contains the position of the beginning of the record
!!
!! SOURCE

subroutine xderiveWRecInit(wff,ierr,me_proc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveWRecInit'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in),optional :: me_proc
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: me
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit
!arrays
#endif

! *************************************************************************

 ierr=0

#if defined HAVE_MPI_IO
 me=-1;if (present(me_proc)) me=me_proc
 if (me==-1.or.me==0) then

!  Write the first word of the record
   posit=wff%offwff;delim_record=0
   call rwRecordMarker(2,posit,delim_record,wff,ierr)

 end if

 wff%off_recs = wff%offwff
 wff%offwff = wff%offwff + wff%nbOct_recMarker
#endif

 RETURN
 ABI_UNUSED((/wff%me,me_proc/))

end subroutine xderiveWRecInit
!!***

!---------------------------------------------------------------------------------

!!****f* m_wffile/xderiveRRecEnd
!! NAME
!!  xderiveRRecEnd
!!
!! FUNCTION
!!  Initializes the end-of-record offset for MPI/IO.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  me_proc= (optional argument) index of current proc
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      ioarr,outxfhist,rwwf,wffreaddatarec,wffreadnpwrec
!!
!! CHILDREN
!!      xnullifyoff
!!
!! NOTES
!!  We assume that:
!!    wff%off_recs contains the position of the beginning of the record
!!
!! SOURCE

subroutine xderiveRRecEnd(wff,ierr)

 use defs_basis

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRRecEnd'
!End of the abilint section

 integer,intent(out) ::  ierr
 type(wffile_type),intent(inout) :: wff

! *************************************************************************

 ierr=0
#if defined HAVE_MPI_IO
!Define offset end of record
 wff%offwff = wff%off_recs + wff%lght_recs + 2*wff%nbOct_recMarker
#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine xderiveRRecEnd
!!***

!-------------------------------------------------------------------------------

!!****f* m_wffile/xderiveRRecInit
!! NAME
!!  xderiveRRecInit
!!
!! FUNCTION
!!  Initializes the record length for MPI/IO.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  me_proc= (optional argument) index of current proc
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      ioarr,outxfhist,rwwf,wffreaddatarec,wffreadnpwrec
!!
!! CHILDREN
!!      xnullifyoff
!!
!! NOTES
!!  We assume that:
!!    wff%offwff contains the position of the beginning of the record
!!
!! SOURCE

subroutine xderiveRRecInit(wff,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRRecInit'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit
#endif

! *************************************************************************

 ierr=0

#if defined HAVE_MPI_IO
 wff%off_recs = wff%offwff

!Read the length of the record
 posit=wff%off_recs
 call rwRecordMarker(1,posit,delim_record,wff,ierr)

 wff%lght_recs = delim_record
 wff%offwff =  wff%offwff + wff%nbOct_recMarker
#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine xderiveRRecInit
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/clsopn
!! NAME 
!! clsopn
!!
!! FUNCTION
!! Close wavefunction file (provided its access is standard F90 IO), then reopen the same.
!! Uses fortran inquire statement to reopen with same characteristics.
!!
!! INPUTS
!!  wff=number of unit to which on which file is already opened.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      nselt3,nstdy3,nstpaw3,optics_paw,optics_vloc,outkss,outwant
!!      partial_dos_fractions,rhofermi3,vtorho,vtorho3,wffile
!!
!! CHILDREN
!!      xnullifyoff
!!
!! SOURCE

subroutine clsopn(wff)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clsopn'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
!scalars
 integer :: ios,unit
 logical :: nmd,od
 character(len=11) :: fm
 character(len=500) :: message
 character(len=fnlen) :: filnam

! *************************************************************************

 if ( ANY(wff%accesswff==(/IO_MODE_FORTRAN_MASTER,IO_MODE_FORTRAN/) ))then

   unit=wff%unwff
   inquire (unit=unit,iostat=ios,opened=od,name=filnam,form=fm,named=nmd)

!  ios is a status specifier.  If an error condition exists,
!  ios is assigned a processor-dependent value > 0.
   if (ios/=0) then
     write(message, '(/,a,/,a,i8,a,i8,/,a,/,a,/,a)' ) &
&     ' clsopn : ERROR -',&
&     '  Attempt to inquire about unit=',unit,&
&     '  indicates error condition iostat=',ios,&
&     '  May be due to temporary problem with file, disks or network.',&
&     '  Action : check whether there might be some external problem,',&
&     '  then resubmit.'
     MSG_ERROR(message)

!    od is a logical variable which is set to true if the specified
!    unit is connected to a file; otherwise it is set to false.
#if !defined FC_HITACHI
   else if (.not.od) then
     write(message, '(/,a,/,a,i8,/,a,/,a,/,a,/,a)' ) &
&     ' clsopn : ERROR -',&
&     '  Tried to inquire about unit',unit,&
&     '  and found it not connected to a file.',&
&     '  May be due to temporary problem with file, disks or network.',&
&     '  Action : check whether there might be some external problem,',&
&     '  then resubmit.'
     MSG_ERROR(message)
#endif

!    nmd is a logical variable assigned the value true if the file
!    has a name; otherwise false.  A scratch file is not named.
   else if (.not.nmd) then

!    No action for the time being. Possibility to debug.

   else

!    May now close the file and then reopen it
!    (file is already opened according to above checks)

#if defined FC_HITACHI
     if (.not.od) then
       write(message, '(/,a,/,a,i8,/,a,/,a,/,a)' ) &
&       ' clsopn : WARNING - (it might be a bug on SR8k sytem)',&
&       '  Tried to inquire about unit',unit,&
&       '  and found it not connected to a file.',&
&       '  May be due to temporary problem with file, disks or network.',&
&       '  Action : disregard this error and continue the process anyway.'
       MSG_WARNING(message)
     end if
#endif
     close (unit=unit)
     open (unit=unit,file=filnam,form=fm,status='old') !VALGRIND complains filnam is just a few thousand bytes inside a block of 8300

   end if

 else if (wff%accesswff == IO_MODE_MPI) then
   call xnullifyOff(wff)
 else if (wff%accesswff == IO_MODE_ETSF) then
!  We do nothing, ETSF access already not being sequential.
 end if

end subroutine clsopn
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/wff_usef90
!! NAME 
!! wff_usef90
!!
!! FUNCTION
!!  1 if a Fortran file is going to be read by this node, 0 otherwise.
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wff_usef90(wff)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wff_usef90'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: wff_usef90
 type(wffile_type),intent(in) :: wff

! *************************************************************************

 wff_usef90=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff ==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) wff_usef90=1

end function wff_usef90
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/wff_ireadf90
!! NAME 
!! wff_ireadf90
!!
!! FUNCTION
!!  1 if a Fortran file is going to be read by this node, 0 otherwise.
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wff_ireadf90(wff)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wff_ireadf90'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: wff_ireadf90
 type(wffile_type),intent(in) :: wff

! *************************************************************************

 wff_ireadf90=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) wff_ireadf90=1

end function wff_ireadf90
!!***

!----------------------------------------------------------------------

END MODULE m_wffile
!!***
