!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_skip_int
!! NAME
!! hdr_skip_int
!!
!! FUNCTION
!! Skip wavefunction or density file header, after having rewound the file.
!! Two instances of the hdr_skip routines are defined :
!!  hdr_skip_int to which only the unit number is given
!!  hdr_skip_wfftype to which a wffil datatype is given
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG,MB,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  unit = number of unit to be read
!!
!! OUTPUT
!!  ierr = error code returned by the MPI calls
!!
!! SIDE EFFECTS
!!
!! NOTES
!! No checking performed, since hdr_skip is assumed to be used only
!! on temporary wavefunction files.
!! This initialize further reading and checking by rwwf
!!
!! PARENTS
!!      m_bse_io
!!
!! CHILDREN
!!      flush_unit,getrecordmarkerlength_wffile,mpi_bcast,mpi_file_read_at
!!      mpi_file_sync,rwrecordmarker
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine hdr_skip_int(unitfi,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_skip_int'
 use interfaces_59_io_mpi, except_this_one => hdr_skip_int
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unitfi
 integer,intent(out) :: ierr

!Local variables-------------------------------
 type(wffile_type) :: wff

! *************************************************************************

!Use default values for wff
 wff%unwff=unitfi
 wff%accesswff=IO_MODE_FORTRAN
 wff%me=0
 wff%master=0
!Then, transmit to hdr_skip_wfftype
 call hdr_skip_wfftype(wff,ierr)
end subroutine hdr_skip_int
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/hdr_skip_wfftype
!! NAME
!! hdr_skip_wfftype
!!
!! FUNCTION
!! Skip wavefunction or density file header, after having rewound the file.
!! Two instances of the hdr_skip routines are defined :
!!  hdr_skip_int to which only the unit number is given
!!  hdr_skip_wfftype to which a wffil datatype is given
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG,MB,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  unit = number of unit to be read
!!
!! OUTPUT
!!  ierr = error code returned by the MPI calls
!!
!! SIDE EFFECTS
!!
!! NOTES
!! No checking performed, since hdr_skip is assumed to be used only
!! on temporary wavefunction files.
!! This initialize further reading and checking by rwwf
!!
!! PARENTS
!!      hdr_io,hdr_io_netcdf,hdr_skip
!!
!! CHILDREN
!!      flush_unit,getrecordmarkerlength_wffile,mpi_bcast,mpi_file_read_at
!!      mpi_file_sync,rwrecordmarker
!!
!! SOURCE
 subroutine hdr_skip_wfftype(wff,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile

 use m_io_tools, only : flush_unit

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_skip_wfftype'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif


!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer, intent(out) :: ierr

!Local variables-------------------------------
 integer :: headform,mu,npsp,unit,usepaw
 integer :: integers(17)
 character(len=6) :: codvsn
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit,positloc
 integer :: statux(MPI_STATUS_SIZE)
#endif

!*************************************************************************

 unit=wff%unwff
 ierr=0

 if( wff%accesswff==IO_MODE_FORTRAN .or. (wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me) ) then

   rewind (unit)

!  Pick off headform from WF file
   read(unit) codvsn,headform                  ! XG040806 This does not work, but I do not understand why ?!
!  read(unit) integers(1),headform             ! This works ...
!  ! MT012408 Because codvsn is char*6 and not an integer !
   if(headform==1   .or. headform==2   .or. &
&   headform==51  .or. headform==52  .or. &
&   headform==101 .or. headform==102       ) headform=22

   if (headform<44) then
     read (unit) integers(1:13),npsp
   else
     read (unit) integers(1:13),npsp,integers(15:17),usepaw
   end if

!  Skip rest of header records
   do mu=1,2+npsp
     read (unit)
   end do
   if ((headform>=44).and.(usepaw==1)) then
     read (unit)
     read (unit)
   end if

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then

   headform=wff%headform
   if(headform==1   .or. headform==2   .or. &
&   headform==51  .or. headform==52  .or. &
&   headform==101 .or. headform==102       ) headform=22

!  Causes all previous writes to be transferred to the storage device
   call flush_unit(wff%unwff)
   call MPI_FILE_SYNC(wff%fhwff,ierr)

!  Check FORTRAN record marker length (only at first call)
   if (wff%nbOct_recMarker<=0) then
     call getRecordMarkerLength_wffile(wff)
   end if

   if (wff%master==wff%me) then

!    Reading the first record of the file -------------------------------------
!    read (unitfi)   codvsn,headform,..............
     posit = 0
     call rwRecordMarker(1,posit,delim_record,wff,ierr)

!    Reading the second record of the file ------------------------------------
!    read(unitfi) bantot, hdr%date, hdr%intxc.................
!    Pick off npsp and usepaw from WF file
     positloc  = posit + wff%nbOct_recMarker + wff%nbOct_int*13
     call MPI_FILE_READ_AT(wff%fhwff,positloc,npsp,1,MPI_INTEGER,statux,ierr)
!    call MPI_FILE_READ_AT_ALL(wff%fhwff,positloc,npsp,1,MPI_INTEGER,statux,ierr)
     if (headform >= 44) then
       positloc = positloc +  wff%nbOct_int*4
       call MPI_FILE_READ_AT(wff%fhwff,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
!      call MPI_FILE_READ_AT_ALL(wff%fhwff,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
     end if
     call rwRecordMarker(1,posit,delim_record,wff,ierr)

!    Reading the rest of the file ---------------------------------------------
     do mu=1,2+npsp
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
     end do
     if ((headform>=44).and.(usepaw==1)) then
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
     end if

     wff%offwff=posit

   end if

   if (wff%spaceComm/=MPI_COMM_SELF) then
     call MPI_BCAST(wff%offwff,1,wff%offset_mpi_type,wff%master,wff%spaceComm,ierr)
   end if
#endif

 end if

end subroutine hdr_skip_wfftype
!!***
