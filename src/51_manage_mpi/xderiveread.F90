!{\src2tex{textfont=tt}}
!!****f* ABINIT/xderiveRead
!! NAME
!!  xderiveRead
!!
!! FUNCTION
!!  Generic subroutines to read wf files.
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (MB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
!!  We use several procedures with the same generic name
!!    xderiveRead  contains
!!               xderiveRead_int        : read integer value
!!               xderiveRead_int1d      : read integer array 1d
!!               xderiveRead_int2d      : read integer array 2d
!!               xderiveRead_int2d_displ: read integer array 2d non-contiguous
!!               xderiveRead_dp         : read double precision value
!!               xderiveRead_dp1d       : read double precision array 1d
!!               xderiveRead_dp2d       : read double precision array 2d
!!               xderiveRead_dp2d_displ : read double precision array 2d non-contiguous
!!               xderiveRead_char       : read character string
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
!!***

!!****f* ABINIT/xderiveRead_int
!! NAME
!!  xderiveRead_int
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: integer scalar.
!!
!! INPUTS
!! (none)
!!
!! OUTPUT
!!  xval= data buffer
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xderiveRead_int(wff,xval,ierr)

 use m_profiling

 use defs_datatypes
 use m_wffile
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRead_int'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 xval=0; ierr=0

#if defined HAVE_MPI_IO
 call MPI_FILE_READ_AT(wff%fhwff,wff%offwff,xval,1,MPI_INTEGER,statux,ierr)
 wff%offwff = wff%offwff + wff%nbOct_int
#endif

 RETURN
 ABI_UNUSED(wff%me)

end subroutine xderiveRead_int
!!***


!!****f* ABINIT/xderiveRead_int1d
!! NAME
!!  xderiveRead_int1d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  xval= data buffer array
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xderiveRead_int1d(wff,xval,n1,spaceComm,ierr)

 use m_profiling

 use defs_datatypes
 use m_wffile
 use m_xmpi

#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRead_int1d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval(:)
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit,nboct,dispoct,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 xval(:)=0 ; ierr=0 ! Initialization, for the compiler
 if(.false.)write(std_out,*)wff%me,n1,spaceComm

#if defined HAVE_MPI_IO
 nboct = wff%nbOct_int * n1
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then
!  Compute offset for local part
!  dispoct = sum (nboct, rank=0..me)
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
     posit = posit+dispoct-nboct
   end if
   call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1,MPI_INTEGER,statux,ierr)

!  get the total number of bits wrote by processors
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
   else
     totoct=nboct
   end if
 else
   ierr = 1
   nboct =0
   totoct = 0
 end if

!new offset
 wff%offwff = wff%offwff + totoct
#endif

end subroutine xderiveRead_int1d
!!***


!!****f* ABINIT/xderiveRead_int2d
!! NAME
!!  xderiveRead_int2d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  xval= data buffer array
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xderiveRead_int2d(wff,xval,n1,n2,spaceComm,ierr)

 use m_profiling

 use defs_datatypes
 use m_wffile
 use m_xmpi

#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRead_int2d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval(:,:)
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,dispoct,nboct,posit,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 xval(:,:)=0 ; ierr=0 ! Initialization, for the compiler
 if(.false.)write(std_out,*)wff%me,n1,n2,spaceComm

#if defined HAVE_MPI_IO
 nboct = wff%nbOct_int * n1 * n2
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then
!  Compute offset for local part
!  dispoct = sum (nboct, rank=0..me)
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
     posit = posit + dispoct - nboct
   end if
   call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)

!  get the total number of bits wrote by processors
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
   else
     totoct=nboct
   end if
 else
   ierr = 1
   nboct =0
   totoct = 0
 end if

!new offset
 wff%offwff=wff%offwff + totoct
#endif

end subroutine xderiveRead_int2d
!!***


!!****f* ABINIT/xderiveRead_dp
!! NAME
!!  xderiveRead_dp
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision scalar.
!!
!! INPUTS
!! (none)
!!
!! OUTPUT
!!  xval= data buffer
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xderiveRead_dp(wff,xval,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRead_dp'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 xval=zero ; ierr=0
 if(.false.)write(std_out,*)wff%me
#if defined HAVE_MPI_IO
 call MPI_FILE_READ_AT(wff%fhwff,wff%offwff,xval,1,MPI_DOUBLE_PRECISION,statux,ierr)
 wff%offwff = wff%offwff + wff%nbOct_dp
#endif

end subroutine xderiveRead_dp
!!***


!!****f* ABINIT/xderiveRead_dp1d
!! NAME
!!  xderiveRead_dp1d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: one-dimensional double precision arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!  xval= data buffer array
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

 subroutine xderiveRead_dp1d(wff,xval,n1,spaceComm,ierr)

 use m_profiling
 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile

#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRead_dp1d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval(:)

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit,nboct,dispoct,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif

!*********************************************************************

 xval(:)=zero ; ierr=0 ! Initialization, for the compiler
 if(.false.)write(std_out,*)wff%me,n1,spaceComm

#if defined HAVE_MPI_IO
 nboct = wff%nbOct_dp * n1
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then
!  Compute offset for local part
!  dispoct = sum (nboct, rank=0..me)
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
     posit = posit + dispoct - nboct
   end if

   call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1,MPI_DOUBLE_PRECISION,statux,ierr)

!  get the total number of bits wrote by processors
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
   else
     totoct=nboct
   end if
 else
   ierr = 1
   nboct =0
   totoct = 0
 end if

!new offset
 wff%offwff=wff%offwff + totoct
#endif

end subroutine xderiveRead_dp1d
!!***


!!****f* ABINIT/xderiveRead_dp2d
!! NAME
!!  xderiveRead_dp2d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!  xval= data buffer array
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xderiveRead_dp2d(wff,xval,n1,n2,spaceComm,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile

#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRead_dp2d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval(:,:)

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,dispoct,nboct,posit,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 xval(:,:)=zero ; ierr=0 ! Initialization, for the compiler
 if(.false.)write(std_out,*)wff%me,n1,n2,spaceComm

#if defined HAVE_MPI_IO
 nboct = wff%nbOct_dp * n1 *n2
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then
!  Compute offset for local part
!  dispoct = sum (nboct, rank=0..me)
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
     posit = posit + dispoct - nboct
   end if
   call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1*n2,MPI_DOUBLE_PRECISION,statux,ierr)

!  get the total number of bits wrote by processors
   if (spaceComm/=MPI_COMM_SELF) then
     call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
   else
     totoct=nboct
   end if
 else
   ierr = 1
   nboct =0
   totoct = 0
 end if

!new offset
 wff%offwff=wff%offwff + totoct
#endif

end subroutine xderiveRead_dp2d
!!***


!!****f* ABINIT/xderiveRead_int2d_displ
!! NAME
!!  xderiveRead_int2d_displ
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  displace= number of elements for the offset
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!  xval= data buffer array
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xderiveRead_int2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRead_int2d_displ'
!End of the abilint section

  implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(out):: xval(:,:)
 integer,intent(in):: displace(:)

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: filetype,i1,i2,ipos,nb,nbval,totsize,wfftempo
!arrays
 integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: buf_val(:),depl(:),depl1(:),depl_val(:),length1(:),type1(:),val(:)
#endif

! *********************************************************************

 xval(:,:)=0 ; ierr=0
 if(.false.)write(std_out,*)wff%me,n1,n2,spaceComm,displace

#if defined HAVE_MPI_IO
 nb=n1*n2
 call xsum_mpi(nb,totsize,spaceComm,ierr)
 ABI_ALLOCATE(depl_val,(0:totsize-1))
 ABI_ALLOCATE(depl,(nb))
 ABI_ALLOCATE(buf_val,(0:totsize-1))
 ABI_ALLOCATE(val,(nb))

!Map displacements
 depl_val(0:totsize-1)=-1
 do i2=1,n2
   do i1=1,n1
     ipos=(displace(i2)-1)*n1 + i1-1
     depl_val(ipos)=ipos
   end do
 end do
!To save time, the location describe by array map must be in increasing order
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     depl(nbval)=depl_val(i1)
   end if
 end do

!Build MPI datatype for view
 ABI_ALLOCATE(length1,(nbval+2))
 ABI_ALLOCATE(depl1,(nbval+2))
 ABI_ALLOCATE(type1,(nbval+2))
 length1(1)=1;depl1(1)=0;type1(1)=MPI_LB
 do i1=2,nbval+1
   length1(i1) = 1
   depl1(i1)= depl(i1-1)*wff%nbOct_int
   type1(i1)= MPI_INTEGER
 end do
 length1(nbval+2)=1;depl1(nbval+2)=totsize*wff%nbOct_int;type1(nbval+2)=MPI_UB
 call MPI_TYPE_STRUCT(nbval+2,length1,depl1,type1,filetype,ierr)
 call MPI_TYPE_COMMIT(filetype,ierr)
 ABI_DEALLOCATE(length1)
 ABI_DEALLOCATE(depl1)
 ABI_DEALLOCATE(type1)

!Write data
 call MPI_FILE_OPEN(spaceComm,wff%fname,MPI_MODE_RDWR,MPI_INFO_NULL,wfftempo,ierr)
 call MPI_FILE_SET_VIEW(wfftempo,wff%offwff,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
 call MPI_FILE_READ_ALL(wfftempo,val,nbval,MPI_INTEGER,statux,ierr)
 call MPI_FILE_CLOSE(wfftempo,ierr)

!Retrieve xval
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     buf_val(i1)=val(nbval)
   end if
 end do
 do i2=1,n2
   do i1=1,n1
     ipos=(displace(i2)-1)*n1 + i1-1
     xval(i1,i2)=buf_val(ipos)
   end do
 end do

!Update offset
 wff%offwff = wff%offwff + totsize*wff%nbOct_int

!Free memory
 call MPI_TYPE_FREE(filetype,ierr)
 ABI_DEALLOCATE(depl)
 ABI_DEALLOCATE(depl_val)
 ABI_DEALLOCATE(buf_val)
 ABI_DEALLOCATE(val)
#endif

end subroutine xderiveRead_int2d_displ
!!***

! ==========================================================================

!!****f* ABINIT/xderiveRead_dp2d_displ
!! NAME
!!  xderiveRead_dp2d_displ
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  displace= number of elements for the offset
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!  xval= data buffer array
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xderiveRead_dp2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveRead_dp2d_displ'
!End of the abilint section

  implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(out):: xval(:,:)
 integer,intent(in):: displace(:)

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: filetype,i1,i2,ipos,nb,nbval,totsize,wfftempo
!arrays
 integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: depl(:),depl1(:),depl_val(:),length1(:),type1(:)
 real(dp), allocatable :: buf_val(:),val(:)
#endif

! *********************************************************************

 xval(:,:)=zero ; ierr=0
 if(.false.)write(std_out,*)wff%me,n1,n2,displace,spaceComm

#if defined HAVE_MPI_IO
 nb=n1*n2
 call xsum_mpi(nb,totsize,spaceComm,ierr)
 ABI_ALLOCATE(depl_val,(0:totsize-1))
 ABI_ALLOCATE(depl,(nb))
 ABI_ALLOCATE(buf_val,(0:totsize-1))
 ABI_ALLOCATE(val,(nb))

!Map displacements
 depl_val(0:totsize-1)=-1
 do i2=1,n2
   do i1=1,n1
     ipos=(displace(i2)-1)*n1 + i1-1
     depl_val(ipos)=ipos
   end do
 end do
!To save time, the location describe by array map must be in increasing order
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     depl(nbval)=depl_val(i1)
   end if
 end do

!Build MPI datatype for view
 ABI_ALLOCATE(length1,(nbval+2))
 ABI_ALLOCATE(depl1,(nbval+2))
 ABI_ALLOCATE(type1,(nbval+2))
 length1(1)=1;depl1(1)=0;type1(1)=MPI_LB
 do i1=2,nbval+1
   length1(i1) = 1
   depl1(i1)= depl(i1-1)*wff%nbOct_dp
   type1(i1)= MPI_DOUBLE_PRECISION
 end do
 length1(nbval+2)=1;depl1(nbval+2)=totsize*wff%nbOct_dp;type1(nbval+2)=MPI_UB
 call MPI_TYPE_STRUCT(nbval+2,length1,depl1,type1,filetype,ierr)
 call MPI_TYPE_COMMIT(filetype,ierr)
 ABI_DEALLOCATE(length1)
 ABI_DEALLOCATE(depl1)
 ABI_DEALLOCATE(type1)

!Write data
 call MPI_FILE_OPEN(spaceComm,wff%fname,MPI_MODE_RDWR,MPI_INFO_NULL,wfftempo,ierr)
 call MPI_FILE_SET_VIEW(wfftempo,wff%offwff,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
 call MPI_FILE_READ_ALL(wfftempo,val,nbval,MPI_DOUBLE_PRECISION,statux,ierr)
 call MPI_FILE_CLOSE(wfftempo,ierr)

!Retrieve xval
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     buf_val(i1)=val(nbval)
   end if
 end do
 do i2=1,n2
   do i1=1,n1
     ipos=(displace(i2)-1)*n1 + i1-1
     xval(i1,i2)=buf_val(ipos)
   end do
 end do

!Update offset
 wff%offwff = wff%offwff + totsize*wff%nbOct_dp

!Free memory
 call MPI_TYPE_FREE(filetype,ierr)
 ABI_DEALLOCATE(depl)
 ABI_DEALLOCATE(depl_val)
 ABI_DEALLOCATE(buf_val)
 ABI_DEALLOCATE(val)
#endif

end subroutine xderiveRead_dp2d_displ
!!***


!!****f* ABINIT/xderiveReadVal_char
!! NAME
!!  xderiveReadVal_char
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: character string.
!!
!! INPUTS
!!  n= number of elements in the array
!!
!! OUTPUT
!!  xval= data buffer array
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xderiveReadVal_char(wff,xval,n,ierr)

 use m_profiling

 use defs_datatypes
 use m_wffile
 use m_xmpi
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveReadVal_char'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n
 integer,intent(out) :: ierr
 character(len=*),intent(out) :: xval

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 xval=' ' ; ierr=0
 if(.false.)write(std_out,*)wff%me,n

#if defined HAVE_MPI_IO
 call MPI_FILE_READ_AT(wff%fhwff,wff%offwff,xval,n,MPI_CHARACTER,statux,ierr)
 wff%offwff = wff%offwff + wff%nbOct_ch * n
#endif

end subroutine xderiveReadVal_char
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xmpi_read_int2d
!! NAME
!!  xmpi_read_int2d
!!
!! FUNCTION
!!  Generic routine to read arrays with MPI I/O.
!!  Target: integer two-dimensional arrays.
!!
!! INPUTS
!!  at_option= 
!!    xmpio_at     ==> Local reading.
!!    xmpio_at_all ==> Collective reading.
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  xval= data buffer array
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_file_read_at_all
!!
!! SOURCE

subroutine xmpi_read_int2d(wff,xval,spaceComm,at_option,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile
 use m_errors      

 use m_fstrings,  only : toupper

#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_read_int2d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spaceComm,at_option
 integer,intent(out) :: ierr
 type(wffile_type),intent(inout) :: wff
!array
 integer,intent(out) :: xval(:,:)

!Local variables-------------------------------
 integer :: n1,n2
#ifdef HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,nboct,posit,totoct
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 ierr=0 
 n1 = SIZE(xval,DIM=1)
 n2 = SIZE(xval,DIM=2)

#ifdef HAVE_MPI_IO
 nboct = wff%nbOct_int * n1 *n2
 posit = wff%offwff
 delim_record = posit - wff%off_recs + wff%lght_recs - wff%nbOct_recMarker

 if (delim_record >= nboct) then

   select case (at_option)
     case (xmpio_at) 
       call MPI_FILE_READ_AT(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)

     case (xmpio_at_all)
       call MPI_FILE_READ_AT_ALL(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)

       case default
       write(msg,('(a,i0)'))" Wrong value for at_option: ",at_option
       MSG_ERROR(msg)
   end select

   totoct=nboct
 else
   write(msg,('(a,2(i0,1x))'))" delim_record < nboct: ",delim_record,nboct
   MSG_WARNING(msg)
   ierr=MPI_ERR_UNKNOWN
   totoct=0
 end if
!
!Increment the offset.
 wff%offwff=wff%offwff + totoct
#endif

 RETURN
 ABI_UNUSED(xval(1,1))
 ABI_UNUSED((/wff%me,spaceComm,at_option/))

end subroutine xmpi_read_int2d
!!***
