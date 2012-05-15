!{\src2tex{textfont=tt}}
!!****f* ABINIT/xderiveWrite
!! NAME
!!  xderiveWrite
!!
!! FUNCTION
!!  Generic subroutines to read/write wf files.
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (MB,MD,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
!!  We use several procedures with the same generic name
!!    xderiveWrite  contains
!!               xderiveWrite_int         : write integer value
!!               xderiveWrite_int1d       : write integer array 1d
!!               xderiveWrite_int2d       : write integer array 2d
!!               xderiveWrite_int2d_displ : write integer array 2d non contiguous
!!               xderiveWrite_dp          : write double precision value
!!               xderiveWrite_dp1d        : write double precision array 1d
!!               xderiveWrite_dp2d        : write double precision array 2d
!!               xderiveWrite_dp2d_displ  : write double precision array 2d non contiguous
!!               xderiveWrite_char        : write character string
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!!***

!!****f* ABINIT/xderiveWrite_int
!! NAME
!!  xderiveWrite_int
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: integer scalar.
!!
!! INPUTS
!!  xval= data buffer
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_int(wff,xval,ierr)

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
#define ABI_FUNC 'xderiveWrite_int'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(out) :: ierr
 integer,intent(in):: xval
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval
#if defined HAVE_MPI_IO
 call MPI_FILE_WRITE_AT(wff%fhwff,wff%offwff,xval,1,MPI_INTEGER,statux,ierr)
 wff%offwff = wff%offwff+wff%nbOct_int
#endif

end subroutine xderiveWrite_int
!!***


!!****f* ABINIT/xderiveWrite_int1d
!! NAME
!!  xderiveWrite_int1d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_int1d(wff,xval,n1,spaceComm,ierr)

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
#define ABI_FUNC 'xderiveWrite_int1d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: xval(:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: dispoct,nboct,posit,totoct
 integer :: statux(MPI_STATUS_SIZE)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,n1,spaceComm,xval
#if defined HAVE_MPI_IO
 nboct = n1*wff%nbOct_int
 posit = wff%offwff

!dispoct = sum (nboct, rank=0..me)
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
   posit = posit+dispoct-nboct
 end if
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1,MPI_INTEGER,statux,ierr)
!gather the bigest offset
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
 else
   totoct=nboct
 end if
 wff%offwff = wff%offwff+totoct

!Disable old code
#endif

end subroutine xderiveWrite_int1d
!!***


!!****f* ABINIT/xderiveWrite_int2d
!! NAME
!!  xderiveWrite_int2d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_int2d(wff,xval,n1,n2,spaceComm,ierr)

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
#define ABI_FUNC 'xderiveWrite_int2d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: dispoct,nboct,posit,totoct
 integer  :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,n1,n2,spaceComm,xval
#if defined HAVE_MPI_IO
 nboct = n1*n2*wff%nbOct_int
 posit = wff%offwff

!dispoct = sum(nboct, rank=0..me)
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
   posit = posit + dispoct-nboct
 end if
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)
!gather the biggest offset
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
 else
   totoct=nboct
 end if
 wff%offwff = wff%offwff + totoct
#endif

end subroutine xderiveWrite_int2d
!!***


!!****f* ABINIT/xderiveWrite_dp
!! NAME
!!  xderiveWrite_dp
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision scalar.
!!
!! INPUTS
!!  xval= data buffer
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_dp(wff,xval,ierr)

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
#define ABI_FUNC 'xderiveWrite_dp'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval
#if defined HAVE_MPI_IO
 call MPI_FILE_WRITE_AT(wff%fhwff,wff%offwff,xval,1,MPI_DOUBLE_PRECISION,statux,ierr)
 wff%offwff = wff%offwff+wff%nbOct_dp
#endif

end subroutine xderiveWrite_dp
!!***


!!****f* ABINIT/xderiveWrite_dp1d
!! NAME
!!  xderiveWrite_dp1d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: one-dimensional double precision arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_dp1d(wff,xval,n1,spaceComm,ierr)

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
#define ABI_FUNC 'xderiveWrite_dp1d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: nboct,dispoct,totoct,posit
 integer  :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,n1,spaceComm,xval
#if defined HAVE_MPI_IO
 nboct = n1*wff%nbOct_dp
 posit = wff%offwff
!dispoct = sum (nboct, rank = 0..me)
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
   posit = posit + dispoct - nboct
 end if
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1,MPI_DOUBLE_PRECISION,statux,ierr)
!Gather the biggest offset
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
 else
   totoct=nboct
 end if
 wff%offwff = wff%offwff + totoct
#endif

end subroutine xderiveWrite_dp1d
!!***

!!****f* ABINIT/xderiveWrite_dp2d
!! NAME
!!  xderiveWrite_dp2d
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_dp2d(wff,xval,n1,n2,spaceComm,ierr)

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
#define ABI_FUNC 'xderiveWrite_dp2d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: nboct,dispoct,totoct,posit
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval,n1,n2,spaceComm

#if defined HAVE_MPI_IO
 nboct = n1*n2*wff%nbOct_dp
 posit = wff%offwff
!dispoct = sum(nboct, rank=0..me)
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_SCAN(nboct,dispoct,1,wff%offset_mpi_type,MPI_SUM,spaceComm,ierr)
   posit = posit+dispoct-nboct
 end if
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1*n2,MPI_DOUBLE_PRECISION,statux,ierr)
 posit = posit+nboct
!gather the biggest offset
 if (spaceComm/=MPI_COMM_SELF) then
   call MPI_ALLREDUCE(dispoct,totoct,1,wff%offset_mpi_type,MPI_MAX,spaceComm,ierr)
 else
   totoct=nboct
 end if
 wff%offwff = wff%offwff+totoct
#endif

end subroutine xderiveWrite_dp2d
!!***


!!****f* ABINIT/xderiveWrite_int2d_displ
!! NAME
!!  xderiveWrite_int2d_displ
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!  displace= number of elements for the offset
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_int2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)

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
#define ABI_FUNC 'xderiveWrite_int2d_displ'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: displace(:),xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: filetype,i1,i2,ipos,nb,nbval,totsize,wfftempo
!arrays
 integer :: statux(MPI_STATUS_SIZE)
 integer, allocatable :: buf_val(:),depl(:),depl1(:),depl_val(:)
 integer, allocatable :: length1(:),type1(:),val(:)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval,n1,n2,spaceComm,displace

#if defined HAVE_MPI_IO
 nb = n1*n2
 call xsum_mpi(nb,totsize,spaceComm,ierr)
 ABI_ALLOCATE(depl_val,(0:totsize-1))
 ABI_ALLOCATE(depl,(nb))
 ABI_ALLOCATE(buf_val,(0:totsize-1))
 ABI_ALLOCATE(val,(nb))

!Map displacements
!Put xval in a buffer at its position
 depl_val(0:totsize-1)=-1
 do i2=1,n2
   do i1=1,n1
!    ipos location of xval(i1,i2) in the array associated with record to be written
     ipos=(displace(i2)-1)*n1 + i1-1
     buf_val(ipos) = xval(i1,i2)
     depl_val(ipos) = ipos
   end do
 end do
!To save time, the location describe by array map must be in increasing order
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     val(nbval)=buf_val(i1)
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
 call MPI_FILE_WRITE_ALL(wfftempo,val,nbval,MPI_INTEGER,statux,ierr)
 call MPI_FILE_CLOSE(wfftempo,ierr)

!Update offset
 wff%offwff = wff%offwff + totsize*wff%nbOct_int

!Free memory
 call MPI_TYPE_FREE(filetype,ierr)
 ABI_DEALLOCATE(depl)
 ABI_DEALLOCATE(depl_val)
 ABI_DEALLOCATE(buf_val)
 ABI_DEALLOCATE(val)
#endif

end subroutine xderiveWrite_int2d_displ
!!***


!!****f* ABINIT/xderiveWrite_dp2d_displ
!! NAME
!!  xderiveWrite_dp2d_displ
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: two-dimensional double precision arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  n2= second dimension of the array
!!  spaceComm= MPI communicator
!!  xval= data buffer array
!!  displace= number of elements for the offset
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_dp2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)

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
#define ABI_FUNC 'xderiveWrite_dp2d_displ'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: displace(:)
 real(dp),intent(in) :: xval(:,:)
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: filetype,i1,i2,ipos,nb,nbval,totsize,wfftempo
!arrays
 integer :: statux(MPI_STATUS_SIZE)
 integer, allocatable :: depl(:),depl1(:),depl_val(:)
 integer, allocatable :: length1(:),type1(:)
 real(dp),allocatable :: buf_val(:),val(:)
#endif

! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval,n1,n2,spaceComm,displace

#if defined HAVE_MPI_IO
 nb = n1*n2
 call xsum_mpi(nb,totsize,spaceComm,ierr)
 ABI_ALLOCATE(depl_val,(0:totsize-1))
 ABI_ALLOCATE(depl,(nb))
 ABI_ALLOCATE(buf_val,(0:totsize-1))
 ABI_ALLOCATE(val,(nb))

!Map displacements
!Put xval in a buffer at its position
 depl_val(0:totsize-1)=-1
 do i2=1,n2
   do i1=1,n1
!    ipos location of xval(i1,i2) in the array associated with record to be written
     ipos=(displace(i2)-1)*n1 + i1-1
     buf_val(ipos) = xval(i1,i2)
     depl_val(ipos) = ipos
   end do
 end do
!To save time, the location describe by array map must be in increasing order
 nbval=0
 do i1=0,totsize-1
   if (depl_val(i1)/=-1) then
     nbval=nbval+1
     val(nbval)=buf_val(i1)
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
 call MPI_FILE_WRITE_ALL(wfftempo,val,nbval,MPI_DOUBLE_PRECISION,statux,ierr)
 call MPI_FILE_CLOSE(wfftempo,ierr)

 wff%offwff = wff%offwff + totsize*wff%nbOct_dp

!Free memory
 call MPI_TYPE_FREE(filetype,ierr)
 ABI_DEALLOCATE(depl)
 ABI_DEALLOCATE(depl_val)
 ABI_DEALLOCATE(buf_val)
 ABI_DEALLOCATE(val)
#endif

end subroutine xderiveWrite_dp2d_displ
!!***


!!****f* ABINIT/xderiveWrite_char
!! NAME
!!  xderiveWrite_char
!!
!! FUNCTION
!!  Generic routine to read/write wf files with MPI I/O.
!!  Target: character string.
!!
!! INPUTS
!!  xval= data buffer array
!!  n= number of elements in the string
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  wff= structured info for reading/writing the wavefunctions
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_write_at
!!
!! SOURCE

subroutine xderiveWrite_char(wff,xval,n,ierr)

 use m_profiling

 use defs_datatypes
 use m_xmpi
 use m_wffile

#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xderiveWrite_char'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: xval

!Local variables-------------------------------
#if defined HAVE_MPI_IO
 integer :: statux(MPI_STATUS_SIZE)
#endif
! *********************************************************************

 ierr=0
 if(.false.)write(std_out,*)wff%me,xval,n

#if defined HAVE_MPI_IO
 call MPI_FILE_WRITE_AT(wff%fhwff,wff%offwff,xval,n,MPI_CHARACTER,statux,ierr)
 wff%offwff = wff%offwff + wff%nbOct_ch * n
#endif

end subroutine xderiveWrite_char
!!***
