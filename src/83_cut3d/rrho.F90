!{\src2tex{textfont=tt}}
!!****f* ABINIT/rrho
!! NAME
!! rrho
!!
!! FUNCTION
!! Reads in the charge in mkdens3D format
!! The file was opened in the calling program, unit number 19.
!! The header was already read in the case of the unformatted file
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! densfileformat=integer flag :
!! 0=formatted (ASCII) or 1=unformatted (BINARY) or 2=unformatted NETCDF-ETSF file 
!! nr1=grid_full size along x
!! nr2=grid_full size along y
!! nr3=grid_full size along z
!! nspden=number of spin polartized densities (1 for non-spin polarized, 2 for spin-polarized)
!! wff=structure type containing the density information
!!
!! OUTPUT
!! grid_full(nr1,nr2,nr3)=grid_full matrix
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,etsf_io_main_get,leave_new,wffclose,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rrho(densfileformat,grid_full,nr1,nr2,nr3,nspden,wff)

 use m_profiling

 use defs_basis
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
 use defs_datatypes
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rrho'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: densfileformat,nr1,nr2,nr3,nspden
 type(wffile_type),intent(inout) :: wff
!arrays
 real(dp),intent(out),target :: grid_full(nr1,nr2,nr3,nspden)

!Local variables--------------------------------------------------------
#if defined HAVE_TRIO_ETSF_IO
 type(etsf_main), target :: main_folder
 logical :: lstat
 character(len=500) :: message
 character(len = etsf_io_low_error_len) :: errmess
 type(etsf_io_low_error) :: error
#endif
!scalars
 integer :: ierr,ir1,ir2,ir3,ispden

! *************************************************************************

 select case (densfileformat)

!  Formatted (only one spin component is allowed)
   case (0)
     do ir3=1,nr3
       do ir2=1,nr2
         do ir1=1,nr1
           read(unit=wff%unwff,fmt=*) grid_full(ir1,ir2,ir3,1)
         end do
       end do
     end do

!    Unformatted, on one record
   case (1)
     do ispden=1,nspden
       read(unit=wff%unwff) grid_full(1:nr1,1:nr2,1:nr3,ispden)
     end do
!    ETSF case
   case (2)
#if defined HAVE_TRIO_ETSF_IO
     main_folder%density%data4D => grid_full
     call etsf_io_main_get(wff%unwff, main_folder, lstat, error)
     if (.not. lstat) then
!      We handle the error
       call etsf_io_low_error_to_str(errmess, error)
       write(message, "(A,A,A,A)") ch10, " rrho: ERROR -", ch10, &
&       errmess(1:min(475, len(errmess)))
       call wrtout(std_out, message, 'COLL')
       call leave_new('COLL')
     end if
#else
     write(std_out,*)&
&     'Error, ETSF_IO support is not compiled. Reconfigure with --enable-etsf-io.'
     stop     
#endif

     case default
     write(std_out,*)&
&     'Error, value for 3D function file format is invalid: ',densfileformat
     stop

 end select

 call wffclose(wff, ierr)

end subroutine rrho
!!***
