!{\src2tex{textfont=tt}}
!!****f* ABINIT/write_var_netcdf
!!
!! NAME
!! write_var_netcdf
!!
!! FUNCTION
!! Write the history into a netcdf dataset
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! arr_int
!! arr_real
!! marr
!! narr
!! typevar
!! varname
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      prttagm,prttagm_images
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_read_dim,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine write_var_netcdf(arr_int,arr_real,marr,narr,ncid,typevar,varname)

 use m_profiling

! define dp,sixth,third,etc...
 use defs_basis
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_var_netcdf'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: narr,marr,ncid
 character(len=*),intent(in) :: varname
 character(len=3),intent(in) :: typevar
!arrays
 integer,intent(in) :: arr_int(marr)
 real(dp),intent(in) :: arr_real(marr)


!Local variables-------------------------------
!scalars
 integer :: var_id,var_type,vardim_id,dimvalue,ncerr
! integer :: ncerr
 logical :: lstat
!arrays

#if defined HAVE_TRIO_ETSF_IO
type(etsf_io_low_error) :: error_data
#endif

!***************************************************************************
!Beginning of executable session
!***************************************************************************

#if defined HAVE_TRIO_ETSF_IO

 if (ncid>0) then

!  write(std_out,*) 'write_var_netcdf 01',lstat
!  ###########################################################
!  ### 01. Put the file in definition mode

   call etsf_io_low_set_define_mode(ncid, lstat, error_data)

!  ncerr = nf90_redef(ncid)
!  if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," nf90_redef")
   
!  write(std_out,*) 'write_var_netcdf 02',lstat
!  ###########################################################
!  ### 02. Define the dimensions

   if (narr==1)then
     
     call etsf_io_low_read_dim(ncid, "one", dimvalue, lstat, vardim_id, error_data)
!    ncerr = nf90_inq_varid(ncid, "one", vardim_id)
!    if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," nf90_inq_varid")
   else

     call etsf_io_low_write_dim(ncid, trim(varname), narr, lstat, vardim_id, error_data)
!    ncerr = nf90_def_dim(ncid,trim(varname),narr,vardim_id)
!    if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," nf90_def_dim")
   end if

!  write(std_out,*) 'write_var_netcdf 03',lstat
!  ###########################################################
!  ### 03. Define the variables

!  if (typevar=='INT') then
!  var_type=NF90_INT
!  else if (typevar=='DPR') then
!  var_type=NF90_DOUBLE
!  end if

   if (typevar=='INT') then
     var_type=etsf_io_low_integer
   else if (typevar=='DPR') then
     var_type=etsf_io_low_double
   end if

   if (narr==1)then
     call etsf_io_low_def_var(ncid, trim(varname), var_type,&
&     (/ "one" /)         , lstat, var_id, error_data)
   else
     call etsf_io_low_def_var(ncid, trim(varname), var_type,&
&     (/ trim(varname) /)  , lstat, var_id, error_data)
   end if

!  ncerr = nf90_def_var(ncid, trim(varname), var_type, vardim_id, var_id)
!  if (ncerr /= NF90_NOERR)  call handle_ncerr(ncerr," define variable "//trim(varname))

!  write(std_out,*) 'write_var_netcdf 04',lstat
!  ###########################################################
!  ### 04. Put the file in data mode

   call etsf_io_low_set_write_mode(ncid, lstat, error_data)

!  ncerr = nf90_enddef(ncid)
!  if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," nf90_enddef")

!  write(std_out,*) 'write_var_netcdf 05',lstat
!  #####################################################################
!  ### Write variables into the dataset

   if (typevar=='INT') then
     
!    call etsf_io_low_write_var(ncid, trim(varname), arr_int(1:narr), lstat,&
!    & start = (/ 1 /) , count = (/ narr /),error_data=error_data)
!    call etsf_io_low_error_handle(error_data)

!    write(std_out,*)  trim(varname), lstat,'NARR',narr,'MARR',marr,arr_int(:)
     ncerr = nf90_put_var(ncid, var_id,&
&     arr_int, start = (/ 1 /) , count = (/ narr /) )

   else if (typevar=='DPR') then

!    call etsf_io_low_write_var(ncid, trim(varname), arr_real(1:narr), lstat,&
!    & start = (/ 1 /) , count = (/ narr /),error_data=error_data)
!    call etsf_io_low_error_handle(error_data)

!    write(std_out,*)  trim(varname), lstat,'NARR',narr,'MARR',marr,arr_real(:)
     ncerr = nf90_put_var(ncid, var_id,&
&     arr_real, start = (/ 1 /) , count = (/ narr /) )
   end if
!  if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable")

 end if

#endif

end subroutine write_var_netcdf
!!***
