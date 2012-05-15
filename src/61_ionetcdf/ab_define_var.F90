!{\src2tex{textfont=tt}}
!!****f* ABINIT/ab_define_var
!!
!! NAME
!! ab_define_var
!!
!! FUNCTION
!! Write the definition of a variable, including units and mnemonics
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
!! ncid = Identifier of the netcdf dataset
!! var_dim_id = Identifier of the Dimensions
!! var_id     = Identifier of the variable
!! var_mnemo  = String of mnemonics
!! var_name   = String with the name of the variable
!! var_type   = NetCDF type of variable (NF90_DOUBLE, etc)
!! var_units  = String of units
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      write_eig,write_md_hist
!!
!! CHILDREN
!!      handle_ncerr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ab_define_var(ncid, var_dim_id, var_id, var_type, var_name, var_mnemo, var_units)

 use m_profiling

! define dp,sixth,third,etc...
 use defs_basis

#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab_define_var'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ncid 
 integer, intent(out) :: var_id
 character(len=*), intent(in) :: var_mnemo,var_units
 character(len=*), intent(in) :: var_name
 integer,intent(in) :: var_type
!arrays
 integer,intent(in) :: var_dim_id(:)

!Local variables-------------------------------
!scalars
integer :: ncerr
!arrays

! *************************************************************************

#if defined HAVE_TRIO_NETCDF

 ncerr = nf90_def_var(ncid, trim(var_name), var_type, var_dim_id, var_id)
 if (ncerr /= NF90_NOERR)  call handle_ncerr(ncerr," define variable "//trim(var_name))

 ncerr = nf90_put_att(ncid, var_id,  "units",trim(var_units))
 if ( ncerr /= NF90_NOERR) call handle_ncerr(ncerr," define attribute for "//trim(var_name))

 ncerr = nf90_put_att(ncid, var_id,  "mnemonics", trim(var_mnemo))
 if ( ncerr /= NF90_NOERR) call handle_ncerr(ncerr," define attribute for "//trim(var_name))

#endif

end subroutine ab_define_var
!!***
