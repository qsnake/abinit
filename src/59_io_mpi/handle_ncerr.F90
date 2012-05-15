!{\src2tex{textfont=tt}}
!!****f* ABINIT/handle_ncerr
!! NAME
!! handle_ncerr
!!
!! FUNCTION
!! Rudimentary Error catching for NetCDF using routines: prints out a string
!! describing what was being done when the error code showed up.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (Mver)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!   ncerr = error code. If different from the preset nf90_noerr, then error
!!   message = subroutine specified message (description of action which
!!    produced the error)
!!
!! OUTPUT
!!  (only writes)
!!
!! PARENTS
!!      ab_define_var,cut3d,gstate,hdr_io_netcdf,ini_wf_netcdf,inwffil,ioarr
!!      loper3,newsp,outwf,outxfhist,read_md_hist,rwwf,uderiv,wffile,wffopen
!!      wffreadnpwrec,write_eig,write_md_hist,wrt_moldyn_netcdf
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine handle_ncerr(ncerr,message)

 use m_profiling

use defs_basis
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'handle_ncerr'
!End of the abilint section

 implicit none

!Arguments --------------------------
 integer         ,intent(in) :: ncerr
 character(len=*),intent(in) :: message

!Local variables --------------------

! *********************************************************************

#if defined HAVE_TRIO_NETCDF
 if (ncerr /= nf90_noerr) then
   write(std_out,*) 'Error in netcdf call while : ', trim(message)
   write(std_out,*) trim(nf90_strerror(ncerr))
   stop
 end if
#else
 write(std_out,*) ' handle_ncerr : Error : NETCDF not defined at compile time. handle_ncerr should not be called.'
 stop
#endif

end subroutine handle_ncerr
!!***
