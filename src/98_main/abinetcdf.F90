!{\src2tex{textfont=tt}}
!!****p* ABINIT/abinetcdf
!! NAME
!! abinetcdf
!!
!! FUNCTION
!! Tests if NetCDF support is working properly.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (JPMinet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  (print all)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program abinetcdf

 use defs_basis
 use m_xmpi, only : xmpi_world
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abinetcdf'
!End of the abilint section

 implicit none
!Arguments -----------------------------------

!Local variables-------------------------------
!scalars
 integer :: EnID,dim_len,i,kxID,ncid,status
 character(len=8) :: att_value,dim_name
!arrays
 real(dp) :: energy(10),energy_read(10)

! *************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

 call random_number(energy)

#if defined HAVE_TRIO_NETCDF
 write(std_out,'(2a)') "We are using NetCDF version ", NF90_INQ_LIBVERS()

!Creating dataset with one variable of dimenssion 10, with one attribute

 status=NF90_CREATE(path="test.nc", cmode=NF90_CLOBBER, ncid=ncid)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_DEF_DIM(ncid, "kx", 10, kxID)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_DEF_VAR(ncid, "En", NF90_DOUBLE, kxID, EnID)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_PUT_ATT(ncid, EnID, "units", "anything")
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_ENDDEF(ncid)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_PUT_VAR(ncid, EnID, energy)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_CLOSE(ncid)
 if (status /= NF90_NOERR) call handle_err(status)

!Reading dataset and comparing with original

 status=NF90_OPEN(path="test.nc", mode=NF90_NOWRITE , ncid=ncid)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_Inquire_Dimension(ncid, 1, dim_name, dim_len)
 if (status /= NF90_NOERR) call handle_err(status)
 if (dim_name /= "kx") call problem("Inquire_dimension")
 if (dim_len /= 10) call problem("Inquire_dimension")

 status=NF90_get_att(ncid, 1,"units",att_value)
 if (status /= NF90_NOERR) call handle_err(status)
 if (att_value /= "anything") call problem("get_att")

 status=NF90_get_var(ncid=ncid, varid=1, values=energy_read)
 if (status /= NF90_NOERR) call handle_err(status)
 do i=1,10
   if (energy_read(i) /= energy(i)) call problem("get var")
 end do

 status=NF90_close(ncid)
 if (status /= NF90_NOERR) call handle_err(status)

 write(std_out,'(a)') "Test OK: NetCDF correctly integrated in ABINIT tree :-)"
#else
 write(std_out,'(a)') "This build of ABINIT does not provide NetCDF support"
#endif

 end program abinetcdf
!!***

!!****f* ABINIT/handle_err
!! NAME
!! handle_err
!!
!! FUNCTION
!!
!! PARENTS
!!      abinetcdf
!!
!! CHILDREN
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine handle_err(status)

 use m_profiling

 use defs_basis
#if defined HAVE_TRIO_NETCDF
  use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'handle_err'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: status

!Local variables-------------------------------

! *************************************************************************

#if defined HAVE_TRIO_NETCDF
   if (status /= NF90_NOERR) then
     write(std_out,'(a)') trim(NF90_STRERROR(status))
     stop "Stopped"
   end if
#endif
 end subroutine handle_err
!!***

!!****f* ABINIT/problem
!! NAME
!! problem
!!
!! FUNCTION
!!
!! PARENTS
!!      abinetcdf
!!
!! CHILDREN
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine problem(msg)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'problem'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=*),intent(in) :: msg

!Local variables-------------------------------

! *************************************************************************

   write(std_out,'(2a)') "Problem with NetCDF function ", msg
   stop
 end subroutine problem
!!***
