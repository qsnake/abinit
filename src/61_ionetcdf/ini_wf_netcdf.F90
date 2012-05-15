!{\src2tex{textfont=tt}}
!!****f* ABINIT/ini_wf_netcdf
!! NAME
!! ini_wf_netcdf
!!
!! FUNCTION
!! Do initialization of additional dimensions and variables in
!! wavefunction files in NetCDF format.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpw = maximal number of plane waves per kpoint
!!  ncid_hdr = id of netcdf file handle
!!  response = 0 for GS case and 1 for RF case
!!
!! OUTPUT
!!  <NetCDF file attached to ncid_hdr is eventually modified>
!!
!! NOTES
!!
!!  Normally mpw and response should not change, so if they already exist
!!  they are not updated to the new input values...
!!
!! PARENTS
!!      inwffil,newsp,outwf,uderiv,vtorho,vtorho3
!!
!! CHILDREN
!!      handle_ncerr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ini_wf_netcdf(mpw,ncid_hdr,response)

 use m_profiling

 use defs_basis
 use defs_datatypes
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ini_wf_netcdf'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpw,ncid_hdr,response

!Local variables-------------------------------
!no_abirules
#if defined HAVE_TRIO_NETCDF
        integer :: ncerr
        integer :: dimg3_id,mband_id,nkpt_id,nspinor_id,nsppol_id
        integer :: complex_id,eigendim_id,mpw_id
        integer :: kg_id,eigen_id,cg_id
        integer :: test_id,mband
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' ini_wf_netcdf: enter'
!write(std_out,*)' ini_wf_netcdf: arguments = ', mpw,ncid_hdr,response
!stop
!ENDDEBUG

#if defined HAVE_TRIO_NETCDF
 ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="eigendim",dimid=test_id)

!if additional dimensions, kg, eigen, and so on are not defined yet... do so
 if (ncerr /= NF90_NOERR) then

!  get ids for needed dimensions that already exist
   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nspinor",dimid=nspinor_id)
   call handle_ncerr(ncerr," inquire nspinor")
   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nsppol",dimid=nsppol_id)
   call handle_ncerr(ncerr," inquire nsppol")
   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nkpt",dimid=nkpt_id)
   call handle_ncerr(ncerr," inquire nkpt")
   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dimg3",dimid=dimg3_id)
   call handle_ncerr(ncerr," inquire dimg3")
   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="mband",dimid=mband_id)
   call handle_ncerr(ncerr," inquire mband")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=mband_id,len=mband)
   call handle_ncerr(ncerr," get mband")


!  Add dimensions and variable ids for eigen,occ,kg, cg
!  
!  definition :   cg(2,mpw*nspinor*mband*mkmem*nsppol)
!  here saved as: cg(2,mpw,nspinor,mband,nkpt,nsppol)
!  definition :   eigen((2*mband)**response*mband*nkpt*nsppol)
!  here saved as: eigen((2*mband)**response*mband,nkpt,nsppol)
!  definition :   kg(3,mpw*mkmem)
!  here saved as: kg(3,mpw,nkpt)
!  
!  go to definition mode
   ncerr = nf90_redef(ncid=ncid_hdr)
   call handle_ncerr(ncerr," return to definition mode ")

   ncerr = nf90_def_dim(ncid=ncid_hdr,name="mpw",len=mpw,dimid=mpw_id)
   call handle_ncerr(ncerr," define mpw")
   ncerr = nf90_def_dim(ncid=ncid_hdr,name="complex",len=2,dimid=complex_id)
   call handle_ncerr(ncerr," define complex")
   ncerr = nf90_def_dim(ncid=ncid_hdr,name="eigendim",len=(2*mband)**response*mband,dimid=eigendim_id)
   call handle_ncerr(ncerr," define eigendim")

   ncerr = nf90_def_var(ncid=ncid_hdr,name="kg",xtype=NF90_INT,&
&   dimids=(/dimg3_id,mpw_id,nkpt_id/),varid=kg_id)
   call handle_ncerr(ncerr," define kg")
!  DEBUG
!  write(std_out,*) 'ini_wf_netcdf : inquiring kg'
!  ncerr = nf90_inq_varid(ncid=ncid_hdr,name="kg",varid=kg_id)
!  call handle_ncerr(ncerr," inquire kg ")
!  ENDDEBUG

   ncerr = nf90_def_var(ncid=ncid_hdr,name="eigen",xtype=NF90_DOUBLE,&
&   dimids=(/eigendim_id,nkpt_id,nsppol_id/),varid=eigen_id)
   call handle_ncerr(ncerr," define eigen")
   ncerr = nf90_def_var(ncid=ncid_hdr,name="cg",xtype=NF90_DOUBLE,&
&   dimids=(/complex_id,mpw_id,nspinor_id,mband_id,nkpt_id,nsppol_id/),varid=cg_id)
   call handle_ncerr(ncerr," define cg")

!  end netCDF definition mode
   ncerr = nf90_enddef(ncid=ncid_hdr)
   call handle_ncerr(ncerr," end definition mode")
!  probably not needed
   ncerr = nf90_sync(ncid=ncid_hdr)
   call handle_ncerr(ncerr," sync file ")

   write(std_out,*) 'ini_wf_netcdf : finished defining variables in netCDF'

 end if
!end test on whether the file is initialized already

#endif
!if NETCDF is undefined, do nothing

end subroutine ini_wf_netcdf
!!***
