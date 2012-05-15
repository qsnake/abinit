!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_etsf_init
!! NAME
!! abi_etsf_init
!!
!! FUNCTION
!!  Create a NetCDF file following the ETSF file format specifications.
!!  It declares the dimensions and set-up the different variables, according
!!  to the itype argument.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2012 ABINIT group (Yann Pouillon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  filapp = character string giving the root to form the name of the GEO file
!!  itype = an integer to define what to put in the output file. This can
!!          be one of the following values (maybe a sum latter):
!!          1 for a density file,
!!          2 for a wavefunction file,
!!          4 for a KSS file,
!!          8 for the exchange potential,
!!         16 for the correlation potential.
!!  kdep= .true. if the data for the array sizes are dependant on k points.
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  Data written in file whose name is filapp//'-etsf.nc'
!!
!! PARENTS
!!      m_io_kss,scfcv,sigma
!!
!! CHILDREN
!!      abi_etsf_dims_init,etsf_io_data_init,etsf_io_low_close
!!      etsf_io_low_error_to_str,etsf_io_low_open_modify,etsf_io_main_def
!!      ini_wf_etsf,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine abi_etsf_init(dtset, filapp, itype, kdep, lmn_size, psps, wfs)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
 use defs_abitypes
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
 use m_abi_etsf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_init'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_61_ionetcdf, except_this_one => abi_etsf_init
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itype
 logical,intent(in) :: kdep
 character(len=fnlen),intent(in) :: filapp
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),intent(in) :: wfs
!arrays
 integer,intent(in) :: lmn_size(psps%npsp)

!Local variables-------------------------------
#if defined HAVE_TRIO_ETSF_IO
 type(etsf_dims) :: dims
 type(etsf_groups_flags) :: flags
 logical :: lstat
 integer :: ncid, var_main
 type(etsf_io_low_error) :: error
 character(len=etsf_io_low_error_len) :: errmess
 character(len=80) :: file_title
#endif
!scalars
 character(len=500) :: message
 character(len=fnlen) :: filetsf

! *************************************************************************

#if defined HAVE_TRIO_ETSF_IO
!Initialize the filename
 filetsf = trim(filapp)//'-etsf.nc'
 write(message, '(a,a,a)' ) ch10,' abi_etsf_init : about to create file ',&
& filetsf
 call wrtout(std_out,message,'COLL')

!Set-up the dimensions
!=====================
 call abi_etsf_dims_init(dims, dtset, itype, psps, wfs)

!Set-up the variables
!====================
!These mandatory values are always written by the hdr_io_etsf() routine.
 flags%geometry  = etsf_geometry_all
 flags%kpoints   = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights
 flags%electrons = etsf_electrons_all - etsf_electrons_x_functional - etsf_electrons_c_functional
 flags%basisdata = etsf_basisdata_basis_set
 if (dtset%usewvl == 0) then
   flags%basisdata = flags%basisdata + etsf_basisdata_kin_cutoff + etsf_basisdata_n_coeff
 end if
!These variables may be written depending on prt<something> input variables.
 if (itype == 1) then
   flags%main      = etsf_main_density
   write(file_title, "(A)") "Density file"
 else if (itype == 2) then
   if (dtset%usewvl == 0) then
     flags%basisdata = flags%basisdata + etsf_basisdata_red_coord_pw
   else
     flags%basisdata = flags%basisdata + etsf_basisdata_coord_grid + &
&     etsf_basisdata_n_coeff_grid
   end if
   flags%main      = etsf_main_wfs_coeff
   write(file_title, "(A)") "Wavefunctions file"
 else if (itype == 4) then
   if (dtset%usewvl == 0) then
     flags%basisdata = flags%basisdata + etsf_basisdata_red_coord_pw
   else
     flags%basisdata = flags%basisdata + etsf_basisdata_coord_grid + &
&     etsf_basisdata_n_coeff_grid
   end if
   flags%main      = etsf_main_wfs_coeff
   flags%gwdata    = etsf_gwdata_all
   write(file_title, "(A)") "KSS file"
 else if (itype == 8) then
   flags%main      = etsf_main_pot_x_only
   write(file_title, "(A)") "Exchange potential file"
 else if (itype == 16) then
   flags%main      = etsf_main_pot_c_only
   write(file_title, "(A)") "Correlation potential file"
 else if (itype == 24) then
   flags%main      = etsf_main_pot_xc
   write(file_title, "(A)") "Exchange-correlation potential file"
 end if

!Actually create the file
!========================
!If the group contains main, we remove it for a while to be sure to
!add it at the end, after ABINIT private variables.
 var_main = flags%main
 flags%main = etsf_main_none
 call etsf_io_data_init(filetsf, flags, dims, file_title, &
& 'File generated by ABINIT with ETSF_IO', lstat, error, &
& overwrite = .true., k_dependent = kdep)
 if (.not.lstat) goto 1000

!We now add the private ABINIT informations when required.
!We open the file to add the private ABINIT variables
 call etsf_io_low_open_modify(ncid, trim(filetsf), lstat, error_data = error)
 if (.not.lstat) goto 1000
!We add the private data
 call ini_wf_etsf(dtset, lmn_size, psps%npsp, psps%ntypat, ncid)

!We now add the main part as last variables in the ETSF file.
!We add the main group
 call etsf_io_main_def(ncid, lstat, error, flags = var_main)
 if (.not.lstat) goto 1000
 
!We close the file.
 call etsf_io_low_close(ncid, lstat, error_data = error)
 if (.not.lstat) goto 1000

 1000 continue
 if (.not. lstat) then
!  We handle the error
   call etsf_io_low_error_to_str(errmess, error)
   write(message, "(A,A,A,A)") ch10, " abi_etsf_init: ERROR -", ch10, &
&   errmess(1:min(475, len(errmess)))
   call wrtout(std_out, message, 'COLL')
   call leave_new('COLL')
 end if

#endif

end subroutine abi_etsf_init
!!***
