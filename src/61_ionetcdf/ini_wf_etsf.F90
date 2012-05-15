!{\src2tex{textfont=tt}}
!!****f* ABINIT/ini_wf_etsf
!! NAME
!! ini_wf_etsf
!!
!! FUNCTION
!! Do initialization of additional dimensions and variables in wavefunction files in ETSF format.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  npsp=number of pseudopotentials.
!!  ntypat=number of types of atoms in cell.
!!  unwff=the unit of the open NetCDF file.
!!
!! SIDE EFFECTS 
!!  New dimensions and variables are added to the initial NetCDF file.
!!
!! PARENTS
!!      abi_etsf_init,pawmkaewf
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_write_dim
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ini_wf_etsf(dtset, lmn_size, npsp, ntypat, unwff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
#if defined HAVE_TRIO_ETSF_IO
  use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ini_wf_etsf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp,ntypat,unwff
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: lmn_size(npsp)

!Local variables-------------------------------
#if defined HAVE_TRIO_ETSF_IO
 type(etsf_io_low_error) :: error_data
 logical :: lstat
 integer :: rhoijdim1
#endif

! *************************************************************************

#if defined HAVE_TRIO_ETSF_IO
!Add the none-ETSF dimensions and variables.
!If dimensions already exist, it will check that definitions are coherent.

!Define dimensions.
 call etsf_io_low_write_dim(unwff, "npsp", npsp, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(unwff, "codvsnlen", 6, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!version 44 add first dimension for rhoij = max(lmn_size)*(max(lmn_size)+1)/2
 rhoijdim1 = maxval(lmn_size)
 rhoijdim1 = rhoijdim1 * (rhoijdim1 + 1) / 2
!impose rhoijdim1 >= 1 : if 0, it defaults to NF90_UNLIMITED
 rhoijdim1 = max(rhoijdim1, 1)
 call etsf_io_low_write_dim(unwff, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(unwff, "psptitlen", 132, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 if (dtset%usewvl == 1) then ! Add the BigDFT private dimensions.
   call etsf_io_low_write_dim(unwff, "number_of_wavelet_resolutions", 2, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if
 
!If variables already exist, it will check that definitions are coherent.

!Define variables.
 call etsf_io_low_def_var(unwff, "date", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "codvsn", etsf_io_low_character, (/ "codvsnlen" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "ecut_eff", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "ecutsm", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "etot", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)
 
 call etsf_io_low_def_var(unwff, "headform", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "fform", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "intxc", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "ixc", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "occopt", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "pertcase", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "pertcase", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "residm", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "stmbias", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "tphysel", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "tsmear", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Version 44 add ecutdg and usepaw
 call etsf_io_low_def_var(unwff, "ecutdg", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "usepaw", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Multi-dimensional variables.
 call etsf_io_low_def_var(unwff, "pspcod", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "pspdat", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "pspso", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "pspxc", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "qptn", etsf_io_low_double, (/ "number_of_reduced_dimensions" /), &
& lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "so_psp", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "symafm", etsf_io_low_integer, (/ "number_of_symmetry_operations" /), &
& lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "title", etsf_io_low_character, (/ pad("psptitlen"), pad("npsp") /), &
& lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 if (npsp /= ntypat) then
   call etsf_io_low_def_var(unwff, "zionpsp", etsf_io_low_double, (/ "npsp" /), lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

 call etsf_io_low_def_var(unwff, "znuclpsp", etsf_io_low_double, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Version 44 add lmn_size and rhoij
 call etsf_io_low_def_var(unwff, "lmn_size", etsf_io_low_integer, (/ "number_of_atom_species" /), &
& lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(unwff, "rhoij", etsf_io_low_double, &
& (/ pad("rhoijdim1"), pad("number_of_components"), pad("number_of_atoms") /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Add BigDFT variables.
 call etsf_io_low_def_var(unwff, "usewvl", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Add the BigDFT private variables.
 if (dtset%usewvl == 1) then
   call etsf_io_low_def_var(unwff, "number_of_wavelets", etsf_io_low_integer, (/ "number_of_wavelet_resolutions" /), &
&   lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

#endif
!if ETSF_IO is undefined, do nothing

end subroutine ini_wf_etsf
!!***
