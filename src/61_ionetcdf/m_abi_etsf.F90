!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abi_etsf
!! NAME
!! m_abi_etsf
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2012 ABINIT group (DCA,YP,MJV,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abi_etsf

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_errors
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 implicit none

 private

#if defined HAVE_TRIO_ETSF_IO
 public :: abi_etsf_dims_init   

CONTAINS  !===========================================================
!!***

!!****f* m_abi_etsf/abi_etsf_dims_init
!! NAME
!! abi_etsf_dims_init
!!
!! FUNCTION
!!  Initialize the structure defining ETSF dimensions.
!!  starting from values stored in the dataset_type, the pseudopotential_type.
!!  and the wave function handler for the BIGDFT part.
!!
!! INPUTS
!!  dtset<type(dataset_type)>=all input variables for this dataset
!!  psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs<wvl_wf_type>=Object to handle wave functions for the BIGDFT part
!!    Presently not used, likely will be needed when ETSF-IO will be generalized to deal with PAW.
!!  itype = an integer to define what to put in the output file. This can
!!          be one of the following values (maybe a sum latter):
!!          1 for a density file,
!!          2 for a wavefunction file,
!!          4 for a KSS file,
!!          8 for the exchange potential,
!!         16 for the correlation potential.
!!
!! OUTPUT
!!  dims=structure with ETSF dimensions.
!!
!! PARENTS
!!      abi_etsf_init,pawmkaewf
!!
!! CHILDREN
!!
!! SOURCE

subroutine abi_etsf_dims_init(dims, dtset, itype, psps, wfs)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_dims_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itype
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),intent(in) :: wfs
 type(etsf_dims),intent(inout) :: dims
! *************************************************************************

!Set-up the dimensions
!=====================
 dims%max_number_of_angular_momenta  = psps%mpsang
!In the case of BigDFT, the number of coefficients are the number of wavelets.
 if (dtset%usewvl == 0) then
  dims%max_number_of_coefficients      = dtset%mpw
  dims%max_number_of_basis_grid_points = etsf_no_dimension
#if defined HAVE_DFT_BIGDFT
 else
  dims%max_number_of_coefficients      = wfs%Glr%wfd%nvctr_c + 7 * wfs%Glr%wfd%nvctr_f
  dims%max_number_of_basis_grid_points = wfs%Glr%wfd%nvctr_c
#endif
 end if
 dims%max_number_of_projectors       = 1
 dims%max_number_of_states           = dtset%mband
 dims%number_of_atoms                = dtset%natom
 dims%number_of_atom_species         = dtset%ntypat
 dims%number_of_components           = dtset%nspden
 if (dtset%usepaw == 1) then
  dims%number_of_grid_points_vector1  = dtset%ngfftdg(1)
  dims%number_of_grid_points_vector2  = dtset%ngfftdg(2)
  dims%number_of_grid_points_vector3  = dtset%ngfftdg(3)
#if defined HAVE_DFT_BIGDFT
 else if (dtset%usewvl == 1) then
!In the case of BigDFT, the grid size is not defined by ngfft.
  dims%number_of_grid_points_vector1  = wfs%Glr%d%n1 * 2
  dims%number_of_grid_points_vector2  = wfs%Glr%d%n2 * 2
  dims%number_of_grid_points_vector3  = wfs%Glr%d%n3 * 2
#endif
 else
  dims%number_of_grid_points_vector1  = dtset%ngfft(1)
  dims%number_of_grid_points_vector2  = dtset%ngfft(2)
  dims%number_of_grid_points_vector3  = dtset%ngfft(3)
 end if
 dims%number_of_kpoints              = dtset%nkpt
 dims%number_of_spinor_components    = dtset%nspinor
 dims%number_of_spins                = dtset%nsppol
 dims%number_of_symmetry_operations  = dtset%nsym
!The density real_or_complex.
 if (iand(itype, 1) /= 0) then
  dims%real_or_complex_density        = 1
 else
  dims%real_or_complex_density        = etsf_no_dimension
 end if
!The coefficient of wavefunctions real_or_complex.
 if (iand(itype, 2) /= 0 .or. iand(itype, 4) /= 0) then
  if (dtset%usewvl == 0) then
   dims%real_or_complex_coefficients= 2 ! used in plane waves
  else
   dims%real_or_complex_coefficients= 1 ! used in wavelets
  end if
 else
  dims%real_or_complex_coefficients   = etsf_no_dimension
 end if
!The gw corrections real_or_complex.
!Todo: Currently not exported.
!if (.false. .and. iand(itype, 4) /= 0) then
 if (iand(itype, 4) /= 0) then
  dims%real_or_complex_gw_corrections = 2 ! used in plane waves
! dims%real_or_complex_gw_corrections = 1 ! used in plane waves
 else
  dims%real_or_complex_gw_corrections = etsf_no_dimension
 end if
!The potential real_or_complex.
 if (iand(itype, 8) /= 0 .or. iand(itype, 16) /= 0) then
  dims%real_or_complex_potential      = 1
 else
  dims%real_or_complex_potential      = etsf_no_dimension
 end if
 dims%real_or_complex_wavefunctions  = etsf_no_dimension

end subroutine abi_etsf_dims_init
!!***

#endif

END MODULE m_abi_etsf
!!***
