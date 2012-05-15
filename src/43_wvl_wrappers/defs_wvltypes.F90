!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_wvltypes
!! NAME
!! defs_wvltypes
!!
!! FUNCTION
!! This module contains definitions of all structured datatypes for the
!! wavelet part of the ABINIT package.
!!
!! List of datatypes :
!! * wvl_projectors_type : structure to store projectors for wavelets calculations.
!! * wvl_wf_type : structure to store wavefunctions for wavelets calculations.
!! * wvl_data : container for all required wavelets data.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module defs_wvltypes

 use m_profiling

 use defs_basis
#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only : atoms_data, orbitals_data, locreg_descriptors, &
      & nonlocal_psp_descriptors, GPU_pointers, rho_descriptors, orthon_data, &
      & SIC_data, diis_objects, communications_arrays
#endif

 implicit none
!!***

 !!****t* defs_wvltypes/wvl_projectors_type
 !! NAME
 !! wvl_projectors_type
 !!
 !! FUNCTION
 !! This type constains all physical data associated to
 !! projectors of a given system (Z and atomic positions).
 !!
 !! SOURCE
 type wvl_projectors_type
    ! These arrays are compacted arrays. One should use a nonlocal_psp_descriptors
    ! object to read these arrays and grep values on real space grid.
#if defined HAVE_DFT_BIGDFT
    type(nonlocal_psp_descriptors) :: keys
#endif

    ! Data for projectors.
    !  size (%nprojel).
    real(dp), pointer :: proj(:)

 end type wvl_projectors_type
!!***

 !!****t* defs_wvltypes/wvl_wf_type
 !! NAME
 !! wvl_wf_type
 !!
 !! FUNCTION
 !! This type constains all physical data associated to
 !! wavefunctions of a given system (Z and atomic positions).
 !!
 !! SOURCE
 type wvl_wf_type
#if defined HAVE_DFT_BIGDFT
    ! This stores the repartition of all orbitals between,
    ! spinor, k-points, spin...
    type(orbitals_data) :: orbs
    ! Contains the information needed for describing completely a
    ! wavefunction localisation region
    type(locreg_descriptors) :: Glr

    ! Some GPU internals.
    type(GPU_pointers) :: GPU

    ! Some internal pointers to deal with DIIS wfs optimisation.
    type(diis_objects) :: diis

    ! Parameters to distribute wavefunctions among processors
    type(communications_arrays) :: comms
#endif

    ! wavefunctions, size (mvctrp * mbandp)
    real(dp), pointer :: psi(:)
    ! wavefunction gradients, size (mvctrp * mbandp)
    real(dp), pointer :: hpsi(:)

    ! Temporary wavefunction storage when several proc are used.
    real(dp), pointer :: psit(:)
 end type wvl_wf_type
!!***

!!****t* defs_wvltypes/wvl_internal_type
!! NAME
!! wvl_internal_type
!!
!! FUNCTION
!! This type is a gathering for all internal variables wavelets required. It is
!! included in the datatypes strutcture.
!!
!! NOTES
!! If you modify this datatype, the copy, init/free routines are in wvl_utils.F90.
!!
!! SOURCE
type wvl_internal_type
  real(dp) :: h(3)
  ! The hgrid values in each direction, after the application of the
  ! boundary conditions. In free boundary conditions, the three values are equal.
  real(dp) :: shift(3)
  ! Shift applied by BigDFT on the atomic position (in cartesian coordinates).

#if defined HAVE_DFT_BIGDFT
  type(locreg_descriptors) :: Glr
  ! Contains the description of the global localisation region.

  type(atoms_data) :: atoms
  ! A copy of the current dtset values.

  type(rho_descriptors) :: rhodsc
  ! Internals to handle the rho array.

  type(orthon_data) :: orthpar
  ! Some parameters to orthogonalise.

  type(SIC_data) :: SIC
  ! Self-interaction correction coefficients.
#endif

  character(len=4) :: exctxpar
  ! parallelisation scheme of the exact exchange operator
  !   BC (Blocking Collective)
  !   OP2P (Overlap Point-to-Point)  
  real(dp) :: exctxfac
  ! Exact exchange coefficient to mix.  
end type wvl_internal_type
!!***

!!****t* defs_wvltypes/wvl_data
!! NAME
!! wvl_data
!!
!! FUNCTION
!! This type is a container to limit the number of arguments in
!! ABINIT, it should not contains attributes other than types
!! created for wavelets.
!!
!! SOURCE
type wvl_data
   ! The data associated to projectors (computed from pseudo)
   type(wvl_projectors_type) :: projectors
   ! The data associated to the wavefunctions
   type(wvl_wf_type) :: wfs
   ! The pointers to internal BigDFT data structures
   type(wvl_internal_type) :: descr
end type wvl_data
!!***
end module defs_wvltypes
