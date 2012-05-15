!!****m* ABINIT/defs_scfcvargs
!! NAME
!! defs_scfcvargs
!!
!! FUNCTION
!! This module contains definitions of high-level structured datatypes for the
!! ABINIT package, related to the routine scfcv.F90.
!!
!! If you are sure a new high-level structured datatype is needed,
!! write it here, and DOCUMENT it properly (not all datastructure here are
!! well documented, it is a shame ...).
!! Do not forget : you will likely be the major winner if you document
!! properly.
!! Proper documentation of a structured datatype means :
!!  (1) Mention it in the list just below
!!  (2) Describe it in the NOTES section
!!  (3) Put it in alphabetical order in the the main section of this module
!!  (4) Document each of its records, except if they are described elsewhere
!!      (this exception is typically the case of the dataset associated with
!!      input variables, for which there is a help file)
!!  (5) Declare variables on separated lines in order to reduce the occurence 
!!      of bzr conflicts.
!!
!! List of datatypes :
!! * scfcvargs   : a container for all the arguments except rprimd and xred
!!                 of the scfcv routine
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module defs_scfcvargs

 use m_profiling

 use defs_basis

! Contains: type(MPI_type)
!           type(datafiles_type)
!           type(dataset_type)
 use defs_abitypes

! Contains: type(macro_uj_type)
!           type(pawrhoij_type)
!           type(pawrad_type)
!           type(pawtab_type)
!           type(hdr_type)
!           type(pawang_type)
!           type(pawfgr_type)
!           type(pseudopotential_type)
!           type(results_gs_type)
 use defs_datatypes

! Contains: type(electronpositron_type)
 use m_electronpositron, only : electronpositron_type

!           type(efield_type)
 use m_efield

 use defs_scftypes

! Contains: type(recursion_type)
! use defs_rectypes

! Contains: type(wvl_data)
 !use defs_wvltypes

! Contains: type(paw_dmft_type)
 use m_paw_dmft, only: paw_dmft_type

! Contains: type(m_wffile)
 use m_wffile

!#if defined HAVE_DFT_BIGDFT
! use BigDFT_API, only : atoms_data
!#endif

 implicit none

!Structures
!!***

!!****t* defs_scftypes/ab_scfcvargs
!! NAME
!! ab_scfcvargs
!!
!! FUNCTION
!! This datatype has the purpouse of store all the
!! arguments needed by the routine scfcv 
!!
!! NOTES
!!
!!
!! SOURCE

type ab_scfcvargs

!scalars
 type(ab_scfcv_args_in),pointer :: ab_scfcv_in
 type(ab_scfcv_args_inout),pointer :: ab_scfcv_inout
 type(dataset_type),pointer :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(paw_dmft_type),pointer :: paw_dmft
 type(wffile_type),pointer :: wffnew,wffnow
! TODO: check if we need to add efield object here
! type(efield_type), pointer :: dtefield
!arrays

end type ab_scfcvargs

end module defs_scfcvargs
!!***
