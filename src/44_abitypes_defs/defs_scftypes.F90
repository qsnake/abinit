!!****m* ABINIT/defs_scftypes
!! NAME
!! defs_scftypes
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
!!  (5) Declare variables on separated lines in order to reduce the occurence of bzr conflicts.
!!
!! List of datatypes :
!! * ab_scfcv_args_in   : a container for all the input arguments of the
!!                   scfcv routine
!! * ab_scfcv_args_inout: a container for all the input/output arguments
!!                   of the scfcv routine
!!
!!                   These variables are excluded:
!!
!!                   type(electronpositron_type) :: electronpositron
!!                     (Definition too high 56>50)
!!                   type(wffile_type) :: wffnew,wffnow
!!                     (Definition too high 51>50)
!!                   type(paw_dmft_type) :: paw_dmft
!!                     (Definition too high 66>50)
!!
!!                   real(dp), pointer :: xred(:,:),xred_old(:,:)
!!                     (General argument)
!!                   real(dp), pointer :: rprimd(:,:)
!!                     (General argument)
!!
!!                  type(dataset_type) :: dtset
!!                     (Too basic, used inside mover)
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


module defs_scftypes

 use m_profiling

 use defs_basis
 use m_results_gs , only : results_gs_type

! Contains: type(MPI_type)
!           type(datafiles_type)
!           type(dataset_type)
 use defs_abitypes
 use m_scf_history, only : scf_history_type

! Contains: type(efield_type)
 use m_efield

! Contains: type(macro_uj_type)
!           type(pawrhoij_type)
!           type(pawrad_type)
!           type(pawtab_type)
!           type(hdr_type)
!           type(pawang_type)
!           type(pawfgr_type)
!           type(pseudopotential_type)
!           type(results_gs_type)
!           type(scf_history_type)
 use defs_datatypes

! Contains: type(electronpositron_type)
! use m_electronpositron, only : electronpositron_type

! Contains: type(recursion_type)
 use defs_rectypes

! Contains: type(wvl_data)
 use defs_wvltypes

! Contains: type(paw_dmft_type)
! use m_paw_dmft, only: paw_dmft_type

! Contains: type(m_wffile)
! use m_wffile

!#if defined HAVE_DFT_BIGDFT
! use BigDFT_API, only : atoms_data
!#endif

 implicit none

!Structures
!!***

!!****t* defs_scftypes/ab_scfcv_args_in
!! NAME
!! ab_scfcv_args_in
!!
!! FUNCTION
!! This datatype has the purpouse of store all the input only
!! arguments needed by the routine scfcv
!!
!! NOTES
!!
!!
!! SOURCE

type ab_scfcv_args_in

!scalars
 integer,pointer :: iapp,mcg,ndtpawuj,pwind_alloc
 real(dp),pointer :: cpus,ecore
 real(dp),pointer :: fatvshift
 type(pawang_type),pointer :: pawang
 type(pseudopotential_type),pointer :: psps

!arrays
 integer,pointer :: atindx(:),atindx1(:)
 integer,pointer :: indsym(:,:,:)
!no_abirules
 integer, pointer :: kg(:,:)
 integer, pointer :: nattyp(:),npwarr(:),pwind(:,:,:)
 real(dp), pointer :: phnons(:,:,:)
 real(dp), pointer :: pwnsfac(:,:)
 real(dp), pointer :: ylm(:,:)
 real(dp), pointer :: ylmgr(:,:,:)
 type(pawrad_type), pointer :: pawrad(:)
 type(pawtab_type), pointer :: pawtab(:)

end type ab_scfcv_args_in
!!***

!!****t* defs_scftypes/ab_scfcv_args_inout
!! NAME
!! ab_scfcv_args_inout
!!
!! FUNCTION
!! This datatype has the purpouse of store most of the
!! ioput/output arguments needed by the routine scfcv
!!
!! NOTES
!!
!!
!! SOURCE

type ab_scfcv_args_inout

!scalars
  integer,pointer :: initialized,nfftf
  type(MPI_type),pointer :: mpi_enreg
  type(datafiles_type),pointer :: dtfil
!  type(dataset_type),pointer :: dtset
  type(efield_type),pointer :: dtefield
!  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),pointer :: hdr
  type(pawfgr_type),pointer :: pawfgr
  type(recursion_type),pointer :: rec_set
  type(results_gs_type),pointer :: results_gs
  type(scf_history_type),pointer :: scf_history
!  type(wffile_type),pointer :: wffnew,wffnow
  type(wvl_data),pointer :: wvl
!  type(paw_dmft_type), pointer :: paw_dmft

!arrays
  integer, pointer :: irrzon(:,:,:)
  integer, pointer :: symrec(:,:,:)
  real(dp), pointer :: cg(:,:)
  real(dp), pointer :: eigen(:)
  real(dp), pointer :: occ(:)
!  real(dp), pointer :: rprimd(:,:)
  real(dp), pointer :: rhog(:,:),rhor(:,:)
  real(dp), pointer :: taug(:,:),taur(:,:)
  real(dp), pointer :: resid(:)
!  real(dp), pointer :: xred(:,:),xred_old(:,:)
  type(macro_uj_type),pointer :: dtpawuj(:)
  type(pawrhoij_type), pointer :: pawrhoij(:)

end type ab_scfcv_args_inout
!!***

end module defs_scftypes
!!***
