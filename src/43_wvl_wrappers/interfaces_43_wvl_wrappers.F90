!!****m* ABINIT/interfaces_43_wvl_wrappers
!! NAME
!! interfaces_43_wvl_wrappers
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/43_wvl_wrappers
!!
!! COPYRIGHT
!! Copyright (C) 2010-2011 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_43_wvl_wrappers

 implicit none

interface
 subroutine wvl_descr_atoms_set(acell, icoulomb, natom, ntypat, typat, wvl)
  use defs_basis
  use defs_wvltypes
  implicit none
  integer, intent(in) :: icoulomb
  integer, intent(in) :: natom
  integer, intent(in) :: ntypat
  type(wvl_internal_type), intent(inout) :: wvl
  real(dp), intent(in) :: acell(3)
  integer, intent(in) :: typat(natom)
 end subroutine wvl_descr_atoms_set
end interface

interface
 subroutine wvl_descr_free(wvl)
  use defs_wvltypes
  implicit none
  type(wvl_internal_type), intent(inout) :: wvl
 end subroutine wvl_descr_free
end interface

interface
 subroutine wvl_descr_psp_set(nsppol, psps, wvl)
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: nsppol
  type(pseudopotential_type), intent(in) :: psps
  type(wvl_internal_type), intent(inout) :: wvl
 end subroutine wvl_descr_psp_set
end interface

interface
 subroutine wvl_projectors_free(proj)
  use defs_wvltypes
  implicit none
  type(wvl_projectors_type),intent(inout) :: proj
 end subroutine wvl_projectors_free
end interface

interface
 subroutine wvl_projectors_set(me, natom, proj, psps, rprimd, wfs,&  
  &  wvl, wvl_frmult, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: me
  integer, intent(in) :: natom
  type(wvl_projectors_type),intent(out) :: proj
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_wf_type),intent(in) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: wvl_frmult
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine wvl_projectors_set
end interface

interface
 subroutine wvl_wfs_free(wfs)
  use defs_wvltypes
  implicit none
  type(wvl_wf_type),intent(inout) :: wfs
 end subroutine wvl_wfs_free
end interface

interface
 subroutine wvl_wfs_lr_copy(wfs, wvl)
  use defs_wvltypes
  implicit none
  type(wvl_wf_type), intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
 end subroutine wvl_wfs_lr_copy
end interface

interface
 subroutine wvl_wfs_set(fixmom, kpt, me, natom, nband, nkpt, nproc, nspinor,&  
  &  nsppol, nwfshist, occ, psps, rprimd, wfs, wtk, wvl, wvl_crmult, wvl_frmult, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: me
  integer, intent(in) :: natom
  integer, intent(in) :: nband
  integer, intent(in) :: nkpt
  integer, intent(in) :: nproc
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  integer, intent(in) :: nwfshist
  real(dp), intent(in) :: fixmom
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_wf_type),intent(out) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: wvl_crmult
  real(dp), intent(in) :: wvl_frmult
  real(dp), intent(in) :: kpt(3,nkpt)
  real(dp), intent(in) :: occ(:)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: wtk(nkpt)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine wvl_wfs_set
end interface

interface
 subroutine derfcf(derfc_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derfc_yy
  real(dp),intent(in) :: yy
 end subroutine derfcf
end interface

interface
 subroutine derf_ab(derf_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derf_yy
  real(dp),intent(in) :: yy
 end subroutine derf_ab
end interface

end module interfaces_43_wvl_wrappers
!!***
