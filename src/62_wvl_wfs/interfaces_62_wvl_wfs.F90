!!****m* ABINIT/interfaces_62_wvl_wfs
!! NAME
!! interfaces_62_wvl_wfs
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/62_wvl_wfs
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

module interfaces_62_wvl_wfs

 implicit none

interface
 subroutine wvl_nl_gradient(grnl, mpi_enreg, natom, rprimd, wvl, xcart)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: natom
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(inout) :: grnl(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine wvl_nl_gradient
end interface

interface
 subroutine wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)
  use defs_basis
  use defs_abitypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer, intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(in) :: hdr
  type(hdr_type), intent(in) :: hdr0
  type(mpi_type), intent(in) :: mpi_enreg
  type(wffile_type),intent(in) :: wff
  type(wvl_wf_type), intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(in) :: xred(3, dtset%natom)
 end subroutine wvl_read
end interface

interface
 subroutine wvl_write(dtset, eigen, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)
  use defs_basis
  use defs_abitypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer, intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(in) :: mpi_enreg
  type(wffile_type),intent(in) :: wff
  type(wvl_wf_type), intent(in) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in), target :: eigen(dtset%mband)
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(inout) :: xred(3, dtset%natom)
 end subroutine wvl_write
end interface

interface
 subroutine wvl_setngfft(ixc, mgfft, mpi_enreg, natom, nfft, ngfft, nsppol, psps, rprimd,&  
  &  wvl, wvl_crmult, wvl_frmult, xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: ixc
  integer, intent(out) :: mgfft
  integer, intent(in) :: natom
  integer, intent(out) :: nfft
  integer, intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(inout) :: wvl
  real(dp), intent(in) :: wvl_crmult
  real(dp), intent(in) :: wvl_frmult
  integer, intent(out) :: ngfft(13)
  real(dp), intent(inout) :: rprimd(3,3)
  real(dp), intent(inout) :: xred(3, natom)
 end subroutine wvl_setngfft
end interface

interface
 subroutine wvl_tail_corrections(dtset, energies, etotal, mpi_enreg, psps,&  
  &  vtrial, wvl, xcart)
  use defs_basis
  use m_energies
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(in),target :: vtrial(dtset%nfft)
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine wvl_tail_corrections
end interface

end module interfaces_62_wvl_wfs
!!***
