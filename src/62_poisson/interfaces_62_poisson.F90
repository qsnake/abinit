!!****m* ABINIT/interfaces_62_poisson
!! NAME
!! interfaces_62_poisson
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/62_poisson
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

module interfaces_62_poisson

 implicit none

interface
 subroutine Psolver_hartree(dtset, enhartr, mpi_enreg, rhor, rprimd, vhartr, wvl)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  real(dp), intent(out) :: enhartr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  real(dp),intent(in) :: rhor(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vhartr(dtset%nfft)
 end subroutine Psolver_hartree
end interface

interface
 subroutine PSolver_kernel(dtset, iaction, kernel_array, mpi_enreg, rprimd, wvl)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: iaction
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), pointer :: kernel_array(:)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine PSolver_kernel
end interface

interface
 subroutine Psolver_rhohxc(dtset, enhartr, enxc, envxc, mpi_enreg, rhor, rprimd,&  
  &  vhartr, vxc, vxcavg, wvl)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  real(dp), intent(out) :: enhartr
  real(dp), intent(out) :: envxc
  real(dp), intent(out) :: enxc
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp), intent(out) :: vxcavg
  type(wvl_internal_type), intent(in) :: wvl
  real(dp),intent(inout) :: rhor(dtset%nfft, dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vhartr(dtset%nfft)
  real(dp),intent(out) :: vxc(dtset%nfft, dtset%nspden)
 end subroutine Psolver_rhohxc
end interface

end module interfaces_62_poisson
!!***
