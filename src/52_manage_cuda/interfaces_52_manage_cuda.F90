!!****m* ABINIT/interfaces_52_manage_cuda
!! NAME
!! interfaces_52_manage_cuda
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/52_manage_cuda
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

module interfaces_52_manage_cuda

 implicit none

interface
 subroutine alloc_hamilt_gpu(atindx1,dtset,gprimd,mpi_enreg,nattyp,option,psps,use_gpu_cuda)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: option
  integer,intent(in) :: use_gpu_cuda
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: nattyp(dtset%ntypat)
 end subroutine alloc_hamilt_gpu
end interface

interface
 subroutine dealloc_hamilt_gpu(option,use_gpu_cuda)
  implicit none
  integer,intent(in) :: option
  integer,intent(in) :: use_gpu_cuda
 end subroutine dealloc_hamilt_gpu
end interface

end module interfaces_52_manage_cuda
!!***
