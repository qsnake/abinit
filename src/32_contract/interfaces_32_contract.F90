!!****m* ABINIT/interfaces_32_contract
!! NAME
!! interfaces_32_contract
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/32_contract
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

module interfaces_32_contract

 implicit none

interface
 subroutine contract_dp_ge_val(calling_routine,dp_name,dp_value,limit)
  use defs_basis
  implicit none
  character(len=*),intent(in) :: calling_routine
  character(len=*),intent(in) :: dp_name
  real(dp),intent(in) :: dp_value
  real(dp),intent(in) :: limit
 end subroutine contract_dp_ge_val
end interface

interface
 subroutine contract_int_ge_val(calling_routine,integer_name,integer_value,limit)
  implicit none
  integer,intent(in) :: integer_value
  integer,intent(in) :: limit
  character(len=*),intent(in) :: calling_routine
  character(len=*),intent(in) :: integer_name
 end subroutine contract_int_ge_val
end interface

interface
 subroutine contract_int_le_val(calling_routine,integer_name,integer_value,limit)
  implicit none
  integer,intent(in) :: integer_value
  integer,intent(in) :: limit
  character(len=*),intent(in) :: calling_routine
  character(len=*),intent(in) :: integer_name
 end subroutine contract_int_le_val
end interface

interface
 subroutine contract_int_list(calling_routine,integer_name,integer_value,list,nlist)
  implicit none
  integer,intent(in) :: integer_value
  integer,intent(in) :: nlist
  character(len=*),intent(in) :: calling_routine
  character(len=*),intent(in) :: integer_name
  integer,intent(in) :: list(nlist)
 end subroutine contract_int_list
end interface

end module interfaces_32_contract
!!***
