!!****m* ABINIT/interfaces_44_abitypes_defs
!! NAME
!! interfaces_44_abitypes_defs
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/44_abitypes_defs
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

module interfaces_44_abitypes_defs

 implicit none

interface
 subroutine cprj_alloc(cprj,ncpgr,nlmn)
  use defs_datatypes
  implicit none
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nlmn(:)
  type(cprj_type),intent(inout) :: cprj(:,:)
 end subroutine cprj_alloc
end interface

interface
 subroutine cprj_free(cprj)
  use defs_datatypes
  implicit none
  type(cprj_type),intent(inout) :: cprj(:,:)
 end subroutine cprj_free
end interface

interface
 subroutine cprj_nullify(cprj)
  use defs_datatypes
  implicit none
  type(cprj_type),intent(inout) :: cprj(:,:)
 end subroutine cprj_nullify
end interface

interface
 subroutine cprj_set_zero(cprj)
  use defs_datatypes
  implicit none
  type(cprj_type),intent(inout) :: cprj(:,:)
 end subroutine cprj_set_zero
end interface

interface
 subroutine cprj_copy(cprj_in,cprj_out,&  
  &  icpgr) ! optional argument
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: icpgr
  type(cprj_type),intent(in) :: cprj_in(:,:)
  type(cprj_type),intent(inout) :: cprj_out(:,:)
 end subroutine cprj_copy
end interface

interface
 subroutine cprj_axpby(alpha,beta,cprjx,cprjy)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(in) :: alpha
  real(dp),intent(in) :: beta
  type(cprj_type),intent(in) :: cprjx(:,:)
  type(cprj_type),intent(inout) :: cprjy(:,:)
 end subroutine cprj_axpby
end interface

interface
 subroutine cprj_zaxpby(alpha,beta,cprjx,cprjy)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(in) :: alpha(2)
  real(dp),intent(in) :: beta(2)
  type(cprj_type),intent(in) :: cprjx(:,:)
  type(cprj_type),intent(inout) :: cprjy(:,:)
 end subroutine cprj_zaxpby
end interface

interface
 subroutine cprj_lincom(alpha,cprj_in,cprj_out,nn)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nn
  real(dp),intent(in) :: alpha(2,nn)
  type(cprj_type),intent(in) :: cprj_in(:,:)
  type(cprj_type),intent(inout) :: cprj_out(:,:)
 end subroutine cprj_lincom
end interface

interface
 function paw_overlap(cprj1,cprj2,typat,pawtab,spinor_comm) result(onsite)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: spinor_comm
  integer,intent(in) :: typat(:)
  type(cprj_type),intent(in) :: cprj1(:,:)
  type(cprj_type),intent(in) :: cprj2(:,:)
  real(dp) :: onsite(2)
  type(pawtab_type),intent(in) :: pawtab(:)
 end function paw_overlap
end interface

interface
 subroutine cprj_output(cprj)
  use defs_datatypes
  implicit none
  type(cprj_type),intent(in) :: cprj(:,:)
 end subroutine cprj_output
end interface

interface
 subroutine rhoij_alloc(cplex,nlmn,nspden,nspinor,nsppol,pawrhoij,typat,&  ! Mandatory arguments
  &  mpi_enreg,ngrhoij,nlmnmix,use_rhoij_,use_rhoijres) ! Optional arguments
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in),optional :: ngrhoij
  integer,intent(in),optional :: nlmnmix
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in),optional :: use_rhoij_
  integer,intent(in),optional :: use_rhoijres
  type(mpi_type),intent(inout),optional :: mpi_enreg
  integer,intent(in) :: nlmn(:)
  integer,intent(in) :: typat(:)
  type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 end subroutine rhoij_alloc
end interface

interface
 subroutine rhoij_free(pawrhoij)
  use defs_datatypes
  implicit none
  type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 end subroutine rhoij_free
end interface

interface
 subroutine rhoij_nullify(pawrhoij)
  use defs_datatypes
  implicit none
  type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 end subroutine rhoij_nullify
end interface

interface
 subroutine rhoij_copy(pawrhoij_in,pawrhoij_out,&  
  &  keep_cplex,keep_nspden,mpi_enreg) ! optional arguments
  use defs_abitypes
  use defs_datatypes
  implicit none
  logical,intent(in),optional :: keep_cplex
  logical,intent(in),optional :: keep_nspden
  type(mpi_type),intent(inout),optional :: mpi_enreg
  type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
  type(pawrhoij_type),intent(inout) :: pawrhoij_out(:)
 end subroutine rhoij_copy
end interface

interface
 subroutine rhoij_allgather(mpi_enreg,pawrhoij_in,pawrhoij_gathered)
  use defs_abitypes
  use defs_datatypes
  implicit none
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawrhoij_type),intent(inout) :: pawrhoij_gathered(:)
  type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
 end subroutine rhoij_allgather
end interface

interface
 subroutine rhoij_io(pawrhoij,unitfi,nsppol_in,nspinor_in,nspden_in,nlmn_type,typat,headform,rdwr_mode,form,natinc)
  use defs_datatypes
  implicit none
  integer,intent(in) :: headform
  integer,optional,intent(in) :: natinc
  integer,intent(in) :: nspden_in
  integer,intent(in) :: nspinor_in
  integer,intent(in) :: nsppol_in
  integer,intent(in) :: unitfi
  character(len=*),optional,intent(in) :: form
  character(len=*),intent(in) :: rdwr_mode
  integer,intent(in) :: nlmn_type(:)
  integer,intent(in) :: typat(:)
  type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 end subroutine rhoij_io
end interface

interface
 subroutine rhoij_unpack(rhoij)
  use defs_datatypes
  implicit none
  type(pawrhoij_type),intent(inout) :: rhoij(:)
 end subroutine rhoij_unpack
end interface

interface
 subroutine rhoij_init_unpacked(rhoij)
  use defs_datatypes
  implicit none
  type(pawrhoij_type),intent(inout) :: rhoij(:)
 end subroutine rhoij_init_unpacked
end interface

interface
 subroutine rhoij_destroy_unpacked(rhoij)
  use defs_datatypes
  implicit none
  type(pawrhoij_type),intent(inout) :: rhoij(:)
 end subroutine rhoij_destroy_unpacked
end interface

interface
 subroutine rhoij_mpi_sum(pawrhoij,comm1,comm2)
  use defs_datatypes
  implicit none
  integer,intent(in) :: comm1
  integer,optional,intent(in) :: comm2
  type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 end subroutine rhoij_mpi_sum
end interface

end module interfaces_44_abitypes_defs
!!***
