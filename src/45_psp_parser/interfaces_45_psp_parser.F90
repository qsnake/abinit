!!****m* ABINIT/interfaces_45_psp_parser
!! NAME
!! interfaces_45_psp_parser
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/45_psp_parser
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

module interfaces_45_psp_parser

 implicit none

interface
 subroutine inpspheads(filnam,npsp,pspheads)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: npsp
  character(len=fnlen), intent(in) :: filnam(npsp)
  type(pspheader_type),intent(out) :: pspheads(npsp)
 end subroutine inpspheads
end interface

interface
 subroutine psxml2ab( psxml, znucl, zion, pspcod, pspxc, lmax, iwrite )
  use defs_basis
  use m_xml_pseudo_types
  implicit none
  integer,intent(in) :: iwrite
  integer,intent(out) :: lmax
  integer,intent(out) :: pspcod
  integer,intent(out) :: pspxc
  type(pseudo_t),intent(in) :: psxml
  real(dp),intent(out) :: zion
  real(dp),intent(out) :: znucl
 end subroutine psxml2ab
end interface

interface
 subroutine upfheader2abi (filpsp,&  
  &  znucl, zion, pspxc,&  
  &  lmax_, n1xccc, nproj_l, nprojso_l)
  use defs_basis
  implicit none
  integer, intent(out) :: lmax_
  integer, intent(inout) :: n1xccc
  integer, intent(out) :: pspxc
  character(len=fnlen), intent(in) :: filpsp
  real(dp), intent(out) :: zion
  real(dp), intent(out) :: znucl
  integer, intent(out) :: nproj_l(0:3)
  integer, intent(out) :: nprojso_l(1:3)
 end subroutine upfheader2abi
end interface

interface
 subroutine upfxc2abi(dft, pspxc)
  implicit none
  integer, intent(out) :: pspxc
  character(len=20), intent(in) :: dft
 end subroutine upfxc2abi
end interface

end module interfaces_45_psp_parser
!!***
