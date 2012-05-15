!!****m* ABINIT/interfaces_11_qespresso_ext
!! NAME
!! interfaces_11_qespresso_ext
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/11_qespresso_ext
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

module interfaces_11_qespresso_ext

 implicit none

interface
 subroutine read_pseudo (is, iunps)
  implicit none
  integer :: is
  integer :: iunps
 end subroutine read_pseudo
end interface

interface
 subroutine scan_begin (iunps, string, rew)
  implicit none
  integer :: iunps
  logical :: rew
  character (len=*) :: string
 end subroutine scan_begin
end interface

interface
 subroutine scan_end (iunps, string)
  implicit none
  integer :: iunps
  character (len=*) :: string
 end subroutine scan_end
end interface

interface
 subroutine read_pseudo_header (is, iunps)
  implicit none
  integer :: is
  integer :: iunps
 end subroutine read_pseudo_header
end interface

interface
 subroutine read_pseudo_local (is, iunps)
  implicit none
  integer :: is
  integer :: iunps
 end subroutine read_pseudo_local
end interface

interface
 subroutine read_pseudo_mesh (is, iunps)
  implicit none
  integer :: is
  integer :: iunps
 end subroutine read_pseudo_mesh
end interface

interface
 subroutine read_pseudo_nl (is, iunps)
  implicit none
  integer :: is
  integer :: iunps
 end subroutine read_pseudo_nl
end interface

interface
 subroutine read_pseudo_nlcc (is, iunps)
  implicit none
  integer :: is
  integer :: iunps
 end subroutine read_pseudo_nlcc
end interface

interface
 subroutine read_pseudo_pswfc (is, iunps)
  implicit none
  integer :: is
  integer :: iunps
 end subroutine read_pseudo_pswfc
end interface

interface
 subroutine read_pseudo_rhoatom (is, iunps)
  implicit none
  integer :: is
  integer :: iunps
 end subroutine read_pseudo_rhoatom
end interface

end module interfaces_11_qespresso_ext
!!***
