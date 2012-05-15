!!****m* ABINIT/interfaces_42_nlstrain
!! NAME
!! interfaces_42_nlstrain
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/42_nlstrain
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

module interfaces_42_nlstrain

 implicit none

interface
 subroutine contistr01(istr,rank,gm,gprimd,eisnl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: bb(2,((rank+2)*(rank+3))/2)
  real(dp),intent(out) :: eisnl(3)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contistr01
end interface

interface
 subroutine contistr03(istr,rank,gm,gprimd,eisnl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: bb(2,((rank+4)*(rank+5))/2)
  real(dp),intent(out) :: eisnl(3)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contistr03
end interface

interface
 subroutine contistr10(istr,rank,gm,gprimd,eisnl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+2)*(rank+3))/2)
  real(dp),intent(in) :: bb(2,((rank+1)*(rank+2))/2)
  real(dp),intent(out) :: eisnl(3)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contistr10
end interface

interface
 subroutine contistr12(istr,rank,gm,gprimd,eisnl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+2)*(rank+3))/2)
  real(dp),intent(in) :: bb(2,((rank+3)*(rank+4))/2)
  real(dp),intent(out) :: eisnl(3)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contistr12
end interface

interface
 subroutine contistr21(istr,rank,gm,gprimd,eisnl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+3)*(rank+4))/2)
  real(dp),intent(in) :: bb(2,((rank+2)*(rank+3))/2)
  real(dp),intent(out) :: eisnl(3)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contistr21
end interface

interface
 subroutine contistr30(istr,rank,gm,gprimd,eisnl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+4)*(rank+5))/2)
  real(dp),intent(in) :: bb(2,((rank+1)*(rank+2))/2)
  real(dp),intent(out) :: eisnl(3)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contistr30
end interface

interface
 subroutine contstr21(istr1,istr2,rank,gm,gprimd,e2nl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr1
  integer,intent(in) :: istr2
  integer,intent(in) :: rank
  real(dp),intent(out) :: e2nl
  real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: bb(2,((rank+6)*(rank+5))/2)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contstr21
end interface

interface
 subroutine contstr22(istr1,istr2,rank,gm,gprimd,e2nl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr1
  integer,intent(in) :: istr2
  integer,intent(in) :: rank
  real(dp),intent(out) :: e2nl
  real(dp),intent(in) :: aa(2,((rank+5)*(rank+6))/2)
  real(dp),intent(in) :: bb(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contstr22
end interface

interface
 subroutine contstr23(istr1,istr2,rank,gm,gprimd,e2nl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr1
  integer,intent(in) :: istr2
  integer,intent(in) :: rank
  real(dp),intent(out) :: e2nl
  real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: bb(2,((rank+3)*(rank+4))/2)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contstr23
end interface

interface
 subroutine contstr24(istr1,istr2,rank,gm,gprimd,e2nl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr1
  integer,intent(in) :: istr2
  integer,intent(in) :: rank
  real(dp),intent(out) :: e2nl
  real(dp),intent(in) :: aa(2,((rank+3)*(rank+4))/2)
  real(dp),intent(in) :: bb(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contstr24
end interface

interface
 subroutine contstr25(istr1,istr2,rank,gm,gprimd,e2nl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr1
  integer,intent(in) :: istr2
  integer,intent(in) :: rank
  real(dp),intent(out) :: e2nl
  real(dp),intent(in) :: aa(2,((rank+3)*(rank+4))/2)
  real(dp),intent(in) :: bb(2,((rank+3)*(rank+4))/2)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contstr25
end interface

interface
 subroutine contstr25a(istr1,istr2,rank,gm,gprimd,e2nl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr1
  integer,intent(in) :: istr2
  integer,intent(in) :: rank
  real(dp),intent(out) :: e2nl
  real(dp),intent(in) :: aa(2,((rank+3)*(rank+4))/2)
  real(dp),intent(in) :: bb(2,((rank+3)*(rank+4))/2)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contstr25a
end interface

interface
 subroutine contstr26(istr1,istr2,rank,gm,gprimd,e2nl,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr1
  integer,intent(in) :: istr2
  integer,intent(in) :: rank
  real(dp),intent(out) :: e2nl
  real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: bb(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: gm(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine contstr26
end interface

end module interfaces_42_nlstrain
!!***
