!!****m* ABINIT/interfaces_62_cg_noabirule
!! NAME
!! interfaces_62_cg_noabirule
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/62_cg_noabirule
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

module interfaces_62_cg_noabirule

 implicit none

interface
 subroutine bracketing (nv1,nv2,dp_dum_v2dp,v,grad,a,x,b,fa,fx,fb)
  use defs_basis
  implicit none
  integer,intent(in) :: nv1
  integer,intent(in) :: nv2
  real(dp),intent(inout) :: a
  real(dp),intent(out) :: b
  real(dp),intent(out) :: fa
  real(dp),intent(out) :: fb
  real(dp),intent(out) :: fx
  real(dp),intent(inout) :: x
  interface
   function dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)::dp_dum_v2dp
   end function dp_dum_v2dp
  end interface
  real(dp),intent(inout) :: grad(nv1,nv2)
  real(dp),intent(inout) :: v(nv1,nv2)
 end subroutine bracketing
end interface

interface
 function brent(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,itmax,v,grad,ax,xx,bx,tol,xmin)
  use defs_basis
  implicit none
  integer,intent(in) :: itmax
  integer,intent(in) :: nv1
  integer,intent(in) :: nv2
  real(dp),intent(in) :: ax
  real(dp) :: brent
  real(dp),intent(in) :: bx
  real(dp),intent(in) :: tol
  real(dp),intent(out) :: xmin
  real(dp),intent(in) :: xx
  interface
   function dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)::dp_dum_v2dp
   end function dp_dum_v2dp
  end interface
  interface
   function v2dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)            :: v2dp_dum_v2dp(nv1,nv2)
   end function v2dp_dum_v2dp
  end interface
  interface
   subroutine sub_dum_dp_v2dp_v2dp(nv1,nv2,inarg1, inarg2, ioarg3)
    use defs_basis
    use defs_datatypes
    integer, intent(in)    :: nv1,nv2
    real(dp),intent(in)    :: inarg1
    real(dp),intent(inout) :: inarg2(nv1,nv2)
    real(dp),intent(inout):: ioarg3(nv1,nv2)
   end subroutine sub_dum_dp_v2dp_v2dp
  end interface
  real(dp),intent(inout) :: grad(nv1,nv2)
  real(dp),intent(inout) :: v(nv1,nv2)
 end function brent
end interface

interface
 subroutine cgpr(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,dtol,itmax,v,fmin,delta)
  use defs_basis
  implicit none
  integer,intent(in) :: itmax
  integer,intent(in) :: nv1
  integer,intent(in) :: nv2
  real(dp),intent(out) :: delta
  real(dp),intent(in) :: dtol
  real(dp),intent(out) :: fmin
  interface
   function dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)::dp_dum_v2dp
   end function dp_dum_v2dp
  end interface
  interface
   function v2dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)            :: v2dp_dum_v2dp(nv1,nv2)
   end function v2dp_dum_v2dp
  end interface
  interface
   subroutine sub_dum_dp_v2dp_v2dp(nv1,nv2,inarg1, inarg2, ioarg3)
    use defs_basis
    use defs_datatypes
    integer, intent(in)    :: nv1,nv2
    real(dp),intent(in)    :: inarg1
    real(dp),intent(inout) :: inarg2(nv1,nv2)
    real(dp),intent(inout):: ioarg3(nv1,nv2)
   end subroutine sub_dum_dp_v2dp_v2dp
  end interface
  real(dp),intent(inout) :: v(nv1,nv2)
 end subroutine cgpr
end interface

interface
 function dotproduct(nv1,nv2,v1,v2)
  use defs_basis
  implicit none
  integer,intent(in) :: nv1
  integer,intent(in) :: nv2
  real(dp) :: dotproduct
  real(dp),intent(in) :: v1(nv1,nv2)
  real(dp),intent(in) :: v2(nv1,nv2)
 end function dotproduct
end interface

interface
 subroutine linmin(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,v,grad,fmin)
  use defs_basis
  implicit none
  integer,intent(in) :: nv1
  integer,intent(in) :: nv2
  real(dp),intent(out) :: fmin
  interface
   function dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)::dp_dum_v2dp
   end function dp_dum_v2dp
  end interface
  interface
   function v2dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)            :: v2dp_dum_v2dp(nv1,nv2)
   end function v2dp_dum_v2dp
  end interface
  interface
   subroutine sub_dum_dp_v2dp_v2dp(nv1,nv2,inarg1, inarg2, ioarg3)
    use defs_basis
    use defs_datatypes
    integer, intent(in)    :: nv1,nv2
    real(dp),intent(in)    :: inarg1
    real(dp),intent(inout) :: inarg2(nv1,nv2)
    real(dp),intent(inout):: ioarg3(nv1,nv2)
   end subroutine sub_dum_dp_v2dp_v2dp
  end interface
  real(dp),intent(inout) :: grad(nv1,nv2)
  real(dp),intent(inout) :: v(nv1,nv2)
 end subroutine linmin
end interface

end module interfaces_62_cg_noabirule
!!***
