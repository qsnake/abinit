!!****m* ABINIT/interfaces_28_numeric_noabirule
!! NAME
!! interfaces_28_numeric_noabirule
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/28_numeric_noabirule
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

module interfaces_28_numeric_noabirule

 implicit none

interface
 function ass_leg_pol(l,m,xarg)
  implicit none
  integer, intent(in) :: l
  integer, intent(in) :: m
  double precision :: ass_leg_pol
  double precision, intent(in) :: xarg
 end function ass_leg_pol
end interface

interface
 subroutine coeffs_gausslegint(xmin,xmax,x,weights,n)
  implicit none
  integer :: n
  double precision :: xmax
  double precision :: xmin
  double precision :: weights(n)
  double precision :: x(n)
 end subroutine coeffs_gausslegint
end interface

interface
 subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result )
  use defs_basis
  implicit none
  integer, intent(in) :: ntab
  real(dp), intent(in) :: a
  real(dp), intent(in) :: b
  real(dp), intent(out) :: result
  real(dp), intent(inout) :: e(ntab)
  real(dp), intent(in) :: ftab(ntab)
  real(dp), intent(inout) :: work(ntab)
  real(dp), intent(in) :: xtab(ntab)
  real(dp), intent(inout) :: y(3,ntab)
 end subroutine cspint
end interface

interface
 subroutine dzgedi(a,lda,n,ipvt,det,work,job)
  implicit none
  integer :: job
  integer :: lda
  integer :: n
  real*8 :: det(2,2)
  real*8 :: a(2,lda,n)
  integer :: ipvt(n)
  real*8 :: work(2,n)
 end subroutine dzgedi
end interface

interface
 subroutine dzgefa(a,lda,n,ipvt,info)
  implicit none
  integer :: info
  integer :: lda
  integer :: n
  real*8 :: a(2,lda,n)
  integer :: ipvt(n)
 end subroutine dzgefa
end interface

interface
 subroutine GAMMA_FUNCTION(X,GA)
  use defs_basis
  implicit none
  real(dp),intent(out) :: ga
  real(dp),intent(in) :: x
 end subroutine GAMMA_FUNCTION
end interface

interface
 function interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function interp
end interface

interface
 function dinterp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: dinterp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function dinterp
end interface

interface
 function taylor_interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: taylor_interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function taylor_interp
end interface

interface
 function dtaylor_interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: dtaylor_interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function dtaylor_interp
end interface

interface
 subroutine calculate_taylor_c(n,z,f,z0,c)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: z0
  complex(gwpc) :: c(n)
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end subroutine calculate_taylor_c
end interface

interface
 SUBROUTINE INTRPL(L,X,Y,N,U,V,dv,dv2,ideriv)
  implicit none
  integer, parameter :: NQQ=12000
  integer :: L
  integer :: N
  integer :: ideriv
  double precision :: DV(NQQ)
  double precision :: DV2(NQQ)
  double precision :: U(N)
  double precision :: V(N)
  double precision :: X(L)
  double precision :: Y(L)
 end subroutine INTRPL
end interface

interface
 SUBROUTINE CALJY0(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALJY0
end interface

interface
 DOUBLE PRECISION FUNCTION BESJ0(X)
 implicit none
 double precision :: X
end function BESJ0
end interface

interface
 DOUBLE PRECISION FUNCTION BESY0(X)
 implicit none
 double precision :: X
end function BESY0
end interface

interface
 SUBROUTINE CALJY1(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALJY1
end interface

interface
 DOUBLE PRECISION FUNCTION BESJ1(X)
 implicit none
 double precision :: X
end function BESJ1
end interface

interface
 DOUBLE PRECISION FUNCTION BESY1(X)
 implicit none
 double precision :: X
end function BESY1
end interface

interface
 subroutine jacobi(a,n,np,d,v,nrot)
  implicit none
  integer :: n
  integer :: np
  integer :: nrot
  real*8 :: a(np,np)
  real*8 :: d(np)
  real*8 :: v(np,np)
 end subroutine jacobi
end interface

interface
 SUBROUTINE CALCK0(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALCK0
end interface

interface
 DOUBLE PRECISION FUNCTION BESK0(X)
 implicit none
 double precision :: X
end function BESK0
end interface

interface
 DOUBLE PRECISION FUNCTION BESEK0(X)
 implicit none
 double precision :: X
end function BESEK0
end interface

interface
 SUBROUTINE CALCK1(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALCK1
end interface

interface
 DOUBLE PRECISION FUNCTION BESK1(X)
 implicit none
 double precision :: X
end function BESK1
end interface

interface
 DOUBLE PRECISION  FUNCTION BESEK1(X)
 implicit none
 double precision :: X
end function BESEK1
end interface

interface
 subroutine gen_oh(code, num, x, y, z, w, a, b, v)
  implicit none
  integer :: code
  integer :: num
  double precision :: a
  double precision :: b
  double precision :: v
  double precision :: w(*)
  double precision :: x(*)
  double precision :: y(*)
  double precision :: z(*)
 end subroutine gen_oh
end interface

interface
 SUBROUTINE LD0006(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 6)
  double precision :: X( 6)
  double precision :: Y( 6)
  double precision :: Z( 6)
 end subroutine LD0006
end interface

interface
 SUBROUTINE LD0014(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 14)
  double precision :: X( 14)
  double precision :: Y( 14)
  double precision :: Z( 14)
 end subroutine LD0014
end interface

interface
 SUBROUTINE LD0026(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 26)
  double precision :: X( 26)
  double precision :: Y( 26)
  double precision :: Z( 26)
 end subroutine LD0026
end interface

interface
 SUBROUTINE LD0038(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 38)
  double precision :: X( 38)
  double precision :: Y( 38)
  double precision :: Z( 38)
 end subroutine LD0038
end interface

interface
 SUBROUTINE LD0050(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 50)
  double precision :: X( 50)
  double precision :: Y( 50)
  double precision :: Z( 50)
 end subroutine LD0050
end interface

interface
 SUBROUTINE LD0074(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 74)
  double precision :: X( 74)
  double precision :: Y( 74)
  double precision :: Z( 74)
 end subroutine LD0074
end interface

interface
 SUBROUTINE LD0086(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 86)
  double precision :: X( 86)
  double precision :: Y( 86)
  double precision :: Z( 86)
 end subroutine LD0086
end interface

interface
 SUBROUTINE LD0110(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 110)
  double precision :: X( 110)
  double precision :: Y( 110)
  double precision :: Z( 110)
 end subroutine LD0110
end interface

interface
 SUBROUTINE LD0146(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 146)
  double precision :: X( 146)
  double precision :: Y( 146)
  double precision :: Z( 146)
 end subroutine LD0146
end interface

interface
 SUBROUTINE LD0170(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 170)
  double precision :: X( 170)
  double precision :: Y( 170)
  double precision :: Z( 170)
 end subroutine LD0170
end interface

interface
 SUBROUTINE LD0194(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 194)
  double precision :: X( 194)
  double precision :: Y( 194)
  double precision :: Z( 194)
 end subroutine LD0194
end interface

interface
 SUBROUTINE LD0230(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 230)
  double precision :: X( 230)
  double precision :: Y( 230)
  double precision :: Z( 230)
 end subroutine LD0230
end interface

interface
 SUBROUTINE LD0266(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 266)
  double precision :: X( 266)
  double precision :: Y( 266)
  double precision :: Z( 266)
 end subroutine LD0266
end interface

interface
 SUBROUTINE LD0302(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 302)
  double precision :: X( 302)
  double precision :: Y( 302)
  double precision :: Z( 302)
 end subroutine LD0302
end interface

interface
 SUBROUTINE LD0350(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 350)
  double precision :: X( 350)
  double precision :: Y( 350)
  double precision :: Z( 350)
 end subroutine LD0350
end interface

interface
 SUBROUTINE LD0434(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 434)
  double precision :: X( 434)
  double precision :: Y( 434)
  double precision :: Z( 434)
 end subroutine LD0434
end interface

interface
 SUBROUTINE LD0590(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 590)
  double precision :: X( 590)
  double precision :: Y( 590)
  double precision :: Z( 590)
 end subroutine LD0590
end interface

interface
 SUBROUTINE LD0770(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 770)
  double precision :: X( 770)
  double precision :: Y( 770)
  double precision :: Z( 770)
 end subroutine LD0770
end interface

interface
 SUBROUTINE LD0974(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W( 974)
  double precision :: X( 974)
  double precision :: Y( 974)
  double precision :: Z( 974)
 end subroutine LD0974
end interface

interface
 SUBROUTINE LD1202(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(1202)
  double precision :: X(1202)
  double precision :: Y(1202)
  double precision :: Z(1202)
 end subroutine LD1202
end interface

interface
 SUBROUTINE LD1454(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(1454)
  double precision :: X(1454)
  double precision :: Y(1454)
  double precision :: Z(1454)
 end subroutine LD1454
end interface

interface
 SUBROUTINE LD1730(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(1730)
  double precision :: X(1730)
  double precision :: Y(1730)
  double precision :: Z(1730)
 end subroutine LD1730
end interface

interface
 SUBROUTINE LD2030(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(2030)
  double precision :: X(2030)
  double precision :: Y(2030)
  double precision :: Z(2030)
 end subroutine LD2030
end interface

interface
 SUBROUTINE LD2354(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(2354)
  double precision :: X(2354)
  double precision :: Y(2354)
  double precision :: Z(2354)
 end subroutine LD2354
end interface

interface
 SUBROUTINE LD2702(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(2702)
  double precision :: X(2702)
  double precision :: Y(2702)
  double precision :: Z(2702)
 end subroutine LD2702
end interface

interface
 SUBROUTINE LD3074(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(3074)
  double precision :: X(3074)
  double precision :: Y(3074)
  double precision :: Z(3074)
 end subroutine LD3074
end interface

interface
 SUBROUTINE LD3470(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(3470)
  double precision :: X(3470)
  double precision :: Y(3470)
  double precision :: Z(3470)
 end subroutine LD3470
end interface

interface
 SUBROUTINE LD3890(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(3890)
  double precision :: X(3890)
  double precision :: Y(3890)
  double precision :: Z(3890)
 end subroutine LD3890
end interface

interface
 SUBROUTINE LD4334(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(4334)
  double precision :: X(4334)
  double precision :: Y(4334)
  double precision :: Z(4334)
 end subroutine LD4334
end interface

interface
 SUBROUTINE LD4802(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(4802)
  double precision :: X(4802)
  double precision :: Y(4802)
  double precision :: Z(4802)
 end subroutine LD4802
end interface

interface
 SUBROUTINE LD5294(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(5294)
  double precision :: X(5294)
  double precision :: Y(5294)
  double precision :: Z(5294)
 end subroutine LD5294
end interface

interface
 SUBROUTINE LD5810(X,Y,Z,W,N)
  implicit none
  integer :: N
  double precision :: W(5810)
  double precision :: X(5810)
  double precision :: Y(5810)
  double precision :: Z(5810)
 end subroutine LD5810
end interface

interface
 SUBROUTINE nfourier(rindata,coutdata,iflag,Iwmax,L,Beta)
  implicit none
  integer :: Iwmax
  integer :: L
  integer :: iflag
  double precision :: Beta
  complex*16 :: coutdata(Iwmax+1)
  double precision :: rindata(L)
 end subroutine nfourier
end interface

interface
 SUBROUTINE invfourier(cindata,routdata,Iwmax,L,iflag,beta)
  implicit none
  integer :: Iwmax
  integer :: L
  integer :: iflag
  double precision :: beta
  complex*16 :: cindata(0:Iwmax)
  complex*16 :: routdata(L)
 end subroutine invfourier
end interface

interface
 SUBROUTINE ludcmp(a,n,np,indx,id,info)
  implicit none
  integer :: id
  integer :: info
  integer :: n
  integer :: np
  real*8 :: a(np,np)
  integer :: indx(n)
 end subroutine ludcmp
end interface

interface
 SUBROUTINE lubksb(a,n,np,indx,b)
  implicit none
  integer :: n
  integer :: np
  real*8 :: a(np,np)
  real*8 :: b(n)
  integer :: indx(n)
 end subroutine lubksb
end interface

interface
 subroutine polyn_coeff(n,x,y,coeff)
  implicit none
  integer :: n
  double precision :: coeff(n)
  double precision :: x(n)
  double precision :: y(n)
 end subroutine polyn_coeff
end interface

interface
 subroutine smooth(a,mesh,it)
  implicit none
  integer, intent(in) :: it
  integer, intent(in) :: mesh
  real*8, intent(inout) :: a(mesh)
 end subroutine smooth
end interface

interface
 subroutine sort_dp(n,list,iperm,tol)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: tol
  integer, intent(inout) :: iperm(n)
  double precision, intent(inout) :: list(n)
 end subroutine sort_dp
end interface

interface
 subroutine sort_int(n,list,iperm)
  implicit none
  integer :: n
  integer :: iperm(n)
  integer :: list(n)
 end subroutine sort_int
end interface


interface
 function uniformrandom(seed) 
  implicit none
  integer :: seed
  double precision :: uniformrandom
 end function uniformrandom
end interface

interface
 subroutine zgedi(a,lda,n,ipvt,det,work,job)
  implicit none
  integer :: job
  integer :: lda
  integer :: n
  complex*16 :: det(2)
  integer :: ipvt(*)
  complex*16 :: work(*)
  complex*16 :: a(lda,*)
 end subroutine zgedi
end interface

interface
 subroutine zgefa(a,lda,n,ipvt,info)
  implicit none
  integer :: info
  integer :: lda
  integer :: n
  integer :: ipvt(*)
  complex*16 :: a(lda,*)
 end subroutine zgefa
end interface

end module interfaces_28_numeric_noabirule
!!***
