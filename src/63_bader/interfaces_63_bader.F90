!!****m* ABINIT/interfaces_63_bader
!! NAME
!! interfaces_63_bader
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/63_bader
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

module interfaces_63_bader

 implicit none

interface
 subroutine addout(aim_dtset)
  use defs_abitypes
  implicit none
  type(aim_dataset_type),intent(in) :: aim_dtset
 end subroutine addout
end interface

interface
 subroutine adini(aim_dtset,inpstr,lenstr)
  use defs_abitypes
  implicit none
  integer,intent(in) :: lenstr
  type(aim_dataset_type), intent(out) :: aim_dtset
  character(len=*),intent(in) :: inpstr
 end subroutine adini
end interface

interface
 subroutine aim_follow(aim_dtset,vv,npmax,srch,iatinit,iposinit,iat,ipos,nstep)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(out) :: iat
  integer,intent(in) :: iatinit
  integer,intent(out) :: ipos
  integer,intent(in) :: iposinit
  integer,intent(in) :: npmax
  integer,intent(out) :: nstep
  type(aim_dataset_type),intent(in) :: aim_dtset
  logical,intent(inout) :: srch
  real(dp),intent(inout) :: vv(3)
 end subroutine aim_follow
end interface

interface
 subroutine consist(aim_dtset,tstngr,tstvpt)
  use defs_abitypes
  implicit none
  integer,intent(in) :: tstngr
  integer,intent(in) :: tstvpt
  type(aim_dataset_type),intent(in) :: aim_dtset
 end subroutine consist
end interface

interface
 subroutine cpdrv(aim_dtset)
  use defs_abitypes
  implicit none
  type(aim_dataset_type),intent(in) :: aim_dtset
 end subroutine cpdrv
end interface

interface
 subroutine critic(aim_dtset,vv,ev,zz,dmax,ires,sort)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(out) :: ires
  integer,intent(in) :: sort
  type(aim_dataset_type), intent(in) :: aim_dtset
  real(dp),intent(in) :: dmax
  real(dp),intent(out) :: ev(3)
  real(dp),intent(inout) :: vv(3)
  real(dp),intent(out) :: zz(3,3)
 end subroutine critic
end interface

interface
 subroutine ordr(aa,dd,nn,cff)
  use defs_basis
  implicit none
  integer,intent(in) :: cff
  integer,intent(in) :: nn
  real(dp),intent(inout) :: aa(nn)
  real(dp),intent(inout) :: dd(nn,nn)
 end subroutine ordr
end interface

interface
 subroutine  critics(aim_dtset,inxat,stwo,sthree,sfour,dstmax)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: inxat
  type(aim_dataset_type), intent(in) :: aim_dtset
  real(dp),intent(in) :: dstmax
  logical,intent(in) :: sfour
  logical,intent(in) :: sthree
  logical,intent(in) :: stwo
 end subroutine critics
end interface

interface
 subroutine defad(aim_dtset)
  use defs_abitypes
  implicit none
  type(aim_dataset_type),intent(out) :: aim_dtset
 end subroutine defad
end interface

interface
 subroutine drvaim(aim_dtset)
  use defs_abitypes
  implicit none
  type(aim_dataset_type),intent(in) :: aim_dtset
 end subroutine drvaim
end interface

interface
 subroutine evspln(aa,bb,indx,dxx,ndim,fld,sdfd,kod,val,der,dder)
  use defs_basis
  implicit none
  integer,intent(in) :: indx
  integer,intent(in) :: kod
  integer,intent(in) :: ndim
  real(dp),intent(in) :: aa
  real(dp),intent(in) :: bb
  real(dp),intent(out) :: dder
  real(dp),intent(out) :: der
  real(dp),intent(in) :: dxx
  real(dp),intent(out) :: val
  real(dp),intent(in) :: fld(ndim)
  real(dp),intent(in) :: sdfd(ndim)
 end subroutine evspln
end interface

interface
 subroutine graph(unts,untg)
  implicit none
  integer,intent(in) :: untg
  integer,intent(in) :: unts
 end subroutine graph
end interface

interface
 subroutine initaim(aim_dtset)
  use defs_abitypes
  implicit none
  type(aim_dataset_type),intent(in) :: aim_dtset
 end subroutine initaim
end interface

interface
 subroutine inpar(instr,lenstr)
  implicit none
  integer,intent(out) :: lenstr
  character(len=*),intent(out) :: instr
 end subroutine inpar
end interface

interface
 subroutine inspln(idir,snn,tnn)
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: snn
  integer,intent(in) :: tnn
 end subroutine inspln
end interface

interface
 subroutine integrho(aim_dtset)
  use defs_abitypes
  implicit none
  type(aim_dataset_type),intent(in) :: aim_dtset
 end subroutine integrho
end interface

interface
 subroutine integvol()
  implicit none
 end subroutine integvol
end interface

interface
 subroutine onestep(vv,chg,grho,hh,np,npmax,deltar)
  use defs_basis
  implicit none
  integer,intent(out) :: np
  integer,intent(in) :: npmax
  real(dp),intent(out) :: chg
  real(dp),intent(out) :: deltar
  real(dp),intent(in) :: hh
  real(dp),intent(out) :: grho(3)
  real(dp),intent(inout) :: vv(3)
 end subroutine onestep
end interface

interface
 subroutine plint()
  implicit none
 end subroutine plint
end interface

interface
 subroutine rsurf(aim_dtset,rr,grho,theta,phi,rr0,iatinit,npmax,srch)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iatinit
  integer,intent(in) :: npmax
  type(aim_dataset_type),intent(in) :: aim_dtset
  real(dp),intent(in) :: phi
  real(dp),intent(out) :: rr
  real(dp),intent(in) :: rr0
  logical,intent(in) :: srch
  real(dp),intent(in) :: theta
  real(dp),intent(out) :: grho(3)
 end subroutine rsurf
end interface

interface
 subroutine surf(aim_dtset)
  use defs_abitypes
  implicit none
  type(aim_dataset_type) :: aim_dtset
 end subroutine surf
end interface

interface
 subroutine vgh_rho(vv,rho,grho,hrho,rdmin,iat,ipos,chs)
  use defs_basis
  implicit none
  integer,intent(in) :: chs
  integer,intent(inout) :: iat
  integer,intent(inout) :: ipos
  real(dp),intent(out) :: rdmin
  real(dp),intent(out) :: rho
  real(dp),intent(out) :: grho(3)
  real(dp),intent(out) :: hrho(3,3)
  real(dp),intent(in) :: vv(3)
 end subroutine vgh_rho
end interface

interface
 function vnorm(vv,dir)
  use defs_basis
  implicit none
  integer,intent(in) :: dir
  real(dp) :: vnorm
  real(dp),intent(in) :: vv(3)
 end function vnorm
end interface

interface
 function vec_prod(uu,vv)
  use defs_basis
  implicit none
  real(dp),intent(in) :: uu(3)
  real(dp) :: vec_prod(3)
  real(dp),intent(in) :: vv(3)
 end function vec_prod
end interface

interface
 subroutine mprod(aa,bb,cc)
  use defs_basis
  implicit none
  real(dp),intent(in) :: aa(3,3)
  real(dp),intent(in) :: bb(3,3)
  real(dp),intent(out) :: cc(3,3)
 end subroutine mprod
end interface

interface
 subroutine bschg1(vv,dir)
  use defs_basis
  implicit none
  integer,intent(in) :: dir
  real(dp),intent(inout) :: vv(3)
 end subroutine bschg1
end interface

interface
 subroutine bschg2(aa,dir)
  use defs_basis
  implicit none
  integer,intent(in) :: dir
  real(dp),intent(inout) :: aa(3,3)
 end subroutine bschg2
end interface

end module interfaces_63_bader
!!***
