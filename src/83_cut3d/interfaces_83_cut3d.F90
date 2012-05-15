!!****m* ABINIT/interfaces_83_cut3d
!! NAME
!! interfaces_83_cut3d
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/83_cut3d
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

module interfaces_83_cut3d

 implicit none

interface
 subroutine hirsh(grid_den,natom,nrx,nry,nrz,ntypat,rprimd,xcart,typat,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrx
  integer,intent(in) :: nry
  integer,intent(in) :: nrz
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: grid_den(nrx,nry,nrz)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xcart(3,natom)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine hirsh
end interface

interface
 subroutine kptindex(shiftx,shifty,shiftz,ikpt_x,ikpt_y,ikpt_z,nkx,nky,nkz,&  
  &  ikpt_x3,ikpt_y3,ikpt_z3)
  use defs_basis
  implicit none
  integer,intent(in) :: ikpt_x
  integer,intent(out) :: ikpt_x3
  integer,intent(in) :: ikpt_y
  integer,intent(out) :: ikpt_y3
  integer,intent(in) :: ikpt_z
  integer,intent(out) :: ikpt_z3
  integer,intent(in) :: nkx
  integer,intent(in) :: nky
  integer,intent(in) :: nkz
  real(dp),intent(in) :: shiftx
  real(dp),intent(in) :: shifty
  real(dp),intent(in) :: shiftz
 end subroutine kptindex
end interface

interface
 subroutine lineint(gridtt,gridux,griddy,gridmz,nr1,nr2,nr3,nspden,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspden
  real(dp),intent(in) :: griddy(nr1,nr2,nr3)
  real(dp),intent(in) :: gridmz(nr1,nr2,nr3)
  real(dp),intent(in) :: gridtt(nr1,nr2,nr3)
  real(dp),intent(in) :: gridux(nr1,nr2,nr3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine lineint
end interface

interface
 subroutine localorb_s(chr_inputfname,ecut,exchn2n3d,headform,istwfk,kpt,natom,nband,nkpt,&  
  npwarr,nr1,nr2,nr3,nspinor,nsppol,ntypat,paral_kgb,rprimd,tau,typat,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: headform
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  character(len=50),intent(in) :: chr_inputfname
  real(dp),intent(in) :: ecut
  integer,intent(in) :: istwfk(nkpt)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: tau(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine localorb_s
end interface

interface
 subroutine normalize(v)
  use defs_basis
  implicit none
  real(dp),intent(inout) :: v(3)
 end subroutine normalize
end interface

interface
 subroutine overlap_wf(cwave1,e_kpt,exchn2n3d,csppol,cbandpick,ckpt,&  
  &  ecut,headform,istwfk,kpt,nband,nbands,nkpt,npwarr,&  
  &  nr1,nr2,nr3,nspinor,nsppol,paral_kgb,rprimd,coverlap)
  use defs_basis
  implicit none
  integer,intent(in) :: cbandpick
  integer,intent(in) :: ckpt
  integer,intent(in) :: csppol
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: headform
  integer,intent(in) :: nbands
  integer,intent(in) :: nkpt
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: paral_kgb
  complex(dp),intent(out) :: coverlap
  real(dp),intent(in) :: ecut
  complex(dp),intent(in) :: cwave1(nr1,nr2,nr3)
  real(dp),intent(out) :: e_kpt(nbands)
  integer,intent(in) :: istwfk(nkpt)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine overlap_wf
end interface

interface
 subroutine planeint(gridtt,gridux,griddy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,tau)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspden
  real(dp),intent(in) :: griddy(nr1,nr2,nr3)
  real(dp),intent(in) :: gridmz(nr1,nr2,nr3)
  real(dp),intent(in) :: gridtt(nr1,nr2,nr3)
  real(dp),intent(in) :: gridux(nr1,nr2,nr3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: tau(3,natom)
 end subroutine planeint
end interface

interface
 subroutine pointint(gridt,gridu,gridd,gridm,nr1,nr2,nr3,nspden,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspden
  real(dp),intent(in) :: gridd(nr1,nr2,nr3)
  real(dp),intent(in) :: gridm(nr1,nr2,nr3)
  real(dp),intent(in) :: gridt(nr1,nr2,nr3)
  real(dp),intent(in) :: gridu(nr1,nr2,nr3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine pointint
end interface

interface
 subroutine recip(x1,hkl,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: hkl(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: x1(3)
 end subroutine recip
end interface

interface
 subroutine reduce(r,rcart,rprimd)
  use defs_basis
  implicit none
  real(dp),intent(out) :: r(3)
  real(dp),intent(in) :: rcart(3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine reduce
end interface

interface
 subroutine rrho(densfileformat,grid_full,nr1,nr2,nr3,nspden,wff)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(in) :: densfileformat
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspden
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out),target :: grid_full(nr1,nr2,nr3,nspden)
 end subroutine rrho
end interface

interface
 subroutine rtau(filnam,tau,nat,ntypat)
  use defs_basis
  implicit none
  integer,intent(in) :: nat
  integer,intent(in) :: ntypat
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(out) :: tau(3,nat)
 end subroutine rtau
end interface

interface
 subroutine vdot(x1,x2,x3)
  use defs_basis
  implicit none
  real(dp),intent(in) :: x1(3)
  real(dp),intent(in) :: x2(3)
  real(dp),intent(out) :: x3(3)
 end subroutine vdot
end interface

interface
 subroutine volumeint(gridtt,gridux,griddy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,tau)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspden
  real(dp),intent(in) :: griddy(nr1,nr2,nr3)
  real(dp),intent(in) :: gridmz(nr1,nr2,nr3)
  real(dp),intent(in) :: gridtt(nr1,nr2,nr3)
  real(dp),intent(in) :: gridux(nr1,nr2,nr3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: tau(3,natom)
 end subroutine volumeint
end interface

interface
 subroutine wffile(ecut,exchn2n3d,headform,istwfk,kpt,natom,nband,nkpt,npwarr,&  
  &  nr1,nr2,nr3,nspinor,nsppol,ntypat,paral_kgb,rprimd,tau,typat,wff,znucl)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: headform
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: ecut
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: istwfk(nkpt)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: tau(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine wffile
end interface

interface
 subroutine wfread(cwave0,e_kpt,exchn2n3d,csppol,cbandpick,ckpt,&  
  &  ecut,headform,istwfk,kpt,nband,nbands,nkpt,npwarr,&  
  &  nr1,nr2,nr3,nspinor,nsppol,paral_kgb,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: cbandpick
  integer,intent(in) :: ckpt
  integer,intent(in) :: csppol
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: headform
  integer,intent(in) :: nbands
  integer,intent(in) :: nkpt
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: ecut
  complex(dp),intent(out) :: cwave0(nr1,nr2,nr3)
  real(dp),intent(out) :: e_kpt(nbands)
  integer,intent(in) :: istwfk(nkpt)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine wfread
end interface

end module interfaces_83_cut3d
!!***
