!!****m* ABINIT/interfaces_65_nonlocal
!! NAME
!! interfaces_65_nonlocal
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/65_nonlocal
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

module interfaces_65_nonlocal

 implicit none

interface
 subroutine cont13(rank1,rank3,rank2)
  use defs_basis
  implicit none
  real(dp),intent(in) :: rank1(2,3)
  real(dp),intent(out) :: rank2(6)
  real(dp),intent(in) :: rank3(2,10)
 end subroutine cont13
end interface

interface
 subroutine cont22(gxa,gmet,rank2)
  use defs_basis
  implicit none
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gxa(2,6)
  real(dp),intent(out) :: rank2(6)
 end subroutine cont22
end interface

interface
 subroutine cont22cso(gxa1,gxa2,gmet,rank2c)
  use defs_basis
  implicit none
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gxa1(2,6)
  real(dp),intent(in) :: gxa2(2,6)
  real(dp),intent(out) :: rank2c(2,6)
 end subroutine cont22cso
end interface

interface
 subroutine cont22so(gxa1,gxa2,amet,rank2)
  use defs_basis
  implicit none
  real(dp),intent(in) :: amet(2,3,3)
  real(dp),intent(in) :: gxa1(2,6)
  real(dp),intent(in) :: gxa2(2,6)
  real(dp),intent(out) :: rank2(6)
 end subroutine cont22so
end interface

interface
 subroutine cont24(gxa,rank4,rank2)
  use defs_basis
  implicit none
  real(dp),intent(in) :: gxa(2,6)
  real(dp),intent(out) :: rank2(6)
  real(dp),intent(in) :: rank4(2,15)
 end subroutine cont24
end interface

interface
 subroutine cont3(gxa,gmet,rank2)
  use defs_basis
  implicit none
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gxa(2,10)
  real(dp),intent(out) :: rank2(6)
 end subroutine cont3
end interface

interface
 subroutine cont33cso(gxa1,gxa2,gmet,rank2c)
  use defs_basis
  implicit none
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gxa1(2,10)
  real(dp),intent(in) :: gxa2(2,10)
  real(dp),intent(out) :: rank2c(2,6)
 end subroutine cont33cso
end interface

interface
 subroutine cont33so(gxa1,gxa2,gmet,amet,rank2)
  use defs_basis
  implicit none
  real(dp),intent(in) :: amet(2,3,3)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gxa1(2,10)
  real(dp),intent(in) :: gxa2(2,10)
  real(dp),intent(out) :: rank2(6)
 end subroutine cont33so
end interface

interface
 subroutine cont35(gxa,rank5,rank2)
  use defs_basis
  implicit none
  real(dp),intent(in) :: gxa(2,10)
  real(dp),intent(out) :: rank2(6)
  real(dp),intent(in) :: rank5(2,21)
 end subroutine cont35
end interface

interface
 subroutine ctocprj(atindx,cg,choice,cprj,gmet,gprimd,iatom,idir,&  
  &  iorder_cprj,istwfk,kg,kpt,mband,mcg,mcprj,mgfft,mkmem,mpi_enreg,mpsang,&  
  &  mpw,natom,nattyp,nband,ncprj,ngfft,nkpt,nloalg,npwarr,nspinor,&  
  &  nsppol,ntypat,paral_kgb,ph1d,psps,rmet,typat,ucvol,uncp,unkg,&  
  &  unylm,useylmgr,wffnow,xred,ylm,ylmgr)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: iatom
  integer,intent(in) :: idir
  integer,intent(in) :: iorder_cprj
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: ncprj
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: uncp
  integer,intent(in) :: unkg
  integer,intent(in) :: unylm
  integer,intent(in) :: useylmgr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in),target :: atindx(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(cprj_type),intent(out) :: cprj(ncprj,mcprj)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in),target :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in),target :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*useylmgr)
 end subroutine ctocprj
end interface

interface
 subroutine ddkten(compact,idir,rank,temp,tmpfac)
  use defs_basis
  implicit none
  integer,intent(in) :: compact
  integer,intent(in) :: idir
  integer,intent(in) :: rank
  real(dp),intent(inout) :: temp(2,(rank*(rank+1))/2)
  real(dp),intent(inout) :: tmpfac(2,((rank+1)*(rank+2))/2)
 end subroutine ddkten
end interface

interface
 subroutine getcprj(choice,cpopt,cwavef,cwaveprj,dimekb1,dimekb2,dimffnl,ekb,ffnl,&  
  &  idir,indlmn,istwf_k,kg_k,kpg,kpoint,lmnmax,matblk,mgfft,mpi_enreg,&  
  &  natom,nattyp,ngfft,nkpg,nloalg,npw_k,nspinor,ntypat,&  
  &  phkxred,ph1d,ph3d,ucvol,usepaw,useylm)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cpopt
  integer,intent(in) :: dimekb1
  integer,intent(in) :: dimekb2
  integer,intent(in) :: dimffnl
  integer,intent(in) :: idir
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: natom
  integer,intent(in) :: nkpg
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  integer,intent(in) :: useylm
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  real(dp),intent(in) :: cwavef(2,npw_k*nspinor)
  type(cprj_type),intent(out) :: cwaveprj(natom,nspinor)
  real(dp),intent(in) :: ekb(dimekb1,dimekb2)
  real(dp),intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kpg(npw_k,nkpg)
  real(dp),intent(in) :: kpoint(3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(inout) :: ph3d(2,npw_k,matblk)
  real(dp),intent(in) :: phkxred(2,natom)
 end subroutine getcprj
end interface

interface
 subroutine metcon(rank,gmet,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2)
  real(dp),intent(out) :: bb(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: gmet(3,3)
 end subroutine metcon
end interface

interface
 subroutine metcon_so(rank,gmet,amet,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: amet(3,3)
  real(dp),intent(out) :: bb(2,((rank+1)*(rank+2))/2)
  real(dp),intent(in) :: gmet(3,3)
 end subroutine metcon_so
end interface

interface
 subroutine metric_so(amet,gprimd,pauli)
  use defs_basis
  implicit none
  real(dp),intent(out) :: amet(2,3,3,2,2)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: pauli(2,2,2,3)
 end subroutine metric_so
end interface

interface
 subroutine metstr(istr,rank,iterm,gmet,gprimd,aa,bb)
  use defs_basis
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: iterm
  integer,intent(in) :: rank
  real(dp),intent(in) :: aa(2,((rank+3)*(rank+4))/2)
  real(dp),intent(out) :: bb(2,((rank+3)*(rank+4))/2)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine metstr
end interface

interface
 subroutine mkffkg(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,&  
  &  kg_k,kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,&  
  &  npw,ntens,ntypat,parity)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: idir
  integer,intent(in) :: ipw1
  integer,intent(in) :: ispinor
  integer,intent(in) :: itypat
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mblkpw
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: nffkg
  integer,intent(in) :: nffnl
  integer,intent(in) :: nincpw
  integer,intent(in) :: nkpg
  integer,intent(in) :: nlang
  integer,intent(in) :: npw
  integer,intent(in) :: ntens
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: ffkg(mblkpw,nffkg)
  real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(out) :: kpgx(mblkpw,ntens)
  real(dp),intent(in) :: kpt(3)
  integer,intent(out) :: parity(nffkg)
 end subroutine mkffkg
end interface

interface
 subroutine mkffkg3(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,&  
  &  kg_k,kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,&  
  &  npw,ntens,ntypat,parity)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: idir
  integer,intent(in) :: ipw1
  integer,intent(in) :: ispinor
  integer,intent(in) :: itypat
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mblkpw
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: nffkg
  integer,intent(in) :: nffnl
  integer,intent(in) :: nincpw
  integer,intent(in) :: nkpg
  integer,intent(in) :: nlang
  integer,intent(in) :: npw
  integer,intent(in) :: ntens
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: ffkg(nffkg,mblkpw)
  real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(out) :: kpgx(mblkpw,ntens)
  real(dp),intent(in) :: kpt(3)
  integer,intent(out) :: parity(nffkg)
 end subroutine mkffkg3
end interface

interface
 subroutine mkffnl(dimekb,dimffnl,ekb,ffnl,ffspl,gmet,gprimd,ider,idir,indlmn,&  
  &  kg,kpg,kpt,lmnmax,lnmax,mpsang,mqgrid,nkpg,npw,ntypat,pspso,&  
  &  qgrid,rmet,usepaw,useylm,ylm,ylm_gr)
  use defs_basis
  implicit none
  integer,intent(in) :: dimekb
  integer,intent(in) :: dimffnl
  integer,intent(in) :: ider
  integer,intent(in) :: idir
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: nkpg
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  integer,intent(in) :: useylm
  real(dp),intent(in) :: ekb(dimekb,ntypat*(1-usepaw))
  real(dp),intent(out) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  real(dp),intent(in) :: ffspl(mqgrid,2,lnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpg(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  integer,intent(in) :: pspso(ntypat)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: ylm(npw,mpsang*mpsang*useylm)
  real(dp),intent(in) :: ylm_gr(npw,3+6*(ider/2),mpsang*mpsang*useylm)
 end subroutine mkffnl
end interface

interface
 subroutine mkkpg(kg,kpg,kpt,nkpg,npw)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpg
  integer,intent(in) :: npw
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(out) :: kpg(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
 end subroutine mkkpg
end interface

interface
 subroutine nonlop(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimffnlin,dimffnlout,&  
  &  enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&  
  &  kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&  
  &  mpi_enreg,mpsang,mpssoang,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,&  
  &  nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&  
  &  phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&  
  &  tim_nonlop,ucvol,useylm,vectin,vectout,&  
  &  use_gpu_cuda) !optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cpopt
  integer,intent(in) :: dimenl1
  integer,intent(in) :: dimenl2
  integer,intent(in) :: dimffnlin
  integer,intent(in) :: dimffnlout
  integer,intent(in) :: idir
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: natom
  integer,intent(in) :: nkpgin
  integer,intent(in) :: nkpgout
  integer,intent(in) :: nnlout
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: nspinor
  integer,intent(in) :: nspinortot
  integer,intent(in) :: ntypat
  integer,intent(in) :: only_SO
  integer,intent(in) :: paw_opt
  integer,intent(in) :: signs
  integer,intent(in) :: tim_nonlop
  integer,optional,intent(in) :: use_gpu_cuda
  integer,intent(in) :: useylm
  real(dp),intent(in) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5))
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
  real(dp),intent(out) :: enlout(nnlout)
  real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kgin(3,npwin)
  integer,intent(in) :: kgout(3,npwout)
  real(dp),intent(in) :: kpgin(npwin,nkpgin*useylm)
  real(dp),intent(in) :: kpgout(npwout,nkpgout*useylm)
  real(dp),intent(in) :: kptin(3)
  real(dp),intent(in) :: kptout(3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(inout) :: ph3din(2,npwin,matblk)
  real(dp),intent(inout) :: ph3dout(2,npwout,matblk)
  real(dp),intent(in) :: phkxredin(2,natom)
  real(dp),intent(in) :: phkxredout(2,natom)
  real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp),intent(out) :: svectout(2,npwout*nspinor*(paw_opt/3))
  real(dp),intent(inout) :: vectin(2,npwin*nspinor)
  real(dp),intent(out) :: vectout(2,npwout*nspinor)
 end subroutine nonlop
end interface

interface
 subroutine nonlop_gpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimffnlin,dimffnlout,&  
  &  enl,enlout,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&  
  &  kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&  
  &  mpi_enreg,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&  
  &  npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&  
  &  ph3din,ph3dout,signs,sij,svectout,ucvol,vectin,vectout)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cpopt
  integer,intent(in) :: dimenl1
  integer,intent(in) :: dimenl2
  integer,intent(in) :: dimffnlin
  integer,intent(in) :: dimffnlout
  integer,intent(in) :: idir
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: natom
  integer,intent(in) :: nkpgin
  integer,intent(in) :: nkpgout
  integer,intent(in) :: nnlout
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: nspinor
  integer,intent(in) :: nspinortot
  integer,intent(in) :: ntypat
  integer,intent(in) :: paw_opt
  integer,intent(in) :: signs
  real(dp),intent(in) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5))
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
  real(dp),intent(out) :: enlout(nnlout)
  real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kgin(3,npwin)
  integer,intent(in) :: kgout(3,npwout)
  real(dp),intent(in) :: kpgin(npwin,nkpgin)
  real(dp),intent(in) :: kpgout(npwout,nkpgout)
  real(dp),intent(in) :: kptin(3)
  real(dp),intent(in) :: kptout(3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(inout) :: ph3din(2,npwin,matblk)
  real(dp),intent(inout) :: ph3dout(2,npwout,matblk)
  real(dp),intent(in) :: phkxredin(2,natom)
  real(dp),intent(in) :: phkxredout(2,natom)
  real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp),intent(out),target :: svectout(2,npwout*nspinor*(paw_opt/3))
  real(dp),intent(inout) :: vectin(2,npwin*nspinor)
  real(dp),intent(out),target :: vectout(2,npwout*nspinor)
 end subroutine nonlop_gpu
end interface

interface
 subroutine nonlop_pl(choice,dimekb1,dimekb2,dimffnlin,dimffnlout,ekb,enlout,&  
  &  ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,kgin,kgout,kpgin,kpgout,&  
  &  kptin,kptout,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,&  
  &  natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,npwin,npwout,nspinor,nspinortot,&  
  &  ntypat,only_SO,phkxredin,phkxredout,ph1d,ph3din,ph3dout,signs,&  
  &  ucvol,vectin,vectout)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: dimekb1
  integer,intent(in) :: dimekb2
  integer,intent(in) :: dimffnlin
  integer,intent(in) :: dimffnlout
  integer,intent(in) :: idir
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: natom
  integer,intent(in) :: nkpgin
  integer,intent(in) :: nkpgout
  integer,intent(in) :: nnlout
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: nspinor
  integer,intent(in) :: nspinortot
  integer,intent(in) :: ntypat
  integer,intent(in) :: only_SO
  integer,intent(in) :: signs
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  real(dp),intent(in) :: ekb(dimekb1,dimekb2,nspinortot**2)
  real(dp),intent(out) :: enlout(nnlout)
  real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kgin(3,npwin)
  integer,intent(in) :: kgout(3,npwout)
  real(dp),intent(in) :: kpgin(npwin,nkpgin)
  real(dp),intent(in) :: kpgout(npwout,nkpgout)
  real(dp),intent(in) :: kptin(3)
  real(dp),intent(in) :: kptout(3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(inout) :: ph3din(2,npwin,matblk)
  real(dp),intent(inout) :: ph3dout(2,npwout,matblk)
  real(dp),intent(in) :: phkxredin(2,natom)
  real(dp),intent(in) :: phkxredout(2,natom)
  real(dp),intent(inout) :: vectin(2,nspinor*npwin)
  real(dp),intent(out) :: vectout(2,nspinor*npwout)
 end subroutine nonlop_pl
end interface

interface
 subroutine nonlop_ylm(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimffnlin,dimffnlout,&  
  &  enl,enlout,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&  
  &  kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&  
  &  mpi_enreg,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&  
  &  npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&  
  &  ph3din,ph3dout,signs,sij,svectout,ucvol,vectin,vectout,cprjin_left)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cpopt
  integer,intent(in) :: dimenl1
  integer,intent(in) :: dimenl2
  integer,intent(in) :: dimffnlin
  integer,intent(in) :: dimffnlout
  integer,intent(in) :: idir
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: natom
  integer,intent(in) :: nkpgin
  integer,intent(in) :: nkpgout
  integer,intent(in) :: nnlout
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: nspinor
  integer,intent(in) :: nspinortot
  integer,intent(in) :: ntypat
  integer,intent(in) :: paw_opt
  integer,intent(in) :: signs
  real(dp),intent(in) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5))
  type(cprj_type),optional,intent(in) :: cprjin_left(natom,nspinor)
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
  real(dp),intent(out) :: enlout(nnlout)
  real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kgin(3,npwin)
  integer,intent(in) :: kgout(3,npwout)
  real(dp),intent(in) :: kpgin(npwin,nkpgin)
  real(dp),intent(in) :: kpgout(npwout,nkpgout)
  real(dp),intent(in) :: kptin(3)
  real(dp),intent(in) :: kptout(3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(inout) :: ph3din(2,npwin,matblk)
  real(dp),intent(inout) :: ph3dout(2,npwout,matblk)
  real(dp),intent(in) :: phkxredin(2,natom)
  real(dp),intent(in) :: phkxredout(2,natom)
  real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp),intent(out) :: svectout(2,npwout*nspinor*(paw_opt/3))
  real(dp),intent(inout) :: vectin(2,npwin*nspinor)
  real(dp),intent(out) :: vectout(2,npwout*nspinor)
 end subroutine nonlop_ylm
end interface

interface
 subroutine opernl2(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&  
  &  ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat,&  
  &  jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&  
  &  mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw,&  
  &  ntypat,ph3d,sign,vect)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: ia3
  integer,intent(in) :: idir
  integer,intent(in) :: ispinor
  integer,intent(in) :: istwf_k
  integer,intent(in) :: itypat
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mincat
  integer,intent(in) :: mlang1
  integer,intent(in) :: mlang3
  integer,intent(in) :: mlang4
  integer,intent(in) :: mlang5
  integer,intent(in) :: mlang6
  integer,intent(in) :: mproj
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: nffnl
  integer,intent(in) :: nincat
  integer,intent(in) :: nkpg
  integer,intent(in) :: nlang
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  integer,intent(in) :: sign
  integer,intent(in) :: nloalg(5)
  real(dp),intent(out) :: d2gxdis(2,mlang5,mincat,mproj)
  real(dp),intent(out) :: d2gxds2(2,mlang6,mincat,mproj)
  real(dp),intent(out) :: dgxdis(2,mlang1,mincat,mproj)
  real(dp),intent(inout) :: dgxds(2,mlang4,mincat,mproj)
  real(dp),intent(inout) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj)
  real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(inout) :: gxa(2,mlang3,mincat,mproj)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: jproj(nlang)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: ph3d(2,npw,matblk)
  real(dp),intent(inout) :: vect(2,npw)
 end subroutine opernl2
end interface

interface
 subroutine opernl3(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&  
  &  ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat,&  
  &  jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&  
  &  mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw,&  
  &  ntypat,ph3d,sign,vect)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: ia3
  integer,intent(in) :: idir
  integer,intent(in) :: ispinor
  integer,intent(in) :: istwf_k
  integer,intent(in) :: itypat
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mincat
  integer,intent(in) :: mlang1
  integer,intent(in) :: mlang3
  integer,intent(in) :: mlang4
  integer,intent(in) :: mlang5
  integer,intent(in) :: mlang6
  integer,intent(in) :: mproj
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: nffnl
  integer,intent(in) :: nincat
  integer,intent(in) :: nkpg
  integer,intent(in) :: nlang
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  integer,intent(in) :: sign
  integer,intent(in) :: nloalg(5)
  real(dp),intent(out) :: d2gxdis(2,mlang5,mincat,mproj)
  real(dp),intent(out) :: d2gxds2(2,mlang6,mincat,mproj)
  real(dp),intent(out) :: dgxdis(2,mlang1,mincat,mproj)
  real(dp),intent(inout) :: dgxds(2,mlang4,mincat,mproj)
  real(dp),intent(inout) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj)
  real(dp),intent(in) :: ffnl(1,npw,nffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(inout) :: gxa(2,mlang3,mincat,mproj)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: jproj(nlang)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: ph3d(2,npw,matblk)
  real(dp),intent(inout) :: vect(2,npw)
 end subroutine opernl3
end interface

interface
 subroutine opernl4a(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&  
  &  ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat,&  
  &  jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&  
  &  mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw,&  
  &  ntypat,ph3d,vect)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: ia3
  integer,intent(in) :: idir
  integer,intent(in) :: ispinor
  integer,intent(in) :: istwf_k
  integer,intent(in) :: itypat
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mincat
  integer,intent(in) :: mlang1
  integer,intent(in) :: mlang3
  integer,intent(in) :: mlang4
  integer,intent(in) :: mlang5
  integer,intent(in) :: mlang6
  integer,intent(in) :: mproj
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: nffnl
  integer,intent(in) :: nincat
  integer,intent(in) :: nkpg
  integer,intent(in) :: nlang
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  integer,intent(in) :: nloalg(5)
  real(dp),intent(out) :: d2gxdis(2,mlang5,mincat,mproj)
  real(dp),intent(out) :: d2gxds2(2,mlang6,mincat,mproj)
  real(dp),intent(out) :: dgxdis(2,mlang1,mincat,mproj)
  real(dp),intent(out) :: dgxds(2,mlang4,mincat,mproj)
  real(dp),intent(out) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj)
  real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(out) :: gxa(2,mlang3,mincat,mproj)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: jproj(nlang)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: ph3d(2,npw,matblk)
  real(dp),intent(in) :: vect(2,npw)
 end subroutine opernl4a
end interface

interface
 subroutine opernl4b(choice,dgxds,dgxdt,ffnl,gmet,gxa,&  
  &  ia3,idir,indlmn,ispinor,itypat,jproj,kg_k,kpg_k,kpt,&  
  &  lmnmax,matblk,mincat,mlang3,mlang4,mproj,ndgxdt,nffnl,nincat,&  
  &  nkpg,nlang,nloalg,npw,ntypat,ph3d,vect)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: ia3
  integer,intent(in) :: idir
  integer,intent(in) :: ispinor
  integer,intent(in) :: itypat
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mincat
  integer,intent(in) :: mlang3
  integer,intent(in) :: mlang4
  integer,intent(in) :: mproj
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: nffnl
  integer,intent(in) :: nincat
  integer,intent(in) :: nkpg
  integer,intent(in) :: nlang
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  integer,intent(in) :: nloalg(5)
  real(dp),intent(in) :: dgxds(2,mlang4,mincat,mproj)
  real(dp),intent(in) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj)
  real(dp),intent(in) :: ffnl(1,npw,nffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gxa(2,mlang3,mincat,mproj)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: jproj(nlang)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: ph3d(2,npw,matblk)
  real(dp),intent(out) :: vect(2,npw)
 end subroutine opernl4b
end interface

interface
 subroutine opernla_ylm(choice,cplex,dimffnl,d2gxdt,dgxdt,ffnl,gx,ia3,idir,&  
  &  indlmn,istwf_k,kpg,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpg,nlmn,&  
  &  nloalg,npw,nspinor,ph3d,signs,ucvol,vect)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cplex
  integer,intent(in) :: dimffnl
  integer,intent(in) :: ia3
  integer,intent(in) :: idir
  integer,intent(in) :: istwf_k
  integer,intent(in) :: matblk
  integer,intent(in) :: nd2gxdt
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: nincat
  integer,intent(in) :: nkpg
  integer,intent(in) :: nlmn
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: signs
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: nloalg(5)
  real(dp),intent(out) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
  real(dp),intent(out) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
  real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn)
  real(dp),intent(out) :: gx(cplex,nlmn,nincat,nspinor)
  integer,intent(in) :: indlmn(6,nlmn)
  real(dp),intent(in) :: kpg(npw,nkpg)
  real(dp),intent(in) :: ph3d(2,npw,matblk)
  real(dp),intent(in) :: vect(2,npw*nspinor)
 end subroutine opernla_ylm
end interface

interface
 subroutine opernlb_ylm(choice,cplex,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,&  
  &  ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,&  
  &  nspinor,paw_opt,ph3d,svect,ucvol,vect)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cplex
  integer,intent(in) :: cplex_fac
  integer,intent(in) :: dimffnl
  integer,intent(in) :: ia3
  integer,intent(in) :: idir
  integer,intent(in) :: matblk
  integer,intent(in) :: ndgxdtfac
  integer,intent(in) :: nincat
  integer,intent(in) :: nkpg
  integer,intent(in) :: nlmn
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: paw_opt
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: nloalg(5)
  real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
  real(dp),intent(in) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat*(paw_opt/3),nspinor)
  real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn)
  real(dp),intent(in) :: gxfac(cplex_fac,nlmn,nincat,nspinor)
  real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
  integer,intent(in) :: indlmn(6,nlmn)
  real(dp),intent(in) :: kpg(npw,nkpg)
  real(dp),intent(in) :: ph3d(2,npw,matblk)
  real(dp),intent(inout) :: svect(2,npw*nspinor*(paw_opt/3))
  real(dp),intent(inout) :: vect(2,npw*nspinor)
 end subroutine opernlb_ylm
end interface

interface
 subroutine opernlc_ylm(atindx1,cplex,cplex_enl,cplex_fac,dgxdt,dgxdtfac,dgxdtfac_sij,&  
  &  dimenl1,dimenl2,enl,gx,gxfac,gxfac_sij,iatm,indlmn,itypat,&  
  &  lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,nincat,nlmn,&  
  &  nspinor,nspinortot,optder,paw_opt,sij)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: cplex_enl
  integer,intent(in) :: cplex_fac
  integer,intent(in) :: dimenl1
  integer,intent(in) :: dimenl2
  integer,intent(in) :: iatm
  integer,intent(in) :: itypat
  integer,intent(in) :: natom
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: ndgxdtfac
  integer,intent(in) :: nincat
  integer,intent(inout) :: nlmn
  integer,intent(in) :: nspinor
  integer,intent(in) :: nspinortot
  integer,intent(in) :: optder
  integer,intent(in) :: paw_opt
  real(dp) :: lambda
  type(mpi_type) , intent(in) :: mpi_enreg
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
  real(dp),intent(out) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
  real(dp),intent(out) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
  real(dp),intent(inout) :: gx(cplex,nlmn,nincat,nspinor)
  real(dp),intent(out) :: gxfac(cplex_fac,nlmn,nincat,nspinor)
  real(dp),intent(out) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
  integer,intent(in) :: indlmn(6,nlmn)
  real(dp),intent(in) :: sij(((paw_opt+1)/3)*nlmn*(nlmn+1)/2)
 end subroutine opernlc_ylm
end interface

interface
 subroutine opernld_ylm(choice,cplex,cplex_fac,d2gxdt,dgxdt,dgxdtfac,enlk,enlout,fnlk,&  
  &  gx,gxfac,gxfac_sij,ia3,natom,nd2gxdt,ndgxdt,ndgxdtfac,nincat,&  
  &  nlmn,nnlout,nspinor,paw_opt,strnlk,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cplex
  integer,intent(in) :: cplex_fac
  integer,intent(in) :: ia3
  integer,intent(in) :: natom
  integer,intent(in) :: nd2gxdt
  integer,intent(in) :: ndgxdt
  integer,intent(in) :: ndgxdtfac
  integer,intent(in) :: nincat
  integer,intent(in) :: nlmn
  integer,intent(in) :: nnlout
  integer,intent(in) :: nspinor
  integer,intent(in) :: paw_opt
  real(dp),intent(out) :: enlk
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
  real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
  real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
  real(dp),intent(inout) :: enlout(nnlout)
  real(dp),intent(inout) :: fnlk(3*natom)
  real(dp),intent(in) :: gx(cplex,nlmn,nincat,nspinor)
  real(dp),intent(in) :: gxfac(cplex_fac,nlmn,nincat,nspinor)
  real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
  real(dp),intent(inout) :: strnlk(6)
 end subroutine opernld_ylm
end interface

interface
 subroutine ph1d3d(iatom,jatom,kg_k,matblk,natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)
  use defs_basis
  implicit none
  integer,intent(in) :: iatom
  integer,intent(in) :: jatom
  integer,intent(in) :: matblk
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: npw_k
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
  real(dp),intent(out) :: ph3d(2,npw_k,matblk)
  real(dp),intent(in) :: phkxred(2,natom)
 end subroutine ph1d3d
end interface

interface
 subroutine scalewf_nonlop(istwf_k,mpi_enreg,npw,option,vect)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: option
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(inout) :: vect(2,npw)
 end subroutine scalewf_nonlop
end interface

interface
 subroutine strsocv(red,gprimd,cart)
  use defs_basis
  implicit none
  real(dp),intent(out) :: cart(6)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: red(6,3)
 end subroutine strsocv
end interface

interface
 subroutine trace2(gxa,gmet,trace)
  use defs_basis
  implicit none
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gxa(2,6)
  real(dp),intent(out) :: trace(2)
 end subroutine trace2
end interface

end module interfaces_65_nonlocal
!!***
