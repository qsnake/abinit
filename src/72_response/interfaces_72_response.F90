!!****m* ABINIT/interfaces_72_response
!! NAME
!! interfaces_72_response
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/72_response
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

module interfaces_72_response

 implicit none

interface
 subroutine accrho3(counter,cplex,cwave0,cwave1,cwavef,cwaveprj0,cwaveprj1,&  
  &  dimcprj,dimffnlk,ffnlk,dimphkxred,eloc0_k,filstat,gbound,&  
  &  gs_hamkq,iband,idir,ipert,isppol,kg_k,kg1_k,kpg_k,kpt,kptopt,lmnmax,matblk,&  
  &  mgfft,mpi_enreg,natom,nband_k,ncpgr,nkpg,npw_k,npw1_k,nspinor,ntypat,n4,n5,n6,&  
  &  occ_k,option,paral_kgb,pawrhoij1,ph3d,phkxred,prtvol,rhoaug1,tim_fourwf,&  
  &  usecprj,vlocal,wf_corrected,wtk_k)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: counter
  integer,intent(in) :: cplex
  integer,intent(in) :: dimffnlk
  integer,intent(in) :: dimphkxred
  integer,intent(in) :: iband
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: kptopt
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nkpg
  integer,intent(in) :: npw1_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: tim_fourwf
  integer,intent(in) :: usecprj
  integer,intent(in) :: wf_corrected
  real(dp),intent(out) :: eloc0_k
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(in) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: wtk_k
  real(dp),intent(in),target :: cwave0(2,npw_k*nspinor)
  real(dp),intent(in),target :: cwave1(2,npw1_k*nspinor)
  real(dp),intent(in),target :: cwavef(2,npw1_k*nspinor)
  type(cprj_type),intent(in) :: cwaveprj0(natom,nspinor*usecprj)
  type(cprj_type),intent(in) :: cwaveprj1(natom,nspinor*gs_hamkq%usepaw)
  integer,intent(in) :: dimcprj(natom*gs_hamkq%usepaw)
  real(dp),intent(in) :: ffnlk(npw_k,dimffnlk,lmnmax,1+gs_hamkq%usepaw*(ntypat-1))
  integer,intent(in) :: gbound(2*mgfft+8,2)
  integer,intent(in) :: kg1_k(3,npw1_k)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kpg_k(npw_k,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: occ_k(nband_k)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(natom*gs_hamkq%usepaw)
  real(dp),intent(inout) :: ph3d(2,npw1_k,matblk)
  real(dp),intent(in) :: phkxred(2,dimphkxred)
  real(dp),intent(inout) :: rhoaug1(cplex*n4,n5,n6)
  real(dp),intent(in) :: vlocal(n4,n5,n6*(option/2))
 end subroutine accrho3
end interface

interface
 subroutine asria_calc(asr,d2asr,d2cart,mpert,natom)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  real(dp),intent(out) :: d2asr(2,3,natom,3,natom)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
 end subroutine asria_calc
end interface

interface
 subroutine asria_corr(asr,d2asr,d2cart,mpert,natom)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  real(dp),intent(in) :: d2asr(2,3,natom,3,natom)
  real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)
 end subroutine asria_corr
end interface

interface
 subroutine asrprs(asr,asrflag,rotinv,uinvers,vtinvers,singular,d2cart,mpert,natom,xcart)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: asrflag
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: rotinv
  real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: singular(1:3*natom*(3*natom-1)/2)
  real(dp),intent(inout) :: uinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)
  real(dp),intent(inout) :: vtinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine asrprs
end interface

interface
 subroutine  bec3(cg,cg1,dtefield,natom,d2lo,idirpert,ipert,mband,mkmem,&  
  &  mpw,mpw1,mpert,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: idirpert
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
  real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine bec3
end interface

interface
 subroutine cart29(blkflg,blkval,carflg,d2cart,&  
  &  gprimd,iblok,mpert,natom,nblok,ntypat,rprimd,typat,ucvol,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  integer,intent(out) :: carflg(3,mpert,3,mpert)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine cart29
end interface

interface
 subroutine cart39(flg1,flg2,gprimd,ipert,natom,rprimd,vec1,vec2)
  use defs_basis
  implicit none
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: flg1(3)
  integer,intent(out) :: flg2(3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vec1(3)
  real(dp),intent(out) :: vec2(3)
 end subroutine cart39
end interface

interface
 subroutine cgwf3(band,berryopt,cgq,cplex,cwavef,cwave0,cwaveprj,cwaveprj0,dcwavef,&  
  &  dimcprj,dimekb,dime1kb,dimffnlk,dimffnl1,dimphkxred,dkinpw,eig0nk,eig0_kq,eig1_k,&  
  &  ekb_typ,e1kbfr,e1kbsc,ffnlk,ffnlkq,ffnl1,filstat,gbound,ghc,gh1c_n,grad_berry,gsc,gscq,&  
  &  gs_hamkq,gvnlc,gvnl1,icgq,idir,indlmn_typ,&  
  &  ipert,igscq,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,&  
  &  kpt,lmnmax,matblk,mcgq,mgfft,mgscq,mpi_enreg,mpsang,mpssoang,mpw1,natom,nband,nbdbuf,&  
  &  nkpg,nkpg1,nline,npw,npw1,nspinor,ntypat,n4,n5,n6,opt_gvnl1,ortalg,paral_kgb,ph3d,phkxred,prtvol,&  
  &  quit,resid,sciss,sij_typ,tolwfr,usecprj,usedcwavef,usee1kb,vlocal,vlocal1,wfoptalg)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: band
  integer,intent(in) :: berryopt
  integer,intent(in) :: cplex
  integer,intent(in) :: dime1kb
  integer,intent(in) :: dimekb
  integer,intent(in) :: dimffnl1
  integer,intent(in) :: dimffnlk
  integer,intent(in) :: dimphkxred
  integer,intent(in) :: icgq
  integer,intent(in) :: idir
  integer,intent(in) :: igscq
  integer,intent(in) :: ipert
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgscq
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw1
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nbdbuf
  integer,intent(in) :: nkpg
  integer,intent(in) :: nkpg1
  integer,intent(in) :: nline
  integer,intent(in) :: npw
  integer,intent(in) :: npw1
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: opt_gvnl1
  integer,intent(in) :: ortalg
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: quit
  integer,intent(in) :: usecprj
  integer,intent(in) :: usedcwavef
  integer,intent(in) :: usee1kb
  integer,intent(in) :: wfoptalg
  real(dp),intent(in) :: eig0nk
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(in) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: resid
  real(dp),intent(in) :: sciss
  real(dp),intent(in) :: tolwfr
  real(dp),intent(in) :: cgq(2,mcgq)
  real(dp),intent(inout) :: cwave0(2,npw*nspinor)
  real(dp),intent(inout) :: cwavef(2,npw1*nspinor)
  type(cprj_type),intent(inout) :: cwaveprj(natom,nspinor)
  type(cprj_type),intent(inout) :: cwaveprj0(natom,nspinor*usecprj)
  real(dp),intent(inout) :: dcwavef(2,npw1*nspinor*((usedcwavef+1)/2))
  integer,intent(in) :: dimcprj(natom*gs_hamkq%usepaw)
  real(dp),intent(in) :: dkinpw(npw)
  real(dp),intent(in) :: e1kbfr(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
  real(dp),intent(in) :: e1kbsc(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
  real(dp),intent(in) :: eig0_kq(nband)
  real(dp),intent(out) :: eig1_k(2*nband**2)
  real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
  real(dp),intent(in) :: ffnl1(npw1,dimffnl1,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlk(npw,dimffnlk,lmnmax,1+gs_hamkq%usepaw*(ntypat-1))
  real(dp),intent(in) :: ffnlkq(npw1,dimffnl1,lmnmax,1)
  integer,intent(in) :: gbound(2*mgfft+8,2)
  real(dp),intent(out) :: gh1c_n(2,npw1*nspinor)
  real(dp),intent(out) :: ghc(2,npw1*nspinor)
  real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband)
  real(dp),intent(out) :: gsc(2,npw1*nspinor*gs_hamkq%usepaw)
  real(dp),intent(in) :: gscq(2,mgscq)
  real(dp),intent(out) :: gvnl1(2,npw1*nspinor)
  real(dp),intent(out) :: gvnlc(2,npw1*nspinor)
  integer,intent(in) :: indlmn_typ(6,lmnmax,1)
  integer,intent(in) :: kg1_k(3,npw1)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw1(npw1)
  real(dp),intent(in) :: kpg1_k(npw1,nkpg1)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(inout) :: ph3d(2,npw1,matblk)
  real(dp),intent(in) :: phkxred(2,dimphkxred)
  real(dp),intent(in) :: sij_typ(dimekb,gs_hamkq%usepaw)
  real(dp),intent(inout) :: vlocal(n4,n5,n6)
  real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6)
 end subroutine cgwf3
end interface

interface
 subroutine chkph3(carflg,idir,mpert,natom)
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: carflg(3,mpert,3,mpert)
 end subroutine chkph3
end interface

interface
 subroutine chneu9(chneut,d2cart,mpert,natom,ntypat,selectz,typat,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: chneut
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: selectz
  real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine chneu9
end interface

interface
 subroutine corrmetalwf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,edocc,eig1,fermie1,ghc,iband,&  
  &  ibgq,icgq,istwf_k,mcgq,mcprjq,mpi_enreg,natom,nband,npw1,nspinor,occ,rocceig,timcount,&  
  &  usepaw,wf_corrected)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: ibgq
  integer,intent(in) :: icgq
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcgq
  integer,intent(in) :: mcprjq
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: npw1
  integer,intent(in) :: nspinor
  integer,intent(in) :: timcount
  integer,intent(in) :: usepaw
  integer,intent(out) :: wf_corrected
  real(dp),intent(in) :: fermie1
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cgq(2,mcgq)
  type(cprj_type),intent(in) :: cprjq(natom,mcprjq)
  real(dp),intent(out) :: cwave1(2,npw1*nspinor)
  real(dp),intent(in) :: cwavef(2,npw1*nspinor)
  type(cprj_type),intent(in) :: cwaveprj(natom,nspinor*usepaw)
  type(cprj_type),intent(out) :: cwaveprj1(natom,nspinor*usepaw)
  real(dp),intent(out) :: edocc(nband)
  real(dp),intent(in) :: eig1(2*nband**2)
  real(dp),intent(in) :: ghc(2,npw1*nspinor)
  real(dp),intent(in) :: occ(nband)
  real(dp),intent(in) :: rocceig(nband,nband)
 end subroutine corrmetalwf1
end interface

interface
 subroutine d2kindstr2(cwavef,ecut,ecutsm,effmass,ekinout,gmet,gprimd,&  
  &  istwfk,kg_k,kpt,npw,nspinor)
  use defs_basis
  implicit none
  integer,intent(in) :: istwfk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  real(dp),intent(in) :: cwavef(2,npw*nspinor)
  real(dp),intent(out) :: ekinout(36)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpt(3)
 end subroutine d2kindstr2
end interface

interface
 subroutine d2sym3(blkflg,d2,indsym,mpert,natom,nsym,&  
  &  qpt,symq,symrec,symrel,timrev)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  integer,intent(inout) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(inout) :: d2(2,3,mpert,3,mpert)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: qpt(3)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine d2sym3
end interface

interface
 subroutine d3output(blkflg,d3,mband,mpert,nkpt,unddb)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: nkpt
  integer,intent(in) :: unddb
  integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert)
 end subroutine d3output
end interface

interface
 subroutine d3sym(blkflg,d3,indsym,mpert,natom,nsym,&  
  &  symrec,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert)
  real(dp),intent(inout) :: d3(2,3,mpert,3,mpert,3,mpert)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine d3sym
end interface

interface
 subroutine  die3(cg,cg1,dtefield,d2lo,idirpert,ipert,mband,mkmem,&  
  &  mpw,mpw1,mpert,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: idirpert
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
  real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine die3
end interface

interface
 subroutine dyfnl3(atindx1,becfrnl,cg,cprj,dimcprj,dyfrnl,dyfr_cplex,dyfr_nondiag,eigen,gsqcut,indsym,&  
  &  istwfk,kg,kptns,kptopt,mband,mgfft,mgfftf,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nband,&  
  &  nfftf,ngfft,ngfftf,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,nsym,ntypat,occ,&  
  &  paral_kgb,paw_ij,pawang,pawbec,pawprtvol,pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,ph1df,psps,&  
  &  qphon,rprimd,symafm,symrec,typat,unkg,unpaw,unylm,usecprj,usexcnhat,wfftgs,vtrial,&  
  &  vxc,wtk,xred,ylm)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: dyfr_nondiag
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgfftf
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: pawbec
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: unkg
  integer,intent(in) :: unpaw
  integer,intent(in) :: unylm
  integer,intent(in) :: usecprj
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wfftgs
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(out) :: becfrnl(3,natom,3*pawbec)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
  integer,intent(in) :: dimcprj(natom*psps%usepaw)
  real(dp),intent(out) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*natom)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(ntypat)
  real(dp),target,intent(in) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 end subroutine dyfnl3
end interface

interface
 subroutine dyfro3(atindx1,dyfrnl,dyfrlo,dyfrwf,dyfrxc,dyfr_cplex,dyfr_nondiag,&  
  &  gmet,gprimd,gsqcut,indsym,mgfft,mpi_enreg,mqgrid,natom,nattyp,&  
  &  nfft,ngfft,nspden,nsym,ntypat,n1xccc,n3xccc,paral_kgb,pawtab,ph1d,qgrid,&  
  &  qphon,rhog,rprimd,symq,symrec,typat,ucvol,usepaw,vlspl,vxc,&  
  &  xcccrc,xccc1d,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: dyfr_nondiag
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(out) :: dyfrlo(3,3,natom)
  real(dp),intent(inout) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(out) :: dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(out) :: dyfrxc(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine dyfro3
end interface

interface
 subroutine dyout3(becfrnl,berryopt,blkflg,carflg,ddboun,ddkfil,dyew,dyfrlo,dyfrnl,&  
  &  dyfrx1,dyfrx2,dyfr_cplex,dyfr_nondiag,d2cart,d2cart_bbb,&  
  &  d2eig0,d2k0,d2lo,d2loc0,d2matr,d2nl,d2nl0,d2nl1,d2ovl,d2vn,&  
  &  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&  
  &  iout,mband,mpert,natom,ntypat,&  
  &  outd2,pawbec,prtbbb,prtvol,qphon,qzero,typat,rfdir,rfpert,rfphon,rfstrs,usepaw,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: ddboun
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: dyfr_nondiag
  integer,intent(in) :: iout
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: outd2
  integer,intent(in) :: pawbec
  integer,intent(in) :: prtbbb
  integer,intent(in) :: prtvol
  integer,intent(in) :: qzero
  integer,intent(in) :: rfphon
  integer,intent(in) :: rfstrs
  integer,intent(in) :: usepaw
  integer,intent(in) :: ddkfil(3)
  integer,intent(in) :: rfdir(3)
  real(dp),intent(in) :: becfrnl(3,natom,3*pawbec)
  integer,intent(in) :: blkflg(3,mpert,3,mpert)
  integer,intent(in) :: carflg(3,mpert,3,mpert)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(in) :: d2eig0(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2k0(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2loc0(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2matr(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2nl0(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2nl1(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2ovl(2,3,mpert,3,mpert*usepaw)
  real(dp),intent(in) :: d2vn(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: dyfrlo(3,3,natom)
  real(dp),intent(in) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(in) :: dyfrx1(2,3,natom,3,natom)
  real(dp),intent(in) :: dyfrx2(3,3,natom)
  real(dp),intent(in) :: eltcore(6,6)
  real(dp),intent(in) :: elteew(6+3*natom,6)
  real(dp),intent(in) :: eltfrhar(6,6)
  real(dp),intent(in) :: eltfrkin(6,6)
  real(dp),intent(in) :: eltfrloc(6+3*natom,6)
  real(dp),intent(in) :: eltfrnl(6+3*natom,6)
  real(dp),intent(in) :: eltfrxc(6+3*natom,6)
  real(dp),intent(in) :: qphon(3)
  integer,intent(in) :: rfpert(mpert)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine dyout3
end interface

interface
 subroutine dyxc13(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,kxc,mgfft,mpert,mpi_enreg,mqgrid,&  
  &  natom,nfft,ngfft,nkxc,nspden,ntypat,n1xccc,paral_kgb,pawtab,&  
  &  ph1d,qgrid,qphon,rfdir,rfpert,rprimd,timrev,typat,ucvol,usepaw,xcccrc,xccc1d,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpert
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: timrev
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: rfdir(3)
  integer,intent(in) :: atindx(natom)
  integer,intent(out) :: blkflgfrx1(3,natom,3,natom)
  real(dp),intent(out) :: dyfrx1(2,3,natom,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: qphon(3)
  integer,intent(in) :: rfpert(mpert)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine dyxc13
end interface

interface
 subroutine  ebp3(cg,cg1,dtefield,eberry,mband,mkmem,&  
  &  mpw,mpw1,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(out) :: eberry
  real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
 end subroutine ebp3
end interface

interface
 subroutine  edie3(cg,cg1,dtefield,eberry,idir_efield,mband,mkmem,&  
  &  mpw,mpw1,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: idir_efield
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(out) :: eberry
  real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine edie3
end interface

interface
 subroutine eig1fixed(band,cwave0,dimekb,dimffnlk,dimffnl1,dkinpw,eig1_k,ekb_typ,&  
  &  ffnlk,ffnlkq,ffnl1,gs_hamkq,gvnl1,idir,indlmn_typ,&  
  &  ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lmnmax,matblk,mgfft,mpi_enreg,&  
  &  mpsang,mpssoang,natom,nband,nkpg,nkpg1,npw,npw1,nspinor,ntypat,ph3d,prtvol)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: band
  integer,intent(in) :: dimekb
  integer,intent(in) :: dimffnl1
  integer,intent(in) :: dimffnlk
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nkpg
  integer,intent(in) :: nkpg1
  integer,intent(in) :: npw
  integer,intent(in) :: npw1
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(gs_hamiltonian_type),intent(in) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cwave0(2,npw*nspinor)
  real(dp),intent(in) :: dkinpw(npw)
  real(dp),intent(out) :: eig1_k(2*nband**2)
  real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
  real(dp),intent(in) :: ffnl1(npw1,dimffnl1,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlk(npw,dimffnlk,lmnmax,1)
  real(dp),intent(in) :: ffnlkq(npw1,dimffnl1,lmnmax,1)
  real(dp),intent(out) :: gvnl1(2,npw1*nspinor)
  integer,intent(in) :: indlmn_typ(6,lmnmax,1)
  integer,intent(in) :: kg1_k(3,npw1)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw1(npw1)
  real(dp),intent(in) :: kpg1_k(npw1,nkpg1)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(inout) :: ph3d(2,npw1,matblk)
 end subroutine eig1fixed
end interface

interface
 subroutine eig2tot(bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0,eigenq,eigen1,eig2nkq,elph2_imagden,esmear,&  
  &  gh0c1_pert,gh1c_pert,ieig2rf,istwfk_pert,mband,mk1mem,mpert,npert,mpi_enreg,mpw1,&  
  &  nkpt_rbz,npwar1,nspinor,nsppol,smdelta,eigbrd)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdeigrf
  integer,intent(in) :: dim_eig2nkq
  integer,intent(in) :: dim_eig2rf
  integer,intent(in) :: ieig2rf
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: npert
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: smdelta
  real(dp),intent(in) :: elph2_imagden
  real(dp),intent(in) :: esmear
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
  integer,intent(in) :: clflg(3,mpert)
  real(dp),intent(out) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert*dim_eig2nkq)
  real(dp),intent(out),optional :: eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)
  real(dp),intent(in) :: eigen0(nkpt_rbz*mband*nsppol)
  real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
  real(dp),intent(in) :: eigenq(nkpt_rbz*mband*nsppol)
  real(dp),intent(in) :: gh0c1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
  real(dp),intent(in) :: gh1c_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
  integer,intent(in) :: istwfk_pert(nkpt_rbz,3,mpert)
  integer,intent(in) :: npwar1(nkpt_rbz,mpert)
 end subroutine eig2tot
end interface

interface
 subroutine eigen_meandege(eigen0,eigenresp,eigenresp_mean,mband,nband,nkpt,nsppol,option)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: option
  real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigenresp((3-option)*mband**(3-option)*nkpt*nsppol)
  real(dp),intent(out) :: eigenresp_mean(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine eigen_meandege
end interface

interface
 subroutine elph2_fanddw(dim_eig2nkq,displ,eig2nkq,eigen_corr,gprimd,mband,natom,nkpt,nsppol,option,phfrq)
  use defs_basis
  implicit none
  integer,intent(in) :: dim_eig2nkq
  integer,intent(in) :: mband
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: option
  real(dp),intent(in) :: displ(2*3*natom*3*natom)
  real(dp),intent(in) :: eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq)
  real(dp),intent(out) :: eigen_corr(mband*nkpt*nsppol)
  real(dp) :: gprimd(3,3)
  real(dp),intent(in) :: phfrq(3*natom)
 end subroutine elph2_fanddw
end interface

interface
 subroutine eltfrhar3(eltfrhar,rprimd,gsqcut,mpi_enreg,nfft,ngfft,rhog)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: eltfrhar(6,6)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine eltfrhar3
end interface

interface
 subroutine eltfrkin3(cg,eltfrkin,ecut,ecutsm,effmass,&  
  &  istwfk,kg,kptns,mband,mgfft,mkmem,mpi_enreg,&  
  &  mpw,nband,nkpt,ngfft,npwarr,nspinor,nsppol,occ,&  
  &  rprimd,unkg,wfftgs,wtk)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: unkg
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wfftgs
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(out) :: eltfrkin(6,6)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine eltfrkin3
end interface

interface
 subroutine eltfrloc3(atindx,eltfrloc,gmet,gprimd,gsqcut,mgfft,&  
  &  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,ph1d,qgrid,rhog,vlspl)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(out) :: eltfrloc(6+3*natom,6)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
 end subroutine eltfrloc3
end interface

interface
 subroutine eltfrnl3(atindx,atindx1,cg,eltfrnl,istwfk,&  
  &  kg,kptns,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nband,&  
  &  nkpt,ngfft,nloalg,npwarr,nspinor,nsppol,ntypat,occ,ph1d,&  
  &  psps,rprimd,unkg,wfftgs,unylm,useylm,wtk,xred,ylm,ylmgr)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: unkg
  integer,intent(in) :: unylm
  integer,intent(in) :: useylm
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wfftgs
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(out) :: eltfrnl(6+3*natom,6)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*useylm)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,9,mpsang*mpsang*useylm)
 end subroutine eltfrnl3
end interface

interface
 subroutine eltfrxc3(eltfrxc,enxc,kxc,mpi_enreg,natom,&  
  &  nfft,ngfft,nkxc,nspden,ntypat,n1xccc,n3xccc,paral_kgb,rhor,rprimd,&  
  &  typat,vxc,xcccrc,xccc1d,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: enxc
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: eltfrxc(6+3*natom,6)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine eltfrxc3
end interface

interface
 subroutine eltxccore(eltfrxc,is2_in,natom,nfft,ntypat,&  
  &  n1,n1xccc,n2,n3,rprimd,typat,ucvol,vxc_core,vxc10_core,vxc1is_core,&  
  &  xcccrc,xccc1d,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: is2_in
  integer,intent(in) :: n1
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(inout) :: eltfrxc(6+3*natom,6)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vxc10_core(nfft)
  real(dp),intent(in) :: vxc1is_core(nfft)
  real(dp),intent(in) :: vxc_core(nfft)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine eltxccore
end interface

interface
 subroutine etot3(berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,&  
  &  efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&  
  &  enl0,enl1,epaw1,etotal,evar,exc1,ipert,natom,optene)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: optene
  real(dp),intent(out) :: deltae
  real(dp),intent(in) :: eberry
  real(dp),intent(in) :: edocc
  real(dp),intent(in) :: eeig0
  real(dp),intent(in) :: eew
  real(dp),intent(in) :: efrhar
  real(dp),intent(in) :: efrkin
  real(dp),intent(in) :: efrloc
  real(dp),intent(in) :: efrnl
  real(dp),intent(in) :: efrx1
  real(dp),intent(in) :: efrx2
  real(dp),intent(in) :: ehart1
  real(dp),intent(in) :: eii
  real(dp),intent(in) :: ek0
  real(dp),intent(in) :: ek1
  real(dp),intent(inout) :: elast
  real(dp),intent(in) :: eloc0
  real(dp),intent(in) :: elpsp1
  real(dp),intent(in) :: enl0
  real(dp),intent(in) :: enl1
  real(dp),intent(in) :: epaw1
  real(dp),intent(out) :: etotal
  real(dp),intent(out) :: evar
  real(dp),intent(in) :: exc1
 end subroutine etot3
end interface

interface
 subroutine ewald3(dyew,gmet,natom,qphon,rmet,sumg0,typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: sumg0
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(*)
 end subroutine ewald3
end interface

interface
 subroutine ewald4(elteew,gmet,gprimd,natom,ntypat,rmet,rprimd,&  
  &  typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: elteew(6+3*natom,6)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine ewald4
end interface

interface
 subroutine gath3(becfrnl,berryopt,blkflg,carflg,dyew,dyfrwf,dyfrx1,&  
  &  dyfr_cplex,dyfr_nondiag,d2bbb,d2cart,d2cart_bbb,d2matr,d2nfr,&  
  &  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&  
  &  gprimd,mband,mpert,natom,ntypat,outd2,pawbec,prtbbb,&  
  &  rfasr,rfpert,rprimd,typat,ucvol,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: dyfr_nondiag
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: outd2
  integer,intent(in) :: pawbec
  integer,intent(in) :: prtbbb
  integer,intent(in) :: rfasr
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: becfrnl(3,natom,3*pawbec)
  integer,intent(inout) :: blkflg(3,mpert,3,mpert)
  integer,intent(out) :: carflg(3,mpert,3,mpert)
  real(dp),intent(in) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(out) :: d2matr(2,3,mpert,3,mpert)
  real(dp),intent(in) :: d2nfr(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
  real(dp),intent(in) :: dyfrx1(2,3,natom,3,natom)
  real(dp),intent(in) :: eltcore(6,6)
  real(dp),intent(in) :: elteew(6+3*natom,6)
  real(dp),intent(in) :: eltfrhar(6,6)
  real(dp),intent(in) :: eltfrkin(6,6)
  real(dp),intent(in) :: eltfrloc(6+3*natom,6)
  real(dp),intent(in) :: eltfrnl(6+3*natom,6)
  real(dp),intent(in) :: eltfrxc(6+3*natom,6)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: rfpert(mpert)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine gath3
end interface

interface
 subroutine gaugetransfo(cg_k,cwavef,cwavef_d,eig_k,eig1_k,iband,nband_k,&  
  &  mband,npw_k,npw1_k,nspinor,nsppol,occ_k)
  use defs_basis
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: mband
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw1_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: cg_k(2,npw_k*nspinor*nband_k)
  real(dp),intent(in) :: cwavef(2,npw1_k*nspinor)
  real(dp),intent(out) :: cwavef_d(2,npw1_k*nspinor)
  real(dp),intent(in) :: eig1_k(2*nsppol*mband**2)
  real(dp),intent(in) :: eig_k(mband*nsppol)
  real(dp),intent(in) :: occ_k(nband_k)
 end subroutine gaugetransfo
end interface

interface
 subroutine gbefd3(cg,cg1,dtefield,grad_berry,idir_efield,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt,&  
  &  npwarr,npwar1,nspinor,&  
  &  nsppol,qmat,pwindall,rprimd)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: idir_efield
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: grad_berry(2,mpw1,dtefield%nband_occ)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine gbefd3
end interface

interface
 subroutine getgh1c(berryopt,cplex,cwave,cwaveprj,dimcprj,dimekb,dime1kb,dimffnlk,dimffnlkq,dimffnl1,dimphkxred,dkinpw,&  
  &  ekb_typ,e1kbfr,e1kbsc,ffnlk,ffnlkq,ffnl1,filstat,gbound,gh1c,grad_berry,gs1c,gs_hamkq,gvnl1,idir,&  
  &  indlmn_typ,ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&  
  &  mpssoang,natom,nkpg,nkpg1,npw,npw1,nspinor,ntypat,n4,n5,n6,optlocal,optnl,opt_gvnl1,paral_kgb,ph3d,&  
  &  phkxred,prtvol,sij_opt,sij_typ,tim_getgh1c,usecprj,usee1kb,usevnl,vlocal1,wfraug)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: cplex
  integer,intent(in) :: dime1kb
  integer,intent(in) :: dimekb
  integer,intent(in) :: dimffnl1
  integer,intent(in) :: dimffnlk
  integer,intent(in) :: dimffnlkq
  integer,intent(in) :: dimphkxred
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nkpg
  integer,intent(in) :: nkpg1
  integer,intent(in) :: npw
  integer,intent(in) :: npw1
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: opt_gvnl1
  integer,intent(in) :: optlocal
  integer,intent(in) :: optnl
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: sij_opt
  integer,intent(in) :: tim_getgh1c
  integer,intent(in) :: usecprj
  integer,intent(in) :: usee1kb
  integer,intent(in) :: usevnl
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(in) :: gs_hamkq
  real(dp),intent(in) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cwave(2,npw*nspinor)
  type(cprj_type),intent(inout),target :: cwaveprj(natom,nspinor*usecprj)
  integer,intent(in) :: dimcprj(natom*gs_hamkq%usepaw)
  real(dp),intent(in) :: dkinpw(npw)
  real(dp),intent(in) :: e1kbfr(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
  real(dp),intent(in) :: e1kbsc(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
  real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
  real(dp),intent(in) :: ffnl1(npw1,dimffnl1,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlk(npw,dimffnlk,lmnmax,1+gs_hamkq%usepaw*(ntypat-1))
  real(dp),intent(in) :: ffnlkq(npw1,dimffnlkq,lmnmax,1)
  integer,intent(in) :: gbound(2*mgfft+8,2)
  real(dp),intent(out) :: gh1c(2,npw1*nspinor)
  real(dp),intent(in) :: grad_berry(2,npw1*nspinor*(berryopt/4))
  real(dp),intent(out) :: gs1c(2,npw1*nspinor*((sij_opt+1)/2))
  real(dp),intent(inout),target :: gvnl1(2,npw1*nspinor*usevnl)
  integer,intent(in) :: indlmn_typ(6,lmnmax,1)
  integer,intent(in) :: kg1_k(3,npw1)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw1(npw1)
  real(dp),intent(in) :: kpg1_k(npw1,nkpg1)
  real(dp),intent(in) :: kpg_k(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(inout) :: ph3d(2,npw1,matblk)
  real(dp),intent(in) :: phkxred(2,dimphkxred)
  real(dp),intent(in) :: sij_typ(dimekb,gs_hamkq%usepaw)
  real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6*optlocal)
  real(dp),intent(inout) :: wfraug(2,n4,n5,n6*optlocal)
 end subroutine getgh1c
end interface

interface
 subroutine getshell(gmet,kneigh,kg_neigh,kptindex,kptopt,kptrlatt,kpt2,&  
  &  kpt3,mkmem,mkmem_max,mpi_enreg,mvwtk,&  
  &  nkpt2,nkpt3,nneigh,nshiftk,rmet,rprimd,shiftk,wtk2)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mkmem
  integer,intent(out) :: mkmem_max
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nkpt3
  integer,intent(out) :: nneigh
  integer,intent(inout) :: nshiftk
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(inout) :: kptrlatt(3,3)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(out) :: kg_neigh(30,nkpt2,3)
  integer,intent(out) :: kneigh(30,nkpt2)
  real(dp),intent(in) :: kpt2(3,nkpt2)
  real(dp),intent(out) :: kpt3(3,nkpt3)
  integer,intent(out) :: kptindex(2,nkpt3)
  real(dp),intent(out) :: mvwtk(30,nkpt2)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: shiftk(3,nshiftk)
  real(dp),intent(in) :: wtk2(nkpt2)
 end subroutine getshell
end interface

interface
 subroutine gradberry3(cg,cg1,dtefield,grad_berry,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt,&  
  &  npwarr,npwar1,nspinor,nsppol,qmat,pwindall)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: grad_berry(2,mpw1,dtefield%nband_occ)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
 end subroutine gradberry3
end interface

interface
 subroutine hartrestr(gmet,gprimd,gsqcut,idir,ipert,mpi_enreg,natom,nfft,ngfft,&  
  &  paral_kgb,rhog,vhartr1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: vhartr1(nfft)
 end subroutine hartrestr
end interface

interface
 subroutine  initberry3(dtefield,dtfil,dtset,gmet,kg,kg1,mband,mkmem,mpi_enreg,&  
  &  mpw,mpw1,nkpt,npwarr,npwar1,nsppol,occ,pwindall,rprimd)
  use defs_basis
  use m_efield
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(efield_type),intent(out) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mkmem)
  integer,intent(in) :: npwar1(nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  integer,intent(out) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine initberry3
end interface

interface
 subroutine inpgkk(bantot1,eigen1,filegkk,hdr1,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bantot1
  character(len=fnlen),intent(in) :: filegkk
  type(hdr_type), intent(out) :: hdr1
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(out) :: eigen1(bantot1)
 end subroutine inpgkk
end interface

interface
 subroutine inprep8 (dimekb,filnam,lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,&  
  &  ntypat,unddb,usepaw,vrsddb)
  use defs_basis
  implicit none
  integer,intent(out) :: dimekb
  integer,intent(out) :: lmnmax
  integer,intent(out) :: mband
  integer,intent(out) :: mblktyp
  integer, intent(out) :: msym
  integer,intent(out) :: natom
  integer,intent(out) :: nblok
  integer,intent(out) :: nkpt
  integer,intent(out) :: ntypat
  integer,intent(in) :: unddb
  integer,intent(out) :: usepaw
  integer,intent(in) :: vrsddb
  character(len=fnlen),intent(in) :: filnam
 end subroutine inprep8
end interface

interface
 subroutine insy3(gprimd,idir,indsym,ipert,natom,nsym,nsym1,&  
  &  rfmeth,symafm,symaf1,symq,symrec,symrel,symrl1,syuse,tnons,tnons1)
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(out) :: nsym1
  integer,intent(in) :: rfmeth
  integer,intent(in) :: syuse
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: symaf1(nsym)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  integer,intent(out) :: symrl1(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(out) :: tnons1(3,nsym)
 end subroutine insy3
end interface

interface
 subroutine ioddb8_in(filnam,matom,mband,&  
  &  mkpt,msym,mtypat,unddb,vrsddb,&  
  &  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&  
  &  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&  
  &  pawecutdg,rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&  
  &  typat,usepaw,wtk,xred,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(out) :: intxc
  integer,intent(out) :: iscf
  integer,intent(out) :: ixc
  integer,intent(in) :: matom
  integer,intent(in) :: mband
  integer,intent(in) :: mkpt
  integer,intent(in) :: msym
  integer,intent(in) :: mtypat
  integer,intent(out) :: natom
  integer,intent(out) :: nkpt
  integer,intent(out) :: nspden
  integer,intent(out) :: nspinor
  integer,intent(out) :: nsppol
  integer,intent(out) :: nsym
  integer,intent(out) :: ntypat
  integer,intent(out) :: occopt
  integer,intent(in) :: unddb
  integer,intent(in) :: usepaw
  integer,intent(in) :: vrsddb
  real(dp),intent(out) :: dilatmx
  real(dp),intent(out) :: ecut
  real(dp),intent(out) :: ecutsm
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(out) :: kptnrm
  real(dp),intent(out) :: pawecutdg
  real(dp),intent(out) :: sciss
  real(dp),intent(out) :: tolwfr
  real(dp),intent(out) :: tphysel
  real(dp),intent(out) :: tsmear
  integer,intent(out) :: ngfft(18)
  real(dp),intent(out) :: acell(3)
  real(dp),intent(out) :: amu(mtypat)
  real(dp),intent(out) :: kpt(3,mkpt)
  integer,intent(out) :: nband(mkpt)
  real(dp),intent(out) :: occ(mband*mkpt)
  real(dp),intent(out) :: rprim(3,3)
  real(dp),intent(out) :: spinat(3,matom)
  integer,intent(out) :: symafm(msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
  integer,intent(out) :: typat(matom)
  real(dp),intent(out) :: wtk(mkpt)
  real(dp),intent(out) :: xred(3,matom)
  real(dp),intent(out) :: zion(mtypat)
  real(dp),intent(out) :: znucl(mtypat)
 end subroutine ioddb8_in
end interface

interface
 subroutine ioddb8_out (dscrpt,filnam,matom,mband,&  
  &  mkpt,msym,mtypat,unddb,vrsddb,&  
  &  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&  
  &  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&  
  &  pawecutdg,rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&  
  &  typat,usepaw,wtk,xred,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: intxc
  integer,intent(in) :: iscf
  integer,intent(in) :: ixc
  integer,intent(in) :: matom
  integer,intent(in) :: mband
  integer,intent(in) :: mkpt
  integer,intent(in) :: msym
  integer,intent(in) :: mtypat
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: unddb
  integer,intent(in) :: usepaw
  integer,intent(in) :: vrsddb
  real(dp),intent(in) :: dilatmx
  character(len=fnlen),intent(in) :: dscrpt
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(in) :: kptnrm
  real(dp),intent(in) :: pawecutdg
  real(dp),intent(in) :: sciss
  real(dp),intent(in) :: tolwfr
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(mtypat)
  real(dp),intent(in) :: kpt(3,mkpt)
  integer,intent(in) :: nband(mkpt)
  real(dp),intent(in) :: occ(mband*mkpt)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: spinat(3,matom)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(in) :: tnons(3,msym)
  integer,intent(in) :: typat(matom)
  real(dp),intent(in) :: wtk(mkpt)
  real(dp),intent(in) :: xred(3,matom)
  real(dp),intent(in) :: zion(mtypat)
  real(dp),intent(in) :: znucl(mtypat)
 end subroutine ioddb8_out
end interface

interface
 subroutine irred_perts(Dtset,out_unt)
  use defs_abitypes
  implicit none
  integer,intent(in) :: out_unt
  type(dataset_type),intent(in) :: Dtset
 end subroutine irred_perts
end interface

interface
 subroutine kpg3(dkinpw,ecut,ecutsm,effmass,gmet,idir,kg,kpt,npw)
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: npw
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  real(dp),intent(out) :: dkinpw(npw)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpt(3)
 end subroutine kpg3
end interface

interface
 subroutine kpgstr(dkinpw,ecut,ecutsm,effmass,gmet,gprimd,istr,kg,kpt,npw)
  use defs_basis
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: npw
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  real(dp),intent(out) :: dkinpw(npw)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpt(3)
 end subroutine kpgstr
end interface

interface
 subroutine mkcor3(cplex,idir,ipert,natom,ntypat,n1,n1xccc,&  
  &  n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xccc3d1,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: n1
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(out) :: xccc3d1(cplex*n1*n2*n3)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mkcor3
end interface

interface
 subroutine mkrho3(cg,cg1,cplex,gprimd,irrzon,istwfk_rbz,&  
  &  kg,kg1,mband,mgfft,mkmem,mk1mem,mpi_enreg,mpw,mpw1,nband_rbz,&  
  &  nfft,ngfft,nkpt_rbz,npwarr,npwar1,nspden,nspinor,nsppol,nsym,&  
  &  occ_rbz,paral_kgb,phnons,rhog1,rhor1,rprimd,symafm,symrel,&  
  &  ucvol,unkg,unkg1,wffnow,wfftgs,wtk_rbz)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: unkg
  integer,intent(in) :: unkg1
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: rhog1(2,nfft)
  real(dp),intent(out) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
 end subroutine mkrho3
end interface

interface
 subroutine mkvxcstr3(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,&  
  &  nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,vxc1,xccc3d1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 end subroutine mkvxcstr3
end interface

interface
 subroutine mkvxcstrgga3(cplex,dgprimdds,gprimd,kxc,mpi_enreg,nfft,ngfft,&  
  &  nkxc,nspden,paral_kgb,qphon,rhor1tmp,vxc1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dgprimdds(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhor1tmp(cplex*nfft,2)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
 end subroutine mkvxcstrgga3
end interface

interface
 subroutine newfermie1(cplex,fermie1,fe1fixed,istep,&  
  &  mpi_enreg,nfft,nfftot,nspden,occopt,&  
  &  prtvol,rhorfermi,ucvol,vtrial1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: istep
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: fe1fixed
  real(dp),intent(inout) :: fermie1
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: rhorfermi(nfft,nspden)
  real(dp),intent(in) :: vtrial1(cplex*nfft,nspden)
 end subroutine newfermie1
end interface

interface
 subroutine newvtr3(cplex,dbl_nnsclo,dielar,dtset,etotal,ffttomix,&  
  &  initialized,iscf,ispmix,istep,mix,mixtofft,&  
  &  mpi_enreg,nfft,nfftmix,ngfft,ngfftmix,npawmix,pawrhoij,&  
  &  qphon,rhor,rprimd,usepaw,vresid,vtrial)
  use defs_basis
  use m_ab6_mixing
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(out) :: dbl_nnsclo
  integer,intent(in) :: initialized
  integer,intent(in) :: iscf
  integer,intent(in) :: ispmix
  integer,intent(in) :: istep
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftmix
  integer,intent(in) :: npawmix
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  type(ab6_mixing_object), intent(inout) :: mix
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftmix(18)
  real(dp),intent(in) :: dielar(7)
  integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
  integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft))
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*usepaw)
  real(dp),intent(in) :: qphon(3)
  real(dp), intent(in), target :: rhor(cplex*nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vresid(cplex*nfft,dtset%nspden)
  real(dp),intent(inout) :: vtrial(cplex*nfft,dtset%nspden)
 end subroutine newvtr3
end interface

interface
 subroutine nlopt(blkflg,carflg,d3,d3cart,gprimd,mpert,natom,rprimd,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
  integer,intent(out) :: carflg(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert)
  real(dp),intent(out) :: d3cart(2,3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine nlopt
end interface

interface
 subroutine nselt3(atindx,atindx1,blkflg,cg,cg1,cplex,&  
  &  d2bbb,d2lo,d2nl,ecut,ecutsm,effmass,&  
  &  gmet,gprimd,gsqcut,idir,&  
  &  ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mband,mgfft,&  
  &  mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&  
  &  natom,nattyp,nband_rbz,nfft,ngfft,&  
  &  nkpt_rbz,nkxc,nloalg,npwarr,npwar1,nspden,nspinor,nsppol,&  
  &  nsym1,ntypat,occ_rbz,&  
  &  paral_kgb, ph1d,prtbbb,psps,qphon,rhog,&  
  &  rhor,rhor1,rmet,rprimd,symrc1,type,ucvol,&  
  &  unkg,unkg1,&  
  &  wffnow,wfftgs,unylm,unylm1,wtk_rbz,&  
  &  xred,ylm,ylm1,ylmgr,ylmgr1)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtbbb
  integer,intent(in) :: unkg
  integer,intent(in) :: unkg1
  integer,intent(in) :: unylm
  integer,intent(in) :: unylm1
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  integer,intent(out) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: type(natom)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang)
 end subroutine nselt3
end interface

interface
 subroutine nstdy3(atindx,blkflg,cg,cg1,cplex,dtfil,dtset,d2bbb,d2lo,d2nl,eigen0,eigen1,&  
  &  gmet,gsqcut,idir,indkpt1,indsy1,ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,&  
  &  mpert,mpi_enreg,mpw,mpw1,nattyp,nband_rbz,nfft,ngfft,nkpt,nkpt_rbz,nkxc,&  
  &  npwarr,npwar1,nspden,nsppol,nsym1,occ_rbz,ph1d,psps,rhor1,rmet,rprimd,&  
  &  symrc1,ucvol,wffnow,wfftgs,wtk_rbz,xred,ylm,ylm1)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(out) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*dtset%mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*dtset%mk1mem*nsppol)
  real(dp),intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*dtset%prtbbb)
  real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: indkpt1(nkpt_rbz)
  integer,intent(in) :: indsy1(4,nsym1,dtset%natom)
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*dtset%mkmem)
  integer,intent(in) :: kg1(3,mpw1*dtset%mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(dtset%ntypat)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrc1(3,3,nsym1)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*dtset%mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine nstdy3
end interface

interface
 subroutine nstpaw3(blkflg,cg,cgq,cg1,cplex,cprj,cprjq,dimcprj,docckqde,doccde_rbz,dtfil,dtset,d2lo,d2nl,d2ovl,&  
  &  eigenq,eigen0,eigen1,eovl1,gmet,gprimd,gsqcut,idir,indkpt1,indsy1,ipert,irrzon1,istwfk_rbz,&  
  &  kg,kg1,kpt_rbz,kxc,mgfftf,mpert,mpi_enreg,mpw,mpw1,nband_rbz,ncpgr,nfftf,ngfftf,nhat1,&  
  &  nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,nspinor,nsppol,nsym1,n3xccc,occkq,occ_rbz,&  
  &  paw_an,paw_an1,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,&  
  &  pawrhoij1,pawtab,phnons1,ph1d,ph1df,psps,rhor1,rmet,rprimd,symaf1,symrc1,symrl1,&  
  &  ucvol,usecprj,usexcnhat,useylmgr1,vhartr1,vpsp1,vtrial,vtrial1,vxc,&  
  &  wffnow,wfftgs,wfftkq,wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr1)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mgfftf
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: n3xccc
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  integer,intent(in) :: usecprj
  integer,intent(in) :: usexcnhat
  integer,intent(in) :: useylmgr1
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eovl1
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in) :: pawang1
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  type(wffile_type),intent(inout) :: wfftkq
  integer,intent(in) :: ngfftf(18)
  integer,intent(out) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*nspinor*dtset%mband*dtset%mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*dtset%mband*dtset%mk1mem*nsppol)
  real(dp),intent(in) :: cgq(2,mpw1*nspinor*dtset%mband*dtset%mkqmem*nsppol)
  type(cprj_type),intent(in) :: cprj(dtset%natom,nspinor*dtset%mband*dtset%mkmem*nsppol*usecprj)
  type(cprj_type),intent(in) :: cprjq(dtset%natom,nspinor*dtset%mband*dtset%mkqmem*nsppol*usecprj)
  real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2ovl(2,3,mpert,3,mpert)
  integer,intent(in) :: dimcprj(dtset%natom*psps%usepaw)
  real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indkpt1(nkpt_rbz)
  integer,intent(in) :: indsy1(4,nsym1,dtset%natom)
  integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*dtset%mkmem)
  integer,intent(in) :: kg1(3,mpw1*dtset%mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  real(dp),intent(in) :: kxc(nfftf,nkxc)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  real(dp),intent(in) :: nhat1(cplex*nfftf,nspden*psps%usepaw)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
  type(paw_an_type),intent(in) :: paw_an(dtset%natom*psps%usepaw)
  type(paw_an_type),intent(inout) :: paw_an1(dtset%natom*psps%usepaw)
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij1(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij1(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
  real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: rhor1(cplex*nfftf,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: symrl1(3,3,nsym1)
  real(dp),intent(in) :: vhartr1(cplex*nfftf)
  real(dp),target,intent(in) :: vpsp1(cplex*nfftf)
  real(dp),target,intent(in) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vtrial1(cplex*nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),target,intent(in) :: xccc3d1(cplex*n3xccc)
  real(dp),intent(in) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*dtset%mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr1(mpw1*dtset%mk1mem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 end subroutine nstpaw3
end interface

interface
 subroutine nstwf3(cg,cg1,ddkfil,dtfil,dtset,d2bbb_k,d2nl_k,eig_k,eig1_k,gs_hamkq,&  
  &  icg,icg1,idir,ikpt,ipert,isppol,istwf_k,kg_k,kg1_k,kpt,mpert,&  
  &  mpi_enreg,mpsang,mpw,mpw1,nband_k,npw_k,npw1_k,nsppol,&  
  &  occ_k,psps,rmet,wffddk,wffnow,wfftgs,wtk_k,ylm,ylm1)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: idir
  integer,intent(in) :: ikpt
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mpert
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(inout) :: nband_k
  integer,intent(inout) :: npw1_k
  integer,intent(inout) :: npw_k
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(in) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  real(dp),intent(in) :: wtk_k
  integer,intent(in) :: ddkfil(3)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*dtset%mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*dtset%mk1mem*nsppol)
  real(dp),intent(out) :: d2bbb_k(2,3,dtset%mband,dtset%mband*dtset%prtbbb)
  real(dp),intent(out) :: d2nl_k(2,3,mpert)
  real(dp),intent(inout) :: eig1_k(2*nsppol*dtset%mband**2)
  real(dp),intent(in) :: eig_k(dtset%mband*nsppol)
  integer,intent(in) :: kg1_k(3,npw1_k)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(in) :: rmet(3,3)
  type(wffile_type),intent(inout) :: wffddk(3)
  real(dp),intent(in) :: ylm(npw_k,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(npw1_k,mpsang*mpsang*psps%useylm)
 end subroutine nstwf3
end interface

interface
 subroutine nstwf4(atindx,atindx1,cg,cg1,d2nl_k,ecut,ecutsm,effmass,&  
  &  gmet,gprimd,icg,icg1,ikpt,isppol,istwf_k,kg_k,kg1_k,kpt,mband,&  
  &  mgfft,mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&  
  &  natom,nattyp,nband_k,ngfft,nloalg,npw_k,npw1_k,nspinor,nsppol,&  
  &  ntypat,occ_k,ph1d,psps,rmet,ucvol,wffnow,wfftgs,&  
  &  wtk_k,xred,ylm,ylmgr)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(inout) :: nband_k
  integer,intent(inout) :: npw1_k
  integer,intent(inout) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  real(dp),intent(in) :: wtk_k
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: d2nl_k(2,3,mpert)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg1_k(3,npw1_k)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kpt(3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(npw_k,mpsang*mpsang)
  real(dp),intent(in) :: ylmgr(npw_k,3,mpsang*mpsang)
 end subroutine nstwf4
end interface

interface
 subroutine outbsd(bdeigrf,dtset,eig2nkq,mpert,nkpt_rbz,unitout)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdeigrf
  integer,intent(in) :: mpert
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: unitout
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: eig2nkq(2,dtset%mband*dtset%nsppol,nkpt_rbz,3,mpert,3,mpert)
 end subroutine outbsd
end interface

interface
 subroutine outgkk(bantot0,bantot1,outfile,eigen0,eigen1,hdr0,hdr1,mpi_enreg,phasecg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bantot0
  integer,intent(in) :: bantot1
  type(hdr_type),intent(inout) :: hdr0
  type(hdr_type),intent(inout) :: hdr1
  type(mpi_type),intent(inout) :: mpi_enreg
  character(len=fnlen),intent(in) :: outfile
  real(dp),intent(in) :: eigen0(bantot0)
  real(dp),intent(in) :: eigen1(2*bantot1)
  real(dp),intent(in) :: phasecg(2,bantot1)
 end subroutine outgkk
end interface

interface
 subroutine phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&  
  &  mpert,msym,natom,nsym,ntypat,phfrq,qphnrm,qphon,rprimd,&  
  &  symdynmat,symrel,typat,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: symdynmat
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(out) :: displ(2*3*natom*3*natom)
  real(dp),intent(out) :: eigval(3*natom)
  real(dp),intent(out) :: eigvec(2*3*natom*3*natom)
  integer,intent(in) :: indsym(4,msym*natom)
  real(dp),intent(out) :: phfrq(3*natom)
  real(dp),intent(inout) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine phfrq3
end interface

interface
 subroutine prtene3(berryopt,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&  
  &  ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,eovl1,epaw1,exc1,iout,ipert,natom,usepaw)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: iout
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: eberry
  real(dp),intent(in) :: edocc
  real(dp),intent(in) :: eeig0
  real(dp),intent(in) :: eew
  real(dp),intent(in) :: efrhar
  real(dp),intent(in) :: efrkin
  real(dp),intent(in) :: efrloc
  real(dp),intent(in) :: efrnl
  real(dp),intent(in) :: efrx1
  real(dp),intent(in) :: efrx2
  real(dp),intent(in) :: ehart01
  real(dp),intent(in) :: ehart1
  real(dp),intent(in) :: eii
  real(dp),intent(in) :: ek0
  real(dp),intent(in) :: ek1
  real(dp),intent(in) :: eloc0
  real(dp),intent(in) :: elpsp1
  real(dp),intent(in) :: enl0
  real(dp),intent(in) :: enl1
  real(dp),intent(in) :: eovl1
  real(dp),intent(in) :: epaw1
  real(dp),intent(in) :: exc1
 end subroutine prtene3
end interface

interface
 subroutine prtph3(displ,eivec,enunit,iodyn,iout,natom,phfrq,qphnrm,qphon)
  use defs_basis
  implicit none
  integer,intent(in) :: eivec
  integer,intent(in) :: enunit
  integer,intent(in) :: iodyn
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: qphon(3)
 end subroutine prtph3
end interface

interface
 subroutine psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&  
  &  nblok,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: dimekb
  integer,intent(inout) :: fullinit
  integer,intent(in) :: lmnmax
  integer,intent(inout) :: nblok
  integer,intent(in) :: ntypat
  integer,intent(in) :: nunit
  integer,intent(in) :: usepaw
  integer,intent(in) :: useylm
  integer,intent(in) :: vrsddb
  real(dp),intent(inout) :: ekb(dimekb,ntypat)
  integer,intent(inout) :: indlmn(6,lmnmax,ntypat)
  type(pawtab_type),intent(inout) :: pawtab(ntypat*usepaw)
  integer,intent(in) :: pspso(ntypat)
 end subroutine psddb8
end interface

interface
 subroutine q0dy3_apply(natom,dyewq0,dyew)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(inout) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: dyewq0(3,3,natom)
 end subroutine q0dy3_apply
end interface

interface
 subroutine q0dy3_calc(natom,dyewq0,dyew,option)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: option
  real(dp),intent(in) :: dyew(2,3,natom,3,natom)
  real(dp),intent(out) :: dyewq0(3,3,natom)
 end subroutine q0dy3_calc
end interface

interface
 subroutine  qmatrix(cg,dtefield,qmat,mpw,mpw1,mkmem,mband,npwarr,nkpt,nspinor,nsppol,pwindall)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(efield_type),intent(in) :: dtefield
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(out) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3)
 end subroutine qmatrix
end interface

interface
 subroutine read_blok8(ddb_blk,iblok,&  
  &  mband,mpert,msize,nkpt,nunit,&  
  &  blkval2,kpt) !optional
  use defs_basis
  use m_ddb_blk
  implicit none
  integer, intent(in) :: iblok
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: nkpt
  integer,intent(in) :: nunit
  type(ddb_blk_type), pointer :: ddb_blk
  real(dp),intent(out),optional :: blkval2(2,msize,mband,nkpt)
  real(dp),intent(out),optional :: kpt(3,nkpt)
 end subroutine read_blok8
end interface

interface
 subroutine redgr (frin,frredgr,mpi_enreg,nfft,ngfft,paral_kgb)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  type(mpi_type) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: frin(nfft)
  real(dp),intent(out) :: frredgr(nfft,3)
 end subroutine redgr
end interface

interface
 subroutine resp3dte(atindx,atindx1,cg,cg1,cg3,cplex,dtfil,dtset,d3lo,&  
  &  gmet,gprimd,i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,&  
  &  kg,mband,mgfft,mkmem,mk1mem,&  
  &  mpert,mpi_enreg,mpsang,mpw,natom,nfft,nkpt,nspden,nspinor,nsppol,&  
  &  npwarr,occ,ph1d,psps,rmet,ucvol,vtrial1,xred,ylm)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: i1dir
  integer,intent(in) :: i1pert
  integer,intent(in) :: i2dir
  integer,intent(in) :: i2pert
  integer,intent(in) :: i3dir
  integer,intent(in) :: i3pert
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: d3lo(2,3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 end subroutine resp3dte
end interface

interface
 subroutine rhofermi3(atindx,atindx1,cg,cgq,cplex,&  
  &  doccde_rbz,docckqde,dtfil,dtset,&  
  &  edocc,eeig0,eigenq,eigen0,eigen1,&  
  &  fe1fixed,gmet,gprimd,idir,&  
  &  ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,mband,mgfft,&  
  &  mkmem,mkqmem,mk1mem,mpi_enreg,mpsang,mpw,mpw1,&  
  &  natom,nattyp,nband_rbz,nfft,nkpt_rbz,npwarr,npwar1,nspden,&  
  &  nsppol,nsym1,ntypat,occkq,occ_rbz,phnons1,&  
  &  ph1d,prtvol,psps,rhorfermi,rmet,rprimd,symaf1,symrl1,ucvol,&  
  &  wfftgs,wfftkq,wtk_rbz,xred,ylm,ylm1,ylmgr1)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkqmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: edocc
  real(dp),intent(inout) :: eeig0
  real(dp),intent(out) :: fe1fixed
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wfftgs
  type(wffile_type),intent(inout) :: wfftkq
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*mband*mkqmem*nsppol)
  real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
  real(dp),intent(out) :: eigen1(2*mband*mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  integer,intent(in) :: npwar1(nkpt_rbz,2)
  integer,intent(in) :: npwarr(nkpt_rbz,2)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: occkq(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: rhorfermi(cplex*nfft,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrl1(3,3,nsym1)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm)
 end subroutine rhofermi3
end interface

interface
 subroutine rhotov3(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&  
  &  kxc,mpi_enreg,natom,nfft,ngfft,nhat1,nkxc,nspden,n3xccc,&  
  &  optene,optres,paral_kgb,qphon,rhog,rhog1,rhor,rhor1,&  
  &  rprimd,ucvol,usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,xccc3d1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(out) :: ehart01
  real(dp),intent(out) :: ehart1
  real(dp),intent(out) :: elpsp1
  real(dp),intent(out) :: exc1
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vres2
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat1(cplex*nfft,nspden*usepaw)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhog1(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vhartr1(cplex*nfft)
  real(dp),intent(in) :: vpsp1(cplex*nfft)
  real(dp),intent(out) :: vresid1(cplex*nfft,nspden)
  real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 end subroutine rhotov3
end interface

interface
 subroutine scfcv3(atindx,atindx1,blkflg,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cpus,&  
  &  dimcprj,dim_eig2rf,doccde_rbz,docckqde,dtfil,dtset,&  
  &  d2bbb,d2lo,d2nl,d2ovl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&  
  &  ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek0,ek1,eloc0,elpsp1,&  
  &  enl0,enl1,eovl1,epaw1,etotal,exc1,fermie,gh0c1_set,gh1c_set,hdr,idir,indkpt1,&  
  &  indsy1,initialized,ipert,irrzon1,istwfk_rbz,&  
  &  kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&  
  &  mpert,mpi_enreg,mpsang,mpw,mpw1,nattyp,nband_rbz,ncpgr,&  
  &  nfftf,ngfftf,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&  
  &  nsym1,n3xccc,occkq,occ_rbz,&  
  &  paw_an,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,pawrhoij1,pawtab,&  
  &  pertcase,phnons1,ph1d,ph1df,&  
  &  prtbbb,psps,qphon,resid,residm,rhog,rhog1,&  
  &  rhor,rhor1,rprimd,symaf1,symrc1,symrl1,&  
  &  usecprj,useylmgr,useylmgr1,wffddk,wffnew,wffnow,wfftgs,wfftkq,vpsp1,vtrial,vxc,&  
  &  wtk_rbz,wvl,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: dim_eig2rf
  integer,intent(in) :: idir
  integer,intent(inout) :: initialized
  integer,intent(in) :: ipert
  integer,intent(in) :: mgfftf
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkqmem
  integer,intent(in) :: mpert
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: n3xccc
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym1
  integer,intent(in) :: pertcase
  integer,intent(in) :: prtbbb
  integer,intent(in) :: usecprj
  integer,intent(in) :: useylmgr
  integer,intent(in) :: useylmgr1
  real(dp),intent(in) :: cpus
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eberry
  real(dp),intent(out) :: edocc
  real(dp),intent(out) :: eeig0
  real(dp),intent(in) :: eew
  real(dp),intent(in) :: efrhar
  real(dp),intent(in) :: efrkin
  real(dp),intent(in) :: efrloc
  real(dp),intent(in) :: efrnl
  real(dp),intent(in) :: efrx1
  real(dp),intent(in) :: efrx2
  real(dp),intent(out) :: ehart01
  real(dp),intent(out) :: ehart1
  real(dp),intent(in) :: eii
  real(dp),intent(out) :: ek0
  real(dp),intent(out) :: ek1
  real(dp),intent(out) :: eloc0
  real(dp),intent(out) :: elpsp1
  real(dp),intent(out) :: enl0
  real(dp),intent(out) :: enl1
  real(dp),intent(out) :: eovl1
  real(dp),intent(out) :: epaw1
  real(dp),intent(out) :: etotal
  real(dp),intent(out) :: exc1
  real(dp),intent(inout) :: fermie
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in) :: pawang1
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(out) :: residm
  type(wffile_type),intent(inout) :: wffddk
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  type(wffile_type),intent(inout) :: wfftkq
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  integer,intent(out) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*dtset%nsppol)
  real(dp),intent(inout) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol)
  real(dp),intent(out) :: cg1_active(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
  real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol)
  type(cprj_type),intent(in) :: cprj(dtset%natom,dtset%nspinor*dtset%mband*mkmem*dtset%nsppol*usecprj)
  type(cprj_type),intent(in) :: cprjq(dtset%natom,dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol*usecprj)
  real(dp),intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
  real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(out) :: d2ovl(2,3,mpert,3,mpert*psps%usepaw)
  integer,intent(in) :: dimcprj(dtset%natom*psps%usepaw)
  real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(out) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(out) :: gh0c1_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
  real(dp),intent(out) :: gh1c_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
  integer,intent(in) :: indkpt1(nkpt_rbz)
  integer,intent(in) :: indsy1(4,nsym1,dtset%natom)
  integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  real(dp),intent(in) :: kxc(nfftf,nkxc)
  integer,intent(in) :: nattyp(psps%ntypat)
  integer,intent(in) :: nband_rbz(nkpt_rbz*dtset%nsppol)
  integer,intent(in) :: npwar1(nkpt_rbz)
  integer,intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
  real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom*psps%usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
  real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out) :: resid(dtset%mband*nkpt_rbz*nspden)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhog1(2,nfftf)
  real(dp),intent(in) :: rhor(nfftf,nspden)
  real(dp),intent(inout) :: rhor1(cplex*nfftf,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: symrl1(3,3,nsym1)
  real(dp),intent(in) :: vpsp1(cplex*nfftf)
  real(dp),target,intent(inout) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm*useylmgr)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm*useylmgr1)
 end subroutine scfcv3
end interface

interface
 subroutine smeared_delta(eigen0,eigenq,esmear,mband,smdelta,smdfunc)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: smdelta
  real(dp),intent(in) :: esmear
  real(dp),intent(in) :: eigen0(mband)
  real(dp),intent(in) :: eigenq(mband)
  real(dp),intent(out) :: smdfunc(mband,mband)
 end subroutine smeared_delta
end interface

interface
 subroutine sydy3(cplex,dyfrow,indsym,natom,nondiag,nsym,qphon,sdyfro,symq,symrec)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: natom
  integer,intent(in) :: nondiag
  integer,intent(in) :: nsym
  real(dp),intent(in) :: dyfrow(cplex,3,3,natom,1+(natom-1)*nondiag)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out) :: sdyfro(cplex,3,3,natom,1+(natom-1)*nondiag)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine sydy3
end interface

interface
 subroutine sygra3(natom,desym,deunsy,indsym,ipert,nsym,qpt,symrec)
  use defs_basis
  implicit none
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(out) :: desym(2,3,natom)
  real(dp),intent(in) :: deunsy(2,3,natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: qpt(3)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine sygra3
end interface

interface
 subroutine symdyma(dmati,indsym,natom,nsym,qptn,rprimd,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(inout) :: dmati(2*3*natom*3*natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: qptn(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symdyma
end interface

interface
 subroutine symph3(iout,acell,eigvec,indsym,natom,nsym,phfrq,rprim,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: eigvec(2*3*natom*3*natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: rprim(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symph3
end interface

interface
 subroutine syper3(indsym,mpert,natom,nsym,pertsy,rfdir,rfpert,symq,symrec,symrel)
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: rfdir(3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: pertsy(3,mpert)
  integer,intent(in) :: rfpert(mpert)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine syper3
end interface

interface
 subroutine sytens(indsym,mpert,natom,nsym,rfpert,symrec,symrel)
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(inout) :: rfpert(3,mpert,3,mpert,3,mpert)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine sytens
end interface

interface
 subroutine vloca3(atindx,cplex,gmet,gsqcut,idir,ipert,&  
  &  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&  
  &  ntypat,n1,n2,n3,paral_kgb,ph1d,qgrid,qphon,ucvol,vlspl,vpsp1,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(out) :: vpsp1(cplex*nfft)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine vloca3
end interface

interface
 subroutine vlocalstr(gmet,gprimd,gsqcut,istr,mgfft,mpi_enreg,&  
  &  mqgrid,natom,nattyp,nfft,ngfft,ntypat,paral_kgb,ph1d,qgrid,&  
  &  ucvol,vlspl,vpsp1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istr
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(out) :: vpsp1(nfft)
 end subroutine vlocalstr
end interface

interface
 subroutine vtorho3(atindx,atindx1,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,cpus,dbl_nnsclo,&  
  &  dimcprj,dim_eig2rf,doccde_rbz,docckqde,dtefield,dtfil,dtset,&  
  &  edocc,eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,enl0,enl1,&  
  &  fermie1,gh0c1_set,gh1c_set,gmet,gprimd,hdr,idir,indsy1,&  
  &  ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,mband,&  
  &  mkmem,mkqmem,mk1mem,mpi_enreg,mpsang,mpw,mpw1,&  
  &  natom,nband_rbz,ncpgr,nfftf,nhat1,nkpt_rbz,npwarr,npwar1,nres2,nspden,&  
  &  nsppol,nsym1,ntypat,nvresid1,occkq,occ_rbz,optres,&  
  &  paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,pawrhoij1,pawtab,&  
  &  phnons1,ph1d,prtvol,psps,pwindall,qmat,resid,residm,rhog1,rhor1,rmet,rprimd,symaf1,symrc1,symrl1,ucvol,&  
  &  usecprj,useylmgr1,wffddk,wffnew,wffnow,wfftgs,wfftkq,vtrial,vtrial1,wtk_rbz,xred,ylm,ylm1,ylmgr1)
  use defs_basis
  use defs_abitypes
  use m_efield
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: dbl_nnsclo
  integer,intent(in) :: dim_eig2rf
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mband
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkqmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt_rbz
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym1
  integer,intent(in) :: ntypat
  integer,intent(in) :: optres
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  integer,intent(in) :: useylmgr1
  real(dp),intent(in) :: cpus
  type(efield_type) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: edocc
  real(dp),intent(out) :: eeig0
  real(dp),intent(out) :: ek0
  real(dp),intent(out) :: ek1
  real(dp),intent(out) :: eloc0
  real(dp),intent(out) :: enl0
  real(dp),intent(out) :: enl1
  real(dp),intent(in) :: fermie1
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: nres2
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in) :: pawang1
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(out) :: residm
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffddk
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  type(wffile_type),intent(inout) :: wfftkq
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)
  real(dp),intent(inout) :: cg1(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: cg1_active(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*mband*mkqmem*nsppol)
  type(cprj_type),intent(in) :: cprj(natom,dtset%nspinor*dtset%mband*mkmem *dtset%nsppol*usecprj)
  type(cprj_type),intent(inout) :: cprj1(natom,dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*usecprj)
  type(cprj_type),intent(in) :: cprjq(natom,dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol*usecprj)
  integer,intent(in) :: dimcprj(natom*psps%usepaw)
  real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
  real(dp),intent(out) :: eigen1(2*mband*mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol)
  real(dp),intent(out) :: gh0c1_set(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(out) :: gh1c_set(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsy1(4,nsym1,natom)
  integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk_rbz(nkpt_rbz)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg1(3,mpw1*mk1mem)
  real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
  integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol)
  real(dp), intent(out) :: nhat1(cplex*nfftf,dtset%nspden*psps%usepaw)
  integer,intent(in) :: npwar1(nkpt_rbz,2)
  integer,intent(in) :: npwarr(nkpt_rbz,2)
  real(dp),intent(inout) :: nvresid1(cplex*nfftf,nspden)
  real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
  real(dp),intent(in) :: occkq(mband*nkpt_rbz*nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
  type(paw_ij_type),intent(in) :: paw_ij1(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
  real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
  real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt_rbz,2,3)
  real(dp),intent(out) :: resid(mband*nkpt_rbz*nsppol)
  real(dp),intent(out) :: rhog1(2,nfftf)
  real(dp),intent(inout) :: rhor1(cplex*nfftf,nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symaf1(nsym1)
  integer,intent(in) :: symrc1(3,3,nsym1)
  integer,intent(in) :: symrl1(3,3,nsym1)
  real(dp),intent(inout) :: vtrial(nfftf,nspden)
  real(dp),intent(inout) :: vtrial1(cplex*nfftf,nspden)
  real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm*useylmgr1)
 end subroutine vtorho3
end interface

interface
 subroutine vtowfk3(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,cpus,&  
  &  dimcprj,dimekb,dime1kb,dim_eig2rf,dimffnlk,dimffnl1,dimphkxred,dkinpw,dtfil,dtset,&  
  &  edocc_k,eeig0_k,eig0_k,eig0_kq,eig1_k,ekb_typ,e1kbfr,e1kbsc,&  
  &  ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k,&  
  &  fermie1,ffnlk,ffnlkq,ffnl1,gbound,gh0c1_set,gh1c_set,grad_berry,gs_hamkq,&  
  &  ibg,ibgq,ibg1,icg,icgq,icg1,idir,ikpt,indlmn_typ,ipert,&  
  &  isppol,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lmnmax,matblk,mband,mcgq,mcprjq,mgfft,mkmem,mk1mem,&  
  &  mpi_enreg,mpsang,mpssoang,mpw,mpw1,natom,nband_k,ncpgr,&  
  &  nkpg,nkpg1,nkpt,nnsclo_now,npw_k,npw1_k,nspinor,nsppol,&  
  &  ntypat,n4,n5,n6,occ_k,pawrhoij1,ph3d,phkxred,prtvol,psps,resid_k,rhoaug1,rocceig,&  
  &  sij_typ,usecprj,usee1kb,wffddk,wffnew,wffnow,wfftgs,vlocal,vlocal1,wtk_k)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: dim_eig2rf
  integer,intent(in) :: dime1kb
  integer,intent(in) :: dimekb
  integer,intent(in) :: dimffnl1
  integer,intent(in) :: dimffnlk
  integer,intent(in) :: dimphkxred
  integer,intent(in) :: ibg
  integer,intent(in) :: ibg1
  integer,intent(in) :: ibgq
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: icgq
  integer,intent(in) :: idir
  integer,intent(in) :: ikpt
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mband
  integer,intent(in) :: mcgq
  integer,intent(in) :: mcprjq
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw
  integer,intent(in) :: mpw1
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(inout) :: nband_k
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nkpg
  integer,intent(in) :: nkpg1
  integer,intent(in) :: nkpt
  integer,intent(in) :: nnsclo_now
  integer,intent(inout) :: npw1_k
  integer,intent(inout) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  integer,intent(in) :: usee1kb
  real(dp),intent(in) :: cpus
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie1
  type(gs_hamiltonian_type),intent(in) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffddk
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wffile_type),intent(inout) :: wfftgs
  real(dp),intent(in) :: wtk_k
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(inout) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
  real(dp),intent(out) :: cg1_active(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(in) :: cgq(2,mcgq)
  type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
  type(cprj_type),intent(out) :: cprj1(natom,nspinor*mband*mk1mem*nsppol*usecprj)
  type(cprj_type),intent(in) :: cprjq(natom,mcprjq)
  integer,intent(in) :: dimcprj(natom*gs_hamkq%usepaw)
  real(dp),intent(in) :: dkinpw(npw_k)
  real(dp),intent(in) :: e1kbfr(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
  real(dp),intent(in) :: e1kbsc(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
  real(dp),intent(out) :: edocc_k(nband_k)
  real(dp),intent(out) :: eeig0_k(nband_k)
  real(dp),intent(in) :: eig0_k(nband_k)
  real(dp),intent(in) :: eig0_kq(nband_k)
  real(dp),intent(out) :: eig1_k(2*nband_k**2)
  real(dp),intent(out) :: ek0_k(nband_k)
  real(dp),intent(out) :: ek1_k(nband_k)
  real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
  real(dp),intent(out) :: eloc0_k(nband_k)
  real(dp),intent(out) :: enl0_k(nband_k)
  real(dp),intent(out) :: enl1_k(nband_k)
  real(dp),intent(in) :: ffnl1(npw1_k,dimffnl1,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlk(npw_k,dimffnlk,lmnmax,1+gs_hamkq%usepaw*(ntypat-1))
  real(dp),intent(in) :: ffnlkq(npw1_k,dimffnl1,lmnmax,1)
  integer,intent(in) :: gbound(2*mgfft+8,2)
  real(dp),intent(out) :: gh0c1_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(out) :: gh1c_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
  real(dp),intent(in) :: grad_berry(2,mpw1,nband_k)
  integer,intent(in) :: indlmn_typ(6,lmnmax,1)
  integer,intent(in) :: kg1_k(3,npw1_k)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kinpw1(npw1_k)
  real(dp),intent(in) :: kpg1_k(npw1_k,nkpg1)
  real(dp),intent(in) :: kpg_k(npw_k,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: occ_k(nband_k)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(natom*gs_hamkq%usepaw)
  real(dp),intent(inout) :: ph3d(2,npw1_k,matblk)
  real(dp),intent(in) :: phkxred(2,dimphkxred)
  real(dp),intent(out) :: resid_k(nband_k)
  real(dp),intent(inout) :: rhoaug1(cplex*n4,n5,n6)
  real(dp),intent(in) :: rocceig(nband_k,nband_k)
  real(dp),intent(inout) :: sij_typ(dimekb,gs_hamkq%usepaw)
  real(dp),intent(inout) :: vlocal(n4,n5,n6)
  real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6)
 end subroutine vtowfk3
end interface

interface
 subroutine wfkfermi3(cg,cgq,cplex,dimekb,dimffnlk,dimffnl1,dkinpw,dtfil,dtset,&  
  &  eig1_k,ekb_typ,fe1fixed_k,fe1norm_k,ffnlk,ffnlkq,ffnl1,gbound,gs_hamkq,&  
  &  icg,icgq,idir,ikpt,indlmn_typ,ipert,isppol,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,&  
  &  lmnmax,matblk,mband,mcgq,mgfft,mkmem,mpi_enreg,&  
  &  mpsang,mpssoang,mpw,natom,nband_k,nkpg,nkpg1,&  
  &  npw_k,npw1_k,nsppol,ntypat,n4,n5,n6,occ_k,ph3d,prtvol,&  
  &  rhoaug1,rocceig,wfftgs,wtk_k)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: dimekb
  integer,intent(in) :: dimffnl1
  integer,intent(in) :: dimffnlk
  integer,intent(in) :: icg
  integer,intent(in) :: icgq
  integer,intent(in) :: idir
  integer,intent(in) :: ikpt
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mband
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(inout) :: nband_k
  integer,intent(in) :: nkpg
  integer,intent(in) :: nkpg1
  integer,intent(in) :: npw1_k
  integer,intent(inout) :: npw_k
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(in) :: gs_hamkq
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wfftgs
  real(dp),intent(in) :: wtk_k
  real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: cgq(2,mcgq)
  real(dp),intent(in) :: dkinpw(npw_k)
  real(dp),intent(out) :: eig1_k(2*nband_k**2)
  real(dp),intent(in) :: ekb_typ(dimekb,1,dtset%nspinor**2)
  real(dp),intent(out) :: fe1fixed_k(nband_k)
  real(dp),intent(out) :: fe1norm_k(nband_k)
  real(dp),intent(in) :: ffnl1(npw1_k,dimffnl1,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlk(npw_k,dimffnlk,lmnmax,1)
  real(dp),intent(in) :: ffnlkq(npw1_k,dimffnl1,lmnmax,1)
  integer,intent(in) :: gbound(2*mgfft+8,2)
  integer,intent(in) :: indlmn_typ(6,lmnmax,1)
  integer,intent(in) :: kg1_k(3,npw1_k)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kinpw1(npw1_k)
  real(dp),intent(in) :: kpg1_k(npw1_k,nkpg1)
  real(dp),intent(in) :: kpg_k(npw_k,nkpg)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(inout) :: ph3d(2,npw1_k,matblk)
  real(dp),intent(inout) :: rhoaug1(cplex*n4,n5,n6)
  real(dp),intent(in) :: rocceig(nband_k,nband_k)
 end subroutine wfkfermi3
end interface

interface
 subroutine wings3(carflg,d2cart,mpert)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(inout) :: carflg(3,mpert,3,mpert)
  real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)
 end subroutine wings3
end interface

interface
 subroutine write_blok8(ddb_blk,iblok,&  
  &  choice,mband,mpert,msize,nkpt,nunit,&  
  &  blkval2,kpt) !optional
  use defs_basis
  use m_ddb_blk
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: iblok
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: nkpt
  integer,intent(in) :: nunit
  type(ddb_blk_type), pointer :: ddb_blk
  real(dp),intent(in),optional :: blkval2(2,msize,mband,nkpt)
  real(dp),intent(in),optional :: kpt(3,nkpt)
 end subroutine write_blok8
end interface

interface
 subroutine wrtloctens(blkflg,d2bbb,d2nl,mband,mpert,natom,prtbbb,rprimd,usepaw)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: prtbbb
  integer,intent(in) :: usepaw
  integer,intent(in) :: blkflg(3,mpert,3,mpert)
  real(dp),intent(inout) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
  real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine wrtloctens
end interface

end module interfaces_72_response
!!***
