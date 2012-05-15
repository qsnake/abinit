!!****m* ABINIT/interfaces_67_common
!! NAME
!! interfaces_67_common
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/67_common
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

module interfaces_67_common

 implicit none

interface
 subroutine berry_linemin(bcut,chc,cg,cgq,costh,cprj_k,detovc,dphase_aux1,dhc,dhd,&  
  &  dimffnl,dimlmn,dimlmn_srt,direc,dtefield,ffnl,gen_eigenpb,&  
  &  gs_hamk,hel,iband,ikpt,iline,isppol,kg_k,lmnmax,matblk,&  
  &  mband,mcg,mcgq,mgfft,mkgq,mpi_enreg,mpw,&  
  &  natom,nkpt,nloalg,npw,npwarr,nspinor,ntypat,pwind,ph3d,&  
  &  phase_init,phase_end,pwind_alloc,pwnsfac,pwnsfacq,sinth,&  
  &  spaceComm_distrb,thetam,xnorm)
  use defs_basis
  use defs_abitypes
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimffnl
  integer,intent(in) :: iband
  integer,intent(in) :: ikpt
  integer,intent(in) :: iline
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkgq
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: spaceComm_distrb
  real(dp),intent(in) :: chc
  real(dp),intent(out) :: costh
  real(dp),intent(in) :: dhc
  real(dp),intent(in) :: dhd
  type(efield_type),intent(in) :: dtefield
  logical,intent(in) :: gen_eigenpb
  type(gs_hamiltonian_type),intent(in) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: sinth
  real(dp),intent(out) :: thetam
  real(dp),intent(in) :: xnorm
  integer,intent(out) :: hel(2,3)
  integer,intent(in) :: nloalg(5)
  real(dp),intent(out) :: bcut(2,3)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: cgq(2,mcgq)
  type(cprj_type),intent(in) :: cprj_k(natom,dtefield%nband_occ*gs_hamk%usepaw)
  real(dp),intent(in) :: detovc(2,2,3)
  integer,intent(in) :: dimlmn(natom)
  integer,intent(in) :: dimlmn_srt(natom)
  real(dp),intent(inout) :: direc(2,npw*nspinor)
  real(dp),intent(inout) :: dphase_aux1(3)
  real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  integer,intent(in) :: kg_k(3,npw)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(inout) :: ph3d(2,npw,matblk)
  real(dp),intent(out) :: phase_end(3)
  real(dp),intent(inout) :: phase_init(3)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: pwnsfacq(2,mkgq)
 end subroutine berry_linemin
end interface

interface
 subroutine berryphase(atindx1,bdberry,cg,gprimd,istwfk,kberry,kg,kpt_,&  
  &  kptopt,kptrlatt,mband,mcg,&  
  &  mkmem,mpi_enreg,mpw,natom,nattyp,nband,nberry,npwarr,nspinor,nsppol,ntypat,&  
  &  nkpt_,rprimd,ucvol,unkg,wffnow,xred,zion)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nberry
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: unkg
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: bdberry(4)
  integer,intent(in) :: kptrlatt(3,3)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: gprimd(1:3,1:3)
  integer,intent(in) :: istwfk(nkpt_)
  integer,intent(in) :: kberry(3,nberry)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt_(3,nkpt_)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt_*nsppol)
  integer,intent(in) :: npwarr(nkpt_)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine berryphase
end interface

interface
 subroutine berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,&  
  &  gprimd,hdr,indlmn,kg,lmnmax,mband,mcg,mcprj,&  
  &  mkmem,mpi_enreg,mpw,natom,npwarr,nsppol,ntypat,&  
  &  nkpt,option,pawrhoij,pawtab,pel,pelev,pion,pwind,&  
  &  pwind_alloc,pwnsfac,&  
  &  rprimd,typat,ucvol,unit_out,usecprj,usepaw,wffnow,xred,zion)
  use defs_basis
  use defs_abitypes
  use m_efield
  use defs_datatypes
  use m_wffile
  implicit none
  integer, intent(in) :: lmnmax
  integer, intent(in) :: mband
  integer, intent(in) :: mcg
  integer, intent(in) :: mcprj
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: natom
  integer, intent(in) :: nkpt
  integer, intent(in) :: nsppol
  integer, intent(in) :: ntypat
  integer, intent(in) :: option
  integer, intent(in) :: pwind_alloc
  integer, intent(in) :: unit_out
  integer, intent(in) :: usecprj
  integer, intent(in) :: usepaw
  type(efield_type), intent(inout) :: dtefield
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(in) :: ucvol
  type(wffile_type), intent(inout) :: wffnow
  integer, intent(in) :: atindx1(natom)
  real(dp), intent(in) :: cg(2,mcg)
  type(cprj_type),intent(in) :: cprj(natom,mcprj*usecprj)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: indlmn(6,lmnmax,ntypat)
  integer, intent(in) :: kg(3,mpw*mkmem)
  integer, intent(in) :: npwarr(nkpt)
  type(pawrhoij_type), intent(in) :: pawrhoij(natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp), intent(out) :: pel(3)
  real(dp), intent(out) :: pelev(3)
  real(dp), intent(out) :: pion(3)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: typat(natom)
  real(dp), intent(inout) :: xred(3,natom)
  real(dp), intent(in) :: zion(ntypat)
 end subroutine berryphase_new
end interface

interface
 subroutine calc_cs(cg,dtefield,gprimd,kg,kptns,mcg,mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,&  
  &  occopt,usepaw)
  use defs_basis
  use m_efield
  use defs_abitypes
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: occopt
  integer,intent(in) :: usepaw
  type(efield_type),intent(in) :: dtefield
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: npwarr(nkpt)
 end subroutine calc_cs
end interface

interface
 subroutine calc_efg(gprimd,natom,nfft,ngfft,nspden,ntypat,paral_kgb,&  
  &  paw_an,pawang,pawrad,pawrhoij,pawtab,&  
  &  ptcharge,prtefg,quadmom,rhor,rprimd,typat,ucvol,usepaw,xred,zion)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtefg
  integer,intent(in) :: usepaw
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  type(paw_an_type),intent(in) :: paw_an(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ptcharge(ntypat)
  real(dp),intent(in) :: quadmom(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine calc_efg
end interface

interface
 subroutine calc_fc(natom,ntypat,pawrad,pawrhoij,pawtab,psps,typat)
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pseudopotential_type),intent(in) :: psps
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine calc_fc
end interface

interface
 subroutine cgwf(berryopt,cg,cgq,chkexit,cpus,dimffnl,dphase_k,dtefield,&  
  &  ffnl,filnam_ds1,filstat,&  
  &  gsc,gs_hamk,icg,igsc,ikpt,inonsc,&  
  &  isppol,kg_k,kinpw,lmnmax,matblk,mband,&  
  &  mcg,mcgq,mgfft,mgsc,mkgq,mpi_enreg,mpsang,&  
  &  mpssoang,mpw,natom,nband,nbdblock,nkpt,nline,nloalg,npw,npwarr,&  
  &  nspinor,nsppol,ntypat,nvloc,n4,n5,n6,ortalg,&  
  &  paral_kgb,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,&  
  &  pwnsfacq,quit,resid,subham,subovl,subvnl,tolwfr,&  
  &  use_subovl,vlocal,wfoptalg,zshift,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: chkexit
  integer,intent(in) :: dimffnl
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: ikpt
  integer,intent(in) :: inonsc
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgsc
  integer,intent(in) :: mkgq
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nbdblock
  integer,intent(in) :: nkpt
  integer,intent(in) :: nline
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: nvloc
  integer,intent(in) :: ortalg
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: pwind_alloc
  integer,intent(inout) :: quit
  integer,intent(in) :: use_subovl
  integer,intent(in) :: wfoptalg
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  character(len=fnlen),intent(in) :: filnam_ds1
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: tolwfr
  integer,intent(in) :: nloalg(5)
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: cgq(2,mcgq)
  real(dp),intent(out) :: dphase_k(3)
  real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  real(dp),intent(inout) :: gsc(2,mgsc)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw(npw)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(inout) :: ph3d(2,npw,matblk)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: pwnsfacq(2,mkgq)
  real(dp),intent(out) :: resid(nband)
  real(dp),intent(out) :: subham(nband*(nband+1))
  real(dp),intent(out) :: subovl(nband*(nband+1)*use_subovl)
  real(dp),intent(out) :: subvnl(nband*(nband+1)*(1-gs_hamk%usepaw))
  real(dp),intent(inout) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
  real(dp),intent(in) :: zshift(nband)
 end subroutine cgwf
end interface

interface
 subroutine clnup1(acell,dtset,eigen,fermie,&  
  &  fnameabo_dos,fnameabo_eig,fred,&  
  &  mpi_enreg,nfft,ngfft,occ,prtfor,&  
  &  resid,rhor,rprimd,vxcavg,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: prtfor
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: fnameabo_dos
  character(len=fnlen),intent(in) :: fnameabo_eig
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: fred(3,dtset%natom)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine clnup1
end interface

interface
 subroutine clnup2(n1xccc,fred,gresid,grewtn,grxc,iscf,natom,prtfor,prtstr,&  
  &  prtvol,start,strten,synlgr,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iscf
  integer,intent(in) :: n1xccc
  integer,intent(in) :: natom
  integer,intent(in) :: prtfor
  integer,intent(in) :: prtstr
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: gresid(3,natom)
  real(dp),intent(in) :: grewtn(3,natom)
  real(dp),intent(in) :: grxc(3,natom)
  real(dp),intent(in) :: start(3,natom)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(in) :: synlgr(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine clnup2
end interface

interface
 subroutine conducti_nc(filnam,filnam_out,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen) :: filnam
  character(len=fnlen) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine conducti_nc
end interface

interface
 subroutine conducti_paw(filnam,filnam_out,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen) :: filnam
  character(len=fnlen) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine conducti_paw
end interface

interface
 subroutine conducti_paw_core(filnam,filnam_out,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen) :: filnam
  character(len=fnlen) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine conducti_paw_core
end interface

interface
 subroutine constrf(diffor,fcart,forold,fred,iatfix,ionmov,maxfor,natom,&  
  &  nconeq,prtvol,rprimd,wtatcon,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: ionmov
  integer,intent(in) :: natom
  integer,intent(in) :: nconeq
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: diffor
  real(dp),intent(out) :: maxfor
  real(dp),intent(inout) :: fcart(3,natom)
  real(dp),intent(inout) :: forold(3,natom)
  real(dp),intent(out) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtatcon(3,natom,nconeq)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine constrf
end interface

interface
 subroutine dens_in_sph(cmax,cg,gmet,istwfk,kg_k,natom,ngfft,mpi_enreg,npw_k,&  
  &  paral_kgb,ph1d,rmax,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwfk
  integer,intent(in) :: natom
  integer,intent(in) :: npw_k
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: cg(2,npw_k)
  real(dp),intent(out) :: cmax(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: ph1d(2,(2*ngfft(1)+1+2*ngfft(2)+1+2*ngfft(3)+1)*natom)
  real(dp),intent(in) :: rmax(natom)
 end subroutine dens_in_sph
end interface

interface
 subroutine dielmt(dielinv,gmet,kg_diel,&  
  &  npwdiel,nspden,occopt,prtvol,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dielmt
end interface

interface
 subroutine dielmt2(gmet,kg_diel,&  
  &  npwdiel,nspden,occopt,susmat,&  
  &  dieldiag,dtset,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: dieldiag(2,npwdiel,nspden)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dielmt2
end interface

interface
 subroutine dieltcel(dielinv,gmet,kg_diel,kxc,&  
  &  nfft,ngfft,nkxc,npwdiel,nspden,occopt,option,paral_kgb,prtvol,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dieltcel
end interface

interface
 subroutine emispec_paw(filnam,filnam_out,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen) :: filnam
  character(len=fnlen) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine emispec_paw
end interface

interface
 subroutine energy(atindx,atindx1,cg,compch_fft,dtfil,dtset,electronpositron,&  
  &  energies,eigen,etotal,gsqcut,indsym,irrzon,kg,mcg,mpi_enreg,nattyp,nfftf,ngfftf,nhat,&  
  &  nhatgr,nhatgrdim,npwarr,n3xccc,occ,optene,paw_ij,pawang,pawfgr,&  
  &  pawfgrtab,pawrhoij,pawtab,phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,&  
  &  taug,taur,usexcnhat,vhartr,vtrial,vpsp,vxc,wffnow,wfs,wvl,xccc3d,xred,ylm,&  
  &  vxctau) ! optional arguments
  use defs_wvltypes
  use m_energies
  use defs_abitypes
  use m_wffile
  use defs_basis
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: optene
  integer,intent(in) :: usexcnhat
  real(dp),intent(out) :: compch_fft
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_wf_type),intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  integer, intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  real(dp), intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type), intent(in) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(out) :: strsxc(6)
  integer, intent(in) :: symrec(3,3,dtset%nsym)
  real(dp), intent(inout) :: taug(2,nfftf*dtset%usekden)
  real(dp), intent(inout) :: taur(nfftf,dtset%nspden*dtset%usekden)
  real(dp), intent(out) :: vhartr(nfftf)
  real(dp), intent(in) :: vpsp(nfftf)
  real(dp), intent(out) :: vtrial(nfftf,dtset%nspden)
  real(dp), intent(out) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(out),optional :: vxctau(nfftf,dtset%nspden*dtset%usekden,4)
  real(dp), intent(in) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine energy
end interface

interface
 subroutine etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&  
  &  hel,nkpt,nstr,sdeg,theta)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: chc
  real(dp),intent(in) :: dhc
  real(dp),intent(in) :: dhd
  real(dp),intent(out) :: e0
  real(dp),intent(out) :: e1
  real(dp),intent(in) :: sdeg
  real(dp),intent(in) :: theta
  integer,intent(in) :: hel(2,3)
  integer,intent(in) :: nstr(3)
  real(dp),intent(in) :: bcut(2,3)
  real(dp),intent(in) :: detovc(2,2,3)
  real(dp),intent(in) :: detovd(2,2,3)
  real(dp),intent(in) :: efield_dot(3)
 end subroutine etheta
end interface

interface
 subroutine etotfor(atindx1,deltae,diffor,dtset,efield_dot,&  
  &  elast,electronpositron,energies,&  
  &  etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,&  
  &  grxc,gsqcut,indsym,kxc,mag_cart,maxfor,mgfft,mpi_enreg,nattyp,&  
  &  nfft,ngfft,nhat,nkxc,ntypat,nvresid,n1xccc,n3xccc,optene,optforces,optres,&  
  &  pawang,pawfgrtab,pawrhoij,pawtab,pel,ph1d,pion,psps,rhog,rhor,rprimd,symrec,synlgr,&  
  &  usepaw,usexcnhat,vhartr,vpsp,vxc,wvl,xccc3d,xred)
  use m_energies
  use defs_abitypes
  use m_electronpositron
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: optforces
  integer,intent(in) :: optres
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(out) :: deltae
  real(dp),intent(out) :: diffor
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: elast
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  real(dp),intent(out) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: efield_dot(3)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(out) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),intent(inout) :: grnl(3*dtset%natom)
  real(dp),intent(out) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: mag_cart(3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*psps%usepaw)
  real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine etotfor
end interface

interface
 subroutine evdw_wannier(csix,corrvdw,mwan,natom,nsppol,nwan,vdw_nfrag,&  
  &  vdw_supercell,vdw_typfrag,vdw_xc,rprimd,wann_centres,wann_spreads,xcart)    
  use defs_basis
  implicit none
  integer , intent(in) :: mwan
  integer , intent(in) :: natom
  integer , intent(in) :: nsppol
  integer , intent(in) :: vdw_nfrag
  integer , intent(in) :: vdw_xc
  real(dp), intent(out) :: corrvdw
  integer , intent(in) :: vdw_supercell(3)
  real(dp), intent(out), allocatable :: csix(:,:,:,:)
  integer , intent(in) :: nwan(nsppol)
  real(dp), intent(in) :: rprimd(3,3)
  integer , intent(in) :: vdw_typfrag(natom)
  real(dp), intent(in) :: wann_centres(3,mwan,nsppol)
  real(dp), intent(in) :: wann_spreads(mwan,nsppol)
  real(dp), intent(in) :: xcart(3,natom)
 end subroutine evdw_wannier
end interface

interface
 subroutine getFu(sn,sl,rn,rl,fu) ! sn-->spread(n), sl-->spread(l), rn --> rc(n), rl --> rc(l) 
  use defs_basis
  implicit none
  real(dp),intent(out) :: fu
  real(dp),intent(in) :: rl
  real(dp),intent(in) :: rn
  real(dp),intent(in) :: sl
  real(dp),intent(in) :: sn
 end subroutine getFu
end interface

interface
 subroutine order_wannier(mwan,natom,nwan,nsppol,ord,vdw_typfrag,wanncent,xcart)
  use defs_basis
  implicit none
  integer, intent(in) :: mwan
  integer, intent(in) :: natom
  integer, intent(in) :: nsppol
  integer, intent(in) :: nwan(nsppol)
  integer, intent(inout) :: ord(mwan,nsppol)
  integer, intent(in) :: vdw_typfrag(natom)
  real(dp),intent(in) :: wanncent(3,mwan,nsppol)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine order_wannier
end interface

interface
 subroutine ovlp_wann(mwan,nwan,nsppol,ord,wanncent,wannspr,xi)
  use defs_basis
  implicit none
  integer, intent(in) :: mwan
  integer, intent(in) :: nsppol
  integer, intent(in) :: nwan(nsppol)
  integer, intent(in) :: ord(mwan,nsppol)
  real(dp),intent(in) :: wanncent(3,mwan,nsppol)
  real(dp),intent(in) :: wannspr(mwan,nsppol)
  real(dp), intent(out) :: xi(mwan,nsppol)
 end subroutine ovlp_wann
end interface

interface
 subroutine ewald(eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: eew
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(out) :: grewtn(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine ewald
end interface

interface
 subroutine ewald2(gmet,natom,ntypat,rmet,rprimd,stress,&  
  &  typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: stress(6)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine ewald2
end interface

interface
 subroutine extraprho(atindx,atindx1,cg,dtset,gmet,gprimd,gsqcut,istep,&  
  &  kg,mcg,mgfft,mpi_enreg,mqgrid,nattyp,nfft,ngfft,npwarr,ntypat,pawrhoij,&  
  &  pawtab,ph1d,psps,qgrid,rhor,rprimd,scf_history,ucvol,usepaw,&  
  &  xred_new,xred_old,ylm,zion,znucl)
  use m_scf_history
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(scf_history_type),intent(inout) :: scf_history
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: npwarr(dtset%nkpt)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred_new(3,dtset%natom)
  real(dp),intent(in) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine extraprho
end interface

interface
 subroutine extrapwf(atindx,atindx1,cg,dtset,istep,kg,mcg,mgfft,mpi_enreg,&  
  &  nattyp,ngfft,npwarr,ntypat,pawtab,psps,scf_history,usepaw,xred_old,ylm)
  use m_scf_history
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(scf_history_type),intent(inout) :: scf_history
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,mcg)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: npwarr(dtset%nkpt)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine extrapwf
end interface

interface
 subroutine fconv(fcart,iatfix,iexit,itime,natom,ntime,&  
  &  optcell,strfact,strtarget,strten,tolmxf)
  use defs_basis
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: natom
  integer,intent(in) :: ntime
  integer,intent(in) :: optcell
  real(dp),intent(in) :: strfact
  real(dp),intent(in) :: tolmxf
  real(dp),intent(in) :: fcart(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
 end subroutine fconv
end interface

interface
 subroutine filterpot(cplex,gmet,gsqcut,nfft,ngfft,nspden,paral_kgb,qphon,vpot)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(inout) :: vpot(cplex*nfft,nspden)
 end subroutine filterpot
end interface

interface
 subroutine fixsym(iatfix,indsym,natom,nsym)
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: iatfix(3,natom)
  integer,intent(in) :: indsym(4,nsym,natom)
 end subroutine fixsym
end interface

interface
 subroutine forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&  
  &  grhf,grnl,grxc,gsqcut,indsym,&  
  &  maxfor,mgfft,mpi_enreg,n1xccc,n3xccc,&  
  &  nattyp,nfft,ngfft,ntypat,&  
  &  pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,&  
  &  vresid,vxc,wvl,xred,&  
  &  electronpositron) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_electronpositron
  use defs_wvltypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: diffor
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(in) :: gsqcut
  real(dp),intent(out) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(inout) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),intent(in) :: grnl(3*dtset%natom)
  real(dp),intent(out) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(inout) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine forces
end interface

interface
 subroutine forstr(atindx1,cg,diffor,dtefield,dtset,eigen,electronpositron,energies,favg,fcart,&  
  &  forold,fred,gresid,grewtn,grhf,grxc,gsqcut,indsym,&  
  &  kg,kxc,maxfor,mcg,mgfftf,mpi_enreg,n3xccc,nattyp,&  
  &  nfftf,ngfftf,nhat,nkxc,npwarr,&  
  &  ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,pawfgr,&  
  &  pawfgrtab,pawrhoij,pawtab,pel,ph1d,ph1df,pion,psps,rhog,rhor,rprimd,stress_needed,&  
  &  strsxc,strten,symrec,synlgr,ucvol,unkg,unylm,usexcnhat,vhartr,vpsp,&  
  &  vxc,wffnow,wvl,xccc3d,xred,ylm,ylmgr)
  use defs_wvltypes
  use m_energies
  use defs_abitypes
  use m_wffile
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftf
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optfor
  integer,intent(in) :: optres
  integer,intent(in) :: stress_needed
  integer,intent(in) :: unkg
  integer,intent(in) :: unylm
  integer,intent(in) :: usexcnhat
  real(dp),intent(inout) :: diffor
  type(efield_type), intent(in) :: dtefield
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(in) :: energies
  real(dp),intent(in) :: gsqcut
  real(dp),intent(inout) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(out) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),intent(out) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp),intent(in) :: kxc(dtset%nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strsxc(6)
  real(dp),intent(out) :: strten(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),intent(in) :: vhartr(nfftf)
  real(dp),intent(in) :: vpsp(nfftf)
  real(dp),intent(in) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine forstr
end interface

interface
 subroutine forstrnps(atindx1,cg,dtefield,ecut,ecutsm,effmass,eigen,electronpositron,&  
  &  grnl,istwfk,kg,kinstr,npsstr,kpt,mband,mcg,mgfft,mkmem,mpi_enreg,mpsang,&  
  &  mpw,natom,nattyp,nband,ngfft,nkpt,nloalg,npwarr,nspinor,nsppol,nsym,&  
  &  ntypat,occ,optfor,paw_ij,pawtab,ph1d,psps,rprimd,&  
  &  stress_needed,symrec,typat,unkg,unylm,use_gpu_cuda,wffnow,wtk,xred,ylm,ylmgr)
  use defs_abitypes
  use m_wffile
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: optfor
  integer,intent(in) :: stress_needed
  integer,intent(in) :: unkg
  integer,intent(in) :: unylm
  integer,intent(in) :: use_gpu_cuda
  type(efield_type), intent(in) :: dtefield
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  type(electronpositron_type),pointer :: electronpositron
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(out) :: grnl(3*natom)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(out) :: kinstr(6)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(out) :: npsstr(6)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 end subroutine forstrnps
end interface

interface
 subroutine fresid(dtset,gresid,mpi_enreg,nfft,ngfft,ntypat,option,&  
  &  pawtab,rhor,rprimd,ucvol,work,xred_new,xred_old,znucl)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: work(nfft,dtset%nspden)
  real(dp),intent(in) :: xred_new(3,dtset%natom)
  real(dp),intent(in) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine fresid
end interface

interface
 subroutine fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,mpi_enreg,mqgrid,nattyp,nfft,&  
  &  ngfft,ntypat,pawtab,ph1d,qgrid,ucvol,usepaw,vresid,zion,znucl)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine fresidrsp
end interface

interface
 subroutine getcgqphase(dtset, timrev, cg,  mcg,  cgq, mcgq, mpi_enreg,&  
  &  nkpt_rbz, npwarr, npwar1, phasecg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: mcg
  integer, intent(in) :: mcgq
  integer, intent(in) :: nkpt_rbz
  integer, intent(in) :: timrev
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(in) :: mpi_enreg
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: cgq(2,mcgq)
  integer, intent(in) :: npwar1(nkpt_rbz)
  integer, intent(in) :: npwarr(nkpt_rbz)
  real(dp),intent(out) :: phasecg(2, dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)
 end subroutine getcgqphase
end interface

interface
 subroutine gipaw_j_dia_bare(jdia,nfft,ngfft,nhat,nspden,rhor,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: jdia(3,3,nfft)
  real(dp),intent(in) :: nhat(nfft,nspden)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine gipaw_j_dia_bare
end interface

interface
 subroutine hartre1(cplex,gmet,gsqcut,nfft,ngfft,paral_kgb,qphon,rhog,vhartr,&  
  &  sum,rcut,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  real(dp),intent(in) :: rcut
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: sum(3)
  real(dp),intent(out) :: vhartr(cplex*nfft)
 end subroutine hartre1
end interface

interface
 subroutine initberry(dtefield,dtset,gmet,gprimd,kg,mband,&  
  &  mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,nsppol,&  
  &  nsym,ntypat,occ,pawang,pawrad,pawtab,psps,&  
  &  pwind,pwind_alloc,pwnsfac,&  
  &  rprimd,symrec,typat,usepaw,xred)
  use defs_basis
  use m_efield
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(out) :: pwind_alloc
  integer,intent(in) :: usepaw
  type(efield_type),intent(out) :: dtefield
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,pointer :: pwind(:,:,:)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),pointer :: pwnsfac(:,:)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine initberry
end interface

interface
 subroutine initmv(cgindex,dtfil,dtset,gmet,kg,kneigh,kg_neigh,kptindex,&  
  &  kpt3,mband,mkmem,mkmem_max,mpi_enreg,mpw,nband,nkpt2,&  
  &  nkpt3,nneigh,npwarr,nsppol,occ,pwind)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkmem_max
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nkpt3
  integer,intent(in) :: nneigh
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(out) :: cgindex(nkpt2,nsppol)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg_neigh(30,nkpt2,3)
  integer,intent(in) :: kneigh(30,nkpt2)
  real(dp),intent(in) :: kpt3(3,nkpt3)
  integer,intent(in) :: kptindex(2,nkpt3)
  integer,intent(in) :: nband(nkpt2*nsppol)
  integer,intent(in) :: npwarr(nkpt2)
  real(dp),intent(in) :: occ(mband*nkpt2*nsppol)
  integer,intent(out) :: pwind(mpw,nneigh,mkmem)
 end subroutine initmv
end interface

interface
 subroutine initro(atindx,densty,gmet,gsqcut,izero,mgfft,mpi_enreg,mqgrid,natom,nattyp,&  
  &  nfft,ngfft,nspden,ntypat,paral_kgb,pawtab,ph1d,qgrid,rhog,rhor,spinat,ucvol,usepaw,zion,znucl)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: izero
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: gsqcut
  type(mpi_type) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(in) :: densty(ntypat,4)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: rhog(2,nfft)
  real(dp),intent(out) :: rhor(nfft,nspden)
  real(dp),intent(in) :: spinat(3,natom)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine initro
end interface

interface
 subroutine ionion_realSpace(dtset, eew, grewtn, rprimd, xred, zion)
  use defs_basis
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eew
  real(dp),intent(out) :: grewtn(3,dtset%natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: zion(dtset%ntypat)
 end subroutine ionion_realSpace
end interface

interface
 subroutine jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,nspden,&  
  &  option,paral_kgb,slabwsrad,rhog,rhor,rprimd,vjell,slabzstart,slabzend)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: slabwsrad
  real(dp),intent(in) :: slabzend
  real(dp),intent(in) :: slabzstart
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,min(option,nspden))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vjell(nfft)
 end subroutine jellium
end interface

interface
 subroutine jvec_to_B(cshield,gcart,jvec,natom,nfft,ngfft,paral_kgb,rprimd,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: cshield(3,3,natom)
  real(dp),intent(in) :: gcart(ngfft(1),ngfft(2),ngfft(3),3)
  real(dp),intent(in) :: jvec(3,3,nfft)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine jvec_to_B
end interface

interface
 subroutine ks_ddiago(Diago_ctl,nband_k,nfftc,mgfftc,ngfftc,natom,&  
  &  typat,nfftf,nspinor,nspden,nsppol,ntypat,Pawtab,Pawfgr,Paw_ij,&  
  &  Psps,rprimd,vtrial,xred,onband_diago,eig_ene,eig_vec,Cprj_k,comm,ierr,&  
  &  Electronpositron) ! Optional arguments
  use m_hamiltonian
  use defs_datatypes
  use m_electronpositron
  use defs_basis
  implicit none
  integer,intent(in) :: comm
  integer,intent(out) :: ierr
  integer,intent(in) :: mgfftc
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: nfftc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(out) :: onband_diago
  type(ddiago_ctl_type),intent(in) :: Diago_ctl
  type(electronpositron_type),optional,pointer :: Electronpositron
  type(pawfgr_type),intent(in) :: Pawfgr
  type(pseudopotential_type),intent(in) :: Psps
  integer,intent(in) :: ngfftc(18)
  type(cprj_type),pointer :: Cprj_k(:,:)
  type(paw_ij_type),intent(in) :: Paw_ij(natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),pointer :: eig_ene(:)
  real(dp),pointer :: eig_vec(:,:,:)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine ks_ddiago
end interface

interface
 subroutine linemin(bcut,chc,costh,detovc,detovd,dhc,dhd,dphase_aux1,&  
  &  efield_dot,iline,nkpt,nstr,hel,phase_end,phase_init,sdeg,sinth,thetam)
  use defs_basis
  implicit none
  integer,intent(in) :: iline
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: chc
  real(dp),intent(out) :: costh
  real(dp),intent(in) :: dhc
  real(dp),intent(in) :: dhd
  real(dp),intent(in) :: sdeg
  real(dp),intent(out) :: sinth
  real(dp),intent(out) :: thetam
  integer,intent(out) :: hel(2,3)
  integer,intent(in) :: nstr(3)
  real(dp),intent(out) :: bcut(2,3)
  real(dp),intent(in) :: detovc(2,2,3)
  real(dp),intent(in) :: detovd(2,2,3)
  real(dp),intent(inout) :: dphase_aux1(3)
  real(dp),intent(in) :: efield_dot(3)
  real(dp),intent(out) :: phase_end(3)
  real(dp),intent(inout) :: phase_init(3)
 end subroutine linemin
end interface

interface
 subroutine linopt(nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,occv,evalv,efermi,pmat,&  
  v1,v2,nmesh,de,sc,brod,fnam)
  use defs_basis
  implicit none
  integer, intent(in) :: nkpt
  integer, intent(in) :: nmesh
  integer, intent(in) :: nspin
  integer, intent(in) :: nstval
  integer, intent(in) :: nsymcrys
  integer, intent(in) :: v1
  integer, intent(in) :: v2
  real(dp), intent(in) :: brod
  real(dp), intent(in) :: de
  real(dp), intent(in) :: efermi
  character(256), intent(in) :: fnam
  real(dp), intent(in) :: omega
  real(dp), intent(in) :: sc
  real(dp), intent(in) :: evalv(nstval,nspin,nkpt)
  real(dp), intent(in) :: occv(nstval,nspin,nkpt)
  complex(dpc), intent(in) :: pmat(nstval,nstval,nkpt,3,nspin)
  real(dp), intent(in) :: symcrys(3,3,nsymcrys)
  real(dp), intent(in) :: wkpt(nkpt)
 end subroutine linopt
end interface

interface
 subroutine mag_out(dtefield,mpi_enreg,nkpt,rprimd,wtk)
  use defs_basis
  use m_efield
  use defs_abitypes
  implicit none
  integer,intent(in) :: nkpt
  type(efield_type),intent(in) :: dtefield
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine mag_out
end interface

interface
 subroutine magcart(dtefield,mpi_enreg,nkpt,rprimd,wtk)
  use defs_basis
  use m_efield
  use defs_abitypes
  implicit none
  integer,intent(in) :: nkpt
  type(efield_type),intent(inout) :: dtefield
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine magcart
end interface

interface
 subroutine make_efg_el(efg,gcart,natom,nfft,ngfft,nspden,paral_kgb,rhor,rprimd,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: efg(3,3,natom)
  real(dp),intent(in) :: gcart(ngfft(1),ngfft(2),ngfft(3),3)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine make_efg_el
end interface

interface
 subroutine make_efg_ion(efg,natom,ntypat,rprimd,typat,ucvol,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp) :: ucvol
  real(dp),intent(out) :: efg(3,3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine make_efg_ion
end interface

interface
 subroutine make_grad_berry(cg,cgq,cprj_k,detovc,dimffnl,dimlmn,dimlmn_srt,direc,dtefield,ffnl,grad_berry,&  
  &  gs_hamk,iband,icg,ikpt,isppol,kg_k,lmnmax,matblk,mband,&  
  &  mcg,mcgq,mgfft,mkgq,mpi_enreg,mpsang,mpssoang,mpw,natom,nkpt,&  
  &  nloalg,npw,npwarr,nspinor,nsppol,ntypat,pwind,ph3d,pwind_alloc,&  
  &  pwnsfac,pwnsfacq)
  use defs_basis
  use defs_abitypes
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimffnl
  integer,intent(in) :: iband
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkgq
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: pwind_alloc
  type(efield_type),intent(inout) :: dtefield
  type(gs_hamiltonian_type),intent(in) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nloalg(5)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: cgq(2,mcgq)
  type(cprj_type),intent(in) :: cprj_k(natom,dtefield%nband_occ*gs_hamk%usepaw)
  real(dp),intent(out) :: detovc(2,2,3)
  integer,intent(in) :: dimlmn(natom)
  integer,intent(in) :: dimlmn_srt(natom)
  real(dp),intent(inout) :: direc(2,npw*nspinor)
  real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  real(dp),intent(out) :: grad_berry(2,npw*nspinor)
  integer,intent(in) :: kg_k(3,npw)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(inout) :: ph3d(2,npw,matblk)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: pwnsfacq(2,mkgq)
 end subroutine make_grad_berry
end interface

interface
 subroutine mklocl(dtset, dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&  
  &  mpi_enreg,natom,nattyp,nfft,ngfft,nspden,ntypat,option,ph1d,psps,qprtrb,&  
  &  rhog,rhor,rprimd,ucvol,vprtrb,vpsp,wvl,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: qprtrb(3)
  real(dp),intent(out) :: dyfrlo(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grtn(3,natom)
  real(dp),intent(out) :: lpsstr(6)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vprtrb(2)
  real(dp),intent(out) :: vpsp(nfft)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine mklocl
end interface

interface
 subroutine mklocl_realspace(dtset, grtn, mpi_enreg, natom, nattyp, nfft, ngfft,&  
  &  nspden, ntypat, option,  psps, rhog, rhor,&  
  &  rprimd, ucvol, vpsp, wvl, xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: grtn(3,natom)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vpsp(nfft)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine mklocl_realspace
end interface

interface
 subroutine createIonicPotential_new(geocode,iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,gridcart,&  
  hxh,hyh,hzh,n1i,n2i,n3i,pkernel,pot_ion,spaceworld)
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: n1i
  integer, intent(in) :: n2i
  integer, intent(in) :: n3i
  integer, intent(in) :: nat
  integer, intent(in) :: nproc
  integer, intent(in) :: ntypes
  integer, intent(in) :: spaceworld
  character(len=1), intent(in) :: geocode
  real(kind=8), intent(in) :: hxh
  real(kind=8), intent(in) :: hyh
  real(kind=8), intent(in) :: hzh
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), dimension(*), intent(inout) :: pot_ion
  real(kind=8), dimension(3,n1i*n2i*n3i), intent(in) :: gridcart
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
 end subroutine createIonicPotential_new
end interface

interface
 subroutine ind_positions_(periodic,i,n,j,go)
  implicit none
  integer, intent(in) :: i
  integer, intent(out) :: j
  integer, intent(in) :: n
  logical, intent(out) :: go
  logical, intent(in) :: periodic
 end subroutine ind_positions_
end interface

interface
 subroutine local_forces_new(geocode,iproc,ntypes,nat,iatype,rxyz,gridcart,psppar,nelpsp,hxh,hyh,hzh,&  
  n1,n2,n3,rho,pot,floc)
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: nat
  integer, intent(in) :: ntypes
  character(len=1), intent(in) :: geocode
  real(kind=8), intent(in) :: hxh
  real(kind=8), intent(in) :: hyh
  real(kind=8), intent(in) :: hzh
  real(kind=8), dimension(*), intent(in) :: pot
  real(kind=8), dimension(*), intent(in) :: rho
  real(kind=8), dimension(3,nat), intent(out) :: floc
  real(kind=8), dimension(3,n1*n2*n3), intent(in) :: gridcart
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
 end subroutine local_forces_new
end interface

interface
 subroutine mklocl_recipspace(dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&  
  &  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,option,paral_kgb,ph1d,qgrid,qprtrb,&  
  &  rhog,ucvol,vlspl,vprtrb,vpsp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: qprtrb(3)
  real(dp),intent(out) :: dyfrlo(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grtn(3,natom)
  real(dp),intent(out) :: lpsstr(6)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: vprtrb(2)
  real(dp),intent(out) :: vpsp(nfft)
 end subroutine mklocl_recipspace
end interface

interface
 subroutine mklocl_wavelets(dtset, grtn, mpi_enreg, option, rhor, rprimd,&  
  &  vpsp, wvl, xcart)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  real(dp),intent(inout) :: grtn(3,dtset%natom)
  real(dp),intent(in) :: rhor(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vpsp(dtset%nfft)
  real(dp),intent(inout) :: xcart(3,dtset%natom)
 end subroutine mklocl_wavelets
end interface

interface
 subroutine mkresi(cg,dimffnl,eig_k,ffnl,filstat,gs_hamk,icg,ikpt,isppol,kg_k,kinpw,lmnmax,&  
  &  matblk,mcg,mgfft,mpi_enreg,mpsang,mpssoang,natom,nband,npw,nspinor,&  
  &  ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,resid_k,usepaw,vlocal)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimffnl
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: nvloc
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: usepaw
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(in) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(out) :: eig_k(nband)
  real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw(npw)
  real(dp),intent(inout) :: ph3d(2,npw,matblk)
  real(dp),intent(out) :: resid_k(nband)
  real(dp),intent(inout) :: vlocal(n4,n5,n6,nvloc)
 end subroutine mkresi
end interface

interface
 subroutine mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&  
  &  rhog,rhor,rprimd,tim_mkrho,ucvol,unkg,wffnow,wfs,wvl,&  
  &  option) !optional
  use defs_basis
  use defs_wvltypes
  use defs_abitypes
  use m_wffile
  use m_paw_dmft
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in),optional :: option
  integer,intent(in) :: tim_mkrho
  integer,intent(in) :: unkg
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(in) :: paw_dmft
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_wf_type),intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2, &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(out) :: rhog(2,dtset%nfft)
  real(dp), intent(out) :: rhor(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
 end subroutine mkrho
end interface

interface
 subroutine mksubham(cg,ghc,ghc_block,gsc,gvnlc,gvnlc_block,iblock,icg,igsc,ikpt,isppol,istwf_k,&  
  &  isubh,isubo,mcg,mgsc,mpi_enreg,nband_k,nbdblock,npw_k,&  
  &  nspinor,subham,subovl,subvnl,use_subovl,use_vnl,wfoptalg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iblock
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: istwf_k
  integer,intent(inout) :: isubh
  integer,intent(inout) :: isubo
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: nbdblock
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: use_subovl
  integer,intent(in) :: use_vnl
  integer,intent(in) :: wfoptalg
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(inout) :: ghc(2,npw_k*nspinor)
  real(dp),intent(in) :: ghc_block(2,npw_k*nspinor,nbdblock)
  real(dp),intent(in) :: gsc(2,mgsc)
  real(dp),intent(inout) :: gvnlc(2,npw_k*nspinor)
  real(dp),intent(in) :: gvnlc_block(2,npw_k*nspinor,nbdblock*use_vnl)
  real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
  real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp),intent(inout) :: subvnl(nband_k*(nband_k+1)*use_vnl)
 end subroutine mksubham
end interface

interface
 subroutine mlwfovlp(atindx1,cg,cprj,dtset,dtfil,eigen,gprimd,hdr,kg,&  
  &  mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpw,natom,&  
  &  nattyp,nfft,ngfft,nkpt,npwarr,nsppol,ntypat,&  
  &  pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfftc
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(in) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(cprj_type) :: cprj(natom,mcprj)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer :: kg(3,mpw*mkmem)
  integer :: nattyp(ntypat)
  integer :: npwarr(nkpt)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mlwfovlp
end interface

interface
 subroutine mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,gprimd,just_augmentation,kg,&  
  &  lproj,max_num_bands,mband,mkmem,mpi_enreg,mpw,mwan,natom,nattyp,&  
  &  nkpt,npwarr,nspinor,&  
  &  nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,proj_radial,&  
  &  proj_site,proj_x,proj_z,proj_zona,psps,spin,ucvol)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: lproj
  integer,intent(in) :: max_num_bands
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: mwan
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: spin
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp) :: ucvol
  complex(dpc),intent(out) :: A_matrix(max_num_bands,mwan,nkpt,nsppol)
  logical,intent(in) :: band_in(mband,nsppol)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  logical,intent(in) :: just_augmentation(mwan,nsppol)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer :: nattyp(ntypat)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: num_bands(nsppol)
  integer,intent(in) :: nwan(nsppol)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  integer,intent(in) :: proj_l(mband,nsppol)
  integer,intent(in) :: proj_m(mband,nsppol)
  integer,intent(inout) :: proj_radial(mband,nsppol)
  real(dp),intent(in) :: proj_site(3,mband,nsppol)
  real(dp),intent(in) :: proj_x(3,mband,nsppol)
  real(dp),intent(in) :: proj_z(3,mband,nsppol)
  real(dp),intent(in) :: proj_zona(mband,nsppol)
 end subroutine mlwfovlp_proj
end interface

interface
 subroutine mlwfovlp_projpaw(A_paw,band_in,cprj,just_augmentation,max_num_bands,mband,mkmem,&  
  &  mwan,natom,nband,nkpt,&  
  &  nspinor,nsppol,ntypat,nwan,pawrad,pawtab,&  
  &  proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,&  
  &  rprimd,spin,typat,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: max_num_bands
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mwan
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: spin
  type(pseudopotential_type),intent(in) :: psps
  complex(dpc),intent(out) :: A_paw(max_num_bands,mwan,nkpt,nsppol)
  logical,intent(in) :: band_in(mband,nsppol)
  type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  logical,intent(in) :: just_augmentation(mwan,nsppol)
  integer,intent(in) :: nband(nsppol*nkpt)
  integer,intent(in) :: nwan(nsppol)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: proj_l(mband,nsppol)
  integer,intent(in) :: proj_m(mband,nsppol)
  integer,intent(in) :: proj_radial(mband,nsppol)
  real(dp),intent(in) :: proj_site(3,mband,nsppol)
  real(dp),intent(in) :: proj_x(3,mband,nsppol)
  real(dp),intent(in) :: proj_z(3,mband,nsppol)
  real(dp),intent(in) :: proj_zona(mband,nsppol)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mlwfovlp_projpaw
end interface

interface
 subroutine mlwfovlp_pw(cg,cm1,g1,iwav,kg,mband,&  
  &  mkmem,mpi_enreg,mpw,nfft,ngfft,nkpt,nntot,&  
  &  npwarr,nspinor,nsppol,ovikp,seed_name,spin)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nntot
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: spin
  type(mpi_type),intent(in) :: mpi_enreg
  character(len=fnlen) :: seed_name
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(out) :: cm1(2,mband,mband,nntot,nkpt,nsppol)
  integer,intent(in) :: g1(3,nkpt,nntot)
  integer,intent(in) :: iwav(mband,nkpt,nsppol)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  integer,intent(in) :: ovikp(nkpt,nntot)
 end subroutine mlwfovlp_pw
end interface

interface
 subroutine mlwfovlp_radial(alpha,lmax,lmax2,radial,rvalue,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: lmax2
  integer,intent(in) :: rvalue
  real(dp),intent(in) :: alpha
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: radial(lmax2)
 end subroutine mlwfovlp_radial
end interface

interface
 subroutine mlwfovlp_seedname(fname_w90,filew90_win,filew90_wout,filew90_amn,&  
  &  filew90_ramn,filew90_mmn,filew90_eig,nsppol,seed_name)
  use defs_basis
  implicit none
  integer,intent(in) :: nsppol
  character(len=fnlen),intent(in) :: fname_w90
  character(len=fnlen),intent(out) :: filew90_amn(nsppol)
  character(len=fnlen),intent(out) :: filew90_eig(nsppol)
  character(len=fnlen),intent(out) :: filew90_mmn(nsppol)
  character(len=fnlen),intent(out) :: filew90_ramn(nsppol)
  character(len=fnlen),intent(out) :: filew90_win(nsppol)
  character(len=fnlen),intent(out) :: filew90_wout(nsppol)
  character(len=fnlen),intent(out) :: seed_name(nsppol)
 end subroutine mlwfovlp_seedname
end interface

interface
 subroutine mlwfovlp_setup(atom_symbols,band_in,dtset,filew90_win,gamma_only,&  
  &  g1,lwanniersetup,mband,natom,nband_inc,nkpt,&  
  &  nntot,num_bands,num_nnmax,nsppol,nwan,ovikp,&  
  &  proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,&  
  &  real_lattice,recip_lattice,rprimd,seed_name,spin,spinors,xcart,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: lwanniersetup
  integer,intent(in) :: mband
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(out) :: nntot
  integer,intent(in) :: nsppol
  integer,intent(in) :: num_nnmax
  integer,intent(in) :: spin
  type(dataset_type),intent(in) :: dtset
  logical,intent(in) :: gamma_only
  logical,intent(in) :: spinors
  character(len=3),intent(out) :: atom_symbols(natom)
  logical,intent(out) :: band_in(mband,nsppol)
  character(len=fnlen),intent(in) :: filew90_win(nsppol)
  integer,intent(out) :: g1(3,nkpt,num_nnmax)
  integer,intent(out) :: nband_inc(nsppol)
  integer,intent(out) :: num_bands(nsppol)
  integer,intent(out) :: nwan(nsppol)
  integer,intent(out) :: ovikp(nkpt,num_nnmax)
  integer,intent(out) :: proj_l(mband,nsppol)
  integer,intent(out) :: proj_m(mband,nsppol)
  integer,intent(out) :: proj_radial(mband,nsppol)
  real(dp),intent(out) :: proj_site(3,mband,nsppol)
  real(dp),intent(out) :: proj_x(3,mband,nsppol)
  real(dp),intent(out) :: proj_z(3,mband,nsppol)
  real(dp),intent(out) :: proj_zona(mband,nsppol)
  real(dp),intent(in) :: real_lattice(3,3)
  real(dp),intent(in) :: recip_lattice(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  character(len=fnlen),intent(in) :: seed_name(nsppol)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mlwfovlp_setup
end interface

interface
 subroutine mlwfovlp_ylmfac(ylmc_fac,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)
  use defs_basis
  implicit none
  integer, intent(in) :: lmax
  integer, intent(in) :: lmax2
  integer, intent(in) :: mband
  integer, intent(in) :: nwan
  integer,intent(in) :: proj_l(mband)
  integer,intent(in) :: proj_m(mband)
  real(dp),intent(in) :: proj_x(3,mband)
  real(dp),intent(in) :: proj_z(3,mband)
  complex(dp),intent(out) :: ylmc_fac(lmax2,nwan)
 end subroutine mlwfovlp_ylmfac
end interface

interface
 subroutine mlwfovlp_ylmfar(ylmr_fac,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)
  use defs_basis
  implicit none
  integer, intent(in) :: lmax
  integer, intent(in) :: lmax2
  integer, intent(in) :: mband
  integer, intent(in) :: nwan
  integer,intent(in) :: proj_l(mband)
  integer,intent(in) :: proj_m(mband)
  real(dp),intent(in) :: proj_x(3,mband)
  real(dp),intent(in) :: proj_z(3,mband)
  real(dp),intent(out) :: ylmr_fac(lmax2,nwan)
 end subroutine mlwfovlp_ylmfar
end interface

interface
 subroutine moddiel(cplex,dielar,mpi_enreg,nfft,ngfft,nspden,optreal,optres,paral_kgb,qphon,rprimd,vresid,vrespc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: optreal
  integer,intent(in) :: optres
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vresid(cplex*nfft,nspden)
  real(dp),intent(out) :: vrespc(cplex*nfft,nspden)
 end subroutine moddiel
end interface

interface
 subroutine msig(fcti,npti,xi,filnam_out_sig)
  use defs_basis
  implicit none
  integer,intent(in) :: npti
  character(len=fnlen),intent(in) :: filnam_out_sig
  real(dp),intent(in) :: fcti(npti)
  real(dp),intent(in) :: xi(npti)
 end subroutine msig
end interface

interface
 subroutine newkpt(ceksp2,cg,debug,ecut1,ecut2,ecut2_eff,eigen,exchn2n3d,fill,&  
  &  formeig,gmet1,gmet2,headform1,indkk,iout,ireadwf,istwfk1,istwfk2,&  
  &  kg2,kptns1,kptns2,mband2,mcg,mkmem1,mkmem2,&  
  &  mpi_enreg1,mpi_enreg2,mpw1,mpw2,&  
  &  nband1,nband2,ngfft,nkpt1,nkpt2,npwarr1,npwarr2,nspinor1,nspinor2,&  
  &  nsppol1,nsppol2,nsym,occ,optorth,prtvol,restart,rprimd,sppoldbl,&  
  &  symrel,tnons,unkg2,wffinp,wffout)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: ceksp2
  integer,intent(in) :: debug
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: fill
  integer,intent(in) :: formeig
  integer,intent(in) :: headform1
  integer,intent(in) :: iout
  integer,intent(in) :: ireadwf
  integer,intent(in) :: mband2
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem1
  integer,intent(in) :: mkmem2
  integer,intent(in) :: mpw1
  integer,intent(in) :: mpw2
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nspinor1
  integer,intent(in) :: nspinor2
  integer,intent(in) :: nsppol1
  integer,intent(in) :: nsppol2
  integer,intent(in) :: nsym
  integer,intent(in) :: optorth
  integer,intent(in) :: prtvol
  integer,intent(in) :: restart
  integer,intent(in) :: sppoldbl
  integer,intent(in) :: unkg2
  real(dp),intent(in) :: ecut1
  real(dp),intent(in) :: ecut2
  real(dp),intent(in) :: ecut2_eff
  type(mpi_type),intent(inout) :: mpi_enreg1
  type(mpi_type),intent(inout) :: mpi_enreg2
  type(wffile_type),intent(inout) :: wffinp
  type(wffile_type),intent(inout) :: wffout
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: cg(2,mcg)
  real(dp),intent(out) :: eigen(mband2*(2*mband2)**formeig*nkpt2*nsppol2)
  real(dp),intent(in) :: gmet1(3,3)
  real(dp),intent(in) :: gmet2(3,3)
  integer,intent(in) :: indkk(nkpt2*sppoldbl,6)
  integer,intent(in) :: istwfk1(nkpt1)
  integer,intent(in) :: istwfk2(nkpt2)
  integer,intent(in) :: kg2(3,mpw2*mkmem2)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  integer,intent(in) :: nband1(nkpt1*nsppol1)
  integer,intent(in) :: nband2(nkpt2*nsppol2)
  integer,intent(in) :: npwarr1(nkpt1)
  integer,intent(in) :: npwarr2(nkpt2)
  real(dp),intent(out) :: occ(mband2*nkpt2*nsppol2)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine newkpt
end interface

interface
 subroutine nlinopt(nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,evalv,efermi,&  
  pmat,v1,v2,v3,nmesh,de,sc,brod,tol,fnam)
  use defs_basis
  implicit none
  integer, intent(in) :: nkpt
  integer, intent(in) :: nmesh
  integer, intent(in) :: nspin
  integer, intent(in) :: nstval
  integer, intent(in) :: nsymcrys
  integer, intent(in) :: v1
  integer, intent(in) :: v2
  integer, intent(in) :: v3
  real(dp), intent(in) :: brod
  real(dp), intent(in) :: de
  real(dp), intent(in) :: efermi
  character(256), intent(in) :: fnam
  real(dp), intent(in) :: omega
  real(dp), intent(in) :: sc
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: evalv(nstval,nspin,nkpt)
  complex(dpc), intent(inout) :: pmat(nstval,nstval,nkpt,3,nspin)
  real(dp), intent(in) :: symcrys(3,3,nsymcrys)
  real(dp), intent(in) :: wkpt(nkpt)
 end subroutine nlinopt
end interface

interface
 subroutine nres2vres(dtset,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,nhat,&  
  &  nkxc,nresid,n3xccc,optnc,optxc,pawang,pawfgrtab,pawrhoij,pawtab,&  
  &  rhor,rprimd,usepaw,usexcnhat,vresid,wvl,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: izero
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: optnc
  integer,intent(in) :: optxc
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(in) :: nresid(nfft,dtset%nspden)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%ntypat*usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vresid(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine nres2vres
end interface

interface
 subroutine odamix(deltae,dtset,efield_dot,elast,energies,etotal,&  
  &  gprimd,gsqcut,kxc,mag_cart,mpi_enreg,nfft,ngfft,nhat,&  
  &  nkxc,ntypat,nvresid,n3xccc,optres,paw_ij,&  
  &  paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,pel,&  
  &  pion,psps,rhog,rhor,rprimd,strsxc,taug,taur,ucvol,usepaw,&  
  &  usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred,&  
  &  vxctau) ! optional argument
  use defs_basis
  use m_energies
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optres
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(out) :: deltae
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: elast
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vxcavg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: efield_dot(3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: mag_cart(3)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(inout) :: taug(2,nfft*dtset%usekden)
  real(dp),intent(inout) :: taur(nfft,dtset%nspden*dtset%usekden)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout),optional :: vxctau(nfft,dtset%nspden*dtset%usekden,4)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine odamix
end interface

interface
 subroutine optics_vloc(cg,dtfil,dtset,eigen0,gprimd,hdr,kg,mband,mcg,mkmem,mpi_enreg,mpw,&  
  &  nkpt,npwarr,nsppol,wffnow)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wffnow
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
 end subroutine optics_vloc
end interface

interface
 subroutine outwf(cg,dtset,eigen,filnam,hdr,kg,kptns,mband,mcg,mkmem,&  
  &  mpi_enreg,mpw,mxfh,natom,nband,nkpt,npwarr,&  
  &  nsppol,nstep,nxfh,occ,resid,response,unwff2,&  
  &  wffnow,wfs,wvl,xfhist)
  use defs_basis
  use defs_abitypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer, intent(in) :: mband
  integer, intent(in) :: mcg
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: mxfh
  integer, intent(in) :: natom
  integer, intent(in) :: nkpt
  integer, intent(in) :: nsppol
  integer, intent(in) :: nstep
  integer, intent(in) :: nxfh
  integer, intent(in) :: response
  integer, intent(in) :: unwff2
  type(dataset_type), intent(in) :: dtset
  character(len=fnlen), intent(in) :: filnam
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  type(wffile_type), intent(inout) :: wffnow
  type(wvl_wf_type), intent(in) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(inout) :: cg(2,mcg)
  real(dp), intent(in) :: eigen((2*mband)**response*mband*nkpt*nsppol)
  integer, intent(in) :: kg(3,mpw*mkmem)
  real(dp), intent(in) :: kptns(3,nkpt)
  integer, intent(in) :: nband(nkpt*nsppol)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ(mband*nkpt*nsppol)
  real(dp), intent(in) :: resid(mband*nkpt*nsppol)
  real(dp), intent(in) :: xfhist(3,natom+4,2,mxfh)
 end subroutine outwf
end interface

interface
 subroutine partial_dos_fractions(cg,dos_fractions,dos_fractions_m,dtfil,&  
  &  dtset,hdr,mbesslang,mcg,mpi_enreg,m_dos_flag,ndosfraction,partial_dos,wffnow)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: mcg
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: partial_dos
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wffnow
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(out) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(out) :: dos_fractions_m(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)
 end subroutine partial_dos_fractions
end interface

interface
 subroutine poslifetime(dtset,electronpositron,gprimd,mpi_enreg,n3xccc,nfft,ngfft,nzlmopt,&  
  &  paw_an,pawang,pawrad,pawrhoij,pawtab,rhor,ucvol,xccc3d)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nzlmopt
  type(dataset_type), intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type), intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  type(paw_an_type),intent(in) :: paw_an(dtset%natom*dtset%usepaw)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*dtset%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
 end subroutine poslifetime
end interface

interface
 subroutine prteigrs(eigen,enunit,fermie,fname_eig,iout,iscf,kptns,kptopt,mband,nband,&  
  &  nkpt,nnsclo_now,nsppol,occ,occopt,option,prteig,prtvol,resid,tolwfr,vxcavg,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: enunit
  integer,intent(in) :: iout
  integer,intent(in) :: iscf
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nnsclo_now
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: prteig
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: fname_eig
  real(dp),intent(in) :: tolwfr
  real(dp),intent(in) :: vxcavg
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: resid(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine prteigrs
end interface

interface
 subroutine prtene(dtset,energies,iout,usepaw)
  use m_energies
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  type(energies_type),intent(in) :: energies
 end subroutine prtene
end interface

interface
 subroutine prtimg(dynimage,imagealgo_str,imgmov,iout,mpi_enreg,nimage,nimage_tot,&  
  &  prt_all_images,prtvolimg,resimg)
  use m_results_img
  use defs_abitypes
  implicit none
  integer,intent(in) :: imgmov
  integer,intent(in) :: iout
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: prtvolimg
  character(len=60),intent(in) :: imagealgo_str
  type(mpi_type),intent(inout) :: mpi_enreg
  logical,intent(in) :: prt_all_images
  integer,intent(in) :: dynimage(nimage_tot)
  type(results_img_type),target,intent(inout) :: resimg(nimage)
 end subroutine prtimg
end interface

interface
 subroutine prtrhomxmn(iout,mpi_enreg,nfft,ngfft,nspden,option,rhor,optrhor,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in),optional :: optrhor
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in),optional :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(nfft,nspden)
 end subroutine prtrhomxmn
end interface

interface
 subroutine prtxf(fred,iatfix,iout,iwfrc,natom,rprimd,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: iwfrc
  integer,intent(in) :: natom
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prtxf
end interface

interface
 subroutine prtxvf(fcart,fred,iatfix,iout,natom,prtvel,vel,xcart,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: prtvel
  real(dp),intent(in) :: fcart(3,natom)
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: vel(3,natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prtxvf
end interface

interface
 subroutine rhophi(cx,phi,rho)
  use defs_basis
  implicit none
  real(dp),intent(out) :: phi
  real(dp),intent(out) :: rho
  real(dp),intent(in) :: cx(2)
 end subroutine rhophi
end interface

interface
 subroutine rhotov(dtset,energies,gprimd,gsqcut,kxc,mpi_enreg,nfft,ngfft,&  
  &  nhat,nhatgr,nhatgrdim,nkxc,vresidnew,n3xccc,optene,optres,optxc,&  
  &  rhog,rhor,rprimd,strsxc,ucvol,usepaw,usexcnhat,&  
  &  vhartr,vnew_mean,vpsp,vres_mean,vres2,vtrial,vxcavg,vxc,wvl,xccc3d,&  
  &  electronpositron,taug,taur,vxctau) ! optional argument
  use defs_basis
  use m_energies
  use defs_abitypes
  use m_electronpositron
  use defs_wvltypes
  implicit none
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vres2
  real(dp),intent(out) :: vxcavg
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*usepaw)
  real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
  real(dp),intent(in),optional :: taur(nfft,dtset%nspden*dtset%usekden)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(out) :: vnew_mean(dtset%nspden)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(out) :: vres_mean(dtset%nspden)
  real(dp),intent(out) :: vresidnew(nfft,dtset%nspden)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(out),optional :: vxctau(nfft,dtset%nspden*dtset%usekden,4)
  real(dp),intent(inout) :: xccc3d(n3xccc)
 end subroutine rhotov
end interface

interface
 subroutine scprqt(choice,cpus,deltae,diffor,dtset,&  
  &  eigen,etotal,favg,fcart,fermie,fname_eig,filnam1,initGS,&  
  &  iscf,istep,kpt,maxfor,moved_atm_inside,mpi_enreg,&  
  &  nband,nkpt,nstep,occ,optres,&  
  &  prtfor,prtxml,quit,res2,resid,residm,response,tollist,usepaw,&  
  &  vxcavg,wtk,xred,&  
  &  electronpositron) ! optional argument)
  use defs_basis
  use defs_abitypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: initGS
  integer,intent(in) :: iscf
  integer,intent(in) :: istep
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: nkpt
  integer,intent(in) :: nstep
  integer,intent(in) :: optres
  integer,intent(in) :: prtfor
  integer,intent(in) :: prtxml
  integer,intent(out) :: quit
  integer,intent(in) :: response
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: cpus
  real(dp),intent(in) :: deltae
  real(dp),intent(in) :: diffor
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: filnam1
  character(len=fnlen),intent(in) :: fname_eig
  real(dp),intent(in) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: res2
  real(dp),intent(in) :: residm
  real(dp),intent(in) :: vxcavg
  real(dp),intent(in) :: eigen(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: favg(3)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt*dtset%nsppol)
  real(dp),intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: resid(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: tollist(12)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine scprqt
end interface

interface
 subroutine set_twind(dtefield)
  use m_efield
  implicit none
  type(efield_type), intent(inout) :: dtefield
 end subroutine set_twind
end interface

interface
 subroutine setup1(acell,amass,bantot,dtset,ecut_eff,ecutc_eff,gmet,&  
  &  gprimd,gsqcut_eff,gsqcutc_eff,natom,ngfft,ngfftc,nkpt,nsppol,&  
  &  response,rmet,rprim,rprimd,ucvol,usepaw)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(out) :: bantot
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: response
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecut_eff
  real(dp),intent(in) :: ecutc_eff
  real(dp),intent(out) :: gsqcut_eff
  real(dp),intent(out) :: gsqcutc_eff
  real(dp),intent(out) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftc(18)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(out) :: amass(natom)
  real(dp),intent(out) :: gmet(3,3)
  real(dp),intent(out) :: gprimd(3,3)
  real(dp),intent(out) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rprimd(3,3)
 end subroutine setup1
end interface

interface
 subroutine setup2(dtset,&  
  &  npwtot,start,wfs,xred)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(wvl_wf_type),intent(in) :: wfs
  integer,intent(in) :: npwtot(dtset%nkpt)
  real(dp),intent(out) :: start(3,dtset%natom)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine setup2
end interface

interface
 subroutine setup_positron(atindx,atindx1,cg,dtefield,dtfil,dtset,ecore,eigen,etotal,electronpositron,&  
  &  energies,forces_needed,fred,gmet,gprimd,grewtn,gsqcut,hdr,ifirst_gs,indsym,istep,istep_mix,kg,&  
  &  kxc,maxfor,mcg,mgfft,mpi_enreg,n3xccc,nattyp,nfft,ngfft,nhat,nkxc,npwarr,nvresid,occ,optres,&  
  &  paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,pel,ph1d,ph1dc,pion,psps,rhog,rhor,&  
  &  rprimd,stress_needed,strsxc,symrec,ucvol,usexcnhat,vhartr,vpsp,vxc,wffnow,&  
  &  xccc3d,xred,ylm,ylmgr)
  use m_energies
  use defs_abitypes
  use m_wffile
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: forces_needed
  integer,intent(in) :: ifirst_gs
  integer,intent(in) :: istep
  integer,intent(inout) :: istep_mix
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: optres
  integer,intent(in) :: stress_needed
  integer,intent(in) :: usexcnhat
  type(efield_type),intent(in) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(hdr_type),intent(inout) :: hdr
  real(dp),intent(in) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type), intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: fred(3,dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(dtset%natom)
  real(dp),intent(inout) :: nhat(nfft,dtset%nspden*dtset%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*dtset%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*dtset%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*dtset%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: ph1dc(2,(3*(2*dtset%mgfft+1)*dtset%natom)*dtset%usepaw)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strsxc(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vpsp(nfft)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine setup_positron
end interface

interface
 subroutine setvtr(atindx1,dtset,energies,gmet,gprimd,grewtn,gsqcut,&  
  &  istep,kxc,mgfft,moved_atm_inside,moved_rhor,mpi_enreg,&  
  &  nattyp,nfft,ngfft,nhat,nhatgr,nhatgrdim,nkxc,ntypat,n1xccc,n3xccc,&  
  &  optene,pawtab,ph1d,psps,rhog,rhor,rmet,rprimd,strsxc,&  
  &  ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,xccc3d,xred,&  
  &  electronpositron,taug,taur,vxctau) ! optional argument
  use m_energies
  use defs_abitypes
  use defs_wvltypes
  use defs_basis
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(inout) :: moved_atm_inside
  integer,intent(inout) :: moved_rhor
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: usexcnhat
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: vxcavg
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfft,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim)
  type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: strsxc(6)
  real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
  real(dp),intent(inout),optional :: taur(nfft,dtset%nspden*dtset%usekden)
  real(dp),intent(inout) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfft,dtset%nspden)
  real(dp),intent(out),optional :: vxctau(nfft,dtset%nspden*dtset%usekden,4)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine setvtr
end interface

interface
 subroutine smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&  
  &  mcg_k,mcg_q,mcg1_k,minbd,mpw,nband_occ,npw_k1,npw_k2,nspinor,&  
  &  pwind_k,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_k,smat_k_paw,usepaw)
  use defs_basis
  implicit none
  integer,intent(in) :: ddkflag
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: itrs
  integer,intent(in) :: job
  integer,intent(in) :: maxbd
  integer,intent(in) :: mcg1_k
  integer,intent(in) :: mcg_k
  integer,intent(in) :: mcg_q
  integer,intent(in) :: minbd
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_occ
  integer,intent(in) :: npw_k1
  integer,intent(in) :: npw_k2
  integer,intent(in) :: nspinor
  integer,intent(in) :: shiftbd
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: cg(2,mcg_k)
  real(dp),intent(out) :: cg1_k(2,mcg1_k)
  real(dp),intent(in) :: cgq(2,mcg_q)
  real(dp),intent(out) :: dtm_k(2)
  integer,intent(in) :: pwind_k(mpw)
  real(dp),intent(in) :: pwnsfac_k(4,mpw)
  integer,intent(inout) :: sflag_k(nband_occ)
  real(dp),intent(out) :: smat_inv(2,nband_occ,nband_occ)
  real(dp),intent(inout) :: smat_k(2,nband_occ,nband_occ)
  real(dp),intent(in) :: smat_k_paw(2,usepaw*nband_occ,usepaw*nband_occ)
 end subroutine smatrix
end interface

interface
 subroutine spin_current(cg,dtfil,dtset,gprimd,hdr,kg,mcg,mpi_enreg,psps)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mcg
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 end subroutine spin_current
end interface

interface
 subroutine stress(atindx1,berryopt,eei,efield,ehart,eii,gsqcut,kinstr,&  
  &  mgfft,mpi_enreg,mqgrid,n1xccc,n3xccc,natom,nattyp,&  
  &  nfft,ngfft,nlstr,nspden,nsym,ntypat,paral_kgb,pawtab,pel,pion,ph1d,&  
  &  prtvol,qgrid,rhog,rprimd,strten,strsxc,symrec,typat,usepaw,vlspl,&  
  &  vxc,xccc1d,xccc3d,xcccrc,xred,zion,&  
  &  electronpositron) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: berryopt
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
  integer,intent(in) :: prtvol
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: ehart
  real(dp),intent(in) :: eii
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: efield(3)
  real(dp),intent(in) :: kinstr(6)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nlstr(6)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strsxc(6)
  real(dp),intent(out) :: strten(6)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc*(1-usepaw),6,ntypat)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine stress
end interface

interface
 subroutine strhar(ehart,gprimd,gsqcut,harstr,mpi_enreg,nfft,ngfft,rhog,ucvol,&  
  &  rhog2) ! optional argument
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: ehart
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: harstr(6)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(in),optional :: rhog2(2,nfft)
 end subroutine strhar
end interface

interface
 subroutine sygrad(fred,natom,dedt,nsym,symrec,indsym)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp),intent(in) :: dedt(3,natom)
  real(dp),intent(out) :: fred(3,natom)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine sygrad
end interface

interface
 subroutine symrhg(cplex,gprimd,irrzon,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym,paral_kgb,&  
  &  phnons,rhog,rhor,rprimd,symafm,symrel)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzon(nfftot**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: phnons(2,nfftot**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symrhg
end interface

interface
 subroutine testsusmat(compute,dielop,dielstrt,dtset,istep)
  use defs_abitypes
  implicit none
  integer,intent(in) :: dielop
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  logical,intent(out) :: compute
  type(dataset_type),intent(in) :: dtset
 end subroutine testsusmat
end interface

interface
 subroutine uderiv(bdberry,cg,gprimd,hdr,istwfk,kberry,kg,kpt_,kptopt,kptrlatt,&  
  &  mband,mcg,mkmem,mpi_enreg,mpw,natom,nband,nberry,npwarr,nspinor,nsppol,nkpt_,&  
  &  unddk,unkg,wffnow,fnameabo_1wf)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nberry
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: unddk
  integer,intent(in) :: unkg
  character(len=fnlen),intent(in) :: fnameabo_1wf
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: bdberry(4)
  integer,intent(in) :: kberry(3,20)
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: gprimd(1:3,1:3)
  integer,intent(in) :: istwfk(nkpt_)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt_(3,nkpt_)
  integer,intent(in) :: nband(nkpt_*nsppol)
  integer,intent(in) :: npwarr(nkpt_)
 end subroutine uderiv
end interface

interface
 subroutine update_mmat(berryopt,cg,cgq,dimffnl,dtefield,ffnl,filstat,&  
  &  gmet,gprimd,gs_hamk,icg,ikpt,kg,kinpw,lmnmax,matblk,&  
  &  mband,mcg,mcgq,mgfft,mkgq,mkmem,mpi_enreg,mpsang,mpssoang,mpw,&  
  &  natom,nkpg,nkpt,npw_k,npwarr,nspinor,ntypat,&  
  &  nvloc,n4,n5,n6,pawtab,psps,pwind,pwind_alloc,&  
  &  paral_kgb,ph3d,prtvol,pwnsfac,pwnsfacq,rmet,ucvol,vlocal,ylm,ylmgr)
  use defs_basis
  use defs_abitypes
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: dimffnl
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcgq
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkgq
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: nkpg
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: nvloc
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: pwind_alloc
  type(efield_type),intent(inout) :: dtefield
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type), intent(in) :: gs_hamk
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: cgq(2,mcgq)
  real(dp),intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer, intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kinpw(npw_k)
  integer, intent(in) :: npwarr(nkpt)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(inout) :: ph3d(2,npw_k,matblk)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: pwnsfacq(2,mkgq)
  real(dp),intent(in) :: rmet(3,3)
  real(dp), intent(inout) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(in) :: ylm(mpw*mkmem,mpsang*mpsang)
  real(dp), intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang)
 end subroutine update_mmat
end interface

interface
 subroutine vso_realspace_local(dtset,hdr,position_op,psps,vso_realspace)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: position_op(3,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3))
  real(dp),intent(out) :: vso_realspace(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3), &
  &         dtset%nspinor,dtset%nspinor,3)
 end subroutine vso_realspace_local
end interface

interface
 subroutine vso_realspace_nonlop(atindx,atindx1,dtfil,dtset,gmet,gprimd,hdr,kg,&  
  &  mpi_enreg,nattyp,ph1d,position_op,psps,rmet,ucvol,vso_realspace_nl,ylm,ylmgr)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(dtset%ntypat)
  real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: position_op(3,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3))
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(out) :: vso_realspace_nl(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3), &
  &         dtset%nspinor,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine vso_realspace_nonlop
end interface

interface
 subroutine vtorhotf(dtfil,dtset,ek,enl,entropy,fermie,gprimd,grnl,&  
  &  irrzon,mpi_enreg,natom,nfft,nspden,nsppol,nsym,phnons,rhog,rhor,rprimd,ucvol,vtrial)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(out) :: entropy
  real(dp),intent(out) :: fermie
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grnl(3*natom)
  integer,intent(in) :: irrzon((dtset%ngfft(1)*dtset%ngfft(1)*dtset%ngfft(1))**(1-1/nsym), &
  &         2,(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(1)*dtset%ngfft(1))**(1-1/nsym), &
  &         (nspden/nsppol)-3*(nspden/4))
  real(dp),intent(inout) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vtrial(nfft,nspden)
 end subroutine vtorhotf
end interface

interface
 function zfermim12(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermim12
 end function zfermim12
end interface

interface
 function zfermi12(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi12
 end function zfermi12
end interface

interface
 function zfermi1(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi1
 end function zfermi1
end interface

interface
 function zfermi32(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi32
 end function zfermi32
end interface

interface
 function zfermi2(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi2
 end function zfermi2
end interface

interface
 function zfermi52(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi52
 end function zfermi52
end interface

interface
 function zfermi3(xx)
  use defs_basis
  implicit none
  real(dp), intent(in) :: xx
  real(dp) :: zfermi3
 end function zfermi3
end interface

interface
 function ifermim12(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermim12
 end function ifermim12
end interface

interface
 function ifermi12(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi12
 end function ifermi12
end interface

interface
 function ifermi32(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi32
 end function ifermi32
end interface

interface
 function ifermi52(ff)
  use defs_basis
  implicit none
  real(dp), intent(in) :: ff
  real(dp) :: ifermi52
 end function ifermi52
end interface

interface
 function fp12a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fp12a1
  real(dp) :: x
 end function fp12a1
end interface

interface
 function fp32a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fp32a1
  real(dp) :: x
 end function fp32a1
end interface

interface
 function xp12a1 (y)
  use defs_basis
  implicit none
  real(dp) :: xp12a1
  real(dp) :: y
 end function xp12a1
end interface

interface
 function fm12a1 (x)
  use defs_basis
  implicit none
  real(dp) :: fm12a1
  real(dp) :: x
 end function fm12a1
end interface

interface
 subroutine fm12a1t (cktf,rtnewt,tsmear,vtrial,rhor_middx,rhor_mid,&  
  &  nfft)
  use defs_basis
  implicit none
  integer :: nfft
  real(dp) :: cktf
  real(dp) :: rtnewt
  real(dp) :: tsmear
  real(dp) :: rhor_mid(nfft)
  real(dp) :: rhor_middx(nfft)
  real(dp) :: vtrial(nfft)
 end subroutine fm12a1t
end interface

interface
 subroutine waveformat(cg,cg_disk,cg_index,cg_new,dk,ii,ikpt,&  
  &  ikpt_,isgn,isppol,jj,jkpt,jkpt_,kg_kpt,kpt,kg_jl,maxband,mband,mcg,mcg_disk,&  
  &  minband,mkmem,mpw,nkpt,nkpt_,npwarr,nsppol,nspinor,shift_g_2,tr)
  use defs_basis
  implicit none
  integer,intent(in) :: ii
  integer,intent(in) :: ikpt
  integer,intent(in) :: ikpt_
  integer,intent(in) :: isgn
  integer,intent(in) :: isppol
  integer,intent(in) :: jj
  integer,intent(in) :: jkpt
  integer,intent(in) :: jkpt_
  integer,intent(in) :: maxband
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcg_disk
  integer,intent(in) :: minband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: cg_disk(2,mcg_disk,2)
  integer,intent(in) :: cg_index(mband,nkpt_,nsppol)
  real(dp),intent(out) :: cg_new(2,mpw,maxband)
  real(dp),intent(in) :: dk(3)
  integer,intent(in) :: kg_jl(3,mpw,2)
  integer,intent(in) :: kg_kpt(3,mpw*nspinor,nkpt_)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: npwarr(nkpt_)
  logical,intent(in) :: shift_g_2(nkpt,nkpt)
  real(dp),intent(in) :: tr(2)
 end subroutine waveformat
end interface

interface
 subroutine wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wfs, wvl)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  type(wvl_wf_type),intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2, &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp),intent(inout) :: rhor(dtset%nfft,dtset%nspden)
 end subroutine wvl_mkrho
end interface

interface
 subroutine wvl_newvtr(dtset, mpi_enreg, nele, offset, vhartr, vpsp, vtrial, vxc, wvl)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer,intent(out) :: nele
  integer,intent(out) :: offset
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  real(dp),intent(in) :: vhartr(dtset%nfft)
  real(dp),intent(in) :: vpsp(dtset%nfft)
  real(dp),intent(out) :: vtrial(dtset%nfft*dtset%nspden)
  real(dp),intent(in) :: vxc(dtset%nfft*dtset%nspden)
 end subroutine wvl_newvtr
end interface

end module interfaces_67_common
!!***
