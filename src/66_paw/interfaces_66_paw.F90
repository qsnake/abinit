!!****m* ABINIT/interfaces_66_paw
!! NAME
!! interfaces_66_paw
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/66_paw
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

module interfaces_66_paw

 implicit none

interface
 subroutine Lij(dtefield,ntypat,pawrad,pawtab,psps)
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  type(pseudopotential_type),intent(in) :: psps
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
 end subroutine Lij
end interface

interface
 subroutine atomden(MPI_enreg,natom,ntypat,typat,ngrid,r_vec_grid,rho,a,b,c,atom_pos,&  
  &  natomgr,natomgrmax,atomrgrid,density,prtvol,calctype)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natomgrmax
  integer,intent(in) :: ngrid
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=7),intent(in) :: calctype
  real(dp),intent(in) :: a(3)
  real(dp),intent(in) :: atom_pos(3,natom)
  real(dp),intent(in) :: atomrgrid(natomgrmax,ntypat)
  real(dp),intent(in) :: b(3)
  real(dp),intent(in) :: c(3)
  real(dp),intent(in) :: density(natomgrmax,ntypat)
  integer,intent(in) :: natomgr(ntypat)
  real(dp),intent(in) :: r_vec_grid(3,ngrid)
  real(dp),intent(inout) :: rho(ngrid)
  integer,intent(in) :: typat(natom)
 end subroutine atomden
end interface

interface
 subroutine chkpawovlp(natom,ntypat,pawovlp,pawtab,rmet,typat,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp) :: pawovlp
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rmet(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine chkpawovlp
end interface

interface
 subroutine denfgr(atindx1,gmet,spaceComm_in,natom,nattyp,ngfft,nhat,nspinor,nsppol,nspden,ntypat,&  
  &  pawfgr,pawrad,pawrhoij,pawtab,prtvol,psps,rhor,rhor_paw,rhor_n_one,rhor_nt_one,&  
  &  rprimd,typat,ucvol,xred,abs_n_tilde_nt_diff,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: spaceComm_in
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),optional,intent(out) :: abs_n_tilde_nt_diff(nspden)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(pawfgr%nfft,nspden)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rhor(pawfgr%nfft,nspden)
  real(dp),intent(out) :: rhor_n_one(pawfgr%nfft,nspden)
  real(dp),intent(out) :: rhor_nt_one(pawfgr%nfft,nspden)
  real(dp),intent(out) :: rhor_paw(pawfgr%nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),optional,intent(in) :: znucl(ntypat)
 end subroutine denfgr
end interface

interface
 subroutine expibi(dtefield,gprimd,natom,rprimd,xred)
  use defs_basis
  use m_efield
  implicit none
  integer,intent(in) :: natom
  type(efield_type),intent(inout) :: dtefield
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine expibi
end interface

interface
 subroutine expibr(dtefield,gprimd,natom,pawfgrtab,xred)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  type(efield_type),intent(inout) :: dtefield
  real(dp),intent(in) :: gprimd(3,3)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine expibr
end interface

interface
 subroutine fourier_interpol(cplex,nspden,optin,optout,nfft_in,ngfft_in,nfft_out,ngfft_out,&  
  &  paral_kgb,MPI_enreg,rhor_in,rhor_out,rhog_in,rhog_out)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft_in
  integer,intent(in) :: nfft_out
  integer,intent(in) :: nspden
  integer,intent(in) :: optin
  integer,intent(in) :: optout
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: MPI_enreg
  integer,intent(in) :: ngfft_in(18)
  integer,intent(in) :: ngfft_out(18)
  real(dp),intent(inout) :: rhog_in(2,nfft_in)
  real(dp),intent(out) :: rhog_out(2,nfft_out)
  real(dp),intent(inout) :: rhor_in(cplex*nfft_in,nspden)
  real(dp),intent(out) :: rhor_out(cplex*nfft_out,nspden)
 end subroutine fourier_interpol
end interface

interface
 subroutine initang(pawang)
  use defs_datatypes
  implicit none
  type(pawang_type),intent(inout) :: pawang
 end subroutine initang
end interface

interface
 subroutine initrhoij(cplex,indlmn,lexexch,lmnmax,lpawu,mpi_enreg,natom,natom_paw,&  
  &  nspden,nspinor,nsppol,ntypat,pawrhoij,pawtab,spinat,typat,&  
  &  ngrhoij,nlmnmix,use_rhoij_,use_rhoijres) ! Optional arguments
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: natom_paw
  integer,intent(in),optional :: ngrhoij
  integer,intent(in),optional :: nlmnmix
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in),optional :: use_rhoij_
  integer,intent(in),optional :: use_rhoijres
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: lexexch(ntypat)
  integer,intent(in) :: lpawu(ntypat)
  type(pawrhoij_type),intent(out) :: pawrhoij(natom_paw)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(in) :: typat(natom)
 end subroutine initrhoij
end interface

interface
 subroutine int_ang(ang_phipphj,mpsang)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(out) :: ang_phipphj(mpsang**2,mpsang**2,8)
 end subroutine int_ang
end interface

interface
 subroutine linear_optics_paw(filnam,filnam_out,mpi_enreg_seq)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen),intent(in) :: filnam
  character(len=fnlen),intent(in) :: filnam_out
  type(mpi_type),intent(inout) :: mpi_enreg_seq
 end subroutine linear_optics_paw
end interface

interface
 subroutine mag_loc_k(atindx,atindx1,cg,cprj,dimffnl,dtefield,ffnl,gmet,gprimd,icg,&  
  &  ikpt,indlmn,istwfk_k,kg_k,kpg_k,kpt,lmnmax,matblk,mcg,mgfft,mpi_enreg,&  
  &  mpsang,mpssoang,natom,nattyp,nband_k,ngfft,nkpg,nloalg,npw_k,&  
  &  nspinor,ntypat,pawtab,phkxred,ph1d,ph3d,ucvol)
  use defs_basis
  use defs_abitypes
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimffnl
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: istwfk_k
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: nkpg
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(cprj_type),intent(in) :: cprj(natom,nband_k)
  real(dp),intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kpg_k(npw_k,nkpg)
  real(dp),intent(in) :: kpt(3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(inout) :: ph3d(2,npw_k,matblk)
  real(dp),intent(in) :: phkxred(2,natom)
 end subroutine mag_loc_k
end interface

interface
 subroutine make_efg_onsite(efg,natom,ntypat,paw_an,pawang,pawrhoij,pawrad,pawtab,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(out) :: efg(3,3,natom)
  type(paw_an_type),intent(in) :: paw_an(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine make_efg_onsite
end interface

interface
 subroutine make_fc_paw(fc,natom,ntypat,pawrhoij,pawrad,pawtab,psps,typat)
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pseudopotential_type),intent(in) :: psps
  type(nuclear_type),intent(out) :: fc(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine make_fc_paw
end interface

interface
 subroutine nhatgrid(atindx1,gmet,mpi_enreg,natom,natom_tot,nattyp,ngfft,ntypat,&  
  &  optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,ucvol,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natom_tot
  integer,intent(in) :: ntypat
  integer,intent(in) :: optcut
  integer,intent(in) :: optgr0
  integer,intent(in) :: optgr1
  integer,intent(in) :: optgr2
  integer,intent(in) :: optrad
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine nhatgrid
end interface

interface
 subroutine optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,eigen0,gprimd,hdr,indlmn,kg,lmnmax,&  
  &  mband,mcg,mcprj,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nsppol,&  
  &  pawrad,pawtab,wffnow)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cg(2,mcg)
  type(cprj_type) :: cprj(natom,mcprj)
  integer,intent(in) :: dimcprj(natom)
  real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,dtset%ntypat)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: npwarr(nkpt)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
 end subroutine optics_paw
end interface

interface
 subroutine optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,eigen0,hdr,indlmn,lmnmax,&  
  &  mband,mcprj,mkmem,mpi_enreg,mpsang,natom,nkpt,nsppol,pawrad,pawtab)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mband
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: atindx1(natom)
  type(cprj_type) :: cprj(natom,mcprj)
  integer,intent(in) :: dimcprj(natom)
  real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
  integer,intent(in) :: indlmn(6,lmnmax,dtset%ntypat)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
 end subroutine optics_paw_core
end interface

interface
 subroutine partial_dos_fractions_paw(atindx1,cprj,dimcprj,dos_fractions,dos_fractions_m,&  
  &  dos_fractions_paw1,dos_fractions_pawt1,&  
  &  dtfil,dtset,fatbands_flag,indlmn,lmnmax,mbesslang,mcprj,mkmem,&  
  &  mpi_enreg,m_dos_flag,ndosfraction,paw_dos_flag,pawrad,pawtab)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: fatbands_flag
  integer,intent(in) :: lmnmax
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: paw_dos_flag
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: atindx1(dtset%natom)
  type(cprj_type) :: cprj(dtset%natom,mcprj)
  integer,intent(in) :: dimcprj(dtset%natom)
  real(dp),intent(inout) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(inout) :: dos_fractions_m(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*mbesslang*max(m_dos_flag,fatbands_flag))
  real(dp),intent(out) :: dos_fractions_paw1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(out) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  integer,intent(in) :: indlmn(6,lmnmax,dtset%ntypat)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
 end subroutine partial_dos_fractions_paw
end interface

interface
 subroutine paw_mknewh0(nsppol,nspden,nfftf,pawspnorb,pawprtvol,Cryst,Psps,&  
  &  Pawtab,Paw_an,Paw_ij,Pawang,Pawfgrtab,vxc,vxc_val,vtrial)
  use defs_basis
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawspnorb
  type(crystal_structure),intent(in) :: Cryst
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(paw_an_type),intent(in) :: Paw_an(Cryst%natom)
  type(paw_ij_type),intent(inout) :: Paw_ij(Cryst%natom)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
  real(dp),intent(in) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
  real(dp),intent(in) :: vxc_val(nfftf,nspden)
 end subroutine paw_mknewh0
end interface

interface
 subroutine paw_symcprj(ik_bz,nspinor,nband_k,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_bz) 
  use m_bz_mesh
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: nband_k
  integer,intent(in) :: nspinor
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(cprj_type),intent(inout) :: Cprj_bz(Cryst%natom,nspinor*nband_k)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
 end subroutine paw_symcprj
end interface

interface
 subroutine paw_symcprj_op(ik_bz,nspinor,nband_k,Cryst,Kmesh,Psps,Pawtab,Pawang,in_Cprj,out_Cprj) 
  use m_bz_mesh
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: nband_k
  integer,intent(in) :: nspinor
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
  type(cprj_type),intent(in) :: in_Cprj(Cryst%natom,nspinor*nband_k)
  type(cprj_type),intent(out) :: out_Cprj(Cryst%natom,nspinor*nband_k)
 end subroutine paw_symcprj_op
end interface

interface
 subroutine pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj1,ipert,isppol,natom,&  
  &  nspinor,occ_k,option,pawrhoij,usetimerev,wtk_k)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ipert
  integer,intent(in) :: isppol
  integer,intent(in) :: natom
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  real(dp),intent(in) :: occ_k
  logical,intent(in) :: usetimerev
  real(dp),intent(in) :: wtk_k
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(in) :: cwaveprj(natom,nspinor)
  type(cprj_type),intent(in) :: cwaveprj1(natom,nspinor)
  type(pawrhoij_type),intent(inout) :: pawrhoij(natom)
 end subroutine pawaccrhoij
end interface

interface
 subroutine pawalloc(dtset,idtset,mpsang,mqgrid_vl,npsp,option,paw_size,paw_size_old,&  
  &  pawang,pawrad,pawtab,pspheads)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: idtset
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid_vl
  integer,intent(in) :: npsp
  integer,intent(in) :: option
  integer,intent(in) :: paw_size
  integer,intent(in) :: paw_size_old
  type(dataset_type),intent(in) :: dtset
  type(pawang_type),intent(inout) :: pawang
  type(pawrad_type),intent(inout) :: pawrad(paw_size)
  type(pawtab_type),intent(inout) :: pawtab(paw_size)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine pawalloc
end interface

interface
 subroutine pawdenpot(compch_sph,epaw,epawdc,ipert,ixc,mpi_enreg,&  
  &  natom,natom_tot,nspden,ntypat,nzlmopt,option,paral_kgb,paw_an,paw_an0,&  
  &  paw_ij,pawang,pawprtvol,pawrad,pawrhoij,pawspnorb,pawtab,pawxcdev,spnorbscl,xclevel,xc_denpos,znucl,&  
  &  electronpositron) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: ipert
  integer,intent(in) :: ixc
  integer,intent(in) :: natom
  integer,intent(in) :: natom_tot
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: nzlmopt
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: xclevel
  real(dp),intent(out) :: compch_sph
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(out) :: epaw
  real(dp),intent(out) :: epawdc
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp), intent(in) :: spnorbscl
  real(dp), intent(in) :: xc_denpos
  type(paw_an_type),intent(inout) :: paw_an(natom)
  type(paw_an_type),intent(inout) :: paw_an0(natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp) :: znucl(ntypat)
 end subroutine pawdenpot
end interface

interface
 subroutine pawdensities(compch_sph,cplex,iatom,lmselectin,lmselectout,lm_size,nhat1,nspden,nzlmopt,&  
  &  opt_compch,opt_dens,opt_l,opt_print,pawang,pawprtvol,pawrad,pawrhoij,pawtab,rho1,trho1,&  
  &  one_over_rad2) ! optional
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: iatom
  integer,intent(in) :: lm_size
  integer,intent(in) :: nspden
  integer,intent(in) :: nzlmopt
  integer,intent(in) :: opt_compch
  integer,intent(in) :: opt_dens
  integer,intent(in) :: opt_l
  integer,intent(in) :: opt_print
  integer,intent(in) :: pawprtvol
  real(dp),intent(out) :: compch_sph
  type(pawang_type),intent(in) :: pawang
  type(pawrad_type),intent(in) :: pawrad
  type(pawrhoij_type),intent(in) :: pawrhoij
  type(pawtab_type),intent(in) :: pawtab
  logical,intent(in) :: lmselectin(lm_size)
  logical,intent(out) :: lmselectout(lm_size)
  real(dp),intent(out) :: nhat1(cplex*pawrad%mesh_size,lm_size,nspden*(1-((opt_dens+1)/2)))
  real(dp),intent(in),target,optional :: one_over_rad2(pawrad%mesh_size)
  real(dp),intent(out) :: rho1(cplex*pawrad%mesh_size,lm_size,nspden)
  real(dp),intent(out) :: trho1(cplex*pawrad%mesh_size,lm_size,nspden*(1-(opt_dens/2)))
 end subroutine pawdensities
end interface

interface
 subroutine pawdij(cplex,dtset,enunit,fatvshift,gprimd,ipert,mpi_enreg,natom,natom_tot,nfft,ngfft,&  
  &  nspden,ntypat,paral_kgb,paw_an,paw_ij,pawang,pawfgrtab,pawprtvol,pawrad,&  
  &  pawspnorb,pawtab,pawxcdev,qphon,typat,ucvol,vtrial,vxc,xred,&  
  &  electronpositron) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: enunit
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: natom_tot
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  real(dp),intent(in) :: fatvshift
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  type(paw_an_type),intent(in) :: paw_an(natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: qphon(3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vtrial(cplex*nfft,nspden)
  real(dp),intent(in) :: vxc(cplex*nfft,nspden)
  real(dp),intent(in) :: xred(3,natom_tot)
 end subroutine pawdij
end interface

interface
 subroutine pawdijhartree(cplex,iatom,natom,ntypat,paw_ij,pawrhoij,pawtab)
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: iatom
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
 end subroutine pawdijhartree
end interface

interface
 subroutine pawdijso(iatom,itypat,natom,ntypat,paw_an,paw_ij,pawang,pawrad,pawtab,pawxcdev,spnorbscl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iatom
  integer,intent(in) :: itypat
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawxcdev
  type(pawang_type),intent(in) :: pawang
  real(dp), intent(in) :: spnorbscl
  type(paw_an_type),intent(in) :: paw_an(natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
 end subroutine pawdijso
end interface

interface
 subroutine pawexpiqr(gprimd,pawfgrtab,qphon,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  type(pawfgrtab_type),intent(inout) :: pawfgrtab
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: xred(3)
 end subroutine pawexpiqr
end interface

interface
 subroutine pawfrnhat(cplex,gprimd,idir,ipert,mpi_enreg,natom,natom_tot,nfft,ngfft,nspden,ntypat,&  
  &  optfr,paw_ij1,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,qphon,rprimd,&  
  &  ucvol,vpsp1,vtrial,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: natom_tot
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: optfr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  type(paw_ij_type),intent(inout) :: paw_ij1(natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vpsp1(cplex*nfft*(optfr+1)/2)
  real(dp),intent(in) :: vtrial(nfft,nspden*(optfr+1)/2)
  real(dp),intent(in) :: xred(3,natom_tot)
 end subroutine pawfrnhat
end interface

interface
 subroutine pawfrnhat_recipspace(atindx,cplex,gmet,gsqcut,idir,ipert,mgfft,mpi_enreg,&  
  &  mqgrid,natom,nattyp,nfft,nfftot,ngfft,nhatfr,nspden,ntypat,paral_kgb,&  
  &  paw_ij1,pawang,pawrhoij,pawtab,ph1d,qgrid,qphon,typat,ucvol,usepaw,vtrial,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: mgfft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(out) :: nhatfr(cplex*nfft,nspden)
  type(paw_ij_type),intent(inout) :: paw_ij1(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: qphon(3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vtrial(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine pawfrnhat_recipspace
end interface

interface
 subroutine pawgrnl(atindx1,dimnhat,dimvtrial,dyfrnl,dyfr_cplex,grnl,gsqcut,mgfft,mpi_enreg,natom,natom_tot,&  
  &  nattyp,nfft,ngfft,nhat,nlstr,nspden,nsym,ntypat,optgr,optgr2,optstr,paral_kgb,&  
  &  pawang,pawfgrtab,pawrhoij,pawtab,ph1d,psps,qphon,rprimd,symrec,typat,vtrial,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimnhat
  integer,intent(in) :: dimvtrial
  integer,intent(in) :: dyfr_cplex
  integer,intent(in) :: mgfft
  integer,intent(in) :: natom
  integer,intent(in) :: natom_tot
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: optgr
  integer,intent(in) :: optgr2
  integer,intent(in) :: optstr
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom_tot)
  real(dp),intent(inout) :: dyfrnl(dyfr_cplex,3,3,natom_tot,natom_tot*optgr2)
  real(dp),intent(inout) :: grnl(3*natom_tot*optgr)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfft,dimnhat)
  real(dp),intent(inout) :: nlstr(6*optstr)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom_tot)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom_tot)
  real(dp),intent(in) :: vtrial(nfft,dimvtrial)
  real(dp),intent(in) :: xred(3,natom_tot)
 end subroutine pawgrnl
end interface

interface
 subroutine pawgylm(gylm,gylmgr,gylmgr2,lm_size,nfgd,optgr0,optgr1,optgr2,pawtab,rfgd,rfgd_allocated)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lm_size
  integer,intent(in) :: nfgd
  integer,intent(in) :: optgr0
  integer,intent(in) :: optgr1
  integer,intent(in) :: optgr2
  integer,intent(in) :: rfgd_allocated
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(out) :: gylm(nfgd,optgr0*lm_size)
  real(dp),intent(out) :: gylmgr(3,nfgd,optgr1*lm_size)
  real(dp),intent(out) :: gylmgr2(6,nfgd,optgr2*lm_size)
  real(dp),intent(in) :: rfgd(3,nfgd)
 end subroutine pawgylm
end interface

interface
 subroutine pawgylmg(gprimd,gylmg,kg,kpg,kpt,lmax,nkpg,npw,ntypat,pawtab,ylm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: nkpg
  integer,intent(in) :: npw
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gylmg(npw,lmax**2,ntypat)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpg(npw,nkpg)
  real(dp),intent(in) :: kpt(3)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ylm(npw,lmax**2)
 end subroutine pawgylmg
end interface

interface
 subroutine pawinit(ecutshp_eff,indlmn,lcutdens,lmix,lmnmax,mpsang,nphi,nsym,ntheta,ntypat,&  
  &  pawang,pawrad,pawspnorb,pawtab,pawxcdev)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lcutdens
  integer,intent(in) :: lmix
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: nphi
  integer,intent(in) :: nsym
  integer,intent(in) :: ntheta
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawxcdev
  real(dp) :: ecutshp_eff
  type(pawang_type),intent(inout) :: pawang
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(inout) :: pawtab(ntypat)
 end subroutine pawinit
end interface

interface
 subroutine pawlsylm(pawang)
  use defs_datatypes
  implicit none
  type(pawang_type),intent(inout) :: pawang
 end subroutine pawlsylm
end interface

interface
 subroutine pawmkaewf(Dtset,natom,mpw,mband,mcg,mcprj,nkpt,mkmem,nsppol,ntypat,nband,istwfk,npwarr,kpt,&  
  &  paral_kgb,ngfftf,kg,dimcprj,Pawfgrtab,Pawrad,Pawtab,gmet,rprimd,ucvol,&  
  &  Psps,Hdr,Dtfil,typat,eigen,occ,cg,Cprj,Wffnow,MPI_enreg,ierr,pseudo_norms,set_k,set_band)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in),optional :: set_band
  integer,intent(in),optional :: set_k
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(hdr_type),intent(inout) :: Hdr
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pseudopotential_type),intent(in) :: Psps
  type(wffile_type),intent(inout) :: Wffnow
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftf(18)
  type(cprj_type),intent(in) :: Cprj(natom,mcprj)
  type(pawfgrtab_type),intent(in) :: Pawfgrtab(natom)
  type(pawrad_type),intent(in) :: Pawrad(ntypat)
  type(pawtab_type),intent(in) :: Pawtab(ntypat)
  real(dp),intent(in) :: cg(2,mcg)
  integer,intent(in) :: dimcprj(natom)
  real(dp),target,intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),target,intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),optional,intent(out) :: pseudo_norms(nsppol,nkpt,mband)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
 end subroutine pawmkaewf
end interface

interface
 subroutine pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,&  
  &  natom,natom_tot,nfft,ngfft,nhatgrdim,nspden,ntypat,paral_kgb,pawang,pawfgrtab,&  
  &  pawgrnhat,pawnhat,pawrhoij,pawrhoij0,pawtab,qphon,rprimd,ucvol,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ider
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: izero
  integer,intent(in) :: natom
  integer,intent(in) :: natom_tot
  integer,intent(in) :: nfft
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  real(dp),intent(out) :: compch_fft
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  real(dp),intent(out) :: pawgrnhat(cplex*nfft,nspden,3*nhatgrdim)
  real(dp),intent(out) :: pawnhat(cplex*nfft,nspden)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij0(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom_tot)
 end subroutine pawmknhat
end interface

interface
 subroutine pawmknhat_psipsi(cprj1,cprj2,ider,izero,mpi_enreg,natom,nfft,ngfft,nhat12_grdim,&  
  &  nspinor,ntypat,typat,paral_kgb,pawang,pawfgrtab,grnhat12,nhat12,pawtab)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ider
  integer,intent(in) :: izero
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nhat12_grdim
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: ngfft(18)
  type(cprj_type),intent(in) :: Cprj1(natom,nspinor)
  type(cprj_type),intent(in) :: Cprj2(natom,nspinor)
  real(dp),intent(out) :: grnhat12(2,nfft,nspinor**2,3*nhat12_grdim)
  real(dp),intent(out) :: nhat12(2,nfft,nspinor**2)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawmknhat_psipsi
end interface

interface
 subroutine pawmkrho(compch_fft,cplex,gprimd,idir,indlmn,indsym,ipert,lmnmax,mpi_enreg,&  
  &  natom,nspden,nsym,ntypat,paral_kgb,pawang,pawfgr,pawfgrtab,pawprtvol,pawrhoij,&  
  &  pawtab,qphon,rhopsg,rhopsr,rhor,rprimd,symafm,symrec,typat,ucvol,xred,&  
  &  pawang_sym,pawnhat,pawrhoij0,rhog) ! optional arguments
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: pawprtvol
  real(dp),intent(out) :: compch_fft
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawang_type),intent(in),optional :: pawang_sym
  type(pawfgr_type),intent(in) :: pawfgr
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  real(dp),intent(out),target,optional :: pawnhat(cplex*pawfgr%nfft,nspden)
  type(pawrhoij_type),intent(inout),target :: pawrhoij(natom)
  type(pawrhoij_type),intent(in),target,optional :: pawrhoij0(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out),optional :: rhog(2,pawfgr%nfft)
  real(dp),intent(inout) :: rhopsg(2,pawfgr%nfftc)
  real(dp),intent(inout) :: rhopsr(cplex*pawfgr%nfftc,nspden)
  real(dp),intent(out) :: rhor(cplex*pawfgr%nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine pawmkrho
end interface

interface
 subroutine pawmkrhoij(atindx1,cprj,dimcprj,istwfk,kptopt,mband,mband_cprj,mcprj,mkmem,mpi_enreg,&  
  &  natom,nband,nkpt,nspinor,nsppol,occ,paral_kgb,paw_dmft,&  
  &  pawprtvol,pawrhoij,unpaw,wtk)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_paw_dmft
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: mband
  integer,intent(in) :: mband_cprj
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: unpaw
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type),intent(in) :: paw_dmft
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),target,intent(in) :: cprj(natom,mcprj)
  integer,intent(in) :: dimcprj(natom)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawrhoij_type),intent(inout) :: pawrhoij(natom)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine pawmkrhoij
end interface

interface
 subroutine pawnabla_init(mpsang,lmnmax,ntypat,indlmn,pawrad,pawtab)
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: ntypat
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(inout) :: pawtab(ntypat)
 end subroutine pawnabla_init
end interface

interface
 subroutine pawnstd2e(epawnst,ipert1,ipert2,mpi_enreg,natom,natom_tot,ntypat,nzlmopt1,nzlmopt2,&  
  &  paw_an0,paw_an1,paw_ij1,pawang,pawprtvol,pawrad,pawrhoij1,pawrhoij2,&  
  &  pawtab,pawxcdev,xclevel)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipert1
  integer,intent(in) :: ipert2
  integer,intent(in) :: natom
  integer,intent(in) :: natom_tot
  integer,intent(in) :: ntypat
  integer,intent(in) :: nzlmopt1
  integer,intent(in) :: nzlmopt2
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: pawxcdev
  integer,intent(in) :: xclevel
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(out) :: epawnst(2)
  type(paw_an_type),intent(in) :: paw_an0(natom)
  type(paw_an_type),intent(inout) :: paw_an1(natom)
  type(paw_ij_type),intent(inout) :: paw_ij1(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij1(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij2(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
 end subroutine pawnstd2e
end interface

interface
 subroutine pawpolev(natom,ntypat,pawrhoij,pawtab,pelev,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(out) :: pelev(3)
  integer,intent(in) :: typat(natom)
 end subroutine pawpolev
end interface

interface
 subroutine pawprt(dtset,indlmn,lmnmax,paw_ij,pawrhoij,pawtab,&  
  &  electronpositron) ! optional argument
  use defs_abitypes
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: lmnmax
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer,optional :: electronpositron
  integer,intent(in) :: indlmn(6,lmnmax,dtset%ntypat)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
 end subroutine pawprt
end interface

interface
 subroutine pawpupot(ispden,paw_ij,pawprtvol,pawtab,vpawu,VUKS)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ispden
  integer,intent(in) :: pawprtvol
  real(dp),intent(inout) :: VUKS
  type(paw_ij_type),intent(in) :: paw_ij
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(inout) :: vpawu(paw_ij%cplex_dij,pawtab%lpawu*2+1,pawtab%lpawu*2+1)
 end subroutine pawpupot
end interface

interface
 subroutine pawpuxinit(dmatpuopt,exchmix,jpawu,llexexch,llpawu,indlmn,lmnmax,ntypat,pawang,&  
  &  pawprtvol,pawrad,pawtab,upawu,use_dmft,useexexch,usepawu)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dmatpuopt
  integer,intent(in) :: lmnmax
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: use_dmft
  integer,intent(in) :: useexexch
  integer,intent(in) :: usepawu
  real(dp),intent(in) :: exchmix
  type(pawang_type), intent(in) :: pawang
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  real(dp),intent(in) :: jpawu(ntypat)
  integer,intent(in) :: llexexch(ntypat)
  integer,intent(in) :: llpawu(ntypat)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(inout) :: pawtab(ntypat)
  real(dp),intent(in) :: upawu(ntypat)
 end subroutine pawpuxinit
end interface

interface
 subroutine pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,iband1,iband2,ispinor1,ispinor2,istwf_k,kg_diel,&  
  &  lmax_diel,mgfftdiel,mpi_enreg,natom,nband,ndiel4,ndiel5,ndiel6,&  
  &  ngfftdiel,npwdiel,nspinor,ntypat,optreal,paral_kgb,&  
  &  pawang,pawtab,ph3d_diel,typat,wfprod,wfraug)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband1
  integer,intent(in) :: iband2
  integer,intent(in) :: ispinor1
  integer,intent(in) :: ispinor2
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmax_diel
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: optreal
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx1(natom)
  type(cprj_type),intent(in) :: cprj_k(natom,nspinor*nband)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat)
  integer,intent(in) :: kg_diel(3,npwdiel)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: wfprod(2,npwdiel*(1-optreal))
  real(dp),intent(inout) :: wfraug(2,ndiel4,ndiel5,ndiel6*optreal)
 end subroutine pawsushat
end interface

interface
 subroutine pawtwdij(dtefield,gprimd,natom,nfftf,nspden,ntypat,&  
  &  paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,typat,&  
  &  vtrial,vxc)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  type(pawang_type), intent(in) :: pawang
  real(dp),intent(in) :: gprimd(3,3)
  type(paw_an_type),intent(in) :: paw_an(natom)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
  type(pawrad_type), intent(in) :: pawrad(ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vtrial(nfftf,nspden)
  real(dp),intent(in) :: vxc(nfftf,nspden)
 end subroutine pawtwdij
end interface

interface
 subroutine pawtwdij_1(dtefield,gprimd,natom,ntypat,pawrad,pawtab,psps,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawtwdij_1
end interface

interface
 subroutine pawtwdij_2a(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawtwdij_2a
end interface

interface
 subroutine pawtwdij_2b(dtefield,gprimd,natom,ntypat,pawrad,pawtab,psps,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawtwdij_2b
end interface

interface
 subroutine pawtwdij_2c(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawtwdij_2c
end interface

interface
 subroutine pawtwdij_2d(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawtwdij_2d
end interface

interface
 subroutine pawtwdij_2e(dtefield,gprimd,natom,ntypat,pawrad,pawtab,psps,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawtwdij_2e
end interface

interface
 subroutine pawtwdij_2f(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine pawtwdij_2f
end interface

interface
 subroutine pawuenergy(iatom,eldaumdc,eldaumdcdc,pawprtvol,pawtab,paw_ij&  
  &  ,e_ee,e_dc,e_dcdc,dmft_dc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, optional, intent(in) :: dmft_dc
  integer,intent(in) :: iatom
  integer,intent(in) :: pawprtvol
  real(dp), optional, intent(inout) :: e_dc
  real(dp), optional, intent(inout) :: e_dcdc
  real(dp), optional, intent(inout) :: e_ee
  real(dp),intent(inout) :: eldaumdc
  real(dp),intent(inout) :: eldaumdcdc
  type(paw_ij_type),intent(in) :: paw_ij
  type(pawtab_type),intent(in) :: pawtab
 end subroutine pawuenergy
end interface

interface
 subroutine pawuj_det(dtpawuj,ndtpawuj,ujdet_filename,ures)
  use defs_basis
  use defs_datatypes
  implicit none
  integer :: ndtpawuj
  character(len=*),intent(in) :: ujdet_filename
  real(dp),intent(out) :: ures
  type(macro_uj_type),intent(in) :: dtpawuj(0:ndtpawuj)
 end subroutine pawuj_det
end interface

interface
 subroutine pawuj_ini(dtpawuj,ndtset)
  use defs_datatypes
  implicit none
  integer :: ndtset
  type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtset)
 end subroutine pawuj_ini
end interface

interface
 subroutine pawuj_red(dtset,dtpawuj,fatvshift,natom,ntypat,paw_ij,pawrad,pawtab,ndtpawuj)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndtpawuj
  integer,intent(in) :: ntypat
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fatvshift
  type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)
  type(paw_ij_type),intent(in) :: paw_ij(natom)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
 end subroutine pawuj_red
end interface

interface
 subroutine linvmat(inmat,oumat,nat,nam,option,gam,prtvol)
  use defs_basis
  implicit none
  integer,intent(in) :: nat
  integer,intent(in),optional :: option
  integer,intent(in),optional :: prtvol
  real(dp),intent(in) :: gam
  character(len=500),intent(in) :: nam
  real(dp),intent(in) :: inmat(nat,nat)
  real(dp),intent(inout) :: oumat(:,:)
 end subroutine linvmat
end interface

interface
 subroutine lprtmat(commnt,chan,prtvol,mmat,nat)
  use defs_basis
  implicit none
  integer,intent(in) :: chan
  integer,intent(in) :: nat
  integer,intent(in) :: prtvol
  character(len=500),intent(in) :: commnt
  real(dp),intent(in) :: mmat(nat,nat)
 end subroutine lprtmat
end interface

interface
 subroutine lcalcu(magv,natom,rprimd,xred,chi,chi0,pawujat,ures,prtvol,gam,opt)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in),optional :: opt
  integer,intent(in),optional :: pawujat
  integer,intent(in),optional :: prtvol
  real(dp),intent(in),optional :: gam
  real(dp),intent(out) :: ures
  real(dp),intent(in) :: chi(natom)
  real(dp),intent(in) :: chi0(natom)
  integer,intent(in) :: magv(natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine lcalcu
end interface

interface
 subroutine pawuj_nullify(dtpawuj)
  use defs_datatypes
  implicit none
  type(macro_uj_type),intent(out) :: dtpawuj
 end subroutine pawuj_nullify
end interface

interface
 subroutine pawuj_free(dtpawuj)
  use defs_datatypes
  implicit none
  type(macro_uj_type),intent(inout) :: dtpawuj
 end subroutine pawuj_free
end interface

interface
 subroutine pawxenergy(eexex,pawprtvol,pawrhoij,pawtab)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: pawprtvol
  real(dp),intent(inout) :: eexex
  type(pawrhoij_type),intent(in) :: pawrhoij
  type(pawtab_type),intent(in) :: pawtab
 end subroutine pawxenergy
end interface

interface
 subroutine pawxpot(pawprtvol,pawtab,paw_ij,pawrhoij)
  use defs_datatypes
  implicit none
  integer,intent(in) :: pawprtvol
  type(paw_ij_type),intent(inout) :: paw_ij
  type(pawrhoij_type),intent(in) :: pawrhoij
  type(pawtab_type),intent(in) :: pawtab
 end subroutine pawxpot
end interface

interface
 subroutine prtfatbands(dos_fractions_m,dtset,fildata,fermie,eigen,&  
  &  mbesslang,m_dos_flag,ndosfraction,pawfatbnd,pawtab)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: pawfatbnd
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: fildata
  real(dp),intent(in) :: dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)
 end subroutine prtfatbands
end interface

interface
 subroutine qijb_bk(dtefield,gprimd,natom,ntypat,pawang,pawrad,pawtab,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine qijb_bk
end interface

interface
 subroutine qijb_kk(dtefield,gprimd,natom,ntypat,&  
  &  pawang,pawrad,pawtab,typat)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(efield_type),intent(inout) :: dtefield
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine qijb_kk
end interface

interface
 subroutine read_atomden(MPI_enreg,natom,nspden,ntypat,pawfgr,&  
  &  rhor_paw,typat,rprimd,xred,prtvol,file_prefix)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  type(mpi_type),intent(in) :: MPI_enreg
  character(len=7), intent(in) :: file_prefix
  type(pawfgr_type),intent(in) :: pawfgr
  real(dp),intent(inout) :: rhor_paw(pawfgr%nfft,nspden)
  real(dp), intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp), intent(in) :: xred(3,natom)
 end subroutine read_atomden
end interface

interface
 subroutine setnoccmmp(compute_dmat,dimdmat,dmatpawu,dmatudiag,impose_dmat,indsym,natom,natpawu,&  
  &  nspinor,nsppol,nsym,ntypat,paw_ij,pawang,pawprtvol,pawrhoij,pawtab,&  
  &  spinat,symafm,useexexch,usepawu)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: compute_dmat
  integer,intent(in) :: dimdmat
  integer,intent(in) :: dmatudiag
  integer,intent(in) :: impose_dmat
  integer,intent(in) :: natom
  integer,intent(in) :: natpawu
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  integer,intent(in) :: useexexch
  integer,intent(in) :: usepawu
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: dmatpawu(dimdmat,dimdmat,nspinor*nsppol,natpawu*impose_dmat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(in) :: symafm(nsym)
 end subroutine setnoccmmp
end interface

interface
 subroutine setrhoijpbe0(dtset,indlmn,initialized,istep,istep_mix,lmnmax,natom,ntypat,pawrhoij,pawtab,typat)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: initialized
  integer,intent(in) :: istep
  integer,intent(inout) :: istep_mix
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(dataset_type),intent(in) :: dtset
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(natom)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  integer,intent(in) :: typat(natom)
 end subroutine setrhoijpbe0
end interface

interface
 subroutine setsymrhoij(gprimd,lmax,nsym,pawprtvol,rprimd,sym,zarot)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: nsym
  integer,intent(in) :: pawprtvol
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: sym(3,3,nsym)
  real(dp),intent(out) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)
 end subroutine setsymrhoij
end interface

interface
 subroutine simple_j_dia(jdia,natom,nfft,pawfgrtab)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  real(dp),intent(out) :: jdia(3,3,nfft)
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
 end subroutine simple_j_dia
end interface

interface
 subroutine smatrix_k_paw(cprj_k,cprj_kb,dtefield,kdir,kfor,natom,smat_k_paw,typat,bdir,bfor)
  use defs_basis
  use m_efield
  use defs_datatypes
  implicit none
  integer,optional,intent(in) :: bdir
  integer,optional,intent(in) :: bfor
  integer,intent(in) :: kdir
  integer,intent(in) :: kfor
  integer,intent(in) :: natom
  type(efield_type),intent(in) :: dtefield
  type(cprj_type),intent(in) :: cprj_k(natom,dtefield%nband_occ)
  type(cprj_type),intent(in) :: cprj_kb(natom,dtefield%nband_occ)
  real(dp),intent(out) :: smat_k_paw(2,dtefield%nband_occ,dtefield%nband_occ)
  integer,intent(in) :: typat(natom)
 end subroutine smatrix_k_paw
end interface

interface
 subroutine smatrix_paw(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&  
  &  mcg_k,mcg_q,mcg1_k,minbd,mpw,nband_occ,npw_k1,npw_k2,nspinor,&  
  &  pwind_k,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_k)
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
  real(dp),intent(in) :: cg(2,mcg_k)
  real(dp),intent(out) :: cg1_k(2,mcg1_k)
  real(dp),intent(in) :: cgq(2,mcg_q)
  real(dp),intent(out) :: dtm_k(2)
  integer,intent(in) :: pwind_k(mpw)
  real(dp),intent(in) :: pwnsfac_k(4,mpw)
  integer,intent(inout) :: sflag_k(nband_occ)
  real(dp),intent(out) :: smat_inv(2,nband_occ,nband_occ)
  real(dp),intent(inout) :: smat_k(2,nband_occ,nband_occ)
 end subroutine smatrix_paw
end interface

interface
 subroutine smatrix_pawinit(atindx1,cm2,cprj,ikpt1,ikpt2,isppol,&  
  &  g1,gprimd,kpt,mband,mbandw,mkmem,mpi_enreg,&  
  &  natom,nband,nkpt,nspinor,nsppol,ntypat,pawang,pawrad,pawtab,rprimd,&  
  &  seed_name,typat,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ikpt1
  integer,intent(in) :: ikpt2
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mbandw
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  character(len=fnlen) :: seed_name
  integer,intent(in) :: g1(3)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cm2(2,mbandw,mbandw)
  type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: nband(nsppol*nkpt)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine smatrix_pawinit
end interface

interface
 subroutine spline_paw_fncs(dphi,dtphi,nnl,npts,pawrad,pawtab,points,phi,tphi)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nnl
  integer,intent(in) :: npts
  type(pawrad_type),intent(in) :: pawrad
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(out) :: dphi(npts,nnl)
  real(dp),intent(out) :: dtphi(npts,nnl)
  real(dp),intent(out) :: phi(npts,nnl)
  real(dp),intent(in) :: points(npts)
  real(dp),intent(out) :: tphi(npts,nnl)
 end subroutine spline_paw_fncs
end interface

interface
 subroutine sym_cprj_kn(cprj_fkn,cprj_ikn,cprj_sym,dimlmn,iband,indlmn,&  
  &  isym,itim,kpt,lmax,lmnmax,natom,nband,nspinor,nsym,ntypat,&  
  &  typat,zarot)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iband
  integer,intent(in) :: isym
  integer,intent(in) :: itim
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(cprj_type),intent(out) :: cprj_fkn(natom,nband*nspinor)
  type(cprj_type),intent(in) :: cprj_ikn(natom,nband*nspinor)
  integer,intent(in) :: cprj_sym(4,nsym,natom)
  integer,intent(in) :: dimlmn(natom)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  real(dp),intent(in) :: kpt(3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)
 end subroutine sym_cprj_kn
end interface

interface
 subroutine symdij(gprimd,indlmn,indsym,ipert,lmnmax,natom,nsym,ntypat,option_dij,&  
  &  paw_ij,pawang,pawprtvol,rprimd,symafm,symrec,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipert
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: option_dij
  integer,intent(in) :: pawprtvol
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine symdij
end interface

interface
 subroutine symdij_all(gprimd,indlmn,indsym,ipert,lmnmax,natom,nsym,ntypat,&  
  &  paw_ij,pawang,pawprtvol,rprimd,symafm,symrec,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipert
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawprtvol
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine symdij_all
end interface

interface
 subroutine symrhoij(choice,gprimd,indlmn,indsym,ipert,lmnmax,natom,nrhoij,nsym,ntypat,optrhoij,&  
  &  pawang,pawprtvol,pawrhoij,rprimd,symafm,symrec,typat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: ipert
  integer,intent(in) :: lmnmax
  integer,intent(in) :: natom
  integer,intent(in) :: nrhoij
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: optrhoij
  integer,intent(in) :: pawprtvol
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  integer,intent(in) :: indsym(4,nsym,natom)
  type(pawrhoij_type),intent(inout) :: pawrhoij(nrhoij)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: typat(natom)
 end subroutine symrhoij
end interface

interface
 subroutine transgrid(cplex,mpi_enreg,nspden,optgrid,optin,optout,paral_kgb,pawfgr,rhog,rhogf,rhor,rhorf)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nspden
  integer,intent(in) :: optgrid
  integer,intent(in) :: optin
  integer,intent(in) :: optout
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawfgr_type),intent(in) :: pawfgr
  real(dp),intent(inout) :: rhog(2,pawfgr%nfftc)
  real(dp),intent(inout) :: rhogf(2,pawfgr%nfft)
  real(dp),intent(inout) :: rhor(cplex*pawfgr%nfftc,nspden)
  real(dp),intent(inout) :: rhorf(cplex*pawfgr%nfft,nspden)
 end subroutine transgrid
end interface

end module interfaces_66_paw
!!***
