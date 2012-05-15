!!****m* ABINIT/interfaces_70_gw
!! NAME
!! interfaces_70_gw
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/70_gw
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

module interfaces_70_gw

 implicit none

interface
 subroutine accumulate_chi0_q0(ik_bz,isym_kbz,itim_kbz,gwcomp,gw_eet,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&  
  &  chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2,chi0_head,chi0_lwing,chi0_uwing)
  use m_gsphere
  use m_bz_mesh
  use defs_basis
  use m_gwdefs
  use m_crystal
  implicit none
  integer,intent(in) :: gw_eet
  integer,intent(in) :: gwcomp
  integer,intent(in) :: ik_bz
  integer,intent(in) :: isym_kbz
  integer,intent(in) :: itim_kbz
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nspinor
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: deltaf_b1b2
  complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  complex(dpc),intent(inout) :: chi0_head(3,3,Ep%nomega)
  complex(dpc),intent(inout) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
  complex(dpc),intent(inout) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)
  complex(dpc),intent(in) :: green_enhigh_w(Ep%nomega)
  complex(dpc),intent(in) :: green_w(Ep%nomega)
  complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
  complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 end subroutine accumulate_chi0_q0
end interface

interface
 function q0limit(ii,qlwl,nspinor,rhotwx,b1,b2,b3)
  use defs_basis
  implicit none
  integer,intent(in) :: ii
  integer,intent(in) :: nspinor
  complex(gwpc) :: q0limit
  real(dp),intent(in) :: b1(3)
  real(dp),intent(in) :: b2(3)
  real(dp),intent(in) :: b3(3)
  real(dp),intent(in) :: qlwl(3)
  complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 end function q0limit
end interface

interface
 subroutine accumulate_sfchi0_q0(ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&  
  &  factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,sf_chi0,sf_head,sf_lwing,sf_uwing)
  use m_bz_mesh
  use defs_basis
  use m_gsphere
  use m_crystal
  implicit none
  integer,intent(in) :: ikbz
  integer,intent(in) :: iomegal
  integer,intent(in) :: iomegar
  integer,intent(in) :: isym_kbz
  integer,intent(in) :: itim_kbz
  integer,intent(in) :: my_wl
  integer,intent(in) :: my_wr
  integer,intent(in) :: nomegasf
  integer,intent(in) :: npwe
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nspinor
  integer,intent(in) :: symchi
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: factocc
  real(dp),intent(in) :: wl
  real(dp),intent(in) :: wr
  complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
  complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
  complex(gwpc),intent(inout) :: sf_chi0(npwe,npwe,my_wl:my_wr)
  complex(dpc),intent(inout) :: sf_head(3,3,my_wl:my_wr)
  complex(dpc),intent(inout) :: sf_lwing(npwe,my_wl:my_wr,3)
  complex(dpc),intent(inout) :: sf_uwing(npwe,my_wl:my_wr,3)
 end subroutine accumulate_sfchi0_q0
end interface

interface
 subroutine assemblychi0sf(ik_bz,nspinor,symchi,Ltg_q,npwepG0,npwe,rhotwg,Gsph_epsG0,&  
  &  factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,nomegasf,chi0sf)
  use m_bz_mesh
  use m_gsphere
  use defs_basis
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: iomegal
  integer,intent(in) :: iomegar
  integer,intent(in) :: my_wl
  integer,intent(in) :: my_wr
  integer,intent(in) :: nomegasf
  integer,intent(in) :: npwe
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nspinor
  integer,intent(in) :: symchi
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: factocc
  real(dp),intent(in) :: wl
  real(dp),intent(in) :: wr
  complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)
  complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 end subroutine assemblychi0sf
end interface

interface
 subroutine calc_ffm(epsm1,nq,npw,nomega,omega,gprimd,qq,gvec,nfmidm)
  use defs_basis
  implicit none
  integer,intent(in) :: nfmidm
  integer,intent(in) :: nomega
  integer,intent(in) :: npw
  integer,intent(in) :: nq
  complex(gwpc),intent(in) :: epsm1(npw,npw,nomega,nq)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npw)
  complex(dpc),intent(in) :: omega(nomega)
  real(dp),intent(in) :: qq(3,nq)
 end subroutine calc_ffm
end interface

interface
 subroutine calc_rpa_functional(gwrpacorr,iqcalc,iq,Ep,Pvc,Qmesh,Dtfil,gmet,chi0,spaceComm,ec_rpa)
  use m_vcoul
  use m_bz_mesh
  use defs_abitypes
  use m_gwdefs
  use defs_basis
  implicit none
  integer,intent(in) :: gwrpacorr
  integer,intent(in) :: iq
  integer,intent(in) :: iqcalc
  integer,intent(in) :: spaceComm
  type(datafiles_type),intent(in) :: Dtfil
  type(epsilonm1_parameters),intent(in) :: Ep
  type(vcoul_t),intent(in) :: Pvc
  type(bz_mesh_type),intent(in) :: Qmesh
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  real(dp),intent(inout) :: ec_rpa(gwrpacorr)
  real(dp),intent(in) :: gmet(3,3)
 end subroutine calc_rpa_functional
end interface

interface
 subroutine calc_sig_ppm_comp(npwc,nomega,rhotwgp,botsq,otq,omegame0i_io,zcut,theta_mu_minus_e0i,ket,ppmodel,npwx,npwc1,npwc2)
  use defs_basis
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  integer,intent(in) :: npwx
  integer,intent(in) :: ppmodel
  real(dp),intent(in) :: omegame0i_io
  real(dp),intent(in) :: theta_mu_minus_e0i
  real(dp),intent(in) :: zcut
  complex(gwpc),intent(in) :: botsq(npwc,npwc1)
  complex(gwpc),intent(inout) :: ket(npwc,nomega)
  complex(gwpc),intent(in) :: otq(npwc,npwc2)
  complex(gwpc),intent(in) :: rhotwgp(npwx)
 end subroutine calc_sig_ppm_comp
end interface

interface
 subroutine calc_delta_ppm(Sigp,ptwsq,niter)
  use defs_basis
  use m_gwdefs
  implicit none
  integer,intent(in) :: niter
  type(sigma_parameters),intent(in) :: Sigp
  complex(gwpc),intent(inout) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 end subroutine calc_delta_ppm
end interface

interface
 subroutine calc_delta_ppm_sc(Sigp,nomega,otq,omegame0k,omegame0lumo,npwc2,&  
  &  qbzpg,ikbz,jkbz,ptwsq,delta,niter)
  use defs_basis
  use m_gwdefs
  implicit none
  integer :: ikbz
  integer :: jkbz
  integer,intent(in) :: niter
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc2
  type(sigma_parameters),intent(in) :: Sigp
  complex(gwpc),intent(out) :: delta(Sigp%npwc,Sigp%npwc,nomega)
  real(dp),intent(in) :: omegame0k(nomega)
  real(dp),intent(in) :: omegame0lumo(nomega)
  complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
  complex(gwpc),intent(in) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
  real(dp),intent(in) :: qbzpg(Sigp%npwc)
 end subroutine calc_delta_ppm_sc
end interface

interface
 subroutine calc_sig_ppm_delta(npwc,nomega,rhotwgp,botsq,otq,omegame0i,zcut,theta_mu_minus_e0i,&  
  &  ket,npwx,npwc1,npwc2,omega4sd,e0,delta)
  use defs_basis
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  integer,intent(in) :: npwx
  real(dp),intent(in) :: e0
  real(dp),intent(in) :: theta_mu_minus_e0i
  real(dp),intent(in) :: zcut
  complex(gwpc),intent(in) :: botsq(npwc,npwc1)
  complex(gwpc),intent(in) :: delta(npwc,npwc,nomega)
  complex(gwpc),intent(inout) :: ket(npwc,nomega)
  complex(dpc),intent(in) :: omega4sd(nomega)
  real(dp),intent(in) :: omegame0i(nomega)
  complex(gwpc),intent(in) :: otq(npwc,npwc2)
  complex(gwpc),intent(in) :: rhotwgp(npwx)
 end subroutine calc_sig_ppm_delta
end interface

interface
 subroutine calc_sig_ppm_delta_clos(npwc,nomega,ikbz,jkbz,qbzpg,botsq,otq,omegame0k,omegame0lumo,zcut,ptwsq,&  
  &  ket,npwc1,npwc2,gw_eet_scale,niter)
  use defs_basis
  implicit none
  integer,intent(in) :: ikbz
  integer,intent(in) :: jkbz
  integer,intent(in) :: niter
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  real(dp),intent(in) :: gw_eet_scale
  real(dp),intent(in) :: zcut
  complex(gwpc),intent(in) :: botsq(npwc,npwc1)
  complex(dpc),intent(inout) :: ket(nomega)
  real(dp),intent(in) :: omegame0k(nomega)
  real(dp),intent(in) :: omegame0lumo(nomega)
  complex(gwpc),intent(in) :: otq(npwc,npwc2)
  complex(gwpc),intent(in) :: ptwsq(npwc,npwc,niter+1)
  real(dp),intent(in) :: qbzpg(npwc)
 end subroutine calc_sig_ppm_delta_clos
end interface

interface
 subroutine gw_eet_sigma(Sigp,Sr,Dtset,Cryst,Wfs,Kmesh,Qmesh,Gsph_Max,Gsph_c,Psps,Vcp,QP_BSt,PPm,&  
  &  isppol,iq_bz,ik_bz,jk_bz,ik_ibz,jk_ibz,itim_q,isym_q,iq_ibz,tabr_ki,&  
  &  tabr_kj,spinrot_ki,spinrot_kj,ph_mkit,ph_mkjt,nfftot_gw,ngfft_gw,&  
  &  use_padfft,igfftcg0,gw_gbound,gw_mgfft,ib1,ib2,nomega_tot,nomega_sigc,&  
  &  fact_sp,nspinor,botsq,otq,sigcme_tmp,sigc,nbhomo,tim_fourdp,wtqp,wtqm,&  
  &  MPI_enreg,extrapolar_distrb,can_symmetrize)
  use m_vcoul
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use m_sigma_results
  use defs_abitypes
  use m_ppmodel
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: gw_mgfft
  integer,intent(in) :: ib1
  integer,intent(in) :: ib2
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: iq_bz
  integer,intent(in) :: iq_ibz
  integer,intent(in) :: isppol
  integer,intent(in) :: isym_q
  integer,intent(in) :: itim_q
  integer,intent(in) :: jk_bz
  integer,intent(in) :: jk_ibz
  integer,intent(out) :: nbhomo
  integer,intent(in) :: nfftot_gw
  integer,intent(in) :: nomega_sigc
  integer,intent(in) :: nomega_tot
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  integer,intent(in) :: wtqm
  integer,intent(in) :: wtqp
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(gvectors_type),intent(in) :: Gsph_Max
  type(gvectors_type),intent(in) :: Gsph_c
  type(bz_mesh_type),intent(in) :: Kmesh
  type(mpi_type),intent(inout) :: MPI_enreg
  type(ppmodel_type),intent(inout) :: PPm
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(bz_mesh_type),intent(in) :: Qmesh
  type(sigma_parameters),intent(in) :: Sigp
  type(sigma_results),intent(in) :: Sr
  type(vcoul_t),intent(in) :: Vcp
  type(wfs_descriptor),intent(inout) :: Wfs
  real(dp),intent(in) :: fact_sp
  complex(dpc),intent(in) :: ph_mkit
  complex(dpc),intent(in) :: ph_mkjt
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(in) :: botsq(Sigp%npwc,PPm%dm2_botsq)
  logical,intent(in) :: can_symmetrize(Wfs%nsppol)
  integer,intent(in) :: extrapolar_distrb(ib1:ib2,ib1:ib2,Kmesh%nbz,Wfs%nsppol)
  integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
  integer,intent(in) :: igfftcg0(Sigp%npwc)
  complex(gwpc),intent(in) :: otq(Sigp%npwc,PPm%dm2_otq)
  complex(dpc),intent(inout) :: sigc(2,nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nsppol*Sigp%nsig_ab)
  complex(dpc),intent(inout) :: sigcme_tmp(nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nsppol*Sigp%nsig_ab)
  real(dp),intent(in) :: spinrot_ki(4)
  real(dp),intent(in) :: spinrot_kj(4)
  integer,intent(in) :: tabr_ki(nfftot_gw)
  integer,intent(in) :: tabr_kj(nfftot_gw)
 end subroutine gw_eet_sigma
end interface

interface
 subroutine gw_eet_chi0(Ep,Dtset,Cryst,Wfs,Kmesh,Gsph_epsG0,Gsph_wfn,Psps,Ltg_q,nbvw,qpoint,&  
  &  nfftot_gw,ngfft_gw,use_padfft,igfftepsG0,gw_gbound,gw_mgfft,is,&  
  &  ik_bz,ik_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k,&  
  &  ikmq_ibz,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,&  
  &  qp_energy,chi0,spin_fact,qp_occ,nspinor,tim_fourdp,bbp_ks_distrb,nbmax)
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use defs_abitypes
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: dim_rtwg
  integer,intent(in) :: gw_mgfft
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: isym_k
  integer,intent(in) :: itim_k
  integer,intent(in) :: itim_kmq
  integer,intent(out) :: nbmax
  integer,intent(in) :: nbvw
  integer,intent(in) :: nfftot_gw
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(gvectors_type),intent(in) :: Gsph_wfn
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfs
  complex(dpc),intent(in) :: ph_mkmqt
  complex(dpc),intent(in) :: ph_mkt
  real(dp),intent(in) :: spin_fact
  integer,intent(in) :: ngfft_gw(18)
  integer,intent(in) :: bbp_ks_distrb(nbvw,Kmesh%nbz,Wfs%nsppol)
  complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
  integer,intent(in) :: igfftepsG0(Ep%npwepG0)
  real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: qpoint(3)
  real(dp),intent(in) :: spinrot_k(4)
  real(dp),intent(in) :: spinrot_kmq(4)
  integer,intent(in) :: tabr_k(nfftot_gw)
  integer,intent(in) :: tabr_kmq(nfftot_gw)
 end subroutine gw_eet_chi0
end interface

interface
 subroutine drho_tw_g(paral_kgb,nspinor,npwvec,nr,ngfft,map2sphere,use_padfft,igfftg0,gbound,&  
  &  wfn1,i1,ktabr1,ktabp1,wfn2,dim_rtwg,rhotwg,tim_fourdp,MPI_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: dim_rtwg
  integer,intent(in) :: i1
  integer,intent(in) :: map2sphere
  integer,intent(in) :: npwvec
  integer,intent(in) :: nr
  integer,intent(in) :: nspinor
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  type(mpi_type),intent(inout) :: MPI_enreg
  complex(dpc),intent(in) :: ktabp1
  integer,intent(in) :: gbound(:,:)
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: igfftg0(npwvec*map2sphere)
  integer,intent(in) :: ktabr1(nr)
  complex(gwpc),intent(out) :: rhotwg(npwvec*dim_rtwg)
  complex(gwpc),intent(in) :: wfn1(nr*nspinor)
  complex(gwpc),intent(in) :: wfn2(nr*nspinor)
 end subroutine drho_tw_g
end interface

interface
 subroutine calc_dwfwfg(MPI_enreg,paral_kgb,tim_fourdp,ktabr_k,ktabi_k,nfftot,ngfft_gw,ph_mkt,wfr_jb,wfr_kb,wfg2_jk)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: ktabi_k
  integer,intent(in) :: nfftot
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: MPI_enreg
  complex(dpc),intent(in) :: ph_mkt
  integer,intent(in) :: ngfft_gw(18)
  integer,intent(in) :: ktabr_k(nfftot)
  complex(gwpc),intent(out) :: wfg2_jk(nfftot)
  complex(gwpc),intent(in) :: wfr_jb(nfftot)
  complex(gwpc),intent(in) :: wfr_kb(nfftot)
 end subroutine calc_dwfwfg
end interface

interface
 subroutine calc_ddwfwfg(MPI_enreg,paral_kgb,tim_fourdp,ktabi_k,nfftot,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: ktabi_k
  integer,intent(in) :: nfftot
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: MPI_enreg
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(out) :: wfg2_jk(nfftot)
  complex(gwpc),intent(in) :: wfr_jb(nfftot)
  complex(gwpc),intent(in) :: wfr_kb(nfftot)
 end subroutine calc_ddwfwfg
end interface

interface
 subroutine fft4eet_0(ik_bz,Ep,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax,is,nfftot_gw,ngfft_gw,&  
  &  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,itim_k,tabr_k,ph_mkt,&  
  &  spinrot_k,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,nspinor,tim_fourdp,&  
  &  wfr1,wfr2,ibv,frhorho,spin_fact,qp_occ,qp_energy,chi0,igstart)
  use m_gsphere
  use m_bz_mesh
  use defs_basis
  use m_gwdefs
  use m_wfs
  implicit none
  integer,intent(in) :: dim_rtwg
  integer,intent(in) :: gw_mgfft
  integer,intent(in) :: ibv
  integer,intent(in) :: igstart
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: itim_k
  integer,intent(in) :: itim_kmq
  integer,intent(in) :: nbmax
  integer,intent(in) :: nfftot_gw
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(wfs_descriptor),intent(inout) :: Wfs
  complex(dpc),intent(in) :: ph_mkmqt
  complex(dpc),intent(in) :: ph_mkt
  real(dp),intent(in) :: spin_fact
  integer,intent(in) :: nbhomo(2)
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
  integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
  integer,intent(in) :: igfftepsG0(Ep%npwepG0)
  real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: spinrot_k(4)
  real(dp),intent(in) :: spinrot_kmq(4)
  integer,intent(in) :: tabr_k(nfftot_gw)
  integer,intent(in) :: tabr_kmq(nfftot_gw)
  complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
  complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)
 end subroutine fft4eet_0
end interface

interface
 subroutine fft4eet(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax,&  
  &  is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0,&  
  &  gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k,&  
  &  itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,grottbm1,&  
  &  nspinor,tim_fourdp,wfr1,wfr2,ibv,qplg,kplqg,&  
  &  niter,frhorho,frhoj,fjj,spin_fact,qp_occ,qp_energy,chi0,igstart)
  use m_wfs
  use m_bz_mesh
  use m_crystal
  use defs_basis
  use m_gsphere
  use m_gwdefs
  implicit none
  integer,intent(in) :: dim_rtwg
  integer,intent(in) :: gw_mgfft
  integer,intent(in) :: ibv
  integer,intent(in) :: igstart
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: isym_k
  integer,intent(in) :: itim_k
  integer,intent(in) :: itim_kmq
  integer,intent(in) :: nbmax
  integer,intent(in) :: nfftot_gw
  integer,intent(in) :: niter
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(wfs_descriptor),intent(inout) :: Wfs
  complex(dpc),intent(in) :: ph_mkmqt
  complex(dpc),intent(in) :: ph_mkt
  real(dp),intent(in) :: spin_fact
  integer,intent(in) :: nbhomo(2)
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(gwpc),intent(out) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))
  complex(gwpc),intent(out) :: frhoj(Ep%npwe,Ep%npwe)
  complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
  integer,intent(in) :: grottbm1(Ep%npwvec,2,Cryst%nsym)
  integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
  integer,intent(in) :: igfftepsG0(Ep%npwepG0)
  real(dp),intent(in) :: kplqg(Ep%npwe)
  real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: qplg(Ep%npwe,3)
  real(dp),intent(in) :: spinrot_k(4)
  real(dp),intent(in) :: spinrot_kmq(4)
  integer,intent(in) :: tabr_k(nfftot_gw)
  integer,intent(in) :: tabr_kmq(nfftot_gw)
  complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
  complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)
 end subroutine fft4eet
end interface

interface
 subroutine fft4eet_kb(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,Psps,nbhomo,nbmax,&  
  &  is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0,&  
  &  gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k,&  
  &  itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,grottbm1,&  
  &  nspinor,tim_fourdp,fnlloc,fnlmax,fnlkpr,mtwk,mtwkp,wfr1,wfr2,ibv,qplg,kplqg,&  
  &  niter,frhorho,frhoj,fjj,spin_fact,qp_occ,qp_energy,chi0,igstart)
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: dim_rtwg
  integer,intent(in) :: gw_mgfft
  integer,intent(in) :: ibv
  integer,intent(in) :: igstart
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: isym_k
  integer,intent(in) :: itim_k
  integer,intent(in) :: itim_kmq
  integer,intent(in) :: nbmax
  integer,intent(in) :: nfftot_gw
  integer,intent(in) :: niter
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfs
  complex(dpc),intent(in) :: ph_mkmqt
  complex(dpc),intent(in) :: ph_mkt
  real(dp),intent(in) :: spin_fact
  integer,intent(in) :: nbhomo(2)
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(gwpc),intent(out) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))
  complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  integer,intent(in) :: fnlloc(Cryst%ntypat,2)
  integer,intent(in) :: fnlmax(Cryst%ntypat)
  complex(gwpc),intent(out) :: frhoj(Ep%npwe,Ep%npwe)
  complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
  integer,intent(in) :: grottbm1(Ep%npwvec,2,Cryst%nsym)
  integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
  integer,intent(in) :: igfftepsG0(Ep%npwepG0)
  real(dp),intent(in) :: kplqg(Ep%npwe)
  complex(gwpc),intent(in) :: mtwk(Wfs%nfftot*nspinor,nbhomo(1))
  complex(gwpc),intent(in) :: mtwkp(Wfs%nfftot*nspinor,nbmax)
  real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: qplg(Ep%npwe,3)
  real(dp),intent(in) :: spinrot_k(4)
  real(dp),intent(in) :: spinrot_kmq(4)
  integer,intent(in) :: tabr_k(nfftot_gw)
  integer,intent(in) :: tabr_kmq(nfftot_gw)
  complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
  complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)
 end subroutine fft4eet_kb
end interface

interface
 subroutine calc_eet_prep(Ep,Cryst,Wfs,Kmesh,Psps,is,nbhomo,nbmax,ik_ibz,ikmq_ibz,nspinor,&  
  &  fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)
  use m_wfs
  use m_bz_mesh
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: nbhomo
  integer,intent(in) :: nbmax
  integer,intent(in) :: nspinor
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_parameters),intent(in) :: Ep
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(in) :: Wfs
  complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  integer,intent(in) :: fnlloc(Cryst%ntypat,2)
  integer,intent(in) :: fnlmax(Cryst%ntypat)
  complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbhomo)
  complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,nbmax)
 end subroutine calc_eet_prep
end interface

interface
 subroutine fft4eet_sig(Sigp,Dtset,Cryst,Wfs,Gsph_c,Sr,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw,&  
  &  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k,&  
  &  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,grottbm1,&  
  &  nspinor,tim_fourdp,MPI_enreg,wfr1,vc_sqrt_qbz,i_sz,kb,qplg,kplqg,niter,&  
  &  ptwsq,ik_bz,ikmq_bz,npwc1,npwc2,botsq,otq,sigmac)
  use defs_basis
  use m_wfs
  use m_sigma_results
  use defs_abitypes
  use m_crystal
  use m_gsphere
  use m_gwdefs
  implicit none
  integer,intent(in) :: gw_mgfft
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_bz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: isym_kmq
  integer,intent(in) :: itim_k
  integer,intent(in) :: itim_kmq
  integer,intent(in) :: kb
  integer,intent(in) :: nbhomo
  integer,intent(in) :: nbmax
  integer,intent(in) :: nfftot_gw
  integer,intent(in) :: niter
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(gvectors_type),intent(in) :: Gsph_c
  type(mpi_type),intent(inout) :: MPI_enreg
  type(sigma_parameters),intent(in) :: Sigp
  type(sigma_results),intent(in) :: Sr
  type(wfs_descriptor),intent(inout) :: Wfs
  real(dp),intent(in) :: i_sz
  complex(dpc),intent(in) :: ph_mkmqt
  complex(dpc),intent(in) :: ph_mkt
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)
  integer,intent(in) :: grottbm1(Sigp%npwvec,2,Cryst%nsym)
  integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
  integer,intent(in) :: igfftepsG0(Sigp%npwc)
  real(dp),intent(in) :: kplqg(Sigp%npwc)
  complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
  complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
  real(dp),intent(in) :: qplg(Sigp%npwc,3)
  complex(dpc),intent(inout) :: sigmac(nomega)
  real(dp),intent(in) :: spinrot_k(4)
  real(dp),intent(in) :: spinrot_kmq(4)
  integer,intent(in) :: tabr_k(nfftot_gw)
  integer,intent(in) :: tabr_kmq(nfftot_gw)
  complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
  complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 end subroutine fft4eet_sig
end interface

interface
 subroutine fft4eet_sig_sc(Sigp,Dtset,Cryst,Wfs,Gsph_c,Sr,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw,&  
  &  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k,&  
  &  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,grottbm1,&  
  &  nspinor,tim_fourdp,MPI_enreg,wfr,vc_sqrt_qbz,i_sz,jb,kb,qplg,kplqg,niter,&  
  &  ptwsq,ik_bz,ikmq_bz,npwc1,npwc2,botsq,otq,sigmac)
  use defs_basis
  use m_wfs
  use m_sigma_results
  use defs_abitypes
  use m_crystal
  use m_gsphere
  use m_gwdefs
  implicit none
  integer,intent(in) :: gw_mgfft
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_bz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: isym_kmq
  integer,intent(in) :: itim_k
  integer,intent(in) :: itim_kmq
  integer,intent(in) :: jb
  integer,intent(in) :: kb
  integer,intent(in) :: nbhomo
  integer,intent(in) :: nbmax
  integer,intent(in) :: nfftot_gw
  integer,intent(in) :: niter
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(gvectors_type),intent(in) :: Gsph_c
  type(mpi_type),intent(inout) :: MPI_enreg
  type(sigma_parameters),intent(in) :: Sigp
  type(sigma_results),intent(in) :: Sr
  type(wfs_descriptor),intent(inout) :: Wfs
  real(dp),intent(in) :: i_sz
  complex(dpc),intent(in) :: ph_mkmqt
  complex(dpc),intent(in) :: ph_mkt
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)
  integer,intent(in) :: grottbm1(Sigp%npwvec,2,Cryst%nsym)
  integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
  integer,intent(in) :: igfftepsG0(Sigp%npwc)
  real(dp),intent(in) :: kplqg(Sigp%npwc)
  complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
  complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
  real(dp),intent(in) :: qplg(Sigp%npwc,3)
  complex(dpc),intent(inout) :: sigmac(nomega)
  real(dp),intent(in) :: spinrot_k(4)
  real(dp),intent(in) :: spinrot_kmq(4)
  integer,intent(in) :: tabr_k(nfftot_gw)
  integer,intent(in) :: tabr_kmq(nfftot_gw)
  complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
  complex(gwpc),intent(in) :: wfr(Wfs%nfftot*nspinor,nbmax)
 end subroutine fft4eet_sig_sc
end interface

interface
 subroutine fft4eet_sig_kb(Sigp,Dtset,Cryst,Wfs,Kmesh,Gsph_c,Psps,Sr,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw,&  
  &  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k,&  
  &  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,grottbm1,&  
  &  nspinor,tim_fourdp,MPI_enreg,fnlloc,fnlmax,fnlkr,mtwk,mtwkp,wfr1,&  
  &  vc_sqrt_qbz,i_sz,kb,qplg,kplqg,niter,ptwsq,ik_bz,ikmq_bz,&  
  &  npwc1,npwc2,botsq,otq,sigmac)
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use m_sigma_results
  use defs_abitypes
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: gw_mgfft
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_bz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: isym_kmq
  integer,intent(in) :: itim_k
  integer,intent(in) :: itim_kmq
  integer,intent(in) :: kb
  integer,intent(in) :: nbhomo
  integer,intent(in) :: nbmax
  integer,intent(in) :: nfftot_gw
  integer,intent(in) :: niter
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: use_padfft
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(gvectors_type),intent(in) :: Gsph_c
  type(bz_mesh_type),intent(in) :: Kmesh
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pseudopotential_type),intent(in) :: Psps
  type(sigma_parameters),intent(in) :: Sigp
  type(sigma_results),intent(in) :: Sr
  type(wfs_descriptor),intent(inout) :: Wfs
  real(dp),intent(in) :: i_sz
  complex(dpc),intent(in) :: ph_mkmqt
  complex(dpc),intent(in) :: ph_mkt
  integer,intent(in) :: ngfft_gw(18)
  complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)
  complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  integer,intent(in) :: fnlloc(Cryst%ntypat,2)
  integer,intent(in) :: fnlmax(Cryst%ntypat)
  integer,intent(in) :: grottbm1(Sigp%npwvec,2,Cryst%nsym)
  integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
  integer,intent(in) :: igfftepsG0(Sigp%npwc)
  real(dp),intent(in) :: kplqg(Sigp%npwc)
  complex(gwpc),intent(in) :: mtwk(Wfs%nfftot*nspinor,nbmax)
  complex(gwpc),intent(in) :: mtwkp(Wfs%nfftot*nspinor)
  complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
  complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
  real(dp),intent(in) :: qplg(Sigp%npwc,3)
  complex(dpc),intent(inout) :: sigmac(nomega)
  real(dp),intent(in) :: spinrot_k(4)
  real(dp),intent(in) :: spinrot_kmq(4)
  integer,intent(in) :: tabr_k(nfftot_gw)
  integer,intent(in) :: tabr_kmq(nfftot_gw)
  complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
  complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 end subroutine fft4eet_sig_kb
end interface

interface
 subroutine calc_eet_sig_prep(Sigp,Cryst,Wfs,Kmesh,Psps,is,nbmax,ib1,ib2,ik_ibz,&  
  &  jk_ibz,nspinor,fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)
  use m_wfs
  use m_bz_mesh
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: ib1
  integer,intent(in) :: ib2
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: is
  integer,intent(in) :: jk_ibz
  integer,intent(in) :: nbmax
  integer,intent(in) :: nspinor
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(sigma_parameters),intent(in) :: Sigp
  type(wfs_descriptor),intent(in) :: Wfs
  complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  integer,intent(in) :: fnlloc(Cryst%ntypat,2)
  integer,intent(in) :: fnlmax(Cryst%ntypat)
  complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbmax)
  complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,ib1:ib2)
 end subroutine calc_eet_sig_prep
end interface

interface
 subroutine check_delta_sigma(qpgsq,qpgpsq,delta,omegame0k,omegame0lumo,ig,igp,gw_eet_scale,niter)
  use defs_basis
  implicit none
  integer,intent(in) :: ig
  integer,intent(in) :: igp
  integer,intent(in) :: niter
  complex(gwpc),intent(inout) :: delta
  real(dp),intent(in) :: gw_eet_scale
  real(dp),intent(in) :: omegame0k
  real(dp),intent(in) :: omegame0lumo
  real(dp),intent(in) :: qpgpsq
  real(dp),intent(in) :: qpgsq
 end subroutine check_delta_sigma
end interface

interface
 subroutine gw_eet_sigma_vkb(Sigp,Cryst,Wfs,Kmesh,Psps,isppol,ik_ibz,jk_ibz,ib1,ib2,nspinor,&  
  &  tim_fourdp,MPI_enreg,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)
  use m_wfs
  use m_bz_mesh
  use defs_abitypes
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: ib1
  integer,intent(in) :: ib2
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: isppol
  integer,intent(in) :: jk_ibz
  integer,intent(in) :: nbmax
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pseudopotential_type),intent(in) :: Psps
  type(sigma_parameters),intent(in) :: Sigp
  type(wfs_descriptor),intent(inout) :: Wfs
  complex(gwpc),intent(out) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  complex(gwpc),intent(out) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  integer,intent(out) :: fnlloc(Cryst%ntypat,2)
  integer,intent(out) :: fnlmax(Cryst%ntypat)
  complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbmax)
  complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,ib1:ib2)
 end subroutine gw_eet_sigma_vkb
end interface

interface
 subroutine gw_eet_chi0_vkb(Ep,Cryst,Wfs,Kmesh,Psps,is,ik_ibz,ikmq_ibz,nspinor,tim_fourdp,&  
  &  nbhomo,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)
  use m_wfs
  use m_bz_mesh
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: nbmax
  integer,intent(in) :: nspinor
  integer,intent(in) :: tim_fourdp
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_parameters),intent(in) :: Ep
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfs
  integer,intent(in) :: nbhomo(2)
  complex(gwpc),intent(out) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  complex(gwpc),intent(out) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
  integer,intent(out) :: fnlloc(Cryst%ntypat,2)
  integer,intent(out) :: fnlmax(Cryst%ntypat)
  complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbhomo(1))
  complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,nbmax)
 end subroutine gw_eet_chi0_vkb
end interface

interface
 subroutine calc_sigc_cd(npwc,npwx,nspinor,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&  
  &  omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket,npoles_missing,calc_poles)
  use defs_basis
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegae
  integer,intent(in) :: nomegaei
  integer,intent(in) :: nomegaer
  integer,intent(inout) :: npoles_missing
  integer,intent(in) :: npwc
  integer,intent(in) :: npwx
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: theta_mu_minus_e0i
  logical, intent(in), optional :: calc_poles(nomega)
  complex(gwp) :: epsm1q(npwc,npwc,nomegae)
  complex(gwpc),intent(inout) :: ket(nspinor*npwc,nomega)
  complex(dpc),intent(in) :: omega(nomegae)
  real(dp),intent(in) :: omegame0i(nomega)
  complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 end subroutine calc_sigc_cd
end interface

interface
 subroutine calc_sigc_me(ikcalc,nomega_sigc,minbnd,maxbnd,Dtset,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_Max,Gsph_c,Vcp,Kmesh,Qmesh,&  
  &  Ltg_k,PPm,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,&  
  &  gwc_ngfft,rho_ngfft,rho_nfftot,rhor,use_aerhor,aepaw_rhor,sigcme_tmp)
  use m_vcoul
  use defs_basis
  use m_paw_toolbox
  use m_bz_mesh
  use m_bands_sym
  use m_wfs
  use m_sigma_results
  use defs_abitypes
  use m_ppmodel
  use m_crystal
  use m_gsphere
  use m_paw_pwij
  use defs_datatypes
  use m_gwdefs
  use m_screening
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: maxbnd
  integer,intent(in) :: minbnd
  integer,intent(in) :: nomega_sigc
  integer,intent(in) :: rho_nfftot
  integer,intent(in) :: use_aerhor
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_results),intent(inout) :: Er
  type(gvectors_type),intent(in) :: Gsph_Max
  type(gvectors_type),intent(in) :: Gsph_c
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_k
  type(ppmodel_type),intent(inout) :: PPm
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(bz_mesh_type),intent(in) :: Qmesh
  type(sigma_parameters),intent(in) :: Sigp
  type(sigma_results),intent(in) :: Sr
  type(vcoul_t),intent(in) :: Vcp
  type(wfs_descriptor),intent(inout) :: Wfd
  type(wfs_descriptor),intent(inout) :: Wfdf
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: gwc_ngfft(18)
  integer,intent(in) :: rho_ngfft(18)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  type(bands_symmetries),intent(in) :: QP_sym(Wfd%nsppol)
  real(dp),intent(in) :: aepaw_rhor(rho_nfftot,Wfd%nspden*use_aerhor)
  real(dp),intent(in) :: rhor(rho_nfftot,Wfd%nspden)
  complex(dpc),intent(out) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd, &
  &         minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 end subroutine calc_sigc_me
end interface

interface
 subroutine calc_sigc_pole_cd(npwc,npwx,nspinor,ncoeff,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&  
  &  omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket,npoles_missing,calc_poles)
  use defs_basis
  implicit none
  integer,intent(in) :: ncoeff
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegae
  integer,intent(in) :: nomegaei
  integer,intent(in) :: nomegaer
  integer,intent(inout) :: npoles_missing
  integer,intent(in) :: npwc
  integer,intent(in) :: npwx
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: theta_mu_minus_e0i
  logical, intent(in), optional :: calc_poles(nomega)
  real(gwp) :: epsm1q(npwc,npwc,ncoeff)
  complex(gwpc),intent(inout) :: ket(nspinor*npwc,nomega)
  complex(dpc),intent(in) :: omega(nomegae)
  real(dp),intent(in) :: omegame0i(nomega)
  complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 end subroutine calc_sigc_pole_cd
end interface

interface
 subroutine calc_sigx_me(ikcalc,minbnd,maxbnd,Cryst,QP_BSt,Sigp,Gsph_x,Vcp,Kmesh,Qmesh,&  
  &  Ltg_k,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,gwx_ngfft,ngfftf,&  
  &  prtvol,pawcross,sigxme_tmp)
  use m_vcoul
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use m_paw_toolbox
  use m_bands_sym
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  use m_paw_pwij
  implicit none
  integer,intent(in) :: ikcalc
  integer,intent(in) :: maxbnd
  integer,intent(in) :: minbnd
  integer,intent(in) :: pawcross
  integer,intent(in) :: prtvol
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(in) :: Gsph_x
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_k
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(bz_mesh_type),intent(in) :: Qmesh
  type(sigma_parameters),intent(in) :: Sigp
  type(vcoul_t),intent(in) :: Vcp
  type(wfs_descriptor),intent(inout) :: Wfd
  type(wfs_descriptor),intent(inout) :: Wfdf
  integer,intent(in) :: gwx_ngfft(18)
  integer,intent(in) :: ngfftf(18)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  type(bands_symmetries),intent(in) :: QP_sym(Wfd%nsppol)
  complex(dpc),intent(out) :: sigxme_tmp(minbnd:maxbnd,minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 end subroutine calc_sigx_me
end interface

interface
 subroutine calc_vhxc_me(Wfd,Mflags,Mels,Cryst,Dtset,gsqcutf_eff,nfftf,ngfftf,&  
  &  vtrial,vhartr,vxc,Psps,Pawtab,Paw_an,Pawang,Pawfgrtab,Paw_ij,dijexc_core,&  
  &  rhor,rhog,usexcnhat,nhat,nhatgr,nhatgrdim,kstab,&  
  &  taug,taur) ! optional arguments
  use m_melemts
  use m_wfs
  use defs_abitypes
  use m_crystal
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: usexcnhat
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(melements_type),intent(out) :: Mels
  type(melements_flags_type),intent(in) :: Mflags
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfd
  real(dp),intent(in) :: gsqcutf_eff
  integer,intent(in) :: ngfftf(18)
  type(paw_an_type),intent(in) :: Paw_an(Cryst%natom)
  type(paw_ij_type),intent(inout) :: Paw_ij(Cryst%natom)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
  real(dp),intent(in) :: dijexc_core(:,:,:)
  integer,intent(in) :: kstab(2,Wfd%nkibz,Wfd%nsppol)
  real(dp),intent(in) :: nhat(nfftf,Wfd%nspden*Wfd%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,Wfd%nspden,3*nhatgrdim)
  real(dp),intent(in) :: rhog(2,nfftf)
  real(dp),intent(in) :: rhor(nfftf,Wfd%nspden)
  real(dp),intent(in),optional :: taug(2,nfftf*Dtset%usekden)
  real(dp),intent(in),optional :: taur(nfftf,Wfd%nspden*Dtset%usekden)
  real(dp),intent(in) :: vhartr(nfftf)
  real(dp),intent(in) :: vtrial(nfftf,Wfd%nspden)
  real(dp),intent(in) :: vxc(nfftf,Wfd%nspden)
 end subroutine calc_vhxc_me
end interface

interface
 subroutine update_cprj(natom,nkibz,nbnds,nsppol,nspinor,m_lda_to_qp,dimlmn,Cprj_ibz)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nbnds
  integer,intent(in) :: nkibz
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(cprj_type),intent(inout) :: Cprj_ibz(natom,nspinor*nbnds*nkibz*nsppol)
  integer,intent(in) :: dimlmn(natom)
  complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 end subroutine update_cprj
end interface

interface
 subroutine cchi0(use_tr,Dtset,Cryst,qpoint,Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,Gsph_wfn,&  
  &  Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,nbvw,ngfft_gw,nfftot_gw,ngfftf,nfftf_tot,&  
  &  chi0,ktabr,ktabrf,Ltg_q,chi0_sumrule,Wfd,Wfdf)
  use m_wfs
  use m_bz_mesh
  use m_paw_toolbox
  use defs_abitypes
  use m_gsphere
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nbvw
  integer,intent(in) :: nfftf_tot
  integer,intent(in) :: nfftot_gw
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(gvectors_type),intent(in) :: Gsph_wfn
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(wfs_descriptor),intent(inout) :: Wfd
  type(wfs_descriptor),intent(inout) :: Wfdf
  logical,intent(in) :: use_tr
  integer,intent(in) :: ngfft_gw(18)
  integer,intent(in) :: ngfftf(18)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  real(dp),intent(out) :: chi0_sumrule(Ep%npwe)
  integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz)
  integer,intent(in) :: ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
  real(dp),intent(in) :: qpoint(3)
 end subroutine cchi0
end interface

interface
 subroutine calc_delta_chi0(Ep,frhorho,frhoj,fjj,niter)
  use defs_basis
  use m_gwdefs
  implicit none
  integer,intent(in) :: niter
  type(epsilonm1_parameters),intent(in) :: Ep
  complex(gwpc),intent(inout) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))
  complex(gwpc),intent(inout) :: frhoj(Ep%npwe,Ep%npwe)
  complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 end subroutine calc_delta_chi0
end interface

interface
 subroutine calc_chi0_delta_clos(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,niter,epsv,epslumo,qpgsq,frhorho,frhoj,fjj,chi0)
  use m_gsphere
  use defs_abitypes
  use defs_basis
  use m_gwdefs
  use m_bz_mesh
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: niter
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: epslumo
  real(dp),intent(in) :: epsv
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(gwpc),intent(in) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))
  complex(gwpc),intent(in) :: frhoj(Ep%npwe,Ep%npwe)
  complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
  real(dp),intent(in) :: qpgsq(Ep%npwe)
 end subroutine calc_chi0_delta_clos
end interface

interface
 subroutine calc_chi0_delta0(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,epsv,epslumo,qpgsq,frhorho,chi0)
  use m_gsphere
  use defs_abitypes
  use defs_basis
  use m_gwdefs
  use m_bz_mesh
  implicit none
  integer,intent(in) :: ik_bz
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: epslumo
  real(dp),intent(in) :: epsv
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
  real(dp),intent(in) :: qpgsq(Ep%npwe)
 end subroutine calc_chi0_delta0
end interface

interface
 subroutine calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg,spin_fact,qp_occ,qp_energy,&  
  &  ibv,ibc,ik_ibz,ikmq_ibz,is,chi0)
  use m_bz_mesh
  use defs_basis
  use m_gwdefs
  use m_gsphere
  implicit none
  integer,intent(in) :: ibc
  integer,intent(in) :: ibv
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: spin_fact
  complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
  complex(gwpc),intent(in) :: rhotwg(Ep%npwepG0)
 end subroutine calc_corr_chi0
end interface

interface
 subroutine calc_corr_sig(Sigp,Sr,nomega,nspinor,npwc1,npwc2,botsq,otq,rhotwg,is,ibv,kb,&  
  &  ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,i_sz,vc_sqrt_qbz,sigmac)
  use m_sigma_results
  use m_gwdefs
  use defs_basis
  implicit none
  integer,intent(in) :: ibv
  integer,intent(in) :: ik_bz
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_bz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: is
  integer,intent(in) :: kb
  integer,intent(in) :: nomega
  integer,intent(in) :: npwc1
  integer,intent(in) :: npwc2
  integer,intent(in) :: nspinor
  type(sigma_parameters),intent(in) :: Sigp
  type(sigma_results),intent(in) :: Sr
  real(dp),intent(in) :: i_sz
  complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)
  complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
  complex(gwpc),intent(in) :: rhotwg(Sigp%npwc)
  complex(dpc),intent(inout) :: sigmac(nomega)
  complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 end subroutine calc_corr_sig
end interface

interface
 subroutine calc_delta0(Dtset,Ep,qpgsq,delta)
  use defs_basis
  use defs_abitypes
  use m_gwdefs
  implicit none
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  complex(gwpc),intent(out) :: delta(Ep%npwe,Ep%npwe,Ep%nomega)
  real(dp), intent(in) :: qpgsq(Ep%npwe)
 end subroutine calc_delta0
end interface

interface
 subroutine calc_chi0_delta0_bis(ik_bz,Ep,Gsph_epsG0,Ltg_q,paux,delta,chi0)
  use m_gsphere
  use defs_basis
  use m_gwdefs
  use m_bz_mesh
  implicit none
  integer,intent(in) :: ik_bz
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  complex(gwpc),intent(out) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(gwpc),intent(in) :: delta(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(gwpc),intent(in) :: paux(Ep%npwe*(Ep%npwe+1)/2)
 end subroutine calc_chi0_delta0_bis
end interface

interface
 subroutine check_delta(Ep,qpgsq,delta,epsv,epslumo,ig,igp)
  use defs_basis
  use m_gwdefs
  implicit none
  integer,intent(in) :: ig
  integer,intent(in) :: igp
  type(epsilonm1_parameters),intent(in) :: Ep
  complex(gwpc),intent(inout) :: delta
  real(dp),intent(in) :: epslumo
  real(dp),intent(in) :: epsv
  real(dp),intent(in) :: qpgsq(Ep%npwe)
 end subroutine check_delta
end interface

interface
 subroutine cchi0q0(use_tr,Dtset,Cryst,Ep,Psps,Kmesh,QP_BSt,KS_BSt,Gsph_epsG0,Gsph_wfn,&  
  &  Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,Pawfgrtab,Paw_onsite,ktabr,ktabrf,nbvw,ngfft_gw,&  
  &  nfftot_gw,ngfftf,nfftf_tot,chi0,chi0_head,chi0_lwing,chi0_uwing,Ltg_q,chi0_sumrule,Wfd,Wfdf)
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use m_paw_toolbox
  use defs_abitypes
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nbvw
  integer,intent(in) :: nfftf_tot
  integer,intent(in) :: nfftot_gw
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(gvectors_type),intent(in) :: Gsph_wfn
  type(bandstructure_type),intent(in) :: KS_BSt
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(wfs_descriptor),intent(inout) :: Wfd
  type(wfs_descriptor),intent(inout) :: Wfdf
  logical,intent(in) :: use_tr
  integer,intent(in) :: ngfft_gw(18)
  integer,intent(in) :: ngfftf(18)
  type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Psps%usepaw)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  complex(gwpc),intent(out) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
  complex(dpc),intent(out) :: chi0_head(3,3,Ep%nomega)
  complex(dpc),intent(out) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
  real(dp),intent(out) :: chi0_sumrule(Ep%npwe)
  complex(dpc),intent(out) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)
  integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz)
  integer,intent(in) :: ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
 end subroutine cchi0q0
end interface

interface
 subroutine chi0q0_intraband(Wfd,Cryst,Ep,Psps,BSt,Gsph_epsG0,Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,use_tr,usepawu,&  
  &  ngfft_gw,chi0,chi0_head,chi0_lwing,chi0_uwing)
  use defs_basis
  use m_wfs
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  use m_paw_pwij
  implicit none
  integer,intent(in) :: usepawu
  type(bandstructure_type),intent(in) :: BSt
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfd
  logical,intent(in) :: use_tr
  integer,intent(in) :: ngfft_gw(18)
  type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Psps%usepaw)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  complex(dpc),intent(out) :: chi0_head(3,3,Ep%nomega)
  complex(dpc),intent(out) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
  complex(dpc),intent(out) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)
 end subroutine chi0q0_intraband
end interface

interface
 subroutine check_completeness(use_tr,Dtset,Cryst,qpoint,Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,&  
  &  Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,ngfft_gw,nfftot_gw,ngfftf,nfftf_tot,deltaI,ktabr,ktabrf,Ltg_q,Wfd,Wfdf)
  use m_wfs
  use m_bz_mesh
  use m_paw_toolbox
  use defs_abitypes
  use m_gsphere
  use m_crystal
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nfftf_tot
  integer,intent(in) :: nfftot_gw
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(wfs_descriptor),intent(inout) :: Wfd
  type(wfs_descriptor),intent(inout) :: Wfdf
  logical,intent(in) :: use_tr
  integer,intent(in) :: ngfft_gw(18)
  integer,intent(in) :: ngfftf(18)
  type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  complex(gwpc),intent(out) :: deltaI(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,1)
  integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz)
  integer,intent(in) :: ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
  real(dp),intent(in) :: qpoint(3)
 end subroutine check_completeness
end interface

interface
 subroutine cohsex_me(ikcalc,nomega_sigc,minbnd,maxbnd,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_c,Vcp,Kmesh,Qmesh,&  
  &  Ltg_k,Pawtab,Pawang,Paw_pwff,Psps,Wfd,QP_sym,gwc_ngfft,accesswff,prtvol,sigcme_tmp)
  use m_vcoul
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use m_sigma_results
  use m_bands_sym
  use m_crystal
  use m_gsphere
  use m_paw_pwij
  use defs_datatypes
  use m_gwdefs
  use m_screening
  implicit none
  integer,intent(in) :: accesswff
  integer,intent(in) :: ikcalc
  integer,intent(in) :: maxbnd
  integer,intent(in) :: minbnd
  integer,intent(in) :: nomega_sigc
  integer,intent(in) :: prtvol
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_results),intent(inout) :: Er
  type(gvectors_type),intent(in) :: Gsph_c
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_k
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(bz_mesh_type),intent(in) :: Qmesh
  type(sigma_parameters),intent(in) :: Sigp
  type(sigma_results),intent(in) :: Sr
  type(vcoul_t),intent(in) :: Vcp
  type(wfs_descriptor),intent(inout) :: Wfd
  type(pawang_type),intent(in) :: pawang
  integer,intent(in) :: gwc_ngfft(18)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  type(bands_symmetries),intent(in) :: QP_sym(Wfd%nsppol)
  complex(dpc),intent(out) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd, &
  &         minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 end subroutine cohsex_me
end interface

interface
 subroutine cutoff_m_elem(ep,kmesh,gvec,Wf,energy,z0,wdth,occ,direction,gprimd)
  use m_bz_mesh
  use defs_basis
  use m_gwdefs
  use m_wfs
  implicit none
  integer,intent(in) :: direction
  type(wfs_descriptor),optional,intent(in) :: Wf
  type(epsilonm1_parameters),intent(in) :: ep
  type(bz_mesh_type),target,intent(in) :: kmesh
  real(dp),intent(in) :: wdth
  real(dp),intent(in) :: z0
  real(dp),intent(in) :: energy(ep%nbnds,ep%nkibz,ep%nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,ep%npwvec)
  real(dp),intent(in) :: occ(ep%nbnds,ep%nkibz,ep%nsppol)
 end subroutine cutoff_m_elem
end interface

interface
 subroutine matrixelements(npwwfn,wfg1,wfg2,gvec,kpoint,res)
  use defs_basis
  implicit none
  integer,intent(in) :: npwwfn
  integer,intent(in) :: gvec(3,npwwfn)
  real(dp),intent(in) :: kpoint(3)
  complex(gwpc),intent(out) :: res(3)
  complex(gwpc),intent(in) :: wfg1(npwwfn)
  complex(gwpc),intent(in) :: wfg2(npwwfn)
 end subroutine matrixelements
end interface

interface
 subroutine matrixelements_cutoff(npwwfn,wfg1,wfg2,gvec,kpoint,z0,wdth,direction,res)
  use defs_basis
  implicit none
  integer,intent(in) :: direction
  integer,intent(in) :: npwwfn
  real(dp), intent(in) :: wdth
  real(dp), intent(in) :: z0
  integer,intent(in) :: gvec(3,npwwfn)
  real(dp),intent(in) :: kpoint(3)
  complex(gwpc),intent(out) :: res(3)
  complex(gwpc),intent(in) :: wfg1(npwwfn)
  complex(gwpc),intent(in) :: wfg2(npwwfn)
 end subroutine matrixelements_cutoff
end interface

interface
 subroutine check_zarot(npwvec,Cryst,ngfft,gvec,psps,pawang,grottb,grottbm1)
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: npwvec
  type(crystal_structure),intent(in) :: Cryst
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: grottb(npwvec,Cryst%timrev,Cryst%nsym)
  integer,intent(in) :: grottbm1(npwvec,Cryst%timrev,Cryst%nsym)
  integer,intent(in) :: gvec(3,npwvec)
 end subroutine check_zarot
end interface

interface
 subroutine paw_check_symcprj(Wfd,ik_bz,band,spin,sym_mode,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_bz)
  use m_bz_mesh
  use m_crystal
  use defs_datatypes
  use m_wfs
  implicit none
  integer,intent(in) :: band
  integer,intent(in) :: ik_bz
  integer,intent(in) :: spin
  integer,intent(in) :: sym_mode
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfd
  type(cprj_type),intent(out) :: Cprj_bz(Cryst%natom,Wfd%nspinor)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 end subroutine paw_check_symcprj
end interface

interface
 function dotproductqrc(r,c,b1,b2,b3)
  use defs_basis
  implicit none
  complex(gwpc) :: dotproductqrc
  real(dp),intent(in) :: b1(3)
  real(dp),intent(in) :: b2(3)
  real(dp),intent(in) :: b3(3)
  complex(gwpc),intent(in) :: c(3)
  real(dp),intent(in) :: r(3)
 end function dotproductqrc
end interface

interface
 function pdtqrc(R,C,b1,b2,b3)
  use defs_basis
  implicit none
  complex(dpc) :: pdtqrc
  complex(dpc),intent(in) :: C(3)
  real(dp),intent(in) :: R(3)
  real(dp),intent(in) :: b1(3)
  real(dp),intent(in) :: b2(3)
  real(dp),intent(in) :: b3(3)
 end function pdtqrc
end interface

interface
 subroutine fsumrule(nomega,omega,eps,omegaplasma,method)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(in) :: nomega
  real(dp),intent(in) :: omegaplasma
  real(dp),intent(in) :: eps(nomega)
  real(dp),intent(in) :: omega(nomega)
 end subroutine fsumrule
end interface

interface
 subroutine make_transitions(Wfd,chi0alg,nbnds,nbvw,nsppol,symchi,timrev,TOL_DELTA_OCC,&  
  &  max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,gw_energy,occ,qpoint,bbp_ks_distrb)
  use defs_basis
  use m_bz_mesh
  use m_wfs
  implicit none
  integer,intent(in) :: chi0alg
  integer,intent(in) :: nbnds
  integer,intent(in) :: nbvw
  integer,intent(in) :: nsppol
  integer,intent(in) :: symchi
  integer,intent(in) :: timrev
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: TOL_DELTA_OCC
  type(wfs_descriptor),intent(in) :: Wfd
  real(dp),intent(out) :: max_rest
  real(dp),intent(out) :: min_rest
  real(dp),intent(out) :: my_max_rest
  real(dp),intent(out) :: my_min_rest
  integer,intent(in) :: bbp_ks_distrb(Wfd%mband,Wfd%mband,Kmesh%nbz,Wfd%nsppol)
  real(dp),intent(in) :: gw_energy(nbnds,Kmesh%nibz,nsppol)
  real(dp),intent(in) :: occ(nbnds,Kmesh%nibz,nsppol)
  real(dp),intent(in) :: qpoint(3)
 end subroutine make_transitions
end interface

interface
 subroutine get_rhor(fname,accesswff,nspden,nfft_asked,ngfft_asked,paral_kgb,MPI_enreg,rhor_out,get_pawden)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: accesswff
  integer,intent(in) :: nfft_asked
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: MPI_enreg
  character(len=fnlen),intent(in) :: fname
  logical,intent(in), optional :: get_pawden
  integer,intent(in) :: ngfft_asked(18)
  real(dp),intent(out) :: rhor_out(nfft_asked,nspden)
 end subroutine get_rhor
end interface

interface
 subroutine sigma_distribution(Wfd,Kmesh,Ltg_kgw,Qmesh,nsppol,can_symmetrize,kptgw,mg0,my_nbks,proc_distrb,got,bks_mask,global)
  use m_bz_mesh
  use defs_basis
  use m_wfs
  implicit none
  integer,intent(out) :: my_nbks
  integer,intent(in) :: nsppol
  type(bz_mesh_type),intent(in) :: Kmesh
  type(little_group),intent(in) :: Ltg_kgw
  type(bz_mesh_type),intent(in) :: Qmesh
  type(wfs_descriptor),intent(inout) :: Wfd
  logical,optional,intent(in) :: global
  integer,intent(in) :: mg0(3)
  logical,optional,intent(in) :: bks_mask(Wfd%mband,Kmesh%nbz,nsppol)
  logical,intent(in) :: can_symmetrize(Wfd%nsppol)
  integer,optional,intent(inout) :: got(Wfd%nproc)
  real(dp),intent(in) :: kptgw(3)
  integer,intent(out) :: proc_distrb(Wfd%mband,Kmesh%nbz,nsppol)
 end subroutine sigma_distribution
end interface

interface
 subroutine chi0_bbp_mask(Ep,use_tr,QP_BSt,mband,ikmq_ibz,ik_ibz,spin,spin_fact,bbp_mask)
  use defs_basis
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: ikmq_ibz
  integer,intent(in) :: mband
  integer,intent(in) :: spin
  type(epsilonm1_parameters),intent(in) :: Ep
  type(bandstructure_type),intent(in) :: QP_BSt
  real(dp),intent(in) :: spin_fact
  logical,intent(in) :: use_tr
  logical,intent(out) :: bbp_mask(mband,mband)
 end subroutine chi0_bbp_mask
end interface

interface
 subroutine completechi0_deltapart(ik_bz,qzero,symchi,npwe,npwvec,nomega,nspinor,&  
  &  nfftot,ngfft,igfft0,Gsph_FFT,Ltg_q,green_enhigh_w,wfwfg,chi0)
  use m_gsphere
  use defs_basis
  use m_bz_mesh
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: nfftot
  integer,intent(in) :: nomega
  integer,intent(in) :: npwe
  integer,intent(in) :: npwvec
  integer,intent(in) :: nspinor
  integer,intent(in) :: symchi
  type(gvectors_type),intent(in) :: Gsph_FFT
  type(little_group),intent(in) :: Ltg_q
  logical,intent(in) :: qzero
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)
  complex(dpc),intent(in) :: green_enhigh_w(nomega)
  integer,intent(in) :: igfft0(npwvec)
  complex(gwpc),intent(in) :: wfwfg(nfftot*nspinor**2)
 end subroutine completechi0_deltapart
end interface

interface
 subroutine output_chi0sumrule(qeq0,iq,npwe,omegaplasma,chi0sumrule,epsm1_w0,vc_sqrt)
  use defs_basis
  implicit none
  integer,intent(in) :: iq
  integer,intent(in) :: npwe
  real(dp),intent(in) :: omegaplasma
  logical,intent(in) :: qeq0
  real(dp),intent(inout) :: chi0sumrule(npwe)
  complex(gwpc),intent(in) :: epsm1_w0(npwe,npwe)
  complex(gwpc),intent(in) :: vc_sqrt(npwe)
 end subroutine output_chi0sumrule
end interface

interface
 subroutine accumulate_chi0sumrule(ik_bz,symchi,npwe,factor,delta_ene,&  
  &  Ltg_q,Gsph_epsG0,npwepG0,rhotwg,chi0sumrule)
  use defs_basis
  use m_gsphere
  use m_bz_mesh
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: npwe
  integer,intent(in) :: npwepG0
  integer,intent(in) :: symchi
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: delta_ene
  real(dp),intent(in) :: factor
  real(dp),intent(inout) :: chi0sumrule(npwe)
  complex(gwpc),intent(in) :: rhotwg(npwepG0)
 end subroutine accumulate_chi0sumrule
end interface

interface
 subroutine mlwfovlp_qp(cg,Cprj_BZ,dtset,dtfil,eigen,mband,mcg,mcprj,mkmem,mpw,natom,&  
  &  nkpt,npwarr,nspden,nsppol,ntypat,Hdr,Pawtab,rprimd,MPI_enreg)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  type(hdr_type),intent(in) :: Hdr
  type(mpi_type),intent(inout) :: MPI_enreg
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(cprj_type),target,intent(inout) :: Cprj_BZ(natom,mcprj)
  type(pawtab_type),intent(in) :: Pawtab(ntypat*Dtset%usepaw)
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine mlwfovlp_qp
end interface

interface
 subroutine calc_coh(nspinor,nsig_ab,nfftot,ngfft,npwc,gvec,wfg2_jk,epsm1q_o,vc_sqrt,i_sz,iqibz,same_band,sigcohme)
  use defs_basis
  implicit none
  integer,intent(in) :: iqibz
  integer,intent(in) :: nfftot
  integer,intent(in) :: npwc
  integer,intent(in) :: nsig_ab
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: i_sz
  logical,intent(in) :: same_band
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(in) :: epsm1q_o(npwc,npwc)
  integer,intent(in) :: gvec(3,npwc)
  complex(gwpc),intent(out) :: sigcohme(nsig_ab)
  complex(gwpc),intent(in) :: vc_sqrt(npwc)
  complex(gwpc),intent(in) :: wfg2_jk(nsig_ab*nfftot)
 end subroutine calc_coh
end interface

interface
 subroutine calc_coh_comp(iqibz,i_sz,same_band,nspinor,nsig_ab,ediff,npwc,gvec,&  
  &  ngfft,nfftot,wfg2_jk,vc_sqrt,botsq,otq,sigcohme)
  use defs_basis
  implicit none
  integer,intent(in) :: iqibz
  integer,intent(in) :: nfftot
  integer,intent(in) :: npwc
  integer,intent(in) :: nsig_ab
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: ediff
  real(dp),intent(in) :: i_sz
  logical,intent(in) :: same_band
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(in) :: botsq(npwc,npwc)
  integer,intent(in) :: gvec(3,npwc)
  complex(gwpc),intent(in) :: otq(npwc,npwc)
  complex(gwpc),intent(out) :: sigcohme(nsig_ab)
  complex(gwpc),intent(in) :: vc_sqrt(npwc)
  complex(gwpc),intent(in) :: wfg2_jk(nsig_ab*nfftot)
 end subroutine calc_coh_comp
end interface

interface
 subroutine paw_qpscgw(Wfd,nscf,nfftf,ngfftf,Dtset,Cryst,Kmesh,Psps,QP_BSt,&  
  &  Pawang,Pawrad,Pawtab,Pawfgrtab,prev_Pawrhoij,&  
  &  QP_pawrhoij,QP_paw_ij,QP_paw_an,QP_energies,qp_nhat,nhatgrdim,qp_nhatgr,qp_compch_sph,qp_compch_fft,MPI_enreg)
  use m_wfs
  use m_bz_mesh
  use m_energies
  use defs_abitypes
  use m_crystal
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nscf
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(bz_mesh_type),intent(in) :: Kmesh
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(energies_type),intent(inout) :: QP_energies
  type(wfs_descriptor),intent(inout) :: Wfd
  real(dp),intent(out) :: qp_compch_fft
  real(dp),intent(out) :: qp_compch_sph
  integer,intent(in) :: ngfftf(18)
  type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)
  type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  type(paw_an_type),intent(inout) :: QP_paw_an(Cryst%natom)
  type(paw_ij_type),intent(out) :: QP_paw_ij(Cryst%natom)
  type(pawrhoij_type),intent(out) :: QP_pawrhoij(Cryst%natom)
  type(pawrhoij_type),intent(inout) :: prev_Pawrhoij(Cryst%natom)
  real(dp),intent(out) :: qp_nhat(nfftf,Dtset%nspden)
  real(dp),intent(out) :: qp_nhatgr(nfftf,Dtset%nspden,3*nhatgrdim)
 end subroutine paw_qpscgw
end interface

interface
 subroutine print_psps(psps,unit,prtvol,mode_paral)
  use defs_datatypes
  implicit none
  integer,intent(in),optional :: prtvol
  integer,intent(in),optional :: unit
  character(len=4),intent(in),optional :: mode_paral
  type(pseudopotential_type),intent(in) :: psps
 end subroutine print_psps
end interface

interface
 subroutine plot_psps(psps,root_filename)
  use defs_basis
  use defs_datatypes
  implicit none
  type(pseudopotential_type),intent(in) :: psps
  character(len=fnlen),intent(in),optional :: root_filename
 end subroutine plot_psps
end interface

interface
 subroutine q0fit(nq,q,gvec,nomega,omega,npwvec,chi0,qcut,metal,&  
  &  nop,op,ninv,gprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: ninv
  integer,intent(in) :: nomega
  integer,intent(in) :: nop
  integer,intent(in) :: npwvec
  integer,intent(in) :: nq
  logical,intent(in) :: metal
  real(dp),intent(in) :: qcut
  complex(gwpc),intent(inout) :: chi0(nq,npwvec,npwvec,nomega)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: gvec(3,npwvec)
  complex(dpc),intent(in) :: omega(nomega)
  real(dp),intent(in) :: op(3,3,nop)
  real(dp),intent(in) :: q(3,nq)
 end subroutine q0fit
end interface

interface
 subroutine setup_screening(codvsn,acell,rprim,ngfftf,ikss_fname,Dtset,Psps,Pawtab,&  
  &  ngfft_gw,Hdr_kss,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Ltg_q,Gsph_epsG0,Gsph_wfn,Vcp,Ep,comm)
  use m_vcoul
  use m_bz_mesh
  use m_crystal
  use defs_abitypes
  use defs_basis
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  implicit none
  integer,intent(in) :: comm
  type(crystal_structure),intent(out) :: Cryst
  type(dataset_type),intent(inout) :: Dtset
  type(epsilonm1_parameters),intent(out) :: Ep
  type(gvectors_type),intent(out) :: Gsph_epsG0
  type(gvectors_type),intent(out) :: Gsph_wfn
  type(hdr_type),intent(out) :: Hdr_kss
  type(hdr_type),intent(out) :: Hdr_out
  type(bandstructure_type),intent(out) :: KS_BSt
  type(bz_mesh_type),intent(out) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(bz_mesh_type),intent(out) :: Qmesh
  type(vcoul_t),intent(out) :: Vcp
  character(len=6),intent(in) :: codvsn
  character(len=fnlen),intent(in) :: ikss_fname
  integer,intent(out) :: ngfft_gw(18)
  integer,intent(in) :: ngfftf(18)
  type(little_group),pointer :: Ltg_q(:)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine setup_screening
end interface

interface
 subroutine setup_sigma(codvsn,acell,rprim,ngfftf,Dtset,Dtfil,Psps,Pawtab,&  
  &  gwx_ngfft,gwc_ngfft,Hdr_kss,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Gsph_Max,Gsph_c,Vcp,Er,Sigp,comm)
  use m_vcoul
  use m_bz_mesh
  use m_crystal
  use defs_abitypes
  use defs_basis
  use m_gsphere
  use defs_datatypes
  use m_gwdefs
  use m_screening
  implicit none
  integer,intent(in) :: comm
  type(crystal_structure),intent(out) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(epsilonm1_results),intent(out) :: Er
  type(gvectors_type),intent(out) :: Gsph_Max
  type(gvectors_type),intent(out) :: Gsph_c
  type(hdr_type),intent(out) :: Hdr_kss
  type(hdr_type),intent(out) :: Hdr_out
  type(bandstructure_type),intent(out) :: KS_BSt
  type(bz_mesh_type),intent(out) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(bz_mesh_type),intent(out) :: Qmesh
  type(sigma_parameters),intent(out) :: Sigp
  type(vcoul_t),intent(out) :: Vcp
  character(len=6),intent(in) :: codvsn
  integer,intent(out) :: gwc_ngfft(18)
  integer,intent(out) :: gwx_ngfft(18)
  integer,intent(in) :: ngfftf(18)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine setup_sigma
end interface

interface
 subroutine sigma_tables(Sigp,Kmesh,Bnd_sym)
  use m_bz_mesh
  use m_bands_sym
  use m_gwdefs
  implicit none
  type(bz_mesh_type),intent(in) :: Kmesh
  type(sigma_parameters),intent(inout) :: Sigp
  type(bands_symmetries),optional,intent(in) :: Bnd_sym(Kmesh%nibz,Sigp%nsppol)
 end subroutine sigma_tables
end interface

interface
 subroutine interpolate_sigmak(Cryst,Kmesh,kptrlatt,nshiftk,shiftk,lbd,ubd,isigk_ij,onkpt,okpt,osigk_ij,ierr)
  use m_bz_mesh
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: lbd
  integer,intent(in) :: nshiftk
  integer,intent(in) :: onkpt
  integer,intent(in) :: ubd
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  integer,intent(in) :: kptrlatt(3,3)
  complex(dpc),intent(in) :: isigk_ij(lbd:ubd,lbd:ubd,Kmesh%nibz)
  real(dp),intent(in) :: okpt(3,onkpt)
  complex(dpc),intent(out) :: osigk_ij(lbd:ubd,lbd:ubd,onkpt)
  real(dp),intent(in) :: shiftk(3,nshiftk)
 end subroutine interpolate_sigmak
end interface

interface
 subroutine approxdelta(nomegasf,omegasf,egwdiff_re,smear,iomegal,iomegar,wl,wr,spmeth)
  use defs_basis
  implicit none
  integer,intent(out) :: iomegal
  integer,intent(out) :: iomegar
  integer,intent(in) :: nomegasf
  integer,intent(in) :: spmeth
  real(dp),intent(in) :: egwdiff_re
  real(dp),intent(in) :: smear
  real(dp),intent(out) :: wl
  real(dp),intent(out) :: wr
  real(dp),intent(in) :: omegasf(nomegasf)
 end subroutine approxdelta
end interface

interface
 subroutine calc_kkweight(ne,omegae,nsp,omegasp,delta,omegamax,kkw)
  use defs_basis
  implicit none
  integer,intent(in) :: ne
  integer,intent(in) :: nsp
  real(dp),intent(in) :: delta
  real(dp),intent(in) :: omegamax
  complex(dpc),intent(out) :: kkw(nsp,ne)
  complex(dpc),intent(in) :: omegae(ne)
  real(dp),intent(in) :: omegasp(nsp)
 end subroutine calc_kkweight
end interface

interface
 subroutine setup_spectral(nomega,omega,nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&  
  &  method,zcut,omegaplasma,my_wl,my_wr,kkweight)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(out) :: my_wl
  integer,intent(out) :: my_wr
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegasf
  real(dp),intent(in) :: max_rest
  real(dp),intent(in) :: min_rest
  real(dp),intent(in) :: my_max_rest
  real(dp),intent(in) :: my_min_rest
  real(dp),intent(in) :: omegaplasma
  real(dp),intent(in) :: zcut
  complex(dpc),intent(out) :: kkweight(nomegasf,nomega)
  complex(dpc),intent(in) :: omega(nomega)
  real(dp),intent(out) :: omegasf(nomegasf)
 end subroutine setup_spectral
end interface

interface
 subroutine hilbert_transform(npwe,nomega,nomegasf,my_wl,my_wr,kkweight,sf_chi0,chi0,spmeth)
  use defs_basis
  implicit none
  integer,intent(in) :: my_wl
  integer,intent(in) :: my_wr
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegasf
  integer,intent(in) :: npwe
  integer,intent(in) :: spmeth
  complex(gwpc), intent(inout) :: chi0(npwe,npwe,nomega)
  complex(dpc),intent(in) :: kkweight(nomegasf,nomega)
  complex(gwpc), intent(inout) :: sf_chi0(npwe,npwe,my_wl:my_wr)
 end subroutine hilbert_transform
end interface

interface
 subroutine hilbert_transform_headwings(npwe,nomega,nomegasf,my_wl,my_wr,kkweight,&  
  &  sf_lwing,sf_uwing,sf_head,chi0_lwing,chi0_uwing,chi0_head,spmeth)
  use defs_basis
  implicit none
  integer,intent(in) :: my_wl
  integer,intent(in) :: my_wr
  integer,intent(in) :: nomega
  integer,intent(in) :: nomegasf
  integer,intent(in) :: npwe
  integer,intent(in) :: spmeth
  complex(dpc), intent(inout) :: chi0_head(3,3,nomega)
  complex(dpc), intent(inout) :: chi0_lwing(npwe,nomega,3)
  complex(dpc), intent(inout) :: chi0_uwing(npwe,nomega,3)
  complex(dpc),intent(in) :: kkweight(nomegasf,nomega)
  complex(dpc), intent(inout) :: sf_head(3,3,my_wl:my_wr)
  complex(dpc), intent(inout) :: sf_lwing(npwe,my_wl:my_wr,3)
  complex(dpc), intent(inout) :: sf_uwing(npwe,my_wl:my_wr,3)
 end subroutine hilbert_transform_headwings
end interface

interface
 subroutine symmetrize_afm_chi0(Cryst,Gsph,Ltg_q,npwe,nomega,chi0,chi0_head,chi0_lwing,chi0_uwing)
  use m_gsphere
  use defs_basis
  use m_bz_mesh
  use m_crystal
  implicit none
  integer,intent(in) :: nomega
  integer,intent(in) :: npwe
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(in) :: Gsph
  type(little_group),intent(in) :: Ltg_q
  complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)
  complex(dpc),optional,intent(inout) :: chi0_head(3,3,nomega)
  complex(dpc),optional,intent(inout) :: chi0_lwing(npwe,nomega,3)
  complex(dpc),optional,intent(inout) :: chi0_uwing(npwe,nomega,3)
 end subroutine symmetrize_afm_chi0
end interface

interface
 subroutine assemblychi0q0_sym(nqlwl,qlwl,ik_bz,isym_kbz,itim_kbz,gwcomp,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&  
  &  chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2,lwing,uwing)
  use m_gsphere
  use defs_basis
  use m_bz_mesh
  use m_gwdefs
  use m_crystal
  implicit none
  integer,intent(in) :: gwcomp
  integer,intent(in) :: ik_bz
  integer,intent(in) :: isym_kbz
  integer,intent(in) :: itim_kbz
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nqlwl
  integer,intent(in) :: nspinor
  type(crystal_structure),intent(in) :: Cryst
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: deltaf_b1b2
  complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  complex(dpc),intent(in) :: green_enhigh_w(Ep%nomega)
  complex(dpc),intent(in) :: green_w(Ep%nomega)
  complex(dpc),intent(inout) :: lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
  real(dp),intent(in) :: qlwl(3,nqlwl)
  complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
  complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
  complex(dpc),intent(inout) :: uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)
 end subroutine assemblychi0q0_sym
end interface

interface
 subroutine assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,npwepG0,rhotwg,Gsph_epsG0,chi0)
  use m_gsphere
  use m_bz_mesh
  use m_gwdefs
  use defs_basis
  implicit none
  integer,intent(in) :: ik_bz
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nspinor
  type(epsilonm1_parameters),intent(in) :: Ep
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
  complex(dpc),intent(in) :: green_w(Ep%nomega)
  complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 end subroutine assemblychi0_sym
end interface

interface
 subroutine mkrhotwg_sigma(ii,nspinor,npw,rhotwg,rhotwg_I)
  use defs_basis
  implicit none
  integer,intent(in) :: ii
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  complex(gwpc),intent(in) :: rhotwg(npw*nspinor**2)
  complex(gwpc),intent(out) :: rhotwg_I(npw)
 end subroutine mkrhotwg_sigma
end interface

interface
 subroutine assemblychi0sfq0(nqlwl,qlwl,ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&  
  &  factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,chi0sf,lwing_sf,uwing_sf)
  use m_gsphere
  use defs_basis
  use m_bz_mesh
  use m_crystal
  implicit none
  integer,intent(in) :: ikbz
  integer,intent(in) :: iomegal
  integer,intent(in) :: iomegar
  integer,intent(in) :: isym_kbz
  integer,intent(in) :: itim_kbz
  integer,intent(in) :: my_wl
  integer,intent(in) :: my_wr
  integer,intent(in) :: nomegasf
  integer,intent(in) :: npwe
  integer,intent(in) :: npwepG0
  integer,intent(in) :: nqlwl
  integer,intent(in) :: nspinor
  integer,intent(in) :: symchi
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(in) :: Gsph_epsG0
  type(little_group),intent(in) :: Ltg_q
  real(dp),intent(in) :: factocc
  real(dp),intent(in) :: wl
  real(dp),intent(in) :: wr
  complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)
  complex(dpc),intent(inout) :: lwing_sf(npwe,my_wl:my_wr,3)
  real(dp),intent(in) :: qlwl(3,nqlwl)
  complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
  complex(gwpc),intent(in) :: rhotwx(3)
  complex(dpc),intent(inout) :: uwing_sf(npwe,my_wl:my_wr,3)
 end subroutine assemblychi0sfq0
end interface

interface
 subroutine wfk_read_ene(wfk_fname,accesswff,anbnds,energies_p,Hdr,prtvol,comm) 
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: accesswff
  integer,intent(inout) :: anbnds
  integer,intent(in) :: comm
  integer,intent(in) :: prtvol
  type(hdr_type),intent(out) :: Hdr
  character(len=fnlen),intent(in) :: wfk_fname
  real(dp),pointer :: energies_p(:,:,:)
 end subroutine wfk_read_ene
end interface

interface
 subroutine write_deltaI(unt,iqibz,deltaI,npwe)
  use defs_basis
  implicit none
  integer , intent(in) :: iqibz
  integer , intent(in) :: npwe
  integer , intent(in) :: unt
  complex(gwpc),intent(in) :: deltaI(npwe,npwe,1)
 end subroutine write_deltaI
end interface

end module interfaces_70_gw
!!***
