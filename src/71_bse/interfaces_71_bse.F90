!!****m* ABINIT/interfaces_71_bse
!! NAME
!! interfaces_71_bse
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/71_bse
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

module interfaces_71_bse

 implicit none

interface
 subroutine calc_optical_mels(Wfd,Kmesh,KS_Bst,Cryst,Psps,Pawtab,Hur,inclvkb,minb,maxb,nkbz,qpoint,opt_cvk)
  use m_wfs
  use m_bz_mesh
  use m_crystal
  use m_paw_commutator
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: inclvkb
  integer,intent(in) :: maxb
  integer,intent(in) :: minb
  integer,intent(in) :: nkbz
  type(crystal_structure),intent(in) :: Cryst
  type(bandstructure_type),intent(in) :: KS_Bst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfd
  type(hur_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
  complex(dpc),intent(out) :: opt_cvk(minb:maxb,minb:maxb,nkbz,Wfd%nsppol)
  real(dp),intent(in) :: qpoint(3)
 end subroutine calc_optical_mels
end interface

interface
 subroutine exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_Max,Gsph_c,Vcp,Wfd,W,Hdr_bse,&  
  &  nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,rhxtwg_q0,is_resonant,fname)
  use m_vcoul
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use defs_datatypes
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  use m_gsphere
  use m_screen
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nfftot_osc
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(in) :: Gsph_Max
  type(gvectors_type),intent(in) :: Gsph_c
  type(hdr_type),intent(inout) :: Hdr_bse
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(bz_mesh_type),intent(in) :: Qmesh
  type(vcoul_t),intent(in) :: Vcp
  type(screen_t),intent(inout) :: W
  type(wfs_descriptor),intent(inout) :: Wfd
  character(len=fnlen),intent(in) :: fname
  logical,intent(in) :: is_resonant
  integer,intent(in) :: ngfft_osc(18)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  integer,intent(in) :: ktabr(nfftot_osc,BSp%nkbz)
  complex(gwpc),intent(in) :: rhxtwg_q0(BSp%npweps,BSp%lomo:BSp%nbnds, &
  &         BSp%lomo:BSp%nbnds,Wfd%nkibz,Wfd%nsppol)
 end subroutine exc_build_block
end interface

interface
 subroutine exc_build_ham(BSp,BS_files,Cryst,Kmesh,Qmesh,ktabr,Gsph_Max,Gsph_c,Vcp,&  
  &  Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff)
  use m_vcoul
  use m_wfs
  use m_bz_mesh
  use defs_datatypes
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  use m_gsphere
  use m_screen
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nfftot_osc
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  type(gvectors_type),intent(in) :: Gsph_Max
  type(gvectors_type),intent(in) :: Gsph_c
  type(hdr_type),intent(inout) :: Hdr_bse
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(bz_mesh_type),intent(in) :: Qmesh
  type(vcoul_t),intent(in) :: Vcp
  type(screen_t),intent(inout) :: W
  type(wfs_descriptor),intent(inout) :: Wfd
  integer,intent(in) :: ngfft_osc(18)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  integer,intent(in) :: ktabr(nfftot_osc,BSp%nkbz)
 end subroutine exc_build_ham
end interface

interface
 subroutine exc_den(BSp,BS_files,ngfft,nfftot,Kmesh,ktabr,Wfd)
  use m_bz_mesh
  use m_bs_defs
  use m_wfs
  implicit none
  integer,intent(in) :: nfftot
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(bz_mesh_type),intent(in) :: Kmesh
  type(wfs_descriptor),intent(inout) :: Wfd
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ktabr(nfftot,BSp%nkbz)
 end subroutine exc_den
end interface

interface
 subroutine exc_diago_resonant(Bsp,BS_files,Hdr_bse,prtvol,comm)
  use defs_abitypes
  use m_bs_defs
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: prtvol
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(hdr_type),intent(in) :: Hdr_bse
 end subroutine exc_diago_resonant
end interface

interface
 subroutine exc_print_eig(BSp,bseig_fname,gw_gap,exc_gap)
  use defs_basis
  use m_bs_defs
  implicit none
  type(excparam),intent(in) :: BSp
  character(len=*),intent(in) :: bseig_fname
  complex(dpc),intent(out) :: exc_gap
  complex(dpc),intent(in) :: gw_gap
 end subroutine exc_print_eig
end interface

interface
 subroutine exc_diago_coupling(Bsp,BS_files,Hdr_bse,prtvol,comm)
  use defs_abitypes
  use m_bs_defs
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: prtvol
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(hdr_type),intent(in) :: Hdr_bse
 end subroutine exc_diago_coupling
end interface

interface
 subroutine exc_diago_coupling_hegv(Bsp,BS_files,Hdr_bse,prtvol,comm)
  use defs_abitypes
  use m_bs_defs
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: prtvol
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(hdr_type),intent(in) :: Hdr_bse
 end subroutine exc_diago_coupling_hegv
end interface

interface
 subroutine exc_diago_driver(Wfd,Bsp,BS_files,KS_BSt,QP_BSt,Cryst,Kmesh,Psps,&  
  &  Pawtab,Hur,Hdr_bse,drude_plsmf)
  use m_wfs
  use m_bz_mesh
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  use m_paw_commutator
  use defs_basis
  use defs_datatypes
  implicit none
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  type(hdr_type),intent(in) :: Hdr_bse
  type(bandstructure_type),intent(in) :: KS_BSt
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(wfs_descriptor),intent(inout) :: Wfd
  real(dp),intent(in) :: drude_plsmf
  type(hur_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 end subroutine exc_diago_driver
end interface

interface
 subroutine exc_interp_ham(BSp,BS_files,Dtset,Cryst,Kmesh,Qmesh,KS_BSt,ktabr,Gsph_Max,Gsph_c,Wfd,Hdr_bse,&  
  &  nfftot_osc,ngfft_osc,Psps,Pawtab,KS_pawrhoij,KS_Paw_ij,Pawang,Pawrad,Pawfgr,Paw_pwff,ngfftc,ngfftf,nfftf,ks_aerhor,ks_vtrial)
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in) :: nfftot_osc
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  type(dataset_type),intent(in) :: Dtset
  type(gvectors_type),intent(in) :: Gsph_Max
  type(gvectors_type),intent(in) :: Gsph_c
  type(hdr_type),intent(inout) :: Hdr_bse
  type(bandstructure_type),intent(in) :: KS_BSt
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pawfgr_type),intent(in) :: Pawfgr
  type(pseudopotential_type),intent(in) :: Psps
  type(bz_mesh_type),intent(in) :: Qmesh
  type(wfs_descriptor),intent(inout) :: Wfd
  integer,intent(in) :: ngfft_osc(18)
  integer,intent(in) :: ngfftc(18)
  integer,intent(in) :: ngfftf(18)
  type(pawrhoij_type),intent(in) :: KS_Pawrhoij(Cryst%natom)
  type(paw_ij_type),intent(in) :: KS_paw_ij(Cryst%natom)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)
  type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  real(dp),intent(in) :: ks_aerhor(nfftf,Wfd%nspden)
  real(dp),intent(in) :: ks_vtrial(nfftf,Wfd%nspden)
  integer,intent(in) :: ktabr(nfftot_osc,BSp%nkbz)
 end subroutine exc_interp_ham
end interface

interface
 subroutine exc_iterative_diago(BSp,BS_files,Hdr_bse,prtvol,comm)
  use defs_abitypes
  use m_bs_defs
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: prtvol
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(hdr_type),intent(in) :: Hdr_bse
 end subroutine exc_iterative_diago
end interface

interface
 subroutine exc_plot(Bsp,Bs_files,Wfd,Kmesh,Cryst,Psps,Pawtab,Pawrad,paw_add_onsite,spin_opt,which_fixed,eh_rcoord,nrcell,ngfftf)
  use m_wfs
  use m_bz_mesh
  use m_bs_defs
  use m_crystal
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: spin_opt
  integer,intent(in) :: which_fixed
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: Bsp
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfd
  logical,intent(in) :: paw_add_onsite
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: nrcell(3)
  type(pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
  real(dp),intent(in) :: eh_rcoord(3)
 end subroutine exc_plot
end interface

interface
 subroutine build_spectra(BSp,BS_files,Cryst,Kmesh,KS_BSt,QP_BSt,Psps,Pawtab,Wfd,Hur,drude_plsmf,comm)
  use m_wfs
  use m_bz_mesh
  use m_bs_defs
  use m_crystal
  use m_paw_commutator
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: comm
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  type(bandstructure_type),intent(in) :: KS_BSt
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_BSt
  type(wfs_descriptor),intent(inout) :: Wfd
  real(dp),intent(in) :: drude_plsmf
  type(hur_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 end subroutine build_spectra
end interface

interface
 subroutine exc_write_data(BSp,BS_files,what,eps,dos)
  use defs_basis
  use m_bs_defs
  implicit none
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  character(len=*),intent(in) :: what
  real(dp),optional,intent(in) :: dos(BSp%nomega)
  complex(dpc),intent(in) :: eps(BSp%nomega,BSp%nq)
 end subroutine exc_write_data
end interface

interface
 subroutine exc_eps_rpa(nbnds,lomo,homo,Kmesh,Bst,nq,nsppol,opt_cvk,ucvol,broad,nomega,omega,eps_rpa,dos)
  use m_bz_mesh
  use defs_datatypes
  use defs_basis
  implicit none
  integer,intent(in) :: homo
  integer,intent(in) :: lomo
  integer,intent(in) :: nbnds
  integer,intent(in) :: nomega
  integer,intent(in) :: nq
  integer,intent(in) :: nsppol
  type(bandstructure_type),intent(in) :: BSt
  type(bz_mesh_type),intent(in) :: Kmesh
  real(dp),intent(in) :: broad
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: dos(nomega)
  complex(dpc),intent(out) :: eps_rpa(nomega,nq)
  complex(dpc),intent(in) :: omega(nomega)
  complex(dpc),intent(in) :: opt_cvk(lomo:nbnds,lomo:nbnds,Kmesh%nbz,nsppol,nq)
 end subroutine exc_eps_rpa
end interface

interface
 subroutine exc_eps_resonant(Bsp,BS_files,minb,maxb,nkbz,nsppol,opt_cvk,ucvol,nomega,omega,eps_exc,dos_exc)
  use defs_basis
  use m_bs_defs
  implicit none
  integer,intent(in) :: maxb
  integer,intent(in) :: minb
  integer,intent(in) :: nkbz
  integer,intent(in) :: nomega
  integer,intent(in) :: nsppol
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: dos_exc(nomega)
  complex(dpc),intent(out) :: eps_exc(nomega,BSp%nq)
  complex(dpc),intent(in) :: omega(nomega)
  complex(dpc),intent(in) :: opt_cvk(minb:maxb,minb:maxb,nkbz,nsppol,BSp%nq)
 end subroutine exc_eps_resonant
end interface

interface
 subroutine exc_eps_coupling(Bsp,BS_files,minb,maxb,nkbz,nsppol,opt_cvk,ucvol,nomega,omega,eps_exc,dos_exc)
  use defs_basis
  use m_bs_defs
  implicit none
  integer,intent(in) :: maxb
  integer,intent(in) :: minb
  integer,intent(in) :: nkbz
  integer,intent(in) :: nomega
  integer,intent(in) :: nsppol
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: dos_exc(nomega)
  complex(dpc),intent(out) :: eps_exc(nomega,BSp%nq)
  complex(dpc),intent(in) :: omega(nomega)
  complex(dpc),intent(in) :: opt_cvk(minb:maxb,minb:maxb,nkbz,nsppol,BSp%nq)
 end subroutine exc_eps_coupling
end interface

interface
 subroutine exc_write_tensor(BSp,BS_files,what,tensor)
  use defs_basis
  use m_bs_defs
  implicit none
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  character(len=*),intent(in) :: what
  complex(dpc),intent(in) :: tensor(BSp%nomega,6)
 end subroutine exc_write_tensor
end interface

interface
 subroutine exc_haydock_driver(BSp,BS_files,Cryst,Kmesh,Hdr_bse,KS_BSt,QP_Bst,Wfd,Psps,Pawtab,Hur)
  use m_wfs
  use m_bz_mesh
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  use m_paw_commutator
  use defs_datatypes
  implicit none
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  type(hdr_type),intent(in) :: Hdr_bse
  type(bandstructure_type),intent(in) :: KS_BSt
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(in) :: QP_Bst
  type(wfs_descriptor),intent(inout) :: Wfd
  type(hur_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 end subroutine exc_haydock_driver
end interface

interface
 subroutine haydock_herm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hmat,nkets,kets,green,comm)
  use defs_basis
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: hsize
  integer,intent(in) :: my_t1
  integer,intent(in) :: my_t2
  integer,intent(in) :: nkets
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  type(hdr_type),intent(in) :: Hdr_bse
  complex(dp),intent(out) :: green(BSp%nomega,nkets)
  complex(dpc),intent(in) :: hmat(hsize,my_t1:my_t2)
  complex(dpc),intent(in) :: kets(hsize,nkets)
 end subroutine haydock_herm
end interface

interface
 subroutine haydock_herm_algo(niter_done,niter_max,nomega,omega,tol_iter,check,hsize,my_t1,my_t2,hmat,&  
  &  factor,term_type,aa,bb,phi_nm1,phi_n,green,inn,is_converged,comm)
  use defs_basis
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: hsize
  integer,intent(out) :: inn
  integer,intent(in) :: my_t1
  integer,intent(in) :: my_t2
  integer,intent(in) :: niter_done
  integer,intent(in) :: niter_max
  integer,intent(in) :: nomega
  integer,intent(in) :: term_type
  complex(dpc),intent(in) :: factor
  logical,intent(out) :: is_converged
  real(dp),intent(in) :: tol_iter
  logical,intent(in) :: check(2)
  complex(dpc),intent(inout) :: aa(niter_max)
  real(dp),intent(inout) :: bb(niter_max)
  complex(dpc),intent(out) :: green(nomega)
  complex(dpc),intent(in) :: hmat(hsize,my_t1:my_t2)
  complex(dpc),intent(in) :: omega(nomega)
  complex(dpc),intent(inout) :: phi_n(my_t2-my_t1+1)
  complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
 end subroutine haydock_herm_algo
end interface

interface
 subroutine haydock_restart(BSp,restart_file,ftype,iq_search,hsize,niter_file,aa_file,bb_file,phi_nm1_file,phi_n_file,comm)
  use defs_basis
  use m_bs_defs
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: ftype
  integer,intent(in) :: hsize
  integer,intent(in) :: iq_search
  integer,intent(out) :: niter_file
  type(excparam),intent(in) :: BSp
  character(len=*),intent(in) :: restart_file
  complex(dpc),pointer :: aa_file(:)
  real(dp),pointer :: bb_file(:)
  complex(dpc),pointer :: phi_n_file(:)
  complex(dpc),pointer :: phi_nm1_file(:)
 end subroutine haydock_restart
end interface

interface
 subroutine haydock_mdf_to_tensor(BSp,Cryst,eps,tensor_cart,tensor_red)
  use defs_basis
  use m_bs_defs
  use m_crystal
  implicit none
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  complex(dpc),intent(in) :: eps(BSp%nomega,BSp%nq)
  complex(dpc),intent(out) :: tensor_cart(BSp%nomega,6)
  complex(dpc),intent(out) :: tensor_red(BSp%nomega,6)
 end subroutine haydock_mdf_to_tensor
end interface

interface
 subroutine haydock_psherm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hreso,hcoup,nkets,kets,green,comm)
  use defs_basis
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: hsize
  integer,intent(in) :: my_t1
  integer,intent(in) :: my_t2
  integer,intent(in) :: nkets
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_structure),intent(in) :: Cryst
  type(hdr_type),intent(in) :: Hdr_bse
  complex(dp),intent(out) :: green(BSp%nomega,BSp%nq)
  complex(dpc),intent(in) :: hcoup(hsize,my_t1:my_t2)
  complex(dpc),intent(in) :: hreso(hsize,my_t1:my_t2)
  complex(dpc),intent(in) :: kets(hsize,nkets)
 end subroutine haydock_psherm
end interface

interface
 subroutine haydock_psherm_optalgo(niter_done,niter_tot,nomega,omega,tol_iter,check,hsize,my_t1,my_t2,hreso,hcoup,&  
  &  factor,term_type,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,green,inn,is_converged,comm)
  use defs_basis
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: hsize
  integer,intent(out) :: inn
  integer,intent(in) :: my_t1
  integer,intent(in) :: my_t2
  integer,intent(in) :: niter_done
  integer,intent(in) :: niter_tot
  integer,intent(in) :: nomega
  integer,intent(in) :: term_type
  complex(dpc),intent(in) :: factor
  logical,intent(out) :: is_converged
  real(dp),intent(in) :: ket0_hbar_norm
  real(dp),intent(in) :: tol_iter
  logical,intent(in) :: check(2)
  complex(dpc),intent(inout) :: aa(niter_tot)
  real(dp),intent(inout) :: bb(niter_tot+1)
  complex(dpc),intent(inout) :: cc(niter_tot+1)
  complex(dpc),intent(out) :: green(nomega)
  complex(dpc),intent(in) :: hcoup(hsize,my_t1:my_t2)
  complex(dpc),intent(in) :: hreso(hsize,my_t1:my_t2)
  complex(dpc),intent(in) :: ket0(my_t2-my_t1+1)
  complex(dpc),intent(in) :: omega(nomega)
  complex(dpc),intent(inout) :: phi_n(my_t2-my_t1+1)
  complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
  complex(dpc),intent(inout) :: phi_np1(my_t2-my_t1+1)
 end subroutine haydock_psherm_optalgo
end interface

interface
 subroutine setup_bse(codvsn,acell,rprim,ngfftf,ngfft_osc,Dtset,Dtfil,BS_files,Psps,Pawtab,BSp,&  
  &  Cryst,Kmesh,Qmesh,KS_BSt,QP_bst,Hdr_kss,Gsph_Max,Gsph_c,Vcp,Hdr_bse,w_fname,comm,Wvl)
  use m_vcoul
  use m_bz_mesh
  use defs_abitypes
  use m_bs_defs
  use m_gsphere
  use m_crystal
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: comm
  type(excfiles),intent(out) :: BS_files
  type(excparam),intent(inout) :: Bsp
  type(crystal_structure),intent(out) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(gvectors_type),intent(out) :: Gsph_Max
  type(gvectors_type),intent(out) :: Gsph_c
  type(hdr_type),intent(out) :: Hdr_bse
  type(hdr_type),intent(out) :: Hdr_kss
  type(bandstructure_type),intent(out) :: KS_BSt
  type(bz_mesh_type),intent(out) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(bandstructure_type),intent(out) :: QP_Bst
  type(bz_mesh_type),intent(out) :: Qmesh
  type(vcoul_t),intent(out) :: Vcp
  type(wvl_internal_type), intent(in) :: Wvl
  character(len=6),intent(in) :: codvsn
  character(len=fnlen),intent(out) :: w_fname
  integer,intent(out) :: ngfft_osc(18)
  integer,intent(in) :: ngfftf(18)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine setup_bse
end interface

interface
 subroutine check_kramerskronig(n,o,eps)
  use defs_basis
  implicit none
  integer,intent(in) :: n
  complex(dpc),intent(in) :: eps(n)
  real(dp),intent(in) :: o(n)
 end subroutine check_kramerskronig
end interface

interface
 subroutine check_fsumrule(n,o,e2,omegaplasma)
  use defs_basis
  implicit none
  integer,intent(in) :: n
  real(dp),intent(in) :: omegaplasma
  real(dp),intent(in) :: e2(n)
  real(dp),intent(in) :: o(n)
 end subroutine check_fsumrule
end interface

end module interfaces_71_bse
!!***
