!{\src2tex{textfont=tt}}
!!****f* ABINIT/exc_interp_ham
!! NAME
!!  exc_interp_ham
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  BSp<excparam>=The parameters for the Bethe-Salpeter calculation. 
!!  BS_files<excfiles>=File names internally used in the BS code.
!!  nsppol=Number of independent spin polarizations.
!!  Cryst<crystal_structure>=Info on the crystalline structure.
!!  Kmesh<BZ_mesh_type>=The list of k-points in the BZ, IBZ and symmetry tables.
!!  Qmesh<BZ_mesh_type>=The list of q-points for epsilon^{-1} and related symmetry tables. 
!!  ktabr(nfftot_osc,BSp%nkbz)=The FFT index of $(R^{-1}(r-\tau))$ where R is symmetry needed to obtains 
!!    the k-points from the irreducible image.  Used to symmetrize u_Sk where S = \transpose R^{-1}
!!  Gsph_Max<gvectors_type>=Info on the G-sphere used to describe wavefunctions and W (the largest one is actually stored).  
!!  Gsph_c<gvectors_type>=Info on the G-sphere used to describe the correlation part.
!!  nfftot_osc=Total Number of FFT points used for the oscillator matrix elements.
!!  ngfft_osc(18)=Info on the FFT algorithm used to calculate the oscillator matrix elements.
!!  Psps<Pseudopotential_type>=Variables related to pseudopotentials
!!  Pawtab(Psps%ntypat)<pawtab_type>=PAW tabulated starting data.
!!  Pawang<pawang_type>=PAW angular mesh and related data.
!!  Paw_pwff(Cryst%ntypat*Wfd%usepaw)<Paw_pwff_type>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  Wfd<wfs_descriptor>=Handler for the wavefunctions.
!!  nfftf=Number of points on the fine FFT mesh used for ks_vtrial
!!  ngfftf(18)=Info on the fine FFT mesh.
!!  ks_vtrial(nfftf,nspden)=Local part of the Kohn-Sham Hamiltonian.
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      bst_plot_bands,bstruct_clean,destroy_bz_mesh_type,get_dos,ks_intp_free
!!      make_mesh,make_path,shexc_free,shexc_init,shexc_solve,shirley_interp
!!      wfd_bloch_to_shirley,wfd_destroy
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine exc_interp_ham(BSp,BS_files,Dtset,Cryst,Kmesh,Qmesh,KS_BSt,ktabr,Gsph_Max,Gsph_c,Wfd,Hdr_bse,&
& nfftot_osc,ngfft_osc,Psps,Pawtab,KS_pawrhoij,KS_Paw_ij,Pawang,Pawrad,Pawfgr,Paw_pwff,ngfftc,ngfftf,nfftf,ks_aerhor,ks_vtrial)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_bs_defs
 use m_xmpi
 use m_errors
 use m_shirley
 use m_shexc

 use m_gwdefs,            only : czero_gw, GW_Q0_DEFAULT
 use m_crystal,           only : crystal_structure
 use m_gsphere,           only : gvectors_type, gsph_fft_tabs
 use m_ebands,            only : bstruct_clean, bstruct_init, pack_eneocc, bst_plot_bands, get_dos
 use m_vcoul,             only : vcoul_t, vcoul_init, vcoul_free, vcoul_nullify
 use m_bz_mesh,           only : bz_mesh_type, get_BZ_item, make_path, make_mesh, destroy_bz_mesh_type, print_BZ_mesh, find_qmesh
 use m_paw_pwij,          only : paw_pwff_type, paw_pwij_type, init_paw_pwij, destroy_paw_pwij, paw_rho_tw_g
 use m_wfs,               only : wfs_descriptor, wfd_get_ur, wfd_get_cprj, wfd_change_ngfft, wfd_ihave_ur, &
&                                wfd_distribute_bbp, wfd_destroy, wfd_nullify
 use m_commutator_vkbr,   only : kb_potential, nullify_kb_potential, destroy_kb_potential, init_kb_potential,nc_ihr_comm
 use m_oscillators,       only : rho_tw_g
 use m_screen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_interp_ham'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_osc,nfftf
 type(dataset_type),intent(in) :: Dtset
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(BZ_mesh_type),intent(in) :: Kmesh,Qmesh
 type(Bandstructure_type),intent(in) :: KS_BSt
 type(crystal_structure),intent(in) :: Cryst
 type(Gvectors_type),intent(in) :: Gsph_Max,Gsph_c
 type(Pseudopotential_type),intent(in) :: Psps
 type(Hdr_type),intent(inout) :: Hdr_bse
 type(pawang_type),intent(in) :: Pawang
 type(wfs_descriptor),intent(inout) :: Wfd
 type(Pawfgr_type),intent(in) :: Pawfgr
!arrays
 integer,intent(in) :: ngfft_osc(18),ngfftc(18),ngfftf(18)
 integer,intent(in) :: ktabr(nfftot_osc,BSp%nkbz)
 real(dp),intent(in) :: ks_vtrial(nfftf,Wfd%nspden)
 real(dp),intent(in) :: ks_aerhor(nfftf,Wfd%nspden)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)
 type(Pawrhoij_type),intent(in) :: KS_Pawrhoij(Cryst%natom)
 type(Paw_ij_type),intent(in) :: KS_paw_ij(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0,map2sphere=1,dim_rtwg1=1,brav1=1,nshiftk1=1,k1=1
 integer,parameter :: istwfk1=1,nomega1=1
 integer :: nproc,nspinor,nsppol,comm,prtvol
 !integer :: use_padfft !,iw !band,
 integer :: my_rank,usepaw,natom
 integer :: kpath_nbounds,intp_mband,intp_nk,kptopt !kpath_ntot,nband_k,
 !integer :: ik_ibz,ik_bz,isym_k,itim_k !,eh1,eh2
 real(dp) :: sh_coverage,eh_coverage
 integer :: spline_opt,min_bsize
 !integer :: sh_size !,iqlwl
 !integer :: iq_bz,iq_ibz,isym_q,itim_q
 !integer :: spin !,ikpt !ierr,
 !real(dp) :: q0vol,fcc_const
 logical :: do_resonant,do_coupling,use_coulomb_term,use_exchange_term
 !character(len=fnlen) :: fname
! character(len=500) :: msg
 !type(Bandstructure_type) :: BSt_intp
 type(wfs_descriptor) :: Wsh
 type(BZ_mesh_type) :: New_Kmesh
!arrays
 integer :: kptrlatt(3,3),bidx(4) !,nsh(Wfd%nsppol)
 !integer :: got(Wfd%nproc) !,bbp_distrb(Wfd%mband,Wfd%mband)
 integer,allocatable :: intp_nband(:,:) !,task_distrib(:,:,:,:)
 integer,allocatable :: kpath_ndiv(:)
 real(dp),allocatable :: kpath_bounds(:,:) !,intp_wtk(:)
 real(dp) :: shiftk(3,1) !,kbz(3),spinrot_k(4),
 !real(dp) :: dos(BSp%nomega)
 !complex(dpc) :: eps_rpa(BSp%nomega,BSp%nq)
 real(dp),pointer :: intp_kpt(:,:)
 !real(dp),allocatable :: intp_ene(:,:,:)
 !complex(gwpc) :: ihrc(3,Wfd%nspinor**2)
 !complex(dpc),allocatable :: ks_ihrc(:,:,:)
 !complex(dpc),allocatable :: wij(:,:,:)
 !complex(gwpc),allocatable :: vphi(:)
 !complex(gwpc),allocatable :: vc_q0(:)
 !logical :: bbp_mask(Wfd%mband,Wfd%mband)
 !type(Cprj_type),allocatable :: Cp1(:,:),Cp2(:,:)
 !type(Paw_pwij_type),allocatable :: Pwij_q0(:)
 !type(ksintp_t),allocatable :: KS_intp(:,:)
 !type(Cprj_type),allocatable :: Cp_v(:,:),Cp_c(:,:)
 type(shexc_t) :: SHexc

!************************************************************************

 MSG_WARNING("Calling experimental code shirley interpolation!")

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(nfftot_osc==PRODUCT(ngfft_osc(1:3)),"mismatch in FFT size")

 ABI_UNUSED(Qmesh%nbz)
 ABI_UNUSED(ktabr(1,1))
 ABI_UNUSED(Paw_pwff(1)%nq_spl)
 ABI_UNUSED(Hdr_bse%intxc)

 if (BSp%have_complex_ene) then
   MSG_ERROR("Complex energies are not supported yet")
 end if

 my_rank = Wfd%my_rank
 nproc   = Wfd%nproc
 comm    = Wfd%comm

 nspinor = Wfd%nspinor
 nsppol  = Wfd%nsppol
 usepaw  = Wfd%usepaw
 natom   = Wfd%natom

 do_resonant = (BS_files%in_hreso == BSE_NOFILE)
 do_coupling = (BS_files%in_hcoup == BSE_NOFILE)

 if (.not. do_resonant .and. .not. do_coupling) then
   MSG_COMMENT("Skipping calculation of both resonant and coupling block")
   RETURN
 end if

 if (.FALSE.) then ! Interpolate the band structure.
   kpath_nbounds = 7
   ABI_ALLOCATE(kpath_bounds,(3,kpath_nbounds))
   ABI_ALLOCATE(kpath_ndiv,(kpath_nbounds-1))

   kpath_bounds(:,1) = (/zero, zero, zero/)   ! #Gamma
   kpath_bounds(:,2) = (/half, zero, zero/)   ! #M
   kpath_bounds(:,3) = (/one/3, one/3, zero/) ! #K
   kpath_bounds(:,4) = (/one/3, one/3, half/) ! #H
   kpath_bounds(:,5) = (/half, zero, half/)   ! #L
   kpath_bounds(:,6) = (/zero, zero, half/)   ! #A
   kpath_bounds(:,7) = (/zero, zero, zero/)   ! #Gamma
   call make_path(kpath_nbounds,kpath_bounds,Cryst%gmet,"G",5,kpath_ndiv,intp_nk,intp_kpt)
   ABI_DEALLOCATE(kpath_bounds)
   ABI_DEALLOCATE(kpath_ndiv)

 else 
   kptrlatt(:,:) = 0
   kptrlatt(1,1) = 2 !4  ! 20 !
   kptrlatt(2,2) = 2 !4  ! 20 !
   kptrlatt(3,3) = 2 !4  ! 20 !
   shiftk(:,1)  = (/zero,zero,zero/)

   kptopt = 1
   call make_mesh(New_Kmesh,Cryst,kptopt,kptrlatt,nshiftk1,shiftk,want_tetra=0)
   !call print_BZ_mesh(New_Kmesh,header="New Kmesh")

   intp_nk  = New_Kmesh%nibz
   ABI_ALLOCATE(intp_kpt,(3,intp_nk))
   intp_kpt = New_Kmesh%ibz
 end if

 ABI_ALLOCATE(intp_nband,(intp_nk,nsppol))
 intp_nband=10
 intp_mband = MAXVAL( MAXVAL(intp_nband, DIM=1))

 !sh_coverage=0.9999
 !sh_coverage=0.999
 !sh_coverage=0.99
 sh_coverage=0.9
 !sh_coverage=one
 min_bsize = intp_mband

 call wfd_bloch_to_shirley(Wfd,Cryst,Kmesh,KS_Bst,Psps,Pawtab,Pawang,Pawrad,min_bsize,sh_coverage,Wsh)

 spline_opt=0
#if 0
 ! KS_intp stores the interpolated eigenvalues and the transformation Optimal set --> KS states.
 ABI_ALLOCATE(KS_intp,(intp_nk,nsppol))
 ABI_ALLOCATE(intp_ene,(intp_mband,intp_nk,nsppol))

 ABI_ALLOCATE(intp_wtk,(intp_nk))
 intp_wtk=New_Kmesh%wt ! FIXME Be careful here.

 call shirley_interp(Wsh,"Vectors",Dtset,Cryst,Psps,Pawtab,Pawfgr,Pawang,Pawrad,&
&  KS_pawrhoij,KS_Paw_ij,ngfftc,ngfftf,nfftf,ks_vtrial,spline_opt,&
&  intp_nband,intp_mband,intp_nk,intp_kpt,intp_ene,KS_intp,intp_wtk=intp_wtk,BSt_intp=BSt_intp) 

 call ks_intp_free(KS_intp)
 ABI_DEALLOCATE(KS_intp)

 ABI_DEALLOCATE(intp_wtk)

 call bst_plot_bands(BSt_intp,Cryst%gmet,"interpolated",ierr)

 ! FIXME Very Bad scaling wrt nkpt!
 call get_dos(BSt_intp,New_Kmesh,1,"test_dos",0.1/Ha_eV,0.05/Ha_eV)

 call destroy_BZ_mesh_type(New_Kmesh)
 call bstruct_clean(BSt_intp)
                                          
 ABI_DEALLOCATE(intp_kpt)
 ABI_DEALLOCATE(intp_ene)
#endif

 ABI_DEALLOCATE(intp_nband)
 !
 ! ========================================================
 ! === Setup of the basis set for the eh representation ===
 ! ========================================================

 eh_coverage = 0.9
 !eh_coverage=0.999
 !eh_coverage=0.9999

 use_coulomb_term  = .TRUE.
 use_exchange_term = .TRUE. 
 prtvol=0

 bidx = (/Bsp%lomo,Bsp%homo,Bsp%lumo,Bsp%humo/)

 call shexc_init(SHexc,eh_coverage,Dtset,Wsh,Cryst,Gsph_Max,Gsph_c,Psps,Pawtab,KS_pawrhoij,KS_Paw_ij,Pawang,Pawrad,Pawfgr,&
  &  kptopt,kptrlatt,nshiftk1,shiftk,nspinor,nsppol,Wfd%nspden,prtvol,use_coulomb_term,use_exchange_term,bidx,Bsp%q,&
  &  nfftf,ngfftf,ngfftc,ks_aerhor,ks_vtrial)

 call shexc_solve(SHexc,Bsp,Cryst)

 ! Master node writes final results on file.
 !call exc_write_data(BSp,BS_files,"RPA_NLF_MDF",eps_rpa,dos=dos)

 call shexc_free(Shexc)

 ! * Free memory.
 !
 ! Deallocate Shirley basis set.
 call wfd_destroy(Wsh)

 MSG_ERROR("Interpolation completed")

end subroutine exc_interp_ham
!!***
