!{\src2tex{textfont=tt}}
!!****f* ABINIT/sigma
!! NAME
!! sigma
!!
!! FUNCTION
!! Calculate the matrix elements of the self-energy operator.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, MT, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! codvsn=code version
!! Dtfil<type(datafiles_type)>=variables related to files
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Pawang<type(pawang_type)>=paw angular mesh and related data
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data
!! Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!! Psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in sigma, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!   the call to pspini. The next time the code enters screening, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=dimensionless real space primitive translations
!! xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  Converged=.TRUE. if degw are within the user-specified tolerance.
!!  Output is written on the main abinit output file. Some results are stored in external files
!!
!! PARENTS
!!      driver,gw_driver
!!
!! NOTES
!!
!! ON THE USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ... Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! CHILDREN
!!      abi_etsf_electrons_put,abi_etsf_geo_put,abi_etsf_init,bstruct_clean
!!      calc_sigc_me,calc_sigx_me,calc_vhxc_me,check_sym_ug,chkpawovlp
!!      classify_bands,cohsex_me,copy_bandstructure,cprj_alloc,cprj_free
!!      cutoff_density,denfgr,destroy_bands_symmetries,destroy_bz_mesh_type
!!      destroy_crystal,destroy_epsilonm1_results,destroy_gsphere
!!      destroy_little_group,destroy_melements,destroy_paw_an,destroy_paw_ij
!!      destroy_paw_pwaves_lmn,destroy_paw_pwff,destroy_sigma_parameters
!!      destroy_sigma_results,energies_init,etsf_dump_qp,fourdp,get_gftt
!!      get_rhor,getem1_from_ppm,getph,hdr_clean,init_gsphere,init_paw_an
!!      init_paw_ij,init_paw_pwaves_lmn,init_paw_pwff,init_pawfgr
!!      init_sigma_results,initmpi_seq,ioarr,metric,mkdump_er,mkrdim,nhatgrid
!!      nullify_bands_symmetries,nullify_little_group,nullify_paw_an
!!      nullify_paw_ij,paw_check_symcprj,paw_dijhf,paw_mkdijexc_core,paw_qpscgw
!!      pawdenpot,pawdij,pawfgrtab_free,pawfgrtab_init,pawfgrtab_print,pawinit
!!      pawmknhat,pawprt,pawpuxinit,ppm_free,ppm_init,print_melements
!!      print_ngfft,print_paw_ij,print_pawtab,print_psps,prtrhomxmn,pspini,rdgw
!!      rdqps,recalculate_epsm1_freq_grid,reportgap,reset_mflags,rhoij_alloc
!!      rhoij_copy,rhoij_free,setsymrhoij,setup_little_group,setup_ppmodel
!!      setup_sigma,setvtr,show_qp,sigma_tables,solve_dyson,split_work2,symdij
!!      symdij_all,test_charge,timab,update_occ,updt_m_lda_to_qp,vcoul_free
!!      wfd_change_ngfft,wfd_copy,wfd_destroy,wfd_get_cprj,wfd_init,wfd_mkrho
!!      wfd_plot_ur,wfd_print,wfd_read_kss,wfd_read_wfk,wfd_reset_ur_cprj
!!      wfd_rotate,wfd_test_ortho,write_sigma_results
!!      write_sigma_results_header,wrqps,wrtout,xbarrier_mpi,xc_kernel
!!      xc_kernel_ada,zero_melements
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine sigma(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,xred,converged)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_timer
 use m_errors
#if defined HAVE_DFT_LIBXC
 use libxc_functionals
#endif

 use m_numeric_tools, only : imax_loc, iseven
 use m_blas,          only : xdotc
 use m_io_tools,      only : get_unit
 use m_header,        only : hdr_clean
 use m_geometry,      only : normv
 use m_crystal,       only : destroy_crystal, crystal_structure
 use m_ebands,        only : update_occ, copy_bandstructure, reportgap, get_valence_idx, get_bandenergy, &
&                            bstruct_clean, bstruct_init, bst_print_fs
 use m_energies,      only : energies_type, energies_init
 use m_bz_mesh,       only : bz_mesh_type, findqg0, destroy_bz_mesh_type, little_group, setup_little_group, &
&                            nullify_little_group, destroy_little_group, make_path
 use m_gsphere,       only : gvectors_type, destroy_gsphere, init_gsphere
 use m_vcoul,         only : vcoul_t, vcoul_free, cutoff_density
 use m_qparticles,    only : wrqps, rdqps, rdgw, show_QP, updt_m_lda_to_qp
 use m_fft_mesh,      only : get_gftt, print_ngfft
 use m_io_kss,        only : wfd_read_kss
 use m_screening,     only : mkdump_er, destroy_epsilonm1_results, epsilonm1_results, recalculate_epsm1_freq_grid
 use m_ppmodel,       only : ppm_init, ppm_free, setup_ppmodel, getem1_from_PPm, ppmodel_type
 use m_wfs,           only : wfd_init, wfd_destroy, wfd_reset_ur_cprj, wfd_print, wfs_descriptor,&
&                            wfd_rotate, wfd_get_cprj, wfd_iam_master, wfd_change_ngfft, wfd_plot_ur, wfd_test_ortho,&
&                            wfd_read_wfk, check_sym_ug, wfd_copy
 use m_paw_dmft,      only : paw_dmft_type
 use m_paw_toolbox,   only : nullify_paw_ij, init_paw_ij, destroy_paw_ij, init_pawfgr, &
&                            pawfgrtab_free, pawfgrtab_init, pawfgrtab_print,&
&                            nullify_paw_an, init_paw_an, destroy_paw_an, print_pawtab, print_paw_ij,&
&                            paw_pwaves_lmn_t,init_paw_pwaves_lmn,destroy_paw_pwaves_lmn,nullify_paw_pwaves_lmn
 use m_paw_pwij,      only : paw_pwff_type, init_paw_pwff, destroy_paw_pwff
 use m_paw_slater,    only : paw_mkdijexc_core, paw_dijhf
 use m_melemts,       only : reset_mflags, print_melements, destroy_melements, melements_flags_type, melements_type,&
&                            zero_melements
 use m_sigma_results, only : write_sigma_results_header, write_sigma_results, init_sigma_results, destroy_sigma_results,&
&                            etsf_dump_qp, sigma_results
 use m_dyson_solver,  only : solve_dyson
 use m_bands_sym,     only : bands_symmetries, nullify_bands_symmetries, destroy_bands_symmetries, bsym_failed

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_61_ionetcdf
 use interfaces_62_iowfdenpot
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_69_wfdesc
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(out) :: converged
 character(len=6),intent(in) :: codvsn
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pawang_type),intent(inout) :: Pawang
 type(Pseudopotential_type),intent(inout) :: Psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,Dtset%natom)
 type(Pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=40,tim_fourdp=5
 integer,parameter :: nshiftk=1
 integer,save :: nsym_old=-1
 integer :: accessfil,approx_type,b1gw,b2gw,choice,cplex,cplex_dij,band,my_nspins
 integer :: dim_kxcg,fformr,has_dijU,has_dijso,iab,bmin,bmax,irr_idx1,irr_idx2
 integer :: iat,ib,ib1,ib2,id_required,ider,idir,ii,ik,ierr
 integer :: ik_bz,ikcalc,ik_ibz,ikxc,ipert
 integer :: isp,is_idx,istat,istep,itypat,izero,jj
 integer :: first_band,last_band
 integer :: ks_iv,lmn2_size_max,master,mband
 integer :: mgfftf,mod10,moved_atm_inside!,mgfft
 integer :: moved_rhor,my_maxb,my_minb,n3xccc
 integer :: nbsc,ndij,nfftf,nfftf_tot,gwc_nfft,gwc_nfftot,gwx_nfft,gwx_nfftot
 integer :: nhatgrdim,nkxc,nkxc1,nprocs,nscf,nspden_rhoij,nzlmopt,optene
 integer :: optcut,optgr0,optgr1,optgr2,option,option_test,option_dij,optrad,optrhoij,psp_gencond
 integer :: my_rank,rdwr,rdwrpaw,rhoxsp_method,comm,use_aerhor
 integer :: use_umklp,usexcnhat,ios
 integer :: ioe0j,spin,io,jb,nomega_sigc
 integer :: iomega,ppm_unt
 !integer :: jb_qp,ib_ks,ks_irr
 real(dp) :: compch_fft,compch_sph,diecut_eff_dum,r_s,rhoav,alpha,opt_ecut
 real(dp) :: drude_plsmf,my_plsmf,dummy,ecore,ecut_eff,ecutdg_eff,ehartree
 real(dp) :: exchange_energy,gsqcutc_eff,gsqcutf_eff,nelect,norm,oldefermi
 real(dp) :: ucvol,vxcavg,vxcavg_qp
 real(dp) :: gwc_gsq,gwx_gsq,gw_gsq
 complex(dpc) :: max_degw,cdummy
 logical :: use_paw_aeur,only_one_kpt,pawden_exists,dbg_mode,pole_screening
 character(len=500) :: msg,sigma_type
 character(len=fnlen) :: wfk_fname,pawden_fname
 type(BZ_mesh_type) :: Kmesh,Qmesh
 type(Bandstructure_type) :: KS_BSt,QP_BSt
 type(vcoul_t) :: Vcp
 type(Crystal_structure) :: Cryst
 type(Energies_type) :: KS_energies,QP_energies
 type(Epsilonm1_results) :: Er
 type(Gvectors_type) :: Gsph_Max,Gsph_x,Gsph_c
 type(Hdr_type) :: Hdr_kss,Hdr_sigma
 type(melements_flags_type) :: KS_mflags,QP_mflags
 type(melements_type) :: KS_me,QP_me
 type(MPI_type) :: MPI_enreg_seq 
 !type(MPI_type) :: Denfgr_MPI
 type(paw_dmft_type) :: Paw_dmft
 type(Pawfgr_type) :: Pawfgr
 type(PPmodel_type) :: PPm
 type(Sigma_parameters) :: Sigp
 type(Sigma_results) :: Sr
 type(wfs_descriptor) :: Wfd,Wfdf
 type(Wvl_data) :: Wvl
!arrays
 integer,save :: paw_gencond(6)=(/-1,-1,-1,-1,-1,-1/)
 integer :: gwc_ngfft(18),ngfftc(18),ngfftf(18),gwx_ngfft(18),my_spins(Dtset%nsppol)
 integer,allocatable :: nq_spl(:)
 integer,allocatable :: tmp_gfft(:,:)
 integer,allocatable :: ks_vbik(:,:),nband(:,:)  !istart(:),istop(:)
 integer,allocatable :: l_size_atm(:),nlmn_type(:),qp_vbik(:,:)
 integer,allocatable :: tmp_kstab(:,:,:)
 integer,allocatable :: ks_irreptab(:,:,:),qp_irreptab(:,:,:) 
 real(dp),parameter ::  k0(3)=zero
 !real(dp) :: shiftk(3,nshiftk)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),strsxc(6),tsec(2)
 real(dp),allocatable :: grewtn(:,:),qmax(:)
 real(dp),allocatable :: ks_nhat(:,:),ks_nhatgr(:,:,:),ks_rhog(:,:)
 real(dp),allocatable :: ks_rhor(:,:),ks_vhartr(:),ks_vtrial(:,:),ks_vxc(:,:)
 real(dp),allocatable :: ks_taug(:,:),ks_taur(:,:)
 real(dp),allocatable :: kxc(:,:),qp_kxc(:,:),ph1d(:,:),ph1df(:,:)
 real(dp),allocatable :: prev_rhor(:,:),prev_taur(:,:),qp_nhat(:,:)
 real(dp),allocatable :: qp_nhatgr(:,:,:),qp_rhog(:,:),qp_rhor_paw(:,:)
 real(dp),allocatable :: qp_rhor_n_one(:,:),qp_rhor_nt_one(:,:)
 real(dp),allocatable :: qp_rhor(:,:),qp_vhartr(:),qp_vtrial(:,:),qp_vxc(:,:)
 real(dp),allocatable :: qp_taur(:,:),qp_taug(:,:)
 real(dp),allocatable :: vpsp(:),xccc3d(:),xred_dummy(:,:)
 real(dp),allocatable :: dijexc_core(:,:,:)
 real(dp),allocatable :: dij_hf(:,:,:)
 real(dp),allocatable :: igwene(:,:,:)
 complex(dpc),allocatable :: omega(:),em1_ppm(:,:,:) 
 real(dp),allocatable :: ks_aepaw_rhor(:,:) !,ks_n_one_rhor(:,:),ks_nt_one_rhor(:,:)
 complex(dpc) :: ovlp(2)
 complex(dpc),allocatable :: ctmp(:,:),hbare(:,:,:,:) 
 complex(dpc),target,allocatable :: sigxme(:,:,:,:),sigcme(:,:,:,:,:)
 complex(dpc),allocatable :: hlda(:,:,:,:),htmp(:,:,:,:)
 complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:)
 complex(dpc),allocatable :: uks2qp(:,:)
 complex(gwpc),allocatable :: kxcg(:,:),fxc_ADA(:,:,:)
 complex(gwpc),pointer :: ug1(:)
 complex(dpc),pointer :: sigxme_p(:,:,:),sigcme_p(:,:,:,:)
 real(dp),pointer :: qp_ene(:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(Bands_Symmetries),target,allocatable :: KS_sym(:,:)
 type(Bands_Symmetries),pointer :: QP_sym(:,:)
 type(Cprj_type),allocatable :: Cp1(:,:),Cp2(:,:)
 type(Little_group),allocatable :: Ltg_k(:)
 type(Paw_an_type),allocatable :: KS_paw_an(:),QP_paw_an(:)
 type(Paw_ij_type),allocatable :: KS_paw_ij(:),QP_paw_ij(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(Pawrhoij_type),allocatable :: KS_Pawrhoij(:),Pawrhoij_dum(:)
 type(Pawrhoij_type),allocatable :: QP_pawrhoij(:),prev_Pawrhoij(:)
 type(Paw_pwff_type),allocatable :: Paw_pwff(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

!************************************************************************

 DBG_ENTER('COLL')

 ABI_TIMER_START("")
 call timab(401,1,tsec) ! sigma(Total)

 ABI_TIMER_START("init1")
 call timab(402,1,tsec) ! sigma(Init1)

 write(msg,'(7a)')&
& ' SIGMA: Calculation of the GW corrections ',ch10,ch10,&
& ' Based on a program developped by R.W. Godby, V. Olevano, G. Onida, and L. Reining.',ch10,&
& ' Incorporated in ABINIT by V. Olevano, G.-M. Rignanese, and M. Torrent.'
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
#if defined HAVE_GW_DPC
 if (gwpc/=8) then
   write(msg,'(6a)')ch10,&
&   ' Number of bytes for double precision complex /=8 ',ch10,&
&   ' Cannot continue due to kind mismatch in BLAS library ',ch10,&
&   ' Some BLAS interfaces are not generated by abilint '
   MSG_ERROR(msg)
 end if
 write(msg,'(a,i2,a)')'.Using double precision arithmetic ; gwpc = ',gwpc,ch10
#else
 write(msg,'(a,i2,a)')'.Using single precision arithmetic ; gwpc = ',gwpc,ch10
#endif
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
!
!=== Initialize MPI variables, and parallelization level ===
!* gwpara 0--> sequential run, 1--> parallelism over k-points, 2--> parallelism over bands.
!* In case of gwpara==1 memory is not parallelized.
!* If gwpara==2, bands are divided among processors but each proc has all the states where GW corrections are required.
 comm = xmpi_world
 my_rank = xcomm_rank(comm)
 nprocs  = xcomm_size(comm)
 master=0

!Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 converged = .FALSE.
!
!accesswff defines the format of the output.
!1--> Plain Fortran file
!2--> Set all outputs to netcdf format (not implemented)
!3--> Set all outputs to ETSF format
 accessfil=0
 if (Dtset%accesswff==IO_MODE_NETCDF) accessfil=1
 if (Dtset%accesswff==IO_MODE_ETSF  ) accessfil=3
 if (Dtset%accesswff==IO_MODE_MPI   ) accessfil=4

!=== Some variables need to be initialized/nullify at start ===
 call energies_init(KS_energies)
 usexcnhat=0
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)
!
!=== Define FFT grid(s) sizes ===
!* Be careful! This mesh is only used for densities, potentials and the matrix elements of v_Hxc. It is NOT the
!(usually coarser) GW FFT mesh employed for the oscillator matrix elements that is defined in setmesh.F90.
!See also NOTES in the comments at the beginning of this file.
!NOTE: This mesh is defined in invars2m using ecutwfn, in GW Dtset%ecut is forced to be equal to Dtset%ecutwfn.

 call init_pawfgr(Dtset,Pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
& gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=gmet,k0=k0)

 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')
 nfftf_tot=PRODUCT(ngfftf(1:3))
!
!=== Open and read pseudopotential files ===

 call pspini(Dtset,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,level,Pawrad,Pawtab,Psps,rprimd)
 if (psp_gencond==1) call print_psps(Psps,std_out,0,'COLL')

 call timab(402,2,tsec) ! Init1
!
!===============================================
!==== Initialize Sigp, Er and basic objects ====
!===============================================
!* Sigp is completetly initialized here.
!* Er is only initialized with dimensions, (SCR|SUSC) file is read in mkdump_Er
 ABI_TIMER_START("setup")
 call timab(403,1,tsec) ! setup_sigma

 call setup_sigma(codvsn,acell,rprim,ngfftf,Dtset,Dtfil,Psps,Pawtab,&
& gwx_ngfft,gwc_ngfft,Hdr_kss,Hdr_sigma,Cryst,Kmesh,Qmesh,KS_BSt,Gsph_Max,Gsph_c,Vcp,Er,Sigp,comm)

 ABI_TIMER_STOP("setup")
 call timab(403,2,tsec) ! setup_sigma
 call timab(402,1,tsec) ! Init1

!XG090617 Please, do not remove this write, unless you have checked
!that the code executes correctly on max+g95 (especially, Tv5#70).
!It is one more a silly write, perhaps needed because the compiler does not treat correctly
!non-nullified pointers.
 if (sigma_needs_w(Sigp)) write(std_out,*)' screening after setup_sigma : Er%Hscr%headform=',Er%Hscr%headform
!END XG090617

 pole_screening = .FALSE.
 if (Er%fform==2002) then
   pole_screening = .TRUE.
   MSG_WARNING(' EXPERIMENTAL - Using a pole-fit screening!')
 end if

 call print_ngfft(gwc_ngfft,header='FFT mesh for oscillator strengths used for Sigma_c')
 call print_ngfft(gwx_ngfft,header='FFT mesh for oscillator strengths used for Sigma_x')

 mod10=MOD(Sigp%gwcalctyp,10)
 b1gw=Sigp%minbdgw
 b2gw=Sigp%maxbdgw

 gwc_nfftot=PRODUCT(gwc_ngfft(1:3))
 gwc_nfft  =gwc_nfftot  !no FFT //

 gwx_nfftot=PRODUCT(gwx_ngfft(1:3))
 gwx_nfft  =gwx_nfftot  !no FFT //
!
!TRYING TO RECREATE AN "ABINIT ENVIRONMENT"
 KS_energies%e_corepsp=ecore/Cryst%ucvol

!=== Calculate KS occupation numbers and ks_vbk(nkibz,nsppol) ====
!* ks_vbk gives the (valence|last Fermi band) index for each k and spin.
!* fixmom is passed to fermi.F90 to fix the problem with newocc in case of magnetic metals
 ABI_ALLOCATE(ks_vbik,(KS_BSt%nkpt,KS_BSt%nsppol))
 ABI_ALLOCATE(qp_vbik,(KS_BSt%nkpt,KS_BSt%nsppol))

 call update_occ(KS_BSt,Dtset%fixmom,prtvol=0)
 ks_vbik(:,:) = get_valence_idx(KS_BSt)
!
!============================
!==== PAW initialization ====
!============================

 if (Dtset%usepaw==1) then
   call chkpawovlp(Cryst%natom,Cryst%ntypat,Dtset%pawovlp,Pawtab,Cryst%rmet,Cryst%typat,Cryst%xred)

   ABI_ALLOCATE(nlmn_type,(Cryst%ntypat))
   do itypat=1,Cryst%ntypat
     nlmn_type(itypat)=Pawtab(itypat)%lmn_size
   end do

   cplex_dij=Dtset%nspinor; cplex=1; ndij=1

   ABI_ALLOCATE(KS_Pawrhoij,(Cryst%natom))
   nspden_rhoij=Dtset%nspden; if (Dtset%pawspnorb>0.and.Dtset%nspinor==2) nspden_rhoij=4

   call rhoij_alloc(Dtset%pawcpxocc,nlmn_type,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,KS_Pawrhoij,Cryst%typat)
   ABI_DEALLOCATE(nlmn_type)

!  === Initialize values for several basic arrays ===
!  TODO Check pawxcdev>2 since gaunt coefficients are allocated with different size
   if (psp_gencond==1.or.&
&   paw_gencond(1)/=Dtset%pawlcutd .or.paw_gencond(2)/=Dtset%pawlmix  .or.&
&   paw_gencond(3)/=Dtset%pawnphi  .or.paw_gencond(4)/=Dtset%pawntheta.or.&
&   paw_gencond(5)/=Dtset%pawspnorb.or.paw_gencond(6)/=Dtset%pawxcdev) then

     call timab(553,1,tsec)
     diecut_eff_dum=ABS(Dtset%diecut)*Dtset%dilatmx**2

     call pawinit(diecut_eff_dum,Psps%indlmn,Dtset%pawlcutd,Dtset%pawlmix,Psps%lmnmax,Psps%mpsang,&
&     Dtset%pawnphi,Cryst%nsym,Dtset%pawntheta,Cryst%ntypat,Pawang,Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev)

     paw_gencond(1)=Dtset%pawlcutd; paw_gencond(2)=Dtset%pawlmix
     paw_gencond(3)=Dtset%pawnphi; paw_gencond(4)=Dtset%pawntheta
     paw_gencond(5)=Dtset%pawspnorb; paw_gencond(6)=Dtset%pawxcdev
     call timab(553,2,tsec)
   else
     if (Pawtab(1)%has_kij  ==1) Pawtab(1:Cryst%ntypat)%has_kij  =2
     if (Pawtab(1)%has_nabla==1) Pawtab(1:Cryst%ntypat)%has_nabla=2
   end if
   Psps%n1xccc=MAXVAL(Pawtab(1:Cryst%ntypat)%usetcore)

!  Initialize optional flags in Pawtab to zero
!  (Cannot be done in Pawinit since the routine is called only if some pars. are changed)
   Pawtab(:)%has_nabla = 0
   Pawtab(:)%usepawu   = 0
   Pawtab(:)%useexexch = 0
   Pawtab(:)%exchmix   =zero

!  if (psp_gencond==1) then !.or. nsym_old/=Cryst%nsym) then
   call setsymrhoij(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,Cryst%rprimd,Cryst%symrec,Pawang%zarot)
   nsym_old=Cryst%nsym
!  end if

!  === Initialize and compute data for LDA+U ===
   Paw_dmft%use_dmft=Dtset%usedmft
   if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
     call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,&
&     Psps%indlmn,Psps%lmnmax,Cryst%ntypat,Pawang,Dtset%pawprtvol,Pawrad,Pawtab,Dtset%upawu,&
&     Dtset%usedmft,Dtset%useexexch,Dtset%usepawu)
   end if

   ABI_CHECK(Dtset%useexexch==0,"LEXX not yet implemented in GW")
   ABI_CHECK(Paw_dmft%use_dmft==0,"DMFT + GW not available")

   if (my_rank==master) call print_pawtab(Pawtab)
!  3-

!  Optionally read core orbitals from file and calculate $ \<\phi_i|Sigma_x^\core|\phi_j\> $ for the HF decoupling.
   if (Sigp%use_sigxcore==1) then
     lmn2_size_max=MAXVAL(Pawtab(:)%lmn2_size)
     ABI_ALLOCATE(dijexc_core,(cplex_dij*lmn2_size_max,ndij,Cryst%ntypat))

     call paw_mkdijexc_core(ndij,cplex_dij,lmn2_size_max,Cryst,Psps,Pawtab,Pawrad,dijexc_core,Dtset%prtvol)
   end if ! HF decoupling
!  
!  === Get Pawrhoij from the header of the KSS file ===
   call rhoij_copy(Hdr_kss%pawrhoij,KS_Pawrhoij,MPI_enreg=MPI_enreg_seq)

!  === Re-symmetrize symrhoij ===
!  this call leads to a SIGFAULT, likely some pointer is not initialized correctly
   choice=1; optrhoij=1; ipert=0; idir=0
!  call symrhoij(choice,Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,Psps%lmnmax,Cryst%natom,Cryst%natom,Cryst%nsym,&
!  &  Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,KS_Pawrhoij,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)

!  === Evaluate form factor of radial part of phi.phj-tphi.tphj ===
   rhoxsp_method=1 ! Arnaud-Alouani
!  rhoxsp_method=2 ! Shiskin-Kresse

!  The q-grid must contain the FFT mesh used for sigma_c and the G-sphere for the exchange part.
!  We use the FFT mesh for sigma_c since COHSEX and the extrapolar method require oscillator
!  strengths on the FFT mesh.
   ABI_ALLOCATE(tmp_gfft,(3,gwc_nfftot))
   call get_gftt(gwc_ngfft,k0,gmet,gwc_gsq,tmp_gfft)
   ABI_DEALLOCATE(tmp_gfft)

   gwx_gsq =  Dtset%ecutsigx/(two*pi**2)
!  allocate(tmp_gfft(3,gwx_nfftot)); q0=zero
!  call get_gftt(gwx_ngfft,q0,gmet,gwx_gsq,tmp_gfft)
!  deallocate(tmp_gfft)
   gw_gsq = MAX(gwx_gsq,gwc_gsq)

!  * Set up q-grid, make qmax 20% larger than largest expected.
   ABI_ALLOCATE(nq_spl,(Psps%ntypat))
   ABI_ALLOCATE(qmax,(Psps%ntypat))
   qmax = SQRT(gw_gsq)*1.2d0 ! qmax=Psps%qgrid_ff(Psps%mqgrid_ff)
   nq_spl = Psps%mqgrid_ff
!  write(std_out,*)"using nq_spl",nq_spl,"qmax=",qmax
   ABI_ALLOCATE(Paw_pwff,(Psps%ntypat))

   call init_paw_pwff(Paw_pwff,rhoxsp_method,nq_spl,qmax,gmet,Pawrad,Pawtab,Psps)

   ABI_DEALLOCATE(nq_spl)
   ABI_DEALLOCATE(qmax)
!  
!  === Variables/arrays related to the fine FFT grid ===
   ABI_ALLOCATE(ks_nhat,(nfftf,Dtset%nspden))
   ks_nhat=zero
   ABI_ALLOCATE(Pawfgrtab,(Cryst%natom))
   ABI_ALLOCATE(l_size_atm,(Cryst%natom))
   do iat=1,Cryst%natom
     l_size_atm(iat)=Pawtab(Cryst%typat(iat))%l_size
   end do
   cplex=1
   call pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,Dtset%nspden)

   ABI_DEALLOCATE(l_size_atm)
   compch_fft=greatest_real
   usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
!  * 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
!  * 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   write(msg,'(a,i2)')' sigma : using usexcnhat = ',usexcnhat
   call wrtout(std_out,msg,'COLL')
!  
!  === Identify parts of the rectangular grid where the density has to be calculated ===

   optcut=0; optgr0=Dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-Dtset%pawstgylm
   if (Dtset%pawcross==1) optrad=1
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(Cryst%atindx1,gmet,MPI_enreg_seq,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%ucvol,Cryst%xred)

   if (Dtset%pawcross==1) then
     ABI_ALLOCATE(Paw_onsite,(Cryst%natom))
     call init_paw_pwaves_lmn(Paw_onsite,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%rprimd,Cryst%xcart,&
&     Psps,Pawtab,Pawrad,Pawfgrtab)
!    call nullify_paw_pwaves_lmn(Paw_onsite)
   end if

   call pawfgrtab_print(Pawfgrtab,unit=std_out,prtvol=Dtset%pawprtvol)
 end if !End of PAW Initialization
!
!Allocate these arrays anyway, since they are passed to subroutines.
 if (.not.allocated(ks_nhat))      then
   ABI_ALLOCATE(ks_nhat,(nfftf,0))
 end if
 if (.not.allocated(dijexc_core))  then
   ABI_ALLOCATE(dijexc_core,(1,1,0))
 end if
!
!==================================================
!==== Read KS band structure from the KSS file ====
!==================================================
!
!* Initialize Wfd, allocate wavefunctions and precalculate tables to do the FFT using the coarse gwc_ngfft.
 mband=Sigp%nbnds
 ABI_ALLOCATE(bks_mask,(mband,Kmesh%nibz,Dtset%nsppol))
 bks_mask=.FALSE.
 ABI_ALLOCATE(keep_ur ,(mband,Kmesh%nibz,Sigp%nsppol))
 keep_ur=.FALSE.

 ABI_ALLOCATE(nband,(Kmesh%nibz,Sigp%nsppol))
 nband=mband

 my_nspins=Dtset%nsppol; my_spins=(/(isp,isp=1,Dtset%nsppol)/)
 my_minb=1; my_maxb=Sigp%nbnds

 select case (Dtset%gwpara)
   case (1)
     call wrtout(std_out,' sigma: parallelized over transitions','COLL')

     if (Dtset%nsppol==2.and.iseven(nprocs)) then ! Distribute spins.
       my_nspins=1
       my_spins(1)=1; if (my_rank+1>nprocs/2) my_spins(1)=2
     end if

     do isp=1,my_nspins
       spin = my_spins(isp)
       bks_mask(my_minb:my_maxb,:,spin)=.TRUE. 
       if (MODULO(Dtset%gwmem,10)==1) keep_ur(my_minb:my_maxb,:,spin)=.TRUE.
     end do

   case (2)

#if 0
     call wrtout(std_out,'sigma: loop over bands done in parallel ','COLL')
     ABI_ALLOCATE(istart,(nprocs))
     ABI_ALLOCATE(istop,(nprocs))
     call split_work2(Sigp%nbnds,nprocs,istart,istop,msg,ierr)
     if (ierr/=0) then 
       MSG_WARNING(msg)
     end if

     my_minb=istart(my_rank+1); my_maxb=istop(my_rank+1)
     ABI_DEALLOCATE(istart)
     ABI_DEALLOCATE(istop)
     if (my_minb>my_maxb) then
       write(msg,'(3a,2(i6,a),a)')&
&       ' One or more processors has zero number of bands ',ch10,&
&       ' my_minb = ',my_minb,' my_maxb = ',my_maxb,ch10,&
&       ' This is a waste, decrease the number of processors '
       MSG_ERROR(msg)
     end if
     write(msg,'(4(a,i4))')' treating ',my_maxb-my_minb+1,' bands from ',my_minb,' up to ',my_maxb,' by rank ',my_rank
     call wrtout(std_out,msg,'PERS')
     bks_mask(my_minb:my_maxb,:,:)=.TRUE. ! bands are block-distributed.
#endif

!    alternating planes of bands
     do band=1,mband 
       if (MODULO(band,nprocs)==my_rank) then 
         bks_mask(band,:,:)=.TRUE.
         if (MODULO(Dtset%gwmem,10)==1) keep_ur(band,:,:)=.TRUE.
       end if
     end do

!    bks_mask=.TRUE.; if (MODULO(Dtset%gwmem,10)==1) keep_ur=.TRUE.
!    if (nprocs>1) then ! Memory distribution over spins.
     if (Dtset%nsppol==2.and.iseven(nprocs)) then ! Distribute bands and spins.
       bks_mask=.FALSE.; keep_ur=.FALSE.
       do isp=1,my_nspins 
         spin = my_spins(isp)
         do band=1,mband 
           if (MODULO(band,nprocs/2)==my_rank-(spin-1)*(nprocs/2)) bks_mask(band,:,spin)=.TRUE.
         end do
       end do
     end if

     case default 
     MSG_ERROR("Wrong value for gwpara")
 end select
!
!Then each node owns the wavefunctions where GW corrections are required.
 do isp=1,my_nspins 
   spin = my_spins(isp)
   do ikcalc=1,Sigp%nkptgw
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Irred k-point for GW
     ii=Sigp%minbnd(ikcalc,spin); jj=Sigp%maxbnd(ikcalc,spin)
     bks_mask(ii:jj,ik_ibz,spin) = .TRUE.
     if (MODULO(Dtset%gwmem,10)==1) keep_ur(ii:jj,ik_ibz,spin)=.TRUE.
   end do
 end do

 opt_ecut=zero
 if (gw_uses_wfk_file) opt_ecut=Dtset%ecutwfn
!if (gw_uses_wfk_file) opt_ecut=Hdr_kss%ecut_eff*1.2
 opt_ecut=zero

 call wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,Sigp%npwwfn,mband,nband,Kmesh%nibz,Sigp%nsppol,bks_mask,&
& Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_kss%istwfk,Kmesh%ibz,gwc_ngfft,&
& Gsph_Max%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)

 ABI_DEALLOCATE(bks_mask)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(keep_ur)

 ABI_TIMER_STOP("init1")
 call timab(402,2,tsec) ! sigma(Init1)

 call timab(404,1,tsec) ! rdkss

 if (gw_uses_wfk_file) then
!  if (.FALSE..and.gw_uses_wfk_file) then
   MSG_WARNING("Wfd is init using WFK file")
   wfk_fname = Dtfil%fnameabi_kss; ii=LEN_TRIM(wfk_fname)
   wfk_fname(ii-2:ii) = "WFK"
#if defined HAVE_MPI_IO
   call wfd_read_wfk(Wfd,wfk_fname,IO_MODE_MPI)
#else
   call wfd_read_wfk(Wfd,wfk_fname,Dtset%accesswff)
#endif
 else 
   call wfd_read_kss(Wfd,Dtfil%fnameabi_kss,Sigp%nbnds,Dtset%accesswff,nelect)
 end if

 if (Dtset%pawcross==1) then
   call wfd_copy(Wfd,Wfdf)
   call wfd_change_ngfft(Wfdf,Cryst,Psps,ngfftf)
 end if

 call wfd_test_ortho(Wfd,Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

 call timab(404,2,tsec) ! rdkss

 ABI_TIMER_START("init2")
 call timab(405,1,tsec) ! Init2

 if (.FALSE.) then ! plot KSS wavefunctions. Change bks_mask to select particular states.
   ABI_ALLOCATE(bks_mask,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))
   bks_mask=.FALSE.
!  bks_mask(1:4,1,1)=.TRUE.
!  bks_mask=.TRUE.
   call wfd_plot_ur(Wfd,Cryst,Psps,Pawtab,Pawrad,ngfftf,bks_mask)
   ABI_DEALLOCATE(bks_mask)
 end if

!Debugging section.
 if (.FALSE.) then
!  
   if (.FALSE..and.Wfd%usepaw==1) then
     ABI_ALLOCATE(Cp1,(Wfd%natom,Wfd%nspinor))
     call cprj_alloc(Cp1,0,Wfd%nlmn_atm)
     ABI_ALLOCATE(Cp2,(Wfd%natom,Wfd%nspinor))
     call cprj_alloc(Cp2,0,Wfd%nlmn_atm)

     call wfd_change_ngfft(Wfd,Cryst,Psps,ngfftf) 

     do spin=1,Wfd%nsppol
       do ik_bz=1,Kmesh%nbz
         ik_ibz = Kmesh%tab(ik_bz)
         do band=1,Wfd%nband(ik_ibz,spin)
           call paw_check_symcprj(Wfd,ik_bz,band,spin,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cp1) 
           call paw_check_symcprj(Wfd,ik_bz,band,spin,2,Cryst,Kmesh,Psps,Pawtab,Pawang,Cp2) 

           do iat=1,Cryst%natom
             do isp=1,Wfd%nspinor
               write(789,'(3i2,/,(f8.4))') band,ik_bz,spin,Cp1(iat,isp)%cp
               write(790,'(3i2,/,(f8.4))') band,ik_bz,spin,Cp2(iat,isp)%cp
               write(791,'(3i2,/,(f8.4))') band,ik_bz,spin,Cp1(iat,isp)%cp(1,:)**2 + Cp1(iat,isp)%cp(2,:)**2
               write(792,'(3i2,/,(f8.4))') band,ik_bz,spin,Cp2(iat,isp)%cp(1,:)**2 + Cp2(iat,isp)%cp(2,:)**2
             end do
           end do
         end do
       end do
     end do

     call cprj_free(Cp1)
     ABI_DEALLOCATE(Cp1)
     call cprj_free(Cp2)
     ABI_DEALLOCATE(Cp2)
   end if

   call check_sym_ug(Wfd,Cryst,Kmesh,ierr) 
   ABI_CHECK(ierr==0,"Check sym_ug")
 end if
!
!==============================================================
!==== Find little group of the k-points for GW corrections ====
!==============================================================
!* The little group is calculated only if sys_sigma.
!* If use_umklp==1 then also symmetries requiring an umklapp to preserve k_gw are included.
!
 ABI_ALLOCATE(Ltg_k,(Sigp%nkptgw))
 call nullify_little_group(Ltg_k)
 use_umklp=1
 do ikcalc=1,Sigp%nkptgw
   if (Sigp%symsigma/=0) then
     call setup_little_group(Sigp%kptgw(:,ikcalc),Qmesh,Cryst,use_umklp,Ltg_k(ikcalc),0)
   end if
 end do
!
!=== Compute structure factor phases and large sphere cut-off ===
!WARNING cannot use Dtset%mgfft, this has to be checked better
!mgfft=MAXVAL(ngfftc(:))
!allocate(ph1d(2,3*(2*mgfft+1)*Cryst%natom),ph1df(2,3*(2*mgfftf+1)*Cryst%natom))
!write(std_out,*)' CHECK ',Dtset%mgfftdg,mgfftf
!if (Dtset%mgfftdg/=mgfftf) then
!write(std_out,*)"WARNING Dtset%mgfftf /= mgfftf"
!write(std_out,*)'HACKING Dtset%mgfftf'
!Dtset%mgfftdg=mgfftf
!end if
 ABI_ALLOCATE(ph1d,(2,3*(2*Dtset%mgfft+1)*Cryst%natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*Cryst%natom))

 call getph(Cryst%atindx,Cryst%natom,ngfftc(1),ngfftc(2),ngfftc(3),ph1d,Cryst%xred)

 if (Psps%usepaw==1.and.Pawfgr%usefinegrid==1) then
   call getph(Cryst%atindx,Cryst%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,Cryst%xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if
!
!===================================================================================
!==== Classify the GW wavefunctions according to the irreducible representation ====
!===================================================================================
!* Warning still under development.
!* Only for SCGW.
!bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw

 ABI_ALLOCATE(KS_sym,(Wfd%nkibz,Wfd%nsppol))
 call nullify_bands_symmetries(KS_sym)

 if (Sigp%symsigma==1.and.Sigp%gwcalctyp>=20) then
!  
!  call check_zarot(Sigp%npwvec,Cryst,gwc_ngfft,Gsph_Max%gvec,Psps,Pawang,Gsph_Max%rottb,Gsph_Max%rottbm1)

   use_paw_aeur=.FALSE. ! should pass ngfftf but the dense mesh is not forced to be symmetric
   do spin=1,Wfd%nsppol
     do ikcalc=1,Sigp%nkptgw
       ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
       first_band = Sigp%minbnd(ikcalc,spin)
       last_band  = Sigp%maxbnd(ikcalc,spin)
!      call classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,ngfftf,Cryst,KS_BSt,Pawtab,Pawrad,Pawang,Psps,&
       call classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,Wfd%ngfft,Cryst,KS_BSt,Pawtab,Pawrad,Pawang,Psps,&
&       Dtset%tolsym,KS_sym(ik_ibz,spin))
     end do
   end do
!  
!  Recreate the Sig_ij tables taking advantage of the classification of the bands.
   call sigma_tables(Sigp,Kmesh,KS_sym)
 end if

 ABI_TIMER_STOP("init2")
 call timab(405,2,tsec) ! Init2

 ABI_TIMER_START("make_vhxc")
 call timab(406,1,tsec) ! make_vhxc
!
!===========================
!=== COMPUTE THE DENSITY ===
!===========================
!* Evaluate the planewave part (complete charge in case of NC pseudos).

 ABI_ALLOCATE(ks_rhor,(nfftf,Dtset%nspden))
 ABI_ALLOCATE(ks_taur,(nfftf,Dtset%nspden*Dtset%usekden))
 call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_rhor)
 if(Dtset%usekden==1)call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_taur,optcalc=1)

!TODO this has to be done in a better way, moreover wont work for PAW
 call cutoff_density(ngfftf,Dtset%nspden,Dtset%nsppol,Vcp,ks_rhor)
!
!========================================
!==== Additional computation for PAW ====
!========================================
 if (Dtset%usepaw==1) then
!  
!  Calculate the compensation charge nhat.
   nhatgrdim=0; if (Dtset%xclevel==2) nhatgrdim=usexcnhat*Dtset%pawnhatxc
   cplex=1; ider=2*nhatgrdim; izero=0
   if (nhatgrdim>0)  then
     ABI_ALLOCATE(ks_nhatgr,(nfftf,Dtset%nspden,3*nhatgrdim))
   end if
   if (nhatgrdim==0)  then
     ABI_ALLOCATE(ks_nhatgr,(0,0,0))
   end if

   call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,Cryst%gprimd,MPI_enreg_seq,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Dtset%paral_kgb,Pawang,&
&   Pawfgrtab,ks_nhatgr,ks_nhat,KS_Pawrhoij,KS_Pawrhoij,Pawtab,k0,Cryst%rprimd,Cryst%ucvol,Cryst%xred)

   if (nhatgrdim==0)  then
     ABI_DEALLOCATE(ks_nhatgr)
   end if

!  === Evaluate onsite energies, potentials, densities ===
!  * Initialize variables/arrays related to the PAW spheres.
!  * Initialize also lmselect (index of non-zero LM-moments of densities).
!  TODO call init_paw_ij in scfcv, fix small issues
   ABI_ALLOCATE(KS_paw_ij,(Cryst%natom))
   call nullify_paw_ij(KS_paw_ij)
!  cplex=1; cplex_dij=Dtset%nspinor
   has_dijso=Dtset%pawspnorb; has_dijU=Dtset%usepawu

   call init_paw_ij(KS_paw_ij,cplex,cplex_dij,Dtset%nspinor,Dtset%nsppol,&
&   Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&   has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=1,has_dijxc_val=1,&
&   has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1)

   ABI_ALLOCATE(KS_paw_an,(Cryst%natom))
   call nullify_paw_an(KS_paw_an)
   nkxc1=0
   call init_paw_an(Cryst%natom,Cryst%ntypat,nkxc1,Dtset%nspden,cplex,Dtset%pawxcdev,&
&   Cryst%typat,Pawang,Pawtab,KS_paw_an,has_vxc=1,has_vxcval=1)
!  
!  Calculate onsite vxc with and without core charge.
   nzlmopt=-1; option=0; compch_sph=greatest_real
   call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,ipert,&
&   Dtset%ixc,MPI_enreg_seq,Cryst%natom,Cryst%natom,Dtset%nspden,&
&   Cryst%ntypat,nzlmopt,option,Dtset%paral_kgb,KS_Paw_an,KS_Paw_an,KS_paw_ij,&
&   Pawang,Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,&
&   Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,Dtset%xclevel,Dtset%xc_denpos,Psps%znuclpsp)

   write(std_out,*)" Silly write for XLF"
 end if !PAW

 if (.not.allocated(ks_nhatgr))  then
   ABI_ALLOCATE(ks_nhatgr,(nfftf,Dtset%nspden,0))
 end if

 call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,ks_rhor,Cryst%ucvol,&
& Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)
!
!For PAW, add the compensation charge on the FFT mesh, then get rho(G).
 if (Dtset%usepaw==1) ks_rhor=ks_rhor+ks_nhat

 call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_rhor,ucvol=ucvol)
 if(Dtset%usekden==1)call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_taur,optrhor=1,ucvol=ucvol)

 ABI_ALLOCATE(ks_rhog,(2,nfftf))
 ABI_ALLOCATE(ks_taug,(2,nfftf*Dtset%usekden))
 call fourdp(1,ks_rhog,ks_rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
 if(Dtset%usekden==1)call fourdp(1,ks_taug,ks_taur(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)

!
!The following steps have been gathered in the setvtr routine:
!- get Ewald energy and Ewald forces
!- compute local ionic pseudopotential vpsp
!- eventually compute 3D core electron density xccc3d
!- eventually compute vxc and vhartr
!- set up ks_vtrial
!
!*******************************************************************
!**** NOTE THAT HERE Vxc CONTAINS THE CORE-DENSITY CONTRIBUTION ****
!*******************************************************************

 ABI_ALLOCATE(grewtn,(3,Cryst%natom))
 ABI_ALLOCATE(xred_dummy,(3,Cryst%natom))
 xred_dummy(:,:)=xred(:,:)
 nkxc=0
 if (Dtset%nspden==1) nkxc=2
 if (Dtset%nspden>=2) nkxc=3 ! check GGA and spinor, quite a messy part!!!
#if defined HAVE_DFT_LIBXC
 if (Dtset%ixc<0 .and. libxc_functionals_ismgga()) nkxc=0 ! in case of MGGA, fxc and kxc are not available
!and we dont need them for the screening part (for now ...)
#endif
 if (nkxc/=0)  then
   ABI_ALLOCATE(kxc,(nfftf,nkxc))
 end if

 n3xccc=0; if (Psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))
 ABI_ALLOCATE(ks_vhartr,(nfftf))
 ABI_ALLOCATE(ks_vtrial,(nfftf,Dtset%nspden))
 ABI_ALLOCATE(vpsp,(nfftf))
 ABI_ALLOCATE(ks_vxc,(nfftf,Dtset%nspden))

 optene=4; moved_atm_inside=0; moved_rhor=0; istep=1

 call setvtr(Cryst%atindx1,Dtset,KS_energies,gmet,gprimd,grewtn,gsqcutf_eff,&
& istep,kxc,mgfftf,moved_atm_inside,moved_rhor,MPI_enreg_seq,&
& Cryst%nattyp,nfftf,ngfftf,ks_nhat,ks_nhatgr,nhatgrdim,nkxc,Cryst%ntypat,Psps%n1xccc,n3xccc,&
& optene,Pawtab,ph1df,Psps,ks_rhog,ks_rhor,Cryst%rmet,Cryst%rprimd,strsxc,&
& Cryst%ucvol,usexcnhat,ks_vhartr,vpsp,ks_vtrial,ks_vxc,vxcavg,Wvl%descr,xccc3d,xred_dummy,taug=ks_taug,taur=ks_taur)
!TODO here xred is INOUT due to ionion_realSpace and xredcart!

!============================
!==== Compute KS PAW Dij ====
!============================
 if (Dtset%usepaw==1) then
   write(std_out,*)"Another silly write for XLF"
!  
!  Calculate the unsymmetrized Dij. 
   ipert=0; idir=0
   call pawdij(cplex,Dtset,Dtset%enunit,one,Cryst%gprimd,ipert,MPI_enreg_seq,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,Dtset%nspden,Cryst%ntypat,&
&   Dtset%paral_kgb,KS_paw_an,KS_paw_ij,Pawang,Pawfgrtab,&
&   Dtset%pawprtvol,Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
&   k0,Cryst%typat,Cryst%ucvol,ks_vtrial,ks_vxc,Cryst%xred)
!  
!  Symmetrize KS Dij 
#if 0
   call symdij(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,&
&   Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,0,KS_paw_ij,Pawang,&
&   Dtset%pawprtvol,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
#else
   call symdij_all(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,&
&   Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,KS_paw_ij,Pawang,&
&   Dtset%pawprtvol,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
#endif

!  
!  Output the pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   call pawprt(Dtset,Psps%indlmn,Psps%lmnmax,KS_paw_ij,KS_Pawrhoij,Pawtab)
 end if

 ABI_TIMER_STOP("make_vhxc")
 call timab(406,2,tsec) ! make_vhxc
!
!=== Calculate Vxc(b1,b2,k,s)=<b1,k,s|v_{xc}|b2,k,s>  for all the states included in GW ===
!* ks_vxcvalme is calculated without NLCC, ks_vxcme contains NLCC (if any)
!* This part is parallelized within MPI_COMM_WORD since each node has all GW wavefunctions.
!* ks_vUme is zero unless we are using LDA+U as starting point, see calc_vHxc_braket
!* Note that vH matrix elements are calculated using the true uncutted interaction.
 ABI_TIMER_START("vhxc_me")
 call timab(407,1,tsec) ! vHxc_me

 call reset_mflags(KS_mflags)
 KS_mflags%has_vhartree=1
 KS_mflags%has_vxc     =1
 KS_mflags%has_vxcval  =1
 if (Dtset%usepawu>0     )  KS_mflags%has_vu     =1
 if (Dtset%useexexch>0   )  KS_mflags%has_lexexch=1
 if (Sigp%use_sigxcore==1)  KS_mflags%has_sxcore =1
 if (Sigp%gwcalctyp<10   )  KS_mflags%only_diago =1 ! off-diagonal elements only for SC on wavefunctions.

 if (.FALSE.) then ! quick and dirty hack to test HF contribution.
   MSG_WARNING("testing on-site HF")
   lmn2_size_max=MAXVAL(Pawtab(:)%lmn2_size)
   ABI_ALLOCATE(dij_hf,(cplex_dij*lmn2_size_max,ndij,Cryst%natom))
   call paw_dijhf(ndij,cplex_dij,lmn2_size_max,Cryst,Psps,Pawtab,Pawrad,Pawang,KS_Pawrhoij,dij_hf,Dtset%prtvol)

   do iat=1,Cryst%natom
     itypat = Cryst%typat(iat)
     ii = Pawtab(itypat)%lmn2_size
     KS_Paw_ij(iat)%dijxc(:,:) = dij_hf(1:cplex_dij*ii,:,iat)
   end do
   ABI_DEALLOCATE(dij_hf)

!  option_dij=3
!  call symdij(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,&
!  &   Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
!  &   KS_paw_ij,Pawang,Dtset%pawprtvol,Cryst%rprimd,&
!  &   Cryst%symafm,Cryst%symrec,Cryst%typat)
 end if

 ABI_ALLOCATE(tmp_kstab,(2,Wfd%nkibz,Wfd%nsppol))
 tmp_kstab=0
 do spin=1,Sigp%nsppol
   do ikcalc=1,Sigp%nkptgw ! No spin dependent!
     ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
     tmp_kstab(1,ik_ibz,spin)=Sigp%minbnd(ikcalc,spin)
     tmp_kstab(2,ik_ibz,spin)=Sigp%maxbnd(ikcalc,spin)
   end do
 end do

 call calc_vhxc_me(Wfd,KS_mflags,KS_me,Cryst,Dtset,gsqcutf_eff,nfftf,ngfftf,&
& ks_vtrial,ks_vhartr,ks_vxc,Psps,Pawtab,KS_paw_an,Pawang,Pawfgrtab,KS_paw_ij,dijexc_core,&
& ks_rhor,ks_rhog,usexcnhat,ks_nhat,ks_nhatgr,nhatgrdim,tmp_kstab,taug=ks_taug,taur=ks_taur)
 ABI_DEALLOCATE(tmp_kstab)

!#ifdef DEV_HAVE_SCGW_SYM
!Set KS matrix elements connecting different irreps to zero. Do not touch unknown bands!.
 if (Sigp%gwcalctyp>=20 .and. Sigp%symsigma > 0) then 
   bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw
   ABI_ALLOCATE(ks_irreptab,(bmin:bmax,Kmesh%nibz,Sigp%nsppol))
   ks_irreptab=0
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw
       ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
       first_band = Sigp%minbnd(ikcalc,spin)
       last_band  = Sigp%maxbnd(ikcalc,spin)
       if (.not.bsym_failed(KS_sym(ik_ibz,spin))) then
         ks_irreptab(first_band:last_band,ik_ibz,spin) = KS_sym(ik_ibz,spin)%b2irrep(first_band:last_band)
!        ks_irreptab(bmin:bmax,ik_ibz,spin) = KS_sym(ik_ibz,spin)%b2irrep(bmin:bmax)
       end if
     end do
   end do
   call zero_melements(KS_me,ks_irreptab)
   ABI_DEALLOCATE(ks_irreptab)
 end if
!#endif

 call print_melements(KS_me,header="Matrix elements in the KS basis set",prtvol=Dtset%prtvol)
!
!If possible, calculate the EXX energy from the between the frozen core
!and the valence electrons using KS wavefunctions
!MG: BE careful here, since exchange_energy is meaningful only is all occupied states are calculated.
 if( KS_mflags%has_sxcore ==1 ) then
   exchange_energy=zero
   do spin=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       do ib=b1gw,b2gw
         if (Sigp%nsig_ab==1) then
           exchange_energy = exchange_energy + half*KS_BSt%occ(ib,ik,spin)*Kmesh%wt(ik)*KS_me%sxcore(ib,ib,ik,spin)
         else
           exchange_energy = exchange_energy + half*KS_BSt%occ(ib,ik,spin)*Kmesh%wt(ik)*SUM(KS_me%sxcore(ib,ib,ik,:))
         end if
       end do
     end do
   end do
   write(msg,'(a,2(es16.6,a))')' CORE Exchange energy with KS wavefunctions: ',exchange_energy,' Ha ,',exchange_energy*Ha_eV,' eV'
   call wrtout(std_out,msg,'COLL')
!  call wrtout(ab_out,msg,'COLL')
 end if

 ABI_TIMER_STOP("vhxc_me")
 call timab(407,2,tsec) ! vHxc_me

 ABI_TIMER_START("hqp_init")
 call timab(408,1,tsec) ! hqp_init

!Do not break this coding! When gwcalctyp>10, the order of the bands can be interexchanged after
!the diagonalization. Therefore, we have to correctly assign the matrix elements to the corresponding
!bands and we cannot skip the following even though it looks unuseful.
 if (Sigp%gwcalctyp>=10) then
   call wrtout(std_out,ch10//' *************** KS Energies *******************','COLL')
 end if

!=== QP_BSt stores energies and occ. used for the calculation ===
!* Initialize QP_BSt with KS values.
!* In case of SC update QP_BSt using the QPS file.
 call copy_bandstructure(KS_BSt,QP_BSt)

 ABI_ALLOCATE(qp_rhor,(nfftf,Dtset%nspden))
 ABI_ALLOCATE(qp_taur,(nfftf,Dtset%nspden*Dtset%usekden))
 QP_sym => KS_sym

 if (Sigp%gwcalctyp<10) then  ! one-shot GW, just do a copy of the KS density.
   qp_rhor=ks_rhor
   if(Dtset%usekden==1)qp_taur=ks_taur
   QP_sym => KS_sym
 else                         ! Self-consistent GW.
!  
!  * Read the unitary matrix and the QP energies of the previous step from the QPS file.
   call energies_init(QP_energies)
   QP_energies%e_corepsp=ecore/Cryst%ucvol

!  $ m_lda_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}> $
   ABI_ALLOCATE(m_lda_to_qp,(Sigp%nbnds,Sigp%nbnds,Kmesh%nibz,Sigp%nsppol))
   m_lda_to_qp=czero
   do ib=1,Sigp%nbnds
     m_lda_to_qp(ib,ib,:,:)=cone ! Initialize the QP amplitudes with KS wavefunctions.
   end do

!  * Now read m_lda_to_qp and update the energies in QP_BSt.
!  TODO switch on the renormalization of n in sigma.
   ABI_ALLOCATE(prev_rhor,(nfftf,Dtset%nspden))
   ABI_ALLOCATE(prev_taur,(nfftf,Dtset%nspden*Dtset%usekden))
   ABI_ALLOCATE(prev_Pawrhoij,(Cryst%natom*Psps%usepaw))

   call rdqps(QP_BSt,Dtfil%fnameabi_qps,Dtset%usepaw,Dtset%nspden,1,nscf,&
&   nfftf,ngfftf,Cryst%ucvol,Dtset%paral_kgb,Cryst,Pawtab,MPI_enreg_seq,nbsc,m_lda_to_qp,prev_rhor,prev_Pawrhoij)

!  Find the irreps associated to the QP amplitudes starting from the analogous table for the KS states.
!  bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw
!  allocate(qp_irreptab(bmin:bmax,Kmesh%nibz,Sigp%nsppol))
!  qp_irreptab=0
!  !qp_irreptab=ks_irreptab

!  do jb_qp=bmin,bmax
!  do ib_ks=bmin,bmax
!  if (ABS(m_lda_to_qp(ib_ks,jb_qp,ik_ibz,spin)) > tol12) then ! jb_qp has same the same character as ib_ks.
!  ks_irr = ks_irreptab(ib_ks,ib_ks,ik_ibz,spin)
!  qp_irreptab(jb_qp,jb_qp,ik_ibz,spin) = ks_irr
!  do ii=bmin,bmax
!  if (ks_irr == ks_irreptab(ii,ib_ks,ik_ibz,spin)) then
!  qp_irreptab(jb_qp,ii,ik_ibz,spin) = ks_irr
!  end if
!  end do
!  end if
!  end do
!  end do

   if (nscf==0) prev_rhor=ks_rhor
   if (nscf==0 .and. Dtset%usekden==1) prev_taur=ks_taur

   if (nscf>0.and.Sigp%gwcalctyp>=20.and.wfd_iam_master(Wfd)) then ! Print the unitary transformation on std_out.
     call show_QP(QP_BSt,m_lda_to_qp,fromb=Sigp%minbdgw,tob=Sigp%maxbdgw,unit=std_out,tolmat=0.001_dp)
   end if
!  
!  === Compute QP wfg as linear combination of KS states ===
!  * Wfd%ug is modified inside calc_wf_qp
!  * For PAW, update also the on-site projections.
!  * WARNING the first dimension of MPI_enreg MUST be Kmesh%nibz
!  TODO here we should use nbsc instead of nbnds

   call wfd_rotate(Wfd,Cryst,m_lda_to_qp)
!  === Reinit the storage mode of Wfd as ug have been changed ===
!  * Update also the wavefunctions for GW corrections on each processor
   call wfd_reset_ur_cprj(Wfd)
!  
!  Compute QP occupation numbers.
   write(msg,'(3a)')ch10,' sigma : calculating QP occupation numbers ',ch10
   call wrtout(std_out,msg,'COLL')

   call update_occ(QP_BSt,Dtset%fixmom,prtvol=0)
   qp_vbik(:,:) = get_valence_idx(QP_BSt)

!  #ifdef DEV_HAVE_SCGW_SYM
!  Calculate the irreducible representations of the new QP amplitdues.
   if (Sigp%symsigma==1.and.Sigp%gwcalctyp>=20) then

     ABI_ALLOCATE(QP_sym,(Wfd%nkibz,Wfd%nsppol))
     call nullify_bands_symmetries(QP_sym)
     use_paw_aeur=.FALSE. ! should pass ngfftf but the dense mesh is not forced to be symmetric
     do spin=1,Wfd%nsppol
       do ikcalc=1,Sigp%nkptgw
         ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
!        Quick fix for SCGW+symm TODO fix properly!
         first_band = Sigp%minbnd(ikcalc,spin)
         last_band  = Sigp%maxbnd(ikcalc,spin)
!        first_band = MINVAL(Sigp%minbnd(:,spin))
!        last_band  = MAXVAL(Sigp%maxbnd(:,spin))
!        call classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,ngfftf,Cryst,QP_BSt,Pawtab,Pawrad,Pawang,Psps,&
         call classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,Wfd%ngfft,Cryst,QP_BSt,Pawtab,Pawrad,Pawang,Psps,&
&         Dtset%tolsym,QP_sym(ik_ibz,spin))
       end do
     end do
!    
!    Recreate the Sig_ij tables taking advantage of the classification of the bands.
     call sigma_tables(Sigp,Kmesh,QP_sym)
   end if
!  #endif
!  
!  Compute QP density using the updated wfg.
   call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_rhor)
   if(Dtset%usekden==1) call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_taur,optcalc=1)
!  
!  ========================================
!  ==== QP self-consistent GW with PAW ====
!  ========================================
   if (Dtset%usepaw==1) then

     ABI_ALLOCATE(qp_nhat,(nfftf,Dtset%nspden))
     nhatgrdim=0; if (Dtset%xclevel==2) nhatgrdim=usexcnhat
     ABI_ALLOCATE(qp_nhatgr,(nfftf,Dtset%nspden,3*nhatgrdim))

     ABI_ALLOCATE(QP_pawrhoij,(Cryst%natom))
     ABI_ALLOCATE(QP_paw_ij,(Cryst%natom))
     ABI_ALLOCATE(QP_paw_an,(Cryst%natom))
!    
!    Calculate new QP quantities: nhat, nhatgr, rho_ij, paw_ij, and paw_an.
     call paw_qpscgw(Wfd,nscf,nfftf,ngfftf,Dtset,Cryst,Kmesh,Psps,QP_BSt,&
&     Pawang,Pawrad,Pawtab,Pawfgrtab,prev_Pawrhoij,&
&     QP_pawrhoij,QP_paw_ij,QP_paw_an,QP_energies,qp_nhat,nhatgrdim,qp_nhatgr,compch_sph,compch_fft,MPI_enreg_seq)
   end if

!  Allocate these arrays anyway, since they are passed to subroutines.
   if (.not.allocated(qp_nhat  ))  then
     ABI_ALLOCATE(qp_nhat,(nfftf,0))
   end if
   if (.not.allocated(qp_nhatgr))  then
     ABI_ALLOCATE(qp_nhatgr,(nfftf,Dtset%nspden,0))
   end if

!  here I should renormalize the density
   call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,qp_rhor,Cryst%ucvol,&
&   Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)

   if (Dtset%usepaw==1) qp_rhor(:,:)=qp_rhor(:,:)+qp_nhat(:,:) ! Add the "hat" term.

   call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,qp_rhor,ucvol=ucvol)
   if(Dtset%usekden==1)call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,qp_taur,optrhor=1,ucvol=ucvol)

!  
!  Simple mixing of the PW density to damp oscillations in the Hartree potential.
   if (nscf>0 .and. (ABS(Dtset%rhoqpmix-one)>tol12) ) then
     write(msg,'(2a,f6.3)')ch10,' sigma: mixing QP densities using rhoqpmix= ',Dtset%rhoqpmix
     call wrtout(std_out,msg,'COLL')
     qp_rhor = prev_rhor + Dtset%rhoqpmix*(qp_rhor-prev_rhor)
     if(Dtset%usekden==1)qp_taur = prev_taur + Dtset%rhoqpmix*(qp_taur-prev_taur) ! this is this line that I'm not if it is useful
   end if

   ABI_DEALLOCATE(prev_rhor)
   ABI_DEALLOCATE(prev_taur)
   if (Psps%usepaw==1.and.nscf>0) then
     call rhoij_free(prev_pawrhoij)
   end if
   ABI_DEALLOCATE(prev_pawrhoij)
   istat = ABI_ALLOC_STAT

   ABI_ALLOCATE(qp_rhog,(2,nfftf))
   ABI_ALLOCATE(qp_taug,(2,nfftf*Dtset%usekden))
   call fourdp(1,qp_rhog,qp_rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
   if(Dtset%usekden==1)call fourdp(1,qp_taug,qp_taur(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
!  
!  ===========================================
!  ==== Optional output of the QP density ====
!  ===========================================
   if (Dtset%prtden/=0.and.wfd_iam_master(Wfd)) then
     rdwr=2; fformr=52; rdwrpaw=0
     call ioarr(accessfil,qp_rhor,Dtset,dummy,fformr,Dtfil%fnameabo_qp_den,Hdr_sigma,&
&     MPI_enreg_seq,nfftf,Pawrhoij_dum,rdwr,rdwrpaw,Wvl%descr)
     if (accessfil==3) then ! Complete the missing geometry and electronic information.
       call abi_etsf_geo_put(Dtset,Dtfil%fnameabo_qp_den,Psps)
       call abi_etsf_electrons_put(Dtset,Dtfil%fnameabo_qp_den)
     end if
   end if
!  
!  ===========================================
!  === Optional output of the full QP density
!  ===========================================
   if (Wfd%usepaw==1.and.Dtset%pawprtden==1) then
     call xbarrier_mpi(comm)  ! First make sure that all processors are here
     ABI_ALLOCATE(qp_rhor_paw   ,(nfftf,Wfd%nspden))
     ABI_ALLOCATE(qp_rhor_n_one ,(nfftf,Wfd%nspden))
     ABI_ALLOCATE(qp_rhor_nt_one,(nfftf,Wfd%nspden))

     call denfgr(Cryst%atindx1,Cryst%gmet,comm,Cryst%natom,Cryst%nattyp,ngfftf,qp_nhat,&
&     Wfd%nspinor,Wfd%nsppol,Wfd%nspden,Cryst%ntypat,Pawfgr,Pawrad,QP_pawrhoij,Pawtab,Dtset%prtvol,&
&     Psps,qp_rhor,qp_rhor_paw,qp_rhor_n_one,qp_rhor_nt_one,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

     ABI_DEALLOCATE(qp_rhor_n_one)
     ABI_DEALLOCATE(qp_rhor_nt_one)
     if (Dtset%prtvol>9) then ! Print a normalisation check
       norm = SUM(qp_rhor_paw(:,1))*Cryst%ucvol/PRODUCT(Pawfgr%ngfft(1:3))
       write(msg,'(a,F8.4)') '  QUASIPARTICLE DENSITY CALCULATED - NORM OF DENSITY: ',norm
       call wrtout(std_out,msg,'PERS')
     end if

     if (my_rank==master) then ! Write the density to file
       rdwr=2; fformr=52; rdwrpaw=0
       call ioarr(accessfil,qp_rhor_paw,Dtset,dummy,fformr,Dtfil%fnameabo_qp_pawden,&
&       Hdr_sigma,MPI_enreg_seq,nfftf,Pawrhoij_dum,rdwr,rdwrpaw,Wvl%descr)
       if ( accessfil == 3 ) then
!        Complete the geometry informations with missing values from hdr_io().
         call abi_etsf_geo_put(Dtset,Dtfil%fnameabo_qp_pawden, Psps)
!        Complete the electrons definition with missing values from hdr_io().
         call abi_etsf_electrons_put(Dtset,Dtfil%fnameabo_qp_pawden)
       end if
     end if
     ABI_DEALLOCATE(qp_rhor_paw)
   end if
!  
!  ===========================================
!  === Optional output of the QP amplitudes
!  ===========================================
   if (.FALSE.) then ! plot QP amplitudes. Change bks_mask to select particular states.
     ABI_ALLOCATE(bks_mask,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))
     bks_mask=.FALSE.
!    bks_mask(1,1,1)=.TRUE.
     call wfd_plot_ur(Wfd,Cryst,Psps,Pawtab,Pawrad,ngfftf,bks_mask)
     ABI_DEALLOCATE(bks_mask)
   end if

   nkxc=0
   if (Dtset%nspden==1) nkxc=2
   if (Dtset%nspden>=2) nkxc=3 !check GGA and spinor that is messy !!!
#if defined HAVE_DFT_LIBXC
   if (Dtset%ixc<0 .and. libxc_functionals_ismgga()) nkxc=0 ! in case of MGGA, fxc and kxc are not available
!  and we dont need them for the screening part (for now ...)
#endif
   if (nkxc/=0)  then
     ABI_ALLOCATE(qp_kxc,(nfftf,nkxc))
   end if
!  
!  **** NOTE THAT Vxc CONTAINS THE CORE-DENSITY CONTRIBUTION ****
!  TODO here xred is INOUT due to ionion_realSpace and xredcart, why?
   n3xccc=0; if (Psps%n1xccc/=0) n3xccc=nfftf
   ABI_ALLOCATE(qp_vhartr,(nfftf))
   ABI_ALLOCATE(qp_vtrial,(nfftf,Dtset%nspden))
   ABI_ALLOCATE(qp_vxc,(nfftf,Dtset%nspden))

   optene=4; moved_atm_inside=0; moved_rhor=0; istep=1

   call setvtr(Cryst%atindx1,Dtset,QP_energies,gmet,gprimd,grewtn,gsqcutf_eff,&
&   istep,qp_kxc,mgfftf,moved_atm_inside,moved_rhor,MPI_enreg_seq,&
&   Cryst%nattyp,nfftf,ngfftf,qp_nhat,qp_nhatgr,nhatgrdim,nkxc,Cryst%ntypat,Psps%n1xccc,n3xccc,&
&   optene,Pawtab,ph1df,Psps,qp_rhog,qp_rhor,Cryst%rmet,Cryst%rprimd,strsxc,&
&   Cryst%ucvol,usexcnhat,qp_vhartr,vpsp,qp_vtrial,qp_vxc,vxcavg_qp,Wvl%descr,xccc3d,xred_dummy,taug=qp_taug,taur=qp_taur)

   if (allocated(qp_kxc))  then
     ABI_DEALLOCATE(qp_kxc)
   end if

   if (Dtset%usepaw==1) then
!    
!    === Compute QP Dij ===
     ipert=0; idir=0
     call pawdij(cplex,Dtset,Dtset%enunit,one,Cryst%gprimd,ipert,MPI_enreg_seq,&
&     Cryst%natom,Cryst%natom,nfftf,ngfftf,Dtset%nspden,Cryst%ntypat,&
&     Dtset%paral_kgb,QP_paw_an,QP_paw_ij,Pawang,Pawfgrtab,&
&     Dtset%pawprtvol,Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
&     k0,Cryst%typat,Cryst%ucvol,qp_vtrial,qp_vxc,Cryst%xred)
!    
!    === Symmetrize total Dij ===
     option_dij=0

#if 0
     call symdij(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,&
&     Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
&     QP_paw_ij,Pawang,Dtset%pawprtvol,Cryst%rprimd,&
&     Cryst%symafm,Cryst%symrec,Cryst%typat)
#else
     call symdij_all(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,&
&     Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
&     QP_paw_ij,Pawang,Dtset%pawprtvol,Cryst%rprimd,&
&     Cryst%symafm,Cryst%symrec,Cryst%typat)
#endif

!    
!    Output the QP pseudopotential strengths Dij and the augmentation occupancies Rhoij.
     call pawprt(Dtset,Psps%indlmn,Psps%lmnmax,QP_paw_ij,QP_Pawrhoij,Pawtab)
   end if

   ehartree=half*SUM(qp_rhor(:,1)*qp_vhartr(:))/DBLE(nfftf)*Cryst%ucvol

   write(msg,'(a,80a)')ch10,('-',ii=1,80)
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(5a,f9.4,3a,es21.14,2a,es21.14)')ch10,&
&   ' QP results after the unitary transformation in the KS subspace: ',ch10,ch10,&
&   '  Number of electrons    = ',qp_rhog(1,1)*Cryst%ucvol,ch10,ch10,&
&   '  QP Band energy    [Ha] = ',get_bandenergy(QP_BSt),ch10,&
&   '  QP Hartree energy [Ha] = ',ehartree
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,80a)')ch10,('-',ii=1,80)
   call wrtout(ab_out,msg,'COLL')
!  
!  TODO Since plasmonpole model 2-3-4 depend on the Fourier components of the density
!  in case of self-consistency we might calculate here the ppm coefficients using qp_rhor
 end if ! gwcalctyp>=10
!
!=== KS hamiltonian hlda(b1,b1,k,s)= <b1,k,s|H_s|b1,k,s> ===
 ABI_ALLOCATE(hlda,(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sigp%nsppol*Sigp%nsig_ab))
 hlda=czero

 if (Dtset%nspinor==1) then
   do spin=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       do ib=b1gw,b2gw
         hlda(ib,ib,ik,spin) = KS_BSt%eig(ib,ik,spin)
       end do
     end do
   end do
 else ! Spinorial case
!  * Note that here vxc contains the contribution of the core.
!  * Scale ovlp if orthonormalization is not satisfied as npwwfn might be < npwvec.
!  TODO add spin-orbit case
   if (Wfd%usepaw==1) then
     ABI_ALLOCATE(Cp1,(Wfd%natom,Wfd%nspinor))
     call cprj_alloc(Cp1,0,Wfd%nlmn_atm)
   end if

   do spin=1,Sigp%nsppol
     do ik_ibz=1,Kmesh%nibz
       do ib=b1gw,b2gw

         ug1  => Wfd%Wave(ib,ik_ibz,spin)%ug
         cdummy = xdotc(Wfd%npwwfn*Wfd%nspinor,ug1,1,ug1,1)
         ovlp(1) = REAL(cdummy) 
         ovlp(2) = AIMAG(cdummy)

         if (Psps%usepaw==1) then
           call wfd_get_cprj(Wfd,ib,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
           ovlp = ovlp + paw_overlap(Cp1,Cp1,Cryst%typat,Pawtab)
         end if
!        write(std_out,*)ovlp(1),ovlp(2)
         norm=DBLE(ovlp(1)+ovlp(2))
         ovlp(1)=DBLE(ovlp(1)/norm)
         ovlp(2)=DBLE(ovlp(2)/norm)
!        ovlp(2)=cone-ovlp(1)
         hlda(ib,ib,ik_ibz,1) = KS_BSt%eig(ib,ik_ibz,1)*ovlp(1)-KS_me%vxc(ib,ib,ik_ibz,3)
         hlda(ib,ib,ik_ibz,2) = KS_BSt%eig(ib,ik_ibz,1)*ovlp(2)-KS_me%vxc(ib,ib,ik_ibz,4)
         hlda(ib,ib,ik_ibz,3) = KS_me%vxc(ib,ib,ik_ibz,3)
         hlda(ib,ib,ik_ibz,4) = KS_me%vxc(ib,ib,ik_ibz,4)
       end do
     end do
   end do

   if (Wfd%usepaw==1) then
     call cprj_free(Cp1)
     ABI_DEALLOCATE(Cp1)
   end if
 end if
!
!=== Initialize Sigma results ===
!TODO it is better if we use ragged arrays indexed by the k-point
 call init_sigma_results(Sigp,Kmesh%nibz,Dtset%usepawu,Sr)
!
!=== Setup of the bare Hamiltonian := T + v_{loc} + v_{nl} + v_H ===
!* The representation depends wheter we are updating the wfs or not.
!* ks_vUme is zero unless we are using LDA+U as starting point, see calc_vHxc_braket
!* Note that vH matrix elements are calculated using the true uncutted interaction.

 if (Sigp%gwcalctyp<10) then  ! * For one-shot GW use the KS representation.
   Sr%hhartree=hlda-KS_me%vxcval
!  === Additional goodies for PAW ===
!  * LDA +U Hamiltonian
!  * LEXX.
!  * Core contribution estimated using Fock exchange.
   if (Dtset%usepaw==1) then
     if (Sigp%use_sigxcore==1) Sr%hhartree=hlda - (KS_me%vxc - KS_me%sxcore)
     if (Dtset%usepawu>0) Sr%hhartree=Sr%hhartree-KS_me%vu
     if (Dtset%useexexch>0) then
       MSG_ERROR("useexexch > 0 not implemented")
       Sr%hhartree = Sr%hhartree - KS_me%vlexx
     end if
   end if
 else !  Self-consistent on energies and|or wavefunctions.
!  * For NC get the bare Hamiltonian  $H_{bare}= T+v_{loc}+ v_{nl}$ in the KS representation
!  * For PAW, calculate the matrix elements of h0, store also the new Dij in QP_Paw_ij.
!  * h0 is defined as T+vH[tn+nhat+tnZc] + vxc[tnc] + dij_eff and
!  dij_eff = dij^0 + dij^hartree + dij^xc-dij^xc_val + dijhat - dijhat_val.
!  In the above expression tn, tnhat are QP quantities.
   if (Dtset%usepaw==0) then
     ABI_ALLOCATE(hbare,(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sigp%nsppol*Sigp%nsig_ab))
     hbare=hlda-KS_me%vhartree-KS_me%vxcval
!    
!    * Change basis from KS to QP, hbare is overwritten: A_{QP} = U^\dagger A_{KS} U
     ABI_ALLOCATE(htmp,(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sigp%nsppol*Sigp%nsig_ab))
     ABI_ALLOCATE(ctmp,(b1gw:b2gw,b1gw:b2gw))
     ABI_ALLOCATE(uks2qp,(b1gw:b2gw,b1gw:b2gw))
     htmp=hbare; hbare=czero

     do spin=1,Sigp%nsppol
       do ik=1,Kmesh%nibz
         uks2qp(:,:) = m_lda_to_qp(b1gw:b2gw,b1gw:b2gw,ik,spin)
         do iab=1,Sigp%nsig_ab
           is_idx=spin; if (Sigp%nsig_ab>1) is_idx=iab
           ctmp(:,:)=MATMUL(htmp(:,:,ik,is_idx),uks2qp(:,:))
           hbare(:,:,ik,is_idx)=MATMUL(TRANSPOSE(CONJG(uks2qp)),ctmp)
         end do
       end do
     end do
     ABI_DEALLOCATE(htmp)
     ABI_DEALLOCATE(ctmp)
     ABI_DEALLOCATE(uks2qp)
   end if ! usepaw==0

!  === Calculate the QP matrix elements===
!  * This part is parallelized within MPI_COMM_WORD since each node has all GW wavefunctions.
!  * For PAW, we have to construct the new bare Hamiltonian.
   call wrtout(std_out,ch10//' *************** QP Energies *******************','COLL')

   call reset_mflags(QP_mflags)
!  if (Sigp%gwcalctyp<20) QP_mflags%only_diago=1 ! For e-only, no need of off-diagonal elements.
   QP_mflags%has_vhartree=1
   if (Dtset%usepaw==1)    QP_mflags%has_hbare  =1
!  QP_mflags%has_vxc     =1
!  QP_mflags%has_vxcval  =1
!  if (Sigp%use_sigxcore==1) QP_mflags%has_sxcore =1
!  if (Dtset%usepawu>0)    QP_mflags%has_vu     =1
!  if (Dtset%useexexch>0)  QP_mflags%has_lexexch=1

   ABI_ALLOCATE(tmp_kstab,(2,Wfd%nkibz,Wfd%nsppol))
   tmp_kstab=0
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw ! No spin dependent!
       ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
       tmp_kstab(1,ik_ibz,spin)=Sigp%minbnd(ikcalc,spin)
       tmp_kstab(2,ik_ibz,spin)=Sigp%maxbnd(ikcalc,spin)
     end do
   end do

   call calc_vhxc_me(Wfd,QP_mflags,QP_me,Cryst,Dtset,gsqcutf_eff,nfftf,ngfftf,&
&   qp_vtrial,qp_vhartr,qp_vxc,Psps,Pawtab,QP_paw_an,Pawang,Pawfgrtab,QP_paw_ij,dijexc_core,&
&   qp_rhor,qp_rhog,usexcnhat,qp_nhat,qp_nhatgr,nhatgrdim,tmp_kstab,taug=qp_taug,taur=qp_taur)
   ABI_DEALLOCATE(tmp_kstab)

!  #ifdef DEV_HAVE_SCGW_SYM
   if (Sigp%gwcalctyp>=20 .and. Sigp%symsigma>0) then 
     bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw
     ABI_ALLOCATE(qp_irreptab,(bmin:bmax,Kmesh%nibz,Sigp%nsppol))
     qp_irreptab=0
     do spin=1,Sigp%nsppol
       do ikcalc=1,Sigp%nkptgw
         ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
         first_band = Sigp%minbnd(ikcalc,spin)
         last_band  = Sigp%maxbnd(ikcalc,spin)
         if (.not.bsym_failed(QP_sym(ik_ibz,spin))) then
           qp_irreptab(first_band:last_band,ik_ibz,spin) = QP_sym(ik_ibz,spin)%b2irrep(first_band:last_band)
!          qp_irreptab(bmin:bmax,ik_ibz,spin) = QP_sym(ik_ibz,spin)%b2irrep(bmin:bmax)
         end if
       end do
     end do
     call zero_melements(QP_me,qp_irreptab)
     ABI_DEALLOCATE(qp_irreptab)
   end if
!  #endif

   call print_melements(QP_me,header="Matrix elements in the QP basis set",prtvol=Dtset%prtvol)
!  
!  * Output the QP pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   if (Dtset%usepaw==1) then
     call wrtout(std_out," *** After calc_vHxc_braket *** ",'COLL')
     call print_paw_ij(QP_Paw_ij,std_out,Dtset%pawprtvol,"COLL") ! TODO terminate the implementation of this routine.
     call pawprt(Dtset,Psps%indlmn,Psps%lmnmax,QP_paw_ij,QP_Pawrhoij,Pawtab)
   end if

   if (Dtset%usepaw==0) then
     Sr%hhartree=hbare+QP_me%vhartree
   else
     Sr%hhartree=QP_me%hbare
   end if

!  #ifdef DEV_HAVE_SCGW_SYM 
   if (Sigp%gwcalctyp>=20 .and. Sigp%symsigma > 0) then 
!    bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw
     do spin=1,Sigp%nsppol
       do ik_ibz=1,Kmesh%nibz
         if (.not.bsym_failed(QP_sym(ik_ibz,spin))) then
           bmin=Sigp%minbnd(ik_ibz,spin); bmax=Sigp%minbnd(ik_ibz,spin)
           do ib2=bmin,bmax
             irr_idx2 = QP_sym(ik_ibz,spin)%b2irrep(ib2)
             do ib1=bmin,bmax
               irr_idx1 = QP_sym(ik_ibz,spin)%b2irrep(ib1)
               if (irr_idx1/=irr_idx2 .and. ALL((/irr_idx1,irr_idx2/)/=0) ) Sr%hhartree(ib1,ib2,ik_ibz,spin) = czero
             end do
           end do
         end if
       end do
     end do
   end if
!  #endif

   ABI_DEALLOCATE(qp_rhog)
   ABI_DEALLOCATE(qp_vhartr)
   ABI_DEALLOCATE(qp_vxc)
   ABI_DEALLOCATE(qp_taug)
   ABI_DEALLOCATE(qp_nhat)
   istat = ABI_ALLOC_STAT
   call destroy_melements(QP_me)
 end if ! gwcalctyp<10
!
!=== Free some memory ===
 if (allocated(hbare))  then
   ABI_DEALLOCATE(hbare)
 end if
 if (allocated(hlda ))  then
   ABI_DEALLOCATE(hlda)
 end if
!
!=== Prepare the storage of QP amplitudes and energies ===
!* Initialize with KS wavefunctions and energies.
 Sr%eigvec_qp   =czero
 Sr%en_qp_diago =zero
 do ib=1,Sigp%nbnds
   Sr%en_qp_diago(ib,:,:)=KS_BSt%eig(ib,:,:)
   Sr%eigvec_qp(ib,ib,:,:)=cone
 end do
!
!=== Store <n,k,s|V_xc[n_val]|n,k,s> and <n,k,s|V_U|n,k,s> ===
!* Note that we store the matrix elements of V_xc in the KS basis set, not in the QP basis set
!* Matrix elements of V_U are zero unless we are using LDA+U as starting point
 do ib=b1gw,b2gw
   Sr%vxcme(ib,:,:)=KS_me%vxcval(ib,ib,:,:)
   if (Dtset%usepawu>0) Sr%vUme (ib,:,:)=KS_me%vu(ib,ib,:,:)
 end do
!
!=== Initial guess for the GW energies ===
!* Save the energies of the previous iteration.
 do spin=1,Sigp%nsppol
   do ik=1,Kmesh%nibz
     do ib=1,Sigp%nbnds
       Sr%e0 (ib,ik,spin)=QP_BSt%eig(ib,ik,spin)
       Sr%egw(ib,ik,spin)=QP_BSt%eig(ib,ik,spin)
     end do
     Sr%e0gap(ik,spin)=zero
     ks_iv=ks_vbik(ik,spin)
     if (Sigp%nbnds>=ks_iv+1) Sr%e0gap(ik,spin)=Sr%e0(ks_iv+1,ik,spin)-Sr%e0(ks_iv,ik,spin)
   end do
 end do
!
!=== If required apply a scissor operator or update the energies ===
!TODO check if other Sr entries have to be updated
!moreover this part should be done only in case of semiconductors
!FIXME To me it makes more sense if we apply the scissor to KS_BS but I have to RECHECK csigme
 if (ABS(Sigp%soenergy)>tol6) then
   write(msg,'(6a,f10.5,2a)')ch10,&
&   ' sigma : performing a first self-consistency',ch10,&
&   '  update of the energies in G by a scissor operator',ch10, &
&   '  applying a scissor operator of ',Sigp%soenergy*Ha_eV,' [eV] ',ch10
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   do spin=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       ks_iv=ks_vbik(ik,spin)
       if (Sigp%nbnds>=ks_iv+1) then
         Sr%egw    (ks_iv+1:Sigp%nbnds,ik,spin) = Sr%egw    (ks_iv+1:Sigp%nbnds,ik,spin)+Sigp%soenergy
         QP_BSt%eig(ks_iv+1:Sigp%nbnds,ik,spin) = QP_BSt%eig(ks_iv+1:Sigp%nbnds,ik,spin)+Sigp%soenergy
       end if
     end do
   end do
!  $call apply_scissor(QP_BSt,Sigp%soenergy)
 else if (.FALSE.) then
   write(msg,'(4a)')ch10,&
&   ' sigma : performing a first self-consistency',ch10,&
&   '  update of the energies in G by a previous GW calculation'
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
!  TODO Recheck this part, is not clear to me!
   ABI_ALLOCATE(igwene,(QP_Bst%mband,QP_Bst%nkpt,QP_Bst%nsppol))
   call rdgw(QP_BSt,'__in.gw__',igwene,extrapolate=.TRUE.)
   ABI_DEALLOCATE(igwene)
   Sr%egw=QP_BSt%eig
!  
!  * Recalculate the new fermi level.
   call update_occ(QP_BSt,Dtset%fixmom,prtvol=0)
 end if
!
!In case of AC refer all the energies wrt to the fermi level
!Take care because results from ppmodel cannot be used for AC
!FIXME check ks_energy or qp_energy (in case of SCGW?)

 if (mod10==SIG_GW_AC) then
!  All these quantities will be passed to csigme
!  if I skipped the self-consistent part then here I have to use fermi
   QP_BSt%eig = QP_BSt%eig -QP_BSt%fermie
   Sr%egw = Sr%egw-QP_BSt%fermie
   Sr%e0  = Sr%e0 -QP_BSt%fermie
   oldefermi=QP_BSt%fermie
!  TODO Recheck fermi
!  Clean EVERYTHING in particulare the treatment of E fermi
   QP_BSt%fermie=zero
 end if
!
!=== Setup frequencies around the KS\QP eigenvalues to compute Sigma derivatives (notice the spin) ===
!TODO it is better using an odd Sr%nomega4sd so that the KS\QP eigenvalue is in the middle
 ioe0j=Sr%nomega4sd/2+1
 do spin=1,Sigp%nsppol
   do ikcalc=1,Sigp%nkptgw
     ib1=Sigp%minbnd(ikcalc,spin)
     ib2=Sigp%maxbnd(ikcalc,spin)
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc))
     do jb=ib1,ib2
       do io=1,Sr%nomega4sd
         Sr%omega4sd(jb,ik_ibz,io,spin)=Sr%egw(jb,ik_ibz,spin)+Sigp%deltae*(io-ioe0j)
       end do
     end do
   end do
 end do

 ABI_TIMER_STOP("hqp_init")
 call timab(408,2,tsec) ! hqp_init

 ABI_TIMER_START("getW")
 call timab(409,1,tsec) ! getW

!=== Get epsilon^{-1} either from the _SCR or the _SUSC file and store it in Er%epsm1 ===
!* If Er%mqmem==0, allocate and read a single q-slice inside csigme.
!TODO Er%nomega should be initialized so that only the frequencies really needed are stored in memory

 if (sigma_needs_w(Sigp)) then
   if (Dtset%gwgamma==0) then
     id_required=4; ikxc=0; approx_type=0; option_test=0; dim_kxcg=0
     ABI_ALLOCATE(kxcg,(nfftf_tot,dim_kxcg))

   else if (Dtset%gwgamma==1.OR.Dtset%gwgamma==2) then ! ALDA TDDFT kernel vertex
!    ABI_CHECK(Dtset%usepaw==0,"GWGamma=1 or 2 + PAW not available")
     ABI_CHECK(Er%ID==0,"Er%ID should be 0")

     if (Dtset%usepaw==1) then ! If we have PAW, we need the full density on the fine grid
       ABI_ALLOCATE(ks_aepaw_rhor,(nfftf,Wfd%nspden))
!      write(std_out,*) 'getpawden: ',Dtset%getpawden,' irdpawden: ',Dtset%irdpawden
       if (Dtset%getpawden==0.and.Dtset%irdpawden==0) then
         MSG_ERROR("Must use get/irdpawden to provide a _PAWDEN file!")
       end if
       call wrtout(std_out,' Checking for existence of: '//TRIM(Dtfil%filpawdensin),"COLL")
       inquire(file=TRIM(Dtfil%filpawdensin),iostat=io)
       if (ios/=0) then
         msg = ' File: '//TRIM(Dtfil%filpawdensin)//' exists but iostat returns nonzero value!'
         MSG_ERROR(msg)
       end if
       call get_rhor(Dtfil%filpawdensin,Dtset%accesswff,Wfd%nspden,nfftf_tot, &
&       ngfftf,Dtset%paral_kgb,MPI_enreg_seq,ks_aepaw_rhor,get_pawden=.TRUE.)
     end if ! Dtset%usepaw==1

     id_required=4; ikxc=7; approx_type=1; dim_kxcg=1
     if (Dtset%gwgamma==1) option_test=1 ! TESTELECTRON, vertex in chi0 *and* sigma
     if (Dtset%gwgamma==2) option_test=0 ! TESTPARTICLE, vertex in chi0 only
     ABI_ALLOCATE(kxcg,(nfftf_tot,dim_kxcg))

     dbg_mode=.FALSE.
     if (Dtset%usepaw==1) then ! Use PAW all-electron density
       call xc_kernel(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,ks_aepaw_rhor,&
&       Er%npwe,dim_kxcg,kxcg,Gsph_c%gvec,xmpi_self,dbg_mode=dbg_mode)
     else ! Norm-conserving
       call xc_kernel(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,ks_rhor,&
&       Er%npwe,dim_kxcg,kxcg,Gsph_c%gvec,xmpi_self,dbg_mode=dbg_mode)
     end if

   else if (Dtset%gwgamma==3.OR.Dtset%gwgamma==4) then ! ADA non-local kernel vertex
!    ABI_CHECK(Dtset%usepaw==0,"GWGamma + PAW not available")
     ABI_CHECK(Er%ID==0,"Er%ID should be 0")
     ABI_CHECK(Sigp%nsppol==1,"ADA vertex for GWGamma not available yet for spin-polarised cases")

     if (Dtset%usepaw==1) then ! If we have PAW, we need the full density on the fine grid
       ABI_ALLOCATE(ks_aepaw_rhor,(nfftf,Sigp%nsppol))
!      write(std_out,*) 'getpawden: ',Dtset%getpawden,' irdpawden: ',Dtset%irdpawden
       if (Dtset%getpawden==0.and.Dtset%irdpawden==0) then
         MSG_ERROR("Must use get/irdpawden to provide a _PAWDEN file!")
       end if
       call wrtout(std_out,' Checking for existence of: '//TRIM(Dtfil%filpawdensin),"COLL")
       inquire(file=TRIM(Dtfil%filpawdensin),iostat=io)
       if (ios/=0) then
         msg = ' File: '//TRIM(Dtfil%filpawdensin)//' exists but iostat returns nonzero value!'
         MSG_ERROR(msg)
       end if
       call get_rhor(Dtfil%filpawdensin,Dtset%accesswff,Wfd%nspden,nfftf_tot, &
&       ngfftf,Dtset%paral_kgb,MPI_enreg_seq,ks_aepaw_rhor,get_pawden=.TRUE.)
     end if ! Dtset%usepaw==1

     id_required=4; ikxc=7; approx_type=2;
     if (Dtset%gwgamma==3) option_test=1 ! TESTELECTRON, vertex in chi0 *and* sigma
     if (Dtset%gwgamma==4) option_test=0 ! TESTPARTICLE, vertex in chi0 only
     ABI_ALLOCATE(fxc_ADA,(Er%npwe,Er%npwe,Er%nqibz))
!    Use userrd to set kappa
     if (Dtset%userrd==zero) Dtset%userrd = 2.1_dp
!    Set correct value of kappa (should be scaled with alpha*r_s where)
!    r_s is Wigner-Seitz radius and alpha=(4/(9*Pi))^(1/3)
     rhoav = (drude_plsmf*drude_plsmf)/four_pi
     r_s = (three/(four_pi*rhoav))**third 
     alpha = (four*ninth*piinv)**third
     Dtset%userrd = Dtset%userrd/(alpha*r_s)

     dbg_mode=.TRUE.
     if (Dtset%usepaw==1) then ! Use PAW all-electron density
       call xc_kernel_ADA(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,&
&       ks_aepaw_rhor,Er%npwe,Er%nqibz,Er%qibz,&
&       fxc_ADA,Gsph_c%gvec,xmpi_self,kappa_init=Dtset%userrd,dbg_mode=dbg_mode)
     else ! Norm conserving
       call xc_kernel_ADA(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,&
&       ks_rhor,Er%npwe,Er%nqibz,Er%qibz,&
&       fxc_ADA,Gsph_c%gvec,xmpi_self,kappa_init=Dtset%userrd,dbg_mode=dbg_mode)
     end if

   end if ! Dtset%gwgamma

!  If reconstruction is wanted, the frequency grid might be different
   if ((Dtset%gw_reconst_scr==1).and.(mod10==SIG_GW_CD.or.mod10==SIG_QPGW_CD) &
&   .and.(Dtset%nfreqre/=Er%nomega_r.or.Dtset%nfreqim/=Er%nomega_i &
&   .or.(ABS(Er%omega(Er%nomega_r)-Dtset%freqremax)>tol6))) then
     MSG_COMMENT("Requested freq. grid seems to be different to the one used in screening calculation!")
     my_plsmf = Dtset%ppmfrq
     if (my_plsmf<tol6) my_plsmf = drude_plsmf
     call recalculate_epsm1_freq_grid(Er,Dtset%nfreqre,Dtset%nfreqim,Dtset%freqremin,&
&     Dtset%freqremax,my_plsmf)
     MSG_COMMENT("Frequency grid recalculated for pole-fit reconstruction.")
   end if

   if (Dtset%gwgamma<3) then

!    #ifdef HAVE_MPI_IO
!    MSG_WARNING("Calling mkdump_Er in MPI_IO mode")
!    call mkdump_Er(Er,Vcp,Er%npwe,Gsph_c%gvec,dim_kxcg,kxcg,id_required,&
!    &     approx_type,ikxc,option_test,Dtfil%fnameabo_scr,IO_MODE_MPI,&
!    &     nfftf_tot,ngfftf,comm)
!    #else
     call mkdump_Er(Er,Vcp,Er%npwe,Gsph_c%gvec,dim_kxcg,kxcg,id_required,&
&     approx_type,ikxc,option_test,Dtfil%fnameabo_scr,Dtset%accesswff,&
&     nfftf_tot,ngfftf,comm,reconstruct_scr=Dtset%gw_reconst_scr)
!    #endif
   else
     call mkdump_Er(Er,Vcp,Er%npwe,Gsph_c%gvec,dim_kxcg,kxcg,id_required,&
&     approx_type,ikxc,option_test,Dtfil%fnameabo_scr,Dtset%accesswff,&
&     nfftf_tot,ngfftf,comm,fxc_ADA)

   end if
   if (allocated(kxcg))  then
     ABI_DEALLOCATE(kxcg)
   end if
   if (allocated(fxc_ADA))  then
     ABI_DEALLOCATE(fxc_ADA)
   end if
   if (allocated(ks_aepaw_rhor))  then
     ABI_DEALLOCATE(ks_aepaw_rhor)
   end if
 end if
!
!================================================
!==== Calculate plasmonpole model parameters ====
!================================================
!TODO Maybe its better if we use mqmem as input variable

 use_aerhor=0
 ABI_ALLOCATE(ks_aepaw_rhor,(nfftf,Wfd%nspden*use_aerhor))

 if (sigma_needs_ppm(Sigp)) then
   my_plsmf=drude_plsmf; if (Dtset%ppmfrq>tol6) my_plsmf=Dtset%ppmfrq
   call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,Sigp%ppmodel,my_plsmf)

!  PPm%force_plsmf= force_ppmfrq  ! this line to change the plasme frequency in HL expression.
   
   if (Wfd%usepaw==1 .and. Ppm%userho==1) then
!    * For PAW and ppmodel 2-3-4 we need the AE rho(G) without compensation charge.
!    * It would be possible to calculate rho(G) using Paw_pwff, though. It should be faster but
!    results will depend on the expression used for the matrix elements. This approach is safer.
     use_aerhor=1
     ABI_DEALLOCATE(ks_aepaw_rhor)
     istat = ABI_ALLOC_STAT
     ABI_ALLOCATE(ks_aepaw_rhor,(nfftf,Wfd%nspden))

!    Check if input density file is available, otherwise compute
     pawden_exists = .FALSE.
     pawden_fname = TRIM(Dtfil%filnam_ds(3))//'_PAWDEN'
     call wrtout(std_out,' Checking for existence of: '//TRIM(pawden_fname),"COLL")
     inquire(file=pawden_fname,iostat=ios,exist=pawden_exists)
     if (ios/=0) then
       msg = ' File: '//TRIM(pawden_fname)//' exists but iostat returns nonzero value!'
       MSG_ERROR(msg)
     end if

     if (pawden_exists) then ! Read density from file
       
       call get_rhor(pawden_fname,Dtset%accesswff,Wfd%nspden,nfftf_tot, &
&       ngfftf,Dtset%paral_kgb,MPI_enreg_seq,ks_aepaw_rhor,get_pawden=.TRUE.)

     else ! Have to calculate PAW AW rhor from scratch
!      call nullify_mpi_enreg(Denfgr_MPI) ! Denfgr is run in parallel inside xmpi_world (see xcomm_init)
!      call initmpi_seq(Denfgr_MPI)       ! It is a dirty trick, I know, but it is the simplest solution
!      as an MPI_type has to be provided.
       ABI_ALLOCATE(qp_rhor_n_one,(pawfgr%nfft,Dtset%nspden))
       ABI_ALLOCATE(qp_rhor_nt_one,(pawfgr%nfft,Dtset%nspden))

!      FIXME
       MSG_WARNING(" denfgr in sigma seems to produce wrong results")
       write(std_out,*)" input tilde ks_rhor integrates: ",SUM(ks_rhor(:,1))*Cryst%ucvol/nfftf

       call denfgr(Cryst%atindx1,Cryst%gmet,comm,Cryst%natom,Cryst%nattyp,ngfftf,ks_nhat,Wfd%nspinor,&
&       Wfd%nsppol,Wfd%nspden,Cryst%ntypat,Pawfgr,Pawrad,KS_Pawrhoij,Pawtab,Dtset%prtvol,Psps,ks_rhor,&
&       ks_aepaw_rhor,qp_rhor_n_one,qp_rhor_nt_one,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

       ABI_DEALLOCATE(qp_rhor_n_one)
       ABI_DEALLOCATE(qp_rhor_nt_one)
!      call destroy_mpi_enreg(Denfgr_MPI)
     end if

     write(msg,'(a,f8.4)')' sigma: PAW AE density used for PPmodel integrates to: ',SUM(ks_aepaw_rhor(:,1))*Cryst%ucvol/nfftf
     call wrtout(std_out,msg,'PERS')

     if (Er%mqmem/=0) then ! Calculate ppmodel parameters for all q.
       call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,nfftf,Gsph_c%gvec,ngfftf,ks_aepaw_rhor(:,1))
     end if

   else ! NC or PAW with PPmodel 1.
     if (Er%mqmem/=0) then ! Calculate ppmodel parameters for all q.
       call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,nfftf,Gsph_c%gvec,ngfftf,ks_rhor(:,1))
     end if
   end if ! PAW or NC PPm and/or needs density
   
!  If nfreqre is set, then output the dielectric function for these frequencies to file
   write(msg,'(a,i0,a,f8.2,a)')&
&   ' Checking if PPm em1 is output, nfreqre is: ',Dtset%nfreqre,' freqremax is: ',Dtset%freqremax*Ha_eV,' eV'
   MSG_WARNING(msg)

   if (Dtset%nfreqre>1) then
     ABI_CHECK(Dtset%freqremax>tol6,'freqremax must be set to maximum frequency')
     ABI_ALLOCATE(em1_ppm,(1,1,Dtset%nfreqre))
     ABI_ALLOCATE(omega,(Dtset%nfreqre))
     do iomega=1,Dtset%nfreqre
       omega(iomega) =  CMPLX((Dtset%freqremax/REAL((Dtset%nfreqre-1)))*(iomega-1),zero)
     end do 

     call getem1_from_PPm(PPm,1,1,Er%Hscr%zcut,Dtset%nfreqre,omega,Vcp,em1_ppm,&
&     only_ig1=1,only_ig2=1)
     
     ppm_unt = get_unit()
     open(unit=ppm_unt,file=Dtfil%fnameabo_em1_lf,status='unknown',form='formatted',iostat=ios)
     write(ppm_unt,'(a)')'#'
     write(ppm_unt,'(a)')'# Macroscopic plasmon-pole Dielectric Function without local fields'
     write(ppm_unt,'(a)')'# Omega [eV]    Re epsilon_M       IM eps_M '
     write(ppm_unt,'(a)')'# ppmodel: ',PPm%model
     do iomega=1,Dtset%nfreqre
       write(ppm_unt,'((3es16.8))')REAL(omega(iomega))*Ha_eV,REAL(em1_ppm(1,1,iomega))
     end do
     close(ppm_unt)

     ABI_DEALLOCATE(em1_ppm)
     ABI_DEALLOCATE(omega)
     msg = ' Wrote epsilon-1 for PPm to file: '//TRIM(Dtfil%fnameabo_em1_lf)
     MSG_WARNING(msg)
   end if ! Check for output of eps^-1 for PPm

 end if ! sigma_needs_ppm

 ABI_TIMER_STOP("getW")
 call timab(409,2,tsec) ! getW
!
 if (wfd_iam_master(Wfd)) then ! Write info on the run on ab_out, then open files to store final results.
   call reportgap(KS_BSt,header='KS Band Gaps',unit=ab_out)
   call write_sigma_results_header(Sigp,Er,Cryst,Kmesh,Qmesh)
   open(unt_gw,file=Dtfil%fnameabo_gw,status='unknown',form='formatted',iostat=ios)
   write(unt_gw,*)Sigp%nkptgw,Sigp%nsppol
!  TODO Change this to a standard abinit output file
   open(687,file='_GWDIAG',status='unknown',form='formatted',iostat=ios)
   write(687,*)Sigp%nkptgw,Sigp%nsppol
   open(unt_sig,file=Dtfil%fnameabo_sig,status='unknown',form='formatted',iostat=ios)
   open(unt_sgr,file=Dtfil%fnameabo_sgr,status='unknown',form='formatted',iostat=ios)
   if (mod10==SIG_GW_AC) then ! Sigma along the imaginary axis.
     open(unt_sgm,file=Dtfil%fnameabo_sgm,status='unknown',form='formatted',iostat=ios)
   end if
 end if
!
!=======================================================================
!==== Calculate self-energy and output the results for each k-point ====
!=======================================================================
!* Here it would be possible to calculate the QP correction for the same k-point using a PPmodel
!in the first iteration just to improve the initial guess and CD or AC in the second step. Useful
!if the KS level is far from the expected QP results. Previously it was possible to do such
!calculation by simply listing the same k-point twice in kptgw. Now this trick is not allowed anymore.
!Everything, indeed, should be done in a clean and transparent way inside csigme.

 call wfd_print(Wfd,mode_paral='PERS')

 sigma_type = sigma_type_from_key(mod10)
 call wrtout(std_out,sigma_type,'COLL')
 
 if (Sigp%gwcalctyp<10) then
   msg = " Perturbative Calculation"
   if (Sigp%gwcalctyp==1) msg = " Newton Raphson method "
 else if (Sigp%gwcalctyp<20) then
   msg = " Self-Consistent on Energies only"
 else
   msg = " Self-Consistent on Energies and Wavefunctions"
 end if
 call wrtout(std_out,msg,'COLL')
!
!=================================================
!==== Calculate the matrix elements of Sigma =====
!=================================================

 nomega_sigc=Sr%nomega_r+Sr%nomega4sd; if (mod10==SIG_GW_AC) nomega_sigc=Sr%nomega_i

 ib1=Sigp%minbdgw ! min and max band indeces for GW corrections.
 ib2=Sigp%maxbdgw
 ABI_ALLOCATE(sigxme,(ib1:ib2,ib1:ib2,Sigp%nkptgw,Sigp%nsppol*Sigp%nsig_ab))
 sigxme=czero
 ABI_ALLOCATE(sigcme,(nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nkptgw,Sigp%nsppol*Sigp%nsig_ab))
 sigcme=czero
!
!==========================================================
!==== Exchange part using the dense gwx_ngfft FFT mesh ====
!==========================================================
!TODO : load distribution is not optimal if band parallelism is used.
!but this problem was also affecting the old implementation.
 call timab(421,1,tsec) ! calc_sigx_me 

 only_one_kpt=(Kmesh%nbz==1)
 call init_gsphere(Gsph_x,only_one_kpt,Cryst,Sigp%npwx,gvec=Gsph_Max%gvec)

 do ikcalc=1,Sigp%nkptgw
   ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW
   ib1=MINVAL(Sigp%minbnd(ikcalc,:)) ! min and max band indeces for GW corrections (for this k-point)
   ib2=MAXVAL(Sigp%maxbnd(ikcalc,:)) 

   sigxme_p => sigxme(ib1:ib2,ib1:ib2,ikcalc,:)

   call calc_sigx_me(ikcalc,ib1,ib2,Cryst,QP_bst,Sigp,Gsph_x,Vcp,Kmesh,Qmesh,Ltg_k(ikcalc),Pawtab,Pawang,Paw_pwff,&
&   Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym(ik_ibz,:),gwx_ngfft,ngfftf,Dtset%prtvol,Dtset%pawcross,sigxme_p)
 end do

 call destroy_gsphere(Gsph_x)

!for the time being, do not remove this barrier!
 call xbarrier_mpi(Wfd%comm)

 call timab(421,2,tsec) ! calc_sigx_me 
!
!==========================================================
!==== Correlation part using the coarse gwc_ngfft mesh ====
!==========================================================

 if (mod10/=SIG_HF) then 
   do ikcalc=1,Sigp%nkptgw
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW
     ib1=MINVAL(Sigp%minbnd(ikcalc,:)) ! min and max band indeces for GW corrections (for this k-point)
     ib2=MAXVAL(Sigp%maxbnd(ikcalc,:)) 
     
     sigcme_p => sigcme(:,ib1:ib2,ib1:ib2,ikcalc,:)
     
     if ( ANY(mod10==(/SIG_SEX, SIG_COHSEX/)) ) then ! Calculate static COHSEX or SEX using the coarse gwc_ngfft mesh.
       call timab(423,1,tsec) ! cohsex_me
       call cohsex_me(ikcalc,nomega_sigc,ib1,ib2,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_c,Vcp,Kmesh,Qmesh,&
&       Ltg_k(ikcalc),Pawtab,Pawang,Paw_pwff,Psps,Wfd,QP_sym(ik_ibz,:),&
&       gwc_ngfft,Dtset%accesswff,Dtset%prtvol,sigcme_p)
       call timab(423,2,tsec) ! cohsex_me

     else
       call timab(424,1,tsec) ! calc_sigc_me
       call calc_sigc_me(ikcalc,nomega_sigc,ib1,ib2,Dtset,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_Max,Gsph_c,Vcp,Kmesh,Qmesh,&
&       Ltg_k(ikcalc),PPm,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym(ik_ibz,:),&
&       gwc_ngfft,ngfftf,nfftf,ks_rhor,use_aerhor,ks_aepaw_rhor,sigcme_p)
       call timab(424,2,tsec) ! calc_sigc_me
     end if

   end do
 end if

 call xbarrier_mpi(Wfd%comm)
!
!=====================================================
!==== Solve Dyson equation storing results in Sr% ====
!=====================================================
!* Use perturbative approach or AC to find QP corrections.
!* If qp-GW, diagonalize also H0+Sigma in the KS basis set to get the
!new QP amplitudes and energies (Sr%eigvec_qp and Sr%en_qp_diago.
!TODO AC with spinor not implemented yet.
!TODO Diagonalization of Sigma+hhartre with AC is wrong.
!
 call timab(425,1,tsec) ! solve_dyson
 do ikcalc=1,Sigp%nkptgw
   ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW
   ib1=MINVAL(Sigp%minbnd(ikcalc,:)) ! min and max band indeces for GW corrections (for this k-point)
   ib2=MAXVAL(Sigp%maxbnd(ikcalc,:)) 
   
   sigcme_p => sigcme(:,ib1:ib2,ib1:ib2,ikcalc,:)
   sigxme_p => sigxme(ib1:ib2,ib1:ib2,ikcalc,:)

   qp_ene => QP_BSt%eig

   call solve_dyson(ikcalc,ib1,ib2,nomega_sigc,Sigp,Kmesh,sigxme_p,sigcme_p,qp_ene,Sr,Dtset%prtvol,Dtfil,Wfd%comm)
!  
!  === Calculate direct gap for each spin and print out final results ===
!  * We use the valence index of the KS system because we still do not know
!  all the QP corrections. Ideally one should use the QP valence index
   do spin=1,Sigp%nsppol
     if ( Sigp%maxbnd(ikcalc,spin) >= ks_vbik(ik_ibz,spin)+1 .and. &
&     Sigp%minbnd(ikcalc,spin) <= ks_vbik(ik_ibz,spin)  ) then
       ks_iv=ks_vbik(ik_ibz,spin)
       Sr%egwgap (ik_ibz,spin)=  Sr%egw(ks_iv+1,ik_ibz,spin) -  Sr%egw(ks_iv,ik_ibz,spin)
       Sr%degwgap(ik_ibz,spin)= Sr%degw(ks_iv+1,ik_ibz,spin) - Sr%degw(ks_iv,ik_ibz,spin)
     else ! The "gap" cannot be computed
       Sr%e0gap  (ik_ibz,spin)=zero
       Sr%egwgap (ik_ibz,spin)=zero
       Sr%degwgap(ik_ibz,spin)=zero
     end if
   end do

   if (wfd_iam_master(Wfd)) call write_sigma_results(ikcalc,ik_ibz,Sigp,Sr,KS_BSt)
 end do !ikcalc

 ABI_DEALLOCATE(sigxme)
 ABI_DEALLOCATE(sigcme)

 call timab(425,2,tsec) ! solve_dyson
 call timab(426,1,tsec) ! finalize
!
!=== Update the energies in QP_BSt ===
!* If QPSCGW, use diagonalized eigenvalues otherwise perturbative results.
 if (Sigp%gwcalctyp>=10) then
   do ib=1,Sigp%nbnds
     QP_BSt%eig(ib,:,:)=Sr%en_qp_diago(ib,:,:)
   end do
 else
   QP_BSt%eig=Sr%egw
 end if

!================================================================================
!==== This part is done only if all k-points in the IBZ have been calculated ====
!================================================================================
 if (Sigp%nkptgw==Kmesh%nibz) then

   call update_occ(QP_BSt,Dtset%fixmom,prtvol=Dtset%prtvol) ! Recalculate new occupations and Fermi level.
   qp_vbik(:,:) = get_valence_idx(QP_BSt)

   write(msg,'(2a,3x,2(es16.6,a))')ch10,' New Fermi energy : ',QP_BSt%fermie,' Ha ,',QP_BSt%fermie*Ha_eV,' eV'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
!  
!  === If all k-points and all occupied bands are calculated, output EXX ===
!  FIXME here be careful about the check on  ks_vbik in case of metals
!  if (my_rank==master.and.Sigp%nkptgw==Kmesh%nibz.and.ALL(Sigp%minbnd(:)==1).and.ALL(Sigp%maxbnd(:)>=MAXVAL(nbv(:)))) then
!  if (ALL(Sigp%minbnd==1).and. ALL(Sigp%maxbnd>=MAXVAL(MAXVAL(ks_vbik(:,:),DIM=1))) ) then
   if (ALL(Sigp%minbnd==1).and. ALL(Sigp%maxbnd>=ks_vbik) ) then  ! FIXME here the two arrays use a different indexing.
     exchange_energy=zero
     do spin=1,Sigp%nsppol
       do ik=1,Kmesh%nibz
         do ib=b1gw,b2gw
           if (Sigp%nsig_ab==1) then
             exchange_energy = exchange_energy + half*QP_BSt%occ(ib,ik,spin)*Kmesh%wt(ik)*Sr%sigxme(ib,ik,spin)
           else
             exchange_energy = exchange_energy + half*QP_BSt%occ(ib,ik,spin)*Kmesh%wt(ik)*SUM(Sr%sigxme(ib,ik,:))
           end if
         end do
       end do
     end do
     write(msg,'(a,2(es16.6,a))')' New Exchange energy : ',exchange_energy,' Ha ,',exchange_energy*Ha_eV,' eV'
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if
!  
!  * Report the QP gaps (Fundamental and Optical)
   call reportgap(QP_BSt,header='QP Band Gaps',unit=ab_out)

!  if ((.FALSE. .and. Dtset%prtwf/=0)) then
!  call wfd_write_wfk(Wfd,Cryst,Psps,Pawtab,QP_Pawrhoij,BSt,Wvl,accesswff,wfk_fname,codvsn,&
!  &  etot,residm,intxc,ixc,pawcpxocc,pawecutdg,stmbias,qptn,so_psp,ngfftdg)
!  end if
 end if ! Sigp%nkptgw==Kmesh%nibz
!
!=== Write SCF data in case of self-consistent calculation ===
!* Save Sr%en_qp_diago, Sr%eigvec_qp and m_lda_to_qp in the _QPS file.
!* Note that in the first iteration qp_rhor contains KS rhor, then the mixed rhor.
 if (Sigp%gwcalctyp>=10) then

   if (wfd_iam_master(Wfd)) then
     call updt_m_lda_to_qp(Sigp,Kmesh,nscf,Sr,m_lda_to_qp) ! Calculate the new m_lda_to_qp
     call wrqps(Dtfil%fnameabo_qps,Sigp,Cryst,Kmesh,Psps,Pawtab,QP_Pawrhoij,&
&     Dtset%nspden,nscf,nfftf,ngfftf,Sr,QP_BSt,m_lda_to_qp,qp_rhor)
   end if
!  
!  === Report the MAX variation for each kptgw and spin ===
   call wrtout(ab_out,ch10//' Convergence of QP corrections ','COLL')
   converged=.TRUE.
   do spin=1,Sigp%nsppol
     write(msg,'(a,i2,a)')' >>>>> For spin ',spin,' <<<<< '
     call wrtout(ab_out,msg,'COLL')
     do ikcalc=1,Sigp%nkptgw
       ib1    = Sigp%minbnd(ikcalc,spin)
       ib2    = Sigp%maxbnd(ikcalc,spin)
       ik_bz  = Sigp%kptgw2bz(ikcalc)
       ik_ibz = Kmesh%tab(ik_bz)
       ii     = imax_loc( ABS(Sr%degw(ib1:ib2,ik_ibz,spin)) )
       max_degw = Sr%degw(ii,ik_ibz,spin)
       write(msg,('(a,i3,a,2f8.3,a,i3)'))'   kptgw no:',ikcalc,'; Maximum DeltaE = (',max_degw*Ha_eV,') for band index:',ii
       call wrtout(ab_out,msg,'COLL')
       converged = (converged .and. ABS(max_degw) < Dtset%gw_toldfeig)
     end do
   end do
 end if
!
!==========================================
!==== Dump GW results to a NETCDF file ====
!==========================================
!WARNING complex real GW corrections is not activated in abi_etsf_init. Here solve problem with v5/t63
 if (wfd_iam_master(Wfd).and.Dtset%accesswff==IO_MODE_ETSF) then
   call abi_etsf_init(Dtset,Dtfil%fnameabo_qp_eig,4,.FALSE.,Hdr_sigma%lmn_size,Psps,Wvl%wfs)
   call etsf_dump_QP(Sr,QP_BSt,KS_BSt,Hdr_sigma,Cryst,Dtfil%fnameabo_qp_eig)
 end if

!----------------------------- END OF THE CALCULATION ------------------------
!
!=====================
!==== Close Files ====
!=====================
 if (wfd_iam_master(Wfd)) then
   close(unt_gw )
   close(687)
   close(unt_sig)
   close(unt_sgr)
   if (mod10==SIG_GW_AC) close(unt_sgm)
 end if

 ABI_DEALLOCATE(ks_vbik)
 ABI_DEALLOCATE(qp_vbik)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(qp_rhor)
 ABI_DEALLOCATE(ks_rhor)
 ABI_DEALLOCATE(ks_rhog)
 ABI_DEALLOCATE(qp_taur)
 ABI_DEALLOCATE(ks_taur)
 ABI_DEALLOCATE(ks_taug)
 ABI_DEALLOCATE(ks_vhartr)
 ABI_DEALLOCATE(ks_vtrial)
 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(ks_vxc)
 ABI_DEALLOCATE(xccc3d)
 ABI_DEALLOCATE(grewtn)
 ABI_DEALLOCATE(xred_dummy)
!if (allocated(qp_irreptab)) deallocate(qp_irreptab)

 if (allocated(kxc))  then
   ABI_DEALLOCATE(kxc)
 end if
 if (allocated(m_lda_to_qp))  then
   ABI_DEALLOCATE(m_lda_to_qp)
 end if
 if (allocated(qp_vtrial  ))  then
   ABI_DEALLOCATE(qp_vtrial)
 end if
!
!===============================================
!==== Free arrays and local data structures ====
!===============================================
 ABI_DEALLOCATE(Pawfgr%fintocoa)
 ABI_DEALLOCATE(Pawfgr%coatofin)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ks_nhat)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ks_nhatgr)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(qp_nhatgr)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(dijexc_core)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ks_aepaw_rhor)
 istat = ABI_ALLOC_STAT

 if (Dtset%usepaw==1) then ! Deallocation for PAW.
   call rhoij_free(KS_Pawrhoij)
   ABI_DEALLOCATE(KS_Pawrhoij)
   call pawfgrtab_free(Pawfgrtab)
   ABI_DEALLOCATE(Pawfgrtab)
   call destroy_paw_ij(KS_paw_ij)
   ABI_DEALLOCATE(KS_paw_ij)
   call destroy_paw_an(KS_paw_an)
   ABI_DEALLOCATE(KS_paw_an)
   call destroy_paw_pwff(Paw_pwff)
   ABI_DEALLOCATE(Paw_pwff)
   if (Sigp%gwcalctyp>=10) then
     call rhoij_free(QP_pawrhoij)
     ABI_DEALLOCATE(QP_pawrhoij)
     call destroy_paw_ij(QP_paw_ij)
     ABI_DEALLOCATE(QP_paw_ij)
     call destroy_paw_an(QP_paw_an)
     ABI_DEALLOCATE(QP_paw_an)
   end if
   if (Dtset%pawcross==1) then
     call destroy_paw_pwaves_lmn(Paw_onsite)
     ABI_DEALLOCATE(Paw_onsite)
     Wfdf%bks_comm = xmpi_comm_null
     call wfd_destroy(Wfdf)
   end if
 end if
!
 call wfd_destroy(Wfd)
 call destroy_little_group(Ltg_k)
 ABI_DEALLOCATE(Ltg_k)
 call destroy_BZ_mesh_type(Kmesh)
 call destroy_BZ_mesh_type(Qmesh)
 call destroy_gsphere(Gsph_Max)
 call destroy_gsphere(Gsph_c)
 call vcoul_free(Vcp)
 call destroy_crystal(Cryst)
 call destroy_Sigma_results(Sr)
 call destroy_Epsilonm1_results(Er)
 call ppm_free(PPm)
 call hdr_clean(Hdr_sigma)
 call hdr_clean(Hdr_kss)
 call bstruct_clean(KS_BSt)
 call bstruct_clean(QP_BSt)
 call destroy_melements(KS_me)
 
 call destroy_bands_symmetries(KS_sym)
 ABI_DEALLOCATE(KS_sym)

 if (Sigp%symsigma==1.and.Sigp%gwcalctyp>=20) then
   call destroy_bands_symmetries(QP_sym)
   ABI_DEALLOCATE(QP_sym)
 end if

 call destroy_sigma_parameters(Sigp)

 call timab(426,2,tsec) ! finalize

 ABI_TIMER_STOP("")
 call timab(401,2,tsec)

 DBG_EXIT('COLL')

end subroutine sigma
!!***
