!{\src2tex{textfont=tt}}
!!****f* ABINIT/screening
!! NAME
!! screening
!!
!! FUNCTION
!! Calculate screening and dielectric functions
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (GMR, VO, LR, RWG, MT, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! codvsn=code version
!! Dtfil<datafiles_type)>=variables related to file names and unit numbers.
!! Pawang<pawang_type)>=paw angular mesh and related data
!! Pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! Pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  Before entering the first time in screening, a significant part of Psps has been initialized:
!!  the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid, ntypat,n1xccc,usepaw,useylm,
!!  and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!  the call to pspini. The next time the code enters screening, Psps might be identical to the
!!  one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!! Output is written on the main output file.
!! The symmetrical inverse dielectric matrix is stored in the _SCR file
!!
!! SIDE EFFECTS
!!  Dtset<type(dataset_type)>=all input variables for this dataset
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      driver,gw_driver
!!
!! CHILDREN
!!      apply_scissor,bstruct_clean,calc_rpa_functional,cchi0,cchi0q0
!!      check_completeness,chi0q0_intraband,chkpawovlp,coeffs_gausslegint
!!      copy_bandstructure,cutoff_density,cutoff_m_elem,destroy_bz_mesh_type
!!      destroy_crystal,destroy_epsilonm1_parameters,destroy_gpairs_type
!!      destroy_gsphere,destroy_little_group,destroy_paw_an,destroy_paw_ij
!!      destroy_paw_pwaves_lmn,destroy_paw_pwff,destroy_spectra,dump_spectra
!!      energies_init,fourdp,free_scrhdr,get_gftt,getph,hdr_clean
!!      init_gpairs_type,init_paw_an,init_paw_ij,init_paw_pwaves_lmn
!!      init_paw_pwff,init_pawfgr,init_scrhdr,initmpi_seq,isfile
!!      make_epsm1_driver,metric,mkem1_q0,mkrdim,nhatgrid,nullify_gpairs_type
!!      nullify_paw_an,nullify_paw_ij,outeps,output_chi0sumrule,pawdenpot
!!      pawdij,pawfgrtab_free,pawfgrtab_init,pawinit,pawmknhat,pawnabla_init
!!      pawprt,pawpuxinit,print_arr,print_ngfft,print_pawtab,print_psps
!!      prtrhomxmn,pspini,rdgw,rdqps,repr_dielconst,rhoij_alloc,rhoij_copy
!!      rhoij_free,rotate_fft_mesh,scr_hdr_io,setsymrhoij,setup_screening
!!      setvtr,symdij,symdij_all,test_charge,timab,update_occ,vcoul_free
!!      wfd_change_ngfft,wfd_copy,wfd_destroy,wfd_init,wfd_mkrho,wfd_print
!!      wfd_read_kss,wfd_read_wfk,wfd_rotate,wfd_test_ortho,write_deltai
!!      write_screening,wrtout,xc_kernel,xc_kernel_ada,xmpi_split_work2_i4b
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine screening(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_timer
 use m_errors
#if defined HAVE_DFT_LIBXC
 use libxc_functionals
#endif

 use m_energies,      only : energies_type, energies_init
 use m_numeric_tools, only : print_arr, iseven
 use m_geometry,      only : normv, vdotw
 use m_gwdefs,        only : GW_TOLQ0, GW_TOLQ, destroy_epsilonm1_parameters, epsilonm1_parameters, GW_Q0_DEFAULT,&
&                            gw_uses_wfk_file
 use m_header,        only : hdr_clean
 use m_crystal,       only : destroy_crystal, crystal_structure
 use m_ebands,        only : update_occ, copy_bandstructure, get_valence_idx, get_occupied, apply_scissor, &
&                            bstruct_clean, bst_ismetal
 use m_bz_mesh,       only : bz_mesh_type, destroy_bz_mesh_type, little_group, destroy_little_group
 use m_gsphere,       only : destroy_gsphere, nullify_Gpairs_type, destroy_Gpairs_type, init_Gpairs_type,&
&                            gvectors_type, gpairs_type
 use m_vcoul,         only : vcoul_t, vcoul_free, cutoff_density
 use m_qparticles,    only : rdqps, rdgw, show_QP
 use m_screening,     only : make_epsm1_driver, outeps, mkem1_q0
 use m_io_screening,  only : init_ScrHdr, scr_hdr_io, write_screening, free_scrhdr, scrhdr_type
 use m_spectra,       only : spectra_type, dump_spectra, repr_dielconst, destroy_spectra, W_EM_LF, W_EM_NLF, W_EELF
 use m_fft_mesh,      only : rotate_FFT_mesh, cigfft, get_gftt, print_ngfft
 use m_io_kss,        only : wfd_read_kss
 use m_wfs,           only : wfd_init, wfd_destroy,  wfd_nullify, wfd_print, wfs_descriptor, wfd_rotate, wfd_test_ortho,&
&                            wfd_read_wfk, wfd_test_ortho, wfd_copy, wfd_change_ngfft
 use m_paw_toolbox,   only : nullify_paw_ij, init_paw_ij, destroy_paw_ij, init_pawfgr, pawfgrtab_free, pawfgrtab_init,&
&                            nullify_paw_an, init_paw_an, destroy_paw_an, print_pawtab, paw_pwaves_lmn_t,init_paw_pwaves_lmn,&
&                            nullify_paw_pwaves_lmn,destroy_paw_pwaves_lmn
 use m_paw_pwij,      only : paw_pwff_type, init_paw_pwff, destroy_paw_pwff
 use m_io_tools,      only : get_unit, flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screening'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_69_wfdesc
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=6),intent(in) :: codvsn
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pawang_type),intent(inout) :: Pawang
 type(Pseudopotential_type),intent(inout) :: Psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 type(Pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Dtset%usepaw)
 type(Pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Dtset%usepaw)

!Local variables ------------------------------
 character(len=4) :: ctype='RPA '
!scalars
 integer,parameter :: level=30,tim_fourdp=4,NOMEGA_PRINTED=15,localrdwf0=0
 integer,save :: nsym_old=-1
 integer :: spin,my_nspins,isp,ik_ibz !,iatom !
 integer :: choice,cplex,cplex_dij,dim_kxcg,dim_wing
 integer :: fileID,fform_chi0,fform_em1,iat,band,ider,idir,ierr
 integer :: ifft,ii,ikbz,ikxc,initialized,iomega,ios,ipert
 integer :: iqibz,iqcalc,is_qeq0,istat,isym,itypat,izero,ifirst,ilast
 integer :: label,master,mgfftf,mgfftgw
 integer :: moved_atm_inside,moved_rhor,my_maxb,my_minb !,mkmem_
 integer :: nbcw,nbsc,nbvw,nkxc,nkxc1,n3xccc,optene,istep
 integer :: nfftf,nfftf_tot,nfftgw,nfftgw_tot,nhatgrdim,nprocs,nspden_rhoij
 integer :: nscf,nzlmopt,mband
 integer :: optcut,optgr0,optgr1,optgr2,option,approx_type,option_test,optgrad
 integer :: optrad,optrhoij,psp_gencond,my_rank,rdwr
 integer :: rhoxsp_method,comm,test_type,tordering,unt_em1,unt_susc,unt_delI,usexcnhat
 real(dp) :: compch_fft,compch_sph,diecut_eff_dum,domegareal,e0
 real(dp) :: ecore,ecut_eff,ecutdg_eff,nelect
 real(dp) :: gsqcutc_eff,gsqcutf_eff,omegaplasma,ucvol,vxcavg
 real(dp) :: gw_gsq,length,r_s,alpha,rhoav,opt_ecut,factor
 logical :: found,iscompatibleFFT,use_tr,is_first_qcalc,add_chi0_intraband !,ltest
 logical :: update_energies=.FALSE.
 character(len=10) :: string
 character(len=500) :: msg
 character(len=80) :: bar
 type(Bandstructure_type) :: KS_BSt,QP_BSt
 type(BZ_mesh_type) :: Kmesh,Qmesh
 type(vcoul_t) :: Vcp
 type(Crystal_structure) :: Cryst
 type(Epsilonm1_parameters) :: Ep
 type(Energies_type) :: KS_energies
 type(Gpairs_type) :: Gpairs_q
 type(Gvectors_type) :: Gsph_epsG0,Gsph_wfn
 type(Hdr_type) :: Hdr_kss,Hdr_local
 type(MPI_type) :: MPI_enreg_seq
 type(Pawfgr_type) :: Pawfgr
 type(ScrHdr_type) :: Hem1,Hchi0
 type(wfs_descriptor) :: Wfd,Wfdf
 type(spectra_type) :: Spectra
 type(wvl_internal_type) :: wvl
!arrays
 integer,save :: paw_gencond(6)=(/-1,-1,-1,-1,-1,-1/)
 integer :: ibocc(Dtset%nsppol),ngfft_gw(18),ngfftc(18),ngfftf(18),my_spins(Dtset%nsppol)
 integer,allocatable :: irottb(:,:)
 integer,allocatable :: istart(:),istop(:),ktabr(:,:),ktabrf(:,:),l_size_atm(:)
 integer,allocatable :: ks_vbik(:,:),ks_occ_idx(:,:),qp_vbik(:,:),nlmn(:),nband(:,:)
 integer,allocatable :: nq_spl(:)
 integer,allocatable :: gw_gfft(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),k0(3),qtmp(3),rmet(3,3),rprimd(3,3),tsec(2),strsxc(6)
 real(dp),allocatable :: igwene(:,:,:)
 real(dp),allocatable :: chi0_sumrule(:)
 real(dp),allocatable :: ec_rpa(:)
 real(dp),allocatable :: nhat(:,:),nhatgr(:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph1df(:,:)
 real(dp),allocatable :: rhog(:,:),rhor(:,:),rhor_p(:,:),rhor_kernel(:,:)
 real(dp),allocatable :: taug(:,:),taur(:,:)
 real(dp),allocatable :: z(:),zw(:)
 real(dp),allocatable :: grewtn(:,:),xred_dummy(:,:),kxc(:,:),qmax(:)
 real(dp),allocatable :: ks_vhartr(:),vpsp(:),ks_vtrial(:,:),ks_vxc(:,:),xccc3d(:)
 complex(dpc) :: wng(3),em1_00
 complex(gwpc),allocatable :: arr_99(:,:)
 complex(gwpc),allocatable :: kxcg(:,:),fxc_ADA(:,:,:)
 complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:)
 complex(dpc),allocatable :: chi0_lwing(:,:,:),chi0_uwing(:,:,:),chi0_head(:,:,:),wtest(:),eps_head(:,:,:)
 complex(dpc),allocatable :: chi0intra_lwing(:,:,:),chi0intra_uwing(:,:,:),chi0intra_head(:,:,:)
 complex(gwpc),allocatable,target :: chi0(:,:,:),chi0intra(:,:,:),deltaI(:,:,:)
 complex(gwpc),pointer :: epsm1(:,:,:),vc_sqrt(:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 character(len=80) :: title(2)
 character(len=fnlen) :: fname,gw_fname,wfk_fname
 type(Little_group),pointer :: Ltg_q(:)
 type(Paw_an_type),allocatable :: Paw_an(:)
 type(Paw_ij_type),allocatable :: Paw_ij(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:),prev_Pawrhoij(:)
 type(Paw_pwff_type),allocatable :: Paw_pwff(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

#undef HAVE_GW_CUTOFF
#if defined HAVE_GW_CUTOFF
 ! * Variables added for cutoffed matrix elements
 integer :: direction
 real(dp) :: width,z0
#endif

!************************************************************************

 DBG_ENTER("COLL")

 ABI_TIMER_START("")
 call timab(301,1,tsec) ! overall time

 ABI_TIMER_START("init")
 call timab(302,1,tsec) ! screening(init

 write(msg,'(6a)')&
& ' SCREENING: Calculation of the susceptibility and dielectric matrices ',ch10,ch10,&
& ' Based on a program developped by R.W. Godby, V. Olevano, G. Onida, and L. Reining.',ch10,&
& ' Incorporated in ABINIT by V. Olevano, G.-M. Rignanese, and M. Torrent.'
 call wrtout(ab_out, msg,'COLL')
 call wrtout(std_out,msg,'COLL')

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
 call wrtout(ab_out, msg,'COLL')
!
!=== Initialize MPI variables, and parallelization level ===
!* gwpara: 0--> sequential run, 1--> parallelism over k-points, 2--> parallelism over bands.
!* If gwpara==2, each node has both fully and partially occupied states while conduction bands are divided

 comm = xmpi_world
 my_rank = xcomm_rank(comm)
 nprocs  = xcomm_size(comm)
 master = 0

!* Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

!=== Some variables need to be initialized/nullify at start ===
 call energies_init(KS_energies)
 usexcnhat=0

!Nullify the pointers in the data types.
 call nullify_Gpairs_type(Gpairs_q)

 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)

!
!=== Define FFT grid(s) sizes ===
!* Be careful! This mesh is only used for densities and potentials. It is NOT the (usually coarser)
!GW FFT mesh employed for the oscillator matrix elements that is defined in setmesh.F90.
!See also NOTES in the comments at the beginning of this file.
!* NOTE: The mesh is defined in invars2m using ecutwfn, in GW Dtset%ecut is forced to be equal to Dtset%ecutwfn.

 k0(:)=zero
 call init_pawfgr(Dtset,Pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
& gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=gmet,k0=k0)

 call print_ngfft(ngfftf,'Dense FFT mesh used for densities and potentials')
 nfftf_tot=PRODUCT(ngfftf(1:3))

!=============================================
!==== Open and read pseudopotential files ====
!=============================================
 call pspini(Dtset,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,level,Pawrad,Pawtab,Psps,rprimd)
 if (psp_gencond==1) call print_psps(Psps,std_out,0,'COLL')

!=== Initialize dimensions and basic objects ===
 call setup_screening(codvsn,acell,rprim,ngfftf,Dtfil%fnameabi_kss,Dtset,Psps,Pawtab,&
& ngfft_gw,Hdr_kss,Hdr_local,Cryst,Kmesh,Qmesh,KS_BSt,Ltg_q,Gsph_epsG0,Gsph_wfn,Vcp,Ep,comm)

 ABI_TIMER_STOP("init")
 call timab(302,2,tsec) ! screening(init)

 call print_ngfft(ngfft_gw,'FFT mesh used for oscillator strengths')

 nfftgw_tot=PRODUCT(ngfft_gw(1:3))
 mgfftgw   =MAXVAL (ngfft_gw(1:3))
 nfftgw    =nfftgw_tot ! no FFT //

!TRYING TO RECREATE AN "ABINIT ENVIRONMENT"
 KS_energies%e_corepsp=ecore/Cryst%ucvol
!
!==========================
!=== PAW initialization ===
!==========================
 if (Dtset%usepaw==1) then
   ABI_TIMER_START("paw_init")
   call timab(315,1,tsec) ! screening(pawin

   call chkpawovlp(Cryst%natom,Cryst%ntypat,Dtset%pawovlp,Pawtab,Cryst%rmet,Cryst%typat,Cryst%xred)

   ABI_ALLOCATE(Pawrhoij,(Cryst%natom))
   ABI_ALLOCATE(nlmn,(Cryst%ntypat))
   do itypat=1,Cryst%ntypat
     nlmn(itypat)=Pawtab(itypat)%lmn_size
   end do
   nspden_rhoij=Dtset%nspden; if (Dtset%pawspnorb>0.and.Dtset%nspinor==2) nspden_rhoij=4
   call rhoij_alloc(Dtset%pawcpxocc,nlmn,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,&
&   Pawrhoij,Cryst%typat,MPI_enreg=MPI_enreg_seq)
   ABI_DEALLOCATE(nlmn)
!  
!  ==== Initialize values for several basic arrays stored in Pawinit ====
!  TODO Check pawxcdev>2 since gaunt coefficients are allocated with different sizes.
   if (psp_gencond==1.or.&
&   paw_gencond(1)/=Dtset%pawlcutd .or.paw_gencond(2)/=Dtset%pawlmix  .or.&
&   paw_gencond(3)/=Dtset%pawnphi  .or.paw_gencond(4)/=Dtset%pawntheta.or.&
&   paw_gencond(5)/=Dtset%pawspnorb.or.paw_gencond(6)/=Dtset%pawxcdev) then

     diecut_eff_dum=ABS(Dtset%diecut)*Dtset%dilatmx**2

     call pawinit(diecut_eff_dum,Psps%indlmn,Dtset%pawlcutd,Dtset%pawlmix,Psps%lmnmax,Psps%mpsang,&
&     Dtset%pawnphi,Cryst%nsym,Dtset%pawntheta,Cryst%ntypat,Pawang,Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev)

     paw_gencond(1)=Dtset%pawlcutd;  paw_gencond(2)=Dtset%pawlmix
     paw_gencond(3)=Dtset%pawnphi;   paw_gencond(4)=Dtset%pawntheta
     paw_gencond(5)=Dtset%pawspnorb; paw_gencond(6)=Dtset%pawxcdev
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
!  
!  * Evaluate <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for the long wavelength limit.
!  TODO solve problem with memory leak and clean this part as well as the associated flag
   call pawnabla_init(Psps%mpsang,Psps%lmnmax,Cryst%ntypat,Psps%indlmn,Pawrad,Pawtab)

!  FIXME see above. Note that here nsym_old comes from KSS.
!  if (psp_gencond==1.or.nsym_old/=Cryst%nsym) then
   call setsymrhoij(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,rprimd,Cryst%symrec,Pawang%zarot)
   nsym_old=Cryst%nsym
!  end if

!  * Initialize and compute data for LDA+U.
!  paw_dmft%use_dmft=dtset%usedmft
   if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
     call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,&
&     Psps%indlmn,Psps%lmnmax,Cryst%ntypat,Pawang,Dtset%pawprtvol,Pawrad,Pawtab,Dtset%upawu,&
&     Dtset%usedmft,Dtset%useexexch,Dtset%usepawu)
   end if

   ABI_CHECK(Dtset%useexexch==0,"LEXX not yet implemented in GW")
   ABI_CHECK(Dtset%usedmft==0,"DMFT + GW not available")

   if (my_rank==master) call print_pawtab(Pawtab)

!  === Eventually open temporary storage file ===
!  FIXME also mkmem_ not yet defined
!  if (mkmem_==0) then
!  open(Dtfil%unpaw,file=dtfil%fnametmp_paw,form='unformatted',status='unknown')
!  rewind(unit=Dtfil%unpaw)
!  end if

!  === Get Pawrhoij from the header of the KSS file ===
   call rhoij_copy(Hdr_kss%Pawrhoij,Pawrhoij,MPI_enreg=MPI_enreg_seq)

!  === Re-symmetrize symrhoij ===
!  this call leads to a SIGFAULT, likely some pointer is not initialized correctly
   choice=1; optrhoij=1; ipert=0; idir=0
!  call symrhoij(choice,Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,Psps%lmnmax,Cryst%natom,Cryst%natom,Cryst%nsym,&
!  &  Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,Pawrhoij,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
!  
!  === Evaluate form factors for the radial part of phi.phj-tphi.tphj ===
   rhoxsp_method=2
   if (Dtset%userie==1) rhoxsp_method=1 ! Arnaud-Alouani
   if (Dtset%userie==2) rhoxsp_method=2 ! Shiskin-Kresse

   ABI_ALLOCATE(gw_gfft,(3,nfftgw_tot))
   call get_gftt(ngfft_gw,(/zero,zero,zero/),gmet,gw_gsq,gw_gfft)
   ABI_DEALLOCATE(gw_gfft)

!  Set up q grids, make qmax 20% larger than largest expected:
   ABI_ALLOCATE(nq_spl,(Psps%ntypat))
   ABI_ALLOCATE(qmax,(Psps%ntypat))
   nq_spl = Psps%mqgrid_ff
   qmax = SQRT(gw_gsq)*1.2d0  !qmax = Psps%qgrid_ff(Psps%mqgrid_ff)
   ABI_ALLOCATE(Paw_pwff,(Psps%ntypat))

   call init_paw_pwff(Paw_pwff,rhoxsp_method,nq_spl,qmax,gmet,Pawrad,Pawtab,Psps)
   ABI_DEALLOCATE(nq_spl)
   ABI_DEALLOCATE(qmax)
!  
!  === Variables/arrays related to the fine FFT grid ===
   ABI_ALLOCATE(nhat,(nfftf,Dtset%nspden))
   nhat=zero; cplex=1
   ABI_ALLOCATE(Pawfgrtab,(Cryst%natom))
   ABI_ALLOCATE(l_size_atm,(Cryst%natom))
   do iat=1,Cryst%natom
     l_size_atm(iat)=Pawtab(Cryst%typat(iat))%l_size
   end do
   call pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,Dtset%nspden)
   ABI_DEALLOCATE(l_size_atm)
   compch_fft=greatest_real
   usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
!  * 0 --> Vloc in atomic data is Vbare    (Blochl s formulation)
!  * 1 --> Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   write(msg,'(a,i3)')' screening : using usexcnhat = ',usexcnhat
   call wrtout(std_out,msg,'COLL')
!  
!  Identify parts of the rectangular grid where the density has to be calculated.
   optcut=0;optgr0=Dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-Dtset%pawstgylm
   if (Dtset%pawcross==1) optrad=1
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(Cryst%atindx1,gmet,MPI_enreg_seq,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%ucvol,Cryst%xred)

   if (Dtset%pawcross==1) then
     ABI_ALLOCATE(Paw_onsite,(Cryst%natom))
     optgrad=1
     call init_paw_pwaves_lmn(Paw_onsite,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%rprimd,Cryst%xcart,&
&     Psps,Pawtab,Pawrad,Pawfgrtab,optgrad)
!    call nullify_paw_pwaves_lmn(Paw_onsite)
   end if
   ABI_TIMER_STOP("paw_init")
   call timab(315,2,tsec) ! screening(pawin
 end if ! End of PAW initialization.

 call timab(316,1,tsec) ! screening(wfs
 ABI_TIMER_START("wfd_init")

!Allocate these arrays anyway, since they are passed to subroutines.
 if (.not.allocated(nhat))  then
   ABI_ALLOCATE(nhat,(nfftf,0))
 end if
!
!=====================================================
!=== Prepare the distribution of the wavefunctions ===
!=====================================================
!* valence and partially occupied are replicate on each node  while conduction bands are MPI distributed.
!This method is mandatory if gwpara==2 and/or we are using awtr==1 or the spectral method.
!* If awtr==1, we evaluate chi0 taking advantage of time-reversal (speed-up~2)
!* Useful indeces:
!nbvw = Max. number of fully/partially occupied states over spin
!nbcw = Max. number of unoccupied states considering the spin
!TODO: Here for semiconducting systems we have to be sure that each processor has all the
!states considered in the SCGW, moreover nbsc<nbvw
!TODO in case of SCGW vale and conduction has to be recalculated to avoid errors
!if a metal becomes semiconductor or viceversa.
!TODO: Ideally nbvw should include only the states v such that the transition
!c-->v is taken into account in cchi0 (see GW_TOLDOCC). In the present implementation
!Wfd_val contains all the states whose occupation is less than tol8. Due the the long
!tail of the smearing function it happens that a large number of states are allocated
!on each node. This facilitate the calculation of the density but it is a waste of memory
!and CPU time in cchi0.
!Solution: change wfd_mkrho to perform the calculation in parallel also if gwpara==2
!Wfd_val should contain only the occupied states contributing to chi0

 ABI_ALLOCATE(ks_occ_idx,(KS_BSt%nkpt,KS_BSt%nsppol))
 ABI_ALLOCATE(ks_vbik   ,(KS_BSt%nkpt,KS_BSt%nsppol))
 ABI_ALLOCATE(qp_vbik   ,(KS_BSt%nkpt,KS_BSt%nsppol))

 call update_occ(KS_BSt,Dtset%fixmom,prtvol=0)
 ks_occ_idx = get_occupied(KS_BSt,tol8) ! tol8 to be consistent when the density
 ks_vbik    = get_valence_idx(KS_BSt)   ! is calculated locally using Wfd_val without MPI.

 ibocc(:)=MAXVAL(ks_occ_idx(:,:),DIM=1) ! Max occupied band index for each spin.
 ABI_DEALLOCATE(ks_occ_idx)

 use_tr=.FALSE.; nbvw=0
 if (Dtset%gwpara==2.or.Ep%awtr==1.or.Dtset%spmeth>0) then
!  if (.TRUE.) then
   use_tr=.TRUE.
   nbvw=MAXVAL(ibocc)
   nbcw=Ep%nbnds-nbvw
   write(msg,'(4a,i5,2a,i5,2a,i5)')ch10,&
&   ' screening : taking advantage of time-reversal symmetry ',ch10,&
&   ' Maximum band index for partially occupied states nbvw = ',nbvw,ch10,&
&   ' Remaining bands to be divided among processors   nbcw = ',nbcw,ch10,&
&   ' Number of bands treated by each node ~',nbcw/nprocs
   if (nprocs>1.or.Ep%awtr==1.or.Dtset%spmeth>0) call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'PERS')
   if (Cryst%timrev/=2) then
     write(msg,'(a)')'Time-reversal cannot be used since Cryst%timrev/=2'
     MSG_ERROR(msg)
   end if
 end if

 mband=Ep%nbnds
 ABI_ALLOCATE(nband,(Kmesh%nibz,Dtset%nsppol))
 nband=mband
 ABI_ALLOCATE(bks_mask,(mband,Kmesh%nibz,Dtset%nsppol))
 bks_mask=.FALSE.
 ABI_ALLOCATE(keep_ur,(mband,Kmesh%nibz,Dtset%nsppol))
 keep_ur=.FALSE.

 my_nspins=Dtset%nsppol; my_spins=(/(isp,isp=1,Dtset%nsppol)/)
 my_minb=1; my_maxb=Ep%nbnds

 select case (Dtset%gwpara) 
   case (1)
     call wrtout(std_out,'Parallel over transitions. NO memory reduction except for the spin.','COLL')

     if (Dtset%nsppol==2.and.iseven(nprocs)) then ! Distribute spins
       my_nspins=1
       my_spins(1)=1; if (my_rank+1>nprocs/2) my_spins(1)=2
     end if
     
     do isp=1,my_nspins
       spin = my_spins(isp)
       bks_mask(1:Ep%nbnds,:,spin) = .TRUE. 
       if (MODULO(Dtset%gwmem,10)==1) keep_ur(my_minb:my_maxb,:,spin)=.TRUE.
     end do

   case (2)
     call wrtout(std_out,'Distributing bands across the nodes','COLL')
     ABI_ALLOCATE(istart,(nprocs))
     ABI_ALLOCATE(istop,(nprocs))

!    
!    it is now meaningless to distinguish gwcomp=0 or 1
!    since the workload is well balanced later on
!    if (Ep%gwcomp==0) then ! No completeness trick, each proc has fully and partially occupied states.
     call xmpi_split_work2_i4b(nbcw,nprocs,istart,istop,msg,ierr)
     if (ierr/=0) then 
       MSG_WARNING(msg)
     end if

     my_minb=nbvw+istart(my_rank+1)
     my_maxb=nbvw+istop (my_rank+1)

!    else ! Use completeness trick, divide entire set of states.
!    call xmpi_split_work2_i4b(Ep%nbnds,nprocs,istart,istop,msg,ierr)
!    if (ierr/=0) then 
!    MSG_WARNING(msg)
!    end if
!    
!    my_minb=istart(my_rank+1)
!    my_maxb=istop (my_rank+1)
!    end if
     
     if (my_maxb-my_minb+1<=0) then
       write(msg,'(3a,2(i4,a),a)')&
&       ' One or more processors has zero number of bands ',ch10,&
&       ' my_minb= ',my_minb,' my_maxb= ',my_maxb,ch10,&
&       ' This is a waste, decrease the number of processors. '
!      MSG_PERS_WARNING(msg)
       MSG_PERS_ERROR(msg)
     end if
     ABI_DEALLOCATE(istart)
     ABI_DEALLOCATE(istop)

     bks_mask(my_minb:my_maxb,:,:)=.TRUE.
     if (MODULO(Dtset%gwmem,10)==1) keep_ur(my_minb:my_maxb,:,:)=.TRUE.

!    Memory distribution over spins.
#if 0
     bks_mask=.FALSE.; my_spin=0 
     if (Dtset%nsppol==2.and.iseven(nprocs)) then
       my_spin=1; if (1+my_rank>nprocs/2) my_spin=2
       do band=1,mband 
         if (MODULO(band,nprocs/2)==my_rank-(my_spin-1)*(nprocs/2)) bks_mask(band,:,my_spin)=.TRUE.
       end do
     else 
       do band=1,mband 
         if (MODULO(band,nprocs)==my_rank) bks_mask(band,:,:)=.TRUE.
       end do
     end if
#endif

!    This is to have the occupied states on each node.
     do isp=1,my_nspins
       spin = my_spins(isp)
       bks_mask(1:nbvw,:,spin) = .TRUE.
       if (MODULO(Dtset%gwmem,10)==1) keep_ur(1:nbvw,:,spin)=.TRUE.
     end do

     case default 
     MSG_ERROR("Wrong value for gwpara")
 end select

!=== Initialize the Wf_info object ===
!* Allocate %ug and if required, also %ur.

 opt_ecut=zero
 if (gw_uses_wfk_file) opt_ecut=Dtset%ecutwfn
!if (gw_uses_wfk_file) opt_ecut=Hdr_kss%ecut_eff*1.2

 call wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,Ep%npwwfn,mband,nband,Ep%nkibz,Dtset%nsppol,bks_mask,&
& Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_kss%istwfk,Kmesh%ibz,ngfft_gw,&
& Gsph_wfn%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)

 ABI_DEALLOCATE(bks_mask)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(keep_ur)

 call wfd_print(Wfd,mode_paral='PERS')
!FIXME: Rewrite the treatment of use_tr branches in cchi0 ...
!Use a different nbvw for each spin.
!Now use_tr means that one can use time-reversal symmetry.

!==================================================
!==== Read KS band structure from the KSS file ====
!==================================================
!
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
   call wfd_read_kss(Wfd,Dtfil%fnameabi_kss,Ep%nbnds,Dtset%accesswff,nelect)
 end if

 if (Dtset%pawcross==1) then
   call wfd_copy(Wfd,Wfdf)
   call wfd_change_ngfft(Wfdf,Cryst,Psps,ngfftf)
 end if

 call wfd_test_ortho(Wfd,Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

 call timab(316,2,tsec) ! screening(wfs  
 ABI_TIMER_STOP("wfd_init")

 call timab(319,1,tsec) ! screening(1)

 if (Cryst%nsym/=Dtset%nsym .and. Dtset%usepaw==1) then 
   MSG_ERROR('Cryst%nsym/=Dtset%nsym, check pawinit and symrhoij')
 end if
!
!=== Get the FFT index of $ (R^{-1}(r-\tau)) $ ===
!* S= $\transpose R^{-1}$ and k_BZ = S k_IBZ
!* irottb is the FFT index of $ R^{-1} (r-\tau) $ used to symmetrize u_Sk.
 ABI_ALLOCATE(irottb,(nfftgw,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,ngfft_gw,irottb,iscompatibleFFT)

 ABI_ALLOCATE(ktabr,(nfftgw,Kmesh%nbz))
 do ikbz=1,Kmesh%nbz
   isym=Kmesh%tabo(ikbz)
   do ifft=1,nfftgw
     ktabr(ifft,ikbz)=irottb(ifft,isym)
   end do
 end do
 ABI_DEALLOCATE(irottb)

 if (Dtset%usepaw==1 .and. Dtset%pawcross==1) then
   ABI_ALLOCATE(irottb,(nfftf,Cryst%nsym))
   call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,ngfftf,irottb,iscompatibleFFT)

   ABI_ALLOCATE(ktabrf,(nfftf,Kmesh%nbz))
   do ikbz=1,Kmesh%nbz
     isym=Kmesh%tabo(ikbz)
     do ifft=1,nfftf
       ktabrf(ifft,ikbz)=irottb(ifft,isym)
     end do
   end do
   ABI_DEALLOCATE(irottb)
 end if

!=== Compute structure factor phases and large sphere cut-off ===
!WARNING cannot use Dtset%mgfft, this has to be checked better
!mgfft=MAXVAL(ngfftc(:))
!allocate(ph1d(2,3*(2*mgfft+1)*Cryst%natom),ph1df(2,3*(2*mgfftf+1)*Cryst%natom))
 write(std_out,*)' CHECK ',Dtset%mgfftdg,mgfftf
 if (Dtset%mgfftdg/=mgfftf) then
   write(std_out,*)"WARNING Dtset%mgfftf /= mgfftf"
!  write(std_out,*)'HACKING Dtset%mgfftf'
!  Dtset%mgfftdg=mgfftf
 end if
 ABI_ALLOCATE(ph1d,(2,3*(2*Dtset%mgfft+1)*Cryst%natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*Cryst%natom))
 call getph(Cryst%atindx,Cryst%natom,ngfftc(1),ngfftc(2),ngfftc(3),ph1d,Cryst%xred)

 if (Psps%usepaw==1.and.Pawfgr%usefinegrid==1) then
   call getph(Cryst%atindx,Cryst%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,Cryst%xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

!=== Initialize QP_BSt using KS bands ===
!* In case of SCGW, update QP_BSt using the QPS file.
 call copy_bandstructure(KS_BSt,QP_BSt)

 call timab(319,2,tsec) ! screening(1)
!
!============================
!==== Self-consistent GW ====
!============================

 if (Ep%gwcalctyp>=10) then
   ABI_TIMER_START("KS-->QP")
   call timab(304,1,tsec) ! KS => QP; [wfrg]
!  
!  Initialize with KS eigenvalues and eigenfunctions.
   ABI_ALLOCATE(m_lda_to_qp,(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol))
   m_lda_to_qp = czero
   do spin=1,Wfd%nsppol
     do ik_ibz=1,Wfd%nkibz
       do band=1,Wfd%nband(ik_ibz,spin)
         m_lda_to_qp(band,band,:,:) = cone
       end do
     end do
   end do
!  
!  Read Unitary transformation and QP energies.
!  TODO switch on the renormalization of n in screening, QPS should report bdgw
   ABI_ALLOCATE(rhor_p,(nfftf,Dtset%nspden))
   ABI_ALLOCATE(prev_Pawrhoij,(Cryst%natom*Psps%usepaw))

   call rdqps(QP_BSt,Dtfil%fnameabi_qps,Dtset%usepaw,Dtset%nspden,1,nscf,&
&   nfftf,ngfftf,Cryst%ucvol,Dtset%paral_kgb,Cryst,Pawtab,MPI_enreg_seq,nbsc,m_lda_to_qp,rhor_p,prev_Pawrhoij)

   ABI_DEALLOCATE(rhor_p)
   ABI_DEALLOCATE(prev_Pawrhoij)

!  FIXME this is to preserve the old implementation for the head and the wings in ccchi0q0
!  But has to be rationalized
   KS_BSt%eig=QP_BSt%eig
!  
!  Calculate new occ. factors and fermi level.
   call update_occ(QP_BSt,Dtset%fixmom)
   qp_vbik(:,:) = get_valence_idx(QP_BSt)
!  
!  === Update only the wfg treated with GW ===
!  For PAW update and re-symmetrize cprj in the full BZ, TODO add rotation in spinor space
   if (nscf/=0) then 
     call wfd_rotate(Wfd,Cryst,m_lda_to_qp)
   end if

   ABI_DEALLOCATE(m_lda_to_qp)
   call timab(304,2,tsec)
   ABI_TIMER_STOP("KS-->QP")
 end if ! gwcalctyp>=10

 call timab(305,1,tsec) ! screening(densit
!
!=== In case update the eigenvalues ===
!* Either use a scissor operator or an external GW file.
 gw_fname = "__in.gw__"
 inquire(file=gw_fname,exist=update_energies)

 if (ABS(Ep%soenergy)>tol6) then
   write(msg,'(5a,f7.3,a)')&
&   ' screening : performing a first self-consistency',ch10,&
&   ' update of the energies in W by a scissor operator',ch10,&
&   ' applying a scissor operator of [eV] : ',Ep%soenergy*Ha_eV,ch10
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   call apply_scissor(QP_BSt,Ep%soenergy)
 else if (update_energies) then
   write(msg,'(4a)')&
&   ' screening : performing a first self-consistency',ch10,&
&   ' update of the energies in W by a previous GW calculation via GW file: ',TRIM(gw_fname)
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   ABI_ALLOCATE(igwene,(QP_Bst%mband,QP_Bst%nkpt,QP_Bst%nsppol))
   call rdgw(QP_BSt,gw_fname,igwene,extrapolate=.FALSE.)
!  call rdgw(QP_BSt,gw_fname,igwene,extrapolate=.TRUE.)
   ABI_DEALLOCATE(igwene)
   call update_occ(QP_BSt,Dtset%fixmom)
 end if

#if defined HAVE_GW_CUTOFF
!Here we should check that 0<= z0 and d<= 1 MG It wont work in case of gwpara==2
 if (.FALSE.) then
   z0 = dtset%userra
   width  = dtset%userrb
   direction = dtset%useria
   call cutoff_m_elem(Ep,Kmesh,Wfd%gvec,Wfd,KS_BSt%eig,z0,width,KS_BSt%occ,direction,gprimd)
   return
 end if
#endif
!
!========================
!=== COMPUTE DENSITY ====
!========================
!* Evaluate PW part (complete charge in case of NC pseudos)
!TODO this part has to be rewritten. If I decrease the tol on the occupations
!I have to code some MPI stuff also if use_tr==.TRUE.

 ABI_ALLOCATE(rhor,(nfftf,Dtset%nspden))
 ABI_ALLOCATE(taur,(nfftf,Dtset%nspden*Dtset%usekden))
 if (Dtset%gwgamma>0) then
   ABI_ALLOCATE(rhor_kernel,(nfftf,Dtset%nspden))
   rhor_kernel=zero
 end if

 call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,rhor)
 if(Dtset%usekden==1)call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,taur,optcalc=1)

!TODO this has to be done in a better way, moreover wont work for PAW
 call cutoff_density(ngfftf,Dtset%nspden,Wfd%nsppol,Vcp,rhor)

 call timab(305,2,tsec) ! screening(densit
!
 if (Dtset%usepaw==1) then ! Additional computation for PAW.
   call timab(320,1,tsec) ! screening(paw
!  
!  Add the compensation charge to the PW density.
   nhatgrdim=0; if (Dtset%xclevel==2) nhatgrdim=usexcnhat*Dtset%pawnhatxc
   cplex=1; ider=2*nhatgrdim; izero=0
   if (nhatgrdim>0)  then
     ABI_ALLOCATE(nhatgr,(nfftf,Dtset%nspden,3))
   end if
   call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,Cryst%gprimd,MPI_enreg_seq,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Dtset%paral_kgb,Pawang,&
&   Pawfgrtab,nhatgr,nhat,Pawrhoij,Pawrhoij,Pawtab,k0,Cryst%rprimd,Cryst%ucvol,Cryst%xred)
!  
!  === Evaluate onsite energies, potentials, densities ===
!  * Initialize variables/arrays related to the PAW spheres.
!  * Initialize also lmselect (index of non-zero LM-moments of densities).
!  TODO call init_paw_ij in scfcv, fix small issues
   cplex=1; cplex_dij=Dtset%nspinor
   ABI_ALLOCATE(Paw_ij,(Cryst%natom))
   call nullify_paw_ij(Paw_ij)

   call init_paw_ij(Paw_ij,cplex,cplex_dij,Dtset%nspinor,Wfd%nsppol,&
&   Wfd%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&   has_dij=1,has_dijhartree=1,has_exexch_pot=1,has_pawu_occ=1)

   ABI_ALLOCATE(Paw_an,(Cryst%natom))
   call nullify_paw_an(Paw_an)
   nkxc1=0
   call init_paw_an(Cryst%natom,Cryst%ntypat,nkxc1,Dtset%nspden,cplex,Dtset%pawxcdev,&
&   Cryst%typat,Pawang,Pawtab,Paw_an,has_vxc=1,has_vxcval=0)

   nzlmopt=-1; option=0; compch_sph=greatest_real
   call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,ipert,Dtset%ixc,MPI_enreg_seq,&
&   Cryst%natom,Cryst%natom,Dtset%nspden,Cryst%ntypat,nzlmopt,option,Dtset%paral_kgb,Paw_an,Paw_an,&
&   Paw_ij,Pawang,Dtset%pawprtvol,Pawrad,Pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,&
&   Dtset%xclevel,Dtset%xc_denpos,Psps%znuclpsp)
   call timab(320,2,tsec) ! screening(paw
 end if ! usepaw

 ABI_TIMER_START("screening2")
 call timab(321,1,tsec) ! screening(2)

 if (.not.allocated(nhatgr))  then
   ABI_ALLOCATE(nhatgr,(nfftf,Dtset%nspden,0))
 end if

 call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,rhor,ucvol,&
& Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,omegaplasma)

!For PAW, add the compensation charge the FFT mesh, then get rho(G).
 if (Dtset%usepaw==1) rhor(:,:)=rhor(:,:)+nhat(:,:)

 call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,rhor,ucvol=ucvol)
 if(Dtset%usekden==1)call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,taur,ucvol=ucvol,optrhor=1)

 if (Dtset%gwgamma>0) rhor_kernel = rhor

 ABI_ALLOCATE(rhog,(2,nfftf))
 ABI_ALLOCATE(taug,(2,nfftf*Dtset%usekden))
 call fourdp(1,rhog,rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
 if(Dtset%usekden==1)call fourdp(1,taug,taur(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)

!The following steps have been gathered in the setvtr routine:
!- get Ewald energy and Ewald forces
!- compute local ionic pseudopotential vpsp
!- eventually compute 3D core electron density xccc3d
!- eventually compute vxc and vhartr
!- set up ks_vtrial
!**************************************************************
!**** NOTE THAT Vxc CONTAINS THE CORE-DENSITY CONTRIBUTION ****
!**************************************************************

 ABI_ALLOCATE(grewtn,(3,Cryst%natom))
 ABI_ALLOCATE(xred_dummy,(3,Cryst%natom))
 xred_dummy=Cryst%xred
 nkxc=0
 if (Dtset%nspden==1) nkxc=2
 if (Dtset%nspden>=2) nkxc=3 ! check GGA and spinor that is messy !!!
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

 optene=4; moved_atm_inside=0; moved_rhor=0; initialized=1; istep=1
 call setvtr(Cryst%atindx1,Dtset,KS_energies,Cryst%gmet,Cryst%gprimd,grewtn,gsqcutf_eff,istep,kxc,mgfftf,&
& moved_atm_inside,moved_rhor,MPI_enreg_seq,Cryst%nattyp,nfftf,ngfftf,nhat,nhatgr,nhatgrdim,nkxc,Cryst%ntypat,&
& Psps%n1xccc,n3xccc,optene,Pawtab,ph1df,Psps,rhog,rhor,Cryst%rmet,Cryst%rprimd,strsxc,Cryst%ucvol,usexcnhat,&
& ks_vhartr,vpsp,ks_vtrial,ks_vxc,vxcavg,wvl,xccc3d,xred_dummy,taug=taug,taur=taur)
!TODO here xred is INOUT due to ionion_realSpace and xredcart!

 ABI_DEALLOCATE(grewtn)
 ABI_DEALLOCATE(xred_dummy)
 ABI_DEALLOCATE(xccc3d)
 istat = ABI_ALLOC_STAT
!
!============================
!==== Compute KS PAW Dij ====
!============================
 if (Dtset%usepaw==1) then
!  
!  Calculate unsymmetrized Dij. 
   cplex=1; ipert=0; idir=0
   call pawdij(cplex,Dtset,Dtset%enunit,one,Cryst%gprimd,ipert,MPI_enreg_seq,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,Dtset%nspden,Cryst%ntypat,&
&   Dtset%paral_kgb,Paw_an,Paw_ij,Pawang,Pawfgrtab,Dtset%pawprtvol,&
&   Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,k0,&
&   Cryst%typat,Cryst%ucvol,ks_vtrial,ks_vxc,Cryst%xred)
!  
!  Symmetrize KS Dij 
#if 0
   call symdij(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,0,Paw_ij,Pawang,&
&   Dtset%pawprtvol,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
#else
   call symdij_all(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert,Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,Paw_ij,Pawang,&
&   Dtset%pawprtvol,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
#endif
!  
!  Output of the pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   call pawprt(Dtset,Psps%indlmn,Psps%lmnmax,Paw_ij,Pawrhoij,Pawtab)
 end if
!
!=== Calculate the frequency mesh ===
!* First omega is always zero without broadening.
!FIXME what about metals? I think we should add eta, this means we need to know if the system is metallic, for example using occopt
!MS Modified to account for non-zero starting frequency (19-11-2010)
!MS Modified for tangent grid (07-01-2011)
 ABI_ALLOCATE(Ep%omega,(Ep%nomega))
 Ep%omega(1)=CMPLX(Ep%omegaermin,zero,kind=dpc)

 if (Ep%nomegaer>1) then ! Avoid division by zero.
   if (Dtset%cd_use_tangrid==0) then
     domegareal=(Ep%omegaermax-Ep%omegaermin)/(Ep%nomegaer-1)
     do iomega=2,Ep%nomegaer
       Ep%omega(iomega)=CMPLX(Ep%omegaermin+(iomega-1)*domegareal,zero,kind=dpc)
     end do
   else ! We are using tangent transformed grid
     MSG_WARNING(' EXPERIMENTAL - Using tangent transform grid for contour deformation.')
     Ep%omegaermax = Dtset%cd_max_freq
     Ep%omegaermin = zero
     ifirst=1; ilast=Ep%nomegaer
     if (Dtset%cd_subset_freq(1)/=0) then ! Only a subset of frequencies is being calculated
       ifirst=Dtset%cd_subset_freq(1); ilast=Dtset%cd_subset_freq(2)
     end if 
     factor = Dtset%cd_halfway_freq/TAN(pi*quarter)
!    *Important*-here nfreqre is used because the step is set by the original grid
     domegareal=(ATAN(Ep%omegaermax/factor)*two*piinv)/(Dtset%nfreqre-1) ! Stepsize in transformed variable
     do iomega=1,Ep%nomegaer
       Ep%omega(iomega)=CMPLX(factor*TAN((iomega+ifirst-2)*domegareal*pi*half),zero,kind=dpc)
     end do
     Ep%omegaermin = REAL(Ep%omega(1))
     Ep%omegaermax = REAL(Ep%omega(Ep%nomegaer))
   end if
 end if

 if (Ep%plasmon_pole_model.and.Ep%nomega==2) then
   e0=Dtset%ppmfrq; if (e0<0.1d-4) e0=omegaplasma
   Ep%omega(2)=CMPLX(zero,e0,kind=dpc)
 end if
!
!=== For AC, use Gauss-Legendre quadrature method ===
!* Replace $ \int_0^\infty dx f(x) $ with $ \int_0^1 dz f(1/z - 1)/z^2 $.
!* Note that the grid is not log as required by CD, so we cannot use the same SCR file.
 if (Ep%analytic_continuation) then
   ABI_ALLOCATE(z,(Ep%nomegaei))
   ABI_ALLOCATE(zw,(Ep%nomegaei))
   call coeffs_gausslegint(zero,one,z,zw,Ep%nomegaei)
   do iomega=1,Ep%nomegaei
     Ep%omega(Ep%nomegaer+iomega)=CMPLX(zero,one/z(iomega)-one,kind=dpc)
   end do
   ABI_DEALLOCATE(z)
   ABI_DEALLOCATE(zw)
 else if (Ep%contour_deformation.and.(Dtset%cd_custom_imfrqs/=0)) then
   Ep%omega(Ep%nomegaer+1)=CMPLX(zero,Dtset%cd_imfrqs(1))
   do iomega=2,Ep%nomegaei
     if (Dtset%cd_imfrqs(iomega)<=Dtset%cd_imfrqs(iomega-1)) then
       MSG_ERROR(' Specified imaginary frequencies need to be strictly increasing!')
     end if
     Ep%omega(Ep%nomegaer+iomega)=CMPLX(zero,Dtset%cd_imfrqs(iomega))
   end do
 else if (Ep%contour_deformation.and.(Ep%nomegaei/=0)) then
   e0=Dtset%ppmfrq; if (e0<0.1d-4) e0=omegaplasma
   do iomega=1,Ep%nomegaei
     Ep%omega(Ep%nomegaer+iomega)=CMPLX(zero,e0/three*(EXP(two/(Ep%nomegaei+1)*LOG(four)*iomega)-one),kind=dpc)
   end do
 end if

 if (Dtset%cd_full_grid/=0) then ! Full grid will be calculated
!  Grid values are added after the last imaginary freq.
   do ios=1,Ep%nomegaei
     do iomega=2,Ep%nomegaer
       Ep%omega(Ep%nomegaer+Ep%nomegaei+(ios-1)*(Ep%nomegaer-1)+(iomega-1)) = &
&       CMPLX(REAL(Ep%omega(iomega)),AIMAG(Ep%omega(Ep%nomegaer+ios)))
     end do
   end do
 end if
!
!Report frequency mesh for chi0.
 write(msg,'(2a)')ch10,' calculating chi0 at frequencies [eV] :'
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 do iomega=1,Ep%nomega
   write(msg,'(i3,2es16.6)')iomega,Ep%omega(iomega)*Ha_eV
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end do

!* Allocate chi0, wings and array for chi0_sumrule check.
 ABI_ALLOCATE(chi0_sumrule,(Ep%npwe))

 write(msg,'(a,f16.1,a)')' Memory required for chi0 matrix= ',2.0*gwpc*Ep%npwe**2*Ep%nI*Ep%nJ*Ep%nomega*b2Mb," [Mb]."
 call wrtout(std_out,msg,'COLL')

 ABI_ALLOCATE(chi0,(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"Out of memory in chi0")
 if (Dtset%chkgwcomp==1)  then
   ABI_ALLOCATE(deltaI,(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,1))
 end if

!=== Open file to write independent matrix elements of \epsilon^-1 ===
!if (Dtset%prtvol==10) then
 if (my_rank==master) then
   fname=dtfil%fnameabo_em1
   call isfile(fname,'new')
   open(dtfil%unem1ggp,file=fname,status='unknown',form='formatted')
 end if
!end if
!
!============================== END OF THE INITIALIZATION PART ===========================
!
!======================================================================
!==== Loop over q-points. Calculate \epsilon^{-1} and save on disc ====
!======================================================================
 ABI_TIMER_STOP("screening2")
 call timab(321,2,tsec) ! screening(2)

 ABI_TIMER_START("q-loop")

 do iqibz=1,Qmesh%nibz
!  
   call timab(306,1,tsec)
   is_first_qcalc=(iqibz==1)
!  
!  Selective q-point calculation.
   found=.FALSE.; label=iqibz
   if (Ep%nqcalc/=Ep%nqibz) then
     do iqcalc=1,Ep%nqcalc
       qtmp(:)=Qmesh%ibz(:,iqibz)-Ep%qcalc(:,iqcalc)
       found=(normv(qtmp,gmet,'G')<GW_TOLQ)
       if (found) then
         label=iqcalc; EXIT !iqcalc
       end if
     end do
     if (.not.found) CYCLE !iqibz
     qtmp(:)=Ep%qcalc(:,1)-Qmesh%ibz(:,iqibz)
     is_first_qcalc=(normv(qtmp,gmet,'G')<GW_TOLQ)
   end if

   bar=REPEAT('-',80)
   write(msg,'(4a,1x,a,i2,a,f9.6,2(",",f9.6),3a)')ch10,ch10,bar,ch10,&
&   ' q-point number ',label,'        q = (',(Qmesh%ibz(ii,iqibz),ii=1,3),') [r.l.u.]',ch10,bar
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   is_qeq0=0
   if (normv(Qmesh%ibz(:,iqibz),gmet,'G')<GW_TOLQ0) is_qeq0=1

!  === Find the independent set of G-Gp pairs for this q-point. ===
!  Useful to write the independent matrix elements of epsilon^-1 or to speed-up the KK transform
!  In the long wavelength limit we set q==0, because we still can use symmetryes for the Body.
   if (Dtset%prtvol>=10) then
     qtmp(:)=Qmesh%ibz(:,iqibz); if (normv(qtmp,gmet,'G')<GW_TOLQ0) qtmp(:)=zero
     call init_Gpairs_type(Gpairs_q,qtmp,Gsph_epsG0,Cryst)
   end if
   call timab(306,2,tsec)

   if (is_qeq0==1) then  ! Special treatment of the long wavelenght limit.
     call timab(307,1,tsec)

     ABI_ALLOCATE(chi0_lwing,(Ep%npwe*Ep%nI,Ep%nomega,3))
     ABI_ALLOCATE(chi0_uwing,(Ep%npwe*Ep%nJ,Ep%nomega,3))
     ABI_ALLOCATE(chi0_head,(3,3,Ep%nomega))

     call cchi0q0(use_tr,Dtset,Cryst,Ep,Psps,Kmesh,QP_BSt,KS_BSt,Gsph_epsG0,Gsph_wfn,&
&     Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,Pawfgrtab,Paw_onsite,ktabr,ktabrf,nbvw,ngfft_gw,nfftgw,&
&     ngfftf,nfftf_tot,chi0,chi0_head,chi0_lwing,chi0_uwing,Ltg_q(iqibz),chi0_sumrule,Wfd,Wfdf)
!    
!    Add the intraband term if required and metallic occupation scheme is used.
!    add_chi0_intraband=.TRUE.
     add_chi0_intraband=.FALSE.
     if (add_chi0_intraband .and. bst_ismetal(QP_BSt)) then

       ABI_ALLOCATE(chi0intra,(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega))
       istat = ABI_ALLOC_STAT
       ABI_CHECK(istat==0,"Out of memory in chi0intra")

       ABI_ALLOCATE(chi0intra_lwing,(Ep%npwe*Ep%nI,Ep%nomega,3))
       ABI_ALLOCATE(chi0intra_uwing,(Ep%npwe*Ep%nJ,Ep%nomega,3))
       ABI_ALLOCATE(chi0intra_head,(3,3,Ep%nomega))

!      $call wfk_read_ene(wfk_fname,Dtset%accesswff,localrdwf,inbnds,energies_p,occ_p,Hdr,prtvol,comm,kg_p) 
!      $call wfd_read_wfk(Wfd,wfk_fname,Dtset%accesswff)

       call chi0q0_intraband(Wfd,Cryst,Ep,Psps,QP_BSt,Gsph_epsG0,Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,use_tr,Dtset%usepawu,&
&       ngfft_gw,chi0intra,chi0intra_head,chi0intra_lwing,chi0intra_uwing)

       call wrtout(std_out,"Head of chi0 and chi0_intra","COLL")
       do iomega=1,Ep%nomega
         write(std_out,*)Ep%omega(iomega)*Ha_eV,REAL(chi0(1,1,iomega)),REAL(chi0intra(1,1,iomega))
         write(std_out,*)Ep%omega(iomega)*Ha_eV,AIMAG(chi0(1,1,iomega)),AIMAG(chi0intra(1,1,iomega))
       end do

       chi0       = chi0       + chi0intra   
       chi0_head  = chi0_head  + chi0intra_head  
       chi0_lwing = chi0_lwing + chi0intra_lwing
       chi0_uwing = chi0_uwing + chi0intra_uwing

       ABI_DEALLOCATE(chi0intra)
       ABI_DEALLOCATE(chi0intra_lwing)
       ABI_DEALLOCATE(chi0intra_uwing)
       ABI_DEALLOCATE(chi0intra_head)
     end if

     if (.FALSE.) then
       write(std_out,*)"head of chi0"
       do iomega=1,Ep%nomega
         call print_arr(chi0_head(:,:,iomega),max_r=3,max_c=3)
       end do

       ABI_ALLOCATE(wtest,(Ep%npwe))
       do iomega=1,Ep%nomega
         write(std_out,*)"iomega=",iomega

         write(std_out,*)"symmetrized e_00 via tensor"
         wng = MATMUL(chi0_head(:,:,iomega),GW_Q0_DEFAULT)
         write(std_out,*) one - vdotw(GW_Q0_DEFAULT,wng,Cryst%gmet,"G") * Vcp%vcqlwl_sqrt(1,1)*Vcp%vcqlwl_sqrt(1,1)

         write(std_out,*)"symmetrized e_0G via tensor"
         do ii=1,Ep%npwe
           wng = chi0_uwing(ii,iomega,:)
           wtest(ii) = - vdotw(GW_Q0_DEFAULT,wng,Cryst%gmet,"G") * Vcp%vcqlwl_sqrt(1,1) * Vcp%vcqlwl_sqrt(ii,1)
         end do
         call print_arr(wtest,max_r=9,unit=std_out)

         write(std_out,*)"symmetrized e_G0 via tensor"
         do ii=1,Ep%npwe
           wng = chi0_lwing(ii,iomega,:)
           wtest(ii) = - vdotw(GW_Q0_DEFAULT,wng,Cryst%gmet,"G") * Vcp%vcqlwl_sqrt(1,1) * Vcp%vcqlwl_sqrt(ii,1)
         end do
         call print_arr(wtest,max_r=9,unit=std_out)
       end do

       ABI_ALLOCATE(eps_head,(3,3,Ep%nomega))
       call mkem1_q0(Ep%npwe,1,1,Ep%nomega,Cryst,Vcp,Gsph_epsG0%gvec,chi0_head,chi0_lwing,chi0_uwing,chi0,eps_head)
       ABI_DEALLOCATE(eps_head)

       length = normv(GW_Q0_DEFAULT,Cryst%gmet,"G")
       do iomega=1,Ep%nomega
         write(std_out,*)"iomega=",iomega

         em1_00 = one / vdotw(GW_Q0_DEFAULT/length, MATMUL(chi0_head(:,:,iomega),GW_Q0_DEFAULT/length),Cryst%gmet,"G")
         write(std_out,*)"e^1_{00} from tensor"
         write(std_out,*) em1_00

         write(std_out,*)"symmetrized e^-1_0G via tensor"
         do ii=1,Ep%npwe
           wng = chi0_uwing(ii,iomega,:)
           wtest(ii) = em1_00*vdotw(GW_Q0_DEFAULT/length,wng,Cryst%gmet,"G")
         end do
         wtest(1) = em1_00
         call print_arr(wtest,max_r=9,unit=std_out)

         write(std_out,*)"symmetrized e^-1_G0 via tensor"
         do ii=1,Ep%npwe
           wng = chi0_lwing(ii,iomega,:)
           wtest(ii) = em1_00*vdotw(GW_Q0_DEFAULT/length,wng,Cryst%gmet,"G")
         end do
         wtest(1) = em1_00
         call print_arr(wtest,max_r=9,unit=std_out)
       end do !iomega
       ABI_DEALLOCATE(wtest)
     end if

     if (Dtset%chkgwcomp==1) then
       call check_completeness(use_tr,Dtset,Cryst,Qmesh%ibz(:,iqibz),Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,&
&       Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,ngfft_gw,nfftgw,ngfftf,nfftf_tot,deltaI,ktabr,ktabrf,Ltg_q(iqibz),Wfd,Wfdf)
     end if

     call timab(307,2,tsec)

   else ! Calculate cchi0 for q/=0.

     call timab(308,1,tsec)

     call cchi0(use_tr,Dtset,Cryst,Qmesh%ibz(:,iqibz),Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,Gsph_wfn,&
&     Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,nbvw,ngfft_gw,nfftgw,ngfftf,nfftf_tot,chi0,ktabr,ktabrf,&
&     Ltg_q(iqibz),chi0_sumrule,Wfd,Wfdf)

     if (Dtset%chkgwcomp==1) then
       call check_completeness(use_tr,Dtset,Cryst,Qmesh%ibz(:,iqibz),Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,&
&       Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,ngfft_gw,nfftgw,ngfftf,nfftf_tot,deltaI,ktabr,ktabrf,Ltg_q(iqibz),Wfd,Wfdf)
     end if

     call timab(308,2,tsec)

   end if
!  
!  ==== Print chi0(q,G,Gp,omega), then calculate epsilon and epsilon^-1 for this q ====
!  * Only master works but this part could be parallelized over frequencies
   ABI_TIMER_START("build_W")
   call timab(309,1,tsec)

   do iomega=1,MIN(Ep%nomega,NOMEGA_PRINTED)
     write(msg,'(1x,a,i4,a,2f9.4,a)')' chi0(G,G'') at the ',iomega,' th omega',Ep%omega(iomega)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     write(msg,'(1x,a,i3,a,i4,a)')' chi0(q =',iqibz, ', omega =',iomega,', G,G'')'
     if (Ep%nqcalc/=Ep%nqibz) write(msg,'(a,i3,a,i4,a)')'  chi0(q=',iqcalc,', omega=',iomega,', G,G'')'
     call wrtout(std_out,msg,'COLL')
!    arr99 is needed to avoid the update of all the tests. Now chi0 is divided by ucvol inside (cchi0|cchi0q0).
!    TODO should be removed but GW tests have to be updated.
     ii = MIN(9,Ep%npwe)
     ABI_ALLOCATE(arr_99,(ii,ii))
     arr_99 = chi0(1:ii,1:ii,iomega)*ucvol
     call print_arr(arr_99,max_r=2,unit=ab_out)
     call print_arr(arr_99,unit=std_out)
     ABI_DEALLOCATE(arr_99)
!    call print_arr(chi0(:,:,iomega),max_r=2,unit=ab_out)
!    call print_arr(chi0(:,:,iomega),unit=std_out)
   end do

   if (Ep%nomega>NOMEGA_PRINTED) then
     write(msg,'(a,i3,a)')' No. of calculated frequencies > ',NOMEGA_PRINTED,', stop printing '
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if
!  
!  === Write chi0 on _SUSC file ===
!  * Master creates and write the header if this is the first q-point calculated.
   if (my_rank==master) then
     title(1)='CHI0 file: chi0'
     title(2)=' '
     if (is_qeq0==1) then
       string='0'; if (Dtset%usepaw==0.and.Ep%inclvkb/=0) call int2char(Ep%inclvkb,string)
       title(1)=title(1)(1:21)//', calculated using inclvkb = '//string
     end if
     unt_susc=Dtfil%unchi0
!    
!    * Open file and write header for polarizability files.
     if (is_first_qcalc) then
       open(unit=unt_susc,file=dtfil%fnameabo_sus,status='unknown',form='unformatted',iostat=ios)
       msg = ' Opening file '//TRIM(dtfil%fnameabo_sus)//' as new-unformatted'
       ABI_CHECK(ios==0,msg)
       fileID=1; ikxc=0; test_type=0; tordering=1
       call init_ScrHdr(fileID,ikxc,test_type,tordering,title,Ep%npwe,Gsph_epsG0%gvec,Ep,Hdr_local,Hchi0)
       rdwr=2; fform_chi0=1102 ! Use the new format
       call scr_hdr_io(fform_chi0,rdwr,unt_susc,xmpi_self,0,Dtset%accesswff,Hchi0)
       call free_scrhdr(Hchi0)
     end if
     call write_screening(unt_susc,Dtset%accesswff,Ep%npwe,Ep%nomega,chi0)
   end if
!  
!  === Write deltaI on _DELI file ===
   if (Dtset%chkgwcomp==1 .and. my_rank==master) then
     unt_delI=Dtfil%undelI
!    
!    * Open file and write header for deltaI files.
     if (is_first_qcalc) then
       open(unit=unt_delI,file=dtfil%fnameabo_delI,status='unknown',form='formatted',iostat=ios)
       msg = ' Opening file '//TRIM(dtfil%fnameabo_delI)//' as formatted'
       ABI_CHECK(ios==0,msg)
       title(1)='deltaI file: incompleteness matrix'
       title(2)=' '
       write(unt_delI,'(a)') title(1)
     end if
     call write_deltaI(unt_delI,iqibz,deltaI,Ep%npwe)
   end if

!  Calculate the RPA functional correlation energy if the polarizability on a
!  Gauss-Legendre mesh along imaginary axis is available
   if ( Ep%analytic_continuation .and. Dtset%gwrpacorr>0 ) then
     if( is_first_qcalc ) then
       ABI_ALLOCATE(ec_rpa,(Dtset%gwrpacorr))
       ec_rpa(:)=zero
     end if
     call calc_rpa_functional(Dtset%gwrpacorr,label,iqibz,Ep,Vcp,Qmesh,Dtfil,gmet,chi0,comm,ec_rpa)
     if(label==Ep%nqcalc)  then
       ABI_DEALLOCATE(ec_rpa)
     end if
   end if
!  
!  ==========================================================
!  === Calculate RPA \tilde\epsilon^{-1} overwriting chi0 ===
!  ==========================================================
   approx_type=0 !RPA
   option_test=0 !TESTPARTICLE
   dim_wing=0; if (is_qeq0==1) dim_wing=3

   if (dim_wing==0) then
     if (.not.allocated(chi0_lwing))  then
       ABI_ALLOCATE(chi0_lwing,(Ep%npwe*Ep%nI,Ep%nomega,dim_wing))
     end if
     if (.not.allocated(chi0_uwing))  then
       ABI_ALLOCATE(chi0_uwing,(Ep%npwe*Ep%nJ,Ep%nomega,dim_wing))
     end if
     if (.not.allocated(chi0_head ))  then
       ABI_ALLOCATE(chi0_head,(dim_wing,dim_wing,Ep%nomega))
     end if
   end if

#if 0
!  Using the random q for the optical limit is one of the reasons
!  why sigma breaks the initial energy degeneracies.
   chi0_lwing=czero
   chi0_uwing=czero
   chi0_head=czero
#endif

!  If the vertex is being included for the spectrum, calculate the kernel now and pass it on
   if (Dtset%gwgamma==0) then
     approx_type=0; option_test=0; dim_kxcg=0
     ABI_ALLOCATE(kxcg,(Ep%npwe,dim_kxcg))

   else if (Dtset%gwgamma==1.OR.Dtset%gwgamma==2) then ! ALDA TDDFT kernel vertex
     ABI_CHECK(Dtset%usepaw==0,"GWGamma + PAW not available")
     MSG_WARNING('EXPERIMENTAL: Kernel is being added to screening, the SCR file will be non-standard!!')
     ikxc=7; approx_type=1; dim_kxcg=1
     if (Dtset%gwgamma==1) option_test=1 ! TESTELECTRON, vertex in chi0 *and* sigma
     if (Dtset%gwgamma==2) option_test=0 ! TESTPARTICLE, vertex in chi0 only
     ABI_ALLOCATE(kxcg,(Ep%npwe,dim_kxcg))
     call xc_kernel(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,rhor_kernel,&
&     Ep%npwe,dim_kxcg,kxcg,Gsph_epsG0%gvec,xmpi_self)

   else if (Dtset%gwgamma==3.OR.Dtset%gwgamma==4) then ! ADA non-local kernel vertex
     ABI_CHECK(Wfd%usepaw==0,"ADA vertex + PAW not available")
     ABI_CHECK(Wfd%nsppol==1,"ADA vertex for GWGamma not available yet for spin-polarised cases")
     MSG_WARNING('EXPERIMENTAL: Kernel is being added to screening, the SCR file will be non-standard!!')
     ikxc=7; approx_type=2
     if (Dtset%gwgamma==3) option_test=1 ! TESTELECTRON, vertex in chi0 *and* sigma
     if (Dtset%gwgamma==4) option_test=0 ! TESTPARTICLE, vertex in chi0 only
     ABI_ALLOCATE(fxc_ADA,(Ep%npwe,Ep%npwe,Ep%nqibz))
!    Use userrd to set kappa
     if (Dtset%userrd==zero) Dtset%userrd = 2.1_dp
!    Set correct value of kappa (should be scaled with alpha*r_s where)
!    r_s is Wigner-Seitz radius and alpha=(4/(9*Pi))^(1/3)
     rhoav = (omegaplasma*omegaplasma)/four_pi
     r_s = (three/(four_pi*rhoav))**third 
     alpha = (four*ninth*piinv)**third
     Dtset%userrd = Dtset%userrd*alpha*r_s

     call xc_kernel_ADA(Dtset,Cryst,ikxc,ngfftf,nfftf,Wfd%nspden,&
&     rhor_kernel,Ep%npwe,Ep%nqibz,Ep%qibz,&
&     fxc_ADA,Gsph_epsG0%gvec,xmpi_self,kappa_init=Dtset%userrd)
   end if

   if (allocated(rhor_kernel))  then
     ABI_DEALLOCATE(rhor_kernel)
   end if

   if (approx_type<2) then
     call make_epsm1_driver(iqibz,dim_wing,Ep%npwe,Ep%nI,Ep%nJ,Ep%nomega,Ep%omega,&
&     approx_type,option_test,Vcp,nfftf_tot,ngfftf,dim_kxcg,kxcg,Gsph_epsG0%gvec,&
&     chi0_head,chi0_lwing,chi0_uwing,chi0,Spectra,comm)
   else
     call make_epsm1_driver(iqibz,dim_wing,Ep%npwe,Ep%nI,Ep%nJ,Ep%nomega,Ep%omega,&
&     approx_type,option_test,Vcp,nfftf_tot,ngfftf,dim_kxcg,kxcg,Gsph_epsG0%gvec,&
&     chi0_head,chi0_lwing,chi0_uwing,chi0,Spectra,comm,fxc_ADA(:,:,iqibz))
   end if

   ABI_DEALLOCATE(chi0_lwing)
   ABI_DEALLOCATE(chi0_uwing)
   ABI_DEALLOCATE(chi0_head)

   if (my_rank==master .and. is_qeq0==1) then
     call repr_dielconst(Spectra,msg)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')

     if (Ep%nomegaer>2) then
       call dump_spectra(Spectra,W_EELF  ,Dtfil%fnameabo_eelf)
       call dump_spectra(Spectra,W_EM_LF ,Dtfil%fnameabo_em1_lf)
       call dump_spectra(Spectra,W_EM_NLF,Dtfil%fnameabo_em1_nlf)
     end if
   end if ! master and is_qeq0==1

   call destroy_spectra(Spectra)
   if (allocated(kxcg))  then
     ABI_DEALLOCATE(kxcg)
   end if
   if (allocated(fxc_ADA))  then
     ABI_DEALLOCATE(fxc_ADA)
   end if

   epsm1   => chi0
   vc_sqrt => Vcp%vc_sqrt(:,iqibz)  ! Contains vc^{1/2}(q,G), complex-valued due to a possible cutoff
!  
!  Output the sum rule evaluation.
   call output_chi0sumrule((is_qeq0==1),iqibz,Ep%npwe,omegaplasma,chi0_sumrule,epsm1(:,:,1),vc_sqrt)
!  
!  Write heads and wings on the main output.
   if (is_qeq0==1) then
     write(msg,'(1x,2a)')' Heads and wings of the symmetrical epsilon^-1(G,G'') ',ch10
     call wrtout(ab_out,msg,'COLL')
     do iomega=1,Ep%nomega
       write(msg,'(2x,a,i4,a,2f9.4,a)')&
&       ' Upper and lower wings at the ',iomega,' th omega',Ep%omega(iomega)*Ha_eV,' [eV]'
       call wrtout(ab_out,msg,'COLL')
       call print_arr(epsm1(1,:,iomega),max_r=9,unit=ab_out)
       call print_arr(epsm1(:,1,iomega),max_r=9,unit=ab_out)
       call wrtout(ab_out,ch10,'COLL')
     end do
   end if

   call timab(309,2,tsec)
   ABI_TIMER_STOP("build_W")


   call timab(310,1,tsec) ! wrscr
   ABI_TIMER_START("write_W")

   if (my_rank==master) then

     if (Dtset%prtvol>=10) then ! Write the independent matrix elements of \tilde epsilon^-1.
       write(msg,'(a,3(f10.6),a)')' Symmetrical epsilon^-1(G,G'') at q = ( ',Qmesh%ibz(:,iqibz),' ) [r.l.u.]'
       call outeps(Ep%npwe,Ep%nomega,Ep%omega,epsm1,Gsph_epsG0,Gpairs_q,msg,dtfil%unem1ggp,Dtset%prtvol)
     end if
!    
!    === Write the symmetrical epsilon^-1 on file ===
!    This might be parallelized but I have to use xsum_mpi in cchi0 and cchi0q0
     title(1)='SCR file: epsilon^-1'
     if (is_qeq0==1) then
       string='0'; if (Dtset%usepaw==0.and.Ep%inclvkb/=0) call int2char(Ep%inclvkb,string)
       title(1)=title(1)(1:21)//', calculated using inclvkb = '//string
     end if
     title(2)='TESTPARTICLE'
     ctype='RPA'
     title(2)(14:17)=ctype !this has to be modified
     unt_em1=Dtfil%unscr

     if (is_first_qcalc) then
!      === Open file and write the header for the SCR file ===
!      * Here we write the RPA approximation for \tilde\epsilon^{-1}
       open(unit=unt_em1,file=dtfil%fnameabo_scr,status='unknown',form='unformatted',iostat=ios)
       msg = ' Opening file '//TRIM(dtfil%fnameabo_scr)//' as new-unformatted '
       ABI_CHECK(ios==0,msg)
       fileID=4; ikxc=0; test_type=0; tordering=1
       call init_ScrHdr(fileID,ikxc,test_type,tordering,title,Ep%npwe,Gsph_epsG0%gvec,Ep,Hdr_local,Hem1)
!      here we still use the old fform
       rdwr=2; fform_em1=1002
       call scr_hdr_io(fform_em1,rdwr,unt_em1,xmpi_self,0,Dtset%accesswff,Hem1)
       call free_scrhdr(Hem1)
     end if

     call write_screening(unt_em1,Dtset%accesswff,Ep%npwe,Ep%nomega,epsm1)

   end if ! my_rank==master

   ABI_TIMER_STOP("write_W")
   call timab(310,2,tsec)
 end do ! Loop over q-points

 ABI_TIMER_STOP("q-loop")
!
!----------------------------- END OF THE CALCULATION ------------------------
!
!=== Close Files ===
!if (Dtset%prtvol==10 .and. my_rank==0) close(dtfil%unem1ggp)
 if (my_rank==master) then
   close(dtfil%unem1ggp)
   close(unt_em1)
   close(unt_susc)
   if (Dtset%chkgwcomp==1) close(unt_delI)
 end if
!
!=====================
!==== Free memory ====
!=====================
 ABI_DEALLOCATE(chi0_sumrule)
 ABI_DEALLOCATE(chi0)
 if (Dtset%chkgwcomp==1)  then
   ABI_DEALLOCATE(deltaI)
 end if

 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(ks_vbik)
 ABI_DEALLOCATE(qp_vbik)
 ABI_DEALLOCATE(ktabr)
 ABI_DEALLOCATE(taur)
 ABI_DEALLOCATE(taug)
 ABI_DEALLOCATE(ks_vhartr)
 ABI_DEALLOCATE(ks_vtrial)
 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(ks_vxc)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)

 ABI_DEALLOCATE(Pawfgr%fintocoa)
 ABI_DEALLOCATE(Pawfgr%coatofin)
 ABI_DEALLOCATE(nhatgr)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(nhat)
 istat = ABI_ALLOC_STAT

 if (Dtset%usepaw==1) then ! Optional deallocation for PAW.
   call rhoij_free(Pawrhoij)
   ABI_DEALLOCATE(Pawrhoij)
   call pawfgrtab_free(Pawfgrtab)
   ABI_DEALLOCATE(Pawfgrtab)
   call destroy_paw_ij(Paw_ij)
   ABI_DEALLOCATE(Paw_ij)
   call destroy_paw_an(Paw_an)
   ABI_DEALLOCATE(Paw_an)
   call destroy_paw_pwff(Paw_pwff)
   ABI_DEALLOCATE(Paw_pwff)
   if (Dtset%pawcross==1) then
     call destroy_paw_pwaves_lmn(Paw_onsite)
     ABI_DEALLOCATE(Paw_onsite)
     ABI_DEALLOCATE(ktabrf)
     Wfdf%bks_comm = xmpi_comm_null
     call wfd_destroy(Wfdf)
   end if
 end if

 call wfd_destroy(Wfd)
 call destroy_BZ_mesh_type(Kmesh)
 call destroy_BZ_mesh_type(Qmesh)
 call destroy_crystal(Cryst)
 call destroy_gsphere(Gsph_epsG0)
 call destroy_gsphere(Gsph_wfn)
 call destroy_Gpairs_type(Gpairs_q)
 call vcoul_free(Vcp)
 call destroy_epsilonm1_parameters(Ep)
 call hdr_clean(Hdr_kss)
 call hdr_clean(Hdr_local)
 call bstruct_clean(KS_BSt)
 call bstruct_clean(QP_BSt)
 call destroy_little_group(Ltg_q)
 ABI_DEALLOCATE(Ltg_q)

 ABI_TIMER_STOP("")
 call timab(301,2,tsec)

 DBG_EXIT("COLL")

end subroutine screening
!!***
