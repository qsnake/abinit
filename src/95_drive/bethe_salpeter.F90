!{\src2tex{textfont=tt}}
!!****f* ABINIT/bethe_salpeter
!! NAME
!!  bethe_salpeter
!!
!! FUNCTION
!!  Main routine to calculate dielectric properties by solving the Bethe-Salpeter equation in
!!  Frequency-Reciprocal space on a transition (electron-hole) basis set.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (M.Giantomassi, L. Reining, V. Olevano, F. Sottile, S. Albrecht, G. Onida)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=Length scales of primitive translations (bohr)
!! codvsn=Code version
!! Dtfil<datafiles_type>=Variables related to files.
!! Dtset<dataset_type>=All input variables for this dataset.
!! Pawang<pawang_type)>=PAW angular mesh and related data.
!! Pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! Pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! Psps<pseudopotential_type>=Variables related to pseudopotentials.
!!   Before entering the first time in the routine, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!   the call to pspini. The next time the code enters bethe_salpeter, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=Dimensionless real space primitive translations.
!! xred(3,natom)=Reduced atomic coordinates.
!!
!! Input files used during the calculation.
!!  KSS        : Kohn Sham electronic structure file.
!!  SCR (SUSC) : Files containing the symmetrized epsilon^-1 or the irreducible RPA polarizability,
!!               respectively. Used to construct the screening W.
!!  GW file    : Optional file with the GW QP corrections.
!!
!! OUTPUT
!!  Output is written on the main output file and on the following external files:
!!   * _RPA_NLF_MDF: macroscopic RPA dielectric function without non-local field effects.
!!   * _GW_NLF_MDF: macroscopic RPA dielectric function without non-local field effects calculated
!!                 with GW energies or the scissors operator.
!!   * _EXC_MDF: macroscopic dielectric function with excitonic effects obtained by solving the
!!              Bethe-Salpeter problem at different level of sophistication.
!!
!! PARENTS
!!      driver
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
!!      bstruct_clean,chkpawovlp,denfgr,destroy_bs_parameters
!!      destroy_bz_mesh_type,destroy_crystal,destroy_gsphere,destroy_hur
!!      destroy_paw_an,destroy_paw_ij,destroy_paw_pwff,energies_init
!!      exc_build_ham,exc_den,exc_diago_driver,exc_haydock_driver
!!      exc_interp_ham,fourdp,get_gftt,getph,hdr_clean,init_paw_an,init_paw_ij
!!      init_paw_pwff,init_pawfgr,initmpi_seq,make_hur_commutator,metric,mkrdim
!!      nhatgrid,nullify_hur,nullify_paw_an,nullify_paw_ij,pawdenpot,pawdij
!!      pawfgrtab_free,pawfgrtab_init,pawinit,pawmknhat,pawnabla_init,pawprt
!!      pawpuxinit,print_ngfft,print_pawtab,print_psps,prtrhomxmn,pspini,rdqps
!!      rhoij_alloc,rhoij_copy,rhoij_free,rotate_fft_mesh,screen_free
!!      screen_init,screen_nullify,setsymrhoij,setup_bse,setvtr,symdij
!!      test_charge,timab,update_occ,vcoul_free,wfd_destroy,wfd_init,wfd_mkrho
!!      wfd_print,wfd_read_kss,wfd_read_wfk,wfd_reset_ur_cprj,wfd_rotate
!!      wfd_test_ortho,wfd_wave_free,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine bethe_salpeter(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,xred)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_errors
 use m_timer

 use m_header,          only : hdr_clean
 use m_fft_mesh,        only : rotate_FFT_mesh, get_gftt, print_ngfft
 use m_crystal,         only : crystal_structure, destroy_crystal
 use m_bz_mesh,         only : bz_mesh_type, destroy_bz_mesh_type
 use m_ebands,          only : update_occ, reportgap, bstruct_clean, copy_bandstructure, update_occ, get_valence_idx
 use m_gsphere,         only : gvectors_type, destroy_gsphere
 use m_io_kss,          only : wfd_read_kss
 use m_vcoul,           only : vcoul_t, vcoul_free, cutoff_density
 use m_qparticles,      only : rdqps !, show_QP , rdgw
 use m_paw_dmft,        only : paw_dmft_type
 use m_paw_toolbox,     only : nullify_paw_ij, init_paw_ij, destroy_paw_ij, init_pawfgr, pawfgrtab_free, pawfgrtab_init,&
&                              nullify_paw_an, init_paw_an, destroy_paw_an, print_pawtab, print_paw_ij
 use m_paw_commutator,  only : HUr_commutator, destroy_Hur, nullify_Hur, make_Hur_commutator
 use m_paw_pwij,        only : paw_pwff_type, init_paw_pwff, destroy_paw_pwff
 use m_wfs,             only : wfd_init, wfd_destroy, wfd_print, wfs_descriptor, wfd_test_ortho, &
&                              wfd_read_wfk, wfd_wave_free, wfd_rotate, wfd_reset_ur_cprj
 use m_gwdefs,          only : gw_uses_wfk_file
 use m_energies,        only : energies_type, energies_init

 use m_screen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bethe_salpeter'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_69_wfdesc
 use interfaces_71_bse
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=6),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: Dtfil
 type(dataset_type),intent(inout) :: Dtset
 type(pawang_type),intent(inout) :: Pawang
 type(pseudopotential_type),intent(inout) :: Psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,Dtset%natom)
 type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0,level=40,ipert0=0,idir0=0,cplex1=1
 integer,save :: nsym_old=-1
 integer :: ii,band,spin,ik_ibz,mqmem
 integer :: has_dijU,has_dijso
 integer :: ik_bz,mband
 integer :: accessfil,choice,cplex_dij
 integer :: iat,ider !,ierr
 integer :: usexcnhat,nfft_osc,mgfft_osc
 integer :: istat,isym,itypat,izero
 integer :: optcut,optgr0,optgr1,optgr2,option,optrad,optrhoij,psp_gencond
 integer :: nhatgrdim,nkxc1,nprocs,nspden_rhoij,nzlmopt,ifft
 integer :: my_rank,rhoxsp_method,comm,master
 integer :: mgfftf,spin_opt,which_fixed
 integer :: my_spin,nscf,nbsc,nkxc,n3xccc
 integer :: ndij,nfftf,nfftf_tot,nfftot_osc,my_minb,my_maxb
 integer :: optene,moved_atm_inside,moved_rhor,initialized,istep
 real(dp) :: nelect_kss,ucvol,drude_plsmf,ecore,ecut_eff,ecutdg_eff,opt_ecut
 real(dp) :: gsqcutc_eff,gsqcutf_eff
 real(dp) :: compch_fft,compch_sph,diecut_eff_dum,gsq_osc
 real(dp) :: vxcavg !,vxcavg_qp
 logical :: iscompatibleFFT,paw_add_onsite
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname,w_fname
 type(Pawfgr_type) :: Pawfgr
 type(excfiles) :: BS_files
 type(excparam) :: BSp
 type(paw_dmft_type) :: Paw_dmft
 type(MPI_type) :: MPI_enreg_seq
 type(Crystal_structure) :: Cryst
 type(BZ_mesh_type) :: Kmesh,Qmesh
 type(Gvectors_type) :: Gsph_Max,Gsph_c
 type(Hdr_type) :: Hdr_kss,Hdr_bse
 type(Bandstructure_type) :: KS_BSt,QP_BSt
 type(Energies_type) :: KS_energies
 type(vcoul_t) :: Vcp
 type(wfs_descriptor) :: Wfd
 type(screen_t) :: W
 type(screen_info_t) :: W_info
 type(wvl_internal_type) :: wvl
!arrays
 integer,save :: paw_gencond(6)=(/-1,-1,-1,-1,-1,-1/)
 integer :: ngfft_osc(18),ngfftc(18),ngfftf(18),nrcell(3)
 integer,allocatable :: ktabr(:,:),l_size_atm(:),nlmn_type(:)
 integer,allocatable :: nband(:,:),nq_spl(:),irottb(:,:)
 integer,allocatable :: qp_vbik(:,:)
 integer,allocatable :: gfft_osc(:,:)
 real(dp),parameter :: k0(3)=zero
 real(dp) :: tsec(2),gmet(3,3),gprimd(3,3),qphon(3),rmet(3,3),rprimd(3,3),eh_rcoord(3),strsxc(6)
 real(dp),allocatable :: ph1df(:,:),prev_rhor(:,:),ph1d(:,:)!,qp_nhat(:,:)
 real(dp),allocatable :: ks_nhat(:,:),ks_nhatgr(:,:,:),ks_rhog(:,:),ks_rhor(:,:),qp_aerhor(:,:)
 real(dp),allocatable :: qp_rhor(:,:),qp_rhog(:,:) !,qp_vhartr(:),qp_vtrial(:,:),qp_vxc(:,:)
 real(dp),allocatable :: grewtn(:,:),qmax(:)
 real(dp),allocatable :: vpsp(:),xccc3d(:),xred_dummy(:,:)
 real(dp),allocatable :: ks_vhartr(:),ks_vtrial(:,:),ks_vxc(:,:)
 real(dp),allocatable :: kxc(:,:) !,qp_kxc(:,:)
 complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(Pawrhoij_type),allocatable :: KS_Pawrhoij(:)
 type(Pawrhoij_type),allocatable :: prev_Pawrhoij(:) !QP_pawrhoij(:),
 type(Paw_pwff_type),allocatable :: Paw_pwff(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(HUr_commutator),allocatable :: Hur(:)
 type(Paw_ij_type),allocatable :: KS_paw_ij(:)
 type(Paw_an_type),allocatable :: KS_paw_an(:)

!************************************************************************

 DBG_ENTER('COLL')

 ABI_TIMER_START("")
 call timab(650,1,tsec) ! bse(Total)

 ABI_TIMER_START("init1")
 call timab(651,1,tsec) ! bse(Init1)

 write(msg,'(8a)')&
& ' Exciton: Calculation of dielectric properties by solving the Bethe-Salpeter equation ',ch10,&
& ' in frequency domain and reciprocal space on a transitions basis set. ',ch10,&
& ' Based on a program developed by L. Reining, V. Olevano, F. Sottile, ',ch10,&
& ' S. Albrecht, and G. Onida. Incorporated in ABINIT by M. Giantomassi. ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

#ifdef HAVE_GW_DPC
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

 comm = xmpi_world
 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master=0
!
!* Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
!
!* accesswff defines the format of the output.
!1--> Plain Fortran file
!2--> Set all outputs to netcdf format (not implemented)
!3--> Set all outputs to ETSF format
!
 accessfil=0
 if (Dtset%accesswff==IO_MODE_NETCDF) accessfil=1
 if (Dtset%accesswff==IO_MODE_ETSF  ) accessfil=3
 if (Dtset%accesswff==IO_MODE_MPI   ) accessfil=4

!===================================================
!=== Initialize names for input and output files ===
!===================================================

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

!TODO Recheck getng, should use same trick as that used in screening and sigma.

 call init_pawfgr(Dtset,Pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
& gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=gmet,k0=k0)

 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')
 nfftf_tot=PRODUCT(ngfftf(1:3))
!
!===========================================
!=== Open and read pseudopotential files ===
!===========================================
 call pspini(Dtset,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,level,Pawrad,Pawtab,Psps,rprimd)
 if (psp_gencond==1) call print_psps(Psps,std_out,0,'COLL')

!=== Initialization of basic objects including the BSp structure that defines the parameters of the run ===
!call timab(652,1,tsec) ! setup_bse

 call setup_bse(codvsn,acell,rprim,ngfftf,ngfft_osc,Dtset,Dtfil,BS_files,Psps,Pawtab,BSp,&
& Cryst,Kmesh,Qmesh,KS_BSt,QP_BSt,Hdr_kss,Gsph_Max,Gsph_c,Vcp,Hdr_bse,w_fname,comm,wvl)

!call timab(652,2,tsec) ! setup_bse

 nfftot_osc=PRODUCT(ngfft_osc(1:3))
 nfft_osc  =nfftot_osc  !no FFT //
 mgfft_osc =MAXVAL(ngfft_osc(1:3))

 call print_ngfft(ngfft_osc,header='FFT mesh used for oscillator strengths')

!TRYING TO RECREATE AN "ABINIT ENVIRONMENT"
 KS_energies%e_corepsp=ecore/Cryst%ucvol
!
!============================
!==== PAW initialization ====
!============================
 if (Dtset%usepaw==1) then

   call chkpawovlp(Cryst%natom,Cryst%ntypat,Dtset%pawovlp,Pawtab,Cryst%rmet,Cryst%typat,xred)

   ABI_ALLOCATE(nlmn_type,(Cryst%ntypat))
   do itypat=1,Cryst%ntypat
     nlmn_type(itypat)=Pawtab(itypat)%lmn_size
   end do

   cplex_dij=Dtset%nspinor; ndij=1

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

     paw_gencond(1)=Dtset%pawlcutd ; paw_gencond(2)=Dtset%pawlmix
     paw_gencond(3)=Dtset%pawnphi  ; paw_gencond(4)=Dtset%pawntheta
     paw_gencond(5)=Dtset%pawspnorb; paw_gencond(6)=Dtset%pawxcdev
     call timab(553,2,tsec)
   else
     if (Pawtab(1)%has_kij  ==1) Pawtab(1:Cryst%ntypat)%has_kij  =2
     if (Pawtab(1)%has_nabla==1) Pawtab(1:Cryst%ntypat)%has_nabla=2
   end if
   Psps%n1xccc=MAXVAL(Pawtab(1:Cryst%ntypat)%usetcore)

!  Initialize optional flags in Pawtab to zero
!  (Cannot be done in Pawinit since the routine is called only if some parameters are changed)
   Pawtab(:)%has_nabla = 0
   Pawtab(:)%usepawu   = 0
   Pawtab(:)%useexexch = 0
   Pawtab(:)%exchmix   =zero

!  * Evaluate <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for the long wavelength limit.
!  TODO solve problem with memory leak and clean this part as well as the associated flag
   call pawnabla_init(Psps%mpsang,Psps%lmnmax,Cryst%ntypat,Psps%indlmn,Pawrad,Pawtab)

!  if (psp_gencond==1) then !.or. nsym_old/=Cryst%nsym) then
   call setsymrhoij(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,Cryst%rprimd,Cryst%symrec,Pawang%zarot)
   nsym_old=Cryst%nsym
!  end if

!  === Initialize and compute data for LDA+U ===
   if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
     Paw_dmft%use_dmft=Dtset%usedmft
     call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,&
&     Psps%indlmn,Psps%lmnmax,Cryst%ntypat,Pawang,Dtset%pawprtvol,Pawrad,Pawtab,Dtset%upawu,&
&     Dtset%usedmft,Dtset%useexexch,Dtset%usepawu)
     MSG_ERROR("BS equation with LDA+U not completely coded")
   end if
   ABI_CHECK(Dtset%usedmft==0,"DMFT + BSE not allowed")
   ABI_CHECK(Dtset%useexexch==0,"LEXX + BSE not allowed")

   if (my_rank==master) call print_pawtab(Pawtab)

!  === Get Pawrhoij from the header of the KSS file ===
   call rhoij_copy(Hdr_kss%pawrhoij,KS_Pawrhoij)

!  === Re-symmetrize symrhoij ===
!  this call leads to a SIGFAULT, likely some pointer is not initialized correctly
   choice=1; optrhoij=1
!  call symrhoij(choice,Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert0,Psps%lmnmax,Cryst%natom,Cryst%natom,Cryst%nsym,&
!  &  Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,KS_Pawrhoij,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)

!  === Evaluate form factor of radial part of phi.phj-tphi.tphj ===
   rhoxsp_method=1 ! Arnaud-Alouani
!  if (Dtset%userie==1) rhoxsp_method=1 ! Arnaud-Alouani
!  if (Dtset%userie==2) rhoxsp_method=2 ! Shiskin-Kresse

   ABI_ALLOCATE(gfft_osc,(3,nfftot_osc))
   call get_gftt(ngfft_osc,k0,gmet,gsq_osc,gfft_osc)
   ABI_DEALLOCATE(gfft_osc)

!  * Set up q grids, make qmax 20% larger than largest expected:
   ABI_ALLOCATE(nq_spl,(Psps%ntypat))
   ABI_ALLOCATE(qmax,(Psps%ntypat))
   nq_spl = Psps%mqgrid_ff
   qmax = SQRT(gsq_osc)*1.2d0 ! qmax=Psps%qgrid_ff(Psps%mqgrid_ff)
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
   call pawfgrtab_init(Pawfgrtab,cplex1,l_size_atm,Dtset%nspden)
   ABI_DEALLOCATE(l_size_atm)
   compch_fft=greatest_real
   usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
!  * 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
!  * 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   write(msg,'(a,i2)')' bethe_salpeter : using usexcnhat = ',usexcnhat
   call wrtout(std_out,msg,'COLL')
!  
!  === Identify parts of the rectangular grid where the density has to be calculated ===
   optcut=0; optgr0=Dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-Dtset%pawstgylm
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(Cryst%atindx1,gmet,MPI_enreg_seq,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
   optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%ucvol,Cryst%xred)
 end if !End of PAW Initialization

!Allocate these arrays anyway, since they are passed to subroutines.
 if (.not.allocated(ks_nhat))  then
   ABI_ALLOCATE(ks_nhat,(nfftf,0))
 end if

!==================================================
!==== Read KS band structure from the KSS file ====
!==================================================

!* Initialize wave function handler, allocate wavefunctions.
 my_minb=1; my_maxb=BSp%nbnds; my_spin=0 ! nsppol==2 not coded anyway.
 mband=BSp%nbnds
 ABI_ALLOCATE(nband,(Kmesh%nibz,Dtset%nsppol))
 nband=mband

!At present, no memory distribution, each node has the full set of states.
 ABI_ALLOCATE(bks_mask,(mband,Kmesh%nibz,Dtset%nsppol))
 bks_mask=.TRUE.

 ABI_ALLOCATE(keep_ur,(mband,Kmesh%nibz,Dtset%nsppol))
 keep_ur=.FALSE.
 if (MODULO(Dtset%gwmem,10)==1) then
   do spin=1,Dtset%nsppol
     if (spin==my_spin.or.my_spin==0) keep_ur(:,:,spin)=.TRUE.
   end do
 end if

 opt_ecut=zero !; if (gw_uses_wfk_file) opt_ecut=Dtset%ecutwfn

 call wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,BSp%npwwfn,mband,nband,Kmesh%nibz,Dtset%nsppol,bks_mask,&
& Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_kss%istwfk,Kmesh%ibz,ngfft_osc,&
& Gsph_Max%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)

 ABI_DEALLOCATE(bks_mask)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(keep_ur)

 call wfd_print(Wfd,mode_paral='PERS')

 ABI_TIMER_STOP("init1")
 call timab(651,2,tsec) ! bse(Init1)

 call timab(653,1,tsec) ! bse(rdkss)

 my_minb=1; my_maxb=BSp%nbnds

 if (gw_uses_wfk_file) then
!  if (.FALSE..and.gw_uses_wfk_file) then
   MSG_ERROR("Wfd is init using WFK file")
   wfk_fname = Dtfil%fnameabi_kss; ii=LEN_TRIM(wfk_fname)
   wfk_fname(ii-2:ii) = "WFK"
!  #ifdef HAVE_MPI_IO
!  call wfd_read_wfk(Wfd,wfk_fname,IO_MODE_MPI)
!  #else
   call wfd_read_wfk(Wfd,wfk_fname,Dtset%accesswff)
!  #endif
 else
   call wfd_read_kss(Wfd,Dtfil%fnameabi_kss,Bsp%nbnds,Dtset%accesswff,nelect_kss)
 end if

 call wfd_test_ortho(Wfd,Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

 call timab(653,2,tsec) ! bse(rdkss)

 call timab(655,1,tsec) ! bse(mkrho)

!=== Calculate the FFT index of $(R^{-1}(r-\tau))$ ===
!* S=\transpose R^{-1} and k_BZ = S k_IBZ
!* irottb is the FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.
 ABI_ALLOCATE(irottb,(nfftot_osc,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,ngfft_osc,irottb,iscompatibleFFT)

 ABI_ALLOCATE(ktabr,(nfftot_osc,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym=Kmesh%tabo(ik_bz)
   do ifft=1,nfftot_osc
     ktabr(ifft,ik_bz)=irottb(ifft,isym)
   end do
 end do
 ABI_DEALLOCATE(irottb)
!
!===========================
!=== COMPUTE THE DENSITY ===
!===========================
!* Evaluate Planewave part (complete charge in case of NC pseudos).
!
 ABI_ALLOCATE(ks_rhor,(nfftf,Wfd%nspden))
 call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_rhor)

!TODO this has to be done in a better way, moreover wont work for PAW
!Check Vcp!
!! call cutoff_density(ngfftf,Dtset%nspden,Dtset%nsppol,Vcp,ks_rhor,MPI_enreg_seq)
!
!=== Additional computation for PAW ===
 nhatgrdim=0
 if (Dtset%usepaw==1) then
!  
!  * Calculate the compensation charge nhat.
   if (Dtset%xclevel==2) nhatgrdim=usexcnhat*Dtset%pawnhatxc
   ider=2*nhatgrdim; izero=0; qphon(:)=zero
   if (nhatgrdim>0)  then
     ABI_ALLOCATE(ks_nhatgr,(nfftf,Dtset%nspden,3))
   end if

   call pawmknhat(compch_fft,cplex1,ider,idir0,ipert0,izero,Cryst%gprimd,MPI_enreg_seq,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Dtset%paral_kgb,Pawang,&
&   Pawfgrtab,ks_nhatgr,ks_nhat,KS_Pawrhoij,KS_Pawrhoij,Pawtab,qphon,Cryst%rprimd,Cryst%ucvol,Cryst%xred)
!  
!  === Evaluate onsite energies, potentials, densities ===
!  * Initialize variables/arrays related to the PAW spheres.
!  * Initialize also lmselect (index of non-zero LM-moments of densities).
!  TODO call init_paw_ij in scfcv and respfn, fix small issues
   ABI_ALLOCATE(KS_paw_ij,(Cryst%natom))
   call nullify_paw_ij(KS_paw_ij)

   cplex_dij=Dtset%nspinor
   has_dijso=Dtset%pawspnorb
   has_dijU=Dtset%usepawu
   has_dijso=Dtset%pawspnorb
   has_dijU=Dtset%usepawu

   call init_paw_ij(KS_paw_ij,cplex1,cplex_dij,Dtset%nspinor,Dtset%nsppol,&
&   Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&   has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=0,has_dijxc_val=0,&
&   has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1)

   ABI_ALLOCATE(KS_paw_an,(Cryst%natom))
   call nullify_paw_an(KS_paw_an)

   nkxc1=0
   call init_paw_an(Cryst%natom,Cryst%ntypat,nkxc1,Dtset%nspden,cplex1,Dtset%pawxcdev,&
&   Cryst%typat,Pawang,Pawtab,KS_paw_an,has_vxc=1,has_vxcval=0)
!  
!  === Calculate onsite vxc with and without core charge ===
   nzlmopt=-1; option=0; compch_sph=greatest_real

   call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,ipert0,&
&   Dtset%ixc,MPI_enreg_seq,Cryst%natom,Cryst%natom,Dtset%nspden,&
&   Cryst%ntypat,nzlmopt,option,Dtset%paral_kgb,KS_Paw_an,KS_Paw_an,KS_paw_ij,&
&   Pawang,Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,&
&   Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,Dtset%xclevel,Dtset%xc_denpos,Psps%znuclpsp)
 end if !PAW

 if (.not.allocated(ks_nhatgr))  then
   ABI_ALLOCATE(ks_nhatgr,(nfftf,Dtset%nspden,0))
 end if

 call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,ks_rhor,Cryst%ucvol,&
& Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)
!
!=== For PAW, add the compensation charge on the FFT mesh, then get rho(G) ===
 if (Dtset%usepaw==1) ks_rhor(:,:)=ks_rhor(:,:)+ks_nhat(:,:)
 call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_rhor,ucvol=ucvol)

 ABI_ALLOCATE(ks_rhog,(2,nfftf))

 call fourdp(1,ks_rhog,ks_rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp0)

 call timab(655,2,tsec) ! bse(mkrho)
!
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
 xred_dummy(:,:)=Cryst%xred
 nkxc=0
!if (Wfd%nspden==1) nkxc=2
!if (Wfd%nspden>=2) nkxc=3 ! check GGA and spinor, quite a messy part!!!
 ABI_ALLOCATE(kxc,(nfftf,nkxc))

 n3xccc=0; if (Psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))
 ABI_ALLOCATE(ks_vhartr,(nfftf))
 ABI_ALLOCATE(ks_vtrial,(nfftf,Wfd%nspden))
 ABI_ALLOCATE(vpsp,(nfftf))
 ABI_ALLOCATE(ks_vxc,(nfftf,Wfd%nspden))

 optene=4; moved_atm_inside=0; moved_rhor=0; initialized=1; istep=1
!
!=== Compute structure factor phases and large sphere cut-off ===
!WARNING cannot use Dtset%mgfft, this has to be checked better
!mgfft=MAXVAL(ngfftc(:))
!allocate(ph1d(2,3*(2*mgfft+1)*Cryst%natom),ph1df(2,3*(2*mgfftf+1)*Cryst%natom))
 write(std_out,*)' CHECK ',Dtset%mgfftdg,mgfftf
 if (Dtset%mgfftdg/=mgfftf) then
!  write(std_out,*)"WARNING Dtset%mgfftf /= mgfftf"
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

 ABI_DEALLOCATE(ph1d)

 call setvtr(Cryst%atindx1,Dtset,KS_energies,Cryst%gmet,Cryst%gprimd,grewtn,gsqcutf_eff,&
& istep,kxc,mgfftf,moved_atm_inside,moved_rhor,MPI_enreg_seq,&
& Cryst%nattyp,nfftf,ngfftf,ks_nhat,ks_nhatgr,nhatgrdim,nkxc,Cryst%ntypat,Psps%n1xccc,n3xccc,&
& optene,Pawtab,ph1df,Psps,ks_rhog,ks_rhor,Cryst%rmet,Cryst%rprimd,strsxc,&
& Cryst%ucvol,usexcnhat,ks_vhartr,vpsp,ks_vtrial,ks_vxc,vxcavg,wvl,xccc3d,xred_dummy)
!TODO here xred is INOUT due to ionion_realSpace and xredcart!

 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(vpsp)

!============================
!==== Compute KS PAW Dij ====
!============================
 if (Wfd%usepaw==1) then
   write(std_out,*)"Another silly write for XLF"
!  
!  Calculate the unsymmetrized Dij. 
   call pawdij(cplex1,Dtset,Dtset%enunit,one,Cryst%gprimd,ipert0,MPI_enreg_seq,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,Dtset%nspden,Cryst%ntypat,&
&   Dtset%paral_kgb,KS_paw_an,KS_paw_ij,Pawang,Pawfgrtab,&
&   Dtset%pawprtvol,Pawrad,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
&   k0,Cryst%typat,Cryst%ucvol,ks_vtrial,ks_vxc,Cryst%xred)
!  
!  Symmetrize KS Dij 
   call symdij(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert0,&
&   Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,0,KS_paw_ij,Pawang,&
&   Dtset%pawprtvol,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
!  
!  Output the pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   call pawprt(Dtset,Psps%indlmn,Psps%lmnmax,KS_paw_ij,KS_Pawrhoij,Pawtab)
 end if

 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(xccc3d)
 ABI_DEALLOCATE(grewtn)
 ABI_DEALLOCATE(xred_dummy)

!=== QP_BSt stores energies and occ. used for the calculation ===
!* Initialize QP_BSt with KS values.
!* In case of SC update QP_BSt using the QPS file.
 ABI_ALLOCATE(qp_rhor,(nfftf,Dtset%nspden))
 qp_rhor=ks_rhor

 ABI_ALLOCATE(qp_aerhor,(nfftf,Dtset%nspden))
 qp_aerhor=ks_rhor
 if (Wfd%usepaw==1 .and. BSp%mdlf_type /=0) then ! Have to calculate the AE rhor
   MSG_ERROR("qp_aerhor with PAW not coded")
 end if

#if 0
 ABI_ALLOCATE(qp_rhor_paw   ,(nfftf,Wfd%nspden))
 ABI_ALLOCATE(qp_rhor_n_one ,(nfftf,Wfd%nspden))
 ABI_ALLOCATE(qp_rhor_nt_one,(nfftf,Wfd%nspden))

 ABI_ALLOCATE(qp_nhat,(nfftf,Wfd%nspden))
 qp_nhat = ks_nhat

 call denfgr(Cryst%atindx1,Cryst%gmet,Wfd%comm,Cryst%natom,Cryst%nattyp,ngfftf,qp_nhat,&
& Wfd%nspinor,Wfd%nsppol,Wfd%nspden,Cryst%ntypat,Pawfgr,Pawrad,KS_pawrhoij,Pawtab,Dtset%prtvol,&
& Psps,qp_rhor,qp_rhor_paw,qp_rhor_n_one,qp_rhor_nt_one,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

 ABI_DEALLOCATE(qp_nhat)
 ABI_DEALLOCATE(qp_rhor_n_one)
 ABI_DEALLOCATE(qp_rhor_nt_one)
 norm = SUM(qp_rhor_paw(:,1))*Cryst%ucvol/PRODUCT(Pawfgr%ngfft(1:3))
 write(msg,'(a,F8.4)') '  QUASIPARTICLE DENSITY CALCULATED - NORM OF DENSITY: ',norm
 call wrtout(std_out,msg,'PERS')

 ABI_DEALLOCATE(qp_rhor_n_one)
 ABI_DEALLOCATE(qp_rhor_nt_one)
 ABI_DEALLOCATE(qp_rhor_paw)
#endif

!call copy_bandstructure(KS_BSt,QP_BSt)

 if (.FALSE.) then
!  $ m_lda_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}> $
   ABI_ALLOCATE(m_lda_to_qp,(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol))
   m_lda_to_qp=czero
   do spin=1,Wfd%nsppol
     do ik_ibz=1,Wfd%nkibz
       do band=1,Wfd%nband(ik_ibz,spin)
         m_lda_to_qp(band,band,ik_ibz,spin)=cone ! Initialize the QP amplitudes with KS wavefunctions.
       end do
     end do
   end do
!  
!  * Now read m_lda_to_qp and update the energies in QP_BSt.
!  TODO switch on the renormalization of n in sigma.
   ABI_ALLOCATE(prev_rhor,(nfftf,Wfd%nspden))
   ABI_ALLOCATE(prev_Pawrhoij,(Cryst%natom*Wfd%usepaw))

   call rdqps(QP_BSt,Dtfil%fnameabi_qps,Wfd%usepaw,Wfd%nspden,1,nscf,&
&   nfftf,ngfftf,Cryst%ucvol,Wfd%paral_kgb,Cryst,Pawtab,MPI_enreg_seq,nbsc,m_lda_to_qp,prev_rhor,prev_Pawrhoij)

   ABI_DEALLOCATE(prev_rhor)
   if (Psps%usepaw==1.and.nscf>0) then
     call rhoij_free(prev_pawrhoij)
   end if
   ABI_DEALLOCATE(prev_pawrhoij)
   istat = ABI_ALLOC_STAT
!  
!  if (nscf>0.and.wfd_iam_master(Wfd)) then ! Print the unitary transformation on std_out.
!  call show_QP(QP_BSt,m_lda_to_qp,fromb=Sigp%minbdgw,tob=Sigp%maxbdgw,unit=std_out,tolmat=0.001_dp)
!  end if
!  
!  === Compute QP wfg as linear combination of KS states ===
!  * Wfd%ug is modified inside calc_wf_qp
!  * For PAW, update also the on-site projections.
!  * WARNING the first dimension of MPI_enreg MUST be Kmesh%nibz
!  TODO here we should use nbsc instead of nbnds

   call wfd_rotate(Wfd,Cryst,m_lda_to_qp)

   ABI_DEALLOCATE(m_lda_to_qp)
!  
!  === Reinit the storage mode of Wfd as ug have been changed ===
!  * Update also the wavefunctions for GW corrections on each processor
   call wfd_reset_ur_cprj(Wfd)

   call wfd_test_ortho(Wfd,Cryst,Pawtab,unit=ab_out,mode_paral="COLL")
!  
!  Compute QP occupation numbers.
   write(msg,'(3a)')ch10,' bethe_salpeter : calculating QP occupation numbers ',ch10
   call wrtout(std_out,msg,'COLL')

   call update_occ(QP_BSt,Dtset%fixmom,prtvol=0)
   ABI_ALLOCATE(qp_vbik,(QP_BSt%nkpt,QP_BSt%nsppol))
   qp_vbik(:,:) = get_valence_idx(QP_BSt)
   ABI_DEALLOCATE(qp_vbik)

   call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_rhor)
 end if

 ABI_ALLOCATE(qp_rhog,(2,nfftf))
 call fourdp(1,qp_rhog,qp_rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Wfd%paral_kgb,0)

 if (.FALSE.) then
!  if (.TRUE.) then
   call exc_interp_ham(BSp,BS_files,Dtset,Cryst,Kmesh,Qmesh,KS_BSt,ktabr,Gsph_Max,Gsph_c,Wfd,Hdr_bse,&
&   nfftot_osc,ngfft_osc,Psps,Pawtab,KS_pawrhoij,KS_Paw_ij,Pawang,Pawrad,Pawfgr,Paw_pwff,ngfftc,ngfftf,nfftf,qp_aerhor,ks_vtrial)
 end if

!
!States up to lomo-1 are useless now since only the states close to the gap are 
!needed to construct the EXC Hamiltonian. Here we deallocate the wavefunctions 
!to make room for the excitonic Hamiltonian that is allocated in exc_build_ham.
!and for the screening that is allocated below.
!Hereafter bands from 1 up to lomo-1 should not be accessed!
 ABI_ALLOCATE(bks_mask,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))
 bks_mask=.FALSE.
 if (Bsp%lomo>1) bks_mask(1:Bsp%lomo-1,:,:)=.TRUE.
 call wfd_wave_free(Wfd,"All",bks_mask)
 ABI_DEALLOCATE(bks_mask)
!
!================================================================
!Build the screened interaction W in the irreducible wedge.
!* W(q,G1,G2) = vc^{1/2} (q,G1) e^{-1}(q,G1,G2) vc^{1/2) (q,G2)
!* Use Coulomb term for q-->0,
!* Only the first small Q is used, shall we average if nqlwl>1?
!================================================================
!TODO clean this part and add an option to retrieve a single frequency to save memory.
 ABI_TIMER_START("readmake_W")
 call timab(654,1,tsec) ! bse(rdmkeps^-1)

 call screen_nullify(W)
 if (BSp%use_coulomb_term) then !  Init W.
   W_info%mat_type = MAT_INV_EPSILON
!  W_info%mat_type = MAT_W
!  W_info%vtx_family
!  W_info%ixc
!  W_info%use_ada
!  W_info%use_mdf = MDL_BECHSTEDT
!  W_info%use_ppm
!  W_info%vtx_test
!  W_info%wint_method
!  
!  W_info%ada_kappa
!  W_info%eps_inf = 23
!  W_info%drude_plsmf

   mqmem = Qmesh%nibz !TODO out-of-core solution not implemented yet.
   call screen_init(W,W_Info,Dtset,Cryst,Qmesh,Gsph_c,Vcp,w_fname,mqmem,Dtset%npweps,&
&   Dtset%accesswff,ngfftf,nfftf_tot,Wfd%nsppol,Wfd%nspden,qp_aerhor,Wfd%prtvol,Wfd%comm)
 end if
 call timab(654,2,tsec) ! bse(rdmkeps^-1)
 ABI_TIMER_STOP("readmake_W")
!
!=============================================
!==== Build of the excitonic Hamiltonian =====
!=============================================
!
 ABI_TIMER_START("mkexcham")
 call timab(656,1,tsec) ! bse(mkexcham)

 call exc_build_ham(BSp,BS_files,Cryst,Kmesh,Qmesh,ktabr,Gsph_Max,Gsph_c,Vcp,&
& Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff)
!
!Free Em1 to make room for the full excitonic Hamiltonian.
 call screen_free(W)

 ABI_TIMER_STOP("mkexcham")
 call timab(656,2,tsec) ! bse(mkexcham)
!
!=========================================
!==== Macroscopic dielectric function ====
!=========================================
 ABI_TIMER_START("mkexceps")
 call timab(657,1,tsec) ! bse(mkexceps)
!
!First deallocate the internal %ur buffers to make room for the excitonic Hamiltonian.
 call timab(658,1,tsec) ! bse(wfd_wave_free)
 ABI_ALLOCATE(bks_mask,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))
 bks_mask=.TRUE.
 call wfd_wave_free(Wfd,"Real_space",bks_mask)
 ABI_DEALLOCATE(bks_mask)
 call timab(658,2,tsec) ! bse(wfd_wave_free)
!
!Compute the commutator [r,Vu] (PAW only).
 call timab(659,1,tsec) ! bse(make_Hur_commutator)
 ABI_ALLOCATE(HUr,(Cryst%natom*Wfd%usepaw))
 call nullify_Hur(HUr)
 if (Dtset%usepawu/=0) then !TODO here I need KS_Paw_ij
   MSG_ERROR("Commutator for LDA+U not tested")
   call make_Hur_commutator(Wfd%nsppol,Wfd%pawprtvol,Cryst,Psps,Pawtab,Pawang,Pawrad,KS_Paw_ij,Hur)
 end if
 call timab(659,2,tsec) ! bse(make_Hur_commutator)

 select case (BSp%algorithm)
   case (BSE_ALGO_DDIAGO, BSE_ALGO_CG)
     call timab(660,1,tsec) ! bse(exc_diago_driver)
     call exc_diago_driver(Wfd,Bsp,BS_files,KS_BSt,QP_BSt,Cryst,Kmesh,Psps,&
&     Pawtab,Hur,Hdr_bse,drude_plsmf)
     call timab(660,2,tsec) ! bse(exc_diago_driver)

     if (.FALSE.) then ! Calculate electron-hole excited state density. Not tested at all.
       call exc_den(BSp,BS_files,ngfftf,nfftf,Kmesh,ktabr,Wfd)
     end if

     if (.FALSE.) then
       paw_add_onsite=.FALSE.; spin_opt=1; which_fixed=1; eh_rcoord=(/zero,zero,zero/); nrcell=(/2,2,2/)
!      call exc_plot(Bsp,Bs_files,Wfd,Kmesh,Cryst,Psps,Pawtab,Pawrad,paw_add_onsite,spin_opt,which_fixed,eh_rcoord,nrcell,ngfftf)
     end if
!    
   case (BSE_ALGO_Haydock)
     call timab(661,1,tsec) ! bse(exc_haydock_driver)
     call exc_haydock_driver(BSp,BS_files,Cryst,Kmesh,Hdr_bse,KS_BSt,QP_BSt,Wfd,Psps,Pawtab,Hur)
     call timab(661,2,tsec) ! bse(exc_haydock_driver)
!    
     case default
     write(msg,'(a,i0)')" Wrong BSE algorithm: ",BSp%algorithm
     MSG_ERROR(msg)
 end select

 call timab(657,2,tsec) ! bse(mkexceps)

 ABI_TIMER_STOP("mkexceps")
!
!=====================
!==== Free memory ====
!=====================
 ABI_DEALLOCATE(ktabr)
 ABI_DEALLOCATE(Pawfgr%fintocoa)
 ABI_DEALLOCATE(Pawfgr%coatofin)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ks_vhartr)
 ABI_DEALLOCATE(ks_vtrial)
 ABI_DEALLOCATE(ks_vxc)
 ABI_DEALLOCATE(ks_nhat)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ks_nhatgr)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ks_rhog)
 ABI_DEALLOCATE(ks_rhor)
 ABI_DEALLOCATE(qp_rhog)
 ABI_DEALLOCATE(qp_rhor)
 ABI_DEALLOCATE(qp_aerhor)
!
!* Destroy local data structures.
 call destroy_crystal(Cryst)
 call destroy_gsphere(Gsph_Max)
 call destroy_gsphere(Gsph_c)
 call destroy_bz_mesh_type(Kmesh)
 call destroy_bz_mesh_type(Qmesh)
 call hdr_clean(Hdr_kss)
 call hdr_clean(Hdr_bse)
 call bstruct_clean(KS_BSt)
 call bstruct_clean(QP_BSt)
 call vcoul_free(Vcp)
 call destroy_Hur(Hur)
 ABI_DEALLOCATE(Hur)
 istat = ABI_ALLOC_STAT
 call destroy_bs_parameters(BSp)
 call wfd_destroy(Wfd)

 if (Dtset%usepaw==1) then ! Optional deallocation for PAW.
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
 end if

 call timab(650,2,tsec) ! bse(Total)
 ABI_TIMER_STOP("")

 DBG_EXIT('COLL')

end subroutine bethe_salpeter
!!***
