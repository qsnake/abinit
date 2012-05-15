!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_sigc_me
!! NAME
!! calc_sigc_me
!!
!! FUNCTION
!! Calculate diagonal and off-diagonal matrix elements of the self-energy operator.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (FB, GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! minbnd, maxbnd= min and Max band index for GW correction (for this k-point)
!! Dtset <type(dataset_type)>=all input variables in this dataset
!!    %accesswff
!!    %paral_kgb
!!    %nspinor=Number of spinorial components.
!!    %gwcomp=If 1 use an extrapolar approximation to accelerate convergence.
!!    %gwencomp=Extrapolar energy.
!! Er <Epsilonm1_results> (see the definition of this structured datatype)
!!    %mqmem=if 0 use out-of-core method in which a single q-slice of espilon is read inside the loop over k
!!    %nomega_i=Number of imaginary frequencies.
!!    %nomega_r=Number of real frequencies.
!!    %nomega=Total number of frequencies.
!! Gsph_c<Gvectors_type>= info on G-sphere for Sigma_c
!! Gsph_Max<Gvectors_type>= info on biggest G-sphere
!!    %nsym=number of symmetry operations
!!    %rottb(Sigp%npwvec,timrev,nsym)=index of (IS) G where I is the identity or the inversion
!!      operation and G is one of the npwvec vectors in reciprocal space
!!    %timrev=2 if time-reversal symmetry is used, 1 otherwise
!!    %gvec(3,Sigp%npwvec)=integer coordinates of each plane wave in reciprocal space
!! ikcalc=index in the array Sigp%kptgw2bz of the k-point where GW corrections are calculated
!! Ltg_k datatype containing information on the little group
!! Kmesh <BZ_mesh_type>
!!    %nbz=Number of points in the BZ
!!    %nibz=Number of points in IBZ
!!    %kibz(3,nibz)=k-point coordinates, irreducible Brillouin zone
!!    %kbz(3,nbz)=k-point coordinates, full Brillouin zone
!!    %ktab(nbz)= table giving for each k-point in the BZ (kBZ), the corresponding
!!    %ktabi(nbz)= for each k-point in the BZ defines whether inversion has to be considered
!!    %ktabp(nbz)= phase factor associated to tnons
!! gwc_ngfft(18)=Information about 3D FFT for the oscillator strengths used for the correlation part,
!! Vcp <vcoul_t datatype> containing information on the cutoff technique
!!    %vc_sqrt(npwc,nqibz)= square-root of the coulombian potential for q-points in the IBZ
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang <type(pawang_type)>=paw angular mesh and related data
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    %usepaw=1 for PAW, 0 for NC pseudopotentials.
!! Qmesh <bz_mesh_type> : datatype gathering information of the q-mesh used
!!    %ibz=q points where $\tilde\epsilon^{-1}$ has been computed
!!    %bz(3,nqbz)=coordinates of all q-points in BZ
!! Sigp <sigma_parameters> (see the definition of this structured datatype)
!!    %npwvec= Max betwee npweps and npwwfn used to dimension arrays
!! Cryst<Crystal_structure>=Info on unit cell and symmetries
!!    %natom=number of atoms in unit cell
!!    %ucvol=unit cell volume
!!    %nsym=number of symmetry operations
!!    %typat(natom)=type of each atom
!!  much slower but it requires less memory
!! PPm<PPmodel_type>= Datatype gathering information on the Plasmonpole technique (see also ppm_get_qbz).
!!    %model=type of ppmodel
!!    %npwc=number of G-vectors for the correlation part.
!!    %dm2_otq =size of second dimension of otq array
!!    %dm2_bots=size of second dimension of botsq arrays
!!    %dm_eig  =size of second dimension of eig arrays
!! QP_BSt<Bandstructure_type>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!  Paw_pwff<Paw_pwff_type>=Form factor used to calculate the onsite mat. elements of a plane wave.
!! QP_sym(Wfd%nsppol)<bands_symmetries>=Datatype collecting data on the irreducible representaions of the
!!    little group of kcalc in the KS representation as well as the symmetry of the bdgw_k states.
!!  Sr=sigma_results (see the definition of this structured datatype)
!!  use_aerhor=1 is aepaw_rhor is used, 0 otherwise.
!!  aepaw_rhor(rho_nfftot,Wfd%nspden*use_aerhor)=AE PAW density used to generate PPmodel paramenters if mqmem==0
!!
!! OUTPUT
!!
!! NOTES
!!  1) The treatment of the divergence of Gygi+Baldereschi (PRB 1986) is included.
!!  2) The calculation of energy derivative is based on finite elements.
!!  3) On the symmetrization of Sigma matrix elements ***/
!!        If  Sk = k+G0 then  M_G(k, Sq)= e^{-i( Sq+G).t} M_{ S^-1(G}   (k,q)
!!        If -Sk = k+G0 then  M_G(k,-Sq)= e^{-i(-Sq+G).t} M_{-S^-1(G)}^*(k,q)
!!
!!     Notice the absence of G0 in the expression. Moreover, when we sum over the little group, it turns out
!!     that there is a cancellation of the phase factor associated to the non-symmorphic operations due to a
!!     similar term coming from the symmetrization of \epsilon^{-1}. Mind however that the nonsymmorphic phase
!!     has to be considered when epsilon^-1 is reconstructed starting from the q-points in the IBZ.
!!
!!  4) The unitary transformation relating wavefunctions
!!     at symmetric k-points should be taken into account during the symmetrization
!!     of the oscillator matrix elements. In case of G_oW_o and GW_o calculations, however,
!!     it is possible to make an invariant by just including all the degenerate states and
!!     averaging the final results over the degenerate subset. 
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      calc_coh_comp,calc_sig_ppm,calc_sig_ppm_comp,calc_sigc_cd,calc_wfwfg
!!      coeffs_gausslegint,cprj_alloc,cprj_copy,cprj_free,destroy_paw_pwij
!!      epsm1_pole_symmetrizer,epsm1_pole_symmetrizer_inplace,epsm1_symmetrizer
!!      epsm1_symmetrizer_inplace,findqg0,flush_unit,get_bz_item,get_epsm1
!!      get_gftt,get_pole_epsm1,gsph_fft_tabs,gw_eet_sigma,init_paw_pwij
!!      initmpi_seq,paw_cross_rho_tw_g,paw_rho_tw_g,paw_symcprj,ppm_get_qbz
!!      print_little_group,rho_tw_g,rotate_fft_mesh,setup_ppmodel
!!      sigma_distribution,symmetrize_me,timab,wfd_barrier,wfd_change_ngfft
!!      wfd_get_cprj,wfd_get_ur,wfd_paw_get_aeur,wrtout,xmax_mpi,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine calc_sigc_me(ikcalc,nomega_sigc,minbnd,maxbnd,Dtset,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_Max,Gsph_c,Vcp,Kmesh,Qmesh,&
& Ltg_k,PPm,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,&
& gwc_ngfft,rho_ngfft,rho_nfftot,rhor,use_aerhor,aepaw_rhor,sigcme_tmp)

 use m_profiling

 use defs_basis
 use m_gwdefs !,        only : czero_gw, cone_gw, j_gw, sigma_parameters, sigma_type_from_key, sigma_is_herm
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_defs_ptgroups
#ifdef HAVE_CLIB
 use m_clib
#endif
 use m_errors

 use m_blas,          only : xdotc, xgemv
 use m_numeric_tools, only : hermitianize, imin_loc
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_structure
 use m_bz_mesh,       only : bz_mesh_type, get_BZ_item, findqg0, little_group, print_little_group
 use m_gsphere,       only : gvectors_type, gsph_fft_tabs
 use m_fft_mesh,      only : get_gftt, rotate_fft_mesh, cigfft
 use m_vcoul,         only : vcoul_t
 use m_paw_pwij,      only : paw_pwff_type, paw_pwij_type, init_paw_pwij, destroy_paw_pwij, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_wfs,           only : wfd_get_ur, wfs_descriptor, wfd_get_cprj, wfd_barrier, wfd_change_ngfft, wfd_paw_get_aeur
 use m_paw_toolbox,   only : paw_pwaves_lmn_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_screening,     only : epsilonm1_results, epsm1_symmetrizer, epsm1_symmetrizer_inplace, get_epsm1, &
 &                           epsm1_pole_symmetrizer, epsm1_pole_symmetrizer_inplace, get_pole_epsm1
 use m_ppmodel,       only : setup_ppmodel, ppm_get_qbz, ppmodel_type, calc_sig_ppm
 use m_sigma_results, only : sigma_results 
 use m_bands_sym,     only : bands_symmetries, symmetrize_me, bsym_failed
 use m_io_tools,      only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sigc_me'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_66_paw
 use interfaces_70_gw, except_this_one => calc_sigc_me
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,rho_nfftot,nomega_sigc,minbnd,maxbnd
 integer,intent(in) :: use_aerhor
 type(Crystal_structure),intent(in) :: Cryst
 type(Bandstructure_type),intent(in) :: QP_BSt
 type(BZ_mesh_type),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_results),intent(inout) :: Er
 type(Gvectors_type),intent(in) :: Gsph_Max
 type(Gvectors_type),intent(in) :: Gsph_c
 type(Little_group),intent(in) :: Ltg_k
 type(PPmodel_type),intent(inout) :: PPm
 type(Pseudopotential_type),intent(in) :: Psps
 type(pawang_type),intent(in) :: pawang
 type(Sigma_parameters),intent(in) :: Sigp
 type(Sigma_results),intent(in) :: Sr
 type(wfs_descriptor),intent(inout) :: Wfd,Wfdf
!arrays
 integer,intent(in) :: gwc_ngfft(18),rho_ngfft(18)
 real(dp),intent(in) :: rhor(rho_nfftot,Wfd%nspden)
 real(dp),intent(in) :: aepaw_rhor(rho_nfftot,Wfd%nspden*use_aerhor)
 !complex(dpc),intent(out) :: sigcme_tmp(:,:,:,:)
 complex(dpc),intent(out) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd,minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 !complex(dpc),pointer :: sigcme_tmp(:,:,:,:)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(Paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(bands_symmetries),intent(in) :: QP_sym(Wfd%nsppol)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp2=2
 integer :: npw_k !npw_kgw,
 integer :: iab,ib,ib1,ib2,ierr,ig,iggp,igp,ii,iik,itim_q,i1,i2,npls 
 integer :: ik_bz,ik_ibz,io,iiw,isym_q,iq_bz,iq_ibz,isppol,istat,isym,jb,is_idx,spin
 integer :: band,band1,band2,idle,rank
 integer :: jik,jk_bz,jk_ibz,kb,nspinor
 integer :: nomega_tot,nq_summed,ispinor,ibsp,dimcprj_gw
 integer :: spad,spadc,spadc1,spadc2,irow,my_nbks
 integer :: comm,ndegs,wtqm,wtqp,mod10
 integer :: isym_kgw,isym_ki,gwc_mgfft,use_padfft,gwc_fftalga,gwc_nfftot,nfftf,mgfftf,use_padfftf
 integer :: nbmax,nbopt,iwc,ifft
 real(dp) :: e0i,fact_sp,theta_mu_minus_e0i,tol_empty,z2,en_high,norm,gw_gsq
 real(dp) :: w_localmax,w_max
 complex(dpc) :: ctmp,omegame0i2_ac,omegame0i_ac,scprod,ph_mkgwt,ph_mkt !,drotp !,drotm
 logical :: iscompatibleFFT,q_is_gamma,pole_screening
 character(len=500) :: msg,sigma_type
 complex(gwpc),allocatable :: botsq(:,:),otq(:,:),eig(:,:)
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: g0(3),spinor_padc(2,4)
 integer,pointer :: kg_k(:,:) 
 integer,allocatable :: proc_distrb(:,:,:),extrapolar_distrb(:,:,:,:)
 integer,allocatable :: degtab(:,:,:)
 integer,allocatable :: igfftcg0(:),gw_gfft(:,:),gw_gbound(:,:),irottb(:,:),ktabr(:,:)
 integer,allocatable :: igfftfcg0(:),gboundf(:,:),ktabrf(:,:)
 integer :: got(Wfd%nproc),npoles_missing(minbnd:maxbnd)
 real(dp) :: ksum(3),kgw(3),kgw_m_ksum(3),omegap(Er%nomega_i),omegap2(Er%nomega_i),q0(3),tsec(2),qbz(3)
 real(dp) :: gl_knots(Er%nomega_i),gl_wts(Er%nomega_i)
 real(dp) :: spinrot_kbz(4),spinrot_kgw(4),w_maxval(minbnd:maxbnd)
 real(dp),pointer :: qp_ene(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: omegame0i(:)
 !real(dp),allocatable :: ks_rhor_paw(:,:)
 real(gwpc), allocatable :: epsm1_pole_qbz(:,:,:) !,epsm1_pole_trcc_qbz(:,:,:)
 complex(gwpc) :: sigcohme(Sigp%nsig_ab)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:),rhotwg(:),rhotwgp(:)
 complex(gwpc),allocatable :: botsq_conjg_transp(:,:),ac_epsm1cqwz2(:,:,:)
 complex(gwpc),allocatable :: epsm1_qbz(:,:,:),epsm1_trcc_qbz(:,:,:)
 !complex(gwpc),allocatable :: epsm1_herm_qbz(:,:,:)
 complex(gwpc),allocatable :: ac_integr(:,:,:),sigc_ket(:,:),ket1(:,:),ket2(:,:)
 complex(gwpc),allocatable :: herm_sigc_ket(:,:),aherm_sigc_ket(:,:)
 complex(gwpc),allocatable :: rhotwg_ki(:,:)
 complex(gwpc),allocatable :: sigcme2(:,:),sigcme_3(:),sigcme_new(:),sigctmp(:,:)
 complex(gwpc),allocatable :: wfr_bdgw(:,:),wfr_sum(:),wf1swf2_g(:)
 complex(gwpc),allocatable :: ur_ae_sum(:),ur_ae_onsite_sum(:),ur_ps_onsite_sum(:)
 complex(gwpc),allocatable :: ur_ae_bdgw(:,:),ur_ae_onsite_bdgw(:,:),ur_ps_onsite_bdgw(:,:)
 complex(gwpc),allocatable :: otq_transp(:,:)
 complex(gwpc),pointer :: cg_jb(:),cg_sum(:)
 complex(dpc) :: ovlp(2)
 complex(dpc),allocatable :: sym_cme(:,:,:),sigc(:,:,:,:,:)
 logical :: rank_mask(Wfd%nproc),can_symmetrize(Wfd%nsppol)
!logical :: me_calc_poles(Sr%nomega_r+Sr%nomega4sd)
 type(sigijtab_t),pointer :: Sigcij_tab(:)
 type(Cprj_type),allocatable :: Cprj_kgw(:,:),Cprj_ksum(:,:)
 type(Paw_pwij_type),allocatable :: Pwij_qg(:),Pwij_fft(:)

!************************************************************************

 DBG_ENTER("COLL")
 !
 ! === Initial check ===
 ABI_CHECK(Sr%nomega_r==Sigp%nomegasr,"")
 ABI_CHECK(Sr%nomega4sd==Sigp%nomegasrd,"")

 ABI_CHECK(Sigp%npwc==Gsph_c%ng,"")
 ABI_CHECK(Sigp%npwvec==Gsph_Max%ng,"")

 !ABI_CHECK( ALL(gwc_ngfft(1:3) == Wfd%ngfft(1:3))," cannot change FFT on-the-fly yet")

 call timab(431,1,tsec) ! calc_sigc_me
 call timab(432,1,tsec) ! Init
 !
 ! === Initialize MPI variables === 
 comm=Wfd%comm
                                                     
 ! * Fake MPI_type for sequential part.
 call initmpi_seq(MPI_enreg_seq) 

 qp_ene => QP_BSt%eig(:,:,:)
 qp_occ => QP_BSt%occ(:,:,:)
 !
 ! Index of the GW point in the BZ array, its image in IBZ and time-reversal ===
 jk_bz=Sigp%kptgw2bz(ikcalc)
 call get_BZ_item(Kmesh,jk_bz,kgw,jk_ibz,isym_kgw,jik,ph_mkgwt)
 !$call get_IBZ_item(Kmesh,jk_ibz,kibz,wtk)
 spinrot_kgw=Cryst%spinrot(:,isym_kgw)
 !
 ib1=minbnd
 ib2=maxbnd

 write(msg,'(2a,3f8.3,2a,2(i3,a))')ch10,&
&  ' Calculating <nk|Sigma_c(omega)|nk> at k = ',kgw(:),ch10,&
&  ' bands n = from ',ib1,' to ',ib2,ch10
 call wrtout(std_out,msg,'COLL')

 w_maxval = zero

 if ( ANY(gwc_ngfft(1:3) /= Wfd%ngfft(1:3)) ) call wfd_change_ngfft(Wfd,Cryst,Psps,gwc_ngfft) 
 gwc_mgfft   = MAXVAL(gwc_ngfft(1:3))
 gwc_fftalga = gwc_ngfft(7)/100 !; gwc_fftalgc=MOD(gwc_ngfft(7),10)

 if (Dtset%pawcross==1) then
   mgfftf = MAXVAL(rho_ngfft(1:3)) 
 end if

 can_symmetrize = .FALSE.
 if (Sigp%symsigma>0) then
   can_symmetrize = .TRUE.
   if (Sigp%gwcalctyp >= 20) then
    do spin=1,Wfd%nsppol
      can_symmetrize(spin) = .not.bsym_failed(QP_sym(spin))
      if (.not.can_symmetrize(spin)) then
        write(msg,'(a,i0,4a)')" Symmetrization cannot be performed for spin: ",spin,ch10,&
&         " band classification encountered the following problem: ",ch10,TRIM(QP_sym(spin)%err_msg)
        MSG_WARNING(msg)
      end if
    end do
   end if
   ABI_CHECK(Wfd%nspinor==1,'Symmetrization with nspinor=2 not implemented')
 end if

 ABI_UNUSED(Pawang%l_max)
 !
 ! Print type of calculation.
 mod10=MOD(Sigp%gwcalctyp,10); sigma_type = sigma_type_from_key(mod10)
 call wrtout(std_out,sigma_type,'COLL')
 !
 ! === Set up logical flags for Sigma calculation ===
 if (mod10==SIG_GW_AC.and.Sigp%gwcalctyp/=1) then
   MSG_ERROR("not implemented")
 end if

 ! === Initialize some values ===
 nspinor = Wfd%nspinor
 spinor_padc(:,:)=RESHAPE((/0,0,Sigp%npwc,Sigp%npwc,0,Sigp%npwc,Sigp%npwc,0/),(/2,4/))
 npoles_missing=0
 !
 ! === Normalization of theta_mu_minus_e0i ===
 ! * If nsppol==2, qp_occ $\in [0,1]$
 SELECT CASE (Wfd%nsppol)
 CASE (1)
   fact_sp=half; tol_empty=0.01   ! below this value the state is assumed empty
   if (Wfd%nspinor==2) then
    fact_sp=one; tol_empty=0.005  ! below this value the state is assumed empty
   end if
 CASE (2)
   fact_sp=one; tol_empty=0.005   ! to be consistent and obtain similar results if a metallic
 CASE DEFAULT                     ! spin unpolarized system is treated using nsppol==2
   MSG_BUG('Wrong nsppol')
 END SELECT
 !
 ! Allocate arrays used to accumulate the matrix elements of \Sigma_c over
 ! k-points and bands. Note that for AC requires only the imaginary frequencies
 !nomega_sigc=Sr%nomega_r+Sr%nomega4sd
 !if (mod10==SIG_GW_AC) nomega_sigc=Sr%nomega_i
 !
 ! === Define the G-G0 shifts for the FFT of the oscillators ===
 ! * Sigp%mG0 gives the MAX G0 component to account for umklapp.
 ! * Note the size MAX(Sigp%npwx,Sigp%npwc).
 !
 !
 ! === Precalculate the FFT index of $(R^{-1}(r-\tau))$ ===
 ! * S=\transpose R^{-1} and k_BZ = S k_IBZ
 ! * irottb is the FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.
 gwc_nfftot = PRODUCT(gwc_ngfft(1:3))
 ABI_ALLOCATE(irottb,(gwc_nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,gwc_ngfft,irottb,iscompatibleFFT)
 if (.not.iscompatibleFFT) then
   msg = "FFT mesh is not compatible with symmetries. Results might be affected by large errors!"
   MSG_WARNING(msg)
 end if

 ABI_ALLOCATE(ktabr,(gwc_nfftot,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym=Kmesh%tabo(ik_bz)
   do ifft=1,gwc_nfftot
     ktabr(ifft,ik_bz)=irottb(ifft,isym)
   end do
 end do
 ABI_DEALLOCATE(irottb)
 if (Psps%usepaw==1 .and. Dtset%pawcross==1) then
   nfftf = PRODUCT(rho_ngfft(1:3))
   ABI_ALLOCATE(irottb,(nfftf,Cryst%nsym))
   call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,rho_ngfft,irottb,iscompatibleFFT)

   ABI_ALLOCATE(ktabrf,(nfftf,Kmesh%nbz))
   do ik_bz=1,Kmesh%nbz
     isym=Kmesh%tabo(ik_bz)
     do ifft=1,nfftf
       ktabrf(ifft,ik_bz)=irottb(ifft,isym)
     end do
   end do
   ABI_DEALLOCATE(irottb)
 end if

 Sigcij_tab => Sigp%Sigcij_tab(ikcalc,1:Wfd%nsppol)

 got=0
 ABI_ALLOCATE(proc_distrb,(Wfd%mband,Kmesh%nbz,Wfd%nsppol))
 call sigma_distribution(Wfd,Kmesh,Ltg_k,Qmesh,Wfd%nsppol,can_symmetrize,kgw,Sigp%mg0,my_nbks,proc_distrb,got,global=.TRUE.)

 write(msg,'(a,i0,a)')" Will sum ",my_nbks," (b,k,s) states in Sigma_c."
 call wrtout(std_out,msg,'PERS')

 if (Sigp%gwcomp==1) then
   en_high=MAXVAL(qp_ene(Sigp%nbnds,:,:)) + Sigp%gwencomp
   write(msg,'(6a,e11.4,a)')ch10,&
&    ' Using the extrapolar approximation to accelerate convergence',ch10,&
&    ' with respect to the number of bands included',ch10,&
&    ' with extrapolar energy: ',en_high*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
   ABI_ALLOCATE(wf1swf2_g,(nspinor*gwc_nfftot))
 endif

 if (Sigp%gwcomp==1.or.dtset%gw_eet/=-1) then
   ! Setup of MPI table for extrapolar contributions. 
   ABI_ALLOCATE(extrapolar_distrb,(ib1:ib2,ib1:ib2,Kmesh%nbz,Wfd%nsppol))
   extrapolar_distrb = xmpi_undefined_rank

   do spin=1,Wfd%nsppol
     do ik_bz=1,Kmesh%nbz
        if (ANY(proc_distrb(:,ik_bz,spin) /= xmpi_undefined_rank) ) then ! This BZ point will be calculated.
           rank_mask = .FALSE. ! The set of node that will treat (k,s).
           do band=1,Wfd%mband
             rank = proc_distrb(band,ik_bz,spin)
             if (rank /= xmpi_undefined_rank) rank_mask(rank+1)=.TRUE.
           end do
           do band2=ib1,ib2
             do irow=1,Sigcij_tab(spin)%col(band2)%size1   ! Looping over the non-zero elements of sigma_ij.
               band1 = Sigcij_tab(spin)%col(band2)%bidx(irow)
               idle = imin_loc(got,mask=rank_mask)
               got(idle) = got(idle)+1
               extrapolar_distrb(band1,band2,ik_bz,spin) = idle-1
             end do
           end do
        end if
     end do
   end do

   write(msg,'(a,i0,a)')" Will treat ",COUNT(extrapolar_distrb==Wfd%my_rank)," extrapolar terms."
   call wrtout(std_out,msg,'PERS')
 end if

 ABI_ALLOCATE(rhotwg_ki,(nspinor*Sigp%npwc,minbnd:maxbnd))
 rhotwg_ki=czero_gw
 ABI_ALLOCATE(rhotwg   ,(nspinor*Sigp%npwc))
 ABI_ALLOCATE(rhotwgp  ,(nspinor*Sigp%npwc))
 ABI_ALLOCATE(vc_sqrt_qbz,(Sigp%npwc))
 
 if (Er%mqmem==0) then ! Use out-of-core solution for epsilon.
   MSG_COMMENT('Reading q-slices from file. Slower but less memory.')
 end if                                                                                !
 !
 ! === Additional allocations for PAW ===
 if (Psps%usepaw==1) then
   ABI_ALLOCATE(Cprj_ksum,(Cryst%natom,nspinor))
   call cprj_alloc(Cprj_ksum,0,Wfd%nlmn_atm)
   !
   ! For the extrapolar method we need the onsite terms of the PW in the FT mesh.
   ! * gw_gfft is the set of plane waves in the FFT Box for the oscillators.
   if (Sigp%gwcomp==1) then
     ABI_ALLOCATE(gw_gfft,(3,gwc_nfftot))
     q0=zero
     call get_gftt(gwc_ngfft,q0,Cryst%gmet,gw_gsq,gw_gfft)
     ABI_ALLOCATE(Pwij_fft,(Psps%ntypat))
     call init_paw_pwij(Pwij_fft,gwc_nfftot,(/zero,zero,zero/),gw_gfft,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
   end if
 end if ! usepaw==1
 !
 !
 if (mod10==SIG_GW_AC) then ! Calculate Gauss-Legendre quadrature knots and weights for analytic continuation
   call coeffs_gausslegint(zero,one,gl_knots,gl_wts,Er%nomega_i)
   
   do io=1,Er%nomega_i ! First frequencies are always real
     if (ABS(AIMAG(one*Er%omega(Er%nomega_r+io))-(one/gl_knots(io)-one)) > 0.0001) then
      write(msg,'(3a)')&
&       ' Frequencies in the SCR file are not compatible with the analytic continuation. ',ch10,&
&       ' Verify the frequencies in the SCR file. '
      MSG_WARNING(msg)
      if (Wfd%my_rank==Wfd%master) write(std_out,*)AIMAG(Er%omega(Er%nomega_r+io)),(one/gl_knots(io)-one)
      MSG_ERROR("Cannot continue!")
     end if
   end do
   !
   ! * To calculate \int_0^\infty domegap f(omegap), we calculate \int_0^1 dz f(1/z-1)/z^2.
   omegap(:)=one/gl_knots(:)-one
   omegap2(:)=omegap(:)*omegap(:)
   ABI_ALLOCATE(ac_epsm1cqwz2,(Sigp%npwc,Sigp%npwc,Er%nomega_i))
   istat = ABI_ALLOC_STAT
   ABI_ALLOCATE(ac_integr,(Sigp%npwc,Sigp%npwc,Sr%nomega_i))
   istat = ABI_ALLOC_STAT
 end if
 !
 ! === Calculate total number of frequencies and allocate related arrays ===
 ! * sigcme2 is used to accumulate the diagonal matrix elements over k-points and
 !   GW bands, used only in case of ppmodel 3 and 4 (TODO save memory)
 nomega_tot=Sr%nomega_r+Sr%nomega4sd
 ABI_ALLOCATE(sigcme2,(nomega_tot,ib1:ib2))
 ABI_ALLOCATE(sigcme_3,(nomega_tot))
 sigcme2=czero_gw; sigcme_3=czero_gw

 ABI_ALLOCATE(sigctmp,(nomega_sigc,Sigp%nsig_ab))
 sigctmp=czero_gw
 ABI_ALLOCATE(sigc_ket,(nspinor*Sigp%npwc,nomega_sigc))

 ! Arrays storing the contribution given by the Hermitian/anti-Hermitian part of \Sigma_c
 ABI_ALLOCATE(aherm_sigc_ket,(Sigp%npwc*nspinor,nomega_sigc))
 ABI_ALLOCATE( herm_sigc_ket,(Sigp%npwc*nspinor,nomega_sigc))

 sigcme_tmp=czero

 ABI_ALLOCATE(sigc,(2,nomega_sigc,ib1:ib2,ib1:ib2,Wfd%nsppol*Sigp%nsig_ab))
 sigc=czero

 !FIXME This quantities are only used for model GW if I am not wrong
 ABI_ALLOCATE(ket1,(Sigp%npwc*nspinor,nomega_tot))
 ABI_ALLOCATE(ket2,(Sigp%npwc*nspinor,nomega_tot))
 ABI_ALLOCATE(omegame0i,(nomega_tot))
 !
 ! Here we divide the states where the QP energies are required into complexes. Note however that this approach is not
 ! based on group theory, and it might lead to spurious results in case of accidental degeneracies.
 !
 nq_summed=Kmesh%nbz
 if (Sigp%symsigma>0) then
   call print_little_group(Ltg_k,std_out,Dtset%prtvol,'COLL')
   nq_summed=SUM(Ltg_k%ibzq(:))
   !
   ! === Find number of complexes and number of bands in each complex ===
   ! The tolerance is a little bit arbitrary (0.001 eV)
   ! It could be reduced, in particular in case of nearly accidental degeneracies
   ABI_ALLOCATE(degtab,(ib1:ib2,ib1:ib2,Wfd%nsppol))
   degtab=0
   do isppol=1,Wfd%nsppol
     do ib=ib1,ib2 
       do jb=ib1,ib2 
        if (ABS(qp_ene(ib,jk_ibz,isppol)-qp_ene(jb,jk_ibz,isppol))<0.001/Ha_ev) then
          degtab(ib,jb,isppol)=1
        end if
       end do
     end do
   end do
!   if (ANY(degtab/=0)) then ! If two states do not belong to the same complex => matrix elements of v_xc differ
!     write(msg,'(a,3f8.3,a)')' Degenerate states at k-point = ( ',kgw(:),' ).'
!     call wrtout(std_out,msg,'COLL')
!     do isppol=1,Wfd%nsppol
!       do ib=ib1,ib2 
!         do jb=ib+1,ib2 
!           if (degtab(ib,jb,isppol)==1) then
!             write(msg,'(a,i2,a,i4,a,i4)')' (spin ',isppol,')',ib,' <====> ',jb
!             call wrtout(std_out,msg,'COLL')
!             if (ABS(Sr%vxcme(ib,jk_ibz,isppol)-Sr%vxcme(jb,jk_ibz,isppol))>ABS(tol6*Sr%vxcme(jb,jk_ibz,isppol))) then 
!               write(msg,'(7a)')&
!&                ' It seems that an accidental degeneracy is occurring at this k-point ',ch10,&
!&                ' In this case, using symsigma=1 might lead to spurious results as the algorithm ',ch10,&
!&                ' will treat these states as degenerate, and it won''t be able to remove the degeneracy. ',ch10,&
!&                ' In order to avoid this deficiency, run the calculation using symsigma=0'
!               MSG_WARNING(msg)
!             end if
!           end if
!         end do
!       end do
!     end do
!   end if
 end if !symsigma

#ifdef HAVE_CLIB
 !%call clib_progress_bar(-1,Kmesh%nbz)
#endif

 write(msg,'(2a,i6,a)')ch10,' calculation status ( ',nq_summed,' to be completed):'
 call wrtout(std_out,msg,'COLL')

 ! Check for pole screening
 pole_screening = .FALSE.
 if (Er%fform==2002) pole_screening = .TRUE.

 ! Here we have a problem in case of CD since epsm1q might be huge
 ! TODO if single q (ex molecule) dont allocate epsm1q, avoid waste of memory
 if ( ANY( mod10==(/SIG_GW_AC, SIG_GW_CD, SIG_QPGW_CD/) )) then
   if (.not.(mod10==SIG_GW_CD.and.Er%mqmem==0)) then
     if (pole_screening) then
       ABI_ALLOCATE(epsm1_pole_qbz,(Sigp%npwc,Sigp%npwc,Er%ncoeff))
       istat = ABI_ALLOC_STAT
     end if
     ABI_ALLOCATE(epsm1_qbz,(Sigp%npwc,Sigp%npwc,Er%nomega))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,"out-of-memory in epsm1_qbz")
   end if
 end if

 ! TODO epsm1_trcc_qbz is needed for SIG_GW_CD with symmetries since 
 ! the Hermitian and the anti-Hermitian part have to be symmetrized 
 ! in a different way.
 !if (mod10==SIG_QPGW_CD) then 
 if (mod10==SIG_QPGW_CD.and.(.not.pole_screening)) then 
   ABI_ALLOCATE(epsm1_trcc_qbz,(Sigp%npwc,Sigp%npwc,Er%nomega))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out-of-memory in epsm1_trcc_qbz")
 end if

 ABI_ALLOCATE(igfftcg0,(Gsph_Max%ng))
 
 ABI_ALLOCATE(wfr_sum,(gwc_nfftot*nspinor))
 if (Dtset%pawcross==1) then
   ABI_ALLOCATE(igfftfcg0,(Gsph_c%ng))
   ABI_ALLOCATE(ur_ae_sum,(nfftf*nspinor))
   ABI_ALLOCATE(ur_ae_onsite_sum,(nfftf*nspinor))
   ABI_ALLOCATE(ur_ps_onsite_sum,(nfftf*nspinor))
 end if
 call timab(432,2,tsec) ! Init
 !
 ! ==========================================
 ! ==== Fat loop over k_i in the full BZ ====
 ! ==========================================
 do isppol=1,Wfd%nsppol

   if (ALL(proc_distrb(:,:,isppol)/=Wfd%my_rank).and.dtset%gw_eet==-1) CYCLE

   call timab(433,1,tsec) ! Init spin

   ! === Load wavefunctions for GW corrections ===
   ABI_ALLOCATE(wfr_bdgw,(gwc_nfftot*nspinor,ib1:ib2))
   do jb=ib1,ib2
     call wfd_get_ur(Wfd,jb,jk_ibz,isppol,wfr_sum)
     wfr_bdgw(:,jb)=wfr_sum
   end do

   if (Wfd%usepaw==1) then ! * Load cprj for GW states, note the indexing.
     dimcprj_gw=nspinor*(ib2-ib1+1)
     ABI_ALLOCATE(Cprj_kgw,(Cryst%natom,ib1:ib1+dimcprj_gw-1))
     call cprj_alloc(Cprj_kgw,0,Wfd%nlmn_atm)
     ibsp=ib1  
     do jb=ib1,ib2
       call wfd_get_cprj(Wfd,jb,jk_ibz,isppol,Cryst,Cprj_ksum,sorted=.FALSE.)
       call paw_symcprj(jk_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_ksum) 
       call cprj_copy(Cprj_ksum,Cprj_kgw(:,ibsp:ibsp+(nspinor-1)))
       ibsp=ibsp+nspinor
     end do
     if (Dtset%pawcross==1) then
       ABI_ALLOCATE(ur_ae_bdgw,(nfftf*nspinor,ib1:ib2))
       istat = ABI_ALLOC_STAT
       ABI_ALLOCATE(ur_ae_onsite_bdgw,(nfftf*nspinor,ib1:ib2))
       istat = ABI_ALLOC_STAT
       ABI_ALLOCATE(ur_ps_onsite_bdgw,(nfftf*nspinor,ib1:ib2))
       istat = ABI_ALLOC_STAT
       do jb=ib1,ib2
         call wfd_paw_get_aeur(Wfdf,jb,jk_ibz,isppol,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,&
&          ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         ur_ae_bdgw(:,jb)=ur_ae_sum
         ur_ae_onsite_bdgw(:,jb)=ur_ae_onsite_sum
         ur_ps_onsite_bdgw(:,jb)=ur_ps_onsite_sum
       end do
     end if
   end if

   call timab(433,2,tsec) ! Init spin

 do ik_bz=1,Kmesh%nbz
   !
   ! * Completeness parts should be performed only once for each k-point
   !done_once(:,:,:)=.FALSE.
   !
   ! === Parallelization over k-points and spin ===
   ! * For the spin there is another check in the inner loop
   if (ALL(proc_distrb(:,ik_bz,isppol)/=Wfd%my_rank).and.dtset%gw_eet==-1) CYCLE

   call timab(434,1,tsec) ! initq
   !
   ! * Find the corresponding irreducible k-point
   call get_BZ_item(Kmesh,ik_bz,ksum,ik_ibz,isym_ki,iik,ph_mkt)
   spinrot_kbz(:)=Cryst%spinrot(:,isym_ki)
   npw_k =  Wfd%Kdata(ik_ibz)%npw
   kg_k  => Wfd%Kdata(ik_ibz)%kg_k

   ! * Identify q and G0 where q+G0=k_GW-k_i
   kgw_m_ksum=kgw-ksum
   call findqg0(iq_bz,g0,kgw_m_ksum,Qmesh%nbz,Qmesh%bz,Sigp%mG0)

   ! === If symsigma, symmetrize the matrix elements ===
   ! * Sum only q"s in IBZ_k. In this case elements are weighted
   !   according to wtqp and wtqm. wtqm is for time-reversal.
   wtqp=1; wtqm=0
   !if (Sigp%symsigma>0) then
   if (can_symmetrize(isppol)) then
     if (Ltg_k%ibzq(iq_bz)/=1) CYCLE
     wtqp=0; wtqm=0
     do isym=1,Ltg_k%nsym_sg
       wtqp=wtqp+Ltg_k%wtksym(1,isym,iq_bz)
       wtqm=wtqm+Ltg_k%wtksym(2,isym,iq_bz)
     end do
   end if

#ifdef HAVE_CLIB
   !%call clib_progress_bar(ik_bz,Kmesh%nbz)
#else
   !%write(msg,'(2(a,i4),a,i3)')' csigme : ik_bz ',ik_bz,'/',Kmesh%nbz,' done by processor ',Wfd%my_rank
   !%call wrtout(std_out,msg,'PERS')
#endif
   !
   ! === Find the corresponding irred q-point ===
   call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
   q_is_gamma = (normv(qbz,Cryst%gmet,"G") < GW_TOL_W0)

   !q_is_gamma = (normv(qbz,Cryst%gmet,"G") < 0.7)
   !if (iq_ibz/=2.and.iq_ibz/=1) CYCLE
   !if (ANY(qbz<=-(half-tol16)) .or. ANY(qbz>(half+tol16))) CYCLE
   !if (q_is_gamma) then; write(std_out,*)"skipping q=Gamma"; CYCLE; end if
   !
   ! Tables for the FFT of the oscillators.
   !  a) FFT index of the G-G0.
   !  b) gw_gbound table for the zero-padded FFT performed in rhotwg. 
   ABI_ALLOCATE(gw_gbound,(2*gwc_mgfft+8,2))
   call gsph_fft_tabs(Gsph_c,g0,gwc_mgfft,gwc_ngfft,use_padfft,gw_gbound,igfftcg0)
   if ( ANY(gwc_fftalga == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
   if (use_padfft==0) then 
     ABI_DEALLOCATE(gw_gbound)
     ABI_ALLOCATE(gw_gbound,(2*gwc_mgfft+8,2*use_padfft))
   end if
   if (Dtset%pawcross==1) then
     ABI_ALLOCATE(gboundf,(2*mgfftf+8,2))
     call gsph_fft_tabs(Gsph_c,g0,mgfftf,rho_ngfft,use_padfftf,gboundf,igfftfcg0)
     if ( ANY(gwc_fftalga == (/2,4/)) ) use_padfftf=0
     if (use_padfftf==0) then 
       ABI_DEALLOCATE(gboundf)
       ABI_ALLOCATE(gboundf,(2*mgfftf+8,2*use_padfftf))
     end if
   end if
   !
   ! === Evaluate oscillator matrix elements ===
   ! * $ <phj/r|e^{-i(q+G)}|phi/r> - <tphj/r|e^{-i(q+G)}|tphi/r> $ in packed form.
   if (Psps%usepaw==1) then
     ABI_ALLOCATE(Pwij_qg,(Psps%ntypat))
     q0 = qbz !;if (q_is_gamma) q0 = (/0.00001_dp,0.00001_dp,0.00001_dp/) ! GW_Q0_DEFAULT
     call init_paw_pwij(Pwij_qg,Sigp%npwc,q0,Gsph_c%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
   end if

   if (Er%mqmem==0) then ! Read q-slice of epsilon^{-1}|chi0 in Er%epsm1(:,:,:,1) (much slower but less memory).
     if (pole_screening) then
       call get_pole_epsm1(Er,Vcp,0,0,Dtset%accesswff,xmpi_self,iqibzA=iq_ibz,&
&       reconstruct_scr=Dtset%gw_reconst_scr)
       if (Dtset%gw_reconst_scr==1) pole_screening=.FALSE.
     else
       call get_epsm1(Er,Vcp,0,0,Dtset%accesswff,xmpi_self,iqibzA=iq_ibz)
     end if
     if (sigma_needs_ppm(Sigp)) then
       if (Wfd%usepaw==1.and.PPm%userho==1) then ! Use PAW AE rhor.
         call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,rho_nfftot,Gsph_c%gvec,&
&          rho_ngfft,aepaw_rhor(:,1),iqiA=iq_ibz)
       else 
         call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,rho_nfftot,Gsph_c%gvec,&
&          rho_ngfft,rhor(:,1),iqiA=iq_ibz)
       end if
     end if
   end if
   !
   ! === Symmetrize PPM parameters and epsm1 (q_IBZ --> q_BZ) ===
   ! * NOTE: We are not considering umklapp with G0/=0. In this case,
   !   indeed the equation is different since we have to use G-G0.
   !   A check, however, is performed in sigma.
   ! * If gwcomp==1 and mod10=1,2,9, one needs both to set up botsq and epsm1_q
   !
   if (sigma_needs_ppm(Sigp)) then
     !if(isym_q==1.and.itim_q==1) then ! If identity, there is no need for rotating the G's.
     !  iq_curr=(iq_ibz-1)*MIN(Er%mqmem,1)+1
     !  botsq%datum => PPm%bigomegatwsq(:,:,iq_curr)
     !  otq%datum   => PPm%omegatw(:,:,iq_curr)
     !  eig%datum   => PPm%eigpot(:,:,iq_curr)
     !else
     !dt_shape(1,:)=1; dt_shape(2,1)=PPm%npwc;   dt_shape(2,2)=PPm%dm2_botsq; call allocate_jpt(botsq,dt_shape,istat)
     ABI_ALLOCATE(botsq,(PPm%npwc,PPm%dm2_botsq))
     !dt_shape(1,:)=1; dt_shape(2,1)=PPm%npwc;   dt_shape(2,2)=PPm%dm2_otq;   call allocate_jpt(otq,dt_shape,istat)
     ABI_ALLOCATE(otq,(PPm%npwc,PPm%dm2_otq))
     !dt_shape(1,:)=1; dt_shape(2,1)=PPm%dm_eig; dt_shape(2,2)=PPm%dm_eig;    call allocate_jpt(eig,dt_shape,istat)
     ABI_ALLOCATE(eig,(PPm%dm_eig,PPm%dm_eig))
     call ppm_get_qbz(PPm,Gsph_c,Qmesh,iq_bz,botsq,otq,eig)
     !end if
   end if

   if ( ANY(mod10==(/SIG_GW_AC,SIG_GW_CD,SIG_QPGW_CD/) )) then
     ! === Numerical integration or model GW with contour deformation or Analytic Continuation ===
     !  TODO In case of AC we should symmetrize only the imaginary frequencies
     if (mod10==SIG_GW_CD.and.Er%mqmem==0) then ! Do in-place symmetrisation
       if (pole_screening) then
         call Epsm1_pole_symmetrizer_inplace(iq_bz,Er%nomega,Sigp%npwc,Er,Gsph_c,Qmesh)
       else
        call Epsm1_symmetrizer_inplace(iq_bz,Er%nomega,Sigp%npwc,Er,Gsph_c,Qmesh,.TRUE.)
       end if
     else
       if (pole_screening) then
         call Epsm1_pole_symmetrizer(iq_bz,Er%ncoeff,Sigp%npwc,Er,Gsph_c,Qmesh,epsm1_pole_qbz)
       else
        call Epsm1_symmetrizer(iq_bz,Er%nomega,Sigp%npwc,Er,Gsph_c,Qmesh,.TRUE.,epsm1_qbz)
       end if
     end if
!     if (pole_screening) then ! Reconstruct epsm1
!       do ig2=1,Sigp%npwc
!         do ig1=1,Sigp%npwc
!           call re_and_im_screening_with_phase(Er%omega,epsm1_qbz(ig1,ig2,:),Er%nomega, &
!&           epsm1_pole_qbz(ig1,ig2,:),Er%ncoeff)
!         end do
!       end do
!     end if 

     !write(std_out,*) ch10,' epsm1_qbz(:,:,1): ',epsm1_qbz(:,:,1)

     if (mod10==SIG_GW_AC) then ! Prepare the first term: \Sum w_i 1/z_i^2 f(1/z_i-1)..
       ! * Note that the first frequencies are always real, skip them.
       ! * Memory is not optimized.
       do iiw=1,Er%nomega_i
         z2=gl_knots(iiw)*gl_knots(iiw)
         ac_epsm1cqwz2(:,:,iiw)= gl_wts(iiw)*epsm1_qbz(:,:,Er%nomega_r+iiw)/z2
       end do
     end if

     !if (mod10==SIG_QPGW_CD) then ! For model GW we need transpose(conjg(epsm1_qbz)) ===
     if (mod10==SIG_QPGW_CD) then 
       do io=1,Er%nomega
        epsm1_trcc_qbz(:,:,io)=TRANSPOSE(CONJG(epsm1_qbz(:,:,io)))
       end do
     end if
   end if !gwcalctyp
   !
   ! === Get Fourier components of the Coulombian interaction in the BZ ===
   ! * In 3D systems, neglecting umklapp,  vc(Sq,sG)=vc(q,G)=4pi/|q+G|
   ! * The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
   do ig=1,Sigp%npwc
     vc_sqrt_qbz(Gsph_c%rottb(ig,itim_q,isym_q))=Vcp%vc_sqrt(ig,iq_ibz)
   end do

   call timab(434,2,tsec) ! initq

   ! TODO check gw_eet_sigma as now igfftg0 and oscillatorm matrix elements are 
   ! given on the ecuteps sphere.

   nbmax=Sigp%nbnds
   if (dtset%gw_eet/=-1) then ! MG Where are the GW wavefunctions? It won't work in parallel
     call timab(435,1,tsec) ! initq
     call gw_eet_sigma(Sigp,Sr,Dtset,Cryst,Wfd,Kmesh,Qmesh,Gsph_Max,Gsph_c,Psps,Vcp,QP_BSt,PPm, &
&                      isppol,iq_bz,ik_bz,jk_bz,ik_ibz,jk_ibz,itim_q,isym_q,iq_ibz,&
&                      ktabr(:,ik_bz),ktabr(:,jk_bz),spinrot_kbz,spinrot_kgw,ph_mkt,ph_mkgwt, &
&                      gwc_nfftot,gwc_ngfft,use_padfft,igfftcg0,gw_gbound,gwc_mgfft,ib1,ib2, &
&                      nomega_tot,nomega_sigc,fact_sp,nspinor,botsq,otq,sigcme_tmp,sigc,nbopt, &
&                      tim_fourdp2,wtqp,wtqm,MPI_enreg_seq,extrapolar_distrb,can_symmetrize)
     nbmax=nbopt
     call timab(435,2,tsec) ! initq
   end if
   !
   ! === Sum over bands ===
   !do isppol=1,Wfd%nsppol
     call timab(445,1,tsec) ! loop

     do ib=1,nbmax

       call timab(436,1,tsec) ! (1)
       !
       ! === Parallelism over spin ===
       ! * This processor has this k-point but what about spin?
       if (proc_distrb(ib,ik_bz,isppol)/=Wfd%my_rank) CYCLE
       !
       ! * Skip empty state ib for HF, SEX, and COHSEX.
       !if (qp_occ(ib,ik_ibz,isppol)<tol_empty.and.(ANY(mod10==(/SIG_HF/))) ) CYCLE

       call wfd_get_ur(Wfd,ib,ik_ibz,isppol,wfr_sum)

       if (Psps%usepaw==1) then ! Load cprj for point ksum, this spin or spinor and *THIS* band.
         ! TODO MG I could avoid doing this but I have to exchange spin and bands ???
         ! For sure there is a better way to do this!
         call wfd_get_cprj(Wfd,ib,ik_ibz,isppol,Cryst,Cprj_ksum,sorted=.FALSE.)
         call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_ksum) 
         if (Dtset%pawcross==1) then
           call wfd_paw_get_aeur(Wfdf,ib,ik_ibz,isppol,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,&
&              ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         end if
       end if

       if (mod10==SIG_GW_AC) then ! Calculate integral over omegap with Gauss-Legendre quadrature.
         ac_integr(:,:,:)=czero_gw
         ! * -1/pi \int domegap epsm1c*(omega-e0i) / ( (omega-e0i)^2 + omegap^2)
         ! * Note that energies are calculated wrt the Fermi level.
         do io=1,Sr%nomega_i
           omegame0i_ac  = Sr%omega_i(io)-qp_ene(ib,ik_ibz,isppol)
           omegame0i2_ac = omegame0i_ac*omegame0i_ac
           do iiw=1,Er%nomega_i
             do iggp=0,Sigp%npwc*Sigp%npwc-1
               ig=iggp/Sigp%npwc+1
               igp= iggp-(ig-1)*Sigp%npwc+1 ! \int domegap epsm1c/((omega-e0i)^2 + omegap^2)
               ac_integr(ig,igp,io)= ac_integr(ig,igp,io) + ac_epsm1cqwz2(ig,igp,iiw)/(omegame0i2_ac + omegap2(iiw))
             end do
           end do
           ac_integr(:,:,io)=ac_integr(:,:,io)*omegame0i_ac
         end do
         ac_integr(:,:,:)=-ac_integr(:,:,:)*piinv
       end if

       call timab(436,2,tsec) ! (1)
       call timab(437,1,tsec) ! rho_tw_g

       do jb=ib1,ib2 ! Get all <k-q,ib,s|e^{-i(q+G).r}|s,jb,k>, at once ===

         call rho_tw_g(Wfd%paral_kgb,nspinor,Sigp%npwc,gwc_nfftot,gwc_ngfft,1,use_padfft,igfftcg0,gw_gbound,&
&          wfr_sum       ,iik,ktabr(:,ik_bz),ph_mkt  ,spinrot_kbz,  &
&          wfr_bdgw(:,jb),jik,ktabr(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&          nspinor,rhotwg_ki(:,jb),tim_fourdp2,Wfd%MPI_enreg)

         if (Psps%usepaw==1) then ! Add on-site contribution, projectors are already in BZ !TODO Recheck this!
           i2=jb; if (nspinor==2) i2=(2*jb-1)
           spad=(nspinor-1)
           call paw_rho_tw_g(Sigp%npwc,nspinor,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_c%gvec,&
&            Cprj_ksum(:,:),Cprj_kgw(:,i2:i2+spad),Pwij_qg,rhotwg_ki(:,jb))
           if (Dtset%pawcross==1) then ! Add paw cross term
             call paw_cross_rho_tw_g(Wfdf%paral_kgb,nspinor,Sigp%npwc,nfftf,rho_ngfft,1,use_padfftf,igfftfcg0,gboundf,&
&             ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum,iik,ktabrf(:,ik_bz),ph_mkt,spinrot_kbz,&
&             ur_ae_bdgw(:,jb),ur_ae_onsite_bdgw(:,jb),ur_ps_onsite_bdgw(:,jb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&             nspinor,rhotwg_ki(:,jb),tim_fourdp2,Wfdf%MPI_enreg)
           end if
         end if
         !
         ! === Multiply by the square root of the Coulomb term ===
         ! * In 3-D systems, the factor sqrt(4pi) is included)
         do ii=1,nspinor
           spad=(ii-1)*Sigp%npwc
           rhotwg_ki(spad+1:spad+Sigp%npwc,jb) = rhotwg_ki(spad+1:spad+Sigp%npwc,jb)*vc_sqrt_qbz(1:Sigp%npwc)
         end do
         !
         ! === Treat analytically the case q --> 0 ===
         ! * The oscillator is evaluated at q=O as it is considered constant in the small cube around Gamma
         !   while the Colulomb term is integrated out out out out.
         ! * In the scalar case we have nonzero contribution only if ib==jb
         ! * For nspinor==2 evalute <ib,up|jb,up> and <ib,dwn|jb,dwn>,
         !   impose orthonormalization since npwwfn might be < npwvec.
         if (ik_bz==jk_bz) then
           if (nspinor==1) then
             rhotwg_ki(1,jb)=czero_gw
             if (ib==jb) rhotwg_ki(1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)

           else ! TODO Recheck this, moreover it wont work if k-centered G-spheres are used.!
             cg_sum  => Wfd%Wave(ib,ik_ibz,isppol)%ug
             cg_jb   => Wfd%Wave(jb,jk_ibz,isppol)%ug
             ctmp = xdotc(Wfd%npwwfn*nspinor,cg_sum,1,cg_jb,1) 
             ovlp(1) = REAL(ctmp) 
             ovlp(2) = AIMAG(ctmp) 
             if (Psps%usepaw==1) then
               i2=(2*jb-1)
               ovlp = ovlp + paw_overlap(Cprj_ksum,Cprj_kgw(:,i2:i2+1),Cryst%typat,Pawtab)
             end if
             !ovlp(2) = -ovlp(1)
             !if (ib==jb) ovlp(2)=cone_gw-ovlp(1)
             if (ib==jb) then
               norm=DBLE(ovlp(1)+ovlp(2))
               ovlp(1)=DBLE(ovlp(1)/norm)
               ovlp(2)=DBLE(ovlp(2)/norm)
             else
               scprod=ovlp(1)+ovlp(2)
               ovlp(1)=ovlp(1)-scprod*half
               ovlp(2)=ovlp(2)-scprod*half
             end if
             rhotwg_ki(1          ,jb) = CMPLX(SQRT(Vcp%i_sz),0.0_gwp)*ovlp(1)
             rhotwg_ki(Sigp%npwc+1,jb) = CMPLX(SQRT(Vcp%i_sz),0.0_gwp)*ovlp(2)
           end if
         end if
       end do !jb  Got all matrix elements from minbnd up to maxbnd.

       theta_mu_minus_e0i=fact_sp*qp_occ(ib,ik_ibz,isppol)

       ! Starting point to evaluate the derivative of Sigma and the Spectral function
       e0i=qp_ene(ib,ik_ibz,isppol)

       ! Frequencies for the spectral function, e0i=qp_ene(ib,ik_ibz,isppol)
       ! FIXME the interval is not centered on eoi ! WHY?
       if (Sr%nomega_r>0) omegame0i(1:Sr%nomega_r)=DBLE(Sr%omega_r(1:Sr%nomega_r))-e0i

       call timab(437,2,tsec) ! rho_tw_g

       do kb=ib1,ib2

         call timab(438,1,tsec) ! (2)
         !
         ! Get frequencies $\omega$-\epsilon_in$ to evaluate  $d\Sigma/dE$, note the spin
         ! subtract e_KS since we have stored e_KS+ Delta \omega in Sr%omega4sd, not required for AC
         do io=Sr%nomega_r+1,nomega_tot
           omegame0i(io)=DBLE(Sr%omega4sd(kb,jk_ibz,io-Sr%nomega_r,isppol))-e0i
         end do
         !
         ! === Get the ket \Sigma|\phi_{k,kb}> according to the method ===
         rhotwgp(:)=rhotwg_ki(:,kb)
         sigc_ket  = czero_gw
         ket1      = czero_gw
         ket2      = czero_gw

         !aherm_sigc_ket = czero_gw 
         ! herm_sigc_ket = czero_gw 

         SELECT CASE (mod10)

         CASE (SIG_GW_PPM) ! GW WITH Plasmon-Pole Model.
           !
           ! * Note that ppmodel 3 or 4 work only in case of standard perturbative approach!
           !   Moreover, for ppmodel 3 and 4, spinorial case is not allowed
           call calc_sig_ppm(PPm,nspinor,Sigp%npwc,nomega_tot,rhotwgp,botsq,otq,&
&           omegame0i,Sigp%zcut,theta_mu_minus_e0i,eig,Sigp%npwc,sigc_ket,sigcme_3)

           if (PPm%model==3.or.PPm%model==4) then
             sigcme2(:,kb)=sigcme2(:,kb) + (wtqp+wtqm)*DBLE(sigcme_3(:)) + (wtqp-wtqm)*j_gw*AIMAG(sigcme_3(:))
           end if

         CASE (SIG_GW_AC) ! GW with Analytic continuation.
           !
           ! * Evaluate \sum_Gp integr_GGp(omegasi) rhotw_Gp TODO this part can be optimized
           do io=1,Sr%nomega_i
             do ispinor=1,nspinor
               spadc=(ispinor-1)*Sigp%npwc

               do ig=1,Sigp%npwc
                 ctmp=czero
                 do igp=1,Sigp%npwc
                   ctmp=ctmp+ac_integr(ig,igp,io)*rhotwgp(igp+spadc)
                 end do
                 sigc_ket(ig+spadc,io)=ctmp
               end do

             end do !ispinor
           end do !io

         CASE (SIG_GW_CD) ! GW with contour deformation.
           ! Check if pole contributions need to be summed (this avoids unnecessary
           ! splint calls and saves time)
           !me_calc_poles = .TRUE.
           do io=1,nomega_tot
             if (omegame0i(io)>=zero.AND.(ABS(one-theta_mu_minus_e0i)>zero)) then
               !me_calc_poles(io) = .TRUE.
               if ( w_maxval(kb) < ABS(omegame0i(io)) ) w_maxval(kb) = ABS(omegame0i(io))
             else if (omegame0i(io)<zero.AND.(ABS(theta_mu_minus_e0i)>zero)) then
               !me_calc_poles(io) = .TRUE.
               if ( w_maxval(kb) < ABS(omegame0i(io)) ) w_maxval(kb) = ABS(omegame0i(io))
             end if 
           end do
           !
           ! Check memory saving
           if (Er%mqmem==0) then
             call calc_sigc_cd(Sigp%npwc,Sigp%npwc,nspinor,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&              Er%omega,Er%epsm1(:,:,:,1),omegame0i,theta_mu_minus_e0i,sigc_ket,npoles_missing(kb))
           else
             call calc_sigc_cd(Sigp%npwc,Sigp%npwc,nspinor,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&              Er%omega,epsm1_qbz,omegame0i,theta_mu_minus_e0i,sigc_ket,npoles_missing(kb))
           end if
           !else
           !  if (Er%mqmem==0) then
           !    call calc_sigc_pole_cd(Sigp%npwc,Sigp%npwc,nspinor,Er%ncoeff,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
!&                Er%omega,Er%epsm1_pole(:,:,:,1),omegame0i,theta_mu_minus_e0i,sigc_ket,npoles_missing(kb),calc_poles=me_calc_poles)
           !  else
           !    call calc_sigc_pole_cd(Sigp%npwc,Sigp%npwc,nspinor,Er%ncoeff,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
!&                Er%omega,epsm1_pole_qbz,omegame0i,theta_mu_minus_e0i,sigc_ket,npoles_missing(kb),calc_poles=me_calc_poles)
           !  end if

#if 0
           if (wtqm/=0) then
             call calc_sigc_cd(Sigp%npwc,Sigp%npwc,nspinor,,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&              Er%omega,epsm1_trcc_qbz,omegame0i,theta_mu_minus_e0i,aherm_sigc_ket,npoles_missing(kb))

              herm_sigc_ket = half*(sigc_ket + aherm_sigc_ket)
             aherm_sigc_ket = half*(sigc_ket - aherm_sigc_ket)
           else
             herm_sigc_ket  = sigc_ket
             aherm_sigc_ket = czero_gw
           end if
#endif

         CASE (SIG_QPGW_PPM) ! MODEL GW calculation WITH PPm  TODO Spinor not tested.
           !
           ! * Calculate \Sigma(E_k) |k> to obtain <j|\Sigma(E_k)|k>
           ABI_ALLOCATE(sigcme_new,(nomega_tot))

           call calc_sig_ppm(PPm,nspinor,Sigp%npwc,nomega_tot,rhotwgp,botsq,otq,&
&            omegame0i,Sigp%zcut,theta_mu_minus_e0i,eig,Sigp%npwc,ket1,sigcme_new)

           if (Sigp%gwcalctyp==28) then
             if (PPm%model/=1.and.PPm%model/=2) then ! This is needed to have Sigp%npwc=PPm%dm2_botsq=PPm%dm2_otq
               write(msg,'(3a)')&
&                ' For the time being, gwcalctyp=28 cannot be used with ppmodel=3,4.',ch10,&
&                ' Use another Plasmon Pole Model when gwcalctyp=28.'
               MSG_ERROR(msg)
             end if
             ABI_ALLOCATE(botsq_conjg_transp,(PPm%dm2_botsq,Sigp%npwc))
             botsq_conjg_transp=TRANSPOSE(botsq) ! Keep these two lines separated, otherwise gfortran messes up
             botsq_conjg_transp=CONJG(botsq_conjg_transp)
             ABI_ALLOCATE(otq_transp,(PPm%dm2_otq,PPm%npwc))
             istat = ABI_ALLOC_STAT
             otq_transp=TRANSPOSE(otq)

             call calc_sig_ppm(PPm,nspinor,Sigp%npwc,nomega_tot,rhotwgp,botsq_conjg_transp,otq_transp,&
&              omegame0i,Sigp%zcut,theta_mu_minus_e0i,eig,Sigp%npwc,ket2,sigcme_3)

             ABI_DEALLOCATE(botsq_conjg_transp)
             ABI_DEALLOCATE(otq_transp)
             sigc_ket= half*(ket1+ket2)
           else
             sigc_ket= ket1
           end if

           ABI_DEALLOCATE(sigcme_new)

         CASE (SIG_QPGW_CD) ! MODEL GW with numerical integration.
           ! Check if pole contributions need to be summed (this avoids unnecessary
           ! splint calls and saves time)
           !me_calc_poles = .TRUE.
           do io=1,nomega_tot
             if (omegame0i(io)>=zero.AND.(ABS(one-theta_mu_minus_e0i)>zero)) then
               !me_calc_poles(io) = .TRUE.
               if ( w_maxval(kb) < ABS(omegame0i(io)) ) w_maxval(kb) = ABS(omegame0i(io))
             else if (omegame0i(io)<zero.AND.(ABS(theta_mu_minus_e0i)>zero)) then
               !me_calc_poles(io) = .TRUE.
               if ( w_maxval(kb) < ABS(omegame0i(io)) ) w_maxval(kb) = ABS(omegame0i(io))
             end if 
           end do
           !
           ! * Calculate \Sigma(E_k)|k> to obtain <j|\Sigma(E_k)|k>
           call calc_sigc_cd(Sigp%npwc,Sigp%npwc,nspinor,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&            Er%omega,epsm1_qbz,omegame0i,theta_mu_minus_e0i,ket1,npoles_missing(kb))

           if (Sigp%gwcalctyp==29) then ! Calculate \Sigma^*(E_k)|k> to obtain <k|\Sigma(E_k)|j>^*
             call calc_sigc_cd(Sigp%npwc,Sigp%npwc,nspinor,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&              Er%omega,epsm1_trcc_qbz,omegame0i,theta_mu_minus_e0i,ket2,npoles_missing(kb))
             sigc_ket = half*(ket1+ket2)
           else
             sigc_ket = ket1
           end if

         CASE DEFAULT
           write(msg,'(a,i0)')" unsupported value for mod10: ",mod10
           MSG_ERROR(msg)
         END SELECT
         
         if (Sigp%gwcomp==1) then ! TODO spinor not implemented
           call calc_sig_ppm_comp(Sigp%npwc,nomega_tot,rhotwgp,botsq,otq,DBLE(Sr%egw(kb,jk_ibz,isppol)-en_high),&
&            Sigp%zcut,theta_mu_minus_e0i,sigc_ket,PPm%model,Sigp%npwc,PPm%dm2_botsq,PPm%dm2_otq)
         end if
 
         call timab(438,2,tsec) ! 
         call timab(439,1,tsec) ! sigma_me
         !
         ! Loop over the non-zero row elements of this column.
         ! 1) If gwcalctyp<20 : only diagonal elements since QP==KS.
         ! 2) If gwcalctyp>=20:
         !     * Only off-diagonal elements connecting states with same character.
         do irow=1,Sigcij_tab(isppol)%col(kb)%size1

           jb = Sigcij_tab(isppol)%col(kb)%bidx(irow)
           rhotwg=rhotwg_ki(:,jb)
           !
           ! === Calculate <\phi_j|\Sigma_c|\phi_k> ===
           ! * Different freqs according to method (AC or Perturbative), see nomega_sigc.

           do iab=1,Sigp%nsig_ab
             spadc1=spinor_padc(1,iab)
             spadc2=spinor_padc(2,iab)
             do io=1,nomega_sigc
               sigctmp(io,iab)=XDOTC(Sigp%npwc,rhotwg(spadc1+1:),1,sigc_ket(spadc2+1:,io),1)
             end do
           end do
           !
           
           if (Sigp%gwcomp==1) then ! Evaluate Extrapolar term TODO this does not work with spinor
             
             !if ( (.not.done_once(jb,kb,isppol)) .and. (gwpara/=2 .or. Wfd%my_rank==modulo(jb+kb*(ib2-ib1+1),Wfd%nproc) ) ) then
             !  done_once(jb,kb,isppol)=.TRUE.
             if ( extrapolar_distrb(jb,kb,ik_bz,isppol) == Wfd%my_rank ) then
               extrapolar_distrb(jb,kb,ik_bz,isppol) = xmpi_undefined_rank ! Do it once as it does not depend on the ib index summed over.
#if 1
               call calc_wfwfg(Wfd%MPI_enreg,Wfd%paral_kgb,tim_fourdp2,ktabr(:,jk_ibz),jik,& ! why jk_ibz?
&                gwc_nfftot,gwc_ngfft,wfr_bdgw(:,jb),wfr_bdgw(:,kb),wf1swf2_g)
#else
               call calc_wfwfg(Wfd%MPI_enreg,Wfd%paral_kgb,tim_fourdp2,ktabr(:,jk_bz),jik,&
&                gwc_nfftot,gwc_ngfft,wfr_bdgw(:,jb),wfr_bdgw(:,kb),wf1swf2_g)
#endif

               if (Psps%usepaw==1) then
                 i1=jb
                 i2=kb
                 if (nspinor==2) then
                   i1=(2*jb-1)
                   i2=(2*kb-1)
                 end if
                 spad=(nspinor-1)
                 call paw_rho_tw_g(gwc_nfftot,Sigp%nsig_ab,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,&
&                  gw_gfft,Cprj_kgw(:,i1:i1+spad),Cprj_kgw(:,i2:i2+spad),Pwij_fft,wf1swf2_g)

                 if (Dtset%pawcross==1) then ! Add paw cross term
                   call paw_cross_rho_tw_g(Wfdf%paral_kgb,nspinor,Sigp%npwc,nfftf,rho_ngfft,1,use_padfftf,igfftfcg0,gboundf,&
&                   ur_ae_bdgw(:,jb),ur_ae_onsite_bdgw(:,jb),ur_ps_onsite_bdgw(:,jb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&                   ur_ae_bdgw(:,kb),ur_ae_onsite_bdgw(:,kb),ur_ps_onsite_bdgw(:,kb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&                   nspinor,wf1swf2_g,tim_fourdp2,Wfdf%MPI_enreg)
                 end if
               end if
               !
               ! === The static contribution from completeness relation is calculated once ===
               call calc_coh_comp(iq_ibz,Vcp%i_sz,(jb==kb),nspinor,Sigp%nsig_ab,DBLE(Sr%egw(kb,jk_ibz,isppol)-en_high),&
&                Sigp%npwc,Gsph_c%gvec,gwc_ngfft,gwc_nfftot,wf1swf2_g,vc_sqrt_qbz,botsq,otq,sigcohme)

               do io=1,nomega_sigc
                 sigctmp(io,:) = sigctmp(io,:)+sigcohme(:)
               end do
             end if ! gwcomp==1
           end if ! gwcom==1
           !
           ! === Accumulate and, in case, symmetrize matrix elements of Sigma_c ===
           do iab=1,Sigp%nsig_ab
             is_idx=isppol; if (nspinor==2) is_idx=iab

             sigcme_tmp(:,jb,kb,is_idx)=sigcme_tmp(:,jb,kb,is_idx) + &
&              (wtqp+wtqm)*DBLE(sigctmp(:,iab)) + (wtqp-wtqm)*j_gw*AIMAG(sigctmp(:,iab))

             sigc(1,:,jb,kb,is_idx)=sigc(1,:,jb,kb,is_idx) + wtqp*      sigctmp(:,iab)
             sigc(2,:,jb,kb,is_idx)=sigc(2,:,jb,kb,is_idx) + wtqm*CONJG(sigctmp(:,iab))
             ! TODO this should be the contribution coming from the anti-hermitian part.
           end do
         end do !jb used to calculate matrix elements of $\Sigma$

         ! shaltaf (030406): this has to be done in a clean way later. TODO does not work with spinor.
         if (mod10==SIG_GW_PPM.and.(PPm%model==3.or.PPm%model==4)) then
           sigcme_tmp(:,kb,kb,isppol)= sigcme2(:,kb)
           sigc(1,:,kb,kb,isppol)= sigcme2(:,kb)
           sigc(2,:,kb,kb,isppol)= czero
         end if

         call timab(439,2,tsec) ! csigme(SigC)

       end do !kb to calculate matrix elements of $\Sigma$
     end do !ib

     call timab(445,2,tsec) ! csigme(SigC)

   !end do !isppol
   !
   ! Deallocate k-dependent quantities.
   ABI_DEALLOCATE(gw_gbound)
   istat = ABI_ALLOC_STAT
   if (Dtset%pawcross==1) then
     ABI_DEALLOCATE(gboundf)
     istat = ABI_ALLOC_STAT
   end if

   if (sigma_needs_ppm(Sigp)) then
     !call destroy_jpt(botsq,istat)
     ABI_DEALLOCATE(botsq)
     istat = ABI_ALLOC_STAT
     !call destroy_jpt(otq,istat)
     ABI_DEALLOCATE(otq)
     istat = ABI_ALLOC_STAT
     !call destroy_jpt(eig,istat)
     ABI_DEALLOCATE(eig)
     istat = ABI_ALLOC_STAT
   end if
   if (Psps%usepaw==1) then
     call destroy_paw_pwij(Pwij_qg)
     ABI_DEALLOCATE(Pwij_qg)
   end if
 end do !ik_bz 

   ABI_DEALLOCATE(wfr_bdgw)
   if (Wfd%usepaw==1) then
     call cprj_free(Cprj_kgw )
     ABI_DEALLOCATE(Cprj_kgw)
     if (Dtset%pawcross==1) then
       ABI_DEALLOCATE(ur_ae_bdgw)
       ABI_DEALLOCATE(ur_ae_onsite_bdgw)
       ABI_DEALLOCATE(ur_ps_onsite_bdgw)
     end if
   end if
 end do !isppol

 ABI_DEALLOCATE(sigcme2)
 ABI_DEALLOCATE(sigcme_3)
 ABI_DEALLOCATE(igfftcg0)
 if (Dtset%pawcross==1) then
   ABI_DEALLOCATE(igfftfcg0)
 end if
 !
 call timab(440,1,tsec) ! wfd_barrier

 ! === Gather contributions from all the CPUs ===
 call wfd_barrier(Wfd)
 call timab(440,2,tsec) ! wfd_barrier
 call timab(441,1,tsec) ! xsum_mpi

 call xsum_mpi(sigcme_tmp,comm,ierr)
 call xsum_mpi(sigc,comm,ierr)

 call timab(441,2,tsec) ! xsum_mpi
 !
 ! === Multiply by constants ===
 ! * For 3D systems sqrt(4pi) is included in vc_sqrt_qbz ===
 sigcme_tmp = sigcme_tmp /(Cryst%ucvol*Kmesh%nbz)
 sigc       = sigc       /(Cryst%ucvol*Kmesh%nbz)
 !
 ! === If we have summed over the IBZ_q now we have to average over complexes ===
 ! * Presently only diagonal terms are considered
 ! * TODO QP-SCGW required a more involved approach, there is a check in sigma
 ! * TODO it does not work if nspinor==2.
 call timab(442,1,tsec) ! final ops
 do spin=1,Wfd%nsppol

   if (can_symmetrize(spin)) then
     if (mod10==SIG_GW_AC) then ! FIXME here there is a problem in case of AC with symmetries
       ABI_ALLOCATE(sym_cme,(Sr%nomega_i,ib1:ib2,ib1:ib2))
     else
       ABI_ALLOCATE(sym_cme,(nomega_tot,ib1:ib2,ib1:ib2))
     end if
     sym_cme=czero
     !
     ! === Average over degenerate diagonal elements ===
     ! NOTE: frequencies for \Sigma_c(\omega) should be equal to avoid spurious results.
     ! another good reason to use a strict criterion for the tolerance on eigenvalues.
     do ib=ib1,ib2 
       ndegs=0
       do jb=ib1,ib2 
         if (degtab(ib,jb,spin)==1) then
           sym_cme(:,ib,ib)=sym_cme(:,ib,ib)+SUM(sigc(:,:,jb,jb,spin),DIM=1)
         end if
         ndegs=ndegs+degtab(ib,jb,spin)
       end do
       sym_cme(:,ib,ib)=sym_cme(:,ib,ib)/ndegs
     end do

     if (Sigp%gwcalctyp >= 20) then 
       do iwc=1,nomega_sigc
         call symmetrize_me(QP_sym(spin),ib1,ib2,sigc(:,iwc,:,:,spin),sym_cme(iwc,:,:))
       end do
     end if
     !
     ! ==== Copy symmetrized values ====
     do ib=ib1,ib2 
       do jb=ib1,ib2 
         !if (mod10==SIG_GW_AC.and.average_real) CYCLE ! this is to check another scheme in case of AC
         sigcme_tmp(:,ib,jb,spin)=sym_cme(:,ib,jb)
       end do
     end do
     ABI_DEALLOCATE(sym_cme)
   end if
 end do
 !
 ! Reconstruct the full sigma matrix from the upper triangle (only for HF, SEX and COHSEX)
 !if (Sigp%gwcalctyp>=20 .and. sigma_is_herm(Sigp) ) then
 !  ABI_CHECK(nspinor==1,"cannot hermitianize non-collinear sigma!")
 !  do isppol=1,Wfd%nsppol
 !    do io=1,nomega_sigc
 !      call hermitianize(sigcme_tmp(io,:,:,isppol),"Upper")
 !    end do
 !  end do
 !end if

 ! GW with contour deformation: check on the number of poles not included.
 if (ANY( mod10 == (/SIG_GW_CD,SIG_QPGW_CD/)) ) then
   call xsum_mpi(npoles_missing,comm,ierr)
   npls = SUM(npoles_missing)
   if (npls>0) then
     write(msg,'(a,i0)')" Total number of missing poles for contour deformation method: ",npls
     MSG_WARNING(msg)
     do band=minbnd,maxbnd
       npls = npoles_missing(band)
       if (npls>0) then
         write(msg,'(a,2(i0,a))')" For band ",band," there are ",npls," missing poles"
         call wrtout(std_out,msg,"COLL")
       end if
     end do
   end if
   ! Print data on the maximum value needed for the screening along the real axis
   w_localmax = MAXVAL(w_maxval)
   call xmax_mpi(w_localmax,w_max,comm,ierr)
   write(msg,'(a,F12.5,a)') ' Max omega value used in W(omega): ',w_max*Ha_eV,' eV'
   call wrtout(std_out,msg,"COLL")
   call flush_unit(std_out)
 end if 
 call timab(442,2,tsec) ! final ops
 !
 ! ===========================
 ! ==== Deallocate memory ====
 ! ===========================
 if (Psps%usepaw==1) then
   if (allocated(gw_gfft))   then
     ABI_DEALLOCATE(gw_gfft)
   end if
   call cprj_free(Cprj_ksum)
   ABI_DEALLOCATE(Cprj_ksum)
   if (allocated(Pwij_fft)) then
     call destroy_paw_pwij(Pwij_fft)
     ABI_DEALLOCATE(Pwij_fft)
   end if
   if (Dtset%pawcross==1) then
     ABI_DEALLOCATE(ur_ae_sum)
     ABI_DEALLOCATE(ur_ae_onsite_sum)
     ABI_DEALLOCATE(ur_ps_onsite_sum)
     ABI_DEALLOCATE(ktabrf)
   end if
 end if

 ABI_DEALLOCATE(wfr_sum)
 ABI_DEALLOCATE(ktabr)
 ABI_DEALLOCATE(sigc_ket)
 ABI_DEALLOCATE(ket1)
 ABI_DEALLOCATE(ket2)
 ABI_DEALLOCATE(rhotwg_ki)
 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(rhotwgp)
 ABI_DEALLOCATE(vc_sqrt_qbz)
 ABI_DEALLOCATE(omegame0i)
 ABI_DEALLOCATE(sigctmp)
 ABI_DEALLOCATE(sigc)

 if (allocated(epsm1_qbz          ))   then
   ABI_DEALLOCATE(epsm1_qbz)
 end if
 if (allocated(epsm1_trcc_qbz     ))   then
   ABI_DEALLOCATE(epsm1_trcc_qbz)
 end if
 if (allocated(epsm1_pole_qbz     ))   then
   ABI_DEALLOCATE(epsm1_pole_qbz)
 end if
! if (allocated(epsm1_pole_trcc_qbz))  deallocate(epsm1_pole_trcc_qbz)
 if (allocated(degtab             ))   then
   ABI_DEALLOCATE(degtab)
 end if
 if (allocated(ac_epsm1cqwz2      ))   then
   ABI_DEALLOCATE(ac_epsm1cqwz2)
 end if
 if (allocated(ac_integr          ))   then
   ABI_DEALLOCATE(ac_integr)
 end if
 if (allocated(aherm_sigc_ket     ))   then
   ABI_DEALLOCATE(aherm_sigc_ket)
 end if
 if (allocated(herm_sigc_ket      ))   then
   ABI_DEALLOCATE(herm_sigc_ket)
 end if

 if (Sigp%gwcomp==1) then
   ABI_DEALLOCATE(wf1swf2_g)
 endif
 if (Sigp%gwcomp==1.or.dtset%gw_eet/=-1) then
   ABI_DEALLOCATE(extrapolar_distrb)
 end if

 ABI_DEALLOCATE(proc_distrb)

 call timab(431,2,tsec)

 DBG_EXIT("COLL")

end subroutine calc_sigc_me
!!***
