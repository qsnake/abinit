!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_sigx_me
!! NAME
!! calc_sigx_me
!!
!! FUNCTION
!! Calculate diagonal and off-diagonal matrix elements of the exchange part of the self-energy operator.
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
!! Gsph_x<Gvectors_type>= Info on the G-sphere used for Sigma_x
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
!! gwx_ngfft(18)=Information about 3D FFT for the oscillator strengths, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! gwx_nfftot=number of points of the FFT grid for GW wavefunctions
!! Vcp <vcoul_t datatype> containing information on the cutoff technique
!!    %vc_sqrt(npwx,nqibz)= square-root of the coulombian potential for q-points in the IBZ
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
!! QP_BSt<Bandstructure_type>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!  Paw_pwff<Paw_pwff_type>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  QP_sym(%nsppol)<bands_symmetries>=Datatype collecting data on the irreducible representaions of the
!!    little group of kcalc in the KS representation as well as the symmetry of the bdgw_k states.
!! prtvol=Flags governing verbosity level.
!!
!! OUTPUT
!!  sigxme_tmp(minbnd:maxbnd,minbnd:maxbnd,%nsppol*Sigp%nsig_ab)=Matrix elements of Sigma_x.
!!
!! NOTES
!!  1) The treatment of the divergence of Gygi+Baldereschi (PRB 1986) is included.
!!
!!  2) On the symmetrization of Sigma matrix elements
!!     If  Sk = k+G0 then  M_G(k, Sq)= e^{-i( Sq+G).t} M_{ S^-1(G}   (k,q)
!!     If -Sk = k+G0 then  M_G(k,-Sq)= e^{-i(-Sq+G).t} M_{-S^-1(G)}^*(k,q)
!!
!! Notice the absence of G0 in the expression. Moreover, when we sum over the little group, it turns out
!! that there is a cancellation of the phase factor associated to the non-symmorphic operations due to a
!! similar term coming from the symmetrization of \epsilon^{-1}. Mind however that the nonsymmorphic phase
!! has to be considered when epsilon^-1 is reconstructed starting from the q-points in the IBZ.
!!
!!  3) the unitary transformation relating wavefunctions
!!     at symmetric k-points should be taken into account during the symmetrization
!!     of the oscillator matrix elements. In case of G_oW_o and GW_o calculations, however,
!!     it is possible to make an invariant by just including all the degenerate states and
!!     averaging the final results over the degenerate subset. Here we divide the states
!!     where the QP energies are required into complexes. Note however that this approach is not
!!     based on group theory, and it might lead to spurious results in case of accidental degeneracies.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,destroy_paw_pwij,findqg0,get_bz_item
!!      gsph_fft_tabs,hermitianize,init_paw_pwij,initmpi_seq,paw_cross_rho_tw_g
!!      paw_rho_tw_g,paw_symcprj,pawmknhat_psipsi,print_little_group,rho_tw_g
!!      rotate_fft_mesh,sigma_distribution,symmetrize_me,timab,wfd_change_ngfft
!!      wfd_get_cprj,wfd_get_ur,wfd_paw_get_aeur,wrtout,xbarrier_mpi,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_sigx_me(ikcalc,minbnd,maxbnd,Cryst,QP_BSt,Sigp,Gsph_x,Vcp,Kmesh,Qmesh,&
& Ltg_k,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,gwx_ngfft,ngfftf,&
& prtvol,pawcross,sigxme_tmp)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_gwdefs!,        only : czero_gw, cone_gw, j_gw, sigma_parameters
 use m_xmpi
 use m_defs_ptgroups
#ifdef HAVE_CLIB
 use m_clib
#endif
 use m_errors

 use m_blas,          only : xdotc, xgemv
 use m_numeric_tools, only : hermitianize
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_structure
 use m_fft_mesh,      only : rotate_FFT_mesh, cigfft 
 use m_bz_mesh,       only : bz_mesh_type, get_BZ_item, findqg0, little_group, print_little_group
 use m_gsphere,       only : gvectors_type, gsph_fft_tabs
 use m_vcoul,         only : vcoul_t
 use m_paw_pwij,      only : paw_pwff_type, paw_pwij_type, init_paw_pwij, destroy_paw_pwij, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_paw_toolbox,   only : pawfgrtab_init, pawfgrtab_free, pawfgrtab_print,paw_pwaves_lmn_t
 use m_wfs,           only : wfs_descriptor, wfd_get_ur, wfd_get_cprj, wfd_change_ngfft, wfd_paw_get_aeur
 use m_oscillators,   only : rho_tw_g
 use m_bands_sym,     only : bands_symmetries, symmetrize_me, bsym_failed
 use m_ptgroups,      only : sum_irreps

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sigx_me'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_66_paw
 use interfaces_70_gw, except_this_one => calc_sigx_me
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,prtvol,minbnd,maxbnd,pawcross
 type(Crystal_structure),intent(in) :: Cryst
 type(Bandstructure_type),intent(in) :: QP_BSt
 type(BZ_mesh_type),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(Gvectors_type),intent(in) :: Gsph_x
 type(Little_group),intent(in) :: Ltg_k
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(pawang_type),intent(in) :: Pawang
 type(wfs_descriptor),intent(inout) :: Wfd,Wfdf
!arrays
 integer,intent(in) :: gwx_ngfft(18),ngfftf(18)
 !complex(dpc),pointer :: sigxme_tmp(:,:,:)
 complex(dpc),intent(out) :: sigxme_tmp(minbnd:maxbnd,minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(Paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(bands_symmetries),intent(in) :: QP_sym(Wfd%nsppol)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=2
 integer,parameter :: use_pawnhat=0,ider0=0 
 integer :: izero,iab,ib_sum,ib,ib1,ib2,ierr,ig,ig_rot,ii,iik,itim_q,i2
 integer :: ik_bz,ik_ibz,isym_q,iq_bz,iq_ibz,spin,istat,isym,jb,is_idx,iatom
 integer :: jik,jk_bz,jk_ibz,kb,nspinor,nsppol,ifft
 integer :: nq_summed,ibsp,dimcprj_gw
 integer :: spad,spadx1,spadx2,irow
 integer :: comm,ndegs,wtqm,wtqp,my_nbks
 integer :: isym_kgw,isym_ki,gwx_mgfft,use_padfft,use_padfftf,gwx_fftalga,gwx_nfftot,nfftf,mgfftf
 integer :: nhat12_grdim
 real(dp) :: fact_sp,theta_mu_minus_esum,tol_empty,norm
 complex(dpc) :: ctmp,scprod,ph_mkgwt,ph_mkt
 complex(gwpc) :: gwpc_sigxme
 logical :: iscompatibleFFT,q_is_gamma
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: g0(3),spinor_padx(2,4) 
 integer,pointer :: igfftxg0(:),igfftfxg0(:)
 integer,allocatable :: degtab(:,:,:),l_size_atm(:)
 integer,allocatable :: gwx_gfft(:,:),gwx_gbound(:,:),gboundf(:,:)
 integer,allocatable ::  ktabr(:,:),irottb(:,:),ktabrf(:,:)
 integer,allocatable :: proc_distrb(:,:,:)
 real(dp) :: ksum(3),kgw(3),kgw_m_ksum(3),qbz(3),q0(3),tsec(2)
 real(dp) :: spinrot_kbz(4),spinrot_kgw(4)
 real(dp),pointer :: qp_ene(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: nhat12(:,:,:),grnhat12(:,:,:,:)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:),rhotwg(:),rhotwgp(:)
 complex(gwpc),allocatable :: rhotwg_ki(:,:)
 complex(gwpc),allocatable :: wfr_bdgw(:,:),wfr_sum(:)
 complex(gwpc),allocatable :: ur_ae_sum(:),ur_ae_onsite_sum(:),ur_ps_onsite_sum(:)
 complex(gwpc),allocatable :: ur_ae_bdgw(:,:),ur_ae_onsite_bdgw(:,:),ur_ps_onsite_bdgw(:,:)
 complex(gwpc),pointer :: cg_jb(:),cg_sum(:)
 complex(dpc) :: ovlp(2)
 complex(dpc),allocatable :: sym_sigx(:,:),sigx(:,:,:,:)
 logical :: can_symmetrize(Wfd%nsppol)
 logical,allocatable :: bks_mask(:,:,:)
 type(sigijtab_t),pointer :: Sigxij_tab(:)
 type(Cprj_type),allocatable :: Cprj_kgw(:,:),Cprj_ksum(:,:)
 type(Paw_pwij_type),allocatable :: Pwij_qg(:),Pwij_fft(:)

!************************************************************************

 DBG_ENTER("COLL")
 !
 ! === Initial check ===
 ABI_CHECK(Sigp%npwx==Gsph_x%ng,"")

 call timab(430,1,tsec) ! csigme (SigX)
 !
 ! === Initialize MPI variables === 
 comm = Wfd%comm
 call initmpi_seq(MPI_enreg_seq) ! * Fake MPI_type for sequential part.
 !
 ! === Initialize some values ===
 nspinor = Wfd%nspinor
 nsppol  = Wfd%nsppol
 spinor_padx(:,:)=RESHAPE((/0,0,Sigp%npwx,Sigp%npwx,0,Sigp%npwx,Sigp%npwx,0/),(/2,4/))
                                                                                      
 qp_ene => QP_BSt%eig(:,:,:)
 qp_occ => QP_BSt%occ(:,:,:)

 ib1=minbnd
 ib2=maxbnd
 !
 ! === Index of the GW point in the BZ array, its image in IBZ and time-reversal ===
 jk_bz=Sigp%kptgw2bz(ikcalc)
 call get_BZ_item(Kmesh,jk_bz,kgw,jk_ibz,isym_kgw,jik,ph_mkgwt)
 !$call get_IBZ_item(Kmesh,jk_ibz,kibz,wtk)
 spinrot_kgw(:)=Cryst%spinrot(:,isym_kgw)
 !
 write(msg,'(2a,3f8.3,2a,2(i3,a))')ch10,&
&  ' Calculating <nk|Sigma_x|nk> at k= ',kgw,ch10,&
&  ' bands n = from ',ib1,' to ',ib2,ch10
 call wrtout(std_out,msg,'COLL')

 if (ANY(gwx_ngfft(1:3) /= Wfd%ngfft(1:3)) ) call wfd_change_ngfft(Wfd,Cryst,Psps,gwx_ngfft) 
 gwx_mgfft   = MAXVAL(gwx_ngfft(1:3)) 
 gwx_fftalga = gwx_ngfft(7)/100 !; gw_fftalgc=MOD(gwx_ngfft(7),10)

 if (pawcross==1) then
   mgfftf = MAXVAL(ngfftf(1:3)) 
 end if

 can_symmetrize = .FALSE.
 if (Sigp%symsigma>0) then
   can_symmetrize = .TRUE.
   if (Sigp%gwcalctyp >= 20) then
    do spin=1,Wfd%nsppol
      can_symmetrize(spin) = .not.bsym_failed(QP_sym(spin))
      if (.not.can_symmetrize(spin)) then
        write(msg,'(a,i0,4a)')&
&         " Symmetrization cannot be performed for spin: ",spin,ch10,&
&         " band classification encountered the following problem: ",ch10,TRIM(QP_sym(spin)%err_msg)
        MSG_WARNING(msg)
      end if
    end do
   end if
   ABI_CHECK(nspinor==1,'Symmetrization with nspinor=2 not implemented')
 end if

 ABI_ALLOCATE(rhotwg_ki,(Sigp%npwx*nspinor,minbnd:maxbnd))
 rhotwg_ki=czero_gw
 ABI_ALLOCATE(rhotwg   ,(Sigp%npwx*nspinor))
 ABI_ALLOCATE(rhotwgp  ,(Sigp%npwx*nspinor))
 ABI_ALLOCATE(vc_sqrt_qbz,(Sigp%npwx))
 !
 ! === Normalization of theta_mu_minus_esum ===
 ! * If nsppol==2, qp_occ $\in [0,1]$
 SELECT CASE (nsppol)
 CASE (1)
   fact_sp=half; tol_empty=0.01   ! below this value the state is assumed empty
   if (Sigp%nspinor==2) then
    fact_sp=one; tol_empty=0.005  ! below this value the state is assumed empty
   end if
 CASE (2)
   fact_sp=one; tol_empty=0.005 ! to be consistent and obtain similar results if a metallic
 CASE DEFAULT                    ! spin unpolarized system is treated using nsppol==2
   MSG_BUG('Wrong nsppol')
 END SELECT

 ! Table for \Sigmax_ij matrix elements.
 Sigxij_tab => Sigp%Sigxij_tab(ikcalc,1:nsppol)

 ! Remove empty states from the list of states that will be distributed.
 ABI_ALLOCATE(bks_mask,(Wfd%mband,Kmesh%nbz,nsppol))
 bks_mask=.FALSE.
 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz
     ik_ibz = Kmesh%tab(ik_bz)
     do ib_sum=1,Sigp%nbnds
       bks_mask(ib_sum,ik_bz,spin) = (qp_occ(ib_sum,ik_ibz,spin)>=tol_empty) 
     end do
   end do
 end do

 ABI_ALLOCATE(proc_distrb,(Wfd%mband,Kmesh%nbz,nsppol))
 call sigma_distribution(Wfd,Kmesh,Ltg_k,Qmesh,nsppol,can_symmetrize,kgw,Sigp%mg0,my_nbks,proc_distrb,bks_mask=bks_mask)

 ABI_DEALLOCATE(bks_mask)

 write(msg,'(2(a,i4),a)')" Will sum ",my_nbks ," (b,k,s) occupied states in Sigma_x."
 call wrtout(std_out,msg,'PERS')
 !
 ! The index of G-G0 in the FFT mesh the oscillators ===
 ! * Sigp%mG0 gives the MAX G0 component to account for umklapp.
 ! * Note the size MAX(Sigp%npwx,Sigp%npwc).
 ABI_ALLOCATE(igfftxg0,(Gsph_x%ng))
 !
 ! === Precalculate the FFT index of $ R^{-1}(r-\tau) $ ===
 ! * S=\transpose R^{-1} and k_BZ = S k_IBZ
 ! * irottb is the FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.
 gwx_nfftot = PRODUCT(gwx_ngfft(1:3))
 ABI_ALLOCATE(irottb,(gwx_nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,gwx_ngfft,irottb,iscompatibleFFT)
 if (.not.iscompatibleFFT) then
   msg = "FFT mesh is not compatible with symmetries. Results might be affected by large errors!"
   MSG_WARNING(msg)
 end if

 ABI_ALLOCATE(ktabr,(gwx_nfftot,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym=Kmesh%tabo(ik_bz)
   do ifft=1,gwx_nfftot
     ktabr(ifft,ik_bz)=irottb(ifft,isym)
   end do
 end do
 ABI_DEALLOCATE(irottb)
 if (Psps%usepaw==1 .and. pawcross==1) then
   nfftf = PRODUCT(ngfftf(1:3))
   ABI_ALLOCATE(irottb,(nfftf,Cryst%nsym))
   call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,ngfftf,irottb,iscompatibleFFT)

   ABI_ALLOCATE(ktabrf,(nfftf,Kmesh%nbz))
   do ik_bz=1,Kmesh%nbz
     isym=Kmesh%tabo(ik_bz)
     do ifft=1,nfftf
       ktabrf(ifft,ik_bz)=irottb(ifft,isym)
     end do
   end do
   ABI_DEALLOCATE(irottb)
 end if
 !
 ! === Additional allocations for PAW ===
 if (Psps%usepaw==1) then
   ABI_ALLOCATE(Cprj_ksum,(Cryst%natom,nspinor))
   call cprj_alloc(Cprj_ksum,0,Wfd%nlmn_atm)

   nhat12_grdim=0
   if (use_pawnhat==1) then ! Compensation charge for \phi_a^*\phi_b
     call wrtout(std_out,"Using nhat12","COLL")
     ABI_ALLOCATE(nhat12  ,(2,gwx_nfftot,nspinor**2))
     ABI_ALLOCATE(grnhat12,(2,gwx_nfftot,nspinor**2,3*nhat12_grdim))
     ABI_ALLOCATE(l_size_atm,(Cryst%natom))
     do iatom=1,Cryst%natom
       l_size_atm(iatom)=pawtab(Cryst%typat(iatom))%lcut_size ! here be careful due to pawlcutd
     end do
     ABI_DEALLOCATE(l_size_atm)
   end if 
 end if ! usepaw==1
 !
 !allocate(sigxme_tmp(ib1:ib2,ib1:ib2,nsppol*Sigp%nsig_ab))
 ABI_ALLOCATE(sigx,(2,ib1:ib2,ib1:ib2,nsppol*Sigp%nsig_ab))
 sigxme_tmp = czero; sigx = czero     
 
 nq_summed=Kmesh%nbz
 if (Sigp%symsigma>0) then
   call print_little_group(Ltg_k,std_out,prtvol,'COLL')
   nq_summed=SUM(Ltg_k%ibzq(:))
   !
   ! === Find number of complexes and number of bands in each complex ===
   ! The tolerance is a little bit arbitrary (0.001 eV)
   ! It could be reduced, in particular in case of nearly accidental degeneracies
   ABI_ALLOCATE(degtab,(ib1:ib2,ib1:ib2,nsppol))
   degtab=0
   do spin=1,nsppol
     do ib=ib1,ib2 
       do jb=ib1,ib2 
        if (ABS(qp_ene(ib,jk_ibz,spin)-qp_ene(jb,jk_ibz,spin))<0.001/Ha_ev) degtab(ib,jb,spin)=1
       end do
     end do
   end do
!   if (ANY(degtab/=0)) then ! If two states do not belong to the same complex => matrix elements of v_xc differ
!     write(msg,'(a,3f8.3,a)')' Degenerate states at k-point = ( ',kgw(:),' ).'
!     call wrtout(std_out,msg,'COLL')
!     do spin=1,nsppol
!       do ib=ib1,ib2 
!         do jb=ib+1,ib2 
!           if (degtab(ib,jb,spin)==1) then
!             write(msg,'(a,i2,a,i4,a,i4)')' (spin ',spin,')',ib,' <====> ',jb
!             call wrtout(std_out,msg,'COLL')
!             if (ABS(Sr%vxcme(ib,jk_ibz,spin)-Sr%vxcme(jb,jk_ibz,spin))>ABS(tol6*Sr%vxcme(jb,jk_ibz,spin))) then 
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

 write(msg,'(2a,i6,a)')ch10,' calc_sigx_me: calculation status ( ',nq_summed,' to be completed):'
 call wrtout(std_out,msg,'COLL')

 ABI_ALLOCATE(wfr_sum,(gwx_nfftot*nspinor))
 if (pawcross==1) then
   ABI_ALLOCATE(ur_ae_sum,(nfftf*nspinor))
   ABI_ALLOCATE(ur_ae_onsite_sum,(nfftf*nspinor))
   ABI_ALLOCATE(ur_ps_onsite_sum,(nfftf*nspinor))
 end if
 !
 ! =======================================
 ! ==== Begin loop over k_i in the BZ ====
 ! =======================================

 do spin=1,nsppol

   if (ALL(proc_distrb(:,:,spin)/=Wfd%my_rank)) CYCLE
   !
   ! * Load wavefunctions for GW corrections.
   ABI_ALLOCATE(wfr_bdgw,(gwx_nfftot*nspinor,ib1:ib2))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out of memory in wfr_bdgw")
   do jb=ib1,ib2
     call wfd_get_ur(Wfd,jb,jk_ibz,spin,wfr_sum)
     wfr_bdgw(:,jb)=wfr_sum 
   end do

   if (Wfd%usepaw==1) then ! * Load cprj for GW states, note the indexing.
     dimcprj_gw=nspinor*(ib2-ib1+1)
     ABI_ALLOCATE(Cprj_kgw,(Cryst%natom,ib1:ib1+dimcprj_gw-1))
     call cprj_alloc(Cprj_kgw,0,Wfd%nlmn_atm)
     ibsp=ib1  
     do jb=ib1,ib2
       call wfd_get_cprj(Wfd,jb,jk_ibz,spin,Cryst,Cprj_ksum,sorted=.FALSE.)
       call paw_symcprj(jk_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_ksum) 
       call cprj_copy(Cprj_ksum,Cprj_kgw(:,ibsp:ibsp+(nspinor-1)))
       ibsp=ibsp+nspinor
     end do
     if (pawcross==1) then
       ABI_ALLOCATE(ur_ae_bdgw,(nfftf*nspinor,ib1:ib2))
       istat = ABI_ALLOC_STAT
       ABI_ALLOCATE(ur_ae_onsite_bdgw,(nfftf*nspinor,ib1:ib2))
       istat = ABI_ALLOC_STAT
       ABI_ALLOCATE(ur_ps_onsite_bdgw,(nfftf*nspinor,ib1:ib2))
       istat = ABI_ALLOC_STAT
       do jb=ib1,ib2
         call wfd_paw_get_aeur(Wfdf,jb,jk_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,&
&          ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         ur_ae_bdgw(:,jb)=ur_ae_sum
         ur_ae_onsite_bdgw(:,jb)=ur_ae_onsite_sum
         ur_ps_onsite_bdgw(:,jb)=ur_ps_onsite_sum
       end do
     end if
   end if

   do ik_bz=1,Kmesh%nbz
     !
     ! === Parallelization over k-points and spin ===
     ! * For the spin there is another check in the inner loop
     if (ALL(proc_distrb(:,ik_bz,spin)/=Wfd%my_rank)) CYCLE
     !
     ! * Find the corresponding irreducible k-point
     call get_BZ_item(Kmesh,ik_bz,ksum,ik_ibz,isym_ki,iik,ph_mkt)
     spinrot_kbz(:)=Cryst%spinrot(:,isym_ki)

     ! * Identify q and G0 where q+G0=k_GW-k_i
     kgw_m_ksum=kgw-ksum
     call findqg0(iq_bz,g0,kgw_m_ksum,Qmesh%nbz,Qmesh%bz,Sigp%mG0)

     ! === Symmetrize the matrix elements ===
     ! * Sum only q"s in IBZ_k. In this case elements are weighted
     !   according to wtqp and wtqm. wtqm is for time-reversal.
     wtqp=1; wtqm=0
     if (can_symmetrize(spin)) then
       if (Ltg_k%ibzq(iq_bz)/=1) CYCLE
       wtqp=0; wtqm=0
       do isym=1,Ltg_k%nsym_sg
         wtqp=wtqp+Ltg_k%wtksym(1,isym,iq_bz)
         wtqm=wtqm+Ltg_k%wtksym(2,isym,iq_bz)
       end do
     end if

     !%call clib_progress_bar(ik_bz,Kmesh%nbz)
     write(msg,'(2(a,i4),a,i3)')' calc_sigx_me : ik_bz ',ik_bz,'/',Kmesh%nbz,' done by processor ',Wfd%my_rank
     call wrtout(std_out,msg,'PERS')
     !
     ! * Find the corresponding irreducible q-point.
     call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
     q_is_gamma = (normv(qbz,Cryst%gmet,"G") < GW_TOL_W0)
     !
     ! Tables for the FFT of the oscillators.
     !  a) FFT index of the G-G0.
     !  b) gwx_gbound table for the zero-padded FFT performed in rhotwg. 
     ABI_ALLOCATE(gwx_gbound,(2*gwx_mgfft+8,2))
     call gsph_fft_tabs(Gsph_x,g0,gwx_mgfft,gwx_ngfft,use_padfft,gwx_gbound,igfftxg0)
     if ( ANY(gwx_fftalga == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
     if (use_padfft==0) then 
       ABI_DEALLOCATE(gwx_gbound)
       ABI_ALLOCATE(gwx_gbound,(2*gwx_mgfft+8,2*use_padfft))
     end if
     if (pawcross==1) then
       ABI_ALLOCATE(gboundf,(2*mgfftf+8,2))
       ABI_ALLOCATE(igfftfxg0,(Gsph_x%ng))
       call gsph_fft_tabs(Gsph_x,g0,mgfftf,ngfftf,use_padfftf,gboundf,igfftfxg0)
       if ( ANY(gwx_fftalga == (/2,4/)) ) use_padfftf=0
       if (use_padfftf==0) then 
         ABI_DEALLOCATE(gboundf)
         ABI_ALLOCATE(gboundf,(2*mgfftf+8,2*use_padfftf))
       end if
     end if
     !
     ! === Evaluate oscillator matrix elements ===
     ! * $ <phj/r|e^{-i(q+G)}|phi/r> - <tphj/r|e^{-i(q+G)}|tphi/r> $ in packed form.
     if (Psps%usepaw==1.and.use_pawnhat==0) then
       q0 = qbz !;if (q_is_gamma) q0 = (/0.00001_dp,0.00001_dp,0.00001_dp/) ! GW_Q0_DEFAULT
       ABI_ALLOCATE(Pwij_qg,(Psps%ntypat))
       call init_paw_pwij(Pwij_qg,Sigp%npwx,q0,Gsph_x%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
     end if
     !
     ! === Get Fourier components of the Coulombian interaction in the BZ ===
     ! * In 3D systems, neglecting umklapp,  vc(Sq,sG)=vc(q,G)=4pi/|q+G|
     ! * The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
     do ig=1,Sigp%npwx
       ig_rot = Gsph_x%rottb(ig,itim_q,isym_q)
       vc_sqrt_qbz(ig_rot)=Vcp%vc_sqrt(ig,iq_ibz)
     end do
     !
     ! === Sum over bands ===
     do ib_sum=1,Sigp%nbnds
       !
       ! === Parallelism over spin ===
       ! * This processor has this k-point but what about spin?
       if (proc_distrb(ib_sum,ik_bz,spin)/=Wfd%my_rank) CYCLE
       !
       ! * Skip empty states.
       if (qp_occ(ib_sum,ik_ibz,spin)<tol_empty) CYCLE

       call wfd_get_ur(Wfd,ib_sum,ik_ibz,spin,wfr_sum)

       if (Psps%usepaw==1) then ! Load cprj for point ksum, this spin or spinor and *THIS* band.
         ! TODO MG I could avoid doing this but I have to exchange spin and bands ???
         ! For sure there is a better way to do this!
         call wfd_get_cprj(Wfd,ib_sum,ik_ibz,spin,Cryst,Cprj_ksum,sorted=.FALSE.)
         call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_ksum) 
         if (pawcross==1) then
           call wfd_paw_get_aeur(Wfdf,ib_sum,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,&
&              ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         end if
       end if

       do jb=ib1,ib2 ! Get all <k-q,ib_sum,s|e^{-i(q+G).r}|s,jb,k>

         if (Psps%usepaw==1.and.use_pawnhat==1) then 
           i2=jb; if (nspinor==2) i2=(2*jb-1)
           spad=(nspinor-1)

           izero=0
           call pawmknhat_psipsi(Cprj_ksum,Cprj_kgw(:,i2:i2+spad),ider0,izero,MPi_enreg_seq,Cryst%natom,gwx_nfftot,gwx_ngfft,&
&          nhat12_grdim,nspinor,Cryst%ntypat,Cryst%typat,Wfd%paral_kgb,Pawang,Pawfgrtab,grnhat12,nhat12,pawtab)

#if 1
           msg = "reinstate optional Argument in rho_tw_g but mind inca slave!"
           MSG_ERROR(msg)
#else
           call rho_tw_g(Wfd%paral_kgb,nspinor,Sigp%npwx,gwx_nfftot,gwx_ngfft,1,use_padfft,igfftxg0,gwx_gbound,&
&            wfr_sum       ,iik,ktabr(:,ik_bz),ph_mkt  ,spinrot_kbz,&
&            wfr_bdgw(:,jb),jik,ktabr(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&            nspinor,rhotwg_ki(:,jb),tim_fourdp,Wfd%MPI_enreg,nhat12=nhat12)
#endif

         else
           call rho_tw_g(Wfd%paral_kgb,nspinor,Sigp%npwx,gwx_nfftot,gwx_ngfft,1,use_padfft,igfftxg0,gwx_gbound,&
&            wfr_sum       ,iik,ktabr(:,ik_bz),ph_mkt  ,spinrot_kbz,&
&            wfr_bdgw(:,jb),jik,ktabr(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&            nspinor,rhotwg_ki(:,jb),tim_fourdp,Wfd%MPI_enreg)

           if (Psps%usepaw==1.and.use_pawnhat==0) then ! Add on-site contribution, projectors are already in BZ.
             i2=jb; if (nspinor==2) i2=(2*jb-1)
             spad=(nspinor-1)
             call paw_rho_tw_g(Sigp%npwx,nspinor,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_x%gvec,&
&              Cprj_ksum(:,:),Cprj_kgw(:,i2:i2+spad),Pwij_qg,rhotwg_ki(:,jb))
           end if
           if (Psps%usepaw==1.and.pawcross==1) then ! Add paw cross term
             call paw_cross_rho_tw_g(Wfdf%paral_kgb,nspinor,Sigp%npwx,nfftf,ngfftf,1,use_padfftf,igfftfxg0,gboundf,&
&             ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum,iik,ktabrf(:,ik_bz),ph_mkt,spinrot_kbz,&
&             ur_ae_bdgw(:,jb),ur_ae_onsite_bdgw(:,jb),ur_ps_onsite_bdgw(:,jb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&             nspinor,rhotwg_ki(:,jb),tim_fourdp,Wfdf%MPI_enreg)
           end if
         end if
         !
         ! === Multiply by the square root of the Coulomb term ===
         ! * In 3-D systems, the factor sqrt(4pi) is included)
         do ii=1,nspinor
           spad=(ii-1)*Sigp%npwx
           rhotwg_ki(spad+1:spad+Sigp%npwx,jb)=rhotwg_ki(spad+1:spad+Sigp%npwx,jb)*vc_sqrt_qbz(1:Sigp%npwx)
         end do
         !
         ! === Treat analytically the case q --> 0 ===
         ! * The oscillator is evaluated at q=O as it is considered constant in the small cube around Gamma
         !   while the Colulomb term is integrated out.
         ! * In the scalar case we have nonzero contribution only if ib_sum==jb
         ! * For nspinor==2 evalute <ib_sum,up|jb,up> and <ib_sum,dwn|jb,dwn>,
         !   impose orthonormalization since npwwfn might be < npwvec.
         if (ik_bz==jk_bz) then
           if (nspinor==1) then
             rhotwg_ki(1,jb)=czero_gw
             if (ib_sum==jb) rhotwg_ki(1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)
           else
             ! TODO Recheck this!
             cg_sum  => Wfd%Wave(ib,ik_ibz,spin)%ug
             cg_jb   => Wfd%Wave(jb,jk_ibz,spin)%ug

             ctmp = xdotc(Wfd%npwwfn*Wfd%nspinor,cg_sum,1,cg_jb,1) 
             ovlp(1) = REAL(ctmp) 
             ovlp(2) = AIMAG(ctmp) 

             if (Psps%usepaw==1) then
               i2=(2*jb-1)
               ovlp = ovlp + paw_overlap(Cprj_ksum,Cprj_kgw(:,i2:i2+1),Cryst%typat,Pawtab)
             end if
             !ovlp(2) = -ovlp(1)
             !if (ib_sum==jb) ovlp(2)=cone_gw-ovlp(1)
             if (ib_sum==jb) then
               norm=DBLE(ovlp(1)+ovlp(2))
               ovlp(1)=DBLE(ovlp(1)/norm)
               ovlp(2)=DBLE(ovlp(2)/norm)
             else
               scprod=ovlp(1)+ovlp(2)
               ovlp(1)=ovlp(1)-scprod*half
               ovlp(2)=ovlp(2)-scprod*half
             end if
             rhotwg_ki(1          ,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)*ovlp(1)
             rhotwg_ki(Sigp%npwx+1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)*ovlp(2)
           end if
         end if
       end do !jb Got all matrix elements from minbnd up to maxbnd.

       theta_mu_minus_esum=fact_sp*qp_occ(ib_sum,ik_ibz,spin)

       do kb=ib1,ib2
         !
         ! * The ket \Sigma|\phi_{k,kb}>.
         rhotwgp(:)=rhotwg_ki(:,kb)
         !
         ! Loop over the non-zero row elements of this column.
         ! 1) If gwcalctyp<20 : only diagonal elements since QP==KS.
         ! 2) If gwcalctyp>=20:
         !     * Only off-diagonal elements connecting states with same character.
         !     * Only the upper triangle if HF, SEX, or COHSEX. 
         do irow=1,Sigxij_tab(spin)%col(kb)%size1

           jb = Sigxij_tab(spin)%col(kb)%bidx(irow)
           rhotwg=rhotwg_ki(:,jb)

           ! === Calculate bare exchange <\phi_j|\Sigma_x|\phi_k> ===
           ! * Do the scalar product only if ib_sum is occupied.
           if (theta_mu_minus_esum/fact_sp >= tol_empty) then
             do iab=1,Sigp%nsig_ab
               spadx1=spinor_padx(1,iab)
               spadx2=spinor_padx(2,iab)

               gwpc_sigxme=-XDOTC(Sigp%npwx,rhotwg(spadx1+1:),1,rhotwgp(spadx2+1:),1)*theta_mu_minus_esum
               !
               ! === Accumulate and symmetrize Sigma_x ===
               ! * -wtqm comes from time-reversal (exchange of band indeces)
               is_idx=spin; if (nspinor==2) is_idx=iab

               sigxme_tmp(jb,kb,is_idx) = sigxme_tmp(jb,kb,is_idx) + & 
&                (wtqp+wtqm)*DBLE(gwpc_sigxme) + (wtqp-wtqm)*j_gw*AIMAG(gwpc_sigxme)

               sigx(1,jb,kb,is_idx) = sigx(1,jb,kb,is_idx) + wtqp *      gwpc_sigxme
               sigx(2,jb,kb,is_idx) = sigx(2,jb,kb,is_idx) + wtqm *CONJG(gwpc_sigxme)

             end do
           end if
         end do ! jb used to calculate matrix elements of $\Sigma_x$

       end do !kb to calculate matrix elements of $\Sigma_x$
     end do !ib_sum
     !
     ! Deallocate k-dependent quantities.
     ABI_DEALLOCATE(gwx_gbound)
     istat = ABI_ALLOC_STAT
     if (pawcross==1) then
       ABI_DEALLOCATE(gboundf)
       istat = ABI_ALLOC_STAT
     end if

     if (Psps%usepaw==1.and.use_pawnhat==0) then
       call destroy_paw_pwij(Pwij_qg)
       ABI_DEALLOCATE(Pwij_qg)
     end if
   end do !ik_bz Got all diagonal (off-diagonal) matrix elements.

   ABI_DEALLOCATE(wfr_bdgw)
   if (Wfd%usepaw==1) then
     call cprj_free(Cprj_kgw )
     ABI_DEALLOCATE(Cprj_kgw)
     if (pawcross==1) then
       ABI_DEALLOCATE(ur_ae_bdgw)
       ABI_DEALLOCATE(ur_ae_onsite_bdgw)
       ABI_DEALLOCATE(ur_ps_onsite_bdgw)
     end if
   end if
 end do !spin

 ABI_DEALLOCATE(igfftxg0)
 if (pawcross==1) then
   ABI_DEALLOCATE(igfftfxg0)
 end if
 !
 ! Gather contributions from all the CPUs.
 call xbarrier_mpi(comm)

 call xsum_mpi(sigxme_tmp,comm,ierr)
 call xsum_mpi(sigx,comm,ierr)
 !
 ! === Multiply by constants ===
 ! * For 3D systems sqrt(4pi) is included in vc_sqrt_qbz ===
 sigxme_tmp = (one/(Cryst%ucvol*Kmesh%nbz)) * sigxme_tmp
 sigx       = (one/(Cryst%ucvol*Kmesh%nbz)) * sigx
 !
 ! === If we have summed over the IBZ_q now we have to average over complexes ===
 ! * Presently only diagonal terms are considered
 ! * TODO QP-SCGW required a more involved approach, there is a check in sigma
 ! * TODO it does not work if spinor==2.

 do spin=1,nsppol
   if (can_symmetrize(spin)) then
     ABI_ALLOCATE(sym_sigx,(ib1:ib2,ib1:ib2))
     sym_sigx=czero
     !
     ! === Average over degenerate diagonal elements ===
     ! NOTE: frequencies for \Sigma_c(\omega) should be equal to avoid spurious results.
     ! another good reason to use a strict criterion for the tollerance on eigenvalues.
     do ib=ib1,ib2 
       ndegs=0
       do jb=ib1,ib2 
         if (degtab(ib,jb,spin)==1) then
           sym_sigx(ib,ib)=sym_sigx(ib,ib) +SUM(sigx(:,jb,jb,spin))
         end if
         ndegs=ndegs+degtab(ib,jb,spin)
       end do
       sym_sigx(ib,ib)=sym_sigx(ib,ib)/ndegs
     end do

     if (Sigp%gwcalctyp >= 20) then 
       call symmetrize_me(QP_sym(spin),ib1,ib2,sigx(:,:,:,spin),sym_sigx)
     end if
     !
     ! ==== Copy symmetrized values ====
     do ib=ib1,ib2
       do jb=ib1,ib2 
         sigxme_tmp(ib,jb,spin)=sym_sigx(ib,jb)
       end do
     end do

     ABI_DEALLOCATE(sym_sigx)
   end if
 end do
 !
 if (Sigp%gwcalctyp>=20) then ! Reconstruct the full sigma_x matrix from the upper triangle.
   ABI_CHECK(nspinor==1,"cannot hermitianize non-collinear sigma!")
   do spin=1,nsppol
     call hermitianize(sigxme_tmp(:,:,spin),"Upper")
   end do
 end if
 !
 ! ===========================
 ! ==== Deallocate memory ====
 ! ===========================
 if (Psps%usepaw==1) then
   if (allocated(gwx_gfft))  then
     ABI_DEALLOCATE(gwx_gfft)
   end if
   call cprj_free(Cprj_ksum)
   ABI_DEALLOCATE(Cprj_ksum)
   if (allocated(Pwij_fft)) then
     call destroy_paw_pwij(Pwij_fft)
     ABI_DEALLOCATE(Pwij_fft)
   end if
   if (use_pawnhat==1) then
     ABI_DEALLOCATE(nhat12)
     ABI_DEALLOCATE(grnhat12)
     istat = ABI_ALLOC_STAT
   end if
   if (pawcross==1) then
     ABI_DEALLOCATE(ur_ae_sum)
     ABI_DEALLOCATE(ur_ae_onsite_sum)
     ABI_DEALLOCATE(ur_ps_onsite_sum)
     ABI_DEALLOCATE(ktabrf)
   end if
 end if

 ABI_DEALLOCATE(wfr_sum)
 ABI_DEALLOCATE(rhotwg_ki)
 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(rhotwgp)
 ABI_DEALLOCATE(vc_sqrt_qbz)
 ABI_DEALLOCATE(ktabr)
 ABI_DEALLOCATE(sigx)
 ABI_DEALLOCATE(proc_distrb)

 if (allocated(degtab))   then
   ABI_DEALLOCATE(degtab)
 end if

 call timab(430,2,tsec) ! csigme (SigX)

 DBG_EXIT("COLL")

end subroutine calc_sigx_me
!!***
