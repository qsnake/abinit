!{\src2tex{textfont=tt}}
!!****f* ABINIT/check_completeness
!! NAME
!! check_completeness
!!
!! FUNCTION
!! Checks the completeness of the wavefunctions basis 
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf, GKA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! use_tr=If .TRUE. valence states are stored in Wfs_val and only resonant transitions are calculated
!!  (time reversal is assumed)
!! Dtset <type(dataset_type)>=all input variables in this dataset
!! Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!!    %natom=number of atoms
!!    %nsym=number of symmetries
!!    %xred(3,natom)=reduced coordinated of atoms
!!    %typat(natom)=type of each atom
!!    %rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!    %timrev= 2 if time reversal can be used, 1 otherwise
!! qpoint(3)=reciprocal space coordinates of the q wavevector
!! Ep<type(epsilonm1_parameters_type)>= Parameters related to the calculation of the inverse dielectric matrix.
!!    %nbnds=number of bands summed over
!!    %npwe=number of planewaves for the irreducible polarizability X^0_GGp
!!    %npwvec=maximum number of G vectors 
!!     used to define the dimension of some arrays e.g igfft
!!    %nsppol=1 for unpolarized, 2 for spin-polarized
!!    %spmeth=1 if we use the spectral method, 0 for standard Adler-Wiser expression
!!    %spsmear=gaussian broadening used to approximate the delta distribution
!!    %zcut=small imaginary shift to avoid poles in X^0
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Kmesh <bz_mesh_type>= datatype gathering parameters related to the k-point sampling
!!    %nibz=number of k-points in the IBZ
!!    %nbz=number of k-points in the BZ
!!    %bz(3,nbz)=reduced coordinates for k-points in the full Brillouin zone
!!    %ibz(3,nibz)=reduced coordinates for k-points in the irreducible wedge
!!    %tab(nbz)=mapping between a kpt in the BZ (array bz) and the irred point in the array ibz
!!    %tabi(nbz)= -1 if inversion is needed to obtain this particular kpt in the BZ, 1 means identity
!!    %tabo(nbz)= for each point in the BZ, the index of the symmetry operation S in reciprocal
!!      space which rotates k_IBZ onto \pm k_BZ (depending on tabi)
!!    %tabp(nbz)= For each k_BZ, it gives the phase factors associated to non-symmorphic operations, i.e
!!      e^{-i 2 \pi k_IBZ \cdot R{^-1}t} == e{-i 2\pi k_BZ cdot t} where :
!!      \transpose R{-1}=S and (S k_IBZ) = \pm k_BZ (depending on ktabi)
!!    %tabr(nfftot,nbz) For each point r on the real mesh and for each k-point in the BZ, tabr
!!      gives the index of (R^-1 (r-t)) in the FFT array where R=\transpose S^{-1} and k_BZ=S k_IBZ.
!!      t is the fractional translation associated to R
!! Gsph_epsG0<gvectors_type data type> The G-sphere used to describe deltaI/eps. (including umklapp G0 vectors)
!!    %ng=number of G vectors for deltaI
!!    %rottbm1(ng,2,nsym)=contains the index (IS^{-1}) G
!!    %phmGt(Ep%npwe,nsym)=phase factors e^{-iG \cdot t} needed to symmetrize oscillator matrix elements and epsilon
!!    %gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!!    %gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!! ngfft_gw(18)= array containing all the information for 3D FFT for the oscillator strengths (see input variable)
!! nfftot_gw=Total number of points in the GW FFT grid
!! Ltg_q<Little group>=Data type gathering information on the little group of the q-points.
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang<pawang_type> angular mesh discretization and related data:
!! QP_BSt<Bandstructure_type>=Quasiparticle energies and occupations (for the moment real quantities)
!!   %mband=MAX number of bands over k-points and spin (==Ep%nbnds)
!!   %occ(mband,nkpt,nsppol)=QP occupation numbers, for each k point in IBZ, and each band
!!   %eig(mband,nkpt,nsppol)=GW energies, for self-consistency purposes
!!  Paw_pwff<Paw_pwff_type>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  Wfd<wfs_descriptor>=Object used to access the wavefunctions
!!
!! OUTPUT
!!  deltaI(Ep%npwe,Ep%npwe)=incompleteness matrix at wavevector qpoint
!!  proportional to Sum_{k,i} <k,i | e^{i(G'-G)r} | k,i > - \sum_j M_{kij}(q+G) M^*_{kij}(q+G')
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      assemblychi0_sym,calc_wfwfg,chi0_bbp_mask,completechi0_deltapart
!!      cprj_alloc,cprj_free,destroy_gsphere,destroy_paw_pwij,get_bz_diff
!!      get_bz_item,get_gftt,gsph_fft_tabs,gsph_in_fftbox,init_paw_pwij
!!      nullify_gsphere,paw_cross_rho_tw_g,paw_rho_tw_g,paw_symcprj
!!      print_little_group,rho_tw_g,symmetrize_afm_chi0,timab,wfd_barrier
!!      wfd_change_ngfft,wfd_distribute_kb_kpbp,wfd_get_cprj,wfd_get_ur
!!      wfd_paw_get_aeur,wrtout,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine check_completeness(use_tr,Dtset,Cryst,qpoint,Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,&
& Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,ngfft_gw,nfftot_gw,ngfftf,nfftf_tot,deltaI,ktabr,ktabrf,Ltg_q,Wfd,Wfdf)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_commutator_vkbr
 use m_errors

 use m_gwdefs,        only : GW_TOL_DOCC, GW_TOL_W0, czero_gw, epsilonm1_parameters, g0g0w
 use m_numeric_tools, only : imin_loc
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_structure
 use m_bz_mesh,       only : bz_mesh_type, get_BZ_item, get_BZ_diff, little_group, print_little_group
 use m_gsphere,       only : gvectors_type, gsph_fft_tabs, nullify_gsphere, destroy_gsphere, gsph_in_fftbox
 use m_fft_mesh,      only : get_gftt
 use m_paw_pwij,      only : paw_pwff_type, paw_pwij_type, init_paw_pwij, destroy_paw_pwij, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_paw_toolbox,   only : paw_pwaves_lmn_t
 !use m_wfs,           only : wfd_get_ur, wfs_descriptor, wfd_distribute_kb_kpbp, wfd_get_cprj, wfd_barrier, wfd_change_ngfft,&
!&                            wfd_paw_get_aeur, wfd_copy, wfd_destroy
 use m_wfs
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_io_tools,      only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_completeness'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_66_paw
 use interfaces_70_gw, except_this_one => check_completeness
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_gw,nfftf_tot
 logical,intent(in) :: use_tr
 type(bandstructure_type),intent(in) :: QP_BSt
 type(bz_mesh_type),intent(in) :: Kmesh
 type(crystal_structure),intent(in) :: Cryst
 type(dataset_type),intent(in) :: Dtset
 type(epsilonm1_parameters),intent(in) :: Ep
 type(gvectors_type),intent(in) :: Gsph_epsG0
 type(little_group),intent(in) :: Ltg_q
 type(pawang_type),intent(in) :: Pawang
 type(pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wfd,Wfdf
!arrays
 integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz),ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
 integer,intent(in) :: ngfft_gw(18),ngfftf(18)
 real(dp),intent(in) :: qpoint(3)
 complex(gwpc),intent(out) :: deltaI(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,1)
 type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)


!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=1,two_poles=2,one_pole=1,fake_nomega=1
 integer :: dim_rtwg,ib1,ib2,ierr
 integer :: ig1,ig2,ik_bz,ik_ibz,ikmq_bz,ikmq_ibz
 integer :: is,istat,isym_k
 integer :: isym_kmq,itim_k,itim_kmq
 integer :: nfound,nkpt_summed,nspinor,nsppol
 integer :: comm,gw_mgfft,use_padfft,gw_fftalga,mgfftf,use_padfftf
 integer :: my_nbbp,my_nbbpks,spin
 integer :: nbmax,dummy,local_std_out
 real(dp) :: spin_fact,weight
 real(dp) :: gw_gsq
 complex(dpc) :: ph_mkmqt,ph_mkt
 complex(gwpc) :: local_czero_gw
 logical :: qzero,isirred_k,isirred_kmq
 character(len=500) :: msg,allup
 type(Gvectors_type) :: Gsph_FFT
 type(Epsilonm1_parameters) :: Epp
!arrays
 integer :: G0(3),umklp_k(3),umklp_kmq(3)
 integer :: wtk_ltg(Kmesh%nbz),got(Wfd%nproc)
 integer,allocatable :: tabr_k(:),tabr_kmq(:),tabrf_k(:),tabrf_kmq(:)
 integer,allocatable :: igfftepsG0(:),gspfft_igfft(:),igfftepsG0f(:)
 integer,allocatable :: gw_gfft(:,:),gw_gbound(:,:),dummy_gbound(:,:),gboundf(:,:)
 integer,allocatable :: bbp_ks_distrb(:,:,:,:)
 real(dp) :: kbz(3),kmq_bz(3),spinrot_k(4),spinrot_kmq(4),q0(3),tsec(2)
 real(dp),pointer :: qp_occ(:,:,:)
 complex(dpc),allocatable :: kkweight(:,:)
 complex(dpc) :: green_w(1),green_enhigh_w(1)
 complex(gwpc),allocatable :: rhotwg(:),wfr1(:)
 complex(gwpc),allocatable :: wfr2(:),wfwfg(:)
 complex(gwpc),allocatable :: ur_ae1(:),ur_ae_onsite1(:),ur_ps_onsite1(:)
 complex(gwpc),allocatable :: ur_ae2(:),ur_ae_onsite2(:),ur_ps_onsite2(:)
 logical,allocatable :: bbp_mask(:,:)
 type(cprj_type),allocatable :: Cprj1_kmq(:,:),Cprj2_k(:,:)
 type(paw_pwij_type),allocatable :: Pwij(:),Pwij_fft(:)
!************************************************************************

 DBG_ENTER("COLL")

 call timab(331,1,tsec) ! cchi0

 if ( ANY(ngfft_gw(1:3) /= Wfd%ngfft(1:3)) ) call wfd_change_ngfft(Wfd,Cryst,Psps,ngfft_gw) 
 gw_mgfft = MAXVAL(ngfft_gw(1:3)) 
 gw_fftalga = ngfft_gw(7)/100 !; gw_fftalgc=MOD(ngfft_gw(7),10)

 if (Dtset%pawcross==1) then
   mgfftf = MAXVAL(ngfftf(1:3)) 
 end if
 !
 ! === Initialize MPI variables ===
 comm = Wfd%comm
 !
 ! == Copy some values ===
 nsppol  = Wfd%nsppol
 nspinor = Wfd%nspinor
 dim_rtwg=1; if (nspinor==2) dim_rtwg=4    !can reduce size depending on Ep%nI and Ep%nj

 qp_occ    => QP_BSt%occ(:,:,:)
 Epp = Ep
 Epp%nomega=1
 Epp%gwcomp=1
 !
 ! === Initialize the completeness correction  ===
 call nullify_gsphere(Gsph_FFT)

 ! Init the largest G-sphere contained in the FFT box for the wavefunctions.
 call gsph_in_fftbox(Gsph_FFT,Cryst,Wfd%ngfft)

 !call print_gsphere(Gsph_FFT,unit=std_out,prtvol=10)

 ABI_ALLOCATE(gspfft_igfft,(Gsph_FFT%ng))
 ABI_ALLOCATE(dummy_gbound,(2*gw_mgfft+8,2))
 !
 ! Mapping between G-sphere and FFT box.
 call gsph_fft_tabs(Gsph_FFT,(/0,0,0/),Wfd%mgfft,Wfd%ngfft,dummy,dummy_gbound,gspfft_igfft)
 ABI_DEALLOCATE(dummy_gbound)
 !
 if (Psps%usepaw==1) then  ! * Prepare the onsite contributions on the GW FFT mesh.   
   ABI_ALLOCATE(gw_gfft,(3,nfftot_gw))
   q0=zero
   call get_gftt(ngfft_gw,q0,Cryst%gmet,gw_gsq,gw_gfft) ! Get the set of plane waves in the FFT Box.
   ABI_ALLOCATE(Pwij_fft,(Psps%ntypat))
   call init_paw_pwij(Pwij_fft,nfftot_gw,(/zero,zero,zero/),gw_gfft,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
 end if
 !
 !
 ! === Setup weights (2 for spin unpolarized sistem, 1 for polarized) ===
 ! * spin_fact is used to normalize the occupation factors to one. Consider also the AFM case.
 SELECT CASE (nsppol)
 CASE (1)
   weight=two/Kmesh%nbz; spin_fact=half
   if (Wfd%nspden==2) then
    weight=one/Kmesh%nbz; spin_fact=half
   end if
   if (nspinor==2) then
    weight=one/Kmesh%nbz; spin_fact=one
   end if
 CASE (2)
   weight=one/Kmesh%nbz; spin_fact=one
 CASE DEFAULT
   MSG_BUG("Wrong nsppol")
 END SELECT
 !
 ! === Weight for points in the IBZ_q ===
 wtk_ltg(:)=1
 if (Ep%symchi==1) then
   do ik_bz=1,Ltg_q%nbz
     wtk_ltg(ik_bz)=0
     if (Ltg_q%ibzq(ik_bz)/=1) CYCLE ! Only k-points in the IBZ_q.
     wtk_ltg(ik_bz)=SUM(Ltg_q%wtksym(:,:,ik_bz))
   end do
 end if

 write(msg,'(a,i2)')' Using symmetries to sum only over the IBZ_q  = ',Ep%symchi
 call wrtout(std_out,msg,'COLL')

 if (use_tr) then
   ! Special care has to be taken in metals and/or spin dependent systems
   ! as Wfs_val might contain unoccupied states.
   write(msg,'(a)')' Using faster algorithm based on time reversal symmetry. '
   call wrtout(std_out,msg,'COLL')
 end if

 ! TODO this table can be calculated for each k-point
 my_nbbpks=0; allup="All"; got=0
 ABI_ALLOCATE(bbp_ks_distrb,(Wfd%mband,Wfd%mband,Kmesh%nbz,nsppol))
 ABI_ALLOCATE(bbp_mask,(Wfd%mband,Wfd%mband))

 do spin=1,nsppol      
   do ik_bz=1,Kmesh%nbz
     if (Ep%symchi==1) then
       if (Ltg_q%ibzq(ik_bz)/=1) CYCLE  ! Only IBZ_q
     end if
     !
     ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)
     !
     ! * Get index of k-q in the BZ, stop if not found as the weight=one/nkbz is not correct.
     call get_BZ_diff(Kmesh,kbz,qpoint,ikmq_bz,g0,nfound)
     ABI_CHECK(nfound==1,"Check kmesh")
     !
     ! * Get ikmq_ibz, non-symmorphic phase, ph_mkmqt, and symmetries from ikmq_bz.
     call get_BZ_item(Kmesh,ikmq_bz,kmq_bz,ikmq_ibz,isym_kmq,itim_kmq)

     call chi0_bbp_mask(Epp,use_tr,QP_BSt,Wfd%mband,ikmq_ibz,ik_ibz,spin,spin_fact,bbp_mask)

     call wfd_distribute_kb_kpbp(Wfd,ikmq_ibz,ik_ibz,spin,allup,my_nbbp,bbp_ks_distrb(:,:,ik_bz,spin),got,bbp_mask) 
     my_nbbpks = my_nbbpks + my_nbbp
   end do
 end do

 ABI_DEALLOCATE(bbp_mask)

 write(msg,'(a,i0,a)')" Will sum ",my_nbbpks," (b,b',k,s) states in deltaI."
 call wrtout(std_out,msg,'PERS')

 got=0 
 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz
     if (Ep%symchi==1) then
       if (Ltg_q%ibzq(ik_bz)/=1) CYCLE  ! Only IBZ_q
     end if
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)
   enddo
 enddo

 if (Psps%usepaw==1) then
   ABI_ALLOCATE(Pwij,(Psps%ntypat))
   call init_paw_pwij(Pwij,Epp%npwepG0,qpoint,Gsph_epsG0%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
 end if

 call wrtout(std_out,' Calculating deltaI(q,G,G")','COLL')

 nkpt_summed=Kmesh%nbz
 if (Ep%symchi==1) then
   nkpt_summed=Ltg_q%nibz_ltg
   call print_little_group(Ltg_q,std_out,Dtset%prtvol,'COLL')
 end if

 write(msg,'(a,i6,a)')' Calculation status : ',nkpt_summed,' to be completed '
 call wrtout(std_out,msg,'COLL')

 !
 ! ============================================
 ! === Begin big fat loop over transitions ===
 ! ============================================
 deltaI=czero_gw

 ! === Loop on spin to calculate trace $\chi_{up,up}+\chi_{down,down}$ ===
 ! * Only $\chi_{up,up} for AFM.
 do is=1,nsppol

   if (ALL(bbp_ks_distrb(:,:,:,is) /= Wfd%my_rank)) CYCLE

   local_std_out = std_out
   local_czero_gw = czero_gw
   ! Allocation of arrays that are private to loop
   ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))
   green_w(1)=cone
   green_enhigh_w(1)=cone
   if (Psps%usepaw==1) then
     ABI_ALLOCATE(Cprj2_k  ,(Cryst%natom,nspinor))
     call cprj_alloc(Cprj2_k,  0,Wfd%nlmn_atm)
     ABI_ALLOCATE(Cprj1_kmq,(Cryst%natom,nspinor))
     call cprj_alloc(Cprj1_kmq,0,Wfd%nlmn_atm)
     if (Dtset%pawcross==1) then
       ABI_ALLOCATE(ur_ae1,(nfftf_tot*Wfd%nspinor))
       ABI_ALLOCATE(ur_ae_onsite1,(nfftf_tot*Wfd%nspinor))
       ABI_ALLOCATE(ur_ps_onsite1,(nfftf_tot*Wfd%nspinor))
       ABI_ALLOCATE(ur_ae2,(nfftf_tot*Wfd%nspinor))
       ABI_ALLOCATE(ur_ae_onsite2,(nfftf_tot*Wfd%nspinor))
       ABI_ALLOCATE(ur_ps_onsite2,(nfftf_tot*Wfd%nspinor))
       ABI_ALLOCATE(igfftepsG0f,(Ep%npwepG0))
       ABI_ALLOCATE(tabrf_k,(nfftf_tot))
       ABI_ALLOCATE(tabrf_kmq,(nfftf_tot))
     end if
   end if
   ABI_ALLOCATE(rhotwg,(Ep%npwepG0*nspinor**2))
   ABI_ALLOCATE(tabr_k,(nfftot_gw))
   ABI_ALLOCATE(tabr_kmq,(nfftot_gw))
   ABI_ALLOCATE(wfr1,(Wfd%nfft*nspinor))
   ABI_ALLOCATE(wfr2,(Wfd%nfft*nspinor))
   ABI_ALLOCATE(igfftepsG0,(Ep%npwepG0))

   do ik_bz=1,Kmesh%nbz ! Loop over k-points in the BZ.

     if (Ep%symchi==1) then
       if (Ltg_q%ibzq(ik_bz)/=1) CYCLE  ! Only IBZ_q
     end if

     if (ALL(bbp_ks_distrb(:,:,ik_bz,is) /= Wfd%my_rank)) CYCLE

     write(msg,'(2(a,i4),a,i2,a,i3)')' ik = ',ik_bz,' / ',Kmesh%nbz,' is = ',is,' done by processor ',Wfd%my_rank
     call wrtout(std_out,msg,'PERS')

     ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt,umklp_k,isirred_k)

     ! * Get index of k-q in the BZ, stop if not found as the weight=one/nkbz is not correct.
     call get_BZ_diff(Kmesh,kbz,qpoint,ikmq_bz,G0,nfound); if (nfound==0) call leave_new('COLL')

     ! * Get ikmq_ibz, non-symmorphic phase, ph_mkmqt, and symmetries from ikmq_bz.
     call get_BZ_item(Kmesh,ikmq_bz,kmq_bz,ikmq_ibz,isym_kmq,itim_kmq,ph_mkmqt,umklp_kmq,isirred_kmq)

     ! * Copy tables for rotated FFT points
     tabr_k(:)  =ktabr(:,ik_bz)
     spinrot_k(:)=Cryst%spinrot(:,isym_k)

     tabr_kmq(:)=ktabr(:,ikmq_bz)
     spinrot_kmq(:)=Cryst%spinrot(:,isym_kmq)

     if (Dtset%pawcross==1) then
       tabrf_k(:)  =ktabrf(:,ik_bz)
       tabrf_kmq(:)=ktabrf(:,ikmq_bz)
     end if

     !
     ! Tables for the FFT of the oscillators.
     !  a) FFT index of G-G0.
     !  b) gw_gbound table for the zero-padded FFT performed in rhotwg. 
     ABI_ALLOCATE(gw_gbound,(2*gw_mgfft+8,2))
     call gsph_fft_tabs(Gsph_epsG0,g0,gw_mgfft,ngfft_gw,use_padfft,gw_gbound,igfftepsG0)
     if ( ANY(gw_fftalga == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
     if (use_padfft==0) then 
       ABI_DEALLOCATE(gw_gbound)
       ABI_ALLOCATE(gw_gbound,(2*gw_mgfft+8,2*use_padfft))
     end if
     if (Dtset%pawcross==1) then
        ABI_ALLOCATE(gboundf,(2*mgfftf+8,2))
       call gsph_fft_tabs(Gsph_epsG0,g0,mgfftf,ngfftf,use_padfftf,gboundf,igfftepsG0f)
       if ( ANY(gw_fftalga == (/2,4/)) ) use_padfftf=0
       if (use_padfftf==0) then 
         ABI_DEALLOCATE(gboundf)
         ABI_ALLOCATE(gboundf,(2*mgfftf+8,2*use_padfftf))
       end if
     end if

     nbmax=Ep%nbnds

     do ib1=1,nbmax ! Loop over all states.

       if (ALL(bbp_ks_distrb(ib1,:,ik_bz,is) /= Wfd%my_rank)) CYCLE

       call wfd_get_ur(Wfd,ib1,ikmq_ibz,is,wfr1)

       if (Psps%usepaw==1) then 
         call wfd_get_cprj(Wfd,ib1,ikmq_ibz,is,Cryst,Cprj1_kmq,sorted=.FALSE.)
         call paw_symcprj(ikmq_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj1_kmq) 
         if (Dtset%pawcross==1) then
           call wfd_paw_get_aeur(Wfdf,ib1,ikmq_ibz,is,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae1,ur_ae_onsite1,ur_ps_onsite1)
         end if
       end if

       do ib2=1,nbmax ! Loop over "valence" states.

         if (bbp_ks_distrb(ib1,ib2,ik_bz,is) /= Wfd%my_rank) CYCLE

         if (qp_occ(ib2,ik_ibz,is) < GW_TOL_DOCC) CYCLE

         call wfd_get_ur(Wfd,ib2,ik_ibz,is,wfr2)

         if (Psps%usepaw==1) then 
           call wfd_get_cprj(Wfd,ib2,ik_ibz,is,Cryst,Cprj2_k,sorted=.FALSE.)
           call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj2_k) 
           if (Dtset%pawcross==1) then
             call wfd_paw_get_aeur(Wfdf,ib2,ik_ibz,is,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae2,ur_ae_onsite2,ur_ps_onsite2)
           end if
         end if

         if (ib1==ib2) then
          ! Identity part. Doesnt work for spinor
          call calc_wfwfg(Wfd%MPI_enreg,Wfd%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw,wfr2,wfr2,wfwfg)

          if (Psps%usepaw==1) then
            call paw_rho_tw_g(nfftot_gw,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,gw_gfft,&
&             Cprj2_k,Cprj2_k,Pwij_fft,wfwfg)

            ! Add PAW cross term
            if (Dtset%pawcross==1) then
              call paw_cross_rho_tw_g(Wfd%paral_kgb,nspinor,Ep%npwepG0,nfftf_tot,ngfftf,1,use_padfftf,igfftepsG0f,gboundf,&
&              ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_kmq,tabrf_kmq,ph_mkmqt,spinrot_kmq,&
&              ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k  ,tabrf_k  ,ph_mkt  ,spinrot_k,&
&              dim_rtwg,wfwfg,tim_fourdp,Wfd%MPI_enreg)
            end if
          end if

          qzero=.FALSE.
          call completechi0_deltapart(ik_bz,qzero,Ep%symchi,Ep%npwe,Gsph_FFT%ng,fake_nomega,nspinor,&
&           nfftot_gw,ngfft_gw,gspfft_igfft,gsph_FFT,Ltg_q,green_enhigh_w,wfwfg,deltaI)
         end if

         !
         ! ==== Form rho-twiddle(r)=u^*_{b1,kmq_bz}(r) u_{b2,kbz}(r) and its FFT transform ====
         call rho_tw_g(Wfd%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&          wfr1,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,&
&          wfr2,itim_k  ,tabr_k  ,ph_mkt  ,spinrot_k,&
&          dim_rtwg,rhotwg,tim_fourdp,Wfd%MPI_enreg) 
         
         if (Psps%usepaw==1) then! Add PAW on-site contribution, projectors are already in the BZ.
           call paw_rho_tw_g(Ep%npwepG0,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_epsG0%gvec,&
&           Cprj1_kmq,Cprj2_k,Pwij,rhotwg)

           ! Add PAW cross term
           if (Dtset%pawcross==1) then
             call paw_cross_rho_tw_g(Wfd%paral_kgb,nspinor,Ep%npwepG0,nfftf_tot,ngfftf,1,use_padfftf,igfftepsG0f,gboundf,&
&             ur_ae1,ur_ae_onsite1,ur_ps_onsite1,itim_kmq,tabrf_kmq,ph_mkmqt,spinrot_kmq,&
&             ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k  ,tabrf_k  ,ph_mkt  ,spinrot_k,&
&             dim_rtwg,rhotwg,tim_fourdp,Wfd%MPI_enreg)
           end if
         end if

         green_w(1)=-cone
         call assemblychi0_sym(ik_bz,nspinor,Epp,Ltg_q,green_w,Ep%npwepG0,rhotwg,Gsph_epsG0,deltaI)

       end do !ib2
     end do !ib1

     ABI_DEALLOCATE(gw_gbound)
     istat = ABI_ALLOC_STAT
     if (Dtset%pawcross==1) then
       ABI_DEALLOCATE(gboundf)
       istat = ABI_ALLOC_STAT
     end if
   end do !ik_bz
   ! Deallocation of arrays private to loop
   ABI_DEALLOCATE(igfftepsG0)
   ABI_DEALLOCATE(wfr1)
   ABI_DEALLOCATE(wfr2)
   ABI_DEALLOCATE(rhotwg)
   ABI_DEALLOCATE(tabr_k)
   ABI_DEALLOCATE(tabr_kmq)
   ABI_DEALLOCATE(wfwfg)
   if (Psps%usepaw==1) then
     call cprj_free(Cprj2_k  )
     ABI_DEALLOCATE(Cprj2_k)
     call cprj_free(Cprj1_kmq)
     ABI_DEALLOCATE(Cprj1_kmq)
     if (Dtset%pawcross==1) then
       ABI_DEALLOCATE(ur_ae1)
       ABI_DEALLOCATE(ur_ae_onsite1)
       ABI_DEALLOCATE(ur_ps_onsite1)
       ABI_DEALLOCATE(ur_ae2)
       ABI_DEALLOCATE(ur_ae_onsite2)
       ABI_DEALLOCATE(ur_ps_onsite2)
       ABI_DEALLOCATE(tabrf_k)
       ABI_DEALLOCATE(tabrf_kmq)
       ABI_DEALLOCATE(gboundf)
       istat = ABI_ALLOC_STAT
       ABI_DEALLOCATE(igfftepsG0f)
     end if
   end if
 end do !is


 call wfd_barrier(Wfd)
 !
 ! === After big loop over transitions, now MPI ===
 ! * Master took care of the contribution in case of metallic|spin polarized systems.
 ! * Collective sum of the contributions of each node.
   call xsum_mpi(deltaI(:,:,1),comm,ierr)
 !
 ! Divide by the volume
 !deltaI=deltaI*weight/Cryst%ucvol
 deltaI=deltaI*weight

 call wfd_barrier(Wfd)
 !
 ! *************************************************
 ! **** Now each node has deltaI(q,G,Gp,1) ****
 ! *************************************************

 ! Impose Hermiticity (valid only for zero or purely imaginary frequencies)
     do ig2=1,Ep%npwe
       do ig1=1,ig2-1
        deltaI(ig2,ig1,1)=CONJG(deltaI(ig1,ig2,1))
       end do
     end do
 !
 ! === Symmetrize deltaI in case of AFM system ===
 ! * Reconstruct $deltaI{\down,\down}$ from $deltaI{\up,\up}$.
 ! * Works only in case of magnetic group Shubnikov type IV.
 if (Cryst%use_antiferro) then
   call symmetrize_afm_chi0(Cryst,Gsph_epsG0,Ltg_q,Ep%npwe,fake_nomega,deltaI)
 end if
 !
 ! =====================
 ! ==== Free memory ====
 ! =====================
 ABI_DEALLOCATE(bbp_ks_distrb)

 if (allocated(gw_gfft       ))  then
   ABI_DEALLOCATE(gw_gfft)
 end if
 if (allocated(kkweight      ))  then
   ABI_DEALLOCATE(kkweight)
 end if
 if (allocated(gspfft_igfft  ))  then
   ABI_DEALLOCATE(gspfft_igfft)
 end if

 call destroy_gsphere(Gsph_FFT)

 if (Psps%usepaw==1) then ! deallocation for PAW.
   call destroy_paw_pwij(Pwij)
   ABI_DEALLOCATE(Pwij)
   if (allocated(Pwij_fft)) then
     call destroy_paw_pwij(Pwij_fft)
     ABI_DEALLOCATE(Pwij_fft)
   end if
 end if

 call timab(331,2,tsec)

 DBG_EXIT("COLL")

end subroutine check_completeness
!!***
