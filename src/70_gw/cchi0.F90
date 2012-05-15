!{\src2tex{textfont=tt}}
!!****f* ABINIT/cchi0
!! NAME
!! cchi0
!!
!! FUNCTION
!! Main calculation of the independent-particle susceptibility chi0 for qpoint!=0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
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
!!    %nomega=total number of frequencies in X^0 (both real and imaginary)
!!    %nomegasf=number of real frequencies used to sample the imaginary part of X^0 (spectral method)
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
!! Gsph_wfn<gvectors_type data type> The G-sphere used to describe the wavefunctions.
!! Gsph_epsG0<gvectors_type data type> The G-sphere used to describe chi0/eps. (including umklapp G0 vectors)
!!    %ng=number of G vectors for chi0
!!    %rottbm1(ng,2,nsym)=contains the index (IS^{-1}) G
!!    %phmGt(Ep%npwe,nsym)=phase factors e^{-iG \cdot t} needed to symmetrize oscillator matrix elements and epsilon
!!    %gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!!    %gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!! nbvw=number of bands in the arrays wfrv
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
!!  chi0(Ep%npwe,Ep%npwe,Ep%nomega)=independent-particle susceptibility matrix at wavevector qpoint and
!!   each frequeny defined by Ep%omega and Ep%nomega
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      accumulate_chi0sumrule,approxdelta,assemblychi0_sym,assemblychi0sf
!!      calc_wfwfg,chi0_bbp_mask,completechi0_deltapart,cprj_alloc,cprj_free
!!      destroy_gsphere,destroy_paw_pwij,get_bz_diff,get_bz_item,get_gftt
!!      gsph_fft_tabs,gsph_in_fftbox,gw_eet_chi0,hilbert_transform
!!      init_paw_pwij,make_transitions,nullify_gsphere,paw_cross_rho_tw_g
!!      paw_rho_tw_g,paw_symcprj,print_little_group,rho_tw_g,setup_spectral
!!      symmetrize_afm_chi0,timab,wfd_barrier,wfd_change_ngfft
!!      wfd_distribute_kb_kpbp,wfd_get_cprj,wfd_get_ur,wfd_paw_get_aeur,wrtout
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cchi0(use_tr,Dtset,Cryst,qpoint,Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,Gsph_wfn,&
& Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,nbvw,ngfft_gw,nfftot_gw,ngfftf,nfftf_tot,&
& chi0,ktabr,ktabrf,Ltg_q,chi0_sumrule,Wfd,Wfdf)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_commutator_vkbr
 use m_timer
 use m_errors
#ifdef HAVE_CLIB
 use m_clib
#endif

 use m_gwdefs,        only : GW_TOL_DOCC, GW_TOL_W0, czero_gw, epsilonm1_parameters, g0g0w
 use m_numeric_tools, only : imin_loc
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_structure
 use m_bz_mesh,       only : bz_mesh_type, get_BZ_item, get_BZ_diff, little_group, print_little_group
 use m_gsphere,       only : gvectors_type, gsph_fft_tabs, nullify_gsphere, destroy_gsphere, gsph_in_fftbox
 use m_fft_mesh,      only : get_gftt
 use m_paw_pwij,      only : paw_pwff_type, paw_pwij_type, init_paw_pwij, destroy_paw_pwij, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_wfs,           only : wfd_get_ur, wfs_descriptor, wfd_distribute_kb_kpbp, wfd_get_cprj, wfd_barrier, wfd_change_ngfft,&
&                            wfd_paw_get_aeur
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_io_tools,      only : flush_unit, get_unit
 use m_paw_toolbox,   only : paw_pwaves_lmn_t
!$ use omp_lib, only: omp_get_thread_num, omp_get_num_threads

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cchi0'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_66_paw
 use interfaces_70_gw, except_this_one => cchi0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbvw,nfftot_gw,nfftf_tot
 logical,intent(in) :: use_tr
 type(Bandstructure_type),intent(in) :: QP_BSt
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Crystal_structure),intent(in) :: Cryst
 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(Gvectors_type),intent(in) :: Gsph_epsG0,Gsph_wfn
 type(Little_group),intent(in) :: Ltg_q
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wfd,Wfdf
!arrays
 integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz),ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
 integer,intent(in) :: ngfft_gw(18),ngfftf(18)
 real(dp),intent(in) :: qpoint(3)
 real(dp),intent(out) :: chi0_sumrule(Ep%npwe)
 complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(Paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=1,two_poles=2,one_pole=1
 integer :: dim_rtwg,ib,ib1,ib2,idle,ierr
 integer :: ig1,ig2,ik_bz,ik_ibz,ikmq_bz,ikmq_ibz
 integer :: io,iomegal,iomegar,is,istat,isym_k
 integer :: isym_kmq,itim_k,itim_kmq,my_wl,my_wr
 integer :: nfound,nkpt_summed,nspinor,nsppol
 integer :: comm,gw_mgfft,use_padfft,gw_fftalga,mgfftf,use_padfftf
 integer :: my_nbbp,my_nbbpks,spin
 integer :: nbmax,nbopt,dummy
 real(dp) :: deltaeGW_b1kmq_b2k,deltaeGW_enhigh_b2k,deltaf_b1kmq_b2k
 real(dp) :: e_b1_kmq,en_high,f_b1_kmq,factor,max_rest,min_rest,my_max_rest
 real(dp) :: my_min_rest,numerator,spin_fact,weight,wl,wr
 real(dp) :: gw_gsq,memreq
 complex(dpc) :: ph_mkmqt,ph_mkt
 complex(gwpc) :: local_czero_gw
 logical :: qzero,ltest,isirred_k,isirred_kmq
 character(len=500) :: msg,allup
 type(Gvectors_type) :: Gsph_FFT
!arrays
 integer :: G0(3),umklp_k(3),umklp_kmq(3)
 integer :: wtk_ltg(Kmesh%nbz),got(Wfd%nproc)
 integer,allocatable :: tabr_k(:),tabr_kmq(:),tabrf_k(:),tabrf_kmq(:)
 integer,allocatable :: igfftepsG0(:),gspfft_igfft(:),igfftepsG0f(:)
 integer,allocatable :: gw_gfft(:,:),gw_gbound(:,:),dummy_gbound(:,:),gboundf(:,:)
 integer,allocatable :: bbp_ks_distrb(:,:,:,:)
 integer,allocatable :: bbp_ks_distrb_eet(:,:,:)
 real(dp) :: kbz(3),kmq_bz(3),spinrot_k(4),spinrot_kmq(4),q0(3),tsec(2)
 real(dp),pointer :: qp_energy(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: omegasf(:)
 complex(dpc),allocatable :: green_enhigh_w(:),green_w(:),kkweight(:,:)
 complex(gwpc),allocatable :: sf_chi0(:,:,:),rhotwg(:),wfr1(:)
 complex(gwpc),allocatable :: wfr2(:),wfwfg(:)
 complex(gwpc),allocatable :: ur_ae1(:),ur_ae_onsite1(:),ur_ps_onsite1(:)
 complex(gwpc),allocatable :: ur_ae2(:),ur_ae_onsite2(:),ur_ps_onsite2(:)
 logical,allocatable :: bbp_mask(:,:)
 type(Cprj_type),allocatable :: Cprj1_kmq(:,:),Cprj2_k(:,:)
 type(Paw_pwij_type),allocatable :: Pwij(:),Pwij_fft(:)
!************************************************************************

 DBG_ENTER("COLL")

 call timab(331,1,tsec) ! cchi0
 ABI_TIMER_START("")

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

 qp_energy => QP_BSt%eig(:,:,:)
 qp_occ    => QP_BSt%occ(:,:,:)
 !
 ! === Initialize the completeness correction  ===
 call nullify_gsphere(Gsph_FFT)

 if (Ep%gwcomp==1) then
   en_high=MAXVAL(qp_energy(Ep%nbnds,:,:)) + Ep%gwencomp
   write(msg,'(a,f8.2,a)')' Using completeness correction with the energy ',en_high*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
   !
   ! Allocation of wfwfg and green_enhigh_w moved inside openmp loop
   !
   ! Init the largest G-sphere contained in the FFT box for the wavefunctions.
   call gsph_in_fftbox(Gsph_FFT,Cryst,Wfd%ngfft)

   !call print_gsphere(Gsph_FFT,unit=std_out,prtvol=10)

!BEGINDEBUG
   ltest = .TRUE.
   do ib=1,Ep%npwe 
      if (ANY (Gsph_FFT%gvec(:,ib) /= Gsph_wfn%gvec(:,ib)) ) then
        write(std_out,*)ib, Gsph_FFT%gvec(:,ib), Gsph_wfn%gvec(:,ib)
        ltest = .FALSE.
      end if
   end do
   ABI_CHECK(ltest,"Bug in Gsph_FFT")
!ENDDEBUG

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
 end if
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

 write(msg,'(a,i2,2a,i2)')&
&  ' Using spectral method for the imaginary part = ',Ep%spmeth,ch10,&
&  ' Using symmetries to sum only over the IBZ_q  = ',Ep%symchi
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

     call chi0_bbp_mask(Ep,use_tr,QP_BSt,Wfd%mband,ikmq_ibz,ik_ibz,spin,spin_fact,bbp_mask)

     call wfd_distribute_kb_kpbp(Wfd,ikmq_ibz,ik_ibz,spin,allup,my_nbbp,bbp_ks_distrb(:,:,ik_bz,spin),got,bbp_mask) 
     my_nbbpks = my_nbbpks + my_nbbp
   end do
 end do

 ABI_DEALLOCATE(bbp_mask)

 write(msg,'(a,i0,a)')" Will sum ",my_nbbpks," (b,b',k,s) states in chi0."
 call wrtout(std_out,msg,'PERS')

 if (dtset%gw_eet/=-1) then
   !@Arjan: do you really need to reset got here? This array is used to optimize the load distribution 
   ! taking into account the number of tasks that have been already assigned to the different procs
   got=0 
   ABI_ALLOCATE(bbp_ks_distrb_eet,(nbvw,Kmesh%nbz,nsppol))
   do spin=1,nsppol
     do ik_bz=1,Kmesh%nbz
       if (Ep%symchi==1) then
         if (Ltg_q%ibzq(ik_bz)/=1) CYCLE  ! Only IBZ_q
       end if
       call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)
       bbp_ks_distrb_eet(:,ik_bz,spin) = xmpi_undefined_rank
       do ib=1,nbvw
         if (Wfd%nproc==1) then
           bbp_ks_distrb_eet(ib,ik_bz,spin)=0
         else
           idle = imin_loc(got)
           got(idle) = got(idle) + 1
           bbp_ks_distrb_eet(ib,ik_bz,spin)=idle-1
         end if
       end do ! ib
     enddo
   enddo
 endif

 if (Psps%usepaw==1) then
   ABI_ALLOCATE(Pwij,(Psps%ntypat))
   call init_paw_pwij(Pwij,Ep%npwepG0,qpoint,Gsph_epsG0%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
   ! Allocate statements moved to inside openmp loop

 end if


 SELECT CASE (Ep%spmeth)

 CASE (0)
   call wrtout(std_out,' Calculating chi0(q,omega,G,G")','COLL')
   ! Allocation of green_w moved inside openmp loop 

 CASE (1,2)
   call wrtout(std_out,' Calculating Im chi0(q,omega,G,G")','COLL')
   !
   ! Find Max and min resonant transitions for this q, report also treated by this proc.
   call make_transitions(Wfd,1,Ep%nbnds,nbvw,nsppol,Ep%symchi,Cryst%timrev,GW_TOL_DOCC,&
&    max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,qp_energy,qp_occ,qpoint,bbp_ks_distrb)
   !
   ! Calculate frequency dependent weights for Hilbert transform.
   ABI_ALLOCATE(omegasf,(Ep%nomegasf))
   ABI_ALLOCATE(kkweight,(Ep%nomegasf,Ep%nomega))
   !my_wl=1; my_wr=Ep%nomegasf
   call setup_spectral(Ep%nomega,Ep%omega,Ep%nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
&    0,Ep%zcut,zero,my_wl,my_wr,kkweight)

   if (.not.use_tr) then
     MSG_BUG('spectral method requires time-reversal')
   end if

   memreq = 2.0*gwpc*Ep%npwe**2*(my_wr-my_wl+1)*b2Gb
   write(msg,'(a,f10.4,a)')' memory required per spectral point: ',2.0*gwpc*Ep%npwe**2*b2Mb,' [Mb]'
   call wrtout(std_out,msg,'PERS')
   write(msg,'(a,f10.4,a)')' memory required by sf_chi0: ',memreq,' [Gb]'
   call wrtout(std_out,msg,'PERS')
   if (memreq > two) then
     MSG_WARNING(' Memory required for sf_chi0 is larger than 2.0 Gb!')
   end if
   if ((Dtset%userre > zero).and.(memreq > Dtset%userre)) then
     write(msg,'(a,f8.3,a)')' Memory required by sf_chi0 is larger than userre:',Dtset%userre,'[Gb]'
     MSG_PERS_ERROR(msg)
   end if

   ABI_ALLOCATE(sf_chi0,(Ep%npwe,Ep%npwe,my_wl:my_wr))
   istat = ABI_ALLOC_STAT
   if (istat/=0) then
      MSG_PERS_ERROR('out-of-memory n allocation sf_chi0')
   end if

   sf_chi0=czero_gw

 CASE DEFAULT
   MSG_BUG("Wrong spmeth")
 END SELECT

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
 chi0=czero_gw; chi0_sumrule=zero

#ifdef HAVE_CLIB
 !%call clib_progress_bar(-1,Kmesh%nbz)
#endif

 ! === Loop on spin to calculate trace $\chi_{up,up}+\chi_{down,down}$ ===
 ! * Only $\chi_{up,up} for AFM.
 do is=1,nsppol

   if (ALL(bbp_ks_distrb(:,:,:,is) /= Wfd%my_rank)) CYCLE

! Allocation of arrays that are private to loop
   if (Ep%gwcomp==1)  then
     ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))
   end if
   if (Ep%gwcomp==1)  then
     ABI_ALLOCATE(green_enhigh_w,(Ep%nomega))
   end if
   if (Ep%spmeth==0)  then
     ABI_ALLOCATE(green_w,(Ep%nomega))
   end if
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

     !%call clib_progress_bar(ik_bz,Kmesh%nbz)
     write(msg,'(2(a,i4),a,i2,a,i3)')' ik = ',ik_bz,' / ',Kmesh%nbz,' is = ',is,' done by processor ',Wfd%my_rank
     call wrtout(std_out,msg,'PERS')
     !
     ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt,umklp_k,isirred_k)
     call get_BZ_diff(Kmesh,kbz,qpoint,ikmq_bz,G0,nfound); if (nfound==0) call leave_new('COLL')

     ! * Get ikmq_ibz, non-symmorphic phase, ph_mkmqt, and symmetries from ikmq_bz.
     call get_BZ_item(Kmesh,ikmq_bz,kmq_bz,ikmq_ibz,isym_kmq,itim_kmq,ph_mkmqt,umklp_kmq,isirred_kmq)

!BEGIN DEBUG
     !if (ANY( g0 /= -umklp_kmq + umklp_k) ) then
     !if (ANY( g0 /= -umklp_kmq ) ) then
     !  write(msg,'(a,3(1x,3i2))')" g0 /= -umklp_kmq + umklp_k ",g0, umklp_kmq, umklp_k
     !  MSG_ERROR(msg)
     !end if
     !if (ANY( umklp_k /=0) ) then
     !  write(msg,'(a,2(1x,3i2))')" umklp_k /= 0 ",umklp_k,umklp_kmq 
     !  MSG_ERROR(msg)
     !end if
     !g0 = -umklp_k + umklp_kmq
     !g0 = +umklp_k - umklp_kmq
     !if (ANY (ABS(g0) > Ep%mg0) ) then
     !  write(msg,'(a,6(1x,i0))')"  ABS(g0) > Ep%mg0 ",g0,Ep%mg0
     !  MSG_ERROR(msg)
     !end if
!END DEBUG

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
     if (Dtset%gw_EET/=-1) then
       call gw_EET_chi0(Ep,Dtset,Cryst,Wfd,Kmesh,Gsph_epsG0,Gsph_wfn,Psps,Ltg_q,nbvw,qpoint, &
&                       nfftot_gw,ngfft_gw,use_padfft,igfftepsG0,gw_gbound,gw_mgfft,is, &
&                       ik_bz,ik_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k,ikmq_ibz,itim_kmq, &
&                       tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,qp_energy,chi0,spin_fact, &
&                       qp_occ,nspinor,tim_fourdp,bbp_ks_distrb_EET,nbopt)
       nbmax=nbopt
     end if

     do ib1=1,nbmax ! Loop over "conduction" states.

!       if ((Dtset%prtvol>9).AND.(mod(ib1,10)==0)) then
!         write(msg,'(2(a,i0))')' Doing ib1: ',ib1,' of: ',nbmax
!         call wrtout(std_out,msg,'PERS')
!       end if

       if (ALL(bbp_ks_distrb(ib1,:,ik_bz,is) /= Wfd%my_rank)) CYCLE

       call wfd_get_ur(Wfd,ib1,ikmq_ibz,is,wfr1)

       if (Psps%usepaw==1) then 
         call wfd_get_cprj(Wfd,ib1,ikmq_ibz,is,Cryst,Cprj1_kmq,sorted=.FALSE.)
         call paw_symcprj(ikmq_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj1_kmq) 
         if (Dtset%pawcross==1) then
           call wfd_paw_get_aeur(Wfdf,ib1,ikmq_ibz,is,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae1,ur_ae_onsite1,ur_ps_onsite1)
         end if
       end if

       e_b1_kmq=qp_energy(ib1,ikmq_ibz,is)
       f_b1_kmq=   qp_occ(ib1,ikmq_ibz,is)

       do ib2=1,nbmax ! Loop over "valence" states.

         if (bbp_ks_distrb(ib1,ib2,ik_bz,is) /= Wfd%my_rank) CYCLE

         deltaf_b1kmq_b2k=spin_fact*(f_b1_kmq-qp_occ(ib2,ik_ibz,is))

#if 1
         if (Ep%gwcomp==0) then ! Skip negligible transitions.
           if (ABS(deltaf_b1kmq_b2k) < GW_TOL_DOCC) CYCLE
         else
           ! * When the completeness correction is used,
           !   we need to also consider transitions with vanishing deltaf
           !if (qp_occ(ib2,ik_ibz,is) < GW_TOL_DOCC) CYCLE
           !
           ! Rangel This is to compute chi correctly when using the extrapolar method
           if (qp_occ(ib2,ik_ibz,is) < GW_TOL_DOCC .and. (ABS(deltaf_b1kmq_b2k) < GW_TOL_DOCC .or. ib1<ib2)) CYCLE
         end if
#endif

         deltaeGW_b1kmq_b2k=e_b1_kmq-qp_energy(ib2,ik_ibz,is)

         call wfd_get_ur(Wfd,ib2,ik_ibz,is,wfr2)

         if (Psps%usepaw==1) then 
           call wfd_get_cprj(Wfd,ib2,ik_ibz,is,Cryst,Cprj2_k,sorted=.FALSE.)
           call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj2_k) 
           if (Dtset%pawcross==1) then
             call wfd_paw_get_aeur(Wfdf,ib2,ik_ibz,is,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae2,ur_ae_onsite2,ur_ps_onsite2)
           end if
         end if

         SELECT CASE (Ep%spmeth)

         CASE (0) ! Standard Adler-Wiser expression.
          ! * Add the small imaginary of the Time-Ordered RF only for non-zero real omega ! FIXME What about metals?

          if (.not.use_tr) then ! Have to sum over all possible resonant and anti-resonant transitions.
            do io=1,Ep%nomega
              green_w(io) = g0g0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0,one_pole)
            end do
          else 
#if 1
            if (Ep%gwcomp==0) then ! cannot be completely skipped in case of completeness correction
              if (ib1<ib2) CYCLE ! Here we GAIN a factor ~2
            end if
#endif
            do io=1,Ep%nomega
              !Rangel: In metals, the intra-band transitions term does not contain the antiresonant part
              !green_w(io) = g0g0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0)
              if (ib1==ib2) green_w(io) = g0g0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0,one_pole)
              if (ib1/=ib2) green_w(io) = g0g0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0,two_poles)

              if (Ep%gwcomp==1) then ! Calculate the completeness correction
                numerator= -spin_fact*qp_occ(ib2,ik_ibz,is)
                deltaeGW_enhigh_b2k = en_high-qp_energy(ib2,ik_ibz,is)
                
                if (REAL(Ep%omega(io))<GW_TOL_W0) then ! Completeness correction is NOT valid for real frequencies
                  green_enhigh_w(io) = g0g0w(Ep%omega(io),numerator,deltaeGW_enhigh_b2k,Ep%zcut,GW_TOL_W0,two_poles)
                else
                  green_enhigh_w(io) = local_czero_gw
                end if
                !
                !Rangel Correction for metals
                !if (deltaf_b1kmq_b2k<0.d0) then
                if (ib1>=ib2 .and. abs(deltaf_b1kmq_b2k) > GW_TOL_DOCC ) then
                  green_w(io)= green_w(io) - green_enhigh_w(io)
                else ! Disregard green_w, since it is already accounted for through the time-reversal
                  green_w(io)=             - green_enhigh_w(io)
                end if
              end if !gwcomp==1
            end do !io
            !
            if (Ep%gwcomp==1.and.ib1==ib2) then ! Add the "delta part" of the extrapolar method. TODO doesnt work for spinor

              call calc_wfwfg(Wfd%MPI_enreg,Wfd%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw,wfr2,wfr2,wfwfg)

              if (Psps%usepaw==1) then
                call paw_rho_tw_g(nfftot_gw,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,gw_gfft,&
&                 Cprj2_k,Cprj2_k,Pwij_fft,wfwfg)

                ! Add PAW cross term
                if (Dtset%pawcross==1) then
                  call paw_cross_rho_tw_g(Wfdf%paral_kgb,nspinor,Ep%npwepG0,nfftf_tot,ngfftf,1,use_padfftf,igfftepsG0f,gboundf,&
&                  ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_kmq,tabrf_kmq,ph_mkmqt,spinrot_kmq,&
&                  ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k  ,tabrf_k  ,ph_mkt  ,spinrot_k,&
&                  dim_rtwg,wfwfg,tim_fourdp,Wfdf%MPI_enreg)
                end if
              end if

              qzero=.FALSE.
              call completechi0_deltapart(ik_bz,qzero,Ep%symchi,Ep%npwe,Gsph_FFT%ng,Ep%nomega,nspinor,&
&               nfftot_gw,ngfft_gw,gspfft_igfft,gsph_FFT,Ltg_q,green_enhigh_w,wfwfg,chi0)

            end if
          end if ! use_tr

         CASE (1,2) ! Spectral method, WARNING time-reversal here is always assumed!
#if 1
           if (deltaeGW_b1kmq_b2k<0) CYCLE
#endif
           call approxdelta(Ep%nomegasf,omegasf,deltaeGW_b1kmq_b2k,Ep%spsmear,iomegal,iomegar,wl,wr,Ep%spmeth)
         END SELECT
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
             call paw_cross_rho_tw_g(Wfdf%paral_kgb,nspinor,Ep%npwepG0,nfftf_tot,ngfftf,1,use_padfftf,igfftepsG0f,gboundf,&
&             ur_ae1,ur_ae_onsite1,ur_ps_onsite1,itim_kmq,tabrf_kmq,ph_mkmqt,spinrot_kmq,&
&             ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k  ,tabrf_k  ,ph_mkt  ,spinrot_k,&
&             dim_rtwg,rhotwg,tim_fourdp,Wfdf%MPI_enreg)
           end if
         end if

         SELECT CASE (Ep%spmeth)

         CASE (0) ! Adler-Wiser.
           call assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,Ep%npwepG0,rhotwg,Gsph_epsG0,chi0)

         CASE (1,2) ! Spectral method ! TODO Does not work with spinor
           call assemblychi0sf(ik_bz,nspinor,Ep%symchi,Ltg_q,Ep%npwepG0,Ep%npwe,rhotwg,Gsph_epsG0,&
&            deltaf_b1kmq_b2k,my_wl,iomegal,wl,my_wr,iomegar,wr,Ep%nomegasf,sf_chi0)

         CASE DEFAULT
           MSG_BUG("Wrong spmeth")
         END SELECT
         !
         ! === Accumulating the sum rule on chi0 ===
         ! *Eq. (5.284) in G.D. Mahan Many-Particle Physics 3rd edition.
         ! TODO Does not work with spinor

         factor=spin_fact*qp_occ(ib2,ik_ibz,is)
         call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_b1kmq_b2k,&
&          Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0_sumrule)
         !
         ! * Include also the completeness correction in the sum rule
         if (Ep%gwcomp==1) then
           factor=-spin_fact*qp_occ(ib2,ik_ibz,is)
           call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_enhigh_b2k,&
&            Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0_sumrule)
           if (ib1==Ep%nbnds) then
             chi0_sumrule(:)=chi0_sumrule(:) + wtk_ltg(ik_bz)*spin_fact*qp_occ(ib2,ik_ibz,is)*deltaeGW_enhigh_b2k
           end if
         end if

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
 if (allocated(green_w       ))  then
   ABI_DEALLOCATE(green_w)
 end if
 if (allocated(wfwfg         ))  then
   ABI_DEALLOCATE(wfwfg)
 end if
 if (allocated(green_enhigh_w))  then
   ABI_DEALLOCATE(green_enhigh_w)
 end if
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
 SELECT CASE (Ep%spmeth)
 CASE (0) ! Adler-Wiser
   ! * Collective sum of the contributions of each node.
   ! * Looping on frequencies to avoid problems with the size of the MPI packet
   do io=1,Ep%nomega
     call xsum_mpi(chi0(:,:,io),comm,ierr)
   end do
   !
   ! Divide by the volume
   chi0=chi0*weight/Cryst%ucvol

 CASE (1,2) ! Spectral method.

   call hilbert_transform(Ep%npwe,Ep%nomega,Ep%nomegasf,&
&   my_wl,my_wr,kkweight,sf_chi0,chi0,Ep%spmeth)

!  Deallocate here before the xsums
   if (allocated(sf_chi0       ))  then
     ABI_DEALLOCATE(sf_chi0)
   end if

   ! === Collective sum of the contributions ===
   ! * Looping over frequencies to avoid problems with the size of the MPI packet
   do io=1,Ep%nomega
    call xsum_mpi(chi0(:,:,io),comm,ierr)
   end do
   chi0=chi0*weight/Cryst%ucvol

 CASE DEFAULT
   MSG_BUG("Wrong spmeth")

 END SELECT
 !
 ! === Collect the sum rule ===
 ! * The pi factor comes from Im[1/(x-ieta)] = pi delta(x)
 call xsum_mpi(chi0_sumrule,comm,ierr)
 chi0_sumrule=chi0_sumrule*pi*weight/Cryst%ucvol

 call wfd_barrier(Wfd)
 !
 ! *************************************************
 ! **** Now each node has chi0(q,G,Gp,Ep%omega) ****
 ! *************************************************

 ! Impose Hermiticity (valid only for zero or purely imaginary frequencies)
 ! MG what about metals, where we have poles around zero?
 do io=1,Ep%nomega
   if (ABS(REAL(Ep%omega(io)))<0.00001) then
     do ig2=1,Ep%npwe
       do ig1=1,ig2-1
        chi0(ig2,ig1,io)=CONJG(chi0(ig1,ig2,io))
       end do
     end do
   end if
 end do
 !
 ! === Symmetrize chi0 in case of AFM system ===
 ! * Reconstruct $chi0{\down,\down}$ from $chi0{\up,\up}$.
 ! * Works only in case of magnetic group Shubnikov type IV.
 if (Cryst%use_antiferro) then
   call symmetrize_afm_chi0(Cryst,Gsph_epsG0,Ltg_q,Ep%npwe,Ep%nomega,chi0)
 end if
 !
 ! =====================
 ! ==== Free memory ====
 ! =====================
 ABI_DEALLOCATE(bbp_ks_distrb)
 if (dtset%gw_eet/=-1) then
   ABI_DEALLOCATE(bbp_ks_distrb_eet)
 endif

 if (allocated(gw_gfft       ))  then
   ABI_DEALLOCATE(gw_gfft)
 end if
 if (allocated(kkweight      ))  then
   ABI_DEALLOCATE(kkweight)
 end if
 if (allocated(omegasf       ))  then
   ABI_DEALLOCATE(omegasf)
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
 ABI_TIMER_STOP("")

 DBG_EXIT("COLL")

end subroutine cchi0
!!***
