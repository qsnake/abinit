!!****f* ABINIT/chi0q0_intraband
!! NAME
!! chi0q0_intraband      
!!
!! FUNCTION
!! Calculate chi0 in the limit q-->0
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  use_tr=If .TRUE. Wfs_val are allocate and only resonant transitions are evaluated (assumes time reversal symmetry)
!!  Ep= datatype gathering differening parameters related to the calculation of the inverse dielectric matrix
!!  Gsph_epsG0<gvectors_data_type>: Info on the G-sphere used to describe chi0/espilon (including umklapp)
!!    %ng=number of G vectors
!!    %rottbm1(ng,2,nsym)=contains the index (IS^{-1}) G  in the array gvec
!!    %phmGt(ng,nsym)=phase factor e^{-iG.\tau} needed to symmetrize oscillator matrix elements and chi0
!!    %gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!    %gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!!  Ep%nbnds=number of bands
!!  ngfft_gw(18)= array containing all the information for 3D FFT for the oscillator strengths.
!!  Ep%nomega=number of frequencies
!!  Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!!   %natom=number of atoms
!!   %nsym=number of symmetry operations
!!   %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!   %typat(natom)=type of each atom
!!   %xred(3,natom)=reduced coordinated of atoms
!!   %rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!   %timrev=2 if time-reversal symmetry can be used, 1 otherwise
!!  Ep%npwe=number of planewaves for sigma exchange (input variable)
!!  Ep%npwvec=dimension of igfft
!!  Ep%nsppol=1 for unpolarized, 2 for spin-polarized
!!  Ep%omega(Ep%nomega)=frequencies
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!     %mpsang=1+maximum angular momentum for nonlocal pseudopotential
!!  Pawang<pawang_type> angular mesh discretization and related data:
!!  Pawrad(ntypat*usepaw)<Pawrad_type>=paw radial mesh and related data
!!  Paw_ij(natom*usepaw)<Paw_ij_type)>=paw arrays given on (i,j) channels
!!  BSt<Bandstructure_type>=Quasiparticle energies and occupations (for the moment real quantities)
!!    %mband=MAX number of bands over k-points and spin (==Ep%nbnds)
!!    %occ(mband,nkpt,nsppol)=QP occupation numbers, for each k point in IBZ, and each band
!!    %eig(mband,nkpt,nsppol)=GW energies, for self-consistency purposes
!!  Paw_pwff<Paw_pwff_type>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!
!! OUTPUT
!!  chi0(Ep%npwe,Ep%npwe,Ep%nomega)=independent-particle susceptibility matrix for wavevector qq,
!!   and frequencies defined by Ep%omega
!!
!! NOTES
!!  *) The terms "head", "wings" and "body" of chi(G,Gp) refer to
!!     G=Gp=0, either G or Gp=0, and neither=0 respectively
!!
!! TODO
!!  Check npwepG0 before Switching on umklapp
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      assemblychi0_sym,clib_progress_bar,cprj_alloc,cprj_copy,cprj_free
!!      destroy_bz_mesh_type,destroy_hur,destroy_kb_potential
!!      destroy_little_group,destroy_paw_pwij,get_bz_item,getnel,gsph_fft_tabs
!!      init_kb_potential,init_kmesh,init_paw_pwij,make_hur_commutator
!!      nullify_hur,nullify_kb_potential,pack_eneocc,paw_rho_tw_g,paw_symcprj
!!      print_arr,print_little_group,rho_tw_g,rotate_fft_mesh
!!      setup_little_group,symmetrize_afm_chi0,unpack_eneocc,wfd_barrier
!!      wfd_change_ngfft,wfd_distribute_bands,wfd_get_cprj,wfd_get_ur,wrtout
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chi0q0_intraband(Wfd,Cryst,Ep,Psps,BSt,Gsph_epsG0,Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,use_tr,usepawu,&
&  ngfft_gw,chi0,chi0_head,chi0_lwing,chi0_uwing)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
#ifdef HAVE_CLIB
 use m_clib
#endif
 use m_errors

 use m_gwdefs,          only : GW_TOL_DOCC, GW_TOL_W0, czero_gw, epsilonm1_parameters, g0g0w
 use m_geometry,        only : vdotw
 use m_numeric_tools,   only : print_arr
 use m_crystal,         only : crystal_structure
 use m_fft_mesh,        only : rotate_FFT_mesh
 use m_ebands,          only : pack_eneocc, unpack_eneocc
 use m_bz_mesh,         only : bz_mesh_type, init_kmesh, destroy_bz_mesh_type, get_BZ_item, &
&                              little_group, print_little_group, destroy_little_group, setup_little_group
 use m_gsphere,         only : gvectors_type, gsph_fft_tabs
 use m_paw_pwij,        only : paw_pwff_type, paw_pwij_type, init_paw_pwij, destroy_paw_pwij, paw_rho_tw_g
 use m_wfs
 use m_oscillators,     only : rho_tw_g
 use m_paw_commutator,  only : HUr_commutator, destroy_Hur, nullify_Hur, make_Hur_commutator, paw_ihr_comm
 use m_commutator_vkbr, only : kb_potential, nullify_kb_potential, destroy_kb_potential, init_kb_potential, nc_ihr_comm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chi0q0_intraband'
 use interfaces_14_hidewrite
 use interfaces_44_abitypes_defs
 use interfaces_62_occeig
 use interfaces_66_paw
 use interfaces_70_gw, except_this_one => chi0q0_intraband
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: usepawu
 logical,intent(in) :: use_tr
 type(Bandstructure_type),intent(in) :: BSt
 type(Crystal_structure),intent(in) :: Cryst
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawang_type),intent(in) :: Pawang
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfft_gw(18)
 complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 complex(dpc),intent(out) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
 complex(dpc),intent(out) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)
 complex(dpc),intent(out) :: chi0_head(3,3,Ep%nomega)
 type(Pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(Paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Psps%usepaw)
 type(Paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp1=1,two_poles=2,one_pole=1
 integer,parameter :: unitdos0=0,option1=1,NOMEGA_PRINTED=15
 integer :: nqlwl,nband_k,iomega,istwf_k,npw_k,my_nband,lbidx
 integer :: band,itim_k,ik_bz,ik_ibz,io,isym_k,spin,istat,iqlwl!,gw_eet !ig,ig1,ig2,my_nbbp,my_nbbpks
 integer :: nkpt_summed,dim_rtwg,use_padfft,gw_fftalga,ifft
 integer :: kptopt,isym,nsppol,nspinor
 integer :: comm,ierr,gw_mgfft,use_umklp,inclvkb
 real(dp) :: spin_fact,deltaf_b1b2,weight
 real(dp) :: deltaeGW_b1b2,zcut
 real(dp),parameter :: dummy_dosdeltae=HUGE(zero)
 real(dp) :: o_entropy,o_nelect,maxocc
 complex(dpc) :: ph_mkt
 logical :: iscompatibleFFT !,ltest 
 character(len=500) :: msg,msg_tmp !,allup
 type(BZ_mesh_type) :: Kmesh
 type(Little_group) :: Ltg_q
 type(kb_potential) :: KBgrad_k
!arrays
 integer :: my_band_list(Wfd%mband)
 integer,pointer :: kg_k(:,:)
 integer,allocatable :: ktabr(:,:),irottb(:,:)
 !integer :: got(Wfd%nproc)
 integer,allocatable :: tabr_k(:),igffteps0(:),gw_gbound(:,:)
 real(dp),parameter :: q0(3)=(/zero,zero,zero/)
 real(dp) :: kpt(3),dedk(3),kbz(3),spinrot_kbz(4) 
 !real(dp),pointer :: ks_energy(:,:,:),qp_energy(:,:,:),qp_occ(:,:,:)
 real(dp) :: shift_ene(BSt%mband,BSt%nkpt,BSt%nsppol)
 real(dp) :: delta_occ(BSt%mband,BSt%nkpt,BSt%nsppol)
 !real(dp) :: eigen_vec(BSt%bantot)
 real(dp) :: o_doccde(BSt%bantot) 
 real(dp) :: eigen_pdelta_vec(BSt%bantot),eigen_mdelta_vec(BSt%bantot)
 real(dp) :: o_occ_pdelta(BSt%bantot),o_occ_mdelta(BSt%bantot)
 real(dp) :: delta_ene(BSt%mband,BSt%nkpt,BSt%nsppol)
 real(dp) :: test_docc(BSt%mband,BSt%nkpt,BSt%nsppol)
 real(dp),allocatable :: qlwl(:,:)
 complex(gwpc) :: comm_kbbs(3,Wfd%nspinor**2)
 complex(dpc),allocatable :: ihr_comm(:,:,:,:,:)
 complex(gwpc),allocatable :: rhotwg(:)
 complex(dpc) :: green_w(Ep%nomega)
 complex(gwpc),allocatable :: ur1(:)
 complex(gwpc),pointer :: ug(:) 
 logical :: bmask(Wfd%mband)
 type(Cprj_type),allocatable :: Cprj1_bz(:,:),Cprj1_ibz(:,:),Cp_bks(:,:)
 type(Paw_pwij_type),allocatable :: Pwij(:)
 type(HUr_commutator),allocatable :: Hur(:)

!************************************************************************

 DBG_ENTER("COLL")

 nsppol  = Wfd%nsppol
 nspinor = Wfd%nspinor

 gw_mgfft = MAXVAL(ngfft_gw(1:3))
 gw_fftalga = ngfft_gw(7)/100 !; gw_fftalgc=MOD(ngfft_gw(7),10)
 !
 ! Calculate <k,b1|i[H,r]|k',b2>.
 inclvkb=2; if (Wfd%usepaw==1) inclvkb=0
 ABI_ALLOCATE(ihr_comm,(3,nspinor**2,Wfd%mband,Wfd%nkibz,nsppol))
 ihr_comm = czero

 if (Wfd%usepaw==1) then 
   ABI_ALLOCATE(Cp_bks,(Cryst%natom,nspinor))
   call cprj_alloc(Cp_bks,0,Wfd%nlmn_atm)
   ABI_ALLOCATE(HUr,(Cryst%natom))
   call nullify_Hur(HUr)
   if (usepawu/=0) then ! For PAW+LDA+U, precalculate <\phi_i|[Hu,r]|phi_j\>.
     call make_Hur_commutator(nsppol,Wfd%pawprtvol,Cryst,Psps,Pawtab,Pawang,Pawrad,Paw_ij,Hur)
   end if
 end if

 do spin=1,nsppol
   do ik_ibz=1,Wfd%nkibz
     npw_k  =  Wfd%npwarr(ik_ibz)
     nband_k=  Wfd%nband(ik_ibz,spin)
     kpt    =  Wfd%kibz(:,ik_ibz)
     kg_k   => Wfd%Kdata(ik_ibz)%kg_k

     istwf_k = Wfd%istwfk(ik_ibz)
     ABI_CHECK(istwf_k==1,"istwf_k/=1 not coded")
     !
     ! Distribute bands.
     bmask=.FALSE.; bmask(1:nband_k)=.TRUE. ! TODO only bands around EF should be included.
     call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,bmask=bmask)
     if (my_nband==0) CYCLE
     
     call nullify_kb_potential(KBgrad_k)
     if (Wfd%usepaw==0.and.inclvkb/=0) then ! Include term <n,k|[Vnl,iqr]|n"k>' for q->0.
       call init_kb_potential(KBgrad_k,Cryst,Psps,inclvkb,istwf_k,npw_k,kpt,kg_k)
     end if

     do lbidx=1,my_nband
       band=my_band_list(lbidx)
       ug => Wfd%Wave(band,ik_ibz,spin)%ug

       if (Wfd%usepaw==0) then  ! Matrix elements of i[H,r] for NC pseudopotentials.        
         comm_kbbs = nc_ihr_comm(nspinor,npw_k,istwf_k,inclvkb,Kmesh%ibz(:,ik_ibz),KBgrad_k,ug,ug,kg_k) 
       else                     ! Matrix elements of i[H,r] for PAW.
         call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cp_bks,sorted=.FALSE.)
         comm_kbbs = paw_ihr_comm(spin,nspinor,npw_k,istwf_k,Kmesh%ibz(:,ik_ibz),Cryst,Pawtab,ug,ug,kg_k,Cp_bks,Cp_bks,HUr)
       end if

       ihr_comm(:,:,band,ik_ibz,spin) = comm_kbbs
     end do

     call destroy_kb_potential(KBgrad_k) ! Not need anymore as we loop only over IBZ.
   end do
 end do
 !
 ! Gather the commutator on each node.
 call xsum_mpi(ihr_comm,Wfd%comm,ierr)

 if (Wfd%usepaw==1) then 
   call cprj_free(Cp_bks)
   ABI_DEALLOCATE(Cp_bks)
   call destroy_Hur(Hur)
 end if

 nqlwl=1
 ABI_ALLOCATE(qlwl,(3,nqlwl))
 !qlwl = GW_Q0_DEFAULT(3) 
 qlwl(:,1) = (/0.00001_dp, 0.00002_dp, 0.00003_dp/)
 !
 write(msg,'(a,i3,a)')' Q-points for long wave-length limit in chi0q_intraband. # ',nqlwl,ch10
 do iqlwl=1,nqlwl
   write(msg_tmp,'(1x,i5,a,2x,3f12.6,a)') iqlwl,')',qlwl(:,iqlwl),ch10
   msg=TRIM(msg)//msg_tmp
 end do
 call wrtout(std_out,msg,'COLL')
 !
 ! delta_ene =  e_{b,k-q} - e_{b,k} = -q. <b,k| i[H,r] |b,k> + O(q^2).
 delta_ene = zero
 do spin=1,nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)
       dedk = REAL(ihr_comm(:,1,band,ik_ibz,spin))
       delta_ene(band,ik_ibz,spin) = -vdotw(qlwl(:,1),dedk,Cryst%gmet,"G")  
     end do
   end do
 end do

 maxocc=two/(nsppol*nspinor)

 ! Calculate the occupations at f(e+delta/2).
 shift_ene = BSt%eig + half*delta_ene

 call pack_eneocc(BSt%nkpt,BSt%nsppol,BSt%mband,BSt%nband,BSt%bantot,shift_ene,eigen_pdelta_vec)
 
 call getnel(o_doccde,dummy_dosdeltae,eigen_pdelta_vec,o_entropy,BSt%fermie,maxocc,BSt%mband,BSt%nband,&
&  o_nelect,BSt%nkpt,BSt%nsppol,o_occ_pdelta,BSt%occopt,option1,BSt%tphysel,BSt%tsmear,unitdos0,BSt%wtk)
 write(std_out,*)"nelect1: ",o_nelect
 ! 
 ! Calculate the occupations at f(e-delta/2).
 shift_ene = BSt%eig - half*delta_ene

 call pack_eneocc(BSt%nkpt,BSt%nsppol,BSt%mband,BSt%nband,BSt%bantot,shift_ene,eigen_mdelta_vec)

 call getnel(o_doccde,dummy_dosdeltae,eigen_mdelta_vec,o_entropy,BSt%fermie,maxocc,BSt%mband,BSt%nband,&
&  o_nelect,BSt%nkpt,BSt%nsppol,o_occ_mdelta,BSt%occopt,option1,BSt%tphysel,BSt%tsmear,unitdos0,BSt%wtk)
 write(std_out,*)"nelect2: ",o_nelect
 !
 ! f(e-delta/2) - f(e+delta/2).
 o_occ_pdelta = o_occ_mdelta - o_occ_pdelta

 call unpack_eneocc(BSt%nkpt,BSt%nsppol,BSt%mband,BSt%nband,BSt%bantot,o_occ_pdelta,delta_occ)
 !
 ! Expand f(e-delta/2) - f(e+delta/2) up to the first order in the small q.
 do spin=1,nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)
       dedk = REAL(ihr_comm(:,1,band,ik_ibz,spin))
       test_docc(band,ik_ibz,spin) = +vdotw(qlwl(:,1),dedk,Cryst%gmet,"G") * BSt%doccde(band,ik_ibz,spin) 
       write(std_out,'(a,3(i0,1x),1x,3es16.8)')" spin,ik_ibz,band, delta_occ: ",&
&      spin,ik_ibz,band,delta_occ(band,ik_ibz,spin),&
&      test_docc(band,ik_ibz,spin),delta_occ(band,ik_ibz,spin)-test_docc(band,ik_ibz,spin)
     end do
   end do
 end do

! MSG_ERROR("DONE")
! do spin=1,nsppol
!   do ik_ibz=1,Wfd%nkibz
!     nband_k = Wfd%nband(ik_ibz,spin)
!     do band=1,nband_k
!       write(std_out,'(a,3i3,2es14.6)')" spin, band, ik_ibz, delta_ene, delta_occ ",&
!&        spin,band,ik_ibz,delta_ene(band,ik_ibz,spin),delta_occ(band,ik_ibz,spin) 
!     end do
!   end do 
! end do

 ABI_DEALLOCATE(ihr_comm)
 ABI_DEALLOCATE(qlwl)

 if ( ANY(ngfft_gw(1:3) /= Wfd%ngfft(1:3)) ) call wfd_change_ngfft(Wfd,Cryst,Psps,ngfft_gw) 

 ! TODO take into account the case of random k-meshes.
 kptopt=3
 call init_kmesh(Kmesh,Cryst,Wfd%nkibz,Wfd%kibz,kptopt)
 !
 !=== Get the FFT index of $ (R^{-1}(r-\tau)) $ ===
 !* S= $\transpose R^{-1}$ and k_BZ = S k_IBZ
 !* irottb is the FFT index of $ R^{-1} (r-\tau) $ used to symmetrize u_Sk.
 ABI_ALLOCATE(irottb,(Wfd%nfftot,Cryst%nsym))

 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,irottb,iscompatibleFFT)
 ABI_CHECK(iscompatibleFFT,"FFT mesh not compatible with symmetries")

 ABI_ALLOCATE(ktabr,(Wfd%nfftot,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym=Kmesh%tabo(ik_bz)
   do ifft=1,Wfd%nfftot
     ktabr(ifft,ik_bz)=irottb(ifft,isym)
   end do
 end do
 ABI_DEALLOCATE(irottb)
 !
 ! === Setup weight (2 for spin unpolarized systems, 1 for polarized) ===
 ! * spin_fact is used to normalize the occupation factors to one.
 ! * Consider also the AFM case.
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

 use_umklp=0
 call setup_little_group(q0,Kmesh,Cryst,use_umklp,Ltg_q,Ep%npwepG0,gvec=Gsph_epsG0%gvec)

 write(msg,'(a,i2)')' Using symmetries to sum only over the IBZ_q  = ',Ep%symchi
 call wrtout(std_out,msg,'COLL')
 !
 ! === Evaluate oscillator matrix elements btw partial waves. Note that q=Gamma is used.
 if (Psps%usepaw==1) then
   ABI_ALLOCATE(Pwij,(Psps%ntypat))
   call init_paw_pwij(Pwij,Ep%npwepG0,(/zero,zero,zero/),Gsph_epsG0%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)

   ABI_ALLOCATE(Cprj1_bz ,(Cryst%natom,nspinor))
   call cprj_alloc(Cprj1_bz, 0,Wfd%nlmn_atm)
   ABI_ALLOCATE(Cprj1_ibz,(Cryst%natom,nspinor))
   call cprj_alloc(Cprj1_ibz,0,Wfd%nlmn_atm)
 end if

 ABI_ALLOCATE(rhotwg,(Ep%npwe*nspinor**2))
 ABI_ALLOCATE(tabr_k,(Wfd%nfftot))
 ABI_ALLOCATE(ur1,(Wfd%nfft*nspinor))
 !
 ! Tables for the FFT of the oscillators.
 !  a) FFT index of the G sphere (only vertical transitions, unlike cchi0, no need to shift the sphere).
 !  b) gw_gbound table for the zero-padded FFT performed in rhotwg. 
 ABI_ALLOCATE(gw_gbound,(2*gw_mgfft+8,2))
 ABI_ALLOCATE(igffteps0,(Gsph_epsG0%ng))

 call gsph_fft_tabs(Gsph_epsG0,(/0,0,0/),gw_mgfft,ngfft_gw,use_padfft,gw_gbound,igffteps0)
 if ( ANY(gw_fftalga == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
 if (use_padfft==0) then 
   ABI_DEALLOCATE(gw_gbound)
   ABI_ALLOCATE(gw_gbound,(2*gw_mgfft+8,2*use_padfft))
 end if

 nkpt_summed=Kmesh%nbz
 if (Ep%symchi/=0) then
   nkpt_summed=Ltg_q%nibz_ltg
   call print_little_group(Ltg_q,std_out,Wfd%prtvol,'COLL')
 end if

#ifdef HAVE_CLIB
 call clib_progress_bar(-1,Kmesh%nbz)
#else
 write(msg,'(a,i0,a)')' Calculation status ( ',nkpt_summed,' to be completed):'
 call wrtout(std_out,msg,'COLL')
#endif
 !
 ! ============================================
 ! === Begin big fat loop over transitions ====
 ! ============================================
 chi0      =czero_gw
 chi0_head =czero_gw
 chi0_lwing=czero_gw
 chi0_uwing=czero_gw
 dim_rtwg=1; if (nspinor==2) dim_rtwg=4 !can reduce size depending on Ep%nI and Ep%nj

 zcut = Ep%zcut
 zcut = 0.1/Ha_eV
 write(std_out,*)" using zcut ",zcut*Ha_eV," [eV]"

 ! === Loop on spin to calculate $ \chi_{\up,\up} + \chi_{\down,\down} $ ===
 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz ! Loop over k-points in the BZ.
     if (Ep%symchi==1) then
       if (Ltg_q%ibzq(ik_bz)/=1) CYCLE ! Only IBZ_q
     end if
     !
     ! * Get ik_ibz, non-symmorphic phase and symmetries from ik_bz.
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt)
     tabr_k=ktabr(:,ik_bz) ! Table for rotated FFT points
     spinrot_kbz(:)=Cryst%spinrot(:,isym_k)
     nband_k=Wfd%nband(ik_ibz,spin)
     !
     ! Distribute bands.
     bmask=.FALSE.; bmask(1:nband_k)=.TRUE. ! TODO only bands around EF should be included.
     call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,bmask=bmask)
     if (my_nband==0) CYCLE

     write(msg,'(2(a,i4),a,i2,a,i3)')' ik = ',ik_bz,' / ',Kmesh%nbz,' spin = ',spin,' done by processor ',Wfd%my_rank
     call wrtout(std_out,msg,'PERS')

     do lbidx=1,my_nband  ! Loop over bands treated by this node.
       band=my_band_list(lbidx)

       call wfd_get_ur(Wfd,band,ik_ibz,spin,ur1)

       if (Psps%usepaw==1) then 
         call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cprj1_ibz,sorted=.FALSE.)
         call cprj_copy(Cprj1_ibz,Cprj1_bz)
         call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj1_bz) 
       end if

       deltaf_b1b2  = spin_fact*delta_occ(band,ik_ibz,spin)
       deltaeGW_b1b2= delta_ene(band,ik_ibz,spin)

       ! Add small imaginary of the Time-Ordered resp function but only for non-zero real omega  FIXME What about metals?
       if (.not.use_tr) then 
         do io=1,Ep%nomega
           !green_w(io) = g0g0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,zcut,-one,one_pole)
           green_w(io) = g0g0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,zcut,GW_TOL_W0,one_pole)
         end do
       else 
         do io=1,Ep%nomega ! This expression implements time-reversal even when the input k-mesh breaks it.
           !green_w(io) = half * g0g0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,zcut,-one,two_poles)
           green_w(io) = half * g0g0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,zcut,GW_TOL_W0,two_poles)
         end do !io
       end if ! use_tr
       !
       ! FFT of u^*_{b1,k}(r) u_{b2,k}(r).
       call rho_tw_g(Wfd%paral_kgb,nspinor,Ep%npwe,Wfd%nfft,ngfft_gw,1,use_padfft,igffteps0,gw_gbound,&
&        ur1,itim_k,tabr_k,ph_mkt,spinrot_kbz,&
&        ur1,itim_k,tabr_k,ph_mkt,spinrot_kbz,&
&        dim_rtwg,rhotwg,tim_fourdp1,Wfd%MPI_enreg)

       if (Psps%usepaw==1) then  ! Add PAW onsite contribution, projectors are already in the BZ.
         call paw_rho_tw_g(Ep%npwe,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_epsG0%gvec,&
&          Cprj1_bz,Cprj1_bz,Pwij,rhotwg)
       end if

       ! ==== Adler-Wiser expression, to be consistent here we use the KS eigenvalues (?) ====
!       gw_eet=-1
!       call accumulate_chi0_q0(ik_bz,isym_k,itim_k,Ep%gwcomp,gw_eet,nspinor,Ep%npwepG0,Ep,&
!&        Cryst,Ltg_q,Gsph_epsG0,chi0,rhotwx(:,1),rhotwg,green_w,green_enhigh_w,deltaf_b1b2,chi0_head,chi0_lwing,chi0_uwing)

       call assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,Ep%npwepG0,rhotwg,Gsph_epsG0,chi0)
     end do !band

   end do !ik_bz
 end do !spin
 !
 ! Collect body, heads and wings within comm
 comm=Wfd%comm
 do io=1,Ep%nomega
   call xsum_mpi(chi0(:,:,io),comm,ierr)
 end do
 !
 call xsum_mpi(chi0_head,comm,ierr)
 call xsum_mpi(chi0_lwing,comm,ierr)
 call xsum_mpi(chi0_uwing,comm,ierr)
 !
 ! Divide by the volume
 chi0       = chi0       * weight/Cryst%ucvol
 chi0_head  = chi0_head  * weight/Cryst%ucvol
 do io=1,Ep%nomega ! Tensor in the basis of the reciprocal lattice vectors.
   chi0_head(:,:,io) = MATMUL(chi0_head(:,:,io),Cryst%gmet) * (two_pi**2)
 end do
 chi0_lwing = chi0_lwing * weight/Cryst%ucvol
 chi0_uwing = chi0_uwing * weight/Cryst%ucvol
 !
 ! ===============================================
 ! ==== Symmetrize chi0 in case of AFM system ====
 ! ===============================================
 ! * Reconstruct $chi0{\down,\down}$ from $chi0{\up,\up}$.
 ! * Works only in the case of magnetic group Shubnikov type IV.
 if (Cryst%use_antiferro) then
   call symmetrize_afm_chi0(Cryst,Gsph_epsG0,Ltg_q,Ep%npwe,Ep%nomega,chi0,chi0_head,chi0_lwing,chi0_uwing)
 end if
 !
 ! ===================================================
 ! ==== Construct heads and wings from the tensor ====
 ! ===================================================
 !do io=1,Ep%nomega    
 !  do ig=2,Ep%npwe
 !    wng = chi0_uwing(ig,io,:)
 !    chi0(1,ig,io) = vdotw(Ep%qlwl(:,1),wng,Cryst%gmet,"G") 
 !    wng = chi0_lwing(ig,io,:)
 !    chi0(ig,1,io) = vdotw(Ep%qlwl(:,1),wng,Cryst%gmet,"G") 
 !  end do
 !  chq = MATMUL(chi0_head(:,:,io), Ep%qlwl(:,1))
 !  chi0(1,1,io) = vdotw(Ep%qlwl(:,1),chq,Cryst%gmet,"G")  ! Use user-defined small q
 !end do

 call wfd_barrier(Wfd)

 ! Impose Hermiticity (valid only for zero or purely imaginary frequencies)
 ! MG what about metals, where we have poles around zero?
 !if (dtset%gw_eet/=-1) then
 !  do io=1,Ep%nomega
 !    if (ABS(REAL(Ep%omega(io)))<0.00001) then
 !      do ig2=1,Ep%npwe
 !        do ig1=1,ig2-1
 !         chi0(ig2,ig1,io)=CONJG(chi0(ig1,ig2,io))
 !        end do
 !      end do
 !    end if
 !  end do
 !end if

 do iomega=1,MIN(Ep%nomega,NOMEGA_PRINTED)
   write(msg,'(1x,a,i4,a,2f9.4,a)')' chi0_intra(G,G'') at the ',iomega,' th omega',Ep%omega(iomega)*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
   call print_arr(chi0(:,:,iomega),unit=std_out)
 end do
 !
 ! =====================
 ! ==== Free memory ====
 ! =====================
 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(tabr_k)
 ABI_DEALLOCATE(ur1)
 ABI_DEALLOCATE(gw_gbound)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ktabr)
 ABI_DEALLOCATE(igffteps0)

 if (Psps%usepaw==1) then ! deallocation for PAW.
   call cprj_free(Cprj1_bz )
   ABI_DEALLOCATE(Cprj1_bz)
   call cprj_free(Cprj1_ibz)
   ABI_DEALLOCATE(Cprj1_ibz)
   call destroy_paw_pwij(Pwij)
   ABI_DEALLOCATE(Pwij)
 end if

 call destroy_little_group(Ltg_q)
 call destroy_bz_mesh_type(Kmesh)

 DBG_EXIT("COLL")

end subroutine chi0q0_intraband
!!***


