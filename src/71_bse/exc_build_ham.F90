!{\src2tex{textfont=tt}}
!!****f* ABINIT/exc_build_ham
!! NAME
!!  exc_build_ham
!!
!! FUNCTION
!!  Calculate and write the excitonic Hamiltonian on an external binary file (Fortran file open
!!  in random mode) for subsequent treatment in the Bethe-Salpeter code.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2009 EXC group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!! Copyright (C) 2009-2012 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
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
!!  Vcp<vcoul_t>=The Coulomb interaction in reciprocal space. A cutoff can be used
!!  W<screen_t>=Data type gathering info and data for W.
!!  nfftot_osc=Total Number of FFT points used for the oscillator matrix elements.
!!  ngfft_osc(18)=Info on the FFT algorithm used to calculate the oscillator matrix elements.
!!  Psps<Pseudopotential_type>=Variables related to pseudopotentials
!!  Pawtab(Psps%ntypat)<pawtab_type>=PAW tabulated starting data.
!!  Pawang<pawang_type>=PAW angular mesh and related data.
!!  Paw_pwff(Cryst%ntypat*Wfd%usepaw)<Paw_pwff_type>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  Wfd<wfs_descriptor>=Handler for the wavefunctions.
!!
!! OUTPUT
!!  The excitonic Hamiltonian is saved on an external binary file (see below).
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,destroy_paw_pwij,exc_build_block,get_bz_item
!!      gsph_fft_tabs,init_paw_pwij,paw_rho_tw_g,rho_tw_g,timab
!!      wfd_change_ngfft,wfd_get_cprj,wfd_get_ur,wrtout,xmpi_distab,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine exc_build_ham(BSp,BS_files,Cryst,Kmesh,Qmesh,ktabr,Gsph_Max,Gsph_c,Vcp,&
& Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_bs_defs
 use m_xmpi
 use m_errors

 use m_gwdefs,       only : czero_gw
 use m_crystal,      only : crystal_structure
 use m_gsphere,      only : gvectors_type, gsph_fft_tabs
 use m_vcoul,        only : vcoul_t
 use m_bz_mesh,      only : bz_mesh_type, get_BZ_item 
 use m_paw_pwij,     only : paw_pwff_type, paw_pwij_type, init_paw_pwij, destroy_paw_pwij, paw_rho_tw_g
 use m_wfs,          only : wfs_descriptor, wfd_get_ur, wfd_get_cprj, wfd_change_ngfft, wfd_ihave_ur, wfd_distribute_bbp
 use m_oscillators,  only : rho_tw_g
 use m_screen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_build_ham'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_71_bse, except_this_one => exc_build_ham
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_osc
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(screen_t),intent(inout) :: W
 type(BZ_mesh_type),intent(in) :: Kmesh,Qmesh
 type(crystal_structure),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
 type(Gvectors_type),intent(in) :: Gsph_Max,Gsph_c
 type(Pseudopotential_type),intent(in) :: Psps
 type(Hdr_type),intent(inout) :: Hdr_bse
 type(pawang_type),intent(in) :: Pawang
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfft_osc(18)
 integer,intent(in) :: ktabr(nfftot_osc,BSp%nkbz)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(Paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0,map2sphere=1,dim_rtwg1=1
 integer :: nproc,nspinor,nsppol,comm
 integer :: use_padfft
 integer :: my_rank,mgfft_osc
 integer :: fftalga_osc
 integer :: ik_ibz,itim_k,isym_k
 integer :: iq_bz,iq_ibz,isym_q,itim_q,iqbz0 
 integer :: ierr,iv,ic,istat,spin
 !real(dp) :: q0vol,fcc_const
 complex(dpc) :: ph_mkt
 logical :: do_resonant,do_coupling
 character(len=fnlen) :: fname
! character(len=500) :: msg
!arrays
 integer :: got(Wfd%nproc) !,bbp_distrb(Wfd%mband,Wfd%mband)
 integer,allocatable :: igfftg0(:),task_distrib(:,:,:,:)
 integer,allocatable :: gbound(:,:)
 real(dp) :: qbz(3),spinrot_k(4),tsec(2)
 complex(gwpc),allocatable :: rhotwg1(:),rhxtwg_q0(:,:,:,:,:) 
 complex(gwpc),target,allocatable :: ur1(:),ur2(:)
 complex(gwpc),pointer :: ptr_ur1(:),ptr_ur2(:)
 logical :: bbp_mask(Wfd%mband,Wfd%mband)
 type(Cprj_type),allocatable :: Cp1(:,:),Cp2(:,:)
 type(Paw_pwij_type),allocatable :: Pwij_q0(:)

!************************************************************************

 call timab(670,1,tsec)
 call timab(671,1,tsec)

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(nfftot_osc==PRODUCT(ngfft_osc(1:3)),"mismatch in FFT size")

 if (BSp%have_complex_ene) then
   MSG_ERROR("Complex energies are not supported yet")
 end if

 my_rank = Wfd%my_rank
 nproc   = Wfd%nproc
 comm    = Wfd%comm

 nspinor = Wfd%nspinor
 nsppol  = Wfd%nsppol

 do_resonant = (BS_files%in_hreso == BSE_NOFILE)
 do_coupling = (BS_files%in_hcoup == BSE_NOFILE)

 if (.not. do_resonant .and. .not. do_coupling) then
   MSG_COMMENT("Skipping calculation of both resonant and coupling block")
   RETURN
 end if

 if ( ANY(ngfft_osc(1:3) /= Wfd%ngfft(1:3)) ) call wfd_change_ngfft(Wfd,Cryst,Psps,ngfft_osc) 
 mgfft_osc   = MAXVAL(ngfft_osc(1:3))
 fftalga_osc = ngfft_osc(7)/100 !; fftalgc_osc=MOD(ngfft_osc(7),10)

! Analytic integration of 4pi/q^2 over the volume element:
! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255
! (see gwa.pdf, appendix A.4)

! If q=0 and C=V then set up rho-twiddle(G=0) to reflect an
! analytic integration of q**-2 over the volume element:
! <q**-2> = 7.44 V**(-2/3)   (for fcc cell)

! q0vol = (8.0*pi**3) / (Cryst%ucvol*BSp%nkbz)
! fcc_const = SQRT(7.44*q0vol**(-2.0/3.0))
! rtw = (6.0*pi**2/(Cryst%ucvol*BSp%nkbz))**(1./3.)
! Average of (q+q')**-2 integration for head of Coulomb matrix 
! INTRTW(QL) = (2*pi*rtw + pi*(rtw**2/QL-QL)*LOG((QL+rtw)/(QL-rtw)))
! &              * (Cryst%ucvol*BSp%nkbz)/(2*pi)**3. * QL*QL

 if (Wfd%usepaw==1) then
   ABI_ALLOCATE(Cp1,(Wfd%natom,nspinor))
   call cprj_alloc(Cp1,0,Wfd%nlmn_atm)
   ABI_ALLOCATE(Cp2,(Wfd%natom,nspinor))
   call cprj_alloc(Cp2,0,Wfd%nlmn_atm)
 end if

 call wrtout(std_out," Calculating all matrix elements for q=0 to save CPU time","COLL")

 ABI_ALLOCATE(ur1,(nfftot_osc*nspinor))
 ABI_ALLOCATE(ur2,(nfftot_osc*nspinor))

 ! Identify q==0
 iqbz0=0
 do iq_bz=1,Qmesh%nbz
   if (ALL(ABS(Qmesh%bz(:,iq_bz))<tol3)) iqbz0 = iq_bz
 end do
 ABI_CHECK(iqbz0/=0,"q=0 not found in q-point list!")

 ! * Get iq_ibz, and symmetries from iqbz0.
 call get_BZ_item(Qmesh,iqbz0,qbz,iq_ibz,isym_q,itim_q)

 if (Wfd%usepaw==1) then ! Prepare onsite contributions at q==0
   ABI_ALLOCATE(Pwij_q0,(Cryst%ntypat))
   call init_paw_pwij(Pwij_q0,BSp%npweps,Qmesh%bz(:,iqbz0),Gsph_Max%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
 end if
 !
 ! Tables for the FFT of the oscillators.
 !  a) FFT index of the G sphere (only vertical transitions, unlike cchi0, no need to shift the sphere).
 !  b) gbound table for the zero-padded FFT performed in rhotwg. 
 ABI_ALLOCATE(igfftg0,(Gsph_Max%ng))
 ABI_ALLOCATE(gbound,(2*mgfft_osc+8,2))
 call gsph_fft_tabs(Gsph_Max,(/0,0,0/),mgfft_osc,ngfft_osc,use_padfft,gbound,igfftg0)
 if ( ANY(fftalga_osc == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
 if (use_padfft==0) then 
   ABI_DEALLOCATE(gbound)
   ABI_ALLOCATE(gbound,(2*mgfft_osc+8,2*use_padfft))
 end if

 ABI_ALLOCATE(rhotwg1,(BSp%npweps))

 ABI_ALLOCATE(rhxtwg_q0,(BSp%npweps,BSp%lomo:BSp%nbnds,BSp%lomo:BSp%nbnds,Wfd%nkibz,nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out-of-memory in rhxtwg_q0")
 rhxtwg_q0 = czero

 ! Distribute the calculation of the matrix elements.
 ABI_ALLOCATE(task_distrib,(BSp%lomo:BSp%nbnds,BSp%lomo:BSp%nbnds,Wfd%nkibz,nsppol))
 call xmpi_distab(nproc,task_distrib)

 got=0
 bbp_mask=.FALSE.; bbp_mask(Bsp%lomo:Bsp%nbnds,Bsp%lomo:Bsp%nbnds)=.TRUE.

 do spin=1,nsppol
   do ik_ibz=1,Wfd%nkibz ! loop over the k-points in IBZ

    ! Distribute the (b,b') entries.
    !%call wfd_distribute_bbp(Wfd,ik_ibz,spin,"All",my_nbbp,bbp_distrb,got=got,bbp_mask=bbp_mask) 
    !%if ( ALL(bbp_distrb/=my_rank) ) CYCLE

    if ( ALL(task_distrib(:,:,ik_ibz,spin)/= my_rank) ) CYCLE

    itim_k=1; isym_k=1; ph_mkt=cone; spinrot_k=Cryst%spinrot(:,isym_k)

    do iv=BSp%lomo,BSp%nbnds ! Loop over band V
      if ( ALL(task_distrib(:,iv,ik_ibz,spin)/=my_rank) ) CYCLE
      !if ( ALL(bbp_distrb(:,iv)/=my_rank) ) CYCLE

      if (wfd_ihave_ur(Wfd,iv,ik_ibz,spin,how="Stored")) then
        ptr_ur1 =>  Wfd%Wave(iv,ik_ibz,spin)%ur
      else
        call wfd_get_ur(Wfd,iv,ik_ibz,spin,ur1)
        ptr_ur1 =>  ur1 
      end if

      if (Wfd%usepaw==1) then 
        call wfd_get_cprj(Wfd,iv,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
      end if

      do ic=BSp%lomo,BSp%nbnds ! Loop over band C 
        if ( task_distrib(ic,iv,ik_ibz,spin)/=my_rank ) CYCLE
        !if ( ALL(bbp_distrb(ic,iv)/=my_rank) ) CYCLE

        if (wfd_ihave_ur(Wfd,ic,ik_ibz,spin,how="Stored")) then
          ptr_ur2 =>  Wfd%Wave(ic,ik_ibz,spin)%ur
        else
          call wfd_get_ur(Wfd,ic,ik_ibz,spin,ur2)
          ptr_ur2 =>  ur2
        end if

        if (Wfd%usepaw==1) then 
          call wfd_get_cprj(Wfd,ic,ik_ibz,spin,Cryst,Cp2,sorted=.FALSE.)
        end if

        call rho_tw_g(Wfd%paral_kgb,nspinor,BSp%npweps,nfftot_osc,ngfft_osc,map2sphere,use_padfft,igfftg0,gbound,&
&         ptr_ur1,1,ktabr(:,1),ph_mkt,spinrot_k,&
&         ptr_ur2,1,ktabr(:,1),ph_mkt,spinrot_k,&
&         dim_rtwg1,rhotwg1,tim_fourdp0,Wfd%MPI_enreg)

        if (Wfd%usepaw==1) then ! Add PAW onsite contribution.
          call paw_rho_tw_g(Bsp%npweps,dim_rtwg1,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_Max%gvec,&
&           Cp1,Cp2,Pwij_q0,rhotwg1)
        end if

        ! If q=0 treat Exchange and Coulomb-term independently
        if (iv <= BSp%homo .and. ic <= BSp%homo .or. &
&           iv >  BSp%homo .and. ic >  BSp%homo) then

          if (iv/=ic) then !COULOMB term: C/=V: ignore them
            rhotwg1(1) = czero_gw
          else
            ! If q=0 and C=V then set up rho-twiddle(G=0) to reflect an
            ! analytic integration of q**-2 over the volume element:
            ! <q**-2> = 7.44 V**(-2/3)   (for fcc cell)
            !rhotwg1(1) = fcc_const * qpg(1,iqbz0)
            rhotwg1(1) = SQRT(Vcp%i_sz) / Vcp%vcqlwl_sqrt(1,1)
            !if (vcut) rhotwg1(1) = 1.0
          end if

        else
!         At present this term is set to zero
!         EXCHANGE term: limit value. 
!         Set up rho-twiddle(G=0) using small vector q instead of zero and k.BSp perturbation theory (see notes)
          rhotwg1(1) = czero_gw
        end if

        rhxtwg_q0(:,iv,ic,ik_ibz,spin) = rhotwg1(:)
      end do ! ic
    end do ! iv
   end do ! ik_ibz
 end do ! spin
 !
 ! Gather matrix elements on each node.
 call xsum_mpi(rhxtwg_q0,comm,ierr)

 ABI_DEALLOCATE(task_distrib)
 ABI_DEALLOCATE(rhotwg1)
 ABI_DEALLOCATE(igfftg0)
 ABI_DEALLOCATE(gbound)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ur1)
 ABI_DEALLOCATE(ur2)

 if (Wfd%usepaw==1) then ! Optional deallocation for PAW.
   call destroy_paw_pwij(Pwij_q0)
   ABI_DEALLOCATE(Pwij_q0)
   call cprj_free(Cp1)
   ABI_DEALLOCATE(Cp1)
   call cprj_free(Cp2)
   ABI_DEALLOCATE(Cp2)
 end if

 call timab(671,2,tsec)
 call timab(672,1,tsec)

 !
 ! ========================
 ! ==== Resonant Block ====
 ! ========================
 if (do_resonant) then
   fname=BS_files%out_hreso
   call exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_Max,Gsph_c,Vcp,&
&    Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,rhxtwg_q0,.TRUE.,fname)
 end if 

 call timab(672,2,tsec)
 call timab(673,1,tsec)

 !
 ! ========================
 ! ==== Coupling Block ====
 ! ========================
 if (do_coupling.and.BSp%use_coupling>0) then
   fname=BS_files%out_hcoup
   call exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_Max,Gsph_c,Vcp,&
&    Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,rhxtwg_q0,.FALSE.,fname)
 end if 
 !
 ! * Free memory.
 ABI_DEALLOCATE(rhxtwg_q0) 

 call timab(673,2,tsec)
 call timab(670,2,tsec)

end subroutine exc_build_ham
!!***
