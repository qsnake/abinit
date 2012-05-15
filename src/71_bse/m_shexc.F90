!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_shexc
!! NAME
!! m_shexc
!!
!! FUNCTION
!!  This module provides routines to read the Bethe-Salpeter Hamiltonian from file
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_shexc

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_bs_defs
 use m_shirley
 use m_errors 
 use m_xmpi
 use m_blas

 use m_gwdefs,            only : czero_gw, GW_Q0_DEFAULT
 use m_crystal,           only : crystal_structure
 use m_gsphere,           only : gvectors_type, gsph_fft_tabs
 use m_ebands,            only : bstruct_clean, bstruct_init, pack_eneocc, bst_plot_bands, get_dos
 use m_vcoul,             only : vcoul_t, vcoul_init, vcoul_free, vcoul_nullify
 use m_bz_mesh,           only : bz_mesh_type, get_BZ_item, make_path, make_mesh, destroy_bz_mesh_type, print_BZ_mesh, find_qmesh,&
&                                nullify_bz_mesh
 !use m_paw_pwij,          only : paw_pwff_type, paw_pwij_type, init_paw_pwij, destroy_paw_pwij, paw_rho_tw_g
 use m_wfs,               only : wfs_descriptor, wfd_get_ur, wfd_get_cprj, wfd_change_ngfft, wfd_ihave_ur, &
&                                wfd_distribute_bbp, wfd_destroy, wfd_nullify
 use m_commutator_vkbr,   only : kb_potential, nullify_kb_potential, destroy_kb_potential, init_kb_potential, nc_ihr_comm
 !use m_oscillators,       only : rho_tw_g
 use m_screen

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_shexc/shexc_t
!! NAME
!! shexc_t
!! 
!! FUNCTION
!! 
!! SOURCE

 type,public :: shexc_t

   integer :: nkbz
   ! Number of k-points in the BZ.

   integer :: nqbz
   ! Number of q-points in the BZ. q = k-k'

   integer :: nqopt
   ! Number of q-points for the optical spectra.

   integer :: nsppol
   ! Number of independent spin polarizations.

   integer :: nsh
   ! Number of Shirley basis set elements used to expand the periodic part u_{nk}.

   integer :: neh_size
   ! Number of basis set elements used to expand the e-h function u_{nk} u_{n'k'}^*

   integer :: prtvol
   ! Verbosity level.

   integer :: comm
   ! MPI communicator.

   logical :: use_coulomb_term =.TRUE.
   logical :: use_exchange_term=.TRUE.

   type(BZ_mesh_type) :: Kmesh
   type(BZ_mesh_type) :: Qmesh
   type(vcoul_t) :: Vc

   integer,pointer :: nreh(:)   SET2NULL
   ! nreh(nsppol)
   ! Number of resonant electron-hole transitions for each spin.

   !integer,pointer :: vcks2t(:,:,:,:)   SET2NULL
   ! vcks2t(v,c,ik_bz,spin) gives the transition index associated to (v,c,kbz,spin)

   ! Create transition table vcks2t
   !allocate(Bsp%vcks2t(BSp%lomo:BSp%homo,BSp%lumo:BSp%humo,BSp%nkbz,nsppol)); Bsp%vcks2t = 0
   
   complex(gwpc),pointer :: sh_ihrc(:,:,:,:)  SET2NULL
   ! sh_ihrc(sh_size,sh_size,3,nkibz)
   ! <sh1| i[H,r] |sh2>

   real(dp),pointer :: qopt(:,:)
   ! qopt(3,nqopt)
   ! reduced coordinates of the q-points for the optical spectra.

   complex(gwpc),pointer :: sh2eh(:,:)  SET2NULL
   ! sh2eh(nsh**2,neh_size)
   ! <O_i* O_j| EH>

   complex(gwpc),pointer :: wij(:,:,:)  SET2NULL
   ! wij(neh_size,neh_size,nqbz))
   ! <EH_a| W(q) | EH_b>

   complex(gwpc),pointer :: vij(:,:)     SET2NULL
   ! vij(neh_size,neh_size))
   ! <EH_a| /tilde v(q=0) | EH_b>

   type(transition),pointer :: Trans(:,:) SET2NULL
   ! Trans(max_nreh,nsppol)

   type(ksintp_t),pointer :: KS_intp(:,:)  SET2NULL
   ! KS_intp(nkbz,nsppol))

 end type shexc_t

 public :: shexc_nullify
 public :: shexc_free
 public :: shexc_init
 public :: shexc_interp
 public :: shexc_solve

!----------------------------------------------------------------------
!!***

CONTAINS  !====================================================================

!----------------------------------------------------------------------

!!****f* m_shexc/shexc_nullify
!! NAME
!!  shexc_nullify
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_shexc
!!
!! CHILDREN
!!      haydock_herm_algo,shexc_interp,timein,wrtout,xbarrier_mpi,xgemv
!!      xsum_mpi
!!
!! SOURCE

subroutine shexc_nullify(SHexc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shexc_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(shexc_t),intent(inout) :: SHexc

!************************************************************************

 !@shexc_t
 ! integer
 nullify(SHexc%nreh)

 ! real
 nullify(SHexc%qopt)

 ! complex
 nullify(SHexc%sh_ihrc)
 nullify(SHexc%sh2eh)
 nullify(SHexc%wij  )
 nullify(SHexc%vij  )

 ! types
 call nullify_BZ_mesh(SHexc%Qmesh) 
 call nullify_BZ_mesh(SHexc%Kmesh) 
 call vcoul_nullify(SHexc%Vc)

 nullify(SHexc%Trans  )
 nullify(SHexc%KS_intp)

end subroutine shexc_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_shexc/shexc_free
!! NAME
!!  shexc_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      exc_interp_ham
!!
!! CHILDREN
!!      haydock_herm_algo,shexc_interp,timein,wrtout,xbarrier_mpi,xgemv
!!      xsum_mpi
!!
!! SOURCE

subroutine shexc_free(SHexc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shexc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(shexc_t),intent(inout) :: SHexc

!************************************************************************

 !@shexc_t

 ! integer
 if (associated(SHexc%nreh))     then
   ABI_DEALLOCATE(SHexc%nreh)
 end if

 ! real 
 if (associated(SHexc%qopt))     then
   ABI_DEALLOCATE(SHexc%qopt)
 end if
                          
 ! complex
 if (associated(SHexc%sh_ihrc))  then
   ABI_DEALLOCATE(SHexc%sh_ihrc)
 end if
 if (associated(SHexc%sh2eh))    then
   ABI_DEALLOCATE(SHexc%sh2eh)
 end if
 if (associated(SHexc%wij  ))    then
   ABI_DEALLOCATE(SHexc%wij)
 end if
 if (associated(SHexc%vij  ))    then
   ABI_DEALLOCATE(SHexc%vij)
 end if
                                                   
 ! types                 
 call destroy_bz_mesh_type(SHexc%Qmesh) 
 call destroy_bz_mesh_type(SHexc%Kmesh) 
 call vcoul_free(SHexc%Vc)

 if (associated(SHexc%Trans)) then
   ABI_DEALLOCATE(SHexc%Trans)
 end if

 if (associated(SHexc%KS_intp)) then 
   call ks_intp_free(SHexc%KS_intp)
   ABI_DEALLOCATE(SHexc%KS_intp)
 end if

end subroutine shexc_free
!!***

!----------------------------------------------------------------------

!!****f* m_shexc/shexc_init
!! NAME
!!  shexc_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      exc_interp_ham
!!
!! CHILDREN
!!      haydock_herm_algo,shexc_interp,timein,wrtout,xbarrier_mpi,xgemv
!!      xsum_mpi
!!
!! SOURCE

subroutine shexc_init(SHexc,eh_coverage,Dtset,Wsh,Cryst,Gsph_Max,Gsph_c,Psps,Pawtab,KS_pawrhoij,KS_Paw_ij,Pawang,Pawrad,Pawfgr,&
&  kptopt,kptrlatt,nshiftk,shiftk,nspinor,nsppol,nspden,prtvol,use_coulomb_term,use_exchange_term,bidx,qopt,&
&  nfftf,ngfftf,ngfftc,ks_aerhor,ks_vtrial)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shexc_init'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,nshiftk,nspden,nfftf,prtvol
 integer,intent(in) :: nsppol,nspinor
 real(dp),intent(in) :: eh_coverage
 logical,intent(in) :: use_coulomb_term,use_exchange_term
 type(dataset_type),intent(in) :: Dtset
 type(crystal_structure),intent(in) :: Cryst
 type(gvectors_type),intent(in) :: Gsph_Max,Gsph_c
 type(Pseudopotential_type),intent(in) :: Psps
 type(pawang_type),intent(in) :: Pawang
 type(Pawfgr_type),intent(in) :: Pawfgr
 type(shexc_t),intent(inout) :: SHexc
 type(wfs_descriptor),intent(inout) :: Wsh
!array
 integer,intent(in) :: ngfftc(18),ngfftf(18),bidx(:)
 integer,intent(inout) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk),qopt(:,:)
 real(dp),intent(in) :: ks_vtrial(nfftf,nspden)
 real(dp),intent(in) :: ks_aerhor(nfftf,nspden)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawrhoij_type),intent(in) :: KS_Pawrhoij(Cryst%natom)
 type(Paw_ij_type),intent(in) :: KS_paw_ij(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: s1=1,k1=1,nomega1=1,nqlwl1=1,istwfk1=1
 integer :: intp_mband,intp_nk,eh1,eh2,sh1,sh2,ii
 integer :: spline_opt,nkbz,nkibz,nqbz,comm,min_bsize,istat
 integer :: iq_bz,neh_size,npwc,npwx,spin,lomo,homo,lumo,humo,ik_ibz
 integer :: inclvkb,npw_k,usepaw,sh_size !nband_k,
 real(dp) :: eps_inf,ir_cut,uv_cut
 character(len=500) :: msg
 type(kb_potential) :: KBgrad_k
 type(wfs_descriptor) :: Weh
!array
 integer,pointer :: kg_k(:,:)
 integer,allocatable :: intp_nband(:,:)
 real(dp) :: kpoint(3)
 real(dp),pointer :: intp_kpt(:,:)
 real(dp),allocatable :: vc_qlwl(:,:)
 real(dp),allocatable :: intp_ene(:,:,:),occ(:,:,:)
 !complex(dpc),pointer :: obloch(:,:) 
 complex(gwpc) :: ihrc(3,nspinor**2)
 complex(gwpc),pointer :: ug_c(:),ug_v(:) 
 complex(dpc),allocatable :: ir_kibz(:,:,:,:,:) !eps_exc(:,:),
 complex(gwpc),allocatable :: w_qbz(:,:,:),wphi(:),vphi(:)
 complex(gwpc),allocatable :: ehg1(:),ehg2(:)
 complex(gwpc),allocatable :: vc_q0(:)
 complex(dpc),allocatable :: intp_gwene(:,:,:)
 !complex(dpc),allocatable :: opt_cvk_bz(:,:,:,:,:) !,hmat(:,:)

!************************************************************************

 MSG_WARNING("Entering shexc_init")

 !@shexc_t
 call shexc_nullify(SHexc)

 comm       = Wsh%comm
 SHexc%comm = comm

 ! Init basic parameters.
 lomo=bidx(1)
 homo=bidx(2)
 lumo=bidx(3) 
 humo=bidx(4) 

 write(std_out,*)" lomo, homo, lumo, humo ",lomo, homo, lumo, humo
 ABI_CHECK(lomo==1,"lomo /= 1 not coded")
 ABI_CHECK(lumo==homo+1,"lumo /= homo+1 not coded")

 SHexc%use_coulomb_term  = use_coulomb_term
 SHexc%use_exchange_term = use_exchange_term
 SHexc%prtvol            = prtvol
 SHexc%nsppol            = nsppol
 SHexc%nsh               = Wsh%nband(k1,s1)

 usepaw = Psps%usepaw
 !
 ! The dense Kmesh.
 call make_mesh(SHexc%Kmesh,Cryst,kptopt,kptrlatt,nshiftk,shiftk)
 !
 ! The dense Q-mesh.
 call find_qmesh(SHexc%Qmesh,Cryst,SHexc%Kmesh)

 nkbz  = SHexc%Kmesh%nbz
 nkibz = SHexc%Kmesh%nibz
 nqbz  = SHexc%Qmesh%nbz

 SHexc%nkbz = nkbz
 SHexc%nqbz = nqbz

 if (SHexc%prtvol>0) then
   call print_BZ_mesh(SHexc%Kmesh,header="Dense Kmesh",prtvol=prtvol)
   call print_BZ_mesh(SHexc%Qmesh,header="Dense Qmesh",prtvol=prtvol)
 end if
 !
 ! =================================================
 ! ==== Init Coulomb tables on the dense Q-mesh ====
 ! =================================================
 ABI_ALLOCATE(vc_qlwl,(3,nqlwl1))
 vc_qlwl(:,1)= GW_Q0_DEFAULT
 !                                                                                                                                              
 call vcoul_init(SHexc%Vc,Gsph_Max,SHexc%Qmesh,SHexc%Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,&
&  Gsph_Max%ng,nqlwl1,vc_qlwl,Cryst%rprimd,ngfftf,comm)

 ABI_DEALLOCATE(vc_qlwl)
 !
 ! =======================================================
 ! ==== Interpolate the band structure on SHexc%Kmesh ====
 ! =======================================================
 intp_nk = SHexc%Kmesh%nibz
 ABI_ALLOCATE(intp_kpt,(3,intp_nk))
 intp_kpt = SHexc%Kmesh%ibz

 ABI_ALLOCATE(intp_nband,(intp_nk,nsppol))
 intp_nband=humo

 intp_mband = -HUGE(0)
 do spin=1,nsppol
   do ii=1,intp_nk
     intp_mband = MAX(intp_mband, intp_nband(ii,spin))
   end do
 end do

 ABI_ALLOCATE(intp_ene,(intp_mband,intp_nk,nsppol))

 ABI_ALLOCATE(SHexc%KS_intp,(intp_nk,nsppol))

 min_bsize=intp_mband; spline_opt=0

 call shirley_interp(Wsh,"Vectors",Dtset,Cryst,Psps,Pawtab,Pawfgr,Pawang,Pawrad,&
&  KS_pawrhoij,KS_Paw_ij,ngfftc,ngfftf,nfftf,ks_vtrial,spline_opt,&
&  intp_nband,intp_mband,intp_nk,intp_kpt,intp_ene,SHexc%KS_intp) 

 ! ====================================================================
 ! ==== Init the transition tables using the interpolated energies ====
 ! ====================================================================
 ir_cut = -HUGE(one)
 uv_cut = +HUGE(one)

 ABI_ALLOCATE(SHexc%nreh,(nsppol))
 SHexc%nreh=0

 ABI_ALLOCATE(intp_gwene,(intp_mband,intp_nk,nsppol))
 intp_gwene = intp_ene

 ABI_ALLOCATE(occ,(intp_mband,intp_nk,nsppol))
 occ=zero
 occ(1:homo,:,:)=two ! FIXME

 call init_transitions(SHexc%Trans,lomo,humo,ir_cut,uv_cut,nkbz,intp_mband,intp_nk,&
&                      nsppol,nspinor,intp_gwene,occ,SHexc%Kmesh%tab,SHexc%nreh) 

 ABI_DEALLOCATE(intp_gwene)
 ABI_DEALLOCATE(occ)

 do spin=1,nsppol
   write(msg,'(a,i2,a,i0)')" For spin: ",spin,' the number of resonant e-h transitions is: ',SHexc%nreh(spin)
   call wrtout(std_out,msg,"COLL")
 end do
                                                                                                           
 if (ANY(SHexc%nreh/=SHexc%nreh(1))) then
   write(msg,'(a,(i0))')" BSE code does not support ifferent number of transitions for the two spin channels",SHexc%nreh
   MSG_ERROR(msg)
 end if
 !
 ! Create transition table vcks2t
 !allocate(Bsp%vcks2t(BSp%lomo:BSp%homo,BSp%lumo:BSp%humo,BSp%nkbz,Dtset%nsppol)); Bsp%vcks2t = 0
 !do spin=1,Dtset%nsppol
 !  do it=1,BSp%nreh(spin)
 !    BSp%vcks2t(BSp%Trans(it,spin)%v,BSp%Trans(it,spin)%c,BSp%Trans(it,spin)%k,spin) = it
 !  end do 
 !end do
 !
 ! ===================================================================
 ! ==== Calculate <Sh1| e^{-ikr} i [H,r] e^{ikr} |Sh2> in the IBZ ====
 ! ===================================================================
 SHexc%nqopt=SIZE(qopt,DIM=2)
 ABI_ALLOCATE(SHexc%qopt,(3,SHexc%nqopt))
 SHexc%qopt = qopt

 ABI_ALLOCATE(ir_kibz,(3,lomo:humo,lomo:humo,intp_nk,nsppol))
 ir_kibz=czero

 sh_size = Wsh%nband(k1,s1)

 inclvkb = 0 !BSp%inclvkb ! FIXME
 call nullify_kb_potential(KBgrad_k)

 ABI_ALLOCATE(SHexc%sh_ihrc,(sh_size,sh_size,3,nkibz))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory sh_ihrc")
 !
 do ik_ibz=1,nkibz
   npw_k  =  Wsh%npwarr(k1)
   kg_k   => Wsh%Kdata(k1)%kg_k

   kpoint =  SHexc%Kmesh%ibz(:,ik_ibz) 

   if (inclvkb/=0.and.usepaw==0) then ! Prepare term i <n,k|[Vnl,r]|n"k>
     call init_kb_potential(KBgrad_k,Cryst,Psps,inclvkb,istwfk1,npw_k,kpoint,kg_k)
   end if
   !
   do sh2=1,sh_size
     ug_v => Wsh%Wave(sh2,k1,s1)%ug
     !%if (usepaw==1) call wfd_get_cprj(Wsh,sh2,k1,s1,Cryst,Cp_v,sorted=.FALSE.)
     do sh1=1,sh_size
       ug_c => Wsh%Wave(sh1,k1,s1)%ug

       if (usepaw==0) then  ! Calculate matrix elements of i[H,r] for NC pseudopotentials.        
         ihrc = nc_ihr_comm(nspinor,npw_k,istwfk1,inclvkb,kpoint,KBgrad_k,ug_c,ug_v,kg_k) 

       else ! Matrix elements of i[H,r] for PAW.
         !%call wfd_get_cprj(Wsh,sh1,k1,s1,Cryst,Cp_c,sorted=.FALSE.)
         !%ihrc = paw_ihr_comm(spin,nspinor,npw_k,istwf_k,kpoint,Cryst,Pawtab,ug_c,ug_v,kg_k,Cp_c,Cp_v,HUr)
       end if

       SHexc%sh_ihrc(sh1,sh2,:,ik_ibz) = ihrc(:,1)
     end do
   end do

   call destroy_kb_potential(KBgrad_k)
 end do ! ik_ibz

 call destroy_kb_potential(KBgrad_k)
 !
 ! ========================================================
 ! === Setup of the basis set for the eh representation ===
 ! ========================================================
 min_bsize = 0
 call wfd_shirley_to_eh(Wsh,Cryst,Psps,Pawtab,Pawang,Pawrad,min_bsize,eh_coverage,Weh,SHexc%sh2eh)
 !
 ! TODO: Project eh transitions onto the Weh basis set.
 !
 ! Basic dimensions.
 neh_size = Weh%nband(k1,s1)
 npwc     = Gsph_c%ng
 npwx     = MIN(Weh%npwarr(k1),SHexc%Vc%ng)

 SHexc%neh_size = neh_size
 !
 ! =================================================
 ! ==== Calculate W in the optimal eh basis set ====
 ! =================================================
 if (SHexc%use_coulomb_term) then
   !
   ABI_ALLOCATE(SHexc%wij,(neh_size,neh_size,nqbz))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out-of-memory in wij")

   ABI_ALLOCATE(w_qbz,(npwc,npwc,nomega1))
   ABI_ALLOCATE(wphi,(npwc))
   ABI_ALLOCATE(ehg1,(npwc))
   ABI_ALLOCATE(ehg2,(npwc))

   eps_inf = 12.0

   do iq_bz=1,nqbz
     call screen_mdielf(iq_bz,npwc,nomega1,1,eps_inf,Cryst,SHexc%Qmesh,SHexc%Vc,Gsph_c,nspden,nfftf,ngfftf,ks_aerhor,"W",w_qbz,comm)
     !
     do eh2=1,neh_size
       ehg2 = Weh%Wave(eh2,k1,s1)%ug(1:npwc)
       wphi = MATMUL(w_qbz(:,:,1), ehg2)
       do eh1=1,neh_size
         ehg1 = Weh%Wave(eh1,k1,s1)%ug(1:npwc)
         SHexc%wij(eh1,eh2,iq_bz) = DOT_PRODUCT(ehg1, wphi)
       end do
     end do
     !
   end do

   ABI_DEALLOCATE(ehg1)
   ABI_DEALLOCATE(ehg2)
   ABI_DEALLOCATE(w_qbz)
   ABI_DEALLOCATE(wphi)
 end if
 !
 ! =================================================
 ! ==== Calculate v in the optimal eh basis set ====
 ! =================================================
 if (SHexc%use_exchange_term) then
   !
   ABI_ALLOCATE(SHexc%vij,(neh_size,neh_size))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out-of-memory in vij")

   ABI_ALLOCATE(vc_q0,(npwx))
   vc_q0 = SHexc%Vc%vcqlwl_sqrt(1:npwx,1)**2 

   ABI_ALLOCATE(ehg1,(npwx))
   ABI_ALLOCATE(ehg2,(npwx))
   ABI_ALLOCATE(vphi,(npwx))

   do eh2=1,neh_size
     ehg2 = Weh%Wave(eh2,k1,s1)%ug(1:npwx)
     vphi = vc_q0 * ehg2
     vphi(1) = czero_gw
     do eh1=1,neh_size
       ehg1 = Weh%Wave(eh1,k1,s1)%ug(1:npwx)
       SHexc%vij(eh1,eh2) = DOT_PRODUCT(ehg1, vphi)
     end do
   end do
   !
   ABI_DEALLOCATE(vc_q0)
   ABI_DEALLOCATE(vphi)
   ABI_DEALLOCATE(ehg1)
   ABI_DEALLOCATE(ehg2)
 end if
 !
 ! Free memory.
 call wfd_destroy(Weh)
 ABI_DEALLOCATE(intp_ene)

end subroutine shexc_init
!!***

!----------------------------------------------------------------------

!!****f* m_shexc/shexc_get_optmel
!! NAME
!!  shexc_get_optmel
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      haydock_herm_algo,shexc_interp,timein,wrtout,xbarrier_mpi,xgemv
!!      xsum_mpi
!!
!! SOURCE

subroutine shexc_get_optmel(SHexc,Cryst,hsize,kets)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shexc_get_optmel'
 use interfaces_32_util
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize
 type(shexc_t),intent(in) :: SHexc
 type(crystal_structure),intent(in) :: Cryst
!arrays
 complex(dpc),intent(out) :: kets(hsize,SHexc%nqopt)

!Local variables ------------------------------
!scalars
 integer :: nband_k,dir,ib_v,ib_c,lomo,humo,spin,ik_bz,ik_ibz,isym_k,itim_k
 integer :: iqo,nkibz,nsppol,nkbz
 real(dp) :: ediff
 complex(dpc) :: emcvk !ph_mkt,
!arrays
 real(dp) :: qopt(3),kbz(3),mat_dp(3,3),qrot(3),b1_vec(3),b2_vec(3),b3_vec(3) !,kpoint(3),qpoint(3)
 real(dp),pointer :: intp_ene(:)
 complex(dpc),pointer :: obloch(:,:) 
 complex(dpc),allocatable :: ks_ihrc(:,:,:)
 complex(dpc),allocatable :: opt_cvk_bz(:,:,:,:,:) 
 complex(dpc),allocatable :: ir_kibz(:,:,:,:,:) !eps_exc(:,:),

!************************************************************************

 !@shexc_t
 MSG_ERROR("Not implemented error")
 !% lomo   = SHexc%lomo
 !% humo   = SHexc%humo
 nkibz  = SHexc%Kmesh%nibz
 nkbz   = SHexc%Kmesh%nbz
 nsppol = SHexc%nsppol

 kets = czero

 ABI_ALLOCATE(ir_kibz,(3,lomo:humo,lomo:humo,nkibz,nsppol))
 ir_kibz=czero

 do spin=1,nsppol
   do ik_ibz=1,nkibz
     nband_k  =  SHexc%KS_intp(ik_ibz,spin)%nband_k
     obloch   => SHexc%KS_intp(ik_ibz,spin)%obloch
     intp_ene => SHexc%KS_intp(ik_ibz,spin)%ene

     ABI_ALLOCATE(ks_ihrc,(nband_k,nband_k,3))
     do dir=1,3
       ks_ihrc(:,:,dir) = MATMUL( TRANSPOSE(CONJG(obloch)),MATMUL(SHexc%sh_ihrc(:,:,dir,ik_ibz),obloch) )
     end do
     !
     ! Save matrix elements of ir in the IBZ.
     do ib_v=lomo,humo !  Loops over the bands C and V start
       do ib_c=lomo,humo
         ediff = intp_ene(ib_c) - intp_ene(ib_v)
         if (ABS(ediff)<tol16) ediff=tol6  ! Treat a possible degeneracy between v and c.
         ir_kibz(:,ib_c,ib_v,ik_ibz,spin) = ks_ihrc(ib_c,ib_v,:) / ediff
       end do
     end do
     ABI_DEALLOCATE(ks_ihrc)
   end do
 end do
 !
 ! ======================================================
 ! ==== Calculate Fcv(kBZ) in the full Brilouin zone ====
 ! ======================================================
 !
 ! Symmetrization of the matrix elements of the position operator.
 ! <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
 !   where S is one of the symrec operations in reciprocal space, R is the 
 !   corresponding operation in real space, \tau being the associated fractional translations.
 !
 ! q.Mcv( Sk) =  S^{-1}q. Mcv(k) 
 ! q.Mcv(-Sk) = -S^{-1}q. CONJG(Mcv(k)) if time-reversal is used.

 ABI_ALLOCATE(opt_cvk_bz,(lomo:humo,lomo:humo,nkbz,nsppol,SHexc%nqopt))
 opt_cvk_bz=czero
                                                                                               
 b1_vec= two_pi * Cryst%gprimd(:,1)
 b2_vec= two_pi * Cryst%gprimd(:,2)
 b3_vec= two_pi * Cryst%gprimd(:,3)

 do spin=1,nsppol
   !
   do ik_bz=1,nkbz
    !
    ! * Get ik_ibz, and symmetries index from ik_bz.
    call get_BZ_item(SHexc%Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)
    
    mat_dp = DBLE(Cryst%symrec(:,:,isym_k))
    call matrginv(mat_dp,3,3) ! Invert
                                                                                               
    do iqo=1,SHexc%nqopt
      qopt = SHexc%qopt(:,iqo)
      qrot = (3-2*itim_k) * MATMUL(mat_dp,qopt)
                                                                                               
      do ib_v=lomo,humo !  Loops over the bands C and V start
        do ib_c=lomo,humo
          !if (ib_c==ib_v) CYCLE 
          emcvk = pdtqrc(qrot,ir_kibz(:,ib_c,ib_v,ik_ibz,spin),b1_vec,b2_vec,b3_vec)
          if (itim_k==2) emcvk = CONJG(emcvk)
          opt_cvk_bz(ib_c,ib_v,ik_bz,spin,iqo) = emcvk
        end do !ib_c
      end do !ib_v
    end do
                                                                                               
   end do !ik_bz
 end do !spin

 ABI_DEALLOCATE(ir_kibz)
 ABI_DEALLOCATE(opt_cvk_bz)

end subroutine shexc_get_optmel
!!***

!----------------------------------------------------------------------

!!****f* m_shexc/shexc_interp
!! NAME
!!  shexc_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine shexc_interp(SHexc,t1,t2,spin1,spin2,isreso,oelm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shexc_interp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: t1,t2,spin1,spin2
 complex(dpc),intent(out) :: oelm
 logical,intent(in) :: isreso
 type(shexc_t),intent(in) :: SHexc

!Local variables ------------------------------
!scalars
 integer :: ik_bz,ikp_bz,iv,ivp,ic,icp,idx
 integer :: nsh,nsh_sq,iq_bz,neh_size,eh,sh1,sh2
 complex(dpc) :: ene_t,ene_tp
!array
 !integer :: g0(3)
 real(dp) :: kmkp(3)
 complex(gwpc),allocatable :: cbuff(:,:),left(:),right(:)
 complex(dpc),pointer :: mat_kv(:),mat_kpvp(:),mat_kc(:),mat_kpcp(:)

!************************************************************************
 
 ABI_CHECK(isreso," Coupling not coded")

 ene_t   = SHexc%Trans(t1,spin1)%en
 ik_bz   = SHexc%Trans(t1,spin1)%k
 iv      = SHexc%Trans(t1,spin1)%v
 ic      = SHexc%Trans(t1,spin1)%c

 ene_tp = SHexc%Trans(t2,spin2)%en
 ikp_bz = SHexc%Trans(t2,spin2)%k
 ivp    = SHexc%Trans(t2,spin2)%v
 icp    = SHexc%Trans(t2,spin2)%c

 nsh      = SHexc%nsh
 nsh_sq   = nsh**2
 neh_size = SHexc%neh_size
                                          
 ABI_ALLOCATE(cbuff,(nsh_sq,nsh_sq))
 ABI_ALLOCATE(left,(neh_size))
 ABI_ALLOCATE(right,(neh_size))

 oelm = czero_gw
 !
 ! Diagonal term.
 if (isreso .and. t1==t2 .and. spin1==spin2) oelm = oelm + ene_t
 !
 ! Direct Coulomb term.
 if (.FALSE..and.SHexc%use_coulomb_term) then
   !
   ! * Find q = K-KP-G0 in the full BZ.
   kmkp = SHexc%Kmesh%bz(:,ik_bz) - SHexc%Kmesh%bz(:,ikp_bz)
   !%call findqg0(iq_bz,g0,kmkp,SHexc%Qmesh%nbz,SHexc%Qmesh%bz,BSp%mG0)
   iq_bz = 1

   !write(std_out,*)" SHAPE(SHexc%sh2eh) ",SHAPE(SHexc%sh2eh)
   !write(std_out,*)" SHAPE(SHexc%wij) ",SHAPE(SHexc%wij(:,:,iq_bz))

   cbuff = MATMUL( SHexc%sh2eh, MATMUL(SHexc%wij(:,:,iq_bz), CONJG(TRANSPOSE(SHexc%sh2eh))) )

   mat_kv   => SHexc%KS_intp(ik_bz ,spin1)%obloch(:,iv)
   mat_kc   => SHexc%KS_intp(ik_bz ,spin1)%obloch(:,ic)

   mat_kpvp => SHexc%KS_intp(ikp_bz,spin2)%obloch(:,ivp)
   mat_kpcp => SHexc%KS_intp(ikp_bz,spin2)%obloch(:,icp)

   left = czero_gw; right= czero_gw
   do eh=1,neh_size
    !
    idx=0
    do sh2=1,nsh
      do sh1=1,nsh
        idx = idx + 1
        !left(eh) = left(eh) + mat_kc(sh1) * mat_kpcp(sh2) * SHexc%sh2eh(idx,eh) 
        !right(eh) = right(eh) + mat_kc(sh1) * mat_kpcp(sh2) * SHexc%sh2eh(idx,eh) 
      end do
    end do
    !
   end do

   !right = MATMUL(cbuff, right)
   !oelm = oelm + DOT_PRODUCT(left, right)
 end if
 !
 ! Exchange term.
 if (.FALSE..and.SHexc%use_exchange_term) then
   cbuff = MATMUL( SHexc%sh2eh, MATMUL(SHexc%vij, CONJG(TRANSPOSE(SHexc%sh2eh))) )
 end if

 ABI_DEALLOCATE(left)
 ABI_DEALLOCATE(right)
 ABI_DEALLOCATE(cbuff)

end subroutine shexc_interp
!!***

!----------------------------------------------------------------------

!!****f* m_shexc/shexc_solve
!! NAME
!!  shexc_solve
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      exc_interp_ham
!!
!! CHILDREN
!!      haydock_herm_algo,shexc_interp,timein,wrtout,xbarrier_mpi,xgemv
!!      xsum_mpi
!!
!! SOURCE

subroutine shexc_solve(SHexc,Bsp,Cryst)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shexc_solve'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_71_bse
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(shexc_t),intent(inout) :: SHexc
 type(excparam),intent(in) :: Bsp
 type(crystal_structure),intent(in) :: Cryst

!Local variables ------------------------------
!scalars
 integer :: hsize,my_nt,my_t1,my_t2,istat,spin1,spin2,t1,t2,niter_file,niter_max
 integer :: nproc,my_rank,master,nsppol,term_type,comm,iq,nkets,ierr,niter_done
 integer :: inn
 complex(dpc) :: oelm,factor
 real(dp) :: cpu_time,wall_time,cpu0,wall0
 real(dp) :: norm,nfact
 logical :: isreso,is_converged
 character(len=500) :: msg
!arrays
 complex(dp) :: green(BSp%nomega,Bsp%nq)
 complex(dpc),allocatable :: hmat(:,:)
 !real(dp),pointer :: bb_file(:)
 real(dp),allocatable :: bb(:)
 complex(dpc),allocatable :: aa(:),phi_nm1(:),phi_n(:),hphi_n(:)
 !complex(dpc),pointer :: aa_file(:),phi_n_file(:),phi_nm1_file(:)
 complex(dpc),allocatable :: ket0(:)
 logical :: check(2)

!************************************************************************

 !@shexc_t
 !
 ! Interpolate the excitonic Hamiltonian.
 comm  = SHexc%comm
 hsize = SHexc%nreh(1)
 my_t1=1; my_t2=hsize
 spin1=1; spin2=1; isreso=.TRUE.

 call timein(cpu0,wall0)

 ABI_ALLOCATE(hmat,(hsize,my_t1:my_t2))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory hmat")
                                                                                                
 do t2=1,hsize
   do t1=1,hsize
     call shexc_interp(SHexc,t1,t2,spin1,spin2,isreso,oelm)
     hmat(t1,t2) = oelm
   end do
 end do

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0
 write(std_out,*)" Build of Hexc cpu_time, wall_time: ",cpu_time, wall_time

 nproc  = xcomm_size(comm)
 my_rank= xcomm_rank(comm)
 master = 0 
 nsppol = SHexc%nsppol

 my_nt = my_t2-my_t1+1
 ABI_CHECK(my_nt>0,"One of the processors has zero columns")

 write(msg,'(a,i0)')' Haydock algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")

 !
 ! Select the terminator for the continued fraction.
 term_type=0; if (Bsp%hayd_term>0) term_type=1 
 write(msg,'(a,i0)')" Using terminator type: ",term_type
 call wrtout(std_out,msg,"COLL")
 !
 ! Calculate green(w) for the different starting points.
 green=czero
 do iq=1,nkets
   !
   ABI_ALLOCATE(ket0,(hsize))
   !ket0=kets(:,iq)
   !
   niter_file=0
   !
   ! For n>1 we have
   ! 1) a_n = <n|H|n>
   ! 2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
   ! 3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
   !
   ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
   ! a_1 = <1|H|1>
   ! b_1 = || H|1> - a_1|1> ||
   ! |2> = [H|1> - a_1|1>]/b_1
   !
   ABI_ALLOCATE(hphi_n,(hsize))
   ABI_ALLOCATE(phi_nm1,(my_nt))
   ABI_ALLOCATE(phi_n,(my_nt))

   niter_max = niter_file + Bsp%niter
   ABI_ALLOCATE(aa,(niter_max))
   ABI_ALLOCATE(bb,(niter_max))
   aa=czero; bb=zero

   if (niter_file==0) then ! Calculation from scratch.
     phi_nm1=ket0(my_t1:my_t2)   ! Select the slice treated by this node.
     norm = xnrm2(hsize,ket0,1) ! Normalization  
     phi_nm1=phi_nm1/norm      
                                                                                
     ! hphi_n = MATMUL(hmat,phi_nm1)
     call xgemv('N',hsize,my_nt,cone,hmat,hsize,phi_nm1,1,czero,hphi_n,1)
     call xsum_mpi(hphi_n,comm,ierr)

     aa(1)=xdotc(my_nt,phi_nm1,1,hphi_n(my_t1:),1)
     call xsum_mpi(aa(1:1),comm,ierr)

     phi_n = hphi_n(my_t1:my_t2) - aa(1)*phi_nm1

     bb(1) = xdotc(my_nt,phi_n,1,phi_n,1)
     call xsum_mpi(bb(1:1),comm,ierr)
     bb(1) = SQRT(bb(1))

     phi_n = phi_n/bb(1)
     niter_done=1

   else ! Use the previously calculates a and b.
     MSG_ERROR("Not implemented error")
     !niter_done=niter_file
     !aa(1:niter_done) = aa_file
     !bb(1:niter_done) = bb_file
     !phi_nm1=phi_nm1_file(my_t1:my_t2)   ! Select the slice treated by this node.
     !phi_n  =phi_n_file  (my_t1:my_t2)   
   end if


   !if (associated(aa_file     )) deallocate(aa_file)
   !if (associated(bb_file     )) deallocate(bb_file)
   !if (associated(phi_nm1_file)) deallocate(phi_nm1_file)
   !if (associated(phi_n_file  )) deallocate(phi_n_file)

   ! Multiplicative factor (k-point sampling and unit cell volume)  
   ! TODO be careful with the spin here
   ! TODO four_pi comes from the coulomb term 1/|q| is already included in the 
   ! oscillators hence the present approach wont work if a cutoff interaction is used.
   nfact = -four_pi/(Cryst%ucvol*SHexc%Kmesh%nbz)
   if (nsppol==1) nfact=two*nfact 

   factor = nfact*(xnrm2(hsize,ket0,1)**2)

   check = (/.FALSE.,.TRUE./) ! TODO for the time being only the imaginary part is tested.
   ! TODO pass negative frequencies.

   call haydock_herm_algo(niter_done,niter_max,BSp%nomega,BSp%omega,BSp%haydock_tol(1),check,hsize,&
&    my_t1,my_t2,hmat,factor,term_type,aa,bb,phi_nm1,phi_n,green(:,iq),inn,is_converged,comm)
   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed 
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !
   !if (my_rank==master) then ! Open the file and writes basic dimensions and info.
   !  write(out_unt)Bsp%q(:,iq)
   !  write(out_unt)MIN(inn,niter_max)  ! NB if the previous loop completed inn=niter_max+1
   !  do it=1,MIN(inn,niter_max)        ! if we exited then inn is not incremented by one.
   !    write(out_unt)it,aa(it),bb(it)
   !  end do
   !end if
   !
   ! hphi_n is used as workspace to gather |n> and |n+1>.
   !
   !hphi_n = czero
   !hphi_n(my_t1:my_t2) = phi_nm1
   !call xsum_master(hphi_n,master,comm,ierr)
   !if (my_rank==master) write(out_unt)hphi_n ! |n-1>

   !hphi_n = czero
   !hphi_n(my_t1:my_t2) = phi_n
   !call xsum_master(hphi_n,master,comm,ierr)
   !if (my_rank==master) write(out_unt)hphi_n ! |n>

   ABI_DEALLOCATE(hphi_n)
   ABI_DEALLOCATE(phi_nm1)
   ABI_DEALLOCATE(phi_n)
   ABI_DEALLOCATE(aa)
   ABI_DEALLOCATE(bb)
   ABI_DEALLOCATE(ket0)

 end do ! iq

 ABI_DEALLOCATE(hmat)

 call xbarrier_mpi(comm)

end subroutine shexc_solve
!!***

!----------------------------------------------------------------------

END MODULE m_shexc
!!***
