!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_shirley
!! NAME
!!  m_shirley
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

MODULE m_shirley

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_cprj_bspline

 use m_fstrings,       only : starts_with
 use m_io_tools,       only : flush_unit, get_unit
 use m_numeric_tools,  only : print_arr, hermitianize, imax_loc, dst_t, dst_init
 use m_blas,           only : xdotc, xgemm, blas_cholesky_ortho
 use m_abilasi,        only : xheev, xhegv, xheevx, xhegvx
 use m_fft_mesh,       only : fft_check_rotrans, setmesh
 use m_geometry,       only : normv, vdotw 
 use m_gsphere,        only : get_kg
 use m_crystal,        only : crystal_structure
 use m_bz_mesh,        only : bz_mesh_type, get_BZ_item
 use m_ebands,         only : update_occ, bst_plot_bands, get_eneocc_vect, pack_eneocc, bstruct_init
 use m_hamiltonian,    only : ddiago_ctl_type, init_ddiago_ctl, destroy_hamiltonian, init_hamiltonian, finalize_hamiltonian
 use m_paw_pwij,       only : paw_pwff_type, paw_pwff_type, init_paw_pwff, destroy_paw_pwff, paw_pwij_type, init_paw_pwij, & 
&                             destroy_paw_pwij, paw_rho_tw_g
 use m_wfs,            only : wfd_get_ur, wfs_descriptor, wfd_barrier, wfd_get_cprj, wfd_barrier, wfd_change_ngfft, wfd_print,&
&                             wfd_sym_ur, wfd_init, wfd_push_ug, wfd_reset_ur_cprj, fft_ur, wfd_ug2cprj, wfd_ihave_ur,&
&                             wfd_test_ortho, kdata_t, kdata_init, kdata_free

 implicit none

 private 
!!***

! Flags used to define the content of ovlp%mat.
 integer,private,parameter :: TYPE_NONE    = 0
 integer,private,parameter :: TYPE_OVERLAP = 1
 integer,private,parameter :: TYPE_EIGVEC  = 2

!----------------------------------------------------------------------

!!****t* m_shirley/ovlp_t
!! NAME
!!  ovlp_t
!! 
!! FUNCTION
!!  Structure used to store the overlap matrix.
!! 
!! SOURCE

 type,public :: ovlp_t
                                                                                  
  integer :: size=0
  ! The size of the overlap matrix.
                                                                                  
  integer :: mband=0
  ! Maximum number of bands.
                                                                                  
  integer :: nkpt=0
  ! Number of k-points.

  integer :: mat_type = TYPE_NONE

  !%real(dp),pointer :: klist(:,:)  SET2NULL
  ! klist(3,nkpt)
  ! The list of k-points used to construct the overlap

  real(dp) :: min_ene = -HUGE(one)
  ! Min energy included.
                                                                                  
  real(dp) :: max_ene = +HUGE(one)
  ! Max energy included 

  integer,pointer :: bk2idx(:,:)  SET2NULL
  ! bk2idx(mband,nkpt) 
  ! Mapping (b,k) --> i 
  ! 0 if the state has been excluded with the energy window.
                                                                                  
  integer,pointer :: idx2bk(:,:)  SET2NULL
  ! idx2bk(2,size))
  ! Mapping i --> (b,k) for i=1,size. k is always the in the BZ.
                                                                                  
  complex(dpc),pointer :: mat(:,:)  SET2NULL     
  ! mat(size,size) 
  ! Stores O_{ij} = <u_i|u_j> 
  ! The matrix is Hermitian hence one could use a packed matrix to save memory, 
  ! but the Lapack routines for the diagonalization are usually slower.

  real(dp),pointer :: eigene(:)  SET2NULL
  ! eigene(size)
  ! The eigenvalues of the overlap operator.

 end type ovlp_t

 public :: ovlp_nullify
 public :: ovlp_free
 public :: ovlp_init
 public :: ovlp_diago_and_prune
 public :: ovlp_diff
 public :: wfd_shirley_to_eh
!!***

!----------------------------------------------------------------------

!!****t* m_shirley/ksintp_t
!! NAME
!!  ksintp_t
!! 
!! FUNCTION
!! Structure used to store the results of the Shirley interpolation for a given k-point and spin.
!! TODO fix issue with the indexing as we might want to read a subset of bands 
!! for ourinterpolation moreover one might used vl and vu in the direct diagonalization
!! to select the energy window we are interested in.
!! 
!! SOURCE

 type,public :: ksintp_t

   integer :: sh_size=0
   ! Numer of Shirley optimal basis set elements.
                                                 
   integer :: nband_k=0
   ! Number of interpolated Kohn-Sham bands.

   !real(dp) :: kpt(3)=HUGE(zero)
                                                 
   real(dp),pointer :: ene(:)  SET2NULL
   ! ene(nband_k)
   ! Interpolated energies
                                                 
   complex(dpc),pointer :: obloch(:,:)  SET2NULL
   ! obloch(sh_size,nband_k)
   ! Matrix storing <O_i|Bloch_b>
                                                 
 end type ksintp_t

 public :: ks_intp_nullify
 public :: ks_intp_free
 public :: ks_intp_init
 !public :: ksintp_read
 !public :: ksintp_write

 interface ks_intp_nullify
   module procedure ks_intp_nullify_0D
   module procedure ks_intp_nullify_2D
 end interface ks_intp_nullify

 interface ks_intp_free
   module procedure ks_intp_free_0D
   module procedure ks_intp_free_2D
 end interface ks_intp_free
!!***

 public :: wfd_bloch_to_shirley
 public :: shirley_hks
 public :: shirley_interp

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_nullify
!! NAME
!!  ovlp_nullify
!!
!! FUNCTION
!!  Nullify the pointers in the data type.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_nullify(Ovlp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ovlp_t),intent(inout) :: Ovlp

! *************************************************************************

 !@ovlp_t
 Ovlp%mat_type = TYPE_NONE

 ! integer
 nullify(Ovlp%bk2idx)
 nullify(Ovlp%idx2bk)

 ! real
 nullify(Ovlp%eigene)
 !%nullify(Ovlp%klist)
 
 ! complex
 nullify(Ovlp%mat)

end subroutine ovlp_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_free
!! NAME
!!  ovlp_free
!!
!! FUNCTION
!!  Free the memory.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_free(Ovlp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ovlp_t),intent(inout) :: Ovlp

! *************************************************************************

 !@ovlp_t
 ! integer
 if (associated(Ovlp%bk2idx))  then
   ABI_DEALLOCATE(Ovlp%bk2idx)
 end if
 if (associated(Ovlp%idx2bk))  then
   ABI_DEALLOCATE(Ovlp%idx2bk)
 end if
 !
 ! real
 if (associated(Ovlp%eigene))  then
   ABI_DEALLOCATE(Ovlp%eigene)
 end if
 !if (associated(Ovlp%klist)) deallocate(Ovlp%klist)
 !
 ! complex
 if (associated(Ovlp%mat))   then
   ABI_DEALLOCATE(Ovlp%mat)
 end if

end subroutine ovlp_free
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_init
!! NAME
!! ovlp_init
!!
!! FUNCTION
!!  Calculates the upper triangle of the overlap matrix <u_i|u_j> for a given spin and for 
!!  all the possible combinations (k,b|k',b') with k and k' in the full Brillouin zone.
!!  The u functions are the periodic part of the Bloch wavefunctions hence there is 
!!  no selection rule in k-space for the matrix elements. The diagonal matrix elements
!!  equals one provided that the input wavefunctions are correctly normalized.
!!
!! INPUTS
!! ene_win(2,Wfd%nsppol)=Window energies for the two spins.
!! Wfd<wfs_descriptor>=Datatype gathering info on the wavefunctions.   
!! use_sym=logical flag defining whether symmetries have to be used for reconstructing the overlap matrix.
!! Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!! Kmesh<BZ_mesh_type>=Datatype describing the k-point sampling used for the wavefunctions.
!! Bst<bandstructure_type>=Band structure energies.
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang<pawang_type> angular mesh discretization and related data.
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_init(O,ene_win,use_sym,ov_ngfft,Wfd,Cryst,Kmesh,Bst,Psps,Pawtab,Pawang,Pawrad)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_init'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: use_sym
 type(Crystal_structure),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(bandstructure_type),intent(in) :: Bst
 type(wfs_descriptor),intent(inout) :: Wfd
 type(ovlp_t),intent(out) :: O(Wfd%nsppol)
!arrays
 integer,intent(in) :: ov_ngfft(18)
 real(dp),intent(in) :: ene_win(2,Wfd%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: dim_rtwg1=1,option_lob0=0
 integer :: ik1_bz,ik2_bz,ik1_ibz,ik2_ibz,istat,ierr,rhoxsp_method,nqpt,iqpt
 integer :: band1,band2,nband_k,nband_k2,nband_k1
 integer :: npw_k1,npw_k2,spin !,ii
 integer :: with_sym,without_sym,row,col,ovlp_size
 integer :: k1_sym,k2_sym,k1_tim,k2_tim,inv_k2_sym,inv_k1_sym
 integer :: band1_stop,nspinor,lk,rk,iq_found,spin_comm
 integer :: itypat,klmn,klmn_size,nsppol,mband,nkpt !,npw_k,istwf_k
 real(dp),parameter :: tnons_tol=tol8,tol_kdiff=tol6
 real(dp) :: fft_fact,ksq,qsq_max,ene_bks
 complex(dpc) :: blk_ovlp,covlp
 complex(dpc) :: k1_eimkt,k2_eimkt,tnons_fact
 complex(gwpc) :: paw_ovlp(1)
 logical :: k1_isirred,k2_isirred,can_use_sym,take_cnj,fft_isok
 character(len=500) :: msg
!arrays
 !integer :: got(Wfd%nproc)
 integer,allocatable :: tmp_idx2bk(:,:)
 integer :: k1_umklp(3),k2_umklp(3)
 integer :: g_gamma(3)=(/0,0,0/)
 !integer,allocatable :: bbp_ks_distrb(:,:,:,:)
 integer,allocatable :: toinv(:,:),multable(:,:,:)
 integer,allocatable :: nq_spl(:),kmk2q(:,:) 
 real(dp) :: kpoint(3),fft_err(3,Cryst%nsym) !,onsite(2)
 real(dp) :: kk1(3),kk2(3),k2mk1(3) !,ovlp_paw(2) !,r1_tau3(3)
 real(dp),allocatable :: qpts(:,:),qmax(:)
 complex(gwpc),allocatable :: ur1(:),ur2(:)
 complex(gwpc),pointer :: ug2(:) !ug1(:)
 complex(dpc),allocatable :: ovlp_ikfk(:,:,:,:)
 !logical,allocatable :: bbp_mask(:,:)
 logical :: k_needs_tr(2)
 type(cprj_type),pointer :: Cp_k1(:,:),Cp_k2(:,:)
 type(Paw_pwij_type),allocatable :: Pwij(:,:)
 type(Paw_pwff_type),allocatable :: Paw_pwff(:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")

 ABI_UNUSED(Pawrad(1)%mesh_size)

 call wfd_change_ngfft(Wfd,Cryst,Psps,ov_ngfft) 

 call fft_check_rotrans(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,fft_err,fft_isok)
                                                                                                                          
 if (.not.fft_isok) then
   write(msg,'(a,3(i0,1x),a)')" Real space FFT mesh ",Wfd%ngfft(1:3)," is not symmetric. Cannot symmetrize in real space"
   MSG_WARNING(msg)
   !%MSG_ERROR(msg)
 end if

 nspinor = Wfd%nspinor
 nsppol  = Wfd%nsppol
 mband   = Wfd%mband

 nkpt    = Kmesh%nbz

 do spin=1,nsppol
  call ovlp_nullify(O(spin))
 end do

 do spin=1,nsppol
   spin_comm = Wfd%bks_comm(0,0,spin)       
   ABI_CHECK(xcomm_size(spin_comm)==1,"ovlp_init is not parallelized")
   !
   O(spin)%mband   = mband
   O(spin)%nkpt    = nkpt
   O(spin)%min_ene = ene_win(1,spin)
   O(spin)%max_ene = ene_win(2,spin)
   !
   ! The size of the overlap matrix and useful tables.
   ABI_ALLOCATE(O(spin)%bk2idx,(mband,nkpt))
   O(spin)%bk2idx=0
   ABI_ALLOCATE(tmp_idx2bk,(2,mband*nkpt))
   tmp_idx2bk=0

   ovlp_size=0
   do ik2_bz=1,nkpt
     ik2_ibz  = Kmesh%tab(ik2_bz)
     nband_k2 = Wfd%nband(ik2_ibz,spin)
     !
     do band2=1,nband_k2
       ene_bks = Bst%eig(band2,ik2_ibz,spin)
       if ( ene_bks >= O(spin)%min_ene .and. ene_bks <= O(spin)%max_ene ) then
         ovlp_size = ovlp_size + 1
         O(spin)%bk2idx(band2,ik2_bz) = ovlp_size
         tmp_idx2bk(1,ovlp_size) = band2
         tmp_idx2bk(2,ovlp_size) = ik2_bz
       end if
     end do
   end do

   O(spin)%size = ovlp_size
   ABI_ALLOCATE(O(spin)%idx2bk,(2,ovlp_size))
   O(spin)%idx2bk = tmp_idx2bk(:,1:ovlp_size)
   ABI_DEALLOCATE(tmp_idx2bk)
   !
   ! Allocate overlap matrix. Could use packed matrix to save memory, but Lapack call is slower.
   ABI_ALLOCATE(O(spin)%mat,(ovlp_size,ovlp_size))
   istat = ABI_ALLOC_STAT
   if (istat/=0) then
     write(msg,'(a,f8.2,a)')" out of memory in O%mat, requiring :",two*dpc*ovlp_size**2*b2Gb," Gb"
     MSG_ERROR(msg)
   end if
   !O(spin)%mat = -HUGE(one)
   write(msg,'(a,f12.1,a,i0)')" Memory required for the overlap matrix: ",two*dpc*ovlp_size**2*b2Mb," Mb; Matrix size= ",ovlp_size
   call wrtout(std_out,msg,"COLL")
 end do
 !
 ! Calculate the overlap matrix  --------------------------------------------------------------------------

 ABI_ALLOCATE(multable,(4,Cryst%nsym,Cryst%nsym))
 ABI_ALLOCATE(toinv,(4,Cryst%nsym))
 call sg_multable(Cryst%nsym,Cryst%symafm,Cryst%symrel,Cryst%tnons,tnons_tol,ierr,multable=multable,toinv=toinv)
 ABI_CHECK(ierr==0,"Group error, cannot continue")
 !
 ! ======================================================
 ! ====  <u_{ik b1}| u_{fk b2}>, ik in IBZ, fk in BZ ====
 ! ======================================================
 !
 ! For PAW we have to recalculate set projections in the IBZ setting k=0 in the exponential.
 ! TODO: For the time being, these terms are calculated in the BZ, symmetrization will be added afterwards.
 if (Wfd%usepaw==1) then
   !
   ! Find the list of q-points: qq=kk1-kk2 and create table (kk1, kk2) --> qq.
   ! <u_k2|u_k1> = <phi_k2|e^{-i(k2-k1).r}|psi_k1>
   ABI_ALLOCATE(qpts,(3,2*nkpt))
   ABI_ALLOCATE(kmk2q,(nkpt,nkpt))

   nqpt=0; qsq_max=zero
   do ik1_bz=1,nkpt 
     kk1 = Kmesh%bz(:,ik1_bz)
     do ik2_bz=1,nkpt
       kk2 = Kmesh%bz(:,ik2_bz)
       k2mk1 = kk2 - kk1
       ksq = normv(k2mk1,Cryst%gmet,"G") 
       qsq_max = MAX(ksq,qsq_max)

       iq_found=0
       do iqpt=1,nqpt
         if (ALL( ABS(k2mk1 - qpts(:,iqpt)) < tol_kdiff) ) then
           iq_found=iqpt
           EXIT
         end if
       end do

       if (iq_found==0) then 
         nqpt = nqpt+1
         qpts(:,nqpt) = k2mk1
         iq_found = nqpt
       end if
       kmk2q(ik1_bz,ik2_bz) = iq_found
     end do
   end do

   qsq_max = qsq_max/(two_pi**2)
   !
   !  * Set up q-grid for form factors, make qmax 20% larger than the largest expected.
   ABI_ALLOCATE(nq_spl,(Cryst%ntypat))
   ABI_ALLOCATE(qmax,(Cryst%ntypat))
   nq_spl = 3001 !Psps%mqgrid_ff
   qmax = SQRT(qsq_max)*1.2d0
   !write(std_out,*)"using nqpt:",nqpt," nq_spl ",nq_spl," qmax=",qmax

   rhoxsp_method=1 ! Arnaud-Alouani
   !rhoxsp_method=2 ! Shiskin-Kresse

   ABI_ALLOCATE(Paw_pwff,(Psps%ntypat))
   call init_paw_pwff(Paw_pwff,rhoxsp_method,nq_spl,qmax,Cryst%gmet,Pawrad,Pawtab,Psps)
   ABI_DEALLOCATE(nq_spl)
   ABI_DEALLOCATE(qmax)

   ABI_ALLOCATE(Pwij,(Psps%ntypat,nqpt))
   do iqpt=1,nqpt
     call init_paw_pwij(Pwij(:,iqpt),1,qpts(:,iqpt),g_gamma,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
     if (iqpt==1) then
       write(std_out,*)qpts(:,iqpt)
       do itypat=1,Cryst%ntypat
         klmn_size = Pwij(itypat,iqpt)%lmn_size*(Pwij(itypat,iqpt)%lmn_size+1)/2
         do klmn=1,klmn_size
           !write(std_out,*)" mqpgij: ",Pwij(itypat,iqpt)%mqpgij(:,1,klmn),pawtab(itypat)%sij(klmn)
           write(std_out,*)" diff mqpg0ij: ",Pwij(itypat,iqpt)%mqpgij(1,1,klmn)-pawtab(itypat)%sij(klmn)
           Pwij(itypat,iqpt)%mqpgij(1,1,klmn) = pawtab(itypat)%sij(klmn)
           Pwij(itypat,iqpt)%mqpgij(2,1,klmn) = zero
         end do
       end do
     end if
   end do
   call destroy_paw_pwff(Paw_pwff)
   ABI_DEALLOCATE(Paw_pwff)
   ABI_DEALLOCATE(qpts)

   ABI_ALLOCATE(Cp_k1 ,(Cryst%natom,nspinor))
   call cprj_alloc(Cp_k1,0,Wfd%nlmn_atm)
                                                                                                              
   ABI_ALLOCATE(Cp_k2 ,(Cryst%natom,nspinor))
   call cprj_alloc(Cp_k2,0,Wfd%nlmn_atm)
 end if ! PAW
 !
 ! ==============================================================
 ! ==== Reconstruct full <u_{kb}| u_{k'b'}> matrix in the BZ ====
 ! ==============================================================
 ! 1) Symmetrization is done in real space. Easier, especially when k-centered G-sphere are used.
 !     u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz) 
 !               =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
 ! 2) Matrix is Hermitian.
 ! 3) <u_{Sk b}| u_{Sk b'}> are obtained from the previously calculated <u_{kb}| u_{kb'}> table.
 !
 ! {A,a} {B,b} = {AB, a+Ab}
 !
 fft_fact = one/Wfd%nfft
 ABI_ALLOCATE(ur1,(Wfd%nfft*nspinor))
 ABI_ALLOCATE(ur2,(Wfd%nfft*nspinor))

 do spin=1,nsppol
   O%mat_type = TYPE_OVERLAP

   SELECT CASE (use_sym)
   
   CASE (.FALSE.)
     call wrtout(std_out," Using version without symmetries","COLL")
     
     do ik2_bz=1,nkpt
       call get_BZ_item(Kmesh,ik2_bz,kk2,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp,k2_isirred)
       nband_k2 = Wfd%nband(ik2_ibz,spin)
       npw_k2   = Wfd%npwarr(ik2_ibz)
       !
       do band2=1,nband_k2 
         !
         col = O(spin)%bk2idx(band2,ik2_bz)
         ug2 => Wfd%Wave(band2,ik2_ibz,spin)%ug
         call wfd_sym_ur(Wfd,Cryst,Kmesh,band2,ik2_bz,spin,ur2)

         if (Wfd%usepaw==1) then 
           call wfd_get_cprj(Wfd,band2,ik2_ibz,spin,Cryst,Cp_k2,sorted=.FALSE.)
           if (.not.k2_isirred) call paw_symcprj(ik2_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cp_k2) 
         end if 
         !
         do ik1_bz=1,ik2_bz 
           call get_BZ_item(Kmesh,ik1_bz,kk1,ik1_ibz,k1_sym,k1_tim,k1_eimkt,k1_umklp,k1_isirred)
           npw_k1     = Wfd%npwarr(ik1_ibz)
           nband_k1   = Wfd%nband(ik1_ibz,spin)
           band1_stop = nband_k1; if (ik1_bz==ik2_bz) band1_stop = band2 

           if (Wfd%usepaw==1) iqpt = kmk2q(ik1_bz,ik2_bz) 
                                                                                                                  
           do band1=1,band1_stop
             row= O(spin)%bk2idx(band1,ik1_bz)

             !$if (ik1_bz==ik2_bz) then
             !!  blk_ovlp = czero
             !!  if (band1==band2) blk_ovlp = cone
             !$else 

               !if (ik2_bz==ik1_bz .and. k1_isirred) then
               !  ug1 => Wfd%Wave(band1,ik1_ibz,spin)%ug
               !  blk_ovlp = xdotc(npw_k1*nspinor,ug1,1,ug2,1)
               !else 
                 call wfd_sym_ur(Wfd,Cryst,Kmesh,band1,ik1_bz,spin,ur1)
                 blk_ovlp = xdotc(Wfd%nfft*nspinor,ur1,1,ur2,1)/Wfd%nfft
               !end if

               if (Wfd%usepaw==1) then
                 call wfd_get_cprj(Wfd,band1,ik1_ibz,spin,Cryst,Cp_k1,sorted=.FALSE.)
                 if (.not.k1_isirred) call paw_symcprj(ik1_bz,nspinor,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cp_k1) 

                 !if (iqpt==1) then
                 !  onsite = paw_overlap(Cp_k1,Cp_k2,Cryst%typat,Pawtab)
                 !  blk_ovlp = blk_ovlp + DCMPLX(onsite(1),onsite(2))
                 !else 
                 paw_ovlp = czero
                 call paw_rho_tw_g(1,dim_rtwg1,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,g_gamma,&
&                  Cp_k1,Cp_k2,Pwij(:,iqpt),paw_ovlp)
!  FIXME not addinf the onsite term works just fine!
                 !blk_ovlp = blk_ovlp + paw_ovlp(1) 
                 !end if
               end if
             !end if
             !
             ! Save overlap.
             O(spin)%mat(row,col) = blk_ovlp
           end do
         end do
       end do
     end do

   CASE (.TRUE.)  ! FIXME does not work yet. Problems with umklapp symmorphic operations and symmetry tables somewhere.

     call wrtout(std_out," Using version with symmetries","COLL")
     MSG_WARNING(" Using version with symmetries (still under development")

     ABI_CHECK(Wfd%usepaw==0,"PAW not coded yet")
     ABI_ALLOCATE(ovlp_ikfk,(mband,Wfd%nkibz,mband,nkpt))
     ovlp_ikfk=+HUGE(one)
     
     do ik2_bz=1,nkpt
       call get_BZ_item(Kmesh,ik2_bz,kk2,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp)
       nband_k2 = Wfd%nband(ik2_ibz,spin)

       do band2=1,nband_k2 
         call wfd_sym_ur(Wfd,Cryst,Kmesh,band2,ik2_bz,spin,ur2)

         do ik1_ibz=1,Wfd%nkibz
           nband_k = Wfd%nband(ik1_ibz,spin)
           do band1=1,nband_k

             call wfd_get_ur(Wfd,band1,ik1_ibz,spin,ur1)
             covlp = xdotc(Wfd%nfft*nspinor,ur1,1,ur2,1) * fft_fact

             if (Wfd%usepaw==1) then ! TODO Add onsite term after IBZ-->BZ symmetrization
             end if

             ovlp_ikfk(band1,ik1_ibz,band2,ik2_bz) = covlp 
           end do
         end do
         !
       end do
     end do
     !
     ! FIXME: Equations are not completed. non-symmorphic phases are missing!
     ! Let <r|k> indicate the periodic part of the Bloch wavefunction that transforms according to:
     !   u_{Sk} = e^{-iSk.t} u_k(S^{-1}(r-t))
     !   S3 = S1^{-1} S2
     !
     ! 1) <S1 k1 | S2 k2> = <k1| S1^{-1}   S2 k2>  e^{i (S1 k1 - S2 k2).t1}  
     !
     ! 2) <T S1 k1 | T S2 k2> = <S1 k1| S2 k2>*     

     ! <S1 k1   | T S2 k2> = <k1| S1^{-1} T S2 k2>
     ! <T S1 k1 |   S2 k2> = <k2| S2^{-1} T S1 k1>^*
     ! <T S1 k1 | T S2 k2> = <k2| S2^{-1}   S1 k1>     ! Problematic if the BZ mesh is not invariant under inversion.
     !                                                   e.g. randomly shifted k-meshes. In this case one should use kptopt=3

     call wrtout(std_out,"Using version with symmetries","COLL")
     with_sym=0; without_sym=0
     do ik2_bz=1,nkpt

       call get_BZ_item(Kmesh,ik2_bz,kk2,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp)

       !$ik2_ibz = Kmesh%tab(ik2_bz)
       nband_k2 = Wfd%nband(ik2_ibz,spin)

       inv_k2_sym = toinv(1,k2_sym)
       k_needs_tr(2) = (k2_tim==2)

       do ik1_bz=1,ik2_bz !nkpt
         call get_BZ_item(Kmesh,ik1_bz,kk1,ik1_ibz,k1_sym,k1_tim,k1_eimkt,k1_umklp)

         nband_k1 = Wfd%nband(ik1_ibz,spin)

         inv_k1_sym = toinv(1,k1_sym)

         k_needs_tr(1) = (k1_tim==2)

         !sym3 = multable(1,inv_k1_sym,k2_sym)

         rk= Kmesh%rottb(ik2_bz,1,inv_k1_sym)

         can_use_sym = .TRUE.
         !can_use_sym = ( ALL(ABS(Cryst%tnons(:,k2_sym))<tol6) .and. ALL(ABS(Cryst%tnons(:,k1_sym))<tol6) )
         !can_use_sym = can_use_sym .and. ( ALL(k2_umklp == g_gamma) .and. ALL(k1_umklp == g_gamma))

         can_use_sym = can_use_sym .and. ( ALL(k2_umklp == g_gamma) .and. ALL(k1_umklp == g_gamma) &
         & .and. ALL(multable(2:4,inv_k1_sym,k2_sym)==g_gamma) )

         !can_use_sym = ( can_use_sym .and. &
         !& ALL( ABS ( -Kmesh%bz(:,rk) + MATMUL(Cryst%symrec(:,:,inv_k1_sym),Kmesh%bz(:,ik2_bz)) ) < tol6 ) )

         !can_use_sym = ( can_use_sym .and. ALL(ABS(Cryst%tnons(:,k1_sym)) <tol6)  .and. ALL(ABS(Cryst%tnons(:,k2_sym)) <tol6) ) 

         kpoint = kk1 - kk2
         !do ii=1,3 ! Wrap in the first BZ thus enforcing traslational invariance.
         !  call wrap2_pmhalf(kk1(ii)-kk2(ii),kpoint(ii),shift)  ! TODO overloaded interface.
         !end do

         if (ANY (ABS(Cryst%tnons(:,k1_sym)) > tol6) ) then
           !tnons_fact = EXP(j_dpc*two_pi*DOT_PRODUCT(kk1-kk2,Cryst%tnons(:,k1_sym)))
           tnons_fact = EXP(j_dpc*two_pi*DOT_PRODUCT(kpoint,Cryst%tnons(:,k1_sym)))
         else 
           tnons_fact = cone
         end if

         take_cnj=ALL(k_needs_tr) 

         if (ALL(k_needs_tr) .or. ALL(.not.k_needs_tr) ) then
           lk = ik1_ibz
           rk= Kmesh%rottb(ik2_bz,1,inv_k1_sym)
           !rk= Kmesh%rottbm1(ik2_bz,1,k1_sym)
         else 
           MSG_ERROR("Need TR")
           take_cnj=.FALSE.
           !if (k_needs_tr(2)) then
           !lk = ik1_ibz
           !rk=Kmesh%rottb(ik2_bz,1,inv_k1_sym)
           !tnons_fact = CONJG(k1_eimkt) * CONJG(k2_eimkt) * EXP(-j_dpc*two_pi*DOT_PRODUCT(kk2,r1_tau3))
           !end if
         end if

         if (can_use_sym) then

           with_sym=with_sym+1
           do band2=1,nband_k2 
             col=O(spin)%bk2idx(band2,ik2_bz)

             band1_stop = nband_k1; if (ik1_bz==ik2_bz) band1_stop = band2 
             do band1=1,band1_stop
               blk_ovlp = tnons_fact * ovlp_ikfk(band1,lk,band2,rk)
               if (take_cnj) blk_ovlp = DCONJG(blk_ovlp)
               row=O(spin)%bk2idx(band1,ik1_bz)
               if (col>=row) O(spin)%mat(row,col) = blk_ovlp
             end do
           end do

         else 
           without_sym=without_sym+1
           do band2=1,nband_k2 
             col=O(spin)%bk2idx(band2,ik2_bz)

             call wfd_sym_ur(Wfd,Cryst,Kmesh,band2,ik2_bz,spin,ur2)

             band1_stop = nband_k1; if (ik1_bz==ik2_bz) band1_stop = band2 

             do band1=1,band1_stop

               call wfd_sym_ur(Wfd,Cryst,Kmesh,band1,ik1_bz,spin,ur1)

               blk_ovlp = xdotc(Wfd%nfft,ur1,1,ur2,1) * fft_fact
               row=O(spin)%bk2idx(band1,ik1_bz)
               if (col>=row) O(spin)%mat(row,col) = blk_ovlp

             end do
           end do
         end if

       end do
     end do

     write(std_out,*)"with_sym",with_sym," without_sym",without_sym
     ABI_DEALLOCATE(ovlp_ikfk)

   END SELECT
   !
   ! Collect ovlp_mat on each node inside spin_comm.
   call xsum_mpi(O(spin)%mat,spin_comm,ierr)
   !call print_arr(O(spin)%mat,max_r=10,max_c=12,unit=std_out,mode_paral="COLL")
 end do ! spin

 ABI_DEALLOCATE(multable)
 ABI_DEALLOCATE(toinv)
 ABI_DEALLOCATE(ur1)
 ABI_DEALLOCATE(ur2)
 !
 if (Wfd%usepaw==1) then
   ABI_DEALLOCATE(kmk2q)
   call destroy_paw_pwij(Pwij)
   ABI_DEALLOCATE(Pwij)
   call cprj_free(Cp_k1)
   ABI_DEALLOCATE(Cp_k1)
   call cprj_free(Cp_k2)
   ABI_DEALLOCATE(Cp_k2)
 end if

 !call wfd_barrier(Wfd)

 DBG_EXIT("COLL")

end subroutine ovlp_init
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_diago_and_prune
!! NAME
!! ovlp_diago_and_prune
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_diago_and_prune(O,sh_coverage,sh_size,base)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_diago_and_prune'
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: sh_size,base
 real(dp),intent(in) :: sh_coverage
 type(ovlp_t),intent(inout) :: O

!Local variables ------------------------------
!scalars
 integer :: ii
 real(dp) :: sum_eige,trace
 real(dp) :: cpu_time,wall_time,cpu0,wall0
 character(len=500) :: msg

!************************************************************************

 ! @ovlp_t
 call timein(cpu0,wall0)
 !
 ! Diagonalize the overlap matrix.
 ABI_ALLOCATE(O%eigene,(O%size))

 call xheev("Vectors","Upper",O%size,O%mat,O%eigene)
 O%mat_type = TYPE_EIGVEC

 trace = SUM(O%eigene)
 write(msg,'(3(a,f8.2))')" Trace: ",trace," Min eig: ",MINVAL(O%eigene)," Max eig: ",MAXVAL(O%eigene)
 call wrtout(std_out,msg,"COLL")

 !do ii=1,O%size
 !  write(345,*)ii, O%eigene(ii)
 !end do
 !
 ! Select the optimal subspace.
 sum_eige=zero; ii=O%size
 do while ( sum_eige < trace * sh_coverage .and. ii/=1)
   sum_eige = sum_eige + O%eigene(ii)
   ii = ii-1
 end do
 !sh_size=O%size; base=0
 sh_size=O%size-ii; base=ii 

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0
 write(std_out,*)" Ovlp_diago cpu_time, wall_time",cpu_time, wall_time

 write(msg,'(2(a,i0),a)')" The optimal basis set contains (",sh_size,"/",O%size,") elements."
 call wrtout(std_out,msg,"COLL")

end subroutine ovlp_diago_and_prune
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ovlp_diff
!! NAME
!! ovlp_diff
!!
!! FUNCTION
!!  Debugging tool used to compare the overlap matrices calculated  with and without symmetried.
!!
!! INPUTS
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ovlp_diff(O1,O2,Kmesh,tol,header,unit,mode_paral,prtvol) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_diff'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 real(dp),intent(in) :: tol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(ovlp_t),intent(in) :: O1,O2
 type(BZ_mesh_type),intent(in) :: Kmesh

!Local variables ------------------------------
!scalars
 integer :: row,col,my_unt,my_prtvol
 integer :: ik1_bz,ik2_bz,ik1_ibz,ik2_ibz,band1
 integer :: k1_sym,k2_sym,k1_tim,k2_tim
 complex(dpc) :: k1_eimkt,k2_eimkt
 character(len=4) :: my_mode
 character(len=500) :: msg      
!arrays
 integer :: k1_umklp(3),k2_umklp(3)
 real(dp) :: kk1(3),kk2(3)

!************************************************************************

 MSG_ERROR("Not tested")

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
                                                                    
 msg=' ==== Diff O1-O2  ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 ABI_CHECK(O1%size==O1%size,"Different sizes")
 !
 ! Compare the overlap matrices 
 do col=1,O2%size
   do row=1,O1%size
     !
     if (ABS(O1%mat(row,col)-O2%mat(row,col))>tol) then
       ik1_bz = O1%idx2bk(1,row)
       band1  = O1%idx2bk(2,row)
       call get_BZ_item(Kmesh,ik1_bz,kk1,ik1_ibz,k1_sym,k1_tim,k1_eimkt,k1_umklp)
       !
       ik2_bz = O1%idx2bk(1,col)
       band1  = O1%idx2bk(2,col)
       call get_BZ_item(Kmesh,ik2_bz,kk2,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp)

       write(my_unt,'(2i3,2(2x,a,3f7.3),4f8.4)')&
&        row,col," bz1 ",Kmesh%bz(:,ik1_bz)," bz2 ",Kmesh%bz(:,ik2_bz),O1%mat(row,col),O2%mat(row,col)
       write(my_unt,'(a,i3,3i3,2f4.1,3i3)')"k2 ",ik2_bz,ik2_ibz,k2_sym,k2_tim,k2_eimkt,k2_umklp
       write(my_unt,'(a,i3,3i3,2f4.1,3i3)')"k1 ",ik1_bz,ik1_ibz,k1_sym,k1_tim,k1_eimkt,k1_umklp

     end if
   end do
 end do

end subroutine ovlp_diff
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ks_intp_nullify_0D
!! NAME
!!  ks_intp_nullify_0D
!!
!! FUNCTION
!!  Nullify the pointers in the data type.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ks_intp_nullify_0D(Ksintp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ks_intp_nullify_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ksintp_t),intent(inout) :: Ksintp

! *************************************************************************

 !@ksintp_t

 ! real
 nullify(Ksintp%ene)
 
 ! complex
 nullify(Ksintp%obloch)

end subroutine ks_intp_nullify_0D
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ks_intp_nullify_2D
!! NAME
!!  ks_intp_nullify_2D
!!
!! FUNCTION
!!  Nullify the pointers in the data type.
!!
!! PARENTS
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ks_intp_nullify_2D(Ksintp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ks_intp_nullify_2D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ksintp_t),intent(inout) :: Ksintp(:,:)

!Local variables
 integer :: ii,jj

! *************************************************************************

 do jj=LBOUND(Ksintp,DIM=2),UBOUND(Ksintp,DIM=2)
   do ii=LBOUND(Ksintp,DIM=1),UBOUND(Ksintp,DIM=1)
     call ks_intp_nullify_0d(Ksintp(ii,jj))
   end do
 end do

end subroutine ks_intp_nullify_2D
!!***

!----------------------------------------------------------------------


!!****f* m_shirley/ksint_free_0D
!! NAME
!!  ks_intp_free_0D
!!
!! FUNCTION
!!  Free the memory.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ks_intp_free_0D(Ksintp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ks_intp_free_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ksintp_t),intent(inout) :: Ksintp

! *************************************************************************

 !@ksintp_t
 !
 ! real
 if (associated(Ksintp%ene))  then
   ABI_DEALLOCATE(Ksintp%ene)
 end if
 !
 ! complex
 if (associated(Ksintp%obloch))   then
   ABI_DEALLOCATE(Ksintp%obloch)
 end if

end subroutine ks_intp_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ks_intp_free_2D
!! NAME
!!  ks_intp_free_2D
!!
!! FUNCTION
!!  Free the pointers in the data type.
!!
!! PARENTS
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ks_intp_free_2D(Ksintp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ks_intp_free_2D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ksintp_t),intent(inout) :: Ksintp(:,:)

!Local variables
 integer :: ii,jj

! *************************************************************************

 do jj=LBOUND(Ksintp,DIM=2),UBOUND(Ksintp,DIM=2)
   do ii=LBOUND(Ksintp,DIM=1),UBOUND(Ksintp,DIM=1)
     call ks_intp_free_0d(Ksintp(ii,jj))
   end do
 end do

end subroutine ks_intp_free_2D
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/ks_intp_init
!! NAME
!!  ks_intp_init
!!
!! FUNCTION
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine ks_intp_init(Ksintp,sh_size,nband_k,intp_ene,obloch,istat)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ks_intp_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sh_size,nband_k
 integer,intent(out) :: istat
 type(ksintp_t),intent(inout) :: Ksintp
!arrays
 real(dp),intent(in) :: intp_ene(nband_k)
 complex(dpc),intent(in) :: obloch(sh_size,nband_k)

!Local variables ------------------------------
!scalars

! *************************************************************************

 !@ksintp_t
 !
 Ksintp%sh_size = sh_size
 Ksintp%nband_k = nband_k

 ABI_ALLOCATE(Ksintp%ene,(nband_k))
 Ksintp%ene = intp_ene

 ABI_ALLOCATE(Ksintp%obloch,(sh_size,nband_k))
 istat = ABI_ALLOC_STAT
 !ABI_CHECK(istat==0,"out of memory %obloch")
 if (istat/=0) return
 Ksintp%obloch = obloch

end subroutine ks_intp_init
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/wfd_bloch_to_shirley
!! NAME
!! wfd_bloch_to_shirley
!!
!! FUNCTION
!!  Returns a new wavefunction descriptor containing the basis set proposed by Shirley.
!!
!! INPUT
!! Wfd<wfs_descriptor>
!! Cryst<crystal_structure>= data type gathering info on symmetries and unit cell
!! Psps<pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat)<pawtab_type>=paw tabulated starting data
!! Pawang <pawang_type>=angular mesh discretization and related data:
!! sh_coverage
!!
!! OUTPUT
!!  Wsh<wfs_descriptor>
!!
!! PARENTS
!!      exc_interp_ham
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine wfd_bloch_to_shirley(Wfd,Cryst,Kmesh,Bst,Psps,Pawtab,Pawang,Pawrad,min_bsize,sh_coverage,Wsh)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_bloch_to_shirley'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: min_bsize
 real(dp),intent(in) :: sh_coverage
 type(Crystal_structure),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wfd
 type(bandstructure_type),intent(in) :: Bst
 type(wfs_descriptor),intent(out) :: Wsh
!arrays
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: istwf1=1,k1=1
 integer :: ik1_bz,istat !,ierr
 integer :: band1,natom,sh,midx,ig
 integer :: npw_gamma,paral_kgb,usepaw,spin !,row,col,useylm_
 integer :: fft_idx,npw_k,nspinor,nsppol,sh_nkibz,sh_mband
 integer :: method,enforce_sym,ii,ik_ibz,istwf_k,tmp_nfft
 real(dp) :: fft_fact !,norm1 !sqrt_norm1,sqrt_norm2
 real(dp) :: cpu_time,wall_time,cpu0,wall0
 logical :: use_sym,fft_isok
 character(len=500) :: msg
!arrays
 integer :: ov_ngfft(18),trial_ngfft(18)
 integer :: sh_istwfk(1),sh_size(Wfd%nsppol),base(Wfd%nsppol)
 !integer :: got(Wfd%nproc)
 !integer,allocatable :: bbp_ks_distrb(:,:,:,:)
 integer,pointer :: kg_k(:,:),kg_gamma(:,:) !,gbound_k(:,:)
 integer,allocatable :: sh_nband(:,:) 
 integer,pointer :: igfft0(:)
 real(dp) :: ene_win(2,Wfd%nsppol)
 real(dp) :: sh_kibz(3,1),gamma_point(3)=(/zero,zero,zero/)
 real(dp) :: pawovlp(2),fft_err(3,Cryst%nsym)
 complex(dpc) :: cdum
 complex(dpc),allocatable :: dpc_tmp(:)
 complex(gwpc),allocatable :: ur1(:) !,ur2(:)
 complex(gwpc),pointer :: ug1(:)
 complex(gwpc),allocatable :: sh_ug(:),sh_ur(:,:),gwpc_tmp(:)
 logical,allocatable :: sh_keep_ur(:,:,:),sh_bks_mask(:,:,:),kg_mask(:)
 type(cprj_type),allocatable :: Cp_sh1(:,:) 
 type(ovlp_t) :: O(Wfd%nsppol),Osym(Wfd%nsppol)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded.")
 ABI_CHECK(wfd%paral_kgb==0,"paral_kgb/=0 not coded.")
 ABI_CHECK(Wfd%rfft_is_symok,"Real space FFT is not symmetric.")

 ABI_UNUSED(Pawrad(1)%mesh_size)

 nspinor   = Wfd%nspinor
 nsppol    = Wfd%nsppol
 paral_kgb = Wfd%paral_kgb
 usepaw    = Wfd%usepaw
 natom     = Wfd%natom

 ! One should have a single spin-undependent basis set!
 ABI_CHECK(nsppol==1,"nsppol==2 not coded")

 write(msg,'(a,f9.6)')" Transformation Bloch --> Shirley basis set with coverage : ",sh_coverage
 call wrtout(std_out,msg,"COLL")

 ! 1) Get the overlap matrix <u_i|u_j> for this spin. 
 ! TODO
 ! The Overlap can be calculated using a coarse real space FFT mesh provided that it is compatible 
 ! with the space group symmetries. We only have to make sure that ngfft encloses the k-centered G-spheres. 
 method      = 1
 enforce_sym = 1
 ov_ngfft(1:6) = 0
                                                                                                             
 do ik_ibz=1,Wfd%nkibz
   istwf_k = Wfd%istwfk(ik_ibz)
   npw_k   = Wfd%npwarr(ik_ibz)
   kg_k    => Wfd%Kdata(ik_ibz)%kg_k
   ABI_CHECK(istwf_k==1,"istwf_k/=1 not coded")
                                                                                                             
   trial_ngfft(7:) = Wfd%ngfft(7:) ! Have to preserve the previous fftalg.
   call setmesh(Cryst%gmet,kg_k,trial_ngfft,npw_k,1,npw_k,tmp_nfft,method,(/0,0,0/),Cryst,enforce_sym)
   !
   do ii=1,6
     ov_ngfft(ii) = MAX(ov_ngfft(ii),trial_ngfft(ii))
   end do
 end do
                                                                                                             
 ov_ngfft(7:) = Wfd%ngfft(7:) ! Have to preserve the previous fftalg.

 use_sym =.FALSE.
 ene_win(1,:) = smallest_real
 ene_win(2,:) = greatest_real ! (TODO Be careful here in parallel, add the possibility of using an energy window.

 call ovlp_init(O,ene_win,use_sym,ov_ngfft,Wfd,Cryst,Kmesh,Bst,Psps,Pawtab,Pawang,Pawrad)

 if (.FALSE.) then  ! TODO Work in progress
   use_sym =.TRUE.
   call ovlp_init(Osym,ene_win,use_sym,ov_ngfft,Wfd,Cryst,Kmesh,Bst,Psps,Pawtab,Pawang,Pawrad)

   do spin=1,nsppol
     call ovlp_diff(O(spin),Osym(spin),Kmesh,tol6,"Error in matrix elements",unit=std_out) 
     call ovlp_free(Osym(spin))
   end do
   MSG_ERROR("Check done")
 end if
 !
 ! 2) Diagonalize the overlap matrix selecting the optimal subspace: [ base(spin):ovlp_size ]
 do spin=1,nsppol
   !
   call ovlp_diago_and_prune(O(spin),sh_coverage,sh_size(spin),base(spin)) ! In exit O%mat stores the eigenvectors.
   !
   ! Make sure we have enough states.
   if (sh_size(spin)<min_bsize) then
     if (O(spin)%size<min_bsize) then
       write(msg,'(2(a,i0),2a)')&
&        " Overlap size is ",O(spin)%size," whereas min_bsize is ",min_bsize,ch10,&
&        " Decrease the number of bands to be interpolated or increase the number of ab-initio input states."
       MSG_ERROR(msg)
     end if
     sh_size(spin) = min_bsize
     write(msg,'(a,2i0)')" Had to enlarge Shirley subspace since input sh_size < min_bsize: ",sh_size(spin),min_bsize
     MSG_COMMENT(msg)
   end if
   !
 end do

 call timein(cpu0,wall0)
 !
 ! 3) Init a new wavefunction descriptor to store the optimal basis set.
 !    *) Wsh must be allocated here since sh_size is know only after the diagonalization of the overlap.
 !    *) Keep the optimal wavefunctions on each node (if possible) to facilitate the interpolation over the fine k-mesh. 
 !    *) Use Gamma-centered basis set to facilitate the operations in reciprocal space.
 !    *) The new basis is orthogonal, but not normalized since <U_i|U_j> = delta_ij e_i.
 !
 ! The optimal basis set is given on the gamma centered basis set with istwfk==1.
 ! FIXME temporary hacking. There is a bug somewhere in kpgsph
 call get_kg(gamma_point,istwf1,Wfd%ecut,Cryst%gmet,npw_gamma,kg_gamma) 
 !call get_kg(gamma_point,istwf1,14.0_dp,Cryst%gmet,npw_gamma,kg_gamma) 
 !
 ! * Index of the G-sphere in the FFT box.
 ABI_ALLOCATE(igfft0,(npw_gamma))
 ABI_ALLOCATE(kg_mask,(npw_gamma))
 call kgindex(igfft0,kg_gamma,kg_mask,Wfd%MPI_enreg,Wfd%ngfft,npw_gamma)
                                                                      
 ABI_CHECK(ALL(kg_mask),"FFT para not yet implemented")
 ABI_DEALLOCATE(kg_mask)

 sh_istwfk    = istwf1
 sh_nkibz     = 1 
 sh_kibz(:,1) = gamma_point

 ABI_ALLOCATE(sh_nband,(sh_nkibz,nsppol))
 do spin=1,nsppol
   sh_nband(:,spin)=sh_size(spin)
 end do
 ! TODO: BE careful in parallel when nsppol==2. I should recreate the communicators.
 
 sh_mband=MAXVAL(sh_nband)      
 ABI_ALLOCATE(sh_bks_mask,(sh_mband,sh_nkibz,nsppol))
 sh_bks_mask=.TRUE.
 ABI_ALLOCATE(sh_keep_ur ,(sh_mband,sh_nkibz,nsppol))
 sh_keep_ur =.TRUE.

 call wfd_init(Wsh,Cryst,Pawtab,Psps,sh_keep_ur,Wfd%paral_kgb,npw_gamma,sh_mband,sh_nband,sh_nkibz,nsppol,&
&  sh_bks_mask,Wfd%nspden,nspinor,Wfd%ecutsm,Wfd%dilatmx,sh_istwfk,sh_kibz,Wfd%ngfft,kg_gamma,Wfd%nloalg,&
&  Wfd%prtvol,Wfd%pawprtvol,Wfd%comm)

 !call wfd_print(Wsh,header="Shirley wavefunction descriptor")

 ABI_DEALLOCATE(sh_keep_ur)
 ABI_DEALLOCATE(sh_bks_mask)
 ABI_DEALLOCATE(sh_nband)
 !
 ! =====================================================================
 ! ==== Rotate the input wavefunctions to get the optimal basis set ====
 ! =====================================================================
 fft_fact = one/Wfd%nfft
 ABI_ALLOCATE(ur1,(Wfd%nfft*nspinor))

 call fft_check_rotrans(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,fft_err,fft_isok)

 if (.not.fft_isok) then
   write(msg,'(a,3(i0,1x),a)')" Real space FFT mesh ",Wfd%ngfft(1:3)," is not symmetric. Cannot symmetrize in real space"
   !MSG_WARNING(msg)
   MSG_ERROR(msg)
 end if

 do spin=1,nsppol
   !
   write(msg,'(a,f12.1,a)')' Memory needed for storing sh_ur= ',two*gwpc*Wfd%nfft*nspinor*sh_size(spin)*b2Mb,' [Mb]'
   call wrtout(std_out,msg,'PERS')

   ABI_ALLOCATE(sh_ur,(Wfd%nfft*nspinor,sh_size(spin)))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out-of-memory in sh_ur")
   sh_ur = czero

   do midx=1,O(spin)%size ! Loop over the Bloch set in the BZ.
     band1  = O(spin)%idx2bk(1,midx)
     ik1_bz = O(spin)%idx2bk(2,midx)

     call wfd_sym_ur(Wfd,Cryst,Kmesh,band1,ik1_bz,spin,ur1)
     !
     do sh=1,sh_size(spin) ! Construct the new optimal basis set.
       sh_ur(:,sh) = sh_ur(:,sh) + O(spin)%mat(midx,base(spin)+sh) * ur1
     end do
     !
   end do
   !
   ! NC: Normalize the basis set using the eigenvalues of the overlap matrix.
   !if (usepaw==0) then
   do sh=1,sh_size(spin)
#if 1
     sh_ur(:,sh) = sh_ur(:,sh)/SQRT(O(spin)%eigene(base(spin)+sh)) 
#else
     norm1 = xdotc(Wfd%nfft*nspinor,sh_ur(:,sh),1,sh_ur(:,sh),1) * fft_fact 
     sh_ur(:,sh) = sh_ur(:,sh)/SQRT(norm1) 
#endif
     !write(std_out,*)" sh_ur integrates to: ", xdotc(Wfd%nfft*nspinor,sh_ur(:,sh),1,sh_ur(:,sh),1) * fft_fact 
   end do
   !
   ! From the FFT mesh to the G-sphere.
   ABI_ALLOCATE(sh_ug,(npw_gamma*nspinor))
   ABI_ALLOCATE(gwpc_tmp,(Wfd%nfft*nspinor))
                                                                                            
   do sh=1,sh_size(spin) ! Construct new optimal basis set in G-space and save results in Wsh.
     gwpc_tmp = sh_ur(:,sh)

#if 1
     ABI_ALLOCATE(dpc_tmp,(Wfd%nfft*nspinor))
     dpc_tmp  = sh_ur(:,sh)
     call fourdp_c2c_ip(dpc_tmp,-1,Wfd%MPI_enreg,Wfd%nfft*nspinor,Wfd%ngfft,paral_kgb,0)
     !
     do ig=1,npw_gamma       
       fft_idx = igfft0(ig)
       if (fft_idx/=0) then ! G-G0 belong to the FFT mesh.
         if (fft_idx>Wfd%nfft .or.fft_idx<0) then
           MSG_ERROR("fft_idx bug")
         end if
         sh_ug(ig)=dpc_tmp(fft_idx) 
       else                 ! Set this component to zero.
         MSG_ERROR("fft_idx bug")
         sh_ug(ig)=czero
       end if
     end do
     ABI_DEALLOCATE(dpc_tmp)
#else
     ! FIXME does not work anymore.
     gbound_k => Wsh%Kdata(k1)%gbound
     call fft_ur(nspinor,npw_gamma,istwf1,paral_kgb,Wfd%nfft,Wfd%mgfft,Wfd%ngfft,gwpc_tmp,&
&      gbound_k,kg_gamma,sh_ug,Wfd%MPI_enreg)
#endif
     ! NC: Normalize the basis set using the eigenvalues of the overlap matrix.
     !if (usepaw==0) then
     !  sh_ug = sh_ug/SQRT(O(spin)%eigene(base(spin)+sh)) 
     !end if

     call wfd_push_ug(Wsh,sh,k1,spin,Cryst,sh_ug,update_ur=.FALSE.,update_cprj=.FALSE.)
   end do
   !
   ! PAW: Normalize the basis set using the eigenvalues of the overlap matrix.
   if (.FALSE. .and. usepaw==1) then
     call wfd_reset_ur_cprj(Wsh)
     ABI_ALLOCATE(Cp_sh1,(natom,nspinor))
     call cprj_alloc(Cp_sh1,0,Wfd%nlmn_atm)

     do sh=1,sh_size(spin)
       ug1 => Wsh%Wave(sh,k1,spin)%ug
       !cdum = xdotc(npw_gamma*nspinor,ug1,1,ug1,1) * Cryst%ucvol
       cdum = xdotc(npw_gamma*nspinor,ug1,1,ug1,1) 
       !if (istwf_k>1) then
       !  cdum=two*DBLE(cdum)
       !  if (istwf_k==2) cdum=cdum-CONJG(ug1(1))*ug1(1)
       !end if

       call wfd_get_cprj(Wsh,sh,k1,spin,Cryst,Cp_sh1,sorted=.FALSE.)

       pawovlp = paw_overlap(Cp_sh1,Cp_sh1,Cryst%typat,Pawtab)
       !pawovlp = pawovlp * Cryst%ucvol
       cdum = cdum + CMPLX(pawovlp(1),pawovlp(2))

       !norm1 = xdotc(Wfd%nfft*nspinor,sh_ur(:,sh),1,sh_ur(:,sh),1) * fft_fact * Cryst%ucvol
       !sh_ur(:,sh) = sh_ur(:,sh)/SQRT(O(spin)%eigene(base(spin)+sh)) 
       !write(std_out,*)" sh_ur integrates to: ", xdotc(Wfd%nfft*nspinor,sh_ur(:,sh),1,sh_ur(:,sh),1) * fft_fact
       !write(std_out,*)" PAW test: ",DBLE(cdum),O(spin)%eigene(base(spin)+sh)

       sh_ug = Wsh%Wave(sh,k1,spin)%ug/SQRT(DBLE(cdum))
       call wfd_push_ug(Wsh,sh,k1,spin,Cryst,sh_ug,update_ur=.FALSE.,update_cprj=.FALSE.)
     end do

     call cprj_free(Cp_sh1)
     ABI_DEALLOCATE(Cp_sh1)
   end if

   ABI_DEALLOCATE(gwpc_tmp)
   ABI_DEALLOCATE(sh_ur)
   ABI_DEALLOCATE(sh_ug)
   !
   ! Free the overlap eigenvectors.
   call ovlp_free(O(spin))
 end do ! spin
 !
 do spin=1,nsppol ! Just to be on the safe side.
   call ovlp_free(O(spin))
 end do

 ABI_DEALLOCATE(kg_gamma)
 ABI_DEALLOCATE(igfft0)
 ABI_DEALLOCATE(ur1)

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0
 write(std_out,*)" Shirley Rotation: cpu_time, wall_time",cpu_time, wall_time

 DBG_EXIT("COLL")

end subroutine wfd_bloch_to_shirley       
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/shirley_hks
!! NAME
!! shirley_hks
!!
!! FUNCTION
!!
!! INPUTS
!! Wsh: 
!! Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang<pawang_type> angular mesh discretization and related data:
!! Paw_ij(natom)<type(paw_ij_type)>=data structure containing PAW arrays given on (i,j) channels.
!! spline_opt
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shirley_hks(Wsh,kpt,spin,spline_opt,Ham_k,Cryst,Psps,Pawtab,Pawang,Paw_ij,sh_size,vloc_ij,hk_ij,sk_ij)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shirley_hks'
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_56_recipspace
 use interfaces_65_nonlocal
 use interfaces_69_wfdesc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spline_opt,spin,sh_size
 type(Crystal_structure),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wsh
 type(gs_hamiltonian_type),intent(inout) :: Ham_k
!arrays
 real(dp),intent(in) :: kpt(3)
 complex(dpc),intent(in) :: vloc_ij(sh_size,sh_size)
 complex(dpc),intent(out) ::  hk_ij(sh_size,sh_size),sk_ij(sh_size,sh_size*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: EVAL_VNLK=0,USE_BSPLINE=1
 integer,parameter :: idir0=0,ider0=0,nnlout0=0,tim_nonlop0=0,k1=1
 integer :: istat,ii,ig,ll,rr,natom,npw_k
 integer :: optder,mkmem_,nkpg,dimffnl,nspinortot
 integer :: choice,cpopt,matblk,paw_opt,signs,step
 integer :: blc1,blc2,istwf_k,useylm_ 
 integer :: iat,iatom,nnlout,nspinor
 real(dp),parameter :: lambda0=zero
 real(dp) :: shift,arg
 complex(dpc) :: kin_ij,vnl_ij
 !character(len=500) :: msg
 type(kdata_t) :: Kdata
!arrays
 integer :: nloalg(5),cp_idx(2,4),bspl_kdiv(3),bspl_kord(3)
 integer,pointer :: kg_k(:,:)
 integer,allocatable :: nlmn_sort(:)
 real(dp) :: k4intp(3),kptns_(3,1),ylmgr_dum(1),enlout(1),estep(4)
 real(dp),allocatable :: kdotg(:),half_gsq(:)
 real(dp),allocatable :: ylm_k(:,:),dum_ylm_gr_k(:,:,:) 
 real(dp),pointer :: ffnl(:,:,:,:)
 real(dp),allocatable :: kpg_k(:,:)
 real(dp),allocatable :: ph3d(:,:,:)
 real(dp),allocatable :: phkxred(:,:),vnl_psi(:,:),vectin(:,:)
 real(dp),allocatable :: s_psi(:,:)
 complex(gwpc),pointer :: ug1(:),ug2(:)
 complex(gwpc),allocatable :: cvnl_psi(:),cs_psi(:),wsg(:)
 !logical,allocatable :: bbp_mask(:,:)
 type(cprj_type),allocatable :: Cp_blc1(:,:),Cp_blc2(:,:)
 type(cprj_type),allocatable,target :: Cp_left(:,:),Cp_right(:,:)
 type(cprj_bspl_t) :: Cp_bspl

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wsh%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(Wsh%nsppol==1,"Wsh%nsppol must be 1")

 natom   = Cryst%natom
 useylm_ = Psps%useylm
 nloalg  = Wsh%nloalg
 nspinor = Wsh%nspinor
 nspinortot = nspinor

 ! Gamma-centered basis set.
 npw_k   = Wsh%Kdata(k1)%npw
 istwf_k = Wsh%istwfk(k1)
 ABI_CHECK(istwf_k==1,"istwfk/=0 not coded")

 kg_k => Wsh%Kdata(k1)%kg_k

 do ii=1,3 ! First wrap in the first BZ thus enforcing traslational invariance.
   call wrap2_pmhalf(kpt(ii),k4intp(ii),shift)  ! TODO overloaded interface.
 end do

 call finalize_hamiltonian(Ham_k,spin,npw_k,istwf_k,k4intp,Paw_ij)

 call kdata_init(Kdata,Cryst,Psps,k4intp,istwf_k,Wsh%ngfft,Wsh%MPI_enreg,kg_k=kg_k)
 !ffnl => Kdata%fnl_dir0der0

 ABI_ALLOCATE(nlmn_sort,(Cryst%natom))
 iat=0 ! nlmn dims sorted by atom type.

 if (Psps%usepaw==1) then
   nlmn_sort = Wsh%nlmn_sort
 else  ! FIXME here lmnmax == lnmax if useylm_==0
   nlmn_sort = 9
   if (Psps%useylm==1) then 
     nlmn_sort=Psps%lmnmax
   else
     MSG_ERROR("useylm==0 not coded")
   end if
   !do itypat=1,Cryst%ntypat
   !  nlmn_sort(iat+1:iat+Cryst%nattyp(itypat))=Pawtab(itypat)%lmn_size
   !  iat=iat+Cryst%nattyp(itypat)
   !end do
   !write(std_out,*)" hacking nlmn_sort",nlmn_sort
   !write(std_out,*)" Psps%lmnmax is ",Psps%lmnmax
 end if
 !
 ! ============================================================================
 ! ==== Evaluate <p_lmn|e^(ikr)U_i> for each k on the k-grid  and each U_i ====
 ! ============================================================================
 !
 ! Here I assume that the G-sphere is gamma-centered.
 ! Real Spherical Harmonics are always used to apply the non-local part even for NC pseudos.
 ! I did non find any easy way to extract only <p_nl|psi> from nonlop_pl.
 !useylm_=1

 ABI_ALLOCATE(Cp_blc2 ,(natom,nspinor))
 call cprj_alloc(Cp_blc2, 0,nlmn_sort)

 if (spline_opt==USE_BSPLINE) then
   ABI_ALLOCATE(Cp_blc1 ,(natom,nspinor))
   call cprj_alloc(Cp_blc1, 0,nlmn_sort)
   ABI_ALLOCATE(Cp_left ,(natom,nspinor))
   call cprj_alloc(Cp_left, 0,nlmn_sort)
   ABI_ALLOCATE(Cp_right,(natom,nspinor))
   call cprj_alloc(Cp_right,0,nlmn_sort)
 end if
 !
 ! ============================
 ! ==== B-spline for Vnl_k ====
 ! ============================
 ! * Tabulate Vnlk_ij on the homogeneous mesh of k-points. 
 ! * useylm is used.
 if (spline_opt==USE_BSPLINE) then 
   ABI_CHECK(Psps%usepaw==0,"PAW with B-splines not coded")
   ABI_CHECK(Psps%useylm==1,"useylm must be 1")
   MSG_ERROR("cprj_bspl_t must be rewritten mind blc_ug!")

   bspl_kdiv = (/4,4,4/) ! TODO should be input.
   !bspl_kdiv = (/8,8,8/) 
   !bspl_kord = (/1,1,1/) 
   bspl_kord = bspl_kdiv ! TODO should be input.

!   call cprj_bspline_init(Cp_bspl,Cryst,Psps,Pawtab,Pawang,nspinor,Wsh%ngfft,bspl_kdiv,bspl_kord,&
!&    npw_k,sh_size,kg_k,blc_ug,nlmn_sort) 
 end if
 ! =======================
 ! ==== Interpolation ====
 ! =======================
 !
 matblk=nloalg(4); if (nloalg(1)>0) matblk=natom
 !
 ! Prepare the kinetic term.
 ABI_ALLOCATE(kdotg,(npw_k))
 ABI_ALLOCATE(half_gsq,(npw_k))
 ABI_ALLOCATE(wsg,(npw_k))

 ! TODO Add new overloaded interface. effmass option!

 do ig=1,npw_k 
   kdotg(ig)    = two_pi**2 * DOT_PRODUCT(k4intp,MATMUL(Cryst%gmet,kg_k(:,ig)))
   half_gsq(ig) = half * vdotw(one*kg_k(:,ig),one*kg_k(:,ig),Cryst%gmet,"G")  
 end do

 ABI_ALLOCATE(vectin,(2, npw_k*nspinor))
 ABI_ALLOCATE(vnl_psi,(2,npw_k*nspinor))
 ABI_ALLOCATE(cvnl_psi,(npw_k*nspinor))

 ABI_ALLOCATE(s_psi,(2,npw_k*nspinor*Psps%usepaw))
 ABI_ALLOCATE(cs_psi,(npw_k*nspinor*Psps%usepaw))

 ! THIS PART IS NEEDED FOR THE CALL TO opernl although some quantities won't be used.
 ! Now I do things cleanly then we try to pass zero-sized arrays!
 ABI_ALLOCATE(ylm_k,(npw_k,Psps%mpsang**2*useylm_))

 if (useylm_==1) then
   kptns_(:,1)=k4intp; optder=0; mkmem_=1
   ABI_ALLOCATE(dum_ylm_gr_k,(npw_k,3+6*(optder/2),Psps%mpsang**2))

   !  Here mband is not used if paral_compil_kpt=0
   call initylmg(Cryst%gprimd,kg_k,kptns_,mkmem_,Wsh%MPI_enreg,Psps%mpsang,npw_k,(/1/),1,&
&    (/npw_k/),1,optder,Cryst%rprimd,0,0,ylm_k,dum_ylm_gr_k)

   ABI_DEALLOCATE(dum_ylm_gr_k)
   istat = ABI_ALLOC_STAT
 end if
 !
 ! Compute (k+G) vectors (only if useylm_=1)
 nkpg=3*nloalg(5)  
 ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
 if (nkpg>0) call mkkpg(kg_k,kpg_k,k4intp,nkpg,npw_k)
 !
 ! ========================================================
 ! ==== Compute nonlocal form factors ffnl at all (k+G) ====
 ! ========================================================
 !
 dimffnl=1+ider0 ! Derivatives are not needed. 
 ABI_ALLOCATE(ffnl,(npw_k,dimffnl,Psps%lmnmax,Psps%ntypat))

 call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,Cryst%gmet,Cryst%gprimd,ider0,idir0,Psps%indlmn,&
&  kg_k,kpg_k,k4intp,Psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,& 
&  Psps%ntypat,Psps%pspso,Psps%qgrid_ff,Cryst%rmet,Psps%usepaw,useylm_,ylm_k,ylmgr_dum)

 ABI_DEALLOCATE(ylm_k)
 !
 ! Allocate the arrays phkxred and ph3d, compute phkxred and eventually ph3d.
 ABI_ALLOCATE(phkxred,(2,Cryst%natom))
 do iat=1,Cryst%natom
   iatom=Cryst%atindx(iat)
   arg=two_pi*DOT_PRODUCT(k4intp,Cryst%xred(:,iat))
   phkxred(1,iatom)=DCOS(arg)
   phkxred(2,iatom)=DSIN(arg)
 end do

 ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
 istat = ABI_ALLOC_STAT
 if (nloalg(1)>0) then ! Allocation as well as precomputation
   if (Wsh%MPI_enreg%mode_para/='b') then
     call ph1d3d(1,Cryst%natom,kg_k,matblk,Cryst%natom,npw_k,Wsh%ngfft(1),Wsh%ngfft(2),Wsh%ngfft(3),phkxred,Wsh%ph1d,ph3d)
   else 
     MSG_ERROR("Stop not coded")
   end if
 end if
!END BORING CODE NEEDED TO CALL opernl
 !
 ! Calculate the upper triangle of <blc2| H_k |blc1>.
 do blc2=1,sh_size  
   ug2 => Wsh%Wave(blc2,k1,1)%ug

   if (spline_opt==EVAL_VNLK) then ! Calculate <G|Vnl|psi> for this k-point
     call wfd_vnlpsi(Wsh,blc2,k1,spin,npw_k,Cryst,Psps,Ham_k,vnl_psi,s_psi,Kext=Kdata)
     cvnl_psi = DCMPLX(vnl_psi(1,:),vnl_psi(2,:))
     if (Psps%usepaw==1) cs_psi = DCMPLX(s_psi(1,:),s_psi(2,:))

   else if (spline_opt==USE_BSPLINE) then
     call cprj_bspline_eval(Cp_bspl,blc2,k4intp,Cp_blc2)
   end if

   do blc1=1,blc2 ! Upper triangle of the hk_ij matrix
     ug1 => Wsh%Wave(blc1,k1,1)%ug

     ! Treat the matrix elements of the Vnl operator.
     if (spline_opt==EVAL_VNLK) then
       vnl_ij  = xdotc(npw_k*nspinor,ug1,1,cvnl_psi,1)

     else if (spline_opt==USE_BSPLINE) then

       call cprj_bspline_eval(Cp_bspl,blc1,k4intp,Cp_blc1)

       ! nonlop_ylm calculates diagonal matrix elements of Vnl while we need 
       ! complex off-diagonal terms as well. Assuming that the pseudotential strengths are
       ! real we can rewrite A_i^* D_{ij} B_j as.
       !
       ! A_i^* D_ij B_j = [rea_i rb_j + ima_i imb_j] D_ij +i [-ima_i reb_j + rea_i imb_j ] D_ij

       cp_idx = RESHAPE( (/1,1,2,2,2,1,1,2/),(/2,4/) )

       do iat=1,Cryst%natom ! already done in cprj_alloc but oh well.
         Cp_left (iat,1)%cp(2,:)=zero
         Cp_right(iat,1)%cp(2,:)=zero
       end do
       !
       ! Non-local part. This coding assumes that Dij are real Spinor not coded! 
       ! TODO Sk_ij for PAW
       do step=1,4    
         ll = cp_idx(1,step)
         rr = cp_idx(2,step)

         do iat=1,Cryst%natom 
           if (ll==2.and.rr==1) then
             Cp_left (iat,1)%cp(1,:) = -Cp_blc1(iat,1)%cp(ll,:) 
           else
             Cp_left (iat,1)%cp(1,:) = Cp_blc1(iat,1)%cp(ll,:) 
           end if
           Cp_right(iat,1)%cp(1,:) = Cp_blc2(iat,1)%cp(rr,:) 
         end do

         signs  =1     ! get contracted elements (energy, forces, stress, ...)
         choice =1     ! => a non-local energy contribution
         cpopt  =2     ! right <p_lmn|in> are already in memory;
         nnlout =1
         paw_opt=0; if (Psps%usepaw==1) paw_opt=4 ! both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
         istwf_k=1
         !paw_opt=Psps%usepaw  !Norm-conserving Vnl (use of Kleinman-Bylander ener.) PAW nonlocal part of H (use of Dij coeffs)
         !paw_opt=1

         ! Calling hacked version of nonlop_ylm TODO implement other cases and sij as well.
         call nonlop_ylm(Cryst%atindx1,choice,cpopt,Cp_right,Ham_k%dimekb1,Ham_k%dimekb2,dimffnl,dimffnl,&
&          Ham_k%ekb,enlout,ffnl,ffnl,Cryst%gprimd,idir0,Psps%indlmn,istwf_k,&
&          kg_k,kg_k,kpg_k,kpg_k,k4intp,k4intp,lambda0,Psps%lmnmax,matblk,Wsh%mgfft,&
&          Wsh%MPI_enreg,Cryst%natom,Cryst%nattyp,Wsh%ngfft,nkpg,nkpg,nloalg,nnlout,&
&          npw_k,npw_k,nspinor,nspinortot,Cryst%ntypat,paw_opt,phkxred,phkxred,Wsh%ph1d,&
&          ph3d,ph3d,signs,Ham_k%sij,vectin,Cryst%ucvol,vectin,vectin,cprjin_left=Cp_left)

         estep(step) = enlout(1)
       end do ! step
       vnl_ij  = DCMPLX( SUM(estep(1:2)), SUM(estep(3:4)) )

       if (blc1==blc2) then ! Calling hacked version of nonlop_ylm TODO implement other cases and sij as well.
         call nonlop_ylm(Cryst%atindx1,choice,cpopt,Cp_blc1,Ham_k%dimekb1,Ham_k%dimekb2,dimffnl,dimffnl,&
&          Ham_k%ekb,enlout,ffnl,ffnl,Cryst%gprimd,idir0,Psps%indlmn,istwf_k,&
&          kg_k,kg_k,kpg_k,kpg_k,k4intp,k4intp,lambda0,Psps%lmnmax,matblk,Wsh%mgfft,&
&          Wsh%MPI_enreg,Cryst%natom,Cryst%nattyp,Wsh%ngfft,nkpg,nkpg,nloalg,nnlout,&
&          npw_k,npw_k,nspinor,nspinortot,Cryst%ntypat,paw_opt,phkxred,phkxred,Wsh%ph1d,&
&          ph3d,ph3d,signs,Ham_k%sij,vectin,Cryst%ucvol,vectin,vectin)

         vnl_ij= DCMPLX( enlout(1), zero)
       end if

       if(blc1==blc2) write(777,*)blc1,blc2,vnl_ij
     end if
     !
     ! ====================================
     ! ==== Assemble final Hamiltonian ====
     ! ====================================
     !
     ! Kinetic energy.
     wsg = kdotg*ug2
     kin_ij = xdotc(npw_k,ug1,1,wsg,1) 

     wsg = half_gsq*ug2
     kin_ij = kin_ij + xdotc(npw_k,ug1,1,wsg,1) 
     if (blc1==blc2) kin_ij = kin_ij + half * vdotw(k4intp,k4intp,Cryst%gmet,"G")
     !
     ! Total Hamiltonian.
     hk_ij(blc1,blc2) = kin_ij + vloc_ij(blc1,blc2) + vnl_ij
     !
     ! PAW Overlap operator.
     if (Psps%usepaw==1) sk_ij(blc1,blc2) = xdotc(npw_k*nspinor,ug1,1,cs_psi,1)
   end do ! blc1
 end do ! blc2

!DEALLOCATE BORING STUFF
 ABI_DEALLOCATE(ffnl)
 ABI_DEALLOCATE(kpg_k)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(phkxred)
 ABI_DEALLOCATE(ph3d)
!END BORING STUFF

 ABI_DEALLOCATE(kdotg)
 ABI_DEALLOCATE(half_gsq)
 ABI_DEALLOCATE(wsg)
 ABI_DEALLOCATE(vectin)
 ABI_DEALLOCATE(vnl_psi)
 ABI_DEALLOCATE(cvnl_psi)
 ABI_DEALLOCATE(s_psi)
 ABI_DEALLOCATE(cs_psi)
 istat = ABI_ALLOC_STAT

 call cprj_free(Cp_blc2 )
 ABI_DEALLOCATE(Cp_blc2)
 call kdata_free(Kdata)

 if (spline_opt==USE_BSPLINE) then
   call cprj_bspline_free(Cp_bspl)
   call cprj_free(Cp_blc1 )
   ABI_DEALLOCATE(Cp_blc1)
   call cprj_free(Cp_left )
   ABI_DEALLOCATE(Cp_left)
   call cprj_free(Cp_right)
   ABI_DEALLOCATE(Cp_right)
 end if

 ABI_DEALLOCATE(nlmn_sort)

 DBG_EXIT("COLL")

 RETURN

 ABI_UNUSED(Pawang%ngnt)
 ABI_UNUSED(Pawtab(1)%basis_size)

end subroutine shirley_hks       
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/shirley_interp
!! NAME
!! shirley_interp
!!
!! FUNCTION
!!
!! INPUTS
!! jobz
!! Dtset
!! Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang<pawang_type> angular mesh discretization and related data:
!! Pawrhoij
!! Paw_ij(natom)<type(paw_ij_type)>=data structure containing PAW arrays given on (i,j) channels.
!! spline_opt
!! sh_coverage
!! ngfftf(18)=Information about the dense 3D FFT used for vtrial.
!! nfftf=Number of points in the FFT grid in vtrial. Might differ from the FFT mesh used for the wavefunctions.
!! vtrial(nfftf,nspden)= trial potential (Hartree)
!! [intp_wtk]=Used to define the k-point Weights when BSt_intp is present.
!!
!! OUTPUT
!!  [BSt_intp]
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      exc_interp_ham,m_shexc
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine shirley_interp(Wsh,jobz,Dtset,Cryst,Psps,Pawtab,Pawfgr,Pawang,Pawrad,&
&  Pawrhoij,Paw_ij,ngfftc,ngfftf,nfftf,vtrial,spline_opt,&
&  intp_nband,intp_mband,intp_nk,intp_kpt,intp_ene,KS_intp,intp_wtk,BSt_intp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shirley_interp'
 use interfaces_14_hidewrite
 use interfaces_44_abitypes_defs
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,spline_opt,intp_nk,intp_mband
 character(len=*),intent(in) :: jobz
 type(wfs_descriptor),intent(inout) :: Wsh
 type(Crystal_structure),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(dataset_type),intent(in) :: Dtset
 type(Pawfgr_type),intent(in) :: Pawfgr
 type(bandstructure_type),optional,intent(out) :: BSt_intp
!arrays
 integer,intent(in) :: ngfftf(18),ngfftc(18)
 integer,intent(in) :: intp_nband(intp_nk,Wsh%nsppol)
 real(dp),intent(in) :: vtrial(nfftf,Wsh%nspden),intp_kpt(3,intp_nk)
 real(dp),optional,intent(in) :: intp_wtk(intp_nk)
 real(dp),intent(out) :: intp_ene(intp_mband,intp_nk,Wsh%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wsh%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Wsh%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wsh%usepaw)
 type(pawrhoij_type),intent(in) :: Pawrhoij(Cryst%natom*Wsh%usepaw)
 type(ksintp_t),intent(inout) :: KS_intp(intp_nk,Wsh%nsppol)

!Local variables ------------------------------
!scalars
 integer,parameter :: EVAL_VNLK=0,USE_BSPLINE=1
 integer,parameter :: istwf1=1,k1=1
 integer :: ii,istat,ib,jj,ierr 
 integer :: nband_k,ikpt,natom,matblk,nefound
 integer :: sh1,sh2,ldz,prtvol,pawprtvol,istwf_k
 integer :: mgfftc,nfftc,onband_diago,comm,paral_kgb,usepaw,spin
 integer :: npw_k,iat,nspinor,nsppol,nspden,sh_size,intp_bantot
 real(dp) :: fft_fact,ene_fact
 !logical,parameter :: debug_with_diago=.TRUE.
 logical,parameter :: debug_with_diago=.FALSE.
 logical :: want_eigenvectors,do_full_diago
 character(len=500) :: msg,frmt1,frmt2
 type(ddiago_ctl_type) :: Diago_ctl
 type(gs_hamiltonian_type) :: Ham_k
 type(dst_t) :: Dst
!arrays
 !integer :: got(Wsh%nproc)
 !integer,allocatable :: bbp_ks_distrb(:,:,:,:)
 integer :: nloalg(5),bspl_kdiv(3),bspl_kord(3)
 integer,allocatable :: nlmn_sort(:)
 integer,allocatable :: intp_istwfk(:),intp_npwarr(:)
 real(dp) :: kpoint(3) !,kk1(3),kk2(3)
 real(dp),pointer :: diag_ene(:),diag_vec(:,:,:)
 real(dp),allocatable :: enek_ij(:),vnl_psi(:,:),opaw_psi(:,:)
 real(dp),allocatable :: copy_vtrial(:,:)
 real(dp),allocatable :: intp_doccde(:),intp_occ(:),ugly_ene(:)
 complex(gwpc),allocatable :: ur1(:),ur2(:),vloc_psi(:)
 complex(dpc),allocatable :: hk_ij(:,:),sk_ij(:,:),eig_vec(:,:),vloc_ij(:,:)
 character(len=10) :: spin_name(2)
 type(cprj_type),pointer :: diag_Cprj(:,:)
 type(cprj_type),allocatable :: Cp_sh1(:,:) 
 type(cprj_type),allocatable,target :: Cp_left(:,:),Cp_right(:,:) 
 type(cprj_bspl_t) :: Cp_bspl

!************************************************************************

 DBG_ENTER("COLL")

 MSG_WARNING("Calling experimental code shirley_interp!")

 ABI_UNUSED(Pawrhoij(1)%cplex)
 ABI_UNUSED(Pawrad(1)%mesh_size)

 ABI_CHECK(Wsh%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(Wsh%paral_kgb==0,"paral_kgb/=0 not coded")
 ABI_CHECK(Wsh%rfft_is_symok,"FFT not symmetric in real space")

 nspinor   = Wsh%nspinor
 nsppol    = Wsh%nsppol
 nspden    = Wsh%nspden
 nloalg    = Wsh%nloalg
 paral_kgb = Wsh%paral_kgb
 usepaw    = Wsh%usepaw
 natom     = Wsh%natom

 prtvol    = Wsh%prtvol
 pawprtvol = Wsh%pawprtvol

 ! The coarse mes used for the Hamiltonian.
 nfftc  = PRODUCT(ngfftc(1:3)) 
 mgfftc = MAXVAL(ngfftc(1:3))

 if (nsppol==1) spin_name=(/'          ','          '/)
 if (nsppol==2) spin_name=(/'SPIN_UP   ','SPIN_DOWN '/)

 call ks_intp_nullify(KS_intp)

 want_eigenvectors = starts_with(jobz,(/"V"/),csens=.FALSE.)        
 if (want_eigenvectors) then
   call wrtout(std_out," Eigenvectors will be calculated." ,"COLL")
 else 
   call wrtout(std_out," Interpolating energies only. Eigenvectors won't be calculated.","COLL")
 end if

!BEGIN DEBUG: Test on the orthogonalization of the input wavefunctions. 
 call wfd_reset_ur_cprj(Wsh)
 call wfd_test_ortho(Wsh,Cryst,Pawtab,unit=std_out,mode_paral="COLL")
!END DEBUG

 select case (spline_opt) 
 case (EVAL_VNLK)
   call wrtout(std_out," Vnlk_ij will be evaluated on the interpolating k-mesh.","COLL")
   !
 case (USE_BSPLINE)
   call wrtout(std_out," Vnlk_ij will be obtained on the interpolating k-mesh via B-spline functions.","COLL")
   ABI_CHECK(Psps%useylm==1,"B-spline technique requires useylm==1") ! Not completely true but it makes life easier.
   !
   ! FIXME This trick is needed for using the B-spline interpolation when NC is used.
   ABI_ALLOCATE(nlmn_sort,(natom))
   iat=0 ! nlmn dims sorted by atom type.
   if (usepaw==1) then
     nlmn_sort = Wsh%nlmn_sort
   else  ! FIXME here lmnmax == lnmax if useylm_==0
     nlmn_sort = 9
     if (Psps%useylm==1) nlmn_sort=Psps%lmnmax
     MSG_ERROR("useylm==0 not coded")
     write(std_out,*)" Psps%lmnmax is ",Psps%lmnmax
   end if
   !
 case default
   write(msg,'(a,i0)')" Wrong value for spline_opt:",spline_opt
   MSG_ERROR(msg)
 end select
 !
 ! =======================
 ! ==== Interpolation ====
 ! =======================
 !
 ! vtrial might be given on a FFT mesh that is denser than the FFT used for Wsh.
 ! If the two FFTs differ, change the mesh for the wavefunctions.
 ! Another possibility would be moving vtrial from the dense to the coarse mesh.
 call wfd_change_ngfft(Wsh,Cryst,Psps,ngfftf) 

 fft_fact = one/Wsh%nfft
 ABI_ALLOCATE(ur1,(Wsh%nfft*nspinor))
 ABI_ALLOCATE(ur2,(Wsh%nfft*nspinor))

 intp_ene=zero

 do spin=1,nsppol

   sh_size = Wsh%nband(k1,spin)
   !
   ! Precompute <sh1| vloc_spin |sh2> as it is not k-dependent.
   ABI_ALLOCATE(vloc_ij,(sh_size,sh_size))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"Out of memory in vloc_ij")
                                                                                        
   ABI_ALLOCATE(vloc_psi,(Wsh%nfft*nspinor))
   !
   ! The upper triangle of the <sh2|vloc|sh1> matrix.
   ! BE careful here as |\tpsi> is not normalized.
   do sh2=1,sh_size
     call wfd_get_ur(Wsh,sh2,k1,spin,ur2)
     if (nspinor==1) then
       vloc_psi = ur2*vtrial(:,spin)
     else 
       MSG_ERROR("vloc_psi doesn't support nspinor==2")
     end if
     !
     ! Diagonal element.
     vloc_ij(sh2,sh2) = xdotc(Wsh%nfft,ur2,1,vloc_psi,1)*fft_fact 
     !
     ! Upper triangle.                                                                     
     do sh1=1,sh2-1 
       call wfd_get_ur(Wsh,sh1,k1,spin,ur1)
       vloc_ij(sh1,sh2) = xdotc(Wsh%nfft,ur1,1,vloc_psi,1)*fft_fact
     end do
     !
   end do
                                                                                        
   ABI_DEALLOCATE(vloc_psi)
   !
   ! ============================================================================
   ! ==== Evaluate <p_lmn|e^(ikr)U_i> for each k on the k-grid  and each U_i ====
   ! ============================================================================
   !
   ! Here I assume that the G-sphere is gamma-centered.
   ! Real Spherical Harmonics are always used to apply the non-local part even for NC pseudos.
   ! I did non find any easy way to extract only <p_nl|psi> from nonlop_pl.

   if (spline_opt==USE_BSPLINE) then
     ABI_ALLOCATE(Cp_sh1 ,(natom,nspinor))
     call cprj_alloc(Cp_sh1, 0,nlmn_sort)
     ABI_ALLOCATE(Cp_left ,(natom,nspinor))
     ABI_ALLOCATE(Cp_right,(natom,nspinor))
     call cprj_alloc(Cp_left, 0,nlmn_sort)
     call cprj_alloc(Cp_right,0,nlmn_sort)
   end if
   !
   ! ==============================
   ! ==== Prepare NL strengths ====
   ! ==============================
   call init_hamiltonian(Ham_k,Psps,Paw_ij,Pawtab,nspinor,nspden,natom,&
&    Cryst%ntypat,Cryst%typat,Cryst%xred,Wsh%nfft,Wsh%mgfft,Wsh%ngfft,Cryst%rprimd,nloalg)
   !
   ! ============================
   ! ==== B-spline for Vnl_k ====
   ! ============================
   ! * Tabulate Vnlk_ij on the homogeneous mesh of k-points. 
   ! * useylm is used.
   if (spline_opt==USE_BSPLINE) then 
     ABI_CHECK(usepaw==0,"PAW with B-splines not coded")
     ABI_CHECK(Psps%useylm==1,"useylm must be 1")

     bspl_kdiv = (/4,4,4/) ! TODO should be input.
     bspl_kord = bspl_kdiv ! TODO should be input.
     MSG_ERROR("cprj_bspl_init must be rewritten mind sh_ug!")

     ! TODO this part has to be rewritten.
     !     call cprj_bspline_init(Cp_bspl,Cryst,Psps,Pawtab,Pawang,nspinor,Wsh%ngfft,bspl_kdiv,bspl_kord,&
     !&      npw_gamma,sh_size,kg_gamma,sh_ug,nlmn_sort) 
   end if
   !
   ! =============================================
   ! ==== Loop over the interpolated k-points ====
   ! =============================================
   !
   ABI_ALLOCATE(hk_ij,(sh_size,sh_size))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"Out of memory in hk_ij")

   ABI_ALLOCATE(sk_ij,(sh_size,sh_size*usepaw))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"Out of memory in sk_ij")

   matblk=nloalg(4); if (nloalg(1)>0) matblk=natom

   npw_k = Wsh%Kdata(k1)%npw

   ABI_ALLOCATE(vnl_psi,(2,npw_k*nspinor))
   ABI_ALLOCATE(opaw_psi,(2,npw_k*nspinor*usepaw))

   do ikpt=1,intp_nk ! Loop over the K-points for the interpolation (NB: k will be wrapped in ]-1/2,1/2]).

     istwf_k = 1 ! no time-reversal tricks!
     kpoint  = intp_kpt(:,ikpt) 
     nband_k = intp_nband(ikpt,spin) 

     call shirley_hks(Wsh,kpoint,spin,spline_opt,Ham_k,Cryst,Psps,Pawtab,Pawang,Paw_ij,sh_size,vloc_ij,hk_ij,sk_ij)
     !
     ! Diagonalize Hk_ij in the optimal Bloch supspace.
     ABI_ALLOCATE(enek_ij,(sh_size))

     do_full_diago = (nband_k == sh_size)
     if (.not.do_full_diago) then
       ldz=1; if (want_eigenvectors) ldz=sh_size
       ABI_ALLOCATE(eig_vec,(ldz,nband_k))
       istat = ABI_ALLOC_STAT
       if (istat/=0) then ! Try to continue.
         MSG_WARNING("Allocation of eig_vec failed. Full diago will be performed")
         do_full_diago = .TRUE.
       end if
     end if

     if (usepaw==0) then  ! Solve H*v = e*v
       if (do_full_diago) then
         call xheev(jobz,"Upper",sh_size,hk_ij,enek_ij) 
       else
         call xheevx(jobz,"Irange","Upper",sh_size,hk_ij,zero,zero,1,nband_k,-tol8,nefound,enek_ij,eig_vec,ldz)
         if (want_eigenvectors) hk_ij(:,1:ldz) = eig_vec
       end if
     else                     ! Solve H*v = e*S*v
       if (do_full_diago) then
         call xhegv(1,jobz,"Upper",sh_size,hk_ij,sk_ij,enek_ij) 
       else
         call xhegvx(1,jobz,"Irange","Upper",sh_size,hk_ij,sk_ij,zero,zero,1,nband_k,-tol8,nefound,enek_ij,eig_vec,ldz)
         if (want_eigenvectors) hk_ij(:,1:ldz) = eig_vec
       end if
     end if
     !
     ! Store the interpolated eigenvalues and eigenstates in the optimal basis set.
     if (want_eigenvectors) then 
       call ks_intp_init(KS_intp(ikpt,spin),sh_size,nband_k,enek_ij,hk_ij,istat)
       ABI_CHECK(istat==0,"out of memory %obloch")
       ! TODO Write unitary transformation for this spin on file 
     end if

     if (allocated(eig_vec))  then
       ABI_DEALLOCATE(eig_vec)
     end if

     if (.TRUE..or.prtvol>0) then ! Write interpolated energies.
       ene_fact=Ha_eV; frmt1='(i4,4x,9(1x,f7.4))'; frmt2='(8x,9(1x,f7.4))'
       write(msg,'(a,3es16.8,2a)')' Eigenvalues in eV for kpt= ( ',kpoint," ), spin ",spin_name(spin)
       call wrtout(std_out,msg,'COLL')

       write(msg,frmt1)ikpt,(enek_ij(ib)*ene_fact,ib=1,MIN(9,nband_k))
       call wrtout(std_out,msg,'COLL') 

       if (nband_k>9) then
         do jj=10,nband_k,9
           write(msg,frmt2) (enek_ij(ib)*ene_fact,ib=jj,MIN(jj+8,nband_k))
           call wrtout(std_out,msg,'COLL') 
         end do
       end if
       call flush_unit(std_out)
     end if
     !
     ! Save the interpolated energies.
     intp_ene(:,ikpt,spin) = enek_ij(1:nband_k)
     !
     if (debug_with_diago) then ! Compare interpolated energies wrt direct diago results.

       ! FIXME here there is a problem with Wd%ecut and Dtset%ecut.
       !if (ikpt==1) write(std_out,*)" CHECK: Dtset%ecut=",Dtset%ecut," Wsh%ecut= ",Wsh%ecut 

       !call init_ddiago_ctl(Diago_ctl,"No Vectors",spin,nspinor,Wsh%ecut,kpoint,nloalg,Cryst%gmet,&
       !&   nband_k=nband_k,effmass=Dtset%effmass,istwf_k=istwf1,prtvol=prtvol)

       call init_ddiago_ctl(Diago_ctl,"No Vectors",spin,nspinor,Dtset%ecut,kpoint,nloalg,Cryst%gmet,&
       &   nband_k=nband_k,effmass=Dtset%effmass,istwf_k=istwf1,prtvol=prtvol)

       nullify(diag_ene)
       nullify(diag_vec)
       nullify(diag_Cprj)

       ! TODO check this part
       comm = xmpi_self
       ABI_ALLOCATE(copy_vtrial,(nfftf,nspden))
       copy_vtrial = vtrial ! TO avoid intent inout

       call ks_ddiago(Diago_ctl,nband_k,nfftc,mgfftc,ngfftc,natom,&
&        Cryst%typat,nfftf,nspinor,nspden,nsppol,Cryst%ntypat,Pawtab,Pawfgr,Paw_ij,&
&        Psps,Cryst%rprimd,copy_vtrial,Cryst%xred,onband_diago,diag_ene,diag_vec,diag_Cprj,comm,ierr)

       ABI_DEALLOCATE(copy_vtrial)

       ABI_CHECK(ierr==0,"Fatal error. Cannot continue!")

       !% call slk_hgg(Hgg,Tgg,Evec,Diago_ctl,Cryst,nband_k,nfftc,mgfftc,Wsh%ngfft,nfftf,&
       !% & nspinor,nspden,Pawtab,Pawfgr,Paw_ij,Psps,vtrial,onband_diago,diago_ene,comm)

       if (.TRUE..or.prtvol>0) then
         ene_fact=Ha_eV
         write(msg,'(a,3es16.8,2a)')' Diago Eigenvalues in eV for kpt= ( ',kpoint," ), spin ",spin_name(spin)
         call wrtout(std_out,msg,'COLL')

         write(msg,'(i4,4x,9(1x,f7.4))')ikpt,(ene_fact*diag_ene(ii),ii=1,MIN(9,nband_k))
         call wrtout(std_out,msg,'COLL') 

         if (nband_k>9) then
           do jj=10,nband_k,9
             write(msg,'(8x,9(1x,f7.4))')(ene_fact*diag_ene(ii),ii=jj,MIN(jj+8,nband_k))
             call wrtout(std_out,msg,'COLL') 
           end do
         end if
       end if
       
       !call dst_init(Dst,nband_k-2,ABS(diag_ene(1:nband_k-2) - enek_ij(1:nband_k-2))*Ha_eV*1000)
       call dst_init(Dst,nband_k-2,ABS(diag_ene(1:nband_k-2) - enek_ij(1:nband_k-2))*Ha_eV*1000)

       write(std_out,'(a,f7.2,a,i0,a,i0)')&
       & " MAX abs error diag_ene - intp_ene ",MAXVAL(ABS(diag_ene(1:nband_k) - enek_ij(1:nband_k)) )*Ha_eV*1000,&
       & " [meV| for band ",imax_loc(ABS(diag_ene(1:nband_k) - enek_ij(1:nband_k)) ),"/",nband_k

       write(std_out,'(a,4es9.1,a)')" min, Max, mean, stdev ",Dst%min, Dst%max, Dst%mean, Dst%stdev, " [meV]."
       call flush_unit(std_out)

       if (associated(diag_ene))  then
         ABI_DEALLOCATE(diag_ene)
       end if
       if (associated(diag_vec))  then
         ABI_DEALLOCATE(diag_vec)
       end if
       if (associated(diag_Cprj)) then
       end if
     end if

     ABI_DEALLOCATE(enek_ij)
   end do ! intp_kpt
   !
   ! Interpolation completed, deallocate memory.
   ABI_DEALLOCATE(hk_ij)
   ABI_DEALLOCATE(vloc_ij)
   ABI_DEALLOCATE(sk_ij)
   istat = ABI_ALLOC_STAT

   ABI_DEALLOCATE(vnl_psi)
   ABI_DEALLOCATE(opaw_psi)
   call destroy_hamiltonian(Ham_k)       

   if (spline_opt==USE_BSPLINE) then
     call cprj_bspline_free(Cp_bspl)
     call cprj_free(Cp_sh1 )
     ABI_DEALLOCATE(Cp_sh1)
     call cprj_free(Cp_left )
     ABI_DEALLOCATE(Cp_left)
     call cprj_free(Cp_right)
     ABI_DEALLOCATE(Cp_right)
   end if
   !
 end do ! spin
 !
 ! Update the energies in io_Bst.
 call xsum_mpi(intp_ene,Wsh%comm,ierr)
 !
 ! Free memory
 ABI_DEALLOCATE(ur1)
 ABI_DEALLOCATE(ur2)
 if (allocated(nlmn_sort))  then
   ABI_DEALLOCATE(nlmn_sort)
 end if

 ! Optionally create Bandstructure_type for storing the interpolated energies.
 if (PRESENT(BSt_intp)) then
   !
   if (.not.PRESENT(intp_wtk)) then
     MSG_ERROR("intp_wtk must be present when BSt_intp is wanted!") 
   end if
   intp_bantot=SUM(intp_nband)

   ABI_ALLOCATE(intp_doccde,(intp_bantot))
   intp_doccde =zero
   ABI_ALLOCATE(intp_occ,(intp_bantot))
   intp_occ    =zero 
   ABI_ALLOCATE(intp_istwfk,(intp_nk))
   intp_istwfk =1
   ABI_ALLOCATE(intp_npwarr,(intp_nk))
   intp_npwarr =Wsh%npwarr(1) ! Meaningless

   ! Have to reshape intp_ene --> ugly_ene due to the ugly convention used in abinit to store energies.
   ABI_ALLOCATE(ugly_ene,(intp_bantot))
   ugly_ene(:)=HUGE(zero)

   call pack_eneocc(intp_nk,nsppol,intp_mband,intp_nband,intp_bantot,intp_ene,ugly_ene)

   call bstruct_init(intp_bantot,BSt_intp,Dtset%nelect,intp_doccde,ugly_ene,intp_istwfk,intp_kpt,&
&   intp_nband,intp_nk,intp_npwarr,nsppol,nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,intp_occ,intp_wtk)

   ABI_DEALLOCATE(intp_doccde)
   ABI_DEALLOCATE(ugly_ene)
   ABI_DEALLOCATE(intp_occ)
   ABI_DEALLOCATE(intp_istwfk)
   ABI_DEALLOCATE(intp_npwarr)

   !%call update_occ(BSt_intp,Dtset%fixmom,Dtset%stmbias)      ??
 end if

 DBG_EXIT("COLL")

end subroutine shirley_interp       
!!***

!----------------------------------------------------------------------

!!****f* m_shirley/wfd_shirley_to_eh
!! NAME
!! wfd_shirley_to_eh
!!
!! FUNCTION
!!  Returns a new wavefunction descriptor containing the basis set for the e-h manifold.
!!
!! INPUT
!! Wfd<wfs_descriptor>
!! Cryst<crystal_structure>= data type gathering info on symmetries and unit cell
!! Psps<pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat)<pawtab_type>=paw tabulated starting data
!! Pawang <pawang_type>=angular mesh discretization and related data:
!! eh_coverage
!!
!! OUTPUT
!!  Weh<wfs_descriptor>
!!
!! PARENTS
!!      m_shexc
!!
!! CHILDREN
!!      blas_cholesky_ortho,fft_ur,flush_unit,fourdp_c2c_ip,get_kg,kgindex
!!      ovlp_diago_and_prune,ovlp_free,ovlp_nullify,timein,wfd_change_ngfft
!!      wfd_get_ur,wfd_init,wfd_print,wfd_push_ug,wfd_test_ortho,wrtout,xgemm
!!
!! SOURCE

subroutine wfd_shirley_to_eh(Wsh,Cryst,Psps,Pawtab,Pawang,Pawrad,min_bsize,eh_coverage,Weh,sh2eh)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_shirley_to_eh'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: min_bsize
 real(dp),intent(in) :: eh_coverage
 type(Crystal_structure),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wsh
 type(wfs_descriptor),intent(out) :: Weh
!arrays
 !integer,intent(in) :: ov_ngfft(18)
 complex(gwpc),pointer :: sh2eh(:,:)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wsh%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wsh%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: istwf1=1,k1=1,s1=1,nkpt1=1
 integer :: istat,ij !,ierr
 integer :: band1,band2,natom,eh,midx,ig,mband,nstates,ovlp_size,comm
 integer :: npw_gamma,paral_kgb,usepaw !,spin !,row,col,useylm_
 integer :: fft_idx,nspinor,eh_nkibz,eh_mband,npw,col,band3,band4
 integer :: nsppol
 real(dp) :: fft_fact,norm1 !sqrt_norm1,sqrt_norm2
 real(dp) :: cpu_time,wall_time,cpu0,wall0
 character(len=500) :: msg
!arrays
 integer :: ov_ngfft(18)
 integer :: eh_istwfk(1),eh_size,base
 !integer :: got(Wsh%nproc)
 !integer,allocatable :: bbp_ks_distrb(:,:,:,:)
 integer,pointer :: kg_gamma(:,:) !,gbound_k(:,:),kg_k(:,:),
 integer,allocatable :: eh_nband(:,:) 
 integer,pointer :: igfft0(:)
 real(dp) :: eh_kibz(3,1),gamma_point(3)=(/zero,zero,zero/)
 !real(dp) :: pawovlp(2)
 complex(dpc) :: cfft_fact
 complex(dpc),allocatable :: dpc_tmp(:)
 complex(gwpc),allocatable :: ur1(:),ur2(:),ur12(:),cf_ovlp(:,:)
 complex(gwpc),target,allocatable :: ur3(:),ur4(:) !ur1(:),ur2(:),
 !complex(gwpc),pointer :: pt_ur1(:),pt_ur2(:) 
 complex(gwpc),pointer :: pt_ur3(:),pt_ur4(:)
 complex(dpc),allocatable :: ur34_big(:,:)
 !complex(gwpc),pointer :: ug1(:)
 complex(gwpc),allocatable :: ehg(:),eh_ur(:,:),eh_ug(:,:),gwpc_tmp(:)
 logical,allocatable :: eh_keep_ur(:,:,:),eh_bks_mask(:,:,:),kg_mask(:)
 type(ovlp_t) :: Oeh

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wsh%nspinor==1,"nspinor==2 not coded.")
 ABI_CHECK(Wsh%paral_kgb==0,"paral_kgb/=0 not coded.")
 ABI_CHECK(Wsh%rfft_is_symok,"Real space FFT is not symmetric.")
 ABI_CHECK(Wsh%nsppol==1,"nsppol must be 1!")

 ABI_UNUSED(Pawang%l_max)
 ABI_UNUSED(Pawrad(1)%mesh_size)

 nspinor   = Wsh%nspinor
 paral_kgb = Wsh%paral_kgb
 usepaw    = Wsh%usepaw
 natom     = Wsh%natom
 nsppol    = Wsh%nsppol
 mband     = Wsh%mband
 nstates   = Wsh%nband(k1,s1)
 comm      = Wsh%comm
 ABI_CHECK(xcomm_size(comm)==1,"ovlp2_init is not parallelized")

 write(msg,'(a,f9.6)')" Transformation Shirley --> e-h basis set with eh_coverage : ",eh_coverage
 call wrtout(std_out,msg,"COLL")

 ! 1) Get the overlap matrix <S_i^* S_j|S_k^* S_l> for this spin. 
 ! The Overlap is calculated using the coarse real space FFT ov_ngfft
 ov_ngfft = Wsh%ngfft
 call wfd_change_ngfft(Wsh,Cryst,Psps,ov_ngfft)

 call ovlp_nullify(Oeh)

 ovlp_size   = nstates**2
 Oeh%mband   = mband
 Oeh%nkpt    = nkpt1
 Oeh%min_ene = smallest_real
 Oeh%max_ene = greatest_real  
 !
 ! The size of the overlap matrix and useful tables.
 ABI_ALLOCATE(Oeh%bk2idx,(mband,nkpt1))
 Oeh%bk2idx=0
                                                                                                                                 
 Oeh%size = ovlp_size
 ABI_ALLOCATE(Oeh%idx2bk,(2,ovlp_size))
 Oeh%idx2bk = 0
 !
 ! Allocate overlap matrix. Could use packed matrix to save memory, but Lapack call is slower.
 ABI_ALLOCATE(Oeh%mat,(ovlp_size,ovlp_size))
 istat = ABI_ALLOC_STAT
 if (istat/=0) then
   write(msg,'(a,f8.2,a)')" out of memory in Oeh%mat, requiring :",two*dpc*ovlp_size**2*b2Gb," Gb"
   MSG_ERROR(msg)
 end if
 !Oeh%mat = -HUGE(one)
 write(msg,'(a,f12.1,a,i0)')" Memory required for the overlap matrix: ",two*dpc*ovlp_size**2*b2Mb," Mb; Matrix size= ",ovlp_size
 call wrtout(std_out,msg,"COLL")
 call flush_unit(std_out)
 !
 ! Calculate the overlap matrix  --------------------------------------------------------------------------
 ! 1) Symmetrization in k-space is not needed here 
 ! 2) Matrix is Hermitian.
 !
 fft_fact = one/Wsh%nfft
 ABI_ALLOCATE(ur3,(Wsh%nfft*nspinor))
 ABI_ALLOCATE(ur4,(Wsh%nfft*nspinor))

 ! TODO Temporary implementation used to speed up this part.
 ABI_ALLOCATE(ur34_big,(Wsh%nfft,nstates**2))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory ur34_big")

 Oeh%mat_type = TYPE_OVERLAP
 npw     = Wsh%npwarr(k1)

 call timein(cpu0,wall0)

 col = 0
 do band4=1,nstates
   !
   if (wfd_ihave_ur(Wsh,band4,k1,s1,how="Stored")) then
     pt_ur4 =>  Wsh%Wave(band4,k1,s1)%ur
   else 
     call wfd_get_ur(Wsh,band4,k1,s1,ur4)
     pt_ur4 =>  ur4 
   end if
   !
   do band3=1,nstates
     col = col+1
     if (wfd_ihave_ur(Wsh,band3,k1,s1,how="Stored")) then
       pt_ur3 =>  Wsh%Wave(band3,k1,s1)%ur
     else 
       call wfd_get_ur(Wsh,band3,k1,s1,ur3)
       pt_ur3 =>  ur3
     end if
     !if (Wsh%usepaw==1) call wfd_get_cprj(Wsh,band3,k1,s1,Cryst,Cp_k3,sorted=.FALSE.)
     !
     ur34_big(:,col) = CONJG(pt_ur3) * pt_ur4
     Oeh%idx2bk(1,col) = band3
     Oeh%idx2bk(2,col) = band4
   end do
 end do

 cfft_fact = cone/Wsh%nfft
 !Oeh%mat = cfft_fact * MATMUL( CONJG(TRANSPOSE(ur34_big)), ur34_big)
 call xgemm("C","N",nstates**2,nstates**2,Wsh%nfft,cfft_fact,ur34_big,Wsh%nfft,ur34_big,Wsh%nfft,czero,Oeh%mat,nstates**2)

 !call print_arr(Oeh%mat,max_r=10,max_c=12,unit=std_out,mode_paral="COLL")

 ABI_DEALLOCATE(ur3)
 ABI_DEALLOCATE(ur4)

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0
 write(std_out,*)" Ovlp2 build cpu_time, wall_time ",cpu_time, wall_time
 call flush_unit(std_out)

 ! 2) Diagonalize the overlap matrix selecting the optimal subspace: [ base(spin):ovlp_size ]
 call ovlp_diago_and_prune(Oeh,eh_coverage,eh_size,base) ! In exit Oeh%mat stores the eigenvectors.
 !
 ! Make sure we have enough states.
 if (eh_size < min_bsize) then
   if (Oeh%size<min_bsize) then
     write(msg,'(2(a,i0),2a)')&
&      " Overlap size is ",Oeh%size," whereas min_bsize is ",min_bsize,ch10,&
&      " Decrease the number of bands to be interpolated or increase the number of ab-initio input states."
     MSG_ERROR(msg)
   end if
   eh_size = min_bsize
   write(msg,'(a,2i0)')" Had to enlarge Shirley subspace since input eh_size < min_bsize: ",eh_size,min_bsize
   MSG_COMMENT(msg)
 end if

 call timein(cpu0,wall0)
 !
 ! 3) Init a new wavefunction descriptor to store the optimal basis set.
 !    *) Weh must be allocated here since eh_size is know only after the diagonalization of the overlap.
 !    *) Keep the optimal wavefunctions on each node (if possible) to facilitate the interpolation over the fine k-mesh. 
 !    *) Use Gamma-centered basis set to facilitate the operations in reciprocal space.
 !    *) The new basis is orthogonal, but not normalized since <U_i|U_j> = delta_ij e_i.
 !
 ! The optimal basis set is given on the gamma centered basis set with istwfk==1.
 ! FIXME temporary hacking. There is a bug somewhere in kpgsph
 call get_kg(gamma_point,istwf1,Wsh%ecut,Cryst%gmet,npw_gamma,kg_gamma) 
 !call get_kg(gamma_point,istwf1,14.0_dp,Cryst%gmet,npw_gamma,kg_gamma) 
 !
 ! * Index of the G-sphere in the FFT box.
 ABI_ALLOCATE(igfft0,(npw_gamma))
 ABI_ALLOCATE(kg_mask,(npw_gamma))
 call kgindex(igfft0,kg_gamma,kg_mask,Wsh%MPI_enreg,Wsh%ngfft,npw_gamma)
                                                                      
 ABI_CHECK(ALL(kg_mask),"FFT para not yet implemented")
 ABI_DEALLOCATE(kg_mask)

 eh_istwfk    = istwf1
 eh_nkibz     = 1 
 eh_kibz(:,1) = gamma_point

 ! TODO: BE careful in parallel when nsppol==2. I should recreate the communicators.
 ABI_ALLOCATE(eh_nband,(eh_nkibz,nsppol))
 eh_nband = eh_size
 !call wfd_change_ngfft(Wsh,Cryst,Psps,ov_ngfft)
 
 eh_mband=MAXVAL(eh_nband)      
 ABI_ALLOCATE(eh_bks_mask,(eh_mband,eh_nkibz,nsppol))
 eh_bks_mask=.TRUE.
 ABI_ALLOCATE(eh_keep_ur ,(eh_mband,eh_nkibz,nsppol))
 eh_keep_ur =.TRUE.

 call wfd_init(Weh,Cryst,Pawtab,Psps,eh_keep_ur,Wsh%paral_kgb,npw_gamma,eh_mband,eh_nband,eh_nkibz,nsppol,&
&  eh_bks_mask,Wsh%nspden,nspinor,Wsh%ecutsm,Wsh%dilatmx,eh_istwfk,eh_kibz,Wsh%ngfft,kg_gamma,Wsh%nloalg,&
&  Wsh%prtvol,Wsh%pawprtvol,Wsh%comm)

 call wfd_print(Weh,header="Shirley wavefunction descriptor")

 ABI_DEALLOCATE(eh_keep_ur)
 ABI_DEALLOCATE(eh_bks_mask)
 ABI_DEALLOCATE(eh_nband)
 !
 ! =====================================================================
 ! ==== Rotate the input wavefunctions to get the optimal basis set ====
 ! =====================================================================

 fft_fact = one/Wsh%nfft
 ABI_ALLOCATE(ur1,(Wsh%nfft*nspinor))
 ABI_ALLOCATE(ur2,(Wsh%nfft*nspinor))
 ABI_ALLOCATE(ur12,(Wsh%nfft*nspinor))
 !
 write(msg,'(a,f12.1,a)')' Memory needed for storing eh_ur= ',two*gwpc*Wsh%nfft*nspinor*eh_size*b2Mb,' [Mb]'
 call wrtout(std_out,msg,'PERS')

 ABI_ALLOCATE(eh_ur,(Wsh%nfft*nspinor,eh_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out-of-memory in eh_ur")
 eh_ur = czero

 ABI_ALLOCATE(eh_ug,(npw_gamma*nspinor,eh_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out-of-memory in eh_ug")
 eh_ug = czero

 do midx=1,ovlp_size ! Loop over the single particle orbitals.
   band1  = Oeh%idx2bk(1,midx) ! TODO to be removed.
   band2  = Oeh%idx2bk(2,midx)

   call wfd_get_ur(Wsh,band1,k1,s1,ur1)
   call wfd_get_ur(Wsh,band2,k1,s1,ur2)
   ur12 = CONJG(ur1) * ur2   ! Use same convention as the one used in ovlp2_init.
   !
   do eh=1,eh_size ! Construct the new optimal basis set.
     eh_ur(:,eh) = eh_ur(:,eh) + Oeh%mat(midx,base+eh) * ur12
   end do
   !
 end do
 !
 ! NC: Normalize the basis set.
 do eh=1,eh_size
   norm1 = xdotc(Weh%nfft*nspinor,eh_ur(:,eh),1,eh_ur(:,eh),1) * fft_fact 
   eh_ur(:,eh) = eh_ur(:,eh)/SQRT(norm1) 
   !write(std_out,*)" eh_ur integrates to: ", xdotc(Weh%nfft*nspinor,eh_ur(:,eh),1,eh_ur(:,eh),1) * fft_fact 
 end do
 !
 ! From the FFT mesh to the G-sphere.
 ABI_ALLOCATE(ehg,(npw_gamma*nspinor))
 ABI_ALLOCATE(gwpc_tmp,(Wsh%nfft*nspinor))
 !
 ! ============================================================================
 ! ==== Construct new optimal basis set in G-space and save results in Weh ====
 ! ============================================================================
 do eh=1,eh_size 
   gwpc_tmp = eh_ur(:,eh)

#if 1
   ABI_ALLOCATE(dpc_tmp,(Wsh%nfft*nspinor))
   dpc_tmp  = eh_ur(:,eh)
   call fourdp_c2c_ip(dpc_tmp,-1,Wsh%MPI_enreg,Wsh%nfft*nspinor,Wsh%ngfft,paral_kgb,0)
   !
   do ig=1,npw_gamma       
     fft_idx = igfft0(ig)
     if (fft_idx/=0) then ! G-G0 belong to the FFT mesh.
       if (fft_idx>Wsh%nfft .or.fft_idx<0) then
         MSG_ERROR("fft_idx bug")
       end if
       ehg(ig)=dpc_tmp(fft_idx) 
     else                 ! Set this component to zero.
       MSG_ERROR("fft_idx bug")
       ehg(ig)=czero
     end if
   end do
   ABI_DEALLOCATE(dpc_tmp)
#else
   ! FIXME does not work anymore.
   gbound_k => Weh%Kdata(k1)%gbound
   call fft_ur(nspinor,npw_gamma,istwf1,paral_kgb,Wsh%nfft,Wsh%mgfft,Wsh%ngfft,gwpc_tmp,&
&    gbound_k,kg_gamma,ehg,Wsh%MPI_enreg)
#endif
   ! NC: Normalize the basis set using the eigenvalues of the overlap matrix.
   !if (usepaw==0) then
   !  ehg = ehg/SQRT(Oeh%eigene(base+eh)) 
   !end if
   !%call wfd_push_ug(Weh,eh,k1,s1,Cryst,ehg,update_ur=.FALSE.,update_cprj=.FALSE.)
   eh_ug(:,eh) = ehg
 end do
 !
 ! ======================================
 ! ==== Orthonormalize the basis set ====
 ! ======================================
 ABI_ALLOCATE(cf_ovlp,(eh_size,eh_size))

 call blas_cholesky_ortho(npw_gamma,eh_size,eh_ug,cf_ovlp)

 ABI_DEALLOCATE(cf_ovlp)
 !
 ! Push data.
 do eh=1,eh_size
   call wfd_push_ug(Weh,eh,k1,s1,Cryst,eh_ug(:,eh),update_ur=.FALSE.,update_cprj=.FALSE.)
 end do

 ABI_DEALLOCATE(gwpc_tmp)
 ABI_DEALLOCATE(ehg)
 ABI_DEALLOCATE(eh_ur)
 ABI_DEALLOCATE(eh_ug)
 ABI_DEALLOCATE(kg_gamma)
 ABI_DEALLOCATE(igfft0)
 ABI_DEALLOCATE(ur1)
 ABI_DEALLOCATE(ur2)
 ABI_DEALLOCATE(ur12)

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0
 write(std_out,*)" E-H Rotation cpu_time, wall_time ",cpu_time, wall_time
 !
 ! ==================================
 ! ==== sh2eh = <S_i* S_j| EH_a> ====
 ! ==================================
 call timein(cpu0,wall0)
 ABI_ALLOCATE(sh2eh,(ovlp_size,eh_size))
 sh2eh=czero

 call wfd_change_ngfft(Weh,Cryst,Psps,Wsh%ngfft) ! Make sure the two set of wave are on the same mesh.
 ABI_ALLOCATE(ur1,(Weh%nfft*nspinor))

 do eh=1,eh_size
   call wfd_get_ur(Weh,eh,k1,s1,ur1)
   do ij=1,ovlp_size
     sh2eh(ij,eh) = DOT_PRODUCT(ur34_big(:,ij),ur1)
   end do
 end do
 ABI_DEALLOCATE(ur1)

 ABI_DEALLOCATE(ur34_big)

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0
 write(std_out,*)" E-H Projection cpu_time, wall_time ",cpu_time, wall_time
 !
 call wfd_test_ortho(Weh,Cryst,Pawtab,unit=std_out,mode_paral="COLL")
 !
 ! Deallocate memory.
 call ovlp_free(Oeh)

 DBG_EXIT("COLL")

end subroutine wfd_shirley_to_eh       
!!***

!----------------------------------------------------------------------

END MODULE m_shirley
!!***
