!{\src2tex{textfont=tt}}
!****m* ABINIT/m_oscillators
!! NAME
!!  m_oscillators
!!
!! FUNCTION
!!  This module contains procedures to calculate the oscillator matrix elements used in the GW code.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_oscillators

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_gwdefs,    only : czero_gw
 use m_crystal,   only : crystal_structure
 use m_bz_mesh,   only : bz_mesh_type, get_BZ_item, get_BZ_diff
 use m_fft_mesh,  only : cigfft
 use m_gsphere,   only : gvectors_type
 use m_wfs,       only : wfs_descriptor, wfd_get_ur

 implicit none

 private 

!Structures
!!***

!----------------------------------------------------------------------

!!****t* m_oscillators/oscillator_t
!! NAME
!! oscillator_t
!!
!! FUNCTION
!! The oscillator_t structured datatype stores the matrix element <k-q,b1,s2|e^{-i(q+G).r}|k,b2,s2> 
!! for a given (k,q,spin) set of indeces, with k being a point in the first Brillouin zone.
!! 
!!
!! SOURCE

 type,public :: oscillator_t

!scalars
  integer :: isppol        ! The spin of the wavefunctions (1 if noncollinear)
  integer :: nspinor       ! The number of spinorial components.
  integer :: ng            ! The number of G in the matrix elements
  integer :: map2sphere    ! 1 if the matrix elements are given on the G-sphere, 0 if on the FFT box

  integer :: ik_bz         ! The index of k in the BZ
  integer :: ik_ibz
  integer :: isym_k
  integer :: itim_k
  real(dp) :: k_bz(3)
  complex(dpc) :: ph_mkt

  integer :: iq_bz         ! The index of q in the BZ
  !integer :: iq_ibz
  !integer :: isym_q
  !integer :: itim_q
  real(dp) :: q_bz(3)
  !complex(gwpc) :: ph_mqt

  integer :: ikmq_bz       ! The index of the (k-q) point folded in the BZ (Umklapp might be needed).
  integer :: ikmq_ibz
  integer :: isym_kmq
  integer :: itim_kmq
  real(dp) :: kmq_bz(3)
  complex(dpc) :: ph_mkmqt

!arrays
  integer :: b1_idx(2)     ! Initial and final band indeces for the "b1" entry.
  integer :: b2_idx(2)     ! Initial and final band indeces for the "b2" entry.
  integer :: g0(3)         ! Reciprocal coordinates of the G0 vector such that (k-q)=kp+G0 with kp is in the BZ.

  complex(gwpc),pointer :: mg(:,:,:)
  ! mg(ng*nspinor**2, b1_idx(1):b1_idx(2), b2_idx(1):b2_idx(2))
  !  The matrix element <k-q,b1,s1|e^{-i(q+G).r}|k,b2,s2>. 
  !  If nspinor==1 ==> delta_{s1,s2} and the spin index is given by ispin
  !  If nspinor==2 ==> The four possible combination of spinor components are packed in the first dimension.

 end type oscillator_t
!!***

 public :: rho_tw_g 
 public :: calc_wfwfg         ! Calculate the Fourier transform of the product u_{bk}^*(r).u_{b"k}(r) at an arbitrary k in the BZ.
 public :: sym_rhotwgq0       ! Symmetrize the oscillator matrix elements in the BZ in the special case of q=0.
 public :: nullify_oscillator
 public :: init_oscillator
 public :: destroy_oscillator
 public :: calc_pw_oscillator

!!***

!----------------------------------------------------------------------

!!****t* m_oscillators/ww_store_t
!! NAME
!! ww_store_t
!!
!! FUNCTION
!! Simple database used to store the matrix elements <k,b1,s|e^{-iG.r}|k,b2,s>  
!!
!! SOURCE

 type ww_store_t

  integer :: msize
  ! The maximum size of the buffer.

  integer :: nstored
  ! The number of matrix elements stored.

  integer :: nspinor
  integer :: nfft

  integer,pointer :: bbs_indeces(:,:)  SET2NULL        
  ! bbs_indeces(3,msize)
  ! The fist nstored entries give the indeces (b1,b2,s) of the matrix elements that
  ! have been stored.

  complex(gwpc),pointer :: buffer(:,:) SET2NULL
  ! buffer(nspinor*nfft,msize)
  ! The calculated matrix elements on the FFT mesh.

 end type ww_store_t
!!***

CONTAINS  !=========================================================================================================================
!!***

!!****f* m_oscillators/rho_tw_g
!! NAME
!! rho_tw_g
!!
!! FUNCTION
!! Calculate rhotwg(G)=<wfn1|exp(-i(q+G).r)|wfn2>
!!
!! INPUTS
!! dim_rtwg=Define the size of the output array rhotwg
!!   === for nspinor==1 ===
!!    dim_rtwg=1
!!   === for nspinor==2 ===
!!    dim_rtwg=2 if only <up|up>, <dwn|dwn> matrix elements are required
!!    dim_rtwg=4 for <up|up>, <dwn|dwn>, <up|dwn> and <dwn|up>.
!! map2sphere= 1 to retrieve Fourier components indexed according to igfftg0.
!!             0 to retrieve Fourier components indexed according to the FFT box.
!!               NOTE: If map2sphere==0 npwvec must be equal to nr
!! use_padfft= Only compatible with map2sphere 1. 
!!             1 if matrix elements are calculated via zero-padded FFT.
!!             0 R-->G Transform in done on the full FFT box.
!! igfftg0(npwvec*map2sphere)=index of G-G_o in the FFT array for each G in the sphere.
!! i1=1 if kbz1 = Sk1, 2 if kbz1 = -Sk_1 (k_1 is in the IBZ)
!! i2=1 if kbz2 = Sk2, 2 if kbz2 = -Sk_2 (k_2 is in the IBZ)
!! ktabr1(nr),ktabr2(nr)= tables R^-1(r-t) for the two k-points
!! ktabp1,ktabp2 = phase factors for non-simmorphic symmetries e^{-i 2\pi kbz.\tau} 
!! MPI_enreg=Information about MPI parallelization
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npwvec=number of plane waves (in the sphere if map2sphere==1, in the FFT box if map2sphere==1)
!! nr=number of FFT grid points
!! nspinor=number of spinorial components.
!! spinrot1(4),spinrot2(4)=components of the spinor rotation matrix :
!!  spinrot(1)=$\cos\phi/2$
!!  spinrot(2)=$\sin\phi/2 \times u_x$
!!  spinrot(3)=$\sin\phi/2 \times u_y$
!!  spinrot(4)=$\sin\phi/2 \times u_z$
!!   where $\phi$ is the angle of rotation, and
!!   $(u_x,u_y,u_z)$ is the normalized direction of the rotation axis
!! tim_fourdp=Option definig the part of the code that has to be analysed. 
!!            0 if not assigned.
!!            1 if called from inside screening 
!!            2 if called from inside sigma
!! wfn1(nr),wfn2(nr)=the two wavefunctions (periodic part)
!! [nhat12(2,nr,nspinor**2)]=Compensation charge in real space to be added to \Psi_1^*\Psi_2 -- Only for PAW.
!!
!! OUTPUT
!! rhotwg(npwsigx)=density of a pair of states, in reciprocal space
!!
!! PARENTS
!!      calc_sig_ppm_eet,calc_sigc_me,calc_sigx_me,cchi0,cchi0q0
!!      cchi0q0_intraband,check_completeness,cohsex_me,exc_build_block
!!      exc_build_ham,m_fft_prof,m_oscillators,rdm
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine rho_tw_g(paral_kgb,nspinor,npwvec,nr,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
& wfn1,i1,ktabr1,ktabp1,spinrot1,&
& wfn2,i2,ktabr2,ktabp2,spinrot2,&
& dim_rtwg,rhotwg,tim_fourdp,MPI_enreg)
!& dim_rtwg,rhotwg,tim_fourdp,MPI_enreg,nhat12)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rho_tw_g'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: paral_kgb,i1,i2,npwvec,nr,tim_fourdp,nspinor,dim_rtwg,map2sphere,use_padfft
 complex(dpc),intent(in) :: ktabp1,ktabp2
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 !integer,intent(in) :: gbound(2*mgfft+8,2)
 integer,intent(in) :: gbound(:,:)
 integer,intent(in) :: igfftg0(npwvec*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 real(dp),intent(in) :: spinrot1(4),spinrot2(4)
 complex(gwpc),intent(in) :: wfn1(nr*nspinor),wfn2(nr*nspinor)
 complex(gwpc),intent(out) :: rhotwg(npwvec*dim_rtwg)
! real(dp),optional,intent(in) :: nhat12(2,nr,nspinor**2)

!Local variables-------------------------------
!scalars
 integer :: ig,ir,ir1,ir2,igfft,iab,spad1,spad2,spad0,ispinor
 integer :: nx,ny,nz,ldx,ldy,ldz,mgfft
 complex(gwpc) :: u1a,u1b,u2a,u2b
!arrays
 integer :: spinor_pad(2,4) 
 complex(dpc) :: spinrot_mat1(2,2),spinrot_mat2(2,2)
 complex(dpc),allocatable :: cwavef1(:),cwavef2(:)
 complex(dpc),allocatable :: usk(:),uu(:)

! *************************************************************************

 SELECT CASE (nspinor)

 CASE (1) ! Collinear case.
   !
   ! Form rho-twiddle(r)=u_1^*(r,b1,kbz1) u_2(r,b2,kbz2), to account for symmetries:
   ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz) 
   !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
   !
   ABI_ALLOCATE(uu,(nr))
   ABI_ALLOCATE(usk,(nr))

   uu  = wfn1(ktabr1)*ktabp1; if (i1==1) uu  = CONJG(uu)
   usk = wfn2(ktabr2)*ktabp2; if (i2==2) usk = CONJG(usk)
   uu  = uu * usk

   ! Add compensation charge.
   !if (PRESENT(nhat12)) then 
   !  uu = uu + CMPLX(nhat12(1,:,1),nhat12(2,:,1))
   !end if

   SELECT CASE (map2sphere)

   CASE (0) ! Need results on the full FFT box thus cannot use zero-padded FFT.

     call fourdp_c2c_ip(uu,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)
     rhotwg=uu
     !call fourdp_c2c_op(uu,rhotwg,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)

   CASE (1) ! Need results on the G-sphere. Call zero-padded FFT routines if required.

     if (use_padfft==1) then
       nx =ngfft(1); ny =ngfft(2); nz =ngfft(3); mgfft = MAXVAL(ngfft(1:3))
       ldx=nx      ; ldy=ny      ; ldz=nz
       call padded_fourwf_cplx(uu,ngfft,nx,ny,nz,ldx,ldy,ldz,mgfft,-1,gbound)
     else
       call fourdp_c2c_ip(uu,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)
     end if
     !
     ! From the FFT to the G-sphere.
     do ig=1,npwvec
       igfft=igfftg0(ig)
       if (igfft/=0) then ! G-G0 belong to the FFT mesh.
         rhotwg(ig)=uu(igfft) 
       else               ! Set this component to zero.
         rhotwg(ig)=czero_gw 
       end if
     end do

   CASE DEFAULT
     MSG_BUG("Wrong map2sphere")
   END SELECT

   ABI_DEALLOCATE(uu)
   ABI_DEALLOCATE(usk)

   RETURN

 CASE (2) ! Spinorial case.

   MSG_ERROR("Add zero-padded FFT") !TODO

   ABI_ALLOCATE(cwavef1,(nr*nspinor))
   ABI_ALLOCATE(cwavef2,(nr*nspinor))

   ! === Apply Time-reversal if required ===
   ! \psi_{-k}^1 =  (\psi_k^2)^*
   ! \psi_{-k}^2 = -(\psi_k^1)^*
   if (i1==1) then 
     cwavef1(:)=wfn1(:) 
   else if (i1==2) then
     cwavef1(1:nr)     = CONJG(wfn1(nr+1:2*nr))
     cwavef1(nr+1:2*nr)=-CONJG(wfn1(1:nr))
   else 
     MSG_ERROR('Wrong i1 in spinor')
   end if

   if (i2==1) then 
     cwavef2(:)=wfn2(:) 
   else if (i2==2) then
     cwavef2(1:nr)     = CONJG(wfn2(nr+1:2*nr))
     cwavef2(nr+1:2*nr)=-CONJG(wfn2(1:nr))
   else 
     MSG_ERROR('Wrong i2 in spinor')
   end if

   ! === Rotate wavefunctions in r-space ===
   do ispinor=1,nspinor
     spad0=(ispinor-1)*nr
     do ir=1,nr
       ir1=ktabr1(ir) ; ir2=ktabr2(ir)
       cwavef1(ir+spad0) = cwavef1(ir1+spad0)*ktabp1
       cwavef2(ir+spad0) = cwavef2(ir2+spad0)*ktabp2
     end do 
   end do !ispinor

   ! === Rotation in spinor space ===
   !spinrots1=spinrot1(1) ; spinrots2=spinrot2(1)
   !spinrotx1=spinrot1(2) ; spinrotx2=spinrot2(2)
   !spinroty1=spinrot1(3) ; spinroty2=spinrot2(3)
   !spinrotz1=spinrot1(4) ; spinrotz2=spinrot2(4)
   spinrot_mat1(1,1)= spinrot1(1) + j_dpc*spinrot1(4)
   spinrot_mat1(1,2)= spinrot1(3) + j_dpc*spinrot1(2)
   spinrot_mat1(2,1)=-spinrot1(3) + j_dpc*spinrot1(2)
   spinrot_mat1(2,2)= spinrot1(1) - j_dpc*spinrot1(4)

   spinrot_mat2(1,1)= spinrot2(1) + j_dpc*spinrot2(4)
   spinrot_mat2(1,2)= spinrot2(3) + j_dpc*spinrot2(2)
   spinrot_mat2(2,1)=-spinrot2(3) + j_dpc*spinrot2(2)
   spinrot_mat2(2,2)= spinrot2(1) - j_dpc*spinrot2(4)

   do ir=1,nr
     !ar=wavefspinor(1,ir)
     !ai=wavefspinor(2,ir)
     !br=wavefspinor(1,npw2+ir)
     !bi=wavefspinor(2,npw2+ir)
     u1a=cwavef1(ir) 
     u1b=cwavef1(ir+nr) 
     cwavef1(ir)   =spinrot_mat1(1,1)*u1a+spinrot_mat1(1,2)*u1b
     cwavef1(ir+nr)=spinrot_mat1(2,1)*u1a+spinrot_mat1(2,2)*u1b
     u2a=cwavef2(ir) 
     u2b=cwavef2(ir+nr) 
     cwavef2(ir)   =spinrot_mat2(1,1)*u2a+spinrot_mat2(1,2)*u2b
     cwavef2(ir+nr)=spinrot_mat2(2,1)*u2a+spinrot_mat2(2,2)*u2b
     !wavefspinor(1,ir)     = spinrots*ar-spinrotz*ai +spinroty*br-spinrotx*bi
     !wavefspinor(2,ir)     = spinrots*ai+spinrotz*ar +spinroty*bi+spinrotx*br
     !wavefspinor(1,npw2+ir)=-spinroty*ar-spinrotx*ai +spinrots*br+spinrotz*bi
     !wavefspinor(2,npw2+ir)=-spinroty*ai+spinrotx*ar +spinrots*bi-spinrotz*br
   end do

   spinor_pad(:,:)=RESHAPE((/0,0,nr,nr,0,nr,nr,0/),(/2,4/))
   ABI_ALLOCATE(uu,(nr))

   do iab=1,dim_rtwg
     spad1=spinor_pad(1,iab)
     spad2=spinor_pad(2,iab)

     uu = CONJG(cwavef1(spad1+1:spad1+nr)) * cwavef2(spad2+1:spad2+nr)

     ! Add compensation charge.
     !if (PRESENT(nhat12)) then 
     !  uu = uu + CMPLX(nhat12(1,:,iab),nhat12(2,:,iab))
     !end if
     spad0=(iab-1)*npwvec

     SELECT CASE (map2sphere)
     
     CASE (0) ! Need results on the full FFT box thus cannot use zero-padded FFT.
       call fourdp_c2c_ip(uu,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)
       rhotwg(spad0+1:spad0+npwvec)=uu 

     CASE (1) ! Need results on the G-sphere. Call zero-padded FFT routines if required.

       if (use_padfft==1) then
         nx =ngfft(1); ny =ngfft(2); nz =ngfft(3); mgfft = MAXVAL(ngfft(1:3))
         ldx=nx      ; ldy=ny      ; ldz=nz
         call padded_fourwf_cplx(uu,ngfft,nx,ny,nz,ldx,ldy,ldz,mgfft,-1,gbound)
       else 
         call fourdp_c2c_ip(uu,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)
       end if

       do ig=1,npwvec      ! Have to map FFT to G-sphere.
         igfft=igfftg0(ig)
         if (igfft/=0) then ! G-G0 belong to the FFT mesh.
           rhotwg(ig+spad0)=uu(igfft) 
         else               ! Set this component to zero
           rhotwg(ig+spad0)=czero_gw
         end if
       end do

     CASE DEFAULT
       MSG_BUG("Wrong map2sphere")
     END SELECT
   end do !iab

   ABI_DEALLOCATE(cwavef1)
   ABI_DEALLOCATE(cwavef2)
   ABI_DEALLOCATE(uu)

   RETURN

 CASE DEFAULT 
   MSG_BUG('Wrong nspinor')
 END SELECT

end subroutine rho_tw_g
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/calc_wfwfg
!! NAME
!! calc_wfwfg
!!
!! FUNCTION
!!  Calculate the Fourier transform of the product u_{bk}^*(r).u_{b"k}(r) at an arbitrary k in the BZ.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet,calc_sigc_me,cchi0,cchi0q0,check_completeness
!!      cohsex_me
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine calc_wfwfg(MPI_enreg,paral_kgb,tim_fourdp,ktabr_k,ktabi_k,nfftot,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_wfwfg'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ktabi_k,nfftot,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ktabr_k(nfftot),ngfft_gw(18)
 complex(gwpc),intent(in) :: wfr_jb(nfftot),wfr_kb(nfftot)
 complex(gwpc),intent(out) :: wfg2_jk(nfftot)

!Local variables-------------------------------
!arrays
 complex(dpc),allocatable :: wfr2_dpcplx(:)
#if ! defined HAVE_GW_DPC
 complex(dpc),allocatable :: wfg2_dpcplx(:)
#endif

! *************************************************************************

 ! There is no need to take into account phases arising from non-symmorphic
 ! operations since the wavefunctions are evaluated at the same k-point.
 ABI_ALLOCATE(wfr2_dpcplx,(nfftot))

 SELECT CASE (ktabi_k)

 CASE (1)
   wfr2_dpcplx = CONJG(wfr_jb(ktabr_k)) * wfr_kb(ktabr_k)

 CASE (2) ! Conjugate the product if time-reversal is used to reconstruct this k-point
   wfr2_dpcplx = wfr_jb(ktabr_k) * CONJG(wfr_kb(ktabr_k))

 CASE DEFAULT
   MSG_ERROR("Wrong ktabi_k")
 END SELECT

 ! Transform to Fourier space (result in wfg2_jk)
#if defined HAVE_GW_DPC
 call fourdp_c2c_op(wfr2_dpcplx,wfg2_jk,-1,MPI_enreg,nfftot,ngfft_gw,paral_kgb,tim_fourdp)
#else
 ABI_ALLOCATE(wfg2_dpcplx,(nfftot))
 call fourdp_c2c_op(wfr2_dpcplx,wfg2_dpcplx,-1,MPI_enreg,nfftot,ngfft_gw,paral_kgb,tim_fourdp)
 wfg2_jk=wfg2_dpcplx
 ABI_DEALLOCATE(wfg2_dpcplx)
#endif

 ABI_DEALLOCATE(wfr2_dpcplx)

end subroutine calc_wfwfg
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/sym_rhotwgq0   
!! NAME
!!  sym_rhotwgq0
!!
!! FUNCTION
!!  Symmetrization of the oscillator matrix elements <k-q,b1|exp(-i(q+G).r)|k,b2> in the special case of q=0.
!!  The matrix elements in the full BZ is obtained from the matrix elements in the IBZ by
!!  rotating the wavefunctions and taking into account time reversal symmetry.
!!  strictly speaking the symmetrization can be performed only for non-degenerate states.
!!
!! INPUTS
!!  Gsph<gvectors_type>=Info on the G-sphere used to describe wavefunctions and W (the largest one is actually stored).  
!!  npw=Number of G-vectors
!!  dim_rtwg=Number of spin-spin combinations, 1 for collinear spin, 4 is nspinor==2 (TODO NOT CODED)
!!  itim_k=2 if time reversal is used to reconstruct the k in the BZ, 1 otherwise.
!!  isym_k=The index of the symmetry symrec rotains k_IBZ onto k_BZ.
!!  rhxtwg_in(dim_rtwg*npw)=The input matrix elements in the IBZ.
!!
!! OUTPUT
!!  rhxtwg_sym(dim_rtwg*npw)=The symmetrized matrix elements in the BZ.
!!
!! NOTES
!! Let M_{G}(k,q) =<k-q,b1|exp(-i(q+G).r)|k,b2>
!!  At q ==0, supposing non-degenerate bands, one obtains:
!!
!!  1) M_{ SG}( Sk) = e^{-iSG.t} M_{G}   (k)
!!  2) M_{-SG}(-Sk) = e^{+iSG.t} M_{G}^* (k)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npw,rhxtwg_in,Gsph) result(rhxtwg_sym)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sym_rhotwgq0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,dim_rtwg,itim_k,isym_k
 type(Gvectors_type),intent(in) :: Gsph
!arrays
 complex(gwpc),intent(in) :: rhxtwg_in(dim_rtwg*npw) 
 complex(gwpc) :: rhxtwg_sym(dim_rtwg*npw) 

!Local variables ------------------------------
!scalars
 integer :: ig
 character(len=500) :: msg

!************************************************************************

 if (dim_rtwg/=1) &
&  MSG_ERROR("dim_rtwg/=1 not coded")

 SELECT CASE (isym_k)

 CASE (1) ! Fractional translation associated to E is assumed to be (zero,zero,zero).

   SELECT CASE (itim_k) 
   CASE (1) ! Identity, no time-reversal. No symmetrization is needed.
     rhxtwg_sym(:) = rhxtwg_in(:)

   CASE (2) ! Identity + Time-reversal.
     do ig=1,npw
       rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = CONJG(rhxtwg_in(ig))
     end do

   CASE DEFAULT
     write(msg,'(a,i0)')"Wrong value of itim_k: ",itim_k
     MSG_ERROR(msg)
   END SELECT

 CASE DEFAULT ! Rotate wavefunctions.

   SELECT CASE (itim_k) 

   CASE (1) ! no time-reversal, only rotation.
    do ig=1,npw
      rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = rhxtwg_in(ig) * Gsph%phmSGt(ig,isym_k) 
    end do

   CASE (2) ! time-reversal + spatial rotation.
    do ig=1,npw
      rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = CONJG( rhxtwg_in(ig) * Gsph%phmSGt(ig,isym_k) )
    end do

   CASE DEFAULT
     write(msg,'(a,i0)')"Wrong value of itim_k: ",itim_k
     MSG_ERROR(msg)
   END SELECT

 END SELECT 

end function sym_rhotwgq0
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/nullify_oscillator
!! NAME
!!  nullify_oscillator
!!
!! FUNCTION
!!  Initialize the pointers in the structures datatype oscillator_t to NULL
!!
!! INPUTS
!!  Osc<oscillator_t>
!!
!! SIDE EFFECTS
!!  Pointers set to NULL.
!!
!! PARENTS
!!      m_oscillators
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine nullify_oscillator(Osc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_oscillator'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(oscillator_t),intent(inout) :: Osc

! *************************************************************************

 !@oscillator_t
 nullify(Osc%mg)

end subroutine nullify_oscillator
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/init_oscillator
!! NAME
!!  init_oscillator
!!
!! FUNCTION
!!  Create an instance of the oscillator_t structured datatype.
!!
!! INPUTS
!!  isppol=The index of the spin.
!!  ik_ibz=The index of the k-point in the BZ.
!!  iq_ibz=The index of the q-point in the BZ.
!!  Kmesh,Qmesh<BZ_mesh_type>=Structures defining the BZ sampling.
!!  ng=Number of G-vectors in the matrix elements (<= Number of FFT points in the wfs_descriptor)
!!  nspinor=Number of spinorial components.
!!  b1_idx(2)=First and last value of the index associated to the b1 band.
!!  b2_idx(2)=First and last value of the index associated to the b2 band.
!!  [map2sphere]=1 if the G index in the oscillator runs over the G-sphere (DEFAULT)
!!               0 if the oscillator is given on the FFT box in reciprocal space
!!
!! SIDE EFFECTS
!!  Osc<oscillator_t>=Initialized in output with correct dimensions and indeces. Osc%mg is allocated here.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine init_oscillator(Osc,isppol,ik_bz,Kmesh,iq_bz,Qmesh,ng,nspinor,b1_idx,b2_idx,map2sphere)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_oscillator'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,iq_bz,ng,nspinor,isppol
 integer,optional,intent(in) :: map2sphere
 type(oscillator_t),intent(out) :: Osc
 type(BZ_mesh_type),intent(in) :: Kmesh,Qmesh
!arrays
 integer,intent(in) :: b1_idx(2),b2_idx(2)

!Local variables ------------------------------
!scalars
 integer :: istat,nfound
!arrays

! *************************************************************************

 !@oscillator_t
 call nullify_oscillator(Osc)

 Osc%isppol  = isppol
 Osc%nspinor = nspinor
 Osc%ng      = ng
 Osc%map2sphere=1; if (PRESENT(map2sphere)) Osc%map2sphere=map2sphere

 ABI_CHECK(Osc%map2sphere==0.or.Osc%map2sphere==1,"wrong map2sphere")

 ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
 Osc%ik_bz   = ik_bz
 call get_BZ_item(Kmesh,Osc%ik_bz,Osc%k_bz,Osc%ik_ibz,Osc%isym_k,Osc%itim_k,Osc%ph_mkt)
                                                                                              
 ! * Get index of k-q in the BZ, stop if not found as the weight=one/nkbz is not correct.
 Osc%iq_bz = iq_bz
 Osc%q_bz  = Qmesh%bz(:,iq_bz)

 !% call get_BZ_item(Qmesh,Osc%iq_bz,Osc%q_bz,Osc%iq_ibz,Osc%isym_q,Osc%itim_q,Osc%ph_mqt)

 ! Get index of k-q = BZ(k-q) + g0. Note that k-q might fall outside the first BZ.
 call get_BZ_diff(Kmesh,Osc%k_bz,Osc%q_bz,Osc%ikmq_bz,Osc%g0,nfound) 
 ABI_CHECK(nfound==1,"nfound/=1")
                                                                                              
 ! * Get ikmq_ibz, non-symmorphic phase, ph_mkmqt, and symmetries from ikmq_bz.
 call get_BZ_item(Kmesh,Osc%ikmq_bz,Osc%kmq_bz,Osc%ikmq_ibz,Osc%isym_kmq,Osc%itim_kmq,Osc%ph_mkmqt)

 ! Allocate the buffer to store the matrix elements.
 Osc%b1_idx = b1_idx
 Osc%b2_idx = b2_idx

 ABI_ALLOCATE( Osc%mg,(ng*nspinor**2, Osc%b1_idx(1):Osc%b1_idx(2), Osc%b2_idx(1):Osc%b2_idx(2)))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out-of-memory in Osc%mg")

end subroutine init_oscillator
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/destroy_oscillator
!! NAME
!!  destroy_oscillator
!!
!! FUNCTION
!!  Free the memory allocated in an instance of the oscillator_t structured datatype.
!!
!! SIDE EFFECTS
!!  Osc<oscillator_t>=All allocatated memory is released.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine destroy_oscillator(Osc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_oscillator'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(oscillator_t),intent(inout) :: Osc

! *************************************************************************

 !@oscillator_t
 if (associated(Osc%mg))  then
   ABI_DEALLOCATE(Osc%mg)
 end if

end subroutine destroy_oscillator
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/calc_pw_oscillator
!! NAME
!!  calc_pw_oscillator
!!
!! FUNCTION
!!  Calculate the planewave contribution to a set of oscillators at fixed (k,q,spin) for a set of (b1,b2) bands.
!!
!! INPUTS
!!  ik_ibz=The index of the k-point in the BZ.
!!  iq_ibz=The index of the q-point in the BZ.
!!  ng=Number of G-vectors in the matrix elements.
!!  nspinor=Number of spinorial components.
!!  b1_idx(2)=First and last value of the index associated to the b1 band.
!!  b2_idx(2)=First and last value of the index associated to the b2 band.
!!  irrottb(Wfd%nfftot,Cryst%nsym) The FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.
!!
!! SIDE EFFECTS
!!  Osc<oscillator_t>)
!!    %mg(ng*nspinor,b1min:b1max,b2min:b2max)=Calculated from the wavefunctions in the irreducible wedge.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine calc_pw_oscillator(Wfd,Cryst,Osc,irottb,mg0)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_pw_oscillator'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_structure),intent(in) :: Cryst
 type(wfs_descriptor),intent(inout) :: Wfd
 type(oscillator_t),intent(out) :: Osc
!arrays
 integer,intent(in) :: mg0(3)
 integer,target,intent(in) :: irottb(Wfd%nfftot,Cryst%nsym)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=0,paral_kgb=0
 integer :: use_padfft,istat,ierr
 integer :: ig01,ig02,ig03,ik_ibz,ikmq_ibz,isym_k,isym_kmq,isppol,dim_rtwg,ib1,ib2
 character(len=500) :: msg
!arrays
 integer,pointer :: tabr_k(:),tabr_kmq(:)
 integer,allocatable :: gbound(:,:),igmg0_fft(:,:,:,:)
 real(dp) :: spinrot_k(4),spinrot_kmq(4)
 complex(gwpc),allocatable :: wfr1(:),wfr2(:)

! *************************************************************************

 use_padfft=0 ! TODO no padded-FFT for the time being 
 ABI_ALLOCATE(gbound,(1,0))

 if (Osc%map2sphere==0) then
   msg = " map2sphere==0, but Osc%ng is not equal to the number of FFT points"
   ABI_CHECK(Osc%ng == Wfd%nfftot,msg)
 end if

 ik_ibz = Osc%ik_ibz
 isym_k = Osc%isym_k
 tabr_k => irottb(:,isym_k)     ! Table for rotated FFT points
 spinrot_k =  Cryst%spinrot(:,isym_k)
 
 ikmq_ibz = Osc%ikmq_ibz
 isym_kmq = Osc%isym_kmq
 tabr_kmq => irottb(:,isym_kmq) ! Table for rotated FFT points
 spinrot_kmq = Cryst%spinrot(:,isym_kmq)

 ! * Get the G-G0 shift for the FFT of the oscillators.
 ig01 = Osc%g0(1)+mg0(1)+1  
 ig02 = Osc%g0(2)+mg0(2)+1
 ig03 = Osc%g0(3)+mg0(3)+1

 ! TODO this part has to be optimzed 
!  For each G, and each G0 vector, this table gives the FFT grid index of G-G0.
!  Note the if G-G0 falls outside the FFT box, the corresponding index will be set to 0.

 ABI_ALLOCATE(igmg0_fft,(Wfd%npwwfn,2*mg0(1)+1,2*mg0(2)+1,2*mg0(3)+1))

 call cigfft(mg0,Wfd%npwwfn,Wfd%ngfft,Wfd%gvec,igmg0_fft,ierr)

 ABI_ALLOCATE(wfr1,(Wfd%nfftot*Wfd%nspinor))
 ABI_ALLOCATE(wfr2,(Wfd%nfftot*Wfd%nspinor))

 isppol=Osc%isppol; dim_rtwg=1; if (Osc%nspinor==2) dim_rtwg=4

 do ib1=Osc%b1_idx(1),Osc%b1_idx(2)
   call wfd_get_ur(Wfd,ib1,ikmq_ibz,isppol,wfr1)

   do ib2=Osc%b2_idx(1),Osc%b2_idx(2)
     call wfd_get_ur(Wfd,ib2,  ik_ibz,isppol,wfr2)

     call rho_tw_g(paral_kgb,Wfd%nspinor,Osc%ng,Wfd%nfftot,Wfd%ngfft,Osc%map2sphere,use_padfft,igmg0_fft(:,ig01,ig02,ig03),gbound,&
&      wfr1, Osc%itim_kmq, tabr_kmq, Osc%ph_mkmqt, spinrot_kmq,&
&      wfr2, Osc%itim_k  , tabr_k  , Osc%ph_mkt  , spinrot_k,  &
&      dim_rtwg,Osc%mg(:,ib1,ib2),tim_fourdp,Wfd%MPI_enreg) 
   end do
 end do

 ABI_DEALLOCATE(wfr1)
 ABI_DEALLOCATE(wfr2)
 ABI_DEALLOCATE(gbound)
 istat = ABI_ALLOC_STAT

 ABI_DEALLOCATE(igmg0_fft)

end subroutine calc_pw_oscillator
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/wwg_store_init
!! NAME
!!  wwg_store_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine wwg_store_init(ww_store,nspinor,nfft,msize)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wwg_store_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msize,nspinor,nfft
 type(ww_store_t),intent(inout) ::  ww_store
!arrays

!Local variables ------------------------------
!scalars
 integer :: ierr
!arrays

! *************************************************************************

 ww_store%msize   = msize
 ww_store%nstored = 0 

 ww_store%nspinor = nspinor
 ww_store%nfft    = nfft

 ABI_ALLOCATE(ww_store%bbs_indeces,(3,msize))
 ww_store%bbs_indeces=-1
 ABI_ALLOCATE(ww_store%buffer,(nspinor*nfft,msize))
 ierr = ABI_ALLOC_STAT

end subroutine wwg_store_init
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/wwg_store_free
!! NAME
!!  wwg_store_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine wwg_store_free(ww_store)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wwg_store_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ww_store_t),intent(inout) ::  ww_store
!arrays

!Local variables ------------------------------

! *************************************************************************

 ww_store%msize   = 0
 ww_store%nstored = 0 
 ww_store%nspinor = 0 
 ww_store%nfft    = 0 

 if (associated(ww_store%bbs_indeces))  then
   ABI_DEALLOCATE(ww_store%bbs_indeces)
 end if
 if (associated(ww_store%buffer     ))  then
   ABI_DEALLOCATE(ww_store%buffer)
 end if

end subroutine wwg_store_free
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/wwg_store_push
!! NAME
!!  wwg_store_push
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine wwg_store_push(ww_store,bbs,ww)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wwg_store_push'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ww_store_t),intent(inout) ::  ww_store

!arrays
 integer,intent(in) :: bbs(3)
 complex(gwpc),intent(in) :: ww(ww_store%nspinor*ww_store%nfft)

!Local variables ------------------------------
!scalars
 integer :: next
!arrays

! *************************************************************************

 next = ww_store%nstored +1
 if (next <= ww_store%msize) then
   ww_store%bbs_indeces(:,next) = bbs
   ww_store%buffer(:,next) = ww
   ww_store%nstored = next
 end if

end subroutine wwg_store_push
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/wwg_store_pull
!! NAME
!!  wwg_store_pull
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine wwg_store_pull(ww_store,bbs,ww,gotit)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wwg_store_pull'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: gotit
 type(ww_store_t),intent(in) ::  ww_store

!arrays
 integer,intent(in) :: bbs(3)
 complex(gwpc),intent(out) :: ww(ww_store%nspinor*ww_store%nfft)

!Local variables ------------------------------
!scalars
 integer :: ist
!arrays
 !integer :: exc_bbs(3)

! *************************************************************************

 gotit=.FALSE.

 !exc_bbs(2) = bbs(1)
 !exc_bbs(1) = bbs(2)
 !exc_bbs(3) = bbs(3)

 do ist =1,ww_store%nstored
   if ( ALL(ww_store%bbs_indeces(:,ist) == bbs) ) then
     ww = ww_store%buffer(:,ist)
     gotit=.TRUE.; EXIT
   end if
 end do

end subroutine wwg_store_pull
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/ostrength
!! NAME
!! ostrength
!!
!! FUNCTION
!! Calculate rhotwg(G)=<ur1|exp(-iG.r)|ur2> 
!! Note that no symmetrization is done here. This version will be used 
!! when both in cchi0 and sigma, the periodic part will be symmetrized 
!! by using wfd_sym_ur.
!!
!! INPUTS
!! dim_rtwg=Define the size of the output array rhotwg
!!   === for nspinor==1 ===
!!    dim_rtwg=1
!!   === for nspinor==2 ===
!!    dim_rtwg=2 if only <up|up>, <dwn|dwn> matrix elements are required
!!    dim_rtwg=4 for <up|up>, <dwn|dwn>, <up|dwn> and <dwn|up>.
!! use_padfft= 1 if matrix elements are calculated via zero-padded FFT.
!!             0 R-->G Transform in done on the full FFT box.
!! igfft(npw)=index of G in the FFT array for each G in the sphere.
!! MPI_enreg=Information about MPI parallelization
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npw=number of plane waves 
!! nfft=number of FFT grid points
!! nspinor=number of spinorial components.
!! tim_fourdp=Option definig the part of the code that has to be analysed. 
!!            0 if not assigned.
!!            1 if called from inside screening 
!!            2 if called from inside sigma
!! ur1(nfft*nspinor),ur2(nfft*nspinor)=the periodic part of the two wavefunctions. 
!!
!! OUTPUT
!! rhotwg(npwsigx)=density of a pair of states, in reciprocal space
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp_c2c_ip,padded_fourwf_cplx
!!
!! SOURCE

subroutine ostrength(paral_kgb,nspinor,npw,nfft,mgfft,ngfft,igfft,gbound,&
& ur1,ur2,dim_rtwg,rhotwg,tim_fourdp,MPI_enreg,use_padfft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ostrength'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: use_padfft
 integer,intent(in) :: paral_kgb,npw,nfft,mgfft,tim_fourdp
 integer,intent(in):: nspinor,dim_rtwg 
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),igfft(npw),ngfft(18)
 complex(gwpc),intent(in) :: ur1(nfft*nspinor),ur2(nfft*nspinor)
 complex(gwpc),intent(out) :: rhotwg(npw*dim_rtwg)

!Local variables-------------------------------
!scalars
 integer :: ig,fft_idx,nx,ny,nz,ldx,ldy,ldz
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 complex(dpc),allocatable :: uprod(:)

! *************************************************************************

 SELECT CASE (nspinor)

 CASE (1) ! Collinear case.
   !
   ABI_ALLOCATE(uprod,(nfft))
   uprod = CONJG(ur1) * ur2

   if (use_padfft==1) then
     nx =ngfft(1); ny =ngfft(2); nz =ngfft(3)
     ldx=nx      ; ldy=ny      ; ldz=nz
     call padded_fourwf_cplx(uprod,ngfft,nx,ny,nz,ldx,ldy,ldz,mgfft,-1,gbound)
   else
     call fourdp_c2c_ip(uprod,-1,MPI_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
   end if

   ! Map FFT to G-sphere.
   do ig=1,npw      
     fft_idx=igfft(ig)
     if (fft_idx/=0) then ! G-G0 belong to the FFT mesh.
       rhotwg(ig)=uprod(fft_idx) 
     else               ! Set this component to zero.
       rhotwg(ig)=czero
     end if
   end do

   ABI_DEALLOCATE(uprod)

   RETURN

 CASE DEFAULT 
   MSG_BUG('Wrong nspinor')
 END SELECT

end subroutine ostrength
!!***

!----------------------------------------------------------------------

END MODULE m_oscillators
!!***
