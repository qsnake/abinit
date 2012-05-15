!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fftw3
!! NAME
!! m_fftw3
!!
!! FUNCTION
!!  This module provides wrappers for the FFTW3 routines: in-place and out-of-place version.
!!  These procedures are mainly used in the GW part in which wavefunctions and other 
!!  matrix elements are stored in complex arrays instead of real arrays with real and imaginary part
!!  as usually done in abinit. It also provides F90 wrappers for complex to complex FFTW3 routines (if available). 
!!  In this case no conversion between the complex and the real(2,:) storage is done thus leading to an
!!  additional speed-up of the calculation. If FFTW3 is not available, we fall back to fourdp, and the wrappers 
!!  take care of the conversion between the complex and the real storage-mode.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  1) MPI parallelism is not supported 
!!  2) For better performance the FFT divisions should contain small factors  (/2, 3, 5, 7, 11/).
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! It seems that MKL wrappers do not like the advanced interfaces for 
! r2c and c2r transforms although they work fine if the true FFTW3 library is used.
#define HAVE_FFT_FFTW3_MKL

MODULE m_fftw3

 use m_profiling

 use defs_basis
 use m_errors
 use m_timer
#ifdef HAVE_OPENMP
 use omp_lib
#endif

 use m_numeric_tools,   only : imax_loc

 implicit none

 private

 public :: fftw3_cleanup        ! Reset FFTW to the pristine state it was in when you started your program, 
 public :: fftw3_set_nthreads   ! Set the number of threads you want FFTW3 to use when HAVE_FFT_FFTW3_THREADS is defined.
 public :: fftw3_show_nthreads  ! Write info on the status of the threaded FFTW3 version.
 public :: fftw3_r2c_op         ! Real to complex transform (out-of-place version).
 public :: fftw3_c2r_op         ! Complex to real transform (out-of-place version).
 public :: fftw3_c2c_op         ! complex to complex transform (out-of-place version).
 public :: fftw3_c2c_ip         ! complex to complex transform (in-place version).
 public :: fftw3_many_dft_op    ! Driver routine for many out-of-place 3D complex-to-complex FFTs.
 public :: fftw3_many_dft_ip    ! Driver routine for many in-place 3D complex-to-complex FFTs.
 public :: fftw3_fftpad         ! Driver routines for zero-padded FFT of wavefunctions.
 public :: fftw3_fftpad_cplx    ! Driver routines for zero-padded FFT of wavefunctions.
 public :: fftw3_gain_wisdom    ! DO NOT USE! Select an optimal FFTW3 plan for the different transforms.
 public :: fftw3_set_in_place   ! Setup an internal flag to use in-place FFTs wherever possible.
 public :: fftw3_is_in_place    ! .TRUE. if in-place FFTs are used wherever possible.

 public :: fftw3_fftpad_tr      ! Still under development.

! flags copied from fftw3.f
 integer,public,parameter :: FFTW_R2HC=0
 integer,public,parameter :: FFTW_HC2R=1
 integer,public,parameter :: FFTW_DHT=2
 integer,public,parameter :: FFTW_REDFT00=3
 integer,public,parameter :: FFTW_REDFT01=4
 integer,public,parameter :: FFTW_REDFT10=5
 integer,public,parameter :: FFTW_REDFT11=6
 integer,public,parameter :: FFTW_RODFT00=7
 integer,public,parameter :: FFTW_RODFT01=8
 integer,public,parameter :: FFTW_RODFT10=9
 integer,public,parameter :: FFTW_RODFT11=10
 integer,public,parameter :: FFTW_FORWARD=-1
 integer,public,parameter :: FFTW_BACKWARD=+1
 integer,public,parameter :: FFTW_MEASURE=0
 integer,public,parameter :: FFTW_DESTROY_INPUT=1
 integer,public,parameter :: FFTW_UNALIGNED=2
 integer,public,parameter :: FFTW_CONSERVE_MEMORY=4
 integer,public,parameter :: FFTW_EXHAUSTIVE=8
 integer,public,parameter :: FFTW_PRESERVE_INPUT=16
 integer,public,parameter :: FFTW_PATIENT=32
 integer,public,parameter :: FFTW_ESTIMATE=64
 integer,public,parameter :: FFTW_ESTIMATE_PATIENT=128
 integer,public,parameter :: FFTW_BELIEVE_PCOST=256
 integer,public,parameter :: FFTW_NO_DFT_R2HC=512
 integer,public,parameter :: FFTW_NO_NONTHREADED=1024
 integer,public,parameter :: FFTW_NO_BUFFERING=2048
 integer,public,parameter :: FFTW_NO_INDIRECT_OP=4096
 integer,public,parameter :: FFTW_ALLOW_LARGE_GENERIC=8192
 integer,public,parameter :: FFTW_NO_RANK_SPLITS=16384
 integer,public,parameter :: FFTW_NO_VRANK_SPLITS=32768
 integer,public,parameter :: FFTW_NO_VRECURSE=65536
 integer,public,parameter :: FFTW_NO_SIMD=131072
 integer,public,parameter :: FFTW_NO_SLOW=262144
 integer,public,parameter :: FFTW_NO_FIXED_RADIX_LARGE_N=524288
 integer,public,parameter :: FFTW_ALLOW_PRUNING=1048576
 integer,public,parameter :: FFTW_WISDOM_ONLY=2097152
! end flags copied from fftw3.f

! ==========================================================================================
! ==== Variables introduced for the FFTW3 interface in abinit. Not belonging to fftw3.f ====
! ==========================================================================================

 integer,public,parameter :: ABINIT_FFTW_CLEANUP = -1
 ! Flag used to deallocate the plans saved by the wrappers. Cannot be mixed with other flags.

 integer,public,parameter :: NULL_PLAN = 0
 ! MKL wrappers might return NULL_PLAN if a particular FFTW3 feature is not available

 integer,public,parameter :: KIND_FFTW_PLAN = 8
 ! It should be at least integer*@SIZEOF_INT_P@
 ! MKL wrappers requires it to be integer*8, so do _not_ use C_INTPTR_T.

 integer,private,save :: THREADS_INITED = 0
 ! 1 if treads have been initialized. 0 otherwise.

 integer,private,save :: ABINIT_FFTW_NTHREADS = 2
 ! The number of threads used at run-time (read from ABINIT_FFTW_NTHREADS)

 integer,private,parameter :: ABINIT_FFTW_NTHREADS_DEFAULT = 2
 ! The default number of threads (used if get_environment_variable is not available)

 integer,private,parameter :: MPLANES = 6
 ! Up to MPLANES initializations (for different combinations of input parameters)
 ! are stored and re-used if available.
!!***

!----------------------------------------------------------------------

!!****t* m_fftw3/fftw3_plan3_t
!! NAME
!! fftw3_plan3_t     
!! 
!! FUNCTION
!!  Structure storing the pointer to the FFTW plan as well as the options used to generate it. 
!! 
!! SOURCE

 type,private :: fftw3_plan3_t
   integer :: isign=0                           ! Sign of the exponential in the FFT
   integer :: ndat=-1                           ! Number of FFTs associated to the plan
   integer :: flags=-HUGE(0)                    ! FFTW3 flags used to construct the plan.
   integer(KIND_FFTW_PLAN) :: plan=NULL_PLAN    ! FFTW3 plan.
   integer :: nthreads=1                        ! The number of threads associated to the plan.
   integer :: idist=-1 
   integer :: odist=-1
   integer :: istride=-1
   integer :: ostride=-1
   integer :: n(3)=-1                           ! The number of FFT divisions.
   integer :: inembed(3)=-1
   integer :: onembed(3)=-1
   !integer(C_INT) :: alignment(2)              ! The alignment of the arrays used to construct the plan.
 end type fftw3_plan3_t
!!***

! private Variables
 !logical,private,save :: FFTW3_IN_PLACE=.FALSE.
 logical,private,save :: FFTW3_IN_PLACE=.TRUE.
 ! Defines whether FFT is done in place or not. 
 ! TODO should be Initialized from ABINIT_FFTW3_IN_PLACE

CONTAINS  !===========================================================

!!****f* m_fftw3/fftw3_c2c_ip
!! NAME
!!  fftw3_c2c_ip
!!
!! FUNCTION
!! Driver routine for in-place 3D complex-complex FFT.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => ff(R); -1 : ff(R) => ff(G)
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator. 
!!   Defaults to FFTW_ESTIMATE.
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz)=
!!    In input: the complex array to be transformed.
!!    In output: the Fourier transformed in the space specified by isign.
!!
!! PARENTS
!!      fftw3_fourdp,m_fftw3
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_c2c_ip(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff,fftw_flags)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_c2c_ip'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 complex(dpc),intent(inout) :: ff(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFT_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank=3
 integer,save :: iplan2insert = 1
 integer :: my_flags,dist,ii,stride
 integer(KIND_FFTW_PLAN) :: my_plan 
 real(dp) :: nfact
!arrays
 integer :: embed(rank),n(rank)
 type(fftw3_plan3_t),save :: Saved_plans(MPLANES)

! *************************************************************************

 my_flags=FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags=fftw_flags

 if (my_flags == ABINIT_FFTW_CLEANUP) then
  call destroy_plans(Saved_plans); RETURN
 end if

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz /) 

 my_plan = retrieve_plan3(n,ndat,embed,stride,dist,embed,stride,dist,isign,my_flags,Saved_plans)

 if (my_plan == NULL_PLAN) then ! No plan exist for these parameters, so initialize a new one.
   my_plan = zplan_many_dft(rank, n, ndat, ff, embed, stride, dist, ff, embed, stride, dist, isign, my_flags)

   ii = insert_plan3(iplan2insert,my_plan,n,ndat,embed,stride,dist,embed,stride,dist,isign,my_flags,Saved_plans)
   iplan2insert = MOD(iplan2insert, MPLANES) + 1
 end if

 ! Now perform the 3D FFT via FFTW.
 call dfftw_execute_dft(my_plan, ff, ff)

 if (isign==FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
  nfact = one / DBLE(nx*ny*nz)
  call ZDSCAL(ldx*ldy*ldz*ndat, nfact, ff, 1) 
 end if

#else 
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,isign/))
 ABI_UNUSED(ff)
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_c2c_ip
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_c2c_op
!! NAME
!!  fftw3_c2c_op
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-complex FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => gg(R); -1 : ff(R) => gg(G)
!! ff(ldx*ldy*ldz)=The array to be transformed.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator. 
!!   Defaults to FFTW_ESTIMATE.
!!
!! OUTPUT 
!! gg(ldx*ldy*ldz*ndat)=The FFT of ff.
!!
!! PARENTS
!!      fftw3_fourdp,m_fftw3
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_c2c_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff,gg,fftw_flags)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_c2c_op'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,isign,ndat
 integer,optional,intent(in) :: fftw_flags
!arrays
 complex(dpc),intent(in) :: ff(ldx*ldy*ldz*ndat)
 complex(dpc),intent(out) :: gg(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFT_FFTW3 
!Local variables-------------------------------
!scalars
 integer,parameter :: rank=3
 integer,save :: iplan2insert = 1
 integer :: my_flags,dist,ii,stride
 integer(KIND_FFTW_PLAN) :: my_plan
 real(dp) :: nfact
!arrays
 integer :: embed(rank),n(rank)
 type(fftw3_plan3_t),save :: Saved_plans(MPLANES)

! *************************************************************************

 my_flags=FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 if (my_flags == ABINIT_FFTW_CLEANUP) then
  call destroy_plans(Saved_plans); RETURN
 end if

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/) 
 n      = (/nx ,ny ,nz/) 

 my_plan = retrieve_plan3(n,ndat,embed,stride,dist,embed,stride,dist,isign,my_flags,Saved_plans)

 if (my_plan == NULL_PLAN) then ! No plan exist for these parameters, so initialize a new one.
   my_plan = zplan_many_dft(rank, n, ndat, ff, embed, stride, dist, gg,embed,stride, dist, isign, my_flags)

   ii = insert_plan3(iplan2insert,my_plan,n,ndat,embed,stride,dist,embed,stride,dist,isign,my_flags,Saved_plans)
   iplan2insert = MOD(iplan2insert, MPLANES) + 1
 end if

 ! Now perform the 3D FFT via FFTW.
 call dfftw_execute_dft( my_plan, ff, gg)

 if (isign==FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
  nfact = one / DBLE(nx*ny*nz)
  call ZDSCAL( ldx*ldy*ldz*ndat, nfact, gg, 1)  
 end if

#else 
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,isign/))
 ABI_UNUSED(ff)
 ABI_UNUSED(gg)
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_c2c_op
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_r2c_op
!! NAME
!!  fftw3_r2c_op
!!
!! FUNCTION
!! Driver routine for out-of-place 3D real-to-complex FFT of lengths nx, ny, nz. 
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the f array (to avoid cache conflicts).
!! ff(nx*ny*nz*ndat)=The real array to be transformed.
!! ndat=Number of FFTs to be done.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator. 
!!   Defaults to FFTW_ESTIMATE.
!!
!! OUTPUT 
!! gg(2,nx*ny*nz*ndat)=The forward FFT of ff.
!!
!! NOTES
!!  FIXME For the time-being. No augmentation of the mesh to reduce memory conflicts, as MKL crashes 
!!  if the advanced interface is used.
!!
!! PARENTS
!!      fftw3_fourdp,m_fftw3
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_r2c_op(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg,fftw_flags)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_r2c_op'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(in) :: ff(nx*ny*nz*ndat)
 real(dp),intent(out) :: gg(2,nx*ny*nz*ndat)

#ifdef HAVE_FFT_FFTW3 
!Local variables-------------------------------
!scalars
 integer,parameter :: rank=3
 integer,save :: iplan2insert = 1
 integer :: nhp,my_flags,idist,odist,padx,i1,i2,i3,igp,igf,imgf,ii,stride
 integer :: i1inv,i2inv,i3inv,idat,padatf     
 integer(KIND_FFTW_PLAN) :: my_plan
 real(dp) :: nfact
!arrays
 integer :: inembed(rank),onembed(rank),n(rank)
 integer,allocatable :: i1inver(:),i2inver(:),i3inver(:)
 real(dp),allocatable :: gg_hp(:,:),gg_test(:,:)
 type(fftw3_plan3_t),save :: Saved_plans(MPLANES)

! *************************************************************************

 my_flags=FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 if (my_flags == ABINIT_FFTW_CLEANUP) then
  call destroy_plans(Saved_plans); RETURN
 end if

 idist = nx*ny*nz
 odist = nhp
                            
 stride = 1
 n      = (/nx,ny,nz/) 
 inembed= (/nx,ny,nz/) 
 onembed= (/(nx/2+1),ny,nz/)

 my_plan = retrieve_plan3(n,ndat,inembed,stride,idist,onembed,stride,odist,FFTW_FORWARD,my_flags,Saved_plans)

 nhp = (nx/2+1)*ny*nz
 ABI_ALLOCATE(gg_hp,(2,nhp*ndat))

 if (my_plan == NULL_PLAN) then ! No plan exist for these parameters, so initialize a new one.

#ifdef HAVE_FFT_FFTW3_MKL    
   if (ndat/=1) MSG_ERROR("ndat/=1 + MKL not coded")
   call dfftw_plan_dft_r2c_3d ( my_plan, nx, ny, nz, ff, gg_hp, my_flags)
   if (my_plan==NULL_PLAN) then 
     MSG_ERROR("dfftw_plan_dft_r2c_3d returned NULL_PLAN")
   end if

   !fftw_plan fftw_plan_many_dft_r2c(int rank, const int *n, int howmany,
   !  double *in, const int *inembed, int istride, int idist,
   !  fftw_complex *out, const int *onembed, int ostride, int odist, unsigned flags);
#else
   my_plan = dplan_many_dft_r2c(rank, n, ndat, ff, inembed, stride, idist, gg_hp, onembed, stride, odist, my_flags)
#endif

   ii = insert_plan3(iplan2insert,my_plan,n,ndat,inembed,stride,idist,onembed,stride,odist,FFTW_FORWARD,my_flags,Saved_plans)
   iplan2insert = MOD(iplan2insert, MPLANES) + 1
 end if

 ! Now perform the 3D FFT via FFTW. r2c are always FFTW_FORWARD
 call dfftw_execute_dft_r2c(my_plan, ff, gg_hp)

 nfact = one / DBLE(nx*ny*nz)
 call ZDSCAL( nhp*ndat, nfact, gg_hp, 1)  ! FFTW returns not normalized FTs
 
 ! Reconstruct full FFT: Hermitian redundancy: out[i] is the conjugate of out[n-i]
 padx = (nx/2+1)

 ABI_ALLOCATE(i1inver,(padx))
 ABI_ALLOCATE(i2inver,(ny))
 ABI_ALLOCATE(i3inver,(nz))
 i1inver(1)=1
 do i1=2,padx
   i1inver(i1)=nx+2-i1
 end do

 i2inver(1)=1
 do i2=2,ny
   i2inver(i2)=ny+2-i2
 end do

 i3inver(1)=1
 do i3=2,nz
   i3inver(i3)=nz+2-i3
 end do

 igp=0
 do idat=1,ndat
   padatf=(idat-1)*nx*ny*nz
   do i3=1,nz
     i3inv = i3inver(i3)
     do i2=1,ny
       i2inv = i2inver(i2)
       do i1=1,padx
         igp=igp+1
         igf = i1 + (i3-1)*nx*ny + (i2-1)*nx + padatf
         gg(:,igf) =  gg_hp(:,igp)
         i1inv = i1inver(i1)
         if (i1inv/=i1) then
           imgf = i1inv + (i3inv-1)*nx*ny + (i2inv-1)*nx + padatf
           gg(1,imgf) =  gg_hp(1,igp)
           gg(2,imgf) = -gg_hp(2,igp)
         end if
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(i1inver)
 ABI_DEALLOCATE(i2inver)
 ABI_DEALLOCATE(i3inver)

 ABI_DEALLOCATE(gg_hp)

#define DEV_MG_DEBUG_THIS 0
#if DEV_MG_DEBUG_THIS == 1
! Test R2C

 ABI_ALLOCATE(gg_hp,(2,nx*ny*nz*ndat))
 ABI_ALLOCATE(gg_test,(2,nx*ny*nz*ndat))

 gg_hp(1,:) = ff 
 gg_hp(2,:) = zero

 call fftw3_many_dft_op(nx,ny,nz,ldx,ldy,ldz,ndat,FFTW_FORWARD,gg_hp,gg_test)

 !call dfftw_plan_dft_3d ( my_plan, nx, ny, nz, gg_hp, gg_test, FFTW_FORWARD,FFTW_ESTIMATE)
 !if (my_plan==NULL_PLAN) then
 !  MSG_ERROR("dfftw_plan_dft_r2c_3d returned NULL_PLAN")
 !end if
 !call dfftw_execute_dft( my_plan, gg_hp, gg_test)
 !nfact = one / DBLE( nx*ny*nz )
 !call ZDSCAL( nx*ny*nz*ndat, nfact, gg_test, 1)  ! FFTW returns not normalized FTs
 !call dfftw_destroy_plan(my_plan)

 if (ANY(ABS(gg_test-gg)>tol6)) then

  write(std_out,*)MAXVAL(ABS(gg_test(1,:)-gg(1,:)))
  write(std_out,*)MAXVAL(ABS(gg_test(2,:)-gg(2,:)))
  ii = imax_loc(ABS(gg_test(1,:)-gg(1,:)))
  write(std_out,*)gg_test(1,ii),gg(1,ii)
  ii = imax_loc(ABS(gg_test(2,:)-gg(2,:)))
  write(std_out,*)gg_test(2,ii),gg(2,ii)

  if (.TRUE.) then
    write(777,*)"real version",nx,ny,nz,ndat
    do ii=1,nx*ny*nz*ndat
     write(777,*)ii,gg(:,ii)
    end do
                                                
    write(778,*)"complex version",nx,ny,nz,ndat
    do ii=1,nx*ny*nz*ndat
     write(778,*)ii,gg_test(:,ii)
    end do
  end if

  MSG_ERROR("gg_test -gg > tol6")
 end if

 gg = gg_test !hack to run MKL wrappers for fftw3
 ABI_DEALLOCATE(gg_hp)
 ABI_DEALLOCATE(gg_test)
#endif
! /*DEV_MG_DEBUG_THIS = 1*/

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz/))
 ABI_UNUSED(ff)
 ABI_UNUSED(gg(1,1))
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_r2c_op
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_c2r_op
!! NAME
!!  fftw3_c2r_op
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-to-real FFT of lengths nx, ny, nz. 
!!
!! INPUTS
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! ff(2,nx*ny*nz*ndat)=The complex array to be transformed.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator. 
!!   Defaults to FFTW_ESTIMATE.
!!
!! OUTPUT 
!! gg(2,nx*ny*nz*ndat)=The backwards real FFT of ff.
!!
!! NOTES
!!  FIXME For the time-being. No augmentation of the mesh to reduce memory conflicts, as MKL crashes 
!!  if the advanced interface is used.
!!
!! PARENTS
!!      fftw3_fourdp,m_fftw3
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_c2r_op(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg,fftw_flags)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_c2r_op'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(in) :: ff(2,nx*ny*nz*ndat)
 real(dp),intent(out) :: gg(nx*ny*nz*ndat)

#ifdef HAVE_FFT_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank=3
 integer,save :: iplan2insert = 1
 integer :: nhp,my_flags,padx,i2,i3,igp,igf,idat,padatf,padatp,ii
 integer :: idist,odist,stride
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: inembed(rank),onembed(rank),n(rank)
 real(dp),allocatable :: ff_hp(:,:)
 type(fftw3_plan3_t),save :: Saved_plans(MPLANES)

! *************************************************************************

 my_flags=FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 if (my_flags == ABINIT_FFTW_CLEANUP) then
  call destroy_plans(Saved_plans); RETURN
 end if

 stride  = 1
 nhp     = (nx/2+1)*ny*nz
 idist   = nhp
 odist   = nx*ny*nz
 n       = (/nx,ny,nz/)
 inembed = (/(nx/2+1),ny,nz/)
 onembed = (/nx,ny,nz/) ! check this

 my_plan = retrieve_plan3(n,ndat,inembed,stride,idist,onembed,stride,odist,FFTW_BACKWARD,my_flags,Saved_plans)

 ! Fill the Hermitian part: Hermitian redundancy: out[i] is the conjugate of out[n-i]
 ABI_ALLOCATE(ff_hp,(2,nhp*ndat))

 padx = (nx/2+1)
 do idat=1,ndat
   padatf=(idat-1)*nx  *ny*nz
   padatp=(idat-1)*padx*ny*nz
   do i3=1,nz
     do i2=1,ny
       igf = (i3-1)*nx  *ny + (i2-1)*nx    + padatf
       igp = (i3-1)*padx*ny + (i2-1)*padx  + padatp
       ff_hp(:,igp+1:igp+padx) = ff(:,igf+1:igf+padx)
     end do
   end do
 end do

 ! NOTE: The c2r transform destroys its input array even for out-of-place transforms.
 if (my_plan == NULL_PLAN) then ! No plan exist for these parameters, so initialize a new one.

#ifdef HAVE_FFT_FFTW3_MKL
   if (ndat/=1) MSG_ERROR("ndat/=1 + MKL not coded")
   call dfftw_plan_dft_c2r_3d( my_plan, nx, ny, nz, ff_hp, gg, my_flags)
   if (my_plan==NULL_PLAN) then 
     MSG_ERROR("dfftw_plan_dft_c2r_3d returned NULL_PLAN")
   end if
#else
   my_plan = dplan_many_dft_c2r(rank, n, ndat, ff_hp, inembed, stride, idist, gg, onembed, stride, odist, my_flags)
#endif

   ii = insert_plan3(iplan2insert,my_plan,n,ndat,inembed,stride,idist,onembed,stride,odist,FFTW_BACKWARD,my_flags,Saved_plans)
   iplan2insert = MOD(iplan2insert, MPLANES) + 1
 end if

 ! Now perform the 3D FFT via FFTW. c2r are always FFTW_BACKWARD
 call dfftw_execute_dft_c2r( my_plan, ff_hp, gg)

 ABI_DEALLOCATE(ff_hp)

#else 
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz/))
 ABI_UNUSED(ff(1,1))
 ABI_UNUSED(gg(1))
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_c2r_op
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_many_dft_op
!! NAME
!!  fftw3_many_dft_op
!!
!! FUNCTION
!! Driver routine for many out-of-place 3D complex-to-complex FFTs of lengths nx, ny, nz. 
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimension of the fin and fout arrays (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! fin(2*ldx*ldy*ldz*ndat)=The complex array to be transformed.
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r, 
!!   -1 for transforming from r to G.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator. 
!!   Defaults to FFTW_ESTIMATE.
!!
!! OUTPUT 
!! fout(2,ldx*ldy*ldz*ndat)=The Fourier transform of fin.
!!
!! PARENTS
!!      fftw3_fourdp,m_fftw3
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_many_dft_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fin,fout,fftw_flags)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_many_dft_op'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(in) :: fin(2*ldx*ldy*ldz*ndat)
 real(dp),intent(out) :: fout(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_FFT_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank=3
 integer,save :: iplan2insert = 1
 integer :: my_flags,dist,ii,stride
 integer(KIND_FFTW_PLAN) :: my_plan
 real(dp) :: nfact
!arrays
 integer :: embed(rank),n(rank)
 type(fftw3_plan3_t),save :: Saved_plans(MPLANES)

! *************************************************************************

 !TODO this one has to be tested!!!!
 my_flags=FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 if (my_flags == ABINIT_FFTW_CLEANUP) then
  call destroy_plans(Saved_plans); RETURN
 end if

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz /)

 my_plan = retrieve_plan3(n,ndat,embed,stride,dist,embed,stride,dist,isign,my_flags,Saved_plans)

 if (my_plan == NULL_PLAN) then ! No plan exist for these parameters, so initialize a new one.
   my_plan = dplan_many_dft(rank, n, ndat, fin, embed, stride, dist, fout,embed, stride, dist, isign, my_flags)

   ii = insert_plan3(iplan2insert,my_plan,n,ndat,embed,stride,dist,embed,stride,dist,isign,my_flags,Saved_plans)
   iplan2insert = MOD(iplan2insert, MPLANES) + 1
 end if

 ! Now perform the 3D FFT via FFTW.
 call dfftw_execute_dft(my_plan, fin, fout)

 if (isign==FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
  nfact = one / DBLE(nx*ny*nz)
  call ZDSCAL(ldx*ldy*ldz*ndat, nfact, fout, 1) 
 end if

#else 
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,isign/))
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(fin(1))
 ABI_UNUSED(fout(1))
#endif

end subroutine fftw3_many_dft_op
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_many_dft_ip
!! NAME
!!  fftw3_many_dft_ip
!!
!! FUNCTION
!! Driver routine for many in-place 3D complex-to-complex FFTs of lengths nx, ny, nz. 
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimension of the finout array (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r, 
!!   -1 for transforming from r to G.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator. 
!!   Defaults to FFTW_ESTIMATE.
!!
!! OUTPUT 
!! finout(2,ldx*ldy*ldz*ndat)=
!!   In input: The complex array to be transformed.
!!   In output: The FFT results.
!!
!! PARENTS
!!      fftw3_fourdp
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,ndat,isign,finout,fftw_flags)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_many_dft_ip'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(inout) :: finout(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_FFT_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank=3
 integer,save :: iplan2insert = 1
 integer :: my_flags,dist,ii,stride
 integer(KIND_FFTW_PLAN) :: my_plan
 real(dp) :: nfact
!arrays
 integer :: embed(rank),n(rank)
 type(fftw3_plan3_t),save :: Saved_plans(MPLANES)

! *************************************************************************

 !TODO this one has to be tested!!!!
 my_flags=FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 if (my_flags == ABINIT_FFTW_CLEANUP) then
   call destroy_plans(Saved_plans); RETURN
 end if

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz /)

 my_plan = retrieve_plan3(n,ndat,embed,stride,dist,embed,stride,dist,isign,my_flags,Saved_plans)

 if (my_plan == NULL_PLAN) then ! No plan exist for these parameters, so initialize a new one.
   my_plan = dplan_many_dft(rank, n, ndat, finout, embed, stride, dist, finout,embed, stride, dist, isign, my_flags)

   ii = insert_plan3(iplan2insert,my_plan,n,ndat,embed,stride,dist,embed,stride,dist,isign,my_flags,Saved_plans)
   iplan2insert = MOD(iplan2insert, MPLANES) + 1
 end if

 ! Now perform the 3D FFT via FFTW.
 call dfftw_execute_dft(my_plan, finout, finout)

 if (isign==FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
  nfact = one / DBLE(nx*ny*nz)
  call ZDSCAL(ldx*ldy*ldz*ndat, nfact, finout, 1) 
 end if

#else 
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,isign/))
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(finout(1))
#endif

end subroutine fftw3_many_dft_ip
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_gain_wisdom
!! NAME
!!  fftw3_gain_wisdom
!!
!! FUNCTION
!!  Initialize optimal FFTW3 plans for the different FFTs algorithms characterized by 
!!  the number of divisions (nx,ny,nz) and ndat.
!!
!! INPUTS
!! nx,ny,nz=Number of point along the three directions.
!! ndat=Number of FFTs to be done.
!!
!! SIDE EFFECTS
!!  Optimal FFTW3 plans are initialized and stored inside each FFTW3 wrapper so 
!!  that they can be subsequently re-used. Note that both FFTW_ESTIMATE and FFTW_PATIENT
!!  overwrite their input. For this reason a separate initialization routine has to be called.
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_gain_wisdom(nx,ny,nz,ndat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_gain_wisdom'
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ndat
!arrays

#ifdef HAVE_FFT_FFTW3
!Local variables-------------------------------
!scalars
 integer :: wisdom_flag,isign,ii
 real(dp) :: cpu0,wall0,cpu,wall
 character(len=500) :: msg
!arrays
 integer :: signs(2)=(/FFTW_FORWARD,FFTW_BACKWARD/)
 real(dp),allocatable :: fin(:),fout(:),ff(:),gg(:,:)
 complex(dpc),allocatable :: ff_cplx(:),gg_cplx(:) 

! *************************************************************************

 call wrtout(std_out," Gaining FFTW3 wisdom. It might take some time...","COLL")
 call timein(cpu0,wall0)

 !wisdom_flag = FFTW_ESTIMATE
 wisdom_flag = FFTW_MEASURE
 !wisdom_flag = FFTW_PATIENT

 ABI_ALLOCATE(ff_cplx,(nx*ny*nz*ndat))
 ABI_ALLOCATE(gg_cplx,(nx*ny*nz*ndat))

 do ii=1,2
   isign = signs(ii)
   ff_cplx = cone
   call fftw3_c2c_ip(nx,ny,nz,nx,ny,nz,ndat,isign,ff_cplx,        fftw_flags=wisdom_flag)
   ff_cplx = cone
   call fftw3_c2c_op(nx,ny,nz,nx,ny,nz,ndat,isign,ff_cplx,gg_cplx,fftw_flags=wisdom_flag)
 end do

 write(msg,'(a)')" fftw3_c2c done"
 call wrtout(std_out,msg,"COLL")

 ABI_DEALLOCATE(ff_cplx)
 ABI_DEALLOCATE(gg_cplx)

 ABI_ALLOCATE(fin ,(2*nx*ny*nz*ndat))
 ABI_ALLOCATE(fout,(2*nx*ny*nz*ndat))

 do ii=1,2
   isign = signs(ii)
   fin=one
   call fftw3_many_dft_op(nx,ny,nz,nx,ny,nz,ndat,isign,fin,fout,fftw_flags=wisdom_flag)
 end do
 write(msg,'(a)')" fftw3_many_dft_op done"
 call wrtout(std_out,msg,"COLL")

 ABI_DEALLOCATE(fin)
 ABI_DEALLOCATE(fout)

 ABI_ALLOCATE(ff,(nx*ny*nz*ndat))
 ABI_ALLOCATE(gg,(2,nx*ny*nz*ndat))

 ! These calls seem to make the code stuck if we link against MKL
 !% call fftw3_r2c_op(nx,ny,nz,nx,ny,nz,ndat,ff,gg,fftw_flags=wisdom_flag)
 !% call fftw3_c2r_op(nx,ny,nz,nx,ny,nz,ndat,ff,gg,fftw_flags=wisdom_flag)
 !% write(msg,'(a)')" fftw3_c2r_op done"
 !% call wrtout(std_out,msg,"COLL")

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(gg)

 call timein(cpu,wall)

 cpu  = cpu - cpu0
 wall = wall - wall0
 write(msg,'(2(a,f6.3),a)')" Gained wisdom in cpu time= ",cpu," sec;  wall time= ",wall, " sec."
 call wrtout(std_out,msg,"COLL")

#else 
 ABI_UNUSED((/nx,ny,nz,ndat/))
 MSG_ERROR("FFTW3 support not activated")
#endif

end subroutine fftw3_gain_wisdom
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/retrieve_plan3
!! NAME
!!  retrieve_plan3
!!
!! FUNCTION
!!  Returns a previously generated FFTW plan that can be safely re-used to perform the FFT.
!!
!! INPUTS
!!  n(3)=FFT divisions used for the present FFT.
!!  ndat=The number of FFTS to be done.
!!  isign=The sign of the present FFT.
!!  fftw_flags=The flags used for the present plan.
!!  Saved_plans(MPLANES)<fftw3_plan3_t>=Array storing the previously created FFTW plans.
!!
!! OUTPUT
!!  retrieve_plan3=The previously generated FFTW plane compatible with the input arguments.
!!    Set tp NULL_PLAN if no plane can be re-used.
!!
!! NOTES
!!  A previously generated plane can be reused for performing similar FFTS provided that 
!!  the following conditions are fulfilled (page 37 of FFTW3 guide)
!!  
!!  1) Array sizes and strides are the same as that of the arrays used to generate the plan.
!!  2) The input and output arras are the same (in-place) or different (out-of-place)
!!  3) Alignment of the arrays is the same as that of the input/output arrays used to generate the plan.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function retrieve_plan3(n,ndat,inembed,istride,idist,onembed,ostride,odist,isign,fftw_flags,Saved_plans)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'retrieve_plan3'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,istride,ostride,idist,odist,isign
 integer,intent(in) :: fftw_flags
 integer(KIND_FFTW_PLAN) :: retrieve_plan3
!arrays
 integer,intent(in) :: n(3),inembed(3),onembed(3)
 type(fftw3_plan3_t),intent(in) ::  Saved_plans(MPLANES)

!Local variables-------------------------------
!scalars
 integer :: ip

! *************************************************************************

 !@fftw3_plan3_t
 retrieve_plan3=NULL_PLAN; if (MPLANES<=0) RETURN

 ! Check if there is already a plan initialized for this combination of parameters.
 do ip=1,MPLANES
   if ( Saved_plans(ip)%isign    == isign      .and.           &
&       Saved_plans(ip)%ndat     == ndat       .and.           &
&       Saved_plans(ip)%flags    == fftw_flags .and.           &
&       Saved_plans(ip)%nthreads == ABINIT_FFTW_NTHREADS .and. &
&       Saved_plans(ip)%idist    == idist .and.                &
&       Saved_plans(ip)%odist    == odist .and.                &
&       Saved_plans(ip)%istride  == istride .and.              &
&       Saved_plans(ip)%ostride  == ostride .and.              &
&   ALL(Saved_plans(ip)%n        == n) .and.                   & 
&   ALL(Saved_plans(ip)%inembed  == inembed) .and.             & 
&   ALL(Saved_plans(ip)%onembed  == onembed) ) then
     retrieve_plan3=Saved_plans(ip)%plan; EXIT
   end if
 end do

end function retrieve_plan3
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/insert_plan3
!! NAME
!!  insert_plan3
!!
!! FUNCTION
!!  Insert the FFTW plan at the location specified by ip. 
!!  Returns 
!!   0 if not plan has been inserted
!!   1 if an old plan has been inserted the make room for the new one.
!!   2 if an empty slot has been used to save the new plan.
!!
!! INPUTS
!!  ip=The index in the array Saved_plans used to insert the new plan. 
!!  plan=The FFTW plan.
!!  n(3)=FFT divisions used for the present FFT.
!!  ndat=The number of FFTS to be done.
!!  isign=The sign of the present FFT.
!!  fftw_flags=The flags used for the present plan.
!!
!! SIDE EFFECTS
!!  Saved_plans(MPLANES)<fftw3_plan3_t>=Array storing the previously created FFTW plans.
!!    Saved_plans(ip) is overwritten with the new plan.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function insert_plan3(ip,plan,n,ndat,inembed,istride,idist,onembed,ostride,odist,isign,fftw_flags,Saved_plans)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'insert_plan3'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,istride,ostride,idist,odist,isign,ip
 integer,intent(in) :: fftw_flags
 integer(KIND_FFTW_PLAN),intent(in) :: plan
!arrays
 integer,intent(in) :: n(3),inembed(3),onembed(3)
 type(fftw3_plan3_t),intent(inout) ::  Saved_plans(MPLANES)

!Local variables-------------------------------
!scalars
 integer :: insert_plan3

! *************************************************************************

 !@fftw3_plan3_t
 insert_plan3 = 0; if (MPLANES<=0) RETURN

 ! Destroy the plan if it already exists.
 if (Saved_plans(ip)%plan/= NULL_PLAN) then 
#ifdef HAVE_FFT_FFTW3     
   call dfftw_destroy_plan( Saved_plans(ip)%plan )
#else 
   MSG_ERROR("FFTW3 support not activated")
#endif
   insert_plan3 = 1
 else 
   insert_plan3 = 2
 end if

 ! Insert the new plan.
 Saved_plans(ip)%isign     = isign 
 Saved_plans(ip)%ndat      = ndat 
 Saved_plans(ip)%flags     = fftw_flags
 Saved_plans(ip)%plan      = plan
 Saved_plans(ip)%nthreads  = ABINIT_FFTW_NTHREADS
 Saved_plans(ip)%idist     = idist
 Saved_plans(ip)%odist     = odist
 Saved_plans(ip)%istride   = istride
 Saved_plans(ip)%ostride   = ostride
 Saved_plans(ip)%n         = n
 Saved_plans(ip)%inembed   = inembed
 Saved_plans(ip)%onembed   = onembed

end function insert_plan3
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/destroy_plan
!! NAME
!!  destroy_plans
!!
!! FUNCTION
!!   Free the FFTW plans stored in memory and reset the object.
!!
!! SIDE EFFECTS
!!  Saved_plans(:)<fftw3_plan3_t>=Array storing the previously created FFTW plans.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine destroy_plans(Saved_plans)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_plans'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!arrays
 type(fftw3_plan3_t),intent(inout) ::  Saved_plans(:)

!Local variables-------------------------------
!scalars
 integer :: ip

! *************************************************************************

 !@fftw3_plan3_t
 do ip=1,SIZE(Saved_plans)
   if (Saved_plans(ip)%plan/= NULL_PLAN) then  ! Destroy the plan if it already exists.
#ifdef HAVE_FFT_FFTW3     
     call dfftw_destroy_plan( Saved_plans(ip)%plan )
#else 
     MSG_ERROR("FFTW3 support not activated")
#endif
   end if
   ! Reset the plan.
   Saved_plans(ip)%isign     = 0
   Saved_plans(ip)%ndat      = -1 
   Saved_plans(ip)%flags     = -HUGE(0)
   Saved_plans(ip)%plan      = NULL_PLAN
   Saved_plans(ip)%nthreads  = 1
   Saved_plans(ip)%idist     =-1
   Saved_plans(ip)%odist     =-1
   Saved_plans(ip)%istride   =-1
   Saved_plans(ip)%ostride   =-1
   Saved_plans(ip)%n         =-1
   Saved_plans(ip)%inembed   =-1
   Saved_plans(ip)%onembed   =-1
 end do


end subroutine destroy_plans
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_cleanup
!! NAME
!!  fftw3_cleanup
!!
!! FUNCTION
!!  Reset FFTW to the pristine state it was in when you started your program, 
!!  All existing plans become undefined. 
!!
!! NOTES
!!  FFTW planner saves some other persistent data, such as the accumulated wisdom and a list of 
!!  algorithms available in the current configuration. If you want to deallocate all of that and reset 
!!  FFTW to the pristine state it was in when you started your program, you can call fftw3_cleanup();
!!  After calling fftw3_cleanup, all existing plans become undefined, and you should not attempt to 
!!  execute them nor to destroy them. You can however create and execute/destroy new plans, in which case 
!!  FFTW starts accumulating wisdom information again. 
!!  fftw3_cleanup does not deallocate your plans, however. To prevent memory leaks, you must still call 
!!  fftw_destroy_plan before executing fftw3_cleanup
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_cleanup()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_cleanup'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

!Local variables-------------------------------
!arrays
 real(dp) :: fin(0),fout(0),ff(0),gg(0,0)
 complex(dpc) :: ff_cplx(0),gg_cplx(0) 

! *************************************************************************

 call fftw3_c2c_ip(0,0,0,0,0,0,0,1,ff_cplx,        fftw_flags=ABINIT_FFTW_CLEANUP)
 call fftw3_c2c_op(0,0,0,0,0,0,0,1,ff_cplx,gg_cplx,fftw_flags=ABINIT_FFTW_CLEANUP)

 call fftw3_many_dft_op(0,0,0,0,0,0,0,1,fin,fout,fftw_flags=ABINIT_FFTW_CLEANUP)

 call fftw3_r2c_op(0,0,0,0,0,0,0,ff,gg,fftw_flags=ABINIT_FFTW_CLEANUP)
 call fftw3_c2r_op(0,0,0,0,0,0,0,ff,gg,fftw_flags=ABINIT_FFTW_CLEANUP)

#ifdef HAVE_FFT_FFTW3_THREADS 
 if (THREADS_INITED==1) then
   call dfftw_cleanup_threads()
   THREADS_INITED = 0
 end if
#elif defined HAVE_FFT_FFTW3
 call dfftw_cleanup()
#else
 MSG_ERROR("FFTW3 support not activated")
#endif

end subroutine fftw3_cleanup
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_set_nthreads
!! NAME
!!  fftw3_set_nthreads
!!
!! FUNCTION
!!  This function sets the number of threads you want FFTW3 to use (or actually, the maximum number). 
!!  It also performs any one-time initialization required to use FFTW3 threads.
!!  All plans subsequently created with any planner routine will use nthreads threads. 
!!  If you pass an nthreads argument of 1 (the default), threads are disabled for subsequent plans.
!!  It does nothing if HAVE_FFT_FFTW3_THREADS is not defined.
!!
!! INPUTS 
!!  [nthreads]=The number of threads you want FFTW3 to use. It not present, nthreads is 
!!    initialized from the environment variable ABINIT_FFTW_NTHREADS.
!!
!! SIDE EFFECTS
!!  The one-time initialization required to use FFTW3 threads is performed when the routine 
!!  is called for the first time.
!!
!! PARENTS
!!      driver,m_fft_prof
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_set_nthreads(nthreads)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_set_nthreads'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: nthreads

!Local variables ------------------------------
!scalars
#ifdef HAVE_FFT_FFTW3_THREADS 
 integer :: iret,istat
 character(len=500) :: msg,str_nthreads
#endif

! *************************************************************************

#ifdef HAVE_FFT_FFTW3_THREADS 
 if (THREADS_INITED==0) then
   call dfftw_init_threads(iret)
   if (iret==0) then
     MSG_WARNING(" dfftw_init_threads returned 0; threaded FFTW3 is not being used!")
     RETURN
   else
    THREADS_INITED=1
   end if
 end if

 ABINIT_FFTW_NTHREADS = ABINIT_FFTW_NTHREADS_DEFAULT  ! Default value
 if (PRESENT(nthreads)) then
   ABINIT_FFTW_NTHREADS = nthreads
 else
#ifdef HAVE_OPENMP
   ABINIT_FFTW_NTHREADS = omp_get_max_threads()
#else
   MSG_WARNING("Using FFTW3 with threads but HAVE_OPENMP is not defined!")
#endif
!! #ifdef HAVE_FC_GETENV
!!    call get_environment_variable('OMP_NUM_THREADS',VALUE=str_nthreads,STATUS=istat)
!!    if (istat==0) read(str_nthreads,'(i3)') ABINIT_FFTW_NTHREADS
!! #endif
 end if

#if 0
 write(msg,'(a,i0)')" Using threaded FFTW3 library with nthreads: ",ABINIT_FFTW_NTHREADS
 call wrtout(std_out,msg,"COLL")
#endif

 call dfftw_plan_with_nthreads(ABINIT_FFTW_NTHREADS)
 
#else
 if (PRESENT(nthreads)) then
   ABI_UNUSED(nthreads) 
 end if
#endif

end subroutine fftw3_set_nthreads
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_show_nthreads
!! NAME
!!  fftw3_show_nthreads
!!
!! FUNCTION
!!   Write on unit unit with mode mode_paral the number of threads used in FFTW3.
!!
!! INPUTS 
!!  [unit]=Unit number for output, defaults to std_out
!!  [mode_paral]= ("COLL" | "PERS"), default to "COLL"
!!
!! PARENTS
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_show_nthreads(unit,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_show_nthreads'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit
 character(len=4),optional,intent(in) :: mode_paral 

!Local variables ------------------------------
!scalars
 integer :: my_unt
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

#ifdef HAVE_FFT_FFTW3_THREADS 
 if (THREADS_INITED==0) then
   write(msg,'(a)')" WARNING: FFTW3 library supports threads, but THREADS_INITED==0 "
 else
   write(msg,'(a,i0)')" Using threaded FFTW3 library with nthreads: ",ABINIT_FFTW_NTHREADS
 end if
#else
   write(msg,'(a)')" FFTW3 library does not support threads. "
#endif

 call wrtout(my_unt,msg,my_mode)

end subroutine fftw3_show_nthreads
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_set_in_place
!! NAME
!!  fftw3_set_in_place
!!
!! FUNCTION
!!  Defines whether FFTW3 transforms have to be done in-place or out-of-place
!!
!! INPUTS 
!!  lflag=.TRUE. for in-place transforms.
!!
!! SIDE EFFECTS
!!  FFTW3_IN_PLACE set to lflag.
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_set_in_place(lflag)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_set_in_place'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: lflag

! *************************************************************************

  FFTW3_IN_PLACE = lflag

end subroutine fftw3_set_in_place
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_is_in_place
!! NAME
!!  fftw3_is_in_place
!!
!! FUNCTION
!!  Returns .TRUE. if FFTW3 transforms are done in place.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function fftw3_is_in_place()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_is_in_place'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: fftw3_is_in_place

! *************************************************************************

  fftw3_is_in_place = FFTW3_IN_PLACE

end function fftw3_is_in_place
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftpad
!! NAME
!!  fftw3_fftpad
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with FFTW3. 
!!  The 3D ffts are computed only on lines and planes which have non zero elements. 
!!  These lines and planes are defined by the two vectors do_fft_x(ldy*nz) and do_fft_y(nz) 
!!  FFT transform is in-place.
!!  
!! INPUTS
!!   logical dimensions of the fft physical dimensions of the f array sign of the transformation
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!   f(2*ldx*ldy*ldz)=
!!     input: The array with the data to be transformed.
!!     output: The results of the FFT.
!!
!! PARENTS
!!      fftw3_fourwf,m_wfs
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_fftpad(ff,nx,ny,nz,ldx,ldy,ldz,mgfft,isign,gbound)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_fftpad'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 real(dp),intent(inout) :: ff(2*ldx*ldy*ldz)

#ifdef HAVE_FFT_FFTW3
!Local variables-------------------------------
!scalars
 integer :: kk,idir,ip,ii,jj,sidx
 integer :: id_yz,g3_max,g3_min,len3,gg2,gg3,g2,ifft_g3,igb,g2min,g2max
 integer,save :: bw_icurrent = 1
 integer,save :: fw_icurrent = 1
!arrays
 integer,save :: fw_dims(6,MPLANES) = -1
 integer,save :: bw_dims(6,MPLANES) = -1
 integer(KIND_FFTW_PLAN),save :: fw_plan (3,MPLANES) = 0
 integer(KIND_FFTW_PLAN),save :: bw_plan (3,MPLANES) = 0
 integer,allocatable :: do_fft_x(:,:),do_fft_y(:)

! *************************************************************************

 ABI_ALLOCATE(do_fft_x,(ldy,ldz))
 ABI_ALLOCATE(do_fft_y,(nz))
 do_fft_y = 0                             ! we have to recalculate them at each call.
 do_fft_x = 0

 g3_min=gbound(3,2)
 g3_max=gbound(4,2)
 len3=g3_max-g3_min+1

 do gg3=1,len3 ! Loop over the z-planes intersecting the G-sphere.
 !
   if (gg3<=g3_max+1) then
     ifft_g3=gg3
   else 
     ifft_g3=gg3+nz-len3 ! Wrap around for negative gg3.
   end if
   do_fft_y(ifft_g3) = 1
   !
   ! Select the set of y for this z-plane.
   igb=2*gg3+3
   g2min = gbound(igb  ,2) 
   g2max = gbound(igb+1,2)
   do_fft_x(1:g2max+1    ,ifft_g3) = 1 ! Positive g_y.
   do_fft_x(g2min+ny+1:ny,ifft_g3) = 1 ! Negative g_y.
 end do
 !do_fft_y=1; do_fft_x=1
 !write(std_out,*)"do_fft_x",do_fft_x
 !write(std_out,*)"do_fft_y",do_fft_y
 !stop

 SELECT CASE (isign)

 CASE (FFTW_BACKWARD) ! G --> R
   !
   ! First check if there is already a table initialized for this combination of parameters
   ip = -1
   do ii=1,MPLANES
     if ( ALL( (/nx,ny,nz,ldx,ldy,ldz/) == bw_dims(:,ii) )) then
       ip=ii; EXIT
     end if
   end do

   if (ip==-1) then ! no table exist for these parameters initialize a new one.
     !
     ! Backward plans
     if ( bw_plan(1,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(1,bw_icurrent) )
     if ( bw_plan(2,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(2,bw_icurrent) )
     if ( bw_plan(3,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(3,bw_icurrent) )
     !
     ! The prototype for dfftw_plan_many_dft is:
     ! dfftw_plan_many_dft(rank, n, howmany, 
     !   fin,  iembed, istride, idist, 
     !   fout, oembed, ostride, odist, isign, my_flags)
     !     
     bw_plan(1, bw_icurrent) = dplan_many_dft(1, (/nx/), 1,       &   ! Single 1D transform of f(Gx,Gy,Gz) along Gx.
&         ff, (/ldx, ldy, ldz/), 1, ldx,                          &
&         ff, (/ldx, ldy, ldz/), 1, ldx, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_plan(2, bw_icurrent) = dplan_many_dft(1, (/ny/), nx,      &   ! nx 1D transforms of f(x,Gy,Gz) along Gy.
&         ff, (/ldx, ldy, ldz/), ldx, 1,                          &
&         ff, (/ldx, ldy, ldz/), ldx, 1, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_plan(3, bw_icurrent) = dplan_many_dft(1, (/nz/), ldx*ldy, & ! ldx*ldy 1D transforms of f(x,y,Gz) along Gz.
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1,                      & ! Note that we have to visit the entire augmented x-y plane!
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_dims(:,bw_icurrent) = (/nx,ny,nz,ldx,ldy,ldz/)
     ip = bw_icurrent ! Got new plan.
     bw_icurrent = MOD( bw_icurrent, MPLANES ) + 1
   end if
   !
   ! 1) Transform along x.
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       do jj=1,ny
         if (do_fft_x(jj,kk) == 1) then 
           sidx = 1 + 2*(jj-1)*ldx + 2*(kk-1)*ldx*ldy ! Pass the pointer, 2 to account for the imag part.
           call dfftw_execute_dft(bw_plan(1,ip), ff(sidx), ff(sidx) )
         end if
       end do
     end if
   end do
   !
   ! 2) Transform along y.
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       sidx = 1 + 2*(kk-1)*ldx*ldy
       call dfftw_execute_dft( bw_plan(2,ip), ff(sidx), ff(sidx) )
     end if
   end do
   !
   ! 3) Transform along z.
   call dfftw_execute_dft( bw_plan(3,ip), ff, ff)
 
 CASE (FFTW_FORWARD) ! R --> G
   !
   ! First check if there is already a table initialized for this combination of parameters
   ip = -1
   do ii=1,MPLANES
     if( ALL( (/nx,ny,nz,ldx,ldy,ldz/) == fw_dims(:,ii) )) then
       ip=ii; EXIT
     end if
   end do

   if (ip==-1) then ! no table exist for these parameters initialize a new one
     ! The prototype for dfftw_plan_many_dft is:
     ! dfftw_plan_many_dft(rank, n, howmany, 
     !   fin,  iembed, istride, idist, 
     !   fout, oembed, ostride, odist, isign, my_flags)
     !
     ! Forward plans
     if ( fw_plan(1,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(1,fw_icurrent) )
     if ( fw_plan(2,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(2,fw_icurrent) )
     if ( fw_plan(3,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(3,fw_icurrent) )

     fw_plan(1, fw_icurrent) = dplan_many_dft(1, (/nx/), 1,      &
&         ff, (/ldx, ldy, ldz/), 1, ldx,                         &
&         ff, (/ldx, ldy, ldz/), 1, ldx, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_plan(2, fw_icurrent) = dplan_many_dft(1, (/ny/), nx,     &
&         ff, (/ldx, ldy, ldz/), ldx, 1,                         &
&         ff, (/ldx, ldy, ldz/), ldx, 1, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_plan(3, fw_icurrent) = dplan_many_dft(1, (/nz/), ldx*ldy,& ! We have to visit the entire augmented x-y plane!
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1,                     &
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_dims(:,fw_icurrent) = (/nx,ny,nz,ldx,ldy,ldz/)
     ip = fw_icurrent ! Got new plan.
     fw_icurrent = MOD( fw_icurrent, MPLANES ) + 1
   end if
   !
   ! 1) Transform along z.
   call dfftw_execute_dft( fw_plan(3,ip), ff, ff)
   !
   ! 2) Transform along y.
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       sidx = 1 + 2*ldx*ldy*(kk-1)
       call dfftw_execute_dft( fw_plan(2,ip), ff(sidx), ff(sidx) )
     end if
   end do
   !
   ! 3) Transform along x. 
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       do jj=1,ny
         if (do_fft_x(jj,kk) == 1) then 
           sidx = 1 + 2*(jj-1)*ldx + 2*(kk-1)*ldx*ldy ! Pass the pointer, 2 to account for the imag part.
           call dfftw_execute_dft( fw_plan(1,ip), ff(sidx), ff(sidx) )
         end if
       end do
     end if
   end do
   !
   ! 4) Normalize the transform.
   !call DSCAL(2*ldx*ldy*nz, one/(nx*ny*nz), ff, 1)
   call ZDSCAL(ldx*ldy*nz, one/(nx*ny*nz), ff, 1)
 
 CASE DEFAULT 
   MSG_BUG("Wrong isign")
 END SELECT

 ABI_DEALLOCATE(do_fft_x)
 ABI_DEALLOCATE(do_fft_y)

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 ABI_UNUSED(ff(1))
#endif

end subroutine fftw3_fftpad     
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftpad_cplx        
!! NAME
!!  fftw3_fftpad_cplx
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with FFTW3. 
!!  The 3D ffts are computed only on lines and planes which have non zero elements. 
!!  These lines and planes are defined by the two vectors do_fft_x(ldy*nz) and do_fft_y(nz) 
!!  FFT transform is in-place. Target: complex arrays.
!!  
!! INPUTS
!!   logical dimensions of the fft physical dimensions of the f array sign of the transformation
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound.
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!  ff(2*ldx*ldy*ldz)=
!!    input: The array with the data to be transformed.
!!    output: The results of the FFT.
!!
!! PARENTS
!!      fftw3_fourwf,m_wfs
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_fftpad_cplx(ff,nx,ny,nz,ldx,ldy,ldz,mgfft,isign,gbound)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_fftpad_cplx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 complex(dpc),intent(inout) :: ff(ldx*ldy*ldz)

#ifdef HAVE_FFT_FFTW3
!Local variables-------------------------------
!scalars
 integer :: kk,idir,ip,ii,jj,sidx
 integer :: id_yz,g3_max,g3_min,len3,gg2,gg3,g2,ifft_g3,igb,g2min,g2max
 integer,save :: bw_icurrent = 1
 integer,save :: fw_icurrent = 1
!arrays
 integer,save :: fw_dims(6,MPLANES) = -1
 integer,save :: bw_dims(6,MPLANES) = -1
 integer(KIND_FFTW_PLAN),save :: fw_plan (3,MPLANES) = 0
 integer(KIND_FFTW_PLAN),save :: bw_plan (3,MPLANES) = 0
 integer,allocatable :: do_fft_x(:,:),do_fft_y(:)
 !type(fftw3_plan3_t),save :: Splans_z(MPLANES)
 !type(fftw3_plan3_t),save :: Splans_y(MPLANES)
 !type(fftw3_plan3_t),save :: Splans_z(MPLANES)

! *************************************************************************

 ABI_ALLOCATE(do_fft_x,(ldy,ldz))
 ABI_ALLOCATE(do_fft_y,(nz))
 do_fft_y = 0                             ! we have to recalculate them at each call.
 do_fft_x = 0

 g3_min=gbound(3,2)
 g3_max=gbound(4,2)
 len3=g3_max-g3_min+1

 do gg3=1,len3 ! Loop over the z-planes intersecting the G-sphere.

   if (gg3<=g3_max+1) then
     ifft_g3=gg3
   else 
     ifft_g3=gg3+nz-len3 ! Wrap around for negative gg3.
   end if
   do_fft_y(ifft_g3) = 1

   ! Select the set of y for this z-plane.
   igb=2*gg3+3
   g2min = gbound(igb  ,2) 
   g2max = gbound(igb+1,2)
   do_fft_x(1:g2max+1    ,ifft_g3) = 1 ! Positive g_y.
   do_fft_x(g2min+ny+1:ny,ifft_g3) = 1 ! Negative g_y.
 end do
 !do_fft_y=1; do_fft_x=1
 !write(std_out,*)"do_fft_x",do_fft_x
 !write(std_out,*)"do_fft_y",do_fft_y
 !stop

 SELECT CASE (isign)

 CASE (FFTW_BACKWARD) ! G --> R
   !
   ! First check if there is already a table initialized for this combination of parameters.
   ip = -1
   do ii=1,MPLANES
     if ( ALL( (/nx,ny,nz,ldx,ldy,ldz/) == bw_dims(:,ii) )) then
       ip=ii; EXIT
     end if
   end do

   if (ip==-1) then ! no table exists for these parameters. Initialize a new one.
     !
     ! Backward plans
     if ( bw_plan(1,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(1,bw_icurrent) )
     if ( bw_plan(2,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(2,bw_icurrent) )
     if ( bw_plan(3,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(3,bw_icurrent) )
     !
     ! The prototype fo dfftw_plan_many_dft is:
     ! dfftw_plan_many_dft(rank, n, howmany, 
     !   fin,  iembed, istride, idist, 
     !   fout, oembed, ostride, odist, isign, my_flags)
     !
     bw_plan(1, bw_icurrent) = zplan_many_dft(1, (/nx/), 1,       &
&         ff, (/ldx, ldy, ldz/), 1, ldx,                          &
&         ff, (/ldx, ldy, ldz/), 1, ldx, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_plan(2, bw_icurrent) = zplan_many_dft(1, (/ny/), nx,      &
&         ff, (/ldx, ldy, ldz/), ldx, 1,                          &
&         ff, (/ldx, ldy, ldz/), ldx, 1, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_plan(3, bw_icurrent) = zplan_many_dft(1, (/nz/), ldx*ldy, & ! We have to visit the entire augmented x-y plane!
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1,                      &
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_dims(:,bw_icurrent) = (/nx,ny,nz,ldx,ldy,ldz/)
     ip = bw_icurrent ! Got new plan.
     bw_icurrent = MOD( bw_icurrent, MPLANES ) + 1
   end if
   !
   ! 1) Transform along x.
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       do jj=1,ny
         if (do_fft_x(jj,kk) == 1) then 
           sidx = 1+ (jj-1)*ldx + (kk-1)*ldx*ldy ! Pass the pointer.
           call dfftw_execute_dft(bw_plan(1,ip), ff(sidx), ff(sidx) )
         end if
       end do
     end if
   end do
   !
   ! 2) Transform along y.
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       sidx = 1 + ldx*ldy*(kk-1)
       call dfftw_execute_dft( bw_plan(2,ip), ff(sidx), ff(sidx) )
     end if
   end do
   !
   ! 3) Transform along z.
   call dfftw_execute_dft( bw_plan(3,ip), ff, ff)
 
 CASE (FFTW_FORWARD) ! R --> G
   !
   ! First check if there is already a table initialized for this combination of parameters
   ip = -1
   do ii=1,MPLANES
     if( ALL( (/nx,ny,nz,ldx,ldy,ldz/) == fw_dims(:,ii) )) then
       ip=ii; EXIT
     end if
   end do

   if (ip==-1) then ! no table exists for these parameters. Initialize a new one.
     ! The prototype for dfftw_plan_many_dft is:
     ! dfftw_plan_many_dft(n, howmany, 
     !   fin,  iembed, istride, idist, 
     !   fout, oembed, ostride, odist, isign, my_flags)
     !
     ! Forward plans
     if ( fw_plan(1,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(1,fw_icurrent) )
     if ( fw_plan(2,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(2,fw_icurrent) )
     if ( fw_plan(3,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(3,fw_icurrent) )

     fw_plan(1, fw_icurrent) = zplan_many_dft(1, (/nx/), 1,       &
&         ff, (/ldx, ldy, ldz/), 1, ldx,                          &
&         ff, (/ldx, ldy, ldz/), 1, ldx, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_plan(2, fw_icurrent) = zplan_many_dft(1, (/ny/), nx,      &
&         ff, (/ldx, ldy, ldz/), ldx, 1,                          &
&         ff, (/ldx, ldy, ldz/), ldx, 1, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_plan(3, fw_icurrent) = zplan_many_dft(1, (/nz/), ldx*ldy, & ! We have to visit the entire augmented x-y plane!
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1,                      &
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_dims(:,fw_icurrent) = (/nx,ny,nz,ldx,ldy,ldz/)
     ip = fw_icurrent ! Got new plan.
     fw_icurrent = MOD( fw_icurrent, MPLANES ) + 1
   end if
   !
   ! 1) Transform along z.
   call dfftw_execute_dft( fw_plan(3,ip), ff, ff)
   !
   ! 2) Transform along y.
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       sidx = 1 + ldx*ldy*(kk-1)
       call dfftw_execute_dft( fw_plan(2,ip), ff(sidx), ff(sidx) )
     end if
   end do
   !
   ! 3) Transform along x. 
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       do jj=1,ny
         if (do_fft_x(jj,kk) == 1) then 
           sidx = 1 + (jj-1)*ldx + (kk-1)*ldx*ldy ! Pass the pointer.
           call dfftw_execute_dft( fw_plan(1,ip), ff(sidx), ff(sidx) )
         end if
       end do
     end if
   end do
   !
   ! 4) Normalize the transform.
   call ZDSCAL(ldx*ldy*nz, one/(nx*ny*nz), ff, 1)
 
 CASE DEFAULT 
   MSG_BUG("Wrong isign")
 END SELECT

 ABI_DEALLOCATE(do_fft_x)
 ABI_DEALLOCATE(do_fft_y)

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 !if (PRESENT(fftw_flags)) then
 !  ABI_UNUSED(fftw_flags)
 !end if
 ABI_UNUSED(ff(1))
#endif

end subroutine fftw3_fftpad_cplx
!!***

#ifdef HAVE_FFT_FFTW3 

!----------------------------------------------------------------------

!!****f* m_fftw3/dplan_many_dft
!! NAME
!!
!! FUNCTION
!!  
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function dplan_many_dft(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,sign,flags) result(plan)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dplan_many_dft'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride, sign,flags,idist,odist
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 real(dp) :: fin(*),fout(*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************
 
 call dfftw_plan_many_dft(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, sign, flags)

 if (plan==NULL_PLAN) then
   call wrtout(std_out,"dfftw_plan_many_dft returned NULL_PLAN!","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),3(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n= ",n," howmany= ",howmany," sign= ",sign," flags= ",flags,ch10,&
&    " inembed= ",inembed," istride= ",istride," idist=",idist,ch10,    &
&    " onembed= ",onembed," ostride= ",ostride," odist=",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function dplan_many_dft
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/zplan_many_dft
!! NAME
!!
!! FUNCTION
!!  
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!! FIXME  technically it should be intent(inout) since FFTW3 can destroy the input 
!! for particular flags.

function zplan_many_dft(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,sign,flags) result(plan)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'zplan_many_dft'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride, sign,flags,idist,odist
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 complex(dpc) :: fin(*),fout(*) 

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************
 
 call dfftw_plan_many_dft(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, sign, flags)

 if (plan==NULL_PLAN) then ! handle the error
   call wrtout(std_out,"dfftw_plan_many_dft returned NULL_PLAN (complex version)","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),3(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n= ",n," howmany= ",howmany," sign= ",sign," flags= ",flags,ch10,&
&    " inembed= ",inembed," istride= ",istride," idist=",idist,ch10,    &
&    " onembed= ",onembed," ostride= ",ostride," odist=",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function zplan_many_dft
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/dplan_many_dft_r2c
!! NAME
!!
!! FUNCTION
!!  
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!! FIXME  technically it should be intent(inout) since FFTW3 can destroy the input 
!! for particular flags.

function dplan_many_dft_r2c(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,flags) result(plan)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dplan_many_dft_r2c'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride,flags,idist,odist
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 real(dp) :: fin(*),fout(*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************

 call dfftw_plan_many_dft_r2c(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, flags)

 if (plan==NULL_PLAN) then ! handle the error.
   call wrtout(std_out,"dfftw_plan_many_dft_r2c returned NULL_PLAN","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),2(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n= ",n," howmany= ",howmany," flags= ",flags,ch10,&
&    " inembed= ",inembed," istride= ",istride," idist=",idist,ch10,    &
&    " onembed= ",onembed," ostride= ",ostride," odist=",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function dplan_many_dft_r2c
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/dplan_many_dft_c2r
!! NAME
!!
!! FUNCTION
!!  
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function dplan_many_dft_c2r(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,flags) result(plan)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dplan_many_dft_c2r'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride,flags,idist,odist
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 real(dp) :: fin(*),fout(*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************

 call dfftw_plan_many_dft_c2r(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, flags)

 if (plan==NULL_PLAN) then ! handle the error.
   call wrtout(std_out,"dfftw_plan_many_dft_c2r returned NULL_PLAN","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),2(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n= ",n," howmany= ",howmany," flags= ",flags,ch10,&
&    " inembed= ",inembed," istride= ",istride," idist=",idist,ch10,    &
&    " onembed= ",onembed," ostride= ",ostride," odist=",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function dplan_many_dft_c2r
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftpad_tr
!! NAME
!!  fftw3_fftpad_tr
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with FFTW3 taking advantage
!!  of time reversal symmetry.
!!  The 3D ffts are computed only on lines and planes which have non zero elements. 
!!  These lines and planes are defined by the two vectors do_fft_x(ldy*nz) and do_fft_y(nz) 
!!  FFT transform is in-place.
!!  
!! INPUTS
!!   logical dimensions of the fft physical dimensions of the f array sign of the transformation
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!   f(2*ldx*ldy*ldz)=
!!     input: The array with the data to be transformed.
!!     output: The results of the FFT.
!!
!! PARENTS
!!      fftw3_fourwf
!!
!! CHILDREN
!!      dfftw_destroy_plan,dfftw_execute_dft,dfftw_execute_dft_c2r,zdscal
!!
!! SOURCE

subroutine fftw3_fftpad_tr(ff,nx,ny,nz,ldx,ldy,ldz,mgfft,isign,gbound)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_fftpad_tr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 real(dp),intent(inout) :: ff(2*ldx*ldy*ldz)

#ifdef HAVE_FFT_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank1=1
 integer :: ndat,istride,ostride,idist,odist
 integer :: kk,idir,ip,ii,jj,sidx,nhp,padz,igf,igp,ifft
 integer :: id_yz,g3_max,g3_min,len3,gg2,gg3,g2,ifft_g3,igb,g2min,g2max
 integer,save :: bw_icurrent = 1
 integer,save :: fw_icurrent = 1
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer,save :: fw_dims(6,MPLANES) = -1
 integer,save :: bw_dims(6,MPLANES) = -1
 integer(KIND_FFTW_PLAN),save :: fw_plan (3,MPLANES) = 0
 integer(KIND_FFTW_PLAN),save :: bw_plan (3,MPLANES) = 0
 integer :: n(1),inembed(1),onembed(1)
 real(dp),allocatable :: gg(:) !ldx*ldy*ldz) ! logical or physical dims?
 real(dp),allocatable :: ff_hp(:)
 integer,allocatable :: do_fft_x(:,:),do_fft_y(:)
 type(fftw3_plan3_t) ::  Saved_plans(MPLANES)

! *************************************************************************

 !TODO
 MSG_ERROR("FFTW3 with istwf_k>1 is still under development!")

 ABI_ALLOCATE(do_fft_x,(ldy,ldz))
 ABI_ALLOCATE(do_fft_y,(nz))
 do_fft_y = 0                             ! we have to recalculate them at each call.
 do_fft_x = 0

 g3_min=gbound(3,2)
 g3_max=gbound(4,2)
 len3=g3_max-g3_min+1

 do gg3=1,g3_max+1 ! Loop over the z-planes intersecting the G-sphere.

   ifft_g3=gg3
   do_fft_y(ifft_g3) = 1

   ! Select the set of y for this z-plane.
   igb=2*gg3+3
   g2min = gbound(igb  ,2) 
   g2max = gbound(igb+1,2)
   do_fft_x(1:g2max+1    ,ifft_g3) = 1 ! Positive g_y.
   do_fft_x(g2min+ny+1:ny,ifft_g3) = 1 ! Negative g_y.
 end do
 !do_fft_y=1; do_fft_x=1
 !write(std_out,*)"do_fft_x",do_fft_x
 !write(std_out,*)"do_fft_y",do_fft_y
 !stop

 SELECT CASE (isign)

 CASE (FFTW_BACKWARD) ! G --> R
   !
   ! First check if there is already a table initialized for this combination of parameters
   ip = -1
   do ii=1,MPLANES
     if ( ALL( (/nx,ny,nz,ldx,ldy,ldz/) == bw_dims(:,ii) )) then
       ip=ii; EXIT
     end if
   end do

   if (ip==-1) then ! no table exist for these parameters initialize a new one.
     !
     ! Backward plans
     if ( bw_plan(1,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(1,bw_icurrent) )
     if ( bw_plan(2,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(2,bw_icurrent) )
     if ( bw_plan(3,bw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( bw_plan(3,bw_icurrent) )
     !
     ! The prototype for dfftw_plan_many_dft is:
     ! dfftw_plan_many_dft(rank, n, howmany, 
     !   fin,  iembed, istride, idist, 
     !   fout, oembed, ostride, odist, isign, my_flags)
     !     
     bw_plan(1, bw_icurrent) = dplan_many_dft(1, (/nx/), 1,       &   ! Single 1D transform of f(Gx,Gy,Gz) along Gx.
&         ff, (/ldx, ldy, ldz/), 1, ldx,                          &
&         ff, (/ldx, ldy, ldz/), 1, ldx, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_plan(2, bw_icurrent) = dplan_many_dft(1, (/ny/), nx,      &   ! nx 1D transforms of f(x,Gy,Gz) along Gy.
&         ff, (/ldx, ldy, ldz/), ldx, 1,                          &
&         ff, (/ldx, ldy, ldz/), ldx, 1, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_plan(3, bw_icurrent) = dplan_many_dft(1, (/nz/), ldx*ldy, & ! ldx*ldy 1D transforms of f(x,y,Gz) along Gz.
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1,                      & ! Note that we have to visit the entire augmented x-y plane!
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1, FFTW_BACKWARD, FFTW_ESTIMATE)

     bw_dims(:,bw_icurrent) = (/nx,ny,nz,ldx,ldy,ldz/)
     ip = bw_icurrent ! Got new plan.
     bw_icurrent = MOD( bw_icurrent, MPLANES ) + 1
   end if
   !
   ! 1) Transform along x.
   !do kk=1,nz
   do kk=1,g3_max+1
     if (do_fft_y(kk) == 1) then
       do jj=1,ny
         if (do_fft_x(jj,kk) == 1) then 
           sidx = 1 + 2*(jj-1)*ldx + 2*(kk-1)*ldx*ldy ! Pass the pointer, 2 to account for the imag part.
           call dfftw_execute_dft(bw_plan(1,ip), ff(sidx), ff(sidx) )
         end if
       end do
     end if
   end do
   !
   ! 2) Transform along y.
   !do kk=1,nz
   do kk=1,g3_max+1
     if (do_fft_y(kk) == 1) then
       sidx = 1 + 2*(kk-1)*ldx*ldy
       call dfftw_execute_dft( bw_plan(2,ip), ff(sidx), ff(sidx) )
     end if
   end do
   !
   ! Now ff contains f(x,y,Gz) for Gz >= 0
   ! Note that f(x,y,-Gz) = f(z,y,Gz)* hence here we can perform the c2r FFT.
   !RETURN

   ! Fill the Hermitian part: Hermitian redundancy: out[i] is the conjugate of out[n-i]
   padz = (nz/2+1)
   !nhp  = padz*nx*ny
   !allocate(ff_hp(2*nhp))
   !do kk=1,nz/2+1
   !  do jj=1,ny
   !    igf = 2 * ( (kk-1)*ldx*ldy + (jj-1)*nx )
   !    igp = 2 * ( (kk-1)*padz*ny + (jj-1)*padz) 
   !    ff_hp(:,igp+1:igp+padz) = ff(:,igf+1:igf+padz)
   !  end do
   !end do
   !deallocate(ff_hp)

   n = (/nz/)
   !n = (/padz/)
   ndat = ldx*ldy
   inembed  = (/ldx*ldy*ldz/)
   !istride = 1
   !idist   = nx*ny
   istride = ldx*ldy
   idist   = 1

   onembed  = (/ldx*ldy*ldz/)
   ostride = ldx*ldy
   odist   = 1

   ABI_ALLOCATE(gg,(ldx*ldy*ldz))

   my_plan = dplan_many_dft_c2r(rank1, n, ndat, ff, inembed, istride, idist, gg, onembed, ostride, odist, FFTW_ESTIMATE)

   ! Now perform the 3D FFT via FFTW. c2r are always FFTW_BACKWARD
   call dfftw_execute_dft_c2r( my_plan, ff, gg)

   do ifft=1,ldx*ldy*ldz
     ff(2*ifft-1) = gg(ifft)
     ff(2*ifft) = zero
   end do
   ABI_DEALLOCATE(gg)
 
 CASE (FFTW_FORWARD) ! R --> G

   MSG_ERROR("R-->G not coded")
   !
   ! First check if there is already a table initialized for this combination of parameters
   ip = -1
   do ii=1,MPLANES
     if( ALL( (/nx,ny,nz,ldx,ldy,ldz/) == fw_dims(:,ii) )) then
       ip=ii; EXIT
     end if
   end do

   if (ip==-1) then ! no table exist for these parameters initialize a new one
     ! The prototype for dfftw_plan_many_dft is:
     ! dfftw_plan_many_dft(rank, n, howmany, 
     !   fin,  iembed, istride, idist, 
     !   fout, oembed, ostride, odist, isign, my_flags)
     !
     ! Forward plans
     if ( fw_plan(1,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(1,fw_icurrent) )
     if ( fw_plan(2,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(2,fw_icurrent) )
     if ( fw_plan(3,fw_icurrent) /= NULL_PLAN ) call dfftw_destroy_plan( fw_plan(3,fw_icurrent) )

     fw_plan(1, fw_icurrent) = dplan_many_dft(1, (/nx/), 1,      &
&         ff, (/ldx, ldy, ldz/), 1, ldx,                         &
&         ff, (/ldx, ldy, ldz/), 1, ldx, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_plan(2, fw_icurrent) = dplan_many_dft(1, (/ny/), nx,     &
&         ff, (/ldx, ldy, ldz/), ldx, 1,                         &
&         ff, (/ldx, ldy, ldz/), ldx, 1, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_plan(3, fw_icurrent) = dplan_many_dft(1, (/nz/), ldx*ldy,& ! We have to visit the entire augmented x-y plane!
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1,                     &
&         ff, (/ldx, ldy, ldz/), ldx*ldy, 1, FFTW_FORWARD, FFTW_ESTIMATE)

     fw_dims(:,fw_icurrent) = (/nx,ny,nz,ldx,ldy,ldz/)
     ip = fw_icurrent ! Got new plan.
     fw_icurrent = MOD( fw_icurrent, MPLANES ) + 1
   end if
   !
   ! 1) Transform along z.
   call dfftw_execute_dft( fw_plan(3,ip), ff, ff)
   !
   ! 2) Transform along y.
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       sidx = 1 + 2*ldx*ldy*(kk-1)
       call dfftw_execute_dft( fw_plan(2,ip), ff(sidx), ff(sidx) )
     end if
   end do
   !
   ! 3) Transform along x. 
   do kk=1,nz
     if (do_fft_y(kk) == 1) then
       do jj=1,ny
         if (do_fft_x(jj,kk) == 1) then 
           sidx = 1 + 2*(jj-1)*ldx + 2*(kk-1)*ldx*ldy ! Pass the pointer, 2 to account for the imag part.
           call dfftw_execute_dft( fw_plan(1,ip), ff(sidx), ff(sidx) )
         end if
       end do
     end if
   end do
   !
   ! 4) Normalize the transform.
   !call DSCAL(2*ldx*ldy*nz, one/(nx*ny*nz), ff, 1)
   call ZDSCAL(ldx*ldy*nz, one/(nx*ny*nz), ff, 1)
 
 CASE DEFAULT 
   MSG_BUG("Wrong isign")
 END SELECT

 ABI_DEALLOCATE(do_fft_x)
 ABI_DEALLOCATE(do_fft_y)

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 ABI_UNUSED(ff(1))
#endif

end subroutine fftw3_fftpad_tr     
!!***

!----------------------------------------------------------------------


END MODULE m_fftw3
!!***
