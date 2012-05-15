!!****f* ABINIT/fftw3_fourdp
!! NAME
!!  fftw3_fourdp
!!
!! FUNCTION
!! Driver routine for 3D FFT of lengths nx, ny, nz. Mainly used for densities or potentials.
!! FFT Transform is out-of-place
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nx,ny,nz=Number of point along the three directions.
!! isign= +1 : fofg(G) => fofr(R); 
!!        -1 : fofr(R) => fofg(G)
!! fofg(2,nx*ny*nz)=The array to be transformed.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator. Defaults to FFTW_ESTIMATE.
!!
!! OUTPUT 
!! fofr(cplex,nx*ny*nz)=The FFT of fofg
!!
!! NOTES
!!  1) No augmentation of the mesh to reduce memory conflicts, as it is not not supported by FFTW3.
!!
!! PARENTS
!!      fourdp
!!
!! CHILDREN
!!      fftw3_c2c_op,fourdp,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fftw3_fourdp(cplex,nx,ny,nz,isign,fofg,fofr,fftw_flags)

 use m_profiling

 use defs_basis
 use m_fftw3
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_fourdp'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nx,ny,nz,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(inout) :: fofg(2,nx*ny*nz)
 real(dp),intent(inout) :: fofr(cplex*nx*ny*nz)

!Local variables-------------------------------
!scalars
 integer :: my_flags,fidx,ii,jj,kk
 integer :: ldx,ldy,ldz
 logical :: have_cache_conflicts
!arrays
 real(dp),allocatable :: work1(:,:,:,:) !,work2(:,:)

! *************************************************************************

 my_flags=FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 ldx=nx; ldy= ny; ldz=nz
 have_cache_conflicts = ( 2*(nx/2)==nx .or. 2*(ny/2)==ny) 
 if (have_cache_conflicts) then
   ldx = 2*(nx/2)+1
   ldy = 2*(ny/2)+1
   ldz = nz
 end if
 have_cache_conflicts = .FALSE.

!TODO At present, no augmentation to avoid cache conflicts. 

 select case (cplex)
   case (2) ! Complex to Complex.

     select case (isign)

       case (FFTW_BACKWARD)  ! +1

         if (have_cache_conflicts) then ! Have to augment the array.
!          write(std_out,*)"have cache conflicts"

           ABI_ALLOCATE(work1,(2,ldx,ldy,ldz))
!          allocate(work2(2,ldx*ldy*ldz))
           fidx=0
           do kk=1,nz
             do jj=1,ny
               fidx = (jj-1)*nx + (kk-1)*nx*ny
               do ii=1,nx
!                widx = ii + (jj-1)*ldx + (kk-1)*ldx*ldy
!                fidx = fidx+1
!                fidx = ii + (jj-1)*nx + (kk-1)*nx*ny
                 work1(:,ii,jj,kk) = fofg(:,ii+fidx)
               end do
             end do
           end do

           call fftw3_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,1,isign,work1,fftw_flags=my_flags)
!          call fftw3_many_dft_op(nx,ny,nz,ldx,ldy,ldz,1,isign,work1,work2,fftw_flags=my_flags)

           fidx=0
           do kk=1,nz
             do jj=1,ny
               do ii=1,nx
!                widx = ii + (jj-1)*ldx + (kk-1)*ldx*ldy
!                fidx = ii + (jj-1)*nx + (kk-1)*nx*ny
                 fidx = fidx+1
                 fofr(2*fidx-1) = work1(1,ii,jj,kk)
                 fofr(2*fidx  ) = work1(2,ii,jj,kk)
               end do
             end do
           end do
           ABI_DEALLOCATE(work1)
!          deallocate(work2)

         else 
           call fftw3_many_dft_op(nx,ny,nz,nx,ny,nz,1,isign,fofg,fofr,fftw_flags=my_flags)
         end if

       case (FFTW_FORWARD) ! -1
         call fftw3_many_dft_op(nx,ny,nz,nx,ny,nz,1,isign,fofr,fofg,fftw_flags=my_flags)

         case default
         MSG_BUG("Wrong isign")
     end select 

   case (1) 
!    Real case.
     select case (isign)
       case (FFTW_FORWARD) ! -1; R --> G 
         call fftw3_r2c_op(nx,ny,nz,nx,ny,nz,1,fofr,fofg,fftw_flags=my_flags)
       case (FFTW_BACKWARD) ! +1; G --> R
         call fftw3_c2r_op(nx,ny,nz,nx,ny,nz,1,fofg,fofr,fftw_flags=my_flags)
         case default
         MSG_BUG("Wrong isign")
     end select

     case default 
     MSG_BUG(" Wrong value for cplex")
 end select 

end subroutine fftw3_fourdp
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/fourdp_c2c_ip
!! NAME
!!  fourdp_c2c_ip
!!
!! FUNCTION
!!  In-place FFT transform of a complex array in double precision.
!!  It calls FFTW3 if available, fallback to fourdp if FFTW3 support is missing. 
!!
!! INPUTS
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r, 
!!   -1 for transforming from r to G.
!! nfft=(effective) number of FFT grid points (for this processor)
!! paral_kgb=Flag related to the kpoint-band-fft parallelism
!! tim_fourdp=timing code of the calling routine (can be set to 0 if not attributed)
!! mpi_enreg=information about MPI parallelization
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! SIDE EFFECTS
!!  ff(nfft)= Changed in output, filled with the FFT results.
!!
!! NOTES
!!  The FFTW3 version does not support paralkgb==1.
!!
!! PARENTS
!!      calc_sig_ppm_eet,debug_tools,m_fft_prof,m_oscillators,m_paw_pwij
!!      m_shirley
!!
!! CHILDREN
!!      fftw3_c2c_op,fourdp,timab
!!
!! SOURCE

subroutine fourdp_c2c_ip(ff,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

 use m_profiling

 use defs_basis
 use m_fftw3
 use m_errors

 use defs_abitypes,  only : MPI_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourdp_c2c_ip'
 use interfaces_18_timing
 use interfaces_53_ffts, except_this_one => fourdp_c2c_ip
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isign,nfft,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 complex(dpc),intent(inout) :: ff(nfft)

!Local variables ------------------------------
!scalars
 integer,parameter :: cplex=2
 integer :: ifft,idx
 integer :: n1,n2,n3,ldx,ldy,ldz
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: fofg(:,:),fofr(:)

! *************************************************************************

 if (ngfft(7) == 312) then !Call the FFTW wrapper

   call timab(260+tim_fourdp,1,tsec) !Keep track of timing

   if (paral_kgb/=0) MSG_ERROR("FFT + paral_kgb/=0 not coded")

   n1 =ngfft(1); n2=ngfft(2); n3=ngfft(3)

!  if (ANY(ngfft(1:3)/=ngfft(4:6)) then 
!  ldx=ngfft(4); ldy=ngfft(5); ldz=ngfft(6)
!  call fftw3_c2c_op(nx,ny,nz,ldx,ldy,ldz,1,isign,ff,gg)
!  else

   ldx=n1; ldy=n2; ldz=n3 ! No augmentation here.
   call fftw3_c2c_ip(n1,n2,n3,ldx,ldy,ldz,1,isign,ff)

   call timab(260+tim_fourdp,2,tsec)

   RETURN

 else ! Fallback to fourdp.

   ABI_ALLOCATE(fofg,(2,nfft))
   ABI_ALLOCATE(fofr,(cplex*nfft))

   if (isign==-1) then ! R->G
     do ifft=1,nfft
       idx = 2*ifft-1
       fofr(idx  ) = DBLE (ff(ifft))
       fofr(idx+1) = AIMAG(ff(ifft))
     end do
   else  ! G->R 
     fofg(1,:) = DBLE (ff)
     fofg(2,:) = AIMAG(ff)
   end if

   call fourdp(cplex,fofg,fofr,isign,MPI_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

!  TODO Call Goedecker FFT for complex2complex transform.
!  ris=DBLE(isign)
!  call sg_fft(fftcache,ldx,ldy,ldz,n1,n2,n3,work1,work2,ris)

!  Copy back the values.
   if (isign==-1) then ! ff contains ff(G)
     ff(:) = DCMPLX(fofg(1,:),fofg(2,:))
   else ! ff contains ff(r)
     do ifft=1,nfft
       idx = 2*ifft-1
       ff(ifft) = DCMPLX(fofr(idx),fofr(idx+1))
     end do
   end if

   ABI_DEALLOCATE(fofg)
   ABI_DEALLOCATE(fofr)

   RETURN
 end if

end subroutine fourdp_c2c_ip
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/fourdp_c2c_op
!! NAME
!!  fourdp_c2c_op
!!
!! FUNCTION
!!  Out-of-place FFT transform of a complex array in double precision.
!!  It calls FFTW3 if available, fallback to fourdp if FFTW3 support is missing. 
!!
!! INPUTS
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r, 
!!   -1 for transforming from r to G.
!! nfft=(effective) number of FFT grid points (for this processor)
!! paral_kgb=Flag related to the kpoint-band-fft parallelism
!! tim_fourdp=timing code of the calling routine (can be set to 0 if not attributed)
!! mpi_enreg=information about MPI parallelization
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! ff(nfft)=The input array to be transformed.
!!
!! SIDE EFFECTS
!! OUTPUT 
!!  gg(nfft)= The FFT results.
!!
!! NOTES
!!  The FFTW3 version does not support paral_kgb==1.
!!
!! PARENTS
!!      calc_sig_ppm_eet,m_fft_prof,m_oscillators
!!
!! CHILDREN
!!      fftw3_c2c_op,fourdp,timab
!!
!! SOURCE

subroutine fourdp_c2c_op(ff,gg,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

 use m_profiling

 use defs_basis
 use m_fftw3
 use m_errors

 use defs_abitypes,  only : MPI_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourdp_c2c_op'
 use interfaces_18_timing
 use interfaces_53_ffts, except_this_one => fourdp_c2c_op
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isign,nfft,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 complex(dpc),intent(in) :: ff(nfft)
 complex(dpc),intent(out) :: gg(nfft)

!Local variables ------------------------------
!scalars
 integer,parameter :: cplex=2
 integer :: ifft,idx
 integer :: n1,n2,n3,ldx,ldy,ldz
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: fofg(:,:),fofr(:)

! *************************************************************************

 if (ngfft(7)==312) then ! call FFTW3 wrapper

   if (paral_kgb/=0) then 
     MSG_ERROR("FFT + paral_kgb/=0 not coded")
   end if

   call timab(260+tim_fourdp,1,tsec) !Keep track of timing

   n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
   ldx=n1; ldy=n2; ldz=n3 ! no augmentation
!  ldx=ngfft(4); ldy=ngfft(5); ldz=ngfft(6)

   call fftw3_c2c_op(n1,n2,n3,ldx,ldy,ldz,1,isign,ff,gg)

   call timab(260+tim_fourdp,2,tsec)

 else ! Fallback to fourdp.
   ABI_ALLOCATE(fofg,(2,nfft))
   ABI_ALLOCATE(fofr,(cplex*nfft))

   if (isign==-1) then ! R->G
     do ifft=1,nfft
       idx = 2*ifft-1
       fofr(idx  ) = DBLE (ff(ifft))
       fofr(idx+1) = AIMAG(ff(ifft))
     end do
   else  ! G->R 
     fofg(1,:) = DBLE (ff)
     fofg(2,:) = AIMAG(ff)
   end if

   call fourdp(cplex,fofg,fofr,isign,MPI_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

!  Copy back the values.
   if (isign==-1) then ! ff contains ff(G)
     gg(:) = DCMPLX(fofg(1,:),fofg(2,:))
   else ! ff contains ff(r)
     do ifft=1,nfft
       idx = 2*ifft-1
       gg(ifft) = DCMPLX(fofr(idx),fofr(idx+1))
     end do
   end if

   ABI_DEALLOCATE(fofg)
   ABI_DEALLOCATE(fofr)
 end if

end subroutine fourdp_c2c_op
!!***

!----------------------------------------------------------------------
