!{\src2tex{textfont=tt}}
!!****f* ABINIT/fftw3_fourwf
!! NAME
!! fftw3_fourwf
!!
!! FUNCTION
!! Carry out composite Fourier transforms between real and reciprocal (G) space.
!! Wavefunctions, contained in a sphere in reciprocal space,
!! can be FFT to real space. They can also be FFT from real space
!! to a sphere. Also, the density maybe accumulated, and a local potential can be applied.
!!
!! The different options are :
!! - option=0 --> reciprocal to real space and output the result.
!! - option=1 --> reciprocal to real space and accumulate the density.
!! - option=2 --> reciprocal to real space, apply the local potential to the wavefunction
!!                in real space and produce the result in reciprocal space.
!! - option=3 --> real space to reciprocal space.
!!                NOTE that in this case, fftalg=1x1 MUST be used. This may be changed in the future.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex= if 1 , denpot is real, if 2 , denpot is complex
!!    (cplex=2 only allowed for option=2, and istwf_k=1)
!!    not relevant if option=0 or option=3, so cplex=0 can be used to minimize memory
!! fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag)
!! gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!! gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!! istwf_k=option parameter that describes the storage of wfs
!! kg_kin(3,npwin)=reduced planewave coordinates, input
!! kg_kout(3,npwout)=reduced planewave coordinates, output
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=information about MPI parallelization
!! ndat=number of FFT to do in //
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npwin=number of elements in fofgin array (for option 0, 1 and 2)
!! npwout=number of elements in fofgout array (for option 2 and 3)
!! n4,n5,n6=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
!! option= if 0: do direct FFT
!!         if 1: do direct FFT, then sum the density
!!         if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!         if 3: do reverse FFT only
!! paral_kgb=Flag related to the kpoint-band-fft parallelism
!! weight_r=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1)
!!
!! weight_i=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1 and (fftalg=4 and fftalgc/=0))
!! fofginb(2,npwin)=holds second input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag) 
!!                 (for non diagonal occupation)
!! use_ndo = use non diagonal occupations.

!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                fofr(2,n4,n5,n6) contains the output Fourier Transform of fofgin;
!!                no use of denpot, fofgout and npwout.
!! for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input density at input,
!!                and the updated density at output (accumulated);
!!                no use of fofgout and npwout.
!! for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input local potential;
!!                fofgout(2,npwout*ndat) contains the output function;
!! for option==3, fofr(2,n4,n5,n6*ndat) contains the input real space wavefunction;
!!                fofgout(2,npwout*ndat) contains its output Fourier transform;
!!                no use of fofgin and npwin.
!!
!! PARENTS
!!      fourwf
!!
!! CHILDREN
!!      fftw3_fftpad_cplx,sg_fftpad
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fftw3_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&  kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,paral_kgb,weight_r,weight_i,&
&  use_ndo,fofginb) ! optional Arguments.

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_fftw3
 use m_errors
 use m_timer

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw3_fourwf'
 use interfaces_53_ffts, except_this_one => fftw3_fourwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,n4,n5,n6,ndat,npwin,npwout,option,paral_kgb,mgfft
 integer,intent(in),optional :: use_ndo
 real(dp),intent(in) :: weight_r,weight_i
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(inout) :: denpot(cplex*n4,n5,n6),fofgin(2,npwin*ndat)
! real(dp) :: fofginb(2,npwin*ndat)
 real(dp),intent(inout),optional :: fofginb(:,:)
 real(dp),intent(inout) :: fofr(2,n4,n5,n6*ndat)
 real(dp),intent(out) :: fofgout(2,npwout*ndat)

!Local variables-------------------------------
!scalars
 integer :: fftalg,fftalga,fftalgc,i1,i2,i3,idat
 integer :: ig,padat,me_fft,n1,n2,n3,fftcache
 integer :: nfftot,nproc_fft !,idx
 real(dp) :: fim,fre
 character(len=500) :: msg
 logical :: luse_ndo
!arrays
 integer :: shiftg(3),symm(3,3)
!no_abirules
!real(dp),allocatable :: tmp1d(:)

! *************************************************************************

 if (ndat/=1) then
   MSG_ERROR("ndat/=1 not coded yet")
 end if

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3); nfftot=n1*n2*n3
 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=mod(fftalg,10)
 me_fft=ngfft(11); nproc_fft=ngfft(10)
 fftcache=ngfft(8)

 luse_ndo=.false.
 if(present(use_ndo).and.present(fofginb)) then
   if(use_ndo==1) luse_ndo=.true.
 end if
 if(luse_ndo) then
   if((size(fofginb,2)==0).and.(use_ndo==1)) then
     write(msg, '(a,a,a,i4,i5)' )&
&     '  fofginb has a dimension equal to zero and use_ndo==1',ch10,&
&     '  Action : check dimension of fofginb',size(fofginb,2),use_ndo
     MSG_ERROR(msg)
   end if
 end if

 if (luse_ndo) then
   MSG_ERROR(" Non diagonal occupations are not compatible with FFTW3")
 end if

 if ( ALL(option /= (/0,1,2,3/)) ) then
   write(msg, '(a,i0,a)' )&
&   '  The option number',option,' is not allowed. Only option=0, 1, 2 or 3 are allowed presently.'
   MSG_ERROR(msg)
 end if

 if( option==1 .and. cplex/=1 )then
   write(msg,'(a,i0)')&
&   ' With the option number 1, cplex must be 1 but it is cplex=',cplex
   MSG_ERROR(msg)
 end if

 if( option==2 .and. (cplex/=1 .and. cplex/=2) )then
   write(msg,'(a,i0)')&
&   ' With the option number 2, cplex must be 1 or 2, but it is cplex=',cplex
   MSG_ERROR(msg)
 end if

 if (paral_kgb/=0) then 
   write(msg,'(a,i0)')' paral_kgb/=0 not yet compatible with FFTW3 wrappers but paral_kgb=',paral_kgb
   MSG_ERROR(msg)
 end if

 shiftg(:)=0
 symm(:,:)=0; symm(1,1)=1; symm(2,2)=1; symm(3,3)=1

 if (ANY( option == (/0,1,2/) ))  then
   call sphere(fofgin,ndat,npwin,fofr,n1,n2,n3,n4,n5,n6,kg_kin,istwf_k,1,mpi_enreg,shiftg,symm,one)

!  if (istwf_k/=2) then
   if (.TRUE.) then
     call fftw3_fftpad(fofr,n1,n2,n3,n4,n5,n6,mgfft,FFTW_BACKWARD,gboundin)
!    call fftw3_many_dft_ip(n1,n2,n3,n4,n5,n6,ndat,FFTW_BACKWARD,fofr) 
   else 
     call fftw3_fftpad_tr(fofr,n1,n2,n3,n4,n5,n6,mgfft,FFTW_BACKWARD,gboundin)
   end if

   select case (option) 
     case (1)
!      $omp parallel do private(i1,i2,i3,idat,padat) shared(denpot,fofr,n1,n2,n3,n6,ndat,weight_r)
       do idat=1,ndat
         padat = (idat-1)*n6
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               denpot(i1,i2,i3)=denpot(i1,i2,i3)+weight_r*(fofr(1,i1,i2,i3+padat)**2+fofr(2,i1,i2,i3+padat)**2)
             end do
           end do
         end do
       end do
       RETURN

     case (2)
       if(cplex==1)then
!        $omp parallel do private(i1,i2,i3,idat) shared(denpot,fofr,n1,n2,n3,ndat)
         do idat=1,ndat
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 fofr(1,i1,i2,i3+n3*(idat-1))=denpot(i1,i2,i3)*fofr(1,i1,i2,i3+n3*(idat-1))
                 fofr(2,i1,i2,i3+n3*(idat-1))=denpot(i1,i2,i3)*fofr(2,i1,i2,i3+n3*(idat-1))
               end do
             end do
           end do
         end do
       else if (cplex==2)then
!        $omp parallel do private(i1,i2,i3,idat,fre,fim) shared(denpot,fofr,n1,n2,n3,ndat)
         do idat=1,ndat
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 fre=fofr(1,i1,i2,i3+n3*(idat-1))
                 fim=fofr(2,i1,i2,i3+n3*(idat-1))
                 fofr(1,i1,i2,i3+n3*(idat-1))=denpot(2*i1-1,i2,i3)*fre-denpot(2*i1,i2,i3)*fim
                 fofr(2,i1,i2,i3+n3*(idat-1))=denpot(2*i1-1,i2,i3)*fim+denpot(2*i1,i2,i3)*fre
               end do
             end do
           end do
         end do
       end if ! cplex=2
!      
!      The data for option==2 is now in fofr.
       call fftw3_fftpad(fofr,n1,n2,n3,n4,n5,n6,mgfft,FFTW_FORWARD,gboundout)
!      call fftw3_many_dft_ip(n1,n2,n3,n4,n5,n6,ndat,FFTW_FORWARD,fofr) 
       
!      $omp parallel do private(i1,i2,i3,ig,idat) shared(npwout,n1,n2,n3,ndat,fofgout,fofr,kg_kout)
       do ig=1,npwout
         i1=kg_kout(1,ig); if (i1<0) i1=i1+n1; i1=i1+1
         i2=kg_kout(2,ig); if (i2<0) i2=i2+n2; i2=i2+1
         i3=kg_kout(3,ig); if (i3<0) i3=i3+n3; i3=i3+1
         do idat=1,ndat
           fofgout(1,ig+npwout*(idat-1))=fofr(1,i1,i2,i3+n3*(idat-1)) 
           fofgout(2,ig+npwout*(idat-1))=fofr(2,i1,i2,i3+n3*(idat-1)) 
         end do
       end do
   end select

 else if (option==3) then !  The data for option==3 is already in fofr.

#if 1
   call fftw3_fftpad(fofr,n1,n2,n3,n4,n5,n6,mgfft,FFTW_FORWARD,gboundout)

!  $omp parallel do private(i1,i2,i3,ig,idat) SHARED(n1,n2,n3,ndat,fofgout,fofr,kg_kout,npwout)
   do ig=1,npwout
     i1=kg_kout(1,ig); if (i1<0) i1=i1+n1; i1=i1+1
     i2=kg_kout(2,ig); if (i2<0) i2=i2+n2; i2=i2+1
     i3=kg_kout(3,ig); if (i3<0) i3=i3+n3; i3=i3+1
     do idat=1,ndat
       fofgout(1,ig+npwout*(idat-1))=fofr(1,i1,i2,i3+n3*(idat-1)) 
       fofgout(2,ig+npwout*(idat-1))=fofr(2,i1,i2,i3+n3*(idat-1)) 
     end do
   end do
#else
!  MSG_WARNING("Entering debugging section in FFTW3_FOURWF")
!  FIXME: 
!  Passing fofr to the FFTW3 wrappers produces wrong results when option==2
!  very strange since fofr in conformable with a vector and FFTW3 only needs 
!  the base address of the array. To work around the problem we have to change
!  the shape of the input before calling fftw3_fftpad although this leads to a general slow down of the run.

   ABI_ALLOCATE(tmp1d,(2*n4*n5*n6))
   tmp1d = RESHAPE(fofr,(/2*n4*n5*n6/))

   call fftw3_fftpad(tmp1d,n1,n2,n3,n4,n5,n6,mgfft,FFTW_FORWARD,gboundout)
!  call fftw3_many_dft_ip(n1,n2,n3,n4,n5,n6,ndat,FFTW_FORWARD,fofr) 

!  $omp parallel do private(i1,i2,i3,ig,idx,idat) shared(n1,n2,n3,n4,n5,ndat,fofgout,kg_kout,npwout,tmp1d)
   do ig=1,npwout
     i1=kg_kout(1,ig); if (i1<0) i1=i1+n1; i1=i1+1
     i2=kg_kout(2,ig); if (i2<0) i2=i2+n2; i2=i2+1
     i3=kg_kout(3,ig); if (i3<0) i3=i3+n3; i3=i3+1
     idx = i1 + (i2-1)*n4 + (i3-1)*n4*n5
     idx = 2*idx-1
     do idat=1,ndat
       fofgout(1,ig+npwout*(idat-1))=tmp1d(idx)
       fofgout(2,ig+npwout*(idat-1))=tmp1d(idx+1)
     end do
   end do

   ABI_DEALLOCATE(tmp1d)
#endif
 end if

 RETURN

 ABI_UNUSED(weight_i)

end subroutine fftw3_fourwf
!!***

!----------------------------------------------------------------------

!!****f* fftw3_fourwf/padded_fourwf_cplx
!! NAME
!!  padded_fourwf_cplx
!!
!! FUNCTION
!!  Driver routine used to transform COMPLEX wavefunctions using 3D zero-padded FFTs.
!!  
!! INPUTS
!!   ngfft(18)=Info on the 3D FFT.
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!   ff(2*ldx*ldy*ldz)=
!!     input: The array with the data to be transformed.
!!     output: The results of the FFT.
!!
!! PARENTS
!!      calc_sig_ppm_eet,m_fft_prof,m_oscillators,m_paw_pwij
!!
!! CHILDREN
!!      fftw3_fftpad_cplx,sg_fftpad
!!
!! SOURCE

subroutine padded_fourwf_cplx(ff,ngfft,nx,ny,nz,ldx,ldy,ldz,mgfft,isign,gbound)

 use m_profiling

 use defs_basis
 use m_errors

 use m_fftw3,   only : fftw3_fftpad_cplx

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'padded_fourwf_cplx'
 use interfaces_53_ffts, except_this_one => padded_fourwf_cplx
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,mgfft,isign
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: gbound(2*mgfft+8,2)
 complex(dpc),intent(inout) :: ff(ldx*ldy*ldz)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,idx,fftalg,fftalga,fftalgc
 !integer :: cplex,istwf_k,option
 real(dp) :: xnorm
 !real(dp) :: weight=one
 character(len=500) :: msg
!arrays
 !integer :: dum_kg_kout(0,0)
 !real(dp),allocatable ::  arr(:,:,:,:) !intent(inout) :
 real(dp),allocatable :: ftarr(:,:,:,:) !intent(out) :: 
 !real(dp) :: dum_denpot(0,0,0),dum_fofgin(0,0) 
 !real(dp),allocatable :: fofgin(:,:)
 real(dp),allocatable :: fofr(:,:,:,:)

! *************************************************************************
 
 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=MOD(fftalg,10)

 select case (fftalga)

   case (1) ! Goedecker"s routines.
!    msg = "Zero-padded FFT with fftalga==1 in GW part is still under development"
!    MSG_WARNING(msg)
!    
!    dummy code needed due to the explicit interface!
!    TODO: move sg_fftpad to 52_fft_no_abirules such that
!    it is possible to cheat the compiler!!!
     ABI_ALLOCATE(fofr,(2,ldx,ldy,ldz))
     fofr = zero
     do iz=1,nz
       do iy=1,ny
         do ix=1,nx
           idx = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
           fofr(1,ix,iy,iz) = DBLE (ff(idx))
           fofr(2,ix,iy,iz) = AIMAG(ff(idx))
         end do
       end do
     end do

     ABI_ALLOCATE(ftarr,(2,ldx,ldy,ldz))
     call sg_fftpad(ngfft(8),mgfft,ldx,ldy,ldz,nx,ny,nz,fofr,ftarr,DBLE(isign),gbound)
!    
!    Normalize and copy the results.
     xnorm=one; if (isign==-1) xnorm=one/DBLE(nx*ny*nz)
     do iz=1,nz
       do iy=1,ny
         do ix=1,nx
           idx = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
           ff(idx) = DCMPLX(ftarr(1,ix,iy,iz), ftarr(2,ix,iy,iz)) * xnorm
         end do
       end do
     end do
     ABI_DEALLOCATE(fofr)
     ABI_DEALLOCATE(ftarr)

   case (3) ! FFTW3 version.   
     call fftw3_fftpad_cplx(ff,nx,ny,nz,ldx,ldy,ldz,mgfft,isign,gbound)

!    case (4)
!    cplex=0; istwf_k=1; option=3
!    $   call sg_fftrisc(cplex,dum_denpot,fofgin,dum_fofgin,fofr,gbound,gbound,istwf_k,dum_gvec,gvec,&
!    $&    mgfft,ngfft,npwwfn,npwwfn,ldx,ldy,ldz,option,weight)

     case default
     write(msg,'(a,i0,a)')" fftalga=", fftalga," not coded "
     MSG_ERROR(msg)
 end select 

end subroutine padded_fourwf_cplx       
!!***

!----------------------------------------------------------------------
