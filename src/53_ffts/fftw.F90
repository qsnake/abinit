!{\src2tex{textfont=tt}}
!!****f* ABINIT/fftw
!! NAME
!! fftw
!!
!! FUNCTION
!! complex-to-complex FFT using FFTW
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (PT,HM)
!! Initial version 30 Nov 2000 by Pascal Thibaudeau and Herve Mathis
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  n1,n2,n3=physical dimension of the transform. Each of these must be a
!!   product of the prime factors 2,3,5. If two ni s are equal
!!   it is recommended to place them behind each other.
!!  isign=sign of exponential in transform
!!  work1(2,n4,n5,n6)=input complex array with alternating real and imaginary
!!   elements; data resides in 2*n4*n5*n6 of this array, spread out.
!!
!! OUTPUT
!!  work2(2,n4,n5,n6)=working space for transform and contains output
!!
!! NOTES
!! Calculates the discrete Fourier transform
!! This routine is optional, and called only when the FFTW libary is used.
!!
!!  WARNINGS:
!! - tested only on mips-sgi-irix6.5 (GNU `config.guess`) architecture
!! - using this subprogram requires an explicit interface in calling subprogram(s)
!! - temporary test for testing work1 and work2 sizes
!!
!! precompilation options
!! -DHAVE_FFT_FFTW2_THREADS: use of multithreading library
!!    set and export shell environment variable FFTW_NTHREADS to the
!!    desired value (2, 4, etc...) before running abinit
!!
!! TODO
!!
!! PARENTS
!!      ccfft
!!
!! CHILDREN
!!      fftw3d_f77_create_plan,fftw_f77_threads_init,fftwnd_f77_destroy_plan
!!      fftwnd_f77_one,fftwnd_f77_threads_one,getenv,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fftw(n1,n2,n3,isign,work1,work2)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftw'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  ! formal parameters
  ! external references (without known interfaces)
#if defined HAVE_FFT_FFTW2
  external :: fftw3d_f77_create_plan,fftwnd_f77_destroy_plan
#endif
!scalars
 integer,intent(in) :: n1,n2,n3,isign
!arrays
 real(dp),intent(in) :: work1(:,:,:,:)
 real(dp),intent(out) :: work2(:,:,:,:)
!no_abirules
#if defined HAVE_FFT_FFTW2
#ifdef HAVE_FFT_FFTW2_THREADS
  external :: getenv
  external :: fftw_f77_threads_init,fftwnd_f77_threads_one
#else
  external :: fftwnd_f77_one
#endif
#endif

!Local variables-------------------------------
  ! local constants
! name of this subprogram
  character(len=4), parameter :: sub_pgm = "fftw"
  ! (the following constants should be set separately)
#ifdef HAVE_FFT_FFTW2_THREADS
  ! static local variables
  ! non-static local variables
  ! (auxiliary non-static) local variables
#endif
!scalars
 integer,parameter :: FFTW_BACKWARD=1,FFTW_ESTIMATE=0,FFTW_FORWARD=-1
 integer,parameter :: FFTW_IN_PLACE=8
 integer :: err_code,j1,j2,j3,plan
 real(dp) :: xnorm
 character(len=500) :: message
#if defined HAVE_FFT_FFTW2_THREADS
 integer :: nthreads
#endif
#ifdef HAVE_FFT_FFTW2_THREADS
 logical,save :: first=.true.
 character(len=60) :: cthreads
#endif
!arrays
 complex,allocatable :: carr(:,:,:)

! *********************************************************************

!source code

 if(isign/=1 .and. isign/=-1)then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' fftw: BUG -',ch10,&
&   '  isign must be 1 or -1 '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 xnorm=one/(n1*n2*n3)

#ifdef HAVE_FFT_FFTW2_THREADS
 call getenv('FFTW_NTHREADS',cthreads)
 if (first) then
   call fftw_f77_threads_init(err_code)
   if (err_code /= 0) then
     write(message,fmt=2000) ch10,sub_pgm,err_code,ch10
     2000    format(a1," ",a,": error # ",i6, &
     " in executing 'call fftw_f77_threads_init'",a1)
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
   first=.false.
 end if
#endif

#if defined HAVE_FFT_FFTW2
 if (isign == FFTW_FORWARD) then
   call fftw3d_f77_create_plan(plan,n1,n2,n3,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
 else
   call fftw3d_f77_create_plan(plan,n1,n2,n3,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
 end if
#endif

 ABI_ALLOCATE(carr,(n1,n2,n3))
!err_code = ABI_ALLOC_STAT
!if (err_code /= 0) then
!write(message,fmt=1000) ch10,sub_pgm,err_code,ch10
!1000 format(a1," ",a,": error # ",i6, &
!" in executing 'allocate(carr(n1,n2,n3)'",a1)
!call wrtout(std_out,message,'PERS')
!call leave_new('PERS')
!end if

!temporary test
 if ((n1 > size(work1,dim=2)) .or. (n1 > size(work2,dim=2))) then
   write(std_out,*) 'pb 1 in fftw'
   stop
 end if
 if ((n2 > size(work1,dim=3)) .or. (n2 > size(work2,dim=3))) then
   write(std_out,*) 'pb 2 in fftw'
   stop
 end if
 if ((n3 > size(work1,dim=4)) .or. (n3 > size(work2,dim=4))) then
   write(std_out,*) 'pb 3 in fftw'
   stop
 end if

 do j3=1,n3
   do j2=1,n2
     do j1=1,n1
       carr(j1,j2,j3)=cmplx(work1(1,j1,j2,j3),work1(2,j1,j2,j3),kind=dpc)
     end do
   end do
 end do

#if defined HAVE_FFT_FFTW2_THREADS
 if (len_trim(cthreads) /= 0) then
   read(cthreads,'(i3)') nthreads
 else
   nthreads=1
 end if
 call fftwnd_f77_threads_one(nthreads,plan,carr,0)
#else
 call fftwnd_f77_one(plan,carr,0)
#endif

 do j3=1,n3
   do j2=1,n2
     do j1=1,n1
       work2(1,j1,j2,j3)=real(carr(j1,j2,j3))
       work2(2,j1,j2,j3)=aimag(carr(j1,j2,j3))
     end do
   end do
 end do

 ABI_DEALLOCATE(carr)
 err_code = ABI_ALLOC_STAT

#if defined HAVE_FFT_FFTW2
 call fftwnd_f77_destroy_plan(plan)
#endif

end subroutine fftw
!!***
