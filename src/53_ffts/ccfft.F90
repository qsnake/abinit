!{\src2tex{textfont=tt}}
!!****f* ABINIT/ccfft
!! NAME
!! ccfft
!!
!! FUNCTION
!! Carry out complex-to-complex Fourier transforms between real
!! and reciprocal (G) space. Library of such routines.
!! Include machine-dependent F90 routines used with fftalg=200.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (PT, XG, FF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fftalga=govern the choice of the fft routine to be used
!!    if 1: SGoedecker routine
!!    if 2: Machine dependent routine, depending on the precompilation options
!!    if 3: FFTW library routine
!!    if 4: new SGoedecker routine, version 2002
!!          Warning : the second and third dimensions of the Fourier space
!!          array are switched, compared to the usual case
!!  fftcache=size of the cache (kB)
!!  isign= Integer specifying which sign to be used for the transformation.
!!         must be either +1 or -1.
!!  mpi_enreg=informations about MPI parallelization
!!  n1,n2,n3=Actual integer dimensions (see ngfft) for the 3D sequence.
!!           Physical dimension of the transform.
!!  n4,n5,n6=Leading dimensions. Generally, n6 is not different to n3.
!!  ndat=number of FFT to do in //
!!  option= 1 if call from fourwf, 2 if call from other routine
!!  paral_kgb=Flag related to the kpoint-band-fft parallelism
!!  work1(2,n4*n5*n6)=Array to be transformed.
!!
!! OUTPUT
!!  inplace = 0 if result is in work2 ; =1 if result is in work1 (machine dependent)
!!  normalized=0 if the backward (isign=-1) FFT is not normalized, so has
!!                      to be normalized outside of ccfft
!!            =1 otherwise
!!  work2(2,n4*n5*n6)=transformed array in case inplace=0.
!!
!! SIDE EFFECTS
!!  work1(2,n4*n5*n6)=at input, array to be transformed
!!                    at output, transformed array (in case inplace=1)
!!
!! NOTES
!! precompilation definitions :
!!   -D(machine_list) :  (case fftalga=200)
!!      choice of machine-dependent FFT library, if permitted
!!   -DHAVE_FFT_FFTW2   : (case fftalga=300) activate the FFTW lib
!!   -Dnolib  : (case fftalga=200) call SGoedecker routine,
!!      instead of machine-dependent one
!!
!! More about fftalga=200
!! Library routines for the following platforms have been implemented :
!!  Compaq/DEC
!!  HP          (in place FFT)
!!  SGI         (in place FFT)
!!  NEC         (in place FFT)
!! For all the other platforms, or if the CPP directive nolib is
!! activated, one uses the fft routine from S. Goedecker.
!!
!! PARENTS
!!      fourdp,fourwf
!!
!! CHILDREN
!!      back,fftw,forw,leave_new,sg_fft,wrtout,z3dfft,zfc3fb,zfft3d,zfft3di
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
& n1,n2,n3,n4,n5,n6,ndat,option,paral_kgb,work1,work2)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ccfft'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts, except_this_one => ccfft
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalga,fftcache,isign,n1,n2,n3,n4,n5,n6,ndat,option,paral_kgb
 integer,intent(out) :: inplace,normalized
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: work1(2,n4*n5*n6*ndat)
 real(dp),intent(out) :: work2(2,n4*n5*n6*ndat)

!Local variables ------------------------------
!scalars
 integer :: cplex,done,ierr,nd2proc,nd3proc
 real(dp) :: ris
 character(len=500) :: message
!no_abirules
!DEBUG
!integer :: i1,i2,i3
!real(dp) :: work3(2,n4*n5*n6)
!ENDDEBUG

!Interfaces -----------------------------------
!Machine-dependent declarations---------------------------

!---------------------------------------------------------

#if defined HAVE_FFT_ASL
 integer :: ifft,nfftot
 integer :: ifax(60)
 real(dp), allocatable, save :: trigs(:)
#endif

!---------------------------------------------------------


#if defined HAVE_FFT_SGIMATH
 integer,parameter :: nmax=15015
 integer,save :: n1_save,n2_save,n3_save
 integer(i4b) :: j1,j2,j3 ! local integers
 real(dp) :: xnorm ! Used to normalize the backward fft.
 complex(dpc),dimension(nmax),save :: coeffs  ! Array of at least ( (N1+15)+(N2+15)+(N3+15) )
!elements.  On entry it contains the
!Sines/Cosines and factorization of N. COEFFS
!needs to be initialized with a call to zfft3di.
 complex(dpc),dimension(:,:,:),allocatable :: zarray ! Complex work array containing the samples
!of the 3D sequence to be transformed.
 logical(lgt),save :: first=.true. ! logical used to initialize the coeffs calculation.
#endif

!*************************************************************************

!DEBUG
!write(std_out,*)' ccfft : enter, fftalga=',fftalga
!ENDDEBUG

 done=0
 if(fftalga==2)then

!  ---------------------------------------------------------

#if defined HAVE_FFT_MLIB
!  The hp routine make the FFT in place
   inplace=1 ; normalized=1
   call Z3DFFT(work1,n1,n2,n3,n4,n5,isign,ierr)
!  Check return code from Z3DFFT
   if(ierr /= 0)then
     write(message, '(a,a,a,a,i6,a,i6,a)' ) ch10,&
&     ' ccfft : BUG -',ch10,&
&     '  Z3DFFT isign=',isign,' error code =',ierr,'.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
   done=1
#endif

!  ---------------------------------------------------------

#if defined HAVE_FFT_ASL
!  Here library ASL NEC routine
!  Conventions for NEC FFT (see the ASL manual):
!  isign=1 implies a negative sign of the exponent
!  isign=-1 means a positive sign of the exponent
!  Normalization is not included in the ASL routine

   ABI_ALLOCATE(trigs,(2*(n4+n5+n6)))
   ifft = -isign
   nfftot = 2*n4*n5*n6
!  The NEC routine makes the FFT in place ; normalisation not included
   inplace=1 ; normalized=0
   call ZFC3FB(n1,n2,n3,work1,n4,n5,n6,ifft,ifax,trigs,work2,ierr)
!  Check return code from ZFC3FB
   if(ierr /= 0)then
     write(message, '(a,a,a,a,i6,a,i6,a)' )ch10,&
&     ' ccfft : BUG -',ch10,&
&     '  ZFC3FB isign=',ifft,' error code =',ierr,'.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
   ABI_DEALLOCATE(trigs)
   done=1
#endif

!  ---------------------------------------------------------

#if defined HAVE_FFT_SGIMATH

   if (first .or. n1/=n1_save .or. n2/=n2_save .or. n3/=n3_save) then
     if ((n1+n2+n3+45) > nmax) then
       write(message, '(a,a,a,a,a,a)' )ch10,&
&       ' ccfft : BUG -',ch10,&
&       '  n1+n2+n3+45>nmax.',ch10,&
&       '  Action : increase nmax in ccfft file.'
       call wrtout(std_out,message,'PERS')
       call leave_new('PERS')
     end if
     call zfft3di(n1,n2,n3,coeffs)
     first=.false. ; n1_save=n1 ; n2_save=n2 ; n3_save=n3
   end if

   inplace=1 ; normalized=0

   call zfft3d(isign,n1,n2,n3,work1,n4,n5,coeffs)
   done=1

#endif

!  ---------------------------------------------------------

 else if(fftalga==3)then

#if defined HAVE_FFT_FFTW2
   inplace=0 ; normalized=0
   call fftw(n1,n2,n3,isign,work1,work2)
   done=1
#else
!  If fftw is not accessible, one should not use fftalga=3
   write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&   ' ccfft : ERROR -',ch10,&
&   '  The library fftw has not been used to generate',ch10,&
&   '  the executable. Hence, fftalg(A)=3 is not allowed.',ch10,&
&   '  Action : check the value of fftalg in your input file.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
#endif

!  ---------------------------------------------------------

 else if(fftalga<1 .or. fftalga>4)then

   write(message, '(a,a,a,a,a,a,i5,a,a)' )ch10,&
&   ' ccfft : ERROR -',ch10,&
&   '  The allowed values of fftalg(A) are 1, 2, 3, and 4 .',ch10,&
&   '  The actual value of fftalg(A) is',fftalga,ch10,&
&   '  Action : check the value of fftalg in your input file.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')

 end if

!---------------------------------------------------------

!In case fftalg is allowed, but the FFT has not yet been done
 if(done==0)then
   inplace=0 ; normalized=0
   if(fftalga/=4)then
!    Call Stefan Goedecker FFT
     ris=real(isign,kind=dp)
     call sg_fft(fftcache,n4,n5,n6,n1,n2,n3,work1,work2,ris)
     ierr=0
   else if(fftalga==4)then
!    Call new version of Stefan Goedecker FFT
     cplex=2
     nd2proc=((n2-1)/mpi_enreg%nproc_fft) +1
     nd3proc=((n6-1)/mpi_enreg%nproc_fft) +1
     if(isign==1)then ! Fourier to Real space (backward)
!      call back(cplex,ndat,n1,n2,n3,n4,n5,n6,n4,n5,n6,nproc_fft,me_fft,paral_kgb,work1,work2)
       call back(cplex,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc, &
&       option,paral_kgb,work1,work2)
     else ! isign=-1, real space to Fourier (forward)
!      call forw(paral_kgb,cplex,ndat,n1,n2,n3,n4,n5,n6,n4,n5,n6,nproc_fft,me_fft,work1,work2)
       call forw(cplex,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc, &
&       option,paral_kgb,work1,work2)
     end if
!    DEBUG
!    write(std_out,*)' ccfft : i1,i2,i3, work2, work3, isign=',isign
!    write(std_out,*)'         n1,n2,n3,n4,n5,n6=',n1,n2,n3,n4,n5,n6
!    do i3=1,n6
!    do i2=1,n5
!    do i1=1,n4
!    write(std_out,'(3i3,4es16.6)')i1,i2,i3,work2(1:2,i1+(i2-1)*n4+(i3-1)*n4*n5),&
!    &                                       work3(1:2,i1+(i2-1)*n4+(i3-1)*n4*n5)
!    end do
!    end do
!    end do
!    stop
!    ENDDEBUG
   end if

 end if

!DEBUG
!write(std_out,*)' ccfft : exit '
!ENDDEBUG

end subroutine ccfft
!!***
