!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sgfft  
!! NAME
!!  m_sgfft      
!!
!! FUNCTION
!!  This module provides low-level interfaces to Goedecker's FFT library.
!!  Working in progress....
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_sgfft

 use m_profiling

 use defs_basis
 use m_errors

 use defs_fftdata,  only : mg

 implicit none

 private 

 public :: sg_fft_rc    ! Real-Complex version.

CONTAINS  !====================================================================
!!***


!!****f* m_fft_mesh/sg_fft_rc
!! NAME
!! sg_fft_rc
!!
!! FUNCTION
!! Conduct Fourier transform of REAL or COMPLEX function f(r)=fofr defined on
!! fft grid in real space, to create complex f(G)=fofg defined on full fft grid 
!! in reciprocal space, in full storage mode, or the reverse operation.
!! For the reverse operation, the final data is divided by nfftot.
!! REAL case when cplex=1, COMPLEX case when cplex=2. Usually used for density and potentials.
!!
!! There are two different possibilities :
!!  fftalgb=0 means using the complex-to-complex FFT routine,
!!   irrespective of the value of cplex
!!  fftalgb=1 means using a real-to-complex FFT or a complex-to-complex FFT,
!!   depending on the value of cplex.
!!  The only real-to-complex FFT available is from SGoedecker library.
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! isign=sign of Fourier transform exponent: current convention uses
!!  +1 for transforming from G to r 
!!  -1 for transforming from r to G.
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft)=f(G), complex.
!! fofr(cplex*nfft)=input function f(r) (real or complex)
!!
!! PARENTS
!!      fourdp
!!
!! CHILDREN
!!      sg_ctrig,sg_fftx,sg_ffty,sg_fftz
!!
!! SOURCE

subroutine sg_fft_rc(cplex,fofg,fofr,isign,nfft,ngfft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sg_fft_rc'
 use interfaces_52_fft_mpi_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: fofg(2,nfft),fofr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: mfac=11
 integer :: fftalg,fftalga,fftalgb,fftcache,i1,i2,i3,ic1,ic2,ic3,index
 integer :: n1,n1half1,n1halfm,n2,n2half1,n3,n4,n4half1,n5,n5half1,n6
 real(dp) :: ris,xnorm
 character(len=500) :: msg
!arrays
 integer :: aft1(mfac),aft2(mfac),aft3(mfac),bef1(mfac),bef2(mfac),bef3(mfac)
 integer :: ind1(mg),ind2(mg),ind3(mg),now1(mfac),now2(mfac),now3(mfac)
 real(dp) :: trig1(2,mg),trig2(2,mg),trig3(3,mg)
 real(dp),allocatable :: wk2d_a(:,:,:,:),wk2d_b(:,:,:,:),wk2d_c(:,:,:,:)
 real(dp),allocatable :: wk2d_d(:,:,:,:),work1(:,:,:,:),work2(:,:,:,:)

! *************************************************************************

 !DBG_ENTER("COLL")

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)

 fftcache=ngfft(8)
 fftalg  =ngfft(7)
 fftalga =fftalg/100
 fftalgb =mod(fftalg,100)/10

 ris=dble(isign)
 xnorm=1.0d0/dble(n1*n2*n3)

 if (fftalgb/=0 .and. fftalgb/=1) then
   write(msg, '(a,i4,a,a,a,a,a)' )&
&   '  The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   '  The second digit (fftalg(B)) must be 0 or 1.',ch10,&
&   '  Action : change fftalg in your input file.'
   MSG_PERS_BUG(msg)
 end if

 if (fftalgb==1 .and. ALL(fftalga/=(/1,3,4/)) )then
   write(msg,'(a,i4,5a)')&
&   '  The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   '  When fftalg(B) is 1, the allowed values for fftalg(A) are 1 and 4.',ch10,&
&   '  Action : change fftalg in your input file.'
   MSG_PERS_BUG(msg)
 end if

 if (n4<n1.or.n5<n2.or.n6<n3) then
   write(msg,'(a,3i8,a,3i8)')'  Each of n4,n5,n6=',n4,n5,n6,'must be >= n1, n2, n3 =',n1,n2,n3
   MSG_PERS_BUG(msg)
 end if

!---------------------------------------------------------
!Here sophisticated algorithm based on S. Goedecker routines, only for the REAL case.
!Take advantage of the fact that fofr is real, and that fofg has corresponding symmetry properties.

#ifdef DEBUG_MODE
 if (n1>mg .or. n2>mg .or. n3>mg) then 
   write(msg, '(a,3i10,a,a,a,i10,a)' )&
&   '  One of the dimensions n1,n2,n3=',n1,n2,n3,',',ch10,&
&   '  exceeds allowed dimension mg=',mg,'.'
   MSG_PERS_BUG(msg)
 end if
#endif

 n1half1=n1/2+1 ; n1halfm=(n1+1)/2
 n2half1=n2/2+1
!n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
 n4half1=(n1half1/2)*2+1
 n5half1=(n2half1/2)*2+1

!This sophisticated algorithm allows to decrease the memory needs.
 ABI_ALLOCATE(work1,(2,n4,n5half1,n6))
 ABI_ALLOCATE(work2,(2,n4,n5half1,n6))

 if(isign==1)then

!  Compute auxiliary arrays needed for FFTs, here forward FFT
   call sg_ctrig(n1,trig1,aft1,bef1,now1,one,ic1,ind1,mfac,mg)
   call sg_ctrig(n2,trig2,aft2,bef2,now2,one,ic2,ind2,mfac,mg)
   call sg_ctrig(n3,trig3,aft3,bef3,now3,one,ic3,ind3,mfac,mg)

!  Transfer fofg to the expanded fft box (only half of it)

!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n3,work1)
   do i3=1,n3
     do i2=1,n2half1
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
         work1(1,i1,i2,i3)=fofg(1,i1+index)
         work1(2,i1,i2,i3)=fofg(2,i1+index)
       end do
     end do
   end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SHARED(aft3,bef3,ind3,ic3)&
!$OMP&SHARED(now3,n1,n2half1,n4,n5half1,n6,ris,trig3,work1,work2)&
!$OMP&PRIVATE(i2)
   do i2=1,n2half1
     call sg_fftz(mfac,mg,n4,n5half1,n6,n1,i2,i2,work1,work2,&
&     trig3,aft3,now3,bef3,ris,ind3,ic3)
   end do
!$OMP END PARALLEL DO

!  Loop over x-y planes

!$OMP PARALLEL PRIVATE(i1,i2,i3,index,wk2d_a,wk2d_b,wk2d_c,wk2d_d) &
!$OMP&SHARED(aft1,aft2,bef1,bef2,fftcache,fofg,fofr,ic1,ic2,ind1,ind2) &
!$OMP&SHARED(n1,n1half1,n1halfm,n2,n2half1,n3) &
!$OMP&SHARED(n4,n5,now1,now2,ris,trig1,trig2,work2)

   ABI_ALLOCATE(wk2d_a,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_c,(2,2*n1halfm+1,n5,1))
   ABI_ALLOCATE(wk2d_d,(2,2*n1halfm+1,n5,1))

!$OMP DO
   do i3=1,n3

     do i2=1,n2half1
       do i1=1,n1
         wk2d_c(1,i1,i2,1)=work2(1,i1,i2,i3)
         wk2d_c(2,i1,i2,1)=work2(2,i1,i2,i3)
       end do
     end do

     call sg_fftx(fftcache,mfac,mg,2*n1halfm+1,n5,1,n2half1,1,wk2d_c,wk2d_d,&
&     trig1,aft1,now1,bef1,ris,ind1,ic1)

     do i1=1,n1half1-1 ! Compute symmetric and antisymmetric combinations
       wk2d_a(1,i1,1,1)=wk2d_d(1,2*i1-1,1,1)
       wk2d_a(2,i1,1,1)=wk2d_d(1,2*i1  ,1,1)
     end do

     if((2*n1half1-2)/=n1)then  ! If n1 odd, must add last data
       wk2d_a(1,n1half1,1,1)=wk2d_d(1,n1,1,1)
       wk2d_a(2,n1half1,1,1)=0.0d0
     end if

     do i2=2,n2half1
       do i1=1,n1half1-1
         wk2d_a(1,i1,i2,1)     = wk2d_d(1,2*i1-1,i2,1)-wk2d_d(2,2*i1,i2,1)
         wk2d_a(2,i1,i2,1)     = wk2d_d(2,2*i1-1,i2,1)+wk2d_d(1,2*i1,i2,1)
         wk2d_a(1,i1,n2+2-i2,1)= wk2d_d(1,2*i1-1,i2,1)+wk2d_d(2,2*i1,i2,1)
         wk2d_a(2,i1,n2+2-i2,1)=-wk2d_d(2,2*i1-1,i2,1)+wk2d_d(1,2*i1,i2,1)
       end do
       if((2*n1half1-2)/=n1)then
         wk2d_a(1,n1half1,i2,1)     = wk2d_d(1,n1,i2,1)
         wk2d_a(2,n1half1,i2,1)     = wk2d_d(2,n1,i2,1)
         wk2d_a(1,n1half1,n2+2-i2,1)= wk2d_d(1,n1,i2,1)
         wk2d_a(2,n1half1,n2+2-i2,1)=-wk2d_d(2,n1,i2,1)
       end if
     end do

     call sg_ffty(fftcache,mfac,mg,n4,n5,1,1,n1halfm,1,1,wk2d_a,wk2d_b,&
&     trig2,aft2,now2,bef2,ris,ind2,ic2)

     do i2=1,n2  ! Take real part data from expanded box and put it in the original box.
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1half1-1 ! copy data
         fofr(2*i1-1+index)=wk2d_b(1,i1,i2,1)
         fofr(2*i1  +index)=wk2d_b(2,i1,i2,1)
       end do
       if((2*n1half1-2)/=n1)then ! If n1 odd, must add last data
         fofr(n1+index)=wk2d_b(1,n1half1,i2,1)
       end if
     end do

   end do ! loop over x-y planes
!$OMP END DO
   ABI_DEALLOCATE(wk2d_a)
   ABI_DEALLOCATE(wk2d_b)
   ABI_DEALLOCATE(wk2d_c)
   ABI_DEALLOCATE(wk2d_d)

!$OMP END PARALLEL

 else if(isign==-1)then

!  Compute auxiliary arrays needed for FFTs, here backward FFT
   call sg_ctrig(n1,trig1,aft1,bef1,now1,-one,ic1,ind1,mfac,mg)
   call sg_ctrig(n2,trig2,aft2,bef2,now2,-one,ic2,ind2,mfac,mg)
   call sg_ctrig(n3,trig3,aft3,bef3,now3,-one,ic3,ind3,mfac,mg)

!  Treat first x-transform in x-y plane, and multiply
!  by overall normalization factor 1/nfftot

!  Loop over x-y planes

!$OMP PARALLEL PRIVATE(i1,i2,i3,index,wk2d_a,wk2d_b,wk2d_c,wk2d_d) &
!$OMP&SHARED(aft1,aft2,bef1,bef2,fftcache,fofr,ic1,ic2,ind1,ind2) &
!$OMP&SHARED(n1,n1half1,n1halfm,n2,n2half1,n3) &
!$OMP&SHARED(n4,n5,now1,now2,ris,trig1,trig2,work1,xnorm)

   ABI_ALLOCATE(wk2d_a,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_b,(2,n4,n5,1))
   ABI_ALLOCATE(wk2d_c,(2,2*n1halfm+1,n5,1))
   ABI_ALLOCATE(wk2d_d,(2,2*n1halfm+1,n5,1))

!$OMP DO
   do i3=1,n3
     do i2=1,n2
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1half1-1 ! copy and normalize data
         wk2d_a(1,i1,i2,1)=fofr(2*i1-1+index)*xnorm
         wk2d_a(2,i1,i2,1)=fofr(2*i1  +index)*xnorm
       end do

       if((2*n1half1-2)/=n1)then ! If n1 odd, must add last data
         wk2d_a(1,n1half1,i2,1)=fofr(n1+index)*xnorm
         wk2d_a(2,n1half1,i2,1)=zero
       end if
     end do

     call sg_ffty(fftcache,mfac,mg,n4,n5,1,1,n1halfm,1,1,wk2d_a,wk2d_b,&
&     trig2,aft2,now2,bef2,ris,ind2,ic2)

     do i1=1,n1halfm ! Decompose symmetric and antisymmetric parts
       wk2d_c(1,2*i1-1,1,1)=wk2d_b(1,i1,1,1)
       wk2d_c(2,2*i1-1,1,1)=0.0d0
       wk2d_c(1,2*i1,1,1)=wk2d_b(2,i1,1,1)
       wk2d_c(2,2*i1,1,1)=0.0d0
     end do

     do i2=2,n2half1
       do i1=1,n1halfm
         wk2d_c(1,2*i1-1,i2,1)= (wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+2-i2,1))*0.5d0
         wk2d_c(2,2*i1-1,i2,1)= (wk2d_b(2,i1,i2,1)-wk2d_b(2,i1,n2+2-i2,1))*0.5d0
         wk2d_c(1,2*i1,i2,1)  = (wk2d_b(2,i1,i2,1)+wk2d_b(2,i1,n2+2-i2,1))*0.5d0
         wk2d_c(2,2*i1,i2,1)  =-(wk2d_b(1,i1,i2,1)-wk2d_b(1,i1,n2+2-i2,1))*0.5d0
       end do
     end do

     call sg_fftx(fftcache,mfac,mg,2*n1halfm+1,n5,1,n2half1,1,wk2d_c,wk2d_d,&
&     trig1,aft1,now1,bef1,ris,ind1,ic1)

     do i2=1,n2half1
       do i1=1,n1
         work1(1,i1,i2,i3)=wk2d_d(1,i1,i2,1)
         work1(2,i1,i2,i3)=wk2d_d(2,i1,i2,1)
       end do
     end do

   end do
!$OMP END DO
   ABI_DEALLOCATE(wk2d_a)
   ABI_DEALLOCATE(wk2d_b)
   ABI_DEALLOCATE(wk2d_c)
   ABI_DEALLOCATE(wk2d_d)
!$OMP END PARALLEL

!$OMP PARALLEL DO SHARED(aft3,bef3,ind3,ic3)&
!$OMP&SHARED(now3,n1,n2half1,n4,n5half1,n6,ris,trig3,work1,work2)&
!$OMP&PRIVATE(i2)
   do i2=1,n2half1
     call sg_fftz(mfac,mg,n4,n5half1,n6,n1,i2,i2,work1,work2,&
&     trig3,aft3,now3,bef3,ris,ind3,ic3)
   end do
!$OMP END PARALLEL DO

!  Transfer fft output to the original fft box

!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofg,n1,n2,n2half1,n3,work2)
   do i3=1,n3
     do i2=1,n2half1
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
         fofg(1,i1+index)=work2(1,i1,i2,i3)
         fofg(2,i1+index)=work2(2,i1,i2,i3)
       end do
     end do
!    Complete missing values with complex conjugate
!    Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
     if(n2half1>2)then
       do i2=2,n2+1-n2half1
         index=n1*((n2+2-i2)-1)
         if(i3/=1)index=index+n1*n2*((n3+2-i3)-1)
         fofg(1,1+index)= work2(1,1,i2,i3)
         fofg(2,1+index)=-work2(2,1,i2,i3)
         do i1=2,n1
           fofg(1,n1+2-i1+index)= work2(1,i1,i2,i3)
           fofg(2,n1+2-i1+index)=-work2(2,i1,i2,i3)
         end do
       end do
     end if
   end do
!$OMP END PARALLEL DO

 end if ! choice of isign

 ABI_DEALLOCATE(work1)
 ABI_DEALLOCATE(work2)

 !DBG_EXIT("COLL")

end subroutine sg_fft_rc
!!***

END MODULE m_sgfft
!!***
