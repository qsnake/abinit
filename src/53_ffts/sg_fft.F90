!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_fft
!! NAME
!! sg_fft
!!
!! FUNCTION
!! Calculates the discrete Fourier transform
!! ftarr(i1,i2,i3)=exp(ris*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) arr(j1,j2,j3)
!!
!! COPYRIGHT
!! Copyright by Stefan Goedecker, Ithaca, NY USA, July 14, 1993
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  nd1,nd2,nd3=memory dimension of arr and ftarr
!!  n1,n2,n3=physical dimension of the transform
!!  arr(2,nd1,nd2,nd3)=input complex array with alternating real and imaginary
!!  elements; data resides in 2*n1*n2*n3 of this array, spread out.
!!  (see SIDE FFECTS).
!!  ris=(real(dp)) sign of exponential in transform
!!
!! OUTPUT
!!  ftarr(2,nd1,nd2,nd3)=working space for transform and contains output
!!
!! SIDE EFFECTS
!!  arr(2,nd1,nd2,nd3) is modified by sg_fftx,sg_ffty,sg_fftz.
!!
!! NOTES
!!  ndi must always be greater or equal to ni.  Recommended choice for nd1
!!  and nd2 is: ni for ni=odd or ni+1 for ni=even (hence 2*(ni/2)+1);
!!  nd3 should always be n3.  Note that choosing nd1 or nd2 larger than
!!  the recommended value can severely degrade efficiency of this routine.
!!  Avoiding even ndi for nd1 and nd2 avoids cache conflicts on cache machines.
!!  Each of n1,n2,n3 must be a
!!  product of the prime factors 2,3,5. If two ni s are equal
!!  it is recommended to place them behind each other.
!!  The largest any of these may be is set by parameter "mg" below.
!!  This fft is particularly efficient for cache architectures.
!!  Note that the meaning of fftcache has changed from the original
!!  ncache of SG (that was the maximum number of COMPLEX*16 in the cache)
!!
!! TODO
!! Use latex for the equation above
!!
!! PARENTS
!!      ccfft
!!
!! CHILDREN
!!      leave_new,sg_ctrig,sg_fftx,sg_ffty,sg_fftz,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sg_fft(fftcache,nd1,nd2,nd3,n1,n2,n3,arr,ftarr,ris)

 use m_profiling

 use defs_basis
 use defs_fftdata

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sg_fft'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_52_fft_mpi_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftcache,n1,n2,n3,nd1,nd2,nd3
 real(dp),intent(in) :: ris
!arrays
 real(dp),intent(inout) :: arr(2,nd1,nd2,nd3)
 real(dp),intent(out) :: ftarr(2,nd1,nd2,nd3)

!Local variables-------------------------------
!mfac sets maximum number of factors (5, 4, 3, or 2) which may be
!contained within any n1, n2, or n3
!mg sets the maximum 1 dimensional fft length (any one of n1, n2, or n3)
!scalars
 integer,parameter :: mfac=11
 integer :: i2,ic,n1i,n3i
 character(len=500) :: message
!arrays
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp) :: trig(2,mg)

! *************************************************************************

!Check that dimension is not exceeded
 if (n1>mg.or.n2>mg.or.n3>mg) then
   write(message, '(a,a,a,a,3i10,a,i10,a)' ) ch10,&
&   ' sg_fft : BUG -',ch10,&
&   '  one of the dimensions n1,n2,n3=',n1,n2,n3,&
&   '  exceeds allowed dimension mg=',mg,ch10
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

!transform along x direction
 call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 call sg_fftx(fftcache,mfac,mg,nd1,nd2,nd3,n2,n3,&
& arr,ftarr,trig,aft,now,bef,ris,ind,ic)

!transform along y direction
 if (n2/=n1)then
   call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 end if
 n1i=1 ; n3i=1
 call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,&
& ftarr,arr,trig,aft,now,bef,ris,ind,ic)

!DEBUG
!if(abs(ris+one)<tol6)then
!ftarr=arr
!return
!end if
!ENDDEBUG

!transform along z direction
 if (n3/=n2)then
   call sg_ctrig(n3,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 end if

!$OMP PARALLEL DO SHARED(aft,arr,bef,ftarr,ind,ic)&
!$OMP SHARED(nd1,nd2,nd3,now,n1,n2,ris,trig)&
!$OMP PRIVATE(i2)
 do i2=1,n2
   call sg_fftz(mfac,mg,nd1,nd2,nd3,n1,i2,i2,arr,ftarr,&
&   trig,aft,now,bef,ris,ind,ic)
 end do
!$OMP END PARALLEL DO

end subroutine sg_fft
!!***
