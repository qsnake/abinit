!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_fftpad
!! NAME
!! sg_fftpad
!!
!! FUNCTION
!! Fast Fourier transform. This is the zero-padding version of "fft".
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
!!  mgfft=maximum size of 1D FFTs
!!  nd1,nd2,nd3=memory dimension of arr and ftarr
!!  n1,n2,n3=physical dimension of the transform
!!  arr(2,nd1,nd2,nd3)=input complex array with alternating real and imaginary
!!   elements; data resides in 2*n1*n2*n3 of this array, spread out.
!!  ris=(real(dp)) sign of exponential in transform
!!  gbound(2*mgfft+8,2)=sphere boundary info
!!
!! OUTPUT
!!  ftarr(2,nd1,nd2,nd3)=working space for transform and contains output
!!
!! SIDE EFFECTS
!!  arr(2,nd1,nd2,nd3) is modified by sg_fftpx,sg_ffty,sg_fftz.
!!
!! NOTES
!!  mfac sets maximum number of factors (5, 4, 3, or 2) which may be
!!  contained within any n1, n2, or n3
!!  mg sets the maximum 1 dimensional fft length (any one of n1, n2, or n3)
!!  XG: the signification of mg is changed with respect to fft3dp !!!
!!
!! PARENTS
!!      fftw3_fourwf,fourwf,m_wfs
!!
!! CHILDREN
!!      sg_ctrig,sg_fftpx,sg_ffty,sg_fftz
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine sg_fftpad(fftcache,mgfft,nd1,nd2,nd3,n1,n2,n3,arr,ftarr,ris,gbound)

 use m_profiling

 use defs_basis
 use m_errors

 use defs_fftdata,  only : mg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sg_fftpad'
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts, except_this_one => sg_fftpad
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftcache,mgfft,n1,n2,n3,nd1,nd2,nd3
 real(dp),intent(in) :: ris
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 real(dp),intent(inout) :: arr(2,nd1,nd2,nd3)
 real(dp),intent(out) :: ftarr(2,nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
 integer,parameter :: mfac=11
 integer :: g3max,g3min,i2,ic,n1i,n3i,n3p
#ifdef DEBUG_MODE
 character(len=500) :: message
#endif
!arrays
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp) :: trig(2,mg)

! *************************************************************************


#ifdef DEBUG_MODE
!Check that dimension is not exceeded
 if (n1>mg.or.n2>mg.or.n3>mg) then
   write(message, '(a,3i10,a,i10)')&
&   ' one of the dimensions n1,n2,n3=',n1,n2,n3,' exceeds the allowed dimension mg=',mg
   MSG_BUG(message)  
 end if
#endif

 g3min=gbound(3,2)
 g3max=gbound(4,2)

!--------------------------------------------------------------------------

 if (abs(ris-one)<tol12) then

!  Handle G -> r  transform (G sphere to fft box)

!  Transform along x direction
   call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)

!  Zero out the untransformed (0) data part of the work array
!  -- at every (y,z) there are 0 s to be added to the ends of
!  the x data so have to zero whole thing.
   ftarr(:,:,:,:)=0.0d0

!  Note the passing of the relevant part of gbound
   call sg_fftpx(fftcache,mfac,mg,mgfft,nd1,nd2,nd3,n2,n3,&
&   arr,ftarr,trig,aft,now,bef,ris,ind,ic,gbound(3,2))

!  Transform along y direction in two regions of z
   if (n2/=n1)then
     call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
   end if

!  First y transform: z=1..g3max+1
   n3p=g3max+1
   n1i=1 ; n3i=1
   call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3p,ftarr,arr,&
&   trig,aft,now,bef,ris,ind,ic)

!  Zero out the untransformed (0) data part of the work array
!  -- only need to zero specified ranges of z
   arr(:,:,:,n3p+1:g3min+n3)=0.0d0

!  Second y transform: z=g3min+1..0 (wrapped around)
   n3p=-g3min
   if (n3p>0) then
     n3i=1+g3min+n3 ; n1i=1
     call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,ftarr,arr,&
&     trig,aft,now,bef,ris,ind,ic)
   end if

!  Transform along z direction
   if (n3/=n2) then
     call sg_ctrig(n3,trig,aft,bef,now,ris,ic,ind,mfac,mg)
   end if

!  $OMP PARALLEL DO SHARED(aft,arr,bef,ftarr,ind,ic)&
!  $OMP&SHARED(nd1,nd2,nd3,now,n1,n2,ris,trig)&
!  $OMP&PRIVATE(i2)
   do i2=1,n2
     call sg_fftz(mfac,mg,nd1,nd2,nd3,n1,i2,i2,arr,ftarr,&
&     trig,aft,now,bef,ris,ind,ic)
   end do
!  $OMP END PARALLEL DO

 else

!  *************************************************
!  Handle r -> G transform (from fft box to G sphere)

!  Transform along z direction
   call sg_ctrig(n3,trig,aft,bef,now,ris,ic,ind,mfac,mg)

!  $OMP PARALLEL DO SHARED(aft,arr,bef,ftarr,ind,ic)&
!  $OMP&SHARED(nd1,nd2,nd3,now,n1,n2,ris,trig)&
!  $OMP&PRIVATE(i2)
   do i2=1,n2
     call sg_fftz(mfac,mg,nd1,nd2,nd3,n1,i2,i2,arr,ftarr,&
&     trig,aft,now,bef,ris,ind,ic)
   end do
!  $OMP END PARALLEL DO

!  Transform along y direction in two regions of z
   if (n2/=n3) then
     call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
   end if

!  First y transform: z=1..g3max+1
   n3p=g3max+1
   n1i=1 ; n3i=1
   call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3p,ftarr,arr,&
&   trig,aft,now,bef,ris,ind,ic)

!  Second y transform: z=g3min+1..0 (wrapped around)
   n3p=-g3min
   if (n3p>0) then
     n1i=1 ; n3i=1+g3min+n3
     call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,ftarr,arr,&
&     trig,aft,now,bef,ris,ind,ic)
   end if

!  Transform along x direction
   if (n1/=n2) then
     call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)
   end if

!  Zero out the untransformed (0) data part of the work array
!  -- at every (y,z) there are 0 s to be added to the ends of
!  the x data so have to zero whole thing.
   ftarr(:,:,:,:)=0.0d0

!  Note the passing of the relevant part of gbound
   call sg_fftpx(fftcache,mfac,mg,mgfft,nd1,nd2,nd3,n2,n3,&
&   arr,ftarr,trig,aft,now,bef,ris,ind,ic,gbound(3,2))

!  Data is now ready to be extracted from fft box to sphere

 end if


end subroutine sg_fftpad
!!***
