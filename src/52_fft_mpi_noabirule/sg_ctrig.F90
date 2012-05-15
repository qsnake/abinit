!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_ctrig
!! NAME
!! sg_ctrig
!!
!! FUNCTION
!! Precalculates trigonometric expressions and bitreversal key IND (Stefan Goedecker lib).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (SG, XG)
!! Copyright (C) Stefan Goedecker, Ithaca, NY USA, July 14, 1993
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! n=Number of FFT points for 1D FFT.
!! ris  = sign of exponential in transform (should be 1 or -1; real)
!! mfac = maximum number of factors in 1D FFTs
!! mg   = maximum length of 1D FFTs
!!
!! OUTPUT
!! trig(2,mg) TO BE DESCRIBED SB 090902
!! aft(mfac) TO BE DESCRIBED SB 090902
!! bef(mfac) TO BE DESCRIBED SB 090902
!! now(mfac) TO BE DESCRIBED SB 090902
!! ic = number of (radix) factors of x transform length (from ctrig)
!! ind(mg) TO BE DESCRIBED SB 090902
!!
!! NOTES
!! This version of sg_ctrig produces cos and tan instead of sin and cos--
!! this allows for much greater efficiency on the superscalar architecture
!! of ibm rs6000 where floating point multiply and add (FMA) is used.
!!
!! TODO
!! Should describe arguments
!! Should suppress one-letter variables
!!
!! PARENTS
!!      m_sgfft,sg_fft,sg_fftpad,sg_fftrisc,sg_fftrisc_2
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sg_ctrig(n,trig,aft,bef,now,ris,ic,ind,mfac,mg)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sg_ctrig'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mfac,mg,n
 integer,intent(out) :: ic
 real(dp),intent(in) :: ris
!arrays
 integer,intent(out) :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp),intent(out) :: trig(2,mg)

!Local variables-------------------------------
!scalars
 integer,save :: nextmx=4
 integer :: i,ii,inc,irep,j,k,l,next,nh
!no_abirules
!"prime" is the set of radices coded elsewhere for fft
 integer,save :: prime(4)=(/5,4,3,2/)
 real(dp) :: angle,trigc,trigs,twopi
 character(len=500) :: message

! *************************************************************************

!**Note**
!2*Pi must not be defined too accurately here or else
!cos(twopi/2) will be exactly 0 and sin/cos below will be
!infinite; if a small error is left in Pi, then sin/cos will
!be about 10**14 and later cos * (sin/cos) will be 1 to within
!about 10**(-14) and the fft routines will work
!The precision on sgi causes the algorithm to fail if
!twopi is defined as 8.d0*atan(1.0d0).

 twopi=6.2831853071795867d0

 angle=ris*twopi/n
!trig(1,0)=1.d0
!trig(2,0)=0.d0
 if (mod(n,2)==0) then
   nh=n/2
   trig(1,nh)=-1.d0
   trig(2,nh)=0.d0
   do i=1,nh-1
     trigc=cos(i*angle)
     trigs=sin(i*angle)
     trig(1,i)=trigc
     trig(2,i)=trigs/trigc
     trig(1,n-i)=trigc
     trig(2,n-i)=-trigs/trigc
   end do
 else
   nh=(n-1)/2
   do i=1,nh
     trigc=cos(i*angle)
     trigs=sin(i*angle)
     trig(1,i)=trigc
     trig(2,i)=trigs/trigc
     trig(1,n-i)=trigc
     trig(2,n-i)=-trigs/trigc
   end do
 end if

 ic=1
 aft(ic)=1
 bef(ic)=n
 next=1

!An infinite loop, with exit or cycle instructions
 do
   if( (bef(ic)/prime(next))*prime(next)<bef(ic) ) then
     next=next+1
     if (next<=nextmx) then
       cycle
     else
       now(ic)=bef(ic)
       bef(ic)=1
     end if
   else
     now(ic)=prime(next)
     bef(ic)=bef(ic)/prime(next)
   end if
   aft(ic+1)=aft(ic)
   now(ic+1)=now(ic)
   bef(ic+1)=bef(ic)
   ic=ic+1
   if (ic>mfac) then
     write(message, '(4a,i10,2a,i5,a)' ) ch10, &
&     ' sg_ctrig: BUG -',ch10,&
&     '  number of factors ic=',ic,ch10,&
&     '  exceeds dimensioned mfac=',mfac,ch10
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
   if (bef(ic)/=1) then
     aft(ic)=aft(ic)*now(ic)
     cycle
   end if
!  If not cycled, exit
   exit
 end do

 ic=ic-1

!DEBUG
!write(std_out,*) 'now',(now(i),i=1,ic)
!write(std_out,*) 'aft',(aft(i),i=1,ic)
!write(std_out,*) 'bef',(bef(i),i=1,ic)
!ENDDEBUG

 do i=1,n
   ind(i)=1
 end do

 irep=1
 inc=n
 do l=ic,1,-1
   inc=inc/now(l)
   ii=0
   do k=1,1+(n-1)/(now(l)*irep)
     do j=0,now(l)-1
       do i=1,irep
         ii=ii+1
         ind(ii)=ind(ii)+j*inc
       end do
     end do
   end do
   irep=irep*now(l)
 end do

 if (irep/=n) then
   write(message, '(a,a,a,a,i10,a,i10)' ) ch10,&
&   ' sg_ctrig : BUG -',ch10,&
&   '  irep should equal n ; irep=',irep,' n=',n
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 if (inc/=1) then
   write(message, '(a,a,a,a,i10)' ) ch10,&
&   ' sg_ctrig : BUG -',ch10,&
&   '  inc should equal 1 in sg_ctrig; inc=',inc
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

end subroutine sg_ctrig
!!***
