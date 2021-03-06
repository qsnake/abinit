!{\src2tex{textfont=tt}}
!!****f* ABINIT/bound
!! NAME
!! bound
!!
!!
!! FUNCTION
!! For given kpt, ngfft, and gmet,
!!  Find distance**2 to boundary point of fft box nearest to kpt
!!  Find distance**2 to boundary point of fft box farthest to kpt
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kpt(3)=real input k vector (reduced coordinates)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  gmet(3,3)=reciprocal space metric (currently in Bohr**-2)
!!
!! OUTPUT
!!  dsqmax=maximum distance**2 from k to boundary in Bohr**-2.
!!  dsqmin=minimum distance**2 from k to boundary in Bohr**-2.
!!  gbound(3)=coords of G on boundary (correspnding to gsqmin)
!!  plane=which plane min occurs in (1,2, or 3 for G1,etc).
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!! Potential trouble: this routine was written assuming kpt lies inside
!! first Brillouin zone.  No measure is taken to fold input kpt back
!! into first zone.  Given arbitrary kpt, this will cause trouble.
!!
!! PARENTS
!!      getcut,getng
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bound'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: plane
 real(dp),intent(out) :: dsqmax,dsqmin
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(out) :: gbound(3)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer :: i1,i1min,i2,i2min,i3,i3min
 real(dp) :: dsm,dsp,dsq
 character(len=500) :: message

! *************************************************************************

 dsq(i1,i2,i3)=gmet(1,1)*(kpt(1)+dble(i1))**2&
& +gmet(2,2)*(kpt(2)+dble(i2))**2&
& +gmet(3,3)*(kpt(3)+dble(i3))**2&
& +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
& +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
& +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))

!Set plane to impossible value
 plane=0

!look at +/- g1 planes:
 dsqmax=zero
 dsqmin=dsq(ngfft(1)/2,-ngfft(2)/2,-ngfft(3)/2)+0.01_dp
 do i2=-ngfft(2)/2,ngfft(2)/2
   do i3=-ngfft(3)/2,ngfft(3)/2
     dsp = dsq(ngfft(1)/2, i2, i3)
     dsm = dsq( - ngfft(1)/2, i2, i3)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = ngfft(1)/2
       i2min = i2
       i3min = i3
       plane=1
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min =  - ngfft(1)/2
       i2min = i2
       i3min = i3
       plane=1
     end if
   end do
 end do
!
!+/- g2 planes:
 do i1=-ngfft(1)/2,ngfft(1)/2
   do i3=-ngfft(3)/2,ngfft(3)/2
     dsp = dsq(i1,ngfft(2)/2,i3)
     dsm = dsq(i1,-ngfft(2)/2,i3)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = i1
       i2min = ngfft(2)/2
       i3min = i3
       plane=2
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min = i1
       i2min =  - ngfft(2)/2
       i3min = i3
       plane=2
     end if
   end do
 end do
!
!+/- g3 planes:
 do i1=-ngfft(1)/2,ngfft(1)/2
   do i2=-ngfft(2)/2,ngfft(2)/2
     dsp = dsq(i1,i2,ngfft(3)/2)
     dsm = dsq(i1,i2,-ngfft(3)/2)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = i1
       i2min = i2
       i3min = ngfft(3)/2
       plane=3
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min = i1
       i2min = i2
       i3min =  - ngfft(3)/2
       plane=3
     end if
   end do
 end do

 if (plane==0) then
!  Trouble: missed boundary somehow
   write(message, '(a,a,a,a,a,a,3f9.4,a,3i5,a,a,a,a,a)' ) ch10,&
&   '  bound: BUG -',ch10,&
&   '  Trouble finding boundary of G sphere for',ch10,&
&   '  kpt=',kpt(:),' and ng=',ngfft(1:3),ch10,&
&   '  Action : check that kpt lies',&
&   ' reasonably within first Brillouin zone; ',ch10,&
&   '  else code bug, contact ABINIT group.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 gbound(1)=i1min
 gbound(2)=i2min
 gbound(3)=i3min

end subroutine bound
!!***
