!{\src2tex{textfont=tt}}
!!****f* ABINIT/klocal
!! NAME
!! klocal
!!
!! FUNCTION
!! In a planewave basis set, the matrix of a local xc kernel:
!!  $f_{\rm xc}(\vec{r},\vec{r}') = f(\vec{r})\delta(\vec{r}-\vec{r}')$
!! is just:
!!  $f_{\rm xc}(\vec{G},\vec{G}') = f(\vec{G}-\vec{G}')$.
!! This subroutine calculates the matrix of such a local xc kernel
!! given $f(\vec{G})$ on the FFT grid.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, MF, XG, GMR, LSI, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License,see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ispxc = 1 for the up-up spin channel.
!!        = 2 for the up-down (and down-up) spin channels.
!!        = 3 for the down-down spin channel.
!!  ispxc must be 1 if nspden = 1.
!!  kg_diel(3,npwdiel) = reduced planewave coordinates for the kxc matrix.
!!  kxcg(2,nfft) = $f(\vec{G})$ on the FFT grid.
!!  nfft = number of fft grid points.
!!  ngfft(1:3) = integer fft box dimensions, see getng for ngfft(4:8).
!!  npwdiel = number of planewaves for the susceptibility matrix.
!!  nspden = number of spin-density components.
!!  option = 0 do not compute the first row and column of the matrix of the
!!             xc kernel (which we assume to the G = 0 row and column).
!!        /= 0 compute the full matrix of the xc kernel.
!!
!! OUTPUT
!!  kxc(2,npwdiel,nspden,npwdiel,nspden) = the matrix of the xc kernel.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      acfd_dyson,xcacfd
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine klocal(ispxc,kg_diel,kxc,kxcg,nfft,ngfft,npwdiel,nspden,option)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'klocal'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ispxc,nfft,npwdiel,nspden,option
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel),ngfft(18)
 real(dp),intent(in) :: kxcg(2,nfft)
 real(dp),intent(out) :: kxc(2,npwdiel,nspden,npwdiel,nspden)

!Local variables -------------------------------------------------------
!For debbuging purposes:
!real(dp) :: c1,c2,c3
!scalars
 integer :: i1,i2,i3,ifft,ipw1,ipw2,ipwstart,isp1,isp2,j1,j2,j3,k1,k2,k3,n1,n2
 integer :: n3
 logical :: ok
 character(len=500) :: message

!***********************************************************************

!Check input parameters.

 if (nspden > 2) then
   write (message,'(4a)') ch10,&
&   ' klocal: ERROR - ',ch10,&
&   '  klocal does not work yet for nspden > 2.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 isp1 = 1
 isp2 = 1
 ok = .true.
 if (nspden == 1) then
   select case (ispxc)
     case (1)
       isp1 = 1
       isp2 = 1
       case default
       ok = .false.
   end select
 else
   select case (ispxc)
     case (1)
       isp1 = 1
       isp2 = 1
     case (2)
       isp1 = 1
       isp2 = 2
     case (3)
       isp1 = 2
       isp2 = 2
       case default
       ok = .false.
   end select
 end if

 if (.not.ok) then
   write (message,'(4a,i10,a,i1,a)') ch10,&
&   ' klocal: BUG - ',ch10,&
&   '  The input ispxc = ',ispxc,' is not compatible with nspden = ',nspden,'.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if (option == 0) then
   ipwstart = 2
   kxc(:,1,isp1,:,isp2) = 0._dp
   kxc(:,:,isp1,1,isp2) = 0._dp
 else
   ipwstart = 1
 end if

!Calculate the xc matrix.

 n1 = ngfft(1) ; n2 = ngfft(2) ; n3 = ngfft(3)

 do ipw2 = ipwstart,npwdiel

   j1 = kg_diel(1,ipw2) ; j2 = kg_diel(2,ipw2) ; j3 = kg_diel(3,ipw2)

!  Fill the diagonal.

   kxc(:,ipw2,isp1,ipw2,isp2) = kxcg(:,1)

!  Fill the off-diagonal elements.

   do ipw1 = ipw2+1,npwdiel

     i1 = kg_diel(1,ipw1) ; i2 = kg_diel(2,ipw1) ; i3 = kg_diel(3,ipw1)

!    Compute the difference between G vectors.
!    The use of two mod calls handles both i1-j1 >= n1 AND i1-j1 < 0.

     k1 = mod(n1+mod(i1-j1,n1),n1)
     k2 = mod(n2+mod(i2-j2,n2),n2)
     k3 = mod(n3+mod(i3-j3,n3),n3)

     ifft = k1+n1*(k2+n2*k3)+1

     kxc(1,ipw1,isp1,ipw2,isp2) =  kxcg(1,ifft)
     kxc(2,ipw1,isp1,ipw2,isp2) =  kxcg(2,ifft)

     kxc(1,ipw2,isp1,ipw1,isp2) =  kxcg(1,ifft)
     kxc(2,ipw2,isp1,ipw1,isp2) = -kxcg(2,ifft)

   end do

 end do

!If needed, copy the up-down to the down-up spin channel.

 if (ispxc == 2) then
   do ipw2 = 1,npwdiel
     do ipw1 = 1,npwdiel
       kxc(1,ipw2,isp2,ipw1,isp1) =  kxc(1,ipw1,isp1,ipw2,isp2)
       kxc(2,ipw2,isp2,ipw1,isp1) = -kxc(2,ipw1,isp1,ipw2,isp2)
     end do
   end do
 end if

!DEBUG
!See kxc_alda.f, "test kernel" DEBUG section.
!do ipw2 = 1,npwdiel
!j1 = kg_diel(1,ipw2) ; j2 = kg_diel(2,ipw2) ; j3 = kg_diel(3,ipw2)
!do ipw1 = ipw2+1,npwdiel
!i1 = kg_diel(1,ipw1) ; i2 = kg_diel(2,ipw1) ; i3 = kg_diel(3,ipw1)
!k1 = mod(n1+mod(i1-j1,n1),n1)
!k2 = mod(n2+mod(i2-j2,n2),n2)
!k3 = mod(n3+mod(i3-j3,n3),n3)
!ifft = k1+n1*(k2+n2*k3)+1
!c1 = 0._dp ; c2 = 0._dp ; c3 = 0._dp
!if (i1-j1 ==  0) c1 = c1+0.0_dp
!if (i2-j2 ==  0) c2 = c2+0.0_dp
!if (i3-j3 ==  0) c3 = c3+0.0_dp
!if (i1-j1 ==  1) c1 = c1+0.5_dp
!if (i2-j2 ==  2) c2 = c2+0.5_dp
!if (i3-j3 ==  3) c3 = c3+0.5_dp
!if (i1-j1 == -1) c1 = c1+0.5_dp
!if (i2-j2 == -2) c2 = c2+0.5_dp
!if (i3-j3 == -3) c3 = c3+0.5_dp
!if ((abs(kxcg(1,ifft)-c1*c2*c3) > tol10).or.(abs(kxcg(2,ifft)) > tol10)) then
!write (std_out,*) ' i1 i2 i3 ifft: ',i1,i2,i3,ifft
!write (std_out,*) ' exp.: ',c1*c2*c3,' got: ',kxcg(:,ifft)
!end if
!end do
!end do
!ENDDEBUG

end subroutine klocal

!!***
