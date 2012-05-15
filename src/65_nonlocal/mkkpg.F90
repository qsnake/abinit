!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkkpg
!! NAME
!! mkkpg
!!
!! FUNCTION
!! Compute all (k+G) vectors (in reduced coordinates) for given k point.
!! Eventually compute related data
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  kg(3,npw)=integer coords of planewaves in basis sphere
!!  kpt(3)=k point in terms of recip. translations
!!  nkpg=second dimension of array kpg
!!  npw=number of plane waves in reciprocal space
!!
!! OUTPUT
!!  kpg(npw,3)= (k+G) components
!!  === if nkpg==9 ===
!!    kpg(npw,4:9)= [(k+G)_a].[(k+G)_b] quantities
!!
!! PARENTS
!!      calc_cs,ctocprj,debug_tools,dyfnl3,eltfrnl3,forstrnps,getcprj,ks_ddiago
!!      m_cprj_bspline,m_shirley,m_wfs,nonlop_ylm,nstpaw3,nstwf3
!!      prep_bandfft_tabs,resp3dte,rhofermi3,update_mmat,vtorho,vtorho3
!!      wfd_vnlpsi
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkkpg(kg,kpg,kpt,nkpg,npw)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkkpg'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpg,npw
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: kpg(npw,nkpg)

!Local variables-------------------------------
!scalars
 integer :: ipw,mu,mua,mub
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)

! *************************************************************************

 if (nkpg==0) return

!-- Test nkpg --
 if (nkpg/=3.and.nkpg/=9) then
   write(message, '(4a)' ) ch10,&
&   ' mkkpg : BUG -',ch10,&
&   '  Bad value for nkpg !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

!-- Compute (k+G) --
 do mu=1,3
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(two_pi,npw,kpg,kpt,kg)
   do ipw=1,npw
     kpg(ipw,mu)=kpt(mu)+dble(kg(mu,ipw))
   end do
!  $OMP END PARALLEL DO
 end do

!-- Compute [(k+G)_a].[(k+G)_b] --
 if (nkpg==9) then
   do mu=4,9
     mua=alpha(mu-3);mub=beta(mu-3)
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(mua,mub,mu,npw,kpg,kkpg)
     do ipw=1,npw
       kpg(ipw,mu)=kpg(ipw,mua)*kpg(ipw,mub)
     end do
!    $OMP END PARALLEL DO
   end do
 end if

end subroutine mkkpg
!!***
