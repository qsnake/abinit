!{\src2tex{textfont=tt}}
!!****f* ABINIT/onestep
!! NAME
!! onestep
!!
!! FUNCTION
!! Advance one step following the gradient from vv(3).
!! It returns a new point in vv(3) and the value and gradient of the
!! electron density at this point in chg and grho(3)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npmax= maximum number of divisions
!!  hh= determines the initial value of the step (to be multiplied by grho)
!!
!! OUTPUT
!!  chg= value of electron density
!!  deltar= the length of the step thaty was needed
!!  grho(3)= gradient of electron density
!!  np= returns the number of divisions that were needed
!!
!! SIDE EFFECTS
!!  vv(3)=starting and updated point
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      aim_follow
!!
!! CHILDREN
!!      vgh_rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine onestep(vv,chg,grho,hh,np,npmax,deltar)

 use m_profiling

 use defs_basis
 use defs_parameters
 use defs_aimprom

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'onestep'
 use interfaces_63_bader, except_this_one => onestep
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npmax
 integer,intent(out) :: np
 real(dp),intent(in) :: hh
 real(dp),intent(out) :: chg,deltar
!arrays
 real(dp),intent(inout) :: vv(3)
 real(dp),intent(out) :: grho(3)

!Local variables ------------------------------
!scalars
 integer :: iat,ii,ipos,jj
 real(dp) :: dt,rr
!arrays
 real(dp) :: hrho(3,3),pom(3),vinter(3,200),vk(3),vkold(3)

!************************************************************************
 dt=hh
 np=1
 deltar=1._dp
 do ii=1,3
   vk(ii)=vv(ii)
 end do


 do while((np<3).or.((np<=npmax).and.(deltar>aim_deltarmin)))
   np=np*2
   dt=dt*0.5_dp
   do ii=1,3
     vkold(ii)=vk(ii)
   end do
   call vgh_rho(vk,chg,grho,hrho,rr,iat,ipos,0)
   do ii=1,3
     vinter(ii,1)=vv(ii)+dt*grho(ii)
   end do
   do jj=2,np
     call vgh_rho(vinter(1,jj-1),chg,grho,hrho,rr,iat,ipos,0)
     if(jj.eq.2) then
       do ii=1,3
         vinter(ii,2)=vv(ii)+2.0*dt*grho(ii)
       end do
     else
       do ii=1,3
         vinter(ii,jj)=vinter(ii,jj-2)+2.0*dt*grho(ii)
       end do
     end if
   end do

   call vgh_rho(vinter(1,np),chg,grho,hrho,rr,iat,ipos,0)
   do ii=1,3
     vinter(ii,np+1)=vinter(ii,np-1)+dt*grho(ii)
   end do

   deltar=0._dp
   do ii=1,3
     vk(ii)=(vinter(ii,np)+vinter(ii,np+1))*0.5_dp
     deltar=deltar+(vkold(ii)-vk(ii))*(vkold(ii)-vk(ii))
   end do
 end do

 pom(:)=vk(:)-vv(:)
 deltar=vnorm(pom,0)
 do ii=1,3
   vv(ii)=vk(ii)
 end do

 call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)
 if(deb) write(unto,*) ':VKf ',np,vk

end subroutine onestep
!!***
