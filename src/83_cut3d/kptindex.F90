!{\src2tex{textfont=tt}}
!!****f* ABINIT/kptindex
!! NAME
!! kptindex
!!
!! FUNCTION
!! choosing out the wavefunction index for an ascending
!! ordered (eg, -pi/a to +pi/a) arrangement of wave vectors
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! shiftx,shifty,shiftz = shift options
!! ikpt_x,ikpt_y,ikpt_z = kpoint index (eg,1 to nkx)
!! nkx,nky,nkz = size of the kpoint grid
!!
!! OUTPUT
!! ikpt_x3,ikpt_y3,ikpt_z3 = desired wf index
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      localorb_S
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine kptindex(shiftx,shifty,shiftz,ikpt_x,ikpt_y,ikpt_z,nkx,nky,nkz,&
     &                     ikpt_x3,ikpt_y3,ikpt_z3)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kptindex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt_x,ikpt_y,ikpt_z,nkx,nky,nkz
 integer,intent(out) :: ikpt_x3,ikpt_y3,ikpt_z3
 real(dp),intent(in) :: shiftx,shifty,shiftz

!Local variables-------------------------------

! *************************************************************************

!-------------------------
!wave function index :
!--------------------------
 if(mod(nkx,2).eq.0)then
   if(shiftx .eq. 0.5)then
     if(ikpt_x.le.nkx/2) ikpt_x3 = ikpt_x + nkx/2
     if(ikpt_x.gt.nkx/2) ikpt_x3 = ikpt_x - nkx/2
   else
     if(ikpt_x.le.nkx/2-1) ikpt_x3 = ikpt_x + nkx/2 +1
     if(ikpt_x.gt.nkx/2-1) ikpt_x3 = ikpt_x - nkx/2 +1
   end if
 else
   if(ikpt_x.le.(nkx-1)/2) ikpt_x3 = ikpt_x + (nkx+1)/2
   if(ikpt_x.gt.(nkx-1)/2) ikpt_x3 = ikpt_x - (nkx-1)/2
 end if

 if(mod(nky,2).eq.0)then
   if(shifty .eq. 0.5)then
     if(ikpt_y.le.nky/2) ikpt_y3 = ikpt_y + nky/2
     if(ikpt_y.gt.nky/2) ikpt_y3 = ikpt_y - nky/2
   else
     if(ikpt_y.le.nky/2-1) ikpt_y3 = ikpt_y + nky/2 +1
     if(ikpt_y.gt.nky/2-1) ikpt_y3 = ikpt_y - nky/2 +1
   end if
 else
   if(ikpt_y.le.(nky-1)/2) ikpt_y3 = ikpt_y + (nky+1)/2
   if(ikpt_y.gt.(nky-1)/2) ikpt_y3 = ikpt_y - (nky-1)/2
 end if

 if(mod(nkz,2).eq.0)then
   if(shiftz .eq. 0.5)then
     if(ikpt_z.le.nkz/2) ikpt_z3 = ikpt_z + nkz/2
     if(ikpt_z.gt.nkz/2) ikpt_z3 = ikpt_z - nkz/2
   else
     if(ikpt_z.le.nkz/2-1) ikpt_z3 = ikpt_z + nkz/2 +1
     if(ikpt_z.gt.nkz/2-1) ikpt_z3 = ikpt_z - nkz/2 +1
   end if
 else
   if(ikpt_z.le.(nkz-1)/2) ikpt_z3 = ikpt_z + (nkz+1)/2
   if(ikpt_z.gt.(nkz-1)/2) ikpt_z3 = ikpt_z - (nkz-1)/2
 end if
 return

end subroutine kptindex
!!***
