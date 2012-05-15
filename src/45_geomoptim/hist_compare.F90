!{\src2tex{textfont=tt}}
!!****f* ABINIT/hist_compare
!! NAME
!! hist_compare
!!
!! FUNCTION
!! Compare 2 HIST records
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  hist_in <type(ab_movehistory)>
!!  tolerance 
!!
!! OUTPUT
!!  similar= 1 the records are consitent
!!           0 the records are not consitent
!!
!! SIDE EFFECTS
!!  hist_out <type(ab_movehistory)>
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine hist_compare(hist_in,hist_out,natom,similar,tolerance)

 use m_profiling

use defs_basis
use defs_mover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hist_compare'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
integer,intent(out) :: similar
real(dp),intent(in) :: tolerance
type(ab_movehistory),intent(in) :: hist_in
type(ab_movehistory),intent(inout) :: hist_out

!arrays

!Local variables-------------------------------
!scalars
integer :: kk,jj
real(dp) :: maxdiff,diff
real(dp) :: x,y

! ***************************************************************

 similar=1

 write(std_out,*) 'Using values from history, iteration:',hist_in%ihist
 write(std_out,*) 'Differences between present history and values stored'
 write(std_out,*) 'on the previous history.(Relative difference)'

 x=hist_out%histXF(1,1,1,hist_out%ihist)
 y=hist_in%histXF(1,1,1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,natom
   do jj=1,3
     x=hist_out%histXF(jj,kk,1,hist_out%ihist)
     y=hist_in%histXF(jj,kk,1,hist_in%ihist)
     diff=2*abs(x-y)/(abs(x)+abs(y))
     if (diff>maxdiff) maxdiff=diff
   end do
 end do
 write(std_out,'(a,e12.5)') 'xcart:    ',maxdiff
 if (maxdiff>tolerance) similar=0


 x=hist_out%histXF(1,1,2,hist_out%ihist)
 y=hist_in%histXF(1,1,2,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,natom
   do jj=1,3
     x=hist_out%histXF(jj,kk,2,hist_out%ihist)
     y=hist_in%histXF(jj,kk,2,hist_in%ihist)
     diff=2*abs(x-y)/(abs(x)+abs(y))
     if (diff>maxdiff) maxdiff=diff
   end do
 end do
 write(std_out,'(a,e12.5)') 'xred:     ',maxdiff
 if (maxdiff>tolerance) similar=0


 x=hist_out%histR(1,1,hist_out%ihist)
 y=hist_in%histR(1,1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,3
   do jj=1,3
     x=hist_out%histR(jj,kk,hist_out%ihist)
     y=hist_in%histR(jj,kk,hist_in%ihist)
     diff=2*abs(x-y)/(abs(x)+abs(y))
     if (diff>maxdiff) maxdiff=diff
   end do
 end do
 write(std_out,'(a,e12.5)') 'rprimd:   ',maxdiff
 if (maxdiff>tolerance) similar=0


 x=hist_out%histA(1,hist_out%ihist)
 y=hist_in%histA(1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,3
   x=hist_out%histA(kk,hist_out%ihist)
   y=hist_in%histA(kk,hist_in%ihist)
   diff=2*abs(x-y)/(abs(x)+abs(y))
   if (diff>maxdiff) maxdiff=diff
 end do
 write(std_out,'(a,e12.5)') 'acell:    ',maxdiff
 if (maxdiff>tolerance) similar=0

 if (similar==1) then
   
   hist_out%histV(:,:,hist_out%ihist)=   hist_in%histV(:,:,hist_in%ihist)
   hist_out%histE(hist_out%ihist)=       hist_in%histE(hist_in%ihist)
   hist_out%histEk(hist_out%ihist)=      hist_in%histEk(hist_in%ihist)
   hist_out%histT(hist_out%ihist)=       hist_in%histT(hist_in%ihist)
   hist_out%histXF(:,:,3,hist_out%ihist)=hist_in%histXF(:,:,3,hist_in%ihist)
   hist_out%histXF(:,:,4,hist_out%ihist)=hist_in%histXF(:,:,4,hist_in%ihist)
   hist_out%histS(:,hist_out%ihist)=     hist_in%histS(:,hist_in%ihist)

 end if

end subroutine hist_compare
!!***
