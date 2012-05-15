!{\src2tex{textfont=tt}}
!!****f* ABINIT/write_deltaI
!! NAME
!!  write_deltaI
!!
!! FUNCTION
!!  Writes the incompleteness matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2012 ABINIT group (GKA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  unt=The unit number of the file to be written (supposed to be already open)
!!  deltaI(npwe,npwe,1)=The matrix to be written, for a single q-point.
!!
!! OUTPUT
!!  (only writing on file)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine write_deltaI(unt,iqibz,deltaI,npwe)

 use m_profiling
    
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_deltaI'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer , intent(in)  :: unt,npwe,iqibz
!arrays
 complex(gwpc),intent(in) :: deltaI(npwe,npwe,1)

!Local variables-------------------------------
 integer :: ipwe,jpwe
 
! *************************************************************************

!DEBUG
!write (std_out,*) ' write_deltaI : enter'
!ENDDEBUG

 write(unt,'(2a,i3,x,2a)')ch10,'# Incompleteness matrix for q-point',iqibz,ch10,'# ig1, ig2    Re     Im '
 do ipwe=1,npwe
   do jpwe=1,npwe
     write(unt,'(i4,x,i4,x,2(f8.4))' ) ipwe,jpwe,deltaI(ipwe,jpwe,1)
   end do
   write(unt,'(x)' )
 end do

!DEBUG                                           ! to be uncommented, if needed
! if(option/=1 .and. option/=2 )then
!  write(message,'(a,a,a,a,a,a,i6)') ch10,&
!&  ' write_deltaI: BUG -',ch10,&
!&  '  The argument option should be 1 or 2,',ch10,&
!&  '  however, option=',option
!  call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
! endif
! if(sizein<1)then
!  write(message,'(a,a,a,a,a,a,i6)') ch10,&
!&  ' write_deltaI: BUG -',ch10,&
!&  '  The argument sizein should be a positive number,',ch10,&
!&  '  however, sizein=',sizein
!  call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
! endif
!ENDDEBUG


!DEBUG
!write (std_out,*) ' write_deltaI : exit'
!stop
!ENDDEBUG

end subroutine write_deltaI
!!***
