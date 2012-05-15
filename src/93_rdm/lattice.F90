!{\src2tex{textfont=tt}}
!!****f* ABINIT/lattice
!! NAME
!! lattice
!!
!! FUNCTION
!! Calculate reciprocal lattice vectors and cell volumes
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! a1(3), a2(3), a3(3)= real-space lattice vectors in au
!!
!! OUTPUT
!! b1(3), b2(3), b3(3)= reciprocal-space lattice vectors in au**-1
!! ucvol= unit cell volume in au**3
!! bzvol= reciprocal cell volume
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine lattice(a1,a2,a3,b1,b2,b3,ucvol,bzvol)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lattice'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: bzvol,ucvol
!arrays
 real(dp),intent(in) :: a1(3),a2(3),a3(3)
 real(dp),intent(out) :: b1(3),b2(3),b3(3)

!Local variables-------------------------------
!scalars
 integer :: ii,lattype
 real(dp) :: a
 character(len=500) :: message
!arrays
 real(dp) :: aa(3,3)

! *************************************************************************

!calculate Brillouin zone and Unit cell volumes
 ucvol=a1(1)*(a2(2)*a3(3)-a2(3)*a3(2))+&
& a1(2)*(a2(3)*a3(1)-a2(1)*a3(3))+a1(3)*(a2(1)*a3(2)-a2(2)*a3(1))
 bzvol=8*pi*pi*pi/ucvol

!calculate reciprocal-space lattice vectors b1-b3
 b1(1)=2.0*pi*(a2(2)*a3(3)-a2(3)*a3(2))/ucvol
 b1(2)=2.0*pi*(a2(3)*a3(1)-a2(1)*a3(3))/ucvol
 b1(3)=2.0*pi*(a2(1)*a3(2)-a2(2)*a3(1))/ucvol
 b2(1)=2.0*pi*(a3(2)*a1(3)-a3(3)*a1(2))/ucvol
 b2(2)=2.0*pi*(a3(3)*a1(1)-a3(1)*a1(3))/ucvol
 b2(3)=2.0*pi*(a3(1)*a1(2)-a3(2)*a1(1))/ucvol
 b3(1)=2.0*pi*(a1(2)*a2(3)-a1(3)*a2(2))/ucvol
 b3(2)=2.0*pi*(a1(3)*a2(1)-a1(1)*a2(3))/ucvol
 b3(3)=2.0*pi*(a1(1)*a2(2)-a1(2)*a2(1))/ucvol

!find Bravais lattice type
 lattype=1
!fcc
 a=(ucvol*4.0)**(1.0/3.0)
 do ii=1,3
   aa(1,ii)=a1(ii)/(a/2)
   aa(2,ii)=a2(ii)/(a/2)
   aa(3,ii)=a3(ii)/(a/2)
   if(aa(1,ii)<=0 .and. aa(2,ii)<=0 .and. aa(3,ii)<=0) then
     aa(:,ii)=-aa(:,ii)
   end if
 end do
 if(aa(1,1)==1 .and. aa(1,2)==0 .and. aa(1,3)==1 .and.&
& aa(2,1)==0 .and. aa(2,2)==1 .and. aa(2,3)==1 .and.&
& aa(3,1)==1 .and. aa(3,2)==1 .and. aa(3,3)==0) then
   lattype=1
 end if
 if(lattype==1) then
!  write(message,'(a)') ' Bravais lattice: fcc'
!  call wrtout(std_out,message,'COLL')
   write(message,'(a,f5.2,a,3x,f5.2,2a)')' a = ',a,' [au]',a*Bohr_Ang,' Angstrom',ch10
   call wrtout(std_out,message,'COLL')
 end if

!tol10 enhance the portability
 write(message,'(2a,3(5x,3f8.2,a))')ch10,&
& ' real-space lattice vectors in cartesians [au]:',&
& a1+tol10,ch10, a2+tol10,ch10, a3+tol10,ch10
 call wrtout(std_out,message,'COLL')

 write(message,'(3a,3(5x,3f8.4,a))')ch10,&
& ' calculated reciprocal lattice vectors in cartesians [au]:',ch10,&
& b1,ch10, b2,ch10, b3,ch10
 call wrtout(std_out,message,'COLL')

 write(message,'(2a,f8.4,3a,f8.2,2a)')ch10,&
& ' brillouin zone volume = ',bzvol,' [au]',ch10,&
& ' unit cell volume      = ',ucvol,' [au]',ch10
 call wrtout(std_out,message,'COLL')

end subroutine lattice
!!***
