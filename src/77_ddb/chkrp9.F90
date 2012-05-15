!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkrp9
!!
!! NAME
!! chkrp9
!!
!! FUNCTION
!! Check if the rprim used for the definition of the unit cell (in the
!! inputs) are consistent with the rprim used in the routine generating
!! the Big Box needed to generate the interatomic forces.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! brav=bravais lattice (1=simple lattice,2=face centered lattice,
!!  3=centered lattice,4=hexagonal lattice)
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!
!! OUTPUT
!!  (only checking)
!!
!! PARENTS
!!      mkifc9
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkrp9(brav,rprim)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkrp9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav
!arrays
 real(dp),intent(in) :: rprim(3,3)

!Local variables -------------------------
!scalars
 integer :: ii,jj
 character(len=500) :: message

! *********************************************************************

 if (brav==1) then
!  Simple Cubic Lattice
!  No condition in this case !
   continue

 else if (brav==2) then
!  Face Centered Lattice
   do ii=1,3
     do jj=1,3
       if (  ( ii==jj .and. abs(rprim(ii,jj))>tol10)&
&       .or.&
&       ( ii/=jj .and. abs(rprim(ii,jj)-.5_dp)>tol10) ) then
         write(message, '(a,a,a,a,a,a,a,a,a,a,a,a,a)' )&
&         ' chkrp9 : ERROR -',ch10,&
&         '  The input variable rprim does not correspond to the',ch10,&
&         '  fixed rprim to be used with brav=2 and ifcflag=1 :',ch10,&
&         '   0  1/2  1/2',ch10,&
&         '  1/2  0   1/2',ch10,&
&         '  1/2 1/2   0 ',ch10,&
&         '  Action : rebuild your DDB by using the latter rprim.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end do

 else if (brav==3) then
!  Body Centered Cubic Lattice
   do ii=1,3
     do jj=1,3
       if (  ( ii==jj .and. abs(rprim(ii,jj)+.5_dp)>tol10)&
&       .or.&
&       ( ii/=jj .and. abs(rprim(ii,jj)-.5_dp)>tol10) ) then
         write(message, '(a,a,a,a,a,a,a,a,a,a,a,a,a)' )&
&         ' chkrp9 : ERROR -',ch10,&
&         '  The input variable rprim does not correspond to the',ch10,&
&         '  fixed rprim to be used with brav=3 and ifcflag=1 :',ch10,&
&         '  -1/2  1/2  1/2',ch10,&
&         '   1/2 -1/2  1/2',ch10,&
&         '   1/2  1/2 -1/2',ch10,&
&         '  Action : rebuild your DDB by using the latter rprim.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end do

 else if (brav==4) then
!  Hexagonal Lattice
   if (abs(rprim(1,1)-1.0_dp)>tol10 .or. &
&   abs(rprim(3,3)-1.0_dp)>tol10 .or. &
&   abs(rprim(2,1)      )>tol10 .or. &
&   abs(rprim(3,1)      )>tol10 .or. &
&   abs(rprim(1,3)      )>tol10 .or. &
&   abs(rprim(2,3)      )>tol10 .or. &
&   abs(rprim(3,2)      )>tol10 .or. &
&   abs(rprim(1,2)+0.5_dp)>tol10 .or. &
&   abs(rprim(2,2)-0.5_dp*sqrt(3.0_dp))>tol10 ) then
     write(message, '(a,a,a,a,a,a,a,a,a,a,a,a,a)' )&
&     ' chkrp9 : ERROR -',ch10,&
&     '  The input variable rprim does not correspond to the',ch10,&
&     '  fixed rprim to be used with brav=4 and ifcflag=1 :',ch10,&
&     '   1      0      0',ch10,&
&     '  -1/2 sqrt[3]/2 0',ch10,&
&     '   0      0      1',ch10,&
&     '  Action : rebuild your DDB by using the latter rprim.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

 else

   write(message, '(a,a,a,i4,a,a,a,a,a)' )&
&   ' chkrp9 : ERROR -',ch10,&
&   '  The value of brav=',brav,' is not allowed.',ch10,&
&   '  Only  1,2,3 or 4 are allowed.',ch10,&
&   '  Action : change the value of brav in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')

 end if

end subroutine chkrp9
!!***
