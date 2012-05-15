!{\src2tex{textfont=tt}}
!!****f* ABINIT/chki8
!!
!! NAME
!! chki8
!!
!! FUNCTION
!! This small subroutine check the identity of inti and intt,
!! who are integers, and eventually send a message and stop
!! if they are found unequal
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! inti=first integer
!! intt=second integer
!! character(len=6) name=name of the variable in the calling routine, to be echoed
!!
!! OUTPUT
!!  (only checking)
!!
!! PARENTS
!!      cmpar8
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chki8(inti,intt,name)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chki8'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: inti,intt
 character(len=6),intent(in) :: name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *********************************************************************

 if(inti/=intt) then
   write(message, '(a,a,a,a,a,a,a,i10,a,a,a,i10,a,a,a)' )&
&   ' chki8 : ERROR -',ch10,&
&   '  Comparing integers for variable',name,'.',ch10,&
&   '  Value from input DDB is',inti,' and',ch10,&
&   '        from transfer DDB is',intt,'.',ch10,&
&   '  Action : check your DDBs.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

end subroutine chki8
!!***
