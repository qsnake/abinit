!{\src2tex{textfont=tt}}
!!****f* ABINIT/contract_dp_ge_val
!! NAME
!! contract_dp_ge_val
!!
!! FUNCTION
!!  "Design by contract" routine, ensuring that the
!!  target real (double precision) argument is greater or equal to some limiting value.
!!  Might be used to test composite quantities :
!!  dp_name might be "arg1-arg2*arg3", with corresponding input value
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  character(len=*) calling_routine = name of the routine of which the argument is checked
!!  character(len=*) dp_name = name of the real (double precision)à variable that is checked
!!  dp_value = value that is checked
!!  limit = lower admitted limit for the real , to which tol12 is subtracted, to
!!          deal in a portable way with the "equal" possibility
!!
!! OUTPUT
!!  (only checking, writing and stopping)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  "Design by contract" routines are simpler than the routines
!!  that check the input variables. Indeed, the error message
!!  can be much primitive.
!!
!! PARENTS
!!      dotprod_vn,dotprodm_vn
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine contract_dp_ge_val(calling_routine,dp_name,dp_value,limit)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'contract_dp_ge_val'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: dp_value,limit
 character(len=*),intent(in) :: calling_routine,dp_name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)' contract_dp_ge_val : enter '
!ENDDEBUG

 if(dp_value<limit-tol12)then
   write(message,'(10a,es16.6,5a,es16.6,a)') ch10,&
&   ' contract_dp_ge_val: BUG -',ch10,&
&   '  For the routine ',trim(calling_routine),',',ch10,&
&   '  the value of "',trim(dp_name),'" should be greater or equal to ',limit,'.',ch10,&
&   '  However, ',trim(dp_name),' = ',dp_value,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!DEBUG
!write(std_out,*)' contract_dp_ge_val : exit'
!stop
!ENDDEBUG

end subroutine contract_dp_ge_val
!!***
