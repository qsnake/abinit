!{\src2tex{textfont=tt}}
!!****f* ABINIT/vander
!! NAME
!! vander
!!
!! FUNCTION
!!   Returns a cutoff function to use as a local potential for making a Kleinman Bylander 
!!   potential out of a SIESTA XML pseudopotential file.
!!   The famous "Vanderbilt generalized cutoff"
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


function vander(a,x) result(f)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vander'
!End of the abilint section

 implicit none

!Arguments ----------------------------
! Generalized gaussian shape
!scalars
 real(dp) :: f
 real(dp),intent(in) :: a,x

!Local variables ----------------------
!! real(dp), parameter :: exp_range = 40.0_dp
!scalars
 real(dp),parameter :: log10_e=0.4343d0
 real(dp) :: gexp
!no_abirules
 real(dp),parameter :: exp_range=(range(1.0_dp)-1)/log10_e

! *************************************************************************

 gexp = sinh( a * x ) / sinh(a)
 gexp = gexp * gexp

 if (gexp .lt. exp_range) then
   f=exp(-gexp)
 else
   f = 0.0_dp
 end if
end function vander
!!***
