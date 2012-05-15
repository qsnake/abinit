!{\src2tex{textfont=tt}}
!!****f* ABINIT/getlambda
!! NAME
!! getlambda
!!
!! FUNCTION
!! Return the abcissas and weights for the coupling constant integration.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, MF, XG, GMR, LSI, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License,see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  idyson = 1 solve the Dyson equation as linear system.
!!             (Gauss-Legendre mesh for the coupling constant integration).
!!         = 2 solve the Dyson equation as a differential equation
!!             (linear mesh for the coupling constant integration).
!!         = 3 solve the Dyson equation iteratively.
!!             (Gauss-Legendre mesh for the coupling constant integration).
!!  nlambda = number of mesh points.
!!
!! OUTPUT
!!  lambda(nlambda) = abscissas for the coupling constant integration.
!!  weight(nlambda) = weights for the coupling constant integration.
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!!
!! TODO
!!
!! PARENTS
!!      acfd_dyson
!!
!! CHILDREN
!!      coeffs_gausslegint,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine getlambda(idyson,lambda,nlambda,weight)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getlambda'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: idyson,nlambda
!arrays
 real(dp),intent(out) :: lambda(nlambda),weight(nlambda)

!Local variables -------------------------------------------------------
!WARNING : This should be moved to a def_* file !
!scalars
 integer,parameter :: solve_DE=2
 integer :: ii
 real(dp),parameter :: lambda_max=1._dp
 real(dp) :: lambda_step
 character(len=500) :: message

!***********************************************************************

!Check input parameters.

 if (nlambda < 1) then
   write (message,'(4a)') ch10,&
&   ' getlambda: BUG - ',ch10,&
&   '  nlambda (number of points for the coupling constant integration) must be >= 2.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Calculate abscissas and weights.

 lambda(1)       =       0._dp; weight(1)       = 0.5_dp
 lambda(nlambda) = lambda_max; weight(nlambda) = 0.5_dp

 if (nlambda > 2) then

   if (idyson == solve_DE) then

!    Linear mesh for the leap-frog solution of the Dyson equation.

     lambda_step = 1._dp/(nlambda-1)

     if (mod(nlambda,2) == 0) then

!      Trapezoidal rule for an even number of supports.

       do ii = 1,nlambda
         weight(ii) = lambda_step
         lambda(ii) = (ii-1)*lambda_step
       end do

       weight(1) = 0.5_dp*lambda_step
       weight(nlambda) = 0.5_dp*lambda_step

     else

!      Simpson rule for an odd number of supports.

       do ii = 1,nlambda
         lambda(ii) = (ii-1)*lambda_step
         weight(ii) = lambda_step/3._dp
       end do

       do ii = 2,nlambda-1,2
         weight(ii) = 4._dp*weight(ii)
       end do

       do ii = 3,nlambda-1,2
         weight(ii) = 2._dp*weight(ii)
       end do

     end if

   else

!    Gauss-Legendre quadrature.

     call coeffs_gausslegint(0._dp,1._dp,lambda(2),weight(2),nlambda-2)

     weight(1) = 0._dp; weight(nlambda) = 0._dp

   end if

 end if

end subroutine getlambda

!!***
