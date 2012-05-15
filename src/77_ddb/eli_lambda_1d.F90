!{\src2tex{textfont=tt}}
!!****f* ABINIT/eli_lambda_1d
!!
!! NAME
!! eli_lambda_1d
!!
!! FUNCTION
!!  In the solving of the 1D (energy only) Eliashberg equations, calculate
!!  the lambda, which is the e-p coupling strength. See Allen and Mitrovic
!!  Solid State Physics vol 37 ed Ehrenreich Seitz and Turnbull, p.45
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   a2f_1d = 1D alpha2F function
!!   elph_ds = elphon dataset
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!
!! OUTPUT
!!   lambda_1d = coupling constant as a function of frequency
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      eliashberg_1d
!!
!! CHILDREN
!!      simpson_int
!!
!! NOTES
!!  lambda is used at points which are differences of Matsubara freqs,
!!  and hence is tabulated on points going through 0.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eli_lambda_1d (a2f_1d,elph_ds,lambda_1d,nmatsu,tc)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eli_lambda_1d'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
 real(dp),intent(in) :: tc
 type(elph_type),intent(in) :: elph_ds
!arrays
 real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
 real(dp),intent(out) :: lambda_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,iomega
 real(dp) :: nu_matsu,nu_matsu2,omega,domega
!arrays
 real(dp) :: lambda_int(elph_ds%na2f),tmplambda(elph_ds%na2f)

! *********************************************************************
!
!MG: the step should be calculated locally using nomega and the extrema of the spectrum.
!One should not rely on previous calls for the setup of elph_ds%domega
!I will remove elph_ds%domega since mka2f.F90 will become a method of gamma_t
 domega =elph_ds%domega

 do imatsu=-nmatsu,nmatsu
   nu_matsu = (two*imatsu)*pi*tc
   nu_matsu2 = nu_matsu*nu_matsu

   tmplambda(:) = zero
   omega=domega
   do iomega=2,elph_ds%na2f
     tmplambda(iomega) = a2f_1d(iomega) * two * omega / (nu_matsu2 + omega*omega)
     omega=omega+domega
   end do
   call simpson_int(elph_ds%na2f,domega,tmplambda,lambda_int)

   lambda_1d(imatsu) = lambda_int(elph_ds%na2f)
 end do

end subroutine eli_lambda_1d
!!***
