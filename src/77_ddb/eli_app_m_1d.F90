!{\src2tex{textfont=tt}}
!!****f* ABINIT/eli_app_m_1d
!!
!! NAME
!! eli_app_m_1d
!!
!! FUNCTION
!!   Apply the linearized Eliashberg matrix once to the input vector.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   lambda_1d = coupling constant as a function of frequency
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!   z_1d = renormalization Z as a function of frequency
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   delta_1d = imaginary gap function as a function of frequency changed
!!
!! PARENTS
!!      eli_m_iter_1d
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eli_app_m_1d (delta_1d,lambda_1d,nmatsu,z_1d)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eli_app_m_1d'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
!arrays
 real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu),z_1d(-nmatsu:nmatsu)
 real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,jmatsu,miguelflag
 real(dp) :: zfact
!arrays
 real(dp) :: delta_tmp(-nmatsu:nmatsu),freqfact(-nmatsu:nmatsu)

! *********************************************************************

 miguelflag = 0


 do imatsu=-nmatsu,nmatsu
   freqfact(imatsu) = one / abs(two*imatsu+one)
 end do

 delta_tmp(:) = delta_1d(:)

 if (miguelflag == 1) then
   do imatsu=-nmatsu,nmatsu
!    zfact = pi*tc / z_1d(imatsu)
     zfact = one / z_1d(imatsu)

     do jmatsu=max(-nmatsu,-nmatsu+imatsu),min(nmatsu,nmatsu+imatsu)
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + delta_1d(jmatsu) &
&       * lambda_1d(imatsu-jmatsu) &
&       * freqfact(jmatsu)
     end do
     delta_tmp(imatsu) = delta_tmp(imatsu)*zfact
   end do

 else

!  i < 0
   do imatsu=-nmatsu,-1

!    j < 0
     do jmatsu=max(-nmatsu,-nmatsu+imatsu),-1
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(jmatsu)*freqfact(jmatsu) &
&       - lambda_1d(imatsu-jmatsu)*delta_1d(imatsu)*freqfact(imatsu)
     end do
!    j > 0
     do jmatsu=0,min(nmatsu,nmatsu+imatsu)
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(jmatsu)*freqfact(jmatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(imatsu)*freqfact(imatsu)
     end do

   end do

!  i > 0
   do imatsu=0,nmatsu

!    j < 0
     do jmatsu=max(-nmatsu,-nmatsu+imatsu),-1
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(jmatsu)*freqfact(jmatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(imatsu)*freqfact(imatsu)
     end do
!    j > 0
     do jmatsu=0,min(nmatsu,nmatsu+imatsu)
       delta_tmp(imatsu) = delta_tmp(imatsu) &
&       + lambda_1d(imatsu-jmatsu)*delta_1d(jmatsu)*freqfact(jmatsu) &
&       - lambda_1d(imatsu-jmatsu)*delta_1d(imatsu)*freqfact(imatsu)
     end do

   end do

 end if

 delta_1d(:) = delta_tmp(:)

end subroutine eli_app_m_1d
!!***
