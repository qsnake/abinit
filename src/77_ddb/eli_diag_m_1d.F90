!{\src2tex{textfont=tt}}
!!****f* ABINIT/eli_diag_m_1d
!!
!! NAME
!! eli_diag_m_1d
!!
!! FUNCTION
!!  diagonalize M matrix. Heavy and should be avoided for production.
!!  Actually, since M is not symmetrical, diagonalize M^{t} M and
!!  get right-eigenvalues and vectors
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
!!   mustar = Coulomb potential parameter in Eliashberg equation
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!   z_1d = renormalization Z as a function of frequency
!!
!! OUTPUT
!!   maxeigval = estimation for maximum eigenvalue of M
!!
!! SIDE EFFECTS
!!   delta_1d = imaginary gap function as a function of frequency
!!
!! PARENTS
!!      eliashberg_1d
!!
!! CHILDREN
!!      dsyev
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eli_diag_m_1d (delta_1d,lambda_1d,maxeigval,mustar,nmatsu,tc,z_1d)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eli_diag_m_1d'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
 real(dp),intent(in) :: mustar,tc
 real(dp),intent(out) :: maxeigval
!arrays
 real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu),z_1d(-nmatsu:nmatsu)
 real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,info,jmatsu,kmatsu,lwork,tmiguel
 real(dp) :: si,sj,sqtimat,sqtjmat
!arrays
 real(dp) :: mtm_eig(2*nmatsu+1),symm_mtm(-nmatsu:nmatsu,-nmatsu:nmatsu)
 real(dp) :: work(3*(2*nmatsu+1))

! *********************************************************************

 tmiguel = 0

 if (tmiguel == 1) then
   do imatsu=-nmatsu,nmatsu
     do jmatsu=-nmatsu,nmatsu

       symm_mtm(imatsu,jmatsu) = zero
       do kmatsu=max(-nmatsu+imatsu,-nmatsu+jmatsu,-nmatsu),min(nmatsu+imatsu,nmatsu+jmatsu,nmatsu)
         symm_mtm(imatsu,jmatsu) = symm_mtm(imatsu,jmatsu) &
&         +  lambda_1d(kmatsu-imatsu)*lambda_1d(kmatsu-jmatsu) &
&         /  ( z_1d(kmatsu)*z_1d(kmatsu) )
       end do
!      symm_mtm(imatsu,jmatsu) = symm_mtm(imatsu,jmatsu) / ((two*imatsu+one)*(two*jmatsu+one))
       symm_mtm(imatsu,jmatsu) = symm_mtm(imatsu,jmatsu) * pi * tc * pi * tc
     end do
   end do

 else

   symm_mtm(:,:) = -mustar

   si = -one
   do imatsu=-nmatsu,nmatsu
     sqtimat = one / sqrt(two*abs(imatsu)+one)
     if (imatsu == 0) si = one
     sj = -one
     do jmatsu=max(-nmatsu,-nmatsu+imatsu),min(nmatsu,nmatsu+imatsu)
       sqtjmat = one / sqrt(two*abs(jmatsu)+one)
       if (jmatsu == 0) sj = one

       symm_mtm(imatsu,jmatsu) = symm_mtm(imatsu,jmatsu) &
&       + lambda_1d(imatsu-jmatsu)*sqtimat*sqtjmat

       symm_mtm(imatsu,imatsu) = symm_mtm(imatsu,imatsu) &
&       - lambda_1d(imatsu-jmatsu)*si*sj*sqtimat*sqtimat
     end do
   end do

 end if

 lwork = 3*(2*nmatsu+1)
 call DSYEV('V', 'U', 2*nmatsu+1, symm_mtm, 2*nmatsu+1, mtm_eig, work, lwork, info )

 write(std_out,*) 'last eigenvalues = '
 write(std_out,*) mtm_eig(2*nmatsu-9:2*nmatsu+1)

 do imatsu=-nmatsu,nmatsu
   delta_1d(imatsu) = symm_mtm(imatsu,nmatsu)*sqrt(two*abs(imatsu)+one)
 end do

 maxeigval = mtm_eig(2*nmatsu+1)

end subroutine eli_diag_m_1d
!!***
