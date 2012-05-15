!{\src2tex{textfont=tt}}
!!****f* ABINIT/k_rpa
!! NAME
!! k_rpa
!!
!! FUNCTION
!! Return the Hartree kernel:
!!  If option = 0, the bare Hartree kernel:
!!   krpa(ipw) = 4.0*pi/gsq(ipw) if gsq(ipw) /= 0.,
!!   krpa(ipw) = 0.0             if gsq(ipw) == 0. (1 <= ipw <= npw).
!!  If option /= 0, the Hartree kernel with a cut-off in real space beyond rcut_coulomb:
!!   krpa(ipw) = (4.0*pi/gsq(ipw))*(1.0-cos(sqrt(gsq(ipw))*rcut_coulomb)) if gsq(ipw) /= 0.,
!!   krpa(ipw) =  2.0*pi*rcut_coulomb**2                                  if gsq(ipw) == 0.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, MF, XG, GMR, LSI, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License,see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gsq(npw) = the squared norm of the planewaves.
!!  npw = number of planewaves in the gsq array.
!!  option = 0 for the bare Hartree kernel, /=0 for the cut-off Hartree kernel.
!!  rcut_coulomb = real space cut-off radius for the Coulomb interaction in Bohr.
!!
!! OUTPUT
!!  krpa(npw) = the Hartree kernel.
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!!
!! PARENTS
!!      acfd_dyson,acfd_intexact
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine k_rpa(gsq,krpa,npw,option,rcut_coulomb)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k_rpa'
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: npw,option
 real(dp),intent(in) :: rcut_coulomb
!arrays
 real(dp),intent(in) :: gsq(npw)
 real(dp),intent(out) :: krpa(npw)

!Local variables -------------------------------------------------------
!scalars
 integer :: ipw

!***********************************************************************

 if (option == 0) then

!  Compute the bare Hartree kernel.

   do ipw = 1,npw

     if (gsq(ipw) > tol12) then
       krpa(ipw) = four_pi/gsq(ipw)
     else
       krpa(ipw) = 0._dp
     end if

   end do

 else

!  Compute the Hartree kernel with a cut-off in real space beyond rcut_coulomb:

   do ipw = 1,npw

     if (gsq(ipw) > tol12) then
       krpa(ipw) = (four_pi/gsq(ipw))*(1._dp-cos(sqrt(gsq(ipw))*rcut_coulomb))
     else
       krpa(ipw) = two_pi*rcut_coulomb**2
     end if

   end do

 end if

end subroutine k_rpa

!!***
