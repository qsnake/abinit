!{\src2tex{textfont=tt}}
!!****f* ABINIT/pmat_renorm
!! NAME
!! pmat_renorm
!!
!! FUNCTION
!! Renormalize the momentum matrix elements according to the scissor shift which is imposed
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mband= number of bands
!!  nkpt = number of k-points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  efermi = Fermi level
!!  sc = scissor shift for conduction bands
!!  evalv = eigenvalues for ground state
!!
!! OUTPUT
!!  pmat(2,mband,mband,nkpt,3,nsppol) = momentum matrix elements, renormalized by denominator change with scissor shift
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine pmat_renorm(efermi, evalv, mband, nkpt, nsppol, pmat, sc)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pmat_renorm'
!End of the abilint section

 implicit none

!Arguments -----------------------------------------------
!scalars
 integer, intent(in) :: nsppol
 integer, intent(in) :: nkpt
 integer, intent(in) :: mband
 real(dp), intent(in) :: efermi
 real(dp), intent(in) :: sc

!arrays
 real(dp), intent(in) :: evalv(mband,nsppol,nkpt)
!no_abirules
 complex(dpc), intent(inout) :: pmat(mband,mband,nkpt,3,nsppol)


!Local variables -----------------------------------------
!scalars
 integer :: iband1,iband2,ikpt,isppol

 real(dp) :: corec, e1, e2
!arrays

! *************************************************************************

 if (abs(sc) < tol8) then
   write(std_out,*) ' No scissor shift to be applied. Returning to main optic routine.'
   return
 end if

 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband1=1,mband ! valence states
       e1 = evalv(iband1,isppol,ikpt)
       if (e1 > efermi) cycle
       do iband2=1,mband ! conduction states
         e2 = evalv(iband2,isppol,ikpt)
         if (e2 < efermi) cycle
         corec = (e2-e1)/(e2+sc-e1)
         pmat(iband2,iband1,ikpt,:,isppol) = corec * pmat(iband2,iband1,ikpt,:,isppol)
         pmat(iband1,iband2,ikpt,:,isppol) = corec * pmat(iband1,iband2,ikpt,:,isppol)
       end do
     end do
   end do
 end do

end subroutine pmat_renorm
!!***
