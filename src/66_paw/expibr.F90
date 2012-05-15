!{\src2tex{textfont=tt}}
!!****f* ABINIT/expibr
!! NAME
!! expibr
!!
!! FUNCTION
!! Routine which computes exp(-i (b2-b1).r) on the fine grid around
!! each PAW sphere. These
!! quantities are stored in an efield structure and used in PAW Berrys Phase 
!! calculations of magnetization.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3) :: dimensioned primitive translations of reciprocal lattice
!!  natom :: number of atoms in unit cell
!!  pawfgrtab <type(pawfgrtab_type)>= atomic data given on fine rectangular grid
!!  xred(natom,3) :: reduced coordinates of atoms in unit cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtefield :: efield structure containing electric field and Berrys Phase data
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine expibr(dtefield,gprimd,natom,pawfgrtab,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'expibr'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: natom
 type(efield_type),intent(inout) :: dtefield

!arrays
 real(dp),intent(in) :: gprimd(3,3),xred(3,natom)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)

!Local variables---------------------------
!scalars
  integer :: iatom,ic,bdir,bsig,eijkl_ind,kdir,mu,nfgd_max
  real(dp) :: phase,phase_xred
!arrays
  real(dp) :: bb(3),bcart(3)
! *************************************************************************
!detefield%expibr(2,natom,nfgd,6) 
!used for PAW field calculations
!stores the on-site phase factors arising from
!$\exp(i(\sigma_b k_b - \sigma_k k_k)\cdot r)$ 
!where $\sigma = \pm 1$, on the fine grid of points in the PAW sphere
!around each atom. These phases are needed to compute the twisted $\hat{D}_{ij}$
!term in the magnetic field calculations.
!Only the following are computed and saved, in the given order:
!1) -k_1 - k_2
!2) +k_1 - k_2
!3) -k_2 - k_3
!4) +k_2 - k_3
!5) -k_3 - k_1
!6) +k_3 - k_1

 if (dtefield%has_expibr == 0) then
   nfgd_max = 0
   do iatom = 1, natom
     if(nfgd_max < pawfgrtab(iatom)%nfgd) nfgd_max = pawfgrtab(iatom)%nfgd
   end do
   ABI_ALLOCATE(dtefield%expibr,(2,natom,nfgd_max,6))
   dtefield%has_expibr = 1
 end if

 if (dtefield%has_expibr > 0) then
   dtefield%expibr(:,:,:,:) = zero

   do iatom = 1, natom
     do bdir = 1, 3
       do bsig = -1, 1, 2

         eijkl_ind = 2*bdir-1+(bsig+1)/2
         kdir = mod(bdir,3)+1

         bb(:) = bsig*dtefield%dkvecs(:,bdir) - &
&         dtefield%dkvecs(:,kdir)  ! k_bra - k_ket

         do mu=1,3
           bcart(mu)=dot_product(bb(:),gprimd(mu,:))
         end do
         phase_xred = two_pi*dot_product(bb(:),xred(:,iatom))

         do ic = 1, pawfgrtab(iatom)%nfgd
           phase = two_pi*dot_product(bcart(:),pawfgrtab(iatom)%rfgd(:,ic)) + &
&           phase_xred 
           dtefield%expibr(1,iatom,ic,eijkl_ind)=cos(phase)
           dtefield%expibr(2,iatom,ic,eijkl_ind)=sin(phase)
         end do ! end loop over nfgd for this atom
       end do ! end loop over bsig
     end do ! end loop over bdir
   end do ! end loop over natom

   dtefield%has_expibr = 2

 end if

 end subroutine expibr
!!***
