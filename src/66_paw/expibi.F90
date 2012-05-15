!{\src2tex{textfont=tt}}
!!****f* ABINIT/expibi
!! NAME
!! expibi
!!
!! FUNCTION
!! Routine which computes exp(-i (b2-b1).R) at each site. These
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
!!  rprimd(3,3) :: dimensioned primitive translations of real space lattice
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
!!      initberry
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine expibi(dtefield,gprimd,natom,rprimd,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield
 use m_errors

 use m_radmesh,  only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'expibi'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: natom
 type(efield_type),intent(inout) :: dtefield

!arrays
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3),xred(3,natom)

!Local variables---------------------------
!scalars
 integer :: bdir,bsig,iatom,kdir,mu,tind
 real(dp) :: bdotr
!arrays
 real(dp) :: bb(3),bcart(3),xcart(3)
! *************************************************************************

 dtefield%expibi(:,:,:) = zero
!expibi(2,natom,9) 
!used for PAW field calculations
!stores the on-site phase factors arising from
!$\langle\phi_{i,k+\sigma k_b}|\phi_{j,k+\sigma k_k}\rangle$ 
!where $\sigma = \pm 1$. These overlaps arise in various Berry
!phase calculations of electric and magnetic polarization. The on-site
!phase factor is $\exp[i(\sigma_b k_b - \sigma_k k_k)]$. Only the following
!are computed and saved, in the given order:
!1) -k_1 - k_2
!2) +k_1 - k_2
!3) -k_2 - k_3
!4) +k_2 - k_3
!5) -k_3 - k_1
!6) +k_3 - k_1
!7)    0 - k_1
!8)    0 - k_2
!9)    0 - k_3

 do iatom = 1, natom

   do bdir = 1, 3
     do bsig = -1, 1

       if(bsig /= 0) then
         kdir = mod(bdir,3)+1
       else 
         kdir = bdir
       end if
       
       bb(:)=bsig*dtefield%dkvecs(:,bdir)-dtefield%dkvecs(:,kdir)
       
       
!      get cartesian positions of atom in cell
       do mu=1,3
         bcart(mu)=dot_product(bb(:),gprimd(mu,:))
         xcart(mu)=dot_product(rprimd(mu,:),xred(:,iatom))
       end do
       bdotr = dot_product(xcart,bcart)

!      here is exp(i (b-k).R) for the given site

       if(bsig /= 0) then
         tind = 1+2*(bdir-1)+(bsig+1)/2
       else
         tind = kdir+6
       end if 
       dtefield%expibi(1,iatom,tind) = cos(two_pi*bdotr)
       dtefield%expibi(2,iatom,tind) = sin(two_pi*bdotr)

     end do ! end loop over bsig
   end do ! end loop over bdir

 end do ! end loop on natom

 dtefield%has_expibi = 2

 end subroutine expibi
!!***
