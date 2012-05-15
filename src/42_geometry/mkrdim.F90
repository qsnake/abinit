!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkrdim
!! NAME
!! mkrdim
!!
!! FUNCTION
!!  Trivial subroutine to make dimensional real space
!!  primitive translations from length scales acell(3)
!!  and dimensionless translations rprim(3,3).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  acell(3)=unit cell length scales (bohr)
!!  rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!              where: rprimd(i,j)=rprim(i,j)*acell(j)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      anaddb,bethe_salpeter,driver,elphon,ingeo,invars1,invars2m,irred_perts
!!      loper3,m_phdos,mkifc9,mkphbs,pred_bfgs,pred_delocint,pred_diisrelax
!!      pred_isothermal,pred_steepdesc,pred_verlet,predict_steepest
!!      predict_string,randomcellpos,rdddb9,rdm,refineblk
!!      scphon_new_frequencies,screening,setup1,setup_bse,setup_screening
!!      setup_sigma,sigma,suscep,symph3,thm9,thmeig,wvl_setboxgeometry
!!      xfpack_x2vin
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkrdim(acell,rprim,rprimd)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkrdim'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 real(dp),intent(out) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj

! *************************************************************************

 do ii=1,3
   do jj=1,3
     rprimd(ii,jj)=rprim(ii,jj)*acell(jj)
   end do
 end do

end subroutine mkrdim
!!***
