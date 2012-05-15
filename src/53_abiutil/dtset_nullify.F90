!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtset_nullify
!! NAME
!! dtset_nullify
!!
!! FUNCTION
!! Nullify all pointers of a dataset before use.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MF, GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SIDE EFFECTS
!!  dtset <type(dataset_type)>=nullify all pointers.
!!
!! PARENTS
!!      m_ab6_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dtset_nullify(dtset)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtset_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(inout) :: dtset

!Local variables-------------------------------
 
! *************************************************************************
 
 nullify(dtset%acell_orig)
 nullify(dtset%amu)
 nullify(dtset%algalch)
 nullify(dtset%atvshift)
 nullify(dtset%bdgw)
 nullify(dtset%cd_imfrqs)
 nullify(dtset%corecs)
 nullify(dtset%densty)
 nullify(dtset%dmatpawu)
 nullify(dtset%dynimage)
 nullify(dtset%gw_freqsp)
 nullify(dtset%gw_qlwl)
 nullify(dtset%iatfix)
 nullify(dtset%iatsph)
 nullify(dtset%istwfk)
 nullify(dtset%jpawu)
 nullify(dtset%kberry)
 nullify(dtset%kpt)
 nullify(dtset%kptgw)
 nullify(dtset%kptns)
 nullify(dtset%lexexch)
 nullify(dtset%lpawu)
 nullify(dtset%mixalch)
 nullify(dtset%nband)
 nullify(dtset%occ_orig)
 nullify(dtset%ptcharge)
 nullify(dtset%qmass)
 nullify(dtset%qptdm)
 nullify(dtset%quadmom)
 nullify(dtset%rprim_orig)
 nullify(dtset%rprimd_orig)
 nullify(dtset%shiftk)
 nullify(dtset%so_psp)
 nullify(dtset%spinat)
 nullify(dtset%symafm)
 nullify(dtset%symrel)
 nullify(dtset%tnons)
 nullify(dtset%typat)
 nullify(dtset%upawu)
 nullify(dtset%vel_orig)
 nullify(dtset%wtatcon)
 nullify(dtset%wtk)
 nullify(dtset%xred_orig)
 nullify(dtset%ziontypat)
 nullify(dtset%znucl)

end subroutine dtset_nullify
!!***
