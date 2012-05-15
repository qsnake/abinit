!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtsetfree
!! NAME
!! dtsetfree
!!
!! FUNCTION
!! Free a dataset after use.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MF, GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SIDE EFFECTS
!!  dtset <type(dataset_type)>=free all associated pointers.
!!
!! PARENTS
!!      afterscfloop,chkinp,cvxclda,driver,kxc_alda,m_ab6_invars_f90,m_io_kss
!!      xc_kernel,xc_kernel_ADA
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dtsetFree(dtset)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtsetFree'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(inout) :: dtset

! *************************************************************************

 if (associated(dtset%acell_orig))  then
   ABI_DEALLOCATE(dtset%acell_orig)
 end if
 if (associated(dtset%algalch))     then
   ABI_DEALLOCATE(dtset%algalch)
 end if
 if (associated(dtset%atvshift))    then
   ABI_DEALLOCATE(dtset%atvshift)
 end if
 if (associated(dtset%bdgw))        then
   ABI_DEALLOCATE(dtset%bdgw)
 end if
 if (associated(dtset%cd_imfrqs))   then
   ABI_DEALLOCATE(dtset%cd_imfrqs)
 end if
 if (associated(dtset%corecs))      then
   ABI_DEALLOCATE(dtset%corecs)
 end if
 if (associated(dtset%dynimage))    then
   ABI_DEALLOCATE(dtset%dynimage)
 end if
 if (associated(dtset%iatfix))      then
   ABI_DEALLOCATE(dtset%iatfix)
 end if
 if (associated(dtset%iatsph))      then
   ABI_DEALLOCATE(dtset%iatsph)
 end if
 if (associated(dtset%istwfk))      then
   ABI_DEALLOCATE(dtset%istwfk)
 end if
 if (associated(dtset%kberry))      then
   ABI_DEALLOCATE(dtset%kberry)
 end if
 if (associated(dtset%lexexch))     then
   ABI_DEALLOCATE(dtset%lexexch)
 end if
 if (associated(dtset%lpawu))       then
   ABI_DEALLOCATE(dtset%lpawu)
 end if
 if (associated(dtset%nband))       then
   ABI_DEALLOCATE(dtset%nband)
 end if
 if (associated(dtset%qmass))       then
   ABI_DEALLOCATE(dtset%qmass)
 end if
 if (associated(dtset%so_psp))      then
   ABI_DEALLOCATE(dtset%so_psp)
 end if
 if (associated(dtset%symafm))      then
   ABI_DEALLOCATE(dtset%symafm)
 end if
 if (associated(dtset%symrel))      then
   ABI_DEALLOCATE(dtset%symrel)
 end if
 if (associated(dtset%typat))       then
   ABI_DEALLOCATE(dtset%typat)
 end if
 if (associated(dtset%amu))         then
   ABI_DEALLOCATE(dtset%amu)
 end if
 if (associated(dtset%densty))      then
   ABI_DEALLOCATE(dtset%densty)
 end if
 if (associated(dtset%dmatpawu))    then
   ABI_DEALLOCATE(dtset%dmatpawu)
 end if
 if (associated(dtset%gw_freqsp))     then
   ABI_DEALLOCATE(dtset%gw_freqsp)
 end if
 if (associated(dtset%gw_qlwl))     then
   ABI_DEALLOCATE(dtset%gw_qlwl)
 end if
 if (associated(dtset%jpawu))       then
   ABI_DEALLOCATE(dtset%jpawu)
 end if
 if (associated(dtset%kpt))         then
   ABI_DEALLOCATE(dtset%kpt)
 end if
 if (associated(dtset%kptgw))       then
   ABI_DEALLOCATE(dtset%kptgw)
 end if
 if (associated(dtset%kptns))       then
   ABI_DEALLOCATE(dtset%kptns)
 end if
 if (associated(dtset%mixalch))     then
   ABI_DEALLOCATE(dtset%mixalch)
 end if
 if (associated(dtset%occ_orig))    then
   ABI_DEALLOCATE(dtset%occ_orig)
 end if
 if (associated(dtset%pimass))      then
   ABI_DEALLOCATE(dtset%pimass)
 end if
 if (associated(dtset%prtatlist))   then
   ABI_DEALLOCATE(dtset%prtatlist)
 end if
 if (associated(dtset%ptcharge))    then
   ABI_DEALLOCATE(dtset%ptcharge)
 end if
 if (associated(dtset%quadmom))     then
   ABI_DEALLOCATE(dtset%quadmom)
 end if
 if (associated(dtset%qptdm))       then
   ABI_DEALLOCATE(dtset%qptdm)
 end if
 if (associated(dtset%ratsph))      then
   ABI_DEALLOCATE(dtset%ratsph)
 end if
 if (associated(dtset%rprim_orig))  then
   ABI_DEALLOCATE(dtset%rprim_orig)
 end if
 if (associated(dtset%rprimd_orig))ABI_DEALLOCATE(dtset%rprimd_orig)
 if (associated(dtset%shiftk))      then
   ABI_DEALLOCATE(dtset%shiftk)
 end if
 if (associated(dtset%spinat))      then
   ABI_DEALLOCATE(dtset%spinat)
 end if
 if (associated(dtset%tnons))       then
   ABI_DEALLOCATE(dtset%tnons)
 end if
 if (associated(dtset%upawu))       then
   ABI_DEALLOCATE(dtset%upawu)
 end if
 if (associated(dtset%vel_orig))    then
   ABI_DEALLOCATE(dtset%vel_orig)
 end if
 if (associated(dtset%wtatcon))     then
   ABI_DEALLOCATE(dtset%wtatcon)
 end if
 if (associated(dtset%wtk))         then
   ABI_DEALLOCATE(dtset%wtk)
 end if
 if (associated(dtset%xred_orig))   then
   ABI_DEALLOCATE(dtset%xred_orig)
 end if
 if (associated(dtset%ziontypat))   then
   ABI_DEALLOCATE(dtset%ziontypat)
 end if
 if (associated(dtset%znucl))       then
   ABI_DEALLOCATE(dtset%znucl)
 end if

end subroutine dtsetFree
!!***
