!{\src2tex{textfont=tt}}
!!****f* ABINIT/anaddb_dtset_clean
!!
!! NAME
!!   anaddb_dtset_clean
!!
!! FUNCTION
!!   deallocate remaining arrays in the anaddb_dtset datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  anaddb_dtset = anaddb datastructure
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
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


subroutine anaddb_dtset_clean(anaddb_dtset)

 use m_profiling

 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anaddb_dtset_clean'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(anaddb_dataset_type), intent(inout) :: anaddb_dtset

! *************************************************************************
 
 if (associated(anaddb_dtset%atifc))  then
   ABI_DEALLOCATE(anaddb_dtset%atifc)
 end if
 if (associated(anaddb_dtset%iatfix))  then
   ABI_DEALLOCATE(anaddb_dtset%iatfix)
 end if
 if (associated(anaddb_dtset%iatprj_bs))  then
   ABI_DEALLOCATE(anaddb_dtset%iatprj_bs)
 end if
 if (associated(anaddb_dtset%qnrml1))  then
   ABI_DEALLOCATE(anaddb_dtset%qnrml1)
 end if
 if (associated(anaddb_dtset%qnrml2))  then
   ABI_DEALLOCATE(anaddb_dtset%qnrml2)
 end if
 if (associated(anaddb_dtset%qpath))  then
   ABI_DEALLOCATE(anaddb_dtset%qpath)
 end if
 if (associated(anaddb_dtset%qph1l))  then
   ABI_DEALLOCATE(anaddb_dtset%qph1l)
 end if
 if (associated(anaddb_dtset%qph2l))  then
   ABI_DEALLOCATE(anaddb_dtset%qph2l)
 end if
 if (associated(anaddb_dtset%ep_qptlist))  then
   ABI_DEALLOCATE(anaddb_dtset%ep_qptlist)
 end if

end subroutine anaddb_dtset_clean
!!***
