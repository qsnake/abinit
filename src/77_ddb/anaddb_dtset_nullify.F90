!{\src2tex{textfont=tt}}
!!****f* ABINIT/anaddb_dtset_nullify
!!
!! NAME
!!   anaddb_dtset_nullify
!!
!! FUNCTION
!!   nullify all arrays in the anaddb_dtset datastructure
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


subroutine anaddb_dtset_nullify(anaddb_dtset)

 use m_profiling

 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anaddb_dtset_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(anaddb_dataset_type), intent(inout) :: anaddb_dtset

! *************************************************************************
 
 nullify(anaddb_dtset%atifc)
 nullify(anaddb_dtset%iatfix)
 nullify(anaddb_dtset%qnrml1)
 nullify(anaddb_dtset%qnrml2)
 nullify(anaddb_dtset%qpath)
 nullify(anaddb_dtset%qph1l)
 nullify(anaddb_dtset%qph2l)
 nullify(anaddb_dtset%ep_qptlist)

end subroutine anaddb_dtset_nullify
!!***
