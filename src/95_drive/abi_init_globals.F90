!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_init_globals
!! NAME
!! abi_init_globals
!!
!! FUNCTION
!! This function initializes the global variables defined in the abinit modules.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   None
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      m_header_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine abi_init_globals()

 use m_profiling

 use defs_basis
 use m_errors

 use m_header,  only : m_header_init

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_init_globals'
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
 integer :: ierr

! *************************************************************************

!Init the internal database with the list of fforms associated to the different filetypes.
 call m_header_init(ierr)
 ABI_CHECK(ierr==0," Fatal error in m_header_init")

end subroutine abi_init_globals
!!***
