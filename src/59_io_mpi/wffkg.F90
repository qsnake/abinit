!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffKg
!! NAME
!! WffKg
!!
!! FUNCTION
!! Check kgwff to  manage WF file in the MPI/IO case
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  wff <type(wffile_type)> = structured info about the wavefunction file
!!  optkg= if 1 , read or write kg_k ; if 0,do not care about kg_k in rwwf
!!
!! OUTPUT
!!
!! PARENTS
!!      inwffil,kss2wfk,m_wfs,nstdy3,nstpaw3,outwf,tddft,vtorho,vtorho3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine WffKg(wff,optkg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_wffile
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffKg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: optkg

! *********************************************************************

#if defined HAVE_MPI_IO
 if (wff%accesswff == IO_MODE_MPI) wff%kgwff=optkg
#else
 ABI_UNUSED((/wff%accesswff,optkg/))
#endif

end subroutine WffKg
!!***
