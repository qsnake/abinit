!{\src2tex{textfont=tt}}
!!****f* ABINIT/filnam_comm
!! NAME
!! filnam_comm
!!
!! FUNCTION
!! Communicate filenames to all processors
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  nfilnam=size of filnam array (number of file names)
!!
!! SIDE EFFECTS
!!  character(len=fnlen) :: filnam(5)=character strings giving file names
!!
!! PARENTS
!!      iofn1
!!
!! CHILDREN
!!      timab,xbarrier_mpi,xcast_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine filnam_comm(nfilnam,filnam)

 use m_profiling

 use defs_basis
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'filnam_comm'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: nfilnam
 character(len=fnlen), intent(inout) :: filnam(nfilnam)

!Local variables-------------------------------
 integer :: ierr
 real(dp) :: tsec(2)

!*************************************************************************

 call xbarrier_mpi(xmpi_world)
 call timab(48,1,tsec)
 call xcast_mpi(filnam,0,xmpi_world,ierr)
 call timab(48,2,tsec)

end subroutine filnam_comm
!!***
