!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_copy
!! NAME
!!  m_copy
!!
!! FUNCTION
!!  This module provides a generic interface (deep_copy) used to return a deep copy of pointers.
!!  The procedure is useful if data types with several pointers have to be copied.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  * The intent for pointer arguments is not specified since
!!    we have to conform to the F90 specifications. However xval is IN while copy is OUT
!!
!!  * copy is a pointer and is supposed to be *not allocated*.
!!    If the value to be copied points to null(), also the copy will be nullified.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_copy

 use m_profiling

 use defs_basis

 implicit none

 private 

 public :: deep_copy     ! Performs deep copy of two pointers 

 interface deep_copy
  module procedure deep_copy_int0d
  module procedure deep_copy_int1d
  module procedure deep_copy_int2d
  module procedure deep_copy_int3d
  module procedure deep_copy_int4d
  module procedure deep_copy_rdp0d
  module procedure deep_copy_rdp1d
  module procedure deep_copy_rdp2d
  module procedure deep_copy_rdp3d
  module procedure deep_copy_rdp4d
  module procedure deep_copy_csp0d
  module procedure deep_copy_csp1d
  module procedure deep_copy_csp2d
  module procedure deep_copy_csp3d
  module procedure deep_copy_csp4d
  module procedure deep_copy_cdp0d
  module procedure deep_copy_cdp1d
  module procedure deep_copy_cdp2d
  module procedure deep_copy_cdp3d
  module procedure deep_copy_cdp4d
  module procedure deep_copy_log0d
  module procedure deep_copy_log1d
  module procedure deep_copy_log2d
  module procedure deep_copy_log3d
  module procedure deep_copy_log4d
  !module procedure deep_copy_ch1d  !Does not work on XLF, do not use it for the time being.
 end interface deep_copy


CONTAINS  !===========================================================
!!***

!!****f* m_copy/deep_copy_int0d
!! NAME
!! deep_copy_int0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine deep_copy_int0d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_int0d'
!End of the abilint section

 integer,intent(in) :: xval
 integer,intent(out) :: copy
! *********************************************************************

  copy=xval

end subroutine deep_copy_int0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_int1d
!! NAME
!! deep_copy_int1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_int1d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_int1d'
!End of the abilint section

 integer,pointer :: xval(:)
 integer,pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_ALLOCATE(copy,(il:iu))
  copy(:)=xval(:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_int1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_int2d
!! NAME
!! deep_copy_int2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_int2d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_int2d'
!End of the abilint section

 integer,pointer :: xval(:,:)
 integer,pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_int2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_int3d
!! NAME
!! deep_copy_int3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_int3d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_int3d'
!End of the abilint section

 integer,pointer :: xval(:,:,:)
 integer,pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_int3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_int4d
!! NAME
!! deep_copy_int4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_int4d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_int4d'
!End of the abilint section

 integer,pointer :: xval(:,:,:,:)
 integer,pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_int4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp0d
!! NAME
!! deep_copy_rdp0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_rdp0d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_rdp0d'
!End of the abilint section

 real(dp),intent(in) :: xval
 real(dp),intent(out) :: copy
! *********************************************************************
  copy=xval

end subroutine deep_copy_rdp0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp1d
!! NAME
!! deep_copy_rdp1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_rdp1d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_rdp1d'
!End of the abilint section

 real(dp),pointer :: xval(:)
 real(dp),pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_ALLOCATE(copy,(il:iu))
  copy(:)=xval(:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_rdp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp2d
!! NAME
!! deep_copy_rdp2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_rdp2d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_rdp2d'
!End of the abilint section

 real(dp),pointer :: xval(:,:)
 real(dp),pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_rdp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp3d
!! NAME
!! deep_copy_rdp3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_rdp3d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_rdp3d'
!End of the abilint section

 real(dp),pointer :: xval(:,:,:)
 real(dp),pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_rdp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp4d
!! NAME
!! deep_copy_rdp4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_rdp4d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_rdp4d'
!End of the abilint section

 real(dp),pointer :: xval(:,:,:,:)
 real(dp),pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_rdp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp0d
!! NAME
!! deep_copy_csp0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_csp0d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_csp0d'
!End of the abilint section

 complex(spc),intent(in) :: xval
 complex(spc),intent(out) :: copy
! *********************************************************************
  copy=xval

end subroutine deep_copy_csp0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp1d
!! NAME
!! deep_copy_csp1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_csp1d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_csp1d'
!End of the abilint section

 complex(spc),pointer :: xval(:)
 complex(spc),pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_ALLOCATE(copy,(il:iu))
  copy(:)=xval(:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_csp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp2d
!! NAME
!! deep_copy_csp2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_csp2d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_csp2d'
!End of the abilint section

 complex(spc),pointer :: xval(:,:)
 complex(spc),pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_csp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp3d
!! NAME
!! deep_copy_csp3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_csp3d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_csp3d'
!End of the abilint section

 complex(spc),pointer :: xval(:,:,:)
 complex(spc),pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_csp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp4d
!! NAME
!! deep_copy_csp4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_csp4d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_csp4d'
!End of the abilint section

 complex(spc),pointer :: xval(:,:,:,:)
 complex(spc),pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_csp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp0d
!! NAME
!! deep_copy_cdp0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_cdp0d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_cdp0d'
!End of the abilint section

 complex(dpc),intent(in) :: xval
 complex(dpc),intent(out) :: copy
! *********************************************************************
  copy=xval

end subroutine deep_copy_cdp0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp1d
!! NAME
!! deep_copy_cdp1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_cdp1d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_cdp1d'
!End of the abilint section

 complex(dpc),pointer :: xval(:)
 complex(dpc),pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_ALLOCATE(copy,(il:iu))
  copy(:)=xval(:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_cdp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp2d
!! NAME
!! deep_copy_cdp2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_cdp2d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_cdp2d'
!End of the abilint section

 complex(dpc),pointer :: xval(:,:)
 complex(dpc),pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_cdp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp3d
!! NAME
!! deep_copy_cdp3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_cdp3d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_cdp3d'
!End of the abilint section

 complex(dpc),pointer :: xval(:,:,:)
 complex(dpc),pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_cdp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp4d
!! NAME
!! deep_copy_cdp4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_cdp4d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_cdp4d'
!End of the abilint section

 complex(dpc),pointer :: xval(:,:,:,:)
 complex(dpc),pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_ALLOCATE(copy,(il1:iu1,il2:il2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_cdp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_ch1d
!! NAME
!! deep_copy_ch1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!  This routine segfaults on XLF, disabled for the time being
!!  Should test whether passing slen fixes the problem
!!
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_ch1d(xval,copy,slen)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_ch1d'
!End of the abilint section

 integer,intent(in) :: slen
 character(len=slen),pointer :: xval(:)
 character(len=slen),pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_ALLOCATE(copy,(il:iu))
  copy(:)=xval(:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_ch1d
!!***

!!****f* m_copy/deep_copy_log0d
!! NAME
!! deep_copy_log0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine deep_copy_log0d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_log0d'
!End of the abilint section

 logical,intent(in) :: xval
 logical,intent(out) :: copy
! *********************************************************************

  copy=xval

end subroutine deep_copy_log0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_log1d
!! NAME
!! deep_copy_log1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_log1d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_log1d'
!End of the abilint section

 logical,pointer :: xval(:)
 logical,pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_ALLOCATE(copy,(il:iu))
  copy(:)=xval(:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_log1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_log2d
!! NAME
!! deep_copy_log2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_log2d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_log2d'
!End of the abilint section

 logical,pointer :: xval(:,:)
 logical,pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_log2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_log3d
!! NAME
!! deep_copy_log3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_log3d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_log3d'
!End of the abilint section

 logical,pointer :: xval(:,:,:)
 logical,pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_log3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_log4d
!! NAME
!! deep_copy_log4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine deep_copy_log4d(xval,copy)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'deep_copy_log4d'
!End of the abilint section

 logical,pointer :: xval(:,:,:,:)
 logical,pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_ALLOCATE(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else 
  nullify(copy)
 end if

end subroutine deep_copy_log4d

!----------------------------------------------------------------------

END MODULE m_copy
!!***
