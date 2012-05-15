!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_jolly_pointers
!! NAME
!!  m_jolly_pointers
!!
!! FUNCTION
!!  This module defines an extension of the standard F90 pointers whose aim 
!!  is to solve the following ambiguity. F90 pointers can be used in two 
!!  different ways: either to point the data stored in a target or to 
!!  reserve a chunk of the heap. In the later method, a pointer is used in the 
!!  same way as an allocatable variable. There is, however, no general and portable 
!!  way to interrogate a pointer to discern which mode is currently used. 
!!  This indecision might lead to errors at run-time. 
!!  For example deallocating a true pointer might affect the values of the target under gfortran or ifc,  
!!  while sunf90 reports stat/=0.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_jolly_pointers

 use m_profiling

 use defs_basis
 use m_errors

 use defs_datatypes,  only : Jpt_gwpc_2D, Jpt_gwpc_3D

 implicit none

 private

 public :: destroy_jpt
 public :: allocate_jpt
 public :: ispointer
 public :: isallocated

 interface destroy_jpt 
  module procedure destroy_jpt_gwpc_2D
  module procedure destroy_jpt_gwpc_3D
 end interface destroy_jpt

 interface allocate_jpt 
  module procedure allocate_jpt_gwpc_2D
  module procedure allocate_jpt_gwpc_3D
 end interface allocate_jpt

 interface ispointer
  module procedure ispointer_gwpc_2D
  module procedure ispointer_gwpc_3D
 end interface ispointer

 interface isallocated
  module procedure isallocated_gwpc_2D
  module procedure isallocated_gwpc_3D
 end interface isallocated

 integer,private,parameter :: JPT_ISPOINTER  =1 ! The pointer is used to store the address in memory.
 integer,private,parameter :: JPT_ISALLOCATED=2 ! The pointer is used as an allocable array.

CONTAINS  !===========================================================
!!***

!!****f* m_jolly_pointers/destroy_jpt_gwpc_2D
!! NAME
!!  destroy_jpt_gwpc_2D
!!
!! FUNCTION
!!  Deallocate or nullify the pointer defined in a variable of type Jpt_gwpc_2D.
!!
!! INPUTS 
!!  Jpt<Jpt_gwpc_2D>=A variable of type Jpt_gwpc_2D.
!!
!! OUTPUT
!!  istat=Status error.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_jpt_gwpc_2D(Jpt,istat)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_jpt_gwpc_2D'
!End of the abilint section

 type(Jpt_gwpc_2D),intent(inout) :: Jpt
 integer,intent(out) :: istat

! *************************************************************************

 istat=0

 select case (Jpt%stat)

 case (JPT_ISPOINTER)
  nullify(Jpt%datum)

 case (JPT_ISALLOCATED)
  ABI_DEALLOCATE(Jpt%datum)
  istat = ABI_ALLOC_STAT
  nullify(Jpt%datum)
  Jpt%stat  = JPT_ISPOINTER

 case DEFAULT
  MSG_BUG("Wrong value for %stat")
 end select

end subroutine destroy_jpt_gwpc_2D
!!***

!!****f* m_jolly_pointers/destroy_jpt_gwpc_3D
!! NAME
!!  destroy_jpt_gwpc_3D
!!
!! FUNCTION
!!  Deallocate or nullify the pointer defined in a Jpt_gwpc_3D data type.
!!
!! INPUTS 
!!  Jpt<Jpt_gwpc_3D>=A variable of type Jpt_gwpc_3D.
!!
!! OUTPUT
!!  istat=Status error.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_jpt_gwpc_3D(Jpt,istat)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_jpt_gwpc_3D'
!End of the abilint section

 type(Jpt_gwpc_3D),intent(inout) :: Jpt
 integer,intent(out) :: istat

! *************************************************************************

 istat=0

 select case (Jpt%stat)

 case (JPT_ISPOINTER)
  nullify(Jpt%datum)

 case (JPT_ISALLOCATED)
  ABI_DEALLOCATE(Jpt%datum)
  istat = ABI_ALLOC_STAT
  nullify(Jpt%datum)
  Jpt%stat  = JPT_ISPOINTER

 case DEFAULT
  MSG_BUG("Wrong value for %stat")
 end select

end subroutine destroy_jpt_gwpc_3D
!!***

!!****f* m_jolly_pointers/allocate_jpt_gwpc_2D
!! NAME
!!  allocate_jpt_gwpc_2D
!!
!! FUNCTION
!!  Allocate the pointer defined in the Jpt_gwpc_2D data type changing its status. 
!!  Performs also additional checking to make sure the pointer is not already allocated.
!!
!! INPUTS 
!!  Jpt<Jpt_gwpc_2D>=A variable of type Jpt_gwpc_2D.
!!  dt_shape(2,2)=Shape used to allocate the pointer.
!!
!! OUTPUT
!!  istat=Status error.
!!
!! SIDE EFFECTS
!!  The %datum pointed defined in the data type is allocated.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine allocate_jpt_gwpc_2D(Jpt,dt_shape,istat)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'allocate_jpt_gwpc_2D'
!End of the abilint section

 type(Jpt_gwpc_2D),intent(inout) :: Jpt
 integer,intent(in) :: dt_shape(2,2)
 integer,intent(out) :: istat

! *************************************************************************

 if (Jpt%stat==JPT_ISALLOCATED) then
  istat=-1
  MSG_ERROR("%datum should not be allocated")
  !?call destroy_jpt(Jpt,istat)
  !?if (istat/=0) return
 end if

 select case (Jpt%stat)

 case (JPT_ISPOINTER) 
  ABI_ALLOCATE(Jpt%datum ,(dt_shape(1,1):dt_shape(2,1),dt_shape(1,2):dt_shape(2,2) ))
  istat = ABI_ALLOC_STAT
  Jpt%stat=JPT_ISALLOCATED

 case DEFAULT
  MSG_BUG("%datum should be a pointer")
 end select

end subroutine allocate_jpt_gwpc_2D
!!***

!!****f* m_jolly_pointers/allocate_jpt_gwpc_3D
!! NAME
!!  allocate_jpt_gwpc_3D
!!
!! FUNCTION
!!  Allocate the pointer defined in the Jpt_gwpc_3D data type changing its status. 
!!  Performs also additional checking to make sure the pointer is not already allocated.
!!
!! INPUTS 
!!  Jpt<Jpt_gwpc_3D>=A variable of type Jpt_gwpc_3D.
!!  dt_shape(2,3)=Shape used to allocate the pointer.
!!
!! OUTPUT
!!  istat=Status error.
!!
!! SIDE EFFECTS
!!  The %datum pointed defined in the data type is allocated.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine allocate_jpt_gwpc_3D(Jpt,dt_shape,istat)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'allocate_jpt_gwpc_3D'
!End of the abilint section

 type(Jpt_gwpc_3D),intent(inout) :: Jpt
 integer,intent(in) :: dt_shape(2,3)
 integer,intent(out) :: istat

! *************************************************************************

 if (Jpt%stat==JPT_ISALLOCATED) then
  istat=-1
  MSG_ERROR("%datum should not be allocated")
  !?call destroy_jpt(Jpt,istat)
  !?if (istat/=0) return
 end if

 select case (Jpt%stat)

 case (JPT_ISPOINTER) 
  ABI_ALLOCATE(Jpt%datum ,(dt_shape(1,1):dt_shape(2,1),dt_shape(1,2):dt_shape(2,2),dt_shape(1,3):dt_shape(2,3) ))
  istat = ABI_ALLOC_STAT
  Jpt%stat=JPT_ISALLOCATED

 case DEFAULT
  MSG_BUG("%datum should be a pointer")
 end select

end subroutine allocate_jpt_gwpc_3D
!!***

!!****f* m_jolly_pointers/ispointer_gwpc_2D
!! NAME
!!  ispointer_gwpc_2D
!!
!! FUNCTION
!!  Returns .TRUE. if the object of type Jpt_gwpc_2D is used as a true pointer.
!!
!! INPUTS 
!!  Jpt<Jpt_gwpc_2D>=A variable of type Jpt_gwpc_2D.
!!
!! OUTPUT
!!  ans=TRUE if %datum is used as pointer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ispointer_gwpc_2D(Jpt) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ispointer_gwpc_2D'
!End of the abilint section

 type(Jpt_gwpc_2D),intent(in) :: Jpt
 logical :: ans

! *************************************************************************
 ans = (Jpt%stat==JPT_ISPOINTER)

end function ispointer_gwpc_2D
!!***

!!****f* m_jolly_pointers/ispointer_gwpc_3D
!! NAME
!!  ispointer_gwpc_3D
!!
!! FUNCTION
!!  Returns .TRUE. if the object of type Jpt_gwpc_3D is used as a true pointer.
!!
!! INPUTS 
!!  Jpt<Jpt_gwpc_3D>=A variable of type Jpt_gwpc_3D.
!!
!! OUTPUT
!!  ans=TRUE if %datum is used as pointer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ispointer_gwpc_3D(Jpt) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ispointer_gwpc_3D'
!End of the abilint section

 type(Jpt_gwpc_3D),intent(in) :: Jpt
 logical :: ans

! *************************************************************************
 ans = (Jpt%stat==JPT_ISPOINTER)

end function ispointer_gwpc_3D
!!***

!!****f* m_jolly_pointers/isallocated_gwpc_2D
!! NAME
!!  isallocated_gwpc_2D
!!
!! FUNCTION
!!  Returns .TRUE. if the object of type Jpt_gwpc_2D is used as an allocatable variable. 
!!
!! INPUTS 
!!  Jpt<Jpt_gwpc_2D>=A variable of type Jpt_gwpc_2D.
!!
!! OUTPUT
!!  ans=TRUE if %datum is used as pointer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function isallocated_gwpc_2D(Jpt) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isallocated_gwpc_2D'
!End of the abilint section

 type(Jpt_gwpc_2D),intent(in) :: Jpt
 logical :: ans

! *************************************************************************
 ans = (Jpt%stat==JPT_ISALLOCATED)

end function isallocated_gwpc_2D 
!!***

!!****f* m_jolly_pointers/isallocated_gwpc_3D
!! NAME
!!  isallocated_gwpc_3D
!!
!! FUNCTION
!!  Returns .TRUE. if the object of type Jpt_gwpc_3D is used as an allocatable variable. 
!!
!! INPUTS 
!!  Jpt<Jpt_gwpc_3D>=A variable of type Jpt_gwpc_3D.
!!
!! OUTPUT
!!  ans=TRUE if %datum is used as pointer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function isallocated_gwpc_3D(Jpt) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isallocated_gwpc_3D'
!End of the abilint section

 type(Jpt_gwpc_3D),intent(in) :: Jpt
 logical :: ans

! *************************************************************************
 ans = (Jpt%stat==JPT_ISALLOCATED)

end function isallocated_gwpc_3D 

END MODULE m_jolly_pointers
!!***
