!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_clib
!! NAME
!! m_clib
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_clib

 !use m_defs_basis
 !use m_errors

#ifdef HAVE_FC_ISO_C_BINDING
#define USE_MODULE use iso_c_binding
#else
#define USE_MODULE use m_iso_c_binding
#endif

 USE_MODULE

 implicit none

 private 

 type,public :: Mallinfo_t
   integer(C_LONG) :: arena 
   integer(C_LONG) :: hblkhd 
   integer(C_LONG) :: usmblks 
   integer(C_LONG) :: fsmblks 
   integer(C_LONG) :: uordblks 
   integer(C_LONG) :: fordblks
 end type Mallinfo_t

!FIXME the interfaces below have been commented out since abilint
! crashes during the analysis of the file (maybe due to the macro USE_MODULE!)      

! ===================================================
! ==== Fortran-bindings declared in intrinsics.c ====
! ===================================================
! interface 
!   subroutine clib_fflush()
!     !USE_MODULE
!     implicit none
!   end subroutine clib_fflush
! end interface
!
! interface 
!   subroutine clib_getenv(ierr,fname)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_getenv
! end interface
!
!! ===================================================
!! ==== Fortran-bindings declared in fsi_posix.c ====
!! ===================================================
! interface 
!   subroutine clib_mkdir(ierr,fname)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_mkdir
! end interface
!
! interface 
!   subroutine clib_chdir(ierr,fname)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_chdir
! end interface
!
! interface 
!   subroutine clib_rename(ierr,from_fname,to_fname)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: from_fname,to_fname
!   end subroutine clib_rename
! end interface
!
! interface 
!   subroutine clib_remove(ierr,fname)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_remove
! end interface
!
! interface 
!   subroutine clib_getcwd(ierr,fname)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_getcwd
! end interface
!
! interface 
!   subroutine clib_gethname(ierr,fname)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_gethname
! end interface
!
!! =====================================================
!! ==== Fortran-bindings declared in progress_bar.c ====
!! =====================================================
! interface
!   subroutine clib_progress_bar(actual, max)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(in) :: actual
!     integer(C_INT),intent(in) :: max
!   end subroutine clib_progress_bar
! end interface
!
!! =================================================
!! ==== Fortran-bindings declared in mallinfo.c ====
!! =================================================
! interface
!   subroutine clib_mallinfo(arena, hblkhd, usmblks, fsmblks, uordblks, fordblks)
!     !USE_MODULE
!     implicit none
!     integer(C_LONG),intent(out) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
!   end subroutine clib_mallinfo
! end interface
!
!
!! ==================================================
!! ==== Fortran-bindings declared in gnu_tools.c ====
!! ==================================================
!
! interface
!   subroutine clib_mtrace(ierr)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_mtrace
! end interface
!
! interface
!   subroutine clib_muntrace(ierr)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_muntrace
! end interface
!
! interface
!   subroutine clib_mcheck(ierr)
!     !USE_MODULE
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_mcheck
! end interface

CONTAINS  !===========================================================


subroutine fmallinfo(Minfo)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fmallinfo'
!End of the abilint section

 type(Mallinfo_t),intent(out) :: Minfo

!Local variables-------------------------------
 integer(C_LONG) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
! *********************************************************************

  call clib_mallinfo(arena,hblkhd,usmblks,fsmblks,uordblks,fordblks) 

  Minfo%arena    = arena
  Minfo%hblkhd   = hblkhd
  Minfo%usmblks  = usmblks
  Minfo%fsmblks  = fsmblks
  Minfo%uordblks = uordblks
  Minfo%fordblks = fordblks

end subroutine fmallinfo 
!!***

subroutine print_mallinfo(Minfo,unit)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_mallinfo'
!End of the abilint section

 integer,intent(in) :: unit
 type(Mallinfo_t),intent(in) :: Minfo
! *********************************************************************

 write(unit,*)' Total space in arena            : ',Minfo%arena
 write(unit,*)' Space in holding block headers  : ',Minfo%hblkhd
 write(unit,*)' Space in small blocks in use    : ',Minfo%usmblks
 write(unit,*)' Space in free small blocks      : ',Minfo%fsmblks
 write(unit,*)' Space in ordinary blocks in use : ',Minfo%uordblks
 write(unit,*)' Space in free ordinary blocks   : ',Minfo%fordblks
 write(unit,*)' End memory statistics '

end subroutine print_mallinfo
!!***

END MODULE m_clib
!!***
