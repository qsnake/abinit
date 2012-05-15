!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_iso_c_binding
!! NAME
!! m_iso_c_binding
!!
!! FUNCTION
!! This module provides a (very) partial implementation of the intrinsic F2003 module iso_c_bindings 
!! for compilers that do not provide a native implementation. 
!! At present it contains the size of basic C types that are used in m_clib to declare Fortran structures 
!! and Fortran interfaces for the bindings.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE m_iso_c_binding

 implicit none

 private 

! character type
 integer(kind=4),public,parameter :: C_CHAR = SIZEOF_CHAR

! integer types
 integer(kind=4),public,parameter :: C_SHORT     = SIZEOF_SHORT
 integer(kind=4),public,parameter :: C_INT       = SIZEOF_INT
 integer(kind=4),public,parameter :: C_LONG      = SIZEOF_LONG
 integer(kind=4),public,parameter :: C_LONG_LONG = SIZEOF_LONG_LONG

! real types
 integer(kind=4),public,parameter :: C_FLOAT       = SIZEOF_FLOAT
 integer(kind=4),public,parameter :: C_DOUBLE      = SIZEOF_DOUBLE
 integer(kind=4),public,parameter :: C_LONG_DOUBLE = SIZEOF_LONG_DOUBLE
      
! special characters
 character(kind=1,len=1),public,parameter :: C_NULL_CHAR       = ACHAR(0)
 character(kind=1,len=1),public,parameter :: C_ALERT           = ACHAR(7)
 character(kind=1,len=1),public,parameter :: C_BACKSPACE       = ACHAR(8)
 character(kind=1,len=1),public,parameter :: C_FORM_FEED       = ACHAR(12)
 character(kind=1,len=1),public,parameter :: C_NEW_LINE        = ACHAR(10)
 character(kind=1,len=1),public,parameter :: C_CARRIAGE_RETURN = ACHAR(13)
 character(kind=1,len=1),public,parameter :: C_HORIZONTAL_TAB  = ACHAR(9)
 character(kind=1,len=1),public,parameter :: C_VERTICAL_TAB    = ACHAR(11)

END MODULE m_iso_c_binding
!!***
