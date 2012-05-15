!!****m* ABINIT/interfaces_27_toolbox_oop
!! NAME
!! interfaces_27_toolbox_oop
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/27_toolbox_oop
!!
!! COPYRIGHT
!! Copyright (C) 2010-2011 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_27_toolbox_oop

 implicit none

interface
 subroutine int2char(iint,string)
  implicit none
  integer,intent(in) :: iint
  character(len=10),intent(out) :: string
 end subroutine int2char
end interface

interface
 subroutine int2char4(iint,string)
  implicit none
  integer,intent(in) :: iint
  character(len=4),intent(out) :: string
 end subroutine int2char4
end interface

end module interfaces_27_toolbox_oop
!!***
