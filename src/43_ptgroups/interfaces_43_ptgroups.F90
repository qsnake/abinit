!!****m* ABINIT/interfaces_43_ptgroups
!! NAME
!! interfaces_43_ptgroups
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/43_ptgroups
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

module interfaces_43_ptgroups

 implicit none

interface
 subroutine ptg_C1 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C1
end interface

interface
 subroutine ptg_C2 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C2
end interface

interface
 subroutine ptg_C2h (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C2h
end interface

interface
 subroutine ptg_C2v (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C2v
end interface

interface
 subroutine ptg_C3 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C3
end interface

interface
 subroutine ptg_C3h (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C3h
end interface

interface
 subroutine ptg_C3i (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C3i
end interface

interface
 subroutine ptg_C3v (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C3v
end interface

interface
 subroutine ptg_C4 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C4
end interface

interface
 subroutine ptg_C4h (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C4h
end interface

interface
 subroutine ptg_C4v (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C4v
end interface

interface
 subroutine ptg_C6 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C6
end interface

interface
 subroutine ptg_C6h (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C6h
end interface

interface
 subroutine ptg_C6v (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_C6v
end interface

interface
 subroutine ptg_Ci (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_Ci
end interface

interface
 subroutine ptg_Cs (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_Cs
end interface

interface
 subroutine ptg_D2 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D2
end interface

interface
 subroutine ptg_D2d (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D2d
end interface

interface
 subroutine ptg_D2h (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D2h
end interface

interface
 subroutine ptg_D3 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D3
end interface

interface
 subroutine ptg_D3d (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D3d
end interface

interface
 subroutine ptg_D3h (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D3h
end interface

interface
 subroutine ptg_D4 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D4
end interface

interface
 subroutine ptg_D4h (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D4h
end interface

interface
 subroutine ptg_D6 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D6
end interface

interface
 subroutine ptg_D6h (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_D6h
end interface

interface
 subroutine ptg_O (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_O
end interface

interface
 subroutine ptg_Oh (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_Oh
end interface

interface
 subroutine ptg_S4 (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_S4
end interface

interface
 subroutine ptg_T (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_T
end interface

interface
 subroutine ptg_Td (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_Td
end interface

interface
 subroutine ptg_Th (nsym,nclass,sym,class_ids,class_names,Irr)
  use m_defs_ptgroups
  implicit none
  integer,intent(out) :: nclass
  integer,intent(out) :: nsym
  integer,pointer :: class_ids(:,:)
  character(len=5),pointer :: class_names(:)
  integer,pointer :: sym(:,:,:)
  type(irrep_t),pointer :: Irr(:)
 end subroutine ptg_Th
end interface

end module interfaces_43_ptgroups
!!***
