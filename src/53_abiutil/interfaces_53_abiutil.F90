!!****m* ABINIT/interfaces_53_abiutil
!! NAME
!! interfaces_53_abiutil
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/53_abiutil
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

module interfaces_53_abiutil

 implicit none

interface
 subroutine dtset_nullify(dtset)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(inout) :: dtset
 end subroutine dtset_nullify
end interface

interface
 subroutine dtsetCopy(dtout, dtin)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtin
  type(dataset_type),intent(out) :: dtout
 end subroutine dtsetCopy
end interface

interface
 subroutine  tells_sizes(chosen_size,name,index,default_size,actual_size)
  implicit none
  integer,intent(in) :: actual_size
  integer,intent(out) :: chosen_size
  integer,intent(in) :: default_size
  integer,intent(in) :: index
  character(len=12),intent(in) :: name
 end subroutine tells_sizes
end interface

interface
 subroutine dtsetFree(dtset)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(inout) :: dtset
 end subroutine dtsetFree
end interface

interface
 subroutine find_getdtset(dtsets,getvalue,getname,idtset,iget,miximage,mxnimage,ndtset_alloc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: getvalue
  integer, intent(in) :: idtset
  integer, intent(out) :: iget
  integer, intent(in) :: mxnimage
  integer, intent(in) :: ndtset_alloc
  character(len=*),intent(in) :: getname
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  real(dp), intent(out) :: miximage(mxnimage,mxnimage)
 end subroutine find_getdtset
end interface

interface
 subroutine mkfilename(filnam,filnam_out,get,idtset,&  
  &  ird,jdtset_,ndtset,stringfil,stringvar,will_read)
  use defs_basis
  implicit none
  integer,intent(in) :: get
  integer,intent(in) :: idtset
  integer,intent(in) :: ird
  integer,intent(in) :: ndtset
  integer,intent(out) :: will_read
  character(len=fnlen),intent(out) :: filnam_out
  character(len=*),intent(in) :: stringfil
  character(len=*),intent(in) :: stringvar
  character(len=fnlen),intent(in) :: filnam(5)
  integer,intent(in) :: jdtset_(0:ndtset)
 end subroutine mkfilename
end interface

interface
 subroutine orthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: vectsize
  type(mpi_type) :: mpi_enreg
  real(dp) :: blockvectorbx(vectsize,blocksize)
  real(dp) :: blockvectorx(vectsize,blocksize)
  real(dp) :: sqgram(blocksize,blocksize)
 end subroutine orthonormalize
end interface

interface
 subroutine zorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: vectsize
  type(mpi_type) :: mpi_enreg
  complex(dpc) :: blockvectorbx(vectsize,blocksize)
  complex(dpc) :: blockvectorx(vectsize,blocksize)
  complex(dpc) :: sqgram(blocksize,blocksize)
 end subroutine zorthonormalize
end interface

end module interfaces_53_abiutil
!!***
