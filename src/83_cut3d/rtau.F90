!{\src2tex{textfont=tt}}
!!****f* ABINIT/rtau
!! NAME
!! rtau
!!
!! FUNCTION
!! Reads in the atomic positions in xmol format
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! filnam=string containing the filename
!! nat=integer number of atoms
!! ntypat=integer number of atom types
!!
!! OUTPUT
!! tau(3,nat)=atomic positions in 3D cartesian space (from XMOL format)
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rtau(filnam,tau,nat,ntypat)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rtau'
!End of the abilint section

 implicit none

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: nat,ntypat
 character(len=fnlen),intent(in) :: filnam
!arrays
 real(dp),intent(out) :: tau(3,nat)

!Local variables--------------------------------------------------------
!scalars
 integer :: iat,idum,itypat,mtypat
 character(len=2) :: type
!arrays
 character(len=2),allocatable :: atypes(:)

! *************************************************************************

 ABI_ALLOCATE(atypes,(ntypat))
 mtypat = 0                ! start with empty array
 write(std_out,*)
 write(std_out,*) 'PROCESSING POSITION FILE ',trim(filnam)
!
 open(unit=19,file=trim(filnam),form='formatted',status='old')

 read(unit=19,fmt=*,end=919) idum
 read(unit=19,fmt=*,end=919)
 if(idum.ne.nat) then
   write(std_out,*) 'Mismatch between the number of atoms in files ',&
&   ' cut.in and ',trim(filnam)
   stop
 end if
 readlp: do iat=1,nat
   read(unit=19,fmt=*,end=919) type,(tau(idum,iat),idum=1,3)
   if (mtypat == 0) then
     mtypat = 1                ! initialize used entries counter
     atypes(1) = type        ! fill first entry
   else
     do itypat=1,mtypat
       if (type == atypes(itypat)) cycle readlp        ! symbol already seen
     end do
!    Here if atomic symbol "type" is met for the first time
     if (mtypat >= ntypat) then
       write(std_out,*) 'Error, more than ntypat=',ntypat,' atom types encountered'
       stop
     end if
     mtypat = mtypat+1                ! bump counter of used entries
     atypes(mtypat) = type        ! save into next available
   end if
 end do readlp
 close(19)

 ABI_DEALLOCATE(atypes)
 return
 919  write(std_out,*) 'Error, premature end of file encountered on ',trim(filnam)
 stop
end subroutine rtau
!!***
