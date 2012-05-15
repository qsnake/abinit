!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkph3
!! NAME
!! chkph3
!!
!! FUNCTION
!! Check the completeness of the dynamical matrix
!! and eventually send a warning
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  idir = direction of the eventual electric field
!!  mpert =maximum number of ipert
!!  natom=number of atoms in unit cell
!!
!! OUTPUT
!!  eventually send a warning message
!!
!! NOTES
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkph3(carflg,idir,mpert,natom)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkph3'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: idir,mpert,natom
!arrays
 integer,intent(in) :: carflg(3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ipert1,ipert2,send
 character(len=500) :: message

! *********************************************************************

 send=0

!Check the elements of the analytical part of the
!dynamical matrix
 do ipert2=1,natom
   do idir2=1,3
     do ipert1=1,natom
       do idir1=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==0)then
           send=1
         end if
       end do
     end do
   end do
 end do

!If some electric field is present
 if(idir/=0)then

!  Check the dielectric constant
   if(carflg(idir,natom+2,idir,natom+2)==0)then
     send=1
   end if

!  Check the effective charges
   do ipert1=1,natom
     do idir1=1,3
       if(carflg(idir1,ipert1,idir,natom+2)==0)then
         send=1
       end if
     end do
   end do

 end if

!If needed, send the message
 if(send==1)then
   write(message, '(a,a,a,a)' )&
&   ' chkph3 : WARNING -',ch10,&
&   '  The dynamical matrix was incomplete :',&
&   ' phonon frequencies may be wrong ...'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end if

end subroutine chkph3
!!***
