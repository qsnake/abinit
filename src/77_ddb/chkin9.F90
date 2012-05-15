!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkin9
!!
!! NAME
!! chkin9
!!
!! FUNCTION
!! Check the value of some input parameters.
!! Send error message and stop if needed.
!! Also transform the meaning of atifc
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! atifc(natom)=list of the atom ifc to be analysed
!! natifc= number of atom ifc to be analysed
!! natom= number of atoms
!!
!! OUTPUT
!! atifc(natom) =  atifc(ia) equals 1 if the analysis of ifc
!!  has to be done for atom ia; otherwise 0.
!!
!! NOTES
!! Only for one processor (no use of wrtout)
!!
!! PARENTS
!!      rdddb9,thmeig
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkin9(atifc,natifc,natom)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkin9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natifc,natom
!arrays
 integer,intent(inout) :: atifc(natom)

!Local variables -------------------------
!scalars
 integer :: iatifc
 character(len=500) :: message
!arrays
 integer,allocatable :: work(:)

! *********************************************************************

 if(natifc>natom)then
   write(message, '(a,a,a,i6,a,a,a,i6,a,a,a)' )&
&   ' chkin9 : ERROR -',ch10,&
&   '  The number of atom ifc in the input files',natifc,',',ch10,&
&   '  is larger than the number of atoms',natom,'.',ch10,&
&   '  Action : change natifc in the input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(natifc>=1)then
   ABI_ALLOCATE(work,(natom))
   work(:)=0

   do iatifc=1,natifc
     if(atifc(iatifc)<=0.or.atifc(iatifc)>natom)then
       write(message, '(a,a,a,i6,a,a,a,a,a,i8,a,a,a)' )&
&       ' chkin9 : ERROR-',ch10,&
&       '  For iatifc=',iatifc,', the number of the atom ifc to be ',ch10,&
&       '  analysed is not valid : either negative, ',ch10,&
&       '  zero, or larger than natom =',natom,'.',ch10,&
&       '  Action : change atifc in your input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     work(atifc(iatifc))=1
   end do

   atifc(1:natom)=work(:)
   ABI_DEALLOCATE(work)
 end if

end subroutine chkin9
!!***
