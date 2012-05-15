!{\src2tex{textfont=tt}}
!!****f* ABINIT/symkchk
!! NAME
!! symkchk
!!
!! FUNCTION
!! Checks that the set of k points chosen for a response function
!! calculation has the full space group symmetry, modulo time reversal
!! if appropirate.
!! Aborts run with error message if not satisfied
!! Currently used only when strain perturbation is treated
!! Based on symkpt.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (DRH, XG, LSI)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! kptns(3,nkpt)= k vectors in reciprocal space
!! nkpt = number of k-points whose weights are wtk
!! nsym=number of space group symmetries
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! timrev: if 1, the time reversal operation has to be taken into account
!! if 0, no time reversal symmetry.
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symkchk(kptns,nkpt,nsym,symrec,timrev)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symkchk'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nkpt,nsym,timrev
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: kptns(3,nkpt)

!Local variables -------------------------
!scalars
 integer :: identi,ii,ikpt,ikpt2,imatch,isym,jj,tident
 real(dp) :: difk,reduce
 character(len=500) :: message
!arrays
 real(dp) :: ksym(3)

! *********************************************************************

 if(timrev/=1 .and. timrev/=0)then
   write(message, '(a,a,a,a,a,i4,a)' )&
&   ' symkchk : BUG -',ch10,&
&   '  timrev should be 0 or 1, while',ch10,&
&   '  it is equal to ',timrev,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(nsym/=1)then
!  Find the identity symmetry operation
   do isym=1,nsym
     tident=1
     do jj=1,3
       if(symrec(jj,jj,isym)/=1)tident=0
       do ii=1,3
         if( ii/=jj .and.&
&         symrec(ii,jj,isym)/=0)tident=0
       end do
     end do
     if(tident==1)then
       identi=isym
       write(message, '(a,i3)' )' symkchk : found identity, with number',identi
       call wrtout(std_out,message,'COLL')
       exit
     end if
   end do
   if(tident==0)then
     write(message, '(a,a,a)' )&
&     ' symkchk : BUG -',ch10,&
&     '  Did not found the identity operation.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

!Here begins the serious business
!The length sorting, etc. of symkpt have been dropped because the
!computational cost is estimated to be negligible.

 if(nsym>1 .or. timrev==1)then

!  Outer loop over kpts
   do ikpt=1,nkpt-1

!    Loop on the symmetries
!    For each k-point and each symmetry transformation, a matching
!    k-pointpt must be found, modulo time reversal if appropriate
     do isym=1,nsym

!      Get the symmetric of the vector
       do ii=1,3
         ksym(ii)= kptns(1,ikpt)*symrec(ii,1,isym)&
&         +kptns(2,ikpt)*symrec(ii,2,isym)&
&         +kptns(3,ikpt)*symrec(ii,3,isym)
       end do

!      Second loop k-points
       do ikpt2=1,nkpt

!        Test for match of symmetric and any vector (including original)
         imatch=1
         do ii=1,3
           difk= ksym(ii)-kptns(ii,ikpt2)
           reduce=difk-anint(difk)
           if(abs(reduce)>tol8)imatch=0
         end do
         if(imatch==1)exit

!        Test for match with time reversal
         if(timrev==1)then
           imatch=1
           do ii=1,3
             difk= ksym(ii)+kptns(ii,ikpt2)
             reduce=difk-anint(difk)
             if(abs(reduce)>tol8)imatch=0
           end do
           if(imatch==1)exit
         end if

!        End secondary loop over k-points
       end do
       if(imatch/=1)then
         write(message, '(a,a,a,a,a,i4,a,i4,a,a,a,a)' )&
&         ' symkchk : ERROR -',ch10,&
&         '   k-point set must have full space-group symmetry',ch10,&
&         '   there is no match for kpt',ikpt,' transformed by symmetry',isym,ch10,&
&         '   Action : change kptopt to 2 or 3 and/or change or use shiftk',ch10,&
&         '            shiftk = 0 0 0 is always a safe choice.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if

!      End loop on isym
     end do

!    End primary loop over k-points
   end do

   write(message, '(a)' )&
&   ' symkchk : k-point set has full space-group symmetry.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end if

end subroutine symkchk
!!***
