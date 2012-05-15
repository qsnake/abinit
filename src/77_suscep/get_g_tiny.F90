!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_g_tiny
!! NAME
!! get_g_tiny
!!
!! FUNCTION
!! Compute the squares of the G vectors, index them with increasing length, index
!! the npw_tiny shortest G vectors, overall and along the primitive translations.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (MF, XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric (bohr**-2)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr**-1)
!!  kg(3,npwdiel)=reciprocal space translations index array (integer)
!!  npw=number of plane waves
!!  npw_tiny=number of (different) shortest vectors to find
!!
!! OUTPUT
!!  gsq_unsorted(npw)=array with dimensional squares of the G vectors (bohr**-2)
!!  ig_tiny(npw_tiny,3)=index of the shortest G vectors along each primitive translation
!!  igsq_tiny(npw_tiny)=index of the shortest G vectors regardless of direction
!!  index_g(npw)=sorting index for G vectors, with increasing length
!!
!! TODO
!! Check possible redundancies with other Abinit routines.
!!
!! PARENTS
!!      xcacfd
!!
!! CHILDREN
!!      leave_new,sort_dp,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine get_g_tiny(gmet,gprimd,gsq_unsorted,ig_tiny,igsq_tiny,index_g,kg,npw,npw_tiny)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_g_tiny'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments-------------------------------------
!scalars
 integer,intent(in) :: npw,npw_tiny
!arrays
 integer,intent(in) :: kg(3,npw)
 integer,intent(out) :: ig_tiny(npw_tiny,3),igsq_tiny(npw_tiny),index_g(npw)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(out) :: gsq_unsorted(npw)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ipw1
 real(dp),parameter :: diffgsq=1.d-2
 real(dp) :: gred1,gred2,gred3,tpisq
 logical,save :: tfirst=.true.
 character(len=500) :: message
!arrays
 real(dp),allocatable :: gsq(:)

!DEBUG
!integer :: i3
!ENDDEBUG

! *************************************************************************
!DEBUG
!write(std_out,*) '%enter: get_g_tiny'
!ENDDEBUG

!tpisq is (2 Pi) **2:
 tpisq=(two_pi)**2

!Perform allocations
 ABI_ALLOCATE(gsq,(npw))

 if(tfirst) then
   write(std_out,*) '%get_g_tiny: primitive g vectors are:'
   do i1=1,3
     write(std_out,'(i4,3(1x,es15.8))') i1,gprimd(1:3,i1)
   end do
 end if

!Find the shortest vectors along the primitive translations
 ig_tiny(:,:)=0
 do i1=1,npw_tiny
   do ipw1=1,npw
     if(kg(1,ipw1)==i1 .and. kg(2,ipw1)==0 .and. kg(3,ipw1)==0) ig_tiny(i1,1)=ipw1
     if(kg(2,ipw1)==i1 .and. kg(3,ipw1)==0 .and. kg(1,ipw1)==0) ig_tiny(i1,2)=ipw1
     if(kg(3,ipw1)==i1 .and. kg(1,ipw1)==0 .and. kg(2,ipw1)==0) ig_tiny(i1,3)=ipw1
   end do
   if(ig_tiny(i1,1)*ig_tiny(i1,2)*ig_tiny(i1,3)==0) then
     write(message, '(6a,i4,2a,3i4,4a)') ch10,&
&     ' get_g_tiny : ERROR -',ch10,&
&     '  Along one or more primitive translations, could not find ',ch10,&
&     '  the', i1,'-th shortest translation. Present reduced components are:',ch10,&
&     ig_tiny(i1,:), ' (the unfound translations have a 0 entry).',ch10,&
&     '  Action : Increase plane wave cutoff diecut.',ch10
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end do

!DEBUG
!write(std_out,*) '%get_g_tiny: show all primitive translations for all g vectors:'
!do ipw1=1,npw
!i1=kg(1,ipw1); i2=kg(2,ipw1); i3=kg(3,ipw1)
!write(std_out,'(a,1(i4,1x),f12.6)')     ipw1,gsq(ipw1)
!write(std_out,'(1(i4,1x),3(f12.6,1x))') i1,i1*gprimd(1,1),i1*gprimd(2,1),i1*gprimd(3,1)
!write(std_out,'(1(i4,1x),3(f12.6,1x))') i2,i2*gprimd(1,2),i2*gprimd(2,2),i2*gprimd(3,2)
!write(std_out,'(1(i4,1x),3(f12.6,1x))') i3,i3*gprimd(1,3),i3*gprimd(2,3),i3*gprimd(3,3)
!end do
!ENDDEBUG

!Compute the squares of g vectors
 do ipw1=1,npw
   index_g(ipw1)=ipw1
   gred1=dble(kg(1,ipw1))
   gred2=dble(kg(2,ipw1))
   gred3=dble(kg(3,ipw1))
   gsq(ipw1)=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&   +2.0_dp*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&   gmet(2,3)*gred2*gred3) )
 end do

!Sort the G vectors by increasing squares
 gsq_unsorted(:)=gsq(:)
 call sort_dp(npw,gsq,index_g,tol14)

!Find the shortest vectors
 igsq_tiny(1)=2
 do i1=2,npw_tiny
   do ipw1=igsq_tiny(i1-1),npw
     if(abs(gsq(ipw1)-gsq(igsq_tiny(i1-1))) > diffgsq*abs(gsq(igsq_tiny(i1-1)))) then

!      DEBUG
!      write(std_out,*) i1,ipw1,abs(gsq(ipw1)-gsq(igsq_tiny(i1-1))),diffgsq*abs(gsq(igsq_tiny(i1-1)))
!      ENDDEBUG

       exit
     end if
   end do
   igsq_tiny(i1)=max(2,ipw1)
 end do

!Print out what was found
 if(tfirst) then
   write(std_out,'(a)') '%get_g_tiny: show index, sorted index, square of g vectors:'

!  DEBUG
!  do ipw1=1,npw
!  write(std_out,*) ipw1,index_g(ipw1),gsq(ipw1),gsq_unsorted(ipw1)
!  end do
!  ENDDEBUG

   write(std_out,'(a,i4,a)') '%get_g_tiny: index for first npw_tiny=', npw_tiny,&
&   ' along primitive translations:'
   do i1=1,3
     write(std_out,'(1x,a,i4,a,3(1x,i4))') 'xyz=',i1,' ig_tiny=',ig_tiny(:,i1)
     do i2=1,npw_tiny
       write(std_out,'(3x,i4,1x,i4,3x,3(1x,3(1x,i4)))') i2,ig_tiny(i2,i1),kg(1:3,ig_tiny(i2,i1))
     end do
   end do

   write(std_out,'(a)')       '%get_g_tiny: on sorted squares:'
   write(std_out,'(a,es14.7)') ' squares different, if relative difference > diffgsq=',diffgsq
   write(std_out,'(a,i4,a)')   ' index for first npw_tiny=', npw_tiny, ' overall:'
   do i1=1,npw_tiny
     write(std_out,'(1x,a,2(1x,i4),1x,es14.7)')&
&     'igsq_tiny=',i1,igsq_tiny(i1),gsq(igsq_tiny(i1))
   end do
 end if

!Map to unsorted
 igsq_tiny(:)=index_g(igsq_tiny(:))

 if(tfirst) then
   write(std_out,'(a)')       '%get_g_tiny: on unsorted squares:'
   do i1=1,npw_tiny
     write(std_out,'(1x,a,i3,2x,a,3(1x,i4),2x,a,1x,es14.7)')&
&     'ipw=',igsq_tiny(i1),'kg=',kg(:,igsq_tiny(i1)),'gsq=',gsq_unsorted(igsq_tiny(i1))
   end do
 end if

!Perform deallocations
 ABI_DEALLOCATE(gsq)

 if(tfirst) tfirst=.false.
!DEBUG
!write(std_out,*) '%get_g_tiny: done'
!ENDDEBUG

end subroutine get_g_tiny
!!***
