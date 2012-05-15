!{\src2tex{textfont=tt}}
!!****f* ABINIT/asrif9
!!
!! NAME
!! asrif9
!!
!! FUNCTION
!! Imposes the Acoustic Sum Rule to Interatomic Forces
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! asr= Option for the imposition of the ASR
!!  0 => no ASR,
!!  1 => modify "asymmetrically" the diagonal element
!!  2 => modify "symmetrically" the diagonal element
!! natom= Number of atoms in the unit cell
!! nrpt= Number of R points in the Big Box
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!! wghatm(natom,natom,nrpt)= Weight associated to the couple of atoms
!!  and the R vector
!! atmfrc(2,3,natom,3,natom,nrpt)= Interatomic Forces
!!
!! OUTPUT
!! atmfrc(2,3,natom,3,natom,nrpt)= ASR-imposed Interatomic Forces
!!
!! TODO
!! List of ouput should be included.
!!
!! PARENTS
!!      hybrid9,mkifc9
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine asrif9(asr,atmfrc,natom,nrpt,rpt,wghatm)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'asrif9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,natom,nrpt
!arrays
 real(dp),intent(in) :: rpt(3,nrpt),wghatm(natom,natom,nrpt)
 real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: found,ia,ib,irpt,izero,mu,nu
 real(dp) :: sum
 character(len=500) :: message

! *********************************************************************

 if(asr==1.or.asr==2)then

   found=0
!  Search for the R vector which is equal to ( 0 , 0 , 0 )
!  This vector leaves the atom a on itself !
   do irpt=1,nrpt
     if (abs(rpt(1,irpt))<=1.0d-10.and.&
&     abs(rpt(2,irpt))<=1.0d-10.and.&
&     abs(rpt(3,irpt))<=1.0d-10) then
       found=1
       izero=irpt
     end if
     if (found==1) exit
   end do
   if(found==0)then
     write(message, '(a,a,a)' )&
&     ' asrif9.f : BUG -',ch10,&
&     '  Not able to find the vector R=(0,0,0).'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   do mu=1,3
     do nu=1,3
       do ia=1,natom
         sum=0.0_dp
         do ib=1,natom

!          Get the sum of interatomic forces acting
!          on the atom ia,
!          either in a symmetrical manner, or an
!          unsymmetrical one.
           if(asr==1)then
             do irpt=1,nrpt
               sum=sum+wghatm(ia,ib,irpt)*atmfrc(1,mu,ia,nu,ib,irpt)
             end do
           else if(asr==2)then
             do irpt=1,nrpt
               sum=sum+&
&               (wghatm(ia,ib,irpt)*atmfrc(1,mu,ia,nu,ib,irpt)+&
&               wghatm(ia,ib,irpt)*atmfrc(1,nu,ia,mu,ib,irpt))/2
             end do
           end if
         end do

!        Correct the self-interaction in order
!        to fulfill the ASR
         atmfrc(1,mu,ia,nu,ia,izero)=atmfrc(1,mu,ia,nu,ia,izero)-sum
         if(asr==2)then
           atmfrc(1,nu,ia,mu,ia,izero)=atmfrc(1,mu,ia,nu,ia,izero)
         end if

       end do
     end do
   end do
 end if

end subroutine asrif9
!!***
