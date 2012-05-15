!{\src2tex{textfont=tt}}
!!****f* ABINIT/wght9
!!
!! NAME
!! wght9
!!
!! FUNCTION
!! Generates a weight to each R points of the Big Box and for
!! each pair of atoms
!! For each R points included in the space generates by moving
!! the unit cell around each atom; the weight will be one.
!! Border conditions are provided.
!! The R points outside the chosen space will have a 0 weight.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,JCC)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! brav = Bravais lattice (1=S.C.;2=F.C.C.;4=Hex.)
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! natom= Number of atoms in the unit cell
!! ngqpt(6)= Numbers used to sample the Brillouin zone
!! nqpt= Number of q points used in the homogeneous grid
!!  sampling the Brillouin zone
!! nqshft=number of shift vectors in the repeated cell
!! nrpt=Number of R points in the Big Box
!! qshft(3,nqshft)=vectors that will be used to determine
!!  the shifts from (0. 0. 0.)
!! rcan(3,natom)=Atomic position in canonical coordinates
!! rpt(3,nprt)=Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3))
!!
!! OUTPUT
!! wghatm(natom,natom,nrpt)
!!  = Weight associated to the couple of atoms
!!    and the R vector
!!  The vector r(atom2)-r(atom1)+rpt should be inside the moving box
!! ngqpt(6)= can be modified
!!
!! PARENTS
!!      mkifc9
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wght9(brav,gprim,natom,ngqpt,nqpt,nqshft,nrpt,qshft,rcan,rpt,wghatm)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wght9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,natom,nqpt,nqshft,nrpt
!arrays
 integer,intent(inout) :: ngqpt(9)
 real(dp),intent(in) :: gprim(3,3),qshft(3,4),rcan(3,natom),rpt(3,nrpt)
 real(dp),intent(out) :: wghatm(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,ii,iqshft,irpt,jqshft,nbordh,tok
 real(dp) :: factor,sum
 character(len=500) :: message
!arrays
 integer :: nbord(9)
 real(dp) :: rdiff(9),red(3,3)

! *********************************************************************

!First analyze the vectors qshft
 if(nqshft/=1)then
   if(brav==4)then
     write(message, '(a,a,a,a,a,i5,a,a,a)' )&
&     ' wght9 : ERROR -',ch10,&
&     '  For the time being, only nqshft=1',ch10,&
&     '  is allowed with brav=4, while it is nqshft=',nqshft,'.',ch10,&
&     '  Action : in the input file, correct either brav or nqshft.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if(nqshft==2)then
!    Make sure that the q vectors form a BCC lattice
     do ii=1,3
       if(abs(abs(qshft(ii,1)-qshft(ii,2))-.5_dp)>1.d-10)then
         write(message, '(a,a,a,a,a,a,a,a,a)' )&
&         ' wght9 : ERROR -',ch10,&
&         '  The test of the q1shft vectors shows that they',ch10,&
&         '  do not generate a body-centered lattice, which',ch10,&
&         '  is mandatory for nqshft=2.',ch10,&
&         '  Action : change the q1shft vectors in your input file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   else if(nqshft==4)then
!    Make sure that the q vectors form a FCC lattice
     do iqshft=1,3
       do jqshft=iqshft+1,4
         tok=0
         do ii=1,3
!          Test on the presence of a +-0.5 difference
           if(abs(abs(qshft(ii,iqshft)-qshft(ii,jqshft))-.5_dp)&
&           <1.d-10)                       tok=tok+1
!          Test on the presence of a 0 or +-1.0 difference
           if(abs(abs(qshft(ii,iqshft)-qshft(ii,jqshft))-1._dp)&
&           <1.d-10     .or.&
&           abs(qshft(ii,iqshft)-qshft(ii,jqshft)) < 1.d-10)&
&           tok=tok+4
         end do
!        Test 1 should be satisfied twice, and test 2 once
         if(tok/=6)then
           write(message, '(9a)' )&
&           ' wght9 : ERROR -',ch10,&
&           '  The test of the q1shft vectors shows that they',ch10,&
&           '  do not generate a face-centered lattice, which',ch10,&
&           '  is mandatory for nqshft=4.',ch10,&
&           '  Action : change the q1shft vectors in your input file.'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
       end do
     end do

   else

     write(message, '(3a,i4,3a)' )&
&     ' wght9 : ERROR -',ch10,&
&     '  nqshft must be 1, 2 or 4. It is nqshft=',nqshft,'.',ch10,&
&     '  Action : change nqshft in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')

   end if
 end if

 factor=0.5_dp
 if(brav==2 .or. brav==3)then
   factor=0.25_dp
 end if
 if(nqshft/=1)factor=factor*2

!DEBUG
!write(std_out,*)'factor,ngqpt',factor,ngqpt(1:3)
!ENDDEBUG


!Begin the big loop on ia and ib
 do ia=1,natom
   do ib=1,natom

!    Simple Lattice
     if (brav==1) then
!      In this case, it is better to work in reduced coordinates
!      As rcan is in canonical coordinates, => multiplication by gprim
       do ii=1,3
         red(1,ii)=  rcan(1,ia)*gprim(1,ii)&
&         +rcan(2,ia)*gprim(2,ii)&
&         +rcan(3,ia)*gprim(3,ii)
         red(2,ii)=  rcan(1,ib)*gprim(1,ii)&
&         +rcan(2,ib)*gprim(2,ii)&
&         +rcan(3,ib)*gprim(3,ii)
       end do
     end if

     do irpt=1,nrpt

!      Initialization of the weights to 1.0
       wghatm(ia,ib,irpt)=1.0_dp

!      Compute the difference vector

!      Simple Cubic Lattice
       if (brav==1) then
!        Change of rpt to reduced coordinates
         do ii=1,3
           red(3,ii)=  rpt(1,irpt)*gprim(1,ii)&
&           +rpt(2,irpt)*gprim(2,ii)&
&           +rpt(3,irpt)*gprim(3,ii)
           rdiff(ii)=red(2,ii)-red(1,ii)+red(3,ii)
         end do
!        Other lattices
       else
         do ii=1,3
           rdiff(ii)=rcan(ii,ib)-rcan(ii,ia)+rpt(ii,irpt)
         end do
!        DEBUG
!        if(ia==1 .and. ib==10)&
!        &     write(std_out,'(a,i5,a,3es16.6)' )' irpt=',irpt,'  rdiff=',rdiff(1:3)
!        ENDDEBUG
       end if

!      ***************************************************************

!      Assignement of weights

       if(nqshft==1 .and. brav/=4)then

         do ii=1,3
!          If the rpt vector is greater than
!          the allowed space => weight = 0.0
           if (abs(rdiff(ii))-1.0d-10>factor*ngqpt(ii)) then
             wghatm(ia,ib,irpt)=0.0_dp
!            If the point is in a boundary position => weight/2
           else if (abs(abs(rdiff(ii))-factor*ngqpt(ii))&
&             <=1.0d-10) then
!            If the point is in a boundary position => weight/2
             wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
           end if
         end do

!        DEBUG
!        if(ia==1 .and. ib==10)&
!        &     write(std_out,*)' wghatm(ia,ib,irpt)=',wghatm(ia,ib,irpt)
!        ENDDEBUG

       else if(brav==4)then
!        Hexagonal

!        Examination of the X and Y boundaries in order to form an hexagon
!        First generate the relevant boundaries
         rdiff(4)=0.5_dp*( rdiff(1)+sqrt(3.0_dp)*rdiff(2) )
         ngqpt(4)=ngqpt(1)
         rdiff(5)=0.5_dp*( rdiff(1)-sqrt(3.0_dp)*rdiff(2) )
         ngqpt(5)=ngqpt(1)

!        Test the four inequalities
         do ii=1,5
           if(ii/=2)then

             nbord(ii)=0
!            If the rpt vector is greater than
!            the allowed space => weight = 0.0
             if (abs(rdiff(ii))-1.0d-10>factor*ngqpt(ii)) then
               wghatm(ia,ib,irpt)=0.0_dp
             else if (abs(abs(rdiff(ii))-factor*ngqpt(ii))&
&               <=1.0d-10) then
!              If the point is in a boundary position increment nbord(ii)
               nbord(ii)=1
             end if

           end if
         end do

!        Computation of weights
         nbordh=nbord(1)+nbord(4)+nbord(5)
         if (nbordh==1) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
         else if (nbordh==2) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/3
         else if (nbordh/=0) then
           write(message, '(a,a,a)' )&
&           ' wght9 : BUG -',ch10,&
&           '  There is a problem of borders and weights (hex).'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
         if (nbord(3)==1)then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
         end if

!        BCC packing of k-points
       else if(nqshft==2 .and. brav/=4)then

!        First, generate the relevant boundaries
         rdiff(4)= rdiff(1)+rdiff(2)
         rdiff(5)= rdiff(1)-rdiff(2)
         rdiff(6)= rdiff(1)+rdiff(3)
         rdiff(7)= rdiff(1)-rdiff(3)
         rdiff(8)= rdiff(3)+rdiff(2)
         rdiff(9)= rdiff(3)-rdiff(2)
         if(ngqpt(2)/=ngqpt(1) .or. ngqpt(3)/=ngqpt(1))then
           write(message, '(a,a,a,a,a,3i6,a,a,a,a)' )&
&           ' wght9 : ERROR -',ch10,&
&           '  In the BCC case, the three ngqpt numbers ',ch10,&
&           '    ',ngqpt(1),ngqpt(2),ngqpt(3),ch10,&
&           '  should be equal.',ch10,&
&           '  Action : use identical ngqpt(1:3) in your input file.'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
         do ii=4,9
           ngqpt(ii)=ngqpt(1)
         end do

!        Test the relevant inequalities
         nbord(1)=0
         do ii=4,9
!          If the rpt vector is greater than
!          the allowed space => weight = 0.0
           if (abs(rdiff(ii))-1.0d-10>factor*ngqpt(ii)) then
             wghatm(ia,ib,irpt)=0.0_dp
           else if (abs(abs(rdiff(ii))-factor*ngqpt(ii))&
&             <=1.0d-10) then
!            If the point is in a boundary position increment nbord(1)
             nbord(1)=nbord(1)+1
           end if
         end do

!        Computation of weights
         if (nbord(1)==1) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
         else if (nbord(1)==2) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/3
         else if (nbord(1)==3) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/4
         else if (nbord(1)==4) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/6
         else if (nbord(1)/=0) then
           write(message, '(a,a,a)' )&
&           ' wght9 : BUG -',ch10,&
&           '  There is a problem of borders and weights (BCC).'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if

!        FCC packing of k-points
       else if(nqshft==4 .and. brav/=4)then

!        First, generate the relevant boundaries
         rdiff(4)= (rdiff(1)+rdiff(2)+rdiff(3))*2._dp/3._dp
         rdiff(5)= (rdiff(1)-rdiff(2)+rdiff(3))*2._dp/3._dp
         rdiff(6)= (rdiff(1)+rdiff(2)-rdiff(3))*2._dp/3._dp
         rdiff(7)= (rdiff(1)-rdiff(2)-rdiff(3))*2._dp/3._dp
         if(ngqpt(2)/=ngqpt(1) .or. ngqpt(3)/=ngqpt(1))then
           write(message, '(a,a,a,a,a,3i6,a,a,a,a)' )&
&           ' wght9 : ERROR -',ch10,&
&           '  In the FCC case, the three ngqpt numbers ',ch10,&
&           '    ',ngqpt(1),ngqpt(2),ngqpt(3),ch10,&
&           '  should be equal.',ch10,&
&           '  Action : use identical ngqpt(1:3) in your input file.'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
         do ii=4,7
           ngqpt(ii)=ngqpt(1)
         end do

!        Test the relevant inequalities
         nbord(1)=0
         do ii=1,7

!          If the rpt vector is greater than
!          the allowed space => weight = 0.0
           if (abs(rdiff(ii))-1.0d-10>factor*ngqpt(ii)) then
             wghatm(ia,ib,irpt)=0.0_dp
!            If the point is in a boundary position increment nbord(1)
           else if (abs(abs(rdiff(ii))-factor*ngqpt(ii))&
&             <=1.0d-10) then
!            DEBUG
!            write(std_out,*)'ia,ib,irpt,ii,nbord(1)'
!            write(std_out,*)ia,ib,irpt,ii,nbord(1)
!            ENDDEBUG
             nbord(1)=nbord(1)+1
           end if
         end do

!        Computation of weights
         if (nbord(1)==1) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/2
         else if (nbord(1)==2) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/3
         else if (nbord(1)==3) then
           wghatm(ia,ib,irpt)=wghatm(ia,ib,irpt)/4
         else if (nbord(1)/=0 .and. wghatm(ia,ib,irpt)>1.d-10) then
!          Interestingly nbord(1)==4 happens for some points
!          outside of the volume
           write(message, '(a,a,a)' )&
&           ' wght9 : BUG -',ch10,&
&           '  There is a problem of borders and weights (FCC).'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if

       else

         write(message, '(a,a,a,a,a,i6,a)' )&
&         ' wght9 : BUG -',ch10,&
&         '  One should not arrive here ... ',ch10,&
&         '  The value nqshft ',nqshft,' is not available'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')

       end if
!      Assignement of weights is done

     end do

!    End of the double loop on ia and ib
   end do
 end do

!Check the results
 do ia=1,natom
   do ib=1,natom
     sum=0.0_dp
     do irpt=1,nrpt
!      Check if the sum of the weights is equal
!      to the number of q points
       sum=sum+wghatm(ia,ib,irpt)
!      DEBUG
!      write(std_out,'(a,3i5)' )' atom1, atom2, irpt ; rpt ; wghatm ',ia,ib,irpt
!      write(std_out,'(3es16.6,es18.6)' )&
!      &    rpt(1,irpt),rpt(2,irpt),rpt(3,irpt),wghatm(ia,ib,irpt)
!      ENDDEBUG
     end do
     if (abs(sum-nqpt)>1.0d-10) then
       write(message, '(a,a,a,a,a,2i4,a,a,es14.4,a,a,i4,a,a,a,a,a,a)' )&
&       ' wght9 : BUG -',ch10,&
&       '  The sum of the weight is not equal to nqpt.',ch10,&
&       '  atoms :',ia,ib,ch10,&
&       '  The sum of the weights is : ',sum,ch10,&
&       '  The number of q points is : ',nqpt,ch10,&
&       '  You might increase "buffer" in bigbx9.f, and recompile.',ch10,&
&       '  Actually, this can also happen when ngqpt is 0 0 0,',ch10,&
&       '  if brav/=1, in which case you should change brav to 1.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end do
 end do

end subroutine wght9
!!***
