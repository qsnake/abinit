!{\src2tex{textfont=tt}}
!!****f* ABINIT/chneu9
!! NAME
!! chneu9
!!
!!
!! FUNCTION
!! Imposition of the Acoustic sum rule on the Effective charges
!! and suppress the imaginary part of the dynamical matrix
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  chneut=(0 => no ASR, 1 => equal repartition,2 => weighted repartition )
!!  mpert =maximum number of ipert
!!  natom=number of atom
!!  ntypat=number of types of atoms in unit cell
!!  selectz=selection of some parts of the effective charge tensor
!!    attached to one atom.
!!    (0=> no selection, 1=> trace only, 2=> symmetric part only)
!!  typat(natom)=type of the atom
!!  zion(ntypat)=atomic charge for every type of atom
!!
!! OUTPUT
!
!!
!! SIDE EFFECTS
!!  Input/Output
!!  d2cart=matrix of second derivatives of total energy, in cartesian
!!       coordinates
!!
!! NOTES
!
!!
!! PARENTS
!!      anaddb,gath3
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chneu9(chneut,d2cart,mpert,natom,ntypat,selectz,typat,zion)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chneu9'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: chneut,mpert,natom,ntypat,selectz
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: zion(ntypat)
 real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
 character(len=500) :: message
!arrays
 real(dp) :: sum(2)
 real(dp),allocatable :: wghtat(:)

! *********************************************************************

 ABI_ALLOCATE(wghtat,(natom))

!In case of acoustic sum rule imposition, compute the weights on
!each atom.
 if (chneut==1)then

!  The weight is the same for all atom
   do ipert1=1,natom
     wghtat(ipert1)=1./natom
   end do

 else if (chneut==2) then

!  The weight is proportional to the diagonal electronic screening
!  charge of the atom
   sum(1)=0.0_dp
   do ipert1=1,natom
     wghtat(ipert1)=0.0_dp
     do idir1=1,3
       wghtat(ipert1)=wghtat(ipert1)+&
&       d2cart(1,idir1,ipert1,idir1,natom+2)+&
&       d2cart(1,idir1,natom+2,idir1,ipert1)-2*zion(typat(ipert1))
     end do
     sum(1)=sum(1)+wghtat(ipert1)
   end do

!  Normalize the weights to unity
   do ipert1=1,natom
     wghtat(ipert1)=wghtat(ipert1)/sum(1)
   end do
 end if

!Calculation of the violation of the charge neutrality
!and imposition of the charge neutrality condition
 if (chneut/=0)then
   write(message, '(a,a,a,a,a,a,a)' )&
&   ' The violation of the charge neutrality conditions',ch10,&
&   ' by the effective charges is as follows :',ch10,&
&   '    atom        electric field',ch10,&
&   ' displacement     direction   '
   call wrtout(ab_out,message,'COLL')
   do idir1=1,3
     do idir2=1,3
       do ii=1,2
         sum(ii)=0.0_dp
         do ipert1=1,natom
           sum(ii)=sum(ii)+d2cart(ii,idir1,ipert1,idir2,natom+2)
         end do
         do ipert1=1,natom
           d2cart(ii,idir1,ipert1,idir2,natom+2)=&
&           d2cart(ii,idir1,ipert1,idir2,natom+2)-sum(ii)*wghtat(ipert1)
         end do
       end do
       write(message, '(i8,i16,2f16.6)' ) idir1,idir2,sum(1),sum(2)
       call wrtout(ab_out,message,'COLL')
     end do
   end do
   write(message, '(a)' )' '
   call wrtout(ab_out,message,'COLL')

!  The same for the symmetrical part
   do idir1=1,3
     do idir2=1,3
       do ii=1,2
         sum(ii)=0.0_dp
         do ipert2=1,natom
           sum(ii)=sum(ii)+d2cart(ii,idir1,natom+2,idir2,ipert2)
         end do
         do ipert2=1,natom
           d2cart(ii,idir1,natom+2,idir2,ipert2)=&
&           d2cart(ii,idir1,natom+2,idir2,ipert2)-sum(ii)*wghtat(ipert2)
         end do
       end do
     end do
   end do
 end if

!Selection of the trace of the effective charge tensor
!attached to each atom
 if(selectz==1)then
   do ipert1=1,natom
     do ii=1,2
       sum(ii)=0.0_dp
       do idir1=1,3
         sum(ii)=sum(ii)+d2cart(ii,idir1,ipert1,idir1,natom+2)
       end do
       do idir1=1,3
         do idir2=1,3
           d2cart(ii,idir1,ipert1,idir2,natom+2)=0.0_dp
         end do
       end do
       do idir1=1,3
         d2cart(ii,idir1,ipert1,idir1,natom+2)=sum(ii)/3.0_dp
       end do
     end do
   end do
!  Do the same for the symmetrical part of d2cart
   do ipert2=1,natom
     do ii=1,2
       sum(ii)=0.0_dp
       do idir1=1,3
         sum(ii)=sum(ii)+d2cart(ii,idir1,natom+2,idir1,ipert2)
       end do
       do idir1=1,3
         do idir2=1,3
           d2cart(ii,idir1,natom+2,idir2,ipert2)=0.0_dp
         end do
       end do
       do idir1=1,3
         d2cart(ii,idir1,natom+2,idir1,ipert2)=sum(ii)/3.0_dp
       end do
     end do
   end do
 end if

!Selection of the symmetric part of the effective charge tensor
!attached to each atom
 if(selectz==2)then
   do ipert1=1,natom
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           sum(ii)=(d2cart(ii,idir1,ipert1,idir2,natom+2)&
&           +d2cart(ii,idir2,ipert1,idir1,natom+2))/2.0_dp
           d2cart(ii,idir1,ipert1,idir2,natom+2)=sum(ii)
           d2cart(ii,idir2,ipert1,idir1,natom+2)=sum(ii)
         end do
       end do
     end do
   end do
!  Do the same for the symmetrical part of d2cart
   do ipert1=1,natom
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           sum(ii)=(d2cart(ii,idir1,ipert1,idir2,natom+2)&
&           +d2cart(ii,idir2,ipert1,idir1,natom+2))/2.0_dp
           d2cart(ii,idir1,ipert1,idir2,natom+2)=sum(ii)
           d2cart(ii,idir2,ipert1,idir1,natom+2)=sum(ii)
         end do
       end do
     end do
   end do
 end if

!Write the effective charge tensor
 write(message, '(a,a,a,a,a,a,a)' )&
& ' Effective charge tensors after ',ch10,&
& ' imposition of the charge neutrality,',ch10,&
& ' and eventual restriction to some part :',ch10,&
& '   atom    displacement  '
 call wrtout(ab_out,message,'COLL')
 do ipert1=1,natom
   do idir1=1,3
     write(message, '(2i10,3es16.6)' )ipert1,idir1,&
&     (d2cart(1,idir1,ipert1,idir2,natom+2),idir2=1,3)
     call wrtout(ab_out,message,'COLL')
   end do
 end do

!Zero the imaginary part of the dynamical matrix
 write(message, '(a)' )&
& ' Now, the imaginary part of the dynamical matrix is zeroed '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 do ipert1=1,natom
   do ipert2=1,natom
     do idir1=1,3
       do idir2=1,3
         d2cart(2,idir1,ipert1,idir2,ipert2)=0.0_dp
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(wghtat)

end subroutine chneu9
!!***
