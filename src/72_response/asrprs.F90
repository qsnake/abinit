!{\src2tex{textfont=tt}}
!!****f* ABINIT/asrprs
!! NAME
!! asrprs
!!
!! FUNCTION
!! Imposition of the Acoustic sum rule on the InterAtomic Forces Plus Rotational
!! Symmetry
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2012 ABINIT group (NH)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  asr=(3 => 1D systems, all elements are modified to give ASR and 
!!            rotational symmetry)
!!      (4 => 0D systems, all elements are modified to give ASR and 
!!            rotational symmetry)
!!  asrflg=(1 => the correction to enforce asr is computed from
!!           d2cart, but NOT applied;
!!          2 => one uses the previously determined correction)
!!  minvers=previously calculated inverted coefficient matrix
!!  mpert =maximum number of ipert
!!  natom=number of atom
!!  rotinv=(1,2,3 => for linear systems along x,y,z
!!          4 => non-linear molecule 
!!  xcart=cartesian coordinates of the ions
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!! d2cart=matrix of second derivatives of total energy, in cartesian coordinates
!! minvers=inverse of the supermatrix for future application of the corrections
!!
!! NOTES
!!
!! PARENTS
!!      anaddb,mkphbs
!!
!! CHILDREN
!!      dgesvd,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine asrprs(asr,asrflag,rotinv,uinvers,vtinvers,singular,d2cart,mpert,natom,xcart)

 use m_profiling

 use defs_basis
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'asrprs'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: asr,asrflag,mpert,natom,rotinv
!arrays
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(inout) :: singular(1:3*natom*(3*natom-1)/2)
 real(dp),intent(inout) :: uinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)
 real(dp),intent(inout) :: vtinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)

!Local variables-------------------------------
!scalars
 integer :: column,idir1,idir2,ii,info,ipert1,ipert2,jj,n3,row,superdim
 real(dp) :: rcond,test
! real(dp) :: tau ! tau is present but commented out in this routine
 character(len=500) :: message
!arrays
 integer :: check(3,natom,3)
 real(dp) :: tmp(natom,3,3),weightf(1:natom,1:natom)
 real(dp),allocatable :: d2cartold(:,:,:,:,:),d2vecc(:),d2veccnew(:),d2vecr(:)
 real(dp),allocatable :: d2vecrnew(:),superm(:,:),umatrix(:,:),vtmatrix(:)
 real(dp),allocatable :: work(:)

! *********************************************************************
 
!DEBUG
!write (std_out,*) ' asrprs : enter'
!ENDDEBUG

 if(asr/=3 .and. asr/=4)then
   write(message,'(a,a,a,a,a,a,i6)') ch10,&
&   ' asrprs: BUG -',ch10,&
&   '  The argument asr should be 3 or 4,',ch10,&
&   '  however, asr=',asr
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 
 if (asr==3.or.asr==4)then
   write(message, '(a,a)' ) ch10, &
&   ' asrprs : imposition of the ASR for the interatomic forces and rotational invariance'
   call wrtout(std_out,message,'COLL')
 end if
 
 write(message, '(a,i6)' )&
& ' asrflag is', asrflag 
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

!if(sizein<1)then
!write(message,'(a,a,a,a,a,a,i6)') ch10,&
!&  ' asrprs: BUG -',ch10,&
!&  '  The argument sizein should be a positive number,',ch10,&
!&  '  however, sizein=',sizein
!call wrtout(std_out,message,'COLL')
!call leave_new('COLL')
!endif
!ENDDEBUG

!variables for the dimensions of the matrices

!n1=3*natom*(3*natom-1)/2
!n2=9*natom
 n3=3*natom

 superdim=9*natom*(natom-1)/2+n3

 ABI_ALLOCATE(d2vecr,(1:superdim))
 ABI_ALLOCATE(d2vecc,(1:superdim))
 d2vecr=0d0
 d2vecc=0d0

!should be changed set to delta function for debugging
 weightf=1d0
!tau=1d-10
 do ii=1, natom
!  do jj=1, ii-1
!  weightf(ii,jj)= & 
!  &     ((xcart(1,ii)-xcart(1,jj))**2+(xcart(2,ii)-xcart(2,jj))**2+(xcart(3,ii)-xcart(3,jj))**2)**tau
!  enddo
   weightf(ii,ii)=0d0
 end do

!allocate(d2cartold(2,3,natom,3,natom))
 ABI_ALLOCATE(d2cartold,(2,3,mpert,3,mpert))

 d2cartold=d2cart

!setup vector with uncorrected derivatives

 do ipert1=1, natom
   do ipert2=1, ipert1-1
     do idir1=1,3
       do idir2=1,3
         row=n3+9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(idir1-1)+idir2
         if(abs(d2cart(1,idir1,ipert1,idir2,ipert2))<1d-6)then
           d2cart(1,idir1,ipert1,idir2,ipert2)=0d0
         else
           d2vecr(row)=4*weightf(ipert1,ipert2)*d2cart(1,idir1,ipert1,idir2,ipert2)
         end if
         if(abs(d2cart(2,idir1,ipert1,idir2,ipert2))<1d-6) then
           d2cart(2,idir1,ipert1,idir2,ipert2)=0d0
         else
           d2vecc(row)=4*weightf(ipert1,ipert2)*d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end do

 if(asrflag==1) then !calculate the pseudo-inverse of the supermatrix
   ABI_ALLOCATE(superm,(1:superdim,1:superdim))
   
   superm=0d0

!  Setting up the supermatrix containing G, A, D

   do ipert1=1, natom
     do idir1=1, 3
!      Setting up G 
       idir2=mod(idir1,3)+1
       row=3*(ipert1-1)+idir1
       do ipert2=1, ipert1-1
         column=9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(rotinv-1)+idir1
         superm(column,row)=xcart(idir2,ipert2)-xcart(idir2,ipert1)
         column=9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(rotinv-1)+idir2
         superm(column,row)=xcart(idir1,ipert1)-xcart(idir1,ipert2)
       end do
       do ipert2=ipert1+1, natom
         column=9*(ipert2-1)*(ipert2-2)/2+9*(ipert1-1)+3*(idir1-1)+rotinv
         superm(column,row)=xcart(idir2,ipert2)-xcart(idir2,ipert1)
         column=9*(ipert2-1)*(ipert2-2)/2+9*(ipert1-1)+3*(idir2-1)+rotinv
         superm(column,row)=xcart(idir1,ipert1)-xcart(idir1,ipert2)
       end do
     end do
     do idir1=1, 3
!      Setting up D
       idir2=mod(idir1,3)+1
       ii=mod(idir1+1,3)+1
       do ipert2=1, ipert1-1
         row=n3+9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(rotinv-1)+idir1
         column=9*natom*(natom-1)/2+3*(ipert1-1)+idir1
         superm(column,row)=superm(column,row)+xcart(idir2,ipert2)-xcart(idir2,ipert1)
         column=9*natom*(natom-1)/2+3*(ipert1-1)+ii
         superm(column,row)=superm(column,row)+xcart(ii,ipert1)-xcart(ii,ipert2)
         row=n3+9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(idir1-1)+rotinv
         column=9*natom*(natom-1)/2+3*(ipert2-1)+idir1
         superm(column,row)=superm(column,row)+xcart(idir2,ipert1)-xcart(idir2,ipert2)
         column=9*natom*(natom-1)/2+3*(ipert2-1)+ii
         superm(column,row)=superm(column,row)+xcart(ii,ipert2)-xcart(ii,ipert1)
       end do 
!      Setting up A
       do idir2=1, 3
         do ipert2=1, ipert1-1
           column=9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(idir1-1)+idir2
           row=n3+column
           superm(column,row)=4*weightf(ipert1,ipert2)
         end do
       end do
     end do
   end do 

!  calculate the pseudo-inverse of the supermatrix

   ABI_ALLOCATE(work,(1:6*superdim))
   ABI_ALLOCATE(vtmatrix,(1:superdim))
   ABI_ALLOCATE(umatrix,(1:superdim,1:superdim))

!  singular value decomposition of superm

   call dgesvd('A','O',superdim,superdim,superm,superdim,singular,umatrix,superdim, &
&   vtmatrix, 1, work,6*superdim,info)

   ABI_DEALLOCATE(vtmatrix)
   ABI_DEALLOCATE(work)

   write(message, '(a,es16.8,es16.8)' )&
&   ' Largest and smallest values from svd', singular(1), singular(superdim) 
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  Invert U and V**T, orthogonal matrices
   
   do ii=1, superdim
     do jj=1, superdim
       uinvers(ii,jj)=umatrix(jj,ii)
       vtinvers(ii,jj)=superm(jj,ii)
     end do
   end do
   
   ABI_DEALLOCATE(umatrix)
   ABI_DEALLOCATE(superm)
   
   write(message, '(a,a)' )&
&   ' asrprs: done with asrflag 1', ch10 
   call wrtout(std_out,message,'COLL')

 end if !asrflag=1
 
 if(asrflag==2) then 
   
   ABI_ALLOCATE(d2vecrnew,(1:superdim))
   ABI_ALLOCATE(d2veccnew,(1:superdim))

!  Calculate V**T**-1 Sigma**-1 U**-1 *rhs

   do ii=1, superdim
     d2vecrnew(ii)=0d0
     d2veccnew(ii)=0d0
     do jj=1, superdim
       d2vecrnew(ii)=d2vecrnew(ii)+uinvers(ii,jj)*d2vecr(jj)
       d2veccnew(ii)=d2veccnew(ii)+uinvers(ii,jj)*d2vecc(jj)
     end do
   end do

   rcond=1d-10*singular(1)
   do ii=1, superdim
     if(singular(ii)>rcond) then
       d2vecrnew(ii)=d2vecrnew(ii)/singular(ii)
       d2veccnew(ii)=d2veccnew(ii)/singular(ii)
     else
       d2vecrnew(ii)=0d0
       d2veccnew(ii)=0d0
     end if
   end do 

   do ii=1, superdim
     d2vecr(ii)=0d0
     d2vecc(ii)=0d0
     do jj=1, superdim
       d2vecr(ii)=d2vecr(ii)+vtinvers(ii,jj)*d2vecrnew(jj)
       d2vecc(ii)=d2vecc(ii)+vtinvers(ii,jj)*d2veccnew(jj)
     end do
   end do

!  Store vector back into the matrix of 2nd order derivates
   
   do ipert1=1, natom
     do ipert2=1, ipert1-1
       do idir1=1,3
         do idir2=1,3
           row=9*(ipert1-1)*(ipert1-2)/2+9*(ipert2-1)+3*(idir1-1)+idir2
           d2cart(1,idir1,ipert1,idir2,ipert2)=d2vecr(row)
           d2cart(2,idir1,ipert1,idir2,ipert2)=d2vecc(row)
           d2cart(1,idir2,ipert2,idir1,ipert1)=d2vecr(row)
           d2cart(2,idir2,ipert2,idir1,ipert1)=d2vecc(row)
         end do
       end do
     end do
   end do
   
!  write(std_out,*) 'test the transpose symmetry'
!  do ipert1=1, natom
!  do ipert2=ipert1+1, natom
!  do idir1=1, 3
!  do idir2=1, 3
!  test=d2cart(1,idir1,ipert1,idir2,ipert2)-d2cart(1,idir2,ipert2,idir1,ipert1)
!  if(abs(test)>1d-5) then
!  write(std_out,*) idir1,ipert1,idir2,ipert2, test
!  end if
!  end do
!  end do
!  end do
!  enddo


   ABI_DEALLOCATE(d2vecrnew)
   ABI_DEALLOCATE(d2veccnew)
   
   check=0

   do ipert1=1, natom
     do idir1=1, 3
       do idir2=1, 3
         d2cart(1,idir1,ipert1,idir2,ipert1)=0d0
         d2cart(2,idir1,ipert1,idir2,ipert1)=0d0
         tmp(ipert1,idir1,idir2)=0d0
         do ipert2=1, natom
           if(ipert2/=ipert1) then
             tmp(ipert1,idir1,idir2)=tmp(ipert1,idir1,idir2) & 
&             -d2cart(1,idir1,ipert1,idir2,ipert2) & 
&             -d2cart(1,idir2,ipert2,idir1,ipert1)
           end if
         end do
       end do
     end do
   end do 
   
   do ipert1=1, natom
     do idir1=1, 3
       do idir2=1, 3
         d2cart(1,idir1,ipert1,idir2,ipert1)=tmp(ipert1,idir1,idir2)/2
         d2cart(1,idir2,ipert1,idir1,ipert1)=d2cart(1,idir1,ipert1,idir2,ipert1)
       end do
     end do
   end do
   
   write(std_out,*) 'this should all be zero'
   
   do ipert1=1, natom
     do idir1=1, 3
       do idir2=1, 3
         test=0d0
         do ipert2=1, natom
           test=test+d2cart(1,idir1,ipert1,idir2,ipert2)+d2cart(1,idir2,ipert2,idir1,ipert1)
         end do
         write(std_out,'(i3,i3,i3,es11.3)') idir1,ipert1,idir2,test
         
         write(message, '(i3,i3,i3,es11.3)' ) idir1,ipert1,idir2,test 
         call wrtout(ab_out,message,'COLL')

       end do
     end do
   end do

   write(std_out,*) 'these as well'
   do ipert2=1, natom
     do idir1=1, 3
       do idir2=1, 3
         test=0d0
         do ipert1=1, natom
           test=test+d2cart(1,idir1,ipert1,idir2,ipert2)
         end do
         write(std_out,'(i3,i3,i3,i3,es11.3)') idir1,ipert1,idir2,ipert2,test
       end do
     end do
   end do

   write(message, '(a,a)' )&
&   ' asrprs: done with asrflag 2', ch10 
   call wrtout(std_out,message,'COLL')
   
 end if !ends asrflag=2

 ABI_DEALLOCATE(d2vecr)
 ABI_DEALLOCATE(d2vecc)
 ABI_DEALLOCATE(d2cartold)

!DEBUG
!write (std_out,*) ' asrprs : exit'
!stop
!ENDDEBUG

end subroutine asrprs
!!***
