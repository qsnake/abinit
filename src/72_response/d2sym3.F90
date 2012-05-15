!{\src2tex{textfont=tt}}
!!****f* ABINIT/d2sym3
!! NAME
!! d2sym3
!!
!! FUNCTION
!! Given a set of calculated elements of the 2DTE matrix d2,
!! build (nearly) all the other matrix elements that can be build using
!! symmetries.
!! 1. Perform first some completion by
!!    symmetrisation (exchange) over the two defining perturbations
!! 2. For each element, uses every symmetry, and build
!!    the element, in case
!!    EITHER all the needed elements are available,
!!    OR the only missing is itself
!!    OR the perturbation is the electric field, in a diamond
!!    symmetry (the last case was coded rather dirty)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  indsym(4,nsym,natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  nsym=number of space group symmetries
!!  qpt(3)=wavevector of the perturbation
!!  symq(4,2,nsym)= (integer) three first numbers define the G vector ;
!!   fourth number is zero if the q-vector is not preserved,
!!              is 1 otherwise
!!   second index is one without time-reversal symmetry,
!!                two with time-reversal symmetry
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!  timrev=1 if the time-reversal symmetry preserves the wavevector,
!!     modulo a reciprocal lattice vector, timrev=0 otherwise
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  Input/Output
!!  d2(2,3,mpert,3,mpert)= matrix of the 2DTE
!!  blkflg(3,mpert,3,mpert)= ( 1 if the element of the dynamical
!!     matrix has been calculated ; 0 otherwise)
!!
!!
!! NOTES
!! The complete search would be to have the possibility
!!   of a set of missing elements. See notes of July 2, 1994,
!!   in the blue notebook 'computer codes'
!!   The partial solution adopted here takes into
!!   account some mirror symmetries
!!   as well as the tetrahedral symmetry of the diamond lattice
!! On 010331, replaced the loops up to mpert by loops up to
!!   natom+2, because of a crash bug under Windows. However,
!!   the problem lies likely in the use of the indsym array.
!!
!! PARENTS
!!      completeperts,rdddb9,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine d2sym3(blkflg,d2,indsym,mpert,natom,nsym,&
& qpt,symq,symrec,symrel,timrev)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2sym3'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,nsym,timrev
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symq(4,2,nsym),symrec(3,3,nsym)
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(inout) :: d2(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: exch12,found,idir1,idir2,idisy1,idisy2,ipert1,ipert2
 integer :: ipesy1,ipesy2,isign,isym,ithree,itirev,noccur,quit,quit1
! integer :: ii,jj ! these appear below in commented out DEBUG sections
 real(dp) :: arg1,arg2,im,norm,re,sumi,sumr,xi,xr
!arrays
 integer :: sym1(3,3),sym2(3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)' d2sym3 : enter d2sym3 '
!ENDDEBUG

!Here look after exchange of 1 and 2 axis,
!for electric field in diamond symmetry
 exch12=0
 if ( qpt(1)**2+qpt(2)**2+qpt(3)**2 < tol16 ) then

   do isym=1,nsym
     exch12=1
     if(symrel(1,1,isym)/=0)exch12=0
     if(symrel(1,2,isym)/=1)exch12=0
     if(symrel(1,3,isym)/=0)exch12=0
     if(symrel(2,1,isym)/=1)exch12=0
     if(symrel(2,2,isym)/=0)exch12=0
     if(symrel(2,3,isym)/=0)exch12=0
     if(symrel(3,1,isym)/=0)exch12=0
     if(symrel(3,2,isym)/=0)exch12=0
     if(symrel(3,3,isym)/=1)exch12=0
!    DEBUG
!    if(exch12==1)then
!    write(std_out,*)' d2sym3 : found exchange 1 2 =',isym
!    end if
!    ENDDEBUG
     if(exch12==1)exit
   end do

 end if

!Exchange of perturbations

!Consider two cases : either time-reversal symmetry
!conserves the wavevector, or not
 if(timrev==0)then

!  do ipert1=1,mpert  See notes
   do ipert1=1,min(natom+2,mpert)
     do idir1=1,3
!      DEBUG
!      write(std_out,*) 'd2sym3: Im', d2(2,idir1,ipert1,idir1,ipert1)
!      ENDDEBUG

!      Since the matrix is hermitian, the diagonal elements are real
       d2(2,idir1,ipert1,idir1,ipert1)=0.0_dp
!      
!      do ipert2=1,mpert See notes
       do ipert2=1,min(natom+2,mpert)
         do idir2=1,3

!          If an element exists
           if(blkflg(idir1,ipert1,idir2,ipert2)==1)then

!            Either complete the symmetric missing element
             if(blkflg(idir2,ipert2,idir1,ipert1)==0)then
               d2(1,idir2,ipert2,idir1,ipert1)=&
&               d2(1,idir1,ipert1,idir2,ipert2)

               d2(2,idir2,ipert2,idir1,ipert1)=&
&               -d2(2,idir1,ipert1,idir2,ipert2)

               blkflg(idir2,ipert2,idir1,ipert1)=1

!              Or symmetrize (the matrix is hermitian) in case both exists
!              (Note : this opportunity has been disabled for more
!              obvious search for bugs in the code )
!              else
!              sumr=d2(1,idir2,ipert2,idir1,ipert1)+d2(1,idir1,ipert1,idir2,ipert2)
!              sumi=d2(1,idir2,ipert2,idir1,ipert1)-d2(1,idir1,ipert1,idir2,ipert2)
!              d2(1,idir2,ipert2,idir1,ipert1)=sumr/2.0_dp
!              d2(1,idir1,ipert1,idir2,ipert2)=sumr/2.0_dp
!              d2(2,idir2,ipert2,idir1,ipert1)=sumi/2.0_dp
!              d2(2,idir1,ipert1,idir2,ipert2)=-sumi/2.0_dp

             end if
           end if

         end do
       end do

     end do
   end do

!  Here, case with time-reversal symmetry
 else

!  do ipert1=1,mpert See notes
   do ipert1=1,min(natom+2,mpert)
     do idir1=1,3
!      do ipert2=1,mpert See notes
       do ipert2=1,min(natom+2,mpert)
         do idir2=1,3
           d2(2,idir1,ipert1,idir2,ipert2)=0.0_dp

!          If an element exists
           if(blkflg(idir1,ipert1,idir2,ipert2)==1)then

!            Either complete the symmetric missing element
             if(blkflg(idir2,ipert2,idir1,ipert1)==0)then
               d2(1,idir2,ipert2,idir1,ipert1)=&
&               d2(1,idir1,ipert1,idir2,ipert2)
               blkflg(idir2,ipert2,idir1,ipert1)=1

!              Or symmetrize (the matrix is hermitian) in case both exists
!              (Note : this opportunity has been disabled for more
!              obvious search for bugs in the code )
!              else
!              sumr=d2(1,idir2,ipert2,idir1,ipert1)+d2(1,idir1,ipert1,idir2,ipert2)
!              d2(1,idir2,ipert2,idir1,ipert1)=sumr/2.0_dp
!              d2(1,idir1,ipert1,idir2,ipert2)=sumr/2.0_dp
             end if

           end if
         end do
       end do
     end do
   end do

 end if
!Big Big Loop : symmetrize three times, because
!of some cases in which one element is not yet available
!at the first pass, and even at the second one !
 do ithree=1,3

!  Big loop on all elements
!  do ipert1=1,mpert See notes
   do ipert1=1,min(natom+2,mpert)
     do idir1=1,3
!      do ipert2=1,mpert See notes
       do ipert2=1,min(natom+2,mpert)
         do idir2=1,3

!          Will get element (idir1,ipert1,idir2,ipert2)
!          so this element should not yet be present ...
           if(blkflg(idir1,ipert1,idir2,ipert2)/=1)then

             d2(1,idir1,ipert1,idir2,ipert2)=0.0_dp
             d2(2,idir1,ipert1,idir2,ipert2)=0.0_dp

!            Loop on all symmetries, including time-reversal
             quit1=0
             do isym=1,nsym
               do itirev=1,2
                 isign=3-2*itirev

                 if(symq(4,itirev,isym)/=0)then
                   found=1

!                  Here select the symmetric of ipert1
                   if(ipert1<=natom)then
                     ipesy1=indsym(4,isym,ipert1)
                     sym1(:,:)=symrec(:,:,isym)
                   else if(ipert1==(natom+2) .and. qpt(1)**2+qpt(2)**2+qpt(3)**2 < tol16)then
                     ipesy1=ipert1
                     sym1(:,:)=symrel(:,:,isym)
                   else
                     found=0
                   end if

!                  Here select the symmetric of ipert2
                   if(ipert2<=natom)then
                     ipesy2=indsym(4,isym,ipert2)
                     sym2(:,:)=symrec(:,:,isym)
                   else if(ipert2==(natom+2) .and. qpt(1)**2+qpt(2)**2+qpt(3)**2 < tol16)then
                     ipesy2=ipert2
                     sym2(:,:)=symrel(:,:,isym)
                   else
                     found=0
                   end if

!                  Now that a symmetric perturbation has been obtained,
!                  including the expression of the symmetry matrix, see
!                  if the symmetric values are available
                   if( found==1 ) then

                     sumr=0.0_dp
                     sumi=0.0_dp
                     noccur=0
                     quit=0
                     do idisy1=1,3
                       do idisy2=1,3
                         if(sym1(idir1,idisy1)/=0 .and. sym2(idir2,idisy2)/=0 )then
                           if(blkflg(idisy1,ipesy1,idisy2,ipesy2)==1)then
                             sumr=sumr+sym1(idir1,idisy1)*sym2(idir2,idisy2)*&
&                             d2(1,idisy1,ipesy1,idisy2,ipesy2)
                             sumi=sumi+sym1(idir1,idisy1)*sym2(idir2,idisy2)*&
&                             d2(2,idisy1,ipesy1,idisy2,ipesy2)

!                            Here, in case the symmetric of the element
!                            is the element, or the symmetric with
!                            respect to permutation of perturbations
!                            (some more conditions on the time-reversal
!                            symmetry must be fulfilled although)
                           else if(  idisy1==idir1 .and. ipesy1==ipert1&
&                             .and. idisy2==idir2 .and. ipesy2==ipert2&
&                             .and.(isign==1 .or. timrev==1 &
&                             .or. (idir1==idir2 .and. ipert1==ipert2)))&
&                             then
                             noccur=noccur+sym1(idir1,idisy1)*sym2(idir2,idisy2)
                           else if(  idisy1==idir2 .and. ipesy1==ipert2&
&                             .and. idisy2==idir1 .and. ipesy2==ipert1&
&                             .and.(isign==-1 .or. timrev==1&
&                             .or. (idir1==idir2 .and. ipert1==ipert2)))&
&                             then
                             noccur=noccur+sym1(idir1,idisy1)*sym2(idir2,idisy2)

!                            Here, electric field case
                           else if( exch12==1 .and. &
&                             ipert1==natom+2 .and. ipert2==natom+2&
&                             .and.(( idisy1+idir1 ==3 &
&                             .and. idisy2==3 .and. idir2==3)&
&                             .or. ( idisy1+idir2 ==3&
&                             .and. idisy2==3 .and. idir1==3)&
&                             .or. ( idisy2+idir2 ==3&
&                             .and. idisy1==3 .and. idir1==3)&
&                             .or. ( idisy2+idir1 ==3&
&                             .and. idisy1==3 .and. idir2==3)))&
&                             then
                             noccur=noccur+sym1(idir1,idisy1)*sym2(idir2,idisy2)

                           else
!                            Not found
                             found=0
                             quit=1
                             exit
                           end if

                         end if
                       end do
                       if(quit==1)exit
                     end do
                   end if

!                  Now, if still found, put the correct value into array d2
                   if(found==1)then

!                    In case of phonons, need to take into account the
!                    time-reversal symmetry, and the shift back to the unit cell
!                    
!                    XG990712 : I am not sure this must be kept for electric field ...
!                    1) Consider time-reversal symmetry
                     sumi=isign*sumi

                     if(ipert1<=natom .and. ipert2<=natom)then
!                      2) Shift the atoms back to the unit cell.
                       arg1=two_pi*( qpt(1)*indsym(1,isym,ipert1)&
&                       +qpt(2)*indsym(2,isym,ipert1)&
&                       +qpt(3)*indsym(3,isym,ipert1) )
                       arg2=two_pi*( qpt(1)*indsym(1,isym,ipert2)&
&                       +qpt(2)*indsym(2,isym,ipert2)&
&                       +qpt(3)*indsym(3,isym,ipert2) )
                       re=cos(arg1)*cos(arg2)+sin(arg1)*sin(arg2)
!                      XG010117 Must use isign
                       im=isign*(cos(arg2)*sin(arg1)-cos(arg1)*sin(arg2))
                     else
                       re=1.0_dp
                       im=0.0_dp
                     end if

!                    Final check, could still fail if the
!                    element was its own symmetric
                     if( abs(1.0_dp-re*noccur) < 1.0d-6&
&                     .and.  abs(im*noccur)  < 1.0d-6 )then
                       found=0

!                      DEBUG
!                      write(std_out,*)' element is its own symmetric ...'
!                      ENDDEBUG
                     end if

                   end if

                   if(found==1)then

!                    DEBUG
!                    write(std_out,*)' all found !  isym, isign= ',isym,isign
!                    write(std_out,'(9i4)' )((sym1(ii,jj),ii=1,3),jj=1,3)
!                    write(std_out,'(9i4)' )((sym2(ii,jj),ii=1,3),jj=1,3)
!                    write(std_out,*)sumr,sumi
!                    ENDDEBUG

                     if(noccur==0)then
                       d2(1,idir1,ipert1,idir2,ipert2)=re*sumr-im*sumi
                       d2(2,idir1,ipert1,idir2,ipert2)=re*sumi+im*sumr
                     else
!                      See page July 2, 1994 in computer codes notebook
                       xr=re*sumr-im*sumi
                       xi=re*sumi+im*sumr
                       norm=1.0_dp+noccur**2-2.0_dp*re*noccur
                       xr=xr/norm
                       xi=xi/norm
                       d2(1,idir1,ipert1,idir2,ipert2)=&
&                       (1.0_dp-re*noccur)*xr-im*noccur*xi
                       d2(2,idir1,ipert1,idir2,ipert2)=&
&                       (1.0_dp-re*noccur)*xi+im*noccur*xr
                     end if

!                    The element has been constructed !
                     blkflg(idir1,ipert1,idir2,ipert2)=1

!                    Exit loop on symmetry operations
                     quit1=1
                     exit

                   end if

!                  End loop on all symmetries + time-reversal
                 end if
               end do
               if(quit1==1)exit
             end do

!            End big loop on all elements
           end if
         end do
       end do
     end do
   end do

!  End Big Big Loop
 end do

!DEBUG
!write(std_out,*)' d2sym3 : exit d2sym3 '
!write(std_out,*)' Write blkflg '
!do ipert1=9,12
!do idir1=1,3
!do ipert2=9,12
!do idir2=1,3
!write(std_out,*)'idir1,ipert1,idir2,ipert2,blkflg=',idir1,ipert1,idir2,ipert2,blkflg(idir1,ipert1,idir2,ipert2)
!enddo
!enddo
!enddo
!enddo
!write(std_out,*)' d2sym3 : timrev=',timrev
!ENDDEBUG

end subroutine d2sym3
!!***
