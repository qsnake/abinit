!{\src2tex{textfont=tt}}
!!****f* ABINIT/gtblk9
!!
!! NAME
!! gtblk9
!!
!! FUNCTION
!! This routine (get block) finds the block that contains the
!! information on the derivatives of the total energy specified
!! by the parameters rfphon,rfelfd,rfstrs,rftyp and
!! the phonon wavevectors qphon (and their normalisation).
!! In case the DDB does not contain this information, the
!! subroutine returns iblok=0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,MM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ddb_blk = ddb blok datastructure
!!   flg(msize,nblok)=flag for every matrix element (0=> the element
!!     is not in the data block), (1=> the element is in the data blok)
!!   nrm(3,nblok)=normalization factors for the three allowed wavevectors
!!   qpt(3,nblok)=wavevector of the perturbation(s). The elements
!!   typ(nblok)=type of the block. (1=> non-stationary block),
!!    (2=> stationary block), (3=> third order derivative).
!! mpert =maximum number of ipert
!! natom=number of atoms
!! qphon(3,3)=wavevectors for the three possible phonons
!!  (note : only one should be used in case of second
!!  derivative of total energy, because we know
!!  that the second is the opposite of this value)
!! qphnrm(3) =normalisation factors for the three possible phonons
!! qtol =tolerance for the identification of two wavevectors
!! rfphon(4) = 1=> response to phonons (for the four possible
!!  derivatives. Two should be used for a second
!!  derivative of total energy)
!! rfelfd(4) = 1=> d/dk, 2=> electric field only, 3=> both
!!  (see comment on rfphon)
!! rfstrs(4) = 1=> uniaxial stresses, 2=> shear stresses, 3=> both
!!  (see comment on rfphon)
!! rftyp =
!!  0 => total energy
!!  1 => non-stationary formation of the 2nd derivative
!!  2 => stationary formulation of the 2nd derivative
!!  3 => third derivative of total energy
!!  4 => first-order derivatives of total energy
!!
!! OUTPUT
!! iblok= number of the block that corresponds to the specifications
!!
!! NOTES
!!
!! PARENTS
!!      anaddb,mkphbs,refineblk,thmeig
!!
!! CHILDREN
!!      gamma9,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gtblk9(ddb_blk,iblok,mpert,&
& natom,qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

 use m_profiling

 use defs_basis
 use m_ddb_blk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gtblk9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_77_ddb, except_this_one => gtblk9
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,rftyp
 integer,intent(out) :: iblok
 real(dp),intent(in) :: qtol
!arrays
 integer,intent(in) :: rfelfd(4),rfphon(4)
 integer,intent(in) :: rfstrs(4)
 real(dp),intent(inout) :: qphnrm(3),qphon(3,3)
 type(ddb_blk_type), pointer :: ddb_blk

!Local variables -------------------------
!scalars
 integer :: blkgam,ider,idir,idir1,idir2,idir3,ii,index,ipert,ipert1,ipert2
 integer :: ipert3,nder,ok
 character(len=500) :: message

!arrays
 integer :: gamma(3)
 integer,allocatable :: worki(:,:)
 real(dp) :: qpt(3)

! *********************************************************************

 write(message,'(a)' )' gtblk9 : enter gtblk9 '
 call wrtout(std_out,message,'COLL')

!Get the number of derivative
 if(rftyp==1.or.rftyp==2)then
   nder=2
 else if(rftyp==3)then
   nder=3
 else if(rftyp==0)then
   nder=0
 else if(rftyp==4)then
   nder=1
 else
   write(message, '(a,a,a,i5,a)' )&
&   ' gtblk9 : BUG -',ch10,&
&   '  rftyp is equal to ',rftyp,'. The only allowed values are 0, 1, 2, 3 or 4.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!In case of a second-derivative, a second phonon wavevector is
!provided.
 if(nder==2)then
   do ii=1,3
     qphon(ii,2)=-qphon(ii,1)
   end do
   qphnrm(2)=qphnrm(1)
 end if

!In case of a third derivative, the sum of wavevectors to gamma
!is checked
 if (nder == 3) then
   qpt(:) = qphon(:,1)/qphnrm(1) + qphon(:,2)/qphnrm(2) + qphon(:,3)/qphnrm(3)
   call gamma9(gamma(nder),qpt,qphnrm(1),qtol)
   if (gamma(nder) == 0) then
     write(message,'(a,a,a,a,a)')&
&     ' gtblk9 : ERROR -',ch10,&
&     '  the sum of the wavevectors of the third-order energy is ',ch10,&
&     '  not equal to zero'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if



!Check the validity of the requirement
 do ider=1,nder

!  Identifies if qphon is at gamma
   call gamma9(gamma(ider),qphon(1:3,ider),qphnrm(ider),qtol)

   if(gamma(ider)==0)then
     if(rfstrs(ider)/=0.or.rfelfd(ider)/=0)then
       write(message, '(a,a,a,a)' )&
&       ' gtblk9 : BUG -',ch10,&
&       '  Not yet able to handle stresses or electric fields',ch10,&
&       '  with non-zero wavevector.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end if
 end do

!Initialise the perturbation table
 ABI_ALLOCATE(worki,(mpert,4))
 worki(:,1:nder)=0

!Build the perturbation table
 do ider=1,nder
!  First the phonons
   if(rfphon(ider)==1)then
     do ipert=1,natom
       worki(ipert,ider)=1
     end do
   end if
!  Then the d/dk
   if(rfelfd(ider)==1.or.rfelfd(ider)==3)then
     worki(natom+1,ider)=1
   end if
!  Then the electric field
   if(rfelfd(ider)==2.or.rfelfd(ider)==3)then
     worki(natom+2,ider)=1
   end if
!  Then the uniaxial stress
   if(rfstrs(ider)==1.or.rfstrs(ider)==3)then
     worki(natom+3,ider)=1
   end if
!  At last, the shear stress
   if(rfstrs(ider)==2.or.rfstrs(ider)==3)then
     worki(natom+4,ider)=1
   end if
 end do

!Examine every blok :
!"!OCL NOPREEX" to avoid zero div. (which appeared around l.185)
!for VPP Fujitsu machine, inserted by MM 19990722
!OCL NOPREEX
 do iblok=1,ddb_blk%nblok

!  If this variable is still 1 at the end of the examination,
!  the blok is the good one...
   ok=1

!  Check the type
   if(rftyp/=ddb_blk%typ(iblok))then
     ok=0
!    SP   Writing this in log file produce quickly G files.
!    write(message,'(a,a,a,i6,a,i2,a,a,a,i2,a)' )&
!    &     ' gtblk9 : COMMENT -',ch10,&
!    &     '  The blok',iblok,' with type',ddb_blk%typ(iblok),',',ch10,&
!    &     '  does not match the requirement rftyp=',rftyp,'.'
!    call wrtout(std_out,message,'COLL')
   end if

!  Check the wavevector

   if( ok==1 )then
     if (nder == 2) then
       call gamma9(blkgam,ddb_blk%qpt(1:3,iblok),ddb_blk%nrm(1,iblok),qtol)
       if(blkgam/=gamma(1))then
         ok=0
!        write(message,'(a,a,a,i6,a,i2,a,a,a,i2,a)' )&
!        &         ' gtblk9 : COMMENT -',ch10,&
!        &         '  The blok',iblok,' with gamma(1)=',gamma(1),',',ch10,&
!        &         '  does not match the requirement blkgam=',blkgam,'.'
!        call wrtout(std_out,message,'COLL')
       else if(blkgam==0)then
         do idir=1,3
           if( abs( ddb_blk%qpt(idir,iblok)/ddb_blk%nrm(1,iblok) - &
&           qphon(idir,1)/qphnrm(1) )>qtol )then
             ok=0
!            write(message,'(a,a,a,i6,a,a,i2)' )&
!            &             ' gtblk9 : COMMENT -',ch10,&
!            &             '  The blok',iblok,ch10,&
!            &             '  does not match the requirement of q along idir=',idir
!            call wrtout(std_out,message,'COLL')
!            write( message, '(es20.10,es14.4,es20.10,es14.4)' )&
!            &             ddb_blk%qpt(idir,iblok),ddb_blk%nrm(1,iblok),&
!            &             qphon(idir,1),qphnrm(1)
!            call wrtout(std_out,message,'COLL')
!            
           end if
         end do
       end if
     else if (nder == 3) then
       do ider = 1, nder
         do idir=1,3
           if( abs( ddb_blk%qpt(idir+3*(ider-1),iblok)/ddb_blk%nrm(ider,iblok) - &
&           qphon(idir,ider)/qphnrm(ider) )>qtol )then
             ok=0
!            write(message,'(a,a,a,i6,a,a,i2)' )&
!            &             ' gtblk9 : COMMENT -',ch10,&
!            &             '  The blok',iblok,ch10,&
!            &             '  does not match the requirement of q along idir=',idir
!            call wrtout(std_out,message,'COLL')
!            write( message, '(es20.10,es14.4,es20.10,es14.4)' )&
!            &             ddb_blk%qpt(idir,iblok),ddb_blk%nrm(1,iblok),&
!            &             qphon(idir,1),qphnrm(1)
!            call wrtout(std_out,message,'COLL')
           end if   ! qphon
         end do    ! idir
       end do      ! nder
     end if       ! nder
   end if        ! ok

!  Check if there is enough information in this blok
   if( ok==1 )then

     if (nder == 0) then
       if (ddb_blk%flg(1,iblok) /= 1) then
         ok = 0
         write(message,'(a,a,a,i6,a,a,a)' )&
&         ' gtblk9 : COMMENT -',ch10,&
&         '  The blok',iblok,' does not match the requirement',ch10,&
&         '  because it lacks the total energy'
         call wrtout(std_out,message,'COLL')
       end if
     end if

     do ipert1=1,mpert

       if ((nder == 4).and.(worki(ipert1,4) == 1).and.(ok == 1)) then
         do idir1 = 1, 3
           index = 3*(ipert1 - 1) + idir1
           if (ddb_blk%flg(index,iblok) /= 1) then
             ok = 0
!            write(message,'(a,a,a,i6,a,a,a,a,a,3i5)' )&
!            &             ' gtblk9 : COMMENT -',ch10,&
!            &             '  The blok',iblok,' does not match the requirement',ch10,&
!            &             '  because it lacks the element with',ch10,&
!            &             '  idir1,ipert1,index=',&
!            &             idir1,ipert1,index
!            call wrtout(std_out,message,'COLL')
           end if ! ddb_blk%flg
         end do
       end if

       if (worki(ipert1,1)==1 .and. ok==1 )then
         do ipert2=1,mpert
           if (worki(ipert2,2)==1 .and. ok==1 )then
             do idir1=1,3
               do idir2=1,3

                 if (nder == 2) then

                   index=idir1+ &
&                   3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
                   if(ddb_blk%flg(index,iblok)/=1)then
                     ok=0
!                    write(message,'(a,a,a,i6,a,a,a,a,a,5i5)' )&
!                    &                     ' gtblk9 : COMMENT -',ch10,&
!                    &                     '  The blok',iblok,' does not match the requirement',ch10,&
!                    &                     '  because it lacks the element with',ch10,&
!                    &                     '  idir1,ipert1,idir2,ipert2,index=',&
!                    &                     idir1,ipert1,idir2,ipert2,index
!                    call wrtout(std_out,message,'COLL')
                   end if ! ddb_blk%flg

                 else if (nder == 3) then

                   do ipert3 = 1, mpert
                     if (worki(ipert3,3) == 1 .and. ok == 1) then
                       do idir3 = 1, 3
                         index = idir1 + &
&                         3*((ipert1 - 1) + mpert*((idir2 - 1) + &
&                         3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))
                         if (ddb_blk%flg(index,iblok) /= 1) then
                           ok = 0
!                          write(message,'(a,a,a,i6,a,a,a,a,a,7i5)' )&
!                          &                           ' gtblk9 : COMMENT -',ch10,&
!                          &                           '  The blok',iblok,' does not match the requirement',ch10,&
!                          &                           '  because it lacks the element with',ch10,&
!                          &                           '  idir1,ipert1,idir2,ipert2,idir3,ipert3,index=',&
!                          &                           idir1,ipert1,idir2,ipert2,idir3,ipert3,index
!                          call wrtout(std_out,message,'COLL')
                         end if ! ddb_blk%flg
                       end do  ! idir3
                     end if   ! worki(ipert3,3)
                   end do    ! i3pert

                 end if

               end do
             end do
           end if
         end do
       end if
     end do
   end if

!  Now that everything has been checked, eventually end the search
   if(ok==1)exit
 end do

 if(ok==0)then
   write(message, '(a,a,a)' )&
&   ' gtblk9 : ',ch10,&
&   '  Unable to find block corresponding to the following specifications :'
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i3)' )' Type (rfmeth) =',rftyp
   call wrtout(std_out,message,'COLL')
   write(message, '(a)' ) ' ider qphon(3)         qphnrm   rfphon rfelfd rfstrs'
   call wrtout(std_out,message,'COLL')
!  DEBUG
!  write(std_out,*)' nder=',nder
!  ENDDEBUG
   do ider=1,nder
     write(message, '(i4,4f6.2,3i7)' )&
&     ider,(qphon(ii,ider),ii=1,3),qphnrm(ider),&
&     rfphon(ider),rfelfd(ider),rfstrs(ider)
     call wrtout(std_out,message,'COLL')
   end do
   iblok=0
!  write(message,'(a,a,a,a,a)' )&
!  &    ' gtblk9 : ERROR -',ch10,&
!  &    '  See the above explanation.',ch10,&
!  &    '  Action : complete your DDB, since a block is missing.'
!  call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
 end if

 if(ok==1)then
   write(message, '(a,i5,a,a)' )&
&   ' gtblk9 : found blok number ',iblok,' agree with',&
&   ' specifications '
   call wrtout(std_out,message,'COLL')
 end if

 ABI_DEALLOCATE(worki)

end subroutine gtblk9
!!***
