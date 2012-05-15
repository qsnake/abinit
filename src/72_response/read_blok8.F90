!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_blok8
!!
!! NAME
!! read_blok8
!!
!! FUNCTION
!! This routine reads blocks of data in the DDBs.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! mpert =maximum number of ipert
!! msize=maximum size of the arrays flags and values
!! nunit=unit number for the data block file
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! ddb_blk = ddb block datastructure
!! ddb_blk%typ=type of the block:
!!   0 => total energy
!!   1 => second-order energy derivatives, non-stationary block
!!   2 => second-order energy derivatives, stationary block
!!   3 => third-order energy derivatives
!!   4 => first-order energy derivatives: forces, stresses and polarization
!!   5 => second-order eigenvalue derivatives
!! ddb_blk%flg(msize)=flag for every matrix element (0=> the element is
!!  not in the data block), (1=> the element is in the data blok)
!! ddb_blk%qpt(9)=wavevector of the perturbation(s). The elements from
!!  1 to 3 are used if we are dealing with the 2nd derivative of
!!  total energy (only one wavevector), while all elements are
!!  used in case of a third order derivative of total energy
!!  (three wavevector could be present)
!! ddb_blk%nrm(3)=normalization factors for the three allowed wavevectors.
!! ddb_blk%val(2,msize)=real(dp), complex, value of the
!!  matrix elements that are present in the data block
!! blkval2(2,msize,mband,nkpt) = value of the matrix elements
!!  that are present in a block of EIGR2D/EIGI2D
!!
!! NOTES
!! only executed by one processor.
!!
!! PARENTS
!!      mblktyp1,mblktyp5,rdddb9,thmeig
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine read_blok8(ddb_blk,iblok,&
&     mband,mpert,msize,nkpt,nunit,&
&     blkval2,kpt) !optional

 use m_profiling

 use defs_basis
 use m_ddb_blk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_blok8'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mband,mpert,msize,nkpt,nunit
 integer, intent(in) :: iblok
 type(ddb_blk_type), pointer :: ddb_blk
!arrays
 real(dp),intent(out),optional :: kpt(3,nkpt)
 real(dp),intent(out),optional :: blkval2(2,msize,mband,nkpt)

!Local variables -------------------------
!scalars
 integer :: band,iband,idir1,idir2,idir3,ii,ikpt,index,ipert1,ipert2,ipert3
 integer :: nelmts
 real(dp) :: ai,ar
 character(len=32) :: name
 character(len=500) :: message

! *********************************************************************
 
!Zero every flag
 ddb_blk%flg(1:msize, iblok)=0
 if(present(blkval2))blkval2(:,:,:,:)=zero
 if(present(kpt))kpt(:,:)=zero

!Read the block type and number of elements
 read(nunit,*)
 read(nunit, '(a32,12x,i8)' )name,nelmts
 if(name==' 2nd derivatives (non-stat.)  - ')then
   ddb_blk%typ(iblok)=1
 else if(name==' 2nd derivatives (stationary) - ')then
   ddb_blk%typ(iblok)=2
 else if(name==' 3rd derivatives              - ')then
   ddb_blk%typ(iblok)=3
 else if(name==' Total energy                 - ')then
   ddb_blk%typ(iblok)=0
 else if(name==' 1st derivatives              - ')then
   ddb_blk%typ(iblok)=4
 else if(name==' 2nd eigenvalue derivatives   - ')then
   ddb_blk%typ(iblok)=5
 else
   write(message, '(a,a,a,a,a,a,a,a)' )&
&   ' blok8 : ERROR -',ch10,&
&   '  The following string appears in the DDB in place of',&
&   ' the block type description :',ch10,name,ch10,&
&   '  Action : check your DDB.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Read the 2nd derivative block
 if(ddb_blk%typ(iblok)==1.or.ddb_blk%typ(iblok)==2)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
     write(message, '(a,a,a,a,a,i10,a,i10,a,a,a)' )&
&     ' blok8 : ERROR -',ch10,&
&     '  There is not enough space to read a second-derivative block.',&
&     ch10,'  The size provided is only ',msize,' although ',&
&     3*mpert*3*mpert,' is needed.',ch10,&
&     '  Action : increase msize and recompile.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Read the phonon wavevector
   read(nunit, '(4x,3es16.8,f6.1)' )(ddb_blk%qpt(ii,iblok),ii=1,3),ddb_blk%nrm(1,iblok)

!  write(std_out,*)' Blok8 : For graphic purposes, format changed'

!  Read every element
   do ii=1,nelmts
     read(nunit,*)idir1,ipert1,idir2,ipert2,ar,ai
     index=idir1+3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
     ddb_blk%flg(index,iblok)=1
     ddb_blk%val(1,index,iblok)=ar
     ddb_blk%val(2,index,iblok)=ai
   end do

!  Read the 3rd derivative block
 else if(ddb_blk%typ(iblok)==3)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert*3*mpert))then
     write(message, '(a,a,a,a,a,i10,a,i10,a,a,a)' )&
&     ' blok8 : ERROR -',ch10,&
&     '  There is not enough space to read a third-derivative block.',&
&     ch10,'  The size provided is only ',msize,' although ',&
&     3*mpert*3*mpert*3*mpert,' is needed.',ch10,&
&     '  Action : increase msize and recompile.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Read the perturbation wavevectors
   read(nunit,'(4x,3es16.8,f6.1)')(ddb_blk%qpt(ii,iblok),ii=1,3),ddb_blk%nrm(1,iblok)
   read(nunit,'(4x,3es16.8,f6.1)')(ddb_blk%qpt(ii,iblok),ii=4,6),ddb_blk%nrm(2,iblok)
   read(nunit,'(4x,3es16.8,f6.1)')(ddb_blk%qpt(ii,iblok),ii=7,9),ddb_blk%nrm(3,iblok)

!  Read every element
   do ii=1,nelmts
     read(nunit,'(6i4,2d22.14)')&
&     idir1,ipert1,idir2,ipert2,idir3,ipert3,ar,ai
     index=idir1+                                              &
&     3*((ipert1-1)+mpert*((idir2-1)+                 &
&     3*((ipert2-1)+mpert*((idir3-1)+3*(ipert3-1)))))
     ddb_blk%flg(index,iblok)=1
     ddb_blk%val(1,index,iblok)=ar
     ddb_blk%val(2,index,iblok)=ai
   end do

!  Read the total energy
 else if(ddb_blk%typ(iblok)==0)then

!  First check if there is enough space to read it
   if(msize<1)then
     write(message, '(a,a,a,a,a,i10,a,a,a,a)' )&
&     ' blok8 : ERROR -',ch10,&
&     '  There is not enough space to read a total energy block.',&
&     ch10,'  The size provided is only ',msize,' although 1',&
&     ' is needed.',ch10,&
&     '  Action : increase msize and recompile.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Read the total energy
   read(nunit,'(2d22.14)')ar,ai
   ddb_blk%flg(1,iblok)=1
   ddb_blk%val(1,1,iblok)=ar
   ddb_blk%val(2,1,iblok)=ai


!  Read the 1st derivative block
 else if (ddb_blk%typ(iblok) == 4) then

!  First check if there is enough space to read it
   if (msize < (3*mpert)) then
     write(message, '(a,a,a,a,a,i10,a,i10,a,a,a)' )&
&     ' blok8 : ERROR -',ch10,&
&     '  There is not enough space to read a first-derivative block.',&
&     ch10,'  The size provided is only ',msize,' although ',&
&     3*mpert,' is needed.',ch10,&
&     '  Action : increase msize and recompile.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Read every element
   do ii=1,nelmts
     read(nunit,'(2i4,2d22.14)')&
&     idir1,ipert1,ar,ai
     index=idir1 + 3*(ipert1 - 1)
     ddb_blk%flg(index,iblok)=1
     ddb_blk%val(1,index,iblok)=ar
     ddb_blk%val(2,index,iblok)=ai
   end do

!  Read the 2nd eigenvalue derivative block
 else if(ddb_blk%typ(iblok)==5)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
     write(message, '(a,a,a,a,a,i10,a,i10,a,a,a)' )&
&     ' blok8 : ERROR -',ch10,&
&     '  There is not enough space to read a second-derivative block.',&
&     ch10,'  The size provided is only ',msize,' although ',&
&     3*mpert*3*mpert*mband*nkpt,' is needed.',ch10,&
&     '  Action : increase msize and recompile.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Read the phonon wavevector
!  read(nunit, '(5x,54a)' ) message
!  write(std_out,*)'blok9: ',message(49:54)
   read(nunit, '(4x,3es16.8,f6.1)' )(ddb_blk%qpt(ii,iblok),ii=1,3),ddb_blk%nrm(1,iblok)

!  Read the K point and band
   if(present(blkval2).and.present(kpt))then
     do ikpt=1,nkpt
       read(nunit, '(9x,3es16.8)')(kpt(ii,ikpt),ii=1,3)
       do iband=1,mband
         read(nunit, '(6x,i3)') band
!        Read every element
         do ii=1,nelmts
           read(nunit,*)idir1,ipert1,idir2,ipert2,ar,ai
           index=idir1+3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
           ddb_blk%flg(index,iblok)=1
           blkval2(1,index,iband,ikpt)=ar
           blkval2(2,index,iband,ikpt)=ai
         end do !nelmts
       end do  !band
     end do   !kpt
   end if
 end if

end subroutine read_blok8
!!***
