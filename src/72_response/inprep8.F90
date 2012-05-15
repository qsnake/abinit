!{\src2tex{textfont=tt}}
!!****f* ABINIT/inprep8
!!
!! NAME
!! inprep8
!!
!! FUNCTION
!! Open Derivative DataBase, then reads the variables that
!! must be known in order to dimension the arrays before complete reading
!!
!! Note : only one processor read or write the DDB.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! character(len=fnlen) filnam: name of input or output file
!! unddb=unit number for input or output
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.
!!
!! OUTPUT
!! dimekb=dimension of ekb (only used for norm-conserving psps)
!! lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!       =if useylm=0, max number of (l,n)   comp. over all type of psps
!! mband=maximum number of bands
!! mblktyp=largest block type
!! msym=maximum number of symmetries
!! natom=number of atoms
!! nblok=number of bloks in the DDB
!! nkpt=number of k points
!! ntypat=number of atom types
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! PARENTS
!!      anaddb,mblktyp1,mblktyp5,mrgddb
!!
!! CHILDREN
!!      chknm8,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inprep8 (dimekb,filnam,lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,&
& ntypat,unddb,usepaw,vrsddb)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inprep8'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: unddb,vrsddb
 integer, intent(out) :: msym
 integer,intent(out) :: dimekb,lmnmax,mband,mblktyp,natom,nblok,nkpt,ntypat,usepaw
 character(len=fnlen),intent(in) :: filnam

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
! integer,save :: count=0 ! this variable is used in debug sections below
 integer :: bantot,basis_size0,blktyp,ddbvrs,iband,iblok,iekb,ii,ikpt,iline,im,ios,iproj
 integer :: itypat,itypat0,jekb,lmn_size0,mproj,mpsang,nekb,nelmts,nsppol
 integer :: occopt,pspso0
 integer :: nsym
 logical :: ddbvrs_is_current_or_old,testn,testv
 character(len=12) :: string
 character(len=32) :: blkname
 character(len=500) :: message
 character(len=6) :: name_old
 character(len=80) :: rdstring
!arrays
 integer,allocatable :: nband(:)
 character(len=9) :: name(9)

! *********************************************************************

!DEBUG
!write(std_out,*)' inprep8 : enter'
!count=count+1
!write(std_out,*)' count=',count
!ENDDEBUG

!Check inprep8 version number (vrsio8) against mkddb version number
!(vrsddb)
 if (vrsio8/=vrsddb) then
   write(message, '(a,a,a,i10,a,a,i10,a)' )&
&   ' inprep8: BUG -',ch10,&
&   '  The input/output DDB version number=',vrsio8,ch10,&
&   '  is not equal to the DDB version number=',vrsddb,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Open the input derivative database.
 write(std_out,'(a,a)') ' inprep8 : open file ',trim(filnam)
 open (unit=unddb,file=filnam,status='old',form='formatted')

!Check the compatibility of the input DDB with the DDB code
 read (unddb,*)
 read (unddb,*)
 read (unddb, '(20x,i10)' )ddbvrs
 if(ddbvrs/=vrsio8 .and. ddbvrs/=vrsio8_old .and. ddbvrs/=vrsio8_old_old)then
   write(message, '(3a,i10,2a,3(a,i10),a)' )&
&   ' inprep8 : BUG - ',ch10,&
&   '  The input DDB version number=',ddbvrs,' does not agree',ch10,&
&   '  with the allowed code DDB version numbers,',&
&   vrsio8,', ',vrsio8_old,' and ',vrsio8_old_old,' .'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Read the 4 n-integers, also testing the names of data,
!and checking that their value is acceptable.
!This is important to insure that any array has a sufficient
!dimension.
 read (unddb,*)
 read (unddb,*)
 read (unddb,*)
 testn=.true.
 testv=.true.
 ddbvrs_is_current_or_old=(ddbvrs==vrsio8.or.ddbvrs==vrsio8_old)

!1. usepaw
 if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),usepaw
 else
   usepaw=0;name(1)='   usepaw'
 end if
 if(name(1)/='   usepaw')testn=.false.
!2. natom
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(2),natom
 else
   read (unddb, '(1x,a6,i10)' )name_old,natom ; name(2)='   '//name_old
 end if
 if(name(2)/='    natom')testn=.false.
 if(natom<=0)testv=.false.
!3. nkpt
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(3),nkpt
 else
   read (unddb, '(1x,a6,i10)' )name_old,nkpt ; name(3)='   '//name_old
 end if
 if(name(3)/='     nkpt')testn=.false.
 if(nkpt <=0)testv=.false.
!4. nsppol
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(4),nsppol
 else
   read (unddb, '(1x,a6,i10)' )name_old,nsppol ; name(4)='   '//name_old
 end if
 if(name(4)/='   nsppol')testn=.false.
 if(nsppol<=0.or.nsppol>2)testv=.false.
!5. nsym
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(5),nsym
 else
   read (unddb, '(1x,a6,i10)' )name_old,nsym ; name(5)='   '//name_old
 end if
 if(name(5)/='     nsym')testn=.false.
 msym = 192
 if (nsym > msym) msym=nsym
!if(nsym <=0.or.nsym >msym )testv=.false.
!6. ntypat
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(6),ntypat
 else
   read (unddb, '(1x,a6,i10)' )name_old,ntypat ; name(6)='   '//name_old
 end if
 if(name(6)/='   ntypat' .and. name(6)/='    ntype')testn=.false.
 if(ntypat<=0)testv=.false.
!7. occopt
!Before reading nband, the last parameters that define
!the dimension of some array, need to know what is their
!representation, given by occopt
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(7),occopt
 else
   read (unddb, '(1x,a6,i10)' )name_old,occopt ; name(7)='   '//name_old
 end if
 if(name(7)/='   occopt')testn=.false.
 if(occopt<0.or.occopt>8)testv=.false.
!Message if the names or values are not right
 if (.not.testn.or..not.testv) then
   write(message, '(a,a)' )' inprep8 : An error has been found in the',&
&   ' positive n-integers contained in the DDB : '
   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )   '     Expected                      Found     '
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i9,a,a,a,i10)' )&
&   '    natom , larger than',0,'    ',trim(name(2)),' =',natom
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i9,a,a,a,i10)' )&
&   '    nkpt  , larger than',0,'    ',trim(name(3)),' =',nkpt
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i1,a,a,a,i10)' )&
&   '    nsppol, either    1 or     ',2,'    ',trim(name(4)),' =',nsppol
   call wrtout(std_out,message,'COLL')
!  write(message, '(a,i10,a,a,a,i10)' )&
!  &   '    nsym  , lower than',msym,'    ',trim(name(5)),' =',nsym
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i9,a,a,a,i10)' )&
&   '    ntypat , larger than',0,'   ',trim(name(6)),' =',ntypat
   call wrtout(std_out,message,'COLL')
   write(message, '(a,a,a,i10)' )&
&   '    occopt,     equal to 0,1 or 2   ',trim(name(7)),' =',occopt
   call wrtout(std_out,message,'COLL')
   write(message, '(a,a,a)' )&
&   ' inprep8 : ERROR -',ch10,&
&   '  See the error message above.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!One more set of parameters define the dimensions of the
!array : nband. Morever, it depends on occopt and
!nkpt, and has to be
!tested after the test on nkpt is performed.
!8. nband
!DEBUG
!write(std_out,*)' inprep8 : nkpt=',nkpt
!if(count==2)stop
!ENDDEBUG
 ABI_ALLOCATE(nband,(nkpt))
 if(occopt==2)then
   im=12
   do iline=1,(nkpt+11)/12
     if(iline==(nkpt+11)/12)im=nkpt-12*(iline-1)
     if(ddbvrs_is_current_or_old)then
       read (unddb, '(1x,a9,5x,12i5)' )name(1),&
&       (nband((iline-1)*12+ii),ii=1,im)
     else
       read (unddb, '(1x,a6,5x,12i5)' )name_old,&
&       (nband((iline-1)*12+ii),ii=1,im) ; name(1)='   '//name_old
     end if
     if (iline==1) then
       call chknm8(name(1),'    nband')
     else
       call chknm8(name(1),'         ')
     end if
   end do
 else
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,i10)' )name(1),nband(1)
   else
     read (unddb, '(1x,a6,i10)' )name_old,nband(1) ; name(1)='   '//name_old
   end if
   call chknm8(name(1),'    nband')
   if(nkpt>1)then
     do ikpt=2,nkpt
       nband(ikpt)=nband(1)
     end do
   end if
 end if

!Check all nband values, and sum them
 bantot=0
 do ikpt=1,nkpt
   if(nband(ikpt)<0)then
     write(message, '(a,a,a,i4,a,i4,3a)' )&
&     ' inprep8 : ERROR -',ch10,&
&     '  For ikpt = ',ikpt,'  nband = ',nband(ikpt),' is negative.',ch10,&
&     '  Action : correct your DDB.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   bantot=bantot+nband(ikpt)
 end do

 mband=maxval(nband(:))

!Skip the rest of variables
!9. acell
 read (unddb,*)
!10. amu
 do iline=1,(ntypat+2)/3
   read (unddb,*)
 end do
!11. dilatmx
 if(ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!12. ecut
 read (unddb,*)
!12b. pawecutdg (PAW only)
 if(ddbvrs==vrsio8.and.usepaw==1) then
   read (unddb,*)
 end if
!13. ecutsm
 if(ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!14. intxc
 if(ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!15. iscf
 read (unddb,*)
!16. ixc
 read (unddb,*)
!17. kpt
 do iline=1,nkpt
   read (unddb,*)
 end do
!18. kptnrm
 read (unddb,*)
!19. ngfft
 read (unddb,*)
!20. nspden
 if(ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!21. nspinor
 if(ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!22. occ
 if(occopt==2)then
   do iline=1,(bantot+2)/3
     read (unddb,*)
   end do
 else
   write(std_out,*)' inprep8 : nband(1)=',nband(1)
   do iline=1,(nband(1)+2)/3
     read (unddb,'(a80)')rdstring
     write(std_out,*)trim(rdstring)
   end do
 end if
!23. rprim
 do iline=1,3
   read (unddb,*)
 end do
!24. sciss
 read (unddb,*)
!25. spinat
 if(ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb,*)
   end do
 end if
!26. symafm
 if(ddbvrs_is_current_or_old)then
   do iline=1,(nsym+11)/12
     read (unddb,*)
   end do
 end if
!27. symrel
 do iline=1,nsym
   read (unddb,*)
 end do
!28old. xred
 if(.not.ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb,*)
   end do
 end if
!28. tnons
 do iline=1,nsym
   read (unddb,*)
 end do
!29. tolwfr
 if(ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!30. tphysel
 if(ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!31. tsmear
 if(ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!32. type
 do iline=1,(natom+11)/12
   read (unddb,*)
 end do
!33old. tolwfr
 if(.not.ddbvrs_is_current_or_old)then
   read (unddb,*)
 end if
!33. wtk
 do iline=1,(nkpt+2)/3
   read (unddb,*)
 end do
!34. xred
 if(ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb,*)
   end do
 end if
!35. znucl
 if(ddbvrs_is_current_or_old)then
   do iline=1,(ntypat+2)/3
     read (unddb,*)
   end do
 end if
!36. zion
 do iline=1,(ntypat+2)/3
   read (unddb,*)
 end do

 read (unddb,*)

!Now, take care of the pseudopotentials
 read(unddb, '(a12)' )string

 if(string=='  Descriptio')then

   read (unddb,*)
   if (ddbvrs==vrsio8_old.or.ddbvrs==vrsio8_old_old) then
     read (unddb, '(10x,i3,14x,i3,11x,i3)', iostat=ios )dimekb,lmnmax,usepaw
     if(ios/=0)then
       backspace(unddb)
       read (unddb, '(10x,i3,14x,i3)')dimekb,lmnmax
       usepaw=0
     end if
   else if (ddbvrs==vrsio8) then
     read (unddb, '(10x,i3)') usepaw
     if (usepaw==0) then
       read (unddb, '(10x,i3,14x,i3)' ) dimekb,lmnmax
     else
       dimekb=0;lmnmax=0
     end if
   end if
   if (usepaw==0) then
     do itypat=1,ntypat
       read(unddb, '(13x,i4,9x,i3,8x,i4)' )itypat0,pspso0,nekb
       read(unddb,*)
       do iekb=1,nekb
         do jekb=1,nekb,4
           read(unddb,*)
         end do
       end do
     end do
   else
     do itypat=1,ntypat
       read(unddb, '(12x,i4,12x,i3,12x,i5)' )itypat0,basis_size0,lmn_size0
       lmnmax=max(lmnmax,lmn_size0)
       read(unddb,*)
       read(unddb,*)
       read(unddb,'(24x,i3)') nekb
       read(unddb,*)
       do iekb=1,nekb,4
         read(unddb,*)
       end do
     end do
   end if

 else if(string==' Description')then

   if (usepaw==1) stop 'BUG: old DDB pspformat not compatible with PAW 1'
   read (unddb, '(10x,i3,10x,i3)' )mproj,mpsang
   dimekb=mproj*mpsang
   usepaw=0
   do itypat=1,ntypat
     read (unddb,*)
!    For f-electrons, one more line has been written
     do iproj=1,mproj*max(1,(mpsang+2)/3)
       read (unddb,*)
     end do
   end do

 else if(string==' No informat')then

   dimekb=0
   lmnmax=0
   usepaw=0

 else
   write(message, '(a,a,a,a,a,a)' )&
&   ' inprep8 : BUG -',ch10,&
&   '  Error when reading the psp information',ch10,&
&   '  String=',string
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Now, the number of blocks
 read(unddb,*)
 read(unddb,*)
 read(unddb, '(24x,i4)' )nblok

!Now, the type of each blok, in turn
 mblktyp=1
 if(nblok>=1)then
   do iblok=1,nblok

     read(unddb,*)
     read(unddb, '(a32,12x,i8)' )blkname,nelmts
     if(blkname==' 2nd derivatives (non-stat.)  - ')then
       blktyp=1
     else if(blkname==' 2nd derivatives (stationary) - ')then
       blktyp=2
     else if(blkname==' 3rd derivatives              - ')then
       blktyp=3
     else if(blkname==' Total energy                 - ')then
       blktyp=0
     else if(blkname==' 1st derivatives              - ')then
       blktyp=4
     else if(blkname==' 2nd eigenvalue derivatives   - ')then
       blktyp=5
     else
       write(message, '(a,a,a,a,a,a,a,a,a,a,a)' )&
&       ' inprep8 : ERROR -',ch10,&
&       '  The following string appears in the DDB in place of',&
&       ' the block type description :',ch10,blkname,ch10,&
&       '  Action : check your DDB.',ch10,&
&       ' Note: If you did use an abinit version prior to 6.12 to generate your DDB',&
&       ' pay attention to the change:: 2rd derivatives ==> 2nd derivatives'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     if(blktyp==1.or.blktyp==2)then
!      Read the phonon wavevector
       read(unddb,*)
     else if(blktyp==3)then
!      Read the perturbation wavevectors
       read(unddb,*)
       read(unddb,*)
       read(unddb,*)
       mblktyp=3
     else if(blktyp==5)then
       read(unddb,*)
       mblktyp=5
     end if

!    Read every element
     if(blktyp==5)then
       do ikpt=1,nkpt
         read(unddb,*)
         do iband=1,nband(ikpt)
           read(unddb,*)
           do ii=1,nelmts
             read(unddb,*)
           end do
         end do
       end do
     else
       do ii=1,nelmts
         read(unddb,*)
       end do
     end if

   end do
 end if

 ABI_DEALLOCATE(nband)

!Close the DDB
 close(unddb)

end subroutine inprep8
!!***
