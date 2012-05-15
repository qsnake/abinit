!{\src2tex{textfont=tt}}
!!****f* ABINIT/ioddb8_in
!!
!! NAME
!! ioddb8_in
!!
!! FUNCTION
!! Open Derivative DataBase, then
!! reads or write Derivative DataBase preliminary information.
!!
!! Note : only one processor read or write the DDB.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! character(len=fnlen) filnam: name of input file
!! matom=maximum number of atoms
!! mband=maximum number of bands
!! mkpt=maximum number of special points
!! msym=maximum number of symetries
!! mtypat=maximum number of atom types
!! unddb=unit number for input
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.
!!
!! OUTPUT
!! acell(3)=length scales of primitive translations (bohr)
!! amu(mtypat)=mass of the atoms (atomic mass unit)
!! dilatmx=the maximal dilatation factor
!! ecut=kinetic energy planewave cutoff (hartree)
!! ecutsm=smearing energy for plane wave kinetic energy (Ha)
!! intxc=control xc quadrature
!! iscf=parameter controlling scf or non-scf choice
!! ixc=exchange-correlation choice parameter
!! kpt(3,mkpt)=k point set (reduced coordinates)
!! kptnrm=normalisation of k points
!! natom=number of atoms in the unit cell
!! nband(mkpt)=number of bands at each k point, for each polarization
!! ngfft(18)=contain all needed information about 3D FFT,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkpt=number of k points
!! nspden=number of spin-density components
!! nspinor=number of spinorial components of the wavefunctions
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetry elements in space group
!! ntypat=number of atom types
!! occ(mband*mkpt)=occupation number for each band and k
!! occopt=option for occupancies
!! pawecutdg=cut-off for fine "double grid" used in PAW calculations (unused for NCPP)
!! rprim(3,3)=dimensionless primitive translations in real space
!! sciss=scissor shift (Ha)
!! spinat(3,matom)=initial spin of each atom, in unit of hbar/2
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym)=symmetry operations in real space
!! tnons(3,msym)=nonsymmorphic translations for symmetry operations
!! tolwfr=tolerance on largest wf residual
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature) in Hartree
!! typat(matom)=type of each atom
!! usepaw=flag for PAW
!! wtk(mkpt)=weight assigned to each k point
!! xred(3,matom)=reduced atomic coordinates
!! zion(mtypat)=valence charge of each type of atom
!! znucl(mtypat)=atomic number of atom type
!!
!! TODO
!!
!! PARENTS
!!      mblktyp1,mblktyp5,rdddb9,thmeig
!!
!! CHILDREN
!!      chknm8,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ioddb8_in(filnam,matom,mband,&
&  mkpt,msym,mtypat,unddb,vrsddb,&
&  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  pawecutdg,rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
&  typat,usepaw,wtk,xred,zion,znucl)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ioddb8_in'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: matom,mband,mkpt,msym,mtypat,unddb,vrsddb
 integer,intent(out) :: intxc,iscf,ixc,natom,nkpt,nspden,nspinor,nsppol,nsym
 integer,intent(out) :: ntypat,occopt
 integer,intent(in)  :: usepaw
 real(dp),intent(out) :: dilatmx,ecut,ecutsm,pawecutdg,kptnrm,sciss,tolwfr,tphysel
 real(dp),intent(out) :: tsmear
 character(len=fnlen),intent(in) :: filnam
!arrays
 integer,intent(out) :: nband(mkpt),ngfft(18),symafm(msym),symrel(3,3,msym)
 integer,intent(out) :: typat(matom)
 real(dp),intent(out) :: acell(3),amu(mtypat),kpt(3,mkpt),occ(mband*mkpt)
 real(dp),intent(out) :: rprim(3,3),spinat(3,matom),tnons(3,msym),wtk(mkpt)
 real(dp),intent(out) :: xred(3,matom),zion(mtypat),znucl(mtypat)

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
 integer :: bantot,ddbvrs,iband,ii,ij,ikpt,iline,im,usepaw0
 logical :: ddbvrs_is_current_or_old,testn,testv
 character(len=500) :: message
 character(len=6) :: name_old
!arrays
 character(len=9) :: name(9)

! *********************************************************************

!DEBUG
!write(std_out,*)' ioddb8_in : enter'
!write(std_out,*)' ioddb8_in : mtypat=',mtypat

!ENDDEBUG

!Check ioddb8 version number (vrsio8) against mkddb version number
!(vrsddb)
 if (vrsio8/=vrsddb) then
   write(message, '(a,a,a,i10,a,a,i10,a)' )&
&   ' ioddb8_in: BUG -',ch10,&
&   '  The input/output DDB version number=',vrsio8,ch10,&
&   '  is not equal to the DDB version number=',vrsddb,'.'
   call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
 end if

 write(std_out,'(a,a)')&
& ' About to open file ',filnam
!Open the input derivative database.
 open (unit=unddb,file=filnam,status='old',form='formatted')

!Check the compatibility of the input DDB with the DDB code
 read (unddb,*)
 read (unddb,*)
 read (unddb, '(20x,i10)' )ddbvrs
 write(std_out,'(a,i10)')' ddbvrs=',ddbvrs
 if(ddbvrs/=vrsio8 .and. ddbvrs/=vrsio8_old .and. ddbvrs/=vrsio8_old_old)then
   write(message, '(3a,i10,2a,3(a,i10),a)' )&
&   ' ioddb8_in : BUG - ',ch10,&
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
   read (unddb, '(1x,a9,i10)' )name(1),usepaw0
 else
   usepaw0=0;name(1)='   usepaw'
 end if
 if(name(1)/='   usepaw')testn=.false.
 if(usepaw0/=usepaw)testv=.false.
!2. natom
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(2),natom
 else
   read (unddb, '(1x,a6,i10)' )name_old,natom ; name(2)='   '//name_old
 end if
 if(name(2)/='    natom')testn=.false.
 if(natom<=0.or.natom>matom)testv=.false.
!3. nkpt
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(3),nkpt
 else
   read (unddb, '(1x,a6,i10)' )name_old,nkpt ; name(3)='   '//name_old
 end if
 if(name(3)/='     nkpt')testn=.false.
 if(nkpt <=0.or.nkpt >mkpt )testv=.false.
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
 if(nsym <=0.or.nsym >msym )testv=.false.
!6. ntypat
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(6),ntypat
 else
   read (unddb, '(1x,a6,i10)' )name_old,ntypat ; name(6)='   '//name_old
 end if
 if(name(6)/='   ntypat' .and. name(6)/='    ntype')testn=.false.
 if(ntypat<=0.or.ntypat>mtypat)testv=.false.
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
   write(message, '(a,a,a)' )' ioddb8_in : An error has been found in one',ch10,&
&   ' of the positive n-integers contained in the DDB : '
   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )&
&   '               Expected                      Found     '
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    usepaw equal to   ',usepaw,'    ',trim(name(1)),' =',usepaw0
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    natom , lower than',matom+1,'    ',trim(name(2)),' =',natom
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nkpt  , lower than',mkpt+1 ,'    ',trim(name(3)),' =',nkpt
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nsppol, lower than',3      ,'    ',trim(name(4)),' =',nsppol
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nsym  , lower than',msym+1 ,'    ',trim(name(5)),' =',nsym
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    ntypat, lower than',mtypat+1,'   ',trim(name(6)),' =',ntypat
   call wrtout(std_out,message,'COLL')
   write(message, '(a,a,a,i10)' )&
&   '    occopt,  between 0 and 7        ',trim(name(7)),' =',occopt
   call wrtout(std_out,message,'COLL')
   write(message, '(a,a,a)' )&
&   ' ioddb8_in : ERROR -',ch10,&
&   '  See the error message above.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!One more set of parameters define the dimensions of the
!array : nband. Morever, it depends on occopt and
!nkpt, and has to be
!tested after the test on nkpt is performed.
!8. nband
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

!check all nband values, and sum them
 bantot=0
 do ikpt=1,nkpt
   if(nband(ikpt)<0)then
     write(message, '(3a,i4,a,i4,3a)' )&
&     ' ioddb8_in : ERROR -',ch10,&
&     '  For ikpt = ',ikpt,'  nband = ',nband(ikpt),' is negative.',ch10,&
&     '  Action : correct your DDB.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   else if(nband(ikpt)>mband)then
     write(message, '(3a,i4,a,i4,a,a,i4,3a)' )&
&     ' ioddb8_in : ERROR -',ch10,&
&     ' For ikpt = ',ikpt,', nband = ',nband(ikpt),ch10,&
&     ' is larger than mband = ',mband,'.',ch10,&
&     ' Action : recompile the calling code with a larger mband.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   bantot=bantot+nband(ikpt)
 end do

!Read the rest of variables, with check of the names
!9. acell
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,3d22.14)' )name(1),acell
 else
   read (unddb, '(1x,a6,3d22.14)' )name_old,acell ; name(1)='   '//name_old
 end if
 call chknm8(name(1),'    acell')
!9. amu
 im=3
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (amu((iline-1)*3+ii),ii=1,im)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,&
&     (amu((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call chknm8(name(1),'      amu')
   else
     call chknm8(name(1),'         ')
   end if
 end do
!11. dilatmx
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),dilatmx
   call chknm8(name(1),'  dilatmx')
 else
   dilatmx=one
 end if
!12. ecut
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),ecut
 else
   read (unddb, '(1x,a6,d22.14)' )name_old,ecut ; name(1)='   '//name_old
 end if
 call chknm8(name(1),'     ecut')
!12b. pawecutdg (PAW only)
 if(ddbvrs==vrsio8.and.usepaw==1) then
   read (unddb, '(1x,a9,d22.14)' )name(1),pawecutdg
 else
   pawecutdg=ecut;name(1)='pawecutdg'
 end if
 call chknm8(name(1),'pawecutdg')
!13. ecutsm
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),ecutsm
   call chknm8(name(1),'   ecutsm')
 else
   ecutsm=zero
 end if
!14. intxc
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),intxc
   call chknm8(name(1),'    intxc')
 else
   intxc=1
 end if
!15. iscf
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),iscf
 else
   read (unddb, '(1x,a6,i10)' )name_old,iscf ; name(1)='   '//name_old
 end if
 call chknm8(name(1),'     iscf')
!16. ixc
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),ixc
 else
   read (unddb, '(1x,a6,i10)' )name_old,ixc ; name(1)='   '//name_old
 end if
 call chknm8(name(1),'      ixc')
!17. kpt
 do iline=1,nkpt
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (kpt(ii,iline),ii=1,3)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,&
&     (kpt(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call chknm8(name(1),'      kpt')
   else
     call chknm8(name(1),'         ')
   end if
 end do
!18. kptnrm
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),kptnrm
 else
   read (unddb, '(1x,a6,d22.14)' )name_old,kptnrm ; name(1)='   '//name_old
 end if
 call chknm8(name(1),'   kptnrm')
!19. ngfft
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,5x,3i5)' )name(1),ngfft(1:3)
 else
   read (unddb, '(1x,a6,5x,3i5)' )name_old,ngfft(1:3) ; name(1)='   '//name_old
 end if
!For the time being, do not check the validity of the name,
!in order to accept both ng and ngfft
!20. nspden
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),nspden
   call chknm8(name(1),'   nspden')
 else
   nspden=0
 end if
!21. nspinor
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),nspinor
   call chknm8(name(1),'  nspinor')
 else
   nspinor=0
 end if
!22. occ
 if(occopt==2)then
   im=3
   do iline=1,(bantot+2)/3
     if(iline==(bantot+2)/3)im=bantot-3*(iline-1)
     if(ddbvrs_is_current_or_old)then
       read (unddb, '(1x,a9,3d22.14)' )name(1),&
&       (occ((iline-1)*3+ii),ii=1,im)
     else
       read (unddb, '(1x,a6,3d22.14)' )name_old,&
&       (occ((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
     end if
     if (iline==1) then
       call chknm8(name(1),'      occ')
     else
       call chknm8(name(1),'         ')
     end if
   end do
 else
   im=3
   do iline=1,(nband(1)+2)/3
     if(iline==(nband(1)+2)/3)im=nband(1)-3*(iline-1)
     if(ddbvrs_is_current_or_old)then
       read (unddb, '(1x,a9,3d22.14)' )name(1),&
&       (occ((iline-1)*3+ii),ii=1,im)
     else
       read (unddb, '(1x,a6,3d22.14)' )name_old,&
&       (occ((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
     end if
     if (iline==1) then
       call chknm8(name(1),'      occ')
     else
       call chknm8(name(1),'         ')
     end if
   end do
   if(nkpt>1)then
     do ikpt=2,nkpt
       do iband=1,nband(1)
         occ(iband+nband(1)*(ikpt-1))=occ(iband)
       end do
     end do
   end if
 end if
!23. rprim
 do iline=1,3
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (rprim(ii,iline),ii=1,3)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,&
&     (rprim(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call chknm8(name(1),'    rprim')
   else
     call chknm8(name(1),'         ')
   end if
 end do
!24. sciss
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),sciss
 else
   read (unddb, '(1x,a6,d22.14)' )name_old,sciss ; name(1)='   '//name_old
 end if
 call chknm8(name(1),'    sciss')
!25. spinat
 if(ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (spinat(ii,iline),ii=1,3)
     if (iline==1) then
       call chknm8(name(1),'   spinat')
     else
       call chknm8(name(1),'         ')
     end if
   end do
 else
!  spinat is set to zero by default in mrgddb.f
!  spinat(:,1:natom)=zero
 end if
!26. symafm
 if(ddbvrs_is_current_or_old)then
   im=12
   do iline=1,(nsym+11)/12
     if(iline==(nsym+11)/12)im=nsym-12*(iline-1)
     read (unddb, '(1x,a9,5x,12i5)' )name(1),&
&     (symafm((iline-1)*12+ii),ii=1,im)
     if (iline==1) then
       call chknm8(name(1),'   symafm')
     else
       call chknm8(name(1),'         ')
     end if
   end do
 else
!  symafm is set to 1 by default in mrgddb.f
!  symafm(1:nsym)=1
 end if
!27. symrel
 do iline=1,nsym
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,5x,9i5)' )name(1),&
&     ((symrel(ii,ij,iline),ii=1,3),ij=1,3)
   else
     read (unddb, '(1x,a6,5x,9i5)' )name_old,&
&     ((symrel(ii,ij,iline),ii=1,3),ij=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call chknm8(name(1),'   symrel')
   else
     call chknm8(name(1),'         ')
   end if
 end do
!28old. xred
 if(.not.ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb, '(1x,a6,3d22.14)' )name(1),&
&     (xred(ii,iline),ii=1,3)
   end do
!  No check of name, to allow the old tn
 end if
!28. tnons
 do iline=1,nsym
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (tnons(ii,iline),ii=1,3)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,&
&     (tnons(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call chknm8(name(1),'    tnons')
   else
     call chknm8(name(1),'         ')
   end if
 end do
!29. tolwfr
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tolwfr
 end if
!Do not check the name, in order to allow both tolwfr and wftol
!30. tphysel
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tphysel
   call chknm8(name(1),'  tphysel')
 else
   tphysel=zero
 end if
!31. tsmear
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tsmear
   call chknm8(name(1),'   tsmear')
 else
   tsmear=zero
 end if
!32. typat
 im=12
 do iline=1,(natom+11)/12
   if(iline==(natom+11)/12)im=natom-12*(iline-1)
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,5x,12i5)' )name(1),&
&     (typat((iline-1)*12+ii),ii=1,im)
   else
     read (unddb, '(1x,a6,5x,12i5)' )name_old,&
&     (typat((iline-1)*12+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
!    Both type and typat are allowed => no check
!    call chknm8(name(1),'    typat')
   else
     call chknm8(name(1),'         ')
   end if
 end do
!33old. tolwfr
 if(.not.ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a6,d22.14)' )name(1),tolwfr
 end if
!Do not check the name, in order to allow both tolwfr and wftol
!33. wtk
 im=3
 do iline=1,(nkpt+2)/3
   if(iline==(nkpt+2)/3)im=nkpt-3*(iline-1)
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (wtk((iline-1)*3+ii),ii=1,im)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,&
&     (wtk((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call chknm8(name(1),'      wtk')
   else
     call chknm8(name(1),'         ')
   end if
 end do
!34. xred
 if(ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (xred(ii,iline),ii=1,3)
     if (iline==1) then
       call chknm8(name(1),'     xred')
     else
       call chknm8(name(1),'         ')
     end if
   end do
 end if
!35. znucl
 if(ddbvrs_is_current_or_old)then
   im=3
   do iline=1,(ntypat+2)/3
     if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (znucl((iline-1)*3+ii),ii=1,im)
     if (iline==1) then
       call chknm8(name(1),'    znucl')
     else
       call chknm8(name(1),'         ')
     end if
   end do
 else
!  znucl is set to zero by default in mrgddb.f
!  znucl(:)=zero
 end if
!36. zion
 im=3
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (zion((iline-1)*3+ii),ii=1,im)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,&
&     (zion((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
!    Do not check the names, to allow both zion and znucl - the latter for 990527 format
!    call chknm8(name(1),'     zion')
   else
     call chknm8(name(1),'         ')
   end if
 end do

end subroutine ioddb8_in
!!***
