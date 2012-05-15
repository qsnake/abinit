!{\src2tex{textfont=tt}}
!!****f* ABINIT/rdddb9
!!
!! NAME
!! rdddb9
!!
!! FUNCTION
!! This routine reads the derivative database entirely,
!! for use in ppddb9, and performs some checks and symmetrisation
!! At the end, the whole DDB is in central memory, contained in the
!! array ddb_blk%val(2,msize,ddb_blk%nblok).
!! The information on it is contained in the four arrays
!!   ddb_blk%flg(msize,ddb_blk%nblok) : blok flag for each element
!!   ddb_blk%qpt(9,ddb_blk%nblok)  : blok wavevector (unnormalized)
!!   ddb_blk%nrm(3,ddb_blk%nblok)  : blok wavevector normalization
!!   ddb_blk%typ(ddb_blk%nblok)    : blok type
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! atifc(natom) = atifc(ia) equals 1 if the analysis of ifc
!!  has to be done for atom ia; otherwise 0
!! ddbun = unit number for DDB io
!! dimekb=dimension of ekb (for the time being, only for norm-
!!                          conserving psps)
!! iout=unit number for output of formatted data
!! filnam=name of input file
!! lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!       =if useylm=0, max number of (l,n)   comp. over all type of psps
!! mband=maximum number of bands
!! mpert =maximum number of ipert
!! msize=maximum size of data blocks
!! msym =maximum number of symmetry elements in space group
!! natifc = number of atoms for which the analysis of ifc is done
!! natom = number of atoms
!! ntypat=number of atom types
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!! acell(3)=length scales of cell (bohr)
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! ddb_blk : ddb blok datatype
!!   contents: ddb_blk%flg(msize,nblok)= flag of existence for each element of the DDB
!!             ddb_blk%nrm(3,nblok)  : blok wavevector normalization
!!             ddb_blk%qpt(9,nblok)  : blok wavevector (unnormalized)
!!             ddb_blk%typ(nblok)    : blok type
!!             ddb_blk%val(2,msize,nblok)= value of each complex element of the DDB
!!             ddb_blk%nblok= number of bloks in the DDB
!! gmet(3,3)=reciprocal space metric tensor in bohr**-2
!! gprim(3,3)=dimensionless reciprocal space primitive translations
!! indsym(4,nsym,natom)=indirect indexing array for symmetries
!! natom=number of atoms in cell
!! nsym=number of space group symmetries
!! occopt=occupation option
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! rprim(3,3)= primitive translation vectors
!! symq(4,2,nsym)= (integer) three first numbers define the G vector ;
!!   fourth number is zero if the q-vector is not preserved,
!!   second index is about time-reversal symmetry
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! tnons(3,nsym)=fractional nonsymmorphic translations
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! xcart(3,natom)=atomic cartesian coordinates
!! xred(3,natom)=fractional dimensionless atomic coordinates
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      cart29,chkin9,d2sym3,d3sym,destroy_pawtab,ioddb8_in,mati3inv,matr3inv
!!      metric,mkrdim,nlopt,nullify_pawtab,psddb8,read_blok8,symatm,symq3
!!      timein,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rdddb9(acell,atifc,amu,ddb_blk,&
& ddbun,dimekb,filnam,gmet,gprim,indsym,iout,&
& lmnmax,mband,mpert,msize,msym,&
& natifc,natom,nkpt,nsym,ntypat,&
& occopt,rmet,rprim,symq,symrec,symrel,&
& tnons,typat,ucvol,usepaw,xcart,xred,zion)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_ddb_blk

 use m_paw_toolbox, only : destroy_pawtab,nullify_pawtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdddb9'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => rdddb9
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! NOTE: these are used for dimensioning and then re-assigned in ioddb8.
!   This is almost definitely bad practice. In particular
!    it should be indsym(4,msym,natom),
!   and
!    the allocation allocate(kpt(3,nkpt)) is strange
!scalars
 integer,intent(in) :: ddbun,dimekb,iout,lmnmax,mband,mpert,msize,msym,natifc
 integer,intent(inout) :: natom,nkpt,nsym,ntypat,occopt,usepaw
 real(dp),intent(out) :: ucvol
 character(len=fnlen),intent(in) :: filnam
 type(ddb_blk_type), pointer :: ddb_blk
!arrays
 integer,intent(inout) :: atifc(natom)
 integer,intent(out) :: indsym(4,nsym,natom)
 integer,intent(out) :: symq(4,2,*),symrec(3,3,msym),symrel(3,3,msym)
 integer,intent(out) :: typat(natom)
 real(dp),intent(out) :: acell(3),amu(ntypat)
 real(dp),intent(out) :: gmet(3,3),gprim(3,3),rmet(3,3)
 real(dp),intent(out) :: rprim(3,3),tnons(3,msym),xcart(3,natom),xred(3,natom)
 real(dp),intent(out) :: zion(ntypat)

!Local variables -------------------------
!mtyplo=maximum number of type, locally
!scalars
 integer,parameter :: msppol=2,mtyplo=6
 integer :: choice,fullinit,iblok,intxc,iscf,isym,ixc
 integer :: nsize,nspden,nspinor,nsppol,nunit,timrev,useylm,vrsddb
 real(dp) :: dilatmx,ecut,ecutsm,kptnrm,pawecutdg,sciss,tcpui,tolsym8,tolwfr
 real(dp) :: tphysel,tsmear,twalli
 character(len=500) :: message
!arrays
 integer :: ngfft(18)
 integer,allocatable :: car3flg(:,:,:,:,:,:),carflg(:,:,:,:),indlmn(:,:,:)
 integer,allocatable :: nband(:),pspso(:),symafm(:),tmpflg(:,:,:,:,:,:)
 real(dp) :: gprimd(3,3),qpt(3),rprimd(3,3)
 real(dp),allocatable :: d2cart(:,:,:,:,:),d3cart(:,:,:,:,:,:,:),ekb(:,:)
 real(dp),allocatable :: kpt(:,:),occ(:),spinat(:,:),tmpval(:,:,:,:,:,:,:)
 real(dp),allocatable :: wtk(:),znucl(:)
 type(pawtab_type),allocatable :: pawtab(:)

! *********************************************************************

 DBG_ENTER("COLL")

 call timein(tcpui,twalli)

!Read the DDB information

 vrsddb=100401

!The checking of pseudopotentials is not done presently
!so that dimensions are fake
 ABI_ALLOCATE(ekb,(dimekb,ntypat))
 ABI_ALLOCATE(indlmn,(6,lmnmax,ntypat))
 ABI_ALLOCATE(pspso,(ntypat))
 ABI_ALLOCATE(pawtab,(ntypat*usepaw))
 call nullify_pawtab(pawtab)

 ABI_ALLOCATE(kpt,(3,nkpt))
 ABI_ALLOCATE(nband,(nkpt))
 ABI_ALLOCATE(occ,(nkpt*mband*msppol))
 ABI_ALLOCATE(spinat,(3,natom))
 ABI_ALLOCATE(symafm,(msym))
 ABI_ALLOCATE(wtk,(nkpt))
 ABI_ALLOCATE(znucl,(ntypat))

!Open the input derivative database file
!and read the preliminary information
 nunit=ddbun
!Note that in this call, mkpt has been replaced by nkpt,
!mtypat by ntypat, and matom by natom.
 call ioddb8_in(filnam,natom,mband,&
& nkpt,msym,ntypat,nunit,vrsddb,&
& acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
& natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
& pawecutdg,rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
& typat,usepaw,wtk,xred,zion,znucl)

!Compute different matrices in real and reciprocal space, also
!checks whether ucvol is positive.
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

!Obtain reciprocal space primitive transl g from inverse trans of r
!(Unlike in abinit, gprim is used throughout ifc; should be changed, later)
 call matr3inv(rprim,gprim)

!Generate atom positions in cartesian coordinates
 call xredxcart(natom,1,rprimd,xcart,xred)

!Transposed inversion of the symmetry matrices, for use in
!the reciprocal space
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

!SYMATM generates for all the atoms and all the symmetries, the atom
!on which the referenced one is sent and also the translation bringing
!back this atom to the referenced unit cell
 tolsym8=tol8
 call symatm(indsym,natom,nsym,symrec,tnons,tolsym8,typat,xred)

!Read the psp information of the input DDB
 useylm=usepaw;choice=1
 call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
& ddb_blk%nblok,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

!Check the correctness of some input parameters,
!and perform small treatment if needed.
 call chkin9(atifc,natifc,natom)

!Read the blocks from the input database, and close it.
 write(message, '(a,a,a,i5,a)' )ch10,ch10,&
& ' rdddb9 : read ',ddb_blk%nblok,' blocks from the input DDB '
 call wrtout(std_out,message,'COLL')
 nunit=ddbun
 do iblok=1,ddb_blk%nblok
!  index=1+nsize*(iblok-1)
   call read_blok8(ddb_blk,iblok,mband,mpert,msize,nkpt,nunit)

!  Here complete the matrix by symmetrisation of the
!  existing elements
   if(ddb_blk%typ(iblok)==1 .or. ddb_blk%typ(iblok)==2) then
     qpt(1)=ddb_blk%qpt(1,iblok)/ddb_blk%nrm(1,iblok)
     qpt(2)=ddb_blk%qpt(2,iblok)/ddb_blk%nrm(1,iblok)
     qpt(3)=ddb_blk%qpt(3,iblok)/ddb_blk%nrm(1,iblok)

!    Examine the symmetries of the q wavevector
     call symq3(nsym,qpt,symq,symrec,timrev,prtvol=0)

     nsize=3*mpert*3*mpert
     ABI_ALLOCATE(tmpflg,(3,mpert,3,mpert,1,1))
     ABI_ALLOCATE(tmpval,(2,3,mpert,3,mpert,1,1))
     tmpflg(:,:,:,:,1,1) = reshape(ddb_blk%flg(1:nsize,iblok),&
&     shape = (/3,mpert,3,mpert/))
     tmpval(1,:,:,:,:,1,1) = reshape(ddb_blk%val(1,1:nsize,iblok),&
&     shape = (/3,mpert,3,mpert/))
     tmpval(2,:,:,:,:,1,1) = reshape(ddb_blk%val(2,1:nsize,iblok),&
&     shape = (/3,mpert,3,mpert/))

!    Then apply symmetry operations
     call d2sym3(tmpflg,tmpval,indsym,mpert,&
&     natom,nsym,qpt,symq,symrec,symrel,timrev)

!    Transform the dynamical matrix in cartesian coordinates
     ABI_ALLOCATE(carflg,(3,mpert,3,mpert))
     ABI_ALLOCATE(d2cart,(2,3,mpert,3,mpert))
     call cart29(tmpflg,tmpval,carflg,d2cart,&
&     gprimd,1,mpert,natom,1,ntypat,rprimd,typat,ucvol,zion)

     ddb_blk%flg(1:nsize,iblok) = reshape(carflg,shape = (/3*mpert*3*mpert/))
     ddb_blk%val(1,1:nsize,iblok) = reshape(d2cart(1,:,:,:,:),&
&     shape = (/3*mpert*3*mpert/))
     ddb_blk%val(2,1:nsize,iblok) = reshape(d2cart(2,:,:,:,:),&
&     shape = (/3*mpert*3*mpert/))


     ABI_DEALLOCATE(carflg)
     ABI_DEALLOCATE(d2cart)
     ABI_DEALLOCATE(tmpflg)
     ABI_DEALLOCATE(tmpval)

   else if (ddb_blk%typ(iblok) == 3) then

     nsize=3*mpert*3*mpert*3*mpert
     ABI_ALLOCATE(tmpflg,(3,mpert,3,mpert,3,mpert))
     ABI_ALLOCATE(tmpval,(2,3,mpert,3,mpert,3,mpert))

     tmpflg(:,:,:,:,:,:) = reshape(ddb_blk%flg(1:nsize,iblok),&
&     shape = (/3,mpert,3,mpert,3,mpert/))
     tmpval(1,:,:,:,:,:,:) = reshape(ddb_blk%val(1,1:nsize,iblok),&
&     shape = (/3,mpert,3,mpert,3,mpert/))
     tmpval(2,:,:,:,:,:,:) = reshape(ddb_blk%val(2,1:nsize,iblok),&
&     shape = (/3,mpert,3,mpert,3,mpert/))

     call d3sym(tmpflg,tmpval,indsym,mpert,natom,nsym,&
&     symrec,symrel)

     ABI_ALLOCATE(d3cart,(2,3,mpert,3,mpert,3,mpert))
     ABI_ALLOCATE(car3flg,(3,mpert,3,mpert,3,mpert))
     call nlopt(tmpflg,car3flg,tmpval,d3cart,gprimd,mpert,natom,rprimd,ucvol)

     ddb_blk%flg(1:nsize,iblok) = reshape(car3flg, shape = (/3*mpert*3*mpert*3*mpert/))
     ddb_blk%val(1,1:nsize,iblok) = reshape(d3cart(1,:,:,:,:,:,:),&
&     shape = (/3*mpert*3*mpert*3*mpert/))
     ddb_blk%val(2,1:nsize,iblok) = reshape(d3cart(2,:,:,:,:,:,:),&
&     shape = (/3*mpert*3*mpert*3*mpert/))


     ABI_DEALLOCATE(d3cart)
     ABI_DEALLOCATE(car3flg)
     ABI_DEALLOCATE(tmpflg)
     ABI_DEALLOCATE(tmpval)
   end if

 end do

 close(ddbun)

 write(message,'(a)' )' Now the whole DDB is in central memory '
 call wrtout(std_out,message,'COLL')
 call wrtout(iout,message,'COLL')

 ABI_DEALLOCATE(ekb)
 ABI_DEALLOCATE(indlmn)
 ABI_DEALLOCATE(kpt)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(pspso)
 ABI_DEALLOCATE(spinat)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(znucl)
 call destroy_pawtab(pawtab)
 ABI_DEALLOCATE(pawtab)

 DBG_EXIT("COLL")

end subroutine rdddb9
!!***
