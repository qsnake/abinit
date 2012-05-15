!{\src2tex{textfont=tt}}
!!****f* ABINIT/thmeig
!! NAME
!! thmeig
!!
!! FUNCTION
!! This routine calculates the thermal corrections to the eigenvalues.
!! The output is this quantity for the input k point.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (PB, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  filnam = root filename for outputs
!!  filnam5 = name of the eig2 database file
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      asria_corr,carttransf,chkin9,create_ddb_blk,destroy_ddb_blk
!!      destroy_pawtab,destroy_tetra,get_tetra_weight,getkgrid,gtblk9
!!      init_tetra,ioddb8_in,leave_new,mati3inv,matr3inv,metric,mkrdim
!!      nullify_ddb_blk,nullify_pawtab,nullify_tetra,outg2f,outphdos,phfrq3
!!      psddb8,read_blok8,sort_dp,symatm,symfind,symlatt,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine thmeig(g2fsmear,acell,amu,anaddb_dtset,d2asr,&
& filnam,mband,mpert,msize,natom,nkpt,ntemper,&
& ntypat,rprim,telphint,temperinc,&
& tempermin,thmflag,typat,xred,&
& ddb_blk,ddbun,dimekb,filnam5,iout,& !new
& lmnmax,msym,nblok2,nsym,occopt,symrel,tnons,usepaw,zion,&
& symrec,natifc,gmet,gprim,indsym,rmet,atifc,ucvol,xcart) !new

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_paw_toolbox, only : destroy_pawtab,nullify_pawtab
 use m_tetrahedron
 use m_errors
 use m_ddb_blk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'thmeig'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_56_recipspace
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => thmeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mpert,msize,ntemper,telphint,thmflag
 real(dp),intent(in) :: g2fsmear,temperinc,tempermin 
 character(len=fnlen),intent(in) :: filnam
 real(dp),intent(out) :: ucvol !new
 integer,intent(in) :: ddbun,dimekb,iout,lmnmax,msym !new
 integer,intent(in) :: usepaw,natifc !new
 character(len=fnlen),intent(in) :: filnam5 !new
 integer,intent(inout) :: natom,nkpt,nsym,ntypat,occopt,nblok2 !new in ==> inout
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
 ! no intent here - should be inout but may cause problems with ifort 10.0
 type(ddb_blk_type), pointer :: ddb_blk
!arrays
 real(dp),intent(inout) :: d2asr(2,3,natom,3,natom)
 integer,intent(out) :: symrel(3,3,msym)
 integer,intent(out) :: indsym(4,nsym,natom),symrec(3,3,msym) !new
 integer,intent(inout) :: typat(natom),atifc(natom)! new in ==> inout
 real(dp),intent(out) :: zion(ntypat),tnons(3,msym),gmet(3,3) !new
 real(dp),intent(out) :: gprim(3,3),rmet(3,3),xcart(3,natom) !new
 real(dp),intent(inout) :: acell(3),amu(ntypat),rprim(3,3),xred(3,natom)! new in ==> inout



!Local variables-------------------------------
!scalars
 integer :: brav,chksymbreak,found,gqpt,iatom1,iatom2,iband,iblok,iblok2,idir1,idir2,ii,ikpt,ilatt,imod,index
 integer :: iomega,iqpt,iqpt1,iqpt2,iqpt2_previous,iqpt3,iscf_fake,itemper
 integer :: mpert_eig2,msize2,nene,ng2f,nqshft,nsym_new,unit_g2f,nqpt,nqpt_computed,qptopt,rftyp
!integer :: mqpt,nqpt2,option
 integer :: unit_phdos,unitout
 integer :: choice,fullinit,intxc,iscf,isym,ixc !new
 integer :: nspden,nspinor,nsppol,nptsym,nunit,use_inversion,useylm,vrsddb !new
 integer :: ierr
 integer,parameter :: msppol=2 !new
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: bosein,deltaene,det,domega,enemax,enemin,fact2i,fact2r,factr
 real(dp) :: gaussfactor,gaussprefactor,gaussval,invdet,omega,omega_max,omega_min,qnrm,qptrlen
 real(dp) :: rcvol,tmp,tol,vec1i,vec1r,vec2i,vec2r,veci,vecr,xx
 real(dp) :: tphysel,tolwfr,sciss,tsmear,pawecutdg,kptnrm,ecut,ecutsm,dilatmx !new
 real(dp) :: tolsym,tolsym8  !new
 character(len=500) :: message
 character(len=fnlen) :: outfile
 type(ddb_blk_type), pointer :: ddb_blk_eig2

!arrays
 integer :: ngqpt(9),qptrlatt(3,3),rfelfd(4),rfphon(4),rfstrs(4),vacuum(3)
 integer :: bravais(11)
 integer,allocatable :: indqpt(:)
 integer,allocatable :: symafm(:),symafm_new(:),pspso(:),nband(:),indlmn(:,:,:) !new
 integer, allocatable :: carflg_eig2(:,:,:,:)
 integer, allocatable :: ptsymrel(:,:,:),symrel_new(:,:,:)
!integer, allocatable :: symrec_new(:,:,:)
 real(dp) :: deigi(mband,nkpt)
 real(dp) :: deigr(mband,nkpt),diff_qpt(3),dwtermi(mband,nkpt),dwtermr(mband,nkpt)
 real(dp) :: gprimd(3,3),mesh(3,3),multi(mband,nkpt),multr(mband,nkpt)
 real(dp) :: qlatt(3,3),qphnrm(3),qpt_search(3,3)
 real(dp) :: rprimd(3,3),shiftq(3,8),slope(2,mband,nkpt),tempqlatt(3),thmeigen(2,mband,nkpt)
 real(dp) :: zeropoint(2,mband,nkpt)!new
 real(dp),allocatable :: displ(:)
 real(dp),allocatable :: dos_phon(:),dtweightde(:,:),d2cart(:,:)
 real(dp),allocatable :: eigvec(:,:,:,:),eigval(:,:),g2f(:,:,:),intweight(:,:,:)
 real(dp),allocatable :: indtweightde(:,:,:),tmpg2f(:,:,:),tmpphondos(:),total_dos(:),tweight(:,:)
 real(dp),allocatable :: phfreq(:,:),spinat(:,:),wtk(:),occ(:),znucl(:),kpt(:,:),ekb(:,:) !new
 real(dp),allocatable :: blkval2(:,:,:,:),blkval2gqpt(:,:,:,:),kpnt(:,:,:)
 real(dp),allocatable :: dedni(:,:,:,:),dednr(:,:,:,:)
 real(dp),allocatable :: eigen_in(:)
 real(dp),allocatable :: qpt_full(:,:),qptnrm(:)
!real(dp),allocatable :: qpt2(:,:)
 real(dp),allocatable :: spqpt(:,:),tnons_new(:,:)
 real(dp),allocatable :: wghtq(:)
!real(dp),allocatable :: wtq_folded(:)

 type(pawtab_type),allocatable :: pawtab(:) !new
 integer :: ngfft(18) !new

 type(t_tetrahedron) :: tetrahedra
 character(len=80) :: errstr

! *********************************************************************

!DEBUG
 write(std_out,*)'-thmeig : enter '
!call flush(6)
!ENDDEBUG

 write(message,'(83a)') ch10,('=',ii=1,80),ch10,&
& ' Computation of the electron-phonon changes to the electronic eigenenergies '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!=========================================================================
!0) Initializations
!=========================================================================

!Current version. This is needed for ioddb8_in (see mrgddb.F90)
 vrsddb=100401
!The checking of pseudopotentials is not done presently
!so that dimensions are fake
 ABI_ALLOCATE(ekb,(dimekb,ntypat))
 ABI_ALLOCATE(indlmn,(6,lmnmax,ntypat))
 ABI_ALLOCATE(pspso,(ntypat))
 ABI_ALLOCATE(pawtab,(ntypat*usepaw))
 call nullify_pawtab(pawtab)

 ABI_ALLOCATE(ddb_blk_eig2,)
 call nullify_ddb_blk(ddb_blk_eig2)

 call nullify_tetra(tetrahedra)

 ABI_ALLOCATE(kpt,(3,nkpt))
 ABI_ALLOCATE(nband,(nkpt))
 ABI_ALLOCATE(occ,(nkpt*mband*msppol))
 ABI_ALLOCATE(spinat,(3,natom))
 ABI_ALLOCATE(symafm,(msym))
 ABI_ALLOCATE(wtk,(nkpt))
 ABI_ALLOCATE(znucl,(ntypat))

!At present, only atom-type perturbations are allowed for eig2 type matrix elements.
 mpert_eig2=natom
 msize2=3*mpert_eig2*3*mpert_eig2

 call create_ddb_blk(msize2, nblok2, ddb_blk_eig2)

 ABI_ALLOCATE(blkval2,(2,msize2,mband,nkpt))
 ABI_ALLOCATE(blkval2gqpt,(2,msize2,mband,nkpt))

 ABI_ALLOCATE(eigvec,(2,3,natom,3*natom))
 ABI_ALLOCATE(phfreq,(3*natom,ddb_blk%nblok))

!Open Derivative DataBase then r/w Derivative DataBase preliminary information.    

 write(std_out, '(a)' )  '- thmeig: Initialize the second-order electron-phonon file with name :'
 write(std_out, '(a,a)' )'-         ',trim(filnam5)

 nunit=ddbun
 call ioddb8_in(filnam5,natom,mband,&
& nkpt,msym,ntypat,nunit,vrsddb,&
& acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
& natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
& pawecutdg,rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
& typat,usepaw,wtk,xred,zion,znucl)

!DEBUG
!write(std_out,*)"after first ioddb8"
!ENDDEBUG

!Compute different matrices in real and reciprocal space, also
!checks whether ucvol is positive.
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

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
& nblok2,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

!Check the correctness of some input parameters,
!and perform small treatment if needed.
 call chkin9(atifc,natifc,natom)

 blkval2gqpt(:,:,:,:)=zero

 ABI_ALLOCATE(carflg_eig2,(3,mpert_eig2,3,mpert_eig2))
 ABI_ALLOCATE(kpnt,(3,nkpt,1))

!DEBUG
!write(std_out,*)'-thmeig : 1 '
!write(std_out,*)' nblok2=',nblok2
!write(std_out,*)' thmflag=',thmflag
!call flush(6)
!ENDDEBUG

!=========================================================================
!1) Take care of the Gamma point for thmflag=3, 5 or 7
!=========================================================================

 if(thmflag==3 .or. thmflag==5 .or. thmflag==7) then
   found=0
   do iblok2=1,nblok2

!    DEBUG
!    write(std_out,*)'-thmeig : 1a '
!    write(std_out,*)' iblok2=',iblok2
!    call flush(6)
!    ENDDEBUG

     choice=1
     nunit=ddbun

     call read_blok8(ddb_blk_eig2,iblok2,mband,mpert_eig2,msize2,&
&     nkpt,nunit,blkval2(:,:,:,:),kpnt(:,:,1))

     write (200,*) 'blkval2 in thmeig'
     write (200,*) blkval2
!    DEBUG
!    write(std_out,*)'-thmeig : iblok2,ddb_blk_eig2%typ(iblok2)=',iblok2,ddb_blk_eig2%typ(iblok2)
!    call flush(6)
!    ENDDEBUG
     
     qnrm = ddb_blk_eig2%qpt(1,iblok2)*ddb_blk_eig2%qpt(1,iblok2)+ &
&     ddb_blk_eig2%qpt(2,iblok2)*ddb_blk_eig2%qpt(2,iblok2)+ &
&     ddb_blk_eig2%qpt(3,iblok2)*ddb_blk_eig2%qpt(3,iblok2)
     if(qnrm < qtol) then
       blkval2gqpt(:,:,:,:) = blkval2(:,:,:,:)
       gqpt=iblok2
       write(std_out,*)'-thmeig: found Gamma point in EIG2 DDB, blok number ',iblok2
       found=1
       exit
     end if
   end do

   if(found==0)then
     write(message,'(4a,i3,2a)') ch10,&
&     ' thmeig : ERROR -',ch10,&
&     '   Was unable to find the blok for Gamma point in EIG2 DDB file, while thmflag= ',thmflag,ch10,&
&     '  Action : compute the contribution from Gamma, and merge it in your EIG2 DDB file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Put blkval2gqpt in cartesian coordinates
   call carttransf(ddb_blk_eig2%flg,blkval2gqpt,carflg_eig2,gprimd,gqpt,mband,&
&   mpert_eig2,msize2,natom,nblok2,nkpt,rprimd)

 end if

 close(ddbun)

!DEBUG
!write(std_out,*)"blkval2gqpt=",blkval2gqpt(1,1,1,1)
!ENDDEBUG

!DEBUG
!write(std_out,*)'-thmeig : 2 '
!call flush(6)
!ENDDEBUG

!=========================================================================
!2) Calculation of dE(n,k)/dn(Q,j) : consider all q and modes
!=========================================================================

 if(thmflag==3 .or. thmflag==4)then


!  Use the first list of q wavevectors
   nqpt=anaddb_dtset%nph1l
   ABI_ALLOCATE(spqpt,(3,nqpt))
   do iqpt=1,anaddb_dtset%nph1l
     spqpt(:,iqpt)=anaddb_dtset%qph1l(:,iqpt)/anaddb_dtset%qnrml1(iqpt)
   end do
   ABI_ALLOCATE(wghtq,(nqpt))
   wghtq(:)=one/nqpt

 else if(thmflag>=5 .and. thmflag<=8)then

!  Generates the q point grid
   ngqpt(1:3)=anaddb_dtset%ngqpt(1:3)
   nqshft=anaddb_dtset%nqshft
   qptrlatt(:,:)=0
   qptrlatt(1,1)=ngqpt(1)
   qptrlatt(2,2)=ngqpt(2)
   qptrlatt(3,3)=ngqpt(3)

   ABI_ALLOCATE(ptsymrel,(3,3,msym))
   ABI_ALLOCATE(symafm_new,(msym))
   ABI_ALLOCATE(symrel_new,(3,3,msym))
   ABI_ALLOCATE(tnons_new,(3,msym))
   if(thmflag==7 .or. thmflag==8) then
!    Re-generate symmetry operations from the lattice and atomic coordinates
     tolsym=tol8
     call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
     use_inversion=1
     call symfind(0,(/zero,zero,zero/),gprimd,0,msym,natom,0,nptsym,nsym_new,&
&     ptsymrel,spinat,symafm_new,symrel_new,tnons_new,tolsym,typat,use_inversion,xred)
     write(std_out,*)' thmeig : found ',nsym_new,' symmetries ',ch10
     qptopt=1
   else
     nsym_new=1
     symrel_new(:,:,1)=0 ; symrel_new(1,1,1)=1 ; symrel_new(2,2,1)=1 ; symrel_new(3,3,1)=1
     tnons_new(:,1)=zero
     symafm_new(1)=1
     qptopt=3
   end if

   brav=anaddb_dtset%brav

   if(brav/=1)then
     write(message,'(4a)') ch10,&
&     ' thmeig : ERROR -',ch10,&
&     '   The possibility to have brav/=1 for thmeig was disabled.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  ----NEW CODING
!  Prepare to compute the q-point grid in the ZB or IZB
   iscf_fake=5 ! Need the weights
   chksymbreak=0 
   vacuum=0
   shiftq(:,1:nqshft)=anaddb_dtset%q1shft(:,1:nqshft)
!  Compute the final number of q points
   call getkgrid(chksymbreak,0,iscf_fake,spqpt,qptopt,qptrlatt,qptrlen,&
&   nsym_new,0,nqpt,nqshft,nsym_new,rprimd,&
&   shiftq,symafm_new,symrel_new,vacuum,wghtq)
   ABI_ALLOCATE(spqpt,(3,nqpt))
   ABI_ALLOCATE(wghtq,(nqpt))
   call getkgrid(chksymbreak,iout,iscf_fake,spqpt,qptopt,qptrlatt,qptrlen,&
&   nsym_new,nqpt,nqpt_computed,nqshft,nsym_new,rprimd,&
&   shiftq,symafm_new,symrel_new,vacuum,wghtq)


!  -----OLD CODING
!  call chkrp9(brav,rprim)
!  option=1
!  mqpt=ngqpt(1)*ngqpt(2)*ngqpt(3)*nqshft
!  if(brav==2)mqpt=mqpt/2
!  if(brav==3)mqpt=mqpt/4
!  allocate(spqpt(3,mqpt))
!  call smpbz(brav,iout,qptrlatt,mqpt,nqpt,nqshft,option,anaddb_dtset%q1shft,spqpt)
!  allocate(wghtq(nqpt))
!  wghtq(:)=one/nqpt

!  write(std_out,*)' after smpbz, nqpt=',nqpt
!  do iqpt=1,nqpt
!  write(std_out,*)iqpt,spqpt(:,iqpt),wghtq(iqpt)
!  end do

!  if(thmflag==7 .or. thmflag==8)then

!  Fold the q point set inside de IBZ
!  allocate(indqpt(nqpt),symrec_new(3,3,nsym_new),wtq_folded(nqpt),qpt2(3,nqpt))
!  do isym=1,nsym_new
!  call mati3inv(symrel_new(:,:,isym),symrec_new(:,:,isym))
!  end do
!  call symkpt(0,gmet,indqpt,ab_out,spqpt,nqpt,nqpt2,nsym_new,&
!  &     symrec_new,use_inversion,wghtq,wtq_folded)

!  write(std_out,*)' after symkpt, nqpt2=',nqpt2

!  do iqpt=1,nqpt2
!  wghtq(iqpt)=wtq_folded(indqpt(iqpt))
!  qpt2(:,iqpt)=spqpt(:,indqpt(iqpt))
!  end do
!  nqpt=nqpt2
!  spqpt(:,1:nqpt)=qpt2(:,1:nqpt)
!  deallocate(qpt2,wtq_folded,symrec_new)

!  write(std_out,*)ch10,' after symkpt, nqpt=',nqpt
!  do iqpt=1,nqpt
!  write(std_out,*)iqpt,spqpt(:,iqpt),wghtq(iqpt)
!  end do

!  end if ! thmflag=7 or 8

   ABI_DEALLOCATE(ptsymrel)
   ABI_DEALLOCATE(symafm_new)
   ABI_DEALLOCATE(symrel_new)
   ABI_DEALLOCATE(tnons_new)
   
 end if

 write(message,'(a,a)')ch10,' thmeig : list of q wavevectors, with integration weights '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 do iqpt=1,nqpt
   write(message,'(i6,3es16.6,es20.6)')iqpt,spqpt(:,iqpt),wghtq(iqpt)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end do

!DEBUG
!write(std_out,*)'-thmeig : 3 '
!call flush(6)
!ENDDEBUG

 if(.not.allocated(indqpt))allocate(indqpt(nqpt))
 ABI_ALLOCATE(dedni,(mband,nkpt,3*natom,nqpt))
 ABI_ALLOCATE(dednr,(mband,nkpt,3*natom,nqpt))
 ABI_ALLOCATE(eigen_in,(nqpt))
 ABI_ALLOCATE(qpt_full,(3,nqpt))
 ABI_ALLOCATE(qptnrm,(nqpt))

 dednr(:,:,:,:) = zero
 dedni(:,:,:,:) = zero

!Prepare the reading of the EIG2 files
 choice=1
 nunit=ddbun
 
 call ioddb8_in(filnam5,natom,mband,&
& nkpt,msym,ntypat,nunit,vrsddb,&
& acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
& natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
& pawecutdg,rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
& typat,usepaw,wtk,xred,zion,znucl)

!Read the psp information of the input EIG2 file
 useylm=usepaw;choice=1
 call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
& nblok2,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

!iqpt2 will be the index of the q point bloks inside the EIG2 file
 iqpt2=0

!Sum on all phonon wavevectors and modes
 do iqpt=1,nqpt

!  Finding the target wavevector in DDB file
   qpt_search(:,:)=0.0d0
   qpt_search(:,1)=spqpt(:,iqpt)
   qphnrm(:)=one
   rfphon(1:2)=1
!  NOTE : at present, no LO-TO splitting included !!!
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   rftyp=1

   write(std_out,'(a,3es16.6)' )' Looking for spqpt=',qpt_search(:,1)

   call gtblk9(ddb_blk,iblok,mpert,natom,&
&   qpt_search,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

   if(iblok==0) then
     write(message,'(4a,3es16.6,2a)') ch10,&
&     ' thmeig : ERROR -',ch10,&
&     '   Was unable to find in DDB file, the blok for point ',spqpt(:,iqpt),ch10,&
&     '  Action : compute the contribution from this point, and merge it in your DDB file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   ABI_ALLOCATE(d2cart,(2,msize))
!  Copy the dynamical matrix in d2cart
   d2cart(:,1:msize)=ddb_blk%val(:,:,iblok)
   write (100,*) 'dynmat in thmeig'
   write (100,*) d2cart(:,1:msize)

!  Eventually impose the acoustic sum rule based on previously calculated d2asr
   if (anaddb_dtset%asr==1 .or. anaddb_dtset%asr==2 .or. anaddb_dtset%asr==5) then
     call asria_corr(anaddb_dtset%asr,d2asr,d2cart,mpert,natom)
   end if

!  DEBUG 
!  do ii=1,msize
!  write(std_out,*)' thmeig : d2cart(:,ii)=',d2cart(:,ii)
!  enddo
!  ENDDEBUG

!  Calculation of the eigenvectors and eigenvalues
!  of the dynamical matrix
   ABI_ALLOCATE(displ,(2*3*natom*3*natom))
   ABI_ALLOCATE(eigval,(3,natom))
   call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&   mpert,msym,natom,nsym,ntypat,phfreq(:,iqpt),qphnrm(1),spqpt(:,iqpt),rprimd,anaddb_dtset%symdynmat,&
&   symrel,typat,ucvol)
   ABI_DEALLOCATE(displ)
   ABI_DEALLOCATE(eigval)
   ABI_DEALLOCATE(d2cart)

!  DEBUG
!  write(std_out,*)"iqpt=",iqpt,"/",nqpt
!  ENDDEBUG

!  Read the next bloks to find the next q point.
   found=0 ; iqpt2_previous=iqpt2
   do while (iqpt2<nblok2)
     iqpt2=iqpt2+1
     call read_blok8(ddb_blk_eig2,iqpt2,mband,mpert_eig2,msize2,&
&     nkpt,nunit,blkval2(:,:,:,:),kpnt(:,:,1))
     write (300,*) 'blkval2 _bis_ in thmeig'
     write (300,*) blkval2
     diff_qpt(:)=ddb_blk_eig2%qpt(1:3,iqpt2)/ddb_blk_eig2%nrm(1,iqpt2)-spqpt(:,iqpt)
     if(diff_qpt(1)**2+diff_qpt(2)**2+diff_qpt(3)**2 < qtol )then
       found=1
       exit
     end if
   end do

!  Usually, the q points come in the right order. However, this is not always the case...
   if(found==0)then
     
!    If the EIG2 database file has to be read again, close it, then search for the right q point,
!    from the beginning of the file
     close(ddbun)

     call ioddb8_in(filnam5,natom,mband,&
&     nkpt,msym,ntypat,nunit,vrsddb,&
&     acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&     natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&     pawecutdg,rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
&     typat,usepaw,wtk,xred,zion,znucl)

!    Read the psp information of the input DDB
     useylm=usepaw;choice=1
     call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
&     nblok2,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

!    And examine again the EIG2 file. Still, not beyond the previously examined value.
     found=0 
     do iqpt2=1,iqpt2_previous
       call read_blok8(ddb_blk_eig2,iqpt2,mband,mpert_eig2,msize2,&
&       nkpt,nunit,blkval2(:,:,:,:),kpnt(:,:,1)) 
       diff_qpt(:)=ddb_blk_eig2%qpt(1:3,iqpt2)/ddb_blk_eig2%nrm(1,iqpt2)-spqpt(:,iqpt)
       if(diff_qpt(1)**2+diff_qpt(2)**2+diff_qpt(3)**2 < qtol )then
         found=1
         exit
       end if
     end do

     if(found==0)then
       write(message,'(4a,3es16.6,2a)') ch10,&
&       ' thmeig : ERROR -',ch10,&
&       '   Was unable to find in EIG2 DDB file, the blok for point ',spqpt(:,iqpt),ch10,&
&       '  Action : compute the contribution from this point, and merge it in your EIG2 DDB file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

   end if

!  Put blkval2 in cartesian coordinates
   call carttransf(ddb_blk_eig2%flg,blkval2,carflg_eig2,gprimd,iqpt,mband,&
&   mpert_eig2,msize2,natom,nblok2,nkpt,rprimd)

   do imod=1,3*natom

!    Calculate the derivative
     deigr(:,:) = zero
     deigi(:,:) = zero
     dwtermr(:,:)=zero
     dwtermi(:,:)=zero
     index=0
     do iatom1=1,natom
       do idir1=1,3
         do iatom2=1,natom
!          Compute factor for SE term
           if(phfreq(imod,iqpt)<tol6)then
             factr = zero
           else
             factr=one/sqrt(amu(typat(iatom1))*amu(typat(iatom2)))/phfreq(imod,iqpt)/amu_emass
           end if
           
           do idir2=1,3
             index = idir1 + 3*((iatom1 - 1) + natom * ((idir2-1)+3*(iatom2-1)))

!            Compute products of polarization vectors
             vecr = eigvec(1,idir1,iatom1,imod)*eigvec(1,idir2,iatom2,imod)+&
&             eigvec(2,idir1,iatom1,imod)*eigvec(2,idir2,iatom2,imod)
             veci = eigvec(2,idir1,iatom1,imod)*eigvec(1,idir2,iatom2,imod)-&
&             eigvec(1,idir1,iatom1,imod)*eigvec(2,idir2,iatom2,imod)

             vec1r = eigvec(1,idir1,iatom1,imod)*eigvec(1,idir2,iatom1,imod)+&
&             eigvec(2,idir1,iatom1,imod)*eigvec(2,idir2,iatom1,imod)
             vec1i = eigvec(2,idir1,iatom1,imod)*eigvec(1,idir2,iatom1,imod)-&
&             eigvec(1,idir1,iatom1,imod)*eigvec(2,idir2,iatom1,imod)

             vec2r = eigvec(1,idir1,iatom2,imod)*eigvec(1,idir2,iatom2,imod)+&
&             eigvec(2,idir1,iatom2,imod)*eigvec(2,idir2,iatom2,imod)
             vec2i = eigvec(2,idir1,iatom2,imod)*eigvec(1,idir2,iatom2,imod)-&
&             eigvec(1,idir1,iatom2,imod)*eigvec(2,idir2,iatom2,imod)
             
!            Compute factor for DW term
             if(phfreq(imod,iqpt)<tol6)then
               fact2r = zero
               fact2i = zero
             else
               fact2r = -wghtq(iqpt)*(vec1r/amu(typat(iatom1)) + vec2r/amu(typat(iatom2)))/phfreq(imod,iqpt)/&
&               amu_emass/2 !/norm(idir1)/norm(idir2)
               fact2i = -wghtq(iqpt)*(vec1i/amu(typat(iatom1)) + vec2i/amu(typat(iatom2)))/phfreq(imod,iqpt)/&
&               amu_emass/2 !/norm(idir1)/norm(idir2)
             end if

             multr(:,:) =(blkval2(1,index,:,:)*vecr - blkval2(2,index,:,:)*veci) !/(norm(idir1)*norm(idir2))
             multi(:,:) =(blkval2(1,index,:,:)*veci + blkval2(2,index,:,:)*vecr) !/(norm(idir1)*norm(idir2))

!            DEBUG
!            write(std_out,*) 'factr et facti',factr,facti
!            write(std_out,*) 'fact2r et fact2i',fact2r,fact2i
!            write(std_out,*) 'multr et multi', multr, multi
!            ENDDEBUG

!            Debye-Waller Term
             if(thmflag==3 .or. thmflag==5 .or. thmflag==7) then
               dwtermr(1:mband,1:nkpt)=dwtermr(1:mband,1:nkpt)+fact2r*blkval2gqpt(1,index,:,:)-fact2i*blkval2gqpt(2,index,:,:)
               dwtermi(1:mband,1:nkpt)=dwtermi(1:mband,1:nkpt)+fact2r*blkval2gqpt(2,index,:,:)+fact2i*blkval2gqpt(1,index,:,:)
             end if

!            Self-energy Term (Fan)
             deigr(1:mband,1:nkpt) = deigr(1:mband,1:nkpt) + wghtq(iqpt)*factr*multr(1:mband,1:nkpt)
             deigi(1:mband,1:nkpt) = deigi(1:mband,1:nkpt) + wghtq(iqpt)*factr*multi(1:mband,1:nkpt)

           end do !idir2
         end do !iatom2
       end do !idir1
     end do !iatom1
!    Eigenvalue derivative or broadening
     if(thmflag==3 .or. thmflag==5 .or. thmflag==7) then  
       dednr(1:mband,1:nkpt,imod,iqpt) = deigr(1:mband,1:nkpt) + dwtermr(1:mband,1:nkpt)
       dedni(1:mband,1:nkpt,imod,iqpt) = deigi(1:mband,1:nkpt) + dwtermi(1:mband,1:nkpt)
     else if(thmflag==4 .or. thmflag==6 .or. thmflag==8) then
       dednr(1:mband,1:nkpt,imod,iqpt) = pi*deigr(1:mband,1:nkpt) 
       dedni(1:mband,1:nkpt,imod,iqpt) = pi*deigi(1:mband,1:nkpt) 
     end if

   end do ! imod
 end do !iqpt

 close(ddbun)

!DEBUG
!write(std_out,*)'mband, nkpt, nqpt',mband,nkpt,nqpt
!write(std_out,*)'dednr ',dednr(5,50,1,50)
!write(std_out,*)'dednr ',dednr(4,54,1,55)
!write(std_out,*)'dedni ',dedni(5,50,1,50)
!ENDDEBUG


!=============================================================================
!3) Evaluation of the Eliashberg type spectral function 
!and phonon DOS via gaussian broadning
!=============================================================================

 if(telphint==1)then
   ng2f = 500  ! number of frequencies
   omega_min=zero
   omega_max=zero
   do iqpt=1,nqpt
     do imod=1,3*natom
       omega_min = min(omega_min,phfreq(imod,iqpt))
       omega_max = max(omega_max,phfreq(imod,iqpt))
     end do
   end do

   ABI_ALLOCATE(dos_phon,(ng2f))
   ABI_ALLOCATE(g2f,(mband,nkpt,ng2f))
   ABI_ALLOCATE(tmpg2f,(mband,nkpt,ng2f))
   ABI_ALLOCATE(tmpphondos,(ng2f))

   write(std_out,'(a,es13.6)') 'omega_min :', omega_min
   write(std_out,'(a,es13.6)') 'omega_max :', omega_max
   write(std_out,'(a,i8)') 'ng2f :', ng2f
   
   omega_max = omega_max + 0.1 * omega_max
   domega = (omega_max-omega_min)/(ng2f-one)   

   gaussprefactor = sqrt(piinv) / g2fsmear    
   gaussfactor = one / g2fsmear

   g2f(:,:,:) = zero
   dos_phon(:) = zero

   do iqpt=1,nqpt
     do imod=1,3*natom
       omega = omega_min     
       tmpg2f(:,:,:) = zero
       tmpphondos(:) = zero
       do iomega=1,ng2f
         xx = (omega-phfreq(imod,iqpt))*gaussfactor
         gaussval = gaussprefactor*exp(-xx*xx)
         tmpg2f(:,:,iomega) = tmpg2f(:,:,iomega) + gaussval*dednr(:,:,imod,iqpt)
         tmpphondos(iomega) = tmpphondos(iomega) + gaussval
         omega = omega+domega
       end do

       g2f(:,:,:) = g2f(:,:,:) + tmpg2f(:,:,:)
       dos_phon(:) = dos_phon(:) + tmpphondos(:)

     end do !imod
   end do !iqpt

   dos_phon(:) = dos_phon(:) / nqpt
   
!  output the g2f
   unit_g2f = 108
   call outg2f(domega,omega_min,omega_max,filnam,g2f,g2fsmear,kpnt,mband,ng2f,nkpt,nqpt,1,telphint,unit_g2f)

!  output the phonon DOS
   unit_phdos = 108
   call outphdos(domega,dos_phon,omega_min,omega_max,filnam,g2fsmear,ng2f,nqpt,1,telphint,unit_g2f)

   
   ABI_DEALLOCATE(dos_phon)
   ABI_DEALLOCATE(g2f)
   ABI_DEALLOCATE(tmpg2f)
   ABI_DEALLOCATE(tmpphondos)
   
 end if !telphint

!=======================================================================
!4) Evaluation of the Eliashberg type spectral function
!and phonon DOS via improved tetrahedron method 
!=======================================================================

 if(telphint==0)then

!  make dimension-ful rprimd and gprimd for transformation of derivatives to cartesian coordinates.
   call mkrdim(acell,rprim,rprimd)
   call matr3inv(rprimd,gprimd)

!  Q point Grid
   qpt_full(:,:) = ddb_blk%qpt(1:3,:)

!  Trivial Q point index  
   do iqpt=1,nqpt
     indqpt(iqpt)=iqpt
     qptnrm(iqpt)= qpt_full(1,iqpt)*qpt_full(1,iqpt)+qpt_full(2,iqpt)*qpt_full(2,iqpt)+qpt_full(3,iqpt)*qpt_full(3,iqpt)
   end do

!  Build qlatt from scratch (for 5.7)
   tol = 0.1_dp
   ilatt = 0
   call sort_dp(nqpt,qptnrm,indqpt,tol)

   do iqpt1=1,nqpt-2
     mesh(1:3,1) = qpt_full(1:3,indqpt(iqpt1))
     do iqpt2=iqpt1+1,nqpt-1
       mesh(1:3,2)= qpt_full(1:3,indqpt(iqpt2))
       do iqpt3=iqpt2+1,nqpt
         mesh(1:3,3)= qpt_full(1:3,indqpt(iqpt3))
         det = mesh(1,1)*mesh(2,2)*mesh(3,3) + mesh(1,2)*mesh(2,3)*mesh(3,1) + mesh(1,3)*mesh(2,1)*mesh(3,2) &
&         -mesh(3,1)*mesh(2,2)*mesh(1,3) - mesh(3,2)*mesh(2,3)*mesh(1,1) - mesh(3,3)*mesh(2,1)*mesh(1,2)
         invdet = one/det
         if (abs(nint(invdet))==nqpt .and. abs(invdet)-nqpt < tol) then
           ilatt = 1
           qlatt(:,:) = mesh(:,:)
           exit
         end if            
       end do
       if(ilatt==1) exit
     end do
     if(ilatt==1) exit
   end do

!  error message if qlatt not found and stop
   if(ilatt==0) then
     write(message, '(a,a)' ) &
&     ' Could not find homogeneous basis vectors for Q point grid ',ch10
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     stop
   end if

!  test if qlatt is righthanded and possibly fixe it
   if(invdet < 0) then
     tempqlatt(:) = qlatt(:,2)
     qlatt(:,2) = qlatt(:,1)
     qlatt(:,1) = tempqlatt(:)    
   end if

   write(std_out,*) 'qlatt',qlatt

!  test if qlatt generates all Q points  TO DO



!  Get tetrahedra, ie indexes of the full kpoints at their summits
   call init_tetra(indqpt,gprimd,qlatt,qpt_full,nqpt,&
&   tetrahedra, ierr, errstr)
   ABI_CHECK(ierr==0,errstr)
   
   rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&   -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&   +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

!  Calculate weights for phonon DOS
!  Special precautions must be taking for Gamma point
!  because of non-analytic term.
!  Non-analyticity must be taken out and treated separatly.

   nene = 100     !nene=number of energies for DOS
   enemin = minval(phfreq) 
   enemax = maxval(phfreq) 
   deltaene = (enemax-enemin)/dble(nene-1)
!  redefine enemin enemax to be at rounded multiples of deltaene
!  enemin = elph_ds%fermie - dble(ifermi)*deltaene
!  enemax = elph_ds%fermie + dble(nene-ifermi-1)*deltaene

   ABI_ALLOCATE(tweight,(nqpt,nene))
   ABI_ALLOCATE(dtweightde,(nqpt,nene))
   ABI_ALLOCATE(intweight,(3*natom,nqpt,nene))
   ABI_ALLOCATE(indtweightde,(3*natom,nqpt,nene))
   
   do iband=1,3*natom
     eigen_in(:) = phfreq(iband,:)

!    calculate general integration weights at each irred kpoint as in Blochl et al PRB 49 16223
     call get_tetra_weight(eigen_in,enemin,enemax,&
&     one,nene,nqpt,tetrahedra,&
&     tweight,dtweightde)

     intweight(iband,:,:) = tweight(:,:)
     indtweightde(iband,:,:) = dtweightde(:,:)
     
   end do !iband
   
!  intdtweightse(nband,nqpt,nene) represents the weight in each energy bin for every kpt and every band
!  So phonon DOS is calculated (neglecting the non-analyticity contribution for now !!!)

   ABI_ALLOCATE(total_dos,(nene))
   ABI_ALLOCATE(g2f,(mband,nkpt,nene))

   total_dos(:) = zero
   do iband=1,3*natom
     do iqpt=1,nqpt
       total_dos(:) = total_dos + indtweightde(iband,iqpt,:)
     end do
   end do

!  For the g2f function
!  Right now for one electronic band and one K point: dednr(1:mband,1:nkpt,imod,iqpt)
!  Once again must pay close attention to the Gamma point
   g2f(:,:,:) = zero
   do ii=1,mband
     do ikpt=1,nkpt
       do iband=1,3*natom
         do iqpt=1,nqpt
           g2f(ii,ikpt,:) = g2f(ii,ikpt,:) + dednr(ii,ikpt,iband,iqpt) * indtweightde(iband,iqpt,:)
         end do
       end do
     end do
   end do

!  output the g2f
   unit_g2f = 108
   call outg2f(deltaene,enemin,enemax,filnam,g2f,g2fsmear,kpnt,mband,nene,nkpt,nqpt,tetrahedra%ntetra,telphint,unit_g2f)

!  output the phonon DOS
   unit_phdos = 108
   call outphdos(deltaene,total_dos,enemin,enemax,filnam,g2fsmear,nene,nqpt,tetrahedra%ntetra,telphint,unit_g2f)

   ABI_DEALLOCATE(tweight)
   ABI_DEALLOCATE(dtweightde)
   ABI_DEALLOCATE(intweight)
   ABI_DEALLOCATE(indtweightde)
   ABI_DEALLOCATE(total_dos)
   ABI_DEALLOCATE(g2f)
 end if !telphint

!=======================================================================
!5) direct evaluation of thermal corrections
!=======================================================================
 
 slope(:,:,:) = zero
 zeropoint(:,:,:) = zero
!Loop on temperatures
 do itemper= 1, ntemper
   tmp=tempermin+temperinc*float(itemper-1)
   thmeigen(:,:,:) = zero

!  Sum on all phonon wavevectors and modes
   do iqpt=1,nqpt
     do imod=1,3*natom

!      Bose-Einstein distribution 
       if(phfreq(imod,iqpt)<tol6)then
         bosein = zero
       else
         bosein = one/(exp(phfreq(imod,iqpt)/(kb_HaK*tmp))-1) 
       end if

!      Calculate total
       thmeigen(1,1:mband,1:nkpt) = thmeigen(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt,imod,iqpt)*(bosein+half)
       thmeigen(2,1:mband,1:nkpt) = thmeigen(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt,imod,iqpt)*(bosein+half)

       if(itemper==1)then
!        Calculate slope of linear regime
         if(phfreq(imod,iqpt)<tol6)then
           slope(1,1:mband,1:nkpt) = slope(1,1:mband,1:nkpt) 
           slope(2,1:mband,1:nkpt) = slope(2,1:mband,1:nkpt) 
         else
           slope(1,1:mband,1:nkpt) = slope(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt,imod,iqpt)*(kb_HaK/phfreq(imod,iqpt))
           slope(2,1:mband,1:nkpt) = slope(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt,imod,iqpt)*(kb_HaK/phfreq(imod,iqpt))
         end if
!        Calculate zero-point renormalization
         zeropoint(1,1:mband,1:nkpt) = zeropoint(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt,imod,iqpt)*half
         zeropoint(2,1:mband,1:nkpt) = zeropoint(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt,imod,iqpt)*half

!        DEBUG
!        write(std_out,*)' For iqpt,imod=',iqpt,imod
!        write(std_out,'(a,8f12.5)' )'  contribution to ZPM correction of ikpt=20, 1:mband',dednr(1:mband,20,imod,iqpt)*half*Ha_eV
!        ENDDEBUG
       end if
     end do ! imod
   end do !iqpt

!  Output
!  unitout should be attributed in dtset to avoid conflicts
   unitout = 42
   outfile = trim(filnam)//"_TBS"
   
!  open TBS file
   open (unit=unitout,file=outfile,form='formatted',status='unknown')
   write(unitout,'(a)')'thmeig: Thermal Eigenvalue corrections (eV)'
!  Write temperature independent results
   if(itemper==1)then
     write(unitout,'(a)')'Temperature independent results (zero-point renormalization and slope)'
     do ikpt=1,nkpt
       write(unitout,'(a,3es16.8)')' Kpt :', kpnt(:,ikpt,1)
       do iband=1,mband
         write(unitout,'(4d22.14)') Ha_eV*zeropoint(1,iband,ikpt),Ha_eV*zeropoint(2,iband,ikpt),&
&         Ha_eV*slope(1,iband,ikpt),Ha_eV*slope(2,iband,ikpt)
       end do
     end do
     write(unitout,'(a)')'Temperature dependent corrections'
   end if
!  Write result in file for each temperature
   write(unitout,'(a,es10.3,a)')'T :', tmp,' K'
   do ikpt=1,nkpt
     write(unitout,'(a,3es16.8)')' Kpt :', kpnt(:,ikpt,1)
     do iband=1,mband
       write(unitout,'(2d22.14)') Ha_eV*thmeigen(1,iband,ikpt), Ha_eV*thmeigen(2,iband,ikpt)
     end do
   end do
 end do !itemper
 
 close(unitout)

!Write temperature-independent results to the main output file 
 write(iout,'(a)')' '
 write(iout,'(80a)') ('-',ii=1,80)
 write(iout,'(a)')' '
 write(iout,'(a)')' Electron-phonon change of electronic structure.'
 write(iout,'(a)')' The temperature-dependent values are written in the _TBS file.'
 write(iout,'(a)')' Here follows, for each electronic wavevector and band :'
 write(iout,'(a)')'      zero-point renormalisation (Ha) and linear slope (Ha/Kelvin)'
 do ikpt=1,nkpt
   write(iout,'(2a,i6,a,3es16.6)')ch10,' Kpt number ',ikpt,', with reduced coordinates :',kpnt(:,ikpt,1)
   do iband=1,mband
     write(iout,'(i6,2es20.6)') iband,zeropoint(1,iband,ikpt),slope(1,iband,ikpt)
   end do
 end do

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

 ABI_DEALLOCATE(dedni)
 ABI_DEALLOCATE(dednr)
 if(allocated(indqpt))ABI_DEALLOCATE(indqpt)
 ABI_DEALLOCATE(eigen_in)
 ABI_DEALLOCATE(qpt_full)
 ABI_DEALLOCATE(qptnrm)
 ABI_DEALLOCATE(wghtq)
 ABI_DEALLOCATE(spqpt)
 ABI_DEALLOCATE(eigvec)
 ABI_DEALLOCATE(phfreq)

 ABI_DEALLOCATE(blkval2gqpt)
 ABI_DEALLOCATE(kpnt)
 ABI_DEALLOCATE(carflg_eig2)

 call destroy_ddb_blk(ddb_blk_eig2)
 ABI_DEALLOCATE(ddb_blk_eig2)

 call destroy_tetra(tetrahedra)

!DEBUG
!close(ddbun)
!ENDEBUG

end subroutine thmeig
!!***

