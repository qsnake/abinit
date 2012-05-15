!!****f* ABINIT/mkphbs
!!
!! NAME
!! mkphbs
!!
!! FUNCTION
!! Function to calculate the phonon band structure, from the IFC
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (XG,MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3) =length scales by which rprim is to be multiplied
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! atmfrc(2,3,natom,3,natom,nrpt) = Interatomic Forces in real space
!! dielt(3,3)=dielectric tensor
!! dyewq0(3,3,natom)=Ewald part of the dynamical matrix, at q=0.
!! gmet(3,3)= metric tensor in reciprocal space.
!! gprim(3,3)=normalized coordinates in reciprocal space
!! indsym(4,nsym,natom)=label given by subroutine symatm, indicating atom
!!  label which gets rotated into given atom by given symmetry
!!  (first three elements are related primitive translation--
!!  see symatm where this is computed)
!! mpert =maximum number of ipert
!! msym = maximum number of symmetries
!! natom=number of atoms in the unit cell
!! nrpt=number of R points in the Big Box
!! nsym=number of symmetries
!! ntypat=number of atom types
!! rmet(3,3)=metric tensor in real space.
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nprt)=canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!! symrel(3,3,nsym)=symmetry operations
!! tcpui=initial cpu time
!! trans(3,natom)=atomic translations : xred = rcan + trans
!! twalli=initial wall clock time
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume
!! wghatm(natom,natom,nrpt)=weights associated to a pair of atoms and to a R vector
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!!
!! NOTES
!! 
!! PARENTS
!!      anaddb,scphon_interpolate_phonon_and_dos
!!
!! CHILDREN
!!      asria_corr,asrprs,atprj_destroy,atprj_init,atprj_nullify,atprj_print
!!      beginprtscphon,end_sortph,endprtscphon,freeze_displ_allmodes,gtblk9
!!      gtdyn9,make_path,mkrdim,outlwf9,phfrq3,prtph3,prtscphon,prtvsound
!!      sortph,symph3,timein,wrap2_pmhalf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkphbs(acell,amu,anaddb_dtset,atmfrc,ddb_blk,&
& d2asr,dielt,dyewq0,outfile_radix,gmet,gprim,indsym,iodyn,&
& mpert,msize,msym,natom,nrpt,nsym,ntypat,&
& qtol,rmet,rprim,rpt,singular,symrel,tcpui,  &
& trans,twalli,typat,ucvol,uinvers,vtinvers,wghatm,xred,zeff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_io_tools
 use m_bz_mesh
 use m_prtscphon
 use m_atprj
 use m_sortph
 use m_ddb_blk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkphbs'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => mkphbs
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars

 type(anaddb_dataset_type),intent(in) :: anaddb_dtset

 integer,intent(in) :: mpert,msym,natom,nrpt,nsym,ntypat
 integer,intent(in) :: iodyn,msize

 real(dp),intent(in) :: tcpui,twalli,ucvol, qtol
 type(ddb_blk_type), pointer :: ddb_blk

!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym)
 integer,intent(in) :: typat(natom)

 character(len=fnlen),intent(in) :: outfile_radix

 real(dp),intent(in) :: acell(3),amu(ntypat),dielt(3,3),gmet(3,3),gprim(3,3)
 real(dp),intent(in) :: rmet(3,3),rprim(3,3),rpt(3,nrpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt),xred(3,natom)
 real(dp),intent(in) :: zeff(3,3,natom)
 real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt),dyewq0(3,3,natom)

 real(dp),intent(inout) :: singular(1:3*natom*(3*natom-1)/2)
 real(dp),intent(inout) :: uinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)
 real(dp),intent(inout) :: vtinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)

 real(dp),intent(inout) :: d2asr(2,3,natom,3,natom)

!Local variables -------------------------
 integer :: iphl1, iblok
 integer :: rftyp, ii
 integer,parameter :: udispl=19,ufreq=18

 real(dp) :: tcpu, twall
 real(dp) :: res

 integer :: rfphon(4)
 integer :: rfelfd(4)
 integer :: rfstrs(4)
 
 real(dp) :: qphnrm(3), qphon(3), qphon_padded(3,3)
 real(dp) :: d2cart(2,msize)
 real(dp) :: real_qphon(3)

 real(dp) :: rprimd(3,3)
 real(dp) :: xcart(3,natom)

 real(dp) :: displ(2*3*natom*3*natom)
 real(dp) :: eigval(3,natom)

 real(dp), allocatable :: phfrq(:)
 real(dp), allocatable :: eigvec(:,:,:,:,:)

 character(len=fnlen) :: phonon_freq_filename, phonon_vec_filename
 character(len=fnlen) :: tmpfilename
 character(500) :: message

 integer :: nfineqpath
 integer, allocatable :: ndiv(:)
 real(dp), pointer :: fineqpath(:,:)

 type(atprj_type) :: t_atprj

! *********************************************************************

!DEBUG
!write(std_out,*)' mkphbs : enter '
!call flush(6) 
!ENDDEBUG

!setup xcart
 call mkrdim(acell,rprim,rprimd)
 xcart = matmul(rprimd, xred)


 nullify(fineqpath)
 nfineqpath = anaddb_dtset%nph1l
 fineqpath => anaddb_dtset%qph1l

 if(anaddb_dtset%nph1l==0) then
   if (anaddb_dtset%nqpath==0) then
!    if there is nothing to do, return
     return
   else
!    allow override of nph1l with nqpath if the former is not set
     ABI_ALLOCATE(ndiv,(anaddb_dtset%nqpath-1))
     call make_path(anaddb_dtset%nqpath,anaddb_dtset%qpath,gmet,'G',20,ndiv,nfineqpath,fineqpath)
   end if
 end if


 write(message, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
& ch10,' Treat the first list of vectors ',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

!Write some information in the lwf-formatted file
 if (anaddb_dtset%eivec==3) then
   call outlwf9(acell,iodyn,msym,natom,nfineqpath,nsym,ntypat,rprim,symrel,typat,xred)
 end if

 phonon_freq_filename = trim(outfile_radix)//"_PHFRQ"
 phonon_vec_filename = trim(outfile_radix)//"_PHVEC"
 call beginprtscphon(phonon_freq_filename, phonon_vec_filename, anaddb_dtset%outscphon)

 if (anaddb_dtset%natprj_bs > 0) then
   call atprj_nullify(t_atprj)
   call atprj_init(t_atprj, natom, anaddb_dtset%natprj_bs, anaddb_dtset%iatprj_bs, outfile_radix)
 end if

 ABI_ALLOCATE(phfrq,(3*natom))
 ABI_ALLOCATE(eigvec,(2,3,natom,3,natom))

 qphnrm = one

 do iphl1=1, nfineqpath

!  Initialisation of the phonon wavevector
   qphon(:)=fineqpath(:,iphl1)
   if (anaddb_dtset%nph1l /= 0) then
     qphnrm(1) = anaddb_dtset%qnrml1(iphl1)
   end if

!  Generation of the dynamical matrix in cartesian coordinates
   if(anaddb_dtset%ifcflag==1)then

!    Get d2cart using the interatomic forces and the
!    long-range coulomb interaction through Ewald summation
     write(message, '(a)' )' mkphbs    : enter gtdyn9 '
     call wrtout(std_out,message,'COLL')

     call timein(tcpu,twall)
     write(message, '(a,f11.3,a,f11.3,a)' )&
&     '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
     call wrtout(std_out,message,'COLL')
     call gtdyn9(acell,atmfrc,dielt,anaddb_dtset%dipdip,&
&     dyewq0,d2cart,gmet,gprim,mpert,natom,&
&     nrpt,qphnrm(1),qphon,rmet,rprim,rpt,&
&     trans,ucvol,wghatm,xred,zeff)

   else if(anaddb_dtset%ifcflag==0)then

!    Look for the information in the DDB (no interpolation here!)
     rfphon(1:2)=1
     rfelfd(1:2)=0
     rfstrs(1:2)=0
     rftyp=anaddb_dtset%rfmeth
     qphon_padded = zero
     qphon_padded(:,1) = qphon
     call gtblk9(ddb_blk,iblok,mpert,natom,&
&     qphon_padded,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

!    Copy the dynamical matrix in d2cart
     d2cart(:,1:msize)=ddb_blk%val(:,:,iblok)

!    Eventually impose the acoustic sum rule based on previously calculated d2asr
     if (anaddb_dtset%asr==1 .or. anaddb_dtset%asr==2 .or. anaddb_dtset%asr==5) then
       call asria_corr(anaddb_dtset%asr,d2asr,d2cart,mpert,natom)
     end if

!    Impose acoustic sum rule plus rotational symmetry for 0D and 1D systems
     if (anaddb_dtset%asr==3 .or. anaddb_dtset%asr==4) then
       call asrprs(anaddb_dtset%asr,2,3,uinvers,vtinvers,singular,d2cart,mpert,natom,xcart)
     end if
   end if

!  DEBUG 
!  do ii=1,msize
!  write(std_out,*)' mkphbs : d2cart(:,ii)=',d2cart(:,ii)
!  enddo
!  ENDDEBUG

!  Calculation of the eigenvectors and eigenvalues
!  of the dynamical matrix
   write(message, '(a)' )' mkphbs    : enter phfrq3 '
   call wrtout(std_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   write(std_out,*) ' in mkphbs with start phfrq3'
   call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&   mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,&
&   rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol)

   if (anaddb_dtset%freeze_displ > tol10) then
     real_qphon = zero
     if (abs(qphnrm(1)) > tol8) then
       real_qphon = qphon / qphnrm(1)
     end if
     call freeze_displ_allmodes(displ, anaddb_dtset%freeze_displ, natom, outfile_radix, phfrq, &
&     real_qphon, rprimd, typat, xcart)
   end if

!  if requested, output projection of each mode on given atoms
   if (anaddb_dtset%natprj_bs > 0) then
     call atprj_print(t_atprj, iphl1, phfrq, eigvec)
   end if

!  DEBUG
!  write(std_out,*)' mkphbs : before sortph '
!  call flush(6) 
!  ENDDEBUG

!  In case eivec == 4, write output files for band2eps
!  (visualization of phonon band structures)
   if (anaddb_dtset%eivec == 4) then
     tmpfilename = trim(outfile_radix)//"_B2EPS"
     call sortph(eigvec,displ,tmpfilename,natom,phfrq,udispl,ufreq)
   end if

!  DEBUG
!  write(std_out,*)' mkphbs : after sortph '
!  call flush(6)
!  ENDDEBUG

!  Write the phonon frequencies
   write(message, '(a)' )' mkphbs    : enter prtph3 '
   call wrtout(std_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call prtph3(displ,anaddb_dtset%eivec,anaddb_dtset%enunit,iodyn,ab_out,natom,phfrq,qphnrm(1),qphon)

!  write to file - present version of prtph3 is not atomic enough to do this and
!  has lots of other junk etc...
   call prtscphon(eigvec, natom, phfrq, qphon)

!  Determine the symmetries of the phonon mode at Gamma
!  TODO: generalize for other q-point little groups.
   if(sum(abs(qphon(:)))<qtol)then
     call symph3(ab_out,acell,eigvec,indsym,natom,nsym,phfrq,rprim,symrel)
   end if

!  if we have an acoustic mode (small q and acoustic type displacements)
!  extrapolate speed of sound in this direction, and Debye frequency
   call wrap2_pmhalf(qphon(1),real_qphon(1),res)
   call wrap2_pmhalf(qphon(2),real_qphon(2),res)
   call wrap2_pmhalf(qphon(3),real_qphon(3),res)
   if (sqrt(real_qphon(1)**2+real_qphon(2)**2+real_qphon(3)**2) < quarter .and. &
&   sqrt(real_qphon(1)**2+real_qphon(2)**2+real_qphon(3)**2) > tol6) then
     call prtvsound(eigvec, gmet, natom, phfrq, real_qphon, ucvol)
   end if

 end do ! iphl1

!deallocate sortph array
 call end_sortph()

!close unit for phonon frequencies
 call endprtscphon()

 if (anaddb_dtset%natprj_bs > 0) then
   call atprj_destroy(t_atprj)
 end if

 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(eigvec)

!DEBUG
!write(std_out,*)' mkphbs : exit '
!call flush(6)
!ENDDEBUG

end subroutine mkphbs
!!***
