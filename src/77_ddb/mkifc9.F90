!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkifc9
!!
!! NAME
!! mkifc9
!!
!! FUNCTION
!! This routine makes the interatomic force constants,
!! taking into account the dipole-dipole interaction,
!! and also perform some checks.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of cell (bohr)
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! ddb_blk = storage object for ddb information read in from DDB file
!! dielt(3,3)=dielectric tensor
!! gmet(3,3)=reciprocal space metric tensor in bohr**-2
!! gprim(3,3)=dimensionless reciprocal space primitive translations
!! indsym(4,msym*natom)=indirect indexing array for symmetries
!! iout=unit number for output of formatted data
!! mpert =maximum number of ipert
!! msym =maximum number of symmetries
!! natom=number of atoms in cell
!! ngqpt_in = input values of ngqpt
!! nsym=number of space group symmetries
!! ntypat=number of atom types
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! rprim(3,3)= primitive translation vectors
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! tcpui,twalli=initial values of cpu and wall clocktime
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3
!! xred(3,natom)=fractional dimensionless atomic coordinates
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!! dyewq0(3,3,natom)=atomic self-interaction correction to the
!!  dynamical matrix. (only when dipdip=1)
!! rcan(3,natom) = Atomic position in canonical coordinates
!! trans(3,natom) = Atomic translations : xred = rcan + trans
!!
!! ifc_obj%nrpt = number of vectors of the lattice in real space
!! ifc_obj%atmfrc(2,3,natom,3,natom,nrpt)
!!  = Interatomic Forces in real space
!! ifc_obj%rpt(3,nprt) = Canonical coordinates of the R points in the unit cell
!! ifc_obj%wghatm(natom,natom,nrpt)= Weight associated to the couple of atoms and the R vector
!!
!! NOTES
!! 1. Should be executed by one processor only.
!! 2. Work variables have been localized here.
!!
!! PARENTS
!!      anaddb,refineblk
!!
!! CHILDREN
!!      asrif9,bigbx9,canat9,chkrp9,destroy_ifc,dymfz9,ewald9,ftifc_q2r,gtdyn9
!!      hybrid9,mkrdim,nanal9,nullify_ifc,omega_decomp,phfrq3,prtph3
!!      q0dy3_apply,q0dy3_calc,rsiaf9,smpbz,symdm9,timein,wght9,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkifc9(acell,amu,anaddb_dtset,&
& ddb_blk,dielt,dyewq0,gmet,gprim,&
& ifc_obj,indsym,iout,mpert,msym,natom,ngqpt_in,&
& nsym,ntypat,rcan,rmet,rprim,&
& symrec,symrel,tcpui,trans,twalli,typat,&
& ucvol,xred,zeff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_ifc
 use m_ddb_blk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkifc9'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_56_recipspace
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => mkifc9
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,msym,natom,nsym,ntypat
 real(dp),intent(in) :: tcpui,twalli,ucvol
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
 type(ifc_type), intent(out) :: ifc_obj
!arrays
 type(ddb_blk_type), intent(in) :: ddb_blk
 integer,intent(in) :: indsym(4,nsym*natom)
 integer,intent(in) :: symrec(3,3,msym),symrel(3,3,msym),typat(natom)
 integer,intent(in) :: ngqpt_in(3)
 real(dp),intent(in) :: acell(3),amu(ntypat)
 real(dp),intent(in) :: gmet(3,3),gprim(3,3),rmet(3,3),rprim(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(inout) :: dielt(3,3),zeff(3,3,natom)
 real(dp),intent(out) :: dyewq0(3,3,natom)
 real(dp),intent(out) :: rcan(3,natom),trans(3,natom)

!Local variables -------------------------
!scalars
 integer :: asr,brav,choice,dipdip,iqpt,mqpt,mrpt,srlr
 integer :: nqpt,nqshft,option,plus,sumg0
 integer :: irpt, irpt_new
 real(dp) :: qphnrm,tcpu,twall
 character(len=500) :: message
!arrays
 integer :: ngqpt(9),qptrlatt(3,3)
 integer,allocatable :: atifc(:)
 real(dp) :: qpt(3),rprimd(3,3), dummy_rpt(3,1)
 real(dp),allocatable :: dyew(:,:,:,:,:),dynmat(:,:,:,:,:,:),spqpt(:,:)

! tmp variables for phonon interpolation
 real(dp) :: d2cart(2,3,mpert,3,mpert),displ(2*3*natom*3*natom)
 real(dp) :: eigval(3,natom),eigvec(2,3,natom,3,natom),phfrq(3*natom)

 type(ifc_type) :: ifc_tmp
 real(dp),allocatable :: dynmatfull(:,:,:,:,:,:),dynmatsr(:,:,:,:,:,:),dynmatlr(:,:,:,:,:,:) ! for OmegaSRLR

!******************************************************************

 asr=anaddb_dtset%asr
 brav=anaddb_dtset%brav
 dipdip=anaddb_dtset%dipdip
 ngqpt = 0
 ngqpt(1:3)=ngqpt_in(1:3)
 nqshft=anaddb_dtset%nqshft
 srlr=anaddb_dtset%prtsrlr

 ABI_ALLOCATE(atifc,(natom))
 atifc(:)=anaddb_dtset%atifc(:)

 call mkrdim(acell,rprim,rprimd)

 if(dipdip==1)then

   if(asr==1.or.asr==2)then
!    Calculation of the non-analytical part for q=0
     sumg0=0
     qpt(1)=0.0_dp
     qpt(2)=0.0_dp
     qpt(3)=0.0_dp
     ABI_ALLOCATE(dyew,(2,3,natom,3,natom))
     call ewald9(acell,dielt,dyew,gmet,gprim,natom,&
&     qpt,rmet,rprim,sumg0,ucvol,xred,zeff)
     option=asr
     call q0dy3_calc(natom,dyewq0,dyew,option)
     ABI_DEALLOCATE(dyew)
   else if (asr==0) then
     dyewq0(:,:,:)=0.0_dp
   end if
 end if

!Check if the rprim are coherent with the choice used in
!the interatomic forces generation
 call chkrp9(brav,rprim)

!Sample the Brillouin zone
 option=1
 qptrlatt(:,:)=0
 qptrlatt(1,1)=ngqpt(1)
 qptrlatt(2,2)=ngqpt(2)
 qptrlatt(3,3)=ngqpt(3)
 mqpt=ngqpt(1)*ngqpt(2)*ngqpt(3)*nqshft
 if(brav==2)mqpt=mqpt/2
 if(brav==3)mqpt=mqpt/4
 ABI_ALLOCATE(spqpt,(3,mqpt))
 call smpbz(brav,iout,qptrlatt,mqpt,nqpt,nqshft,option,anaddb_dtset%q1shft,spqpt)

 ABI_ALLOCATE(dynmat,(2,3,natom,3,natom,nqpt))

!Find symmetrical dynamical matrices
 call symdm9(ddb_blk%flg,ddb_blk%nrm,ddb_blk%qpt,ddb_blk%typ,ddb_blk%val,&
& dynmat,gprim,indsym,mpert,natom,ddb_blk%nblok,nqpt,nsym,anaddb_dtset%rfmeth,rprim,spqpt,&
& symrec,symrel)

!OmegaSRLR: Store full dynamical matrix for decomposition into short- and long-range parts
 ABI_ALLOCATE(dynmatfull,(2,3,natom,3,natom,nqpt))
 dynmatfull=dynmat

 if(dipdip==1)then
!  Take off the dipole-dipole part of the dynamical matrix
   write(std_out, '(a)' )' mkifc9 : will extract the dipole-dipole part,'
   write(std_out, '(a)' )' using ewald9, q0dy3 and nanal9 for every wavevector.'
   do iqpt=1,nqpt
     write(std_out, '(a,i4)' )'    wavevector number :',iqpt
     qpt(1)=spqpt(1,iqpt)
     qpt(2)=spqpt(2,iqpt)
     qpt(3)=spqpt(3,iqpt)
     sumg0=0

     call timein(tcpu,twall)
     write(message, '(a,f11.3,a,f11.3,a)' )&
&     '-ewald9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
     call wrtout(std_out,message,'COLL')
     ABI_ALLOCATE(dyew,(2,3,natom,3,natom))
     call ewald9(acell,dielt,dyew,gmet,gprim,natom,&
&     qpt,rmet,rprim,sumg0,ucvol,xred,zeff)

     call timein(tcpu,twall)
     write(message, '(a,f11.3,a,f11.3,a)' )&
&     '-q0dy3 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
     call wrtout(std_out,message,'COLL')
     call q0dy3_apply(natom,dyewq0,dyew)
     plus=0

     call timein(tcpu,twall)
     write(message, '(a,f11.3,a,f11.3,a)' )&
&     '-nanal9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
     call wrtout(std_out,message,'COLL')
     call nanal9(dyew,dynmat,iqpt,natom,nqpt,plus)
     ABI_DEALLOCATE(dyew)
   end do
 end if

!OmegaSRLR: Store the short-range dynmat and compute long-range as difference
 ABI_ALLOCATE(dynmatsr,(2,3,natom,3,natom,nqpt))
 ABI_ALLOCATE(dynmatlr,(2,3,natom,3,natom,nqpt))
 dynmatsr=dynmat
 dynmatlr=dynmatfull-dynmatsr


!Now, take care of the remaining part of the ynamical matrix

!Passage to canonical normalized coordinates
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-canat9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 call canat9(brav,natom,rcan,rprim,trans,xred)

!Multiply the dynamical matrix by a phase shift
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-dymfz9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 option=1
 call dymfz9(dynmat,natom,nqpt,gprim,option,spqpt,trans)

!Create the Big Box of R vectors in real space
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-bigbx9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')


!Compute the number of points (cells) in real space
 choice=0
 call bigbx9(brav,choice,1,ngqpt,nqshft,mrpt,rprim,dummy_rpt)

 call nullify_ifc(ifc_tmp)
 ifc_tmp%nrpt = mrpt
 ABI_ALLOCATE(ifc_tmp%rpt,(3,ifc_tmp%nrpt))
 ABI_ALLOCATE(ifc_tmp%wghatm,(natom,natom,ifc_tmp%nrpt))
 choice=1
 call bigbx9(brav,choice,mrpt,ngqpt,nqshft,ifc_tmp%nrpt,rprim,ifc_tmp%rpt)

!Weights associated to these R points and to atomic pairs
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-wght9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

 call wght9(brav,gprim,natom,ngqpt,nqpt,nqshft,ifc_tmp%nrpt,anaddb_dtset%q1shft,rcan,ifc_tmp%rpt,ifc_tmp%wghatm)

!Fourier transformation of the dynamical matrices
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-ftiaf9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

 ABI_ALLOCATE(ifc_tmp%atmfrc,(2,3,natom,3,natom,ifc_tmp%nrpt))

 call ftifc_q2r(ifc_tmp%atmfrc,dynmat,gprim,natom,nqpt,ifc_tmp%nrpt,ifc_tmp%rpt,spqpt)

 ABI_DEALLOCATE(dynmat)

!Eventually impose Acoustic Sum Rule
!to the interatomic forces
 if(asr>0)then
   call asrif9(asr,ifc_tmp%atmfrc,natom,ifc_tmp%nrpt,ifc_tmp%rpt,ifc_tmp%wghatm)
 end if

!*** The interatomic forces have been calculated ! ***
 write(message, '(2a)')ch10,&
& ' The interatomic forces have been obtained '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-hybrid9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

!Eventually artificially modify the IFC
 call hybrid9(acell,asr,ifc_tmp%atmfrc,dielt,dipdip,dyew,dyewq0,&
& gmet,gprim,iout,natom,ifc_tmp%nrpt,rcan,rmet,&
& rprim,ifc_tmp%rpt,ucvol,ifc_tmp%wghatm,xred,zeff)


!Analysis of the real-space interatomic force constants
 if (anaddb_dtset%nsphere/=0 .or. anaddb_dtset%rifcsph>tol10 .or. anaddb_dtset%ifcout/=0) then

   write(std_out,'(/,a)' )' mkifc9 : analysis of real-space IFCs'

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-rsiaf9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')

   call rsiaf9(acell,atifc,ifc_tmp%atmfrc,dielt,dipdip,dyewq0,&
&   gprim,anaddb_dtset%ifcana,anaddb_dtset%ifcout,iout,&
&   natom,ifc_tmp%nrpt,anaddb_dtset%nsphere,rcan,anaddb_dtset%rifcsph,rprim,ifc_tmp%rpt,&
&   tcpui,twalli,ifc_tmp%wghatm,zeff)
!  write(std_out,'(a)' )' mkifc9 : real space analysis of the IFC''s done '
 end if

 ABI_DEALLOCATE(atifc)

!Eventually impose Acoustic Sum Rule
!to the interatomic forces
!(Note : here, after the analysis, in which the range
!of the short-range IFCs may have been changed
!That is why it is asked that asr be negative)
!FIXME: asr < 0 is not tested in abinit suite and does not appear to work
!(I get frequencies of 10^105 Hartree...) Modifying this 12/6/2011
!if(asr<0)then
!asr=-asr

 if (asr > 0 .and. (anaddb_dtset%nsphere/=0 .or. anaddb_dtset%rifcsph>tol10)) then
   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-asrif9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call asrif9(asr,ifc_tmp%atmfrc,natom,ifc_tmp%nrpt,ifc_tmp%rpt,ifc_tmp%wghatm)
 end if

!Only conserve the necessary points in rpt: in the FT algorithm
!the order of the points is unimportant
 ifc_obj%nrpt = 0
 do irpt=1,ifc_tmp%nrpt
   if (sum(ifc_tmp%wghatm(:,:,irpt)) /= 0) then
     ifc_obj%nrpt = ifc_obj%nrpt+1
   end if
 end do

 ABI_ALLOCATE(ifc_obj%atmfrc,(2,3,natom,3,natom,ifc_obj%nrpt))
 ABI_ALLOCATE(ifc_obj%rpt,(3,ifc_obj%nrpt))
 ABI_ALLOCATE(ifc_obj%wghatm,(natom,natom,ifc_obj%nrpt))

 irpt_new = 1
 do irpt = 1, ifc_tmp%nrpt
   if (sum(ifc_tmp%wghatm(:,:,irpt)) /= 0) then
     ifc_obj%atmfrc(:,:,:,:,:,irpt_new) = ifc_tmp%atmfrc(:,:,:,:,:,irpt)
     ifc_obj%rpt(:,irpt_new) = ifc_tmp%rpt(:,irpt)
     ifc_obj%wghatm(:,:,irpt_new) = ifc_tmp%wghatm(:,:,irpt)
     irpt_new = irpt_new + 1
   end if
 end do

 call destroy_ifc(ifc_tmp)


!Check that the starting values are well reproduced.
!(This is to be suppressed in a future version)
 write(std_out, '(a,a)' )' mkifc9 : now check that the starting values ',&
& ' are reproduced after the use of interatomic forces '

!call timein(tcpu,twall)
!write(std_out,1000)tcpu-tcpui,twall-twalli
 do iqpt=1,nqpt
   qpt(1)=spqpt(1,iqpt)
   qpt(2)=spqpt(2,iqpt)
   qpt(3)=spqpt(3,iqpt)
   qphnrm=1.0_dp

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-gtdyn9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')

!  The dynamical matrix d2cart is calculated here :
   call gtdyn9(acell,ifc_obj%atmfrc,dielt,dipdip,&
&   dyewq0,d2cart,gmet,gprim,mpert,natom,&
&   ifc_obj%nrpt,qphnrm,qpt,rmet,rprim,ifc_obj%rpt,&
&   trans,ucvol,ifc_obj%wghatm,xred,zeff)

!  Calculation of the eigenvectors and eigenvalues
!  of the dynamical matrix
   call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&   mpert,msym,natom,nsym,ntypat,phfrq,qphnrm,qpt,rprimd,&
&   anaddb_dtset%symdynmat,symrel,typat,ucvol)

!  OmegaSRLR: Perform decomposition of dynamical matrix
   if(srlr==1)then
     call omega_decomp(amu,natom,ntypat,typat,&
&     dynmatfull,dynmatsr,dynmatlr,iqpt,nqpt,eigvec)
   end if

!  Write the phonon frequencies (this is for checking purposes).
!  Note : these phonon frequencies are not written on unit
!  iout, only on unit 6.
   call prtph3(displ,0,anaddb_dtset%enunit,-1,-1,natom,phfrq,qphnrm,qpt)
 end do

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-end of mkifc9 at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(spqpt)
!OmegaSRLR: deallocate memory used by dynmat decomposition
 ABI_DEALLOCATE(dynmatfull)
 ABI_DEALLOCATE(dynmatsr)
 ABI_DEALLOCATE(dynmatlr)


end subroutine mkifc9
!!***
