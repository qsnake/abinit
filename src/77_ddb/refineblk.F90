!{\src2tex{textfont=tt}}
!!****f* ABINIT/refineblk
!!
!! NAME
!! refineblk
!!
!! FUNCTION
!!  Get a first set of interatomic force constants using a coarse q-point grid,
!!  interpolate onto required ngqpt grid,
!!  and reallocate blkval blocks and arrays to be used in the rest of anaddb with full grid.
!!      Should implement Gaal-Nagy's algorithm in PRB <b>73</b> 014117.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of cell (bohr)
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! dielt(3,3)=dielectric tensor
!! gmet(3,3)=reciprocal space metric tensor in bohr**-2
!! gprim(3,3)=dimensionless reciprocal space primitive translations
!! indsym(4,msym*natom)=indirect indexing array for symmetries
!! iout=unit number for output of formatted data
!! mpert =maximum number of ipert
!! msym =maximum number of symmetries
!! natom=number of atoms in cell
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
!! ddb_blk = datastructure with DDB information: type, qpt, 2DTE...
!!   this info is updated with interpolation + additional info in DDB
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      create_ddb_blk,destroy_ddb_blk,destroy_ifc,getkgrid,gtblk9,gtdyn9
!!      mkifc9,mkrdim
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine refineblk(acell,amu,anaddb_dtset,ddb_blk,&
&     dielt,gmet,gprim,indsym,iout,&
&     mpert,msym,natom,nsym,ntypat,rmet,rprim,&
&     symrec,symrel,tcpui,twalli,typat,&
&     ucvol,xred,zeff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_ifc
 use m_ddb_blk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'refineblk'
 use interfaces_42_geometry
 use interfaces_56_recipspace
 use interfaces_77_ddb, except_this_one => refineblk
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,msym,natom,nsym,ntypat
 integer,intent(in) :: iout
 real(dp),intent(in) :: tcpui,twalli,ucvol
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
!arrays
 integer,intent(in) :: indsym(4,nsym*natom)
 integer,intent(in) :: symrec(3,3,msym),symrel(3,3,msym),typat(natom)
 real(dp),intent(in) :: acell(3),amu(ntypat)
 real(dp),intent(in) :: gmet(3,3),gprim(3,3),rmet(3,3),rprim(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(inout) :: dielt(3,3),zeff(3,3,natom)

! these should already be associated in anaddb
 type(ddb_blk_type), pointer :: ddb_blk

!Local variables -------------------------

 integer :: iblok, nblok, iqpt
 integer :: nqpt_fine, nqpt_coarse
 integer :: option, mqpt, nqpt_fine_ibz
 integer :: nqshft

 integer :: rfphon(4), rfelfd(4), rfstrs(4), rftyp
 integer :: ngqpt_coarse(3)
 integer :: qptrlatt(3,3)
 integer, allocatable :: symafm(:)

 real(dp),parameter :: qtol=2.0d-8 ! copied from anaddb: should be in a module somewhere
 real(dp) :: qptrlen
 real(dp) :: qphnrm(3) = one
 real(dp) :: qphon(3,3)
 real(dp) :: rprimd(3,3)
 real(dp) :: dyewq0(3,3,natom)
 real(dp) :: rcan(3,natom),trans(3,natom)
 real(dp), allocatable :: spqpt_fine_ibz(:,:)
 real(dp), allocatable :: d2cart(:,:)
 real(dp), allocatable :: wtq(:)
 real(dp), allocatable :: q1shft(:,:)

 type(ifc_type) :: ifc_coarse

 type(ddb_blk_type), pointer :: ddb_tmp

! *********************************************************************

!Should implement Gaal-Nagy's algorithm in PRB <b>73</b> 014117.
 ngqpt_coarse(1) = anaddb_dtset%ngqpt(1)/anaddb_dtset%qrefine
 ngqpt_coarse(2) = anaddb_dtset%ngqpt(2)/anaddb_dtset%qrefine
 ngqpt_coarse(3) = anaddb_dtset%ngqpt(3)/anaddb_dtset%qrefine

 ABI_ALLOCATE(d2cart,(2,ddb_blk%msize))

!1) get IFC from coarse grid of qpoints, which is complete
!at this stage ddb_blk can be larger than the nqpt_coarse needed for this call
 call mkifc9(acell,amu,anaddb_dtset,ddb_blk,&
& dielt,dyewq0,gmet,gprim,ifc_coarse,indsym,iout,&
& mpert,msym,natom,ngqpt_coarse,nsym,ntypat,rcan,rmet,rprim,&
& symrec,symrel,tcpui,trans,twalli,typat,&
& ucvol,xred,zeff)

 nqpt_fine = anaddb_dtset%ngqpt(1)*anaddb_dtset%ngqpt(2)*anaddb_dtset%ngqpt(3)
 nqpt_coarse = nqpt_fine / anaddb_dtset%qrefine**3

!create new container that ddb_blk will be pointing to
!nqpt_fine is a maximal value - if we are using symops to reduce q we will have much less
!with nsym = 1 and elfd perturbations, might need to increase it though!
 ABI_ALLOCATE(ddb_tmp,)
 call create_ddb_blk(ddb_blk%msize, nqpt_fine, ddb_tmp)
 
!transfer old dyn mats to new container
 ddb_tmp%flg(:,1:ddb_blk%nblok)   = ddb_blk%flg(:,1:ddb_blk%nblok)
 ddb_tmp%typ(1:ddb_blk%nblok)     = ddb_blk%typ(1:ddb_blk%nblok)
 ddb_tmp%nrm(:,1:ddb_blk%nblok)   = ddb_blk%nrm(:,1:ddb_blk%nblok)
 ddb_tmp%qpt(:,1:ddb_blk%nblok)   = ddb_blk%qpt(:,1:ddb_blk%nblok)
 ddb_tmp%val(:,:,1:ddb_blk%nblok) = ddb_blk%val(:,:,1:ddb_blk%nblok)

!make list of irred q in fine grid.
 option=1
 qptrlatt(:,:)=0
 qptrlatt(1,1)=anaddb_dtset%ngqpt(1)
 qptrlatt(2,2)=anaddb_dtset%ngqpt(2)
 qptrlatt(3,3)=anaddb_dtset%ngqpt(3)
 mqpt=anaddb_dtset%nqshft * nqpt_fine
 if(anaddb_dtset%brav==2)mqpt=mqpt/2
 if(anaddb_dtset%brav==3)mqpt=mqpt/4

 ABI_ALLOCATE(symafm,(nsym))
 symafm = 1
 call mkrdim(acell,rprim,rprimd)

!get nqpt_fine_ibz
 ABI_ALLOCATE(spqpt_fine_ibz,(3,0))
 ABI_ALLOCATE(wtq,(0))
 nqshft = anaddb_dtset%nqshft
 ABI_ALLOCATE(q1shft,(3,nqshft))
 q1shft = anaddb_dtset%q1shft
 call getkgrid(1,iout,7,spqpt_fine_ibz,1,qptrlatt,qptrlen,&
& nsym,0,nqpt_fine_ibz,nqshft,nsym,rprimd,q1shft,&
& symafm,symrel,(/0,0,0/),wtq)
 ABI_DEALLOCATE(wtq)
 ABI_DEALLOCATE(spqpt_fine_ibz)
 ABI_DEALLOCATE(q1shft)

!get all ibz q
 ABI_ALLOCATE(spqpt_fine_ibz,(3,nqpt_fine_ibz))
 ABI_ALLOCATE(wtq,(nqpt_fine_ibz))
 ABI_ALLOCATE(q1shft,(3,nqshft))
 call getkgrid(1,iout,7,spqpt_fine_ibz,1,qptrlatt,qptrlen,&
& nsym,nqpt_fine_ibz,nqpt_fine_ibz,nqshft,nsym,rprimd,q1shft,&
& symafm,symrel,(/0,0,0/),wtq)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(wtq)
 ABI_DEALLOCATE(q1shft)

!defaults to search for phonon pert only
 rfphon(1:2)=1
 rfelfd(:)=0
 rfstrs(:)=0
 rftyp=anaddb_dtset%rfmeth

 nblok = ddb_blk%nblok

!update those we have no explicit data for in blkval, before going to main mkifc9 call in anaddb
 do iqpt = 1, nqpt_fine_ibz

   qphon = zero
   qphon (:,1) = spqpt_fine_ibz(:,iqpt)
!  check if it is already accounted for in ddb_tmp
   call gtblk9(ddb_tmp,iblok,mpert,natom,&
&   qphon,qphnrm,qtol,rfphon,rfstrs,rfelfd,rftyp)

   if (iblok /= 0) then
     write (std_out,*) 'q ', spqpt_fine_ibz(:,iqpt), ' is already in existing ddb'
     continue
   end if
   nblok = nblok+1

   ddb_tmp%flg(:,nblok) = 1 ! all elements will be present
   ddb_tmp%typ(nblok) = 1 ! default to non-stationary phonon - actually interpolated
   ddb_tmp%qpt(:, nblok) = zero
   ddb_tmp%qpt(1:3, nblok) = spqpt_fine_ibz(:,iqpt)
   ddb_tmp%nrm(:, nblok) = one

!  calculate dynamical matrices on fine grid ngqpt
!  Get d2cart using the interatomic forces and the
!  long-range coulomb interaction through Ewald summation
   call gtdyn9(acell,ifc_coarse%atmfrc,dielt,anaddb_dtset%dipdip,&
&   dyewq0,d2cart,gmet,gprim,mpert,natom,&
&   ifc_coarse%nrpt,qphnrm(1),spqpt_fine_ibz(:,iqpt),rmet,rprim,ifc_coarse%rpt,&
&   trans,ucvol,ifc_coarse%wghatm,xred,zeff)

   ddb_tmp%val(:,:, nblok) = d2cart(:,:)

 end do

!inform ddb_tmp that we have added information
 ddb_tmp%nblok = nblok 

!deallocate old space 
 call destroy_ddb_blk (ddb_blk)

!deallocate the pointer - it will be pointed onto ddb_tmp and deallocated at the end of anaddb
 ABI_DEALLOCATE(ddb_blk)
 ddb_blk => ddb_tmp

 ABI_DEALLOCATE(spqpt_fine_ibz)
 ABI_DEALLOCATE(d2cart)

 call destroy_ifc (ifc_coarse)

end subroutine refineblk
!!***
