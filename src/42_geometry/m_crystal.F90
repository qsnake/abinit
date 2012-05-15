!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_crystal
!! NAME
!! m_crystal
!!
!! FUNCTION
!! Module containing the definition of the crystal_structure data type and methods used to handle it. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG, YP)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_crystal

 use m_profiling

 use defs_basis
 use m_errors 

 use m_numeric_tools,  only : set2unit

 implicit none

 private 
!!***

!----------------------------------------------------------------------

!!****t* m_crystal/crystal_structure
!! NAME
!! crystal_structure
!! 
!! FUNCTION
!! Structure defining the unit cell (geometry, atomic positions and symmetry operations in real and reciprocal space)
!! 
!! SOURCE

 type,public :: crystal_structure

!scalars
  !integer :: point_group                    ! Point group
  !integer :: bravais,crystsys               ! Bravais lattice, Crystal system
  !integer :: nptsym                         ! No of point symmetries of the Bravais lattice
  !integer :: bravais(11)                    ! bravais(1)=iholohedry, bravais(2)=center
                                             ! bravais(3:11)=coordinates of rprim in the axes of the conventional
                                             ! bravais lattice (*2 if center/=0)
  !integer,pointer ptsymrel(:,:,:)
  !ptsymrel(3,3,nptsym)
  ! nptsym point-symmetry operations of the Bravais lattice in real space in terms of primitive translations.

  integer :: natom     
  ! Number of atoms

  integer :: nsym   
  ! Number of symmetry operations

  integer :: ntypat   
  ! Number of type of atoms

  !$integer :: ntypalch,ntyppure

  integer :: npsp  
  ! No. of pseudopotentials

  integer :: space_group 
  ! Space group

  integer :: timrev  
  ! TODO BE CAREFUL here, as the convention used in abinit is different.
  ! 1 => do not use time-reversal symmetry.
  ! 2 => take advantage of time-reversal symmetry. 

  real(dp) :: ucvol                          
  ! Real space unit cell volume.

  logical :: use_antiferro
  ! .TRUE. if AFM symmetries are present and used.

  logical :: has_inversion 
  ! .TRUE. if inversion symmetry is present

!arrays
  real(dp) :: angdeg(3)
  ! Angles among rprim (degree).

  real(dp) :: gmet(3,3)
  ! Reciprocal space metric ($\textrm{bohr}^{-2}$).

  real(dp) :: gprimd(3,3)
  ! Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)

  real(dp) :: rmet(3,3)
  ! Metric in real space.

  real(dp) :: rprimd(3,3)
  ! Direct lattice vectors, Bohr units.

  integer,pointer :: indsym(:,:,:)    SET2NULL
  ! indsym(4,nsym,natom)
  ! indirect indexing array for atoms, see symatm.F90.

  integer,pointer :: symafm(:)   SET2NULL
  ! symafm(nsym)
  ! (Anti)Ferromagnetic symmetries.

  integer,pointer :: symrec(:,:,:)   SET2NULL
  ! symrec(3,3,nsym)
  ! Symmetry operation in reciprocal space (reduced coordinates)

  integer,pointer :: symrel(:,:,:)   SET2NULL
  ! symrel(3,3,nsym)
  ! Symmetry operations in direct space (reduced coordinates).

  integer,pointer :: atindx(:)     SET2NULL
  integer,pointer :: atindx1(:)    SET2NULL
  ! atindx(natom), atindx1(natom)
  ! Index tables for atoms useful to treat atoms type after type.

  integer,pointer :: typat(:)    SET2NULL
  integer,pointer :: nattyp(:)   SET2NULL
  ! typat(natom), nattyp(ntypat)
  ! Type of each natom and number of atoms of each type.

  real(dp),pointer :: tnons(:,:)   SET2NULL
  ! tnons(3,nsym)
  ! Fractional translations (reduced coordinates)

  real(dp),pointer :: xcart(:,:)  SET2NULL
  ! xcart(3,natom)
  ! Cartesian coordinates.

  real(dp),pointer :: xred(:,:)   SET2NULL
  ! xred(3,natom)
  ! Reduced coordinates.

  real(dp),pointer :: spinrot(:,:)   SET2NULL
  ! spinrot(4,nsym)
  ! spinor rotation matrices.

! Useful quantities that might be added in the future
  !real(dp),pointer :: amu(:)                ! amu(ntypat)

  real(dp),pointer :: ziontypat(:)   SET2NULL
  ! ziontypat(ntypat)
  ! Charge of the pseudo-ion (No of valence electrons needed to screen exactly the pseudopotential).

  !real(dp),pointer :: znucltypat(:)         ! znucltypat(ntypat) from alchemy

  real(dp),pointer :: znucl(:)  SET2NULL
  ! znucl(npsp)
  ! Nuclear charge for each type of pseudopotential

  character(len=132),pointer :: title(:)   SET2NULL
   ! title(ntypat)
   ! The content of first line read from the psp file

 end type crystal_structure
!!***

!----------------------------------------------------------------------

 public :: init_crystal            ! Main Creation method.
 public :: nullify_crystal         ! Set all pointers to NULL.
 public :: destroy_crystal         ! Free the structure.
 public :: print_crystal           ! Print dimensions and basic info stored in the object
 public :: print_symmetries        ! Helper function to print symmetries in a nice format.
 public :: idx_spatial_inversion   ! Return the index of the spatial inversion, 0 if not present.
 public :: isymmorphic             ! Returns .TRUE. if space group is symmorphic.

 !integer,private,parameter :: SPATIAL_INVERSION = RESHAPE(/-1,0,0,0,-1,0,0,0,-1/),(/3.3/)

CONTAINS  !=========================================================================================================================
!!***

!!****f* m_crystal/init_crystal
!! NAME
!!  init_crystal 
!!
!! FUNCTION
!!  Initialize a crystal_structure data type.
!!  Ideally the routine should work in two different modes:
!!  Either the symmetries are directly supplied or the space group
!!  is determined starting from the definition of the unit cell.
!!  Only the first method is implemented, the second one should be
!!  a wrapper for the symmetry finder library. To implement the 
!!  second case I have to add additional entries in the object
!!  and I have also to pass an object describing the (optional) geometry builder.
!!
!! INPUTS
!!  natom=number of atom
!!  ntypat=number of type of atoms
!!  nsym=number of symmetry operations
!!  rprimd(3,3)=dimensional lattive vector (real space)
!!  typat(natom)=type of each atom
!!  xred(3,natom)=reduced coordinates of each atom
!!  symrel(3,3,nsym) [optional]=symmetry operations in real space
!!  space_group=Space group (0 if not available)
!!  tnons(3,nsym) [optional]=fractional Translations
!!  symafm(nsym) [optional]=  ferromagnetic symmetries
!!  remove_inv [optional]= if .TRUE. the inversion is removed from the set of symmetries
!!  timrev ==2 => take advantage of time-reversal symmetry
!!         ==1 ==> do not use time-reversal symmetry 
!!
!! OUTPUT
!!  Cryst<crystal_structure>= the object completely initialized.
!!
!! TODO
!!  Add additional entries in the class:
!!  1) Info on space and point group (generators?).
!!  2) alchemy
!!  3) masses and nuclear (pseudo&AE) charge
!!  4) forces stresses, velocities.
!!  5) constraints for the relaxation
!!  6) Likely I will need also info on the electric field and berryopt
!!
!! PARENTS
!!      elphon,m_crystal_io,vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_crystal(Cryst,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
& ziontypat,znucl,timrev,use_antiferro,remove_inv,title,&
& symrel,tnons,symafm) ! Optional
    
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_crystal'
 use interfaces_32_util
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,nsym,timrev,space_group,npsp
 type(crystal_structure),intent(inout) :: Cryst
 logical,intent(in) :: remove_inv,use_antiferro
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: xred(3,natom),rprimd(3,3),ziontypat(ntypat),znucl(npsp)
 real(dp),optional,intent(in) :: tnons(3,nsym)
 character(len=*),intent(in) :: title(ntypat) 

!Local variables-------------------------------
!scalars
 integer :: iat,indx,itypat,pinv,isym,nsym_noI
 real(dp) :: tolsym8,ucvol
 character(len=500) :: msg      
!arrays
 integer :: symrec(3,3),inversion(3,3)
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3)
 real(dp) :: spinrot(4)
 integer,pointer :: symrel_noI(:,:,:)
 integer,allocatable :: indsym(:,:,:)
 real(dp),pointer :: tnons_noI(:,:)
! *************************************************************************

 !@crystal_structure
 call nullify_crystal(Cryst)

 Cryst%natom  = natom 
 Cryst%ntypat = ntypat
 Cryst%npsp   = npsp

 Cryst%space_group = space_group

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 
 Cryst%ucvol  = ucvol
 Cryst%rprimd = rprimd 
 Cryst%rmet   = rmet
 Cryst%gmet   = gmet
 Cryst%gprimd = gprimd

 Cryst%angdeg(1)=ACOS(Cryst%rmet(2,3)/SQRT(Cryst%rmet(2,2)*Cryst%rmet(3,3)))/two_pi*360.0d0
 Cryst%angdeg(2)=ACOS(Cryst%rmet(1,3)/SQRT(Cryst%rmet(1,1)*Cryst%rmet(3,3)))/two_pi*360.0d0
 Cryst%angdeg(3)=ACOS(Cryst%rmet(1,2)/SQRT(Cryst%rmet(1,1)*Cryst%rmet(2,2)))/two_pi*360.0d0

 ABI_ALLOCATE(Cryst%typat,(natom))
 ABI_ALLOCATE(Cryst%xred,(3,natom))
 ABI_ALLOCATE(Cryst%xcart,(3,natom))
 ABI_ALLOCATE(Cryst%ziontypat,(ntypat))
 ABI_ALLOCATE(Cryst%znucl,(npsp))

 Cryst%typat     = typat 
 Cryst%xred      = xred 
 Cryst%ziontypat = ziontypat
 Cryst%znucl     = znucl

 call xredxcart(natom,1,rprimd,Cryst%xcart,Cryst%xred)

 ABI_ALLOCATE(Cryst%title,(ntypat))
 Cryst%title=title
 !
 ! === Generate index table of atoms, in order for them to be used type after type ===
 ABI_ALLOCATE(Cryst%atindx,(natom))
 ABI_ALLOCATE(Cryst%atindx1,(natom))
 ABI_ALLOCATE(Cryst%nattyp,(ntypat))

 indx=1
 do itypat=1,ntypat
   Cryst%nattyp(itypat)=0
   do iat=1,natom
     if (Cryst%typat(iat)==itypat) then
       Cryst%atindx (iat )=indx 
       Cryst%atindx1(indx)=iat
       indx=indx+1
       Cryst%nattyp(itypat)=Cryst%nattyp(itypat)+1
     end if
   end do
 end do

 Cryst%timrev = timrev
 call set2unit(inversion) ; inversion=-inversion

 if (PRESENT(symrel).and.PRESENT(tnons).and.PRESENT(symafm)) then 
  if (.not.remove_inv) then
   ! * Just a copy 
   Cryst%nsym= nsym
   ABI_ALLOCATE(Cryst%symrel,(3,3,nsym))
   ABI_ALLOCATE(Cryst%symrec,(3,3,nsym))
   ABI_ALLOCATE(Cryst%tnons,(3,nsym))
   ABI_ALLOCATE(Cryst%symafm,(nsym))
   Cryst%symrel=symrel 
   Cryst%tnons=tnons
   Cryst%symafm=symafm
   Cryst%use_antiferro = use_antiferro
   Cryst%has_inversion=.FALSE.
   do isym=1,nsym
    call mati3inv(symrel(:,:,isym),symrec)
    Cryst%symrec(:,:,isym)=symrec
    if (ALL(symrel(:,:,isym)==inversion)) Cryst%has_inversion=.TRUE. 
   end do
  else 
   ! * Remove inversion, just to be compatible with old GW implementation
   ! TODO should be removed!
   call remove_inversion(nsym,symrel,tnons,nsym_noI,symrel_noI,tnons_noI,pinv)
   Cryst%nsym=nsym_noI
   ABI_ALLOCATE(Cryst%symrel,(3,3,nsym_noI))
   ABI_ALLOCATE(Cryst%symrec,(3,3,nsym_noI))
   ABI_ALLOCATE(Cryst%tnons,(3,nsym_noI))
   ABI_ALLOCATE(Cryst%symafm,(nsym_noI))
   Cryst%symrel=symrel_noI
   Cryst%tnons=tnons_noI
   Cryst%has_inversion=.FALSE.
   if (ANY(symafm==-1)) then 
    msg = ' Solve the problem with inversion before adding ferromagnetic symmetries '
    MSG_BUG(msg)
   end if
   Cryst%symafm=1
   Cryst%use_antiferro=use_antiferro 
   do isym=1,nsym_noI
    call mati3inv(symrel_noI(:,:,isym),symrec)
    Cryst%symrec(:,:,isym)=symrec
   end do
   ABI_DEALLOCATE(symrel_noI)
   ABI_DEALLOCATE(tnons_noI)
  end if

 else
  ! * Find symmetries symrec,symrel,tnons,symafm
  ! TODO This should be a wrapper around the abinit library whose usage is not so straightforward
  MSG_BUG('not yet implemented')
 end if

 ! === Obtain a list of rotated atoms ===
 ! $ R^{-1} (xred(:,iat)-\tau) = xred(:,iat_sym) + R_0 $ 
 ! * indsym(4,  isym,iat) gives iat_sym in the original unit cell.
 ! * indsym(1:3,isym,iat) gives the lattice vector $R_0$.
 ! 
 ABI_ALLOCATE(indsym,(4,Cryst%nsym,natom))
 tolsym8=tol8
 call symatm(indsym,natom,Cryst%nsym,Cryst%symrec,Cryst%tnons,tolsym8,Cryst%typat,Cryst%xred)

 ABI_ALLOCATE(Cryst%indsym,(4,Cryst%nsym,natom))
 Cryst%indsym=indsym  
 ABI_DEALLOCATE(indsym)

 ! === Rotation in spinor space ===
 ABI_ALLOCATE(Cryst%spinrot,(4,Cryst%nsym))
 do isym=1,Cryst%nsym
  call getspinrot(Cryst%rprimd,spinrot,Cryst%symrel(:,:,isym))
  Cryst%spinrot(:,isym)=spinrot(:)
 end do

end subroutine init_crystal
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/nullify_crystal
!! NAME
!! nullify_crystal
!!
!! FUNCTION
!!  Nullify the pointers in the crystal_structure data type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_crystal
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine nullify_crystal(Cryst)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_crystal'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(crystal_structure),intent(inout) :: Cryst
! *********************************************************************

 !@crystal_structure

! integer
 nullify(Cryst%indsym )
 nullify(Cryst%symafm )
 nullify(Cryst%symrec )
 nullify(Cryst%symrel )
 nullify(Cryst%atindx )   
 nullify(Cryst%atindx1)   
 nullify(Cryst%typat  )   
 nullify(Cryst%nattyp )   

! real
 nullify(Cryst%tnons    )
 nullify(Cryst%xcart    )
 nullify(Cryst%xred     )
 nullify(Cryst%ziontypat)
 nullify(Cryst%znucl    )
 nullify(Cryst%spinrot  )

!character
 nullify(Cryst%title)

end subroutine nullify_crystal
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/destroy_crystal
!! NAME
!!  destroy_crystal 
!!
!! FUNCTION
!!  Destroy the dynamic arrays in a Crystal_structure data type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter,elphon,m_gwannier,m_screening,mlwfovlp_qp,mrgscr,rdm
!!      screening,sigma,vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_crystal(Cryst)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_crystal'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(crystal_structure),intent(inout) :: Cryst

! *********************************************************************

 DBG_ENTER("COLL")
 
 !@crystal_structure

!integer
 if (associated(Cryst%indsym ))  then
   ABI_DEALLOCATE(Cryst%indsym)
 end if
 if (associated(Cryst%symafm ))  then
   ABI_DEALLOCATE(Cryst%symafm)
 end if
 if (associated(Cryst%symrec ))  then
   ABI_DEALLOCATE(Cryst%symrec)
 end if
 if (associated(Cryst%symrel ))  then
   ABI_DEALLOCATE(Cryst%symrel)
 end if
 if (associated(Cryst%atindx ))  then
   ABI_DEALLOCATE(Cryst%atindx)
 end if
 if (associated(Cryst%atindx1))  then
   ABI_DEALLOCATE(Cryst%atindx1)
 end if
 if (associated(Cryst%typat  ))  then
   ABI_DEALLOCATE(Cryst%typat)
 end if
 if (associated(Cryst%nattyp ))  then
   ABI_DEALLOCATE(Cryst%nattyp)
 end if

!real
 if (associated(Cryst%tnons    ))  then
   ABI_DEALLOCATE(Cryst%tnons)
 end if
 if (associated(Cryst%xcart    ))  then
   ABI_DEALLOCATE(Cryst%xcart)
 end if
 if (associated(Cryst%xred     ))  then
   ABI_DEALLOCATE(Cryst%xred)
 end if
 if (associated(Cryst%ziontypat))  then
   ABI_DEALLOCATE(Cryst%ziontypat)
 end if
 if (associated(Cryst%znucl    ))  then
   ABI_DEALLOCATE(Cryst%znucl)
 end if
 if (associated(Cryst%spinrot  ))  then
   ABI_DEALLOCATE(Cryst%spinrot)
 end if

!character
 if (associated(Cryst%title))  then
   ABI_DEALLOCATE(Cryst%title)
 end if

 DBG_EXIT("COLL")

end subroutine destroy_crystal
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/print_crystal
!! NAME
!!  print_crystal 
!!
!! FUNCTION
!!  Print the content of crystal_structure data type
!!
!! INPUTS
!!  Cryst<crystal_structure>=The structure.
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS"
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing 
!!
!! PARENTS
!!      rdm,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_crystal(Cryst,header,unit,mode_paral,prtvol) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_crystal'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(crystal_structure),intent(in) :: Cryst

!Local variables-------------------------------
 integer :: my_unt,my_prtvol,nu
 character(len=4) :: my_mode
 character(len=500) :: msg      
! ********************************************************************* 

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the Cryst% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a)')' Real(R)+Recip(G) space primitive vectors, cartesian coordinates (Bohr,Bohr^-1):'
 call wrtout(my_unt,msg,my_mode)
 do nu=1,3
  write(msg,'(1x,a,i1,a,3f11.7,2x,a,i1,a,3f11.7)')&
&  'R(',nu,')=',Cryst%rprimd(:,nu)+tol10,&
&  'G(',nu,')=',Cryst%gprimd(:,nu)+tol10 !tol10 is used to be consistent with metric.F90
  call wrtout(my_unt,msg,my_mode)
 end do

 write(msg,'(a,1p,e15.7,a)')&
& ' Unit cell volume ucvol=',Cryst%ucvol+tol10,' bohr^3'
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a,3es16.8,a)')&
& ' Angles (23,13,12)=',Cryst%angdeg(1:3),' degrees'
 call wrtout(my_unt,msg,my_mode)

 if (Cryst%timrev==1) then 
  write(msg,'(a)')' Time-reversal symmetry is not present '
 else if (Cryst%timrev==2) then 
  write(msg,'(a)')' Time-reversal symmetry is present '
 else 
  MSG_BUG('Wrong value for timrev') 
 end if
 call wrtout(my_unt,msg,my_mode)

 !if (Cryst%use_antiferro) then 
 ! write(msg,'(a)')' System has magnetic symmetries '
 ! call wrtout(my_unt,msg,my_mode)
 !end if

 call print_symmetries(Cryst%nsym,Cryst%symrel,Cryst%tnons,Cryst%symafm,unit=my_unt,mode_paral=my_mode)

end subroutine print_crystal
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/print_symmetries
!! NAME
!! print_symmetries
!!
!! FUNCTION
!!  Helper function to print the set of symmetries.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gensymspgr,hdr_vs_dtset,m_crystal
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_symmetries(nsym,symrel,tnons,symafm,unit,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_symmetries'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,optional,intent(in) :: unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 integer,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
 integer :: my_unt,isym,isymin,isymend,ii,jj
 character(len=500) :: msg      
 character(len=4) :: my_mode
! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 !write(msg,'(2a)')ch10,' Rotations                           Translations     Symafm '
 !do isym=1,nsym
 ! write(msg,'(1x,3(3i3,1x),4x,3(f11.7,1x),6x,i2)')symrel(:,:,isym),tnons(:,isym),symafm(isym)
 ! call wrtout(my_unt,msg,my_mode)
 !end do 

 write(msg,'(2a)')ch10,' Symmetry operations in real space (Rotation tnons AFM)'
 call wrtout(my_unt,msg,my_mode)
 do isymin=1,nsym,4
  isymend=isymin+3
  if (isymend>nsym) isymend=nsym
  do ii=1,3
   write(msg,'(4(3i3,f8.3,i3,3x))')((symrel(ii,jj,isym),jj=1,3),tnons(ii,isym),symafm(isym),isym=isymin,isymend)
   call wrtout(my_unt,msg,my_mode)
  end do
  write(msg,'(a)')ch10
  call wrtout(my_unt,msg,my_mode)
 end do

end subroutine print_symmetries 
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/idx_spatial_inversion
!! NAME
!!  idx_spatial_inversion
!!
!! FUNCTION
!!  Return the index of the spatial inversion, 0 if not present
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function idx_spatial_inversion(Cryst) result(inv_idx)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'idx_spatial_inversion'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: inv_idx
 type(crystal_structure),intent(in) :: Cryst

!Local variables-------------------------------
!scalars
 integer :: isym
!arrays
 integer :: inversion(3,3)

! *************************************************************************

 inversion=RESHAPE((/-1,0,0,0,-1,0,0,0,-1/),(/3,3/))

 inv_idx=0
 do isym=1,Cryst%nsym
   if ( ALL(Cryst%symrel(:,:,isym)==inversion) ) then 
    inv_idx=isym; RETURN
   end if
 end do

end function idx_spatial_inversion
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/isymmorphic
!! NAME
!!  isymmorphic
!!
!! FUNCTION
!!  Returns .TRUE. is space group is symmorphic, i.e. all fractional translations are zero.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function isymmorphic(Cryst) result(ans)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isymmorphic'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: ans
 type(crystal_structure),intent(in) :: Cryst

! *************************************************************************

 ans = ALL (ABS(Cryst%tnons)<tol6) 

end function isymmorphic
!!***

!----------------------------------------------------------------------

END MODULE m_crystal
!!***
