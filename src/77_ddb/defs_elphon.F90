!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_elphon
!!
!! NAME
!! defs_elphon
!!
!! FUNCTION
!! This module contains the datastructures for elphon
!!  the different (huge) matrices will either be allocated and
!!  used, or be written to disk. All combinations should be feasible.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!  Contains two datastructures:
!!   elph_type contains data and dimensions for the kpoints near the
!!    fermi surface and the $g_{k k+q}$ matrix elements
!!   phon_type contains the necessary data to interpolate the
!!    phonon bandstructure and eigenvectors in reciprocal space (ie.
!!    interatomic force constants and corresponding real space grid info).
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

module defs_elphon

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_kptrank,  only : kptrank_type, nullify_kptrank, destroy_kptrank
 use m_crystal,  only : crystal_structure

 implicit none

 private 
!!***

!----------------------------------------------------------------------
!!****t* defs_elphon/elph_kgrid_type
!! NAME
!! elph_kgrid_type
!!
!! FUNCTION
!! elph_kgrid_type contains k-point grid data and dimensions
!!  this is a sub object of elph_type
!!
!! SOURCE

  type,public :: elph_kgrid_type

   integer :: nkpt                            ! number of k-points in full grid
   integer :: nkptirr                         ! number of k-points in irreducible grid

   type(kptrank_type) :: kptrank_t            ! ranking of all kpoints on phonon calculation grid, and inverse rank

   integer, pointer :: irr2full(:)      ! correspondence of irred kpoints to a full one
   integer, pointer :: full2irr(:,:)    ! correspondence of full k to one irred kpoints through sym and timrev
   integer, pointer :: full2full(:,:,:) ! symmetry mapping of kpoints 

   real(dp),pointer :: kpt(:,:)               ! coordinates of the full kpoints from phonon calculation
   real(dp),pointer :: kptirr(:,:)            ! irreducible k-points, for preliminary set up
   real(dp),pointer :: wtk(:,:,:)             ! integration weights (see also gkk_intweight)
   real(dp),pointer :: wtkirr(:)              ! weights for irreducible kpoints, to sum over _whole_ BZ (not just Fermi Surface)
  end type elph_kgrid_type
!!***

!----------------------------------------------------------------------

!!****t* defs_elphon/elph_type
!! NAME
!! elph_type
!!
!! FUNCTION
!! elph_type contains data and dimensions for the kpoints near the
!! fermi surface and the $g_{k k+q}$ matrix elements
!!
!! SOURCE

  type,public :: elph_type
    
   type(elph_kgrid_type) :: k_phon               ! object for k-grid of phonon calculation
   type(elph_kgrid_type) :: k_fine               ! object for fine k-grid for FS integration

   integer :: natom,nbranch,nFSband 
   integer :: minFSband,maxFSband                !Index of lower and upper bands used for the FS integration

   integer :: ngkkband                           !Number of bands kept in final gkk matrix elements:
                                                 !either 1 if sum is performed immediately
                                                 !or = nband if all elements are kept based on flag ep_keepbands

   integer :: nqpt_full                          !number of q in full BZ
   integer :: nqptirred                          !number of irred q-points


   integer :: unita2f,unit_gkk2,unit_gkk_rpt
   integer :: unitgkq                            !units for file output
     
   integer :: gkqwrite                              
   integer :: gkk2write
   integer :: gkk_rptwrite

   integer :: ep_scalprod                        !flag to perform the scalar product
   integer :: symgkq                             !flag to symmetrize gkq matrix elements
   integer :: ep_keepbands                       !flag to sum over bands or not
   integer :: tuniformgrid                       !flag to expect uniform grid of q or not

   
   integer :: na2f                               !dimensions and increments for a2F function
   integer :: nsppol                             ! number of spin polarization channels
   integer :: nspinor                            ! number of spinorial components
   real(dp) :: omega_min,omega_max 
   real(dp) :: a2fsmear,domega
   real(dp) :: nelect                            ! number of electrons per unit cell, eventually with extra charges for carriers in semiconductors.
   real(dp) :: occ_factor                        ! normalization for integrals over FS, for num of spins, spinors, etc...

   real(dp) :: mustar                            !mustar parameter
   real(dp) :: fermie                            !Fermi energy (Ha), either comes from wfk file or from anaddb input file

   character(len=fnlen) :: elph_base_name        !base name for output files

   integer,pointer :: qirredtofull(:)            !mapping between the qpoints found in the GGK file
                                                 !and the array of qpoints generated by the code

   real(dp),pointer :: wtq(:)                    !weight for each qpoint in the full grid spqt
                                                 !if a point is not in the IBZ ==>  wtq=0 
                                                 !MG we can also use indqpt

   real(dp),pointer :: n0(:)                     !DOS at the Fermi level (states/Ha/spin)
   real(dp),pointer :: qpt_full(:,:)             !special q points obtained by the Monkhorst & Pack method,
                                                 !in reduced coordinates


   real(dp),pointer :: gkk_intweight(:,:,:)      ! (nFSband,nkpt_fine,nsppol)
                                                 !integration weights for gkk matrix elements on FS:
                                                 !if ep_keepbands == 0 all are 1
                                                 !if ep_keepbands == 1 then = to wtk_phon in elphon
                                                 !DOES NOT INCLUDE FACTOR OF 1/nkpt_phon

   real(dp),pointer :: gkk_qpt(:,:,:,:,:,:)      ! (2, ngkkband*ngkkband, nbranch*nbranch, nkpt_phon, nsppol, nqptirred)
                                                 !Now gkq contains gkk2 matrices on basic qpts, 
                                                 !summed over bands if ngkkband==1


   real(dp),pointer :: gkk_rpt(:,:,:,:,:,:)      ! (2, ngkkband**2, nbranch**2, nkpt_phon, nsppol, nrpt)
                                                 !For the moment, gkk_rpt in memory is out of the question
   real(dp),pointer :: gkk2(:,:,:,:,:,:)         ! (nbranch, ngkkband,ngkkband, nkpt_phon, nkpt_phon, nsppol)

   real(dp),pointer :: gamma_qpt(:,:,:,:)        !gamma matrices integrated over kpoint coeff
                                                 !  and bands: still depends on qpt
                                                 ! dims= 2, nbranch**2, nsppol, nqpt
   real(dp),pointer :: gamma_rpt(:,:,:,:)
                                                 ! dims= 2, nbranch**2, nsppol, nrpt
!NOTE: choice to put nsppol before or after nqpt is a bit arbitrary
!   abinit uses nband,nkpt,nsppol, but here for convenience nkpt_phon,nsppol,nqpt 
!   as interpolation is on qpt

   real(dp),pointer :: phfrq(:,:)                !phonon frequencies
   real(dp),pointer :: a2f(:,:,:)                !a2f function

   real(dp),pointer :: qgrid_data(:,:,:,:)       !e-ph values calculated over the irreducible part of the q-grid:
                                                 !first entry  =  index of the q-point,
                                                 !second index =  branch index
                                                 !the third slice contains the frequency, the linewidth and lambda(q,nu)
                                                 !for that particular phonon mode
                                                 ! dims= nqptirred,elph_ds%nbranch,nsppol,3 

  end type elph_type

 public :: elph_ds_nullify
 public :: elph_ds_clean
!!***

!----------------------------------------------------------------------

!!****t* defs_elphon/phon_type
!! NAME
!! phon_type
!!
!! FUNCTION
!! phon_type contains the necessary data to interpolate the
!! phonon bandstructure and eigenvectors in reciprocal space 
!! (ie. interatomic force constants and corresponding real space grid info).
!!
!! SOURCE

 type,public :: phon_type

  integer :: dipdip
  ! dipole dipole interaction flag

  integer :: mpert 
  ! maximum number of ipert

  integer :: nsym
  integer :: natom 
  integer :: ntypat 

  integer :: nrpt
  ! Number of real space points used to integrate IFC (for interpolation of dynamical matrices)

  integer :: symdynmat
  ! If equal to 1, the dynamical matrix is symmetrized in phfrq3 before the diagonalization.

  real(dp) :: ucvol

  integer,pointer :: indsym(:,:,:)
  integer,pointer :: symrel(:,:,:)
  integer,pointer :: typat(:)

  real(dp),pointer :: acell(:)
  real(dp),pointer :: amu(:)

  real(dp),pointer :: atmfrc(:,:,:,:,:,:)
  ! atmfrc(2,3,natom,3,natom,nrpt)
  ! inter-atomic force constants.

  real(dp),pointer :: dielt(:,:)
  ! dielt(3,3)
  ! dielectric tensor.

  real(dp),pointer :: dyewq0(:,:,:)
  ! dyewq0(3,3,natom)
  ! atomic self-interaction correction to the dynamical matrix (only when dipdip=1).

  real(dp),pointer :: gprim(:,:)
  real(dp),pointer :: gmet(:,:)

  real(dp),pointer :: xred(:,:)
  real(dp),pointer :: zeff(:,:,:)

  real(dp),pointer :: rcan(:,:)
  ! rcan(3,natom) 
  ! canonical positions of atoms.


  real(dp),pointer :: rmet(:,:)
  real(dp),pointer :: rprim(:,:)
  real(dp),pointer :: rprimd(:,:)

  real(dp),pointer :: rpt(:,:)
  ! rpt(3,nrpt)
  ! Canonical positions of R points.

  real(dp),pointer :: trans(:,:)
  ! trans(3,natom) 
  ! Atomic translations : xred = rcan + trans

  real(dp),pointer :: wghatm(:,:,:)
  ! wghatm(natom,natom,nrpt)
  ! Weight for the pair of atoms and the R vector.

 end type phon_type
 
 public :: setup_phon_ds    ! Creation method.
 public :: clean_phon_ds    ! Destructor method (free memory).
 !%public :: inpphon          ! Interpolate phonon frequencies and eigenvectors at 1 qpt.
!!***

!----------------------------------------------------------------------

!!****t* defs_elphon/elph_tr_type
!! NAME
!! phon_type
!!
!! FUNCTION
!! elph_tr_ds contains the necessary data for the transport properties
!!
!! SOURCE

  type,public :: elph_tr_type

     integer :: ifltransport
     integer :: unitgkq_trin,unitgkq_trout
     integer :: gkqwrite,gkqexist
     integer :: onegkksize
     character(len=fnlen) :: ddkfilename
     real(dp),pointer :: el_veloc(:,:,:,:)        ! nkpt nband 3 nsppol
! the 9 = 3x3 is for the full tensorial transport coefficients
     real(dp),pointer :: eta_trin(:,:,:,:,:)      ! 9 bands**2 nFSkpt nsppol qpt 
     real(dp),pointer :: eta_trout(:,:,:,:,:) 

     real(dp),pointer :: gamma_qpt_tr(:,:,:,:,:)    ! 2 9 branches**2 nsppol qpt
     real(dp),pointer :: gamma_qpt_trin(:,:,:,:,:)  !idem
     real(dp),pointer :: gamma_qpt_trout(:,:,:,:,:) !idem

     real(dp),pointer :: gamma_rpt_tr(:,:,:,:,:)    !idem
     real(dp),pointer :: gamma_rpt_trin(:,:,:,:,:)  !idem
     real(dp),pointer :: gamma_rpt_trout(:,:,:,:,:) !idem

     real(dp),pointer :: a2f_1d_tr(:,:,:)           ! nfreq 9 nsppol
     real(dp),pointer :: a2f_1d_trin(:,:,:)
     real(dp),pointer :: a2f_1d_trout(:,:,:)

     real(dp),pointer ::  FSelecveloc_sq(:,:)       ! 3 nsppol

  end type elph_tr_type

 public :: elph_tr_ds_nullify
 public :: elph_tr_ds_clean
!!***

!----------------------------------------------------------------------

 ! Helper functions 
 ! TODO: Create a new low-level module.
 public :: phdispl_cart2red

!----------------------------------------------------------------------

CONTAINS  !=========================================================================================================================
!!***

!!****f* defs_elphon/setup_phon_ds
!! NAME
!! setup_phon_ds
!!
!! FUNCTION
!! This routine copies scalars and arrays into structure phon_ds
!!  to be passed later to phonon interpolation routines.
!!
!! INPUTS
!!   acell =  input length scales of cell (bohr)
!!   amu = mass of the atoms (atomic mass unit)
!!   atmfrc = inter-atomic force constants from anaddb
!!   dielt = dielectric tensor
!!   dipdip = dipole dipole interaction flag
!!   dyewq0 = atomic self-interaction correction to the dynamical matrix (only when dipdip=1)
!!   gmet = metric in reciprocal space (telphint=1)
!!   gprim =dimensionless basis vectors of reciprocal space
!!   indsym = mapping of atoms btw themselves under symmetry
!!   mpert = maximum number of ipert
!!   natom = number of atoms in cell
!!   nsym = number of space group symmetries
!!   ntypat = number of types of atoms
!!   nrpt =number of real space points used to integrate IFC (for interpolation of dynamical matrices)
!!   symdynmat=If equal to 1, the dynamical matrix is symmetrized in phfrq3 before the diagonalization.
!!   rcan = canonical positions of atoms
!!   rmet = metric tensor in real space (bohr^2)
!!   rprim =  primitive translation vectors (normalized)
!!   rprimd =  primitive translation vectors (dimensional)
!!   rpt = canonical positions of R points in the unit cell
!!   symrel = 3x3 matrices of the group symmetries (real space)
!!   trans = Atomic translations : xred = rcan + trans
!!   typat = type integer for each atom in cell
!!   ucvol = unit cell volume in bohr**3
!!   wghatm = Weight for the pair of atoms and the R vector
!!   xred = fractional dimensionless atomic coordinates
!!   zeff = effective charge on each atom, versus electric field and atomic displacement
!!
!! OUTPUT
!!   phon_ds = data structure for phonon interpolation - filled and allocated
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine setup_phon_ds(phon_ds,Cryst,dipdip,mpert,nrpt,symdynmat,&
&     ucvol,acell,amu,atmfrc,dielt,dyewq0,gprim,gmet,&
&     xred,zeff,rcan,rmet,rprim,rprimd,rpt,trans,wghatm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_phon_ds'
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: dipdip,mpert,nrpt,symdynmat
 real(dp),intent(in) :: ucvol
 type(crystal_structure),intent(in) :: Cryst
 type(phon_type),intent(inout) :: phon_ds
 !arrays
 real(dp),intent(in) :: acell(3),amu(Cryst%ntypat),atmfrc(2,3,Cryst%natom,3,Cryst%natom,nrpt)
 real(dp),intent(in) :: dielt(3,3),dyewq0(3,3,Cryst%natom),gmet(3,3),gprim(3,3)
 real(dp),intent(in) :: xred(3,Cryst%natom),zeff(3,3,Cryst%natom),rcan(3,Cryst%natom)
 real(dp),intent(in) :: rmet(3,3),rprim(3,3),rprimd(3,3),rpt(3,nrpt)
 real(dp),intent(in) :: trans(3,Cryst%natom),wghatm(Cryst%natom,Cryst%natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: natom,nsym

! *************************************************************************

 !@phon_type
 natom = Cryst%natom
 nsym  = Cryst%nsym
 !
 ! Copy basic parameters.
 phon_ds%dipdip = dipdip
 phon_ds%mpert  = mpert
 phon_ds%nsym   = Cryst%nsym
 phon_ds%natom  = natom
 phon_ds%ntypat = Cryst%ntypat
 phon_ds%nrpt   = nrpt
 phon_ds%symdynmat = symdynmat
 phon_ds%ucvol  = ucvol
 !
 ! Allocate arrays.
 ABI_ALLOCATE(phon_ds%indsym,(4,nsym,natom))
 phon_ds%indsym(:,:,:) = Cryst%indsym(:,:,:)

 ABI_ALLOCATE(phon_ds%symrel,(3,3,nsym))
 phon_ds%symrel(:,:,:) = Cryst%symrel(:,:,:)

 ABI_ALLOCATE(phon_ds%typat,(natom))
 phon_ds%typat(:) = Cryst%typat(:)

 ABI_ALLOCATE(phon_ds%acell,(3))
 phon_ds%acell(:) = acell(:)

 ABI_ALLOCATE(phon_ds%amu,(Cryst%ntypat))
 phon_ds%amu(:) = amu(:)

 ABI_ALLOCATE(phon_ds%atmfrc,(2,3,natom,3,natom,nrpt))
 phon_ds%atmfrc = atmfrc

 ABI_ALLOCATE(phon_ds%dielt,(3,3))
 phon_ds%dielt(:,:) = dielt(:,:)

 ABI_ALLOCATE(phon_ds%dyewq0,(3,3,natom))
 phon_ds%dyewq0(:,:,:) = dyewq0(:,:,:)

 ABI_ALLOCATE(phon_ds%gprim,(3,3))
 phon_ds%gprim(:,:) = gprim(:,:)

 ABI_ALLOCATE(phon_ds%gmet,(3,3))
 phon_ds%gmet(:,:) = gmet(:,:)

 ABI_ALLOCATE(phon_ds%xred,(3,natom))
 phon_ds%xred(:,:) = xred(:,:)

 ABI_ALLOCATE(phon_ds%zeff,(3,3,natom))
 phon_ds%zeff(:,:,:) = zeff(:,:,:)

 ABI_ALLOCATE(phon_ds%rcan,(3,natom))
 phon_ds%rcan(:,:) = rcan(:,:)

 ABI_ALLOCATE(phon_ds%rmet,(3,3))
 phon_ds%rmet(:,:) = rmet(:,:)

 ABI_ALLOCATE(phon_ds%rprim,(3,3))
 phon_ds%rprim(:,:) = rprim(:,:)

 ABI_ALLOCATE(phon_ds%rprimd,(3,3))
 phon_ds%rprimd(:,:) = rprimd(:,:)

 ABI_ALLOCATE(phon_ds%rpt,(3,nrpt))
 phon_ds%rpt(:,:) = rpt(:,:)

 ABI_ALLOCATE(phon_ds%trans,(3,natom))
 phon_ds%trans(:,:) = trans(:,:)

 ABI_ALLOCATE(phon_ds%wghatm,(natom,natom,nrpt))
 phon_ds%wghatm(:,:,:) = wghatm(:,:,:)

end subroutine setup_phon_ds
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/clean_phon_ds
!!
!! NAME
!! clean_phon_ds
!!
!! FUNCTION
!! This routine cleans structure phon_ds
!!
!! INPUTS
!!   phon_ds = data structure for phonon interpolation - filled and allocated
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE


subroutine clean_phon_ds(phon_ds)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clean_phon_ds'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(phon_type),intent(inout) :: phon_ds
!Local variables-------------------------------

! *************************************************************************

 !@phon_type
 ABI_DEALLOCATE(phon_ds%indsym)
 ABI_DEALLOCATE(phon_ds%symrel)
 ABI_DEALLOCATE(phon_ds%typat)
 ABI_DEALLOCATE(phon_ds%acell)
 ABI_DEALLOCATE(phon_ds%amu)
 ABI_DEALLOCATE(phon_ds%atmfrc)
 ABI_DEALLOCATE(phon_ds%dielt)
 ABI_DEALLOCATE(phon_ds%dyewq0)
 ABI_DEALLOCATE(phon_ds%gprim)
 ABI_DEALLOCATE(phon_ds%gmet)
 ABI_DEALLOCATE(phon_ds%xred)
 ABI_DEALLOCATE(phon_ds%zeff)
 ABI_DEALLOCATE(phon_ds%rcan)
 ABI_DEALLOCATE(phon_ds%rmet)
 ABI_DEALLOCATE(phon_ds%rprim)
 ABI_DEALLOCATE(phon_ds%rprimd)
 ABI_DEALLOCATE(phon_ds%rpt)
 ABI_DEALLOCATE(phon_ds%trans)
 ABI_DEALLOCATE(phon_ds%wghatm)

end subroutine clean_phon_ds
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/phdispl_cart2red
!! NAME
!!  phdispl_cart2red
!!
!! FUNCTION
!!  Calculates the displacement vectors for all branches in reduced coordinates.
!!  $ displ_red = displ_cart \cdot gprimd $ for each phonon branch.
!!
!! INPUTS
!!  natom=Number of atoms.
!!  gprimd(3,3)Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  displ_cart(2,3*natom,3*natom)=Phonon displacement in Cartesian coordinates.
!!
!! OUTPUT
!!  displ_red(2,3*natom,3*natom)=Phonon displacement in reduded coordinates.
!!
!! PARENTS
!!      m_gamma,mka2f,mkph_linwid,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine phdispl_cart2red(natom,gprimd,displ_cart,displ_red)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdispl_cart2red'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: displ_cart(2,3*natom,3*natom)
 real(dp),intent(out) :: displ_red(2,3*natom,3*natom)

!Local variables-------------------------
!scalars
 integer :: nbranch,jbranch,iatom,idir,ibranch,kdir,k1

! *************************************************************************

 displ_red = zero

 nbranch=3*natom

 do jbranch=1,nbranch
   !
   do iatom=1,natom
     do idir=1,3
       ibranch=idir+3*(iatom-1)
       do kdir=1,3
         k1 = kdir+3*(iatom-1)
         ! WARNING: could be non-transpose of rprimd matrix : to be checked.
         ! 23 june 2004: rprimd becomes gprimd
         ! could be gprim and then multiply by acell...
         ! Nope, checked and ok with gprimd 24 jun 2004
         displ_red(1,ibranch,jbranch) = displ_red(1,ibranch,jbranch) + &
&         gprimd(kdir,idir) * displ_cart(1,k1,jbranch)

         displ_red(2,ibranch,jbranch) = displ_red(2,ibranch,jbranch) + &
&         gprimd(kdir,idir) * displ_cart(2,k1,jbranch)

       end do !kdir
     end do !idir
   end do !iatom
   !
 end do !jbranch

end subroutine phdispl_cart2red
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_ds_clean
!!
!! NAME
!!   elph_ds_clean
!!
!! FUNCTION
!!   deallocate remaining arrays in the elph_ds datastructure
!!
!! INPUTS
!!  elph_ds = elphon datastructure
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine elph_ds_clean(elph_ds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_ds_clean'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_type), intent(inout) :: elph_ds

! *************************************************************************

 !@elph_type
 if (associated(elph_ds%qirredtofull))  then
   ABI_DEALLOCATE(elph_ds%qirredtofull)
 end if
 if (associated(elph_ds%wtq))  then
   ABI_DEALLOCATE(elph_ds%wtq)
 end if
 if (associated(elph_ds%n0))  then
   ABI_DEALLOCATE(elph_ds%n0)
 end if
 if (associated(elph_ds%qpt_full))  then
   ABI_DEALLOCATE(elph_ds%qpt_full)
 end if
 if (associated(elph_ds%gkk_intweight))  then
   ABI_DEALLOCATE(elph_ds%gkk_intweight)
 end if
 if (associated(elph_ds%gkk_qpt))  then
   ABI_DEALLOCATE(elph_ds%gkk_qpt)
 end if
 if (associated(elph_ds%gkk_rpt))  then
   ABI_DEALLOCATE(elph_ds%gkk_rpt)
 end if
 if (associated(elph_ds%gkk2))  then
   ABI_DEALLOCATE(elph_ds%gkk2)
 end if
 if (associated(elph_ds%gamma_qpt))  then
   ABI_DEALLOCATE(elph_ds%gamma_qpt)
 end if
 if (associated(elph_ds%gamma_rpt))  then
   ABI_DEALLOCATE(elph_ds%gamma_rpt)
 end if
 if (associated(elph_ds%phfrq))  then
   ABI_DEALLOCATE(elph_ds%phfrq)
 end if
 if (associated(elph_ds%a2f))  then
   ABI_DEALLOCATE(elph_ds%a2f)
 end if
 if (associated(elph_ds%qgrid_data))  then
   ABI_DEALLOCATE(elph_ds%qgrid_data)
 end if

 if (associated(elph_ds%k_phon%wtk))  then
   ABI_DEALLOCATE(elph_ds%k_phon%wtk)
 end if
 if (associated(elph_ds%k_phon%wtkirr))  then
   ABI_DEALLOCATE(elph_ds%k_phon%wtkirr)
 end if
 if (associated(elph_ds%k_phon%kpt))  then
   ABI_DEALLOCATE(elph_ds%k_phon%kpt)
 end if
 if (associated(elph_ds%k_phon%kptirr))  then
   ABI_DEALLOCATE(elph_ds%k_phon%kptirr)
 end if
 if (associated(elph_ds%k_phon%irr2full))  then
   ABI_DEALLOCATE(elph_ds%k_phon%irr2full)
 end if
 if (associated(elph_ds%k_phon%full2irr))  then
   ABI_DEALLOCATE(elph_ds%k_phon%full2irr)
 end if
 if (associated(elph_ds%k_phon%full2full))  then
   ABI_DEALLOCATE(elph_ds%k_phon%full2full)
 end if

 call destroy_kptrank (elph_ds%k_phon%kptrank_t)

 if (associated(elph_ds%k_fine%wtk))  then
   ABI_DEALLOCATE(elph_ds%k_fine%wtk)
 end if
 if (associated(elph_ds%k_fine%wtkirr))  then
   ABI_DEALLOCATE(elph_ds%k_fine%wtkirr)
 end if
 if (associated(elph_ds%k_fine%kpt))  then
   ABI_DEALLOCATE(elph_ds%k_fine%kpt)
 end if
 if (associated(elph_ds%k_fine%kptirr))  then
   ABI_DEALLOCATE(elph_ds%k_fine%kptirr)
 end if
 if (associated(elph_ds%k_fine%irr2full))  then
   ABI_DEALLOCATE(elph_ds%k_fine%irr2full)
 end if
 if (associated(elph_ds%k_fine%full2irr))  then
   ABI_DEALLOCATE(elph_ds%k_fine%full2irr)
 end if
 if (associated(elph_ds%k_fine%full2full))  then
   ABI_DEALLOCATE(elph_ds%k_fine%full2full)
 end if

 call destroy_kptrank (elph_ds%k_fine%kptrank_t)

end subroutine elph_ds_clean
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_tr_ds_clean
!!
!! NAME
!!   elph_tr_ds_clean
!!
!! FUNCTION
!!   deallocate remaining arrays in the elph_tr_ds datastructure
!!
!! INPUTS
!!  elph_tr_ds = elphon transport datastructure
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine elph_tr_ds_clean(elph_tr_ds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_tr_ds_clean'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 
! *************************************************************************

 !@elph_tr_type
 if (associated(elph_tr_ds%el_veloc))  then
   ABI_DEALLOCATE(elph_tr_ds%el_veloc)
 end if
 if (associated(elph_tr_ds%eta_trin))  then
   ABI_DEALLOCATE(elph_tr_ds%eta_trin)
 end if
 if (associated(elph_tr_ds%eta_trout))  then
   ABI_DEALLOCATE(elph_tr_ds%eta_trout)
 end if
 if (associated(elph_tr_ds%FSelecveloc_sq))  then
   ABI_DEALLOCATE(elph_tr_ds%FSelecveloc_sq)
 end if
 if (associated(elph_tr_ds%gamma_qpt_tr))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_tr)
 end if
 if (associated(elph_tr_ds%gamma_qpt_trin))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_trin)
 end if
 if (associated(elph_tr_ds%gamma_qpt_trout))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_trout)
 end if
 if (associated(elph_tr_ds%gamma_rpt_tr))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_tr)
 end if
 if (associated(elph_tr_ds%gamma_rpt_trin))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_trin)
 end if
 if (associated(elph_tr_ds%gamma_rpt_trout))  then
   ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_trout)
 end if
 if (associated(elph_tr_ds%a2f_1d_tr))  then
   ABI_DEALLOCATE(elph_tr_ds%a2f_1d_tr)
 end if
 if (associated(elph_tr_ds%a2f_1d_trin))  then
   ABI_DEALLOCATE(elph_tr_ds%a2f_1d_trin)
 end if
 if (associated(elph_tr_ds%a2f_1d_trout))  then
   ABI_DEALLOCATE(elph_tr_ds%a2f_1d_trout)
 end if

end subroutine elph_tr_ds_clean
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_ds_nullify
!!
!! NAME
!!   elph_ds_nullify
!!
!! FUNCTION
!!   nullify all arrays in the elph_ds datastructure
!!
!! INPUTS
!!  elph_ds = elphon datastructure
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine elph_ds_nullify(elph_ds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_ds_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_type), intent(inout) :: elph_ds

! *************************************************************************

 !@elph_tr_type
 nullify(elph_ds%qirredtofull)
 nullify(elph_ds%wtq)
 nullify(elph_ds%n0)
 nullify(elph_ds%qpt_full)
 nullify(elph_ds%gkk_intweight)
 nullify(elph_ds%gkk_qpt)
 nullify(elph_ds%gkk_rpt)
 nullify(elph_ds%gkk2)
 nullify(elph_ds%gamma_qpt)
 nullify(elph_ds%gamma_rpt)
 nullify(elph_ds%phfrq)
 nullify(elph_ds%a2f)
 nullify(elph_ds%qgrid_data)

 call elph_k_nullify (elph_ds%k_phon)
 call elph_k_nullify (elph_ds%k_fine)

end subroutine elph_ds_nullify
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_k_nullify
!!
!! NAME
!!   elph_k_nullify
!!
!! FUNCTION
!!   nullify all arrays in the elph_k datastructure
!!
!! INPUTS
!!  elph_k = elphon k-points datastructure
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      defs_elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine elph_k_nullify(elph_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_k_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(inout) :: elph_k

! *************************************************************************

 !@elph_kgrid_type
 nullify(elph_k%irr2full)
 nullify(elph_k%full2irr)
 nullify(elph_k%full2full)
 nullify(elph_k%kpt)
 nullify(elph_k%kptirr)
 nullify(elph_k%wtk)
 nullify(elph_k%wtkirr)
 call nullify_kptrank (elph_k%kptrank_t)

end subroutine elph_k_nullify
!!***

!----------------------------------------------------------------------

!!****f* defs_elphon/elph_tr_ds_nullify
!!
!! NAME
!!   elph_tr_ds_nullify
!!
!! FUNCTION
!!   deallocate remaining arrays in the elph_tr_ds datastructure
!!
!! INPUTS
!!  elph_tr_ds = elphon transport datastructure
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine elph_tr_ds_nullify(elph_tr_ds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph_tr_ds_nullify'
!End of the abilint section

 implicit none


!Arguments ------------------------------------
!scalars
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 
! *************************************************************************

 !@elph_tr_type
 nullify(elph_tr_ds%el_veloc)
 nullify(elph_tr_ds%eta_trin)
 nullify(elph_tr_ds%eta_trout)
 nullify(elph_tr_ds%FSelecveloc_sq)
 nullify(elph_tr_ds%gamma_qpt_tr)
 nullify(elph_tr_ds%gamma_qpt_trin)
 nullify(elph_tr_ds%gamma_qpt_trout)
 nullify(elph_tr_ds%gamma_rpt_tr)
 nullify(elph_tr_ds%gamma_rpt_trin)
 nullify(elph_tr_ds%gamma_rpt_trout)
 nullify(elph_tr_ds%a2f_1d_tr)
 nullify(elph_tr_ds%a2f_1d_trin)
 nullify(elph_tr_ds%a2f_1d_trout)

end subroutine elph_tr_ds_nullify
!!***

!----------------------------------------------------------------------

end module defs_elphon
!!***

