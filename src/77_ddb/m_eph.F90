!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_eph
!! NAME
!!  m_eph
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods 
!!  used to perform electron-phonon calculations. 
!! ***** Working in progress, do not rely on this implementation. type declarations might be changed *****
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG,MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! NOTES
!!
!! Most important requirements of the new e-ph implementation:
!!
!! 1) Reduce the memory requirements:
!!    * Memory distribution via MPI.
!!    * Only states in an energy window close to E_F should be stored
!!    * Only the irreducible k-points for a given q and (ipert,idir) should be stored.
!!      Specialized routines will be used to retrieve the symmetrized matrix elements in the full BZ.
!!    * Only the irreducible pertubations at given (q, idir, ipert) should be stored and 
!!      explicitly calculated at the DFPT level.
!!
!! 2) Methods hiding the ugly details of the implementation.
!!
!! 3) out-of-core and in-core solutions (ETSF-IO support should facilitate the implementation but
!!    final specifications are still missing)
!!
!! 4) The objects should be designed taking int account a possible use of the Wannier interpolation.
!!    Since non-homogeneous k-meshes lead to a dramatic increase in the number of k-points, all the
!!    arrays dimensioned with nkibz and nkbz should be carefully designed to reduce memory.
!!
!! Issues:
!!
!! 1) Several quantities depend on the spin : Fermi surface, E_f, eigen, occ, gkk, set of bands around E_f
!!
!! 2) Spinorial case: each gkk is a two-by-two matrix.
!!
!! 3) Only those symmetries which preserve q, and (iatom, idir) can be used to reconstruct the gkk's in the full BZ.
!!    ==> for given q and (idir,ipert) we have IBZ_q^{iatom,idir}, 
!!
!! 4) To facilitate the implementation, the new E-PH part will make use of several data types already introduced in 
!!    the GW part:
!!     crystal_structure 
!!     band_structure
!!     bz_mesh_type
!!
!!    TODO: Better integration between GW data types and the main abinit code is needed
!!    e.g treatment of k-points, time-reversal, Fermi level for semiconductors...
!!
!! 5) To speed up the k-point search either we use the previous approach (rank and invrank)
!!    or we order the points according to their length. The later approach has the advantage
!!    of being applicable also to non-homogeneous k-meshes. Moreover it is less memory demanding 
!!    as the storage of invrank is not needed anymre. On the other hand, looping over shells  
!!    is expected to be a bit slower than the algorithm based on the k-point rank.
!!
!! 6) Which representation for the gkk? The one presently used i.e. (q,idir,ipert) 
!!    or the phonon representation (q,nu)?
!!    One should check whether (q,nu) leads to some simplification when symmetries are used 
!!    to complete the matrix elements.
!!
!! Questions:
!! 1) At Gamma and border zone, \Delta V_{SCF} is Hermitian thus we can store the upper triangle 
!!    of the (b1,b2) matrix. Is it worth taking into account this possibility?
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

MODULE m_eph

 use m_profiling

 use defs_basis
 use m_errors
 use defs_elphon

 implicit none

 private 
!!***

!!****t* m_eph/fermi_surface_type
!! NAME
!! fermi_surface_type
!!
!! FUNCTION
!!  Stores data and table related to the Fermi surface for a single spin polarization.
!!  Only the set of k-points where at least one band falls within an energy 
!!  window centered at Ef are stored in memory.
!!
!! SOURCE

 type,public :: fermi_surface_type

  integer :: nkbz
  ! Number of k-points on the Fermi surface in the full Brillouin zone (FS BZ).

  integer :: nkibz
  ! Number of k-points on the Fermi surface in the irreducible Brillouin zone (FS IBZ).

  integer :: nband
  ! Number of bands around the Fermi level falling within the energy window [ef-ewidth, ef+ewidth].

  integer :: bstart
  ! Initial band index.

  integer :: nsym
  ! Number of symmetries in the full point group.

  integer :: timrev
  ! 1 or 2 depending whether time-reversal cannot or can be used (respectively). 
  ! TODO use abinit conventions.

  integer :: occopt
  ! Option for metallic occupation.

  real(dp) :: efermi
  ! Fermi level.

  real(dp) :: ewidth
  ! Width of the energy window around Ef

  real(dp) :: tsmear
  ! Value of the broadening for occupation factors.

  real(dp) :: tphysel
  ! Physical temperature.

  real(dp) :: dosef 
  ! Dos at the Fermi level for this spin channel.

  real(dp) :: tol_weight
  ! Tolerance on the weight for FS integration. 
  ! A (k,b) state whose weight is smaller than tolweight won"t be stored.

  integer :: mesh_type
  ! 1 for homogeneous mesh
  ! 2 for random mesh (e.g. for Wannier interpolation)

  integer :: int_opt
  ! Flag defining the technique used to perform the integration of the FS.
  !  1 for gaussian.
  !  2 for tetrahedrons.

  real(dp) :: elphsmear
  ! Standard deviation of the gaussian function used to integrate over the FS. Only used if int_opt=1

  integer :: kptrlatt(3,3)=RESHAPE((/0,0,0,0,0,0,0,0,0/),(/3,3/))
  ! Coordinates of three vectors in real space, expressed in reduced coordinates.
  ! They defines a super-lattice in real space. The k point lattice is the reciprocal of
  ! this super-lattice, eventually shifted by shift. ONLY used if mesh_type == 1

  ! real(dp),pointer :: shift(3) !shift(:,:)
   !  shift(3,nshift)
   !  shift for k-points, usually nshift=1

  integer,pointer :: bz2ibz(:)  SET2NULL
  ! bz2ibz(nkpt) 
  ! For each point in the FS BZ, it gives the index of the symmetrical image in the IBZ.

  !integer,pointer :: symopbz2ibz(:)
  ! symopbz2ibz(nkbz)
  ! Index of the symmetry operation such that IS k_ibz = k_bz
  ! where S is one of the symrec ops and I is the identity or the inversion, depending on symtrbz2ibz.

  !integer,pointer :: symtrbz2ibz(:)
  ! symtrbz2ibz(nkbz)
  ! Index of the symmetry operation such that IS k_ibz = k_bz

  !integer, pointer :: kbz2sh(:)
  ! kbz2sh(nkbz)
  ! For each k-point in the FS BZ, it gives the index of the shell to which it belongs.

  
  !integer,pointer :: shlim(:)
  ! shlim(nksh+1)
  ! Index of the first k-point in each shell, =nkbz+1 for nksh+1

  !integer :: nksh                             ! Number of shells of k-points.

  !real(dp),pointer :: shlen(:)
  ! shlen(nksh)
  ! Radius of each shell.

  real(dp) :: gprimd(3,3)
   ! Dimensional reciprocal space primitive translations (Bohr^-1)

  real(dp),pointer :: kibz(:,:)   SET2NULL
  ! kibz(3,nkibz)
  ! Reduced coordinates of the points in the FS IBZ. Ordered by increasing module.

  real(dp),pointer :: kbz(:,:)   SET2NULL
  ! kbz(3,nkbz)
  ! Reduced coordinates of the points in the FS. Ordered by increasing module.
  ! TODO is it really needed? It might be large. Might be recalculated on-the-fly
  
  real(dp),pointer :: eigen(:,:)   SET2NULL
  ! eigen(nband,nkibz)
  ! Energies in the FS IBZ.

  real(dp),pointer :: occ(:,:)   SET2NULL
  ! occ(nband,nkibz)
  ! Occupation number in the FS IBZ.

  real(dp),pointer :: fs_weight(:,:)  SET2NULL
  ! fs_weight(nband,nkibz)
  ! Weights due to the delta function centered at the Fermi level.
  ! Two methods are available, standard gaussian and tetrahedron methods, depending on int_opt.

  real(dp),pointer :: wtk(:)   SET2NULL
  ! wtk(nkibz)
  ! Weights for each point on the FS IBZ. Normalized to one.

 end type fermi_surface_type
!!***
 
 ! Bound methods:
 public :: nullify_fermi_surface
 public :: destroy_fermi_surface

 !public :: init_fermi_surface
 !public :: wannier_interpolate_fermi_surface
 !public :: get_fs_ibz
 !public :: bxsf_write_fermi_surface

 ! example:
 !type(fermi_surface_type),allocatable :: Fsurf(:)
 !
 !allocate(Fsurf(nsppol))
 !do isppol=1,nsppol
 ! call init_Fermi_surface(Fsurf(isppol),Cryst,Kmesh,Bst,tolweight)
 !end do

 interface nullify_fermi_surface
  module procedure nullify_fermi_surface_0D
  module procedure nullify_fermi_surface_1D
 end interface nullify_fermi_surface

 interface destroy_fermi_surface
  module procedure destroy_fermi_surface_0D
  module procedure destroy_fermi_surface_1D
 end interface destroy_fermi_surface

!----------------------------------------------------------------------

!!****t* m_eph/gkk_type
!! NAME
!! gkk_type
!!
!! FUNCTION
!! Structure used to store the matrix elements at fixed (q, idir, ipert)
!!
!! SOURCE

 type,public :: gkk_type

  integer :: nkibz_gkk
  ! Number of points in the IBZ_q^{idir,ipert}

  integer :: fs_nkbz
  ! Numbe of points on the FS BZ.

  integer :: fs_nband
  ! Number of states around E_f. Equivalent to the value stored in the fermi_surface_type.

  integer :: bstart
  ! First band index. Equivalent to the value stored in the fermi_surface_type.

  integer :: nsym_gkk
  ! Number of symmetries preserving (q,idir,ipert)

  integer :: timrev
  ! 1 or 2 depending whether time-reversal cannot or can be used (respectively). 
  ! TODO use abinit conventions.

  integer,pointer :: symgkk2symrec(:)   SET2NULL
  ! symgkk2symrec(nsym_gkk)
  ! Index in the full array symrec of the symmetries preserving (qpt,idir,iper)

  integer,pointer :: ibzgkk2fsbz(:)    SET2NULL
  ! ibzgkk2fsbz(nkibz_gkk)
  ! Index in the full FS BZ of each point in the IBZ_q^{idir,ipert}
  
  integer,pointer :: fsbz2ibzgkk(:)   SET2NULL
  ! fsbz2ibzgkk(fs_nkbz)
  ! For each point in the full BZ gives the index in my IBZ of the symmetrical point. 

  integer,pointer :: symtrbz2ibz(:)   SET2NULL
  ! symtrbz2ibz(fs_nkbz)
  ! 1 or 2, depending wheter time-reversal has to be used to obtain the point.

  integer,pointer :: symopbz2ibz(:)  SET2NULL
  ! ksymopbz2ibz(fs_nkbz)
  ! Index of the symmetry operation in the set of my symmetries such that I S kibz = kbz
  ! TODO Do we need unklapp vectors

  real(dp),pointer :: gkk(:,:,:,:)  SET2NULL
  ! gkk(2,fs_nband,fs_nband,nkibz_gkk)

  !integer,pointer :: symrec_gkk(:,:,:)
  ! symrec_gkk(3,3,nsym_gkk)
  ! symmetry operations in reciprocal space preserving the external q and the perturbation (idir,ipert)
                                                                                                  
  !integer,pointer :: symafm_gkk(:)
  ! symafm_gkk(nsym_gkk)

 end type gkk_type
!!***

! Bound Methods: 
 public :: nullify_gkk
 public :: destroy_gkk
 !! init_gkk
 !! read_gkk_from_file
 !! get_gkk_full_fsbz        ! complete gkk on the full FS BZ.

 interface nullify_gkk
  module procedure nullify_gkk_0D
  module procedure nullify_gkk_1D
 end interface nullify_gkk
                                            
 interface destroy_gkk
  module procedure destroy_gkk_0D
  module procedure destroy_gkk_1D
 end interface destroy_gkk

!----------------------------------------------------------------------

!!****t* m_eph/gkk_handler_type
!! NAME
!! gkk_handler_type
!!
!! FUNCTION
!! Structure used to store the matrix elements for a given q.
!!
!! SOURCE

 type,public :: gkk_handler_type

  integer :: natom
  ! Number of atoms in the unit cell.

  integer :: nirr_perts
  ! Number of irreducible perturbations (idir,ipert) <= 3*natom.

  character(len=fnlen) :: fname
  ! Name of the file storing the gkk matrix elements (either complete database of single q-point and spin index).
  ! Used in the case of out-of-core solution.

  integer,pointer :: pert_list(:,:)    SET2NULL
  ! pert_list(2,nirr_perts)
  ! gives (idir,ipert) for each irreducible perturbation.

  real(dp) :: qpt(3)
  ! The q-point in reduced coordinates.

  type(gkk_type),pointer :: Gkk(:)   SET2NULL
  ! gkk(nirr_perts)
  ! gkk matrix elements for this q-point 

 end type gkk_handler_type 
!!***

! Bound Methods: 
 public :: nullify_gkk_handler
 public :: destroy_gkk_handler
 !% init_gkk_handler(Gkk,FSurf,Cryst,Cryst,qpt,fname)
 !% get_gammaq
 !% symmetrize_gkk_over_perts

 interface nullify_gkk_handler
  module procedure nullify_gkk_handler_0D
  module procedure nullify_gkk_handler_1D
 end interface nullify_gkk_handler
                                            
 interface destroy_gkk_handler
  module procedure destroy_gkk_handler_0D
  module procedure destroy_gkk_handler_1D
 end interface destroy_gkk_handler

 !type(gkk_handler_type),allocatable :: Gkk(:,:)
 !allocate(Gkk(nqibz,nsppol))
 !call nullify_gkk_handler(Gkk))

 !do isppol=1,nsppol
 ! do iqibz=1,nqibz
 !  if (I_treat(iqibz,isppol)) then
 !   call init_gkk_handler(Gkk(iqibz,isppol),Fsurf(isppol),Cryst,qpt,fname)
 !  end if
 ! end do
 !end do

CONTAINS  !===========================================================
!!***

!!****f* m_eph/nullify_fermi_surface_0D
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the fermi_surface_type structure, by nullifying them
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_fermi_surface_0D(FSurf)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_fermi_surface_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(fermi_surface_type),intent(inout) :: FSurf

! ************************************************************************

 !@fermi_surface_type

 ! integer 
 nullify(Fsurf%bz2ibz      )
 !nullify(Fsurf%symopbz2ibz)
 !nullify(Fsurf%symtrbz2ibz)
 !nullify(Fsurf%kbz2sh)
 !nullify(Fsurf%shlim)
 !nullify(Fsurf%shlen)

 ! real
 nullify(Fsurf%kibz     )
 nullify(Fsurf%kbz     )
 nullify(Fsurf%eigen    )
 nullify(Fsurf%occ      )
 nullify(Fsurf%fs_weight)
 nullify(Fsurf%wtk      )

end subroutine nullify_fermi_surface_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/nullify_fermi_surface_1D
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the fermi_surface_type structure, by nullifying them
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_fermi_surface_1D(FSurf)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_fermi_surface_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(fermi_surface_type),intent(inout) :: FSurf(:)

!Local variables-------------------------------
!scalars
 integer :: isppol
!arrays

! ************************************************************************

 do isppol=1,SIZE(FSurf)
  call nullify_fermi_surface_0D(FSurf(isppol))
 end do

end subroutine nullify_fermi_surface_1D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/destroy_fermi_surface_0D
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the fermi_surface_type structure
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_fermi_surface_0D(FSurf)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_fermi_surface_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(fermi_surface_type),intent(inout) :: FSurf

! ************************************************************************

 !@fermi_surface_type

 ! integer 
 if (associated(Fsurf%bz2ibz      ))  then
   ABI_DEALLOCATE(Fsurf%bz2ibz)
 end if
 !if (associated(Fsurf%symopbz2ibz)) !deallocate(Fsurf%symopbz2ibz)
 !if (associated(Fsurf%symtrbz2ibz)) !deallocate(Fsurf%symtrbz2ibz)
 !if (associated(Fsurf%kbz2sh       )) !deallocate(Fsurf%kbz2sh       )
 !if (associated(Fsurf%shlim      )) !deallocate(Fsurf%shlim      )
 !if (associated(Fsurf%shlen      )) !deallocate(Fsurf%shlen      )

 ! real
 if (associated(Fsurf%kibz     ))  then
   ABI_DEALLOCATE(Fsurf%kibz)
 end if
 if (associated(Fsurf%kbz      ))  then
   ABI_DEALLOCATE(Fsurf%kbz)
 end if
 if (associated(Fsurf%eigen    ))  then
   ABI_DEALLOCATE(Fsurf%eigen)
 end if
 if (associated(Fsurf%occ      ))  then
   ABI_DEALLOCATE(Fsurf%occ)
 end if
 if (associated(Fsurf%fs_weight))  then
   ABI_DEALLOCATE(Fsurf%fs_weight)
 end if
 if (associated(Fsurf%wtk      ))  then
   ABI_DEALLOCATE(Fsurf%wtk)
 end if

end subroutine destroy_fermi_surface_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/destroy_fermi_surface_1D 
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the fermi_surface_type structure
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_fermi_surface_1D(FSurf)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_fermi_surface_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(fermi_surface_type),intent(inout) :: FSurf(:)

!Local variables-------------------------------
!scalars
 integer :: isppol

! ************************************************************************

 do isppol=1,SIZE(FSurf)
  call destroy_fermi_surface_0D(FSurf(isppol))
 end do

end subroutine destroy_fermi_surface_1D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/nullify_gkk_0D
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the gkk_type structure, by nullifying them
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_gkk_0D(Gkk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_gkk_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_type),intent(inout) :: Gkk

! ************************************************************************

 !@gkk_type

 ! integer 
 nullify(Gkk%symgkk2symrec)
 nullify(Gkk%ibzgkk2fsbz)
 nullify(Gkk%fsbz2ibzgkk)
 nullify(Gkk%symtrbz2ibz)
 nullify(Gkk%symopbz2ibz)

 ! real
 nullify(Gkk%gkk)
 !nullify(Gkk%symrec_gkk)
 !nullify(Gkk%symafm_gkk)

end subroutine nullify_gkk_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/nullify_gkk_1D
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the gkk_type structure, by nullifying them
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_gkk_1D(Gkk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_gkk_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_type),intent(inout) :: Gkk(:)

!Local variables-------------------------------
!scalars
 integer :: ii
!arrays

! ************************************************************************

 do ii=1,SIZE(Gkk)
  call nullify_gkk_0D(Gkk(ii))
 end do

end subroutine nullify_gkk_1D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/destroy_gkk_0D
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the gkk_type structure
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_gkk_0D(Gkk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_gkk_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_type),intent(inout) :: Gkk

! ************************************************************************

 !@gkk_type
                            
 ! integer 
 if (associated(Gkk%symgkk2symrec))  then
   ABI_DEALLOCATE(Gkk%symgkk2symrec)
 end if
 if (associated(Gkk%ibzgkk2fsbz  ))  then
   ABI_DEALLOCATE(Gkk%ibzgkk2fsbz)
 end if
 if (associated(Gkk%fsbz2ibzgkk  ))  then
   ABI_DEALLOCATE(Gkk%fsbz2ibzgkk)
 end if
 if (associated(Gkk%symtrbz2ibz  ))  then
   ABI_DEALLOCATE(Gkk%symtrbz2ibz)
 end if
 if (associated(Gkk%symopbz2ibz  ))  then
   ABI_DEALLOCATE(Gkk%symopbz2ibz)
 end if
                            
 ! real
 if (associated(Gkk%gkk        ))   then
   ABI_DEALLOCATE(Gkk%gkk)
 end if
 !if (associated(Gkk%symrec_gkk)  deallocate(Gkk%symrec_gkk)
 !if (associated(Gkk%symafm_gkk)  deallocate(Gkk%symafm_gkk)

end subroutine destroy_gkk_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/destroy_gkk_1D 
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the gkk_type structure
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_gkk_1D(Gkk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_gkk_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_type),intent(inout) :: Gkk(:)

!Local variables-------------------------------
!scalars
 integer :: ii
!arrays

! ************************************************************************

 do ii=1,SIZE(Gkk)
  call destroy_gkk_0D(Gkk(ii))
 end do

end subroutine destroy_gkk_1D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/nullify_gkk_handler_0D
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the gkk_handler_type structure, by nullifying them
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_gkk_handler_0D(Gkk_hdl)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_gkk_handler_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_handler_type),intent(inout) :: Gkk_hdl


!Local variables-------------------------------
!scalars
 integer :: ii

! ************************************************************************

 !@gkk_handler_type

 ! integer 
 nullify(Gkk_hdl%pert_list)

 ! nested structures 
 ! NOTE: This coding works only if the array of substructures is already allocated.
 do ii=1,SIZE(Gkk_hdl%Gkk)
  call nullify_gkk(Gkk_hdl%Gkk(ii))
 end do

end subroutine nullify_gkk_handler_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/nullify_gkk_handler_1D
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the gkk_handler_type structure, by nullifying them
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_gkk_handler_1D(Gkk_hdl)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_gkk_handler_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_handler_type),intent(inout) :: Gkk_hdl(:)

!Local variables-------------------------------
!scalars
 integer :: ii
!arrays

! ************************************************************************

 do ii=1,SIZE(Gkk_hdl)
  call nullify_gkk_handler_0D(Gkk_hdl(ii))
 end do

end subroutine nullify_gkk_handler_1D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/destroy_gkk_handler_0D
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the gkk_handler_type structure
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_gkk_handler_0D(Gkk_hdl)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_gkk_handler_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_handler_type),intent(inout) :: Gkk_hdl

! ************************************************************************

 !@gkk_handler_type
                            
 ! integer 
 if (associated(Gkk_hdl%pert_list))  then
   ABI_DEALLOCATE(Gkk_hdl%pert_list)
 end if

 ! substructures
 call destroy_gkk(Gkk_hdl%Gkk)

end subroutine destroy_gkk_handler_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/destroy_gkk_handler_1D
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the gkk_handler_type structure
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_gkk_handler_1D(Gkk_hdl)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_gkk_handler_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_handler_type),intent(inout) :: Gkk_hdl(:)

!Local variables-------------------------------
!scalars
 integer :: ii
!arrays

! ************************************************************************

 do ii=1,SIZE(Gkk_hdl)
  call destroy_gkk_handler_0D(Gkk_hdl(ii))
 end do

end subroutine destroy_gkk_handler_1D
!!***

!----------------------------------------------------------------------

END MODULE m_eph
!!***
