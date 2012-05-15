!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_efield
!! NAME
!!  m_efield
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods 
!!  used to handle electric fields
!!  Imported object from defs_datatypes
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MJV)
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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_efield

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 implicit none

 private
!!***


!!****t* defs_datatypes/efield_type
!! NAME
!! efield_type
!!
!! FUNCTION
!! First-principles calculations in a finite electric field
!!
!! SOURCE

 type, public :: efield_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer variables
  integer :: berryopt            ! value of berryopt in use
  integer :: fmkmem              ! number of k-points in the FBZ per cpu
  integer :: fmkmem_max          ! max of fmkmem
  integer :: fnkpt               ! number of k-points in the FBZ
  integer :: has_expibi          ! 2 if expibi computed, 1 if only allocated, zero else
  integer :: has_expibr          ! 2 if expibr computed, 1 if only allocated, zero else
  integer :: has_Lij             ! 2 if paw_Lij computed, 1 if only allocated, zero else
  integer :: has_Lijr3           ! 2 if Lijr3 computed, 1 if only allocated, zero else
  integer :: has_rij             ! 2 if paw_rij computed, 1 if only allocated, zero else
  integer :: has_twdij0          ! 2 if twdij0 computed, 1 if only allocated, zero else
  integer :: has_tweijkl         ! 2 if tweijkl computed, 1 if only allocated, zero else
  integer :: has_qijb            ! 2 if paw_qijb computed, 1 if only allocated, zero else
  integer :: lmax
  integer :: lmnmax
  integer :: maxnstr             ! max number of strings along idir=1,2,3
  integer :: maxnkstr            ! max number of k-points per string
  integer :: mhcg                ! 2nd dimension of hcg array
  integer :: mkmem_max           ! max of mkmem
  integer :: natom               ! number of atoms in unit cell
  integer :: nband_occ           ! number of occupied bands
                                 ! this number must be the same for every k
  integer :: nsym
  integer :: usecprj             ! 1 if efield%cprj allocated (see below), 0 else
  integer :: usepaw              ! 1 if a PAW calculation, 0 else

! Integer arrays
  integer :: indhk(6,6)          ! index of phase twist terms for <u_k1|H_k2|u_k3> type structures.
                                 ! in this case there are 24 distinct terms.

  integer :: nstr(3)             ! nstr(idir) = number of strings along idir
  integer :: nkstr(3)            ! nkstr(idir) = number of k-points per string

  integer :: twind(6,6)          ! index of phase twist terms. First index gives bra
                                 ! phase in the form 2*bdir-bfor+1, where bdir is the
                                 ! direction of the bra phase vector (1,2,3) and bfor 
                                 ! indexes whether it is added to (1) or subtracted from (2)
                                 ! the base k vector. The second index gives the same information
                                 ! for the ket side. Equivalently, can use bsig=-1 for subtract,
                                 ! bsig=+1 for add, giving 2*bdir + (bsig+1)/2 - 1. 
                                 ! The value of twind is the index number
                                 ! of terms in qijb_bk and related structures, such that 
                                 ! |twind| is the term number, and if negative, take its
                                 ! complex comjugate. All this is necessary because only
                                 ! six combinations of bra and ket shift vectors are stored
                                 ! (1) -k1 -k2
                                 ! (2) +k1 -k2
                                 ! (3) -k2 -k3
                                 ! (4) +k2 -k3
                                 ! (5) -k3 -k1
                                 ! (6) +k3 -k1
                                 ! all others are obtained by appropriate complex conjugation

! Real(dp) scalars
  real(dp) :: sdeg               ! spin degeneracy: sdeg = 2 if nsppol = 1
                                 !                         1 if nsppol = 2

! Real(dp) arrays
  real(dp) :: bfield(3)          ! bfield vector in real space, atomic units
  real(dp) :: dkvecs(3,3)        ! dkvec(:,idir) = vector between a k-poinit
                                 ! and its nearest neighbour along idir
  real(dp) :: efield_dot(3)      ! reciprocal lattice coordinates of the
                                 ! electric field
  real(dp) :: gmet_str(2,2,3)    ! gmet_str(:,:,idir) is the metric of the metric of 
                                 ! the space of strings of direction idir
  real(dp) :: mag_cart(3)        ! magnetization in cartesian coordinates

! Integer pointers
  integer, pointer :: atom_indsym(:,:,:) ! atom_indsym(4,natom,nsym)
                                         ! this is data on how the symmetries map the atoms in the cell
                                         ! see symatm.F90 for full description
  integer, pointer :: cgindex(:,:)    ! cgindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cg array
  integer, pointer :: cgqindex(:,:,:) ! cgqindex(3,6,nkpt*nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cgq and pwnsfacq
                                      ! arrays
                                      ! (see vtorho.f and initberry.f)
  integer, pointer :: cprjindex(:,:)  ! cprjindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the cprj in the cprj array (used only
                                      ! for PAW calculations)
  integer, pointer :: fkgindex(:)     ! same as kgindex, but defined
                                      ! for the FBZ and intended to use
                                      ! with pwindf
  integer, pointer :: idxkstr(:,:,:)  ! idxkstr(maxnkstr,maxnstr,3)
                                      ! idxkstr(ikstr,istr,idir) index (ikpt) of
                                      ! k-point ikstr on string istr along idir
  integer, pointer :: ikpt_dk(:,:,:)  ! ikpt_dk(nkpt,2,3)
                                      ! ikpt_dp(ikpt,ii,idir) = index of the
                                      ! k-point at k+dk (ii=1) and k-dk (ii=2)
  integer, pointer :: indkk_f2ibz(:,:)   ! indkk_f2ibz(1:dtefield%fnkpt,1:6)
                                         ! information needed to fold a
                                         ! k-point in the FBZ into the IBZ;
                                         ! the second index (1:6)
                                         ! is as described in listkk
  integer, pointer :: i2fbz(:)           ! i2fbz(1:nkpt) gives index of IBZ
                                         ! k-points in the FBZ k-point list

  integer, pointer :: kgindex(:)      ! kgind(nkpt)
                                      ! kgind(ikpt) = ikg

  integer, pointer :: lmn_size(:)        ! lmn_size(ntypat)
  integer, pointer :: lmn2_size(:)       ! lmn2_size(ntypat)

  integer, pointer :: nneigh(:)          ! nneigh(nkpt)
                                         ! for each k-point, nneigh stores
                                         ! the number of its nearest neighbours
                                         ! that are not related by symmetry
  integer, pointer :: sflag(:,:,:,:)  ! sflag(nband_occ,nkpt*nsppol,2,3)
                                      ! sflag = 0 : compute the whole row of
                                      !             smat
                                      ! sflag = 1 : the row is up to date

  integer, pointer :: str_neigh(:,:,:)
  integer, pointer :: strg_neigh(:,:,:,:)
! str_neigh(ineigh, istr, idir) is the index ineigh-th neighbour of the istr-th string in
! the direction idir
! str_neigh(ineigh, istr, :, idir) is a 2-dimensional vector which coordinates are 0 or 1,
! useful only if the k-point mesh isn't a full mesh - if it's a single point, a line or a plane.


! Real(dp) pointers

  real(dp), pointer :: chern_k(:,:,:)
! chern_k(2,nkpt,dir)
! Chern-form at each k-point. Used for magnetization.

! the coordinates of the ineigh-th neighbour of the istr-th string in the direction idir are :
! coord_str(:,str_neigh(ineigh,istr,idir),idir) + real(str_neigh(ineigh, istr, :, idir),dp)
  real(dp),pointer :: coord_str(:,:,:)
! coord_str(1:2,istr,idir) are the coordinate of the istr-th string in the direction idir.

  real(dp),pointer :: emat(:,:,:)
! emat(2,nband_occ,nkpt*nsppol)
! Diagonal H matrix for every k-point. Used for magnetization.
! stores <u_nk|H_k|u_nk>

  real(dp), pointer :: expibi(:,:,:)
! expibi(2,natom,9) 
! used for PAW field calculations
! stores the on-site phase factors arising from
! $\langle\phi_{i,k+\sigma_b k_b}|\phi_{j,k+\sigma_k k_k}\rangle$ 
! where $\sigma = \pm 1$. These overlaps arise in various Berry
! phase calculations of electric and magnetic polarization. The on-site
! phase factor is $\exp[i(\sigma_b k_b - \sigma_k k_k)\cdot I]$ where
! $I$ is the nuclear position. Only the following
! are computed and saved, in the given order:
! 1) -k_1 - k_2
! 2) +k_1 - k_2
! 3) -k_2 - k_3
! 4) +k_2 - k_3
! 5) -k_3 - k_1
! 6) +k_3 - k_1
! 7)    0 - k_1
! 8)    0 - k_2
! 9)    0 - k_3

  real(dp), pointer :: expibr(:,:,:,:)
! expibr(2,natom,nfgd,6) 
! used for PAW field calculations
! stores the on-site phase factors arising from
! $\exp(i(\sigma_b k_b - \sigma_k k_k)\cdot r)$ 
! where $\sigma = \pm 1$, on the fine grid of points in the PAW sphere
! around each atom. These phases are needed to compute the twisted $\hat{D}_{ij}$
! term in the magnetic field calculations.
! Only the following are computed and saved, in the given order:
! 1) -k_1 - k_2
! 2) +k_1 - k_2
! 3) -k_2 - k_3
! 4) +k_2 - k_3
! 5) -k_3 - k_1
! 6) +k_3 - k_1

  real(dp), pointer :: fkptns(:,:)       ! fkptns(3,1:dtefield%fnkpt)
                                         ! k-points in FBZ

  real(dp), pointer :: hcg(:,:,:)
! hcg(2,mhcg,3)
! Script H wavefunctions at each k point, used only for finite field
! magnetization. Stored as
! hcg(2,npw_k*nband*nkpt,idir) like cg

! pointer to on-site angular momentum
  real(dp),pointer :: Lij(:,:,:,:) ! Lij(2,lmn2_size_max,ntypat,3)
! gives <(r-R)xp> at each atom type in each of 3 directions
! these are used only in the PAW case with magnetic field

! pointer to on-site angular momentum
  real(dp),pointer :: Lijr3(:,:,:,:) ! Lijr3(2,lmn2_size_max,ntypat,3)
! gives <(r-R)xp/|r-R|^3> at each atom type in each of 3 directions
! these are used only in the PAW case with magnetic field to get shielding

! pointer to magnetization
! this list gives the magnetization at each k point
! in each direction. Used for berryopt = 5 and -5.
  real(dp),pointer :: mag_k(:,:,:) ! mag_k(2,idir,ikpt)

! pointer to onsite local magnetization
! this list gives the on-site <L> part of the magnetization at each k point
! in each direction. Used for berryopt = 5 and -5.
  real(dp),pointer :: mag_local_k(:,:) ! mag_local_k(idir,ikpt)

  real(dp), pointer :: qijb_bk(:,:,:,:)
! qijb_bk(2,lmnmax,natom,6)
! on-site part of <u_nk+b1|u_mk+b2> matrix elements, relevant for PAW only
! only six combinations of bra and ket shift vectors are stored
! (the shifts are obtained from dkvecs above) in the following order:
! (1) -k1 -k2
! (2) +k1 -k2
! (3) -k2 -k3
! (4) +k2 -k3
! (5) -k3 -k1
! (6) +k3 -k1
! all others are obtained by appropriate complex conjugation

  real(dp), pointer :: qijb_kk(:,:,:,:)
! qijb_kk(2,lmnmax,natom,3)
! on-site part of <u_nk|u_mk+b> matrix elements, relevant for PAW only
! vector b described by idir (1,2,3), forward direction; value for
! reverse direction (ifor = 2 in berryphase_new and cgwf) obtained by
! complex conjugation

! pointer to on-site dipole moment
  real(dp),pointer :: rij(:,:,:) ! rij(lmn2_size_max,natom,3)
 ! gives <r-R> at each atom in each of 3 directions
 ! these are used only in the PAW case with electric field

  real(dp), pointer :: smat(:,:,:,:,:,:)
! smat(2,nband_occ,nband_occ,nkpt*nsppol,2,3)
! Overlap matrix for every k-point. In an electric field calculation,
! smat is updated at every iteration.

  real(dp), pointer :: twh(:,:,:,:,:)
! twh(2,nband_occ,nband_occ,nkpt*nsppol,tind)
! Overlap H matrix for every k-point. Used for magnetization.
! stores <u_nk_b|H_k|u_mk_k>
! tind is an index combining the bra (bdir, bfor) and ket
! (kdir, kfor) directions. Recall that bdir = 1,2,3 are the
! three directions in the cell in recip space, and for each one
! we have bfor = 1,2 for forward, backward. Likewise on the
! ket side for kdir and kfor. Then the ket elements are
! indexed as 2*kdir-kfor+1, and the bra elements as 
! 2*bdir-bfor+1. These are combined as 
! tind = 6*( (2*kdir-kfor+1) - 1) + 2*bdir-bfor+1

  real(dp), pointer :: twdij0(:,:,:,:,:)
! twdij0(2,lmnmax,lmnmax,natom,tind)
! for each atom, on-site Dij0 terms including phase shifts
! these include (Torrent CMS 42, 337 (2008) appendix E)
! (1) kinetic <phi_i|exp(I*b_b.r)(-del^2/2)exp(-I*b_k.r)|phi_j>
! (2b) Hartree n_ZC
! (2e) Hartree n_ZC charge compensation
! needed for magnetic fields
! tind is an index combining the bra (bdir, bfor) and ket
! (kdir, kfor) directions: see twind 

  real(dp), pointer :: twdij(:,:,:,:,:)
! twdij(2,lmnmax,lmnmax,natom,tind)
! for each atom, on-site Dij terms including phase shifts
! needed for magnetic fields
! tind is an index combining the bra (bdir, bfor) and ket
! (kdir, kfor) directions. Recall that bdir = 1,2,3 are the
! three directions in the cell in recip space, and for each one
! we have bfor = 1,2 for forward, backward. Likewise on the
! ket side for kdir and kfor. Then the ket elements are
! indexed as 2*kdir-kfor+1, and the bra elements as 
! 2*bdir-bfor+1. These are combined as 
! 6*( (2*kdir-kfor+1) - 1) + 2*bdir-bfor+1

  real(dp), pointer :: tweijkl(:,:,:,:,:)
! tweijkl(2,lmn2max,lmn2max,natom,6)
! for each atom, on-site e_ijkl terms including phase shifts
! these include (Torrent CMS 42, 337 (2008) appendix E)
! (2a*) on-site Hartree
! (2c*) \hat{n} Hartree
! (2d*) \tilde{n}^1 Hartree and charge compensation
! (2f*) \hat{n} Hartree and charge compensation
! needed for magnetic fields
! only six combinations of bra and ket shift vectors are stored
! (the shifts are obtained from dkvecs above) in the following order:
! (1) -k1 -k2
! (2) +k1 -k2
! (3) -k2 -k3
! (4) +k2 -k3
! (5) -k3 -k1
! (6) +k3 -k1
! all others are obtained by appropriate complex conjugation

  real(dp), pointer :: zarot(:,:,:,:)
   !  zarot(l_size_max,l_size_max,l_max,nsym)
   !  Coeffs of the transformation of real spherical
   !  harmonics under the symmetry operations. These are needed when the
   ! cprj's need to be computed in the full BZ, that is,
   ! in the PAW case with kptopt /= 3.

! pointer to cprj
   type(cprj_type),pointer :: cprj(:,:)
! used with finite efield and PAW


 end type efield_type
!!***

 ! Bound methods:
 public :: nullify_efield
 public :: destroy_efield

contains

!!****f* m_efield/nullify_efield
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the efield_type structure, by nullifying them
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_efield(dtefield)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_efield'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(efield_type),intent(out) :: dtefield

! ************************************************************************

! Integer pointers

  nullify(dtefield%atom_indsym) 
  nullify(dtefield%cgindex)    
  nullify(dtefield%cgqindex) 
  nullify(dtefield%cprjindex)
  nullify(dtefield%fkgindex)
  nullify(dtefield%idxkstr)
  nullify(dtefield%ikpt_dk)
  nullify(dtefield%indkk_f2ibz)
  nullify(dtefield%i2fbz)     
  nullify(dtefield%kgindex)  
  nullify(dtefield%lmn_size)
  nullify(dtefield%lmn2_size)
  nullify(dtefield%nneigh)  
  nullify(dtefield%sflag)  
  nullify(dtefield%str_neigh)
  nullify(dtefield%strg_neigh)

! Real(dp) pointers
  nullify(dtefield%chern_k)
  nullify(dtefield%coord_str)
  nullify(dtefield%emat)
  nullify(dtefield%expibi)
  nullify(dtefield%expibr)
  nullify(dtefield%fkptns) 
  nullify(dtefield%hcg)
  nullify(dtefield%Lij) 
  nullify(dtefield%Lijr3)
  nullify(dtefield%mag_k)
  nullify(dtefield%mag_local_k)
  nullify(dtefield%qijb_bk)
  nullify(dtefield%qijb_kk)
  nullify(dtefield%rij)
  nullify(dtefield%smat)
  nullify(dtefield%twh)
  nullify(dtefield%twdij0)
  nullify(dtefield%twdij)
  nullify(dtefield%tweijkl)
  nullify(dtefield%zarot)

! pointer to cprj
   nullify(dtefield%cprj)

end subroutine nullify_efield
!!***

!!****f* m_efield/destroy_efield
!! NAME
!!
!! FUNCTION
!!   deallocate fields in efield structure 
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_efield(dtefield)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_efield'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(efield_type),intent(out) :: dtefield

! ************************************************************************

! Integer pointers
  if(associated(dtefield%atom_indsym))  then
    ABI_DEALLOCATE(dtefield%atom_indsym)
  end if
  if(associated(dtefield%cgindex))  then
    ABI_DEALLOCATE(dtefield%cgindex)
  end if
  if(associated(dtefield%cgqindex))  then
    ABI_DEALLOCATE(dtefield%cgqindex)
  end if
  if(associated(dtefield%cprjindex))  then
    ABI_DEALLOCATE(dtefield%cprjindex)
  end if
  if(associated(dtefield%fkgindex))  then
    ABI_DEALLOCATE(dtefield%fkgindex)
  end if
  if(associated(dtefield%idxkstr))  then
    ABI_DEALLOCATE(dtefield%idxkstr)
  end if
  if(associated(dtefield%ikpt_dk))  then
    ABI_DEALLOCATE(dtefield%ikpt_dk)
  end if
  if(associated(dtefield%indkk_f2ibz))  then
    ABI_DEALLOCATE(dtefield%indkk_f2ibz)
  end if
  if(associated(dtefield%i2fbz))  then
    ABI_DEALLOCATE(dtefield%i2fbz)
  end if
  if(associated(dtefield%kgindex))  then
    ABI_DEALLOCATE(dtefield%kgindex)
  end if
  if(associated(dtefield%lmn_size))  then
    ABI_DEALLOCATE(dtefield%lmn_size)
  end if
  if(associated(dtefield%lmn2_size))  then
    ABI_DEALLOCATE(dtefield%lmn2_size)
  end if
  if(associated(dtefield%nneigh))  then
    ABI_DEALLOCATE(dtefield%nneigh)
  end if
  if(associated(dtefield%sflag))  then
    ABI_DEALLOCATE(dtefield%sflag)
  end if
  if(associated(dtefield%str_neigh))  then
    ABI_DEALLOCATE(dtefield%str_neigh)
  end if
  if(associated(dtefield%strg_neigh))  then
    ABI_DEALLOCATE(dtefield%strg_neigh)
  end if

! Real(dp) pointers

  if(associated(dtefield%chern_k))  then
    ABI_DEALLOCATE(dtefield%chern_k)
  end if
  if(associated(dtefield%coord_str))  then
    ABI_DEALLOCATE(dtefield%coord_str)
  end if
  if(associated(dtefield%emat))  then
    ABI_DEALLOCATE(dtefield%emat)
  end if
  if(associated(dtefield%expibi))  then
    ABI_DEALLOCATE(dtefield%expibi)
  end if
  if(associated(dtefield%expibr))  then
    ABI_DEALLOCATE(dtefield%expibr)
  end if
  if(associated(dtefield%fkptns))  then
    ABI_DEALLOCATE(dtefield%fkptns)
  end if
  if(associated(dtefield%hcg))  then
    ABI_DEALLOCATE(dtefield%hcg)
  end if
  if(associated(dtefield%Lij))  then
    ABI_DEALLOCATE(dtefield%Lij)
  end if
  if(associated(dtefield%Lijr3))  then
    ABI_DEALLOCATE(dtefield%Lijr3)
  end if
  if(associated(dtefield%mag_k))  then
    ABI_DEALLOCATE(dtefield%mag_k)
  end if
  if(associated(dtefield%mag_local_k))  then
    ABI_DEALLOCATE(dtefield%mag_local_k)
  end if
  if(associated(dtefield%qijb_bk))  then
    ABI_DEALLOCATE(dtefield%qijb_bk)
  end if
  if(associated(dtefield%qijb_kk))  then
    ABI_DEALLOCATE(dtefield%qijb_kk)
  end if
  if(associated(dtefield%rij))  then
    ABI_DEALLOCATE(dtefield%rij)
  end if
  if(associated(dtefield%smat))  then
    ABI_DEALLOCATE(dtefield%smat)
  end if
  if(associated(dtefield%twh))  then
    ABI_DEALLOCATE(dtefield%twh)
  end if
  if(associated(dtefield%twdij0))  then
    ABI_DEALLOCATE(dtefield%twdij0)
  end if
  if(associated(dtefield%twdij))  then
    ABI_DEALLOCATE(dtefield%twdij)
  end if
  if(associated(dtefield%tweijkl))  then
    ABI_DEALLOCATE(dtefield%tweijkl)
  end if
  if(associated(dtefield%zarot))  then
    ABI_DEALLOCATE(dtefield%zarot)
  end if

! pointer to cprj
  if(associated(dtefield%cprj)) then
    call cprj_free(dtefield%cprj)
    ABI_DEALLOCATE(dtefield%cprj)
  end if

end subroutine destroy_efield
!!***

end module m_efield
!!***
