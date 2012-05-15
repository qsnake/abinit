!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_toolbox
!! NAME
!!  m_paw_toolbox
!!
!! FUNCTION
!!  This module contains basic tools to initialize, nullify and free basic
!!  PAW objects.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!! This module should contain the definition of the data types as well.
!! At present these procedures are used only in GW but, hopefully, they will
!! become standard PAW methods.
!!
!! NOTES
!!
!! * Routines tagged with "@type_name" are strongly connected to the definition of the data type.
!!   Strongly connected means that the proper functioning of the implementation relies on the
!!   assumption that the tagged procedure is consistent with the type declaration.
!!   Every time a developer changes the structure "type_name" adding new entries, he/she has to make sure
!!   that all the strongly connected routines are changed accordingly to accommodate the modification of the data type.
!!   Typical examples of strongly connected routines are creation, destruction or reset methods.
!!
!! * gcc 4.3 does not hanlde correctly the case in which there are routines with similar names defined
!!   in the same module. This creates problems at link-time if the module contains procedures
!!   beginning with the same name, e.g init_longname, init_longname_slightly_different.
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

MODULE m_paw_toolbox

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_fstrings,    only : toupper
 use m_crystal,     only : crystal_structure
 use defs_abitypes, only : dataset_type

!for PAW partial waves
 use m_radmesh,  only : deducer0
 use m_splines,  only : spline, splint

 implicit none

 private

 public :: reset_paw_ij_flags
 public :: destroy_paw_ij
 public :: init_paw_ij
 public :: destroy_paw_an
 public :: destroy_pawtab
 public :: init_pawfgr
 public :: nullify_paw_ij
 public :: nullify_paw_an
 public :: nullify_pawtab
 public :: init_paw_an
 public :: print_pawtab
 public :: print_paw_ij     !FIXME not yet tested
 public :: pawfgrtab_free
 public :: pawfgrtab_init
 public :: pawfgrtab_print
 public :: pawfgrtab_nullify
 public :: reset_paw_an_flags
 public :: init_pawang
 public :: destroy_pawang
 public :: get_dimcprj !  Helper function returning the number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom.
!!***

!----------------------------------------------------------------------

!!****t* m_paw_toolbox/paw_pwaves_lmn_t
!! NAME
!! paw_pwaves_lmn_t
!!
!! FUNCTION
!!  Datatype used to store the 3D values of the all-electron and of the pseudized part of the
!!  PAW partial waves on the set of FFT points falling inside the spheres around each atom.
!!  The data is mainly used for plotting the true PAW wavefunctions in real space.
!!
!! SOURCE

 type,public ::  paw_pwaves_lmn_t

  integer :: nfgd

  integer :: lmn_size

  !$integer :: ngfft(18)

  integer,pointer :: r0shift(:,:)  SET2NULL
  ! r0shift(3,nfgd)

  !real(dp),pointer :: phk_atm(:,:)  SET2NULL
  ! phk_atmt(2,nfgd)

  real(dp),pointer :: phi(:,:)    SET2NULL
  ! phi (nfgd,lmn_size)
  ! \phi_{nlm}(ifgd) for each point of the FFT mesh located inside the PAW sphere (see pawfgrtab_type).

  real(dp),pointer :: tphi(:,:)   SET2NULL
  ! tphi (nfgd,lmn_size)
  ! \tphi_{nlm}(ifgd) for each point of the FFT mesh located inside the PAW sphere (see pawfgrtab_type).

  real(dp),pointer :: phi_gr(:,:,:)    SET2NULL
  ! phi_gr (3,nfgd,lmn_size)
  ! gradient, in cartesian coordinates, of \phi_{nlm}(ifgd) for each point of the FFT mesh 
  ! located inside the PAW sphere (see pawfgrtab_type).

  real(dp),pointer :: tphi_gr(:,:,:)   SET2NULL
  ! tphi_gr (3,nfgd,lmn_size)
  ! gradient, in cartesian coordinates, of \tphi_{nlm}(ifgd) for each point of the FFT mesh 
  ! located inside the PAW sphere (see pawfgrtab_type).

 end type paw_pwaves_lmn_t
!!***

public :: init_paw_pwaves_lmn
public :: destroy_paw_pwaves_lmn
public :: nullify_paw_pwaves_lmn

!----------------------------------------------------------------------


CONTAINS  !===========================================================
!!***

!!****f* m_paw_toolbox/reset_paw_ij_flags
!! NAME
!! reset_paw_ij_flags
!!
!! FUNCTION
!!  Set all paw_ij flags set to 0.
!!
!! SIDE EFFECTS
!!  Paw_ij_flags<type(paw_ij_flags_type)>=flags in a paw_ij structure
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine reset_paw_ij_flags(Paw_ij_flags)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'reset_paw_ij_flags'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_ij_flags_type),intent(inout) :: Paw_ij_flags

! *************************************************************************

 ! @Paw_ij_flags_type
 Paw_ij_flags%has_dij       =0
 Paw_ij_flags%has_dijfr     =0
 Paw_ij_flags%has_dijhartree=0
 Paw_ij_flags%has_dijhat    =0
 Paw_ij_flags%has_dijso     =0
 Paw_ij_flags%has_dijU      =0
 Paw_ij_flags%has_dijxc     =0
 Paw_ij_flags%has_dijxc_val =0
 Paw_ij_flags%has_exexch_pot=0
 Paw_ij_flags%has_pawu_occ  =0

end subroutine reset_paw_ij_flags
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/copy_paw_ij_flags
!! NAME
!! copy_paw_ij_flags
!!
!! FUNCTION
!!  Copy a Paw_ij_flags structure.
!!
!! INPUTS
!!  ij_flags_in<type(paw_ij_flags_type)>=input a paw_ij_flags structure
!!
!! OUTPUT
!!  ij_flags_out<type(paw_ij_flags_type)>=output a paw_ij_flags structure
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_paw_ij_flags(ij_flags_in, ij_flags_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_paw_ij_flags'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_ij_flags_type),intent(in)    :: ij_flags_in
 type(Paw_ij_flags_type),intent(inout) :: ij_flags_out

! *************************************************************************

 ! @Paw_ij_flags_type
 ij_flags_out%has_dij        = ij_flags_in%has_dij
 ij_flags_out%has_dijfr      = ij_flags_in%has_dijfr
 ij_flags_out%has_dijhartree = ij_flags_in%has_dijhartree
 ij_flags_out%has_dijhat     = ij_flags_in%has_dijhat
 ij_flags_out%has_dijso      = ij_flags_in%has_dijso
 ij_flags_out%has_dijU       = ij_flags_in%has_dijU
 ij_flags_out%has_dijxc      = ij_flags_in%has_dijxc
 ij_flags_out%has_dijxc_val  = ij_flags_in%has_dijxc_val
 ij_flags_out%has_exexch_pot = ij_flags_in%has_exexch_pot
 ij_flags_out%has_pawu_occ  = ij_flags_in%has_pawu_occ

end subroutine copy_paw_ij_flags
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/destroy_paw_ij
!! NAME
!!  destroy_paw_ij
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a paw_ij structure
!!
!! SIDE EFFECTS
!!  paw_ij(:)<type(paw_ij_type)>=paw arrays given on (i,j) channels
!!
!! PARENTS
!!      bethe_salpeter,dyfnl3,ldau_self,m_energy,nstpaw3,respfn,scfcv,scfcv3
!!      screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_paw_ij(Paw_ij)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_paw_ij'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_ij_type),intent(inout) :: Paw_ij(:)

!Local variables-------------------------------
 integer :: iat,natom

! *************************************************************************

 DBG_ENTER("COLL")

 ! @Paw_ij_type
 natom=SIZE(Paw_ij)
 do iat=1,natom
  if (associated(Paw_ij(iat)%dij       ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%dij)
  end if
  if (associated(Paw_ij(iat)%dijfr     ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%dijfr)
  end if
  if (associated(Paw_ij(iat)%dijhartree))  then
    ABI_DEALLOCATE(Paw_ij(iat)%dijhartree)
  end if
  if (associated(Paw_ij(iat)%dijhat    ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%dijhat)
  end if
  if (associated(Paw_ij(iat)%dijU      ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%dijU)
  end if
  if (associated(Paw_ij(iat)%dijso     ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%dijso)
  end if
  if (associated(Paw_ij(iat)%dijxc     ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%dijxc)
  end if
  if (associated(Paw_ij(iat)%dijxc_val ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%dijxc_val)
  end if
  if (associated(Paw_ij(iat)%noccmmp   ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%noccmmp)
  end if
  if (associated(Paw_ij(iat)%nocctot   ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%nocctot)
  end if
  if (associated(Paw_ij(iat)%vpawx     ))  then
    ABI_DEALLOCATE(Paw_ij(iat)%vpawx)
  end if

  ! === Reset all has_* flags ===
  Paw_ij(iat)%has_dij       =0
  Paw_ij(iat)%has_dijfr     =0
  Paw_ij(iat)%has_dijhartree=0
  Paw_ij(iat)%has_dijhat    =0
  Paw_ij(iat)%has_dijso     =0
  Paw_ij(iat)%has_dijU      =0
  Paw_ij(iat)%has_dijxc     =0
  Paw_ij(iat)%has_dijxc_val =0
  Paw_ij(iat)%has_exexch_pot=0
  Paw_ij(iat)%has_pawu_occ  =0
 end do

 DBG_EXIT("COLL")

end subroutine destroy_Paw_ij
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/init_paw_ij
!! NAME
!! init_paw_ij
!!
!! FUNCTION
!!  Initialize a Paw_ij data type.
!!
!! INPUTS
!!  cplex=1 if all on-site PAW quantities are real, 2 if they are complex
!!  cplex_dij=1 if dij are real, 2 if they are complex
!!  natom=Number of atoms.
!!  ntypat=Number of types of atoms in cell.
!!  nspinor=number of spinor components
!!  nsppol=Number of independent spin polarizations.
!!  nspden=Number of spin-density components
!!  pawspnorb=1 if spin-orbit coupling is activated
!!  typat(natom)=Type of each atom
!!  Pawtab(ntypat)<type(pawtab_type)>=PAW tabulated starting data
!!
!! OPTIONAL INPUTS
!!  has_dij=1 to allocate Paw_ij%dij, 0 otherwise (default)
!!  has_dijhat=1 to allocate Paw_ij%dijhat, 0 otherwise (default)
!!  has_dijxc=1 to allocate Paw_ij%dijxc, 0 otherwise (default)
!!  has_dijxc_val=1 to allocate Paw_ij%dijxc_val, 0 otherwise (default)
!!  has_dijhartree=1 to allocate Paw_ij%dijhartree, 0 otherwise (default)
!!  has_dijso=1 to allocate Paw_ij%dijso, used only if pawspnorb>0. 0 otherwise (default)
!!  has_dijU=1 to allocate Paw_ij%dijU, used only if Pawtab(itypat)%usepawu>0. 0 otherwise (default).
!!  has_exexch_pot=1 to allocate potnetial used in PAW+(local exact exchange) formalism, 0 otherwise (default)
!!  has_pawu_occ=1 to allocate occupations used in PAW+U formalism, 0 otherwise (default)
!!
!! OUTPUT
!!  Paw_ij(natom)<type(paw_ij_type)>=data structure containing PAW arrays given on (i,j) channels.
!!   In output all the basic dimensions are defined and the arrays are allocated
!!   according to the input variables.
!!
!! PARENTS
!!      bethe_salpeter,dyfnl3,ldau_self,m_energy,nstpaw3,paw_qpscgw,respfn
!!      scfcv,scfcv3,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_paw_ij(Paw_ij,cplex,cplex_dij,nspinor,nsppol,nspden,pawspnorb,natom,ntypat,typat,Pawtab,&
&                      has_dij,has_dijfr,has_dijhartree,has_dijhat,& ! Optional
&                      has_dijxc,has_dijxc_val,has_dijso,has_dijU,&  ! Optional
&                      has_exexch_pot,has_pawu_occ)                  ! Optional

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_paw_ij'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_dij,nspinor,nspden,nsppol,natom,ntypat,pawspnorb
 integer,optional,intent(in) :: has_dij,has_dijfr,has_dijhat,has_dijxc,has_dijxc_val
 integer,optional,intent(in) :: has_dijso,has_dijhartree,has_dijU,has_exexch_pot,has_pawu_occ
!arrays
 integer,intent(in) :: typat(natom)
 type(Paw_ij_type),intent(inout) :: Paw_ij(natom)
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)

!Local variables-------------------------------
 integer :: iat,itypat,lmn2_size,ndij

! *************************************************************************

 DBG_ENTER("COLL")

 !allocate(Paw_ij(natom))
 !call nullify_paw_ij(Paw_ij)

 ! @Paw_ij_type
 do iat=1,natom
  itypat=typat(iat)
  lmn2_size              =Pawtab(itypat)%lmn2_size
  Paw_ij(iat)%cplex      =cplex
  !Paw_ij(iat)%cplex_dij  =nspinor
  !Paw_ij(iat)%cplex_dij  =MAX(cplex,1+pawspnorb,nspinor)
  Paw_ij(iat)%cplex_dij  =cplex_dij
  Paw_ij(iat)%nspden     =nspden
  Paw_ij(iat)%nsppol     =nsppol
  Paw_ij(iat)%lmn_size   =Pawtab(itypat)%lmn_size
  Paw_ij(iat)%lmn2_size  =lmn2_size
  Paw_ij(iat)%ndij       =MAX(nspinor**2,nspden)
  !Paw_ij(iat)%lmnmix_sz =  do we need this? It seems it is not used anymore and can be removed

  ndij=Paw_ij(iat)%ndij

  ! ==================================
  ! === Allocations (all optional) ===
  ! ==================================

  ! === Allocation for total Dij ===
  Paw_ij(iat)%has_dij=0
  if (PRESENT(has_dij)) then
    if (has_dij/=0) then
      Paw_ij(iat)%has_dij=1
      ABI_ALLOCATE(Paw_ij(iat)%dij,(cplex_dij*lmn2_size,ndij))
    end if
  end if

  ! === Allocation for total Dij_Hartree ===
  Paw_ij(iat)%has_dijhartree=0
  if (PRESENT(has_dijhartree)) then
    if (has_dijhartree/=0) then
      Paw_ij(iat)%has_dijhartree=1
      ABI_ALLOCATE(Paw_ij(iat)%dijhartree,(cplex*lmn2_size))
    end if
  end if

  ! === Allocation for total Dij_hat ===
  Paw_ij(iat)%has_dijhat=0
  if (PRESENT(has_dijhat)) then
    if (has_dijhat/=0) then
      Paw_ij(iat)%has_dijhat=1
      ABI_ALLOCATE(Paw_ij(iat)%dijhat,(cplex_dij*lmn2_size,ndij))
    end if
  end if

  ! === Allocation for total Dij_XC ===
  Paw_ij(iat)%has_dijxc=0
  if (PRESENT(has_dijxc)) then
    if (has_dijxc/=0) then
      Paw_ij(iat)%has_dijxc=1
      ABI_ALLOCATE(Paw_ij(iat)%dijxc,(cplex_dij*lmn2_size,ndij))
    end if
  end if

  ! === Allocation for total Dij_XC_val ===
  Paw_ij(iat)%has_dijxc_val=0
  if (PRESENT(has_dijxc_val)) then
    if (has_dijxc_val/=0) then
      Paw_ij(iat)%has_dijxc_val=1
      ABI_ALLOCATE(Paw_ij(iat)%dijxc_val,(cplex_dij*lmn2_size,ndij))
    end if
  end if

  ! === Allocation for total Dij_U_val ===
  Paw_ij(iat)%has_dijU=0
  if (PRESENT(has_dijU)) then
    if (has_dijU/=0.and.Pawtab(itypat)%usepawu>0) then
      Paw_ij(iat)%has_dijU=1
      ABI_ALLOCATE(Paw_ij(iat)%dijU,(cplex_dij*lmn2_size,ndij))
    end if
  end if

  ! === Allocation for total Dij_SO ===
  Paw_ij(iat)%has_dijso=0
  if (PRESENT(has_dijso)) then
    if (has_dijso/=0.and.pawspnorb>0) then
      Paw_ij(iat)%has_dijso=1
      ABI_ALLOCATE(Paw_ij(iat)%dijso,(cplex_dij*lmn2_size,ndij))
     end if
  end if

  ! === Allocation for frozen part of 1st-order Dij ===
  Paw_ij(iat)%has_dijfr=0
  if (PRESENT(has_dijfr)) then
    if (has_dijfr/=0) then
      Paw_ij(iat)%has_dijfr=1
      ABI_ALLOCATE(Paw_ij(iat)%dijfr,(cplex_dij*lmn2_size,ndij))
    end if
  end if

  ! === Allocation for PAW+U occupations ===
  Paw_ij(iat)%has_pawu_occ=0
  if (PRESENT(has_pawu_occ)) then
    if (has_pawu_occ/=0.and.Pawtab(itypat)%usepawu>0) then
      Paw_ij(iat)%has_pawu_occ=1
      ABI_ALLOCATE(Paw_ij(iat)%noccmmp,(cplex_dij,2*Pawtab(itypat)%lpawu+1,2*Pawtab(itypat)%lpawu+1,ndij))
      ABI_ALLOCATE(Paw_ij(iat)%nocctot,(ndij))
    end if
  end if

  ! === Allocation for PAW+LEXX potential ===
  Paw_ij(iat)%has_exexch_pot=0
  if (PRESENT(has_exexch_pot)) then
    if (has_exexch_pot/=0.and.Pawtab(itypat)%useexexch>0) then
      Paw_ij(iat)%has_exexch_pot=1
    ! TODO solve issue with first dimension
      ABI_ALLOCATE(Paw_ij(iat)%vpawx,(1,lmn2_size,nspden))
     end if
  end if

 end do !iat

 DBG_EXIT("COLL")

end subroutine init_paw_ij
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/destroy_paw_an
!! NAME
!! destroy_paw_an
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a paw_an structure
!!
!! SIDE EFFECTS
!!  Paw_an(:)<type(Paw_an_type)>=various arrays given on ANgular mesh or ANgular moments
!!
!! SIDE EFFECTS
!!  All associated pointers in Paw_an(:) are deallocated
!!
!! NOTES
!!  vh1 and vht1 are defined in the data structure but never used.
!!  Cannot test for association status since these quantities are
!!  not nullified before entering the calculation
!!
!! PARENTS
!!      bethe_salpeter,respfn,scfcv,scfcv3,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_Paw_an(Paw_an)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_Paw_an'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,natom
!integer :: itypat

! *************************************************************************

 DBG_ENTER("COLL")

 !@Paw_an_type

 natom=SIZE(Paw_an)

! do iat=1,natom
!  itypat=typat(iat)
!  deallocate(Paw_an(iat)%lmselect)
!  !deallocate(Paw_an(iat)%vh1,Paw_an(iat)%vht1)      !TODO nullify these arrays
!  deallocate(paw_an(iat)%vxc1,Paw_an(iat)%vxct1)
!  if (Paw_an(iat)%has_vxcval==1 ) deallocate(Paw_an(iat)%vxc1_val,Paw_an(iat)%vxct1_val)
!  if (Pawtab(itypat)%useexexch>0) deallocate(Paw_an(iat)%vxc_ex)
! end do

 do iat=1,natom
  if (associated(Paw_an(iat)%lmselect ))  then
    ABI_DEALLOCATE(Paw_an(iat)%lmselect)
  end if
  if (associated(Paw_an(iat)%vh1      ))  then
    ABI_DEALLOCATE(Paw_an(iat)%vh1)
  end if
  if (associated(Paw_an(iat)%vht1     ))  then
    ABI_DEALLOCATE(Paw_an(iat)%vht1)
  end if
  if (associated(Paw_an(iat)%vxc1     ))  then
    ABI_DEALLOCATE(Paw_an(iat)%vxc1)
  end if
  if (associated(Paw_an(iat)%vxc1_val ))  then
    ABI_DEALLOCATE(Paw_an(iat)%vxc1_val)
  end if
  if (associated(Paw_an(iat)%vxct1    ))  then
    ABI_DEALLOCATE(Paw_an(iat)%vxct1)
  end if
  if (associated(Paw_an(iat)%vxct1_val))  then
    ABI_DEALLOCATE(Paw_an(iat)%vxct1_val)
  end if
  if (associated(Paw_an(iat)%vxc_ex   ))  then
    ABI_DEALLOCATE(Paw_an(iat)%vxc_ex)
  end if
  if (associated(Paw_an(iat)%kxc1     ))  then
    ABI_DEALLOCATE(Paw_an(iat)%kxc1)
  end if
  if (associated(Paw_an(iat)%kxct1    ))  then
    ABI_DEALLOCATE(Paw_an(iat)%kxct1)
  end if

  if (associated(Paw_an(iat)%vxc_ex   ))  then
    ABI_DEALLOCATE(Paw_an(iat)%vxc_ex)
  end if

  ! === Reset all has_* flags ===
  Paw_an(iat)%has_kxc     =0
  Paw_an(iat)%has_vhartree=0
  Paw_an(iat)%has_vxc     =0
  Paw_an(iat)%has_vxcval  =0
 end do !iat


 DBG_EXIT("COLL")

end subroutine destroy_Paw_an
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/destroy_pawtab
!! NAME
!!  destroy_pawtab
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a pawtab structure
!!
!! SIDE EFFECTS
!!  Pawtab(:)<type(pawtab_type)>=PAW arrays tabulated.
!!  All associated pointers in Pawtab(:) are deallocated
!!
!! PARENTS
!!      mblktyp1,mblktyp5,rdddb9,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_pawtab(Pawtab)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_pawtab'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Pawtab_type),intent(inout) :: Pawtab(:)

!Local variables-------------------------------
 integer :: ityp,ntypat

! *************************************************************************

 !@Pawtab_type
 ntypat=SIZE(Pawtab(:))
 if (ntypat==0) return

 do ityp=1,ntypat

  if (associated(Pawtab(ityp)%indklmn))  then
    ABI_DEALLOCATE(Pawtab(ityp)%indklmn)
  end if
  if (associated(Pawtab(ityp)%klmntomn))  then
    ABI_DEALLOCATE(Pawtab(ityp)%klmntomn)
  end if
  if (associated(Pawtab(ityp)%kmix))  then
    ABI_DEALLOCATE(Pawtab(ityp)%kmix)
  end if
  if (associated(Pawtab(ityp)%lnproju))  then
    ABI_DEALLOCATE(Pawtab(ityp)%lnproju)
  end if
  if (associated(Pawtab(ityp)%coredens))  then
    ABI_DEALLOCATE(Pawtab(ityp)%coredens)
  end if
  if (associated(Pawtab(ityp)%dij0))  then
    ABI_DEALLOCATE(Pawtab(ityp)%dij0)
  end if
  if (associated(Pawtab(ityp)%dltij))  then
    ABI_DEALLOCATE(Pawtab(ityp)%dltij)
  end if
  if (associated(Pawtab(ityp)%dshpfunc))  then
    ABI_DEALLOCATE(Pawtab(ityp)%dshpfunc)
  end if
  if (associated(Pawtab(ityp)%eijkl))  then
    ABI_DEALLOCATE(Pawtab(ityp)%eijkl)
  end if
  if (associated(Pawtab(ityp)%fk))  then
    ABI_DEALLOCATE(Pawtab(ityp)%fk)
  end if
  if (associated(Pawtab(ityp)%gnorm))  then
    ABI_DEALLOCATE(Pawtab(ityp)%gnorm)
  end if
  if (associated(Pawtab(ityp)%kij))  then
    ABI_DEALLOCATE(Pawtab(ityp)%kij)
  end if
  if (associated(Pawtab(ityp)%nabla_ij))  then
    ABI_DEALLOCATE(Pawtab(ityp)%nabla_ij)
  end if
  if (associated(Pawtab(ityp)%phi))  then
    ABI_DEALLOCATE(Pawtab(ityp)%phi)
  end if
  if (associated(Pawtab(ityp)%phiphj))  then
    ABI_DEALLOCATE(Pawtab(ityp)%phiphj)
  end if
  if (associated(Pawtab(ityp)%phiphjint))  then
    ABI_DEALLOCATE(Pawtab(ityp)%phiphjint)
  end if
  if (associated(Pawtab(ityp)%ph0phiint))  then
    ABI_DEALLOCATE(Pawtab(ityp)%ph0phiint)
  end if
  if (associated(Pawtab(ityp)%qgrid_shp))  then
    ABI_DEALLOCATE(Pawtab(ityp)%qgrid_shp)
  end if
  if (associated(Pawtab(ityp)%qijl))  then
    ABI_DEALLOCATE(Pawtab(ityp)%qijl)
  end if
  if (associated(Pawtab(ityp)%rad_for_spline))  then
    ABI_DEALLOCATE(Pawtab(ityp)%rad_for_spline)
  end if
  if (associated(Pawtab(ityp)%rhoij0))  then
    ABI_DEALLOCATE(Pawtab(ityp)%rhoij0)
  end if
  if (associated(Pawtab(ityp)%shape_alpha))  then
    ABI_DEALLOCATE(Pawtab(ityp)%shape_alpha)
  end if
  if (associated(Pawtab(ityp)%shape_q))  then
    ABI_DEALLOCATE(Pawtab(ityp)%shape_q)
  end if
  if (associated(Pawtab(ityp)%shapefunc))  then
    ABI_DEALLOCATE(Pawtab(ityp)%shapefunc)
  end if
  if (associated(Pawtab(ityp)%shapefncg))  then
    ABI_DEALLOCATE(Pawtab(ityp)%shapefncg)
  end if
  if (associated(Pawtab(ityp)%sij))  then
    ABI_DEALLOCATE(Pawtab(ityp)%sij)
  end if
  if (associated(Pawtab(ityp)%tcoredens))  then
    ABI_DEALLOCATE(Pawtab(ityp)%tcoredens)
  end if
  if (associated(Pawtab(ityp)%tcorespl))  then
    ABI_DEALLOCATE(Pawtab(ityp)%tcorespl)
  end if
  if (associated(Pawtab(ityp)%tphi))  then
    ABI_DEALLOCATE(Pawtab(ityp)%tphi)
  end if
  if (associated(Pawtab(ityp)%tphitphj))  then
    ABI_DEALLOCATE(Pawtab(ityp)%tphitphj)
  end if
  if (associated(Pawtab(ityp)%tvalespl))  then
    ABI_DEALLOCATE(Pawtab(ityp)%tvalespl)
  end if
  if (associated(Pawtab(ityp)%Vee))  then
    ABI_DEALLOCATE(Pawtab(ityp)%Vee)
  end if
  if (associated(Pawtab(ityp)%Vex))  then
    ABI_DEALLOCATE(Pawtab(ityp)%Vex)
  end if
  if (associated(Pawtab(ityp)%VHntZC))  then
    ABI_DEALLOCATE(Pawtab(ityp)%VHntZC)
  end if
  if (associated(Pawtab(ityp)%VHnZC))  then
    ABI_DEALLOCATE(Pawtab(ityp)%VHnZC)
  end if
  if (associated(Pawtab(ityp)%zioneff))  then
    ABI_DEALLOCATE(Pawtab(ityp)%zioneff)
  end if

  ! === Reset all has_* flags ===
  Pawtab(ityp)%has_kij  =0
  Pawtab(ityp)%has_nabla=0
  Pawtab(ityp)%has_vhntzc = 0
  Pawtab(ityp)%has_vhnzc = 0
  Pawtab(ityp)%usetcore =0
  Pawtab(ityp)%usetvale =0
 end do

end subroutine destroy_pawtab
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/init_pawfgr
!! NAME
!! init_pawfgr
!!
!! FUNCTION
!!  Initialize a pawfgr_type datatype, reporting also info on the FFT mesh
!!  according to the method used (norm-conserving or PAW)
!!
!! INPUTS
!!  k0(3)=input k vector for k+G sphere
!!  Dtset <type(dataset_type)>=all input variables for this dataset
!!   %dilatmx
!!   %usepaw
!!   %natom
!!   %ngfft
!!   %ngfftdg
!!   %nfft
!!   %mgfft
!!   %mgfftdg
!!   %dilatmx
!!   %pawecutdg
!!   %ecut
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!
!! OUTPUT
!!  ecut_eff=effective energy cutoff (hartree) for coarse planewave basis sphere
!!  ecutdg_eff=effective energy cutoff (hartree) for dense planewave basis sphere
!!  gsqcutc_eff=(PAW) Fourier cutoff on G^2 for "large sphere" of radius double for the coarse FFT grid
!   gsqcutf_eff=Fourier cutoff on G^2 for "large sphere" of radius double for the dense FFT grid
!!  nfftf=(effective) number of FFT grid points (for this proc), for dense FFT mesh
!!  mgfftf=maximum size of 1D FFTs, for dense FFT mesh
!!  ngfftc(18),ngfftf(18)=contain all needed information about 3D FFT, for coarse and dense FFT mesh, resp.
!!                        see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  Pawfgr<pawfgr_type>=For PAW, Fine rectangular GRid parameters and related data
!!
!! PARENTS
!!      bethe_salpeter,gstate,respfn,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_pawfgr(Dtset,Pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
&                      gsqcutc_eff,gsqcutf_eff,gmet,k0) ! optional
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_pawfgr'
 use interfaces_14_hidewrite
 use interfaces_53_ffts
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: nfftf,mgfftf
 real(dp),intent(out) :: ecut_eff,ecutdg_eff
 real(dp),intent(out),optional :: gsqcutf_eff,gsqcutc_eff
 type(dataset_type),intent(in) :: Dtset
 type(Pawfgr_type),intent(out) :: Pawfgr
!arrays
 real(dp),intent(in),optional :: gmet(3,3)
 integer,intent(out) :: ngfftc(18),ngfftf(18)
 real(dp),intent(in),optional :: k0(3)

!Local variables-------------------------------
 integer :: ii,nfftc_tot,nfftf_tot
 real(dp) :: boxcut,boxcutc
 character(len=500) :: msg

!************************************************************************

 DBG_ENTER("COLL")

 !@Pawfgr_type

 if ((present(gsqcutc_eff).or.present(gsqcutf_eff)).and.&
&    ((.not.present(gmet)).or.(.not.present(k0)))) then
   msg='To compute gsqcut[c,f]_eff, both k0 and gmet must be present as argument !'
   MSG_BUG(msg)
 end if

 ngfftc(:)=Dtset%ngfft(:)

 SELECT CASE (Dtset%usepaw)

 CASE (0)
  ! === Norm-conserving pseudopotentials ===
  nfftf=Dtset%nfft ; mgfftf=Dtset%mgfft ; ngfftf(:)=Dtset%ngfft(:)
  Pawfgr%usefinegrid=0 
  ABI_ALLOCATE(Pawfgr%coatofin,(0))
  ABI_ALLOCATE(Pawfgr%fintocoa,(0))
  ecut_eff  =Dtset%ecut*Dtset%dilatmx**2
  ecutdg_eff=ecut_eff

 CASE (1)
  ! == PAW calculation ===
  if (Dtset%pawecutdg>=1.0000001_dp*Dtset%ecut) then
   ! * Use fine FFT grid generated according to pawecutdg.
   nfftf=Dtset%nfftdg ; mgfftf=Dtset%mgfftdg ; ngfftf(:)=Dtset%ngfftdg(:)
   nfftc_tot =ngfftc(1)*ngfftc(2)*ngfftc(3)
   nfftf_tot =ngfftf(1)*ngfftf(2)*ngfftf(3)
   Pawfgr%usefinegrid=1 
   ABI_ALLOCATE(Pawfgr%coatofin,(nfftc_tot))
   ABI_ALLOCATE(Pawfgr%fintocoa,(nfftf_tot))
   call indgrid(Pawfgr%coatofin,Pawfgr%fintocoa,nfftc_tot,nfftf_tot,ngfftc,ngfftf)
  else
   ! * Do not use fine FFT mesh. Simple transfer that can be done in parallel with only local info.
   nfftf=Dtset%nfft ; mgfftf=Dtset%mgfft ; ngfftf(:)=Dtset%ngfft(:)
   Pawfgr%usefinegrid=0 
   ABI_ALLOCATE(Pawfgr%coatofin,(Dtset%nfft))
   ABI_ALLOCATE(Pawfgr%fintocoa,(Dtset%nfft))
   do ii=1,Dtset%nfft
    Pawfgr%coatofin(ii)=ii ; Pawfgr%fintocoa(ii)=ii
   end do
  end if

  ! == Store useful dimensions in Pawfgr ===
  Pawfgr%nfftc=Dtset%nfft ; Pawfgr%mgfftc=Dtset%mgfft ; Pawfgr%ngfftc(:)=Dtset%ngfft(:)
  Pawfgr%nfft=nfftf       ; Pawfgr%mgfft=mgfftf       ; Pawfgr%ngfft (:)=ngfftf(:)
  ecutdg_eff=Dtset%pawecutdg*Dtset%dilatmx**2
  ecut_eff  =Dtset%ecut*Dtset%dilatmx**2

 CASE DEFAULT
  write(msg,'(a,i4)')' Wrong value of usepaw: ',Dtset%usepaw
  MSG_BUG(msg)
 END SELECT
 !
 ! === Get boxcut for given gmet, ngfft, and ecut (center at k0) ===
 !     boxcut=ratio of basis sphere diameter to fft box side
 boxcut=-one
 if (Dtset%usepaw==1) then
   if (present(gsqcutc_eff)) then
     write(msg,'(2a)')ch10,' Coarse grid specifications '!(used for wave-functions):'
     call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL')
     call getcut(boxcutc,ecut_eff,gmet,gsqcutc_eff,Dtset%iboxcut,std_out,k0,ngfftc)
   end if
   if (present(gsqcutf_eff)) then
     write(msg,'(2a)')ch10,' Fine grid specifications (used for densities):'
     call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL')
     call getcut(boxcut,ecutdg_eff,gmet,gsqcutf_eff,Dtset%iboxcut,std_out,k0,ngfftf)
   end if
 else if (present(gsqcutc_eff)) then
   call getcut(boxcut,ecut_eff,gmet,gsqcutc_eff,Dtset%iboxcut,std_out,k0,ngfftc)
   gsqcutf_eff=gsqcutc_eff
 end if
 !
 ! === Check that boxcut>=2 if intxc=1; otherwise intxc must be set=0 ===
 if (boxcut>=zero .and. boxcut<two .and. Dtset%intxc==1) then
   write(msg,'(a,es12.4,5a)')&
&   ' boxcut=',boxcut,' is < 2.0  => intxc must be 0;',ch10,&
&   ' Need larger ngfft to use intxc=1.',ch10,&
&   ' Action : you could increase ngfft, or decrease ecut, or put intxc=0.'
   MSG_ERROR(msg)
 end if

 DBG_EXIT("COLL")

end subroutine init_pawfgr
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/nullify_paw_ij
!! NAME
!!  nullify_paw_ij
!!
!! FUNCTION
!!  Nullify pointers and flags in a paw_ij structure
!!
!! SIDE EFFECTS
!!  Paw_ij(:)<type(paw_ij_type)>=PAW arrays given on (i,j) channels. Nullified in output
!!
!! PARENTS
!!      bethe_salpeter,dyfnl3,ldau_self,m_energy,nstpaw3,paw_qpscgw,respfn
!!      scfcv,scfcv3,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_paw_ij(Paw_ij)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_paw_ij'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_ij_type),intent(inout) :: Paw_ij(:)

!Local variables-------------------------------
 integer :: iat,natom

! *************************************************************************

 !@Paw_ij_type

 natom=SIZE(Paw_ij(:))

 do iat=1,natom
  nullify(Paw_ij(iat)%dij       )
  nullify(Paw_ij(iat)%dijfr     )
  nullify(Paw_ij(iat)%dijhartree)
  nullify(Paw_ij(iat)%dijhat    )
  nullify(Paw_ij(iat)%dijU      )
  nullify(Paw_ij(iat)%dijso     )
  nullify(Paw_ij(iat)%dijxc     )
  nullify(Paw_ij(iat)%dijxc_val )
  nullify(Paw_ij(iat)%noccmmp   )
  nullify(Paw_ij(iat)%nocctot   )
  nullify(Paw_ij(iat)%vpawx     )

  ! === Set all has_* flags to zero ===
  Paw_ij(iat)%has_dij       =0
  Paw_ij(iat)%has_dijfr     =0
  Paw_ij(iat)%has_dijhartree=0
  Paw_ij(iat)%has_dijhat    =0
  Paw_ij(iat)%has_dijso     =0
  Paw_ij(iat)%has_dijU      =0
  Paw_ij(iat)%has_dijxc     =0
  Paw_ij(iat)%has_dijxc_val =0
  Paw_ij(iat)%has_exexch_pot=0
  Paw_ij(iat)%has_pawu_occ  =0
 end do !iat

end subroutine nullify_paw_ij
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/nullify_paw_an
!! NAME
!!  nullify_paw_an
!!
!! FUNCTION
!!  Nullify pointers and flags in a paw_an structure
!!
!! SIDE EFFECTS
!!  Paw_an(:)<type(paw_an_type)>=PAW arrays given on ANgular mesh or ANgular moments.
!!                               Nullified in output
!!
!! PARENTS
!!      bethe_salpeter,paw_qpscgw,respfn,scfcv,scfcv3,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_paw_an(Paw_an)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_paw_an'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,natom

! *************************************************************************

 !@Paw_an_type
 natom=SIZE(Paw_an(:))

 do iat=1,natom
  nullify(Paw_an(iat)%lmselect )
  nullify(Paw_an(iat)%vh1      )
  nullify(Paw_an(iat)%vht1     )
  nullify(Paw_an(iat)%vxc1     )
  nullify(Paw_an(iat)%vxct1    )
  nullify(Paw_an(iat)%vxc1_val )
  nullify(Paw_an(iat)%vxct1_val)
  nullify(Paw_an(iat)%vxc_ex   )
  nullify(Paw_an(iat)%kxc1     )
  nullify(Paw_an(iat)%kxct1    )

  ! === Set all has_* flags to zero ===
  Paw_an(iat)%has_kxc      =0
  Paw_an(iat)%has_vhartree =0
  Paw_an(iat)%has_vxc      =0
  Paw_an(iat)%has_vxcval   =0
 end do

end subroutine nullify_paw_an
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/nullify_pawtab
!! NAME
!!  nullify_pawtab
!!
!! FUNCTION
!!  Nullify pointers and flags in a pawtab structure
!!
!! SIDE EFFECTS
!!  Pawtab(:)<type(pawtab_type)>=PAW arrays tabulated.
!!                               Nullified in output
!!
!! PARENTS
!!      mblktyp1,mblktyp5,rdddb9,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_pawtab(Pawtab)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_pawtab'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Pawtab_type),intent(inout) :: Pawtab(:)

!Local variables-------------------------------
 integer :: ityp,ntypat

! *************************************************************************

 !@Pawtab_type
 ntypat=SIZE(Pawtab(:))
 if (ntypat==0) return

 do ityp=1,ntypat

  nullify(Pawtab(ityp)%indklmn)
  nullify(Pawtab(ityp)%klmntomn)
  nullify(Pawtab(ityp)%kmix)
  nullify(Pawtab(ityp)%lnproju)
  nullify(Pawtab(ityp)%coredens)
  nullify(Pawtab(ityp)%dij0)
  nullify(Pawtab(ityp)%dltij)
  nullify(Pawtab(ityp)%dshpfunc)
  nullify(Pawtab(ityp)%eijkl)
  nullify(Pawtab(ityp)%fk)
  nullify(Pawtab(ityp)%gnorm)
  nullify(Pawtab(ityp)%kij)
  nullify(Pawtab(ityp)%nabla_ij)
  nullify(Pawtab(ityp)%phi)
  nullify(Pawtab(ityp)%phiphj)
  nullify(Pawtab(ityp)%phiphjint)
  nullify(Pawtab(ityp)%ph0phiint)
  nullify(Pawtab(ityp)%qgrid_shp)
  nullify(Pawtab(ityp)%qijl)
  nullify(Pawtab(ityp)%rad_for_spline)
  nullify(Pawtab(ityp)%rhoij0)
  nullify(Pawtab(ityp)%shape_alpha)
  nullify(Pawtab(ityp)%shape_q)
  nullify(Pawtab(ityp)%shapefunc)
  nullify(Pawtab(ityp)%shapefncg)
  nullify(Pawtab(ityp)%sij)
  nullify(Pawtab(ityp)%tcoredens)
  nullify(Pawtab(ityp)%tcorespl)
  nullify(Pawtab(ityp)%tphi)
  nullify(Pawtab(ityp)%tphitphj)
  nullify(Pawtab(ityp)%tvalespl)
  nullify(Pawtab(ityp)%Vee)
  nullify(Pawtab(ityp)%Vex)
  nullify(Pawtab(ityp)%VHntZC)
  nullify(Pawtab(ityp)%VHnZC)
  nullify(Pawtab(ityp)%zioneff)

  ! === Set all has_* flags to zero ===
  Pawtab(ityp)%has_kij  =0
  Pawtab(ityp)%has_nabla=0
  Pawtab(ityp)%has_vhntzc = 0
  Pawtab(ityp)%has_vhnzc = 0
  Pawtab(ityp)%usetcore =0
  Pawtab(ityp)%usetvale =0
 end do

end subroutine nullify_pawtab
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/init_paw_an
!! NAME
!!  init_paw_an
!!
!! FUNCTION
!!  Initialize a paw_an data type.
!!
!! SIDE EFFECTS
!!  Paw_an(:)<type(paw_an_type)>=PAW arrays given on ANgular mesh or ANgular moments.
!!                               Initialized in output
!!
!! PARENTS
!!      bethe_salpeter,paw_qpscgw,respfn,scfcv,scfcv3,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE


subroutine init_paw_an(natom,ntypat,nkxc1,nspden,cplex,pawxcdev,typat,Pawang,Pawtab,Paw_an,&
&                      has_vhartree,has_vxc,has_vxcval,has_kxc) ! Optional

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_paw_an'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nkxc1,ntypat,cplex,nspden,pawxcdev
 integer,optional,intent(in) :: has_vhartree,has_vxc,has_vxcval,has_kxc
!arrays
 integer,intent(in) :: typat(natom)
 type(Pawang_type),intent(in) :: Pawang
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)
 type(Paw_an_type),intent(inout) :: Paw_an(:)

!Local variables-------------------------------
 integer :: iat,itypat,lm_size,v_size

! *************************************************************************

 DBG_ENTER("COLL")

 !@Paw_an_type

 !allocate(Paw_an(natom))
 !call nullify_paw_an(Paw_an)

 do iat=1,natom
  itypat=typat(iat)
  lm_size                =Pawtab(itypat)%lcut_size**2
  Paw_an(iat)%angl_size  =Pawang%angl_size
  Paw_an(iat)%cplex      =cplex
  Paw_an(iat)%lm_size    =lm_size
  Paw_an(iat)%mesh_size  =Pawtab(itypat)%mesh_size
  Paw_an(iat)%nkxc1      =nkxc1
  Paw_an(iat)%nspden     =nspden

  ! === Non-zero LM-moments of "one-center" densities/potentials ===
  ! * Filled in pawdenpot.
  ABI_ALLOCATE(Paw_an(iat)%lmselect,(lm_size))

  v_size=Paw_an(iat)%lm_size ; if (pawxcdev==0) v_size=Paw_an(iat)%angl_size

 ! === XC potential inside the sphere ===
 ! * LM-moments of potential if pawxcdev/=0
 ! * (theta,phi) values of potential if pawxcdev=0
  Paw_an(iat)%has_vxc=0
  if (PRESENT(has_vxc)) then
   if (has_vxc>0) then
    Paw_an(iat)%has_vxc=1
    ABI_ALLOCATE(Paw_an(iat)%vxc1 ,(cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
    ABI_ALLOCATE(Paw_an(iat)%vxct1,(cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
   end if
  end if

  ! ==========================
  ! === Optional arguments ===
  ! ==========================

  ! * XC potential inside PAW spheres generated by valence electrons.
  Paw_an(iat)%has_vxcval=0
  if (PRESENT(has_vxcval)) then
   if (has_vxcval>0) then
    Paw_an(iat)%has_vxcval=1
    ABI_ALLOCATE(Paw_an(iat)%vxc1_val ,(cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
    ABI_ALLOCATE(Paw_an(iat)%vxct1_val,(cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
   end if
  end if

  ! * XC potential for local exact exchange inside the sphere.
  if (Pawtab(itypat)%useexexch>0) then
   ABI_ALLOCATE(Paw_an(iat)%vxc_ex,(cplex*Pawtab(itypat)%mesh_size,v_size,nspden))
  end if

  ! * Hartree potential LM-moments inside the sphere.
  Paw_an(iat)%has_vhartree=0
  if (PRESENT(has_vhartree)) then
   if (has_vhartree>0) then
    Paw_an(iat)%has_vhartree=1
! FIXME what about vht1?
!MG This is the coding PRESENTLY used in pawdenpot but the commented code should be the correct one
!MT: don't agreee for nspden (there is no dependance of vH^(1) with nspden)
    ABI_ALLOCATE(Paw_an(iat)%vh1,(cplex*Pawtab(itypat)%mesh_size,1,1))
    !$allocate(Paw_an(iat)%vh1 (cplex*Pawtab(itypat)%mesh_size,lm_size,nspden))
    !$allocate(Paw_an(iat)%vht1(cplex*Pawtab(itypat)%mesh_size,lm_size,nspden))
   end if
  end if

  ! xc kernels inside the sphere.
  Paw_an(iat)%has_kxc=0
  if (PRESENT(has_kxc)) then
   if (has_kxc>0) then
    Paw_an(iat)%has_kxc=1
    ABI_ALLOCATE(Paw_an(iat)%kxc1 ,(cplex*Pawtab(itypat)%mesh_size,v_size,nkxc1))
    ABI_ALLOCATE(Paw_an(iat)%kxct1,(cplex*Pawtab(itypat)%mesh_size,v_size,nkxc1))
   end if
  end if

 end do !iat

 DBG_EXIT("COLL")

end subroutine init_paw_an
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/print_pawtab
!! NAME
!! print_pawtab
!!
!! FUNCTION
!!  Print out the content of a pawtab datastructure
!!
!! INPUTS
!!  Pawtab<pawtab_type> Only for PAW, TABulated data initialized at start
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_pawtab(Pawtab,header,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_pawtab'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
!arrays
 type(Pawtab_type) :: Pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: ityp,ntypat,my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 write(msg,'(6a)')&
&  ' ==================================== ',ch10,&
&  ' ==== Info on PAW TABulated data ==== ',ch10,&
&  ' ==================================== ',ch10
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 ntypat=SIZE(Pawtab(:))

 do ityp=1,ntypat

  ! Print out integer values (dimensions)
  write(msg,'(a)')'                                 '
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a)')'  ****************************** '
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4,a)')'  **** Atom type ',ityp,' ****   '
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a)')'  ****************************** '
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Number of (n,l) elements ....................... ',Pawtab(ityp)%basis_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Number of (l,m,n) elements ..................... ',Pawtab(ityp)%lmn_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Number of (i,j) elements (packed form) ......... ',Pawtab(ityp)%ij_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Max L+1 leading to non-zero Gaunt .............. ',Pawtab(ityp)%l_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Max L+1 leading to non-zero Gaunt (pawlcutd) ... ',Pawtab(ityp)%lcut_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  lmn2_size ...................................... ',Pawtab(ityp)%lmn2_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  lmnmix_sz ...................................... ',Pawtab(ityp)%lmnmix_sz
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Size of radial mesh ............................ ',Pawtab(ityp)%mesh_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  No of Q-points for tcorespl and tvalespl ....... ',Pawtab(ityp)%mqgrid
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  No of Q-points for the radial shape functions .. ',Pawtab(ityp)%mqgrid_shp
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Radial shape function type ..................... ',Pawtab(ityp)%shape_type
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  shape_lambda ................................... ',Pawtab(ityp)%shape_lambda
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Use pseudized core density ..................... ',Pawtab(ityp)%usetcore
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Use pseudized valence density .................. ',Pawtab(ityp)%usetvale
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Option for the use of hat density in XC terms .. ',Pawtab(ityp)%usexcnhat
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Use LDA+U ...................................... ',Pawtab(ityp)%usepawu
  call wrtout(ab_out,msg,'COLL')
  if (Pawtab(ityp)%usepawu/=0) then
    write(msg,'(a,i4)')'  L on which U is applied ........................ ',Pawtab(ityp)%lpawu
    call wrtout(ab_out,msg,'COLL')
  end if
  write(msg,'(a,i4)')'  Use Local Exact exchange ....................... ',Pawtab(ityp)%useexexch
  call wrtout(ab_out,msg,'COLL')
  if (Pawtab(ityp)%useexexch/=0) then
    write(msg,'(a,i4)')'  L on which local exact-exchange is applied ..... ',Pawtab(ityp)%lexexch
    call wrtout(ab_out,msg,'COLL')
  end if
  if (Pawtab(ityp)%usepawu/=0.or.Pawtab(ityp)%useexexch/=0) then
    write(msg,'(a,i4)')'  Number of (i,j) elements for PAW+U or EXX ..... ',Pawtab(ityp)%ij_proj
    call wrtout(ab_out,msg,'COLL')
    write(msg,'(a,i4)')'  Number of projectors on which U or EXX acts .... ',Pawtab(ityp)%nproju
    call wrtout(ab_out,msg,'COLL')
  end if

  ! "Has" flags
  write(msg,'(a,i4)')'  Has kij   ...................................... ',Pawtab(ityp)%has_kij
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Has nabla ...................................... ',Pawtab(ityp)%has_nabla
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Has vhntzc ..................................... ',Pawtab(ityp)%has_vhntzc
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Has vhnzc ...................................... ',Pawtab(ityp)%has_vhnzc
  call wrtout(ab_out,msg,'COLL')
  !
  ! Real scalars
  write(msg,'(a,es16.8)')'  1/q d(tNcore(q))/dq for q=0 .....................',Pawtab(ityp)%dncdq0
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  1/q d(tNvale(q))/dq for q=0 .....................',Pawtab(ityp)%dnvdq0
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  XC energy for the core density ..................',Pawtab(ityp)%exccore
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  Mixing of exact exchange (PBE0) .................',Pawtab(ityp)%exchmix
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  Radius of the PAW sphere ........................',Pawtab(ityp)%rpaw
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  Compensation charge radius (if >rshp, g(r)=0) ...',Pawtab(ityp)%rshp !(if r>rshp, g(r)=zero)
  call wrtout(ab_out,msg,'COLL')
  if (Pawtab(ityp)%shape_type==2) then
   write(msg,'(a,es16.8)')'  Sigma parameter in gaussian shape function ......',Pawtab(ityp)%shape_sigma !(shape_type=2)
   call wrtout(ab_out,msg,'COLL')
  end if
  if (Pawtab(ityp)%usepawu/=0) then
   write(msg,'(a,es16.8)')'  Value of the U parameter [eV] ...................',Pawtab(ityp)%upawu*Ha_eV
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,es16.8)')'  Value of the J parameter [eV] ...................',Pawtab(ityp)%jpawu*Ha_eV
   call wrtout(ab_out,msg,'COLL')
  end if

 end do ! ityp
 !
 ! The other (huge) arrays are not reported..

end subroutine print_pawtab
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/print_paw_ij
!! NAME
!! print_paw_ij
!!
!! FUNCTION
!!  Print out the content of a paw_ij datastructure
!!
!! INPUTS
!! [unit]=the unit number for output
!! [pawprtvol]=verbosity level
!! [mode_paral]=either "COLL" or "PERS"
!!
!! OUTPUT
!! (Only writing)
!!
!! NOTES
!!  The implementation of the routine is not yet completed.
!!  This implementation does not work if Paw_ij contains a 1st-order Dij (DFPT) at a non-zero q
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_paw_ij(Paw_ij,unit,pawprtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_paw_ij'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: pawprtvol
 integer,optional,intent(in) :: unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 type(Paw_ij_type),intent(in) :: Paw_ij(:)

!Local variables-------------------------------
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
!scalars
 integer :: cplex,cplex_dij,iatom,idij,lmn2_size,lmn_size,natom,nspden,nsploop,nsppol,my_unt
 integer :: opt_sym,tmp_cplex_dij,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer,allocatable :: idum(:)
 real(dp),pointer :: dij2p(:)

! *************************************************************************

 DBG_ENTER("COLL")

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(pawprtvol )) my_prtvol=pawprtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 ! TODO: to be consistent, my_unt and my_mode should be passed to print_ij
 ! moreover the pointers should be nullified when paw_ij is initialized

 natom  = SIZE(Paw_ij)
 nsppol = Paw_ij(1)%nsppol
 nspden = Paw_ij(1)%nspden
 nsploop= nsppol; if (Paw_ij(1)%ndij==4) nsploop=4

 do iatom=1,natom

  lmn_size  = Paw_ij(iatom)%lmn_size
  lmn2_size = Paw_ij(iatom)%lmn2_size
  cplex_dij = Paw_ij(iatom)%cplex_dij
  cplex     = Paw_ij(iatom)%cplex

  ! ====================================
  ! === Loop over density components ===
  ! ====================================
  do idij=1,nsploop

   ! * Print title.
   if (ABS(my_prtvol)>=1) then
    if (iatom==1.or.iatom==natom.or.my_prtvol<0) then
     if (nspden==2.and.nsppol==1) then
      write(msg,'(2a,i3,2a)')ch10,&
&      ' >>>>>>>>>> Atom ',iatom,':',ch10,&
&      ' (antiferromagnetism case: only one spin component)'
     else
      write(msg,'(2a,i3,3a)') ch10,&
&      ' >>>>>>>>>> Atom ',iatom,' (component ',TRIM(dspin(idij+2*(nsploop/4))),'):'
     end if
     call wrtout(my_unt,msg,my_mode)
    end if
   end if

   !if (abs(my_prtvol)>=1) then
   ! if (iatom==1.or.iatom==natom.or.my_prtvol<0) then
   !  write(msg, '(a)') '   ************ Dij atomic (Dij0) ***********'
   !  call wrtout(my_unt,msg,my_mode)
   !  call print_ij(Pawtab(itypat)%dij0,lmn2_size,1,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1)
   ! end if
   !end if

   if (abs(my_prtvol)>=1) then
    if (Paw_ij(iatom)%has_dijhartree==2) then
     if (iatom==1.or.iatom==natom.or.my_prtvol<0) then
      write(msg, '(a)')'   ************** Dij Hartree ***************'
      call wrtout(my_unt,msg,my_mode)
      call print_ij(Paw_ij(iatom)%dijhartree,lmn2_size,cplex,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1)
     end if
    end if
   end if

   if (Paw_ij(iatom)%has_dijxc>0) then
    if ((abs(my_prtvol)>=1).and.(idij<=2.or.nspden==4)) then
     if (iatom==1.or.iatom==natom.or.my_prtvol<0) then
      write(msg,'(a)')'   ****************** Dij_xc + Dijhat_xc ****************'
      call wrtout(my_unt,msg,my_mode)
      if (idij<=nsppol.or.idij==2) then
       opt_sym=2; tmp_cplex_dij=1
       dij2p => Paw_ij(iatom)%dijxc(1:cplex_dij*lmn2_size:cplex_dij,idij)
      else
       opt_sym=1; tmp_cplex_dij=cplex_dij
       dij2p => Paw_ij(iatom)%dijxc(1:cplex_dij*lmn2_size:1,idij)
      end if
      call print_ij(dij2p,lmn2_size,tmp_cplex_dij,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1,opt_sym=opt_sym)
     end if
    end if
   end if

   if (Paw_ij(iatom)%has_dijxc_val>0) then
    if ((abs(my_prtvol)>=1).and.(idij<=2.or.nspden==4)) then
     if (iatom==1.or.iatom==natom.or.my_prtvol<0) then
      write(msg,'(a)')'   ****************** Dij_xcval ****************'
      call wrtout(my_unt,msg,my_mode)
      if (idij<=nsppol.or.idij==2) then
       opt_sym=2; tmp_cplex_dij=1
       dij2p => Paw_ij(iatom)%dijxc_val(1:cplex_dij*lmn2_size:cplex_dij,idij)
      else
       opt_sym=1; tmp_cplex_dij=cplex_dij
       dij2p => Paw_ij(iatom)%dijxc(1:cplex_dij*lmn2_size:1,idij)
      end if
      call print_ij(dij2p,lmn2_size,tmp_cplex_dij,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1,opt_sym=opt_sym)
     end if
    end if
   end if

   if (Paw_ij(iatom)%has_dijhat>0) then
    if ((abs(my_prtvol)>=1).and.(idij<=2.or.nspden==4)) then
     if (iatom==1.or.iatom==natom.or.my_prtvol<0) then
      write(msg,'(a)')'   ************* Dij_hat (Veff_ij) **********'
      call wrtout(my_unt,msg,my_mode)
      !if ((idij<=nsppol.or.idij==2).and.cplex==1)then
      if ((idij<=nsppol.or.idij==2))then
       opt_sym=2; tmp_cplex_dij=1
       dij2p => Paw_ij(iatom)%dijhat(1:cplex_dij*lmn2_size:cplex_dij,idij)
       !call print_ij(dijhat(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1)
      else
        opt_sym=1; tmp_cplex_dij=cplex_dij
        dij2p => Paw_ij(iatom)%dijxc(1:cplex_dij*lmn2_size:1,idij)
       !call print_ij(dijhat,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1,opt_sym=1)
      end if
      call print_ij(dij2p,lmn2_size,tmp_cplex_dij,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1,opt_sym=opt_sym)
     end if
    end if
   end if

   if (Paw_ij(iatom)%has_dijso>0) then
    if (abs(my_prtvol)>=1) then
     if (iatom==1.or.iatom==natom.or.my_prtvol<0) then
      write(msg,'(a)')'   ************** Dij SpinOrbit ************'
      call wrtout(my_unt,msg,my_mode)
      dij2p =>  Paw_ij(iatom)%dijso(:,idij)
      call print_ij(dij2p,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1,opt_sym=3)
     end if
    end if
   end if

   if (Paw_ij(iatom)%has_dijU>0) then
    if ((abs(my_prtvol)>=1).and.(idij<=2.or.nspden==4)) then
     if (iatom==1.or.iatom==natom.or.my_prtvol<0) then
      write(msg,'(a)')'   ************* Dij_LDA+U (dijpawu) **********'
      call wrtout(my_unt,msg,my_mode)
      if (idij<=nsppol.or.idij==2) then
       opt_sym=2; tmp_cplex_dij=1
       dij2p => Paw_ij(iatom)%dijU(1:cplex_dij*lmn_size:cplex_dij,idij)
       !call print_ij(dijpawu(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1)
      else
       opt_sym=1; tmp_cplex_dij=cplex_dij
       dij2p => Paw_ij(iatom)%dijU(1:cplex_dij*lmn2_size:1,idij)
       !call print_ij(dijpawu,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1,opt_sym=1)
      end if
      call print_ij(dij2p,lmn2_size,tmp_cplex_dij,lmn_size,1,-1,idum,0,my_prtvol,idum,-1.d0,1,opt_sym=opt_sym)
     end if
    end if
   end if

   !TODO Dij_Exact_Exchange is not printed because there is no entry in the objects
   ! Add new entries in Paw_ij

  end do !idij
 end do !iat

 write(msg,'(a)')' '
 call wrtout(my_unt,msg,my_mode)

 DBG_ENTER("COLL")

end subroutine print_paw_ij
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/pawfgrtab_free
!! NAME
!! pawfgrtab_free
!!
!! FUNCTION
!!  Free all dynamic memory stored in a pawfgrtab datastructure
!!
!! SIDE EFFECTS
!!  Pawfgrt(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!
!! PARENTS
!!      bethe_salpeter,classify_bands,denfgr,exc_plot,m_wfs,pawmkaewf,respfn
!!      scfcv,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_free(Pawfgrt)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawfgrtab_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Pawfgrtab_type),intent(inout) :: Pawfgrt(:)

!Local variables-------------------------------
!scalars
 integer :: iat,natom

! *************************************************************************

 DBG_ENTER("COLL")

 !@Pawfgrtab_type
 natom=SIZE(Pawfgrt)
 do iat=1,natom
  if (associated(Pawfgrt(iat)%ifftsph))  then
    ABI_DEALLOCATE(Pawfgrt(iat)%ifftsph)
  end if
  if (associated(Pawfgrt(iat)%gylm   ))  then
    ABI_DEALLOCATE(Pawfgrt(iat)%gylm)
  end if
  if (associated(Pawfgrt(iat)%gylmgr ))  then
    ABI_DEALLOCATE(Pawfgrt(iat)%gylmgr)
  end if
  if (associated(Pawfgrt(iat)%gylmgr2))  then
    ABI_DEALLOCATE(Pawfgrt(iat)%gylmgr2)
  end if
  if (associated(Pawfgrt(iat)%nhatfr ))  then
    ABI_DEALLOCATE(Pawfgrt(iat)%nhatfr)
  end if
  if (associated(Pawfgrt(iat)%rfgd   ))  then
    ABI_DEALLOCATE(Pawfgrt(iat)%rfgd)
  end if
  if (associated(Pawfgrt(iat)%expiqr))   then
    ABI_DEALLOCATE(Pawfgrt(iat)%expiqr)
  end if
 end do

 DBG_EXIT("COLL")

end subroutine pawfgrtab_free
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/pawfgrtab_init
!! NAME
!! pawfgrtab_init
!!
!! FUNCTION
!!  Initialize a pawfgrtab datastructure
!!
!! OUTPUT
!!  Pawfgrt(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!
!! PARENTS
!!      bethe_salpeter,classify_bands,denfgr,exc_plot,m_wfs,pawmkaewf,respfn
!!      scfcv,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_init(Pawfgrt,cplex,l_size_atm,nspden)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawfgrtab_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden
!arrays
 integer,intent(in) :: l_size_atm(:)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrt(:)

!Local variables-------------------------------
!scalars
 integer :: iat,natom
 character(len=500) :: msg

! *************************************************************************

 natom=SIZE(Pawfgrt)
 if (natom/=SIZE(l_size_atm)) then
  msg='Sizes of assumed shape arrays do not match'
  MSG_BUG(msg)
 end if

 !@Pawfgrtab_type
 call pawfgrtab_nullify(Pawfgrt)

 do iat=1,natom
  Pawfgrt(iat)%cplex             = cplex
  Pawfgrt(iat)%nspden            = nspden
  Pawfgrt(iat)%l_size            = l_size_atm(iat)
  Pawfgrt(iat)%nfgd              = 0  
  ABI_ALLOCATE(Pawfgrt(iat)%ifftsph,(0))
  Pawfgrt(iat)%gylm_allocated    = 0  
  ABI_ALLOCATE(Pawfgrt(iat)%gylm,(0,0))
  Pawfgrt(iat)%gylmgr_allocated  = 0  
  ABI_ALLOCATE(Pawfgrt(iat)%gylmgr,(0,0,0))
  Pawfgrt(iat)%gylmgr2_allocated = 0  
  ABI_ALLOCATE(Pawfgrt(iat)%gylmgr2,(0,0,0))
  Pawfgrt(iat)%nhatfr_allocated  = 0  
  ABI_ALLOCATE(Pawfgrt(iat)%nhatfr,(0,0))
  Pawfgrt(iat)%rfgd_allocated    = 0  
  ABI_ALLOCATE(Pawfgrt(iat)%rfgd,(0,0))
  Pawfgrt(iat)%expiqr_allocated  = 0  
  ABI_ALLOCATE(Pawfgrt(iat)%expiqr,(0,0))
 end do

end subroutine pawfgrtab_init
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/pawfgrtab_print
!! NAME
!! pawfgrtab_print
!!
!! FUNCTION
!!  Reports basic info on the pawfgrtab datatype.
!!
!! INPUTS
!! Pawfgrt<pawfgrtab_type>=The datatype to be printed
!! [mode_paral]=either "COLL" or "PERS", "COLL" if None.
!! [unit]=Unit number for output, std_out if None.
!! [prtvol]=Verbosity level, lowest if None.
!!
!! OUTPUT
!! (only writing)
!!
!! PARENTS
!!      exc_plot,m_wfs,pawmkaewf,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_print(Pawfgrt,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawfgrtab_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 type(Pawfgrtab_type),intent(inout) :: Pawfgrt(:)

!Local variables-------------------------------
!scalars
 integer :: iat,natom,my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 natom=SIZE(Pawfgrt)

 write(msg,'(3a)')ch10,' === Content of the pawfgrtab datatype === ',ch10
 call wrtout(my_unt,msg,my_mode)
 do iat=1,natom
   write(msg,'(3(2a,i0))')ch10,&
&    ' > For atom number : ',iat,ch10,&
&    '    1+ Max l in Gaunt coefficients ',Pawfgrt(iat)%l_size,ch10,&
&    '    Number of fine FFT points in PAW sphere ',Pawfgrt(iat)%nfgd
   call wrtout(my_unt,msg,my_mode)

   if (my_prtvol>=3) then
     write(msg,'(a,6(a,i2,a))')ch10,&
&      '    rfgd_allocated    : ',Pawfgrt(iat)%rfgd_allocated,ch10,&
&      '    gylm_allocated    : ',Pawfgrt(iat)%gylm_allocated,ch10,&
&      '    gylmgr_allocated  : ',Pawfgrt(iat)%gylmgr_allocated,ch10,&
&      '    gylmgr2_allocated : ',Pawfgrt(iat)%gylmgr2_allocated,ch10,&
&      '    nhatgr_allocated  : ',Pawfgrt(iat)%nhatfr_allocated,ch10,&
&      '    expiqr_allocated  : ',Pawfgrt(iat)%expiqr_allocated,ch10
     call wrtout(my_unt,msg,my_mode)
   end if

!  These huge arrays are not printed out!
!  Pawfgrt(iat)%ifftsph
!  Pawfgrt(iat)%rfgd
!  Pawfgrt(iat)%gylm
!  Pawfgrt(iat)%gylmgr
!  Pawfgrt(iat)%gylmgr2
!  Pawfgrt(ia)%nhatfr
!  Pawfgrt(ia)%expiqr
 end do

end subroutine pawfgrtab_print
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/pawfgrtab_nullify
!! NAME
!! pawfgrtab_nullify
!!
!! FUNCTION
!!  Nullify the pointers in a pawfgrtab datastructure
!!
!! SIDE EFFECTS
!!  Pawfgrt(:) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_paw_toolbox
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgrtab_nullify(Pawfgrt)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawfgrtab_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Pawfgrtab_type),intent(inout) :: Pawfgrt(:)

!Local variables-------------------------------
!scalars
 integer :: iat,natom

! *************************************************************************

 !@Pawfgrtab_type
 natom=SIZE(Pawfgrt)
 do iat=1,natom
  nullify(Pawfgrt(iat)%ifftsph)
  nullify(Pawfgrt(iat)%gylm   )
  nullify(Pawfgrt(iat)%gylmgr )
  nullify(Pawfgrt(iat)%gylmgr2)
  nullify(Pawfgrt(iat)%nhatfr )
  nullify(Pawfgrt(iat)%rfgd   )
  nullify(Pawfgrt(iat)%expiqr )
  Pawfgrt(iat)%nfgd              = 0
  Pawfgrt(iat)%gylm_allocated    = 0
  Pawfgrt(iat)%gylmgr_allocated  = 0
  Pawfgrt(iat)%gylmgr2_allocated = 0
  Pawfgrt(iat)%nhatfr_allocated  = 0
  Pawfgrt(iat)%rfgd_allocated    = 0
  Pawfgrt(iat)%expiqr_allocated  = 0
 end do

end subroutine pawfgrtab_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/reset_paw_an_flags
!! NAME
!! reset_paw_an_flags
!!
!! FUNCTION
!!  Reset the flags in a paw_an datastructure
!!
!! SIDE EFFECTS
!!  Paw_an_flags<type(Paw_an_flags_type)>=flags in a paw_an datastructure
!!                                        All "has" flags set to 0.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine reset_paw_an_flags(Paw_an_flags)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'reset_paw_an_flags'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_an_flags_type),intent(inout) :: Paw_an_flags

! *************************************************************************

 ! @Paw_an_flags_type
 Paw_an_flags%has_kxc     =0
 Paw_an_flags%has_vhartree=0
 Paw_an_flags%has_vxc     =0
 Paw_an_flags%has_vxcval  =0

end subroutine reset_paw_an_flags
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/copy_paw_an_flags
!! NAME
!! copy_paw_an_flags
!!
!! FUNCTION
!!  Copy the flags in a paw_an datastructure
!!
!! INPUTS
!!  an_flags_in<type(Paw_an_flags_type)>=(input) flags in a paw_an datastructure
!!  Copy a Paw_an_flags_type
!!
!! OUTPUT
!!  an_flags_out<type(Paw_an_flags_type)>=(output) flags in a paw_an datastructure
!!  Copy a Paw_an_flags_type
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_paw_an_flags(an_flags_in, an_flags_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_paw_an_flags'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Paw_an_flags_type),intent(in)    :: an_flags_in
 type(Paw_an_flags_type),intent(inout) :: an_flags_out

! *************************************************************************

 ! @Paw_an_flags_type
 an_flags_out%has_kxc      = an_flags_in%has_kxc
 an_flags_out%has_vhartree = an_flags_in%has_vhartree
 an_flags_out%has_vxc      = an_flags_in%has_vxc
 an_flags_out%has_vxcval   = an_flags_in%has_vxcval

end subroutine copy_paw_an_flags
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/destroy_pawang
!! NAME
!! destroy_pawang
!!
!! FUNCTION
!!  Free all dynamic memory and reset all flags stored in a pawang datastructure
!!
!! SIDE EFFECTS
!!  Pawang <type(pawang_type)>=ANGular mesh discretization and related data
!!
!! PARENTS
!!      loper3,pawalloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_pawang(Pawang)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_pawang'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Pawang_type),intent(inout) :: Pawang

!Local variables-------------------------------
 integer :: ierr

! *************************************************************************

 DBG_ENTER("COLL")

 !@Pawang_type
 if (associated(pawang%angwgth))    then
   ABI_DEALLOCATE(pawang%angwgth)
 end if
 if (associated(pawang%anginit))    then
   ABI_DEALLOCATE(pawang%anginit)
 end if
 if (associated(pawang%zarot))      then
   ABI_DEALLOCATE(pawang%zarot)
 end if
 if (associated(pawang%gntselect))  then
   ABI_DEALLOCATE(pawang%gntselect)
   ierr = ABI_ALLOC_STAT
 end if
 if (associated(pawang%realgnt))    then
   ABI_DEALLOCATE(pawang%realgnt)
   ierr = ABI_ALLOC_STAT
 end if
 if (associated(pawang%ylmr))       then
   ABI_DEALLOCATE(pawang%ylmr)
 end if
 if (associated(pawang%ylmrgr))     then
   ABI_DEALLOCATE(pawang%ylmrgr)
 end if
 if (associated(pawang%ls_ylm))     then
   ABI_DEALLOCATE(pawang%ls_ylm)
 end if

 pawang%angl_size =0
 pawang%ylm_size  =0
 pawang%gnt_option=0
 pawang%use_ls_ylm=0

 DBG_EXIT("COLL")

end subroutine destroy_pawang
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/init_pawang
!! NAME
!! init_pawang
!!
!! FUNCTION
!!  Initialize a pawang datastructure
!!
!! INPUTS
!!  gnt_option=flag activated if pawang%gntselect and pawang%realgnt have to be allocated
!!             also determine the size of these pointers
!!  lmax=maximum value of angular momentum l+1
!!  ngnt=number of non-zero Gaunt coefficients
!!  nphi,ntheta=dimensions of paw angular mesh
!!  nsym=number of symetries
!   pawxcdev=choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  use_ls_ylm=flag activated if pawang%ls_ylm has to be allocated
!!  use_ylm=flag activated if pawang%ylmr and pawang%ylmrgr have to be allocated
!!  xclevel=XC functional level
!!
!! OUTPUT
!!  Pawang <type(pawang_type)>=ANGular mesh discretization and related data
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_pawang(Pawang,gnt_option,lmax,ngnt,nphi,nsym,ntheta,pawxcdev,use_ls_ylm,use_ylm,xclevel)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_pawang'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Pawang_type),intent(inout) :: Pawang
 integer,intent(in) :: gnt_option,lmax,ngnt,nphi,nsym,ntheta,pawxcdev,use_ls_ylm,use_ylm,xclevel

!Local variables-------------------------------

! *************************************************************************

 !@Pawang_type

 Pawang%nphi=nphi
 Pawang%ntheta=ntheta
 Pawang%angl_size=0
 if (pawxcdev==0) Pawang%angl_size=ntheta*nphi
 if (Pawang%angl_size>0) then
   ABI_ALLOCATE(Pawang%anginit,(3,Pawang%angl_size))
   ABI_ALLOCATE(Pawang%angwgth,(Pawang%angl_size))
 else
   nullify(Pawang%anginit)
   nullify(Pawang%angwgth)
 end if

 Pawang%l_max=lmax
 Pawang%l_size_max=2*lmax-1
 if (use_ylm>0) then
   if (xclevel==2) Pawang%ylm_size=(Pawang%l_size_max+1)**2
   if (xclevel/=2) Pawang%ylm_size= Pawang%l_size_max   **2
   ABI_ALLOCATE(Pawang%ylmr,(Pawang%ylm_size,Pawang%angl_size))
   if (xclevel==2) then
     ABI_ALLOCATE(Pawang%ylmrgr,(3,Pawang%ylm_size,Pawang%angl_size))
   else
     nullify(Pawang%ylmrgr)
   end if
 else
   Pawang%ylm_size=0
   nullify(Pawang%ylmr)
   nullify(Pawang%ylmrgr)
 end if

 Pawang%ngnt=ngnt
 Pawang%gnt_option=gnt_option
 if (gnt_option>0) then
  if (gnt_option==2) then
    ABI_ALLOCATE(Pawang%gntselect,((2*Pawang%l_size_max-1)**2,((2*Pawang%l_max-1)**2)*((2*Pawang%l_max-1)**2+1)/2))
  else
    ABI_ALLOCATE(Pawang%gntselect,((Pawang%l_size_max)**2,(Pawang%l_max**2)*(Pawang%l_max**2+1)/2))
  end if
  if (ngnt>0) then
    ABI_ALLOCATE(Pawang%realgnt,(ngnt))
  else
    nullify(Pawang%realgnt)
  end if
 else
   nullify(Pawang%gntselect)
   nullify(Pawang%realgnt)
 end if

 Pawang%use_ls_ylm=use_ls_ylm
 if (use_ls_ylm>0) then
  ABI_ALLOCATE(pawang%ls_ylm,(2,Pawang%l_max**2*(Pawang%l_max**2+1)/2,2))
 else
  nullify(pawang%ls_ylm)
 end if

 Pawang%nsym=nsym
 if (nsym>0) then
   ABI_ALLOCATE(Pawang%zarot,(Pawang%l_size_max,Pawang%l_size_max,Pawang%l_max,nsym))
 else
   nullify(Pawang%zarot)
 end if

end subroutine init_pawang
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/get_dimcprj
!! NAME
!! get_dimcprj
!!
!! FUNCTION
!!  Helper function returning the number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom.
!!  Used to initialize the dimensioning array that is passed to the cprj_alloc routines when the
!!  Cprj_type structure is allocated and initialized.
!!
!! INPUTS
!! Pawtab(ntypat)<pawtab_type>=PAW tabulated starting data.
!! Cryst<Crystal_structure>=Structure gathering info on the unit cell and the symmetries.
!! sort_mode(len=*)=String defining the sorting of the atoms in the Cprj arrays.
!!   Two modes are possible:
!!   -- "O[rdered]", if atoms are sorted by atom type.
!!   -- "R[andom]", if atoms are sorted randomly i.e. according the values of typat specified in the input file.
!!
!! OUTPUT
!!  dimcprj(natom)=Number of nlm elements in the <p_{lmn}^i|\psi> matrix elements for i=1,...,natom.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_dimcprj(Pawtab,Cryst,sort_mode,dimcprj)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_dimcprj'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(crystal_structure),intent(in) :: Cryst
 character(len=*),intent(in) :: sort_mode
!arrays
 integer,intent(out) :: dimcprj(Cryst%natom)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)

!Local variables-------------------------------
 integer :: iatom,itypat
 character(len=500) :: msg

! *************************************************************************

 SELECT CASE (toupper(sort_mode(1:1)))

 CASE ("O") ! Ordered by atom-type

  iatom=0
  do itypat=1,Cryst%ntypat
   dimcprj(iatom+1:iatom+Cryst%nattyp(itypat))=Pawtab(itypat)%lmn_size
   iatom=iatom+Cryst%nattyp(itypat)
  end do

 CASE ("R") ! Randomly ordered (typat from input file)

  do iatom=1,Cryst%natom
   itypat=Cryst%typat(iatom)
   dimcprj(iatom)=Pawtab(itypat)%lmn_size
  end do

 CASE DEFAULT
  msg = " Wrong value for sort_mode: "//TRIM(sort_mode)
  MSG_ERROR(msg)
 END SELECT

end subroutine get_dimcprj
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/init_paw_pwave_lmn
!! NAME
!! init_paw_pwave_lmn
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!    screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_paw_pwaves_lmn(Paw_onsite,natom,ntypat,typat,rprimd,xcart,Psps,Pawtab,Pawrad,local_pawfgrtab,optgrad)

 use defs_basis
 use m_radmesh,  only: nderiv_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_paw_pwaves_lmn'
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ntypat,natom
 type(pseudopotential_type),intent(in) :: Psps
 integer,intent(in),optional :: optgrad
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: xcart(3,natom),rprimd(3,3)
 type(pawtab_type),intent(in) :: Pawtab(ntypat)
 type(pawrad_type),intent(in) :: Pawrad(ntypat)
 type(pawfgrtab_type),intent(in) :: local_pawfgrtab(natom)
 type(paw_pwaves_lmn_t),intent(out) :: Paw_onsite(natom)

!Local variables-------------------------------
!scalars
 integer :: itypat,ln_size,lmn_size,mesh_size,inl,iatom
 integer :: nfgd,ifgd,ipsang,option_ylmr,normchoice,ii,jlmn,jl,jm,jlm,jln
 real(dp) :: phj,rR,tphj,ybcbeg,ybcend
!arrays
 integer, allocatable :: iperm(:)
 real(dp) :: yvals(4),red(3),shift(3)
 real(dp),allocatable :: ff(:),nrm(:),nrm_sort(:),phigrd(:,:),tphigrd(:,:),ylm_tmp(:,:),ylm(:,:),ylm_gr(:,:,:)
 real(dp),allocatable :: rsph_red(:,:),rsph_cart(:,:),phigrd_gr(:,:),tphigrd_gr(:,:),gg(:,:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_spline(:)

! *************************************************************************

 ABI_CHECK(all(local_pawfgrtab(:)%rfgd_allocated==1),"init_paw_pwaves_lmn: R vectors not allocated in pawfgrtab")

 ! @paw_pwaves_lmn_t
 !
 ! Prepare the spline. Calculate 2nd derivatives of partial waves for each atom type.
 ABI_ALLOCATE(Paw_spline,(ntypat))

 do itypat=1,ntypat
   ln_size  =Pawtab(itypat)%basis_size
   mesh_size=Pawrad(itypat)%mesh_size

   ABI_ALLOCATE(Paw_spline(itypat)%phi ,(mesh_size,ln_size))
   ABI_ALLOCATE(Paw_spline(itypat)%tphi,(mesh_size,ln_size))

   do inl=1,ln_size ! Calculate 2nd derivatives of %phi and %tphi for each ln component.
     ybcbeg=zero; ybcend=zero
     call spline(Pawrad(itypat)%rad,Pawtab(itypat)%phi(:,inl), mesh_size,ybcbeg,ybcend,Paw_spline(itypat)%phi(:,inl))

     ybcbeg=zero; ybcend=zero
     call spline(Pawrad(itypat)%rad,Pawtab(itypat)%tphi(:,inl),mesh_size,ybcbeg,ybcend,Paw_spline(itypat)%tphi(:,inl))
   end do
 end do
 !
 ! === spline for each atom ===
 ! * FFT points within PAW sphere depend on the atom site.
 do iatom=1,natom
   itypat   = typat(iatom)
   ln_size  = Pawtab(itypat)%basis_size
   lmn_size = Pawtab(itypat)%lmn_size
   mesh_size= Pawrad(itypat)%mesh_size
   nfgd     = local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere

   ! The points in the PAW sphere might belong to a different unit cell, in this case one has to
   ! reconstruct the contribution to the AE psi_k given by the lattice-symmetric atom in another sphere.
   ! Here I wrap rr back into the first unit cell keeping trace of the lattice vector L needed so that rr = rr_first_ucell + L
   ! The contribution to the AE u(r) in the first unit cell has to be multiplied by e^{-ikL}.

   Paw_onsite(iatom)%nfgd     = nfgd
   Paw_onsite(iatom)%lmn_size = lmn_size
   !
   ABI_ALLOCATE(Paw_onsite(iatom)%r0shift,(3,nfgd))

   ABI_ALLOCATE(rsph_red,(3,nfgd))
   ABI_ALLOCATE(rsph_cart,(3,nfgd))
   do ifgd=1,nfgd
     rsph_cart(:,ifgd) = local_pawfgrtab(iatom)%rfgd(:,ifgd) + xcart(:,iatom)
   end do
   call xredxcart(nfgd,-1,rprimd,rsph_cart,rsph_red) ! go to reduced coordinates.
   ABI_DEALLOCATE(rsph_cart)

   do ifgd=1,nfgd
     call wrap2_zero_one(rsph_red(1,ifgd),red(1),shift(1)) ! rr = r_cell + shift
     call wrap2_zero_one(rsph_red(2,ifgd),red(2),shift(2))
     call wrap2_zero_one(rsph_red(3,ifgd),red(3),shift(3))
     Paw_onsite(iatom)%r0shift(:,ifgd) = NINT(shift)
     !if (ANY( ABS(shift) > tol12)) then
       !MSG_WARNING("rmR_red is outside the first unit cell.")
       !write(ab_out,*)rsph_red(:,ifgd),shift
     !end if
   end do
   ABI_DEALLOCATE(rsph_red)
   !
   ! * Obtain |r-R| on fine grid, note that rfgd is given in Cartesian coordinates.
   ABI_ALLOCATE(nrm,(nfgd))
   do ifgd=1,nfgd
     nrm(ifgd) = sqrt(dot_product(local_pawfgrtab(iatom)%rfgd(:,ifgd),local_pawfgrtab(iatom)%rfgd(:,ifgd)))
   end do
   !
   ! * Compute Ylm for each r-R vector.
   ipsang = 1 + (Pawtab(itypat)%l_size-1)/2 ! recall l_size=2*l_max-1 where l_max is shifted by 1.
   ABI_ALLOCATE(ylm_tmp,(ipsang**2,nfgd))
   normchoice = 1 ! Use computed norms of input vectors.
   if (present(optgrad) .and. optgrad==1) then
     option_ylmr=2 ! Compute Ylm(r-R) and its gradient
     ABI_ALLOCATE(ylm_gr,(3,ipsang**2,nfgd))
   else
     option_ylmr= 1 ! To compute Ylm(r-R).
     ABI_ALLOCATE(ylm_gr,(3,3,0))
   end if
   call initylmr(ipsang,normchoice,nfgd,nrm,option_ylmr,local_pawfgrtab(iatom)%rfgd,ylm_tmp,ylm_gr)
   !
   !  Exchange dimensions for better memory access.
   ABI_ALLOCATE(ylm,(nfgd,ipsang**2))
   do ii=1,ipsang**2
     ylm(:,ii) = ylm_tmp(ii,:)
   end do
   ABI_DEALLOCATE(ylm_tmp)
   !
   ! In order to do spline fits, the |r-R| data must be sorted
   ! Here we sort the nrm points, keeping track of which goes where
   ABI_ALLOCATE(nrm_sort,(nfgd))
   nrm_sort = nrm

   ABI_ALLOCATE(iperm,(nfgd))
   do ifgd=1,nfgd
     iperm(ifgd)=ifgd
   end do

   call sort_dp(nfgd,nrm_sort,iperm,tol8)

   ! Now make spline fits of phi and tphi  onto the fine grid around the atom
   ABI_ALLOCATE(phigrd,(nfgd,ln_size))
   ABI_ALLOCATE(tphigrd,(nfgd,ln_size))
   ABI_ALLOCATE(ff,(nfgd))
   if (present(optgrad) .and. optgrad==1) then
     ABI_ALLOCATE(phigrd_gr,(nfgd,ln_size))
     ABI_ALLOCATE(tphigrd_gr,(nfgd,ln_size))
     ABI_ALLOCATE(gg,(mesh_size,1))
   end if

   do inl=1,ln_size
     !
     ! * splint phi onto points and reorder indices.
     call splint(mesh_size,Pawrad(itypat)%rad,Pawtab(itypat)%phi(:,inl),Paw_spline(itypat)%phi(:,inl),nfgd,nrm_sort,ff)
     do ifgd=1,nfgd
       ii=iperm(ifgd)
       phigrd(ii,inl) = ff(ifgd)
     end do
     !
     ! * compute d phi/dr, interpolate onto points and reorder indices.
     if (present(optgrad) .and. optgrad==1) then
       ybcbeg=zero; ybcend=zero
       call nderiv_gen(gg,Pawtab(itypat)%phi(:,inl),1,Pawrad(itypat))
       call spline(Pawrad(itypat)%rad,gg(:,1), mesh_size,ybcbeg,ybcend,Paw_spline(itypat)%phi(:,inl))
       call splint(mesh_size,Pawrad(itypat)%rad,Paw_spline(itypat)%phi(:,inl),Paw_spline(itypat)%phi(:,inl),nfgd,nrm_sort,ff)
       do ifgd=1,nfgd
         ii=iperm(ifgd)
         phigrd_gr(ii,inl)  = ff(ifgd)
       end do
     end if
     !
     ! * compute d tphi/dr, interpolate onto points and reorder indices.
     call splint(mesh_size,Pawrad(itypat)%rad,Pawtab(itypat)%tphi(:,inl),Paw_spline(itypat)%tphi(:,inl),nfgd,nrm_sort,ff)
     do ifgd=1,nfgd
       ii=iperm(ifgd)
       tphigrd(ii,inl) = ff(ifgd)
     end do
     if (present(optgrad) .and. optgrad==1) then
       ybcbeg=zero; ybcend=zero
       call nderiv_gen(gg,Pawtab(itypat)%tphi(:,inl),1,Pawrad(itypat))
       call spline(Pawrad(itypat)%rad,gg(:,1), mesh_size,ybcbeg,ybcend,Paw_spline(itypat)%tphi(:,inl))
       call splint(mesh_size,Pawrad(itypat)%rad,Paw_spline(itypat)%tphi(:,inl),Paw_spline(itypat)%phi(:,inl),nfgd,nrm_sort,ff)
       do ifgd=1,nfgd
         ii=iperm(ifgd)
         tphigrd_gr(ii,inl)  = ff(ifgd)
       end do
     end if
   end do !inl

   ABI_DEALLOCATE(ff)
   if (present(optgrad) .and. optgrad==1) then
     ABI_DEALLOCATE(gg)
   end if
   !
   ! === Calculate AE and PS partial waves inside the sphere ===
   ! * recall that <r|phi>=u(r)*Slm(r^)/r, hence avoid division by zero except for s-waves.
   ABI_ALLOCATE(Paw_onsite(iatom)%phi ,(nfgd,lmn_size))
   ABI_ALLOCATE(Paw_onsite(iatom)%tphi,(nfgd,lmn_size))

   if (present(optgrad) .and. optgrad==1) then
     ABI_ALLOCATE(Paw_onsite(iatom)%phi_gr ,(3,nfgd,lmn_size))
     ABI_ALLOCATE(Paw_onsite(iatom)%tphi_gr,(3,nfgd,lmn_size))
   end if

   do jlmn=1,lmn_size
     jl  = Psps%indlmn(1,jlmn,itypat)
     jm  = Psps%indlmn(2,jlmn,itypat)
     jlm = Psps%indlmn(4,jlmn,itypat)
     jln = Psps%indlmn(5,jlmn,itypat)

     do ifgd=1,nfgd ! loop over fine grid points in current PAW sphere
       !if (nrm(ifgd)>tol16) then
       if (nrm(ifgd)>tol10) then ! tol10 to be consistent with initylmr.
         rR  = nrm(ifgd) ! value of |r-R|
         !write(ab_out,*) 'rR:',rR,' phigrd:',phigrd(ifgd,jln),' tphigrd:',tphigrd(ifgd,jln),' ylm:',ylm(ifgd,jlm)
         phj = phigrd (ifgd,jln)*ylm(ifgd,jlm)/rR
         tphj= tphigrd(ifgd,jln)*ylm(ifgd,jlm)/rR
         Paw_onsite(iatom)%phi (ifgd,jlmn) =  phj
         Paw_onsite(iatom)%tphi(ifgd,jlmn) = tphj

         if (present(optgrad) .and. optgrad==1) then
           Paw_onsite(iatom)%phi_gr (1:3,ifgd,jlmn) = phigrd (ifgd,jln)*ylm_gr(1:3,jlm,ifgd) &
&            + phigrd_gr(ifgd,jln)*local_pawfgrtab(iatom)%rfgd(1:3,ifgd)*ylm(ifgd,jlm)
           Paw_onsite(iatom)%tphi_gr (1:3,ifgd,jlmn) = tphigrd (ifgd,jln)*ylm_gr(1:3,jlm,ifgd) &
&            + tphigrd_gr(ifgd,jln)*local_pawfgrtab(iatom)%rfgd(1:3,ifgd)*ylm(ifgd,jlm)
         end if

       else ! Extrapolate if the point is at the origin
         yvals(1) = zero
         if (jl==0) then
           yvals(2:4) = Pawtab(itypat)%phi(2:4,jln)/Pawrad(itypat)%rad(2:4)
           call deducer0(yvals,4,Pawrad(itypat))
         end if
         Paw_onsite(iatom)%phi(ifgd,jlmn) = yvals(1) * ylm(ifgd,jlm)
         yvals(1) = zero
         if (jl==0) then
           yvals(2:4) = Pawtab(itypat)%tphi(2:4,jln)/Pawrad(itypat)%rad(2:4)
           call deducer0(yvals,4,pawrad(itypat))
         end if
         Paw_onsite(iatom)%tphi(ifgd,jlmn) = yvals(1) * ylm(ifgd,jlm)
         ! The gradient is expected to go to zero at the origin
         if (present(optgrad) .and. optgrad==1) then
           Paw_onsite(iatom)%phi_gr (1:3,ifgd,jlmn) = zero
           Paw_onsite(iatom)%tphi_gr (1:3,ifgd,jlmn) = zero
         end if
       end if

     end do !nfgd
   end do !jlmn

   ABI_DEALLOCATE(nrm)
   ABI_DEALLOCATE(nrm_sort)
   ABI_DEALLOCATE(iperm)
   ABI_DEALLOCATE(phigrd)
   ABI_DEALLOCATE(tphigrd)
   ABI_DEALLOCATE(ylm)
   ABI_DEALLOCATE(ylm_gr)
   if (present(optgrad) .and. optgrad==1) then
     ABI_DEALLOCATE(phigrd_gr)
     ABI_DEALLOCATE(tphigrd_gr)
   end if
 end do !iatom
 !
 !* Free 2nd derivates used for spline.
 call destroy_paw_pwaves_lmn(Paw_spline)
 ABI_DEALLOCATE(Paw_spline)

end subroutine init_paw_pwaves_lmn
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/destroy_paw_pwave_lmn
!! NAME
!! destroy_paw_pwave_lmn
!!
!! FUNCTION
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

subroutine destroy_paw_pwaves_lmn(Paw_onsite)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_paw_pwaves_lmn'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(paw_pwaves_lmn_t),intent(inout) :: Paw_onsite(:)

!Local variables-------------------------------
!scalars
 integer :: iatom

! *************************************************************************

 ! @paw_pwaves_lmn_t
 do iatom=LBOUND(Paw_onsite,DIM=1),UBOUND(Paw_onsite,DIM=1)
   if (associated(Paw_onsite(iatom)%phi ))   then
     ABI_DEALLOCATE(Paw_onsite(iatom)%phi)
   end if
   if (associated(Paw_onsite(iatom)%tphi))   then
     ABI_DEALLOCATE(Paw_onsite(iatom)%tphi)
   end if
   if (associated(Paw_onsite(iatom)%r0shift))   then
     ABI_DEALLOCATE(Paw_onsite(iatom)%r0shift)
   end if
   if (associated(Paw_onsite(iatom)%phi_gr ))   then
     ABI_DEALLOCATE(Paw_onsite(iatom)%phi_gr)
   end if
   if (associated(Paw_onsite(iatom)%tphi_gr))   then
     ABI_DEALLOCATE(Paw_onsite(iatom)%tphi_gr)
   end if
 end do

end subroutine destroy_paw_pwaves_lmn
!!***

!----------------------------------------------------------------------

!!****f* m_paw_toolbox/nullify_paw_pwave_lmn
!! NAME
!! nullify_paw_pwave_lmn
!!
!! FUNCTION
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

subroutine nullify_paw_pwaves_lmn(Paw_onsite)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_paw_pwaves_lmn'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(paw_pwaves_lmn_t),intent(inout) :: Paw_onsite(:)

!Local variables-------------------------------
!scalars
 integer :: iatom

! *************************************************************************

 ! @paw_pwaves_lmn_t
 do iatom=LBOUND(Paw_onsite,DIM=1),UBOUND(Paw_onsite,DIM=1)
   nullify(Paw_onsite(iatom)%phi)
   nullify(Paw_onsite(iatom)%tphi)
   nullify(Paw_onsite(iatom)%r0shift)
   nullify(Paw_onsite(iatom)%phi_gr)
   nullify(Paw_onsite(iatom)%tphi_gr)
 end do

end subroutine nullify_paw_pwaves_lmn
!!***

!----------------------------------------------------------------------

END MODULE m_paw_toolbox
!!***
