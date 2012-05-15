!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_scf_history
!! NAME
!!  m_scf_history
!!
!! FUNCTION
!!  This module provides the definition of the scf_history_type used to store
!!  various arrays obtained from previous SCF cycles (density, positions...).
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_scf_history

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

 implicit none

 private

! public procedures.
 public :: init_scf_history
 public :: destroy_scf_history
 public :: nullify_scf_history
!!***

!!****t* m_scf_history/scf_history_type
!! NAME
!! scf_history_type
!!
!! FUNCTION
!! This structured datatype contains various arrays obtained from
!! previous SCF cycles (density, positions...)
!!
!! SOURCE

 type, public :: scf_history_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: history_size
   ! Number of previous SCF cycles stored in history
   ! If history_size<0, scf_history is not used
   ! If history_size=0, scf_history only contains
   !    current values of data (rhor, taur, pawrhoih, xred)
   ! If history_size>0, scf_history contains
   !    current values of data and also
   !    history_size previous values of these data

  integer :: icall
   ! Number of call for the routine extraprho

  integer :: mcg
   ! Size of cg array

  integer :: mcprj
   ! Size of cprj datsatructure array

  integer :: natom
   ! Number of atoms in cell

  integer :: nfft
   ! Size of FFT grid (for density)

  integer :: nspden
   ! Number of independant spin components for density

  integer :: usecg
   ! usecg=1 if the extrapolation of the wavefunctions is active

  real(dp) :: alpha
   ! alpha mixing coefficient for the prediction of density and wavefunctions

  real(dp) :: beta
   ! beta mixing coefficient for the prediction of density and wavefunctions

! Integer arrays

  integer,pointer :: hindex(:)
   ! hindex(history_size)
   ! Indexes of SCF cycles in the history
   ! hindex(1) is the newest SCF cycle
   ! hindex(history_size) is the oldest SCF cycle

! Real (real(dp)) arrays

   real(dp) :: rprimd(3,3)

   real(dp),pointer :: cg(:,:,:)
    ! cg(2,mcg,history_size)
    ! wavefunction coefficients needed for each SCF cycle of history

   real(dp),pointer :: deltarhor(:,:,:)
    ! deltarhor(nfft,nspden,history_size)
    ! Diference between electronic density (in real space)
    ! and sum of atomic densities at the end of each SCF cycle of history

   real(dp),pointer :: atmrho_last(:)
    ! atmrho_last(nfft)
    ! Sum of atomic densities at the end of the LAST SCF cycle

   real(dp),pointer :: rhor_last(:,:)
    ! rhor_last(nfft,nspden)
    ! Last computed electronic density (in real space)

   real(dp),pointer :: taur_last(:,:)
    ! taur_last(nfft,nspden*usekden)
    ! Last computed kinetic energy density (in real space)

   real(dp),pointer :: xreddiff(:,:,:)
    ! xreddiff(3,natom,history_size)
    ! Difference of reduced coordinates of atoms between a
    ! SCF cycle and the previous

   real(dp),pointer :: xred_last(:,:)
    ! xred_last(3,natom)
    ! Last computed atomic positions (reduced coordinates)

! Structured datatypes arrays

  type(pawrhoij_type), pointer :: pawrhoij(:,:)
    ! pawrhoij(natom,history_size)
    ! PAW only: occupancies matrix at the end of each SCF cycle of history

  type(pawrhoij_type), pointer :: pawrhoij_last(:)
    ! pawrhoij_last(natom)
    ! PAW only: last computed occupancies matrix

  type(cprj_type),pointer :: cprj(:,:,:)
    !cprj(natom,nspinor*mband*mkmem*nsppol,history_size)

 end type scf_history_type
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_scf_history/init_scf_history
!! NAME
!!  init_scf_history
!!
!! FUNCTION
!!  Init all scalars and pointers in a scf_history datastructure
!!  according to scf_history%history_size value which has to be
!!  defined before calling this routine
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=MPI-parallelisation information
!!
!! SIDE EFFECTS
!!  scf_history=<type(scf_history_type)>=scf_history datastructure
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_scf_history(dtset,mpi_enreg,scf_history)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_scf_history'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(inout) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 type(scf_history_type),intent(inout) :: scf_history
!Local variables-------------------------------
!scalars
 integer :: jj,my_nspinor,nfft
!arrays

!************************************************************************

 !@scf_history_type

if (scf_history%history_size<0) then
  call nullify_scf_history(scf_history)
else

  nfft=dtset%nfft
  if (dtset%usepaw==1.and.(dtset%pawecutdg>=1.0000001_dp*dtset%ecut)) nfft=dtset%nfftdg

  if (scf_history%history_size>=0) then
    ABI_ALLOCATE(scf_history%rhor_last,(nfft,dtset%nspden))
    ABI_ALLOCATE(scf_history%taur_last,(nfft,dtset%nspden*dtset%usekden))
    ABI_ALLOCATE(scf_history%xred_last,(3,dtset%natom))
    ABI_ALLOCATE(scf_history%pawrhoij_last,(mpi_enreg%natom*dtset%usepaw))
    if (dtset%usepaw==1) then
      call rhoij_nullify(scf_history%pawrhoij_last)
    end if
  end if

  if (scf_history%history_size>0) then

    scf_history%natom=dtset%natom
    scf_history%nfft=nfft
    scf_history%nspden=dtset%nspden
    scf_history%alpha=zero
    scf_history%beta=zero
    scf_history%usecg=0
    scf_history%icall=0
    scf_history%mcg=0
    scf_history%mcprj=0

    ABI_ALLOCATE(scf_history%hindex,(scf_history%history_size))
    scf_history%hindex(:)=0
    ABI_ALLOCATE(scf_history%deltarhor,(nfft,dtset%nspden,scf_history%history_size))
    ABI_ALLOCATE(scf_history%xreddiff,(3,dtset%natom,scf_history%history_size))
    ABI_ALLOCATE(scf_history%atmrho_last,(nfft))

    if (dtset%usepaw==1) then
      ABI_ALLOCATE(scf_history%pawrhoij,(mpi_enreg%natom,scf_history%history_size))
      do jj=1,scf_history%history_size
        call rhoij_nullify(scf_history%pawrhoij(:,jj))
      end do
    else
      nullify(scf_history%pawrhoij)
    end if

    if (dtset%iextrapwf==1) then
      my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
      scf_history%usecg=1
      scf_history%mcg=my_nspinor*dtset%mpw*dtset%mband*dtset%mkmem*dtset%nsppol
      ABI_ALLOCATE(scf_history%cg,(2,scf_history%mcg,scf_history%history_size))
      if (dtset%usepaw==1) then
        scf_history%mcprj=my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
        ABI_ALLOCATE(scf_history%cprj,(dtset%natom,scf_history%mcprj,scf_history%history_size))
        do jj=1,scf_history%history_size
          call cprj_nullify(scf_history%cprj(:,:,jj))
        end do
      else
        nullify(scf_history%cprj)
      end if
    else
      nullify(scf_history%cg)
      nullify(scf_history%cprj)
    end if

  end if
end if

end subroutine init_scf_history
!!***

!----------------------------------------------------------------------

!!****f* m_scf_history/destroy_scf_history
!! NAME
!!  destroy_scf_history
!!
!! FUNCTION
!!  Clean and destroy a scf_history datastructure
!!
!! SIDE EFFECTS
!!  scf_history(:)=<type(scf_history_type)>=scf_history datastructure
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_scf_history(scf_history)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_scf_history'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(scf_history_type),intent(inout) :: scf_history
!Local variables-------------------------------
!scalars
 integer :: jj

!************************************************************************

 !@scf_history_type

 if (associated(scf_history%pawrhoij_last)) then
   call rhoij_free(scf_history%pawrhoij_last)
   ABI_DEALLOCATE(scf_history%pawrhoij_last)
 end if
 if (associated(scf_history%pawrhoij)) then
   do jj=1,size(scf_history%pawrhoij,2)
     call rhoij_free(scf_history%pawrhoij(:,jj))
   end do
   ABI_DEALLOCATE(scf_history%pawrhoij)
 end if
 if (associated(scf_history%cprj)) then
   do jj=1,size(scf_history%cprj,2)
     call cprj_free(scf_history%cprj(:,:,jj))
   end do
   ABI_DEALLOCATE(scf_history%cprj)
 end if

 if (associated(scf_history%hindex))       then
   ABI_DEALLOCATE(scf_history%hindex)
 end if
 if (associated(scf_history%deltarhor))    then
   ABI_DEALLOCATE(scf_history%deltarhor)
 end if
 if (associated(scf_history%xreddiff))     then
   ABI_DEALLOCATE(scf_history%xreddiff)
 end if
 if (associated(scf_history%atmrho_last))  then
   ABI_DEALLOCATE(scf_history%atmrho_last)
 end if
 if (associated(scf_history%xred_last))    then
   ABI_DEALLOCATE(scf_history%xred_last)
 end if
 if (associated(scf_history%rhor_last))    then
   ABI_DEALLOCATE(scf_history%rhor_last)
 end if
 if (associated(scf_history%taur_last))    then
   ABI_DEALLOCATE(scf_history%taur_last)
 end if
 if (associated(scf_history%cg))           then
   ABI_DEALLOCATE(scf_history%cg)
 end if

 scf_history%history_size=-1
 scf_history%usecg=0
 scf_history%icall=0
 scf_history%mcprj=0
 scf_history%mcg=0

end subroutine destroy_scf_history
!!***

!----------------------------------------------------------------------

!!****f* m_scf_history/nullify_scf_history
!! NAME
!!  nullify_scf_history
!!
!! FUNCTION
!!  Nullify (set to null) an scf_history datastructure
!!
!! SIDE EFFECTS
!!  scf_history(:)=<type(scf_history_type)>=scf_history datastructure
!!
!! PARENTS
!!      gstateimg,m_scf_history
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_scf_history(scf_history)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_scf_history'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(scf_history_type),intent(inout) :: scf_history
!Local variables-------------------------------
!scalars

!************************************************************************

 !@scf_history_type

 nullify(scf_history%hindex)
 nullify(scf_history%deltarhor)
 nullify(scf_history%xreddiff)
 nullify(scf_history%atmrho_last)
 nullify(scf_history%xred_last)
 nullify(scf_history%rhor_last)
 nullify(scf_history%taur_last)
 nullify(scf_history%pawrhoij_last)
 nullify(scf_history%pawrhoij)
 nullify(scf_history%cprj)
 nullify(scf_history%cg)

 scf_history%history_size=-1
 scf_history%usecg=0
 scf_history%icall=0
 scf_history%mcprj=0
 scf_history%mcg=0

end subroutine nullify_scf_history
!!***

END MODULE m_scf_history
!!***
