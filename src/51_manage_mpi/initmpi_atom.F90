!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_atom
!! NAME
!!  initmpi_atom
!!
!! FUNCTION
!!  Initializes the mpi informations for parallelism over atoms (PAW).
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  mpi_enreg= informations about MPI parallelization
!!
!! OUTPUT
!!  mpi_enreg=informations about MPI parallelization:
!!     mpi_enreg%comm_atom=communicator over atoms
!!     mpi_enreg%nproc_atom=number of procs for parallelisation over atoms
!!     mpi_enreg%atom_indx()=indexes of atoms treated by current proc.
!!
!! PARENTS
!!      gstate,respfn
!!
!! CHILDREN
!!      mpi_cart_sub,mpi_comm_size
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_atom(dtset,mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_atom'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: iatom
 logical,parameter :: activate_paral_atom=.false.
 integer,allocatable :: atom_indx_tmp(:)
#if defined HAVE_MPI
 integer :: ierr,nproc,spacecomm
 logical :: keepdim(3)
#endif

! ***********************************************************************

 DBG_ENTER("COLL")

 mpi_enreg%nproc_atom=1
 mpi_enreg%natom=dtset%natom

#if defined HAVE_MPI

!In case of a "ground-state" calculation (no DFPT, no GW),
!and parallelisation over kpt or "kpt-bands-FFT" parallelisation,
!build communicator over atoms
 if ((mpi_enreg%paralbd == 0).and. &
& (dtset%paral_rf == 0)) then

   if ((mpi_enreg%paral_compil_kpt==1).and.(activate_paral_atom)) then

!    Take into account possible parallelization over images
     if (mpi_enreg%paral_img==0) then
       nproc=mpi_enreg%nproc
       spacecomm=mpi_enreg%world_comm
     else
       nproc=mpi_enreg%nproc_one_img
       spacecomm=mpi_enreg%comm_one_img
     end if

     if (dtset%paral_kgb==0) then
!      kpt parallelisation: communicator(atom)=communicator(kpt)=world communicator
       mpi_enreg%comm_atom=spacecomm
       mpi_enreg%nproc_atom=nproc

     else
!      "kpt-bands-FFT" parallelisation: communicator(atom)=communicator(kpt-bands)
       keepdim(1)=.false.;keepdim(2)=.true.;keepdim(3)=.true.
       call MPI_CART_SUB(spacecomm,keepdim,mpi_enreg%comm_atom,ierr)
       call MPI_COMM_SIZE(mpi_enreg%comm_atom,mpi_enreg%nproc_atom,ierr)
     end if

!    Select indexes of atoms for the current proc.
     if (mpi_enreg%nproc_atom>1) then
       ABI_ALLOCATE(atom_indx_tmp,(dtset%natom))
       mpi_enreg%natom=0
       do iatom=1,dtset%natom
         if (mod(iatom-1,mpi_enreg%nproc_atom)==0) then
           mpi_enreg%natom=mpi_enreg%natom+1
           atom_indx_tmp(mpi_enreg%natom)=iatom
         end if
       end do
       ABI_ALLOCATE(mpi_enreg%atom_indx,(mpi_enreg%natom))
       mpi_enreg%atom_indx(1:mpi_enreg%natom)=atom_indx_tmp(1:mpi_enreg%natom)
       ABI_DEALLOCATE(atom_indx_tmp)
     end if

   else
!    Other cases: communicator(atom)=null
     mpi_enreg%comm_atom=MPI_COMM_SELF
   end if
 else
   mpi_enreg%comm_atom=MPI_COMM_SELF

 end if

#endif

 DBG_EXIT("COLL")

end subroutine initmpi_atom
!!***
