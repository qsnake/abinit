!{\src2tex{textfont=tt}}
!!****f* ABINIT/exc_diago_driver
!! NAME
!!  exc_diago_driver
!!
!! FUNCTION
!!  Driver routine for the direct diagonalization of the Hermitian excitonic Hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT and EXC groups (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  neh=Rank of the resonant block of the Hamiltoninan.
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!    %exh=Name of the file storing the excitonic resonant part.
!!
!! OUTPUT
!!  Eigenvalues and eigenvectors are written on file.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      build_spectra,exc_diago_coupling,exc_diago_coupling_hegv
!!      exc_diago_resonant,exc_iterative_diago,exc_print_eig,reportgap
!!      xbarrier_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine exc_diago_driver(Wfd,Bsp,BS_files,KS_BSt,QP_BSt,Cryst,Kmesh,Psps,&
&  Pawtab,Hur,Hdr_bse,drude_plsmf)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use defs_datatypes
 use m_xmpi
 use m_errors

 use defs_abitypes,     only : hdr_type
 use m_crystal,         only : crystal_structure
 use m_bz_mesh,         only : bz_mesh_type
 use m_ebands,          only : ReportGap
 use m_wfs,             only : wfs_descriptor
 use m_paw_commutator,  only : HUr_commutator

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_diago_driver'
 use interfaces_71_bse, except_this_one => exc_diago_driver
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: drude_plsmf
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) ::  BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
 type(Crystal_structure),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Bandstructure_type),intent(in) :: KS_BSt,QP_BSt
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(HUr_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: my_rank,master,comm,prtvol
 complex(dpc) :: exc_gap,gw_gap
 logical :: eval_eigenstates
 character(len=500) :: msg
!arrays
 real(dp) :: gaps(2,QP_BSt%nsppol)

!************************************************************************

 DBG_ENTER("COLL")

 comm    = Wfd%comm
 my_rank = Wfd%my_rank
 master  = Wfd%master
 prtvol  = Wfd%prtvol

 if (BSp%have_complex_ene) then
   MSG_ERROR("Complex energies are not supported yet")
 end if
 !
 ! This trick is needed to restart a CG run, use DDIAGO to calculate the spectra reusing an old BSEIG file.
 eval_eigenstates = &
&   (BS_files%in_eig == BSE_NOFILE) .or. &
&   (Bsp%algorithm == BSE_ALGO_CG)

 if (eval_eigenstates) then
   !
   select case (BSp%algorithm)
   case (BSE_ALGO_DDIAGO)
     if (BSp%use_coupling==0) then
       call exc_diago_resonant(BSp,BS_files,Hdr_bse,prtvol,comm)
     else
       if (Bsp%have_complex_ene) then 
         ! Solve Hv = ev with generic complex matrix. 
         call exc_diago_coupling(BSp,BS_files,Hdr_bse,prtvol,comm)
       else 
         ! Solve generalized eigenvalue problem F Hbar with Hbar Hermitian definitive positive matrix.
         call exc_diago_coupling_hegv(BSp,BS_files,Hdr_bse,prtvol,comm)
       end if
     end if

   case (BSE_ALGO_CG)
     if (BSp%use_coupling==0) then
       call exc_iterative_diago(Bsp,BS_files,Hdr_bse,prtvol,comm)
     else 
       MSG_ERROR("CG + coupling not coded")
     end if

   case default
     write(msg,'(a,i0)')" Wrong value for Bsp%algorithm: ",Bsp%algorithm
     MSG_ERROR(msg)
   end select
   !
   if (my_rank==master) then 
     call ReportGap(QP_BSt,header="QP bands",unit=std_out,gaps=gaps)
     gw_gap = MINVAL(gaps(2,:))
     call exc_print_eig(BSp,BS_files%out_eig,gw_gap,exc_gap)
   end if
   call xbarrier_mpi(comm)
   !
 end if

 call build_spectra(BSp,BS_files,Cryst,Kmesh,KS_BSt,QP_BSt,Psps,Pawtab,Wfd,Hur,drude_plsmf,comm)

 DBG_EXIT("COLL")

end subroutine exc_diago_driver
!!***
