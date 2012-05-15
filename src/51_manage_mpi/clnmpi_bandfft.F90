!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_bandfft
!! NAME
!!  clnmpi_bandfft
!!
!! FUNCTION
!!  Cleans-up the mpi informations for the BAND-FFT parallelism (paral_kgb=1).
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!    mpi_enreg%bandfft_kpt(:)=informations about band-fft parallelisation
!!
!! PARENTS
!!      gstate,mpi_enreg_tools
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine clnmpi_bandfft(mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_bandfft'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: ikpt_this_proc,isppol,ikpt,mkmem,nkpt,nsppol
 character(len=500) :: msg
#endif

! ***********************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_MPI
 if (associated(mpi_enreg%bandfft_kpt).and.associated(mpi_enreg%tab_kpt_distrib)) then
   mkmem =size(mpi_enreg%bandfft_kpt)
   nkpt=size(mpi_enreg%proc_distrb,1)
   nsppol=size(mpi_enreg%proc_distrb,3)
   if (nsppol==0.or.nkpt==0) then
     msg=' mpi_enreg%proc_distrb should be allocated !'
     MSG_BUG(msg)
   end if
   nkpt=size(mpi_enreg%tab_kpt_distrib)
   if (nkpt==0) then
     msg=' mpi_enreg%tab_kpt_distrib should be allocated !'
     MSG_BUG(msg)
   end if
   do isppol=1,nsppol
     do ikpt=1,nkpt
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,:,isppol)-mpi_enreg%me_kpt))/=0) then
         cycle
       end if
       ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)
       if ((ikpt_this_proc>mkmem) .or.(ikpt_this_proc<=0)) then
         msg=' The bandfft tab cannot be deallocated !'
         MSG_BUG(msg)
       end if
       if (associated(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ind_kg_mpi_to_seq)) then
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ind_kg_mpi_to_seq)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ind_kg_mpi_to_seq)
       end if
       if (associated(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather)) then
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather)
       end if
       if (mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag1_is_allocated==1) then
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%gbound)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%gbound)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls)
         mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag1_is_allocated=0
       end if
       if (mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag2_is_allocated==1) then
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ffnl_gather)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ffnl_gather)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kinpw_gather)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kinpw_gather)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kpg_k_gather)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kpg_k_gather)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ph3d_gather)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ph3d_gather)
         mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag2_is_allocated=0
       end if
       if (mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag3_is_allocated==1) then
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls_sym)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls_sym)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym_tot)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym_tot)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls_sym)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls_sym)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym_all)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym_all)
         ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%tab_proc)
         nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%tab_proc)
         mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag3_is_allocated=0
       end if
     end do
   end do
   ABI_DEALLOCATE(mpi_enreg%bandfft_kpt)
   nullify(mpi_enreg%bandfft_kpt)
   ABI_DEALLOCATE(mpi_enreg%tab_kpt_distrib)
   nullify(mpi_enreg%tab_kpt_distrib)
 end if
#endif

 DBG_EXIT("COLL")

end subroutine clnmpi_bandfft
!!***
