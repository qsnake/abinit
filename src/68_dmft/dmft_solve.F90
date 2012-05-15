!{\src2tex{textfont=tt}}
!!****f* ABINIT/dmft_solve
!! NAME
!! dmft_solve
!!
!! FUNCTION
!! Solve the DMFT loop from PAW data.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_structure)>=crystal structure data
!!  istep           =  step of iteration for LDA.
!!  lda_occup <type(oper_type)> = occupations in the correlated orbitals in LDA
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft <type(paw_dmft_type)> =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pawprtvol  = option for printing
!!
!! OUTPUT
!!  paw_dmft <type(paw_dmft_type)> =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      check_fourier_green,compute_energy,compute_green,dc_self,destroy_energy
!!      destroy_green,destroy_hu,destroy_self,diff_oper,dyson,fermi_green
!!      icip_green,impurity_solve,init_energy,init_green,init_hu
!!      initialize_self,integrate_green,leave_new,local_ks_green,new_self
!!      print_self,printocc_green,psichi_renormalization,rw_self
!!      spectral_function,wrtout,xbarrier_mpi,xcomm_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dmft_solve(cryst_struc,istep,lda_occup,mpi_enreg,paw_dmft,pawang,pawtab,pawprtvol)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

 use m_paw_dmft, only: paw_dmft_type
 use m_crystal, only : crystal_structure
 use m_green, only : green_type, destroy_green, icip_green,init_green,&
&                    print_green,printocc_green,&
&                    integrate_green,copy_green,compute_green,check_fourier_green
 use m_oper, only : oper_type,diff_oper
 use m_self, only : self_type,initialize_self,destroy_self,print_self,dc_self,rw_self,new_self
 use m_hu, only : hu_type,init_hu,destroy_hu
 use m_energy, only : energy_type,init_energy,destroy_energy,compute_energy,print_energy
 use m_matlu, only : print_matlu

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dmft_solve'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_68_dmft, except_this_one => dmft_solve
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: istep
 integer, intent(in) :: pawprtvol
 type(MPI_type), intent(inout) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
 type(crystal_structure),intent(in) :: cryst_struc
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(oper_type),intent(in)  :: lda_occup
 type(paw_dmft_type), intent(inout)  :: paw_dmft

!Local variables ------------------------------
!scalars
 integer :: check,idmftloop,ierr,istep_iter,me,spaceComm
! type
 type(green_type) :: green
 type(green_type) :: greenlda
 type(hu_type),allocatable :: hu(:)
 type(green_type) :: weiss
 type(self_type) :: self
 type(self_type) :: self_new
 type(energy_type) :: energies_dmft
 character(len=500) :: message
!************************************************************************
 check=paw_dmft%dmftcheck ! checks enabled
 paw_dmft%dmft_fepr=tol5
 paw_dmft%dmft_lcpr=tol5
 paw_dmft%dmft_chpr=tol6
!paw_dmft%dmft_chpr=20_dp ! total number of electron.
 paw_dmft%dmft_prgn=0

 if(check==1) then
   write(message,'(2a)') ch10,' DMFT Checks are enabled '
 else
   write(message,'(2a)') ch10,' DMFT Checks will not be performed'
 end if
 call wrtout(std_out,message,'COLL')

 DBG_ENTER("COLL")

 if(istep==0) then
   write(message,'(2a)') ch10,' BUG: istep should not be equal to zero'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_kpt)
 call xme_init(mpi_enreg,me) 
 call xbarrier_mpi(spaceComm)

 call initialize_self(self,paw_dmft)
 call init_energy(cryst_struc,energies_dmft)

!===========================================================================
!==  First construct LDA green function (Init, Compute, Integrate, Print)
!===========================================================================
 write(message,'(6a)') ch10,' ========================', &
& ch10,' =====  LDA Calculation',&
& ch10,' ========================'
 call wrtout(std_out,message,'COLL')
 call icip_green("LDA",cryst_struc,greenlda,mpi_enreg,paw_dmft,pawang,1,self)

!== Compare greenlda%occup and lda_occup: check that LDA green function is fine
!----------------------------------------------------------------------
 write(message,'(2a)') ch10,&
& '  == Check lda occ. mat. from green with respect to the direct calc =='
 call wrtout(std_out,message,'COLL')
 call diff_oper("Occup from LDA green function",&
& "LDA occupations",greenlda%occup,lda_occup,1,tol4)
 write(message,'(2a)') ch10,&
& '  ***** => Calculation of Green function is thus correct without self ****'
 call wrtout(std_out,message,'COLL')
 call destroy_green(greenlda)
 
!== Orthonormalize psichi
!----------------------------------------------------------------------
 if(paw_dmft%dmft_solv/=-1) then
   call psichi_renormalization(cryst_struc,paw_dmft,pawang)

!  ===========================================================================
!  ==  re-construct LDA green function with new psichi's
!  ===========================================================================
   write(message,'(6a)')&
&   ch10,' ================================================'&
&   ,ch10,' =====  LDA Calculation with renormalized psichi'&
&   ,ch10,' ================================================'
 end if

 call wrtout(std_out,message,'COLL')
 call icip_green("LDA_renormalized",cryst_struc,greenlda,mpi_enreg,&
& paw_dmft,pawang,pawprtvol,self)
 call compute_energy(cryst_struc,energies_dmft,greenlda,paw_dmft,pawprtvol,pawtab,self,occ_type=" lda")
 if(paw_dmft%dmft_prgn==1) then
   if(paw_dmft%lpsichiortho==1) call local_ks_green(greenlda,paw_dmft,prtopt=1)
 end if

!== define Interaction from input upawu and jpawu
!----------------------------------------------------------------------
 ABI_ALLOCATE(hu,(cryst_struc%ntypat))
 call init_hu(cryst_struc,pawtab,hu)

!== define self from scratch or file and double counting
!----------------------------------------------------------------------
 call destroy_self(self)
!- Self allocated
 call initialize_self(self,paw_dmft)
 call dc_self(greenlda%charge_matlu,cryst_struc,hu,self,paw_dmft%dmft_dc,pawprtvol)
!- Read self or do self=hdc
 call rw_self(self,mpi_enreg,paw_dmft,prtopt=2,opt_rw=1,istep_iter=1000*istep)

 call destroy_green(greenlda)  ! destroy LDA green function
 call print_self(self,"print_dc",paw_dmft,prtopt=2)


!===========================================================================
!==  construct green function with the self-energy.
!===========================================================================
 write(message,'(6a)') &
& ch10,' =======================================================', &
& ch10,' =====  DMFT Calculation with input self-energy ========', &
& ch10,' ======================================================='
 call wrtout(std_out,message,'COLL')
 call icip_green("DMFT_inputself",cryst_struc,green,mpi_enreg,&
& paw_dmft,pawang,pawprtvol,self,opt_self=1)

!== Find fermi level
!---------------------------------------------------------------------
 write(message,'(2a,i3,13x,a)') ch10,'   ===  Compute green function from self-energy'
 call fermi_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,self)

!== define weiss field only for the local quantities (opt_oper=2)
!----------------------------------------------------------------------
 call init_green(weiss,paw_dmft,opt_oper_ksloc=2)

!== Check fourier transforms
!----------------------------------------------------------------------
 if(check==1) then
   call check_fourier_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang)
 end if

 write(message,'(6a)') &
& ch10,' ======================================================'&
& ,ch10,' =====  DMFT Loop starts here                  ========'&
& ,ch10,' ======================================================'
 call wrtout(std_out,message,'COLL')
!=======================================================================
!===  dmft loop  =======================================================
 do idmftloop=1, paw_dmft%dmft_iter 
   paw_dmft%idmftloop=idmftloop
!  =======================================================================
   istep_iter=1000*istep+idmftloop

   write(message,'(2a,i3,13x,a)') ch10,&
&   ' =====  DMFT Loop : ITER number',paw_dmft%idmftloop,'========'
   call wrtout(std_out,message,'COLL')

!  == Dyson Equation G,self -> weiss(w)
!  ---------------------------------------------------------------------
   call dyson(green,paw_dmft,self,weiss,opt_weissself=1)  

!  == Printout local "occupations" from weiss field  (useless)
   if(abs(pawprtvol)>3) then
     call integrate_green(cryst_struc,weiss,mpi_enreg,paw_dmft,&
&     pawang,prtopt=2,opt_ksloc=2)
     call printocc_green(weiss,5,paw_dmft,3,opt_weissgreen=1)
   end if
   
!  ===  Prepare data, solve Impurity problem: weiss(w) -> G(w)
!  ---------------------------------------------------------------------
   call initialize_self(self_new,paw_dmft)

   call impurity_solve(cryst_struc,green,hu,mpi_enreg,&
&   paw_dmft,pawang,pawtab,self,self_new,weiss,pawprtvol) ! weiss-> green, or self if dmft_solv=1

!  ==  Compute double counting from charge from green_solver
!  ---------------------------------------------------------------------
   if (green%has_charge_matlu_solver/=2) green%charge_matlu_solver=green%charge_matlu
   call dc_self(green%charge_matlu_solver,cryst_struc,hu,self_new,paw_dmft%dmft_dc,pawprtvol)

!  ==  Solve dyson equation. G_imp(w), weiss_imp(w) -> Self_imp(w)
!  ---------------------------------------------------------------------
!  if dmft_solv==1, self is computed previously
   if(abs(paw_dmft%dmft_solv)/=1) then
     call dyson(green,paw_dmft,self_new,weiss,opt_weissself=2) 
   end if

!  ==  Possibility if imposing self (opt_rw==3)
!  ---------------------------------------------------------------------
   call rw_self(self_new,mpi_enreg,paw_dmft,prtopt=2,opt_rw=3,istep_iter=istep_iter)

!  print dc just computed before and self computed in before in dyson or
!  impurity_solve
   if(abs(pawprtvol)>=3) then
     write(message,'(2a)') ch10,"  == New self" ; call wrtout(std_out,message,'COLL')
     call print_self(self_new,"print_dc",paw_dmft,2)
     write(message,'(2a)') ch10,"  == Old self" ; call wrtout(std_out,message,'COLL')
     call print_self(self,"print_dc",paw_dmft,2)
   end if

!  ==  Compute Energy with NEW self-energy and edc from green_solver
!  ---------------------------------------------------------------------
   call compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self_new,occ_type="nlda")

!  ==  Mix new and old self_energies and double countings
!  ---------------------------------------------------------------------
   call new_self(self,self_new,paw_dmft,1) ! self,self_new => self
   write(message,'(2a)') ch10,"  == After mixing,"
   call wrtout(std_out,message,'COLL')
   call print_self(self,"print_dc",paw_dmft,2) ! print self and DC
   call destroy_self(self_new)

!  ==  Compute green function self -> G(k) 
!  ---------------------------------------------------------------------
   call compute_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,1,self,opt_self=1)

   call integrate_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,prtopt=2,opt_ksloc=3,opt_diff=1)

   call printocc_green(green,5,paw_dmft,3,chtype="DMFT")
   if(paw_dmft%lpsichiortho==1.and.paw_dmft%dmft_prgn==1) call local_ks_green(green,paw_dmft,prtopt=1)

!  ==  Find fermi level
!  ---------------------------------------------------------------------
   call fermi_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,self)
!  call leave_new('COLL')

!  == Save self on disk
!  ---------------------------------------------------------------------
   call rw_self(self,mpi_enreg,paw_dmft,prtopt=2,opt_rw=2)

!  == Test convergency
!  ---------------------------------------------------------------------
   if(green%ifermie_cv==1.and.self%iself_cv==1.and.green%ichargeloc_cv==1.and.paw_dmft%idmftloop>1) then
     write(message,'(a,8x,a)') ch10,"DMFT Loop is converged !"
     call wrtout(std_out,message,'COLL')
     exit
   end if
!  =======================================================================
!  === end dmft loop  ====================================================
 end do
!=======================================================================

!== Save self on disk
!---------------------------------------------------------------------
 call rw_self(self,mpi_enreg,paw_dmft,prtopt=2,opt_rw=2)

 paw_dmft%idmftloop=0

 write(message,'(2a,13x,a)') ch10,' =====  DMFT Loop :  END          ',&
& '========'
!Do not compute here, because, one want a energy computed after the
!solver (for Hubbard I and LDA+U).
 call wrtout(std_out,message,'COLL')
 call compute_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,1,self,opt_self=1)
 call integrate_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,prtopt=2,opt_ksloc=3)
!call compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,opt=0)
!write(message,'(2a,13x,a)') ch10,' =====  DMFT Loop is finished'
!call wrtout(ab_out,message,'COLL')
 call printocc_green(green,9,paw_dmft,3,chtype="converged DMFT")
 if(paw_dmft%dmft_solv<=2.and.paw_dmft%prtdos>=1) then
   call spectral_function(cryst_struc,green,hu,mpi_enreg,&
&   paw_dmft,pawang,pawtab,self,pawprtvol) ! weiss-> green
 end if
 call destroy_green(weiss)
 call destroy_green(green)
!todo_ab rotate back density matrix into unnormalized basis just for
!printout 
 call destroy_hu(hu,cryst_struc%ntypat)
 call destroy_self(self)
 call destroy_energy(energies_dmft,paw_dmft)
 ABI_DEALLOCATE(hu)

 DBG_EXIT("COLL")

end subroutine dmft_solve
!!***
