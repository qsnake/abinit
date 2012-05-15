!{\src2tex{textfont=tt}}
!!****f* ABINIT/impurity_solve
!! NAME
!! impurity_solve
!!
!! FUNCTION
!! Solve the Impurity problem
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc
!!  istep    =  step of iteration for LDA.
!!  lda_occup
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab <type(pawtab)>
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      copy_green,destroy_green_tau,fourier_green,hubbard_one,init_green_tau
!!      integrate_green,ldau_self,print_green,printocc_green,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine impurity_solve(cryst_struc,green,hu,mpi_enreg,paw_dmft,&
& pawang,pawtab,self_old,self_new,weiss,pawprtvol)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_crystal, only : crystal_structure
 use m_green, only : green_type, fourier_green&
& ,init_green_tau,destroy_green_tau,print_green,printocc_green,integrate_green,copy_green
 use m_paw_dmft, only : paw_dmft_type
 use m_hu, only : hu_type
 use m_self, only : self_type
 use m_energy, only : energy_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'impurity_solve'
 use interfaces_14_hidewrite
 use interfaces_68_dmft, except_this_one => impurity_solve
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_structure),intent(in) :: cryst_struc
 type(green_type), intent(inout) :: weiss
 type(green_type), intent(out) :: green
 type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
 type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(self_type), intent(inout) :: self_new
 type(self_type), intent(inout) :: self_old
 integer, intent(in) :: pawprtvol

!Local variables ------------------------------
 character(len=500) :: message
! integer iatom,il,i_nd,isppol,lpawu,im,Nd,nrat,nsweeptot
! real(dp) :: acc,kx
! real(dp), allocatable :: correl(:,:),g0(:,:),gtmp(:,:)
!scalars
!************************************************************************
!character(len=500) :: message

!=======================================================================
!== Prepare data for solver QMC
!=======================================================================
 if(abs(paw_dmft%dmft_solv)==4) then
!  == Initialize weiss and green functions for fourier transformation
!  -------------------------------------------------------------------
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Initialize Weiss field G_0(tau)'
   call wrtout(std_out,message,'COLL')
   call init_green_tau(weiss,paw_dmft)
   call init_green_tau(green,paw_dmft)
!  in init_solver

!  == Print weiss function G_0(tau=0-) before computation (really useless check)
!  ------------------------------------------------------------------------------
   if(abs(pawprtvol)>3) then
     write(message,'(2a,i3,13x,a)') ch10,'   ===  Check G_0(tau=0-) first'
     call wrtout(std_out,message,'COLL')
     call printocc_green(weiss,6,paw_dmft,3)
   end if

!  == Fourier transform of weiss Field
!  ------------------------------------
!  for fourier of KS green functions
!  call fourier_green(cryst_struc,weiss,mpi_enreg,paw_dmft,pawang,pawtab,1)
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Inverse Fourier Transform w->t of Weiss Field'
   call wrtout(std_out,message,'COLL')
   call fourier_green(cryst_struc,weiss,mpi_enreg,paw_dmft,pawang,opt_ksloc=2,opt_tw=-1)

!  == Print weiss function G_0(tau=0-) 
!  --------------------------------------
   call printocc_green(weiss,6,paw_dmft,3,opt_weissgreen=1)

!  for fourier of KS green functions
!  call fourier_green(cryst_struc,weiss,mpi_enreg,paw_dmft,pawang,pawtab,1)
!  == Print G_0(tau) in files
!  ---------------------------
   if(paw_dmft%dmft_prgn==2) then
     call print_green('weiss',weiss,1,paw_dmft,pawprtvol=1,opt_wt=2)
   end if
 end if
!=======================================================================
!== End preparation of QMC
!=======================================================================

!=======================================================================
!== Solve impurity model   =============================================
!=======================================================================
 write(message,'(2a,i3,13x,a)') ch10,'  ===  Solve impurity model'
 call wrtout(std_out,message,'COLL')
 if(abs(paw_dmft%dmft_solv)==1) then

!  == LDA+U for test -> self
!  -------------------
   call ldau_self(cryst_struc,green,paw_dmft,&
&   pawtab,self_new,opt_ldau=1,prtopt=pawprtvol)
 else if(abs(paw_dmft%dmft_solv)==2) then

!  == Hubbard One -> green
!  -------------------
   call hubbard_one(cryst_struc,green,hu,paw_dmft,&
&   pawang,pawtab,pawprtvol,self_old%hdc,weiss)

 else if(abs(paw_dmft%dmft_solv)==4) then

!  == Nothing
!  -------------------
   call copy_green(weiss,green,opt_tw=1)
!  !  todo:  diagonalize green function
!  !!   call diag_matlu(weiss%oper_tau(1)%matlu,weiss_diag,natom,&
!  !&   prtopt=2,eigvectmatlu=eigvectmatlu)
!  ! take xmu and selfDCqmc into account
!  nsweeptot=10000
!  do iatom=1,cryst_struc%natom
!  lpawu=pawtab(cryst_struc%typat(iatom))%lpawu
!  if(lpawu/=-1) then
!  Nd=2*lpawu+1
!  allocate(correl(Nd,Nd))
!  allocate(g0(paw_dmft%dmftqmc_l,Nd))
!  allocate(gtmp(paw_dmft%dmftqmc_l,Nd))
!  do im=1,lpawu
!  do isppol=1,paw_dmft%nsppol
!  i_nd=im+(isppol-1)*lpawu 
!  do il=1,paw_dmft%dmftqmc_l
!  g0(il,i_nd)=weiss%oper_tau(il)%matlu(iatom)%mat(im,im,isppol,1,1)
!  enddo
!  enddo
!  enddo
!  !      call qmc(g0,gtmp,dmftqmc_l,1,1,dmft_idmftloop,acc,nrat,correl,kx,Nd,hu,nsweeptot) 
!  do im=1,lpawu
!  do isppol=1,paw_dmft%nsppol
!  i_nd=im+(isppol-1)*lpawu 
!  do il=1,paw_dmft%dmftqmc_l
!  green%oper_tau(il)%matlu(iatom)%mat(im,im,isppol,1,1)=gtmp(il,i_nd)
!  enddo
!  enddo
!  enddo
!  
!  deallocate(correl,g0,gtmp)
!  endif
!  enddo
!  !   do il=1,paw_dmft%dmftqmc_l
!  !    call rotate_matlu(green%oper_tau(il)%matlu,eigvectmatlu,natom,3)
!  !   enddo

 else if(abs(paw_dmft%dmft_solv)==0) then

!  == Nothing
!  -------------------
!  weiss%occup%has_operks=0 -> only local part is duplicated
   call copy_green(weiss,green,opt_tw=2)
 end if
!call print_green("invWeiss",cryst_struc,weiss,3,paw_dmft,pawtab,2)

!=======================================================================
!== Treat data from QMC
!=======================================================================
 if(abs(paw_dmft%dmft_solv)==4) then

   if(paw_dmft%dmft_prgn==1) then
     call print_green('DMFT',green,1,paw_dmft,pawprtvol=1,opt_wt=2)
   end if

!  == Print local occupations from G(tau)
!  ---------------------------------------
   call printocc_green(green,6,paw_dmft,3)

!  == Fourier back transform of green function G(tau)->G(iw_n)
!  -------------------------------------------------------------------
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Direct Fourier Transform t->w of Weiss Field'
   call wrtout(std_out,message,'COLL')
   call fourier_green(cryst_struc,green,mpi_enreg,paw_dmft,&
&   pawang,opt_ksloc=2,opt_tw=1)

!  == Back fourier transform of G_0(tau) for compensation (try).
!  -------------------------------------------------------------------
   call fourier_green(cryst_struc,weiss,mpi_enreg,paw_dmft,&
&   pawang,opt_ksloc=2,opt_tw=1)
   call destroy_green_tau(weiss)
   call destroy_green_tau(green)
 end if
!=======================================================================
!== End Treat data for QMC
!=======================================================================

!=======================================================================
!== Integrate green function and printout occupations
!For dmft_solv=-1,0,or 1 , the green function was not computed: it
!cannot be integrated
!=======================================================================
 if(paw_dmft%dmft_solv>=2.and.green%w_type=="imag") then
!  ==  Integrate G(iw_n)
!  ---------------------
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Integrate local part of green function'
   call wrtout(std_out,message,'COLL')
   call integrate_green(cryst_struc,green,mpi_enreg,paw_dmft,&
&   pawang,prtopt=2,opt_ksloc=2,opt_after_solver=1)

!  == Print local occupations from integration of G(iw_n)
!  --------------------------------------------------------
   call printocc_green(green,5,paw_dmft,3)

!  == Print G_loc(w)
!  --------------------------------------------------------
   if(paw_dmft%dmft_prgn==1) then
     call print_green('DMFT',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
   end if
 end if

 if(abs(pawprtvol)>0) then
 end if

end subroutine impurity_solve
!!***
