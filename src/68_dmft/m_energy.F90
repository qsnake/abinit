!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_energy
!! NAME
!!  m_energy
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2012 ABINIT group (BAmadon)
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

MODULE m_energy

 use m_profiling

 use defs_basis
 use defs_datatypes

 implicit none

 private 

 public :: init_energy
 public :: compute_energy
 public :: compute_ldau_energy
 public :: destroy_energy
 public :: nullify_energy
 public :: print_energy

!!***

!!****t* m_energy/energy_type
!! NAME
!!  energy_type
!!
!! FUNCTION
!!  This structured datatype contains interaction matrices for the correlated subspace
!!
!! SOURCE

 type, public :: energy_type ! for each typat

  real(dp) :: eband_lda

  real(dp) :: eband_dmft

  real(dp) :: e_dc_tot
      
  real(dp) :: e_hu_tot

  real(dp) :: e_hu_ldau_tot

  real(dp) :: e_hu_mig_tot

  real(dp) :: edmft

  real(dp) :: natom

  real(dp), pointer :: e_dc(:) ! => vee

  real(dp), pointer :: e_hu(:)

  real(dp), pointer :: e_hu_ldau(:)

  real(dp), pointer :: e_hu_mig(:)

 end type energy_type

!!***

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_energy/init_energy
!! NAME
!! init_energy
!!
!! FUNCTION
!!  Allocate variables used in type energy_type.
!!
!! INPUTS
!!
!! OUTPUTS
!! energies_dmft  = structure of data for dmft of type energy_type
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      destroy_paw_ij,init_paw_ij,nullify_paw_ij,pawuenergy,wrtout
!!
!! SOURCE

subroutine init_energy(cryst_struc,energies_dmft)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_energy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(crystal_structure),intent(in) :: cryst_struc
 type(energy_type), intent(inout) :: energies_dmft
!Local variables ------------------------------------
!************************************************************************

 ABI_ALLOCATE(energies_dmft%e_dc,(cryst_struc%natom))
 ABI_ALLOCATE(energies_dmft%e_hu,(cryst_struc%natom))
 ABI_ALLOCATE(energies_dmft%e_hu_ldau,(cryst_struc%natom))
 ABI_ALLOCATE(energies_dmft%e_hu_mig,(cryst_struc%natom))
 energies_dmft%e_dc=zero
 energies_dmft%e_hu=zero
 energies_dmft%e_hu_ldau=zero
 energies_dmft%e_hu_mig=zero
 energies_dmft%eband_lda=zero
 energies_dmft%eband_dmft=zero
 energies_dmft%e_dc_tot=zero
 energies_dmft%e_hu_tot=zero
 energies_dmft%e_hu_ldau_tot=zero
 energies_dmft%e_hu_mig_tot=zero
 energies_dmft%edmft=zero
 energies_dmft%natom=cryst_struc%natom

end subroutine init_energy
!!***

!!****f* m_energy/nullify_energy
!! NAME
!! nullify_energy
!!
!! FUNCTION
!!  nullify energies_dmft
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_paw_ij,init_paw_ij,nullify_paw_ij,pawuenergy,wrtout
!!
!! SOURCE

subroutine nullify_energy(energies_dmft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_energy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(energy_type),intent(inout) :: energies_dmft
!Local variables-------------------------------

!*********************************************************************

 nullify(energies_dmft%e_dc)
 nullify(energies_dmft%e_hu)
 nullify(energies_dmft%e_hu_ldau)
 nullify(energies_dmft%e_hu_mig)


end subroutine nullify_energy
!!***

!!****f* m_energy/destroy_energy
!! NAME
!! destroy_energy
!!
!! FUNCTION
!!  deallocate energies_dmft
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      destroy_paw_ij,init_paw_ij,nullify_paw_ij,pawuenergy,wrtout
!!
!! SOURCE

subroutine destroy_energy(energies_dmft,paw_dmft)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_energy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(energy_type),intent(inout) :: energies_dmft
 type(paw_dmft_type), intent(out) :: paw_dmft
!Local variables-------------------------------
! *********************************************************************
  paw_dmft%edmft=energies_dmft%edmft
 if ( associated(energies_dmft%e_dc) )   then
   ABI_DEALLOCATE(energies_dmft%e_dc)
 end if
 if ( associated(energies_dmft%e_hu) )   then
   ABI_DEALLOCATE(energies_dmft%e_hu)
 end if
 if ( associated(energies_dmft%e_hu_ldau) )  then
   ABI_DEALLOCATE(energies_dmft%e_hu_ldau)
 end if
 if ( associated(energies_dmft%e_hu_mig) )  then
   ABI_DEALLOCATE(energies_dmft%e_hu_mig)
 end if

end subroutine destroy_energy
!!***

!!****f* m_energy/print_energy
!! NAME
!! print_energy
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_energy
!!
!! CHILDREN
!!      destroy_paw_ij,init_paw_ij,nullify_paw_ij,pawuenergy,wrtout
!!
!! SOURCE

subroutine print_energy(cryst_struc,energies_dmft,pawprtvol,pawtab)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_energy'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(crystal_structure),intent(in) :: cryst_struc
 type(energy_type),intent(in) :: energies_dmft
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 integer, intent(in) :: pawprtvol

!Local variables-------------------------------
 integer :: iatom
 character(len=500) :: message
! *********************************************************************
 if(abs(pawprtvol)>=3) then
   do iatom=1,cryst_struc%natom
     if(pawtab(cryst_struc%typat(iatom))%lpawu/=-1) then
       write(message,'(a,4x,a,i3,a,5x,f12.6)')  &
&       ch10,"For Correlated Atom",iatom,", E_hu =",energies_dmft%e_hu(iatom)
       call wrtout(std_out,message,'COLL')
       write(message,'(26x,a,1x,f12.6)')  &
&       ", E_hu_mig =",energies_dmft%e_hu_mig(iatom)
       call wrtout(std_out,message,'COLL')
       write(message,'(26x,a,f12.6)')  &
&       ", E_hu_ldau =",energies_dmft%e_hu_ldau(iatom)
       call wrtout(std_out,message,'COLL')
       write(message,'(26x,a,f12.6)')  &
&       ", E_dc =",energies_dmft%e_dc(iatom)
       call wrtout(std_out,message,'COLL')
     endif
   enddo
 endif
 write(message,'(a,5x,2a,5x,a,7(a,5x,a,2x,f16.9),a,5x,a)') ch10 &
&      ,"-----------------------------------------------",ch10 &
&      ,"--- Energy in DMFT (in Ha)  ",ch10 &
&      ,"--- E_bandlda (1)  (Ha.) = ",energies_dmft%eband_lda,ch10 &
&      ,"--- E_banddmft(2)  (Ha.) = ",energies_dmft%eband_dmft,ch10 &
&      ,"--- E_hu      (3)  (Ha.) = ",energies_dmft%e_hu_tot,ch10 &
&      ,"--- E_hu_mig  (4)  (Ha.) = ",energies_dmft%e_hu_mig_tot,ch10 &
&      ,"--- E_hu_ldau (5)  (Ha.) = ",energies_dmft%e_hu_ldau_tot,ch10 &
&      ,"--- E_dc      (6)  (Ha.) = ",energies_dmft%e_dc_tot,ch10 &
&      ,"--- edmft=(    3-6)(Ha.) = ",energies_dmft%edmft,ch10 &
!&      ,"--- edmft2=(2-1+3-6)(Ha.)= ",zero,ch10 &
&      ,"-----------------------------------------------" 
 call wrtout(std_out,message,'COLL')
end subroutine print_energy
!!***

!!****f* m_energy/compute_energy
!! NAME
!! compute_energy
!!
!! FUNCTION
!!
!! INPUTS
!!  cryst_struc <type(crystal_structure)>=crystal structure data
!!  green  <type(green_type)>= green function data 
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      destroy_paw_ij,init_paw_ij,nullify_paw_ij,pawuenergy,wrtout
!!
!! SOURCE

subroutine compute_energy(cryst_struc,energies_dmft,green,paw_dmft,pawprtvol,pawtab,self,occ_type)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_crystal, only : crystal_structure
 use m_green, only : green_type
 use m_self, only : self_type
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_energy'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(energy_type),intent(inout) :: energies_dmft
 type(crystal_structure),intent(in) :: cryst_struc
 type(green_type),intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(self_type), intent(in) :: self
 integer, intent(in) :: pawprtvol
 character(len=4), intent(in) :: occ_type
! integer :: prtopt

!Local variables-------------------------------
 integer :: iatom,ib,ifreq,ikpt,im,im1,ispinor,ispinor1,isppol,lpawu
 integer :: natom,ndim,nspinor,nsppol,nwlo
 real(dp) :: beta
 complex(dpc) :: xmig_1,xmig_2,xmig_3,se,gr
 character(len=500) :: message
! *********************************************************************
 write(message,'(2a)') ch10,"  == Compute LDA+DMFT energy terms "
 call wrtout(std_out,message,'COLL')

! Only imaginary frequencies here
 if(green%w_type=="real".or.self%w_type=="real") then
   write(message,'(2a,i3,13x,a)') ch10,'   BUG: compute_energy not implemented for real frequency'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 endif

! == Compute band energy
! -----------------------------------------------------------------------
 energies_dmft%eband_lda=zero
 energies_dmft%eband_dmft=zero
 do isppol=1,paw_dmft%nsppol
   do ikpt=1,paw_dmft%nkpt
     do ib=1,paw_dmft%mbandc
       energies_dmft%eband_dmft=energies_dmft%eband_dmft+ &
&         green%occup%ks(isppol,ikpt,ib,ib)*&
&         paw_dmft%eigen_lda(isppol,ikpt,ib)*paw_dmft%wtk(ikpt)
       energies_dmft%eband_lda=energies_dmft%eband_lda+ &
&         occup_fd(paw_dmft%eigen_lda(isppol,ikpt,ib),paw_dmft%fermie_lda,paw_dmft%temp)*&
&         paw_dmft%eigen_lda(isppol,ikpt,ib)*paw_dmft%wtk(ikpt)
!          write(std_out,*) "isppol,ikpt,ib",isppol,ikpt,ib
!          write(std_out,*) "paw_dmft%eigen_lda",paw_dmft%eigen_lda(isppol,ikpt,ib)
!          write(std_out,*) green%occup%ks(isppol,ikpt,ib,ib)
!          write(std_out,*) occup_fd(paw_dmft%eigen_lda(isppol,ikpt,ib),paw_dmft%fermie,paw_dmft%temp)
     enddo
   enddo
 enddo

 if (occ_type==" lda") then
   if(abs(energies_dmft%eband_lda-energies_dmft%eband_dmft)>tol5) then
     write(message,'(5x,a,a,a,15x,a,f12.6,a,15x,a,5x,f12.5)')  "Warning !:"&
&     ,"Differences between band energy from LDA occupations",ch10&
&     ,"and LDA green function is:",energies_dmft%eband_lda-energies_dmft%eband_dmft,ch10&
&     ,"which is larger than",tol5
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') &
&     "   Action: increase number of frequencies, or reduce the number of high energies_dmft bands"
     call wrtout(std_out,message,'COLL')
   else
     write(message,'(a,a,a,10x,a,f12.6,a,10x,a,5x,f12.5)')  "          "&
&     ,"Differences between band energy from LDA occupations",ch10&
&     ,"and LDA green function is:",energies_dmft%eband_lda-energies_dmft%eband_dmft,ch10&
&     ,"which is smaller than",tol5
     call wrtout(std_out,message,'COLL')
   endif
 endif

! == Compute Correlation energy from Migdal formula
! -----------------------------------------------------------------------
 natom=cryst_struc%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 beta=one/paw_dmft%temp
 nwlo=paw_dmft%dmft_nwlo
! write(std_out,*) "beta",beta

 xmig_1=zero
 xmig_2=zero
 xmig_3=zero
 energies_dmft%e_hu_mig_tot = zero
 do iatom=1,natom
   lpawu=pawtab(cryst_struc%typat(iatom))%lpawu
   if(lpawu/=-1) then
     xmig_1=czero
     xmig_2=czero
     xmig_3=czero
     ndim=2*lpawu+1
     do isppol=1,nsppol
       do ispinor = 1 , nspinor
         do ispinor1 = 1, nspinor
           do im=1,ndim
             do im1=1,ndim
               do ifreq=1,nwlo
!                write(std_out,*) ifreq,xmig_1,imag(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)),&
!&                  green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )
                 xmig_1=xmig_1 + j_dpc/beta*       &
&                 imag(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))* &
&                      green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )* &
&                      paw_dmft%wgt_wlo(ifreq)
                 if(ispinor==ispinor1.and.im==im1) then
                   se=(self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)-  &
&                      self%oper (nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))
                 else
                   se=self%oper (ifreq)%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1)
                 endif          
                 xmig_2=xmig_2 + one/beta*real(se)* &
&                      green%oper(ifreq)%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor )* &
&                      paw_dmft%wgt_wlo(ifreq)
                 if(ispinor==ispinor1.and.im==im1.and.ifreq==1) then
                   xmig_3=xmig_3 + &
&                   real(self%oper(nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))* &
&                         green%occup%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor)/two
                 endif
               enddo
!               if(ispinor==ispinor1.and.im==im1) then
!                 xmig_3=xmig_3 + &
!&                 real(self%oper(nwlo )%matlu(iatom)%mat(im ,im1,isppol,ispinor ,ispinor1))* &
!!&                         green%occup%matlu(iatom)%mat(im1,im ,isppol,ispinor1,ispinor)/two
!               endif
             enddo
           enddo
         enddo
       enddo
     enddo
     energies_dmft%e_hu_mig(iatom)=real(xmig_1+xmig_2+xmig_3)
     energies_dmft%e_hu_mig_tot = energies_dmft%e_hu_mig_tot + energies_dmft%e_hu_mig(iatom)
     if(abs(pawprtvol)>=3) then
       write(message,'(2a,3(a,5x,a,2f12.6))')ch10,&
&         "  Interaction energy: Decomposition of Migdal energy",ch10,&
&         "xmig_1=",xmig_1,ch10,&
&         "xmig_3=",xmig_2,ch10,&
&         "xmig_3=",xmig_3
       call wrtout(std_out,message,'COLL')
     endif
   endif
 enddo


 call compute_ldau_energy(cryst_struc,energies_dmft,green,paw_dmft,pawtab)
 if(abs(paw_dmft%dmft_solv)<=1) then
   energies_dmft%e_hu= energies_dmft%e_hu_ldau
   energies_dmft%e_hu_tot= energies_dmft%e_hu_ldau_tot
 else if(paw_dmft%dmft_solv==2) then
   energies_dmft%e_hu= energies_dmft%e_hu_mig
   energies_dmft%e_hu_tot= energies_dmft%e_hu_mig_tot
 else if(paw_dmft%dmft_solv==4) then
   write(message,'(2a)') ch10,"Warning, energy is not computed"
   call wrtout(std_out,message,'COLL')
 endif
 energies_dmft%edmft=energies_dmft%e_hu_mig_tot-energies_dmft%e_dc_tot

 call print_energy(cryst_struc,energies_dmft,pawprtvol,pawtab)
! write(message,'(2a)') ch10," == The LDA+U self-energy is == "
! call wrtout(std_out,message,'COLL')
! call print_oper(self%oper(1),5,paw_dmft,2)
! a voir: energies_dmft%e_hu_tot = energies_dmft%e_hu_ldau_tot

end subroutine compute_energy
!!***

!!****f* m_energy/compute_ldau_energy
!! NAME
!! compute_ldau_energy
!!
!! FUNCTION
!!  Initialize noccmmp from green%occup and compute LDA+U energy with it
!!
!! INPUTS
!!  cryst_struc <type(crystal_structure)>=crystal structure data
!!  green  <type(green_type)>= green function data 
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!
!! PARENTS
!!      m_energy
!!
!! CHILDREN
!!      destroy_paw_ij,init_paw_ij,nullify_paw_ij,pawuenergy,wrtout
!!
!! SOURCE

subroutine compute_ldau_energy(cryst_struc,energies_dmft,green,paw_dmft,pawtab)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_crystal, only : crystal_structure
 use m_green, only : green_type
 use m_paw_dmft, only : paw_dmft_type
 use m_paw_toolbox, only : init_paw_ij,destroy_paw_ij,nullify_paw_ij

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_ldau_energy'
 use interfaces_14_hidewrite
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(energy_type),intent(inout) :: energies_dmft
 type(crystal_structure),intent(in) :: cryst_struc
 type(green_type),intent(in) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
! integer :: prtopt

!Local variables-------------------------------
 integer :: iatom,idijeff,im,im1,ispinor,ispinor1,isppol,ldim,lpawu
 integer :: nsploop
 character(len=500) :: message
 real(dp) :: eldaumdcdc,eldaumdc,e_ee,e_dc,e_dcdc,xe1,xe2
! arrays
 type(paw_ij_type), allocatable :: paw_ij(:)
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
! *********************************************************************

! - allocations
! -----------------------------------------------------------------------
 ABI_ALLOCATE(paw_ij,(cryst_struc%natom))
 call nullify_paw_ij(paw_ij)
!Should be contained in one of the arguments
 call init_paw_ij(paw_ij,2,2,paw_dmft%nspinor,paw_dmft%nsppol,paw_dmft%nspden, &
&      1,paw_dmft%natom,cryst_struc%ntypat,cryst_struc%typat,pawtab,has_pawu_occ=1)
 nsploop=max(paw_dmft%nsppol,paw_dmft%nspinor**2)
 e_ee=zero
 e_dc=zero
 e_dcdc=zero
 isppol=0
 ispinor=0
 ispinor1=0

! - Loop and call to pawuenergy
! -----------------------------------------------------------------------
 do iatom=1,cryst_struc%natom
   lpawu=pawtab(cryst_struc%typat(iatom))%lpawu
   if(lpawu.ne.-1) then
     ldim=2*lpawu+1
! - Setup nocctot and noccmmp
! -----------------------------------------------------------------------
     paw_ij(iatom)%nocctot(:)=zero ! contains nmmp in the n m representation
! Begin loop over spin/spinors to initialize noccmmp
     do idijeff=1,nsploop
       if(nsploop==2) then
         isppol=spinor_idxs(1,idijeff)
         ispinor=1
         ispinor1=1
       else if(nsploop==4) then
         isppol=1
         ispinor=spinor_idxs(1,idijeff)
         ispinor1=spinor_idxs(2,idijeff)
       else
         write(message,'(2a)') " BUG in ldau_self: nsploop should be equal to 2 or 4"
         call wrtout(std_out,message,'COLL')
       endif
! Initialize noccmmp
       do im1 = 1 , ldim
         do im = 1 ,  ldim
            paw_ij(iatom)%noccmmp(1,im,im1,idijeff)=&
&             real(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
            paw_ij(iatom)%noccmmp(2,im,im1,idijeff)=&
&             imag(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
         enddo
       enddo
! Compute nocctot 
       if(green%has_charge_matlu_solver/=2) then
         do im1=1,ldim
           if(nsploop==4) then
             paw_ij(iatom)%nocctot(idijeff)=paw_ij(iatom)%nocctot(idijeff)+&
&              paw_ij(iatom)%noccmmp(1,im1,im1,idijeff)
           else
             paw_ij(iatom)%nocctot(idijeff)=paw_ij(iatom)%nocctot(idijeff)+&
&              paw_ij(iatom)%noccmmp(1,im1,im1,idijeff)
           end if
         enddo
       else
         if(nsploop==4) then
           paw_ij(iatom)%nocctot(1)=green%charge_matlu_solver(iatom,2) !  total nb of elec for nspinor=2 is (iatom,2) !!
           paw_ij(iatom)%nocctot(2)=zero
           paw_ij(iatom)%nocctot(3)=zero
           paw_ij(iatom)%nocctot(4)=zero
         else
           paw_ij(iatom)%nocctot(1)=green%charge_matlu_solver(iatom,1) !  first spin
           paw_ij(iatom)%nocctot(2)=green%charge_matlu_solver(iatom,2) !  second one
         end if
       endif
     enddo
     paw_ij(iatom)%has_pawu_occ=2
     xe1=e_dc
     xe2=e_ee
     call pawuenergy(iatom,eldaumdc,eldaumdcdc,1,pawtab(cryst_struc%typat(iatom)),&
&     paw_ij(iatom),e_ee,e_dc,e_dcdc,paw_dmft%dmft_dc)
     energies_dmft%e_dc(iatom)=e_dc-xe1
     energies_dmft%e_hu_ldau(iatom)=e_ee-xe2 
   endif ! lpawu/=-1
 enddo

! - gather results
! -----------------------------------------------------------------------
 energies_dmft%e_dc_tot=e_dc ! todo_ab: here or not ?
 energies_dmft%e_hu_ldau_tot=e_ee

! - dealloc
! -----------------------------------------------------------------------
 call destroy_paw_ij(paw_ij)
 ABI_DEALLOCATE(paw_ij)


end subroutine compute_ldau_energy
!!***
!!****f* m_energy/occup_fd
!! NAME
!! occup_fd
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!     
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

 function occup_fd(eig,fermie,temp)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'occup_fd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
! Integrate analytic tail 1/(iw-mu)
 real(dp),intent(in) :: eig,fermie,temp
 real(dp) :: occup_fd
!Local variables-------------------------------
! *********************************************************************

 if((eig-fermie) > zero) then
   occup_fd=exp(-(eig-fermie)/temp)/(one+exp(-(eig-fermie)/temp))
 else
   occup_fd=one/(one+exp((eig-fermie)/temp))
 endif

 end function occup_fd

END MODULE m_energy
!!***
!!***
