!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_matlu
!! NAME
!!  m_matlu
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
!! NOTES
!!  subroutines in this module must never call
!!   a subroutine of m_oper, m_green, m_self
!!   in order to avoid circular dependancies
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_matlu

 use m_profiling

 use defs_basis
 use defs_datatypes

 implicit none

 private 

 public :: init_matlu
 public :: inverse_matlu
 public :: destroy_matlu
 public :: diff_matlu
 public :: add_matlu
 public :: nullify_matlu
 public :: print_matlu
 public :: sym_matlu
 public :: copy_matlu
 public :: gather_matlu
 public :: zero_matlu
 public :: trace_matlu
 public :: diag_matlu
 public :: rotate_matlu
!!***
 

!!****t* m_matlu/matlu_type
!! NAME
!!  matlu_type
!!
!! FUNCTION
!!  This structured datatype contains an matrix for the correlated subspace
!!
!! SOURCE

 type, public :: matlu_type ! for each atom

  integer :: lpawu         
  ! value of the angular momentum for correlated electrons

!  integer :: natom       
   ! number of atoms (given for each atom, not useful..could be changed)
!
!  integer :: mband
!  ! Number of bands
!      
!  integer :: mbandc
!  ! Total number of bands in the Kohn-Sham Basis for PAW+DMFT
!
!  integer :: natpawu         ! Number of correlated atoms 
!
!  integer :: nkpt
!  ! Number of k-point in the IBZ.
!
  integer :: nspinor
  ! Number of spinorial component
!
  integer :: nsppol
  ! Number of polarization 

  complex(dpc), pointer :: mat(:,:,:,:,:)
  ! local quantity
      
 end type matlu_type

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_matlu/init_matlu
!! NAME
!! init_matlu
!!
!! FUNCTION
!!  Allocate variables used in type matlu_type.
!!
!! INPUTS
!!  natom = number of atoms
!!  nspinor = number of spinorial components
!!  nsppol = number of polarisation components
!!  lpawu_natom(natom) = value of lpawu for all atoms
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! OUTPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! PARENTS
!!      datafordmft,hubbard_one,m_green,m_matlu,m_oper
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE

subroutine init_matlu(natom,nspinor,nsppol,lpawu_natom,matlu)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_matlu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: natom,nspinor,nsppol
!type
 integer, intent(in) :: lpawu_natom(natom)
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables ------------------------------------
 integer :: iatom,lpawu

!************************************************************************

! matlu%mband       = mband
! matlu%dmftbandf   = dmftbandf
! matlu%dmftbandi   = dmftbandi
! matlu%nkpt        = nkpt
! matlu%mbandc  = 0
 matlu%nsppol      = nsppol
 matlu%nspinor     = nspinor
 call nullify_matlu(matlu,natom)
 do iatom=1,natom
!  lpawu=pawtab(cryst_struc%typat(iatom))%lpawu
  lpawu=lpawu_natom(iatom)
  matlu(iatom)%lpawu=lpawu
  if(lpawu.ne.-1) then
   ABI_ALLOCATE(matlu(iatom)%mat,(2*lpawu+1,2*lpawu+1,nsppol,nspinor,nspinor))
   matlu(iatom)%mat=czero
  else
   ABI_ALLOCATE(matlu(iatom)%mat,(0,0,nsppol,nspinor,nspinor))
  endif
 enddo

end subroutine init_matlu
!!***

!!****f* m_matlu/nullify_matlu
!! NAME
!! nullify_matlu
!!
!! FUNCTION
!!  nullify matlu
!!
!! INPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! PARENTS
!!      hubbard_one,m_green,m_matlu,m_oper
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE

subroutine nullify_matlu(matlu,natom)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_matlu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(matlu_type),intent(inout) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom

!*********************************************************************

 do iatom=1,natom
  nullify(matlu(iatom)%mat)
 enddo


end subroutine nullify_matlu
!!***

!!****f* m_matlu/zero_matlu
!! NAME
!! zero_matlu
!!
!! FUNCTION
!!  zero_matlu 
!!
!! INPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! PARENTS
!!      m_green
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE

subroutine zero_matlu(matlu,natom)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'zero_matlu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(matlu_type),intent(inout) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom

!*********************************************************************

 do iatom=1,natom
  matlu(iatom)%mat=czero
 enddo


end subroutine zero_matlu
!!***

!!****f* m_matlu/destroy_matlu
!! NAME
!! destroy_matlu
!!
!! FUNCTION
!!  deallocate matlu
!!
!! INPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!
!! PARENTS
!!      datafordmft,hubbard_one,m_green,m_matlu,m_oper
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE

subroutine destroy_matlu(matlu,natom)

 use defs_basis
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_matlu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(matlu_type),intent(inout) :: matlu(natom)

!Local variables-------------------------------
 integer :: iatom

! *********************************************************************

 do iatom=1,natom
  if ( associated(matlu(iatom)%mat) )   then
    ABI_DEALLOCATE(matlu(iatom)%mat)
  end if
 enddo

end subroutine destroy_matlu
!!***

!!****f* m_matlu/copy_matlu
!! NAME
!! copy_matlu
!!
!! FUNCTION
!!  Copy matlu1 into matlu2 
!!
!! INPUTS
!!  maltu1 <type(matlu_type)>= density matrix matlu1 in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  maltu2 <type(matlu_type)>= density matrix matlu2 in the local orbital basis and related variables
!!
!! PARENTS
!!      datafordmft,hubbard_one,m_oper,m_self,spectral_function
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE

subroutine copy_matlu(nmat1,nmat2,natom)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_matlu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in) :: natom
 type(matlu_type),intent(in) :: nmat1(natom)
 type(matlu_type),intent(out) :: nmat2(natom)

!Local variables-------------------------------
 integer :: iatom, lpawu
! *********************************************************************

 do iatom=1,natom
  lpawu=nmat1(iatom)%lpawu
  if(lpawu.ne.-1) then
   nmat2(iatom)%mat=nmat1(iatom)%mat
  endif
 enddo


end subroutine copy_matlu
!!***

!!****f* m_matlu/print_matlu
!! NAME
!! print_matlu
!!
!! FUNCTION
!!
!! INPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom= number of atoms
!!  prtopt= option for printing
!!  opt_diag=   0   print non diagonal matrix (real or complex according to nspinor)
!!             -1   print non diagonal complex matrix
!!            >=1   print diagonal matrix (real or complex according to nspinor)
!!  opt_ab_out=  0  print matrix on std_out
!!             /=0  print matrix on ab_out
!!
!! OUTPUT
!!
!! PARENTS
!!      datafordmft,hubbard_one,m_green,m_matlu,m_oper,m_self
!!      psichi_renormalization
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE

subroutine print_matlu(matlu,natom,prtopt,opt_diag,opt_ab_out)

 use defs_basis
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_matlu'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in):: natom
 type(matlu_type),intent(in) :: matlu(natom)
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: opt_diag,opt_ab_out
!Local variables-------------------------------
 integer :: iatom,ispinor,ispinor1,isppol,m1,m,lpawu,nspinor,nsppol,optdiag,optab_out,arg_out
 character(len=500) :: message
! *********************************************************************
 if(present(opt_diag)) then
   optdiag=opt_diag
 else
   optdiag=0
 endif
 if(present(opt_ab_out)) then
   optab_out=opt_ab_out
 else
   optab_out=0
 endif
 if(optab_out==0) then
   arg_out=std_out
 else
   arg_out=ab_out
 endif
 nspinor=matlu(1)%nspinor
 nsppol=matlu(1)%nsppol

 do iatom = 1 , natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu/=-1) then
     write(message,'(2a,i4)')  ch10,'   -------> For Correlated Atom', iatom
     call wrtout(arg_out,  message,'COLL')
     do isppol = 1 , nsppol
       if(nspinor == 1) then
         write(message,'(a,10x,a,i3,i3)')  ch10,'-- polarization spin component',isppol
         call wrtout(arg_out,  message,'COLL')
       endif
       do ispinor = 1 , nspinor
         do ispinor1 = 1, nspinor
           if(nspinor == 2) then
             write(message,'(a,10x,a,i3,i3)')  ch10,'-- spin components',ispinor,ispinor1
             call wrtout(arg_out,  message,'COLL')
           endif
           if(optdiag<=0) then
             do m1=1, 2*lpawu+1
               if(optdiag==0) then
                 if(nspinor==1.and.abs(prtopt)>0) &
&                  write(message,'(5x,20f10.5)') (REAL(matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
                 if(nspinor==2.and.abs(prtopt)>0) &
&                  write(message,'(5x,14(2f9.5,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
               else if(optdiag==-1) then
                 write(message,'(5x,14(2f9.5,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
               endif
               call wrtout(arg_out,  message,'COLL')
             enddo
           elseif (optdiag>=1) then
             if(nspinor==1.and.abs(prtopt)>0) &
&             write(message,'(5x,20f10.5)') (REAL(matlu(iatom)%mat(m,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
             if(nspinor==2.and.abs(prtopt)>0) &
&             write(message,'(5x,14(2f9.5,2x))')((matlu(iatom)%mat(m,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
!            write(std_out,'(5x,14(2f9.5,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
             call wrtout(arg_out,  message,'COLL')
           endif
         end do ! ispinor1
       end do ! ispinor
     enddo ! isppol
   endif ! lpawu/=1
 enddo ! natom


end subroutine print_matlu
!!***

!{\src2tex{textfont=tt}}
!!****f* m_matlu/sym_matlu
!! NAME
!! sym_matlu
!!
!! FUNCTION
!! Symetrise local quantity.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_structure)>=crystal structure data
!!  gloc(natom) <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  gloc(natom) <type(matlu_type)>= density matrix symetrized in the local orbital basis and related variables
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      datafordmft,hubbard_one,m_green,psichi_renormalization
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
 subroutine sym_matlu(cryst_struc,gloc,pawang)
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sym_matlu'
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_structure),intent(in) :: cryst_struc
 type(pawang_type),intent(in) :: pawang
!arrays
 type(matlu_type),intent(inout) :: gloc(cryst_struc%natom)
!Local variables-------------------------------
!scalars
 integer :: at_indx,iatom,irot,ispinor,ispinor1,isppol,lpawu,m1,m2,m3,m4,mu
 integer :: natom,ndim,nsppol,nspinor,nu
 complex(dpc) :: sumrho,summag(3),rotmag(3),ci
 real(dp) :: zarot2
!arrays
! complex(dpc),allocatable :: glocnm(:,:,:,:,:)
 type(matlu_type),allocatable :: glocnm(:)
! complex(dpc),allocatable :: glocnms(:,:,:,:,:)
 type(matlu_type),allocatable :: glocnms(:)
 type(matlu_type),allocatable :: glocsym(:)
 real(dp),allocatable :: symrec_cart(:,:,:)

! DBG_ENTER("COLL")

 ci=cone
 nspinor=gloc(1)%nspinor
 nsppol=gloc(1)%nsppol
 natom=cryst_struc%natom

 ABI_ALLOCATE(glocnm,(natom))
 ABI_ALLOCATE(glocnms,(natom))
 ABI_ALLOCATE(glocsym,(natom))
 call init_matlu(natom,nspinor,nsppol,gloc(1:natom)%lpawu,glocnm)
 call init_matlu(natom,nspinor,nsppol,gloc(1:natom)%lpawu,glocnms)
 call init_matlu(natom,nspinor,nsppol,gloc(1:natom)%lpawu,glocsym)


!=========  Case nspinor ==1 ========================

 if (nspinor==1) then
  ispinor=1
  ispinor1=1
  do iatom=1,cryst_struc%natom
   do isppol=1,nsppol
    if(gloc(iatom)%lpawu/=-1) then
     lpawu=gloc(iatom)%lpawu
     do m1=1, 2*lpawu+1
      do m2=1, 2*lpawu+1
       do irot=1,cryst_struc%nsym
        at_indx=cryst_struc%indsym(4,irot,iatom)
        do m3=1, 2*lpawu+1
         do m4=1, 2*lpawu+1
          zarot2=pawang%zarot(m3,m1,lpawu+1,irot)*pawang%zarot(m4,m2,lpawu+1,irot) 
          glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)=&
&          glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)&
&          +gloc(at_indx)%mat(m3,m4,isppol,ispinor,ispinor1)*zarot2
         end do  ! m3
        end do  ! m4 
       end do  ! irot
       glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)=&
&       glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)/cryst_struc%nsym
      end do ! m2
     end do ! m1
    endif ! lpawu/=-1
   end do ! isppol
  end do ! iatom
  do iatom=1,cryst_struc%natom
    if(gloc(iatom)%lpawu/=-1) then
      gloc(iatom)%mat=glocsym(iatom)%mat
!!      gloc(iatom)%mat(:,:,1,:,:)=(glocsym(iatom)%mat(:,:,1,:,:) &
!!&      + glocsym(iatom)%mat(:,:,2,:,:))/two
!!      gloc(iatom)%mat(:,:,2,:,:)= gloc(iatom)%mat(:,:,1,:,:)
!!      write(std_out,*) "WARNING: SYM non mag"
!!      write(ab_out,*) "WARNING: SYM non mag"
    endif
  end do ! iatom

!=========  Case nspinor ==2 ========================

 else if (nspinor==2) then

!== Allocate temporary arrays
  do iatom=1,cryst_struc%natom
   if(gloc(iatom)%lpawu/=-1) then
    ndim=2*gloc(iatom)%lpawu+1
    ABI_DEALLOCATE(glocnm(iatom)%mat)
    ABI_DEALLOCATE(glocnms(iatom)%mat)
    ABI_DEALLOCATE(glocsym(iatom)%mat)
    ABI_ALLOCATE(glocnm(iatom)%mat,(ndim,ndim,nsppol,4,1))
    ABI_ALLOCATE(glocnms(iatom)%mat,(ndim,ndim,nsppol,4,1))
    ABI_ALLOCATE(glocsym(iatom)%mat,(ndim,ndim,nsppol,2,2))
   endif
  enddo
  ABI_ALLOCATE(symrec_cart,(3,3,cryst_struc%nsym))

!==  Compute symrec_cart
  do irot=1,cryst_struc%nsym
   call symredcart(cryst_struc%gprimd,cryst_struc%rprimd,symrec_cart(:,:,irot),cryst_struc%symrec(:,:,irot))
  end do

!==  Compute density matrix in density and magnetisation representation
  call chg_repr_matlu(gloc,glocnm,cryst_struc%natom,option=1,prtopt=1)

!==  Do the sum over symetrized density matrix (in n,m repr)
  isppol=1
  do iatom=1,cryst_struc%natom
   if(gloc(iatom)%lpawu/=-1) then
    lpawu=gloc(iatom)%lpawu
    ndim=2*gloc(iatom)%lpawu+1
    do m1=1, 2*lpawu+1
     do m2=1, 2*lpawu+1
      sumrho=czero
      rotmag=czero
      do irot=1,cryst_struc%nsym
       summag=czero
       at_indx=cryst_struc%indsym(4,irot,iatom)
       do m3=1, 2*lpawu+1
        do m4=1, 2*lpawu+1
         zarot2=pawang%zarot(m3,m2,lpawu+1,irot)*pawang%zarot(m4,m1,lpawu+1,irot) 
         sumrho=sumrho +  glocnm(at_indx)%mat(m4,m3,isppol,1,1)  * zarot2 
         do mu=1,3
          summag(mu)=summag(mu) + glocnm(at_indx)%mat(m4,m3,isppol,mu+1,1) * zarot2 
         enddo
        end do ! m3
       end do !m4

!       ==  special case of magnetization 
       do nu=1,3
        do mu=1,3
         rotmag(mu)=rotmag(mu)+symrec_cart(mu,nu,irot)*summag(nu) 
        end do
       end do
!      write(std_out,'(a,3i4,2x,3(2f10.5,2x))') "rotmag",irot,m1,m2,(rotmag(mu),mu=1,3)
      end do ! irot

!       ==  Normalizes sum
      sumrho=sumrho/cryst_struc%nsym 
!        sumrho=glocnm(isppol,1,iatom,m1,m2) ! test without sym
      glocnms(iatom)%mat(m1,m2,isppol,1,1)=sumrho
      do mu=1,3
       rotmag(mu)=rotmag(mu)/cryst_struc%nsym 
!          rotmag(mu)=glocnm(isppol,mu+1,iatom,m1,m2) ! test without sym
       glocnms(iatom)%mat(m1,m2,isppol,mu+1,1)=rotmag(mu)
      enddo
     end do  ! m2
    end do ! m1
   endif ! lpawu/=-1
  end do ! iatom

!==  Compute back density matrix in upup dndn updn dnup representation
  call chg_repr_matlu(glocsym,glocnms,cryst_struc%natom,option=-1,prtopt=1)

!==  Put glocsym into gloc
  do iatom=1,cryst_struc%natom
    if(gloc(iatom)%lpawu/=-1) then
      gloc(iatom)%mat=glocsym(iatom)%mat
    endif
  end do ! iatom

  ABI_DEALLOCATE(symrec_cart)
 endif

 call destroy_matlu(glocnm,cryst_struc%natom)
 call destroy_matlu(glocnms,cryst_struc%natom)
 call destroy_matlu(glocsym,cryst_struc%natom)
 ABI_DEALLOCATE(glocnm)
 ABI_DEALLOCATE(glocnms)
 ABI_DEALLOCATE(glocsym)
!==============end of nspinor==2 case ===========


! DBG_EXIT("COLL")

 end subroutine sym_matlu
!!***

!!****f* m_matlu/inverse_matlu
!! NAME
!! inverse_matlu
!!
!! FUNCTION
!! Inverse local quantity.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to inverse
!!  natom=number of atoms in cell.
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: inverse of input matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_oper
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
 subroutine inverse_matlu(matlu,natom,prtopt)
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inverse_matlu'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 integer, intent(in) :: prtopt
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom,tndim
 integer :: nsppol,nspinor
!scalars
 type(coeff2c_type),allocatable :: gathermatlu(:)
 !************************************************************************


 nspinor=matlu(1)%nspinor
 nsppol=matlu(1)%nsppol
 if(prtopt>0) then
 endif
 ABI_ALLOCATE(gathermatlu,(natom))
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=nsppol*nspinor*(2*matlu(iatom)%lpawu+1)
     ABI_ALLOCATE(gathermatlu(iatom)%value,(tndim,tndim))
     gathermatlu(iatom)%value=czero
   endif
 enddo

 call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=nsppol*nspinor*(2*matlu(iatom)%lpawu+1)
     call matcginv_dpc(gathermatlu(iatom)%value,tndim,tndim)
   endif
 enddo
 call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)

 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     ABI_DEALLOCATE(gathermatlu(iatom)%value)
   endif
 enddo
 ABI_DEALLOCATE(gathermatlu)
 end subroutine inverse_matlu
!!***

!!****f* m_matlu/diff_matlu
!! NAME
!! diff_matlu
!!
!! FUNCTION
!!
!! INPUTS
!!  char1 = character describing matlu1
!!  char2 = character describing matlu2
!!  matlu1(natom) <type(matlu_type)>= density matrix 1 in the local orbital basis and related variables
!!  matlu2(natom) <type(matlu_type)>= density matrix 2 in the local orbital basis and related variables
!!  natom = number of atoms
!!  option =1      if diff> toldiff , stop 
!!          0      print diff and toldiff
!!          else   do not test and do not print
!!  toldiff = maximum value for the difference between matlu1 and matlu2
!!
!! OUTPUT
!!
!! PARENTS
!!      datafordmft,m_green,m_oper
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
subroutine diff_matlu(char1,char2,matlu1,matlu2,natom,option,toldiff,ierr)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'diff_matlu'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 integer,intent(in) :: natom,option
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 character(len=*), intent(in) :: char1,char2
 real(dp),intent(in) :: toldiff
 integer,intent(out), optional :: ierr

!local variables-------------------------------
 integer :: iatom,idiff,ispinor,ispinor1,isppol,m1,m,lpawu,nspinor,nsppol
 real(dp) :: matludiff
 character(len=500) :: message
! *********************************************************************
 nsppol=matlu1(1)%nsppol
 nspinor=matlu1(1)%nspinor

 matludiff=zero
 idiff=0
 do iatom = 1 , natom
  lpawu=matlu1(iatom)%lpawu
  if(lpawu/=-1) then
   do isppol = 1 , nsppol
    do ispinor = 1 , nspinor
     do ispinor1 = 1, nspinor
      do m1 = 1 , 2*lpawu+1
       do m = 1 ,  2*lpawu+1
        idiff=idiff+1
        matludiff=matludiff+ &
&        sqrt( real(matlu1(iatom)%mat(m1,m,isppol,ispinor,ispinor1)      &
&            -      matlu2(iatom)%mat(m1,m,isppol,ispinor,ispinor1))**2  &
&            + imag(matlu1(iatom)%mat(m1,m,isppol,ispinor,ispinor1)      &
&            -      matlu2(iatom)%mat(m1,m,isppol,ispinor,ispinor1))**2  )
!       write(std_out,*) m,m1,matlu1(iatom)%mat(m1,m,isppol,ispinor,ispinor1),matlu2(iatom)%mat(m1,m,isppol,ispinor,ispinor1),matludiff
       enddo
      enddo
     end do ! ispinor1
    end do ! ispinor
   enddo ! isppol
  endif ! lpawu/=1
 enddo ! natom
 matludiff=matludiff/float(idiff)
 if(option==1.or.option==0) then
  if( matludiff < toldiff ) then
   write(message,'(5a,6x,3a,4x,e12.4,a,e12.4)') ch10,&
&   '   ** Differences between ',trim(char1),' and ',ch10,trim(char2),' are small enough:',&
&   ch10,matludiff,' is lower than',toldiff
   call wrtout(std_out,message,'COLL')
   if(present(ierr)) ierr=0
  else
   write(message,'(5a,3x,3a,3x,e12.4,a,e12.4)') ch10,&
&   '   Warning : Differences between ',trim(char1),' and ',ch10,trim(char2),' is too large:',&
&   ch10,matludiff,' is larger than',toldiff
   call wrtout(std_out,message,'COLL')
!   write(message,'(8a,4x,e12.4,a,e12.4)') ch10,"  Matrix for ",trim(char1)
   write(message,'(a,3x,a)') ch10,trim(char1)
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu1,natom,prtopt=1,opt_diag=-1)
   write(message,'(a,3x,a)') ch10,trim(char2)
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu2,natom,prtopt=1,opt_diag=-1)
   if(option==1) call leave_new('COLL')
   if(present(ierr)) ierr=-1
  endif
 endif
 
end subroutine diff_matlu
!!***

!!****f* m_matlu/add_matlu
!! NAME
!! add_matlu
!!
!! FUNCTION
!!
!! INPUTS
!!  maltu1 <type(matlu_type)>= density matrix matlu1 in the local orbital basis and related variables
!!  maltu2 <type(matlu_type)>= density matrix matlu2 in the local orbital basis and related variables
!!  natom = number of atoms
!!  sign_matlu2= 1 add matlu1 and matlu2
!!              -1 substract matlu2 to matlu1
!!
!! OUTPUT
!!  maltu3 <type(matlu_type)>= density matrix matlu3, sum/substract matlu1 and matlu2
!!
!! PARENTS
!!      dyson,m_green
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
subroutine add_matlu(matlu1,matlu2,matlu3,natom,sign_matlu2)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'add_matlu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 integer,intent(in) :: natom,sign_matlu2
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 type(matlu_type), intent(out) :: matlu3(natom)

!local variables-------------------------------
 integer :: iatom
! *********************************************************************

 do iatom = 1 , natom
   matlu3(iatom)%mat=matlu1(iatom)%mat+float(sign_matlu2)*matlu2(iatom)%mat
 enddo ! natom
 
end subroutine add_matlu
!!***

!!****f* m_matlu/chg_repr_matlu
!! NAME
!! chg_repr_matlu
!!
!! FUNCTION
!! Change representation of density matrix (useful for nspinor=2)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  glocspsp(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: density matrix in the spin spin representation
!!  glocnm(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: density matrix in the magnetization representation
!!  natom=number of atoms in cell.
!!  option= 1 glocspsp is input, glocnm is computed
!!  option= -1 glocspsp is computed, glocnm is input
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  glocspsp(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: density matrix in the spin spin representation
!!  glocnm(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: density matrix in the magnetization representation
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_matlu
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
 subroutine chg_repr_matlu(glocspsp,glocnm,natom,option,prtopt)
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chg_repr_matlu'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option,prtopt
!arrays
 type(matlu_type),intent(inout) :: glocspsp(natom)
 type(matlu_type),intent(inout) :: glocnm(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,isppol,lpawu,m1,m2,ndim,nsppol,mu
 complex(dpc) :: ci
 character(len=500) :: message

! DBG_ENTER("COLL")

 ci=j_dpc

!==  Compute density matrix in density magnetisation representation
 if (option==1) then
  nsppol=glocspsp(1)%nsppol
  do isppol=1,nsppol
   do iatom=1,natom
    if(glocspsp(iatom)%lpawu/=-1) then
     ndim=2*glocspsp(iatom)%lpawu+1
     do m1=1,ndim
      do m2=1,ndim
       glocnm(iatom)%mat(m1,m2,isppol,1,1)=glocspsp(iatom)%mat(m1,m2,isppol,1,1) &
&                                         +glocspsp(iatom)%mat(m1,m2,isppol,2,2)
       glocnm(iatom)%mat(m1,m2,isppol,4,1)=glocspsp(iatom)%mat(m1,m2,isppol,1,1) &
&                                         -glocspsp(iatom)%mat(m1,m2,isppol,2,2)
       glocnm(iatom)%mat(m1,m2,isppol,2,1)=glocspsp(iatom)%mat(m1,m2,isppol,1,2) &
&                                         +glocspsp(iatom)%mat(m1,m2,isppol,2,1)
       glocnm(iatom)%mat(m1,m2,isppol,3,1)= &
&                               cmplx((imag(glocspsp(iatom)%mat(m1,m2,isppol,2,1))   &
&                                     -imag(glocspsp(iatom)%mat(m1,m2,isppol,1,2))), &
&                                    (-real(glocspsp(iatom)%mat(m1,m2,isppol,2,1))+  &
&                                      real(glocspsp(iatom)%mat(m1,m2,isppol,1,2))),kind=dp)
      enddo  ! m2
     enddo ! m1
     if(abs(prtopt)>=3) then
      write(message,'(a)') "        -- in n, m repr "
      call wrtout(std_out,  message,'COLL')
      do mu=1,4
       do m1=1,ndim
        write(message,'(8x,(14(2f9.5,2x)))')(glocnm(iatom)%mat(m1,m2,isppol,mu,1),m2=1,ndim)
        call wrtout(std_out,  message,'COLL')
       enddo ! m1
       write(message,'(a)') ch10
       call wrtout(std_out,  message,'COLL')
      enddo ! mu
     endif ! prtopt >3
    endif ! lpawu/=-1
   enddo
  enddo

!==  Compute back density matrix in upup dndn updn dnup representation
 else  if (option==-1) then
  isppol=1
  do iatom=1,natom
   if(glocnm(iatom)%lpawu/=-1) then
    lpawu=glocnm(iatom)%lpawu
    ndim=2*glocnm(iatom)%lpawu+1
    do m1=1, 2*lpawu+1
     do m2=1, 2*lpawu+1
      glocspsp(iatom)%mat(m1,m2,isppol,1,1)=half*(glocnm(iatom)%mat(m1,m2,isppol,1,1)+glocnm(iatom)%mat(m1,m2,isppol,4,1))
      glocspsp(iatom)%mat(m1,m2,isppol,2,2)=half*(glocnm(iatom)%mat(m1,m2,isppol,1,1)-glocnm(iatom)%mat(m1,m2,isppol,4,1))
      glocspsp(iatom)%mat(m1,m2,isppol,1,2)=half*(glocnm(iatom)%mat(m1,m2,isppol,2,1)-ci*glocnm(iatom)%mat(m1,m2,isppol,3,1))
      glocspsp(iatom)%mat(m1,m2,isppol,2,1)=half*(glocnm(iatom)%mat(m1,m2,isppol,2,1)+ci*glocnm(iatom)%mat(m1,m2,isppol,3,1))
     end do  ! m2
    end do ! m1
    if(abs(prtopt)>6) then
     write(message,'(a)') "        -- in spin spin repr "
     call wrtout(std_out,  message,'COLL')
     do mu=1,4
      do m1=1,ndim
       write(message,'(8x,14(2f9.5,2x))')(glocspsp(iatom)%mat(m1,m2,isppol,mu,1),m2=1,ndim)
       call wrtout(std_out,  message,'COLL')
      enddo
      write(message,'(a)') ch10
      call wrtout(std_out,  message,'COLL')
     enddo
    endif !prtopt>3
   endif ! lpawu/=-1
  end do ! iatom
 else  
  write(message,'(2a,i4)')  ch10,"stop in chg_repr_matlu"
  call wrtout(std_out,  message,'COLL')
  call leave_new('COLL')
 endif


! DBG_EXIT("COLL")

 end subroutine chg_repr_matlu
!!***

!!****f* m_matlu/trace_matlu
!! NAME
!! trace_matlu
!!
!! FUNCTION
!!  Compute the trace of the matlu matrix
!!
!! INPUTS
!!  maltu(natom) <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  trace_loc(natom,nsppol+1)= trace for each atoms and each polarization, trace_loc(iatom,nsppol+1) is
!!                             the full trace over polarization also.
!!
!! PARENTS
!!      m_oper
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
 subroutine trace_matlu(matlu,natom,trace_loc)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'trace_matlu'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in) :: natom
 type(matlu_type), intent(in) :: matlu(natom)
 real(dp),intent(inout) :: trace_loc(natom,matlu(1)%nsppol+1)

!local variables-------------------------------
 integer :: iatom,isppol,ispinor,m,lpawu
 integer :: nsppol,nspinor
 character(len=500) :: message
! *********************************************************************
 nsppol=matlu(1)%nsppol
 nspinor=matlu(1)%nspinor

 trace_loc=zero
 do iatom = 1 , natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu/=-1) then
     write(message,'(2a,i4)')  ch10,'   -------> For Correlated Atom', iatom
     call wrtout(std_out,  message,'COLL')
     do isppol = 1 , nsppol
       do ispinor = 1 , nspinor
         do m = 1 ,  2*lpawu+1
           trace_loc(iatom,isppol)=trace_loc(iatom,isppol)+&
&           matlu(iatom)%mat(m,m,isppol,ispinor,ispinor)
         enddo 
       enddo 
      trace_loc(iatom,nsppol+1)=trace_loc(iatom,nsppol+1)+trace_loc(iatom,isppol)
     enddo 
     write(message,'(8x,a,f12.6)')   'Nb of Corr. elec. is:'&
&     ,trace_loc(iatom,nsppol+1)
     call wrtout(std_out,  message,'COLL')
   endif
 enddo

 end subroutine trace_matlu
!!***

!!****f* m_matlu/gather_matlu
!! NAME
!! gather_matlu
!!
!! FUNCTION
!! Create new array from matlu
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gloc(natom) <type(matlu_type)>        = density matrix in the spin spin representation
!!  gatherloc(natom) <type(coeff2c_type)> = density matrix where spin and angular momentum are gathered in the same index
!!  natom=number of atoms in cell.
!!  option= 1 go from gloc to gathergloc
!!  option= -1 go from gathergloc to gloc
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  gloc(natom) <type(matlu_type)>        = density matrix in the spin spin representation
!!  gatherloc(natom) <type(coeff2c_type)> = density matrix where spin and angular momentum are gathered in the same index
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_matlu,psichi_renormalization
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
 subroutine gather_matlu(gloc,gathergloc,natom,option,prtopt)
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gather_matlu'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

! type  matlus_type 
!  SEQUENCE
!  complex(dpc), pointer :: mat(:,:)
! end type matlus_type

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option,prtopt
 type(coeff2c_type), intent(inout) :: gathergloc(natom)
 type(matlu_type),intent(inout) :: gloc(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,ispinor,ispinor1,isppol,isppol1
 integer :: jc1,jc2,ml1,ml2,ndim,nspinor,nsppol,tndim
 complex(dpc) :: ci
 character(len=500) :: message

! DBG_ENTER("COLL")
 nsppol=gloc(1)%nsppol
 nspinor=gloc(1)%nspinor

 ci=j_dpc


 do iatom=1,natom
   if(gloc(iatom)%lpawu.ne.-1) then
!==-------------------------------------

     ndim=2*gloc(iatom)%lpawu+1
     tndim=nsppol*nspinor*ndim

!== Put norm into array "gathergloc"
     jc1=0
     do isppol=1,nsppol
       do ispinor=1,nspinor
         do ml1=1,ndim
           jc1=jc1+1
           jc2=0
           do isppol1=1,nsppol
             do ispinor1=1,nspinor
               do ml2=1,ndim
                 jc2=jc2+1
                 if(option==1) then
                   if(isppol==isppol1) then
                     gathergloc(iatom)%value(jc1,jc2)=gloc(iatom)%mat(ml1,ml2,isppol,ispinor,ispinor1)
                   endif
                 else if(option==-1) then
                   if(isppol==isppol1) then
                     gloc(iatom)%mat(ml1,ml2,isppol,ispinor,ispinor1)=gathergloc(iatom)%value(jc1,jc2)
                   endif
                 endif
               enddo
             enddo ! ispinor1
           enddo ! isppol1
         enddo 
       enddo !ispinor
     enddo ! isppol
   endif
 enddo ! iatom
 if(option==1.and.prtopt==3) then
   do iatom=1,natom
     if(gloc(iatom)%lpawu.ne.-1) then
       tndim=nsppol*nspinor*(2*gloc(iatom)%lpawu+1)
       write(message,'(2a,i5)') ch10,' (gathermatlu:) For atom', iatom
       call wrtout(std_out,message,'COLL')
       do im1=1,tndim
         write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&         (gathergloc(iatom)%value(im1,im2),im2=1,tndim)
         call wrtout(std_out,message,'COLL')
       end do
     endif
   enddo ! iatom
 else if(option==-1.and.prtopt==3) then
   call print_matlu(gloc,natom,prtopt)
 endif



! DBG_EXIT("COLL")

 end subroutine gather_matlu
!!***

!!****f* m_matlu/diag_matlu
!! NAME
!! diag_matlu
!!
!! FUNCTION
!! Diagonalize matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to diagonalize
!!  natom=number of atoms 
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  matlu_diag(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: diagonalized density matrix
!!  eigvectmatlu(natom) <type(coeff2c_type)> = Eigenvectors corresponding to the diagonalization
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      hubbard_one
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
 subroutine diag_matlu(matlu,matlu_diag,natom,prtopt,eigvectmatlu)
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'diag_matlu'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 integer, intent(in) :: prtopt
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
 type(matlu_type),intent(out) :: matlu_diag(natom)
 type(coeff2c_type),optional,intent(out) :: eigvectmatlu(natom,matlu(1)%nsppol)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,info,ispinor,isppol,lwork,tndim
 integer :: nsppol,nspinor
 character(len=500) :: message
!arrays
 type(coeff2c_type),allocatable :: gathermatlu(:)
 real(dp),allocatable :: eig(:),rwork(:)
 complex(dpc),allocatable :: zwork(:)
!************************************************************************

 nsppol=matlu(1)%nsppol
 nspinor=matlu(1)%nspinor
 do isppol=1,nsppol
! ===========================
! Define gathermatlu
! ===========================
   ABI_ALLOCATE(gathermatlu,(natom))
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       tndim=nspinor*(2*matlu(iatom)%lpawu+1)
       ABI_ALLOCATE(gathermatlu(iatom)%value,(tndim,tndim))
       gathermatlu(iatom)%value=czero
     endif
   enddo
   if(nsppol==1.and.nspinor==2) then
     call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
   else if(nsppol==2.and.nspinor==1) then
     do iatom=1,natom
       if(matlu(iatom)%lpawu.ne.-1) then
         do im1=1,tndim
           do im2=1,tndim
             gathermatlu(iatom)%value(im1,im2)=matlu(iatom)%mat(im1,im2,isppol,1,1)
           enddo
         enddo
       endif
     enddo
   endif
! ===========================
! Diagonalize
! ===========================
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       tndim=nspinor*(2*matlu(iatom)%lpawu+1)
       lwork=2*tndim-1
       ABI_ALLOCATE(rwork,(3*tndim-2))
       ABI_ALLOCATE(zwork,(lwork))
       ABI_ALLOCATE(eig,(tndim))
       if(prtopt>=3) then
         do im1=1,tndim
           write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&           (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
           call wrtout(std_out,message,'COLL')
         end do
       endif
!         write(std_out,*)"diag"
       call zheev('v','u',tndim,gathermatlu(iatom)%value,tndim,eig,zwork,lwork,rwork,info)
       if(prtopt>=3) then
         do im1=1,tndim
           write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&           (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
           call wrtout(std_out,message,'COLL')
         end do
       endif
!       write(std_out,*) "eig",eig
! ===========================
! Put eigenvalue in matlu_diag
! ===========================
       do ispinor=1,nspinor
         do im1=1,2*matlu(iatom)%lpawu+1
           matlu_diag(iatom)%mat(im1,im1,isppol,ispinor,ispinor)=eig(im1)
         enddo
       enddo
       ABI_DEALLOCATE(zwork)
       ABI_DEALLOCATE(rwork)
       ABI_DEALLOCATE(eig)
!     endif
!   enddo
! ===========================
! Keep eigenvectors gathermatlu
! ===========================
       if (present(eigvectmatlu)) then
         tndim=nspinor*(2*matlu(iatom)%lpawu+1)
         eigvectmatlu(iatom,isppol)%value(:,:)=gathermatlu(iatom)%value(:,:) 
       endif
     endif
   enddo
! End loop over atoms
! ===========================
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       ABI_DEALLOCATE(gathermatlu(iatom)%value)
     endif
   enddo
   ABI_DEALLOCATE(gathermatlu)
 enddo ! isppol

 end subroutine diag_matlu
!!***

!!****f* m_matlu/rotate_matlu
!! NAME
!! rotate_matlu
!!
!! FUNCTION
!! Rotate matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to rotate
!!  rot_mat(natom) <type(coeff2c_type)> = Rotation matrix (usually from diag_matlu)
!!  natom=number of atoms in cell.
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: Rotated matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      hubbard_one
!!
!! CHILDREN
!!      gather_matlu,zgemm
!!
!! SOURCE
 subroutine rotate_matlu(matlu,rot_mat,natom,prtopt)
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_matlu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 integer, intent(in) :: prtopt
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
 type(coeff2c_type),optional,intent(out) :: rot_mat(natom,matlu(1)%nsppol)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,isppol
 integer :: nsppol,nspinor,tndim
! character(len=500) :: message
!arrays
 type(coeff2c_type),allocatable :: gathermatlu(:)
 complex(dpc),allocatable :: temp_mat(:,:)
!************************************************************************
 if(prtopt==1) then
 endif
 nsppol=matlu(1)%nsppol
 nspinor=matlu(1)%nspinor
 do isppol=1,nsppol
! ===========================
! Define gathermatlu and allocate
! ===========================
   ABI_ALLOCATE(gathermatlu,(natom))
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       tndim=nspinor*(2*matlu(iatom)%lpawu+1)
       ABI_ALLOCATE(gathermatlu(iatom)%value,(tndim,tndim))
       gathermatlu(iatom)%value=czero
     endif
   enddo
   if(nsppol==1.and.nspinor==2) then
     call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
   else if(nsppol==2.and.nspinor==1) then
     do iatom=1,natom
       if(matlu(iatom)%lpawu.ne.-1) then
         do im1=1,tndim
           do im2=1,tndim
             gathermatlu(iatom)%value(im1,im2)=matlu(iatom)%mat(im1,im2,isppol,1,1)
           enddo
         enddo
       endif
     enddo
   endif
! ===========================
! Rotate
! ===========================
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       tndim=nspinor*(2*matlu(iatom)%lpawu+1)
       ABI_ALLOCATE(temp_mat,(tndim,tndim))
       temp_mat(:,:)=czero
       call zgemm('n','c',tndim,tndim,tndim,cone,gathermatlu(iatom)%value   ,tndim,&
&        rot_mat(iatom,isppol)%value,tndim,czero,temp_mat                ,tndim)
       call zgemm('n','n',tndim,tndim,tndim,cone,rot_mat(iatom,isppol)%value,tndim,&
&        temp_mat                   ,tndim,czero,gathermatlu(iatom)%value,tndim)
       ABI_DEALLOCATE(temp_mat)
     endif ! lpawu
   enddo ! iatom

! ===========================
! Put data into matlu(iatom)
! ===========================
   if(nsppol==1.and.nspinor==2) then
     call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)
   else if(nsppol==2.and.nspinor==1) then
     do iatom=1,natom
       if(matlu(iatom)%lpawu.ne.-1) then
         do im1=1,tndim
           do im2=1,tndim
             matlu(iatom)%mat(im1,im2,isppol,1,1)= gathermatlu(iatom)%value(im1,im2)
           enddo
         enddo
       endif
     enddo
   endif ! test nsppol/nspinor
! ===========================
! Deallocations
! ===========================
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       ABI_DEALLOCATE(gathermatlu(iatom)%value)
     endif
   enddo
   ABI_DEALLOCATE(gathermatlu)
 enddo ! isppol

 end subroutine rotate_matlu


END MODULE m_matlu
!!***
