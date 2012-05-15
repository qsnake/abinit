!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_green
!! NAME
!!  m_green
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_green

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi

 use m_oper, only : oper_type
 use m_matlu, only : matlu_type

 implicit none

 private

 public :: init_green
 public :: destroy_green
 public :: init_green_tau
 public :: destroy_green_tau
 public :: nullify_green
 public :: print_green
 public :: printocc_green
 public :: compute_green
 public :: integrate_green
 public :: icip_green
 public :: fourier_green
 public :: check_fourier_green
 public :: copy_green
 public :: int_fct
 public :: fourier_fct
!!***

!!****t* m_green/green_type
!! NAME
!!  green_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! SOURCE

 type, public :: green_type ! for each atom

  integer :: dmft_nwlo
  ! dmft frequencies

  character(len=4) :: w_type
  ! type of frequencies used

  integer :: nw
  ! dmft frequencies

  integer :: dmft_nwli
  ! dmft frequencies

  integer :: dmftqmc_l
  ! number of time slices for QMC

  integer :: ifermie_cv

  integer :: ichargeloc_cv

  real(dp) :: charge_ks
  ! Total charge computed from ks orbitals

  integer :: has_charge_matlu_solver
  ! =0 charge_matlu_solver not allocated
  ! =1 charge_matlu_solver is allocated
  ! =2 charge_matlu_solver is calculated (ie calculation of LOCAL CORRELATED occupations is done  from
  ! solver green function)

  integer :: has_charge_matlu
  ! =2 if calculation of LOCAL CORRELATED occupations is done  from

  integer :: has_charge_matlu_prev
  ! =0 charge_matlu_prev not allocated
  ! =1 charge_matlu_prev is allocated
  ! =2 charge_matlu_prev is calculated (ie calculation of LOCAL CORRELATED occupations is done  from
  ! solver green function)

  real(dp), pointer :: charge_matlu_prev(:,:)
  ! Total charge on correlated orbitals from previous iteration

  real(dp), pointer :: charge_matlu(:,:)
  ! Total charge on correlated orbitals
! todo_ba name of charge_matlu is misleading: should be changed

  real(dp), pointer :: charge_matlu_solver(:,:)
  ! Total charge on correlated orbitals obtained from solver

  real(dp), pointer :: tau(:)
  ! value of time in imaginary space

  real(dp), pointer :: omega(:)
  ! value of frequencies

  type(oper_type), pointer :: oper(:)
  ! green function  in different basis

  type(oper_type), pointer :: oper_tau(:)
  ! green function  in different basis

  type(oper_type) :: occup
  ! occupation in different basis

  type(oper_type) :: occup_tau
  ! occupation in different basis

 end type green_type

!----------------------------------------------------------------------


CONTAINS
!========================================================================================
!!***

!!****f* m_green/init_green
!! NAME
!! init_green
!!
!! FUNCTION
!!  Allocate variables used in type green_type.
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  opt_oper_ksloc (optional) = option for init_oper
!!  wtype = "real" Green function will be computed for real frequencies
!!        = "imag" Green function will be computed for imaginary frequencies
!!
!! OUTPUTS
!! green  = variable of type green_type
!!
!! PARENTS
!!      dmft_solve,dyson,hubbard_one,m_green,spectral_function
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine init_green(green,paw_dmft,opt_oper_ksloc,wtype)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_crystal, only : crystal_structure
 use m_oper, only : init_oper,nullify_oper
 use m_paw_dmft, only: paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_green'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!type
 type(green_type), intent(out) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) :: opt_oper_ksloc
 character(len=4), optional :: wtype
!local variables ------------------------------------
 integer :: ifreq,optoper_ksloc

!************************************************************************
 if(present(opt_oper_ksloc)) then
   optoper_ksloc=opt_oper_ksloc
 else
   optoper_ksloc=3
 endif
 if(present(wtype)) then
   green%w_type=wtype
 else
   green%w_type="imag"
 endif

 if(green%w_type=="imag") then
   green%nw=paw_dmft%dmft_nwlo
   green%omega=>paw_dmft%omega_lo
 else if(green%w_type=="real") then
   green%nw=2*paw_dmft%dmft_nwr
   green%omega=>paw_dmft%omega_r
 endif

 green%dmft_nwlo=paw_dmft%dmft_nwlo
 green%dmft_nwli=paw_dmft%dmft_nwli
 green%charge_ks=zero
 green%has_charge_matlu=0
 nullify(green%charge_matlu)
 ABI_ALLOCATE(green%charge_matlu,(paw_dmft%natom,paw_dmft%nsppol+1))
 green%charge_matlu=zero
 green%has_charge_matlu=1

 green%has_charge_matlu_solver=0
 nullify(green%charge_matlu_solver)
 ABI_ALLOCATE(green%charge_matlu_solver,(paw_dmft%natom,paw_dmft%nsppol+1))
 green%charge_matlu_solver=zero
 green%has_charge_matlu_solver=1

 green%has_charge_matlu_prev=0
 nullify(green%charge_matlu_prev)
 ABI_ALLOCATE(green%charge_matlu_prev,(paw_dmft%natom,paw_dmft%nsppol+1))
 green%charge_matlu_prev=zero
 green%has_charge_matlu_prev=1

 call nullify_oper(green%occup)
 call init_oper(paw_dmft,green%occup,opt_ksloc=optoper_ksloc)

 nullify(green%oper)
 ABI_ALLOCATE(green%oper,(green%nw))
 do ifreq=1,green%nw
  call nullify_oper(green%oper(ifreq))
  call init_oper(paw_dmft,green%oper(ifreq),opt_ksloc=optoper_ksloc)
 enddo
 green%ifermie_cv=0
 green%ichargeloc_cv=0

end subroutine init_green
!!***

!!****f* m_green/init_green_tau
!! NAME
!! init_green_tau
!!
!! FUNCTION
!!  Allocate variables used in type green_type.
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUTS
!!  green  <type(green_type)>= green function data
!!
!! PARENTS
!!      impurity_solve,m_green
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine init_green_tau(green,paw_dmft)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_oper, only : init_oper,nullify_oper
 use m_paw_dmft, only: paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_green_tau'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!type
 type(green_type), intent(inout) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
!local variables ------------------------------------
 integer :: itau

!************************************************************************

 green%dmftqmc_l=paw_dmft%dmftqmc_l
 nullify(green%tau)
 ABI_ALLOCATE(green%tau,(green%dmftqmc_l))
 do itau=1,green%dmftqmc_l
  green%tau(itau)=float(itau-1)/float(green%dmftqmc_l)/paw_dmft%temp
 enddo

 call nullify_oper(green%occup_tau)
 call init_oper(paw_dmft,green%occup_tau)

 nullify(green%oper_tau)
 ABI_ALLOCATE(green%oper_tau,(paw_dmft%dmftqmc_l))
 do itau=1,green%dmftqmc_l
  call nullify_oper(green%oper_tau(itau))
  call init_oper(paw_dmft,green%oper_tau(itau))
 enddo
end subroutine init_green_tau
!!***

!!****f* m_green/nullify_green
!! NAME
!! nullify_green
!!
!! FUNCTION
!!  nullify green
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!  green  <type(green_type)>= green function data
!!
!! PARENTS
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine nullify_green(green)

 use defs_basis
 use m_oper, only : nullify_oper

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_green'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(green_type),intent(inout) :: green
!*********************************************************************

 nullify(green%oper)


end subroutine nullify_green
!!***

!!****f* m_green/destroy_green
!! NAME
!! destroy_green
!!
!! FUNCTION
!!  deallocate green
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!
!! PARENTS
!!      dmft_solve,dyson,hubbard_one,m_green,spectral_function
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine destroy_green(green)

 use defs_basis
 use m_oper, only : destroy_oper,nullify_oper

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_green'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(green_type),intent(inout) :: green

!local variables-------------------------------
 integer :: ifreq

! *********************************************************************
 call destroy_oper(green%occup)
! call nullify_oper(green%occup)
 do ifreq=1,green%nw
  call destroy_oper(green%oper(ifreq))
!  call nullify_oper(green%oper(ifreq))
 enddo
 if ( associated(green%oper))       then
   ABI_DEALLOCATE(green%oper)
 end if
 if ( associated(green%charge_matlu))       then
   ABI_DEALLOCATE(green%charge_matlu)
 end if
 green%has_charge_matlu=0

 if ( associated(green%charge_matlu_prev))       then
   ABI_DEALLOCATE(green%charge_matlu_prev)
 end if
 green%has_charge_matlu_prev=0

 if ( associated(green%charge_matlu_solver))       then
   ABI_DEALLOCATE(green%charge_matlu_solver)
 end if
 green%has_charge_matlu_solver=0
! call nullify_green(green)

end subroutine destroy_green
!!***

!!****f* m_green/destroy_green_tau
!! NAME
!! destroy_green_tau
!!
!! FUNCTION
!!  deallocate green
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!
!! PARENTS
!!      impurity_solve,m_green
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine destroy_green_tau(green)

 use defs_basis
 use m_crystal, only : crystal_structure
 use m_oper, only : destroy_oper,nullify_oper

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_green_tau'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(green_type),intent(inout) :: green

!local variables-------------------------------
 integer :: itau

! *********************************************************************
 call destroy_oper(green%occup_tau)
! call nullify_oper(green%occup_tau)
 do itau=1,green%dmftqmc_l
  call destroy_oper(green%oper_tau(itau))
!  call nullify_oper(green%oper_tau(itau))
 enddo
 if ( associated(green%oper_tau))       then
   ABI_DEALLOCATE(green%oper_tau)
 end if
 if ( associated(green%tau))            then
   ABI_DEALLOCATE(green%tau)
 end if
! nullify(green%oper_tau)
! nullify(green%tau)

end subroutine destroy_green_tau
!!***

!!****f* m_green/copy_green
!! NAME
!! copy_green
!!
!! FUNCTION
!!  copy one data structure green1 into green2
!!
!! INPUTS
!!  green1  <type(green_type)>= green function data
!!  green2  <type(green_type)>= green function data
!!  opt_tw = option to precise which data to copy
!!          1: copy only green%occup_tau and green%oper_tau data
!!          2: copy only green%occup and green%oper  data (frequency)
!!
!! OUTPUT
!!
!! PARENTS
!!      dyson,impurity_solve,m_green,spectral_function
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine copy_green(green1,green2,opt_tw)

 use defs_basis
 use m_crystal, only : crystal_structure
 use m_oper, only : copy_oper

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_green'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(green_type),intent(in) :: green1
 type(green_type),intent(out) :: green2
 integer, intent(in) :: opt_tw

!local variables-------------------------------
 integer :: ifreq,itau
! *********************************************************************

 if(opt_tw==2) then
   call copy_oper(green1%occup,green2%occup)
   do ifreq=1, green1%nw
     call copy_oper(green1%oper(ifreq),green2%oper(ifreq))
   enddo
 else if (opt_tw==1) then
   call copy_oper(green1%occup_tau,green2%occup_tau)
   do itau=1,green1%dmftqmc_l
     call copy_oper(green1%oper_tau(itau),green2%oper_tau(itau))
   enddo
 endif

end subroutine copy_green
!!***

!!****f* m_green/printocc_green
!! NAME
!! printocc_green
!!
!! FUNCTION
!!  print occupations
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!! option= 1 :for G(w)
!!         2 :for G(tau)
!!         3 :for G(tau) and check % G(w)
!!         4
!!         <5: write diagonal part of KS occupation matrix
!!         5: for G(w)
!!         6: for G(tau)
!!         7 :for G(tau) and check % G(w)
!!         >8: write all elements of KS occup. matrix.
!!         9: for G(w)
!!
!! OUTPUT
!!
!! PARENTS
!!      dmft_solve,impurity_solve,m_green
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine printocc_green(green,option,paw_dmft,pawprtvol,opt_weissgreen,chtype)

 use defs_basis
 use m_crystal, only : crystal_structure
 use m_oper, only : print_oper
 use m_matlu, only : diff_matlu,print_matlu
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'printocc_green'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type),intent(in) :: green
 integer,intent(in) :: option,pawprtvol
 integer,optional,intent(in) :: opt_weissgreen
 character(len=*), optional, intent(in) :: chtype

!local variables-------------------------------
 character(len=500) :: message
 integer :: optweissgreen
! *********************************************************************
 if(present(opt_weissgreen)) then
   optweissgreen=opt_weissgreen
 else
   optweissgreen=2
 endif


 if(mod(option,4)==1) then
   if(optweissgreen==2) then
     if(present(chtype)) then
       write(message,'(4a)') ch10,"  == The ",trim(chtype)," occupations (integral of the Green function) are  == "
     else
       write(message,'(2a)') ch10,"  == The occupations (integral of the Green function) are  == "
     endif
   else if(optweissgreen==1) then
     write(message,'(2a)') ch10,"  == The integrals of the Weiss function are  == "
   endif
   call wrtout(std_out,message,'COLL')
   call print_oper(green%occup,option,paw_dmft,pawprtvol)
 endif
 if(mod(option,4)>=2) then
   if(optweissgreen==2) then
     write(message,'(2a)') ch10,"  == The occupations (value of G(tau) for tau=0-) are  == "
   else if(optweissgreen==1) then
     write(message,'(2a)') ch10,"  == Values of G_0(tau) for tau=0- are  == "
   endif
   call wrtout(std_out,message,'COLL')
   call print_oper(green%occup_tau,option,paw_dmft,pawprtvol)
!   write(message,'(2a)') ch10," == check: occupations from Green functions are  == "
!   call wrtout(std_out,message,'COLL')
!   call print_oper(green%occup,1,paw_dmft,pawprtvol)
   if(mod(option,4)>=3) then
     call diff_matlu("Local occup from integral of G(w) ",&
&       "Local occup from G(tau=0-) ",&
&       green%occup%matlu,green%occup_tau%matlu,paw_dmft%natom,1,tol4)
     write(message,'(2a)') ch10,&
&     '  *****  => Calculations of occupations in omega and tau spaces are coherent ****'
     call wrtout(std_out,message,'COLL')
   endif
 endif

 if(present(chtype)) then
   if(paw_dmft%prtvol>=4.and.chtype=="converged DMFT".and.green%occup%has_opermatlu==1) then
     write(message,'(4a)') ch10,"  == The DMFT occupations are  == "
     call wrtout(ab_out,message,'COLL')
     call print_matlu(green%occup%matlu,paw_dmft%natom,pawprtvol,opt_ab_out=1)
   endif
 endif

end subroutine printocc_green
!!***

!!****f* m_green/print_green
!! NAME
!! print_green
!!
!! FUNCTION
!!  print green function
!!
!! INPUTS
!!  char1 = character which describes the type of green function
!!  green  <type(green_type)>= green function data
!!  option=1 print local green function
!!         2 print KS green function
!!         3 print both local and KS green function
!!         4 print spectral function is green%w_type="real"
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawprtvol = printing option
!!  opt_wt=1 print green function as a function of frequency
!!         2 print green function as a function of imaginary time
!!
!! OUTPUT
!!
!! PARENTS
!!      impurity_solve,m_green,spectral_function
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine print_green(char1,green,option,paw_dmft,pawprtvol,opt_wt)

 use defs_basis
 use m_crystal, only : crystal_structure
 use m_oper, only : print_oper
 use m_paw_dmft, only : paw_dmft_type
 use m_io_tools, only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_green'
 use interfaces_14_hidewrite
 use interfaces_27_toolbox_oop
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type),intent(in) :: green
 integer,intent(in) :: option,pawprtvol
 integer,intent(in),optional :: opt_wt
 character(len=*), intent(in) :: char1

!local variables-------------------------------
 integer :: iall,iatom,ib,ifreq,ikpt,im,ispinor,isppol,itau
 integer :: lsub,mbandc,natom,ndim,nkpt,nspinor,nsppol,optwt
 character(len=500) :: message
 integer,allocatable :: unitgreenfunc_arr(:)
 integer,allocatable :: unitgreenloc_arr(:)
 character(len=fnlen) :: tmpfil
 character(len=1) :: tag_is,tag_is2
 character(len=4) :: tag_at
 character(len=3) :: tag_ik
 complex(dpc), allocatable :: sf(:)
! *********************************************************************
 if(present(opt_wt)) then
   optwt=opt_wt
 else
   optwt=1
 endif

 if(pawprtvol>200) then
 endif
 natom=green%oper(1)%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 if(option==1.or.option==3) then
   ABI_ALLOCATE(unitgreenfunc_arr,(natom*nsppol*nspinor))
   iall=0
   do iatom=1,natom
     if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
       ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
       call int2char4(iatom,tag_at)
       do isppol=1,nsppol
         write(tag_is,'(i1)')isppol
         do ispinor=1,nspinor
           iall=iall+1
           write(tag_is2,'(i1)')ispinor
           if(optwt==1) then
             tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-omega_iatom'//tag_at//'_isppol'//tag_is//'_ispinor'//tag_is2
           else
             tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-tau_iatom'//tag_at//'_isppol'//tag_is//'_ispinor'//tag_is2
           endif
           if(iall<=4) then
             write(message,'(3a)') ch10,"  == Print green function on file ",tmpfil
             call wrtout(std_out,message,'COLL')
           elseif(iall==5)  then
             write(message,'(3a)') ch10,"  == following values are printed in files"
             call wrtout(std_out,message,'COLL')
           endif
           unitgreenfunc_arr(iall)=300+iall-1
           open (unit=unitgreenfunc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted')
!           rewind(unitgreenfunc_arr(iall))
!           write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenfunc_arr(iall)
!           call wrtout(std_out,message,'COLL')
           write(message,'(a,a)') ch10,"# New record :"
           call wrtout(unitgreenfunc_arr(iall),message,'COLL')
           if(optwt==1) then
             do ifreq=1,green%nw
               write(message,'(2x,30(e10.3,2x))') &
&              green%omega(ifreq),&
&              (green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor),im=1,ndim)
               call wrtout(unitgreenfunc_arr(iall),message,'COLL')
             enddo
           else
             do itau=1,green%dmftqmc_l
               write(message,'(2x,30(e10.3,2x))') &
&              green%tau(itau),&
&              (green%oper_tau(itau)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor),im=1,ndim)
               call wrtout(unitgreenfunc_arr(iall),message,'COLL')
             enddo
           endif
         enddo
       enddo ! isppol
     endif ! lpawu=/-1
   enddo ! iatom
   ABI_DEALLOCATE(unitgreenfunc_arr)
 endif
 if(option==2.or.option==3) then
   ABI_ALLOCATE(unitgreenloc_arr,(nsppol*nkpt))
   iall=0
   do isppol = 1 , nsppol
     write(tag_is,'(i1)')isppol
     do ikpt = 1, nkpt
       write(tag_ik,'(i3)')ikpt
!      do ib1 = 1, mbandc
       iall=iall+1
       if(optwt==1) then
         tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-omega_isppol'//tag_is//'_ikpt'//trim(adjustl(tag_ik))
       else
         tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-tau_isppol'//tag_is//'_ikpt'//trim(adjustl(tag_ik))
       endif
       if(iall<=4)  then
         write(message,'(3a)') ch10,"  == Print green function on file ",tmpfil
         call wrtout(std_out,message,'COLL')
       elseif(iall==5)  then
         write(message,'(3a)') ch10,"  == following values are printed in files"
         call wrtout(std_out,message,'COLL')
       endif
       unitgreenloc_arr(iall)=400+iall-1
       open (unit=unitgreenloc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted')
!      rewind(unitgreenloc_arr(iall))
!       write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenloc_arr(iall)
!       call wrtout(std_out,message,'COLL')
       write(message,'(a,a)') ch10,"# New record : First 20 bands"
       call wrtout(unitgreenloc_arr(iall),message,'COLL')
!       call flush(std_out)
       do lsub=1,mbandc/20+1
         if(optwt==1) then
           do ifreq=1,green%nw
!             call flush(std_out)
             write(message,'(2x,50(e10.3,2x))') &
&            green%omega(ifreq), &
&            (green%oper(ifreq)%ks(isppol,ikpt,ib,ib),ib=20*(lsub-1)+1,min(20*lsub,mbandc))
             call wrtout(unitgreenloc_arr(iall),message,'COLL')
           enddo
         else
           do itau=1,green%dmftqmc_l
!             call flush(std_out)
             write(message,'(2x,50(e10.3,2x))') &
&            green%tau(itau), &
&            (green%oper_tau(itau)%ks(isppol,ikpt,ib,ib),ib=20*(lsub-1)+1,min(20*lsub,mbandc))
             call wrtout(unitgreenloc_arr(iall),message,'COLL')
           enddo
         endif
         if(20*lsub<mbandc) write(message,'(a,a,i5,a,i5)')    &
&           ch10,"# Same record, Following bands : From ",    &
&           20*(lsub),"  to ",min(20*(lsub+1),mbandc)
         call wrtout(unitgreenloc_arr(iall),message,'COLL')
       enddo
     enddo ! ikpt
   enddo ! isppol
   ABI_DEALLOCATE(unitgreenloc_arr)
 endif

 if(green%w_type=="real".and.option==4) then
   write(message,'(a,a)') ch10,"  == About to print spectral function"
   call wrtout(std_out,message,'COLL')
   tmpfil = trim(paw_dmft%filapp)//'SpFunc-'//trim(char1)
   open (unit=3333,file=trim(tmpfil),status='unknown',form='formatted')
   ABI_ALLOCATE(sf,(green%nw))
   iall=0
   sf=czero
   do isppol = 1 , nsppol
     do ikpt = 1, nkpt
       do ib=1,mbandc
         iall=iall+1
         do ifreq=1,green%nw
           sf(ifreq)=sf(ifreq)+green%oper(ifreq)%ks(isppol,ikpt,ib,ib)*green%oper(1)%wtk(ikpt)
         enddo
       enddo
     enddo
   enddo
   do ifreq=1,green%nw
     write(message,*) green%omega(ifreq),(-imag(sf(ifreq)))/3.141592653589793238_dp
     call wrtout(3333,message,'COLL')
   enddo
   ABI_DEALLOCATE(sf)
 endif
 call flush_unit(3333)

end subroutine print_green
!!***

!!****f* m_green/compute_green
!! NAME
!! compute_green
!!
!! FUNCTION
!! compute green function from LDA and self-energy
!!
!! INPUTS
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  opt_self = optional argument, if =1, upfold self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!
!! PARENTS
!!      dmft_solve,fermi_green,m_green,newton,spectral_function
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine compute_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,prtopt,self,opt_self)

 use defs_basis
 use defs_abitypes
 use m_crystal, only : crystal_structure
 use m_matlu, only : sym_matlu,zero_matlu,add_matlu,print_matlu
 use m_oper, only : inverse_oper,loc_oper,print_oper,upfold_oper,init_oper,destroy_oper
 use m_paw_dmft, only : paw_dmft_type
 use m_self, only : self_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_green'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(crystal_structure),intent(in) :: cryst_struc
 type(green_type),intent(out) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: opt_self

!local variables-------------------------------
 integer :: iatom,ib,ib1,ierr,ifreq,ikpt,is
 integer :: myproc,nproc,nspinor,nsppol,option,optself,spaceComm
 character(len=500) :: message
 real(dp) :: fermilevel
 real(dp), allocatable :: Id(:,:)
 type(oper_type) :: self_minus_hdc_oper
 complex(dpc) :: omega_current
! *********************************************************************
 if(present(opt_self)) then
   optself=opt_self
 else
   optself=0
 endif

 if(prtopt>0)  then
   write(message,'(2a,i3,13x,a)') ch10,'  ===  Compute green function '
   call wrtout(std_out,message,'COLL')
 endif

 if(self%nw/=green%nw)  then
   write(message,'(2a,i3,13x,a)') ch10,' BUG: frequencies for green and self not coherent'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 endif

! Initialise spaceComm, myproc, and nproc
 call xcomm_init(mpi_enreg,spacecomm)
 call xme_init(mpi_enreg,myproc)
 call xproc_init(mpi_enreg,nproc)

! Initialise integers
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol

! Initialise for compiler
 omega_current=czero

! ==============================
! Initialise Identity
! ==============================
 ABI_ALLOCATE(Id,(paw_dmft%mbandc,paw_dmft%mbandc))
 Id=zero
 do ib=1, paw_dmft%mbandc
   Id(ib,ib)=one
 enddo ! ib

 if(prtopt/=0) then
   write(message,'(2a)') ch10,'  == Green function is computed:'
   call wrtout(std_out,message,'COLL')
 endif
 option=1
 fermilevel=paw_dmft%fermie
 if(option==123) then
   fermilevel=2.d0
   write(message,'(2a,e14.3,a)') ch10,&
&  '  Warning (special case for check: fermi level=',fermilevel,')'
   call wrtout(std_out,message,'COLL')
 endif

! ====================================================
! Upfold self-energy and double counting  Self_imp -> self(k)
! ====================================================
! if(optself==1) then
!   do ifreq=1,green%nw
!     call upfold_oper(self%oper(ifreq),paw_dmft,1)
!   enddo ! ifreq
!   call upfold_oper(self%hdc,paw_dmft,1)
! endif

! =================================================================
! Initialize green%oper before calculation (important for xsum_mpi)
! =================================================================
 do ifreq=1,green%nw
  do is=1,nsppol
     do ikpt=1,paw_dmft%nkpt
       do ib=1,paw_dmft%mbandc
         do ib1=1,paw_dmft%mbandc
           green%oper(ifreq)%ks(is,ikpt,ib,ib1)=czero
         enddo
       enddo
     enddo
   enddo
   call zero_matlu(green%oper(ifreq)%matlu,green%oper(ifreq)%natom)
 enddo ! ifreq
 call xbarrier_mpi(spacecomm)

! ================================
! == Compute Green function G(k)
! ================================
 do ifreq=1,green%nw
!   ====================================================
!   First Upfold self-energy and double counting  Self_imp -> self(k)
!   ====================================================
   if(mod(ifreq-1,nproc)==myproc) then
     if(green%w_type=="imag") then
       omega_current=cmplx(zero,green%omega(ifreq),kind=dp)
     else if(green%w_type=="real") then
       omega_current=cmplx(green%omega(ifreq),0.1/27.211_dp,kind=dp)
     endif
     call init_oper(paw_dmft,self_minus_hdc_oper)
     call add_matlu(self%oper(ifreq)%matlu,self%hdc%matlu,self_minus_hdc_oper%matlu,green%oper(ifreq)%natom,-1)
!   write(std_out,*) self_minus_hdc_oper%matlu(1)%mat
     call upfold_oper(self_minus_hdc_oper,paw_dmft,1)
     do is = 1 , paw_dmft%nsppol
       do ib = 1 , paw_dmft%mbandc
         do ib1 = 1 , paw_dmft%mbandc
           do ikpt = 1 , paw_dmft%nkpt
             green%oper(ifreq)%ks(is,ikpt,ib,ib1)=       &
&             ( omega_current     &
&             + fermilevel                               &
&             - paw_dmft%eigen_lda(is,ikpt,ib)) * Id(ib,ib1) &
&             - self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
!&             - (self%oper(ifreq)%ks(is,ikpt,ib,ib1)-self%hdc%ks(is,ikpt,ib,ib1))
!             if(prtopt>5) then
!             if(ikpt==1.and.(ifreq==1.or.ifreq==3).and.ib==23.and.ib1==23) then
!              write(std_out,*) 'omega_current                         ',omega_current
!              write(std_out,*) 'fermilevel                            ',fermilevel
!              write(std_out,*) ' paw_dmft%eigen_lda(is,ikpt,ib)       ', paw_dmft%eigen_lda(is,ikpt,ib),Id(ib,ib1)
!              write(std_out,*) 'self_minus_hdc_oper%ks(is,ikpt,ib,ib1)',self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
!              write(std_out,*) 'green                                 ',green%oper(ifreq)%ks(is,ikpt,ib,ib1)
!             endif
!             endif
           enddo ! ikpt
         enddo ! ib1
       enddo ! ib
     enddo ! is
!     call print_oper(green%oper(ifreq),9,paw_dmft,3)
!     if(ifreq==1.or.ifreq==3) then
!       write(std_out,*) 'green1  ifreq         %ks(1,1,23,23)',ifreq,green%oper(ifreq)%ks(1,1,23,23)
!     endif
     call inverse_oper(green%oper(ifreq),2,prtopt)
!     if(ifreq==1.or.ifreq==3) then
!       write(std_out,*) 'green1afterinversion  %ks(1,1,23,23)',ifreq,green%oper(ifreq)%ks(1,1,23,23)
!     endif
! ================================
! == Compute Local Green function
! ================================
     call loc_oper(green%oper(ifreq),paw_dmft,1)
     call sym_matlu(cryst_struc,green%oper(ifreq)%matlu,pawang)
     call destroy_oper(self_minus_hdc_oper)
   endif ! parallelisation
 enddo ! ifreq

! =============================================
! built total green function (sum over procs).
! =============================================
 call xbarrier_mpi(spacecomm)
 do ifreq=1,green%nw
   call xsum_mpi(green%oper(ifreq)%ks,spacecomm,ierr)
   call xbarrier_mpi(spacecomm)
   do iatom=1,green%oper(ifreq)%natom
     if(green%oper(ifreq)%matlu(iatom)%lpawu.ne.-1) then
       call xsum_mpi(green%oper(ifreq)%matlu(iatom)%mat,spacecomm,ierr)
       call xbarrier_mpi(spacecomm)
     endif
   enddo ! iatom
 enddo ! ifreq
 call xbarrier_mpi(spacecomm)
! write(std_out,*) 'afterxsum sym     %matlu(1)%mat(2,5,1,1,1) 1',green%oper(1)%matlu(1)%mat(2,5,1,1,1)

 if(prtopt/=0) then
   write(message,'(2a)') ch10,&
&   '  == Local Green function has been computed and projected on local orbitals'
   call wrtout(std_out,message,'COLL')
 endif
! useless test
 if(abs(prtopt)>=4) then
   write(message,'(2a)') ch10,' == Green function is now printed for first frequency'
   call wrtout(std_out,message,'COLL')
   call print_oper(green%oper(1),2,paw_dmft,2)
 endif


 ABI_DEALLOCATE(Id)

end subroutine compute_green
!!***

!!****f* m_green/integrate_green
!! NAME
!! integrate_green
!!
!! FUNCTION
!!  integrate green function in the band index basis
!!
!! INPUTS
!!  cryst_struc <type(crystal_structure)>=crystal structure data
!!  green <type(green_type)>=green function  (green%oper(:))
!!  opt_ksloc=   1: do only integrate on the KS basis
!!               2: do the integration on the local basis
!!                 (This can work only if psichi are renormalized!!!)
!!               3: do both calculations and test the consistency of it
!!              -1: do the integration on the KS basis, but only
!!                      compute diagonal part of the band-band density matrix
!!                      in order to compute the total charge for fermi_green
!!  paw_dmft <type(m_paw_dmft)>= paw+dmft data
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!   green%occup = occupations
!!
!! PARENTS
!!      dmft_solve,fermi_green,impurity_solve,m_green,newton
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine integrate_green(cryst_struc,green,mpi_enreg,paw_dmft&
&  ,pawang,prtopt,opt_ksloc,opt_after_solver,opt_diff)

 use defs_basis
 use defs_abitypes
 use m_crystal, only : crystal_structure
 use m_matlu, only : sym_matlu,print_matlu,init_matlu,&
& destroy_matlu,nullify_matlu,diff_matlu,zero_matlu
 use m_paw_dmft, only : paw_dmft_type
 use m_oper, only : loc_oper,trace_oper,init_oper,destroy_oper
! use m_oper, only : loc_oper,trace_oper,upfold_oper,print_oper,identity_oper,init_oper,destroy_oper
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'integrate_green'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(crystal_structure),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: opt_ksloc
 integer, optional, intent(in) :: opt_after_solver
 integer, optional, intent(in) :: opt_diff

!local variables-------------------------------
 integer :: iatom,ib,ib1,icomp_chloc,ierr,ifreq,ikpt,im,im1,ispinor,ispinor1,is
 integer :: mband,mbandc,myproc,natom,ndim,nkpt,nproc,nspinor
 integer :: nsppol,option
 integer :: optksloc,spacecomm,optaftsolv
 complex(dpc) :: integral
 character(len=500) :: message
 complex(dpc), allocatable :: ff(:)
 type(matlu_type), allocatable :: matlu_temp(:)
 type(oper_type) :: occup_temp
 real(dp) :: diff_chloc
! real(dp), allocatable :: charge_loc_old(:,:)
! type(oper_type)  :: oper_c
! *********************************************************************

 DBG_ENTER("COLL")

 if(prtopt>0) then
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Integrate green function'
   call wrtout(std_out,message,'COLL')
 endif
 if(green%w_type=="real") then
   write(message,'(2a,i3,13x,a)') ch10,'   BUG: integrate_green not implemented for real frequency'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 endif

! Initialise spaceComm, myproc, and master
 call xcomm_init(mpi_enreg,spacecomm)
 call xme_init(mpi_enreg,myproc)
 call xproc_init(mpi_enreg,nproc)

! Initialise integers
 mband   = paw_dmft%mband
 mbandc  = paw_dmft%mbandc
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 natom   = paw_dmft%natom

! Initialize green%oper before calculation (important for xsum_mpi)
! allocate(charge_loc_old(paw_dmft%natom,paw_dmft%nsppol+1))
! if(.not.present(opt_diff)) then  ! if integrate_green is called in dmft_solve after calculation of self
!   charge_loc_old=green%charge_matlu
! endif
 icomp_chloc=0

! Choose what to compute
 if(present(opt_ksloc)) then
   optksloc=opt_ksloc
 else
   optksloc=3
 endif
 if(present(opt_after_solver)) then
   optaftsolv=opt_after_solver
 else
   optaftsolv=0
 endif
 if(optaftsolv==1.and.abs(optksloc)/=2) then
    write(message,'(2a)') ch10,&
&     "BUG: integration of ks green function should not be done after call to solver : it has not been computed"
    call wrtout(std_out,message,'COLL')
    call leave_new('COLL')
 endif

! Allocations
 ABI_ALLOCATE(matlu_temp,(natom))
 call init_matlu(natom,nspinor,nsppol,green%oper(1)%matlu(:)%lpawu,matlu_temp)

 ABI_ALLOCATE(ff,(green%nw))

! =================================================
! == Integrate Local Green function ===============
 if(abs(optksloc)/2==1) then ! optksloc=2 or 3
! =================================================
   call zero_matlu(green%occup%matlu,green%occup%natom)

! ==  Calculation of \int{G_{LL'}{\sigma\sigma',s}(R)(i\omega_n)}
   if(paw_dmft%lpsichiortho==1) then
!  - Calculation of frequency sum over positive frequency
     if (nspinor==1) option=1
     if (nspinor==2) option=2
     do iatom=1, natom
       ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         do is = 1 , nsppol
           do ispinor = 1, nspinor
             do ispinor1 = 1, nspinor
               do im=1,ndim
                 do im1=1,ndim
                   do ifreq=1,green%nw
                     ff(ifreq)= &
&                     green%oper(ifreq)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
!                     write(std_out,*) green%omega(ifreq),ff(ifreq)," integrate green fw_lo"
                   enddo
!                   call int_fct(ff,(im==im1).and.(ispinor==ispinor1),&
!&                   option,paw_dmft,integral)
                   call int_fct(ff,(im==im1).and.(ispinor==ispinor1),&
&                   2,paw_dmft,integral)  ! test_1
                   green%occup%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=integral
!                   if(im==2.and.im1==5.and.is==1.and.iatom==1) then
!                     write(std_out,*) " occup        %matlu(1)%mat(2,5,1,1,1)",green%occup%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
!                   endif
                 enddo
               enddo
             enddo ! ispinor1
           enddo ! ispinor
         enddo ! is
         matlu_temp(iatom)%mat=green%occup%matlu(iatom)%mat
       endif ! lpawu=/-1
     enddo ! iatom

!   Print density matrix if prtopt high
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,"  = green%occup%matlu from int(gloc(w))"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif

!  - Symetrise: continue sum over k-point: Full BZ
     call sym_matlu(cryst_struc,green%occup%matlu,pawang)
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,&
&       "  = green%occup%matlu from int(gloc(w)) with symetrization"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif
     call sym_matlu(cryst_struc,matlu_temp,pawang)

!  - Post-treatment for summation over negative and positive frequencies:
!    necessary in the case of nspinor==2 AND nspinor==1, but valid anywhere
!    N(ll'sigmasigma')= (N(ll'sigmasigma')+ N*(l'lsigma'sigma))/2
!    because [G_{LL'}^{sigma,sigma'}(iomega_n)]*= G_{L'L}^{sigma',sigma}(-iomega_n)
     if(nspinor>=1) Then
       do iatom=1, natom
         ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
         if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
           do is = 1 , nsppol
             do ispinor = 1, nspinor
               do ispinor1 = 1, nspinor
                 do im=1,ndim
                   do im1=1,ndim
                     green%occup%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)= &
&                          (matlu_temp(iatom)%mat(im,im1,is,ispinor,ispinor1)+ &
&                           conjg(matlu_temp(iatom)%mat(im1,im,is,ispinor1,ispinor)))/two
                   enddo
                 enddo
               enddo ! ispinor1
             enddo ! ispinor
           enddo ! isppol
           matlu_temp(iatom)%mat=green%occup%matlu(iatom)%mat
         endif
       enddo ! iatom
     endif
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,&
&       "  = green%occup%matlu from int(gloc(w)) symetrized with post-treatment"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif
     if(optaftsolv==0) then
       call trace_oper(green%occup,green%charge_ks,green%charge_matlu,2)
       green%has_charge_matlu=2
       green%has_charge_matlu_solver=0
       icomp_chloc=1
     else if(optaftsolv==1) then
!      This is done only when called from impurity_solver with solver
!      green function.
       call trace_oper(green%occup,green%charge_ks,green%charge_matlu_solver,2)
       green%has_charge_matlu_solver=2
     endif
   else
     write(message,'(a,4x,a,a,a,4x,a)') ch10,&
&     " Local basis is not (yet) orthonormal:",&
&     " local green function is thus not integrated",ch10,&
&     " Local occupations are computed from KS occupations"
     call wrtout(std_out,message,'COLL')
   endif

 endif ! optksloc
! =================================================



! =================================================
! == Integrate Kohn Sham Green function ===========
 if(mod(abs(optksloc),2)==1) then ! optksloc=1 or 3 or -1
   green%occup%ks=czero ! important for xsum_mpi
! =================================================
! ==  Calculation of \int{G_{\nu\nu'}{k,s}(i\omega_n)}
   call init_oper(paw_dmft,occup_temp)
   do is = 1 , nsppol
     do ikpt = 1, nkpt
       if(mod(ikpt-1,nproc)==myproc) then
         do ib = 1, mbandc
           do ib1 = 1, mbandc
             if(optksloc==-1.and.(ib/=ib1)) cycle
             do ifreq=1,green%nw
               ff(ifreq)=green%oper(ifreq)%ks(is,ikpt,ib,ib1)
             enddo
!             call int_fct(ff,ib==ib1,nspinor,paw_dmft,integral) ! here, option==1 even if nspinor==2
!             green%occup%ks(is,ikpt,ib,ib1)=integral
             call int_fct(ff,ib==ib1,2,paw_dmft,integral) ! here, option==1 even if nspinor==2
             green%occup%ks(is,ikpt,ib,ib1)=integral
!        write(std_out,'(a,4i6,e14.5,e14.5)') "ks",is,ikpt,ib,ib1,integral
           enddo ! ib1
         enddo ! ib
       endif
     enddo ! ikpt
   enddo ! isppol
   call xbarrier_mpi(spacecomm)
   call xsum_mpi(green%occup%ks,spacecomm,ierr)
   occup_temp%ks=green%occup%ks
!  - Post-treatment for summation over negative and positive frequencies:
!    necessary in the case of nspinor==2, but valid everywhere
!    N(k,n_1,n_2)= (N(k,n_1,n_2)+ N*(k,n_2,n_1))/2
!    because [G_{k}^{n_1,n_2}(iomega_n)]*= G_{k}^{n_2,n_1}(-iomega_n)
   do is = 1 , nsppol
     do ikpt = 1, nkpt
       do ib = 1, mbandc
         do ib1 = 1, mbandc
           green%occup%ks(is,ikpt,ib,ib1)=&
&            (       occup_temp%ks(is,ikpt,ib,ib1)+ &
&              conjg(occup_temp%ks(is,ikpt,ib1,ib))   )/two
         enddo ! ib1
       enddo ! ib
     enddo ! ikpt
   enddo ! isppol
   call destroy_oper(occup_temp)
   do is = 1 , nsppol
     do ikpt = 1, nkpt
       do ib = 1, mbandc
         do ib1 = 1, mbandc
           paw_dmft%occnd(paw_dmft%include_bands(ib),&
&           paw_dmft%include_bands(ib1),ikpt,is)=green%occup%ks(is,ikpt,ib,ib1)
         enddo
       enddo
     enddo
   enddo

   if(optksloc>0) then
!  - Compute local occupations
     call loc_oper(green%occup,paw_dmft,1)
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,&
&        "  = green%occup%matlu from projection of int(gks(w)) without symetrization"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif

!  - Symetrise: continue sum over k-point: Full BZ
     call sym_matlu(cryst_struc,green%occup%matlu,pawang)
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,&
&        "  = green%occup%matlu from projection of int(gks(w)) with symetrization"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif

!  - If the trace of occup matrix in the LOCAL basis was not done
!  before because lpsichiortho/=1 , do it now
     if(paw_dmft%lpsichiortho/=1) then
       if(abs(prtopt)>0) then
         call trace_oper(green%occup,green%charge_ks,green%charge_matlu,2)
         green%has_charge_matlu=2
         icomp_chloc=1
       endif
     endif
   endif ! optksloc>0: only diagonal elements of G\nunu' are computed

!  - Compute trace over ks density matrix
   call trace_oper(green%occup,green%charge_ks,green%charge_matlu,1)
   if(abs(prtopt)>0) then
     write(message,'(a,a,f12.6)') ch10,&
&    "  ==  Total number of electrons from KS green function is :", green%charge_ks
     call wrtout(std_out,message,'COLL')
     write(message,'(8x,a,f12.6,a)') " (should be",paw_dmft%nelectval,")"
     call wrtout(std_out,message,'COLL')
   endif
 endif ! optksloc
! =================================================


! =================================================
! Tests and compute precision on local charge
! =================================================
!  - Check that if, renormalized psichi are used, occupations matrices
!    obtained directly from local green function or, through kohn sham
!    occupations are the same.
 if(abs(optksloc)==3) then ! optksloc= 3
   if(paw_dmft%lpsichiortho==1) then
     call diff_matlu("Local_projection_of_kohnsham_occupations ",&
&     "Integration_of_local_green_function ",&
&       green%occup%matlu,matlu_temp,natom,1,tol4)
     write(message,'(2a)') ch10,&
&     '  ***** => Calculations of Green function in KS and local spaces are coherent ****'
     call wrtout(std_out,message,'COLL')
   endif
 endif

!!***

! == Precision on charge_matlu (done only if local charge was computed ie not for optksloc=-1)
 if(icomp_chloc==1.and.paw_dmft%idmftloop>=1.and.present(opt_diff)) then ! if the computation was done here.
   if(green%has_charge_matlu_prev==2) then
     diff_chloc=zero
     do iatom=1, natom
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         do is=1,nsppol
           diff_chloc=diff_chloc+&
&           (green%charge_matlu_prev(iatom,is)-green%charge_matlu(iatom,is))**2
         enddo
       endif
     enddo
     if(sqrt(diff_chloc)<paw_dmft%dmft_lcpr) then
       green%ichargeloc_cv=1
       write(message,'(a,8x,a,e9.2,a,8x,a)') ch10,"Change of local charge =<",paw_dmft%dmft_lcpr,&
&       ch10,"DMFT Loop: Local charge is converged"
       write(std_out,*) sum(green%charge_matlu_prev),sum(green%charge_matlu)
       call wrtout(std_out,message,'COLL')
     else
       green%ichargeloc_cv=0
       write(message,'(a,8x,a)') ch10,"DMFT Loop: Local charge is not converged"
       call wrtout(std_out,message,'COLL')
     endif
   endif
   green%charge_matlu_prev=green%charge_matlu
   green%has_charge_matlu_prev=2
 endif


 ABI_DEALLOCATE(ff)
 call destroy_matlu(matlu_temp,natom)
 call nullify_matlu(matlu_temp,natom)
 ABI_DEALLOCATE(matlu_temp)
! deallocate(charge_loc_old)
 DBG_EXIT("COLL")

end subroutine integrate_green

!!      m_green
!!***

!!****f* m_green/icip_green
!! NAME
!!  lda_green
!!
!! FUNCTION
!!  init, compute, integrate and print lda green function
!!
!! INPUTS
!!  char1 = character which precises the type of green function computed.
!!  cryst_struc <type(crystal_structure)>= crystal structure data.
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawprtvol  = option for printing
!!  self <type(self_type)>= variables related to self-energy
!!  opt_self = optional argument, if =1, upfold self-energy
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!   wrtout
!!
!! SOURCE

subroutine icip_green(char1,cryst_struc,green,mpi_enreg,paw_dmft,pawang,pawprtvol,self,opt_self)

 use defs_basis
 use defs_abitypes
 use m_crystal, only : crystal_structure
 use m_matlu, only : sym_matlu
 use m_paw_dmft, only : paw_dmft_type
 use m_oper, only : loc_oper
 use m_self, only : self_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'icip_green'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(crystal_structure),intent(in) :: cryst_struc
 type(MPI_type), intent(in) :: mpi_enreg
 type(green_type),intent(inout) :: green
 type(pawang_type),intent(in) :: pawang
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(self_type), intent(inout) :: self
 integer, intent(in) :: pawprtvol
 character(len=*), intent(in) :: char1
 integer, optional, intent(in) :: opt_self
!Local variables ------------------------------------
 integer :: option,optself,optlocks
 character(len=500) :: message
! *********************************************************************
 if(present(opt_self)) then
   optself=opt_self
 else
   optself=0
 endif

 call init_green(green,paw_dmft)

! useless test ok
! call printocc_green(green,1,paw_dmft,2)
! write(std_out,*)" printocc_green zero finished "

!== Compute green%oper(:)%ks
!== Deduce  green%oper(:)%matlu(:)%mat
 call compute_green(cryst_struc,green,mpi_enreg,paw_dmft,&
& pawang,pawprtvol,self,optself)
 if(paw_dmft%dmft_prgn>=1.and.paw_dmft%dmft_prgn<=2) then
   optlocks=paw_dmft%dmft_prgn*2+1 ! if dmft_prgn==2 => do not print
   if(paw_dmft%lpsichiortho==1)  then
     call print_green(char1,green,optlocks,paw_dmft,pawprtvol)
   endif
 endif

!== Integrate green%oper(:)%ks
!== Integrate green%oper(:)%matlu(:)%mat
 call integrate_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,3,opt_ksloc=3)

!== Print green%oper(:)%ks
!== Print green%oper(:)%matlu(:)%mat
 if(char1=="LDA") then
   option=1
   if(self%oper(1)%matlu(1)%lpawu/=-1) then
     if(abs(real(self%oper(1)%matlu(1)%mat(1,1,1,1,1)))>tol7) then
! todo_ab: generalise this
       write(message,'(a,a,2(e15.4))') ch10,&
&        "ERROR:  In LDA calculation, self should be zero"
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     endif
   endif
 else
   option=5
 endif
 call printocc_green(green,option,paw_dmft,3,chtype=char1)

end subroutine icip_green
!!***

!!****f* m_green/fourier_green
!! NAME
!! fourier_green
!!
!! FUNCTION
!!  integrate green function in the band index basis
!!
!! INPUTS
!!  cryst_struc <type(crystal_structure)>= crystal structure data.
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!
!! PARENTS
!!      impurity_solve,m_green
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine fourier_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang,opt_ksloc,opt_tw)

 use defs_basis
 use defs_abitypes
 use m_crystal, only : crystal_structure
 use m_matlu, only : sym_matlu,init_matlu,destroy_matlu,nullify_matlu,print_matlu
 use m_paw_dmft, only : paw_dmft_type
 use m_oper, only : loc_oper

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourier_green'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(crystal_structure),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer,intent(in) :: opt_ksloc,opt_tw ! fourier on ks or local

!local variables-------------------------------
 integer :: iatom,ib,ib1,ierr,ifreq,ikpt,im,im1,iparal,is,ispinor,ispinor1,itau
 integer :: mband,mbandc,myproc,natom,ndim,nkpt,nproc,nspinor,nsppol,spacecomm!,opt_four
 character(len=500) :: message
! complex(dpc):: ybcbeg,ybcend
! arrays
 complex(dpc), allocatable :: fw(:)
 complex(dpc), allocatable :: ft(:)
 type(green_type) :: green_temp
! *********************************************************************
! ybcbeg=czero
! ybcend=czero

! define spaceComm, myproc, and nproc from world communicator
! and mpi_enreg
 call xcomm_init(mpi_enreg,spacecomm)
 call xme_init(mpi_enreg,myproc)
 call xproc_init(mpi_enreg,nproc)

! Initialise integers
 mband   = paw_dmft%mband
 mbandc  = paw_dmft%mbandc
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 natom   = cryst_struc%natom

! Only imaginary frequencies here
 if(green%w_type=="real") then
   write(message,'(2a,i3,13x,a)') ch10,'   BUG: fourier_green not implemented for real frequency'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 endif

! Initialise temporary green function
 call init_green(green_temp,paw_dmft)
 call init_green_tau(green_temp,paw_dmft)

 !green%oper(:)%matlu(1)%mat(1,1,1,1,1)
 ABI_ALLOCATE(fw,(green%nw))
 ABI_ALLOCATE(ft,(green%dmftqmc_l))

!==============================================
! == Inverse fourier transformation ===========
!==============================================

 if(opt_tw==-1) then

!  = For Kohn Sham green function
!==============================================
   if(opt_ksloc ==1) then

     do is = 1 , nsppol
       do ikpt = 1, nkpt
         do ib = 1, mbandc
           do ib1 = 1, mbandc
             do ifreq=1,green%nw
               fw(ifreq)=green%oper(ifreq)%ks(is,ikpt,ib,ib1)
             enddo
             call fourier_fct(fw,ft,ib==ib1,green%dmftqmc_l,-1,paw_dmft) ! inverse fourier
             do itau=1,green%dmftqmc_l
               green_temp%oper_tau(itau)%ks(is,ikpt,ib,ib1)=ft(itau)
             enddo
             if(ib==ib1) then
               green%occup_tau%ks(is,ikpt,ib,ib1)=ft(1)+one
             else
               green%occup_tau%ks(is,ikpt,ib,ib1)=ft(1)
             endif
           enddo ! ib1
         enddo ! ib
       enddo ! ikpt
     enddo ! isppol

!  = Post-treatment: necessary in the case of nspinor==2, but valid anywhere
!    because G(tau)=G_{nu,nu'}(tau)+[G_{nu,nu'}(tau)]*

     do is = 1 , nsppol
       do ikpt = 1, nkpt
         do ib = 1, mbandc
           do ib1 = 1, mbandc
             do itau=1,green%dmftqmc_l
               green%oper_tau(itau)%ks(is,ikpt,ib,ib1)=&
&                (green_temp%oper_tau(itau)%ks(is,ikpt,ib,ib1)+ &
&                 conjg(green_temp%oper_tau(itau)%ks(is,ikpt,ib1,ib)))/two
               if(ib==ib1) then
                 green%occup_tau%ks(is,ikpt,ib,ib1)=green%oper_tau(1)%ks(is,ikpt,ib,ib1)+one
               else
                 green%occup_tau%ks(is,ikpt,ib,ib1)=green%oper_tau(1)%ks(is,ikpt,ib,ib1)
               endif
             enddo
           enddo ! ib1
         enddo ! ib
       enddo ! ikpt
     enddo ! isppol
     call loc_oper(green%occup_tau,paw_dmft,1)
     call sym_matlu(cryst_struc,green%occup_tau%matlu,pawang)
     write(message,'(a,a,i10,a)') ch10,"  green%occup_tau%matlu from green_occup_tau%ks"
     call wrtout(std_out,message,'COLL')
     call print_matlu(green%occup_tau%matlu,natom,prtopt=3)
   endif

!  = For local green function
!==============================================
   if(opt_ksloc ==2) then

     iparal=0
     do iatom=1, natom
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
         do is = 1 , nsppol
           do ispinor = 1, nspinor
             do ispinor1 = 1, nspinor
               do im=1,ndim
                 do im1=1,ndim
                   iparal=iparal+1
                   if(mod(iparal-1,nproc)==myproc) then
                     do ifreq=1,green%nw
                       fw(ifreq)=green%oper(ifreq)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     enddo
                     ! inverse fourier
!                     write(std_out,'(a)') "fourierbeforeposttreatement,ispinor,ispinor1,is,im,im1"
!                     write(std_out,'(a,5i4,f12.5,f12.5)') "fourier",ispinor,ispinor1,is,im,im1
!                     write(std_out,'(a,e12.5,e12.5)')&
!                     &"green%oper(4)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)"&
!&                     ,green%oper(4)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     call fourier_fct(fw,ft,(im==im1).and.(ispinor==ispinor1),green%dmftqmc_l,-1,paw_dmft)
                     do itau=1,green%dmftqmc_l
                       green_temp%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=ft(itau)
!                     write(std_out,*) itau,green_temp%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     enddo
!                    if((im==im1).and.(ispinor==ispinor1)) then
!                      green%occup_tau%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=ft(1)+one
!                    else
!                      green%occup_tau%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=ft(1)
!                    endif
                   endif ! iparal
                 enddo
               enddo
             enddo ! ispinor1
           enddo ! ispinor
         enddo ! isppol
       endif ! lpawu.ne.-1
     enddo ! iatom

!    Parallelisation must be finished here, because in the post
!    treatment, value from different proc will be mixed.
     call xbarrier_mpi(spacecomm)
     do iatom=1,natom
       do itau=1,green%dmftqmc_l
        call xsum_mpi(green_temp%oper_tau(itau)%matlu(iatom)%mat,spacecomm,ierr)
       enddo
     enddo
     call xbarrier_mpi(spacecomm)

!  = Post-treatment: necessary in the case of nspinor==2, but valid anywhere
!    because G(tau)=G_{LL'}^{sigma,sigma'}(tau)+[G_{L'L}^{sigma',sigma}(tau)]*

     if(nspinor>0) Then
       do iatom=1, natom
         if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
           ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
           do is = 1 , nsppol
             do ispinor = 1, nspinor
               do ispinor1 = 1, nspinor
                 do im=1,ndim
                   do im1=1,ndim
!                   write(std_out,'(a,5i4,f12.5,f12.5)') "fourier -1",ispinor,ispinor1,is,im,im1
                     do itau=1,green%dmftqmc_l
                       green%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=&
&                      ((green_temp%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)+ &
&                      conjg(green_temp%oper_tau(itau)%matlu(iatom)%mat(im1,im,is,ispinor1,ispinor))))/two
!                       write(std_out,*) itau,green%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                       if((im==im1).and.(ispinor==ispinor1)) then
                         green%occup_tau%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=&
&                         green%oper_tau(1)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)+one
                       else
                         green%occup_tau%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=&
&                         green%oper_tau(1)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                       endif ! test diag
                     enddo  ! itau
                   enddo  ! im1
                 enddo ! im
               enddo ! ispinor1
             enddo ! ispinor
           enddo ! isppol
         endif
       enddo ! iatom
     endif ! nspinor>0
   endif ! opt_ksloc=2
 endif ! opt_tw==-1


!==============================================
! == Direct fourier transformation
!==============================================
! todo_ba ft useful only for diagonal elements ...

 if(opt_tw==1) then

!  = For local green function
!==============================================
   if(opt_ksloc ==2) then

     iparal=0
     do iatom=1, natom
!  put to zero (usefull at least for parallelism)
       do ifreq=1,green%nw
         green%oper(ifreq)%matlu(iatom)%mat=czero
       enddo
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
         do is = 1 , nsppol
           do ispinor = 1, nspinor
             do ispinor1 = 1, nspinor
               do im=1,ndim
                 do im1=1,ndim
                   iparal=iparal+1
                   if(mod(iparal-1,nproc)==myproc) then
                     do itau=1,green%dmftqmc_l
                       ft(itau)=green%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     enddo
!                     write(std_out,'(a,5i4,f12.5,f12.5)') "fourier",ispinor,ispinor1,is,im,im1
                     call fourier_fct(fw,ft,(im==im1).and.(ispinor==ispinor1),green%dmftqmc_l,1,paw_dmft)
                     do ifreq=1,green%nw
                       green%oper(ifreq)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=fw(ifreq)
                     enddo
                   endif
                 enddo
               enddo
             enddo ! ispinor1
           enddo ! ispinor
         enddo ! isppol
       endif ! lpawu=/-1
     enddo ! iatom
     call xbarrier_mpi(spacecomm)
     do iatom=1,natom
       do ifreq=1,green%nw
       call xsum_mpi(green%oper(ifreq)%matlu(iatom)%mat,spacecomm,ierr)
       enddo
     enddo
   endif ! opt_ksloc=2
 endif ! opt_tw==-1
 ABI_DEALLOCATE(fw)
 ABI_DEALLOCATE(ft)
 call destroy_green_tau(green_temp)
 call destroy_green(green_temp)

end subroutine fourier_green
!!***


!!****f* m_green/check_fourier_green
!! NAME
!! check_fourier_green
!!
!! FUNCTION
!!  Check fourier transformations
!!
!! INPUTS
!!  cryst_struc <type(crystal_structure)>= crystal structure data.
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE

subroutine check_fourier_green(cryst_struc,green,mpi_enreg,paw_dmft,pawang)

 use defs_basis
 use defs_abitypes
 use m_crystal, only : crystal_structure
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_fourier_green'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(crystal_structure),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(paw_dmft_type), intent(inout) :: paw_dmft

!local variables-------------------------------
 type(green_type) :: green_check
 character(len=500) :: message
! *********************************************************************
! Only imaginary frequencies here
 if(green%w_type=="real") then
   write(message,'(2a,i3,13x,a)') ch10,'   BUG: check_fourier_green not implemented for real frequency'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 endif

 call init_green(green_check,paw_dmft)
 call init_green_tau(green_check,paw_dmft)
 call copy_green(green,green_check,opt_tw=2)

 call fourier_green(cryst_struc,green_check,mpi_enreg,paw_dmft,&
& pawang,opt_ksloc=2,opt_tw=-1)

 write(message,'(3a)') ch10,' === Print (for check by user) of occupation matrix'&
&   ,' after  fourier transform with respect to initial one'
 call wrtout(std_out,message,'COLL')
 call printocc_green(green_check,6,paw_dmft,3)

 call fourier_green(cryst_struc,green_check,mpi_enreg,paw_dmft,&
& pawang,opt_ksloc=2,opt_tw=1)

 call integrate_green(cryst_struc,green_check,mpi_enreg,paw_dmft,&
& pawang,prtopt=2,opt_ksloc=2)

 write(message,'(3a)') ch10,' === Print (for check by user) of occupation matrix'&
&   ,' after double fourier transform with respect to initial one'
 call wrtout(std_out,message,'COLL')
 call printocc_green(green_check,5,paw_dmft,3)

 call destroy_green_tau(green_check)
 call destroy_green(green_check)

end subroutine check_fourier_green
!!***


!!****m* m_green/int_fct
!! NAME
!! int_fct
!!
!! FUNCTION
!! Do integration in matsubara space
!!
!! COPYRIGHT
!! Copyright (C) 2006-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ff= function is frequency space
!!  ldiag    = option according to diagonal or non-diagonal elements
!!  option = nspinor
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!  integral = integral of ff over matsubara frequencies
!!
!! SIDE EFFECTS
!!  ft= function is time space
!!
!! PARENTS
!!      m_green
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE
subroutine int_fct(ff,ldiag,option,paw_dmft,integral)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int_fct'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 logical,intent(in) :: ldiag
 integer,intent(in) :: option
 complex(dpc),intent(out) :: integral
 type(paw_dmft_type), intent(in) :: paw_dmft
 complex(dpc), intent(in) :: ff(paw_dmft%dmft_nwlo)

!local variables-------------------------------
 integer :: ifreq
! *********************************************************************

 integral=czero
 if(ldiag) then

  if(option==1) then ! nspinor==1
    do ifreq=1,paw_dmft%dmft_nwlo
!   write(500,*) paw_dmft%omega_lo(ifreq),real(ff(ifreq)),imag(ff(ifreq))
     integral=integral+2.d0*paw_dmft%temp *                         &
&      real( ff(ifreq)-one / ( j_dpc*paw_dmft%omega_lo(ifreq) ) ) *  &
&      paw_dmft%wgt_wlo(ifreq)
    enddo
    integral=integral+half
  endif

  if(option==2) then ! nspinor==2
    do ifreq=1,paw_dmft%dmft_nwlo
     integral=integral+2.d0*paw_dmft%temp *                         &
&          ( ff(ifreq)-one / ( j_dpc*paw_dmft%omega_lo(ifreq) ) ) *  &
&      paw_dmft%wgt_wlo(ifreq)
    enddo
    integral=integral+half
  endif


 else   ! ldiag

! write(std_out,*) "nondiag"
  if(option==1) then
    do ifreq=1,paw_dmft%dmft_nwlo
     integral=integral+2.d0*paw_dmft%temp *   &
&     real( ff(ifreq) ) *                     &
&     paw_dmft%wgt_wlo(ifreq)
    enddo
  endif
  if(option==2) then
    do ifreq=1,paw_dmft%dmft_nwlo
     integral=integral+2.d0*paw_dmft%temp *   &
&           ff(ifreq)   *                     &
&     paw_dmft%wgt_wlo(ifreq)
    enddo
  endif
 endif  ! ldiag

end subroutine int_fct
!!***

!!****m* m_green/fourier_fct
!! NAME
!! fourier_fct
!!
!! FUNCTION
!! Do fourier transformation from matsubara space to imaginary time
!! (A spline is performed )
!!
!! COPYRIGHT
!! Copyright (C) 2006-2009 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ldiag    = option according to diagonal or non-diagonal elements
!!  opt_four = option for direct (1) or inverse (-1) fourier transform
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  fw= function is frequency space
!!  ft= function is time space
!!
!! PARENTS
!!      local_ks_green,m_green
!!
!! CHILDREN
!!      invfourier,leave_new,nfourier,spline_c,wrtout
!!
!! SOURCE
subroutine fourier_fct(fw,ft,ldiag,ltau,opt_four,paw_dmft)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourier_fct'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 logical,intent(in) :: ldiag
 integer,intent(in) :: ltau,opt_four
 type(paw_dmft_type), intent(in) :: paw_dmft
 complex(dpc), intent(inout) :: fw(paw_dmft%dmft_nwlo)
 complex(dpc), intent(inout) :: ft(ltau)

!local variables-------------------------------
 complex(dpc), allocatable ::  splined_li(:)
 complex(dpc), allocatable ::  tospline_li(:)
! complex(dpc), allocatable :: fw1(:)
 real(dp), allocatable :: ftr(:)
 real(dp) :: beta
 integer :: iflag,itau,iwarn
 character(len=500) :: message
! *********************************************************************
 beta=one/paw_dmft%temp
 iflag=0
 if(ldiag) iflag=1

! == inverse fourier transform
 if(opt_four==-1) then

   ABI_ALLOCATE(splined_li,(paw_dmft%dmft_nwli))
!   allocate(fw1(0:paw_dmft%dmft_nwlo-1))
   call spline_c(paw_dmft%dmft_nwlo,paw_dmft%dmft_nwli,paw_dmft%omega_lo,&
&                 paw_dmft%omega_li,splined_li,fw)
   call invfourier(splined_li,ft,paw_dmft%dmft_nwli-1,ltau,iflag,beta)
!   deallocate(fw1)
   ABI_DEALLOCATE(splined_li)

! == direct fourier transform
 else if(opt_four==1) then

   ABI_ALLOCATE(ftr,(ltau))

   iwarn=0
   do itau=1,ltau
     if(abs(imag(ft(itau)))>tol12) then
       if(ldiag) then
         write(message,'(a,a,2(e15.4))') ch10,&
&          "ERROR:  green function is not real in imaginary time space",ft(itau)
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       else
         iwarn=iwarn+1
         ftr(itau)=real(ft(itau))
       endif
     else
       ftr(itau)=real(ft(itau))
     endif
!       write(std_out,*) itau,ftr(itau)
   enddo

   ABI_ALLOCATE(tospline_li,(paw_dmft%dmft_nwli))
   call nfourier(ftr,tospline_li,iflag,paw_dmft%dmft_nwli-1,ltau,beta)
!   do ifreq=1,paw_dmft%dmft_nwli
!    write(std_out,*) ifreq,tospline_li(ifreq)
!   enddo
   call spline_c(paw_dmft%dmft_nwli,paw_dmft%dmft_nwlo,paw_dmft%omega_li,&
&                 paw_dmft%omega_lo,fw,tospline_li)
   ABI_DEALLOCATE(tospline_li)

   ABI_DEALLOCATE(ftr)
   if(iwarn>0) then
     write(message,'(a,a,2(e15.4))') ch10,&
&     "WARNING:  green function is not real in imaginary time space"
     call wrtout(std_out,message,'COLL')
   endif

 endif

end subroutine fourier_fct

END MODULE m_green
!!***
