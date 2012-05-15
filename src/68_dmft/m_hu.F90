!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hu
!! NAME
!!  m_hu
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

MODULE m_hu

 use m_profiling

 use defs_basis
 use defs_datatypes

 implicit none

 private 

 public :: init_hu
 public :: destroy_hu
! public :: qmc_hu
 public :: nullify_hu
 public :: print_hu


!!***

!!****t* m_hu/hu_type
!! NAME
!!  hu_type
!!
!! FUNCTION
!!  This structured datatype contains interaction matrices for the correlated subspace
!!
!! SOURCE

 type, public :: hu_type ! for each typat

  integer :: lpawu         

  real(dp) :: upawu    ! => upaw

  real(dp) :: jpawu    ! => jpaw

  real(dp), pointer :: vee(:,:,:,:) ! => vee

  real(dp), pointer :: uqmc(:)

  real(dp), pointer :: udens(:,:)

 end type hu_type

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_hu/init_hu
!! NAME
!! init_hu
!!
!! FUNCTION
!!  Allocate variables used in type hu_type.
!!
!! INPUTS
!!
!! OUTPUTS
!! hu  = structure of data for dmft of type hu_type
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_hu(cryst_struc,pawtab,hu)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_hu'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(crystal_structure),intent(in) :: cryst_struc
 type(pawtab_type), target, intent(in)  :: pawtab(cryst_struc%ntypat)
 type(hu_type), intent(inout) :: hu(cryst_struc%ntypat)
!Local variables ------------------------------------
 integer :: itypat,i,ij,ij1,ij2,j,lpawu,ms,ms1,m,m1,ndim
 integer, allocatable :: xij(:,:)
 real(dp) :: xtemp
 character(len=500) :: message
!************************************************************************
 write(message,'(2a)') ch10,"  == Compute Interactions for DMFT"
 call wrtout(std_out,message,'COLL')

 xtemp=zero
 call nullify_hu(hu,cryst_struc%ntypat)

! ====================================
!  Compute hu(iatom)%uqmc from vee
! ====================================
 do itypat=1,cryst_struc%ntypat
   lpawu=pawtab(itypat)%lpawu
   hu(itypat)%lpawu=lpawu
   if(lpawu.ne.-1) then
     hu(itypat)%upawu=pawtab(itypat)%upawu
     hu(itypat)%jpawu=pawtab(itypat)%jpawu
     ndim=2*lpawu+1
     write(message,'(2a,i4)')  ch10,'  -------> For Correlated Species', itypat
     call wrtout(std_out,  message,'COLL')
!     allocate(hu(itypat)%vee(ndim,ndim,ndim,ndim))
     ABI_ALLOCATE(hu(itypat)%uqmc,(ndim*(2*ndim-1)))
     ABI_ALLOCATE(hu(itypat)%udens,(2*ndim,2*ndim))
     ABI_ALLOCATE(xij,(2*ndim,2*ndim))
     hu(itypat)%vee => pawtab(itypat)%vee
     hu(itypat)%udens=zero
     ij=0
     do ms=1,2*ndim-1
         xij(ms,ms)=0
       do ms1=ms+1,2*ndim 
         ij=ij+1
         xij(ms,ms1)=ij
         xij(ms1,ms)=ij
         if(ms<=ndim.and.ms1>ndim) then
           m1 = ms1 - ndim
           m  = ms
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%vee(m,m1,m,m1)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         else if(ms<=ndim.and.ms1<=ndim) then
           m1 = ms1
           m  = ms
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)-hu(itypat)%vee(m,m1,m1,m)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%uqmc(ij)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         else
           m1 = ms1 - ndim
           m  = ms  - ndim
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)-hu(itypat)%vee(m,m1,m1,m)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%uqmc(ij)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         endif
       enddo
     enddo
     xij(2*ndim,2*ndim)=0
     write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation "
     call wrtout(std_out,  message,'COLL')
     write(message,'(6x,14(2x,i5))') (m,m=1,2*ndim)
     call wrtout(std_out,  message,'COLL')
!     xtemp1b=0.d0
! ====================================
!  Print hu(iatom)%uqmc 
! ====================================
     ij1=-10
     ij2=-10
     ij=0
     do i=1,2*ndim
       do j=i+1,2*ndim
         ij=ij+1
         if(j==i+1) ij1=ij
         if(j==2*ndim) ij2=ij
       enddo 
!       write(std_out,*) itypat
!       do m=1,i
!        write(std_out,*) i,m
!        write(std_out,*) xij(i,m)
!        write(std_out,*) ij1,ij2
!       enddo
       if(i==1)               write(message,'(i3,14f7.3)') &
&                              i,xtemp, (hu(itypat)%uqmc(m),m=ij1,ij2)
       if(i/=2*ndim.and.i/=1) write(message,'(i3,14f7.3)') i, &
&        (hu(itypat)%uqmc(xij(i,m)), m=1,i-1),xtemp, (hu(itypat)%uqmc(m),m=ij1,ij2)
       if(i==2*ndim)          write(message,'(i3,14f7.3)') i, &
&                  (hu(itypat)%uqmc(xij(i,m)), m=1,i-1),xtemp
       call wrtout(std_out,  message,'COLL')
     enddo 
       write(message,'(5x,a)') "--------------------------------------------------------"
       call wrtout(std_out,  message,'COLL')
     ABI_DEALLOCATE(xij)
   else
     hu(itypat)%upawu=zero
     hu(itypat)%jpawu=zero
!     allocate(hu(itypat)%vee(0,0,0,0))
   endif
 enddo ! itypat

end subroutine init_hu
!!***

!!****f* m_hu/nullify_hu
!! NAME
!! nullify_hu
!!
!! FUNCTION
!!  nullify hu
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine nullify_hu(hu,ntypat)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_hu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntypat
 type(hu_type),intent(inout) :: hu(ntypat)
!Local variables-------------------------------
 integer :: itypat

!*********************************************************************

 do itypat=1,ntypat
  nullify(hu(itypat)%vee)
  nullify(hu(itypat)%uqmc)
  nullify(hu(itypat)%udens)
 enddo


end subroutine nullify_hu
!!***

!!****f* m_hu/destroy_hu
!! NAME
!! destroy_mh
!!
!! FUNCTION
!!  deallocate hu
!!
!! INPUTS
!!  hu
!!
!! OUTPUT
!!
!! PARENTS
!!  
!!
!! CHILDREN
!!   wrtout
!!
!! SOURCE

subroutine destroy_hu(hu,ntypat)

 use defs_basis
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_hu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntypat
 type(hu_type),intent(inout) :: hu(ntypat)

!Local variables-------------------------------
 integer :: itypat

! *********************************************************************

 do itypat=1,ntypat
!  if ( associated(hu(itypat)%vee) )  deallocate(hu(itypat)%vee)
  if ( associated(hu(itypat)%uqmc) )   then
    ABI_DEALLOCATE(hu(itypat)%uqmc)
  end if
  if ( associated(hu(itypat)%udens) )   then
    ABI_DEALLOCATE(hu(itypat)%udens)
  end if
 enddo

end subroutine destroy_hu
!!***

!!****f* m_hu/print_hu
!! NAME
!! print_hu
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
!!      wrtout
!!
!! SOURCE

subroutine print_hu(hu,ntypat,prtopt)

 use defs_basis
 use m_crystal, only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_hu'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in):: ntypat
 type(hu_type),intent(in) :: hu(ntypat)
 integer :: prtopt

!Local variables-------------------------------
 integer :: itypat,lpawu
 character(len=500) :: message
! *********************************************************************

 if(prtopt>0) then
 endif
 do itypat = 1 , ntypat
  lpawu=hu(itypat)%lpawu
  if(lpawu/=-1) then
   write(message,'(2a,i4)')  ch10,'  -------> For Correlated species'
   call wrtout(std_out,  message,'COLL')
  endif ! lpawu/=1
 enddo ! ntypat


end subroutine print_hu

END MODULE m_hu
!!***
