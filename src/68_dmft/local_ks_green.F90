!{\src2tex{textfont=tt}}
!!****f* ABINIT/local_ks_green
!! NAME
!! local_ks_green
!!
!! FUNCTION
!! Compute the sum over k-point of ks green function.
!! do the fourier transformation and print it
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
!!      fourier_fct,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine local_ks_green(green,paw_dmft,prtopt)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_crystal, only : crystal_structure
 use m_green, only : green_type,fourier_fct
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'local_ks_green'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in)  :: paw_dmft
 integer, intent(in) :: prtopt

!Local variables ------------------------------
 character(len=500) :: message
 integer :: iband,ifreq,ikpt,isppol,itau,lsub,ltau,mbandc,nkpt,nsppol
 character(len=1) :: tag_is
 character(len=fnlen) :: tmpfil
 integer,allocatable :: unitgreenlocks_arr(:)
 real(dp) :: beta
 real(dp), allocatable :: tau(:)
 complex(dpc), allocatable :: loc_ks(:,:,:)
 complex(dpc), allocatable :: loc_ks_tau(:,:,:),fw(:),ft(:)
!scalars
!************************************************************************
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 ltau=128
 ABI_ALLOCATE(tau,(ltau))
 do itau=1,ltau
   tau(itau)=float(itau-1)/float(ltau)/paw_dmft%temp
 end do
 beta=one/paw_dmft%temp

!Only imaginary frequencies here
 if(green%w_type=="real") then
   write(message,'(2a,i3,13x,a)') ch10,'   BUG: compute_energy not implemented for real frequency'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!=========================================
!Compute local band ks green function
!=========================================
 ABI_ALLOCATE(loc_ks,(nsppol,mbandc,paw_dmft%dmft_nwlo))
 loc_ks(:,:,:)=czero
 do isppol=1,nsppol
   do iband=1,mbandc
     do ifreq=1,paw_dmft%dmft_nwlo
       do ikpt=1,nkpt
         loc_ks(isppol,iband,ifreq)=loc_ks(isppol,iband,ifreq)+  &
&         green%oper(ifreq)%ks(isppol,ikpt,iband,iband)*paw_dmft%wtk(ikpt)
       end do
     end do
   end do
 end do

!=========================================
!Compute fourier transformation 
!=========================================

 ABI_ALLOCATE(loc_ks_tau,(nsppol,mbandc,ltau))
 ABI_ALLOCATE(fw,(paw_dmft%dmft_nwlo))
 ABI_ALLOCATE(ft,(ltau))
 loc_ks_tau(:,:,:)=czero
 do isppol=1,nsppol
   do iband=1,mbandc
     do ifreq=1,paw_dmft%dmft_nwlo
       fw(ifreq)=loc_ks(isppol,iband,ifreq)
     end do
     call fourier_fct(fw,ft,.true.,ltau,-1,paw_dmft) ! inverse fourier
     do itau=1,ltau
       loc_ks_tau(isppol,iband,itau)=ft(itau)
     end do
   end do
 end do
 ABI_DEALLOCATE(fw)
 ABI_DEALLOCATE(ft)
 do isppol=1,nsppol
   do iband=1,mbandc
     do itau=1,ltau
       loc_ks_tau(isppol,iband,itau)=(loc_ks_tau(isppol,iband,itau)+conjg(loc_ks_tau(isppol,iband,itau)))/two
     end do
   end do
 end do

!=========================================
!Print out ksloc green function
!=========================================
 if(abs(prtopt)==1) then
   ABI_ALLOCATE(unitgreenlocks_arr,(nsppol))
   do isppol=1,nsppol
     write(tag_is,'(i1)')isppol
     tmpfil = trim(paw_dmft%filapp)//'Gtau_locks_isppol'//tag_is
     write(message,'(3a)') ch10," == Print green function on file ",tmpfil
     call wrtout(std_out,message,'COLL')
     unitgreenlocks_arr(isppol)=500+isppol-1
     open (unit=unitgreenlocks_arr(isppol),file=trim(tmpfil),status='unknown',form='formatted')
     rewind(unitgreenlocks_arr(isppol))
     write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenlocks_arr(isppol)
     write(message,'(a,a)') ch10,"# New record : First 40 bands"
     call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
     do lsub=1,mbandc/40+1
       do itau=1, ltau
         write(message,'(2x,50(e10.3,2x))') tau(itau), &
&         (real(loc_ks_tau(isppol,iband,itau)),iband=40*(lsub-1)+1,min(40*lsub,mbandc))
         call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       end do
       write(message,'(2x,50(e10.3,2x))') beta, &
&       ((-one-real(loc_ks_tau(isppol,iband,1))),iband=40*(lsub-1)+1,min(40*lsub,mbandc))
       call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       if(40*lsub<mbandc) then
         write(message,'(a,a,i5,a,i5)')    &
&         ch10,"# Same record, Following bands : From ",    &
&         40*(lsub),"  to ",min(40*(lsub+1),mbandc)
         call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       end if
     end do
!    call flush(unitgreenlocks_arr(isppol))
   end do
   ABI_DEALLOCATE(unitgreenlocks_arr)
 end if

!Deallocations
 ABI_DEALLOCATE(loc_ks)
 ABI_DEALLOCATE(loc_ks_tau)
 ABI_DEALLOCATE(tau)

end subroutine local_ks_green
!!***
