!{\src2tex{textfont=tt}}
!!****f* ABINIT/mean_fftr
!! NAME
!! mean_fftr
!!
!! FUNCTION
!!  Compute the mean of an arraysp(nfft,nspden), over the
!!  FFT grid, for each component nspden, and return it
!!  in meansp(nspden). Take into account the spread
!!  of the array due to parallelism : the actual number of fft
!!  points is nfftot, but the number of points on this proc
!!  is nfft only.
!!  So : for ispden from 1 to nspden
!!       meansp(ispden) = sum(ifft=1,nfftot) arraysp(ifft,ispden) / nfftot
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  arraysp(nfft,nspden)=the array whose average has to be computed
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=number of FFT points stored by one proc
!!  nfftot=total number of FFT points
!!  nspden=number of spin-density components
!!
!! OUTPUT
!!  meansp(nspden)=mean value for each nspden component
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      fresid,multipoles_fftr,newvtr,pawmknhat,prcref,prcref_PMA,prctfvw1
!!      prctfvw2,prctfw3,psolver_rhohxc,rhohxc,rhohxcpositron,rhotov
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,leave_new,timab,wrtout
!!      xcomm_world,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mean_fftr(arraysp,meansp,mpi_enreg,nfft,nfftot,nspden)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mean_fftr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
#if defined DEBUG_CONTRACT
 use interfaces_32_contract
#endif
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nfftot,nspden
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: arraysp(nfft,nspden)
 real(dp),intent(out) :: meansp(nspden)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ispden,old_paral_level,spaceComm
 real(dp) :: invnfftot,tmean
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
!no_abirules
#if defined DEBUG_CONTRACT
 character(len=9) :: subrnm
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' mean_fftr : enter '
!ENDDEBUG

#if defined DEBUG_CONTRACT
 subrnm='mean_fftr'
 call contract_int_ge_val(subrnm,'nfft',nfft,1)
 call contract_int_ge_val(subrnm,'nfftot',nfftot,1)
 call contract_int_list(subrnm,'nspden',nspden,(/1,2,4/),3)
#endif



 if(nspden/=1 .and. nspden/=2 .and. nspden/=4)then
   write(message,'(a,a,a,a,a,a,i6)') ch10,&
&   ' mean_fftr: BUG -',ch10,&
&   '  The argument nspden should be 1, 2 or 4,',ch10,&
&   '  however, nspden=',nspden
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(nfft<0)then
   write(message,'(a,a,a,a,a,a,i6)') ch10,&
&   ' mean_fftr: BUG -',ch10,&
&   '  The argument nfft should be a positive number,',ch10,&
&   '  however, nfft=',nfft
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 invnfftot=one/(dble(nfftot))

 do ispden=1,nspden
   tmean=zero
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&REDUCTION(+:tmean) &
!  $OMP&SHARED(nfft,arraysp,ispden)
   do ifft=1,nfft
     tmean=tmean+arraysp(ifft,ispden)
   end do ! ifft
!  $OMP END PARALLEL DO
   meansp(ispden)=tmean*invnfftot
 end do ! ispden

!XG030514 : MPIWF The values of meansp(ispden) should
!now be summed accross processors in the same WF group,
!and spread on all procs.

!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then

   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
!  call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
!  The same, use xcomm_init to get the right spaceComm.
   if(mpi_enreg%mode_para=='b')then
     spaceComm=mpi_enreg%comm_fft
     call timab(48,1,tsec)
     call xsum_mpi(meansp,nspden,spaceComm ,ierr)
     call timab(48,2,tsec)
   end if
   mpi_enreg%paral_level=old_paral_level
 else if (associated(mpi_enreg%nscatterarr)) then
   call xcomm_world(mpi_enreg,spaceComm)
   call timab(48,1,tsec)
   call xsum_mpi(meansp,nspden,spaceComm ,ierr)
   call timab(48,2,tsec)
 end if

!DEBUG
!write(std_out,*)' mean_fftr : exit'
!stop
!ENDDEBUG

end subroutine mean_fftr
!!***
