!{\src2tex{textfont=tt}}
!!****f* ABINIT/fftwfn
!! NAME
!! fftwfn
!!
!! FUNCTION
!! Calculate wavefunctions in real space using FFT
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! igfft(npwwfn)=index of each plane wave in FFT grid
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkibz=number of k points
!! npwwfn=number of plane waves
!! nsppol=number of independent spin polarizations 
!! tim_fourdp=4 if called from within screening ; =5 if called from within sigma
!! wfg(npwwfn,my_minb:my_maxb,nkibz,nsppol)=wavefunctions in reciprocal space treated by this processor.
!! my_minb,my_maxb = min and max band treated by this processor
!! mpi_enreg= datatype containing information on parallelism to be passed to fourdp
!!
!! OUTPUT
!!  wfr(ngfft(1)*ngfft(2)*ngfft(3),my_minb:my_maxb,nkibz,nsppol) 
!!   wavefunctions in real space, for each band, k point and spin
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!      fourdp,wrtout,xcomm_init,xmaster_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fftwfn(paral_kgb,npwwfn,my_minb,my_maxb,nkibz,nsppol,wfg,wfr,igfft,ngfft,tim_fourdp,mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftwfn'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_maxb,my_minb,nkibz,npwwfn,nsppol,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: igfft(npwwfn),ngfft(18)
 complex(gwpc),intent(in) :: wfg(npwwfn,my_minb:my_maxb,nkibz,nsppol)
 complex(gwpc),intent(out) :: wfr(ngfft(1)*ngfft(2)*ngfft(3),my_minb:my_maxb,nkibz,nsppol)

!Local variables-------------------------------
 character(len=50),parameter :: sub_name='fftwfn.F90'
!scalars
 integer :: ib,ig,ik,is,istat,master,me,nfftot,spaceComm
 character(len=500) :: message
!arrays
 real(dp),allocatable :: wfg_dp(:,:),wfr_dp(:,:)

! *************************************************************************

 write(message,'(a)')' FFT wavefunctions to real space...'
 call wrtout(std_out,message,'COLL')

 call xcomm_init  (mpi_enreg,spaceComm) 
 call xme_init    (mpi_enreg,me)         
 call xmaster_init(mpi_enreg,master)  

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ABI_ALLOCATE(wfg_dp,(2,nfftot))
 ABI_ALLOCATE(wfr_dp,(2,nfftot))
 istat = ABI_ALLOC_STAT
 if (istat/=0) stop "out-of-memory" 

 do is=1,nsppol
   do ik=1,nkibz
     do ib=my_minb,my_maxb
!      Fill FFT array from PW array
       wfg_dp(:,:)=zero
       do ig=1,npwwfn
         wfg_dp(1,igfft(ig))=real (wfg(ig,ib,ik,is))
         wfg_dp(2,igfft(ig))=aimag(wfg(ig,ib,ik,is))
       end do
!      
!      Take FFT to give wfn in real space
       call fourdp(2,wfg_dp(:,:),wfr_dp(:,:),+1,mpi_enreg,nfftot,ngfft,paral_kgb,tim_fourdp)

       wfr(:,ib,ik,is)=cmplx(wfr_dp(1,:),wfr_dp(2,:),gwpc)

     end do 
   end do 
 end do
 ABI_DEALLOCATE(wfg_dp)
 ABI_DEALLOCATE(wfr_dp)

end subroutine fftwfn
!!***
