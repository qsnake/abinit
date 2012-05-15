!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_newvtr
!! NAME
!! wvl_newvtr
!!
!! FUNCTION
!! Compute new trial potential by summing all components.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=informations about MPI parallelization
!!  vhartr(nfft)=array for holding Hartree potential
!!  vpsp(nfft)=array for holding local psp
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!
!! OUTPUT
!!  nele=number of relevant elements in vtrial per spin component
!!  offset=offset for first relevant element in vtrial
!!  vtrial(nfft,nspden)=new potential
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop,rhotov,scfcv,setvtr
!!
!! CHILDREN
!!      leave_new,wrtout,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_newvtr(dtset, mpi_enreg, nele, offset, vhartr, vpsp, vtrial, vxc, wvl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_newvtr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments-------------------------------
!scalars
 integer,intent(out) :: nele,offset
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 real(dp),intent(in) :: vhartr(dtset%nfft),vpsp(dtset%nfft)
 real(dp),intent(in) :: vxc(dtset%nfft*dtset%nspden)
 real(dp),intent(out) :: vtrial(dtset%nfft*dtset%nspden)

!Local variables-------------------------------
!scalars
 integer :: comm,i,ispden,me
 character(len=500) :: message

! *************************************************************************

 if (dtset%usewvl /= 1) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' wvl_newvtr : BUG -',ch10,&
&   '  dtset%usewvl /= 1 not allowed (use newtr() instead) !'
   call wrtout(std_out,message,'COLL')
   call leave_new('PERS')
 end if

 call xcomm_world(mpi_enreg,comm,myrank=me)

#if defined HAVE_DFT_BIGDFT
 offset = wvl%Glr%d%n1i * wvl%Glr%d%n2i * &
& mpi_enreg%nscatterarr(me, 4)
 nEle   = wvl%Glr%d%n1i * wvl%Glr%d%n2i * &
& mpi_enreg%nscatterarr(me, 2)

!Vtrial has a beginning offset, but no one for spin dimensions
![dpSize(1) * dpSize(2) * nscatterarr(me, 2) * nspden]
!Vpsp   has an offset in GGA, [dpSize(1) * dpSize(2) * nscatterarr(me, 2)]
!Vhartr has an offset in GGA, [dpSize(1) * dpSize(2) * nscatterarr(me, 2)]
!Vxc    has no offset in GGA for each spin dimension and is packed,
![dpSize(1) * dpSize(2) * nscatterarr(me, 2) * nspden]
 do ispden = 1, dtset%nspden, 1
   do i = 1, nEle, 1
     vtrial(i + (ispden - 1) * nEle) = vhartr(i + offset) + &
&     vpsp(i + offset) + vxc(i + (ispden - 1) * nEle)
   end do
 end do
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_newvtr: BUG -',ch10,&
& '  BigDFT is not compiled. Use --enable-bigdft during configure.'
 call wrtout(std_out, message, 'COLL')
 call leave_new('COLL')
#endif

end subroutine wvl_newvtr
!!***
