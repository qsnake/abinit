!{\src2tex{textfont=tt}}
!!****f* ABINIT/bfactor
!! NAME
!! bfactor
!!
!! FUNCTION
!! Calculate the nesting factor
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nkptfull = number of k-points in full grid
!!  kptfull(3,nkptfull) = k-point grid
!!  nqpt = number of qpoints
!!  qpt(3,nqpt) = q-point grid (must be a subgrid of the k grid),
!!                the nesting factor will be calculated for each q point in this array
!!  nkpt = eventually reduced number of k-points
!!  weight(nband,nkpt) =  integration weights for each k-point and band (NOT NORMALISED!!!)
!!  nband = number of bands
!!
!! OUTPUT
!!  nestfactor(nqpt) = array containing the nesting factor values
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Inspired to nmsq_gam_sumfs and mkqptequiv
!!  TODO : better use of symmetries to reduce the computational effort
!! Must be called with kpt = full grid! Reduction by symmetry is not possible for q-dependent quantities (or not easy :)
!!
!! PARENTS
!!      mknesting,outelph
!!
!! CHILDREN
!!      get_rank_1kpt,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine bfactor(nkptfull,kptfull,nqpt,qpt,kptrank_t,nkpt,weight,nband,nestfactor)

 use m_profiling

 use defs_basis
 use m_kptrank
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bfactor'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband,nkptfull,nqpt,nkpt
!arrays
 real(dp),intent(in) :: kptfull(3,nkptfull),qpt(3,nqpt),weight(nband,nkpt)
 real(dp),intent(inout) :: nestfactor(nqpt)
 type(kptrank_type), intent(in) :: kptrank_t

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2,ikplusq_irr,ikpt
 integer :: irank_kpt,ikpt_irr,iqpt,symrank_kpt
 real(dp) :: factor,w1,w2
 character(len=500) :: message
!arrays
 real(dp) :: kptpq(3)

! *************************************************************************

 nestfactor(:)=zero

 do iqpt=1,nqpt
   do ikpt=1,nkptfull
     call get_rank_1kpt (kptfull(:,ikpt),irank_kpt,kptrank_t)
     ikpt_irr = kptrank_t%invrank(irank_kpt)

     kptpq(:) = kptfull(:,ikpt) + qpt(:,iqpt)
     call get_rank_1kpt (kptpq,symrank_kpt,kptrank_t)

     ikplusq_irr = kptrank_t%invrank(symrank_kpt)
     if (ikplusq_irr == -1) then
       write (message,'(4a)')ch10,' bfactor : ERROR- :',ch10,' it looks like no kpoint equiv to k+q !!!'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     do ib1=1,nband
       w1 = weight(ib1,ikpt_irr) !weight for distance from the Fermi surface
       if (w1 < tol6 ) cycle
       do ib2=1,nband
         w2 = weight(ib2,ikplusq_irr) !weight for distance from the Fermi surface
         if (w1 < tol6 ) cycle
         nestfactor(iqpt) = nestfactor(iqpt) + w1*w2
       end do !ib2
     end do !ib1

   end do !ikpt
 end do !iqpt


!need prefactor of (1/nkptfull) for normalisation of integration
 factor=1./nkptfull
 nestfactor(:)=factor*nestfactor(:)

end subroutine bfactor
!!***
