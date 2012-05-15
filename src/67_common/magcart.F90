!{\src2tex{textfont=tt}}
!!****f* ABINIT/magcart
!! NAME
!! magcart
!!
!! FUNCTION
!! This routine outputs the orbital magnetization in cartesian coordinates and
!! reduced coordinates.
!!
!! COPYRIGHT
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!  nkpt = number of kpoints
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  wtk = weights of k points
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtefield%mag_cart total magnetization along cartesian real space axes
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      matr3inv,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine magcart(dtefield,mpi_enreg,nkpt,rprimd,wtk)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_efield

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'magcart'
 use interfaces_32_util
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 type(efield_type),intent(inout) :: dtefield
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: rprimd(3,3),wtk(nkpt)

!Local variables-------------------------------
!scalars
 integer :: idir,ierr,ikpt,jdir,mu
 integer :: spaceComm
 real(dp) :: intfac1,intfac2
!arrays
 real(dp) :: bcart(3),gprimd(3,3),mag_loc(3),mag_bare(2,3)

! ************************************************************************

 call xcomm_init(mpi_enreg,spaceComm)

 call matr3inv(rprimd,gprimd)

 mag_loc(:) = zero
 mag_bare(:,:) = zero
 intfac1 = 1.0/(Sp_Lt)
 intfac2 = 1.0/(2.0*Sp_Lt*two_pi**3)
 
 do idir = 1, 3
   do mu=1,3
     bcart(mu)=dot_product(dtefield%dkvecs(:,idir),gprimd(mu,:))
   end do
   do ikpt = 1, nkpt
     mag_loc(idir) = mag_loc(idir) - wtk(ikpt)*dtefield%mag_local_k(idir,ikpt)
     do jdir = 1, 3
       mag_bare(1,jdir) = mag_bare(1,jdir) - bcart(jdir)*wtk(ikpt)*dtefield%mag_k(2,ikpt,idir)
       mag_bare(2,jdir) = mag_bare(2,jdir) + bcart(jdir)*wtk(ikpt)*dtefield%mag_k(1,ikpt,idir)
     end do
   end do ! end loop over ikpt
 end do ! end loop over idir

 mag_loc(:) = mag_loc(:)*intfac1
 call xsum_mpi(mag_loc,spaceComm,ierr)

 mag_bare(:,:) = mag_bare(:,:)*intfac2
 call xsum_mpi(mag_bare,spaceComm,ierr)
 
 dtefield%mag_cart(:) = mag_loc(:) + mag_bare(1,:)


end subroutine magcart
!!***
