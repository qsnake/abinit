!{\src2tex{textfont=tt}}
!!****f* ABINIT/mag_out
!! NAME
!! mag_out
!!
!! FUNCTION
!! This routine prints the magnetization.
!!
!! COPYRIGHT
!!
!! INPUTS
!!  dtefield <type(efield_type)> = variables related to Berry phase. 
!!  mpi_enreg=informations about MPI parallelization
!!  nkpt = number of kpoints
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  wtk = weights of k points
!!
!! OUTPUT
!!  only printing
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      matr3inv,wrtout,xcomm_world,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mag_out(dtefield,mpi_enreg,nkpt,rprimd,wtk)

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
#define ABI_FUNC 'mag_out'
 use interfaces_14_hidewrite
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
 type(efield_type),intent(in) :: dtefield
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: rprimd(3,3),wtk(nkpt)

!Local variables-------------------------------
!scalars
 integer :: idir,ierr,jdir,ikpt,mu
 integer :: spaceComm
 real(dp) :: intfac1, intfac2
!arrays
 character(len=500) :: message
 real(dp) :: bcart(3),chern(2,3), gprimd(3,3)
 real(dp) :: mag_loc(3),  mag_bare(2,3)

! ************************************************************************

 call xcomm_world(mpi_enreg,spaceComm)

 call matr3inv(rprimd,gprimd)

 mag_loc(:) = zero
 mag_bare(:,:) = zero
 chern(:,:) = zero

 intfac1 = 1.0/(Sp_Lt)
 intfac2 = 1.0/(2.0*Sp_Lt*two_pi**3)

 do idir = 1, 3
   do mu=1,3
     bcart(mu)=dot_product(dtefield%dkvecs(:,idir),gprimd(mu,:))
   end do
   do ikpt = 1, nkpt
     mag_loc(idir) = mag_loc(idir) - wtk(ikpt)*dtefield%mag_local_k(idir,ikpt)
     do jdir = 1, 3
       chern(1,jdir) = chern(1,jdir) - wtk(ikpt)*bcart(jdir)*dtefield%chern_k(2,ikpt,idir)
       chern(2,jdir) = chern(2,jdir) + wtk(ikpt)*bcart(jdir)*dtefield%chern_k(1,ikpt,idir)
       mag_bare(1,jdir) = mag_bare(1,jdir) - bcart(jdir)*wtk(ikpt)*dtefield%mag_k(2,ikpt,idir)
       mag_bare(2,jdir) = mag_bare(2,jdir) + bcart(jdir)*wtk(ikpt)*dtefield%mag_k(1,ikpt,idir)
     end do
   end do ! end loop over ikpt
 end do ! end loop over idir
 
 mag_loc(:) = mag_loc(:)*intfac1
 call xsum_mpi(mag_loc,spaceComm,ierr)

 mag_bare(:,:) = mag_bare(:,:)*intfac2
 call xsum_mpi(mag_bare,spaceComm,ierr)
!mag_tot(:) = mag_loc(:) + mag_mmat(:)

 chern(:,:) = chern(:,:)/two_pi !
 call xsum_mpi(chern,spaceComm,ierr)
!two_pi comes from the definition and the rest on going from the integral to the sum

 write(message,'(2a)')ch10,' ********************************************** '
 call wrtout(ab_out,message,'COLL')
 write(message,'(a)')'  Results of orbital magnetization calculation '
 call wrtout(ab_out,message,'COLL')
 write(message,'(2a)')' ********************************************** ',ch10
 call wrtout(ab_out,message,'COLL')

 write(message,'(a)')' On-site local magnetization (mag_loc): '  
 call wrtout(ab_out,message,'COLL')
 do idir = 1, 3
   write(message,'(a,i4,a,es16.8)')'   dir = ',idir,'  mag_loc = ',mag_loc(idir)
   call wrtout(ab_out,message,'COLL')
 end do 
 
 write(message,'(a)')' Bare magnetization (mag_bare): '  
 call wrtout(ab_out,message,'COLL')
 do idir = 1, 3
   write(message,'(a,i4,a,es16.8,a,es16.8)')'   dir = ',idir,'  mag_bare = ',mag_bare(1,idir),' + i',mag_bare(2,idir)
   call wrtout(ab_out,message,'COLL')
 end do 

 write(message,'(a)')' Total magnetization (mag_cart): '  
 call wrtout(ab_out,message,'COLL')
 do idir = 1, 3
   write(message,'(a,i4,a,es16.8)')'   dir = ',idir,'  mag_cart = ',dtefield%mag_cart(idir)
   call wrtout(ab_out,message,'COLL')
 end do 

 write(message,'(a)')' Chern number: '  
 call wrtout(ab_out,message,'COLL')
 do idir = 1, 3
   write(message,'(a,i4,a,es16.8,a,es16.8)')'   dir = ',idir,'  Chern = ',chern(1,idir),' +i',chern(2,idir)
   call wrtout(ab_out,message,'COLL')
 end do 

 write(message,'(2a)')' ********************************************** ',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine mag_out
!!***
