!{\src2tex{textfont=tt}}
!!****f* ABINIT/cvxclda
!! NAME
!! cvxclda
!!
!! FUNCTION
!! Calculate Vxc on the FFT grid.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, YMN, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!! ixc = choice for the exchange-correlation potential.
!! mpi_enreg = informations about MPI parallelization.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nfftot = total number of points on the FFT grid.
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! rho(nfftot,nsppol) = the charge density on the FFT grid.
!!  (total in first half and spin-up in second half if nsppol=2)
!! rprimd(3,3) = dimensional real space primitive translations (bohr).
!!
!! OUTPUT
!!  vxclda(nfftot,nsppol) = the exchange-correlation potential on the FFT grid.
!!                      (spin up in first half and spin down in second half if nsppol=2)
!!
!! NOTES
!! No xc quadrature
!! No nl core correction
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!      dtsetcopy,dtsetfree,leave_new,rhohxc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine cvxclda(dtset,ixc,mpi_enreg,ngfft,nfftot,nsppol,rho,rprimd,vxclda)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cvxclda'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_53_abiutil
 use interfaces_56_xc
 use interfaces_93_rdm, except_this_one => cvxclda
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,nfftot,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rho(nfftot,nsppol),rprimd(3,3)
 real(dp),intent(out) :: vxclda(nfftot,nsppol)

!Local variables ------------------------------
!scalars
 integer,parameter :: nkxc=0,option=0
 integer :: ir,n3xccc,nk3xc,nspden
 real(dp) :: enxc,gsqcut,vxcavg
 character(len=500) :: message
 type(dataset_type) :: dtGW
!arrays
 real(dp) :: strsxc(6)
 real(dp),allocatable :: dum1(:,:),dum2(:,:,:),kxc(:,:),rhog(:,:),vhartr(:)
 real(dp),allocatable :: xccc3d(:)

!************************************************************************

 write(message,'(a,i3)')' cvxclda: calculating Vxc using ixc = ',ixc
 call wrtout(std_out,message,'COLL')
!
!Form Vxc (in Hartree units)
!
 if (ixc==0) then
!  For backward compatibility
   write (message,'(6a)')ch10,&
&   ' cvxclda: WARNING - ',ch10,&
&   ' ixc = 0 is a relativistic Ceperley-Alder xc functional [PRB 26, p. 4199, (1982)]',ch10,&
&   ' in the GW code. It should be used only for backward compatibility.'
   call wrtout(ab_out,message,'COLL') 
   call wrtout(std_out,message,'COLL')
   if (nsppol==2) then 
     write (message,'(3a)')&
&     ' cvxclda: ERROR- ',ch10,&
&     ' ixc = 0 and nsppol==2 not yet implemented ' 
     call wrtout(std_out,message,'COLL') 
     call leave_new('COLL')
   end if 
   do ir=1,nfftot
     vxclda(ir,1)=vxcca(rho(ir,1))
   end do
 else
!  MG this comes from my private branch (5.3.4)
!  if ((ixc<0).or.(ixc>16)) then 
!  this comes from XG merge-public--5.4.3--patch-2
   if (ixc>=10 .and. dtset%xclevel/=2) then
     write (message,'(4a,i3,a)') ch10,&
&     ' cvxclda: ERROR - ',ch10,&
&     '  ixc = ',ixc,' is not allowed at the present time in the GW code.'
     call wrtout(ab_out,message,'COLL') 
     call leave_new('COLL')
   end if
!  
!  Copy the input variables from the current dataset to a temporary one to tune some parameters
   call dtsetCopy(dtGW,dtset)
   dtGW%intxc=0
   write(message,'(a)')&
&   ' cvxclda: calling rhohxc to calculate Vxc[n_val] (excluding non-linear core corrections)'
   call wrtout(std_out,message,'COLL')
!  
!  Note: one must have nfftot=ngfft1*ngfft2*ngfft3, ie the FFT grid must not 
!  be augmented. This is actually enforced at the present time in setmesh.f
!  
   ABI_ALLOCATE(rhog,(2,nfftot))
   ABI_ALLOCATE(vhartr,(nfftot))
   ABI_ALLOCATE(kxc,(nfftot,nkxc))
!  gsqcut and rhog are zeroed because they are not used by rhohxc if 1<=ixc<=16 and option=0
   gsqcut=zero ; rhog(:,:)=zero
!  
!  TODO this is the 3D core electron density for XC core correction (bohr^-3)
!  should implement the non linear core correction 
   n3xccc=0       
   ABI_ALLOCATE(xccc3d,(n3xccc))
!  
!  nkxc=0  ==> no computation of the exchange-correlation kernel
!  option=0  ==> only exc, vxc, strsxc
!  
   nspden=nsppol
   ABI_ALLOCATE(dum1,(nfftot,0))
   ABI_ALLOCATE(dum2,(nfftot,nspden,0))

!  to be adjusted for the call to rhohxc
   nk3xc=1

   call rhohxc(dtGW,enxc,gsqcut,0,kxc,mpi_enreg,nfftot,ngfft,dum1,0,dum2,0,nkxc,nk3xc,nspden,&
&   n3xccc,option,rhog,rho,rprimd,strsxc,1,vhartr,vxclda,vxcavg,xccc3d)
   ABI_DEALLOCATE(dum1)
   ABI_DEALLOCATE(dum2)

   ABI_DEALLOCATE(rhog)
   ABI_DEALLOCATE(vhartr)
   ABI_DEALLOCATE(xccc3d)
   ABI_DEALLOCATE(kxc)
   call dtsetFree(dtGW)

   write(message,'(a,f8.4,2a,f8.4,2a)')&
&   ' cvxclda: rhohxc returned  Exc[n_val]  = ',enxc,  ' [Ha]',&
&   ' and <Vxc[n_val]> = ',vxcavg,' [Ha]',ch10
   call wrtout(std_out,message,'COLL')
 end if

end subroutine cvxclda
!!***
