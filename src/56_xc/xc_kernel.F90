!{\src2tex{textfont=tt}}
!!****f* ABINIT/xc_kernel
!! NAME
!! xc_kernel
!!
!! FUNCTION
!! Calculate the exchange-correlation kernel in reciprocal space.
!! Require density in real space on the *full* FFT mesh
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (Rhaltaf,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Dtset<dataset_type>=all input variables in this dataset
!! Cryst<Crystal_structure>=Info on the crystal structure.
!! ixc = choice for the exchange-correlation potential.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nfft_tot = Total number of points on the FFT grid.
!! nspden=Number of independent spin densities.
!! rhor(nfft_tot,nspden) = the charge density on the full FFT grid.
!!  (total in first half and spin-up in second half if nspden=2)
!! npw: the size of kernel matrix
!! dim_kxcg=dimension of the kernel.
!! comm=MPI communicator.
!! [dbg_mode]=Optional flag used to switch on the debug mode.
!!
!! OUTPUT
!!  FIXME: Why are we using nfft_tot instead of the G-sphere
!!  kxcg(nfft_tot,dim_kxcg) = the exchange-correlation potential on the FFT grid.
!!  warning: the kernel is not divided by the unit cell volume
!!
!! NOTES
!!  No xc quadrature
!!  No nl core correction
!!
!! PARENTS
!!      m_screen,screening,sigma
!!
!! CHILDREN
!!      dtsetcopy,dtsetfree,fourdp,initmpi_seq,mkvxc3,printxsf,rhohxc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine xc_kernel(Dtset,Cryst,ixc,ngfft,nfft_tot,nspden,rhor,npw,dim_kxcg,kxcg,gvec,comm,dbg_mode)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_crystal

 use m_io_tools,      only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xc_kernel'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
 use interfaces_53_ffts
 use interfaces_56_xc, except_this_one => xc_kernel
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,npw,nfft_tot,nspden,dim_kxcg,comm
 logical,optional,intent(in) :: dbg_mode
 type(Crystal_structure),intent(in) :: Cryst
 type(dataset_type),intent(in) :: Dtset
!arrays
 integer,intent(in) :: gvec(3,npw),ngfft(18)
 real(dp),intent(in) :: rhor(nfft_tot,nspden)
 complex(gwpc),intent(out) :: kxcg(nfft_tot,dim_kxcg)

!Local variables ------------------------------
!scalars
 integer,parameter :: paral_kgb0=0
 integer :: cplex,i1,i2,i3,ig,igp,iq,ir,n3xccc,ngfft1,ngfft2,izero
 integer :: ngfft3,nkxc,option,ikxc,nk3xc,my_rank,master,unt_dmp
 real(dp) :: enxc,expo,gpqx,gpqy,gpqz,gsqcut,vxcavg
 character(len=500) :: msg,fname
 type(dataset_type) :: DtGW
 type(MPI_type) :: MPI_enreg_seq 
!arrays
 real(dp) :: qphon(3),strsxc(6)
 real(dp),allocatable :: dum(:),kxcpw_g(:,:),kxcr(:,:),phas(:,:,:)
 real(dp),allocatable :: rhog(:,:),vhartr(:),kxcpw_r(:,:),vxclda(:,:)
 real(dp),allocatable :: xccc3d(:),xx(:,:)
 real(dp),allocatable :: my_kxcg(:,:)

!************************************************************************

 ABI_CHECK(Dtset%nsppol==1,'nsppol/=1 not coded')
 ABI_CHECK(Dtset%nspden==1,'nsppol/=1 not coded')
 ABI_CHECK(nfft_tot==PRODUCT(ngfft(1:3)),"mismatch in nfftot")

!Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

 my_rank = xcomm_rank(comm)
 master=0

 write(msg,'(a,i3)') ' xc_kernel: calculating exchange-correlation kernel using ixc = ',ixc
 call wrtout(std_out,msg,'COLL')

 call dtsetCopy(DtGW,Dtset)
 DtGW%intxc = 0
 DtGW%ixc   = ixc

!Redefine xclevel.
 DtGW%xclevel=0
 if ( ( 1<=DtGW%ixc .and. DtGW%ixc<=10).or.(30<=DtGW%ixc .and. DtGW%ixc<=39) ) DtGW%xclevel=1 ! LDA
 if ( (11<=DtGW%ixc .and. DtGW%ixc<=19).or.(23<=DtGW%ixc .and. DtGW%ixc<=29) ) DtGW%xclevel=2 ! GGA
 if ( 20<=DtGW%ixc .and. DtGW%ixc<=22 ) DtGW%xclevel=3 ! ixc for TDDFT kernel tests

 if (ALL(DtGW%xclevel/=(/1,2/))) then
   write(msg,'(a,i0)')"Unsupported xclevel = ",DtGW%xclevel
   MSG_ERROR(msg)
 end if
 
 ngfft1=ngfft(1)
 ngfft2=ngfft(2)
 ngfft3=ngfft(3)

 if (ixc>=1.and.ixc<11) then ! LDA case
!  nkxc=3
!  nkxc=1
   nkxc=2*min(nspden,2)-1
 else ! GGA
   nkxc=23
   ABI_CHECK(dtset%xclevel==2,"Functional should be GGA")
   MSG_ERROR("GGA functional not tested")
 end if
 
 ABI_ALLOCATE(kxcr,(nfft_tot,nkxc))

!gsqcut and rhog are zeroed because they are not used by rhohxc if 1<=ixc<=16 and option=0
 gsqcut=zero

 ABI_ALLOCATE(rhog,(2,nfft_tot))
 ABI_ALLOCATE(vhartr,(nfft_tot))
 rhog(:,:)=zero
!MG FIXME this is the 3D core electron density for XC core correction (bohr^-3)
!should implement the non linear core correction 
 n3xccc=0       
 ABI_ALLOCATE(xccc3d,(n3xccc))
 ABI_ALLOCATE(vxclda,(nfft_tot,nspden))

 option=2 ! 2 for Hxc and kxcr (no paramagnetic part if nspden=1)
 qphon =zero

!to be adjusted for the call to rhohxc
 nk3xc=1
 izero=0

!Compute the kernel.
 call rhohxc(DtGW,enxc,gsqcut,izero,kxcr,MPI_enreg_seq,nfft_tot,ngfft,&
& dum,0,dum,0,nkxc,nk3xc,nspden,n3xccc,option,rhog,rhor,Cryst%rprimd,&
& strsxc,1,vhartr,vxclda,vxcavg,xccc3d)

!DEBUG print Kxc
 if (present(dbg_mode)) then
   if (dbg_mode .and. my_rank==master) then
     fname = 'xc_Kxc.xsf'
     unt_dmp = get_unit()
     open(unt_dmp,file=fname,status='unknown',form='formatted')
     call printxsf(ngfft1,ngfft2,ngfft3,kxcr(:,1),Cryst%rprimd,(/zero,zero,zero/),&
&     Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,unt_dmp,0)
     close(unt_dmp)
   end if
 end if
!DEBUG

 ABI_DEALLOCATE(xccc3d)
 ABI_DEALLOCATE(vxclda)

 ABI_ALLOCATE(my_kxcg,(2,nfft_tot))

 do ikxc=1,nkxc
   call fourdp(1,my_kxcg,kxcr(:,ikxc),-1,MPI_enreg_seq,nfft_tot,ngfft,paral_kgb0,0)
   kxcg(:,ikxc)=CMPLX(my_kxcg(1,:),my_kxcg(2,:))
 end do

!write(std_out,*)"kxcr(r=0)",kxcr(1,1)
!write(std_out,*)"my_kxg(G=0)",my_kxcg(:,1)
!write(std_out,*)"SUM kxcr/nfft_tot ",SUM(kxcr(:,1))/nfft_tot
!write(std_out,*)"SUM my_kxg ",SUM(kxcg(:,1))

 ABI_DEALLOCATE(my_kxcg)

!MG this part is never executed, but one should use mkvxc3 for the GGA kernel.
 if (DtGW%xclevel==2) then
   MSG_ERROR("check GGA implementation")
   cplex=2
   ABI_ALLOCATE(phas,(cplex*nfft_tot,npw,nspden))
   ABI_ALLOCATE(kxcpw_r,(cplex*nfft_tot,nspden))
   ABI_ALLOCATE(xx,(3,nfft_tot))
   ABI_ALLOCATE(kxcpw_g,(2,nfft_tot))
   
   kxcg = czero

!  find the coordinates for all r in the FFT grid
   ir=0
   do i3=1,ngfft3
     do i2=1,ngfft2
       do i1=1,ngfft1
         ir=ir+1
         xx(1,ir)=dble((i1-1))/ngfft1
         xx(2,ir)=dble((i2-1))/ngfft2
         xx(3,ir)=dble((i3-1))/ngfft3
       end do
     end do
   end do

   do iq=1,1

!    Calculate at once exp(i(G+q).r), for all possible q,G,r
     do ig=1,npw
       gpqx=dble(gvec(1,ig))
       gpqy=dble(gvec(2,ig))
       gpqz=dble(gvec(3,ig))
       do ir=1,nfft_tot
         expo=gpqx*xx(1,ir)+gpqy*xx(2,ir)+gpqz*xx(3,ir)              
         phas(2*ir-1,ig,1)= cos(two_pi*expo)
         phas(2*ir,ig,1) =  sin(two_pi*expo)
       end do
     end do

!    Calculate $K(G,G'',q)=\frac{1}{nfft_tot}\sum_{r} exp(-i(q+G_{2}).r_{2} kxcr(r_{1}r_{2}) exp(i(q+G_{1}).r_{1} $

     do igp=1,npw

       kxcpw_r(:,:)=zero

       call mkvxc3(cplex,kxcr,MPI_enreg_seq,nfft_tot,ngfft,nkxc,nspden,n3xccc,option,&
&       paral_kgb0,qphon(:),phas(:,igp,:),Cryst%rprimd,kxcpw_r,xccc3d)

!      FFT the first index to --> to G space
       call fourdp(cplex,kxcpw_g(:,:),kxcpw_r(:,1),-1,MPI_enreg_seq,nfft_tot,ngfft,paral_kgb0,0)

!      kxcg(:,igp,iq)=CMPLX(kxcpw_g(1,igfft(:)),kxcpw_g(2,igfft(:)))
!      kxcg(:,igp)=CMPLX(kxcpw_g(1,igfft(:)),kxcpw_g(2,igfft(:)))

     end do ! igp
   end do ! iq

   ABI_DEALLOCATE(phas)
   ABI_DEALLOCATE(kxcpw_r)
   ABI_DEALLOCATE(xx)
   ABI_DEALLOCATE(kxcpw_g)
 end if !xclevel==2

 call dtsetFree(DtGW)
 ABI_DEALLOCATE(kxcr)

end subroutine xc_kernel
!!***
