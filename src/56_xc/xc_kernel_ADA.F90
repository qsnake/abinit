!{\src2tex{textfont=tt}}
!!****f* ABINIT/xc_kernel_ADA
!! NAME
!! xc_kernel_ADA
!!
!! FUNCTION
!! Calculate exchange-correlation kernel in reciprocal space
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MS)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Dtset <type(dataset_type)>=all input variables in this dataset
!! Cryst<Crystal_structure>=Info on the unit cell.
!! ixc = choice for the exchange-correlation potential.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nfft = total number of points on the FFT grid.
!! rhor(nfft,nspden) = the charge density on the FFT grid.
!!  (total in first half and spin-up in second half if nsppol=2)
!! npw: the size of kernel matrix
!! dim_kxcg=dimension of the kernel.
!! comm=MPI communicator.
!! [dbg_mode]=Set it to .TRUE. to switch on the debug mode.
!!
!! OUTPUT
!!  kxcg(nfft,dim_kxcg) = the exchange-correlation potential on the FFT grid.
!!  warning: the kernel is not divided by unit cell volume
!!
!! NOTES
!!  No xc quadrature
!!  No nl core correction
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      dtsetcopy,dtsetfree,fourdp,fourdp_6d,initmpi_seq,printxsf,rhohxc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine xc_kernel_ADA(Dtset,Cryst,ixc,ngfft,nfft,nspden,rhor,&
&                        npw,nqibz,qibz,fxc_ADA,gvec,comm,kappa_init,dbg_mode)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_xmpi
 use m_errors
 use m_crystal

 use m_io_tools,      only : get_unit
 use m_numeric_tools, only : hermitianize
 use m_fft_mesh,      only : g2ifft

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xc_kernel_ADA'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
 use interfaces_53_ffts
 use interfaces_56_xc, except_this_one => xc_kernel_ADA
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,nfft,nspden,npw,comm
 real(dp),intent(in),optional :: kappa_init
 logical,optional,intent(in) :: dbg_mode
 type(dataset_type),intent(in) :: Dtset
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npw),ngfft(18)
 integer,intent(in) :: nqibz
 real(dp),intent(in) :: rhor(nfft,nspden)
 real(dp),intent(in) :: qibz(3,nqibz)
 complex(gwpc),intent(out) :: fxc_ADA(npw,npw,nqibz)

!Local variables ------------------------------
!scalars
 integer,parameter :: paral_kgb0=0
 integer :: i1,i2,i3,ig,igp,ir,irp,n3xccc,ngfft1,ngfft2,izero !,isp
 integer :: ngfft3,nkxc,option,ikxc,ierr,nproc
 integer :: nk3xc,igrid,iqbz,my_rank,master,unt_dmp,gmgp_idx
 real(dp) :: enxc,gsqcut,ucvol !,rs,Kx,Kc
 real(dp) :: vxcavg,kappa,abs_qpg_sq,abs_qpgp_sq
 real(dp) :: difx,dify,difz,inv_kappa_sq
 character(len=500) :: msg,fname
 type(dataset_type) :: DtGW
 type(MPI_type) :: MPI_enreg_seq 
!arrays
 real(dp) :: qpg(3),qpgp(3),qphon(3),strsxc(6),q_point(3)
 real(dp),allocatable :: dum(:),kxcr(:,:)
 real(dp),allocatable :: rhog(:,:),vhartr(:),vxclda(:,:)
 real(dp),allocatable :: xccc3d(:),my_rhor(:,:)
 real(dp),allocatable :: my_kxcg(:,:)
 real(dp),allocatable :: rhotilder(:,:)
 complex(gwpc),allocatable :: my_fxc_ADA_ggpq(:,:,:)
 complex(gwpc),allocatable :: FT_fxc_ADA_ggpq(:,:,:),dummy(:,:)
 real(dp),allocatable :: rvec(:,:),my_fxc_ADA_rrp(:,:)
 real(dp) :: rmrp(3),abs_rmrp
 integer :: n1,n2,n3,ig_idx_fft(npw)

! ************************************************************************

 ABI_CHECK(Dtset%nsppol==1,'nsppol/=1 not coded')
 ABI_CHECK(nspden==1,'nsppol/=1 not coded')
 ABI_CHECK(nfft==PRODUCT(ngfft(1:3)),"mismatch in nfftot")

!Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

 my_rank = xcomm_rank(comm)
 nproc   = xcomm_size(comm)
 master=0

 write(msg,'(a,i3)') ' xc_kernel_ADA: calculating exchange-correlation kernel using ixc = ',ixc
 call wrtout(std_out,msg,'COLL')
 call wrtout(std_out,' xc_kernel_ADA: using smeared density','COLL')

 if (.not.present(kappa_init)) then
   kappa = 2.1_dp
 else
   kappa = kappa_init
 end if
 write(msg,'(a,F10.3)') ' xc_kernel_ADA: inverse smearing length, kappa = ',kappa
 call wrtout(std_out,msg,'COLL')
 inv_kappa_sq = one/(kappa*kappa)

 call dtsetCopy(DtGW,Dtset)
 DtGW%intxc = 0
 DtGW%ixc   = ixc

!Redefine xclevel.
 DtGW%xclevel=0
 if( ( 1<=DtGW%ixc .and. DtGW%ixc<=10).or.(30<=DtGW%ixc .and. DtGW%ixc<=39) )DtGW%xclevel=1 ! LDA
 if( (11<=DtGW%ixc .and. DtGW%ixc<=19).or.(23<=DtGW%ixc .and. DtGW%ixc<=29) )DtGW%xclevel=2 ! GGA
 if( 20<=DtGW%ixc .and. DtGW%ixc<=22 )DtGW%xclevel=3 ! ixc for TDDFT kernel tests

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
   nkxc=2*min(DtGW%nspden,2)-1
 else ! GGA
   nkxc=23
   ABI_CHECK(dtset%xclevel==2,"Functional should be GGA")
   MSG_ERROR("GGA functional not implemented for ADA vertex")
 end if
 
 ABI_ALLOCATE(kxcr,(nfft,nkxc))

!gsqcut and rhog are zeroed because they are not used by rhohxc if 1<=ixc<=16 and option=0
 gsqcut=zero

 ABI_ALLOCATE(rhog,(2,nfft))
 ABI_ALLOCATE(vhartr,(nfft))
 rhog(:,:)=zero
!MG FIXME this is the 3D core electron density for XC core correction (bohr^-3)
!should implement the non linear core correction 
 n3xccc=0       
 ABI_ALLOCATE(xccc3d,(n3xccc))
 ABI_ALLOCATE(vxclda,(nfft,nspden))

 option=2 ! 2 for Hxc and kxcr (no paramagnetic part if nspden=1)
 qphon(:)=0.0

!to be adjusted for the call to rhohxc
 nk3xc=1

!Compute the kernel.
 izero=0

!DEBUG print density
 if (present(dbg_mode)) then
   if (dbg_mode.and.my_rank==master) then
     unt_dmp = get_unit()
     fname =  'xc_ADA_den.xsf'
     open(unt_dmp,file=TRIM(fname),status='unknown',form='formatted')
     call printxsf(ngfft1,ngfft2,ngfft3,rhor(:,1),Cryst%rprimd,(/zero,zero,zero/),&
&     Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,unt_dmp,0)
     close(unt_dmp)
   end if
 end if
!DEBUG

!Calculate the smeared density
 ABI_ALLOCATE(my_rhor,(nfft,nspden))
 ABI_ALLOCATE(rhotilder,(nfft,nspden))
 ucvol = Cryst%ucvol
 my_rhor = rhor
!do isp = 1,nsppol
!call calc_smeared_density(my_rhor(:,isp),1,rhotilder(:,isp),nfft,ngfft,npw,&
!&   gvec,Cryst%gprimd,Cryst%ucvol,MPI_enreg_seq,paral_kgb0,kappa_in=kappa) 
!my_rhor(:,isp) = rhotilder(:,isp)
!end do

!DEBUG print smeared density
 if (present(dbg_mode)) then
   if (dbg_mode.and.my_rank==master) then
     unt_dmp = get_unit()
     fname = 'xc_ADA_smeared_den.xsf'
     open(unt_dmp,file=TRIM(fname),status='unknown',form='formatted')
     call printxsf(ngfft1,ngfft2,ngfft3,my_rhor(:,1),Cryst%rprimd,(/zero,zero,zero/),&
&     Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,unt_dmp,0)
     close(unt_dmp)
   end if
 end if
!DEBUG

 call rhohxc(DtGW,enxc,gsqcut,izero,kxcr,MPI_enreg_seq,nfft,ngfft,&
& dum,0,dum,0,nkxc,nk3xc,nspden,n3xccc,option,rhog,my_rhor,Cryst%rprimd,&
& strsxc,1,vhartr,vxclda,vxcavg,xccc3d)

!Check for extreme (NaN) values
!do ir=1,nfft
!if (isnan(kxcr(ir,1))) kxcr(ir,1) = HUGE(kxcr(ir,1))
!end do

!DEBUG test with direct way of calculating Kxc
!do i1=1,nfft  
!rs = (three/(four_pi*my_rhor(i1,1)))**third
!Kx = 16._dp/27._dp*0.3141592653589793e1_dp*(rs**2)*(-0.4581652_dp)
!
!Kc =  -0.4e1_dp / 0.9e1_dp * 0.3141592654e1_dp * rs ** 4 &
!* (0.207271333333333333333333333333e-1_dp * &
!(-0.177442658629204480000000e3_dp * rs - 0.17565190511219200000000e2_dp &
!* sqrt(rs) - 0.1332650665120000e2_dp * rs ** 2 &
!- 0.51031691247948928000000e2_dp * rs ** (0.3e1_dp / 0.2e1_dp)) &
!* rs ** (-0.3e1_dp / 0.2e1_dp) / (rs + 0.37274400e1_dp * sqrt(rs) &
!+ 0.129352000e2_dp) ** 2 / (-sqrt(rs) - 0.1049800_dp) &
!+ 0.518178333333333333333333333333e-2_dp * rs ** (-0.3e1_dp / 0.2e1_dp) &
!* (0.617071835390850041282140897280e3_dp * sqrt(rs) &
!+ 0.659369347307557491857191871552e5_dp * rs ** 2 + &
!0.700403648491298930017835369562e5_dp * rs ** (0.3e1_dp / 0.2e1_dp) &
!+ 0.398437532951539263722720308167e5_dp * rs ** (0.5e1_dp / 0.2e1_dp) &
!+ 0.368852071032531998953472000000e4_dp * rs ** (0.7e1_dp / 0.2e1_dp) &
!+ 0.5330602660480000e2_dp * rs ** (0.9e1_dp / 0.2e1_dp) &
!+ 0.143783940386264738593799346176e5_dp * rs ** 3 &
!+ 0.124672564145568409213848436081e5_dp * rs &
!+ 0.557398029956167136000000e3_dp * rs ** 4) &
!/ (rs + 0.37274400e1_dp * sqrt(rs) + 0.129352000e2_dp) ** 4 &
!/ (-sqrt(rs) - 0.1049800_dp) ** 2)
!kxcr(i1,1) = Kx + Kc
!end do
!END DEBUG

!DEBUG print Kxc
 if (present(dbg_mode)) then
   if (dbg_mode.and.my_rank==master) then
     unt_dmp = get_unit()
     write(fname,'(a)') 'xc_ADA_Kxc.xsf'
     open(unt_dmp,file=TRIM(fname),status='unknown',form='formatted')
     call printxsf(ngfft1,ngfft2,ngfft3,kxcr(:,1),Cryst%rprimd,(/zero,zero,zero/),&
&     Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,unt_dmp,0)
     close(unt_dmp)
   end if
 end if
!DEBUG

 ABI_DEALLOCATE(xccc3d)
 ABI_DEALLOCATE(vxclda)
 ABI_DEALLOCATE(vhartr)

 ABI_ALLOCATE(my_kxcg,(2,nfft))

 do ikxc=1,nkxc
   call fourdp(1,my_kxcg,kxcr(:,ikxc),-1,MPI_enreg_seq,nfft,ngfft,paral_kgb0,0)
!  kxcg(:,ikxc)=CMPLX(my_kxcg(1,:),my_kxcg(2,:))
 end do
!TODO Check symmetry of kxcg

!set up ADA vertex
 ABI_ALLOCATE(my_fxc_ADA_ggpq,(npw,npw,nqibz))
 my_fxc_ADA_ggpq = czero
!Calculate f_xc(R,R')=(kappa^2/2)K_xc[\tilde{n}](G-G')
!x(1/(kappa^2+|q+G|^2) + 1/(kappa^2+|q+G'|^2
!First get G vectors and indices

 ierr=0
 do iqbz=1,nqibz
   q_point(:) = qibz(:,iqbz)
!  
   do ig=1,npw
     do igp=1,npw
!      Calculate |q+G| and |q+G'|
       qpg(:) = two_pi*MATMUL(Cryst%gprimd,q_point(:)+gvec(:,ig))
       qpgp(:) = two_pi*MATMUL(Cryst%gprimd,q_point(:)+gvec(:,igp))
       abs_qpg_sq = 1.0_dp/(1.0_dp+dot_product(qpg,qpg)*inv_kappa_sq)
       abs_qpgp_sq = 1.0_dp/(1.0_dp+dot_product(qpgp,qpgp)*inv_kappa_sq)
       
       gmgp_idx = g2ifft(gvec(:,ig)-gvec(:,igp),ngfft)
       if (gmgp_idx>0) then
         my_fxc_ADA_ggpq(ig,igp,iqbz) = half*CMPLX(my_kxcg(1,gmgp_idx), my_kxcg(2,gmgp_idx))*(abs_qpg_sq+abs_qpgp_sq)
       else 
         ierr=ierr+1
         my_fxc_ADA_ggpq(ig,igp,iqbz) = czero
       end if
     end do
   end do
   if (ierr/=0) then 
     write(msg,'(a,i4,3a)')&
&     ' Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&     ' Enlarge the FFT mesh to get rid of this problem. '
     MSG_WARNING(msg)
   end if
 end do

 fxc_ADA = my_fxc_ADA_ggpq 

!do iqbz=1,nqibz
!call hermitianize(my_fxc_ADA_ggpq(:,:,iqbz),"All")
!end do


!DEBUG check symmetry
 if (.FALSE.) then
!  do iqbz=1,nkptgw
!  do ig=1,npw
!  do igp=ig,npw
!  if (ABS(REAL(fxc_ADA(ig,igp,iqbz))-REAL(fxc_ADA(igp,ig,iqbz)))>tol15.OR.&
!  ABS(AIMAG(fxc_ADA(ig,igp,iqbz))-AIMAG(-fxc_ADA(igp,ig,iqbz)))>tol15) then
!  write(std_out,*) 'Elements:'
!  write(std_out,*) 'fxc_ADA(ig,igp,iqbz):',ig,igp,iqbz,fxc_ADA(ig,igp,iqbz)
!  write(std_out,*) 'fxc_ADA(igp,ig,iqbz):',igp,ig,iqbz,fxc_ADA(igp,ig,iqbz)
!  MSG_ERROR('fxc_ADA not symmetric')
!  end if
!  end do
!  end do
!  end do
   
!  write(std_out,*)"kxcr(r=0)",kxcr(1,1)
!  write(std_out,*)"my_kxg(G=0)",my_kxcg(:,1)
!  write(std_out,*)"SUM kxcr/nfft ",SUM(kxcr(:,1))/nfft
!  write(std_out,*)"SUM my_kxg ",SUM(kxcg(:,1))

!  DEBUG Check FT to real space
!  The real-space expression is:
!  f_xc(R,R')=(1/2)(kappa^2/(4*Pi))
!  \{K_xc[\tilde{n(R)}]+K_xc[\tilde{n(R')}]\}
!  x exp(-kappa|R-R'|)/|R-R'|
   ABI_ALLOCATE(my_fxc_ADA_rrp,(nfft,nfft))
   ABI_ALLOCATE(FT_fxc_ADA_ggpq,(npw,npw,nqibz))
   ABI_ALLOCATE(rvec,(3,nfft))
   ABI_ALLOCATE(dummy,(nfft,nfft))
   my_fxc_ADA_rrp=zero; FT_fxc_ADA_ggpq=czero; dummy=czero; rvec=zero
   
!  First find coordinates of real-space fft points
   igrid = 0
   ngfft1 = ngfft(1)
   ngfft2 = ngfft(2)
   ngfft3 = ngfft(3)
   do i3=0,ngfft3-1
     difz=dble(i3)/dble(ngfft3)
     do i2=0,ngfft2-1
       dify=dble(i2)/dble(ngfft2)
       do i1=0,ngfft1-1
         difx=dble(i1)/dble(ngfft1)
         igrid = igrid + 1
         rvec(1,igrid)=difx*Cryst%rprimd(1,1)+dify*Cryst%rprimd(1,2)+difz*Cryst%rprimd(1,3)
         rvec(2,igrid)=difx*Cryst%rprimd(2,1)+dify*Cryst%rprimd(2,2)+difz*Cryst%rprimd(2,3)
         rvec(3,igrid)=difx*Cryst%rprimd(3,1)+dify*Cryst%rprimd(3,2)+difz*Cryst%rprimd(3,3)
       end do
     end do
   end do
   if (igrid/=nfft) STOP 'ERROR in xc_kernel_ADA: igrid not equal to nfft'
   
!  Construct kernel in real space
   do ir=1,nfft
     do irp=ir,nfft 
       rmrp(:) = rvec(:,ir)-rvec(:,irp)
       abs_rmrp = sqrt(dot_product(rmrp,rmrp))
       my_fxc_ADA_rrp(ir,irp) = eighth*kappa*kappa*piinv* &
       (kxcr(ir,1)+kxcr(irp,1))* &
       EXP(-kappa*abs_rmrp)/(abs_rmrp+1.e-3_dp)
!      (a small convergence factor is introduced
!      to avoid a singularity)
       my_fxc_ADA_rrp(irp,ir) = my_fxc_ADA_rrp(ir,irp)
     end do
   end do

!  Find FFT index for all G   
   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
!  Use the following indexing (N means ngfft of the adequate direction)
!  0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= kg
!  1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index
   do ig=1,npw
     i1=modulo(gvec(1,ig),n1)
     i2=modulo(gvec(2,ig),n2)
     i3=modulo(gvec(3,ig),n3)
     ig_idx_fft(ig)=i1+1+n1*(i2+n2*i3)
   end do
!  FT kernel to reciprocal space for each q
   do iqbz=1,nqibz
     dummy = CMPLX(my_fxc_ADA_rrp,0.0_dp)
!    Multiply with q-point phase factors exp(-iq.r)*f_xc(r,r')*exp(iq.r')
     do ir=1,nfft
       do irp=1,nfft
!        Calculate q (variables defined for other purposes
!        are being re-used as dummy variables)
         q_point(:) = qibz(:,iqbz)
         qpg(:) = two_pi*MATMUL(Cryst%gprimd,q_point(:))
         abs_qpg_sq = dot_product(qpg(:),rvec(:,ir))
         abs_qpgp_sq = dot_product(qpg(:),rvec(:,irp))
         dummy(ir,irp) = EXP(-j_dpc*abs_qpg_sq)* &
&         dummy(ir,irp)* &
&         EXP(j_dpc*abs_qpgp_sq)
       end do
     end do
     call fourdp_6d(2,dummy,-1,MPI_enreg_seq,nfft,ngfft,paral_kgb0,0)
     do ig=1,npw
       do igp=1,npw
         FT_fxc_ADA_ggpq(ig,igp,iqbz) = dummy(ig_idx_fft(ig),ig_idx_fft(igp))
       end do
     end do

!    Output
     msg=''
     if (iqbz<10) write(msg,'(a,i1,a)') './debug_fxc_ADA_q',iqbz,'.dat'
     if ((iqbz>9).and.(iqbz<100)) write(msg,'(a,i2,a)') './debug_fxc_ADA_q',iqbz,'.dat'
     if ((iqbz>99).and.(iqbz<1000)) write(msg,'(a,i3,a)') './debug_fxc_ADA_q',iqbz,'.dat'
     open(777,file=TRIM(msg),STATUS='REPLACE')
     do igp=1,npw
       do ig=1,npw
         write(777,*) ig,igp,REAL(my_fxc_ADA_ggpq(ig,igp,iqbz)),AIMAG(my_fxc_ADA_ggpq(ig,igp,iqbz)), &
&         REAL(FT_fxc_ADA_ggpq(ig,igp,iqbz)),AIMAG(FT_fxc_ADA_ggpq(ig,igp,iqbz)), &
&         ABS(ABS(my_fxc_ADA_ggpq(ig,igp,iqbz))-ABS(FT_fxc_ADA_ggpq(ig,igp,iqbz)))
       end do
       write(777,*) ''
     end do
     close(777)


   end do ! iqbz

   MSG_ERROR('Stopping in xc_kernel_ADA for debugging')

   ABI_DEALLOCATE(rvec)
   ABI_DEALLOCATE(my_fxc_ADA_rrp)
   ABI_DEALLOCATE(FT_fxc_ADA_ggpq)
   
   if (DtGW%xclevel==2) then
     MSG_ERROR(" GGA not implemented for xc_kernel_ADA")
   end if !xclevel==2

 end if ! Debugging section


 call dtsetFree(DtGW)
 ABI_DEALLOCATE(my_kxcg)
 ABI_DEALLOCATE(my_rhor)
 ABI_DEALLOCATE(rhotilder)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(kxcr)

end subroutine xc_kernel_ADA
!!***
