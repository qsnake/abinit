!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhohxc
!! NAME
!! rhohxc
!!
!! FUNCTION
!! Start from the density or spin-density, and
!! compute Hartree (if option>=1) and xc correlation potential and energies.
!! Eventually compute xc kernel (if option=-2, 2, 3, 10 or 12).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MF, GZ, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | intxc=0 for old quadrature; 1 for new improved quadrature
!!   | ixc= choice of exchange-correlation scheme (see above, and below)
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!! (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  izero=if 1, unbalanced components of Vhartree(g) have to be set to zero
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfftf,nspden*nhatdim)= -PAW only- compensation density
!!  nhatdim= -PAW only- 0 if nhat array is not used ; 1 otherwise
!!  nhatgr(nfftf,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the kxc array. If /=0,
!!   the exchange-correlation kernel must be computed.
!!  nspden=number of spin-density components
!!  n3xccc=dimension of the xccc3d array (0 or nfft or cplx*nfft).
!!  option=0 for xc only (exc, vxc, strsxc),
!!         1 for Hxc (idem + vhartr) ,
!!         2 for Hxc and kxc (no paramagnetic part if nspden=1)
!!        10 for xc  and kxc with only LDA part (d2Exc/drho^2)
!!        12 for Hxc and kxc with only LDA part (d2Exc/drho^2)
!!         3 for Hxc, kxc and k3xc
!!        -2 for Hxc and kxc (with paramagnetic part if nspden=1)
!!  rhog(2,nfft)=electron density in G space
!!  rhor(nfft,nspden)=electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if nspden=2)
!!   (total in first comp. and magnetization in comp. 2 to 4 if nspden=4)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!
!! OUTPUT
!!  enxc=returned exchange and correlation energy (hartree).
!!
!!  === Only if abs(option)=3 ===
!!
!!  k3xc(nfft,nk3xc)=third derivative of the XC energy functional of the density,
!!    at each point of the real space grid (only in the LDA or LSDA)
!!    if nspden==1: return k3xc(:,1)= d3Exc/drho3
!!    if nspden>=2, return  k3xc(:,1)=d3Exc/drho_up drho_up drho_up
!!                          k3xc(:,2)=d3Exc/drho_up drho_up drho_dn
!!                          k3xc(:,3)=d3Exc/drho_up drho_dn drho_dn
!!                          k3xc(:,4)=d3Exc/drho_dn drho_dn drho_dn
!!
!!  === Only if abs(option)=2, -2, 3, 10, 12 ===
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!                 (returned only if nkxc/=0)
!!   allowed if LDAs (dtset%xclevel=1 or option=(10 or 12)) :
!!    if nspden==1: return kxc(:,1)= d2Exc/drho2
!!       that is 1/2 ( d2Exc/drho_up drho_up + d2Exc/drho_up drho_dn )
!!    if nspden==1: also return kxc(:,2)= d2Exc/drho_up drho_dn
!!    if nspden>=2, return  kxc(:,1)=d2Exc/drho_up drho_up
!!                          kxc(:,2)=d2Exc/drho_up drho_dn
!!                          kxc(:,3)=d2Exc/drho_dn drho_dn
!!   allowed also if GGAs (dtset%xclevel=2 and option=2 or 3)
!!    for the time being, treat all cases as spin-polarized, with nkxc=23
!!    kxc(:,1)= d2Ex/drho_up drho_up
!!    kxc(:,2)= d2Ex/drho_dn drho_dn
!!    kxc(:,3)= dEx/d(abs(grad(rho_up))) / abs(grad(rho_up))
!!    kxc(:,4)= dEx/d(abs(grad(rho_dn))) / abs(grad(rho_dn))
!!    kxc(:,5)= d2Ex/d(abs(grad(rho_up))) drho_up / abs(grad(rho_up))
!!    kxc(:,6)= d2Ex/d(abs(grad(rho_dn))) drho_dn / abs(grad(rho_dn))
!!    kxc(:,7)= 1/abs(grad(rho_up)) * d/d(abs(grad(rho_up)) (dEx/d(abs(grad(rho_up))) /abs(grad(rho_up)))
!!    kxc(:,8)= 1/abs(grad(rho_dn)) * d/d(abs(grad(rho_dn)) (dEx/d(abs(grad(rho_dn))) /abs(grad(rho_dn)))
!!    kxc(:,9)= d2Ec/drho_up drho_up
!!    kxc(:,10)=d2Ec/drho_up drho_dn
!!    kxc(:,11)=d2Ec/drho_dn drho_dn
!!    kxc(:,12)=dEc/d(abs(grad(rho))) / abs(grad(rho))
!!    kxc(:,13)=d2Ec/d(abs(grad(rho))) drho_up / abs(grad(rho))
!!    kxc(:,14)=d2Ec/d(abs(grad(rho))) drho_dn / abs(grad(rho))
!!    kxc(:,15)=1/abs(grad(rho)) * d/d(abs(grad(rho)) (dEc/d(abs(grad(rho))) /abs(grad(rho)))
!!    kxc(:,16)=rho_up
!!    kxc(:,17)=rho_dn
!!    kxc(:,18)=gradx(rho_up)
!!    kxc(:,19)=gradx(rho_dn)
!!    kxc(:,20)=grady(rho_up)
!!    kxc(:,21)=grady(rho_dn)
!!    kxc(:,22)=gradz(rho_up)
!!    kxc(:,23)=gradz(rho_dn)
!!
!!  strsxc(6)= contribution of xc to stress tensor (hartree/bohr^3),
!!   given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!   (note: fxc is rho*exc in the following)
!!   Explicitely : strsxc(mu,nu) = (1/N) Sum(i=1,N)
!!    ( delta(mu,nu) * [  exc(i)rhotot(i)
!!               - depsxc_drho(up,i)*rhor(up,i)-depsxc_drho(dn,i)*rhor(dn,i)]
!!     - gradrho(up,mu)*gradrho(up,nu) * depsxc_dgradrho(up,i) / gradrho(up,i)
!!     - gradrho(dn,mu)*gradrho(dn,nu) * depsxc_dgradrho(dn,i) / gradrho(dn,i) )
!!  vhartr(nfft)=Hartree potential (returned if option/=0 and option/=10)
!!  vxc(nfft,nspden)=xc potential
!!    (spin up in first half and spin down in second half if nspden=2)
!!    (v^11, v^22, Re[V^12], Im[V^12] if nspden=4)
!!  vxcavg=<Vxc>=unit cell average of Vxc = (1/ucvol) Int [Vxc(r) d^3 r].
!!  vxctau=(only for meta-GGA): derivative of XC energy density with respect to
!!    kinetic energy density (depsxcdtau). The array vxctau(nfft,nspden,4) contains also
!!    the gradient of vxctau (gvxctau) in vxctau(:,:,2:4)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!
!! NOTES
!! Start from the density, and compute Hartree (if option>=1) and xc correlation potential and energies.
!! Eventually compute xc kernel (if option=-2, 2, 3, 10 or 12).
!! Allows a variety of exchange-correlation functionals
!! according to ixc. Here is a list of allowed values.
!!                                                    subroutine name
!!   <0 means use of libxc
!!    0 means no xc applied (usually for testing)
!! *LDA,LSD
!!    1 means new Teter (4/93) with spin-pol option        xcspol
!!    2 means Perdew-Zunger-Ceperley-Alder                 xcpzca
!!    3 means old Teter (4/91) fit to Ceperley-Alder data  xctetr
!!    4 means Wigner                                       xcwign
!!    5 means Hedin-Lundqvist                              xchelu
!!    6 means "X-alpha" xc                                 xcxalp
!!    7 mean Perdew-Wang 92 LSD fit to Ceperley-Alder data xcpbe
!!    8 mean Perdew-Wang 92 LSD , exchange-only            xcpbe
!!    9 mean Perdew-Wang 92 Ex+Ec_RPA  energy              xcpbe
!!   10 means RPA LSD energy (only the energy !!)          xcpbe
!! *GGA
!!   11 means Perdew-Burke-Ernzerhof GGA functional        xcpbe
!!   12 means x-only Perdew-Burke-Ernzerhof GGA functional xcpbe
!!   13 means LDA (ixc==7), except that the xc potential
!!      is given within the van Leeuwen-Baerends GGA       xclb
!!   14 means revPBE GGA functional                        xcpbe
!!   15 means RPBE GGA functional                          xcpbe
!!   16 means HCTH GGA functional                          xchcth
!!   23 means WC GGA functional                            xcpbe
!!   24 means C09x GGA exchange functional                 xcpbe
!! *Fermi-Amaldi
!!   20 means Fermi-Amaldi correction
!!   21 means Fermi-Amaldi correction with LDA(ixc=1) kernel
!!   22 means Fermi-Amaldi correction with hybrid BPG kernel
!!
!! NOTE: please update echo_xc_name.F90 if you add new functional (apart from libxc)
!!
!! Allow for improved xc quadrature (intxc=1) by using the usual FFT grid
!! as well as another, shifted, grid, and combining both results.
!! Spin-polarization is allowed only with ixc=0, 1, and GGAs until now.
!! Note : this routine has been optimized already. See the end of this routine.
!!
!! To make the variable names easier to understand, a rule notation is tentatively proposed here:
!!   rho ---> means density
!!   tau ---> means kinetic energy density
!!   exc ---> means exchange-correlation energy density per particle
!!   epsxc ---> means rho*exc == exchange-correlation energy density
!!   vxc ---> means exchange-correlation potential
!!   bigexc ---> means exchange-correlation energy E_xc (for the moment it is named "enxc")
!!   m_norm ---> means norm of magnetization
!!
!!   g... --> means gradient of something (e.g. : grho --> means gradient of electron density)
!!   g...2 -> means square norm of gradient of something (e.g. : grho2 -> means square norm of gradient of electron density)
!!   l... --> means laplacian of something (e.g. : lrho --> means laplacian of electron density)
!!   d...d... --> means derivative of something with regards to something else.
!!   (d2...d...d...  ---> means second derivative of ... with regards to ... and to ...) etc...
!!   d... --> without the occurence of the second "d" means that this is an array of several derivative of the same quantity (e.g. : depsxc)
!!
!!   ..._b ----> means a block of the quantity "..." (use in mpi loops which treat the data block by block)
!!   ..._updn -> means that spin up and spin down is available in that array as (..,1) and (..,2). (if nspden >=2 of course).
!!   ..._apn --> in case of positrons are concerned.
!!
!!   for more details about notations please see pdf in /doc/theory/MGGA/
!!
!! PARENTS
!!      afterscfloop,calc_vhxc_me,cvxclda,energy,ftfvw1,kxc_alda,nonlinear
!!      nres2vres,odamix,prcref,prcref_PMA,prctfvw2,respfn,rhotov,scfcv,setvtr
!!      xc_kernel,xc_kernel_ADA
!!
!! CHILDREN
!!      dotprod_vn,drivexc,hartre,leave_new,mean_fftr,metric,mkdenpos,size_dvxc
!!      timab,wrtout,xcden,xcmult,xcomm_init,xcpositron,xcpot,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rhohxc(dtset,enxc,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,nspden,n3xccc,option,rhog,rhor,rprimd, &
& strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,&
& k3xc,electronpositron,taug,taur,vxctau) ! optional argument

 use m_profiling

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_xmpi
#if defined HAVE_DFT_LIBXC
 use libxc_functionals
#endif

 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhohxc'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
 use interfaces_56_xc, except_this_one => rhohxc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,nk3xc,n3xccc,nfft,nhatdim,nhatgrdim,nkxc,nspden,option
 integer,intent(in) :: usexcnhat
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: enxc,vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat(nfft,nspden*nhatdim)
 real(dp),intent(in) :: nhatgr(nfft,nspden,3*nhatgrdim),rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),xccc3d(n3xccc)
 real(dp),intent(out) :: kxc(nfft,nkxc),strsxc(6),vhartr(nfft),vxc(nfft,nspden)
 real(dp),intent(out),optional :: k3xc(1:nfft,1:nk3xc)
 real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
 real(dp),intent(in),optional :: taur(nfft,nspden*dtset%usekden)
 real(dp),intent(out),optional :: vxctau(nfft,nspden*dtset%usekden,4)

!Local variables-------------------------------
!scalars
 integer :: cplex,ierr,ifft,ii,index,inkxc,ipositron,ipts,ishift,ispden,iwarn,iwarnp
 integer :: ixc,jj,mpts,ndvxc,nd2vxc,nfftot,ngr,ngr2,ngrad,ngrad_apn,npts,nspden_apn,nspden_eff
 integer :: nspden_updn,nspgrad,nvxcgrho,old_paral_level,order,spaceComm
 real(dp),parameter :: mot=-one/3.0_dp
 real(dp) :: coeff,divshft,doti,dstrsxc,dvdn,dvdz,epsxc,factor,m_norm_min,nelect,s1,s2,s3
 real(dp) :: strdiag,strsxc1_tot,strsxc2_tot,strsxc3_tot,strsxc4_tot
 real(dp) :: strsxc5_tot,strsxc6_tot,ucvol
!real(dp) :: sum
 logical :: test_libxc,test_nhat,uselibxc,mgga,allow3
 character(len=500) :: message
!arrays
 real(dp) :: gm_norm(3),grho(3),gmet(3,3),gprimd(3,3),qphon(3),rmet(3,3)
 real(dp) :: tsec(2),vxcmean(4)
 real(dp),allocatable :: d2vxc_b(:,:),depsxc(:,:),depsxc_apn(:,:),dvxc_apn(:),dvxc_b(:,:)
 real(dp),allocatable :: exc_b(:),fxc_apn(:),grho2_apn(:),grho2_b_updn(:,:),lrhonow(:,:),lrho_b_updn(:,:)
 real(dp),allocatable :: m_norm(:),nhat_up(:),rho_b_updn(:,:),rho_b(:),rhocorval(:,:),rhonow_apn(:,:,:)
 real(dp),allocatable :: tau_b_updn(:,:),vxc_apn(:,:),vxcgr_apn(:),vxcgrho_b(:,:),vxcrho_b_updn(:,:)
 real(dp),allocatable :: vxc_b_apn(:),vxc_ep(:),vxctau_b_updn(:,:),vxclrho_b_updn(:,:)
 real(dp),allocatable,target :: rhonow(:,:,:)
 real(dp),pointer :: rhonow_ptr(:,:,:)
!integer :: i1,i2,i3,ir,n1,n2,n3
#if defined HAVE_OPENMP
 integer,external :: OMP_GET_NUM_THREADS
#endif

! *************************************************************************

!DEBUG
!write(std_out,*) '+++++++++++++++++'
!write(std_out,*) '-> rhohxc : enter'
!write(std_out,*) 'nfft',nfft
!write(std_out,*) 'nhatdim',nhatdim
!write(std_out,*) 'nhatgrdim',nhatgrdim
!write(std_out,*) 'nkxc',nkxc
!write(std_out,*) 'nspden',nspden
!write(std_out,*) 'nhatdim',nhatdim
!write(std_out,*) 'option',option
!sum=0.0
!do jj=1,nspden
!do kk=1,nfft
!sum=sum+abs(rhor(kk,jj))
!end do
!end do
!write(std_out,*) 'SUM(rhor)',sum
!ENDDEBUG

 if(.false.)write(std_out,*)taug ! just to keep taug as an argument while in development

!Check options
 if(nspden/=1 .and. nspden/=2 .and. nspden/=4)then
   write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&   ' rhohxc :  BUG -',ch10,&
&   '  The only allowed values of nspden are 1, 2, or 4,',ch10,&
&   '  while the argument nspden=',nspden
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(option==3)then
!  LHH,FL,GMR
!  if ((option==3).and.((dtset%ixc/=0) .and. (dtset%ixc/=3).and.(dtset%ixc/=7) .and. (dtset%ixc/=8))) then
   allow3=(dtset%ixc > 0).and.(dtset%ixc /= 3).and.(dtset%ixc /= 7).and.(dtset%ixc /= 8)
#if defined HAVE_DFT_LIBXC
   if(.not.allow3)then
     allow3=(dtset%ixc < 0).and.(libxc_functionals_isgga() .or. libxc_functionals_ismgga())
   end if
#endif
   if(allow3)then
!    LHH,FL,GMR

     write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&     ' rhohxc :  ERROR -',ch10,&
&     '  Third-order xc kernel can only be computed for ixc = 0, 3, 7 or 8,',ch10,&
&     '  while it is found to be',dtset%ixc
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

 if(dtset%icoulomb /= 0)then
   write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&   ' rhohxc : BUG -',ch10,&
&   '  To use non-periodic computation (icoulomb /= 0), ',ch10,&
&   '  use PSolver_rhohxc() instead, while the argument icoulomb=',dtset%icoulomb
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(nspden==4.and.dtset%xclevel==2.and.(abs(option)==2))then
   write(message, '(4a)' ) ch10,&
&   ' rhohxc :  BUG -',ch10,&
&   '  When nspden==4 and GGA, the absolute value of option cannot be 2 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Are we using a functional from libxc? Is the functional a MGGA?
#if defined HAVE_DFT_LIBXC
 uselibxc=(dtset%ixc<0)
 mgga=.false.
 if (uselibxc) then
   mgga=libxc_functionals_ismgga()
   if (mgga .and. dtset%usekden == 0) then
     write(message, '(6a)' )ch10,&
&     ' rhohxc : ERROR -',ch10,&
&     '  The functional is a MGGA, but the kinetic energy density',ch10, &
&     '  is not present. Please set "usekden 1" in the input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
 if(dtset%ixc>=31 .and. dtset%ixc<=34)mgga=.true.
#else
 uselibxc=.false.
 mgga=.false.
 if(dtset%ixc>=31 .and. dtset%ixc<=34)mgga=.true.
#endif

 if(dtset%usekden/=0 .and. dtset%intxc/=0)then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' rhohxc : BUG -',ch10,&
&   '  The use of intxc (/=0) is not implemented in the case of kinetic energy density (usekden/=0), ',ch10,&
&   '  use intxc=0'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!In this routine, hartre, xcden and xcpot are called for real
!densities and potentials, corresponding to zero wavevector
 cplex=1
 qphon(:)=zero
 iwarn=0
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

 if(option/=0.and.option/=10)then
   call hartre(cplex,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qphon,rhog,vhartr)
 end if

!Note : hartre is excluded from the timing
 call timab(81,1,tsec)

 enxc=zero
 epsxc=zero
 vxc(:,:)=zero
 if(present(vxctau))vxctau(:,:,:)=zero
 strsxc(:)=zero
 strsxc1_tot=zero
 strsxc2_tot=zero
 strsxc3_tot=zero
 strsxc4_tot=zero
 strsxc5_tot=zero
 strsxc6_tot=zero

 m_norm_min=EPSILON(0.0_dp)**2

!Some initializations for the electron-positron correlation
 ipositron=0; if (present(electronpositron)) ipositron=electronpositron_calctype(electronpositron)
 if (ipositron==2) then
   electronpositron%e_xc  =zero
   electronpositron%e_xcdc=zero
   nspden_apn=1;ngrad_apn=1;iwarnp=1
   ngrad=1;if(dtset%xclevel==2)ngrad=2
   if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad_apn=2
   if (ngrad_apn>ngrad) then
     write(message, '(4a)' ) ch10,&
&     ' rhohxc :  ERROR -',ch10,&
&     '  GGA for the positron can only be performed with GGA pseudopotentials for the electron !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if (ngrad_apn>1.and.option/=0.and.option/=1.and.option/=10.and.option/=12) then
     write(message, '(4a)' ) ch10,&
&     ' rhohxc :  ERROR -',ch10,&
&     '  You cannot compute full GGA XC kernel for electrons-positron systems !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   ABI_ALLOCATE(depsxc_apn,(nfft,ngrad_apn))
 end if

 if (dtset%ixc==0) then

   vxcavg=zero
   if (nkxc/=0) kxc(:,:)=zero
   if(abs(option)==3)then
     if (nk3xc/=0) k3xc(:,:)=zero
   end if
!  No xc at all is applied (usually for testing)
   write(message, '(a,a,a,a)' ) ch10,&
&   ' rhohxc : WARNING -',ch10,&
&   '  Note that no xc is applied (ixc=0).'
   call wrtout(std_out,message,'COLL')

 else if (dtset%ixc/=20) then

!  ngrad=1 is for LDAs or LSDs, ngrad=2 is for GGAs
   ngrad=1;if(dtset%xclevel==2)ngrad=2
!  ixc 31 to 34 are for mgga test purpose only (fake functionals based on LDA but need the gradients too)
   if(dtset%ixc>=31 .and. dtset%ixc<=34)ngrad=2
!  Test: has a compensation density to be added/substracted (PAW) ?
   test_nhat=((nhatdim==1).and.(usexcnhat==0.or.(ngrad==2.and.nhatgrdim==1)))
!  nspden_updn: 1 for non-polarized, 2 for polarized
   nspden_updn=min(nspden,2)
!  nspden_eff: effective value of nspden used to compute gradients of density:
!  1 for non-polarized system,
!  2 for collinear polarized system or LDA (can be reduced to a collinear system)
!  4 for non-collinear polarized system and GGA
   nspden_eff=nspden_updn;if (nspden==4.and.ngrad==2) nspden_eff=4

!  The different components of depsxc will be
!  for nspden=1,   depsxc(:,1)=d(rho.exc)/d(rho) == (depsxcdrho) == (vxcrho)
!  and if ngrad=2, depsxc(:,2)=1/2*1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|)
!  +1/|grad rho|*d(rho.exc)/d(|grad rho|)
!  == (1/2 * 1/|grho_up| * depsxcd|grho_up|) +  1/|grho| * depsxcd|grho|
!  (vxcgrho=1/|grho| * depsxcd|grho|)
!  (do not forget : |grad rho| /= |grad rho_up| + |grad rho_down|
!  and if mgga,    depsxc(:,3)=d(rho.exc)/d(lapl rho) == (depsxcdlrho) == (vxclrho)
!  
!  for nspden>=2,  depsxc(:,1)=d(rho.exc)/d(rho_up) == (depsxcdrho_up) == (vxcrho_up)
!  depsxc(:,2)=d(rho.exc)/d(rho_down) == (depsxcdrho_dn) == (vxcrho_dn)
!  and if ngrad=2, depsxc(:,3)=1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|) == (1/|grho_up| * depsxcd|grho_up|) == (vxcgrho_up)
!  depsxc(:,4)=1/|grad rho_down|*d(rho.exc)/d(|grad rho_down|) == (1/|grho_dn| * depsxcd|grho_dn|) == (vxcgrho_dn)
!  depsxc(:,5)=1/|grad rho|*d(rho.exc)/d(|grad rho|) == (1/|grho| * depsxcd|grho|) == (vxcgrho)
!  and if mgga,     depsxc(:,6)=d(rho.exc)/d(lapl rho_up) == (depsxcdlrho_up) == (vxclrho_up)
!  depsxc(:,7)=d(rho.exc)/d(lapl rho_dn) == (depsxcdlrho_dn) == (vxclrho_dn)
!  Note: if nspden=4, rho_up=(rho+|m|)/2, rho_down=(rho-|m|)/2

   nspgrad=nspden_updn*ngrad;if(nspden_updn==2.and.ngrad==2)nspgrad=5

   if(mgga)then
     if(nspden==1)then
       nspgrad=nspgrad+1
     else if(nspden==2)then
       nspgrad=nspgrad+2
     else if(nspden==4)then
       write(message, '(4a)' ) ch10,&
&       ' rhohxc :  ERROR -',ch10,&
&       '  meta-GGA XC kernel is not yet implemented for non-collinear magnetism case!'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end if
   ABI_ALLOCATE(depsxc,(nfft,nspgrad))
   depsxc(:,:)=zero

!  Non-collinear magnetism: store norm of magnetization
   if (nspden==4) then
     ABI_ALLOCATE(m_norm,(nfft))
     m_norm=zero
     if ((usexcnhat==1).or.(nhatdim==0)) then
       m_norm(:)=sqrt(rhor(:,2)**2+rhor(:,3)**2+rhor(:,4)**2)
     else
       m_norm(:)=sqrt((rhor(:,2)-nhat(:,2))**2+(rhor(:,3)-nhat(:,3))**2 &
&       +(rhor(:,4)-nhat(:,4))**2)
     end if
   end if

!  rhocorval will contain effective density used to compute gradients:
!  - with core density (if NLCC)
!  - without compensation density (if PAW under certain conditions)
!  - in (up+dn,up) or (n,mx,my,mz) format according to collinearity
!  of polarization and use of gradients (GGA)
   if (n3xccc>0.or.test_nhat.or.nspden_eff/=nspden) then
     ABI_ALLOCATE(rhocorval,(nfft,nspden_eff))
     if (nspden==nspden_eff) then
       rhocorval(:,1:nspden)=rhor(:,1:nspden)
     else if (nspden==4) then
       rhocorval(:,1)=rhor(:,1)
       rhocorval(:,2)=half*(rhor(:,1)+m_norm(:))
     else
       rhocorval=zero
     end if
   end if

!  Add core electron density to effective density
   if (n3xccc>0) then
     rhocorval(:,1)=rhocorval(:,1)+xccc3d(:)
     if(nspden_eff==2) then
       rhocorval(:,2)=rhocorval(:,2)+half*xccc3d(:)
     end if
   end if

!  If PAW, substract compensation density from effective density:
!  - if GGA, because nhat gradients are computed separately
!  - if nhat does not have to be included in XC
   if (test_nhat) then
     if (nspden==nspden_eff) then
       rhocorval(:,1:nspden)=rhocorval(:,1:nspden)-nhat(:,1:nspden)
     else if (nspden==4) then
!      if (usexcnhat==0) then
!      rhocorval(:,1)=rhocorval(:,1)-nhat(:,1)
!      rhocorval(:,2)=rhocorval(:,2)-half*nhat(:,2)
!      else
       ABI_ALLOCATE(nhat_up,(nfft))
       nhat_up(:)=half*(nhat(:,1)+(rhor(:,2)*nhat(:,2) &
&       +rhor(:,3)*nhat(:,3) &
&       +rhor(:,4)*nhat(:,4))/m_norm(ifft))
       rhocorval(:,1)=rhocorval(:,1)-nhat(:,1)
       rhocorval(:,2)=rhocorval(:,2)-nhat_up(:)
!      end if
     end if
   end if

!  rhonow will contain effective density (and gradients if GGA)
!  lrhonow will contain the laplacian if we have a MGGA
   ABI_ALLOCATE(rhonow,(nfft,nspden_eff,ngrad*ngrad))
   if (mgga)  then
     ABI_ALLOCATE(lrhonow,(nfft,nspden_eff))
   end if

!  ====================================================================
!  Loop on unshifted or shifted grids
   do ishift=0,dtset%intxc

!    Set up density on unshifted or shifted grid (will be in rhonow(:,:,1)),
!    as well as the gradient of the density, also on the unshifted
!    or shifted grid (will be in rhonow(:,:,2:4)), if needed.
     if ((n3xccc==0).and.(.not.test_nhat).and.(nspden_eff==nspden)) then
       if (mgga) then
         call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhor,rhonow,lrhonow)
       else
         call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhor,rhonow)
       end if
     else if ((ishift>0).and.(test_nhat)) then
       if (mgga) then
         call xcden(cplex,gprimd,0,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhocorval,rhonow,lrhonow)
       else
         call xcden(cplex,gprimd,0,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhocorval,rhonow)
       end if
     else
       if (mgga) then
         call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhocorval,rhonow,lrhonow)
       else
         call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,dtset%paral_kgb,qphon,rhocorval,rhonow)
       end if
     end if

!    -PAW+GGA: add "exact" gradients of compensation density
     if (test_nhat.and.usexcnhat==1) then
       if (ishift==0) then
         if (nspden==nspden_eff) then
           rhonow(:,1:nspden,1)=rhocorval(:,1:nspden)+nhat(:,1:nspden)
         else if (nspden==4) then
           rhonow(:,1,1)=rhocorval(:,1)+nhat(:,1)
           rhonow(:,2,1)=rhocorval(:,2)+nhat_up(:)
         end if
       else
         if (nspden==nspden_eff) then
           rhocorval(:,1:nspden)=rhocorval(:,1:nspden)+nhat(:,1:nspden)
         else if (nspden==4) then
           rhocorval(:,1)=rhocorval(:,1)+nhat(:,1)
           rhocorval(:,2)=rhocorval(:,2)+nhat_up(:)
         end if
         call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,1,nspden_eff,dtset%paral_kgb,qphon,rhocorval,rhonow)
       end if
       if (ngrad==2.and.nhatgrdim==1.and.nspden==nspden_eff) then
         do ii=1,3
           jj=ii+1
           do ispden=1,nspden
             do ifft=1,nfft
               rhonow(ifft,ispden,jj)=rhonow(ifft,ispden,jj)+nhatgr(ifft,ispden,ii)
             end do
           end do
         end do
       end if
     end if

!    Deallocate temporary arrays
     if (ishift==dtset%intxc) then
       if (n3xccc>0.or.test_nhat.or.nspden_eff/=nspden)  then
         ABI_DEALLOCATE(rhocorval)
       end if
       if (test_nhat.and.nspden/=nspden_eff.and.usexcnhat==1)  then
         ABI_DEALLOCATE(nhat_up)
       end if
     end if

!    In case of non-collinear magnetism, extract up and down density and gradients (if GGA)
     if (nspden==4.and.nspden_eff==nspden) then
       if (ngrad==2) then
         do ifft=1,nfft
           gm_norm(1:3)=zero
           if(m_norm(ifft)>m_norm_min) then
!            if(m_norm(ifft)>rhonow(ifft,1,1)*tol10+tol14) then
             do jj=1,3  ! Compute here nabla(|m|)=(m.nabla(m))/|m| == (g|m| = m/|m| * gm)
               do ii=2,4
                 gm_norm(jj)=gm_norm(jj)+rhonow(ifft,ii,1+jj)*rhonow(ifft,ii,1)
               end do
             end do
             gm_norm(1:3)=gm_norm(1:3)/m_norm(ifft)
           end if
           rhonow(ifft,2,2)=half*(rhonow(ifft,1,2)+gm_norm(1))
           rhonow(ifft,2,3)=half*(rhonow(ifft,1,3)+gm_norm(2))
           rhonow(ifft,2,4)=half*(rhonow(ifft,1,4)+gm_norm(3))
         end do
       end if
       rhonow(:,2,1)=half*(rhonow(:,1,1)+m_norm(:))
     end if

!    Make the density positive everywhere (but do not care about gradients)
     call mkdenpos(iwarn,nfft,nspden_updn,1,rhonow(:,1:nspden_updn,1),dtset%xc_denpos)

!    write(std_out,*) 'rhonow',rhonow

!    Uses a block formulation, in order to save simultaneously
!    CPU time and memory : xc routines
!    are called only once over mpts times, while the amount of allocated
!    space is kept at a low value, even if a lot of different
!    arrays are allocated, for use in different xc functionals.
!    !!!  !$OMP PARALLEL LASTPRIVATE(mpts)
!    $OMP PARALLEL
!    $OMP SINGLE
#if defined HAVE_OPENMP
     mpts=8000 / OMP_GET_NUM_THREADS()
#else
     if (mgga) then
       mpts=nfft
     else
       mpts=4000
     end if
#endif
!    $OMP END SINGLE
!    $OMP END PARALLEL

!    The variable order indicates to which derivative of the energy
!    the computation must be done. Computing exc and vxc needs order=1 .
!    Meaningful values are 1, 2, 3. Lower than 1 is the same as 1, and larger
!    than 3 is the same as 3.
!    order=1 or 2 supported for all LSD and GGA ixc
!    order=3 supported only for ixc=3 and ixc=7
     order=1
     if(option==2.or.option==10.or.option==12)order=2
     if(option==-2)order=-2
     if(option==3)order=3

!    $OMP PARALLEL DO PRIVATE(coeff,grho,dstrsxc,vxcgrho_b,dvxc_b) &
!    $OMP&PRIVATE(d2vxc_b,exc_b,grho2_b_updn,ifft,index,ipts,ispden) &
!    $OMP&PRIVATE(npts,rho_b,rho_b_updn,strdiag) &
!    $OMP&PRIVATE(s1,s2,s3,vxcrho_b_updn) &
!    $OMP&REDUCTION(+:epsxc,enxc,strsxc1_tot,strsxc2_tot,strsxc3_tot) &
!    $OMP&REDUCTION(+:strsxc4_tot,strsxc5_tot,strsxc6_tot) &
!    $OMP&SHARED(depsxc,dtset%ixc,kxc,mpts,nfft,ngrad,nspden_updn,order) &
!    $OMP&SHARED(rhonow,vxc)
     do ifft=1,nfft,mpts
!      npts=mpts
!      npts is the number of points to be treated in this bunch
       npts=min(nfft-ifft+1,mpts)

!      Allocation of mandatory arguments of drivexc
       ABI_ALLOCATE(exc_b,(npts))
       ABI_ALLOCATE(rho_b,(npts))
       ABI_ALLOCATE(rho_b_updn,(npts,nspden_updn))
       ABI_ALLOCATE(vxcrho_b_updn,(npts,nspden_updn))
       vxcrho_b_updn(:,:)=zero
!      Allocation of optional arguments
       call size_dvxc(dtset%ixc,ndvxc,ngr2,nd2vxc,nspden_updn,nvxcgrho,order)
       if (ndvxc/=0)  then
         ABI_ALLOCATE(dvxc_b,(npts,ndvxc))
       end if
       if (nvxcgrho/=0)  then
         ABI_ALLOCATE(vxcgrho_b,(npts,nvxcgrho))
       end if
       ixc=dtset%ixc
       test_libxc=.false.
#if defined HAVE_DFT_LIBXC
       test_libxc= (ixc<0).and.(.not.(libxc_functionals_isgga() .or. libxc_functionals_ismgga()))
#endif
       if ((ixc==3 .or. (ixc>=7 .and. ixc<=15) .or. (ixc>=23 .and. ixc<=24) .or. test_libxc ) .and. order==3) then 
         ABI_ALLOCATE(d2vxc_b,(npts,nd2vxc))
       end if


       if (ngrad == 2)  then
         ABI_ALLOCATE(grho2_b_updn,(npts,ngr2))
       end if
       if (mgga) then
         ABI_ALLOCATE(lrho_b_updn,(npts,nspden_updn))
         ABI_ALLOCATE(vxclrho_b_updn,(npts,nspden_updn))
         ABI_ALLOCATE(tau_b_updn,(npts,nspden_updn))
         ABI_ALLOCATE(vxctau_b_updn,(npts,nspden_updn))
       end if

       do ipts=ifft,ifft+npts-1
!        index=ipts-ifft+1 varies from 1 to npts
         index=ipts-ifft+1
         rho_b(index)=rhonow(ipts,1,1)
         if(nspden_updn==1)then
           rho_b_updn(index,1)=rhonow(ipts,1,1)*half
           if(dtset%xclevel==2 .or. (ixc>=31 .and. ixc<=34))then
             grho2_b_updn(index,1)=quarter*(rhonow(ipts,1,2)**2+rhonow(ipts,1,3)**2+rhonow(ipts,1,4)**2)
           end if
           if (mgga) then
             tau_b_updn(index,1)=taur(ipts,1)*half
             lrho_b_updn(index,1)=lrhonow(ipts,1)*half
           end if
         else
           rho_b_updn(index,1)=rhonow(ipts,2,1)
           rho_b_updn(index,2)=rhonow(ipts,1,1)-rhonow(ipts,2,1)
           if(dtset%xclevel==2 .or. (ixc>=31 .and. ixc<=34))then
             grho2_b_updn(index,1)=rhonow(ipts,2,2)**2+   &
&             rhonow(ipts,2,3)**2+   &
&             rhonow(ipts,2,4)**2
             grho2_b_updn(index,2)=(rhonow(ipts,1,2)-rhonow(ipts,2,2))**2 +   &
&             (rhonow(ipts,1,3)-rhonow(ipts,2,3))**2 +   &
&             (rhonow(ipts,1,4)-rhonow(ipts,2,4))**2
             grho2_b_updn(index,3)=rhonow(ipts,1,2)**2+   &
&             rhonow(ipts,1,3)**2+   &
&             rhonow(ipts,1,4)**2
           end if
           if (mgga) then
             tau_b_updn(index,1)=taur(ipts,2)
             tau_b_updn(index,2)=taur(ipts,1)-taur(ipts,2)
             lrho_b_updn(index,1)=lrhonow(ipts,2)
             lrho_b_updn(index,2)=lrhonow(ipts,1)-lrhonow(ipts,2)
           end if
         end if

       end do

#if defined HAVE_DFT_LIBXC
       if (uselibxc) then
         if (mgga) then
           call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&           grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b,lrho_updn=lrho_b_updn,vxclrho=vxclrho_b_updn,&
&           tau_updn=tau_b_updn,vxctau=vxctau_b_updn)
         else if (libxc_functionals_isgga()) then
!          LHH,FL,GMR
           if (order**2 <= 1) then
             call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&             grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b)
           else
             call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&             grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b,dvxc=dvxc_b)
           end if
!          LHH,FL,GMR
         else
!          LHH,FL,GMR
           if (order**2 <= 1) then
             call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho)
           elseif (order**2 <= 4) then
             call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&             dvxc=dvxc_b)
           else
             call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&             dvxc=dvxc_b,d2vxc=d2vxc_b)
           end if
!          LHH,FL,GMR

         end if
       else
#endif
!        case with gradient
         if (dtset%xclevel==2)then
           if (order**2 <= 1 .or. ixc == 16 .or. ixc==17 .or. ixc==26 .or. ixc==27 ) then
             if (ixc /= 13) then
               call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&               grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b)
             else
               call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&               grho2_updn=grho2_b_updn)
             end if
           else if (order /= 3) then
             if (ixc /= 13) then
               call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&               dvxc=dvxc_b,grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b)
             else
               call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&               dvxc=dvxc_b,grho2_updn=grho2_b_updn)
             end if
           else if (order == 3) then
             if (ixc /= 13) then
               call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&               dvxc=dvxc_b,d2vxc=d2vxc_b,grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b)
             else
               call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&               dvxc=dvxc_b,d2vxc=d2vxc_b,grho2_updn=grho2_b_updn)
             end if
           end if
!          cases without gradient
         else
           if (order**2 <=1) then
             if (ixc>=31 .and. ixc<=34) then !fake mgga functionals for testing purpose only (based on LDA functional)
               call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,     &
&               grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b,lrho_updn=lrho_b_updn,vxclrho=vxclrho_b_updn,&
&               tau_updn=tau_b_updn,vxctau=vxctau_b_updn)
             else
               call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho)
             end if
           else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
             call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,      &
&             dvxc=dvxc_b,d2vxc=d2vxc_b)
           else
             call drivexc(exc_b,ixc,npts,nspden_updn,order,rho_b_updn,vxcrho_b_updn,ndvxc,ngr2,nd2vxc,nvxcgrho,      &
&             dvxc=dvxc_b)
           end if
         end if
#if defined HAVE_DFT_LIBXC
       end if
#endif

!      Accumulate enxc, strsxc and store vxc (and eventually kxc)
       dstrsxc=zero
       do ipts=ifft,ifft+npts-1
         index=ipts-ifft+1
         epsxc=epsxc+rho_b(index)*exc_b(index)  !will be normalized with respect to the volume later to get enxc ("bigexc").
         depsxc(ipts,1)=vxcrho_b_updn(index,1)
         if(nspden_updn==1)then
           strdiag=rho_b(index)*(exc_b(index)-vxcrho_b_updn(index,1))
         else if(nspden_updn==2)then
           depsxc(ipts,2)=vxcrho_b_updn(index,2)
!          Note : this is not the complete Vxc in the GGA case
           strdiag=rho_b(index)*exc_b(index) &
&           -rho_b_updn(index,1)*vxcrho_b_updn(index,1)&
&           -(rho_b(index)-rho_b_updn(index,1))*vxcrho_b_updn(index,2)
         end if
         dstrsxc=dstrsxc+strdiag

!        For GGAs, additional terms appear
!        (the LB functional does not lead to additional terms)
         if(ngrad==2 .and. dtset%ixc/=13)then

!          Treat explicitely spin up, spin down and total spin for spin-polarized
!          Will exit when ispden=1 is finished if non-spin-polarized
           do ispden=1,3

             if(nspden_updn==1 .and. ispden>=2)exit

!            If the norm of the gradient vanishes, then the different terms vanishes,
!            but the inverse of the gradient diverges, so skip the update.
             if(grho2_b_updn(index,ispden) < 1.0d-24) then
               depsxc(ipts,ispden+nspden_updn)=zero
               cycle
             end if

!            Compute the derivative of n.e_xc wrt the
!            spin up, spin down, or total density. In the non-spin-polarized
!            case take the coefficient that will be multiplied by the
!            gradient of the total density
             if(nspden_updn==1)then
!              !              Definition of vxcgrho_b changed in v3.3
               if (nvxcgrho == 3) then
                 coeff=half*vxcgrho_b(index,1) + vxcgrho_b(index,3)
               else
                 coeff=half*vxcgrho_b(index,1)
               end if
             else if(nspden_updn==2)then
               if (nvxcgrho == 3) then
                 coeff=vxcgrho_b(index,ispden)
               else if (ispden /= 3) then
                 coeff=vxcgrho_b(index,ispden)
               else if (ispden == 3) then
                 coeff=zero
               end if
             end if
             depsxc(ipts,ispden+nspden_updn)=coeff

!            Store the gradient of up, down or total density, depending on ispden and nspden, at point ipts
             if(nspden_updn==1)then
               grho(1:3)=rhonow(ipts,1,2:4)
             else if(ispden==1 .and. nspden_updn==2)then
               grho(1:3)=rhonow(ipts,2,2:4)
             else if(ispden==2 .and. nspden_updn==2)then
               grho(1:3)=rhonow(ipts,1,2:4)-rhonow(ipts,2,2:4)
             else if(ispden==3 .and. nspden_updn==2)then
               grho(1:3)=rhonow(ipts,1,2:4)
             end if

!            In case of ixc 31 (mGGA functional fake 1),
!            skip the stress tensor to follow a LDA scheme (see doc/theory/MGGA/report_MGGA.pdf)
             if(ixc==31) cycle

!            Compute the contribution to the stress tensor
             s1=-grho(1)*grho(1)*coeff
             s2=-grho(2)*grho(2)*coeff
             s3=-grho(3)*grho(3)*coeff
!            The contribution of the next line comes from the part of Vxc
!            obtained from the derivative wrt the gradient
             dstrsxc=dstrsxc+s1+s2+s3
             strsxc1_tot=strsxc1_tot+s1
             strsxc2_tot=strsxc2_tot+s2
             strsxc3_tot=strsxc3_tot+s3
             strsxc4_tot=strsxc4_tot-grho(3)*grho(2)*coeff
             strsxc5_tot=strsxc5_tot-grho(3)*grho(1)*coeff
             strsxc6_tot=strsxc6_tot-grho(2)*grho(1)*coeff

           end do
         end if

!        For meta-GGAs, add the laplacian term (vxclrho) and kinetic energy density term (vxctau)
         if(mgga)then
           if(nspden_updn==1)then
             depsxc(ipts,3)   = vxclrho_b_updn(index,1)
             if(present(vxctau))then
               vxctau(ipts,1,1) = vxctau_b_updn(index,1)
             end if
           else if (nspden_updn==2)then
             depsxc(ipts,6)   = vxclrho_b_updn(index,1)
             depsxc(ipts,7)   = vxclrho_b_updn(index,2)
             if(present(vxctau))then
               vxctau(ipts,1,1) = vxctau_b_updn(index,1)
               vxctau(ipts,2,1) = vxctau_b_updn(index,2)
             end if
           end if
         end if

       end do

!      Additional electron-positron correlation terms
       if (ipositron==2) then
!        Compute electron-positron XC energy per unit volume, potentials and derivatives
         ngr=0;if (ngrad_apn==2) ngr=npts
         ABI_ALLOCATE(fxc_apn,(npts))
         ABI_ALLOCATE(vxc_b_apn,(npts))
         ABI_ALLOCATE(vxcgr_apn,(ngr))
         ABI_ALLOCATE(vxc_ep,(npts))
         ABI_ALLOCATE(rhonow_apn,(npts,nspden_apn,1))
         ABI_ALLOCATE(grho2_apn,(ngr))
         rhonow_apn(1:npts,1,1)=electronpositron%rhor_ep(ifft:ifft+npts-1,1)
         if (usexcnhat==0) rhonow_apn(1:npts,1,1)=rhonow_apn(1:npts,1,1)-electronpositron%nhat_ep(ifft:ifft+npts-1,1)
         if (.not.electronpositron%posdensity0_limit) then
           call mkdenpos(iwarnp,npts,nspden_apn,1,rhonow_apn(:,1,1),dtset%xc_denpos)
         end if
         if (ngrad_apn==2.and.ngr2==1) grho2_apn(:)=four*grho2_b_updn(:,1)
         if (ngrad_apn==2.and.ngr2==3) grho2_apn(:)=grho2_b_updn(:,3)
         if (ndvxc==0) then
           call xcpositron(fxc_apn,grho2_apn,electronpositron%ixcpositron,ngr,npts,&
&           electronpositron%posdensity0_limit,rho_b,&
&           rhonow_apn(:,1,1),vxc_b_apn,vxcgr_apn,vxc_ep)
         else
           ABI_ALLOCATE(dvxc_apn,(npts))
           call xcpositron(fxc_apn,grho2_apn,electronpositron%ixcpositron,ngr,npts,&
&           electronpositron%posdensity0_limit,rho_b,&
&           rhonow_apn(:,1,1),vxc_b_apn,vxcgr_apn,vxc_ep,dvxce=dvxc_apn)
         end if
         ABI_DEALLOCATE(vxc_ep)
         ABI_DEALLOCATE(rhonow_apn)
         ABI_DEALLOCATE(grho2_apn)
!        Accumulate electron-positron XC energies
         s1=zero
         do ipts=1,npts
           s1=s1+fxc_apn(ipts)
         end do
         electronpositron%e_xc=electronpositron%e_xc+s1*ucvol/dble(nfftot)
!        Add electron-positron dVxc_el/dRho_el to electron-electron one
         if (ndvxc==1) dvxc_b(:,1)=dvxc_b(:,1)+dvxc_apn(:)
         if (ndvxc==3) then
           dvxc_b(:,1)=dvxc_b(:,1)+four*dvxc_apn(:)
           dvxc_b(:,2)=dvxc_b(:,2)+four*dvxc_apn(:)
           dvxc_b(:,3)=dvxc_b(:,3)+four*dvxc_apn(:)
         end if
         if (ndvxc==15) then
           dvxc_b(:, 9)=dvxc_b(:, 9)+four*dvxc_apn(:)
           dvxc_b(:,10)=dvxc_b(:,10)+four*dvxc_apn(:)
           dvxc_b(:,11)=dvxc_b(:,11)+four*dvxc_apn(:)
!          if (ngr>0) then
!          dvxc_b(:,12)=dvxc_b(:,12)+vxcgr_apn(:)
!          dvxc_b(:,13)=
!          dvxc_b(:,14)=
!          dvxc_b(:,15)=
!          end if
         end if
!        Modify stresses - Compute factors for GGA
         do ipts=ifft,ifft+npts-1
           index=ipts-ifft+1
           depsxc_apn(ipts,1)=vxc_b_apn(index)
           dstrsxc=dstrsxc+fxc_apn(index)-rho_b(index)*vxc_b_apn(index)
           if (ngrad_apn==2) then
             depsxc_apn(ipts,2)=vxcgr_apn(index)
             s1=-grho(1)*grho(1)*vxcgr_apn(index)
             s2=-grho(2)*grho(2)*vxcgr_apn(index)
             s3=-grho(3)*grho(3)*vxcgr_apn(index)
             dstrsxc=dstrsxc+s1+s2+s3
             strsxc1_tot=strsxc1_tot+s1
             strsxc2_tot=strsxc2_tot+s2
             strsxc3_tot=strsxc3_tot+s3
             strsxc4_tot=strsxc4_tot-grho(3)*grho(2)*vxcgr_apn(index)
             strsxc5_tot=strsxc5_tot-grho(3)*grho(1)*vxcgr_apn(index)
             strsxc6_tot=strsxc6_tot-grho(2)*grho(1)*vxcgr_apn(index)
           end if ! GGA
         end do ! ipts
!        Deallocations
         ABI_DEALLOCATE(fxc_apn)
         ABI_DEALLOCATE(vxc_b_apn)
         ABI_DEALLOCATE(vxcgr_apn)
         if (ndvxc>0) ABI_DEALLOCATE(dvxc_apn)
       end if

!      Transfer the xc kernel
       if(nkxc/=0 .and. nkxc/=23)then
         if (ndvxc==15.and.(option==10.or.option==12)) then
           if (nkxc>=3) then
             kxc(ifft:ifft+npts-1,1)=dvxc_b(1:npts,1)+dvxc_b(1:npts,9)
             kxc(ifft:ifft+npts-1,2)=dvxc_b(1:npts,10)
             kxc(ifft:ifft+npts-1,3)=dvxc_b(1:npts,2)+dvxc_b(1:npts,11)
             if (nkxc>3) kxc(ifft:ifft+npts-1,4:nkxc)=zero
           else
             kxc(ifft:ifft+npts-1,1)=half*(dvxc_b(1:npts,1)+dvxc_b(1:npts,9)+dvxc_b(1:npts,10))
             if (nkxc>1) kxc(ifft:ifft+npts-1,2:nkxc)=zero
           end if
         else if (nkxc<=ndvxc) then
           kxc(ifft:ifft+npts-1,1:nkxc)=dvxc_b(1:npts,1:nkxc)
         else
           if (allocated(dvxc_b)) kxc(ifft:ifft+npts-1,1:ndvxc)=dvxc_b(1:npts,1:ndvxc) !if ndvxc=0 dvxc_b is not allocated (PMA)
           kxc(ifft:ifft+npts-1,ndvxc+1:nkxc)=zero
         end if
         if (order==3) k3xc(ifft:ifft+npts-1,1:nd2vxc)=d2vxc_b(1:npts,1:nd2vxc)
       else if(nkxc==23)then
         if (ndvxc==15) then
           kxc(ifft:ifft+npts-1,1:15)=dvxc_b(1:npts,1:15)
         else
           kxc(ifft:ifft+npts-1,1:ndvxc)=dvxc_b(1:npts,1:ndvxc)
           kxc(ifft:ifft+npts-1,ndvxc+1:15)=zero
         end if
         do ispden=1,nspden_updn
           do ii=1,4
             kxc(ifft:ifft+npts-1,13+ispden+2*ii)=rhonow(ifft:ifft+npts-1,ispden,ii)
           end do
         end do
       end if

!      Add the diagonal part to the xc stress
       strsxc1_tot=strsxc1_tot+dstrsxc
       strsxc2_tot=strsxc2_tot+dstrsxc
       strsxc3_tot=strsxc3_tot+dstrsxc

       ABI_DEALLOCATE(exc_b)
       ABI_DEALLOCATE(rho_b)
       ABI_DEALLOCATE(rho_b_updn)
       ABI_DEALLOCATE(vxcrho_b_updn)
       if (allocated(dvxc_b))  then
         ABI_DEALLOCATE(dvxc_b)
       end if
       if (allocated(vxcgrho_b))  then
         ABI_DEALLOCATE(vxcgrho_b)
       end if
       if (allocated(d2vxc_b))  then
         ABI_DEALLOCATE(d2vxc_b)
       end if
       if (allocated(grho2_b_updn))  then
         ABI_DEALLOCATE(grho2_b_updn)
       end if
       if (allocated(vxclrho_b_updn))  then
         ABI_DEALLOCATE(vxclrho_b_updn)
       end if
       if (allocated(lrho_b_updn))  then
         ABI_DEALLOCATE(lrho_b_updn)
       end if
       if (allocated(tau_b_updn))  then
         ABI_DEALLOCATE(tau_b_updn)
       end if
       if (allocated(vxctau_b_updn))  then
         ABI_DEALLOCATE(vxctau_b_updn)
       end if

!      End of the loop on blocks of data
     end do

!    $OMP END PARALLEL DO

     strsxc(1)=strsxc1_tot
     strsxc(2)=strsxc2_tot
     strsxc(3)=strsxc3_tot
     strsxc(4)=strsxc4_tot
     strsxc(5)=strsxc5_tot
     strsxc(6)=strsxc6_tot

!    If GGA, multiply the gradient of the density by the proper
!    local partial derivatives of the XC functional
     if (ipositron==2) then
       ABI_ALLOCATE(rhonow_ptr,(nfft,nspden_eff,ngrad*ngrad))
       rhonow_ptr=rhonow
     else
       rhonow_ptr => rhonow
     end if
     if(ngrad==2 .and. dtset%ixc/=13)then
       call xcmult(depsxc,nfft,ngrad,nspden_eff,nspgrad,rhonow_ptr)
     end if

!    Compute contribution from this grid to vxc, and ADD to existing vxc
     if (nspden/=4) then
       if(present(vxctau))then
         call xcpot(cplex,depsxc,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,nspgrad,&
&         dtset%paral_kgb,qphon,rhonow_ptr,vxc,vxctau=vxctau)
       else
         call xcpot(cplex,depsxc,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,nspgrad,&
&         dtset%paral_kgb,qphon,rhonow_ptr,vxc)
       end if

     else

!      If non-collinear magnetism, restore potential in proper axis before adding it
       ABI_ALLOCATE(vxcrho_b_updn,(nfft,4))
       vxcrho_b_updn=zero
       call xcpot(cplex,depsxc,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,nspgrad,&
&       dtset%paral_kgb,qphon,rhonow_ptr,vxcrho_b_updn)
       if (usexcnhat==1.or.nhatdim==0) then
         do ifft=1,nfft
           dvdn=half*(vxcrho_b_updn(ifft,1)+vxcrho_b_updn(ifft,2))
           if(m_norm(ifft)>m_norm_min) then
!            if(m_norm(ifft)>rhor(ifft,1)*tol10+tol14) then
             dvdz=half*(vxcrho_b_updn(ifft,1)-vxcrho_b_updn(ifft,2))/m_norm(ifft)
             vxc(ifft,1)=vxc(ifft,1)+dvdn+rhor(ifft,4)*dvdz
             vxc(ifft,2)=vxc(ifft,2)+dvdn-rhor(ifft,4)*dvdz
             vxc(ifft,3)=vxc(ifft,3)+rhor(ifft,2)*dvdz
             vxc(ifft,4)=vxc(ifft,4)-rhor(ifft,3)*dvdz
           else
             vxc(ifft,1:2)=vxc(ifft,1:2)+dvdn
           end if
         end do
       else
         do ifft=1,nfft
           dvdn=half*(vxcrho_b_updn(ifft,1)+vxcrho_b_updn(ifft,2))
           if(m_norm(ifft)>m_norm_min) then
!            if(m_norm(ifft)>rhor(ifft,1)*tol10+tol14) then
             dvdz=half*(vxcrho_b_updn(ifft,1)-vxcrho_b_updn(ifft,2))/m_norm(ifft)
             vxc(ifft,1)=vxc(ifft,1)+dvdn+(rhor(ifft,4)-nhat(ifft,4))*dvdz
             vxc(ifft,2)=vxc(ifft,2)+dvdn-(rhor(ifft,4)-nhat(ifft,4))*dvdz
             vxc(ifft,3)=vxc(ifft,3)+(rhor(ifft,2)-nhat(ifft,2))*dvdz
             vxc(ifft,4)=vxc(ifft,4)-(rhor(ifft,3)-nhat(ifft,3))*dvdz
           else
             vxc(ifft,1:2)=vxc(ifft,1:2)+dvdn
           end if
         end do
       end if
       ABI_DEALLOCATE(vxcrho_b_updn)
     end if
     if (ipositron==2)  then
       ABI_DEALLOCATE(rhonow_ptr)
     end if
     nullify(rhonow_ptr)

!    Add electron-positron XC potential to electron-electron one
!    Eventually compute GGA contribution
     if (ipositron==2) then
       ABI_ALLOCATE(rhonow_apn,(nfft,nspden_apn,ngrad_apn**2))
       rhonow_apn(1:nfft,1,1:ngrad_apn**2)=rhonow(1:nfft,1,1:ngrad_apn**2)
       if (ngrad_apn==2) then
         call xcmult(depsxc_apn,nfft,ngrad_apn,nspden_apn,ngrad_apn,rhonow_apn)
       end if
       ABI_ALLOCATE(vxc_apn,(nfft,nspden_apn))
       vxc_apn=zero
       call xcpot(cplex,depsxc_apn,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad_apn,&
&       nspden_apn,ngrad_apn,dtset%paral_kgb,qphon,rhonow_apn,vxc_apn)
       vxc(:,1)=vxc(:,1)+vxc_apn(:,1)
       if (nspden_updn==2) vxc(:,2)=vxc(:,2)+vxc_apn(:,1)
       s1=zero
       do ipts=1,nfft
         s1=s1+vxc_apn(ipts,1)*rhonow(ipts,1,1)
       end do
       electronpositron%e_xcdc=electronpositron%e_xcdc+s1*ucvol/dble(nfftot)
       ABI_DEALLOCATE(rhonow_apn)
       ABI_DEALLOCATE(vxc_apn)
       ABI_DEALLOCATE(depsxc_apn)
     end if

!    End loop on unshifted or shifted grids
   end do

!  Normalize enxc, strsxc and vxc
   divshft=one/dble(dtset%intxc+1)
   strsxc(:)=strsxc(:)/dble(nfftot)*divshft
   if (dtset%usewvl == 0) then
     enxc=epsxc*ucvol/dble(nfftot)*divshft
   else
     enxc = epsxc * (dtset%wvl_hgrid / real(2, dp)) ** 3 * divshft
   end if
   do ispden=1,nspden
     do ifft=1,nfft
       vxc(ifft,ispden)=vxc(ifft,ispden)*divshft
       if(present(vxctau)) vxctau(ifft,ispden,1:4)=vxctau(ifft,ispden,1:4)*divshft
     end do
   end do

!  XG030514 : MPIWF Should reduce strsxc and enxc
!  in the group of WF processors
   if(mpi_enreg%paral_compil_fft==1)then
     old_paral_level=mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
     call timab(48,1,tsec)
     call xsum_mpi(strsxc,spaceComm ,ierr)
     call xsum_mpi(enxc  ,spaceComm ,ierr)
     if (ipositron==2) then
       s1=electronpositron%e_xc;s2=electronpositron%e_xcdc
       call xsum_mpi(s1,spaceComm ,ierr)
       call xsum_mpi(s2,spaceComm ,ierr)
       electronpositron%e_xc=s1;electronpositron%e_xcdc=s2
     end if
     call timab(48,2,tsec)
     mpi_enreg%paral_level=old_paral_level
   end if

!  Compute vxcavg
   call mean_fftr(vxc,vxcmean,mpi_enreg,nfft,nfftot,min(nspden,2))

!  write(std_out,*) 'vxcmean',vxcmean

   if(nspden==1)then
     vxcavg=vxcmean(1)
   else
     vxcavg=half*(vxcmean(1)+vxcmean(2))
   end if

   ABI_DEALLOCATE(depsxc)
   ABI_DEALLOCATE(rhonow)
   if (allocated(lrhonow))  then
     ABI_DEALLOCATE(lrhonow)
   end if
   if (nspden==4)  then
     ABI_DEALLOCATE(m_norm)
   end if

!  TO BE DELETED (unuseful test)
!  else if (dtset%ixc/=20) then
!  
!  !  Not an allowed choice for ixc and nspden
!  write(message, '(a,a,a,a,2i8,a)' ) ch10,&
!  &   ' rhohxc: BUG -',ch10,&
!  &   '  ixc,nspden=',dtset%ixc,nspden,' not supported.'
!  call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')

 end if

!Treat separately the Fermi-Amaldi correction.
 if (dtset%ixc==20 .or. dtset%ixc==21 .or. dtset%ixc==22) then

!  Fermi-Amaldi correction : minus Hartree divided by the
!  number of electrons per unit cell. This is not size consistent, but
!  interesting for isolated systems with a few electrons.
   nelect=ucvol*rhog(1,1)
   factor=-one/nelect
   vxc(:,1)=factor*vhartr(:)
   if(nspden>=2) then
     vxc(:,2)=factor*vhartr(:)
   end if

!  Compute corresponding xc energy and stress as well as vxcavg
   call dotprod_vn(1,rhor,enxc,doti,mpi_enreg,nfft,nfftot,1,1,vxc,ucvol)
   enxc=half*enxc
   strsxc(1:3)=-enxc/ucvol

!  Compute average of vxc (one component only).
   call mean_fftr(vxc,vxcmean,mpi_enreg,nfft,nfftot,1)
   vxcavg = vxcmean(1)
!  For ixc=20, the local exchange-correlation kernel is zero, but the Hartree
!  kernel will be modified in tddft. No other use of kxc should be made with ixc==20
   if(nkxc/=0 .and. dtset%ixc==20) then
     do inkxc=1,nkxc
!      $OMP PARALLEL DO PRIVATE(ifft) SHARED(inkxc,kxc,nfft,nkxc)
       do ifft=1,nfft
         kxc(ifft,inkxc)=zero
       end do
!      $OMP END PARALLEL DO
     end do
   end if
!  For ixc=21 or 22, the LDA (ixc=1) kernel has been computed previously.

 end if

 call timab(81,2,tsec)

!DEBUG
!write(std_out,*)' rhohxc : enxc=',enxc
!if(present(vxctau))write(std_out,*)' rhohxc : max vxctau=',maxval(vxctau(:,1,1))
!if(present(vxctau))write(std_out,*)' rhohxc : min vxctau=',minval(vxctau(:,1,1))
!if(present(vxctau))write(std_out,*)' rhohxc : max gvxctau=',maxval(vxctau(:,1,2)),maxval(vxctau(:,1,3)),maxval(vxctau(:,1,4))
!if(present(vxctau))write(std_out,*)' rhohxc : min gvxctau=',minval(vxctau(:,1,2)),minval(vxctau(:,1,3)),minval(vxctau(:,1,4))
!stop
!write(std_out,*)' rhohxc : strsxc(1:3)=',strsxc(1:3)
!ENDDEBUG

!DEBUG
!write(std_out,'(a)') ' rhohxc :  '
!write(std_out,'(a)')&
!& '   ir              rhor(ir)      vxc(ir)       kxc(ir)     xccc3d(ir) '
!n1=ngfft(1)
!n2=ngfft(2)
!n3=ngfft(3)
!do ir=1,nfft
!! if(ir<=11 .or. mod(ir,301)==0 )then
!i3=(ir-1)/n1/n2
!i2=(ir-1-i3*n1*n2)/n1
!i1=ir-1-i3*n1*n2-i2*n1
!if(n3xccc>0)write(message,'(i5,3i3,a,4es14.6)')ir,i1,i2,i3,' ',&
!&   rhor(ir,1),vxc(ir,1),kxc(ir,1),xccc3d(ir)
!if(n3xccc==0)write(message,'(i5,3i3,a,4es14.6)')ir,i1,i2,i3,' ',&
!&   rhor(ir,1),vxc(ir,1),kxc(ir,1)
!call wrtout(std_out,message,'COLL')
!if(nspden_updn==2)then
!write(message,'(a,2es14.6)')'               ',rhor(ir,2),vxc(ir,2)
!call wrtout(std_out,message,'COLL')
!end if
!! end if
!end do
!stop
!ENDDEBUG

!DEBUG
!sum=0.0
!do jj=1,nspden
!do kk=1,nfft
!sum=sum+abs(vxc(kk,jj))
!end do
!end do
!write(std_out,*) 'SUM(vxc)',sum
!write(std_out,*)' rhohxc : debug, stop at exit '
!enxc=zero
!kxc(:,:)=zero
!strsxc(:)=zero
!vhartr(:)=zero
!vxc(:,:)=zero
!vxcavg=zero
!stop
!ENDDEBUG

!call leave_new("COLL")

end subroutine rhohxc

!An optimization of the xc part of this routine
!and the subroutines called by this routine
!has been performed already.
!The number of exponentiations, logarithms, divisions etc has been minimized.
!The results of tests 11-18, 21-24 are as follows, on a P6 at 200 MHz.
!All tests are with intxc=0, except the last one.
!A further optimization should focus on intxc=1 .
!Time is in microsec/data/call
!Test    time     Note
!11     2.643    non spin-pol Teter xcspol
!12     3.103    id           Perdew Zunger
!13     2.836    id           Teter xctetr
!14     2.670    id           Wigner
!15     3.815    id           Hedin-Lundqvist
!16     2.365    id           X-alpha
!17     4.476    id           Perdew Wang
!18    12.894    id           GGA PBE
!21     4.001    spin-pol     Teter xcspol
!22     6.012    id           Perdew Wang
!23    12.797    id           GGA PBE
!24    25.953    id           GGA PBE with intxc=1
!
!Tests where also lead to determine the time for some basic operations :
!(these are very approximative numbers)
!div     0.22
!sqrt    0.3
!log     0.5
!exp     1.3
!**dblep 3.0
!The routine invcb, that compute x**(-one/3.0d0), scores 1.4  .

!!***
