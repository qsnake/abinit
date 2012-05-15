!{\src2tex{textfont=tt}}
!!****f* ABINIT/prctfw3
!! NAME
!! prctfw3
!!
!! FUNCTION
!! Compute new trial potential by applying the Thomas--Fermi--von Weizsaecker
!! charge mixing scheme (see PRB 64 121101).
!! First step is to compute a localy averaged non local potential
!! Written starting from src/5common/energy.F90
!! This function tries a third algorithm (compare with prctfvw1 and 2) and is not efficient enough
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! old trial potential
!! old output potential
!! new trial potential resulting from the mixing choice
!! old output density
!!
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of wavefunction
!!   operator (ground-state symmetries)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eew=Ewald energy (hartree)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  eii=psp core-core energy
!!  entropy=entropy due to the occupation number smearing (if metal)
!!  epaw=PAW spherical part energy
!!  epawdc=PAW spherical part double-counting energy
!!  gsqcut=G^2 cutoff from gsqcut=ecut/(2 Pi^2)
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimension for number of planewaves
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  dtset%nfft=(effective) number of FFT grid points (for  real(dp), intent(out)  :: resid(mband*nkpt*nsppol)
!!  dtset%nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (dtset%nfftf=dtset%nfft for norm-conserving potential runs)
!!  dtset%ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#dtset%ngfft
!!  dtset%ngfft(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (dtset%ngfft=dtset%ngfft for norm-conserving potential runs)
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  ntypat=number of types of atoms in cell
!!  n3xccc=dimension of the xccc3d array (0 or dtset%nfftf).
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each k point
!!  occopt=option for occupancies
!!  optene=option for the computation of total energy (direct scheme or double-counting scheme)
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  tsmear=smearing energy or temperature (if metal)
!!  vpsp(dtset%nfftf)=local pseudopotential in real space (hartree)
!!  wffnow=structured array giving all information about wavefunction file
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!! modified new trial potential based on the TFvW charge mixing
!!
!! SIDE EFFECTS
!!  rhog(2,dtset%nfftf)=work space for rho(G); save intact on return
!!  rhor(dtset%nfftf,nspden)=work space for rho(r); save intact on return
!!  nspinor should not be modified in the call of rdnpw
!!
!! WARNINGS
!! This is experimental code : input, ouptput, results and any other feature may vary greatly.
!!
!! NOTES
!!
!! PARENTS
!!      prcref_PMA
!!
!! CHILDREN
!!      cgpr,dotprod_vn,fourdp,ftfvw1__end,ftfvw1__init,laplacian,mean_fftr
!!      metric,rhotov
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prctfw3(deltae,dtset,&
     & efermi,etotal,gsqcut,&
     & lavnlr,mpi_enreg,&
     & nhat,nhatgr,nhatgrdim,&
     & nkxc,n3xccc,optene,optxc,&
     & psps,rhor_in,rprimd,&
     & usexcnhat,&
     & vpsp,vresid,vrespc,vtrial,wvl,&
     & xccc3d,xred)

 use m_profiling

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use m_energies,only:energies_type
  use ftfvw1
  use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prctfw3'
 use interfaces_42_geometry
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_62_cg_noabirule
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n3xccc,nhatgrdim,nkxc,optene,optxc,usexcnhat
 real(dp),intent(in) :: etotal,gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
  !no_abirules
  !dtset%nfft**(1-1/nsym) is 1 if nsym==1, and dtset%nfft otherwise
  ! WARNING
  ! BEWARE THERE IS TWO DIFFERENT SIZE DECLARED FOR ARRAY NHAT IN RHOTOV AND RHOHXC
  ! THIS MIGHT RESULT IN A BUG
  real(dp),intent(in)   ::   nhat(dtset%nfft,dtset%nspden*psps%usepaw),nhatgr(dtset%nfft,dtset%nspden,3*nhatgrdim)
  real(dp), intent(in)  :: rhor_in(dtset%nfft,dtset%nspden)
  real(dp), intent(in)   :: rprimd(3,3),xred(3,dtset%natom)
  real(dp), intent(in)   :: lavnlr(dtset%nfft,dtset%nspden)
  !real(dp), intent(inout):: vin_old(dtset%nfft,dtset%nspden),vout_unmixed(dtset%nfft,dtset%nspden),vtrial(dtset%nfft,dtset%nspden)
  real(dp), intent(in)    :: vresid(dtset%nfft,dtset%nspden),vtrial(dtset%nfft,dtset%nspden)
  real(dp), intent(out)   :: vrespc(dtset%nfft,dtset%nspden)
  real(dp), intent(inout):: vpsp(dtset%nfft)
  real(dp),dimension(:),intent(inout)    :: xccc3d(n3xccc)
  real(dp), intent(in) :: deltae,efermi

!Local variables-------------------------------
  !previously intent(out) without any reasonable justifications? 21/08/2006
  !end of previously intent(out)
!scalars
 integer :: count,cplex,ifft,ispden,option
 real(dp),parameter :: alpha32=(3._dp*pi*pi)
 real(dp) :: Z,c1,c2,doti,dummy,dummy2,maxi,mini,smear
 real(dp) :: ucvol,vme,vres2,vxcavg
! character(len=fnlen) :: filapp
 type(energies_type) :: energies
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),qphon(3),rhog(2,dtset%nfft)
 real(dp) :: rhor(dtset%nfft,dtset%nspden),rmet(3,3),strsxc(6),vmean(2)
 real(dp) :: vnew_mean(2),vres_mean(2),vtfw(dtset%nfft,dtset%nspden)
 real(dp) :: vtrialold(dtset%nfft,dtset%nspden)
 real(dp) :: xred2(size(dtset%xred_orig,1),size(dtset%xred_orig,2))
 real(dp),allocatable :: deltaW(:,:),g2cart(:),kxc(:,:),laplacerhor(:,:)
 real(dp),allocatable :: newvout(:,:),newvoutfourier(:,:,:),sqrtrhor(:,:)
 real(dp),allocatable :: vhartr(:),vin_oldfourier(:,:,:),vtrialfourier(:,:,:)
 real(dp),allocatable :: vxc(:,:)

! *************************************************************************

!***********************************************************************************
!Getting the localy averaged non-local potential                                ***
!$Vnl(r) = [\sum_{n,k} f_{n,k} \psi_{n,k}(r) (Vnl(r,r') |\psi_{n,k}(r')>)]/n(r)$***
!**********************************************************************************
 rhor=rhor_in
 qphon=zero
 xred2=xred
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!******************************************************************
!Getting the DeltaW factor                                      **
!******************************************************************
 ABI_ALLOCATE(deltaW,(dtset%nfft,dtset%nspden))
 ABI_ALLOCATE(laplacerhor,(dtset%nfft,dtset%nspden))
 ABI_ALLOCATE(sqrtrhor,(dtset%nfft,dtset%nspden))
 ABI_ALLOCATE(vtrialfourier,(2,dtset%nfft,dtset%nspden))
 ABI_ALLOCATE(newvoutfourier,(2,dtset%nfft,dtset%nspden))
 ABI_ALLOCATE(g2cart,(dtset%nfft))
!Compute the real space laplacian of sqrt(rhor)
 sqrtrhor(:,:)=(rhor(:,:))**half
 call laplacian(gprimd,mpi_enreg,dtset%nfft,dtset%nspden,dtset%ngfft,dtset%paral_kgb,rdfuncr=sqrtrhor,&
& laplacerdfuncr=laplacerhor,g2cart_out=g2cart)
!second step: get deltaW
 do ispden=1,dtset%nspden
   do ifft=1,dtset%nfft
     deltaW(ifft,ispden)=((efermi&
&     - half*(alpha32 * (rhor(ifft,ispden))  )**(two_thirds)&
&     - vtrial(ifft,ispden) - lavnlr(ifft,ispden))& ! this one is the true one
&    * sqrtrhor(ifft,ispden)&
&     + one*half*laplacerhor(ifft,ispden) )
   end do
 end do

!******************************************************************
!Finding the density which minimizes the associated Energy      **
!******************************************************************
!compute the total charge number
 cplex=1;
 option=1;
 call dotprod_vn(cplex,& !complex density/pot
&sqrtrhor,&          !the density
&Z,&  !resulting dorproduct integrated over r
&doti,&          !imaginary part of the integral
&mpi_enreg,&     !
&size(rhor,1),&          !number of localy(cpu) attributed grid point
&dtset%nfft,&        !real total number of grid point
&size(rhor,2),&        !dtset%nspden
&option,&        !1=compute only the real part 2=compute also the imaginary part
&sqrtrhor,&          !the potential
&ucvol)          !cell volume
 Z=real(nint(Z),dp)
!enable the use of the functions eneofrho_tfw and deneofrho_tfw
 call ftfvw1__init(dtset,dtset%intxc,dtset%ixc,psps%usepaw,n3xccc,dtset%ngfft,&
& dtset%nfft,&
& nhat,nhatgr,nhatgrdim,&
& nkxc,dtset%nspden,mpi_enreg,deltaW,gprimd,gsqcut,&
& lavnlr,rhor,rprimd,ucvol,&
& psps%usepaw,usexcnhat,&
& vtrial+vresid,&
& vpsp,vtrial,xccc3d,Z)
!minimizes Etfw with respect to sqrtrhor instead of rhor
 call cgpr(size(rhor,1),size(rhor,2),ftfvw1__e,ftfvw1__de,ftfvw1__newdensity,&
& abs(deltae*real(0.0001,dp)/etotal),55,sqrtrhor,dummy,dummy2)
!free the dynamically allocated memory used by eneofrho_tfw and deneofrho_tfw
 call ftfvw1__end()
!new density from the minimised sqrtrhor
 count=0
 do ifft=1,dtset%nfft
   if (sqrtrhor(ifft,1)<zero) then
     count=count+1
   end if
 end do
 rhor(:,:)=sqrtrhor(:,:)*sqrtrhor(:,:)
!rhor(:,:)=    half*rhor(:,:) +&
!&   (one-half)*rhor_in(:,:)
!******************************************************************
!Production of a new trial potential                            **
!******************************************************************
!write(3330,*) vtrial
 ABI_ALLOCATE(newvout,(dtset%nfft,dtset%nspden))
 ABI_ALLOCATE(vin_oldfourier,(2,dtset%nfft,dtset%nspden))
!step1: make V from rho
 call fourdp(1, rhog(:,:), rhor(:,1),-1,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%paral_kgb,0) !
 ABI_ALLOCATE(vhartr,(dtset%nfft))
 ABI_ALLOCATE(vxc,(dtset%nfft,dtset%nspden))
 ABI_ALLOCATE(kxc,(dtset%nfft,1))
 call rhotov(dtset,energies,gprimd,gsqcut,kxc,mpi_enreg,dtset%nfft,dtset%ngfft,&
& nhat(:,1:dtset%nspden*psps%usepaw),&
& nhatgr,nhatgrdim,  &
& nkxc,vrespc,&
& n3xccc,optene,1,optxc,&
& rhog,rhor,rprimd,strsxc,ucvol,&
& psps%usepaw,usexcnhat,vhartr,vnew_mean,vpsp,&
& vres_mean,vres2,vtfw,vxcavg,vxc,wvl,xccc3d)

 if(allocated(kxc))ABI_DEALLOCATE(kxc)
!shift the newvout to remove the constant part. Don't know if this is OK
!WARNING
 call mean_fftr(vtfw,vmean,mpi_enreg,dtset%nfft,dtset%nfft,min(dtset%nspden,2))
 do ispden=1,min(dtset%nspden,2)
   if(dtset%nspden/=2 .or. &
&   ( dtset%occopt>=3 .and. abs(dtset%fixmom+real(99.99,dp))<real(1.0d-10,dp) ))then
     if(dtset%nspden==1)then
       vme=vmean(1)
     else
       vme=(vmean(1)+vmean(2))*half
     end if
   else
     vme=vmean(ispden)
   end if
   vtfw(:,ispden)=vtfw(:,ispden)-vme
 end do
 if(allocated(vhartr))ABI_DEALLOCATE(vhartr)
 if(allocated(vxc))ABI_DEALLOCATE(vxc)

!step2: mix the new Vout with the old Vout
!mixing based on the kinetic energy at the k points
!change potential to fourier space
 do ispden=1,dtset%nspden
   call fourdp(1, newvoutfourier(:,:,ispden), vtfw(:,ispden),-1,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%paral_kgb,0)
   vtfw(:,ispden)=vtrial(:,ispden)+vresid(:,ispden)
   call fourdp(1, vtrialfourier(:,:,ispden), vtfw(:,ispden),-1,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%paral_kgb,0)
!  call fourdp(1, vin_oldfourier(:,:,ispden), vtrial(:,ispden),-1,dtset%nfft,dtset%ngfft,dtset%paral_kgb,0)     !modified for printing...
 end do
!filtering
 c1=0
 c2=0

!write(2201,*) "beginning of step",istep
 do ifft=1,dtset%nfft
   if(g2cart(ifft) < 2.0) then
!    write(2201,*) g2cart(ifft),vtrialfourier(:,ifft,1),newvoutfourier(:,ifft,1)
   end if
 end do
!write(2201,*) "end of step",istep
 close(2201)

 do ispden=1,dtset%nspden
   do ifft=1,dtset%nfft
!    vtrialfourier(:,ifft,ispden) = newvoutfourier(:,ifft,ispden)*exp(-g2cart(ifft)*two) &
!    &+(one-exp(-g2cart(ifft)*two))*(&
!    &vin_oldfourier(:,ifft,ispden)*(one-exp(0.8_dp*g2cart(ifft)/(g2cart(ifft)+half)))& !! mixing proposed by
!    &+vtrialfourier(:,ifft,ispden)*(exp(0.8_dp*g2cart(ifft)/(g2cart(ifft)+half))))    !! raczowski...
!    & vtrialfourier(:,ifft,ispden)) !no further mixing
!    vtrialfourier(:,ifft,ispden) = newvoutfourier(:,ifft,ispden)*exp(-g2cart(ifft)*12.56637061435917d0) &
!    &+(one-exp(-g2cart(ifft)*12.56637061435917d0))*(&
!    & vtrialfourier(:,ifft,ispden))
!    vtrialfourier(:,ifft,ispden) = newvoutfourier(:,ifft,ispden)*exp(-g2cart(ifft)*36.0d0) &
!    &+(one-exp(-g2cart(ifft)*36.0d0))*(&
!    & vtrialfourier(:,ifft,ispden))
     smear=4.0d0*(exp(-1.0d0*g2cart(ifft))+exp(1.0d0*g2cart(ifft)))**(-2)
     smear=zero
     maxi=0.02
     mini=0.004

     if (g2cart(ifft) <maxi) then
       if (g2cart(ifft) < mini) then
         smear=one-0.d0
         c1=c1+1
       else
         smear=(((g2cart(ifft)-maxi)/(mini-maxi))**2)*(one-0.d0)
         c2=c2+1
       end if
     end if



     vtrialfourier(:,ifft,ispden) = &
!    & (g2cart(ifft)/(g2cart(ifft)+0.5))*&
&     newvoutfourier(:,ifft,ispden)*smear &
&     +(one-smear)* (vtrialfourier(:,ifft,ispden))



   end do
 end do
!change resulting potential to real space
!write(3363,*) vtrialfourier
 vtrialold=vtrial
 do ispden=1,dtset%nspden
!  call fourdp(1,vtrialfourier(:,:,ispden),newvout(:,ispden),1,dtset%nfft,dtset%ngfft,dtset%paral_kgb,0)    !output modified->no unholy effect
   call fourdp(1,vtrialfourier(:,:,ispden),vtfw(:,ispden),1,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%paral_kgb,0)
 end do

!shift the newvout to remove the constant part. Don't know if this is OK
!WARNING
 call mean_fftr(vtfw,vmean,mpi_enreg,dtset%nfft,dtset%nfft,min(dtset%nspden,2))
 do ispden=1,min(dtset%nspden,2)
   if(dtset%nspden/=2 .or. &
&   ( dtset%occopt>=3 .and. abs(dtset%fixmom+real(99.99,dp))<real(1.0d-10,dp) ))then
     if(dtset%nspden==1)then
       vme=vmean(1)
     else
       vme=(vmean(1)+vmean(2))*half
     end if
   else
     vme=vmean(ispden)
   end if
   vtfw(:,ispden)=vtfw(:,ispden)-vme
 end do
 vrespc=vtfw-vtrial

 call laplacian(gprimd,mpi_enreg,dtset%nfft,dtset%nspden,dtset%ngfft,dtset%paral_kgb,rdfuncr=vtfw,&
 laplacerdfuncr=rhor,g2cart_in=g2cart)
 rhor=-rhor*4.0d0*pi

!----------------------------------------------------------------------------------------------------------
!last FREE

 if(allocated(deltaW))  then
   ABI_DEALLOCATE(deltaW)
 end if
 if(allocated(newvout))ABI_DEALLOCATE(newvout)
 if(allocated(newvoutfourier))ABI_DEALLOCATE(newvoutfourier)
 if(allocated(vtrialfourier))ABI_DEALLOCATE(vtrialfourier)
 if(allocated(sqrtrhor))ABI_DEALLOCATE(sqrtrhor)
 if(allocated(laplacerhor))ABI_DEALLOCATE(laplacerhor)
 if(allocated(g2cart))ABI_DEALLOCATE(g2cart)
 if(allocated(vin_oldfourier))ABI_DEALLOCATE(vin_oldfourier)
end subroutine prctfw3
!!***
