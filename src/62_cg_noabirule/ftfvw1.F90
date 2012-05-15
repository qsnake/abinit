!{\src2tex{textfont=tt}}
!!****m* ABINIT/ftfvw1
!! NAME
!! ftfvw1
!!
!! FUNCTION
!! Module providing the ability to compute the energy and its first derivative
!! from a density
!! this routine uses the thomas--fermi--von Weisaecker
!! approximation for the kinetic energy
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! This is a module.
!! It is made of two functions and one init subroutine
!!
!! NOTES
!!
!! PARENTS
!! prctfw
!!
!! CHILDREN
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module ftfvw1

 use m_profiling

  use defs_basis
  use defs_abitypes
  use interfaces_53_spacepar
  use interfaces_56_recipspace  ! THIS IS MANDATORY TO CALL LAPLACIAN

  implicit none

  !! common variables copied from input
  integer,private                  :: intxc,ixc,izero
  integer,private                  :: n3xccc
  integer,private                  :: ngfft(18),nfftf
  integer,private                  :: nhatgrdim,usepaw,usexcnhat
  integer,private                  :: nkxc
  integer,private                  :: nspden
  type(dataset_type),private       :: dtset
  type(MPI_type),private           :: mpi_enreg

  real(dp),allocatable,private     :: deltaW(:,:)
  real(dp),private                 :: gprimd(3,3)
  real(dp),private                 :: gsqcut
  real(dp),allocatable,private     :: lavnlfft(:,:)
  real(dp),allocatable,private     :: nhat(:,:),nhatgr(:,:,:)
  real(dp),private                 :: rprimd(3,3)
  real(dp),private                 :: ucvol
  real(dp),allocatable,private     :: vpion(:,:)
  real(dp),allocatable,private     :: xccc3d(:)
  real(dp),private                 :: Z
  !! common variables computed
  logical,private :: ok=.false.
  integer,private :: nfftot
  real(dp),private:: alpha
  real(dp),allocatable,private     :: rhog(:,:),rhogsp(:,:,:),kxc(:,:)

contains
!!***

!!****f* ftfvw1/ftfvw1__init
!! NAME
!! ftfvw1__init
!!
!! FUNCTION
!! initialisation subroutine
!! Copy every variables required for the energy calculation
!! Allocate the required memory
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      prctfvw1,prctfw3
!!
!! CHILDREN
!!      dotprod_vn,fourdp,laplacian,rhohxc
!!
!! SOURCE
  subroutine ftfvw1__init(dtset_in,intxc_in,ixc_in,izero_in,n3xccc_in,&
       &ngfft_in,nfftf_in,&
       & nhat_in,nhatgr_in,nhatgrdim_in,&
       & nkxc_in,nspden_in,mpi_enreg_in,deltaW_in,gprimd_in,&
       & gsqcut_in,lavnlfft_in,rhor,rprimd_in,ucvol_in,&
       & usepaw_in,usexcnhat_in,&
       & vout_unmixed,vpsp,vtrial_in,xccc3d_in,Z_in )
    !use defs_basis
    !use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ftfvw1__init'
 use interfaces_53_ffts
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    integer,intent(in) :: intxc_in,ixc_in,izero_in
    integer,intent(in) :: n3xccc_in
    integer,intent(in) :: ngfft_in(18),nfftf_in
    integer,intent(in) :: nkxc_in
    integer,intent(in) :: nspden_in
    integer,intent(in) :: nhatgrdim_in,usepaw_in,usexcnhat_in
    type(MPI_type),intent(in) :: mpi_enreg_in
    type(dataset_type),intent(in) :: dtset_in

    real(dp),intent(in)  :: deltaW_in(nfftf_in,nspden_in)
    real(dp),intent(in)  :: gprimd_in(3,3)
    real(dp),intent(in)  :: gsqcut_in
    real(dp),intent(in)  :: lavnlfft_in(nfftf_in,nspden_in)
    real(dp),intent(inout):: rhor(nfftf_in,nspden_in)
    real(dp),intent(in)  :: rprimd_in(3,3)
    real(dp),intent(in)  :: ucvol_in
    real(dp),intent(in)  :: vout_unmixed(nfftf_in,nspden_in)
    real(dp),intent(in)  :: vpsp(nfftf_in)
    real(dp),intent(in)  :: vtrial_in(nfftf_in,nspden_in)
    real(dp),intent(in)  :: xccc3d_in(n3xccc_in)
    real(dp),intent(in)  :: Z_in
    real(dp),intent(in) :: nhat_in(nfftf_in,nspden_in*usepaw_in),nhatgr_in(nfftf_in,nspden_in,3*nhatgrdim_in)

    !Local variables-------------------------------
    integer:: cplex,option,ispden
!   real(dp):: doti,dotr
!   real(dp),allocatable :: vhartr(:),vxc(:,:)
!   real(dp)             :: dummy_strsxc(6)
    real(dp),parameter   :: identity(4)=(/1.0_dp,1.0_dp,0.0_dp,0.0_dp/)
    !!allocation and data transfer
    !!Thought it would have been more logical to use the privates intrinsic of the module as
    !!input variables it seems that it is not possible...
    if(.not.ok) then
       dtset = dtset_in
       intxc=intxc_in
       ixc=ixc_in
       izero=izero_in
       n3xccc=n3xccc_in
       ngfft=ngfft_in
       nfftf=nfftf_in
       nkxc=nkxc_in
       nspden=nspden_in
       ABI_ALLOCATE(deltaW,(size(deltaW_in,1),size(deltaW_in,2)))
       ABI_ALLOCATE(lavnlfft,(size(lavnlfft_in,1),size(lavnlfft_in,2)))
       ABI_ALLOCATE(vpion,(size(vtrial_in,1),size(vtrial_in,2)))
       ABI_ALLOCATE(xccc3d,(size(xccc3d_in)))
       ABI_ALLOCATE(kxc,(1,0))
       ABI_ALLOCATE(rhog,(2,nfftf))
       ABI_ALLOCATE(rhogsp,(2,nfftf,nspden))
       ABI_ALLOCATE(nhat,(nfftf,nspden*usepaw))
       ABI_ALLOCATE(nhatgr,(nfftf,nspden,3*nhatgrdim))
       deltaW=deltaW_in
       gprimd=gprimd_in
       gsqcut=gsqcut_in
       lavnlfft=lavnlfft_in
       usepaw=usepaw_in
       nhatgrdim=nhatgrdim_in
       nhat=nhat_in
       nhatgr=nhatgr_in
       usexcnhat=usexcnhat_in
       rprimd=rprimd_in
       ucvol=ucvol_in
!!! COMPUTE VPION = Vext + Vtrial - Vout_unmixed   !!!
       !allocate(vhartr(nfftf),vxc(nfftf,nspden))
       call fourdp(1, rhog(:,:), rhor(:,1),-1,mpi_enreg,nfftf,ngfft,dtset_in%paral_kgb,0)
       cplex=1
       option=1
       !call rhohxc(doti,gsqcut,izero,kxc,mpi_enreg,nfftf,&
       !     &  ngfft,nkxc,nspden,n3xccc,option,rhog,rhor,rprimd,dummy_strsxc,vhartr,&
       !     &  vxc,dotr,xccc3d)
       ! compute the different part of the lagrangian derivative
       do ispden=1,nspden
          vpion(:,ispden)=vpsp(:)*identity(ispden)
          !vpion(:,ispden)=vout_unmixed(:,ispden)&
          !     &-vxc(:,ispden)&
          !     &-vhartr(:)*identity(ispden)
       end do
       vpion(:,:)=vpion(:,:)+vout_unmixed(:,:)-vtrial_in(:,:)
       !call random_number(vpion(:,:))
       !write(2222,*) vpion
       !write(2223,*) vtrial_in
       !write(2224,*) vout_unmixed
       !write(2225,*) vtrial_in(:,1)&
       !    &-vxc(:,1)&
       !       &-vhartr(:)*identity(1)
       !write(2226,*) vtrial_in-vout_unmixed
!!!!!!!!!!!!!!END OF VION COMPUTATION!!!!!!!!!!!!!!!!!!!!!!!!!
       xccc3d=xccc3d_in
       mpi_enreg=mpi_enreg_in
       gprimd(:,:)=gprimd_in(:,:)
       !! set ok to 1 which allow using eneofrho_tfw
       ok = .true.
       !! alpha is a constant factor used many times in the tfw calculation of energy
       alpha=(3._dp*pi*pi)**two_thirds
       !!total number of grid point (from energy.F90)
       nfftot=ngfft(1)*ngfft(2)*ngfft(3)
       Z=Z_in
       !deallocate(vhartr,vxc)
    end if
  end subroutine ftfvw1__init
!!***
  
!!****f* ftfvw1/ftfvw1__end
!! NAME
!! ftfvw1__end
!!
!! FUNCTION
!! ending subroutine
!! deallocate memory areas
!!  
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      prctfvw1,prctfw3
!!
!! CHILDREN
!!      dotprod_vn,fourdp,laplacian,rhohxc
!!
!! SOURCE
  subroutine ftfvw1__end()
    use defs_basis
    use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ftfvw1__end'
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    !Local variables-------------------------------
    if(ok) then
       !! set ok to false which prevent using eneofrho_tfw
       ok = .false.
       !! free memory
       ABI_DEALLOCATE(deltaW)
       ABI_DEALLOCATE(lavnlfft)
       ABI_DEALLOCATE(vpion)
       ABI_DEALLOCATE(xccc3d)
       ABI_DEALLOCATE(kxc)
       ABI_DEALLOCATE(rhog)
       ABI_DEALLOCATE(rhogsp)
       ABI_DEALLOCATE(nhat)
       ABI_DEALLOCATE(nhatgr)

    end if
  end subroutine ftfvw1__end
!!***

!!****f* ftfvw1/ftfvw1__newdensity
!! NAME
!! ftfvw1__newdensity
!!
!! FUNCTION
!! affectation subroutine
!! do the required renormalisation when providing a new value for
!! the density after application of the gradient
!!  
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      dotprod_vn,fourdp,laplacian,rhohxc
!!
!! SOURCE
  subroutine ftfvw1__newdensity(nv1,nv2,x, grad, sqrtrho)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ftfvw1__newdensity'
 use interfaces_53_spacepar
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    integer,intent(in)    :: nv1,nv2
    real(dp),intent(in)   :: x
    real(dp),intent(inout):: grad(nv1,nv2)
    real(dp),intent(inout):: sqrtrho(nv1,nv2)
    !Local variables-------------------------------
    real(dp):: store(nv1,nv2)
    integer :: cplex,option
    real(dp):: doti,dotr
    cplex=1
    option=1
    store=sqrtrho
    sqrtrho(:,:)=sqrtrho(:,:)+x*grad(:,:)
    call dotprod_vn(cplex,& !complex density/pot
         &sqrtrho,&          !the density
         &dotr,&  !resulting dorproduct integrated over r
         &doti,&          !imaginary part of the integral
         &mpi_enreg,&     !
         &size(sqrtrho,1),&          !number of localy(cpu) attributed grid point
         &nfftot,&        !real total number of grid point
         &size(sqrtrho,2),&        !nspden
         &option,&        !1=compute only the real part 2=compute also the imaginary part
         &sqrtrho,&          !the potential
         &ucvol)          !cell volume
    if(Z*dotr.lt.zero) stop 'newdensity: the density renormalisation is negative...'
    sqrtrho(:,:)=sqrtrho(:,:)*sqrt(Z/dotr)
    !sqrtrho(:,:)=sqrt(sqrtrho(:,:)*sqrtrho(:,:)*Z/dotr)
    grad=sqrtrho-store
  end subroutine ftfvw1__newdensity
!!***

!!****f* ftfvw1/ftfvw1__e
!! NAME
!! ftfvw1__e
!!
!! FUNCTION
!! energy
!!  
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE
  function ftfvw1__e(nv1,nv2,sqrtrho_in)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ftfvw1__e'
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_56_xc
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    integer,intent(in)  :: nv1,nv2
    real(dp),intent(in) ::sqrtrho_in(nv1,nv2)
    real(dp)            :: ftfvw1__e
    !Local variables-------------------------------
    integer ::option,cplex,ispden,nk3xc
    real(dp)::doti,dotr,enxc,vxcavg
    real(dp)::fv(nv1,nv2)
    real(dp)::laplacerhotothehalf(nv1,nv2)
    real(dp)::sqrtrho(nv1,nv2),vxc(nv1,nv2),rho(nv1,nv2)
    real(dp)::vhartr(nv1)
    real(dp)::dummy_strsxc(6)
    real(dp),parameter ::    identity(4)=(/1.0_dp,1.0_dp,0.0_dp,0.0_dp/)

    if(ok) then
       !renormalizing the density
       option=1
       cplex=1
       call dotprod_vn(cplex,& !complex density/pot
            &sqrtrho_in,&          !the sqrt density
            &dotr,&  !resulting dorproduct integrated over r
            &doti,&          !imaginary part of the integral
            &mpi_enreg,&     !
            &size(rho,1),&          !number of localy(cpu) attributed grid point
            &nfftot,&        !real total number of grid point
            &size(rho,2),&        !nspden
            &option,&        !1=compute only the real part 2=compute also the imaginary part
            &sqrtrho_in,&          !the sqrt density
            &ucvol)          !cell volume
       if(Z*dotr.lt.zero)  stop 'eneofrho_tfw: the density renormalisation is negative...'
       rho(:,:)=sqrtrho_in(:,:)*sqrtrho_in(:,:)!*(Z/dotr)
       sqrtrho(:,:)=sqrtrho_in(:,:)!*sqrt(Z/dotr)


       call dotprod_vn(cplex,& !complex density/pot
            &sqrtrho,&          !the sqrt density
            &dotr,&  !resulting dorproduct integrated over r
            &doti,&          !imaginary part of the integral
            &mpi_enreg,&     !
            &size(rho,1),&          !number of localy(cpu) attributed grid point
            &nfftot,&        !real total number of grid point
            &size(rho,2),&        !nspden
            &option,&        !1=compute only the real part 2=compute also the imaginary part
            &sqrtrho,&          !the sqrt density
            &ucvol)          !cell volume


       !sqrtrho(:,:)=sqrt(rho(:,:))

       call laplacian(gprimd,mpi_enreg,nfftf,nspden,ngfft,dtset%paral_kgb,rdfuncr=sqrtrho,laplacerdfuncr=laplacerhotothehalf)
       call fourdp(1, rhog(:,:), rho(:,1),-1,mpi_enreg,nfftf,ngfft,dtset%paral_kgb,0)
       option=1 ; nk3xc=1
       call rhohxc(dtset,enxc,gsqcut,izero,kxc,mpi_enreg,nfftf,ngfft,&
            &  nhat,usepaw,nhatgr,nhatgrdim,&
            &  nkxc,nk3xc,nspden,n3xccc,option,rhog,rho,rprimd,dummy_strsxc,usexcnhat,vhartr,&
            &  vxc,vxcavg,xccc3d)

       do ispden=1,nspden
          fv(:,ispden)= (&
               & sqrtrho(:,ispden)*(&
               & lavnlfft(:,ispden) *one&
               & + vpion(:,ispden)  *one&
               & + vhartr(:)*identity(ispden)*half *one&
               & )+three_fifth*alpha*rho(:,ispden)**(one+sixth)*half *one &
               & -one*half*laplacerhotothehalf(:,ispden) *one &
               & +deltaW(:,ispden)*two *one&
                                !& + merge(zero,sign(real(1e22,dp),(sqrtrho(:,ispden))),(sqrtrho(:,ispden).ge.tol16)) *two &
               & )
       end do

       option=1
       cplex=1
       call dotprod_vn(cplex,& !complex density/pot
            &sqrtrho,&          !the sqrt density
            &ftfvw1__e,&  !resulting dotproduct integrated over r
            &doti,&          !imaginary part of the integral
            &mpi_enreg,&     !
            &size(rho,1),&          !number of localy(cpu) attributed grid point
            &nfftot,&        !real total number of grid point
            &size(rho,2),&        !nspden
            &option,&        !1=compute only the real part 2=compute also the imaginary part
            &fv,&          !the potential
            &ucvol)          !cell volume

       !sum       
       ftfvw1__e = ftfvw1__e + enxc
       !do ifftf=1,nfftf
       ! do ispden=1,nspden
       ! eneofrho_tfw=eneofrho_tfw+&
       !& merge(zero,((sqrtrho(ifftf,ispden)-epsilon)*0.1_dp/epsilon)**2,sqrtrho(ifftf,ispden).ge.epsilon)
       !    merge(zero,1e4_dp*((sqrtrho(ifftf,ispden)))**2,sqrtrho(ifftf,ispden).ge.zero)
       ! end do
       ! end do


       !_________________________________________________________
       !COMPONENT OF THE ENERGY
       !_________________________________________________________
!!!      write(std_err,*) 'Etotal=',eneofrho_tfw
!!!       write(std_err,*) 'Exc=',enxc
!!!       ispden=1
!!!       fv(:,ispden)=sqrtrho(:,ispden)*lavnlfft(:,ispden)
!!!       call dotprod_vn(cplex,& !complex density/pot
!!!            &sqrtrho,&          !the sqrt density
!!!            &dotr,&  !resulting dotproduct integrated over r
!!!            &doti,&          !imaginary part of the integral
!!!            &mpi_enreg,&     !
!!!            &size(rho,1),&          !number of localy(cpu) attributed grid point
!!!            &nfftot,&        !real total number of grid point
!!!            &size(rho,2),&        !nspden
!!!            &option,&        !1=compute only the real part 2=compute also the imaginary part
!!!            &fv,&          !the potential
!!!            &ucvol)          !cell volume
!!!       write(std_err,*) 'Elavnl=',dotr
!!!       fv(:,ispden)=sqrtrho(:,ispden)*vpion(:,ispden)
!!!       call dotprod_vn(cplex,& !complex density/pot
!!!            &sqrtrho,&          !the sqrt density
!!!            &dotr,&  !resulting dotproduct integrated over r
!!!            &doti,&          !imaginary part of the integral
!!!            &mpi_enreg,&     !
!!!            &size(rho,1),&          !number of localy(cpu) attributed grid point
!!!            &nfftot,&        !real total number of grid point
!!!            &size(rho,2),&        !nspden
!!!            &option,&        !1=compute only the real part 2=compute also the imaginary part
!!!            &fv,&          !the potential
!!!            &ucvol)          !cell volume
!!!       write(std_err,*) 'Epion=',dotr
!!!       fv(:,ispden)=sqrtrho(:,ispden)*vhartr(:)*identity(ispden)*half
!!!       call dotprod_vn(cplex,& !complex density/pot
!!!            &sqrtrho,&          !the sqrt density
!!!            &dotr,&  !resulting dotproduct integrated over r
!!!            &doti,&          !imaginary part of the integral
!!!            &mpi_enreg,&     !
!!!            &size(rho,1),&          !number of localy(cpu) attributed grid point
!!!            &nfftot,&        !real total number of grid point
!!!            &size(rho,2),&        !nspden
!!!            &option,&        !1=compute only the real part 2=compute also the imaginary part
!!!            &fv,&          !the potential
!!!            &ucvol)          !cell volume
!!!       write(std_err,*) 'Ehartr=',dotr
!!!       fv(:,ispden)=three_fifth*alpha*rho(:,ispden)**(one+sixth)*half
!!!       call dotprod_vn(cplex,& !complex density/pot
!!!            &sqrtrho,&          !the sqrt density
!!!            &dotr,&  !resulting dotproduct integrated over r
!!!            &doti,&          !imaginary part of the integral
!!!            &mpi_enreg,&     !
!!!            &size(rho,1),&          !number of localy(cpu) attributed grid point
!!!            &nfftot,&        !real total number of grid point
!!!            &size(rho,2),&        !nspden
!!!            &option,&        !1=compute only the real part 2=compute also the imaginary part
!!!            &fv,&          !the potential
!!!            &ucvol)          !cell volume
!!!       write(std_err,*) 'Ekrho=',dotr
!!!       fv(:,ispden)=-half*laplacerhotothehalf(:,ispden)*fifth
!!!       call dotprod_vn(cplex,& !complex density/pot
!!!            &sqrtrho,&          !the sqrt density
!!!            &dotr,&  !resulting dotproduct integrated over r
!!!            &doti,&          !imaginary part of the integral
!!!            &mpi_enreg,&     !
!!!            &size(rho,1),&          !number of localy(cpu) attributed grid point
!!!            &nfftot,&        !real total number of grid point
!!!            &size(rho,2),&        !nspden
!!!            &option,&        !1=compute only the real part 2=compute also the imaginary part
!!!            &fv,&          !the potential
!!!            &ucvol)          !cell volume
!!!       write(std_err,*) 'Eknab=',dotr
!!!       fv(:,ispden)=deltaW(:,ispden)*two
!!!       call dotprod_vn(cplex,& !complex density/pot
!!!            &sqrtrho,&          !the sqrt density
!!!            &dotr,&  !resulting dotproduct integrated over r
!!!            &doti,&          !imaginary part of the integral
!!!            &mpi_enreg,&     !
!!!            &size(rho,1),&          !number of localy(cpu) attributed grid point
!!!            &nfftot,&        !real total number of grid point
!!!            &size(rho,2),&        !nspden
!!!            &option,&        !1=compute only the real part 2=compute also the imaginary part
!!!            &fv,&          !the potential
!!!            &ucvol)          !cell volume
!!!       write(std_err,*) 'Ew=',dotr
!!$
    else
       ftfvw1__e=zero
    end if
  end function ftfvw1__e
!!***

!!****f* ftfvw1/ftfvw1__de
!! NAME
!! ftfvw1__de
!!
!! FUNCTION
!! derivative of the energy
!! actually not the derivative but something allowing minimization
!! at constant density
!! formula from the work of rackowski,canning and wang
!! H*phi - int(phi**2H d3r)phi
!!  
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE
  function ftfvw1__de(nv1,nv2,sqrtrho_in)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ftfvw1__de'
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_56_xc
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    integer,intent(in) :: nv1,nv2
    real(dp),intent(in):: sqrtrho_in(nv1,nv2)
    real(dp)           :: ftfvw1__de(nv1,nV2)
    !Local variables-------------------------------
    integer  :: option,cplex,ispden,nk3xc
    real(dp) :: doti,dotr
    real(dp) :: fv(nv1,nv2)
    real(dp) :: laplacerhotothehalf(nv1,nv2)
    real(dp) :: sqrtrho(nv1,nv2),fh(nv1,nv2),vxc(nv1,nv2),rho(nv1,nv2)
    real(dp) :: vhartr  (nv1)
    real(dp) :: dummy_strsxc(6)
    real(dp),parameter ::    identity(4)=(/1.0_dp,1.0_dp,0.0_dp,0.0_dp/)
    if(ok) then
       !renormalizing the density
       option=1
       cplex=1
       call dotprod_vn(cplex,& !complex density/pot
            &sqrtrho_in,&          !the sqrt density
            &dotr,&  !resulting dorproduct integrated over r
            &doti,&          !imaginary part of the integral
            &mpi_enreg,&     !
            &size(rho,1),&          !number of localy(cpu) attributed grid point
            &nfftot,&        !real total number of grid point
            &size(rho,2),&        !nspden
            &option,&        !1=compute only the real part 2=compute also the imaginary part
            &sqrtrho_in,&          !the sqrt density
            &ucvol)          !cell volume
       if(Z*dotr.le.zero)  stop 'the density renormalisation is negative...'
       rho(:,:)=(sqrtrho_in(:,:)*sqrtrho_in(:,:)*Z/dotr)
       sqrtrho(:,:)=sqrtrho_in(:,:)*sqrt(Z/dotr)
       !sqrtrho(:,:)=sqrt(rho(:,:))
       call fourdp(1, rhog(:,:), rho(:,1),-1,mpi_enreg,nfftf,ngfft,dtset%paral_kgb,0)   !get rhog
       option=1 ; nk3xc=1
       call rhohxc(dtset,doti,gsqcut,izero,kxc,mpi_enreg,nfftf,&
            &  ngfft,&
            &  nhat,usepaw,nhatgr,nhatgrdim,&
            &  nkxc,nk3xc,nspden,n3xccc,option,rhog,rho,rprimd,dummy_strsxc,usexcnhat,vhartr,&
            &  vxc,dotr,xccc3d)      !get vhartr and vxc
       ! compute the different part of the lagrangian derivative
       call laplacian(gprimd,mpi_enreg,nfftf,nspden,ngfft,dtset%paral_kgb,rdfuncr=sqrtrho,laplacerdfuncr=laplacerhotothehalf)
       do ispden=1,nspden
          fv(:,ispden)= lavnlfft(:,ispden)   *one&
               & + vpion(:,ispden)           *one&
               & + vxc(:,ispden)             *one&
               & + vhartr(:)*identity(ispden)*one&
               & + alpha*rho(:,ispden)**(two_thirds)*half *one
       end do
       fh(:,:)=(&
            & fv(:,:)*sqrtrho(:,:)                 &
            & -one*half*laplacerhotothehalf(:,:)   *one &
            & +deltaW(:,:)                         *one&
                                !& + merge(zero,1e4_dp,sqrtrho(:,:).ge.zero) * sqrtrho(:,:)&
            & )


       option=1
       cplex=1
       call dotprod_vn(cplex,&
            &sqrtrho,&          !the density
            &dotr,&  !resulting dorproduct integrated over r
            &doti,&          !imaginary part of the integral
            &mpi_enreg,&     !
            &nfftf,&          !number of localy(cpu) attributed grid point
            &nfftot,&        !real total number of grid point
            &nspden,&        !nspden
            &option,&        !1=compute only the real part 2=compute also the imaginary part
            &fh,&          !the potential
            &ucvol)          !cell volume
       ftfvw1__de(:,:)=(fh(:,:)-(dotr*sqrtrho(:,:)/Z))
    else
       ftfvw1__de = zero
    end if
    !write(std_err,*) 'deneofrho 4'
  end function ftfvw1__de


end module ftfvw1
!!***
