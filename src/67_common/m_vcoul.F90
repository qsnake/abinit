!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_vcoul
!! NAME
!!  m_vcoul
!!
!! FUNCTION
!!  This module contains the definition of the vcoul_t as well
!!  as procedures to calculate the Coulomb interaction in reciprocal 
!!  space taking into account a possible cutoff in real space.
!!  Procedures to deal with the singularity for q-->0 are also provided.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG, FB, GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_vcoul

 use m_profiling

 use defs_basis
 use m_errors 
 use m_xmpi

 use m_gwdefs,          only : GW_TOLQ0
 use m_io_tools,        only : get_unit
 use m_numeric_tools,   only : arth, geop, imin_loc, llsfit_svd, l2norm, OPERATOR(.x.), quadrature
 use m_geometry,        only : normv
 use m_bz_mesh,         only : bz_mesh_type, get_BZ_item
 use m_gsphere,         only : gvectors_type

 implicit none

 private 
!!***

!!****t* m_vcoul/vcoul_t
!! NAME
!!  vcoul_t
!!
!! FUNCTION
!!  This data type contains the square root of the Fourier components of the Coulomb interaction
!!  calculated taking into account a possible cutoff. Info on the particular geometry used for the cutoff
!!  as well as quantities required to deal with the Coulomb divergence are reported as well.
!!
!! SOURCE

 type,public :: vcoul_t

  integer  :: nfft                   ! Number of points in FFT grid
  integer  :: ng                     ! Number of G-vectors
  integer  :: nqibz                  ! Number of irreducible q-points
  integer  :: nqlwl                  ! Number of small q-points around Gamma

  real(dp) :: alpha(3)               ! Lenght of the finite surface
  real(dp) :: rcut                   ! Cutoff radius
  real(dp) :: i_sz                   ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2
  real(dp) :: hcyl                   ! Lenght of the finite cylinder along the periodic dimension
  real(dp) :: ucvol                  ! Volume of unit cell

  character(len=50) :: mode          ! String defining the cutoff mode, possible values are: sphere,cylinder,surface,crystal

  integer :: pdir(3)                 ! 1 if the system is periodic along this direction
  integer  :: ngfft(18)              ! Information of FFT grid

                                     ! 1 if the point in inside the cutoff region 0 otherwise
  real(dp) :: boxcenter(3)           ! Reduced coordinates of the center of the box (input variable)
  real(dp) :: vcutgeo(3)             ! For each reduced direction gives the lenght of the finite system
                                     ! 0 if the system is infinite along that particular direction
                                     ! negative value to indicate that a finite size has to be used
  real(dp) :: rprimd(3,3)

  integer(i1b),pointer :: ctab(:)    SET2NULL
  ! For each point of the FFT grid this table reports

  real(dp),pointer :: qibz(:,:)   SET2NULL
  ! qibz(3,nqibz)
  ! q-points in the IBZ.

  real(dp),pointer :: qlwl(:,:)   SET2NULL
  ! qibz(3,nqlwl)
  ! q-points for the treatment of the Coulomb singularity.

  complex(gwpc),pointer :: vc_sqrt(:,:)  SET2NULL
  ! vc_sqrt(ng,nqibz)
  ! Square root of the Coulomb in reciprocal space.
  ! A cut might be applied.

  complex(gwpc),pointer :: vcqlwl_sqrt(:,:)  SET2NULL
  ! vcqs_sqrt(ng,nqlwl)
  ! Square root of the Coulomb calculated for small q-points

 end type vcoul_t
!!***

 public ::  vcoul_nullify        ! Nullify the data type.
 public ::  vcoul_init           ! Main creation method. 
 public ::  vcoul_plot           ! Plot vc in real and reciprocal space.
 public ::  vcoul_print          ! Report info on the object.
 public ::  cutoff_table         ! Create table defining wheter a FFT point lyes inside the cutoff region.
 public ::  cutoff_density       ! Redefine the density according to the cutoff mode.
 public ::  vcoul_free           ! Destruction method.
 public ::  cutoff_sphere        ! FT of the Coulomb truncated within a sphere.
 public ::  cutoff_cylinder      ! FT of the Coulomb truncated within a cylinder
 public ::  cutoff_surface       ! FT of the Coulomb truncated within a slab.
 public ::  cmod_qpg             ! FT of the long ranged Coulomb interaction.

! local private variables used for the integration needed by the cylindrical case.
 integer :: npts_,ntrial_,qopt_
 real(dp) :: ha_,hb_,r0_
 real(dp) :: qpg_perp_,qpg_para_,qpgx_,qpgy_
 real(dp) :: zz_,xx_,rho_
 real(dp) :: hcyl_,rcut_,accuracy_

CONTAINS  !========================================================================================
!!***

!!****f* m_vcoul/vcoul_init
!! NAME
!! vcoul_init
!!
!! FUNCTION
!! Perform general check and initialize the data type containing information on the cutoff technique
!!
!! INPUTS
!!  Qmesh<Bz_mesh_type>=Info on the q-point sampling.
!!  Kmesh<Bz_mesh_type>=Info on the k-point sampling.
!!  Gsph<Gvectors_type>=Info of the G sphere.
!!    %gmet(3,3)=Metric in reciprocal space.
!!    %gprimd(3,3)=Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!    %gvec=G vectors
!!  rcut=Cutoff radius for the cylinder.
!!  icutcoul=Option of the cutoff technique.
!!  vcutgeo(3)= Info on the orientation and extension of the cutoff region.
!!  ng=Number of G-vectors to be used to describe the Coulomb interaction
!!  nqlwl=Number of point around Gamma for treatment of long-wavelength limit
!!  qlwl(3,nqlwl)= The nqlwl "small" q-points 
!!  ngfft(18)=Information on the (fine) FFT grid used for the density.
!!  rprimd(3,3)=Direct lattice vectors.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Vcp=<vcoul_t>=Datatype gathering information on the Coulomb interaction.
!!
!! PARENTS
!!      m_shexc,mrgscr,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_init(Vcp,Gsph,Qmesh,Kmesh,rcut,icutcoul,vcutgeo,ng,nqlwl,qlwl,rprimd,ngfft,comm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vcoul_init'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,nqlwl,icutcoul,comm
 real(dp),intent(in) :: rcut
 type(Bz_mesh_type),intent(in) :: Kmesh,Qmesh
 type(Gvectors_type),intent(in) :: Gsph
 type(vcoul_t),intent(out) :: Vcp
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rprimd(3,3),qlwl(3,nqlwl),vcutgeo(3)

!Local variables-------------------------------
!scalars
 integer :: ii,iqlwl,iq_bz,iq_ibz,istat,npar,npt
 integer :: opt_cylinder,opt_surface,test,master,rank,nprocs
 real(dp),parameter :: tolq0=1.d-3
 real(dp) :: bz_geometry_factor,bz_plane,check,dx,integ,q0_vol,q0_volsph
 real(dp) :: qbz_norm,step,ucvol
 character(len=500) :: msg
!arrays
 integer :: gamma_pt(3,1)
 integer,pointer :: gvec(:,:)
 real(dp) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),gmet(3,3),gprimd(3,3)
 real(dp) :: qbz_cart(3),rmet(3,3)
 real(dp),allocatable :: cov(:,:),par(:),qfit(:,:),sigma(:),var(:)
 real(dp),allocatable :: vcfit(:,:),vcoul(:,:),vcoul_lwl(:,:),xx(:),yy(:)

! *************************************************************************

 DBG_ENTER("COLL")
 !
 ! === Test if the first q-point is zero === 
 ! FIXME this wont work if nqptdm/=0
 !if (normv(Qmesh%ibz(:,1),gmet,'G')<GW_TOLQ0)) STOP 'vcoul_init, non zero first point '
 rank   = xcomm_rank(comm)
 nprocs = xcomm_size(comm)
 master = 0

 call vcoul_nullify(Vcp)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 gvec => Gsph%gvec(:,:) 
 !
 ! === Save dimension and other useful quantities in Vcp% ===
 Vcp%nfft      = PRODUCT(ngfft(1:3))  ! Number of points in the FFT mesh.
 Vcp%ng        = ng                   ! Number of G-vectors in the Coulomb matrix elements.
 Vcp%nqibz     = Qmesh%nibz           ! Number of irred q-point.
 Vcp%nqlwl     = nqlwl                ! Number of small q-directions to deal with singularity and non Analytic behavior.
 Vcp%rcut      = rcut                 ! Cutoff radius for cylinder.
 Vcp%hcyl      = zero                 ! Length of finite cylinder (Rozzi"s method, default is Beigi).
 Vcp%ucvol     = ucvol                ! Unit cell volume.

 Vcp%rprimd    = rprimd(:,:)          ! Dimensional direct lattice.
 Vcp%boxcenter = zero                 ! boxcenter at the moment is supposed to be at the origin.
 Vcp%vcutgeo   = vcutgeo(:)           ! Info on the orientation and extension of the cutoff region.
 Vcp%ngfft(:)  = ngfft(:)             ! Info on the FFT mesh.
 !
 ! === Define geometry and cutoff radius (if used) ===
 Vcp%mode='NONE'
 if (icutcoul==0) Vcp%mode='SPHERE'
 if (icutcoul==1) Vcp%mode='CYLINDER'
 if (icutcoul==2) Vcp%mode='SURFACE'
 if (icutcoul==3) Vcp%mode='CRYSTAL'
 if (icutcoul==4) Vcp%mode='ERF'
 if (icutcoul==5) Vcp%mode='ERFC'
 if (icutcoul==6) Vcp%mode='AUXILIARY_FUNCTION'
 !
 ! === Calculate Fourier coefficient of Coulomb interaction ===
 ! * For the limit q-->0 we consider ng vectors due to a possible anisotropy in case of a cutoff interaction
 ABI_ALLOCATE(Vcp%qibz,(3,Vcp%nqibz))
 Vcp%qibz = Qmesh%ibz(:,:)
 ABI_ALLOCATE(Vcp%qlwl,(3,Vcp%nqlwl))
 Vcp%qlwl = qlwl(:,:)

 ! ===============================================
 ! == Calculation of the FT of the Coulomb term ==
 ! ===============================================
 ABI_ALLOCATE(vcoul    ,(ng,Vcp%nqibz))
 ABI_ALLOCATE(vcoul_lwl,(ng,Vcp%nqlwl))

 a1=rprimd(:,1); b1=two_pi*gprimd(:,1)
 a2=rprimd(:,2); b2=two_pi*gprimd(:,2)
 a3=rprimd(:,3); b3=two_pi*gprimd(:,3)

 SELECT CASE (TRIM(Vcp%mode)) 

 CASE ('SPHERE')
   ! TODO check that L-d > R_c > d
   ! * A non-positive value of rcut imposes the recipe of Spencer & Alavi, PRB 77, 193110 (2008).
   if (Vcp%rcut<tol12) then
     Vcp%rcut = (ucvol*Kmesh%nbz*3.d0/four_pi)**third
     write(msg,'(2a,2x,f8.4,a)')ch10,&
&     ' Using a calculated value for rcut = ',Vcp%rcut, ' to have the same volume as the BvK crystal '
     call wrtout(std_out,msg,'COLL') 
   end if

   Vcp%vcutgeo=zero
   call cutoff_sphere(Qmesh%nibz,Qmesh%ibz,ng,gvec,gmet,Vcp%rcut,vcoul)

   ! q-points for optical limit.
   call cutoff_sphere(Vcp%nqlwl,Vcp%qlwl,ng,gvec,gmet,Vcp%rcut,vcoul_lwl)
   !
   ! === Treat the limit q--> 0 ===
   ! * The small cube is approximated by a sphere, while vc(q=0)=2piR**2.
   ! * if a single q-point is used, the expression for the volume is exact.
   Vcp%i_sz=two_pi*Vcp%rcut**2
   call vcoul_print(Vcp,unit=ab_out) 

 CASE ('CYLINDER')
   test=COUNT(ABS(Vcp%vcutgeo)>tol6) 
   ABI_CHECK(test==1,'Wrong cutgeo for cylinder')

   ! === Beigi method is the default one, i.e infinite cylinder of radius rcut ===
   ! * Negative values to use Rozzi method with finite cylinder of extent hcyl.
   opt_cylinder=1; Vcp%hcyl=zero; Vcp%pdir(:)=0
   do ii=1,3
     check=Vcp%vcutgeo(ii)
     if (ABS(check)>tol6) then 
       Vcp%pdir(ii)=1
       if (check<zero) then  ! use Rozzi's method.
         Vcp%hcyl=ABS(check)*SQRT(SUM(rprimd(:,ii)**2))
         opt_cylinder=2 
       end if
     end if
   end do

   test=COUNT(Vcp%pdir==1) 
   ABI_CHECK((test==1),'Wrong pdir for cylinder')
   if (Vcp%pdir(3)/=1) then 
     MSG_ERROR("The cylinder must be along the z-axis")
   end if

   call cutoff_cylinder(Qmesh%nibz,Qmesh%ibz,ng,gvec,Vcp%rcut,Vcp%hcyl,Vcp%pdir,&
&    Vcp%boxcenter,rprimd,vcoul,opt_cylinder,comm)

   ! q-points for optical limit.
   call cutoff_cylinder(Vcp%nqlwl,Vcp%qlwl,ng,gvec,Vcp%rcut,Vcp%hcyl,Vcp%pdir,&
&    Vcp%boxcenter,rprimd,vcoul_lwl,opt_cylinder,comm) 
   !
   ! === If Beigi, treat the limit q--> 0 ===
   if (opt_cylinder==1) then 
     npar=8; npt=100 ; gamma_pt=RESHAPE((/0,0,0/),(/3,1/))
     ABI_ALLOCATE(qfit,(3,npt))
     ABI_ALLOCATE(vcfit,(1,npt))
     if (Qmesh%nibz==1) STOP
     qfit(:,:)=zero 
     step=half/(npt*(Qmesh%nibz-1))              ; qfit(3,:)=arth(tol6,step,npt)
     !step=(half/(Qmesh%nibz-1)/tol6)**(one/npt) ; qfit(3,:)=geop(tol6,step,npt)
     call cutoff_cylinder(npt,qfit,1,gamma_pt,Vcp%rcut,Vcp%hcyl,Vcp%pdir,Vcp%boxcenter,rprimd,vcfit,opt_cylinder,comm)

     ABI_ALLOCATE(xx,(npt))
     ABI_ALLOCATE(yy,(npt))
     ABI_ALLOCATE(sigma,(npt))
     ABI_ALLOCATE(par,(npar))
     ABI_ALLOCATE(var,(npar))
     ABI_ALLOCATE(cov,(npar,npar))
     do ii=1,npt 
      xx(ii)=normv(qfit(:,ii),gmet,'G') 
     end do
     sigma=one ; yy(:)=vcfit(1,:)
     !call llsfit_svd(xx,yy,sigma,npar,K0fit,chisq,par,var,cov,info)
     !do ii=1,npt
     ! write(99,*)xx(ii),yy(ii),DOT_PRODUCT(par,K0fit(xx(ii),npar))
     !end do
     bz_plane=l2norm(b1.x.b2)
     !integ=K0fit_int(xx(npt),par,npar)
     !write(std_out,*)' SVD fit : chi-square',chisq 
     !write(std_out,*)' fit-parameters : ',par
     !write(std_out,*)' variance ',var
     !write(std_out,*)' bz_plane ',bz_plane
     !write(std_out,*)' SCD integ ',integ
     ! Here Im assuming homogeneous mesh
     dx=(xx(2)-xx(1))
     integ=yy(2)*dx*3.0/2.0
     do ii=3,npt-2
       integ=integ+yy(ii)*dx
     end do
     integ=integ+yy(npt-1)*dx*3.0/2.0
     write(std_out,*)' simple integral',integ
     q0_volsph=(two_pi)**3/(Kmesh%nbz*ucvol) 
     q0_vol=bz_plane*two*xx(npt)
     write(std_out,*)' q0 sphere : ',q0_volsph,' q0_vol cyl ',q0_vol
     Vcp%i_sz=bz_plane*two*integ/q0_vol
     write(std_out,*)' spherical approximation ',four_pi*7.44*q0_volsph**(-two_thirds) 
     write(std_out,*)' Cylindrical cutoff value ',Vcp%i_sz
     !Vcp%i_sz=four_pi*7.44*q0_vol**(-two_thirds) 
     ABI_DEALLOCATE(xx)
     ABI_DEALLOCATE(yy)
     ABI_DEALLOCATE(sigma)
     ABI_DEALLOCATE(par)
     ABI_DEALLOCATE(var)
     ABI_DEALLOCATE(cov)
   else 
     ! In Rozzi"s method the lim q+G --> 0 is finite.
     Vcp%i_sz=vcoul(1,1)
   end if

   call vcoul_print(Vcp,unit=ab_out )

 CASE ('SURFACE') 
  
   test=COUNT(Vcp%vcutgeo/=zero)
   ABI_CHECK(test==2,"Wrong vcutgeo")
   !
   ! === Default is Beigi"s method ===
   opt_surface=1; Vcp%alpha(:)=zero 
   if (ANY(Vcp%vcutgeo<zero)) opt_surface=2
   Vcp%pdir(:)=zero
   do ii=1,3
     check=Vcp%vcutgeo(ii)
     if (ABS(check)>zero) then ! Use Rozzi"s method with a finite surface along x-y
       Vcp%pdir(ii)=1
       if (check<zero) Vcp%alpha(ii)=normv(check*rprimd(:,ii),rmet,'R')
     end if
   end do

   ! Beigi"s method: the surface must be along x-y and R must be L_Z/2.
   if (opt_surface==1) then
     ABI_CHECK(ALL(Vcp%pdir == (/1,1,0/)),"Surface must be in the x-y plane")
     Vcp%rcut = half*SQRT(DOT_PRODUCT(a3,a3))
   end if

   call cutoff_surface(Qmesh%nibz,Qmesh%ibz,ng,gvec,gprimd,Vcp%rcut,&
&    Vcp%boxcenter,Vcp%pdir,Vcp%alpha,vcoul,opt_surface)

   ! q-points for optical limit.
   call cutoff_surface(Vcp%nqlwl,Vcp%qlwl,ng,gvec,gprimd,Vcp%rcut,&
&    Vcp%boxcenter,Vcp%pdir,Vcp%alpha,vcoul_lwl,opt_surface)
   !
   ! === Treat the limit q--> 0 TODO 

 CASE ('AUXILIARY_FUNCTION')
   !
   ! Numerical integration of the exact-exchange divergence through the
   ! auxiliary function of Carrier et al. PRB 75, 205126 (2007).
   !
   do iq_ibz=1,Vcp%nqibz
     call cmod_qpg(Vcp%nqibz,iq_ibz,Vcp%qibz,ng,gvec,gprimd,vcoul(:,iq_ibz))

     if (iq_ibz==1) then ! The singularity is treated using vcoul_lwl.
       vcoul(1, iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/vcoul(2:,iq_ibz)**2
     else
       vcoul(:,iq_ibz) = four_pi/vcoul(:,iq_ibz)**2
     end if
   end do

   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/vcoul_lwl**2
   !
   ! === Integration of 1/q^2 singularity ===
   ! * We use the auxiliary function from PRB 75, 205126 (2007)
   q0_vol=(two_pi)**3/(Kmesh%nbz*ucvol) ; bz_geometry_factor=zero

   ! It might be useful to perform the integration using the routine quadrature
   ! so that we can control the accuracy of the integral and improve the
   ! numerical stability of the GW results. 
   do iq_bz=1,Qmesh%nbz
     qbz_cart(:)=Qmesh%bz(1,iq_bz)*b1(:)+Qmesh%bz(2,iq_bz)*b2(:)+Qmesh%bz(3,iq_bz)*b3(:)
     qbz_norm=SQRT(SUM(qbz_cart(:)**2))
     if (qbz_norm>tolq0) bz_geometry_factor=bz_geometry_factor-faux(Qmesh%bz(:,iq_bz))
   end do

   bz_geometry_factor = bz_geometry_factor + integratefaux()*Qmesh%nbz

   write(msg,'(2a,2x,f8.4)')ch10,&
&    ' integrate q->0 : numerical BZ geometry factor = ',bz_geometry_factor*q0_vol**(2./3.)
   call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL')

   Vcp%i_sz=four_pi*bz_geometry_factor  ! Final result stored here

 CASE ('CRYSTAL')

   do iq_ibz=1,Vcp%nqibz
     call cmod_qpg(Vcp%nqibz,iq_ibz,Vcp%qibz,ng,gvec,gprimd,vcoul(:,iq_ibz))

     if (iq_ibz==1) then ! The singularity is treated using vcoul_lwl.
       vcoul(1, iq_ibz) = zero 
       vcoul(2:,iq_ibz) = four_pi/vcoul(2:,iq_ibz)**2
     else 
       vcoul(:,iq_ibz) = four_pi/vcoul(:,iq_ibz)**2
     end if
   end do

   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/vcoul_lwl**2
   !
   ! === Integration of 1/q^2 singularity ===
   ! * We use the auxiliary function from PRB 75, 205126 (2007)
   q0_vol=(two_pi)**3/(Kmesh%nbz*ucvol) ; bz_geometry_factor=zero

   ! $$ MG: this is to restore the previous implementation, it will facilitate the merge
   ! Analytic integration of 4pi/q^2 over the volume element:
   ! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
   ! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
   ! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255
   ! (see gwa.pdf, appendix A.4)
   Vcp%i_sz=four_pi*7.44*q0_vol**(-two_thirds)

 CASE ('ERF') 
   ! Modified long-range only Coulomb interaction thanks to the error function:
   ! * Vc = erf(r/rcut)/r
   ! * The singularity is treated using vcoul_lwl.
   do iq_ibz=1,Qmesh%nibz
     call cmod_qpg(Qmesh%nibz,iq_ibz,Qmesh%ibz,ng,gvec,gprimd,vcoul(:,iq_ibz))

     ! * The Fourier transform of the error function reads
     if (iq_ibz==1) then 
       vcoul(1, iq_ibz) = zero
       ! FB: Please be careful with the parenthesis in the following. I think I fixed the bug
       ! already twice.
       vcoul(2:,iq_ibz) = four_pi/(vcoul(2:,iq_ibz)**2) *  EXP( -0.25d0 * (Vcp%rcut*vcoul(2:,iq_ibz))**2 )
     else 
       vcoul(:,iq_ibz)  = four_pi/(vcoul(:, iq_ibz)**2) *  EXP( -0.25d0 * (Vcp%rcut*vcoul(: ,iq_ibz))**2 )
     end if
   end do

   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/(vcoul_lwl**2) *  EXP( -0.25d0 * (Vcp%rcut*vcoul_lwl)**2 )

   !
   ! === Treat 1/q^2 singularity === 
   ! * We use the auxiliary function from PRB 75, 205126 (2007)
   q0_vol=(two_pi)**3/(Kmesh%nbz*ucvol) ; bz_geometry_factor=zero
   do iq_bz=1,Qmesh%nbz
     qbz_cart(:)=Qmesh%bz(1,iq_bz)*b1(:)+Qmesh%bz(2,iq_bz)*b2(:)+Qmesh%bz(3,iq_bz)*b3(:)
     qbz_norm=SQRT(SUM(qbz_cart(:)**2))
     if (qbz_norm>tolq0) bz_geometry_factor=bz_geometry_factor-faux(Qmesh%bz(:,iq_bz))
   end do

   bz_geometry_factor = bz_geometry_factor + integratefaux()*Qmesh%nbz

   write(msg,'(2a,2x,f12.4)')ch10,&
&    ' integrate q->0 : numerical BZ geometry factor = ',bz_geometry_factor*q0_vol**(2./3.)
   call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL')
   Vcp%i_sz=four_pi*bz_geometry_factor  ! Final result stored here

 CASE ('ERFC') 
   ! * Use a modified short-range only Coulomb interaction thanks to the complementar error function:
   !   $ V_c = [1-erf(r/r_{cut})]/r $
   ! * The Fourier transform of the error function reads
   !   vcoul=four_pi/(vcoul**2) * ( 1.d0 - exp( -0.25d0 * (Vcp%rcut*vcoul)**2 ) )
   do iq_ibz=1,Qmesh%nibz
     call cmod_qpg(Qmesh%nibz,iq_ibz,Qmesh%ibz,ng,gvec,gprimd,vcoul(:,iq_ibz))

     if (iq_ibz==1) then 
       vcoul(1 ,iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/(vcoul(2:,iq_ibz)**2) * ( one - EXP( -0.25d0 * (Vcp%rcut*vcoul(2:,iq_ibz))**2 ) )
     else 
       vcoul(:, iq_ibz) = four_pi/(vcoul(:, iq_ibz)**2) * ( one - EXP( -0.25d0 * (Vcp%rcut*vcoul(:, iq_ibz))**2 ) )
     end if
   end do
   !
   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/(vcoul_lwl**2) * ( one - EXP( -0.25d0 * (Vcp%rcut*vcoul_lwl)**2 ) )
   !
   ! === Treat 1/q^2 singularity === 
   ! * There is NO singularity in this case.
   Vcp%i_sz=pi*Vcp%rcut**2 ! Final result stored here

 CASE DEFAULT
   write(msg,'(2a)')' Unknown cutoff mode ',TRIM(Vcp%mode)
   MSG_BUG(msg)
 END SELECT
 !
 ! === Store final results in complex array ===
 ! * Rozzi"s cutoff can give real negative values 

 ABI_ALLOCATE(Vcp%vc_sqrt,(ng,Vcp%nqibz))
 Vcp%vc_sqrt=CMPLX(vcoul,zero) 
 Vcp%vc_sqrt=SQRT(Vcp%vc_sqrt) 
 call vcoul_plot(Vcp,Qmesh,Gsph,ng,vcoul,comm)
 ABI_DEALLOCATE(vcoul)

 ABI_ALLOCATE(Vcp%vcqlwl_sqrt,(ng,Vcp%nqlwl))
 istat = ABI_ALLOC_STAT
 Vcp%vcqlwl_sqrt=CMPLX(vcoul_lwl,zero) 
 Vcp%vcqlwl_sqrt=SQRT(Vcp%vcqlwl_sqrt) 
 ABI_DEALLOCATE(vcoul_lwl)
 !
 ! === Create a table for FFT points inside the cutoff region ===
 if (TRIM(Vcp%mode)/='CRYSTAL') then
   ABI_ALLOCATE(Vcp%ctab,(Vcp%nfft))
   call cutoff_table(Vcp)
 else 
   ABI_ALLOCATE(Vcp%ctab,(1))
 end if

 !Vcp%i_sz=zero

 call vcoul_print(Vcp,unit=std_out)

 DBG_EXIT("COLL")

contains !===============================================================

 real(dp) function integratefaux()

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'integratefaux'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !MG here there is a bug in abilint since integratefaux is chaged to integrate, dont know why!
 !real(dp) :: integrate

!Local variables-------------------------------
  integer,parameter :: nref=3  
  integer,parameter :: nqx=50
  integer :: ierr,iqx1,iqy1,iqz1,iqx2,iqy2,iqz2
  real(dp) :: qq1(3),qq2(3)

  ! nqx is the number of sampling points along each axis for the numerical integration
  ! nref is the area where the mesh is refined

  integratefaux=zero
  do iqx1=1,nqx
   do iqy1=1,nqx
    do iqz1=1,nqx
     if(modulo(iqx1,nprocs)/=rank) cycle
     qq1(1)=DBLE(iqx1)/DBLE(nqx)-half
     qq1(2)=DBLE(iqy1)/DBLE(nqx)-half
     qq1(3)=DBLE(iqz1)/DBLE(nqx)-half
     !
     ! Refine the mesh for the point close to the origin
     if( abs(iqx1-nqx/2)<=nref .and. abs(iqy1-nqx/2)<=nref .and. abs(iqz1-nqx/2)<=nref ) then
      do iqx2=1,nqx
       do iqy2=1,nqx
        do iqz2=1,nqx
         qq2(1)=qq1(1) + (DBLE(iqx2)/DBLE(nqx)-half ) /DBLE(nqx)
         qq2(2)=qq1(2) + (DBLE(iqy2)/DBLE(nqx)-half ) /DBLE(nqx)
         qq2(3)=qq1(3) + (DBLE(iqz2)/DBLE(nqx)-half ) /DBLE(nqx)
         !
         ! Treat the remaining divergence in the origin as if it would be a spherical
         ! integration of 1/q^2
         if( iqx1/=nqx/2 .or. iqy1/=nqx/2 .or. iqz1/=nqx/2 .or. iqx2/=nqx/2 .or. iqy2/=nqx/2 .or. iqz2/=nqx/2 ) then
          integratefaux=integratefaux+ faux(qq2) /DBLE(nqx)**6
         else
          integratefaux=integratefaux + 7.7955* ( (two_pi)**3/ucvol/DBLE(nqx)**6 )**(-2./3.) /DBLE(nqx)**6
         end if
        end do
       end do
      end do
     else
      integratefaux=integratefaux+faux(qq1)/DBLE(nqx)**3
     end if
    end do
   end do
  end do

  call xsum_mpi(integratefaux,comm,ierr)

 end function integratefaux

 function faux(qq)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'faux'
!End of the abilint section

  real(dp),intent(in) :: qq(3)
  real(dp) :: faux

  faux= four*( dot_product(b1,b1) * SIN(two_pi*qq(1)/two)**2 &
&             +dot_product(b2,b2) * SIN(two_pi*qq(2)/two)**2 &
&             +dot_product(b3,b3) * SIN(two_pi*qq(3)/two)**2 &
&            )                                              &
&       +two*( dot_product(b1,b2) * SIN(two_pi*qq(1))*SIN(two_pi*qq(2)) &
&             +dot_product(b2,b3) * SIN(two_pi*qq(2))*SIN(two_pi*qq(3)) &
&             +dot_product(b3,b1) * SIN(two_pi*qq(3))*SIN(two_pi*qq(1)) &
&            )
  faux=(two_pi)**2/faux * exp( -0.25d0*Vcp%rcut**2* sum( ( qq(1)*b1(:)+qq(2)*b2(:)+qq(3)*b3(:) )**2 ) )
 end function faux

end subroutine vcoul_init
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcoul_plot
!! NAME
!! vcoul_plot
!!
!! FUNCTION
!! Plot vccut(q,G) as a function of |q+G|. Calculate also vc in real space.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_vcoul
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_plot(Vcp,Qmesh,Gsph,ng,vc,comm)

 use defs_basis
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vcoul_plot'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,comm
 type(Bz_mesh_type),intent(in) :: Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(Gvectors_type),intent(in) :: Gsph
!arrays
 real(dp),intent(in) :: vc(ng,Qmesh%nibz)

!Local variables-------------------------------
!scalars
 integer :: icount,idx_Sm1G,ierr,ig,igs,ii,iq_bz,iq_ibz,iqg,ir,isym,itim
 integer :: master,my_start,my_stop,nqbz,nqibz,nr,ntasks,rank,unt
 real(dp) :: arg,fact,l1,l2,l3,lmax,step,tmp,vcft,vc_bare
 character(len=500) :: msg
 character(len=fnlen) :: filnam
!arrays
 integer,allocatable :: insort(:)
 integer,pointer :: gvec(:,:)
 real(dp) :: b1(3),b2(3),b3(3),gmet(3,3),gprimd(3,3),qbz(3),qpgc(3)
 real(dp),allocatable :: qpg_mod(:),rr(:,:,:),vcr(:,:),vcr_cut(:,:)

!************************************************************************

 if (TRIM(Vcp%mode)/='CYLINDER')  RETURN

 rank = xcomm_rank(comm)
 master=0

 nqibz=Qmesh%nibz 
 nqbz =Qmesh%nbz

 gmet   =  Gsph%gmet(:,:) 
 gprimd =  Gsph%gprimd(:,:) 
 gvec   => Gsph%gvec(:,:)

 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)

 ! === Compare in Fourier space the true Coulomb with the cutted one ===
 if (rank==master) then 
   ABI_ALLOCATE(insort,(nqibz*ng))
   ABI_ALLOCATE(qpg_mod,(nqibz*ng))
   iqg=1
   do iq_ibz=1,nqibz
     do ig=1,ng
       qpg_mod(iqg)=normv(Qmesh%ibz(:,iq_ibz)+gvec(:,ig),gmet,'g')
       insort(iqg)=iqg; iqg=iqg+1
     end do
   end do
   call sort_dp(nqibz*ng,qpg_mod,insort,tol14)

   filnam='_VCoulFT_'; call isfile(filnam,'new')
   unt=get_unit()
   open(unt,file=filnam,status='new',form='formatted')
   write(unt,'(a,i3,a,i6,a)')&
&    '#   |q+G|       q-point (Tot no.',nqibz,')        Gvec (',ng,')     vc_bare(q,G)    vc_cutoff(q,G) '

   do iqg=1,nqibz*ng
     iq_ibz=(insort(iqg)-1)/ng +1
     ig=(insort(iqg))-(iq_ibz-1)*ng
     vc_bare=zero
     if (qpg_mod(iqg)>tol16) vc_bare=four_pi/qpg_mod(iqg)**2
     write(unt,'(f12.6,2x,3f8.4,2x,3i6,2x,2es14.6)')&
&      qpg_mod(iqg),Qmesh%ibz(:,iq_ibz),gvec(:,ig),vc_bare,vc(ig,iq_ibz)
   end do 

   close(unt)
   ABI_DEALLOCATE(insort)
   ABI_DEALLOCATE(qpg_mod)
 end if ! rank==master

 ! === Fourier transform back to real space just to check cutoff implementation ===
 ntasks=nqbz*ng
 call xmpi_split_work(ntasks,comm,my_start,my_stop,msg,ierr)
 if (ierr/=0) then
   MSG_WARNING(msg)
 end if

 l1=SQRT(SUM(Vcp%rprimd(:,1)**2))
 l2=SQRT(SUM(Vcp%rprimd(:,2)**2))
 l3=SQRT(SUM(Vcp%rprimd(:,3)**2))

 nr=50
 lmax=MAX(l1,l2,l3) ; step=lmax/(nr-1)
 fact=one/(Vcp%ucvol*nqbz)
 !
 ! numb coding 
 ABI_ALLOCATE(rr,(3,nr,3))
 rr=zero
 do ii=1,3
   do ir=1,nr
     rr(ii,ir,ii)=(ir-1)*step
   end do
 end do

 ABI_ALLOCATE(vcr,(nr,3))
 ABI_ALLOCATE(vcr_cut,(nr,3))
 vcr=zero; vcr_cut=zero 

 do iq_bz=1,nqbz
   call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym,itim)
   if (ABS(qbz(1))<0.01) qbz(1)=zero
   if (ABS(qbz(2))<0.01) qbz(2)=zero
   if (ABS(qbz(3))<0.01) qbz(3)=zero
   igs=1 ; if (ALL(qbz(:)==zero)) igs=2
   do ig=igs,ng
     icount=ig+(iq_bz-1)*ng 
     if (icount<my_start.or.icount>my_stop) CYCLE
     idx_Sm1G=Gsph%rottbm1(ig,itim,isym) ! IS{^-1}G
     vcft=vc(idx_Sm1G,iq_ibz)
     qpgc(:)=qbz(:)+gvec(:,ig) ; qpgc(:)=b1(:)*qpgc(1)+b2(:)*qpgc(2)+b3(:)*qpgc(3)
     tmp=SQRT(DOT_PRODUCT(qpgc,qpgc)) ; tmp=tmp**2
     do ii=1,3
       do ir=1,nr
         arg=DOT_PRODUCT(rr(:,ir,ii),qpgc)
         vcr_cut(ir,ii)=vcr_cut(ir,ii) + vcft*COS(arg)
         vcr    (ir,ii)=vcr    (ir,ii) + four_pi/tmp*COS(arg)
       end do
     end do
   end do !ig
 end do !iq_ibz

 call xsum_master(vcr_cut,master,comm,ierr)
 call xsum_master(vcr    ,master,comm,ierr)
 call xbarrier_mpi(comm)

 if (rank==master) then 
   filnam='_VCoulR_' ; call isfile(filnam,'new')
   unt=get_unit()
   open(unt,file=filnam,status='new',form='formatted')
   write(unt,'(a)')'# length  vc_bare(r)   vc_cut(r) '
   do ir=1,nr
     write(unt,'(7es18.6)')(ir-1)*step,(fact*vcr(ir,ii),fact*vcr_cut(ir,ii),ii=1,3)
   end do
   close(unt)
 end if
 ABI_DEALLOCATE(rr)
 ABI_DEALLOCATE(vcr)

end subroutine vcoul_plot
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcoul_print
!! NAME
!! vcoul_print
!!
!! FUNCTION
!!  Print the content of a Coulomb datatype.
!!
!! INPUTS
!!  Vcp<vcoul_t>=The datatype whose content has to be printed.
!!  [unit]=The unit number for output 
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS".
!!
!! PARENTS
!!      m_vcoul
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_print(Vcp,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vcoul_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(vcoul_t),intent(in) :: Vcp

!Local variables-------------------------------
!scalars
 integer :: ii,my_unt,my_prtvol,iqlwl
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt    =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode   ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
 my_prtvol=0       ; if (PRESENT(prtvol    )) my_prtvol=prtvol

 SELECT CASE (Vcp%mode) 

 CASE ('SPHERE')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
&    ' === Spherical cutoff === ',ch10,ch10,&
&    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10,&
&    '  Volume of the sphere .. ',four_pi/three*Vcp%rcut**3,' [Bohr^3] '
!FB: This has no meaning here! &  '  Sphere centered at .... ',Vcp%boxcenter,' (r.l.u) ',ch10
!MG It might be useful if the system is not centered on the origin because in this case the 
!   matrix elements of the Coulomb have to be multiplied by a phase depending on boxcenter.
!   I still have to decide if it is useful to code this possibility and which variable use to 
!   define the center (boxcenter is used in the tddft part). 
   call wrtout(my_unt,msg,my_mode) 

 CASE ('CYLINDER')
   ii=imin_loc(ABS(Vcp%pdir-1))
   write(msg,'(5a,f10.4,3a,i2,2a,3f10.2,a)')ch10,&
&    ' === Cylindrical cutoff === ',ch10,ch10,&
&    '  Cutoff radius ............... ',Vcp%rcut,' [Bohr] ',ch10,&
&    '  Axis parallel to direction... ',ii,ch10,&
&    '  Passing through point ....... ',Vcp%boxcenter,' (r.l.u) '
   call wrtout(my_unt,msg,my_mode) 

   write(msg,'(2a)')'  Infinite length  ....... ',ch10
   if (Vcp%hcyl/=zero) write(msg,'(a,f8.5,2a)')'  Finite length of ....... ',Vcp%hcyl,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode) 

 CASE ('SURFACE') 
   write(msg,'(5a,f10.4,3a,3f10.2,2a)')ch10,&
&    ' === Surface cutoff === ',ch10,ch10,&
&    '  Cutoff radius .................... ',Vcp%rcut,' [Bohr] ',ch10,&
&    '  Central plane passing through .... ',Vcp%boxcenter,' (r.l.u) ',ch10
   call wrtout(my_unt,msg,my_mode) 
   !write(msg,'(a)')'  Infinite length  .......'
   !if (Vcp%hcyl/=zero) write(msg,'(a,f8.5,a)')'  Finite length of .......',Vcp%hcyl,' [Bohr] '
   !call wrtout(my_unt,msg,my_mode) 

 CASE ('AUXILIARY_FUNCTION')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',TRIM(Vcp%mode)
   call wrtout(my_unt,msg,my_mode) 

 CASE ('CRYSTAL')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',TRIM(Vcp%mode)
   call wrtout(my_unt,msg,my_mode) 

 CASE ('ERF')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
&    ' === Error function cutoff === ',ch10,ch10,&
&    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 CASE ('ERFC')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
&    ' === Complement Error function cutoff === ',ch10,ch10,&
&    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 CASE DEFAULT
   write(msg,'(2a)')' Unknown cutoff mode: ',TRIM(Vcp%mode)
   MSG_BUG(msg)
 END SELECT

 if (Vcp%nqlwl>0) then 
   write(msg,'(a,i3)')" q-points for optical limit: ",Vcp%nqlwl
   call wrtout(my_unt,msg,my_mode)
   do iqlwl=1,Vcp%nqlwl
     write(msg,'(1x,i5,a,2x,3f12.6)') iqlwl,')',Vcp%qlwl(:,iqlwl)
     call wrtout(my_unt,msg,my_mode)
   end do
 end if

 !TODO add additional information

end subroutine vcoul_print 
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/cutoff_table
!! NAME
!! cutoff_table
!!
!! FUNCTION
!!  Create a table for each point in the FFT grid: 0 if the point is outside the cutoff region, 1 otherwise 
!!  This table can be used to redefine the density in the isolated system (see cutoff_density)
!!
!! INPUTS
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  Vcp<type vcoul_t> 
!!   The following quantities are used :
!!   %ngfft(1:3)= FFT dimensions
!!   %nr= number of point in the FFT mesh
!!   %boxcenter=reduced coordinates of the center of the box  (input variable
!!  (FIXME kept from tddft.F90, should write something in doc) 
!!   %cutoffmode= string of characters defining the cutoff mode 
!!   %rcut= cutoff radius 
!!
!! OUTPUT
!!  Vcp%ctab(Vcp%nr)
!!  For each point of the FFT grid gives 
!!   1 if the point in inside the cutoff region
!!   0 otherwise
!!
!! PARENTS
!!      m_vcoul
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine cutoff_table(Vcp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cutoff_table'
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(vcoul_t),intent(inout) :: Vcp

!Local variables-------------------------------
!scalars
 integer :: ir,ix,iy,iz,ngfft1,ngfft2,ngfft3
 real(dp) :: ucvol
 character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rr(3),rxy(3),rz(3)

!************************************************************************

 if (TRIM(Vcp%mode)/='CRYSTAL') RETURN

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)
 ngfft1=Vcp%ngfft(1) 
 ngfft2=Vcp%ngfft(2) 
 ngfft3=Vcp%ngfft(3)

 Vcp%ctab(:)=1
 do ix=0,ngfft1-1
   rr(1)=DBLE(ix)/ngfft1
   do iy=0,ngfft2-1
     rr(2)=DBLE(iy)/ngfft2
     do iz=0,ngfft3-1
       rr(3)=DBLE(iz)/ngfft3
       ir=1+ix+iy*ngfft1+iz*ngfft1*ngfft2
       rr(:)=rr(:)-Vcp%boxcenter(:)

       SELECT CASE (TRIM(Vcp%mode))
       CASE ('SPHERE')
         if (normv(rr,rmet,'R')>Vcp%rcut) Vcp%ctab(ir)=0
       CASE ('CYLINDER') ! Infinite cylinder
         rxy(1)=rr(1)
         rxy(2)=rr(2)
         rxy(3)=zero
         if (normv(rxy,rmet,'R')>Vcp%rcut) Vcp%ctab(ir)=0
       CASE ('SURFACE')
         rz(:)=zero
         rz(3)=rr(3)
         if (normv(rz,rmet,'R')>Vcp%rcut) Vcp%ctab(ir)=0
       CASE DEFAULT   
         write(msg,'(2a)')' Undefined cutoff mode ',TRIM(Vcp%mode)
         MSG_BUG(msg)
       END SELECT

     end do 
   end do 
 end do

end subroutine cutoff_table
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/cutoff_density
!! NAME
!! cutoff_density
!!
!! FUNCTION
!!  Modify density in case of calculations with Coulomb cutoff
!!
!! INPUTS
!!  ngfft(18)=information on the FFT grid 
!!  (at the moment only the first three elements are used and are checked
!!  against the values contained in the vcoul_t datatype)
!!  nspden=input variable 
!!  Vcp <type vcoul_t>, only ctab(ngfft(1)*ngfft(2)*ngfft(3)) is used: 
!!   For each point of the FFT grid gives 
!!   1 if the point in inside the cutoff region, 0 otherwise
!!  ucvol=unit cell volume 
!!
!! OUTPUT
!!  Only printing
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine cutoff_density(ngfft,nspden,nsppol,Vcp,rhor)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cutoff_density'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspden,nsppol
 type(vcoul_t),intent(in) :: Vcp
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(ngfft(1)*ngfft(2)*ngfft(3),nspden)

!Local variables-------------------------------
!scalars
 integer :: ir,ix,iy,iz,ngfft1,ngfft2,ngfft3
 real(dp),parameter :: EPS=1.d-2
 logical :: ltest
!arrays
 real(dp) :: nel(nspden),overflow(nspden),rr(3)

! *************************************************************************
 
 ! === Do some check ===
 ltest=(ALL(ngfft(1:3)==Vcp%ngfft(1:3))) 
 ABI_CHECK(ltest,'Inconsistent ngfft')

 if (Vcp%mode=='CRYSTAL') RETURN

 ngfft1=ngfft(1) 
 ngfft2=ngfft(2) 
 ngfft3=ngfft(3)
 overflow(:)=zero ; nel(:)=zero

 do ix=0,ngfft1-1
  rr(1)=DBLE(ix)/ngfft1
  do iy=0,ngfft2-1
   rr(2)=DBLE(iy)/ngfft2
   do iz=0,ngfft3-1
    rr(3)=DBLE(iz)/ngfft3
    ir=1+ix+iy*ngfft1+iz*ngfft1*ngfft2

    nel(:)=nel(:)+rhor(ir,:)
    if (Vcp%ctab(ir)==0) then 
     overflow(:)=overflow+rhor(ir,:)
     !rhor(ir,:)=zero
    end if 

   end do 
  end do 
 end do

 nel(:)=nel(:)*Vcp%ucvol/(ngfft1*ngfft2*ngfft3)
 overflow(:)=overflow(:)*Vcp%ucvol/(ngfft1*ngfft2*ngfft3)

 write(std_out,*)' Number of electrons per unit cell       = ',nel(1:nsppol)
 write(std_out,*)' Charge density outside the cutoff region= ',overflow(1:nsppol)
 !
 ! If overflow is larger that few percent of the total charge warn or stop
 ! since one should increase the size of the supercell or the cutoff radius
 ! if (ANY(overflow(:)>EPS*nel(:)))  then 
 !  write(msg,'(4a,f8.5,3a)')ch10,&
 !&  ' cutoff_density : WARNING - ',ch10,&
 !&  '  More than ',eps,' % of the charge density is outside the cutoff region',ch10,&
 !&  '  You should try to increase the cutoff radius '
 !   call wrtout(std_out,msg,'COLL') 
 ! end if 

 !write(std_out,*)' restoring removed charge '
 ! Restore charge neutrality, spreading the overflow inside the cutoff region
 ! Using this approach if the charge is zero close to 
 ! the boundary then there should be no spurious effect
 !do is=1,nspden
 ! fact=nel(is)/(nel(is)-overflow(is))
 ! rhor(:,is)=rhor(:,is)*fact
 !end do

end subroutine cutoff_density
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcoul_nullify
!! NAME
!! vcoul_nullify
!!
!! FUNCTION
!!  Nullify a vcoul_t type 
!!
!! SIDE EFFECTS
!!  Vcp<vcoul_t>=the datatype with all pointers nullified
!!
!! PARENTS
!!      m_shexc,m_vcoul
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_nullify(Vcp) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vcoul_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(vcoul_t),intent(inout) :: Vcp

! *************************************************************************

 nullify(Vcp%ctab       )
 nullify(Vcp%qibz       )
 nullify(Vcp%qlwl       )
 nullify(Vcp%vc_sqrt    )
 nullify(Vcp%vcqlwl_sqrt)
 
end subroutine vcoul_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcoul_free
!! NAME
!! vcoul_free
!!
!! FUNCTION
!!  Destroy a vcoul_t type 
!!
!! SIDE EFFECTS
!!  Vcp<vcoul_t>=the datatype to be destroyed
!!
!! PARENTS
!!      bethe_salpeter,m_shexc,mrgscr,screening,sigma
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_free(Vcp) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vcoul_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(vcoul_t),intent(inout) :: Vcp

!Local variables ------------------------------
!scalars
 integer :: istat

! *************************************************************************

 DBG_ENTER("COLL")

 if (associated(Vcp%ctab       ))   then
   ABI_DEALLOCATE(Vcp%ctab)
   istat = ABI_ALLOC_STAT
 end if
 if (associated(Vcp%qibz       ))   then
   ABI_DEALLOCATE(Vcp%qibz)
   istat = ABI_ALLOC_STAT
 end if
 if (associated(Vcp%qlwl       ))   then
   ABI_DEALLOCATE(Vcp%qlwl)
   istat = ABI_ALLOC_STAT
 end if
 if (associated(Vcp%vc_sqrt    ))   then
   ABI_DEALLOCATE(Vcp%vc_sqrt)
   istat = ABI_ALLOC_STAT
 end if
 if (associated(Vcp%vcqlwl_sqrt))   then
   ABI_DEALLOCATE(Vcp%vcqlwl_sqrt)
   istat = ABI_ALLOC_STAT
 end if

 DBG_EXIT("COLL")
 
end subroutine vcoul_free
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/cutoff_sphere
!! NAME
!! cutoff_sphere
!!
!! FUNCTION
!!  Calculate the Fourier transform of the Coulomb interaction with a spherical cutoff:  
!!   $ v_{cut}(G)= \frac{4\pi}{|q+G|^2} [ 1-cos(|q+G|*R_cut) ] $  (1)
!!
!! INPUTS
!!  gmet(3,3)=Metric in reciprocal space.
!!  gvec(3,ngvec)=G vectors in reduced coordinates.
!!  rcut=Cutoff radius of the sphere.
!!  ngvec=Number of G vectors
!!  nqpt=Number of q-points 
!!  qpt(3,nqpt)=q-points where the cutoff Coulomb is required.
!!
!! OUTPUT
!!  vc_cut(ngvec,nqpt)=Fourier components of the effective Coulomb interaction.
!!
!! NOTES
!!  For |q|<small and G=0 we use 2pi.R_cut^2, namely we consider the limit q-->0 of Eq. (1) 
!!
!! PARENTS
!!      m_vcoul
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine cutoff_sphere(nqpt,qpt,ngvec,gvec,gmet,rcut,vc_cut)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cutoff_sphere'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngvec,nqpt
 real(dp),intent(in) :: rcut
!arrays
 integer,intent(in) :: gvec(3,ngvec)
 real(dp),intent(in) :: gmet(3,3),qpt(3,nqpt)
 real(dp),intent(out) :: vc_cut(ngvec,nqpt)

!Local variables-------------------------------
!scalars
 integer :: ig,igs,iqpt
 real(dp) :: qpg
 logical :: ltest

!************************************************************************

 ltest=ALL(gvec(:,1)==0)
 ABI_CHECK(ltest,'The first G vector should be Gamma')

 do iqpt=1,nqpt
   igs=1
   if (normv(qpt(:,iqpt),gmet,'G')<GW_TOLQ0) then ! For small q and G=0, use the limit q-->0.
     vc_cut(1,iqpt)=two_pi*rcut**2
     igs=2
   end if
   do ig=igs,ngvec
     qpg=normv(qpt(:,iqpt)+gvec(:,ig),gmet,'G')
     vc_cut(ig,iqpt)=four_pi*(one-COS(rcut*qpg))/qpg**2
   end do
 end do
 
end subroutine cutoff_sphere
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/cutoff_cylinder
!! NAME
!! cutoff_cylinder
!!
!! FUNCTION
!!  Calculate the Fourier components of an effective Coulomb interaction 
!!  zeroed outside a finite cylindrical region.
!!  Two methods are implemented: 
!!   method==1: The interaction in the (say) x-y plane is truncated outside the Wigner-Seitz 
!!              cell centered on the wire in the x-y plane. The interaction has infinite 
!!              extent along the z axis and the Fourier transform is singular only at the Gamma point. 
!!              Only orthorombic Bravais lattice are supported.
!!   method==2: The interaction is truncated outside a cylinder of radius rcut. The cylinder has finite 
!!              extent along z. No singularity occurs. 
!!
!! INPUTS
!!  boxcenter(3)= center of the wire in the x-y axis
!!  gvec(3,ng)=G vectors in reduced coordinates 
!!  ng=number of G vectors
!!  qpt(3,nq)= q-points 
!!  nq=number of q-points 
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!   where: rprimd(i,j)=rprim(i,j)*acell(j)
!!  method=1 for Beigi approach (infinite cylinder with interaction truncated outside the W-S cell)
!!         2 for Rozzi method (finite cylinder)
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  vccut(ng,nq)= Fourier components of the effective Coulomb interaction 
!!
!! PARENTS
!!      m_vcoul
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine cutoff_cylinder(nq,qpt,ng,gvec,rcut,hcyl,pdir,boxcenter,rprimd,vccut,method,comm)

 use defs_basis
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cutoff_cylinder'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,nq,method,comm
 real(dp),intent(in) :: rcut,hcyl 
!arrays
 integer,intent(in) :: gvec(3,ng),pdir(3)
 real(dp),intent(in) :: boxcenter(3),qpt(3,nq),rprimd(3,3)
 real(dp),intent(out) :: vccut(ng,nq)

!Local variables-------------------------------
!scalars
 integer,parameter :: N0=1000
 integer :: icount,ig,igs,iq
 integer :: ntasks,ierr,my_start,my_stop
 real(dp) :: j0,j1,k0,k1,qpg2,qpg_xy,tmp
 real(dp) :: qpg_z,quad,rcut2,hcyl2,c1,c2,ucvol,SMALL
 logical :: q_is_zero
 character(len=500) :: msg
!arrays
 real(dp) :: qpg(3),b1(3),b2(3),b3(3),gmet(3,3),rmet(3,3),gprimd(3,3),qc(3),gcart(3)
 real(dp),allocatable :: qcart(:,:)
!************************************************************************

 ABI_UNUSED(pdir)
 ABI_UNUSED(boxcenter)
 !
 ! ===================================================
 ! === Setup for the quadrature of matrix elements ===
 ! ===================================================
 qopt_    =6         ! Quadrature method, see quadrature routine.
 ntrial_  =30        ! Max number of attempts.
 accuracy_=0.001     ! Fractional accuracy required.
 npts_    =6         ! Initial number of point (only for Gauss-Legendre method).
 SMALL    =GW_TOLQ0  ! Below this value (q+G)_i is treated as zero.
 rcut_    =rcut      ! Radial cutoff, used only if method==2
 hcyl_    =hcyl      ! Lenght of cylinder along z, only if method==2

 write(msg,'(3a,2(a,i5,a),a,f8.5)')ch10,&
&  ' cutoff_cylinder: Info on the quadrature method : ',ch10,&
&  '  Quadrature scheme      = ',qopt_,ch10,&
&  '  Max number of attempts = ',ntrial_,ch10,&
&  '  Fractional accuracy    = ',accuracy_
 call wrtout(std_out,msg,'COLL') 
 !
 ! === From reduced to Cartesian coordinates ===
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)

 ABI_ALLOCATE(qcart,(3,nq))
 do iq=1,nq
   qcart(:,iq)=b1(:)*qpt(1,iq)+b2(:)*qpt(2,iq)+b3(:)*qpt(3,iq)
 end do

 ntasks=nq*ng 
 call xmpi_split_work(ntasks,comm,my_start,my_stop,msg,ierr)
 if (ierr/=0) then
   MSG_WARNING(msg)
 end if
 !
 ! ================================================
 ! === Different approaches according to method ===
 ! ================================================
 vccut(:,:)=zero 

 SELECT CASE (method)

 CASE (1)
   ! === Infinite cylinder, interaction is zeroed outside the Wigner-Seitz cell ===
   ! * Beigi"s expression holds only if the BZ is sampled only along z.
  write(msg,'(2(a,f8.4))')' cutoff_cylinder: Using Beigi''s Infinite cylinder '
  call wrtout(std_out,msg,'COLL') 
  if (ANY(qcart(1:2,:)>SMALL)) then 
    write(std_out,*)' qcart = ',qcart(:,:)
    write(msg,'(5a)')&
&     ' found q-points with non zero components in the X-Y plane. ',ch10,&
&     ' This is not allowed, see Notes in cutoff_cylinder.F90. ',ch10,&
&     ' ACTION: Modify the q-point sampling. '
    MSG_ERROR(msg)
  end if 
  ! * Check if Bravais lattice is orthorombic and parallel to the Cartesian versors.
  !   In this case the intersection of the W-S cell with the x-y plane is a rectangle with -ha_<=x<=ha_ and -hb_<=y<=hb_ 
  if ( (ANY(ABS(rprimd(2:3,  1))>tol6)).or.&
&      (ANY(ABS(rprimd(1:3:2,2))>tol6)).or.&
&      (ANY(ABS(rprimd(1:2,  3))>tol6))    &
&    ) then
    msg = ' Bravais lattice should be orthorombic and parallel to the cartesian versors '
    MSG_ERROR(msg)
  end if

  ha_=half*SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1)))
  hb_=half*SQRT(DOT_PRODUCT(rprimd(:,2),rprimd(:,2)))
  r0_=MIN(ha_,hb_)/N0
  ! 
  ! === For each pair (q,G) evaluate the integral defining the cutoff Coulomb  ===
  ! * Here the code assumes that all q-vectors are non zero and q_xy/=0.
  do iq=1,nq
   igs=1 
   ! * Skip singularity at Gamma, it will be treated "by hand" in csigme.
   q_is_zero = (normv(qpt(:,iq),gmet,'G')<GW_TOLQ0)
   !if (q_is_zero) igs=2
   qc(:)=qcart(:,iq)
   write(msg,'(2(a,i3))')' entering loop iq: ',iq,' with igs = ',igs
   call wrtout(std_out,msg,'COLL')

   do ig=igs,ng
    icount=ig+(iq-1)*ng; if (icount<my_start.or.icount>my_stop) CYCLE

    gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
    qpg(:)=qc(:)+gcart(:)
    qpgx_=qpg(1) ; qpgy_=qpg(2) ; qpg_para_=ABS(qpg(3))
    write(std_out,*)"qpgx_=",qpgx_ 
    write(std_out,*)"qpgy_=",qpgy_ 
    write(std_out,*)"qpg_para=",qpg_para_ 

    ! Avoid singularity in K_0{qpg_para_\rho) by using a small q along the periodic dimension.
    if (q_is_zero.and.qpg_para_<tol6) then 
      qpg_para_ = tol6
      write(std_out,*)"setting qpg_para to=",qpg_para_ 
    end if
    !
    ! * Calculate $ 2\int_{WS} dxdy K_0{qpg_para_\rho) cos(x.qpg_x + y.qpg_y) $
    !   where WS is the Wigner-Seitz cell.
    tmp=zero
    ! Difficult part, integrate on a small cirle of radius r0 using spherical coordinates
    !call quadrature(K0cos_dth_r0,zero,r0_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
    !if (ierr/=0) MSG_ERROR("Accuracy not reached")
    !write(std_out,'(i8,a,es14.6)')ig,' 1 ',quad
    !tmp=tmp+quad
    ! Add region with 0<=x<=r0 and y>=+-(SQRT(r0^2-x^2))since WS is rectangular
    !call quadrature(K0cos_dy_r0,zero,r0_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
    !if (ierr/=0) MSG_ERROR("Accuracy not reached")
    !write(std_out,'(i8,a,es14.6)')ig,' 2 ',quad
    !tmp=tmp+quad
    ! Get the in integral in the rectangle with x>=r0, should be the easiest but sometimes has problems to converge
    !call quadrature(K0cos_dy,r0_,ha_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
    !if (ierr/=0) MSG_ERROR("Accuracy not reached")
    !write(std_out,'(i8,a,es14.6)')ig,' 3 ',quad
    ! 
    ! === More stable method: midpoint integration with Romberg extrapolation ===
    call quadrature(K0cos_dy,zero,ha_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
    !write(std_out,'(i8,a,es14.6)')ig,' 3 ',quad
    if (ierr/=0) then 
     MSG_ERROR("Accuracy not reached")
    end if
    ! === Store final result ===
    ! * Factor two comes from the replacement WS -> (1,4) quadrant thanks to symmetries of the integrad.
    tmp=tmp+quad
    vccut(ig,iq)=two*(tmp*two) 

   end do !ig
  end do !iq

 CASE (2) 
   ! === Finite cylinder of lenght hcyl, from Rozzi et al ===
   ! TODO add check on hcyl value that should be smaller that 1/deltaq 
   if (hcyl_<zero) then
     write(msg,'(a,f8.4)')' Negative value for cylinder lenght hcyl_=',hcyl_
     MSG_BUG(msg)
   end if

   if (ABS(hcyl_)>tol12) then 
     write(msg,'(2(a,f8.4))')' cutoff_cylinder: using finite cylinder of length= ',hcyl_,' rcut= ',rcut_
     call wrtout(std_out,msg,'COLL') 
     hcyl2=hcyl_**2 
     rcut2=rcut_**2

     do iq=1,nq
       write(msg,'(a,i3)')' entering loop iq: ',iq
       call wrtout(std_out,msg,'COLL')
       qc(:)=qcart(:,iq)
       do ig=1,ng
         ! === No singularity occurs in finite cylinder, thus start from 1 ===
         icount=ig+(iq-1)*ng ; if (icount<my_start.or.icount>my_stop) CYCLE
         gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
         qpg(:)=qc(:)+gcart(:)
         qpg_para_=ABS(qpg(3)) ; qpg_perp_=SQRT(qpg(1)**2+qpg(2)**2)

         if (qpg_perp_/=zero.and.qpg_para_/=zero) then 
           ! $ 4\pi\int_0^{R_c} d\rho\rho j_o(qpg_perp_.\rho)\int_0^hcyl dz\cos(qpg_para_*z)/sqrt(\rho^2+z^2) $ 
           call quadrature(F2,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then 
             MSG_ERROR("Accuracy not reached")
           end if

           vccut(ig,iq)=four_pi*quad

         else if (qpg_perp_==zero.and.qpg_para_/=zero) then 
           ! $ \int_0^h sin(qpg_para_.z)/\sqrt(rcut^2+z^2)dz $
           call quadrature(F3,zero,hcyl_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then 
             MSG_ERROR("Accuracy not reached")
           end if

           c1=one/qpg_para_**2-COS(qpg_para_*hcyl_)/qpg_para_**2-hcyl_*SIN(qpg_para_*hcyl_)/qpg_para_
           c2=SIN(qpg_para_*hcyl_)*SQRT(hcyl2+rcut2)
           vccut(ig,iq)=four_pi*c1+four_pi*(c2-quad)/qpg_para_

         else if (qpg_perp_/=zero.and.qpg_para_==zero) then 
           ! $ 4pi\int_0^rcut d\rho \rho J_o(qpg_perp_.\rho) ln((h+\sqrt(h^2+\rho^2))/\rho) $
           call quadrature(F4,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then 
             MSG_ERROR("Accuracy not reached")
           end if

           vccut(ig,iq)=four_pi*quad

         else if (qpg_perp_==zero.and.qpg_para_==zero) then 
           ! Use lim q+G --> 0
           vccut(ig,iq)=two_pi*(-hcyl2+hcyl_*SQRT(hcyl2+rcut2)+rcut2*LOG((hcyl_+SQRT(hcyl_+SQRT(hcyl2+rcut2)))/rcut_))

         else 
           MSG_BUG('You should not be here!')
         end if

       end do !ig
     end do !iq

   else 
     ! === Infinite cylinder ===
     write(msg,'(a)')' cutoff_cylinder : using Rozzi''s method with infinite cylinder '
     call wrtout(std_out,msg,'COLL') 
     do iq=1,nq
       write(msg,'(a,i3)')' entering loop iq: ',iq
       call wrtout(std_out,msg,'COLL')
       qc(:)=qcart(:,iq)
       do ig=1,ng
         icount=ig+(iq-1)*ng ; if (icount<my_start.or.icount>my_stop) CYCLE
         gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
         qpg(:)=qc(:)+gcart(:)
         qpg2  =DOT_PRODUCT(qpg,qpg)
         qpg_z =ABS(qpg(3)) ; qpg_xy=SQRT(qpg(1)**2+qpg(2)**2)
         if (qpg_z>SMALL) then 
           ! === Analytic expression ===
           call CALCK0(qpg_z *rcut_,k0,1)
           call CALJY1(qpg_xy*rcut_,j1,0)
           call CALJY0(qpg_xy*rcut_,j0,0)
           call CALCK1(qpg_z *rcut_,k1,1)
           vccut(iq,ig)=(four_pi/qpg2)*(one+rcut_*qpg_xy*j1*k0-qpg_z*rcut_*j0*k1) 
         else 
           if (qpg_xy>SMALL) then 
             ! === Integrate r*Jo(G_xy r)log(r) from 0 up to rcut_  ===
             call quadrature(F5,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
             if (ierr/=0) then 
               MSG_ERROR("Accuracy not reached")
             end if
             vccut(ig,iq)=-four_pi*quad
           else 
             ! === Analytic expression ===
             vccut(ig,iq)=-pi*rcut_**2*(two*LOG(rcut_)-one)
           end if 
         end if 
       end do !ig
     end do !iq 
   end if !finite/infinite

 CASE DEFAULT
   write(msg,'(a,i3)')' Wrong value for method: ',method 
   MSG_BUG(msg)
 END SELECT 
 !
 ! === Collect vccut on each node ===
 call xbarrier_mpi(comm)
 call xsum_mpi(vccut,comm,ierr)

 ABI_DEALLOCATE(qcart)

end subroutine cutoff_cylinder
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/cutoff_surface
!! NAME
!! cutoff_surface
!!
!! FUNCTION
!!  Calculate the Fourier components of an effective Coulomb interaction 
!!  within a slab of thickness 2*rcut which is symmetric with respect to the xy plane. 
!!  In this implementation rcut=L_z/2 where L_z is the periodicity along z 
!!
!! INPUTS
!!  gprimd(3,3)=Dimensional primitive translations in reciprocal space ($\textrm{bohr}^{-1}$).
!!  gvec(3,ng)=G vectors in reduced coordinates. 
!!  ng=Number of G vectors.
!!  qpt(3,nq)=q-points 
!!  nq=Number of q-points 
!!  gmet(3,3)=Metric in reciprocal space.
!!
!! OUTPUT
!!  vc_cut(ng,nq)=Fourier components of the effective Coulomb interaction.
!!
!! NOTES
!!  The Fourier expression for an interaction truncated along the z-direction 
!!  (i.e non-zero only if |z|<R) is :
!!  
!!  vc(q.G) = 4pi/|q+G|^2 * [ 1 + e^{-((q+G)_xy)*R} * ( (q_z+G_z)/(q+G)_xy * sin((q_z+G_z)R) - 
!!   - cos((q_z+G_Z)R)) ]  (1)
!! 
!!  Equation (1) diverges when q_xy+G_xy --> 0 for any non zero q_z+G_z
!!  However if we choose R=L/2, where L defines the periodicity along z, 
!!  and we limit ourselves to consider q-points such as q_z==0, then 
!!  sin((q_z+G_z)R)=sin(G_Z 2pi/L)=0 for every G.  
!!  Under these assumptions we obtain
!!
!!  v(q,G) = 4pi/|q+G|^2 [1-e^{-(q+G)_xy*L/2}\cos((q_z+G_z)R)]
!! 
!!  which is always finite when G_z /=0 while it diverges as 4piR/(q+G)_xy as (q+G)_xy -->0 
!!  but only in the x-y plane.
!!
!! PARENTS
!!      m_vcoul
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine cutoff_surface(nq,qpt,ng,gvec,gprimd,rcut,boxcenter,pdir,alpha,vc_cut,method)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cutoff_surface'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,ng,nq
 real(dp),intent(in) :: rcut
!arrays
 integer,intent(in) :: gvec(3,ng),pdir(3)
 real(dp),intent(in) :: alpha(3),boxcenter(3),gprimd(3,3),qpt(3,nq)
 real(dp),intent(out) :: vc_cut(ng,nq)

!Local variables-------------------------------
!scalars
 integer :: ig,iq,igs
 real(dp),parameter :: SMALL=tol6
 real(dp) :: qpg2,qpg_para,qpg_perp
 character(len=500) :: msg
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gcart(3),qc(3),qpg(3)
 real(dp),allocatable :: qcart(:,:)

! *************************************************************************

 ABI_UNUSED(pdir) 
 ABI_UNUSED(boxcenter) 

 ! === From reduced to cartesian coordinates ===
 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)
 ABI_ALLOCATE(qcart,(3,nq))
 do iq=1,nq
   qcart(:,iq) = b1*qpt(1,iq) + b2*qpt(2,iq) + b3*qpt(3,iq)
 end do
 !
 ! === Different approaches according to method ===
 vc_cut(:,:)=zero 

 SELECT CASE (method)

 CASE (1) ! Beigi's expression.
   ! * q-points with non-zero component along the z-axis are not allowed if 
   !   the simplified Eq.1 for the Coulomb interaction is used.
   if (ANY(ABS(qcart(3,:))>SMALL)) then 
     write(std_out,*)qcart(:,:) 
     write(msg,'(5a)')&
&      ' Found q-points with non-zero component along non-periodic direction ',ch10,&
&      ' This is not allowed, see Notes in cutoff_surface.F90 ',ch10,&
&      ' ACTION : Modify the q-point sampling '
     MSG_ERROR(msg)
   end if 
   ! 
   ! === Calculate truncated Coulomb interaction for a infinite surface ===
   ! * Here I suppose that all the input q-points are different from zero
   do iq=1,nq
     qc(:)=qcart(:,iq)
     igs=1; if (SQRT(DOT_PRODUCT(qc,qc))<tol16) igs=2 ! avoid (q=0, G=0)
     do ig=igs,ng
       gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
       qpg(:)=qc(:)+gcart(:)
       qpg2  =DOT_PRODUCT(qpg(:),qpg(:))
       qpg_para=SQRT(qpg(1)**2+qpg(2)**2) ; qpg_perp=qpg(3)
       ! if (abs(qpg_perp)<SMALL.and.qpg_para<SMALL) stop 'SMALL in cutoff_surface
       vc_cut(ig,iq)=four_pi/qpg2*(one-EXP(-qpg_para*rcut)*COS(qpg_perp*rcut))  
     end do 
   end do

 CASE (2) ! Rozzi"s method
   MSG_ERROR("Work in progress")
   ABI_UNUSED(alpha) ! just to keep alpha as an argument
   !alpha=?? ; ap1sqrt=SQRT(one+alpha**2)
   do iq=1,nq
     qc(:)=qcart(:,iq)
     do ig=1,ng
       gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
       qpg(:)=qc(:)+gcart(:)
       qpg2  =DOT_PRODUCT(qpg(:),qpg(:))
       qpg_para=SQRT(qpg(1)**2+qpg(2)**2) ; qpg_perp =qpg(3)
       if (qpg_para>SMALL) then 
        vc_cut(ig,iq)=four_pi/qpg2*(one+EXP(-qpg_para*rcut)*(qpg_perp/qpg_para*SIN(qpg_perp*rcut)-COS(qpg_perp*rcut))) 
       else 
         if (ABS(qpg_perp)>SMALL) then 
           vc_cut(ig,iq)=four_pi/qpg_perp**2*(one-COS(qpg_perp*rcut)-qpg_perp*rcut*SIN(qpg_perp*rcut)) !&
!&          + 8*rcut*SIN(qpg_perp*rcut)/qpg_perp*LOG((alpha+ap1sqrt)*(one+ap1sqrt)/alpha) ! contribution due to finite surface
         else 
           vc_cut(ig,iq)=-two_pi*rcut**2
         end if 
       end if 
     end do !ig 
   end do !iq

 CASE DEFAULT 
   write(msg,'(a,i3)')' Wrong value of method: ',method 
   MSG_BUG(msg)
 END SELECT

 ABI_DEALLOCATE(qcart)

end subroutine cutoff_surface
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/cmod_qpg
!! NAME
!! cmod_qpg
!!
!! FUNCTION
!! Set up table of lengths |q+G| for the Coulomb potential.
!!
!! INPUTS
!! gvec(3,npwvec)=Reduced coordinates of the G vectors.
!! gprimd(3,3)=Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!! iq=Index specifying the q point.
!! npwvec=Number of planewaves
!! nq=Number of q points.
!! q=Coordinates of q points.
!!
!! OUTPUT
!! qplusg(npwvec)=Norm of q+G vector
!!
!! PARENTS
!!      calc_ffm,m_ppmodel,m_vcoul
!!
!! CHILDREN
!!      calck0,jbessel,quadrature
!!
!! SOURCE

subroutine cmod_qpg(nq,iq,q,npwvec,gvec,gprimd,qplusg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cmod_qpg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,npwvec,nq
!arrays
 integer,intent(in) :: gvec(3,npwvec)
 real(dp),intent(in) :: gprimd(3,3),q(3,nq)
 real(dp),intent(out) :: qplusg(npwvec)

!Local variables ------------------------------
!scalars
 integer :: ig,ii
!arrays
 real(dp) :: gmet(3,3),gpq(3)

!************************************************************************

 ! Compute reciprocal space metrics
 do ii=1,3
   gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&             gprimd(2,ii)*gprimd(2,:)+&
&             gprimd(3,ii)*gprimd(3,:)
 end do

 if (ALL(ABS(q(:,iq))<1.e-3)) then !FIXME avoid this, everything should be under the control of the programmer.
   ! * Treat q as it were zero except when G=0
   qplusg(1)=two_pi*SQRT(DOT_PRODUCT(q(:,iq),MATMUL(gmet,q(:,iq))))
   do ig=2,npwvec
     gpq(:)=gvec(:,ig)
     qplusg(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
   end do
 else
   do ig=1,npwvec
     gpq(:)=gvec(:,ig)+q(:,iq)
     qplusg(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
   end do
 end if

end subroutine cmod_qpg

!----------------------------------------------------------------------

function F2(xx)
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'F2'
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: F2
!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: intr
!************************************************************************

 zz_=xx
 call quadrature(F1,zero,rcut_,qopt_,intr,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then 
   MSG_ERROR("Accuracy not reached")
 end if

 F2=intr*COS(qpg_para_*xx)

end function F2
!!***

!----------------------------------------------------------------------

function F1(rho) 
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'F1'
 use interfaces_32_util
!End of the abilint section

 real(dp),intent(in) :: rho
 real(dp) :: F1

!Local variables-------------------------------
!scalars
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 !F1(\rho;z)= \rho*j_o(qpg_perp_*\rho)/sqrt(\rho**2+z**2)
 arg=rho*qpg_perp_
 call jbessel(bes,besp,bespp,ll,order,arg)

 if (zz_==zero) then 
   F1=bes
 else 
   F1=bes*rho/SQRT(rho**2+zz_**2)
 end if

end function F1
!!***

!----------------------------------------------------------------------

function F3(xx)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'F3'
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: F3
!************************************************************************

 ! F3(z)=z*\sin(qpg_para_*z)/\sqrt(rcut^2+z^2)
 F3=xx*SIN(qpg_para_*xx)/SQRT(rcut_**2+xx**2)

end function F3
!!***

!----------------------------------------------------------------------

function F4(rho)
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'F4'
 use interfaces_32_util
!End of the abilint section

 real(dp),intent(in) :: rho
 real(dp) :: F4

!Local variables-------------------------------
!scalars
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 ! $F4(rho)=\rho*j_o(qpg_perp_.\rho) \ln((hcyl+\sqrt(rho^2+hcyl^2))/\rho)$
 if (ABS(rho)<tol12) then 
   F4=zero
 else
   arg=rho*qpg_perp_
   call jbessel(bes,besp,bespp,ll,order,arg)
   F4=bes*rho*LOG((hcyl_+SQRT(rho**2+hcyl_**2))/rho)
 end if

end function F4
!!***

!----------------------------------------------------------------------

function F5(rho)
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'F5'
 use interfaces_32_util
!End of the abilint section

 real(dp),intent(in) :: rho
 real(dp) :: F5

!Local variables-------------------------------
!scalars
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 ! $F5(\rho)=\rho*j_o(G_perp\rho)log(\rho)$
 if (rho==0) then 
   F5=zero
 else 
   arg=rho*qpg_perp_
   call jbessel(bes,besp,bespp,ll,order,arg)
   F5=bes*rho*LOG(rho)
 end if

end function F5
!!***

!----------------------------------------------------------------------

function K0cos(yy) 
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'K0cos'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 real(dp),intent(in) :: yy
 real(dp) :: K0cos

!Local variables-------------------------------
!scalars
 real(dp) :: k0,rho,arg
!************************************************************************

 ! K0cos(y)=K0(\rho*|qpg_z|)*COS(x.qpg_x+y*qpg_y)
 rho=SQRT(xx_**2+yy**2) ; arg=qpg_para_*rho
 call CALCK0(arg,k0,1)
 K0cos=k0*COS(qpgx_*xx_+qpgy_*yy)

end function K0cos
!!***

!----------------------------------------------------------------------

function K0cos_dy(xx)
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'K0cos_dy'
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: K0cos_dy
!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: quad
!************************************************************************

 !! K0cos_dy(x)=\int_{-b/2}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
 xx_=xx 
 call quadrature(K0cos,-hb_,+hb_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then 
   MSG_ERROR("Accuracy not reached")
 end if

 K0cos_dy=quad

end function K0cos_dy
!!***

!----------------------------------------------------------------------

function K0cos_dy_r0(xx)
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'K0cos_dy_r0'
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: K0cos_dy_r0
!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: quad,yx
!************************************************************************

 ! $ K0cos_dy_r0(x)= \int_{-b/2}^{-y(x)} K0(|qpg_z|\rho) cos(x.qpg_x+y.qpg_y)dy 
 !                  +\int_{y(x)}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
 ! where y(x)=SQRT(r0^2-x^2) and x<=r0
 !
 xx_=xx; yx=SQRT(r0_**2-xx**2)
 call quadrature(K0cos,-hb_,-yx,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then 
   MSG_ERROR("Accuracy not reached")
 end if
 K0cos_dy_r0=quad

 call quadrature(K0cos,+yx,+hb_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then 
   MSG_ERROR("Accuracy not reached")
 end if

 K0cos_dy_r0=quad+K0cos_dy_r0

end function K0cos_dy_r0
!!***

!----------------------------------------------------------------------

function K0cos_dth_r0(rho)
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'K0cos_dth_r0'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 real(dp),intent(in) :: rho
 real(dp) :: K0cos_dth_r0

!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: quad,arg,k0,tmp
!************************************************************************

 ! $ K0cos_dth_r0(\rho)= 
 ! \int_{0}^{2pi)} K0(|qpg_z|\rho)cos(\rho.cos(\theta).qpg_x+\rho.sin(\theta).qpg_y) d\theta $
 !
 ! where y(x)=SQRT(r0^2-x^2) and x<=r0
 !
 rho_=rho 
 call quadrature(Fcos_th,zero,two_pi,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then 
   MSG_ERROR("Accuracy not reached")
 end if

 arg=qpg_para_*rho_ 
 tmp=zero 
 if (arg>tol6) then
   call CALCK0(arg,k0,1)
   tmp=k0*rho_
 end if
 K0cos_dth_r0=quad*tmp

end function K0cos_dth_r0
!!***

!----------------------------------------------------------------------

function Fcos_th(theta) 
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Fcos_th'
!End of the abilint section

 real(dp),intent(in) :: theta
 real(dp) :: Fcos_th

!Local variables-------------------------------
!scalars
!************************************************************************

 ! $ Fcos_th(\theta)=rho*K0(\rho*|qpg_z|)*COS(\rho.COS(\theta).qpg_x+\rho.SIN/(\theta)*qpg_y) $

 !arg=qpg_para_*rho_ 
 !call CALCK0(arg,k0,1)
 !tmp=k0*rho_
 Fcos_th=COS(rho_*COS(theta)*qpgx_+rho_*SIN(theta)*qpgy_)

end function Fcos_th
!!***

!----------------------------------------------------------------------

!the following functions should be used to deal with the singularity in the Cylindrical cutoff
!TODO Not yet used and indeed are still private 

function K0fit(mq,nn) result(vals)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'K0fit'
!End of the abilint section

 integer,intent(in) :: nn
 real(dp),intent(in) :: mq
 real(dp) :: vals(nn)
!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: mqh
!arrays
 real(dp),parameter :: cc(7)=(/-0.57721566,0.42278420,0.23069756,&
&                               0.03488590,0.00262698,0.00010750,0.00000740/) 
 ! *************************************************************************

  if (nn>8.or.nn<1) STOP 'not implemented'
  ! === Eq 9.8.5 in Abramovitz ===
  vals(1)=-LOG(mq*half)*I0(mq) 
  mqh=mq*half
  do ii=2,nn
   vals(ii)=cc(ii-1)*mqh**(2*(ii-2))
  end do
end function K0fit

function K0fit_int(mq,par,nn) result(integ)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'K0fit_int'
!End of the abilint section

 integer,intent(in) :: nn
 real(dp),intent(in) :: mq
 real(dp) :: integ
 real(dp),intent(in) :: par(nn)
!Local variables-------------------------------
!scalars
 integer :: ii,aa
 real(dp) :: mqh
!arrays
 real(dp),parameter :: cc(7)=(/-0.57721566,0.42278420,0.23069756,&
&                               0.03488590,0.00262698,0.00010750,0.00000740/) 
 ! *************************************************************************

 if (nn>8.or.nn<1) STOP 'not implemented'

 mqh=mq*half
 integ=-par(1)*int_I0ln(mqh)
 ! primitive of polynomial \sum_0^{N/2} cc_{2i} (x/2)^{2*i}
 do ii=2,nn
  aa=(2*(ii-1)+1)
  integ=integ+par(ii)*two*cc(ii-1)*(mqh**aa)/aa
 end do

end function K0fit_int

function I0(xx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'I0'
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: I0
!Local variables-------------------------------
 real(dp) :: tt
  ! *************************************************************************

 ! Eq 9.8.1 of Abramovitz, entering the expansion of K0 -->0
 ! Expansion holds for |x|<3.75, Error<1.6*10D-07 
 tt=xx/3.75
 I0=one+3.5156229*tt**2+3.0899424*tt**4 +1.2067492*tt**6 &
       +0.2659732*tt**8+0.0360768*tt**10+0.0045813*tt**12
end function I0
 
! Primitive of x^m Ln(x) for m/=-1
function int_xmln(xx,mm)  result(res)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int_xmln'
!End of the abilint section

 integer,intent(in) :: mm
 real(dp),intent(in) :: xx
 real(dp) :: res
! *********************************************************************

 if (mm==-1) STOP ' invalid value for mm '
 if (xx<=zero) STOP ' invalid value '

 res= (xx**(mm+1))/(mm+1) * (LOG(xx) - one/(mm+1))

end function int_xmln

! Primitive function of ln(x/2)*I0(x) = sum_0^{N/2} 2^{2s+1} c_{2s} T(x/2,2s) 
! where T(x,s)=\int x^s ln(x)dx 
function int_I0ln(xx) result(res)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int_I0ln'
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: res
!Local variables-------------------------------
 real(dp) :: yy
! *********************************************************************
 
 yy=xx*half
 res =  (       one*2    *int_xmln(yy,0)  &
&        +3.5156229*2**3 *int_xmln(yy,2)  &  
&        +3.0899424*2**5 *int_xmln(yy,4)  & 
&        +1.2067492*2**7 *int_xmln(yy,6)  &
&        +0.2659732*2**9 *int_xmln(yy,8)  & 
&        +0.0360768*2**11*int_xmln(yy,10) &
&        +0.0045813*2**13*int_xmln(yy,12) &
&       )

end function int_I0ln

END MODULE m_vcoul
!!***
