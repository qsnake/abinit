!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_electronpositron
!! NAME
!!  m_electronpositron
!!
!! FUNCTION
!!  This module provides the definition of the electronpositron_type used
!!  used to store data for the electron-positron two-component DFT
!!  as methods to operate on it.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_electronpositron

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_energies

 implicit none

 private

! public constants
 integer,public,parameter :: EP_NOTHING  =-1
 integer,public,parameter :: EP_ELECTRON = 0
 integer,public,parameter :: EP_POSITRON = 1

! public procedures
 public :: init_electronpositron
 public :: destroy_electronpositron
 public :: exchange_electronpositron
 public :: electronpositron_calctype
!!***

!!****t* m_electronpositron/electronpositron_type
!! NAME
!!
!! FUNCTION
!!
!! NOTES
!!
!! SOURCE

 type, public :: electronpositron_type

! Integer scalars
  integer :: calctype        ! type of electron-positron calculation:
                             !   0: no calculation
                             !   1: positron in the electrons potential
                             !   2: electrons in the positron potential
  integer :: particle        ! current particle stored in electronpositron%xxx_ep arrays
                             !                 -1: no particle, 0: electron, 1: positron
  integer :: dimcg           ! Dimension of cg array dimcg=dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
  integer :: dimcprj         ! Dimension of cprj array dimcprj=dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*usecprj
  integer :: dimeigen        ! Dimension of eigen array dimeigen=dtset%mband*dtset%nkpt*dtset%nsppol
  integer :: dimocc          ! Dimension of occ array dimocc=dtset%mband*dtset%nkpt*dtset%nsppol
  integer :: has_pawrhoij_ep ! flag for pawrhoij_ep (0: not allocated, 1: allocated, 2: computed)
  integer :: has_pos_ham     ! flag: 1 if current Hamiltonian in memory (vtrial, vpsp, vhartr, vxc, paw_ij%dij)
!                                    is the positronic hamiltonian, 0 is it is the electronic one
  integer :: ixcpositron     ! XC type for electron-positron correlation
  integer :: istep           ! Current index of TC-DFT SCF step
  integer :: istep_scf       ! Current index of DFT SCF step  in current electron/positron minimization
  integer :: lmmax           ! Max. number of (l,m) moments over all types of atom
  integer :: natom           ! Number of atoms
  integer :: nfft            ! Number of points in FFT grid
  integer :: nspden          ! Number of spin density components
  integer :: nstep           ! Max. number of steps for the TC-DFT SCF cycle

! Logical scalars
  logical :: posdensity0_limit ! True if we are in the zero positron density limit
  logical :: scf_converged     ! True if the SCF cycle is converged for a positronic/electronic GS calculation

! Real(dp) scalars
  real(dp) :: e_hartree      !  Hartree electron-positron interaction energy
  real(dp) :: e_xc           !  XC electron-positron interaction energy
  real(dp) :: e_xcdc         !  Double-counting XC electron-positron interaction energy
  real(dp) :: e_paw          !  PAW electron-positron interaction energy
  real(dp) :: e_pawdc        !  Double-counting PAW electron-positron interaction energy
  real(dp) :: e0             !  Energy only due to particle(s) currently evolving
                                  !   calctype=1, energy due to positron  only
                                  !   calctype=2, energy due to electrons only
  real(dp) :: etotal_prev    !  Total energy of the previous GS calculation
  real(dp) :: lambda         ! Electron-positron annihilation rate
  real(dp) :: lifetime       ! Positron lifetime
  real(dp) :: maxfor_prev    ! Max. force of the previous GS calculation
  real(dp) :: posocc         ! Occupation number for the positron
  real(dp) :: postoldfe      ! Tolerance on total energy for the TC-DFT SCF cycle
  real(dp) :: postoldff      ! Tolerance on max. force for the TC-DFT SCF cycle

! Other scalars
  type(energies_type) :: energies_ep  !  Energies of the previous electronic/positronic SCF step

! Logical pointers
  logical, pointer :: lmselect_ep(:,:)
!  lmselect_ep(lmmax,natom)
!  flags selecting the non-zero LM-moments of on-site densities

! Real(dp) pointers
  real(dp), pointer :: cg_ep(:,:)
!  cg_ep(2,dimcg)
!  if typecalc=1: electronic wavefunctions
!  if typecalc=2: positronic wavefunctions

  real(dp), pointer :: eigen_ep(:)
!  eigen(dimeigen)
!  if typecalc=1: electronic eigen energies
!  if typecalc=2: positronic eigen energies

  real(dp), pointer :: fred_ep(:,:)
!  fred_ep(3,natom)
!  if typecalc=1: forces only due to electrons
!  if typecalc=2: forces only due to positron

  real(dp), pointer :: nhat_ep(:,:)
!  nhat_ep(nfft,nspden)
!  if typecalc=1: electronic compensation charge density in real space
!  if typecalc=2: positronic compensation charge density in real space

  real(dp), pointer :: occ_ep(:)
!  occ(dimocc)
!  if typecalc=1: electronic occupations
!  if typecalc=2: positronic occupations

  real(dp), pointer :: rhor_ep(:,:)
!  rhor_ep(nfft,nspden)
!  if typecalc=1: electronic density in real space
!  if typecalc=2: positronic density in real space

  real(dp), pointer :: stress_ep(:)
!  stress_ep(6)
!  if typecalc=1: stresses only due to electrons
!  if typecalc=2: stresses only due to positron

  real(dp), pointer :: vha_ep(:)
!  vha_ep(nfft)
!  if typecalc=1: electronic Hartree potential
!  if typecalc=2: positronic Hartree potential

! Other pointers
  type(pawrhoij_type), pointer :: pawrhoij_ep(:)
!  pawrhoij_ep(natom)
!  Relevant only if PAW
!  if typecalc=1: electronic PAW occupation matrix associated with rhor_ep
!  if typecalc=2: positronic PAW occupation matrix associated with rhor_ep

  type(cprj_type), pointer :: cprj_ep(:,:)
!  cprj_ep(natom,dimcprj)
!  Relevant only if PAW
!  if typecalc=1: electronic WF projected on nl projectors <p_i|Cnk>
!  if typecalc=2: positronic WF projected on nl projectors <p_i|Cnk>

 end type electronpositron_type


CONTAINS

!===========================================================
!!***

!!****f* m_electronpositron/init_electronpositron
!! NAME
!!  init_electronpositron
!!
!! FUNCTION
!!  Init all scalars and pointers in the structure.
!!
!! INPUTS
!!  ireadwf=if 1, read the wavefunction
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  mpi_enreg=informations about MPI parallelization:
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!!
!! SIDE EFFECTS
!!  electronpositron=<type(electronpositron_type)>=electronpositron datastructure
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,energies_copy,fourdp,rhoij_alloc
!!      rhoij_copy,rhoij_free
!!
!! SOURCE

subroutine init_electronpositron(ireadwf,dtset,electronpositron,mpi_enreg,nfft,pawrhoij,pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_electronpositron'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ireadwf,nfft
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 type(pawrhoij_type), intent(in) :: pawrhoij(dtset%natom*dtset%usepaw)
 type(pawtab_type),intent(in)  :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ii,my_nspinor,ncpgr,optfor,optstr
!arrays
 integer,allocatable :: nlmn(:)

!************************************************************************

 !@electronpositron_type

 if (dtset%positron/=0) then

  ABI_ALLOCATE(electronpositron,)

  electronpositron%calctype=0
  electronpositron%particle=-1

  electronpositron%ixcpositron=dtset%ixcpositron
  electronpositron%natom=dtset%natom
  electronpositron%nfft=nfft
  electronpositron%nspden=dtset%nspden
  electronpositron%istep=0
  electronpositron%istep_scf=0

  electronpositron%posocc=dtset%posocc
  electronpositron%nstep=dtset%posnstep
  electronpositron%postoldfe=dtset%postoldfe
  electronpositron%postoldff=dtset%postoldff
  electronpositron%posdensity0_limit=(dtset%ixcpositron/=2)
  electronpositron%scf_converged=.false.
  electronpositron%has_pos_ham=0

  call energies_init(electronpositron%energies_ep)

  electronpositron%e_hartree  =zero
  electronpositron%e_xc       =zero
  electronpositron%e_xcdc     =zero
  electronpositron%e_paw      =zero
  electronpositron%e_pawdc    =zero
  electronpositron%e0         =zero
  electronpositron%etotal_prev=zero
  electronpositron%maxfor_prev=zero

  electronpositron%lambda=zero
  electronpositron%lifetime=zero

  ABI_ALLOCATE(electronpositron%rhor_ep,(nfft,dtset%nspden))
  ABI_ALLOCATE(electronpositron%vha_ep,(nfft))
  ABI_ALLOCATE(electronpositron%pawrhoij_ep,(dtset%natom*dtset%usepaw))

  if (dtset%usepaw==1) then
   electronpositron%has_pawrhoij_ep=1
   ABI_ALLOCATE(nlmn,(dtset%ntypat))
   do ii=1,dtset%ntypat;nlmn(ii)=pawtab(ii)%lmn_size;end do
   call rhoij_alloc(pawrhoij(1)%cplex,nlmn,pawrhoij(1)%nspden,pawrhoij(1)%nspinor,&
&                   pawrhoij(1)%nsppol,electronpositron%pawrhoij_ep,dtset%typat,&
&                   mpi_enreg=mpi_enreg,ngrhoij=pawrhoij(1)%ngrhoij,nlmnmix=pawrhoij(1)%lmnmix_sz,&
&                   use_rhoij_=pawrhoij(1)%use_rhoij_,use_rhoijres=pawrhoij(1)%use_rhoijres)
   ABI_DEALLOCATE(nlmn)
   electronpositron%lmmax=0
   do ii=1,dtset%ntypat
    electronpositron%lmmax=max(electronpositron%lmmax,pawtab(ii)%lcut_size**2)
   end do
   ABI_ALLOCATE(electronpositron%lmselect_ep,(electronpositron%lmmax,dtset%natom))
   if (maxval(pawtab(1:dtset%ntypat)%usexcnhat)==0) then
     ABI_ALLOCATE(electronpositron%nhat_ep,(nfft,dtset%nspden))
   else
     nullify(electronpositron%nhat_ep)
   end if
  else
   electronpositron%has_pawrhoij_ep=0
   electronpositron%lmmax=0
   nullify(electronpositron%lmselect_ep)
   nullify(electronpositron%nhat_ep)
  end if

  if (dtset%positron<=-10) then
   my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
   electronpositron%dimcg=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
   electronpositron%dimocc=dtset%mband*dtset%nkpt*dtset%nsppol
   electronpositron%dimeigen=dtset%mband*dtset%nkpt*dtset%nsppol
   ABI_ALLOCATE(electronpositron%cg_ep,(2,electronpositron%dimcg))
   ABI_ALLOCATE(electronpositron%eigen_ep,(electronpositron%dimeigen))
   ABI_ALLOCATE(electronpositron%occ_ep,(electronpositron%dimocc))
   electronpositron%dimcprj=0
   if (.false.) then !TEMPORARY: will be activated later
!   if (dtset%usepaw==1.and.dtset%pawusecp>0) then
    electronpositron%dimcprj=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
    ABI_ALLOCATE(electronpositron%cprj_ep,(dtset%natom,electronpositron%dimcprj))
    ABI_ALLOCATE(nlmn,(dtset%natom))
    ncpgr=0
    do ii=1,dtset%natom;nlmn(ii)=pawtab(dtset%typat(ii))%lmn_size;end do
    call cprj_alloc(electronpositron%cprj_ep,ncpgr,nlmn)
    ABI_DEALLOCATE(nlmn)
   else
    ABI_ALLOCATE(electronpositron%cprj_ep,(dtset%natom,electronpositron%dimcprj))
   end if
  else
   electronpositron%dimcg   =0
   electronpositron%dimcprj =0
   electronpositron%dimeigen=0
   electronpositron%dimocc  =0
   nullify(electronpositron%cg_ep)
   nullify(electronpositron%eigen_ep)
   nullify(electronpositron%occ_ep)
   nullify(electronpositron%cprj_ep)
  end if

  optfor=0;optstr=0
  if ((dtset%optforces>0.or.dtset%ionmov/=0.or.abs(dtset%toldff)>tiny(0._dp))) optfor=1
  if (dtset%optstress>0.and.dtset%iscf>0.and.(dtset%nstep>0.or.ireadwf==1)) optstr=1

  if (optfor>0) then
   ABI_ALLOCATE(electronpositron%fred_ep,(3,dtset%natom))
   electronpositron%fred_ep(:,:)=zero
  else
   nullify(electronpositron%fred_ep)
  end if

  if (optstr>0) then
   ABI_ALLOCATE(electronpositron%stress_ep,(6))
   electronpositron%stress_ep(:)=zero
  else
   nullify(electronpositron%stress_ep)
  end if

 else !dtset%positron==0
  nullify(electronpositron)
 end if

end subroutine init_electronpositron
!!***

!----------------------------------------------------------------------

!!****f* m_electronpositron/destroy_electronpositron
!! NAME
!!  destroy_electronpositron
!!
!! FUNCTION
!!  Clean and destroy electronpositron datastructure
!!
!! SIDE EFFECTS
!!  electronpositron=<type(electronpositron_type)>=electronpositron datastructure
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,energies_copy,fourdp,rhoij_alloc
!!      rhoij_copy,rhoij_free
!!
!! SOURCE

subroutine destroy_electronpositron(electronpositron)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_electronpositron'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(electronpositron_type),pointer :: electronpositron

!************************************************************************

 !@electronpositron_type

 if (associated(electronpositron)) then

  if (associated(electronpositron%cg_ep))        then
    ABI_DEALLOCATE(electronpositron%cg_ep)
  end if
  if (associated(electronpositron%eigen_ep))     then
    ABI_DEALLOCATE(electronpositron%eigen_ep)
  end if
  if (associated(electronpositron%occ_ep))       then
    ABI_DEALLOCATE(electronpositron%occ_ep)
  end if
  if (associated(electronpositron%rhor_ep))      then
    ABI_DEALLOCATE(electronpositron%rhor_ep)
  end if
  if (associated(electronpositron%nhat_ep))      then
    ABI_DEALLOCATE(electronpositron%nhat_ep)
  end if
  if (associated(electronpositron%vha_ep))       then
    ABI_DEALLOCATE(electronpositron%vha_ep)
  end if
  if (associated(electronpositron%lmselect_ep))  then
    ABI_DEALLOCATE(electronpositron%lmselect_ep)
  end if
  if (associated(electronpositron%fred_ep))      then
    ABI_DEALLOCATE(electronpositron%fred_ep)
  end if
  if (associated(electronpositron%stress_ep))    then
    ABI_DEALLOCATE(electronpositron%stress_ep)
  end if

  if (electronpositron%has_pawrhoij_ep/=0) then
   call rhoij_free(electronpositron%pawrhoij_ep)
  end if
  if (associated(electronpositron%pawrhoij_ep))  then
    ABI_DEALLOCATE(electronpositron%pawrhoij_ep)
  end if

  if (electronpositron%dimcprj/=0) then
   call cprj_free(electronpositron%cprj_ep)
  end if
  if (associated(electronpositron%cprj_ep))  then
    ABI_DEALLOCATE(electronpositron%cprj_ep)
  end if

  nullify(electronpositron%cg_ep)
  nullify(electronpositron%eigen_ep)
  nullify(electronpositron%occ_ep)
  nullify(electronpositron%rhor_ep)
  nullify(electronpositron%nhat_ep)
  nullify(electronpositron%vha_ep)
  nullify(electronpositron%lmselect_ep)
  nullify(electronpositron%fred_ep)
  nullify(electronpositron%stress_ep)
  nullify(electronpositron%pawrhoij_ep)
  nullify(electronpositron%cprj_ep)

  electronpositron%calctype       =0
  electronpositron%particle       =-1
  electronpositron%dimcg          =0
  electronpositron%dimcprj        =0
  electronpositron%dimeigen       =0
  electronpositron%dimocc         =0
  electronpositron%has_pawrhoij_ep=0
  electronpositron%has_pos_ham    =0
  electronpositron%istep          =0
  electronpositron%istep_scf      =0

  electronpositron%posdensity0_limit=.false.
  electronpositron%scf_converged=.false.

  ABI_DEALLOCATE(electronpositron)

 end if

end subroutine destroy_electronpositron
!!***

!----------------------------------------------------------------------

!!****f* m_electronpositron/exchange_electronpositron
!! NAME
!!  exchange_electronpositron
!!
!! FUNCTION
!!  Invert electron and positron quantities between an electronpositron datastructure
!!  and current evoving variables
!!  Example: exchange electronpositron%rhor_ep and rhor
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=wavefunctions
!!  electronpositron=<type(electronpositron_type)>=electronpositron datastructure
!!  energies <type(energies_type)>=all part of total energy.
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fred(3,natom)=forces in reduced coordinates
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  occ(mband*nkpt*nsppol)=occupation number for each band at each k point
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  rhog(2,nfft)=Fourier transform of total electron/positron density
!!  rhor(nfft,nspden)=total electron/positron density (el/bohr**3)
!!  stress(6)=components of the stress tensor (hartree/bohr^3) for the
!!  vhartr(nfftf)=array for holding Hartree potential
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,energies_copy,fourdp,rhoij_alloc
!!      rhoij_copy,rhoij_free
!!
!! SOURCE

subroutine exchange_electronpositron(cg,dtset,eigen,electronpositron,energies,fred,mcg,mpi_enreg,&
&                                    nfft,ngfft,nhat,npwarr,occ,paw_an,pawrhoij,rhog,rhor,stress,vhartr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exchange_electronpositron'
 use interfaces_44_abitypes_defs
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,nfft
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),npwarr(dtset%nkpt)
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: fred(3,dtset%natom),nhat(nfft,dtset%nspden)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: rhog(2,nfft),rhor(nfft,dtset%nspden)
 real(dp),intent(inout) :: stress(6),vhartr(nfft)
 type(paw_an_type),intent(inout) :: paw_an(dtset%natom*dtset%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: iatom,ib,ibsp,icg,icgb,ifft,ii,ilm,ikpt,ispden,isppol,ispinor
 integer :: my_nspinor,nband_k,npw_k
 logical :: ltmp
 real(dp) :: rtmp
 type(energies_type) :: energies_tmp
!arrays
 integer,allocatable :: nlmn(:)
 real(dp) :: ctmp(2)
 type(cprj_type),allocatable :: cprj(:,:)  ! This is temporary waiting for cprj from SCF cycle...
 type(cprj_type),allocatable :: cprj_tmp(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_tmp(:)

!************************************************************************

 if (associated(electronpositron)) then
  if (electronpositron%particle/=EP_NOTHING) then

!  Type of particle stored
   if (electronpositron%particle==EP_ELECTRON) then
     electronpositron%particle=EP_POSITRON
   else if (electronpositron%particle==EP_POSITRON) then
     electronpositron%particle=EP_ELECTRON
   end if

!  Energies
   ctmp(1)=energies%e_electronpositron
!  ctmp(2)=energies%edc_electronpositron
   call energies_copy(electronpositron%energies_ep,energies_tmp)
   call energies_copy(energies,electronpositron%energies_ep)
   call energies_copy(energies_tmp,energies)
   energies%e_electronpositron=ctmp(1)
!  energies%edc_electronpositron=ctmp(2)
   energies%e0_electronpositron=electronpositron%e0
   electronpositron%e0=electronpositron%energies_ep%e0_electronpositron

!  Density and PAW occupation matrix
   do ispden=1,dtset%nspden
     do ifft=1,nfft
       rtmp=rhor(ifft,ispden)
       rhor(ifft,ispden)=electronpositron%rhor_ep(ifft,ispden)
       electronpositron%rhor_ep(ifft,ispden)=rtmp
     end do
     if (associated(electronpositron%nhat_ep).and.size(nhat,2)>0) then
       do ifft=1,nfft
         rtmp=nhat(ifft,ispden)
         nhat(ifft,ispden)=electronpositron%nhat_ep(ifft,ispden)
         electronpositron%nhat_ep(ifft,ispden)=rtmp
       end do
     end if
   end do
   call fourdp(1,rhog,rhor,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   if (dtset%usepaw==1.and.electronpositron%has_pawrhoij_ep==1) then
    ABI_ALLOCATE(nlmn,(dtset%ntypat))
    ABI_ALLOCATE(pawrhoij_tmp,(dtset%natom))
    do iatom=1,dtset%natom;nlmn(dtset%typat(iatom))=pawrhoij(iatom)%lmn_size;end do
    call rhoij_alloc(pawrhoij(1)%cplex,nlmn,pawrhoij(1)%nspden,pawrhoij(1)%nspinor,&
&                    pawrhoij(1)%nsppol,pawrhoij_tmp,dtset%typat,mpi_enreg=mpi_enreg,&
&                    ngrhoij=pawrhoij(1)%ngrhoij,nlmnmix=pawrhoij(1)%lmnmix_sz,&
&                    use_rhoij_=pawrhoij(1)%use_rhoij_,use_rhoijres=pawrhoij(1)%use_rhoijres)
    ABI_DEALLOCATE(nlmn)
    call rhoij_copy(pawrhoij,pawrhoij_tmp)
    call rhoij_copy(electronpositron%pawrhoij_ep,pawrhoij)
    call rhoij_copy(pawrhoij_tmp,electronpositron%pawrhoij_ep)
    if (pawrhoij_tmp(1)%ngrhoij>0.and.pawrhoij(1)%ngrhoij==0) then
     do iatom=1,dtset%natom
      sz1=pawrhoij_tmp(iatom)%ngrhoij
      sz2=pawrhoij_tmp(iatom)%cplex*pawrhoij_tmp(iatom)%lmn2_size
      sz3=pawrhoij_tmp(iatom)%nspden
      ABI_ALLOCATE(pawrhoij(iatom)%grhoij,(sz1,sz2,sz3))
      pawrhoij(iatom)%grhoij(:,:,:)=pawrhoij_tmp(iatom)%grhoij(:,:,:)
     end do
    end if
    if (pawrhoij_tmp(1)%use_rhoijres>0.and.pawrhoij(1)%use_rhoijres==0) then
     do iatom=1,dtset%natom
      sz1=pawrhoij_tmp(iatom)%cplex*pawrhoij_tmp(iatom)%lmn2_size
      sz2=pawrhoij_tmp(iatom)%nspden
      ABI_ALLOCATE(pawrhoij(iatom)%rhoijres,(sz1,sz2))
      pawrhoij(iatom)%rhoijres(:,:)=pawrhoij_tmp(iatom)%rhoijres(:,:)
     end do
    end if
    if (pawrhoij_tmp(1)%use_rhoij_>0.and.pawrhoij(1)%use_rhoij_==0) then
     do iatom=1,dtset%natom
      sz1=pawrhoij_tmp(iatom)%cplex*pawrhoij_tmp(iatom)%lmn2_size
      sz2=pawrhoij_tmp(iatom)%nspden
      ABI_ALLOCATE(pawrhoij(iatom)%rhoij_,(sz1,sz2))
      pawrhoij(iatom)%rhoij_(:,:)=pawrhoij_tmp(iatom)%rhoij_(:,:)
     end do
    end if
    if (pawrhoij_tmp(1)%lmnmix_sz>0.and.pawrhoij(1)%lmnmix_sz==0) then
     do iatom=1,dtset%natom
      ABI_ALLOCATE(pawrhoij(iatom)%kpawmix,(pawrhoij_tmp(iatom)%lmnmix_sz))
      pawrhoij(iatom)%kpawmix(:)=pawrhoij_tmp(iatom)%kpawmix(:)
     end do
    end if
    call rhoij_free(pawrhoij_tmp)
    ABI_DEALLOCATE(pawrhoij_tmp)
   else
    do iatom=1,dtset%natom
     pawrhoij(iatom)%rhoijp=zero
    end do
   end if

!  Hartree potential
   do ifft=1,nfft
    rtmp=vhartr(ifft)
    vhartr(ifft)=electronpositron%vha_ep(ifft)
    electronpositron%vha_ep(ifft)=rtmp
   end do

!  PAW LM-moment selection flags
   if (dtset%usepaw==1.and.electronpositron%lmmax>0) then
    do iatom=1,dtset%natom
     do ilm=1,paw_an(iatom)%lm_size
      ltmp=electronpositron%lmselect_ep(ilm,iatom)
      electronpositron%lmselect_ep(ilm,iatom)=paw_an(iatom)%lmselect(ilm)
      paw_an(iatom)%lmselect(ilm)=ltmp
     end do
    end do
   else
    do iatom=1,dtset%natom
     paw_an(iatom)%lmselect(:)=.true.
    end do
   end if

!  Wave-functions
   if (electronpositron%dimcg>0) then
    do ii=1,electronpositron%dimcg
     ctmp(1:2)=electronpositron%cg_ep(1:2,ii)
     electronpositron%cg_ep(1:2,ii)=cg(1:2,ii)
     cg(1:2,ii)=ctmp(1:2)
    end do
   else
    icg=0
    my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
    do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
      icgb=icg;ibsp=0
      npw_k=npwarr(ikpt);nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
      do ib=1,nband_k
       cg(:,icgb+1:icgb+my_nspinor*npw_k)=zero
       do ispinor=1,my_nspinor
        ibsp=ibsp+1;if (ibsp<my_nspinor*npw_k) cg(1,icgb+ibsp)=one
       end do
       icgb=icgb+my_nspinor*npw_k
      end do
      if (dtset%mkmem/=0) icg=icg+my_nspinor*npw_k*nband_k
     end do
    end do
   end if
   if (dtset%usepaw==1) then
    if(electronpositron%dimcprj>0) then
     ABI_ALLOCATE(nlmn,(dtset%natom))
     ABI_ALLOCATE(cprj_tmp,(dtset%natom,electronpositron%dimcprj))
     do iatom=1,dtset%natom;nlmn(iatom)=cprj(iatom,1)%nlmn;end do
     call cprj_alloc(cprj_tmp,cprj(1,1)%ncpgr,nlmn)
     ABI_DEALLOCATE(nlmn)
     call cprj_copy(electronpositron%cprj_ep,cprj_tmp)
     call cprj_copy(cprj,electronpositron%cprj_ep)
     call cprj_copy(cprj_tmp,cprj)
     call cprj_free(cprj_tmp)
     ABI_DEALLOCATE(cprj_tmp)
    else
!TO BE ACTIVATED WHEN cprj IS PRESENT
!    call cprj_set_zero(cprj)
    end if
   end if

!  Eigenvalues
   if (electronpositron%dimeigen>0) then
    do ii=1,electronpositron%dimeigen
     rtmp=eigen(ii)
     eigen(ii)=electronpositron%eigen_ep(ii)
     electronpositron%eigen_ep(ii)=rtmp
    end do
   else
    eigen(:)=9.99999_dp
   end if

!  Occupations
   if (electronpositron%dimocc>0) then
    do ii=1,electronpositron%dimocc
     rtmp=occ(ii)
     occ(ii)=electronpositron%occ_ep(ii)
     electronpositron%occ_ep(ii)=rtmp
    end do
   else
    occ(:)=9.99999_dp
   end if

!  Forces
   if (associated(electronpositron%fred_ep)) then
    do iatom=1,dtset%natom
     electronpositron%fred_ep(1:3,iatom)=fred(1:3,iatom)-electronpositron%fred_ep(1:3,iatom)
    end do
   end if

!  Stresses
   if (associated(electronpositron%stress_ep)) then
    electronpositron%stress_ep(1:6)=stress(1:6)-electronpositron%stress_ep(1:6)
   end if

  end if
 end if

end subroutine exchange_electronpositron
!!***

!----------------------------------------------------------------------

!!****f* m_electronpositron/electronpositron_calctype
!! NAME
!!  electronpositron_calctype
!!
!! FUNCTION
!!  Returns the value of the calculation type from an electronpositron
!!  structure (can be eventually unassociated)
!!
!! INPUTS
!!  electronpositron=<type(electronpositron_type)>=electronpositron datastructure
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function electronpositron_calctype(electronpositron)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'electronpositron_calctype'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(electronpositron_type),pointer :: electronpositron

!************************************************************************

 if (associated(electronpositron)) then
  electronpositron_calctype=electronpositron%calctype
 else
  electronpositron_calctype=0
 end if


end function electronpositron_calctype
!!***

!----------------------------------------------------------------------

END MODULE m_electronpositron
!!***
