!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhofermi3
!! NAME
!! rhofermi3
!!
!! FUNCTION
!! This routine computes the fixed contribution to the first-order
!! Fermi energy for metallic occupation and Q=0, as well as the
!! Fermi level charge density needed to compute the remainder of the
!! first-order Fermi energy from the self-consistent local potential
!! at each step in the iteration process.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DRH, DCA, XG, GMR, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  edocc=correction to 2nd-order total energy coming from changes of
!!    occupation
!!  eeig0=0th-order eigenenergies part of 2nd-order total energy
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mkqmem =number of k+q points which can fit in memory (GS data); 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  maximum dimension for q points in grids for nonlocal form factors
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nkpt_rbz=number of k points in the IBZ for this perturbation
!!  mpi_enreg=informations about MPI parallelization
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!   Note that this routine will NOT work with nspden==4 : at least the use of fftpac should
!!   be modified.
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!    perturbation
!!  ntypat=number of types of atoms in unit cell.
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band and k (usually 2)
!!  phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional real space primitive translations
!!  symaf1(nsym1)=(anti)ferromagnetic part of symmetry operations
!!  symrl1(3,3,nsym1)=3x3 matrices of the group symmetries
!!  ucvol=unit cell volume in bohr**3.
!!  wfftgs,wfftkq=struct info for GS wf disk files.
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= spherical harmonics for each G and k+g point
!!  ylmgr1(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics at k+q
!!
!!
!! OUTPUT
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues
!!   (hartree) - only digonal elements computed here
!!  fe1fixed=fixed contribution to the first-order Fermi energy
!!   (nonlocal and kinetic in the case of strain)
!!  rhorfermi(cplex*nfft,nspden)=RF electron density in electrons/bohr**3.
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      clsopn,fftpac,hdr_skip,kpg3,kpgstr,leave_new,leave_test,mkffnl,mkkin
!!      mkkpg,occeig,ph1d3d,rdnpw,rwwf,sphereboundary,status,symrhg,timab
!!      wfkfermi3,wrtout,xcomm_world,xdefineoff,xmaster_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rhofermi3(atindx,atindx1,cg,cgq,cplex,&
& doccde_rbz,docckqde,dtfil,dtset,&
& edocc,eeig0,eigenq,eigen0,eigen1,&
& fe1fixed,gmet,gprimd,idir,&
& ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,mband,mgfft,&
& mkmem,mkqmem,mk1mem,mpi_enreg,mpsang,mpw,mpw1,&
& natom,nattyp,nband_rbz,nfft,nkpt_rbz,npwarr,npwar1,nspden,&
& nsppol,nsym1,ntypat,occkq,occ_rbz,phnons1,&
& ph1d,prtvol,psps,rhorfermi,rmet,rprimd,symaf1,symrl1,ucvol,&
& wfftgs,wfftkq,wtk_rbz,xred,ylm,ylm1,ylmgr1)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhofermi3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_65_nonlocal
 use interfaces_67_common
 use interfaces_72_response, except_this_one => rhofermi3
!End of the abilint section

 implicit none

!Arguments -------------------------------
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
!scalars
 integer,intent(in) :: cplex,idir,ipert,mband,mgfft,mk1mem,mkmem,mkqmem,mpsang
 integer,intent(in) :: mpw,mpw1,natom,nfft,nkpt_rbz,nspden,nsppol,nsym1,ntypat
 integer,intent(in) :: prtvol
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: edocc,eeig0
 real(dp),intent(out) :: fe1fixed
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wfftgs,wfftkq
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom)
 integer,intent(in) :: irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nattyp(ntypat),nband_rbz(nkpt_rbz*nsppol)
 integer,intent(in) :: npwar1(nkpt_rbz,2),npwarr(nkpt_rbz,2),symaf1(nsym1)
 integer,intent(in) :: symrl1(3,3,nsym1)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*mband*mkqmem*nsppol)
 real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz),occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: occkq(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3),wtk_rbz(nkpt_rbz)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(out) :: eigen1(2*mband*mband*nkpt_rbz*nsppol)
 real(dp),intent(out) :: rhorfermi(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=17
! integer,save :: count=0 ! used in a commented debug section below
 integer :: bd2tot_index,bdtot_index,buffer_size,counter,dimffnl1
 integer :: dimffnlk,formeig,ia,iatom,iband
 integer :: icg,icgq,ider,idir0,ierr,iexit,ii,ikg,ikg1,ikpt,ilm,ilmn,index
 integer :: index1,ipw,iscf_mod,isppol,istr,istwf_k
 integer :: itypat,master,matblk,mbd2kpsp,mbdkpsp,mcgq,mcgq_disk,me,muig,n1,n2
 integer :: n3,n4,n5,n6,nband_k,nkpg,nkpg1,npw1_k,npw_k,nsp,nvloc
 integer :: spaceComm,spaceworld,tim_rwwf
 real(dp) :: arg,fe1norm,invfe1norm,wtk_k
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamkq
!arrays
 integer :: pspso_typ(1)
 integer,allocatable :: gbound(:,:),indlmn_typ(:,:,:),kg1_k(:,:),kg_dum(:,:)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kpoint(3),kpq(3),tsec(2)
 real(dp) :: ylmgr_dum(1)
 real(dp),allocatable :: buffer1(:),cgq_disk(:,:),dkinpw(:),doccde_k(:)
 real(dp),allocatable :: doccde_kq(:),eig0_k(:),eig0_kq(:),eig1_k(:)
 real(dp),allocatable :: eigkq_dum(:),ekb_typ(:,:,:),fe1fixed_k(:),fe1norm_k(:)
 real(dp),allocatable :: ffnl1(:,:,:,:),ffnlk(:,:,:,:),ffnlkq(:,:,:,:)
 real(dp),allocatable :: ffspl_typ(:,:,:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:)
 real(dp),allocatable :: occ_k(:),occ_kq(:),occkq_dum(:),ph3d(:,:,:)
 real(dp),allocatable :: rhoaug1(:,:,:),rhogfermi(:,:),rocceig(:,:)
 real(dp),allocatable :: ylm1_k(:,:),ylm_k(:,:),ylmgr1_k(:,:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' rhofermi3 : enter '
!count=count+1
!write(std_out,*)' count=',count
!if(count==13)stop
!write(std_out,*)' rhofermi3 : kg1(:,66)',kg1(:,66)
!write(std_out,*)' xred=',xred
!stop
!ENDDEBUG

 if(psps%usepaw==1)then
   write(message, '(4a)' ) ch10,&
&   ' rhofermi3 :  ERROR -',ch10,&
&   '  Not yet allowed within PAW !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Keep track of total time spent in rhofermi3
 call timab(121,1,tsec)
 call timab(124,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Init mpi_comm
 call xcomm_world(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)
!Init master
 call xmaster_init(mpi_enreg,master)

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' rhofermi3 : enter '
   call wrtout(std_out,message,'COLL')
 end if

 iscf_mod=dtset%iscf

!The value of iscf must be modified if ddk perturbation, see loper3.f
 if(ipert==natom+1)iscf_mod=-3

 edocc=0.0_dp ; eeig0=0.0_dp
 fe1fixed=0.0_dp ; fe1norm = 0.0_dp
 bdtot_index=0
 bd2tot_index=0
 icg=0
 icgq=0
 mbdkpsp=mband*nkpt_rbz*nsppol
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol

 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))

!Initialize rhorfermi
 rhorfermi(:,:)=0.0_dp

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 ABI_ALLOCATE(rhoaug1,(cplex*n4,n5,n6))

 ABI_ALLOCATE(kg_dum,(3,0))

!Prepare GS k wf file for reading if mkmem==0
 if (mkmem==0) then
   call clsopn(wfftgs)
   call hdr_skip(wfftgs,ierr)
!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wfftgs,mpi_enreg,nband_rbz,npwarr,dtset%nspinor,nsppol,nkpt_rbz)
 end if

!Prepare GS k+q wf file for reading if mkqmem==0
 if (mkqmem==0) then
   call clsopn(wfftkq)
   call hdr_skip(wfftkq,ierr)
   mcgq_disk=mpw1*dtset%nspinor*mband
   ABI_ALLOCATE(cgq_disk,(2,mcgq_disk))
!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wfftkq,mpi_enreg,nband_rbz,npwar1,dtset%nspinor,nsppol,nkpt_rbz)
 end if

 nvloc=1

!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 ABI_ALLOCATE(gs_hamkq%atindx,(natom))
 ABI_ALLOCATE(gs_hamkq%atindx1,(natom))
 ABI_ALLOCATE(gs_hamkq%gbound,(2*mgfft+8,2))
 ABI_ALLOCATE(gs_hamkq%indlmn,(6,psps%lmnmax,ntypat))
 ABI_ALLOCATE(gs_hamkq%nattyp,(ntypat))
 ABI_ALLOCATE(gs_hamkq%phkxred,(2,natom))
 ABI_ALLOCATE(gs_hamkq%ph1d,(2,3*(2*mgfft+1)*natom))
 ABI_ALLOCATE(gs_hamkq%pspso,(ntypat))
 ABI_ALLOCATE(gs_hamkq%xred,(3,natom))

!Initialize most components of the Ground-state Hamiltonian ar k+q
 gs_hamkq%atindx(:)  =atindx(:)
 gs_hamkq%atindx1(:) =atindx1(:)
 gs_hamkq%gmet(:,:)  =gmet(:,:)
 gs_hamkq%gprimd(:,:)=gprimd(:,:)
 gs_hamkq%indlmn(:,:,:)=psps%indlmn(:,:,:)
 gs_hamkq%lmnmax     =psps%lmnmax
 gs_hamkq%mgfft      =mgfft
 gs_hamkq%mpsang     =mpsang
 gs_hamkq%mpssoang   =psps%mpssoang
 gs_hamkq%natom      =natom
 gs_hamkq%nattyp(:)  =nattyp(:)
 gs_hamkq%nfft       =nfft
 gs_hamkq%ngfft(:)   =dtset%ngfft(:)
 gs_hamkq%nloalg(:)  =dtset%nloalg(:)
 gs_hamkq%nspinor    =dtset%nspinor
 gs_hamkq%ntypat     =ntypat
 gs_hamkq%nvloc      =nvloc
 gs_hamkq%n4         =n4
 gs_hamkq%n5         =n5
 gs_hamkq%n6         =n6
 gs_hamkq%usepaw     =psps%usepaw
 gs_hamkq%ph1d(:,:)  =ph1d(:,:)
 gs_hamkq%pspso(:)   =psps%pspso(:)
 gs_hamkq%ucvol      =ucvol
 gs_hamkq%useylm     =psps%useylm
 gs_hamkq%use_gpu_cuda=dtset%use_gpu_cuda
 gs_hamkq%xred(:,:)  =xred(:,:)

!Special case of array ekb:
!Not the same meaning for norm-cons. psps and paw
 if (psps%usepaw==1) then
   gs_hamkq%dimekb1=psps%dimekb
   gs_hamkq%dimekb2=natom
   ABI_ALLOCATE(gs_hamkq%ekb,(psps%dimekb,natom,dtset%nspinor**2))
!  ekb will receive Dij coeffs later (cf below)
 else
   gs_hamkq%dimekb1=psps%dimekb
   gs_hamkq%dimekb2=ntypat
   ABI_ALLOCATE(gs_hamkq%ekb,(gs_hamkq%dimekb1,gs_hamkq%dimekb2,dtset%nspinor**2))
   gs_hamkq%ekb(:,:,1)=psps%ekb(:,:)
   if (dtset%nspinor==2) then
     gs_hamkq%ekb(:,:,2)=psps%ekb(:,:)
     gs_hamkq%ekb(:,:,3:4)=zero
   end if
 end if

!LOOP OVER SPINS
 do isppol=1,nsppol

   if (nsppol/=1) then
     write(message,*)' ****  In rhofermi3 for isppol=',isppol
     call wrtout(std_out,message,'COLL')
   end if

!  PAW: transfer Dij coeffs into array ekb     ! TO BE COMPLETED LATER
!  if (psps%usepaw==1) then
!  do ispden=1,dtset%nspinor**2
!  isp=isppol;if (dtset%nspinor==2) isp=ispden
!  do iatom=1,natom
!  dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
!  do ilmn=1,dimdij
!  gs_hamk%ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
!  end do
!  if(dimdij+1<=gs_hamk%dimekb1) gs_hamk%ekb(dimdij+1:gs_hamk%dimekb1,iatom,ispden)=zero
!  end do
!  end do
!  end if
!  end do

!  Rewind kpgsph data file if needed:
   if (mkmem==0) rewind(dtfil%unkg)
   if (mk1mem==0) rewind(dtfil%unkg1)
   if (mkmem==0.and.psps%useylm==1) rewind(dtfil%unylm)
   if (mk1mem==0.and.psps%useylm==1) rewind(dtfil%unylm1)
   ikg=0;ikg1=0

   rhoaug1(:,:,:)=0.0_dp

   call timab(125,1,tsec)

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     counter=100*ikpt+isppol
     call status(counter,dtfil%filstat,iexit,level,'loop ikpt     ')
     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt,1)
     npw1_k=npwar1(ikpt,1)

     if(mpi_enreg%paral_compil_kpt==1)then
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&       -mpi_enreg%me))/=0) then
         eigen1(1+bd2tot_index : 2*nband_k**2+bd2tot_index) = 0.0_dp
         bdtot_index=bdtot_index+nband_k
         bd2tot_index=bd2tot_index+2*nband_k**2

!        Skip the rest of the k-point loop
         cycle
       end if
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     ABI_ALLOCATE(ylm1_k,(npw1_k,mpsang*mpsang*psps%useylm))
     if (ipert==natom+1.or.ipert==natom+3.or.ipert==natom+4) then
       ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,mpsang*mpsang*psps%useylm))
     end if

!    Continue to initialize the Hamiltonian at k+q
     kpoint(:)=kpt_rbz(:,ikpt)
     kpq(:)=kpoint(:)+dtset%qptn(1:3)
     gs_hamkq%istwf_k    =istwf_k
     gs_hamkq%kpoint(:)  =kpq(:)
     gs_hamkq%npw        =npw1_k

     ABI_ALLOCATE(doccde_k,(nband_k))
     ABI_ALLOCATE(doccde_kq,(nband_k))
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(eig0_kq,(nband_k))
     ABI_ALLOCATE(eig1_k,(2*nband_k**2))
     ABI_ALLOCATE(fe1fixed_k,(nband_k))
     ABI_ALLOCATE(fe1norm_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(occ_kq,(nband_k))
     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     ABI_ALLOCATE(rocceig,(nband_k,nband_k))

     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
     eig0_kq(:)=eigenq(1+bdtot_index:nband_k+bdtot_index)
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
     occ_kq(:)=occkq(1+bdtot_index:nband_k+bdtot_index)
     doccde_k(:)=doccde_rbz(1+bdtot_index:nband_k+bdtot_index)
     doccde_kq(:)=docckqde(1+bdtot_index:nband_k+bdtot_index)

!    For each pair of active bands (m,n), generates the ratios
!    rocceig(m,n)=(occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
!    and decide to which band to attribute it.
     call occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,&
&     dtset%occopt,occ_k,occ_kq,rocceig)

     if (mkmem==0) then
!      Read (k+G) basis sphere data (same for each spin)
       call status(counter,dtfil%filstat,iexit,level,'read kg data  ')
       nsp=dtset%nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,dtfil%unkg)
!      Read k+g data
       read (dtfil%unkg) kg_k(1:3,1:npw_k)
       call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

!      Eventually read (k+G) spherical harmonics
       if (psps%useylm==1) then
         read(dtfil%unylm)
         read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
       end if

     else

       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)
       if (psps%useylm==1) then
         do ilm=1,mpsang*mpsang
           ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
         end do
       end if

!      End if for choice governed by mkmem
     end if

     if (mkqmem==0) then
       call status(counter,dtfil%filstat,iexit,level,'read k+q wfs  ')
!      Take care of GS wavefunctions at k+q , nothing to be done for kpg sphere
!      Eigenvalues already in eigen1
       tim_rwwf=17
       ABI_ALLOCATE(eigkq_dum,(mband))
       ABI_ALLOCATE(occkq_dum,(mband))
       call rwwf(cgq_disk,eigkq_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcgq_disk,mpi_enreg,nband_k,nband_k,&
&       npw1_k,dtset%nspinor,occkq_dum,-2,0,tim_rwwf,wfftkq)
       ABI_DEALLOCATE(eigkq_dum)
       ABI_DEALLOCATE(occkq_dum)

!      End if for choice governed by mkqmem
     end if

     wtk_k=wtk_rbz(ikpt)

     kg1_k(:,:) = 0
     if (mk1mem==0) then

!      Read (k+q+G) basis sphere data (same for each spin)
       call status(counter,dtfil%filstat,iexit,level,'read kg1 data ')
       nsp=dtset%nspinor
       call rdnpw(ikpt,isppol,nband_k,npw1_k,nsp,0,dtfil%unkg1)
       read (dtfil%unkg1) kg1_k(1:3,1:npw1_k)
       call sphereboundary(gs_hamkq%gbound,istwf_k,kg1_k,mgfft,npw1_k)

!      Eventually read (k+q+G) spherical harmonics
       if (psps%useylm==1) then
         read(dtfil%unylm1)
         if (ipert==natom+1.or.ipert==natom+3.or.ipert==natom+4) then
           read(dtfil%unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang),&
&           (((ylmgr1_k(muig,ii,ilm),muig=1,npw1_k),ii=1,3),ilm=1,mpsang*mpsang)
         else
           read(dtfil%unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang)
         end if
       end if

     else

       kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
       call sphereboundary(gs_hamkq%gbound,istwf_k,kg1_k,mgfft,npw1_k)
       if (psps%useylm==1) then
         do ilm=1,mpsang*mpsang
           ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
         end do
         if (ipert==natom+1.or.ipert==natom+3.or.ipert==natom+4) then
           do ilm=1,mpsang*mpsang
             do ii=1,3
               ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
             end do
           end do
         end if
       end if

!      End if for choice governed by mk1mem
     end if

!    Set up the ground-state Hamiltonian, and some parts of the 1st order Hamiltonian

!    Note that not all these arrays should be allocated in the general case
!    when wtk_k vanishes
     ABI_ALLOCATE(ekb_typ,(psps%dimekb,1,dtset%nspinor**2))
     ABI_ALLOCATE(indlmn_typ,(6,psps%lmnmax,1))
     ABI_ALLOCATE(dkinpw,(npw_k))
     ABI_ALLOCATE(kinpw1,(npw1_k))
     dimffnlk=1
     ABI_ALLOCATE(ffnlk,(npw_k,dimffnlk,psps%lmnmax,1))

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=0;if(ipert>=1.and.ipert<=natom) nkpg=3*dtset%nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!    Preparation of the kinetic and non-local contributions
     if( ipert>=1 .and. ipert<=natom )then

!      Compute nonlocal form factors ffnlk at (k+G), for the displaced atom only.
       call status(counter,dtfil%filstat,iexit,level,'call mkffnl(0)')
       ider=0;idir0=0
!      Need to transfer the infos relative to the displaced atom: ekb, indlmn, pspso, ffspl
       itypat=dtset%typat(ipert)
!      if (psps%usepaw==1) then
!      ekb_typ(:,1)=paw_ij(iatom)%dij() ! TO UNCOMMENT LATER
!      else
       ekb_typ(:,1,1)=psps%ekb(:,itypat)
       if (dtset%nspinor==2) then
         ekb_typ(:,1,2)=psps%ekb(:,itypat)
         ekb_typ(:,1,3:4)=zero
       end if
!      end if
       indlmn_typ(:,:,1)=psps%indlmn(:,:,itypat)
       pspso_typ(1)=psps%pspso(itypat)
       ABI_ALLOCATE(ffspl_typ,(psps%mqgrid_ff,2,psps%lnmax,1))
       ffspl_typ(:,:,:,1)=psps%ffspl(:,:,:,itypat)
       call mkffnl(psps%dimekb,dimffnlk,ekb_typ(:,1,1),ffnlk,ffspl_typ,&
&       gmet,gprimd,ider,idir0,indlmn_typ,kg_k,kpg_k,kpoint,psps%lmnmax,&
&       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,1,&
&       pspso_typ,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
       ABI_DEALLOCATE(ffspl_typ)

     else if(ipert==natom+1)then

!      This flag is needed for the call to mkffnl later
       ider=1;idir0=idir

       call status(counter,dtfil%filstat,iexit,level,'call kpg3     ')
!      Compute the derivative of the kinetic operator vs k in kinpw
       call kpg3(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,idir,kg_k,kpoint,npw_k)
!      No non-local form factor is to be calculated here in this case

     else if(ipert==natom+2)then

       ider=0;idir0=0
!      dkinpw and ffnlk are not needed ...

!      section for strain perturbation
     else if(ipert==natom+3 .or. ipert==natom+4)then

!      istr is 1,2,...,6 and indicates cartesian strain component
       if(ipert==natom+3) then
         istr=idir
       else
         istr=idir+3
       end if

!      Derivatives needed for strain perturbation
!      This flag is needed for the call to mkffnl later
       ider=1;idir0=-istr

       call status(counter,dtfil%filstat,iexit,level,'call kpgstr   ')
!      Compute the derivative of the kinetic operator vs strain in dkinpw
       call kpgstr(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,gprimd,istr,&
&       kg_k,kpoint,npw_k)

!      endsection for strain perturbation

     end if

!    Compute (1/2) (2 Pi)**2 (k+q+G)**2:
     call status(counter,dtfil%filstat,iexit,level,'call mkkin(1) ')
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg1_k,kinpw1,kpq,npw1_k)

!    Compute (k+q+G) vectors (only if useylm=1)
     nkpg1=0;if(ipert>=1.and.ipert<=natom) nkpg1=3*dtset%nloalg(5)
     ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
     if (nkpg1>0) call mkkpg(kg1_k,kpg1_k,kpq,nkpg1,npw1_k)

!    Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
     call status(counter,dtfil%filstat,iexit,level,'call mkffnl(1)')
     dimffnl1=1+ider;if (ider==1.and.idir0==0) dimffnl1=dimffnl1+2*psps%useylm
     ABI_ALLOCATE(ffnlkq,(npw1_k,dimffnl1,psps%lmnmax,1))
     ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,ntypat))
     call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gmet,gprimd,ider,idir0,&
&     psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
&     npw1_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)

     if( ipert>=1 .and. ipert<=natom )then
!      Copy the part needed for the displaced atom, in ffnlkq.
       do ilmn=1,psps%lmnmax
         do ii=1,dimffnl1
           do ipw=1,npw1_k
             ffnlkq(ipw,ii,ilmn,1)=ffnl1(ipw,ii,ilmn,itypat)
           end do
         end do
       end do
     end if

!    Allocate the arrays phkxred and ph3d, compute phkxred
!    and eventually ph3d.
!    NOTE : in this RF case, uses kpq instead of kpt
     call status(counter,dtfil%filstat,iexit,level,'make phkxred  ')
     do ia=1,natom
       iatom=atindx(ia)
       arg=two_pi*(kpq(1)*xred(1,ia)+kpq(2)*xred(2,ia)+kpq(3)*xred(3,ia))
       gs_hamkq%phkxred(1,iatom)=cos(arg)
       gs_hamkq%phkxred(2,iatom)=sin(arg)
     end do

!    Note : use npw1_k
     call status(counter,dtfil%filstat,iexit,level,'make ph3d     ')
     if(dtset%nloalg(1)<=0)then
!      Here, only the allocation, not the precomputation.
       matblk=dtset%nloalg(4)
       ABI_ALLOCATE(ph3d,(2,npw1_k,matblk))
     else
!      Here, allocation as well as precomputation
       matblk=natom
       ABI_ALLOCATE(ph3d,(2,npw1_k,matblk))
       call ph1d3d(1,natom,kg1_k,matblk,natom,npw1_k,n1,n2,n3,&
&       gs_hamkq%phkxred,ph1d,ph3d)
     end if
     gs_hamkq%matblk=matblk

     call status(counter,dtfil%filstat,iexit,level,'call wfkfermi3  ')

!    Compute fixed contributions to 1st-order Fermi energy and
!    Fermi level charge density
     fe1fixed_k(:)=zero
     fe1norm_k(:)=zero
     if(mkqmem/=0)then
       mcgq=mpw1*dtset%nspinor*mband*mkqmem*nsppol
!      Note that wfkfermi3 is called with kpoint, while kpt is used
!      inside wfkfermi3
       call wfkfermi3(cg,cgq,cplex,psps%dimekb,dimffnlk,dimffnl1,dkinpw,dtfil,dtset,&
&       eig1_k,ekb_typ,fe1fixed_k,fe1norm_k,ffnlk,ffnlkq,ffnl1,gbound,gs_hamkq,&
&       icg,icgq,idir,ikpt,indlmn_typ,ipert,isppol,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpoint,&
&       psps%lmnmax, matblk,mband,mcgq,mgfft,mkmem,mpi_enreg,&
&       psps%mpsang,psps%mpssoang,mpw,natom,nband_k,nkpg,nkpg1,&
&       npw_k,npw1_k,nsppol,ntypat,n4,n5,n6,occ_k,ph3d,prtvol,&
&       rhoaug1,rocceig,wfftgs,wtk_k)
     else if(mkqmem==0)then
       mcgq=mpw1*dtset%nspinor*mband
!      Note that wfkfermi3 is called with kpoint, while kpt is used
!      inside wfkfermi3
       call wfkfermi3(cg,cgq_disk,cplex,psps%dimekb,dimffnlk,dimffnl1,dkinpw,dtfil,dtset,&
&       eig1_k,ekb_typ,fe1fixed_k,fe1norm_k,ffnlk,ffnlkq,ffnl1,gbound,gs_hamkq,&
&       icg,icgq,idir,ikpt,indlmn_typ,ipert,isppol,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpoint,&
&       psps%lmnmax, matblk,mband,mcgq,mgfft,mkmem,mpi_enreg,&
&       psps%mpsang,psps%mpssoang,mpw,natom,nband_k,nkpg,nkpg1,&
&       npw_k,npw1_k,nsppol,ntypat,n4,n5,n6,occ_k,ph3d,prtvol,&
&       rhoaug1,rocceig,wfftgs,wtk_k)
     end if

     ABI_DEALLOCATE(dkinpw)
     ABI_DEALLOCATE(ekb_typ)
     ABI_DEALLOCATE(ffnlk)
     ABI_DEALLOCATE(ffnlkq)
     ABI_DEALLOCATE(ffnl1)
     ABI_DEALLOCATE(indlmn_typ)
     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kinpw1)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(kpg1_k)

     call status(counter,dtfil%filstat,iexit,level,'after wfkfermi3 ')

!    Save eigenvalues (hartree)
     eigen1 (1+bd2tot_index : 2*nband_k**2+bd2tot_index) = eig1_k(:)

!    Accumulate sum over k points for 1st-order Fermi energy components
     do iband=1,nband_k
       fe1fixed=fe1fixed+wtk_k*occ_k(iband)*fe1fixed_k(iband)
       fe1norm=fe1norm+wtk_k*occ_k(iband)*fe1norm_k(iband)
     end do

     ABI_DEALLOCATE(doccde_k)
     ABI_DEALLOCATE(doccde_kq)
     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(eig0_kq)
     ABI_DEALLOCATE(eig1_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(occ_kq)
     ABI_DEALLOCATE(fe1fixed_k)
     ABI_DEALLOCATE(fe1norm_k)
     ABI_DEALLOCATE(rocceig)

!    Keep track of total number of bands (all k points so far, even for
!    k points not treated by me)
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2

!    Shift array memory
     if (mkmem/=0) then
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if
     if (mkqmem/=0) then
       icgq=icgq+npw1_k*dtset%nspinor*nband_k
     end if
     if (mk1mem/=0) then
       ikg1=ikg1+npw1_k
     end if
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     if (ipert==natom+1.or.ipert==natom+3.or.ipert==natom+4)  then
       ABI_DEALLOCATE(ylmgr1_k)
     end if
!    End big k point loop
   end do

   call timab(125,2,tsec)

!  Transfer density on augmented fft grid to normal fft grid in real space
!  Also take into account the spin.
   if(iscf_mod>0) &
&   call fftpac(isppol,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rhorfermi,rhoaug1,1)

!  End loop over spins
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
   call timab(166,1,tsec)
!  BEGIN TF_CHANGES
   call leave_test()
!  END TF_CHANGES
   write(message,*) 'rhofermi3: loop on k-points and spins done in parallel'
   call wrtout(std_out,message,'COLL')
   call timab(166,2,tsec)
 end if

 ABI_DEALLOCATE(gs_hamkq%atindx)
 ABI_DEALLOCATE(gs_hamkq%atindx1)
 ABI_DEALLOCATE(gs_hamkq%ekb)
 ABI_DEALLOCATE(gs_hamkq%gbound)
 ABI_DEALLOCATE(gs_hamkq%indlmn)
 ABI_DEALLOCATE(gs_hamkq%nattyp)
 ABI_DEALLOCATE(gs_hamkq%phkxred)
 ABI_DEALLOCATE(gs_hamkq%pspso)
 ABI_DEALLOCATE(gs_hamkq%ph1d)
 ABI_DEALLOCATE(gs_hamkq%xred)


 ABI_DEALLOCATE(rhoaug1)
 if(mkqmem==0)ABI_DEALLOCATE(cgq_disk)

 call status(0,dtfil%filstat,iexit,level,'after loops   ')

 call timab(124,2,tsec)

 if(mpi_enreg%paral_compil_kpt==1)then

   call timab(129,1,tsec)

   buffer_size=cplex*nfft*nspden+2+mbd2kpsp+mbdkpsp
   ABI_ALLOCATE(buffer1,(buffer_size))

!  Pack rhorfermi,fe1fixed,fe1norm
   index1=cplex*nfft*nspden
   buffer1(1:index1)=reshape(rhorfermi,(/index1/))
   buffer1(index1+1       )=fe1fixed
   buffer1(index1+2       )=fe1norm
   index=index1+2
   bdtot_index=0
   bd2tot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       buffer1(index+1:index+2*nband_k**2)=&
&       eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)
       bd2tot_index=bd2tot_index+2*nband_k**2
       index=index+2*nband_k**2
     end do
   end do
   if(index/=buffer_size)buffer1(index+1:buffer_size)=0.0_dp

!  Build sum of everything
   call timab(48,1,tsec)
   write(message, '(a,i8,a)' ) &
&   ' rhofermi3 : MPI_ALLREDUCE, buffer of size',&
&   8*buffer_size,' bytes'
   call wrtout(std_out,message,'COLL')
!  BEGIN TF_CHANGES
   call xcomm_world(mpi_enreg,spaceworld)
!  END TF_CHANGES
   call xsum_mpi(buffer1,buffer_size,spaceworld,ierr)
   call timab(48,2,tsec)

!  Unpack the final result
   index1=cplex*nfft*nspden
   if(iscf_mod > 0) then
     rhorfermi(:,:)  =reshape(buffer1(1:index1),(/cplex*nfft,nspden/))
   end if
   fe1fixed         =buffer1(index1+1)
   fe1norm          =buffer1(index1+2)
   index=index1+2
   bdtot_index=0
   bd2tot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)=&
&       buffer1(index+1:index+2*nband_k**2)
       bd2tot_index=bd2tot_index+2*nband_k**2
       index=index+2*nband_k**2
     end do
   end do
   ABI_DEALLOCATE(buffer1)
   call timab(129,2,tsec)
 end if ! if kpt parallel

 call timab(127,1,tsec)

!Normalize the fixed part of fermie1
 if (abs(fe1norm) > tol10) then
   invfe1norm=1.0_dp/fe1norm
 else
   invfe1norm = zero
 end if
 fe1fixed=fe1fixed*invfe1norm

!DEBUG
!write(std_out,*)' rhofermi3 : fe1fixed=',fe1fixed
!ENDDEBUG

!Compute rhogfermi, and symmetrize the density

 call status(0,dtfil%filstat,iexit,level,'call symrhg   ')

 ABI_ALLOCATE(rhogfermi,(2,nfft))
 call symrhg(cplex,gprimd,&
& irrzon1,mpi_enreg,nfft,nfft,dtset%ngfft,nspden,nsppol,nsym1,dtset%paral_kgb,phnons1,&
& rhogfermi,rhorfermi,rprimd,symaf1,symrl1)
!To have this work with parallel fft, the dim of irrzon and phnons should be nfftot
!We now have both rho(r) and rho(G), symmetrized, and if nsppol=2
!we also have the spin-up density, symmetrized, in rhorfermi(:,2).

 ABI_DEALLOCATE(rhogfermi)
!Normalize the Fermi level charge density
 rhorfermi(:,:)=invfe1norm*rhorfermi(:,:)

 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)
 ABI_DEALLOCATE(kg_dum)

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
   write(message,'(a1,a,a1,a,i2,a)') ch10,' rhofermi3 : exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(127,2,tsec)
 call timab(121,2,tsec)

end subroutine rhofermi3
!!***
