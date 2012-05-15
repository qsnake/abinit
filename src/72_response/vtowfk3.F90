!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtowfk3
!! NAME
!! vtowfk3
!!
!! FUNCTION
!! This routine compute the partial density at a given k-point,
!! for a given spin-polarization, from a fixed potential (vlocal1).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG, AR, DRH, MB, MVer,XW, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cgq(2,mcgq)=array for planewave coefficients of wavefunctions.
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of
!!    RF wavefunctions at k,q.
!!  cplex=1 if rhoaug1 is real, 2 if rhoaug1 is complex
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,mcprjq)= wave functions at k+q projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  cpus= cpu time limit in seconds
!!  dimcprj(natom*usepaw)=array of dimensions of arrays cprj, cprjq (ordered by atom-type)
!!  dimekb=first dimension of ekb (see ekb_typ below)
!!  dime1kb=first dimension of e1kb (see e1kb below)
!!  dim_eig2rf = dimension for the second order eigenvalues
!!  dimffnlk=second dimension of ffnlk (1+number of derivatives)
!!  dimffnl1=second dimension of ffnl1 and ffnlkq (1+number of derivatives)
!!  dimphkxred=second dimension of phkxred
!!  dkinpw(npw_k)=derivative of the (modified) kinetic energy for
!!    each plane wave at k (Hartree)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig0_k(nband_k)=GS eigenvalues at k (hartree)
!!  eig0_kq(nband_k)=GS eigenvalues at k+Q (hartree)
!!  ekb_typ(dimekb,1,nspinor**2)=
!!    ->Norm conserving : (Real) Kleinman-Bylander energies (hartree)
!!                        for the displaced atom
!!                        for number of basis functions (l,n) (lnmax)
!!                        dimekb=lnmax
!!    ->PAW : (Real, symmetric) Frozen part of Dij coefficients
!!                        to connect projectors
!!                        for the displaced atom
!!                        for number of basis functions (l,m,n) (lmnmax)
!!                        dimekb=lmnmax*(lmnmax+1)/2
!!  e1kbfr(dime1kb,dimekb2,use1ekb*nspinor**2)=frozen part of 1st der. of ekb
!                     for the current pertubation (not depending on VHxc^(1))
!!  e1kbsc(dime1kb,dimekb2,use1ekb*nspinor**2)=self-consistent part of 1st der. of ekb
!!                                 for the current pertubation (depending on VHxc^(1))
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  ffnlk(npw_k,dimffnlk,lmnmax,1+usepaw(ntypat-1))=nonloc form factors at k, for the displaced atom.
!!  ffnlkq(npw1_k,dimffnl1,lmnmax,1)=nonloc form fact at k+q for the displaced atom
!!  ffnl1(npw1_k,dimffnl1,lmnmax,ntypat)=nonloc form factors at k+q
!!  gbound(2*mgfft+8,2)=G sphere boundary
!!  grad_berry(2,mpw1,dtefield%nband_occ) = the gradient of the Berry phase term
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  ibg1=shift to be applied on the location of data in the array cprj1
!!  icg=shift to be applied on the location of data in the array cg
!!  icgq=shift to be applied on the location of data in the array cgq
!!  icg1=shift to be applied on the location of data in the array cg1
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  indlmn_typ(6,lmnmax,1)=indlmn info for the displaced atom
!!  ipert=type of the perturbation
!!  isppol=1 index of current spin component
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kinpw1(npw1_k)=(modified) kinetic energy for each plane wave at k+q (Hartree)
!!  kpg_k(npw_k,nkpg)= (k+G) components at k (only if useylm=1)
!!  kpg1_k(npw1_k,nkpg1)= (k+G) components at k+q (only if useylm=1)
!!  kpt(3)=reduced coordinates of k points.
!!  lmnmax= max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mband=maximum number of bands
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nkpg,nkpg1=second dimensions of kpg_k and kpg1_k (0 if useylm=0)
!!  nkpt=number of k points
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!    (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  n4,n5,n6 used for dimensioning real space arrays
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  phkxred(2,dimphkxred)=phase factors exp(2 pi kpoint.xred) at k
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhoaug1(cplex*n4,n5,n6)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output)
!!  rocceig(nband_k,nband_k)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!    if this
!!   ratio has been attributed to the band n (second argument), zero otherwise
!!  sij_typ(dimekb,usepaw)=-PAW only- overlap matrix components for the current perturbation
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  usee1kb=1 if ekb derivatives (e1kbsc, e1kbfr) exist
!!  wffddk=struct info for wf ddk file.
!!  wffnew=struct info for OUTPUT 1st-order wf file
!!  wffnow=struct info for INPUT 1st-order wf file
!!  wfftgs=struct info for GS wf disk files.
!!  vlocal(n4,n5,n6)= GS local potential in real space, on the augmented fft grid
!!  vlocal1(cplex*n4,n5,n6)= RF local pot. in real space, on the augm. fft grid
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q. They are orthogonalized to the occupied states.
!!  cg1_active(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)=pw coefficients of RF
!!    wavefunctions at k,q. They are orthogonalized to the active. Only needed for ieigrf/=0
!!  edocc_k(nband_k)=correction to 2nd-order total energy coming
!!      from changes of occupation
!!  eeig0_k(nband_k)=zero-order eigenvalues contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  eig1_k(2*nband_k**2)=first-order eigenvalues (hartree)
!!  ek0_k(nband_k)=0-order kinetic energy contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  ek1_k(nband_k)=1st-order kinetic energy contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  eloc0_k(nband_k)=zero-order local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  enl0_k(nband_k)=zero-order non-local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  enl1_k(nband_k)=first-order non-local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  gh1c_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(1)}|nK>
!!  gh0c1_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(0)}|\Psi^{(1)}>
!!      The wavefunction is orthogonal to the active space (for metals). It is not
!!      coherent with cg1.
!!  resid_k(nband_k)=residuals for each band over all k points,
!!  rhoaug1(cplex*n4,n5,n6)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output).
!!  ==== if (gs_hamkq%usepaw==1) ====
!!    cprj1(natom,nspinor*mband*mk1mem*nsppol*usecprj)=
!!              1st-order wave functions at k,q projected with non-local projectors:
!!                       cprj1=<p_i|C1nk,q> where p_i is a non-local projector
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!                                            (cumulative, so input as well as output)
!!
!! TODO
!!
!! PARENTS
!!      vtorho3
!!
!! CHILDREN
!!      accrho3,cgwf3,chkexi,corrmetalwf1,cprj_alloc,cprj_diskskip,cprj_free
!!      cprj_get,cprj_put,dotprod_g,getgsc,leave_new,matrixelmt_g,meanvalue_g
!!      sqnorm_g,status,timab,wffreaddatarec,wffreadnpwrec,wffreadskiprec
!!      wffwritedatarec,wffwritenpwrec,wrtout,xcomm_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vtowfk3(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,cpus,&
& dimcprj,dimekb,dime1kb,dim_eig2rf,dimffnlk,dimffnl1,dimphkxred,dkinpw,dtfil,dtset,&
& edocc_k,eeig0_k,eig0_k,eig0_kq,eig1_k,ekb_typ,e1kbfr,e1kbsc,&
& ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k,&
& fermie1,ffnlk,ffnlkq,ffnl1,gbound,gh0c1_set,gh1c_set,grad_berry,gs_hamkq,&
& ibg,ibgq,ibg1,icg,icgq,icg1,idir,ikpt,indlmn_typ,ipert,&
& isppol,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lmnmax,matblk,mband,mcgq,mcprjq,mgfft,mkmem,mk1mem,&
& mpi_enreg,mpsang,mpssoang,mpw,mpw1,natom,nband_k,ncpgr,&
& nkpg,nkpg1,nkpt,nnsclo_now,npw_k,npw1_k,nspinor,nsppol,&
& ntypat,n4,n5,n6,occ_k,pawrhoij1,ph3d,phkxred,prtvol,psps,resid_k,rhoaug1,rocceig,&
& sij_typ,usecprj,usee1kb,wffddk,wffnew,wffnow,wfftgs,vlocal,vlocal1,wtk_k)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vtowfk3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
 use interfaces_59_io_mpi
 use interfaces_66_wfs
 use interfaces_72_response, except_this_one => vtowfk3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,dime1kb,dimekb,dim_eig2rf
 integer,intent(in) :: dimffnl1,dimffnlk,dimphkxred,ibg
 integer,intent(in) :: ibg1,ibgq,icg,icg1,icgq,idir,ikpt,ipert,isppol
 integer,intent(in) :: lmnmax,matblk,mband,mcgq,mcprjq,mgfft,mk1mem,mkmem
 integer,intent(in) :: mpsang,mpssoang,mpw,mpw1,n4,n5,n6,natom,ncpgr,nkpg,nkpg1
 integer,intent(in) :: nkpt,nnsclo_now,nspinor,nsppol,ntypat,prtvol
 integer,intent(in) :: usecprj,usee1kb
 integer,intent(inout) :: nband_k,npw1_k,npw_k
 real(dp),intent(in) :: cpus,fermie1,wtk_k
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffddk,wffnew,wffnow,wfftgs
!arrays
 integer,intent(in) :: dimcprj(natom*gs_hamkq%usepaw),gbound(2*mgfft+8,2)
 integer,intent(in) :: indlmn_typ(6,lmnmax,1),kg1_k(3,npw1_k),kg_k(3,npw_k)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),cgq(2,mcgq)
 real(dp),intent(in) :: dkinpw(npw_k)
 real(dp),intent(in) :: e1kbfr(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
 real(dp),intent(in) :: e1kbsc(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
 real(dp),intent(in) :: eig0_k(nband_k),eig0_kq(nband_k)
 real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
 real(dp),intent(in) :: ffnl1(npw1_k,dimffnl1,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlk(npw_k,dimffnlk,lmnmax,1+gs_hamkq%usepaw*(ntypat-1))
 real(dp),intent(in) :: ffnlkq(npw1_k,dimffnl1,lmnmax,1)
 real(dp),intent(in) :: grad_berry(2,mpw1,nband_k),kinpw1(npw1_k)
 real(dp),intent(in) :: kpg1_k(npw1_k,nkpg1),kpg_k(npw_k,nkpg),kpt(3)
 real(dp),intent(in) :: occ_k(nband_k),phkxred(2,dimphkxred)
 real(dp),intent(in) :: rocceig(nband_k,nband_k)
 real(dp),intent(inout) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(inout) :: ph3d(2,npw1_k,matblk),rhoaug1(cplex*n4,n5,n6)
 real(dp),intent(inout) :: sij_typ(dimekb,gs_hamkq%usepaw),vlocal(n4,n5,n6)
 real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6)
 real(dp),intent(out) :: cg1_active(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh1c_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh0c1_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(out) :: edocc_k(nband_k),eeig0_k(nband_k),eig1_k(2*nband_k**2)
 real(dp),intent(out) :: ek0_k(nband_k),ek1_k(nband_k),eloc0_k(nband_k)
 real(dp),intent(out) :: enl0_k(nband_k),enl1_k(nband_k)
 real(dp),intent(out) :: resid_k(nband_k)
 type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
 type(cprj_type),intent(in) :: cprjq(natom,mcprjq)
 type(cprj_type),intent(out) :: cprj1(natom,nspinor*mband*mk1mem*nsppol*usecprj)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(natom*gs_hamkq%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=14
 integer,save :: nskip=0
! integer,save :: count=0 ! used in commented section below
 integer :: accesswff,counter,iband,ierr,iexit,igs,igscq,ii,dim_dcwf,inonsc
 integer :: iorder_cprj,iorder_cprj1,ipw,iscf_mod,ispinor,me,mgscq,nkpt_max
 integer :: nspinor0,openexit,option,opt_gvnl1,quit,spaceComm,test_ddk
 integer :: tim_fourwf,tocceig,usedcwavef
 real(dp) :: aa,ai,ar,eig0nk,resid
 real(dp) :: residk,scprod
 character(len=500) :: message
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwave0(:,:),cwave1(:,:),cwavef(:,:)
 real(dp),allocatable :: dcwavef(:,:),gh1c_n(:,:),gh0c1(:,:)
 real(dp),allocatable :: gsc(:,:),gscq(:,:),gvnl1(:,:),gvnlc(:,:)
 type(cprj_type),allocatable :: cwaveprj(:,:),cwaveprj0(:,:),cwaveprj1(:,:)
 type(cprj_type),allocatable :: cwaveprj_tmp(:,:)

! *********************************************************************

!Keep track of total time spent in vtowfk3
 call timab(128,1,tsec)

 nkpt_max=50
 if(mpi_enreg%paral_compil_kpt==1)nkpt_max=-1

!DEBUG
!write(std_out,*)' vtowfk3: enter '
!write(std_out,*)' vtowfk3: ikpt=',ikpt
!count=count+1
!write(std_out,*)' count=',count
!if(count==27)stop
!write(std_out,*)' vtowfk3: prtvol,wtk_k,npw_k,npw1_k,ipert'
!write(std_out,*)prtvol,wtk_k,npw_k,npw1_k,ipert
!if(ikpt==4)stop
!write(std_out,*)' vtowfk3 : cg1(:,1)=',cg1(:,1)
!write(std_out,*)' nband_k,natom,npw_k',nband_k,natom,npw_k
!stop
!ENDDEBUG

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,'vtowfk3 : enter'
   call wrtout(std_out,message,'PERS')
 end if

 quit=0
!Init me
 call xme_init(mpi_enreg,me)
!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
 accesswff=dtset%accesswff

 iscf_mod=dtset%iscf

!The value of iscf must be modified if ddk perturbation, see loper3.f
 if(ipert==natom+1) iscf_mod=-3

 ABI_ALLOCATE(gh0c1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gvnlc,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gvnl1,(2,npw1_k*nspinor))
 if (gs_hamkq%usepaw==1)  then
   ABI_ALLOCATE(gsc,(2,npw1_k*nspinor))
 end if

 if(prtvol>2 .or. ikpt<=nkpt_max)then
   write(message, '(a,a,i5,2x,a,3f9.5,2x,a)' ) ch10,&
&   ' Non-SCF iterations; k pt #',ikpt,'k=',kpt(:),'band residuals:'
   call wrtout(std_out,message,'PERS')
 end if

 ABI_ALLOCATE(cwave0,(2,npw_k*nspinor))
 ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))
 ABI_ALLOCATE(cwave1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gh1c_n,(2,npw1_k*nspinor))
!Read the npw and kg records of wf files
!NOTE : it should be possible to use rwwf in the present routine
 call status(0,dtfil%filstat,iexit,level,'before WffRead')
 test_ddk=0
 if( ipert==natom+2 .and. &
& sum( (dtset%qptn(1:3))**2 ) < 1.0d-7 .and. (dtset%berryopt .ne. 4) )then
   test_ddk=1
!  Read npw record
   call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor0,wffddk)
!  Skip k+G record
   call WffReadSkipRec(ierr,1,wffddk)
 end if
 if( ipert==natom+5 .and. &
& sum( (dtset%qptn(1:3))**2 ) < 1.0d-7 .and. (dtset%berryopt .ne. 4) )then
   test_ddk=1
!  Read npw record
!  JWZ, 20-Aug-08
!  original code: call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor,wffddk)
!  problem here is that nspinor is intent(in), but in WffReadNpwRec nspinor is intent(out).
!  I think this should read nspinor0, like the others
!  
   call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor0,wffddk)
!  
!  Skip k+G record
   call WffReadSkipRec(ierr,1,wffddk)
 end if
 if(mkmem==0)then
   call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_k,nspinor0,wfftgs)
!  Skip k+G and eigenvalue records in wfftgs (already in eigen0)
   call WffReadSkipRec(ierr,2,wfftgs)
 end if
 if(mk1mem==0)then
   call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor0,wffnow)
!  Skip k+G record
   call WffReadSkipRec(ierr,1,wffnow)
!  Initialize writing for this k point
   call WffWriteNpwRec(ierr,nband_k,npw1_k,nspinor,wffnew)
   ABI_ALLOCATE(kg_dum,(3,npw1_k))
   kg_dum(:,:) = kg1_k(:,:)
   call WffWriteDataRec(kg_dum,ierr,3,npw1_k,wffnew)
   ABI_DEALLOCATE(kg_dum)
 end if

!Additional stuff for PAW
 if (gs_hamkq%usepaw==1) then
!  1-Compute all <g|S|Cnk+q>
   igscq=0
   mgscq=mpw1*nspinor*mband
   ABI_ALLOCATE(gscq,(2,mgscq))
   call getgsc(cgq,cprjq,dimcprj,dimffnl1,ffnl1,gs_hamkq,gscq,ibgq,icgq,igscq,ikpt,isppol,&
&   kg1_k,psps%lmnmax,matblk,mcgq,mcprjq,mgfft,mgscq,mpi_enreg,psps%mpsang,&
&   psps%mpssoang,natom,nband_k,nkpt,npw1_k,dtset%nspinor,ntypat,ph3d)
!  2-Initialize additional scalars/arrays
   iorder_cprj=0;iorder_cprj1=0
   dim_dcwf=npw1_k*nspinor;if (ipert==natom+2.or.ipert==natom+5) dim_dcwf=0
   ABI_ALLOCATE(dcwavef,(2,dim_dcwf))
   if (usecprj==1) then
     ABI_ALLOCATE(cwaveprj0,(natom,nspinor))
     call cprj_alloc(cwaveprj0,1,dimcprj)
   end if
   ABI_ALLOCATE(cwaveprj,(natom,nspinor))
   ABI_ALLOCATE(cwaveprj1,(natom,nspinor))
   call cprj_alloc(cwaveprj ,0,dimcprj)
   call cprj_alloc(cwaveprj1,0,dimcprj)
 else
   igscq=0;mgscq=0;dim_dcwf=0
 end if

 call timab(139,1,tsec)

!Loop over bands

 do iband=1,nband_k

   if(mpi_enreg%paral_compil_kpt==1)then

     if( (mpi_enreg%proc_distrb(ikpt, iband,isppol) /= me )   ) then

       if(test_ddk==1)then
!        Skip the eigenvalue and the wf records of this band
         call WffReadSkipRec(ierr,2,wffddk)
       end if
       if(mkmem==0)then
         call WffReadSkipRec(ierr,1,wfftgs)
       end if
       if(mk1mem==0)then
         call WffReadSkipRec(ierr,2,wffnow)
!        Fill these records with zeroes (so that they can be read without I/O error)
         call WffWriteDataRec( (/ (zero*dble(ii),ii=1,2*nband_k) /) ,ierr,2*nband_k,wffnew)
         call WffWriteDataRec( (/ (zero*dble(ii),ii=1,2*npw1_k*nspinor) /) ,ierr,2*npw1_k*nspinor,wffnew)
       end if
!      Skip PAW projected WFs (cprj) or write zeros (cprj1)
       if (gs_hamkq%usepaw==1.and.usecprj==1) then
         if(mkmem==0)then
           call cprj_diskskip(mkmem,ncpgr,1,dtfil%unpaw)
         end if
         if(mk1mem==0)then
           ABI_ALLOCATE(cwaveprj_tmp,(natom,nspinor))
           call cprj_alloc(cwaveprj_tmp,0,dimcprj)
           call cprj_put(gs_hamkq%atindx,cwaveprj_tmp,cprj1,natom,iband,ibg1,ikpt,iorder_cprj1,isppol,&
&           mband,mk1mem,mpi_enreg,natom,1,nband_k,dimcprj,nspinor,nsppol,0,dtfil%unpaw1)
           call cprj_free(cwaveprj_tmp)
           ABI_DEALLOCATE(cwaveprj_tmp)
         end if
       end if
       cycle
     end if
   end if ! paral

!  Read ground-state wavefunctions
   if(mkmem/=0)then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(cg,cwave0,iband,icg,npw_k,nspinor)
     do ipw=1,npw_k*nspinor
       cwave0(1,ipw)=cg(1,ipw+(iband-1)*npw_k*nspinor+icg)
       cwave0(2,ipw)=cg(2,ipw+(iband-1)*npw_k*nspinor+icg)
     end do
!    $OMP END PARALLEL DO
   else
     call timab(288,1,tsec)
     call WffReadDataRec(cwave0,ierr,2,npw_k*nspinor,wfftgs)
     call timab(288,2,tsec)
   end if
!  Read PAW ground state projected WF (cprj)
   if (gs_hamkq%usepaw==1.and.usecprj==1) then
     call cprj_get(gs_hamkq%atindx1,cwaveprj0,cprj,natom,iband,ibg,ikpt,iorder_cprj,&
&     isppol,mband,mkmem,mpi_enreg,natom,1,nband_k,nspinor,nsppol,dtfil%unpaw,&
&     icpgr=idir,ncpgr=ncpgr)
   end if

!  Read first-order wavefunctions
   if(mk1mem/=0)then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(cg1,cwavef,iband,icg1,npw1_k,nspinor)
     do ipw=1,npw1_k*nspinor
       cwavef(1,ipw)=cg1(1,ipw+(iband-1)*npw1_k*nspinor+icg1)
       cwavef(2,ipw)=cg1(2,ipw+(iband-1)*npw1_k*nspinor+icg1)
     end do
!    $OMP END PARALLEL DO
   else
     call timab(288,1,tsec)
!    Skip the eigenvalue line
     call WffReadSkipRec(ierr,1,wffnow)
     call WffReadDataRec(cwavef,ierr,2,npw1_k*nspinor,wffnow)
     call timab(288,2,tsec)
   end if
!  Read PAW projected 1st-order WF (cprj)
!  Unuseful for the time being (will be recomputed in cgwf3)
!  if (gs_hamkq%usepaw==1.and.usecprj==1) then
!  call cprj_get(gs_hamkq%atindx1,cwaveprj,cprj1,natom,iband,ibg1,ikpt,iorder_cprj1,&
!  &      isppol,mband,mk1mem,mpi_enreg,natom,1,nband_k,nspinor,nsppol,dtfil%unpaw1)
!  end if

!  Filter the wavefunctions for large modified kinetic energy
!  The GS wavefunctions should already be non-zero
   do ispinor=1,nspinor
     igs=(ispinor-1)*npw1_k
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(cwavef,igs,kinpw1,npw1_k)
     do ipw=1+igs,npw1_k+igs
       if(kinpw1(ipw-igs)>huge(zero)*1.d-11)then
         cwavef(1,ipw)=zero
         cwavef(2,ipw)=zero
       end if
     end do
!    $OMP END PARALLEL DO
   end do

   if(prtvol>=10)then
     call status(0,dtfil%filstat,iexit,level,'after wf read ')
   end if

!  If electric field, the derivative of the wf should be read,
!  and multiplied by i.
   if(test_ddk==1)then
!    Skip the eigenvalue record
     call WffReadSkipRec(ierr,1,wffddk)
!    Read gvnl1
     call WffReadDataRec(gvnl1,ierr,2,npw1_k*nspinor,wffddk)
!    Multiplication by -i
!    MVeithen 021212 : use + i instead,
!    See X. Gonze, Phys. Rev. B 55, 10337 (1997) Eq. (79)
!    the operator used to compute the first-order derivative
!    of the wavefunctions with respect to an electric field
!    is $+i \frac{d}{dk}$
!    This change will affect the computation of the 2dtes from non
!    stationary expressions, see nstdy3.f and nstwf3.f

     do ipw=1,npw1_k*nspinor
!      aa=gvnl1(1,ipw)
!      gvnl1(1,ipw)=gvnl1(2,ipw)
!      gvnl1(2,ipw)=-aa
       aa=gvnl1(1,ipw)
       gvnl1(1,ipw)=-gvnl1(2,ipw)
       gvnl1(2,ipw)=aa
     end do
   end if

!  Unlike in GS calculations, the inonsc loop is inside the band loop
!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!  (often 1 for SCF calculation, =nstep for non-SCF calculations)
   do inonsc=1,nnsclo_now

     counter=100*iband+inonsc
!    Because in this loop, the CPU time matters, the writing
!    in the STATUS file is usually inhibited
     if(prtvol>=10)then
       call status(counter,dtfil%filstat,iexit,level,'loop iband    ')
     end if

!    Not too often, check whether the run must be stopped.
!    If so, iexit will be non-zero.
!    Note that when the number of bands becomes large, the check
!    must be done more often, because treating one band takes also longer ...
!    Only do this in the sequential mode
     if(mpi_enreg%paral_compil_kpt==0)then
       if(iband==1 .or. (nband_k>=16 .and. mod(iband,8)==1) &
&       .or. (nband_k>=32 .and. mod(iband,4)==1) &
&       .or. (nband_k>=64 .and. mod(iband,2)==1) &
&       .or. (nband_k>=128)                        )then
         openexit=1 ; if(dtset%chkexit<=1) openexit=0
         call chkexi(cpus,dtfil%filnam_ds(1),iexit,6,mpi_enreg,openexit)
         if(iexit/=0)quit=1
       end if
     end if

     if(prtvol>=10)then
       call status(counter,dtfil%filstat,iexit,level,'call cgwf3    ')
     end if

!    Note that the following translation occurs in the called routine :
!    iband->band, nband_k->nband, npw_k->npw, npw1_k->npw1
     eig0nk=eig0_k(iband)
     usedcwavef=gs_hamkq%usepaw;if (dim_dcwf==0) usedcwavef=0
     if (inonsc==1) usedcwavef=2*usedcwavef
     opt_gvnl1=0;if (ipert==natom+2) opt_gvnl1=1
     if (ipert==natom+2.and.gs_hamkq%usepaw==1.and.inonsc==1) opt_gvnl1=2
     call cgwf3(iband,dtset%berryopt,cgq,cplex,cwavef,cwave0,cwaveprj,cwaveprj0,dcwavef,&
&     dimcprj,dimekb,dime1kb,dimffnlk,dimffnl1,dimphkxred,dkinpw,eig0nk,eig0_kq,eig1_k,&
&     ekb_typ,e1kbfr,e1kbsc,ffnlk,ffnlkq,ffnl1,dtfil%filstat,gbound,gh0c1,gh1c_n,grad_berry,&
&     gsc,gscq,gs_hamkq,gvnlc,gvnl1,icgq,idir,indlmn_typ,ipert,igscq,&
&     kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lmnmax,matblk,mcgq,mgfft,mgscq,mpi_enreg,&
&     mpsang,mpssoang,mpw1,natom,nband_k,dtset%nbdbuf,nkpg,nkpg1,dtset%nline,&
&     npw_k,npw1_k,nspinor,ntypat,n4,n5,n6,opt_gvnl1,dtset%ortalg,dtset%paral_kgb,ph3d,phkxred,prtvol,&
&     quit,resid,dtset%sciss,sij_typ,dtset%tolwfr,usecprj,usedcwavef,usee1kb,vlocal,&
&     vlocal1,dtset%wfoptalg)

     resid_k(iband)=resid
     if(prtvol>=10)then
       call status(counter,dtfil%filstat,iexit,level,'after cgwf    ')
     end if

!    At this stage, the 1st order function cwavef is orthogonal to cgq (unlike
!    when it is input to cgwf3). Here, restore the "active space" content
!    of the first-order wavefunction, to give cwave1.
     call corrmetalwf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,edocc_k,eig1_k,fermie1,gh0c1,&
&     iband,ibgq,icgq,gs_hamkq%istwf_k,mcgq,mcprjq,mpi_enreg,natom,nband_k,npw1_k,nspinor,&
&     occ_k,rocceig,0,gs_hamkq%usepaw,tocceig)

     if ( abs(occ_k(iband)) <= tol8 ) then

       ek0_k(iband)=zero
       ek1_k(iband)=zero
       eeig0_k(iband)=zero
       enl0_k(iband)=zero
       enl1_k(iband)=zero
       eloc0_k(iband)=zero
       nskip=nskip+1

     else

!      Compute the 0-order kinetic operator contribution (with cwavef)
       call meanvalue_g(ar,kinpw1,0,gs_hamkq%istwf_k,mpi_enreg,npw1_k,nspinor,cwavef,cwavef,0)
!      There is an additional factor of 2 with respect to the bare matrix element
       ek0_k(iband)=two*ar

!      Compute the 1-order kinetic operator contribution (with cwave1 and cwave0), if needed.
!      Note that this is called only for ddk or strain, so that npw1_k=npw_k
       if(ipert==natom+1 .or. ipert==natom+3 .or. ipert==natom+4)then
         call matrixelmt_g(ai,ar,dkinpw,gs_hamkq%istwf_k,mpi_enreg,0,npw_k,nspinor,cwave1,cwave0)
!        There is an additional factor of 4 with respect to the bare matrix element
         ek1_k(iband)=four*ar
       end if

!      Compute eigenvalue part of total energy (with cwavef)
       if (gs_hamkq%usepaw==1) then
         call dotprod_g(scprod,ai,gs_hamkq%istwf_k,mpi_enreg,npw1_k*nspinor,1,cwavef,gsc)
       else
         call sqnorm_g(scprod,gs_hamkq%istwf_k,mpi_enreg,npw1_k*nspinor,cwavef)
       end if
       eeig0_k(iband)=-two*(eig0_k(iband)- (dtset%sciss) )*scprod

!      Compute nonlocal psp contributions to nonlocal energy:
!      <G|Vnl|C1nk(perp)> is contained in gvnlc (with cwavef)
       call dotprod_g(scprod,ai,gs_hamkq%istwf_k,mpi_enreg,npw1_k*nspinor,1,cwavef,gvnlc)
       enl0_k(iband)=two*scprod

!      <G|Vnl1|Cnk> is contained in gvnl1 (with cwave1)
       call dotprod_g(scprod,ai,gs_hamkq%istwf_k,mpi_enreg,npw1_k*nspinor,1,cwave1,gvnl1)
       enl1_k(iband)=four*scprod

!      Removal of the 1st-order kinetic energy from the 1st-order non-local part.
       if(ipert==natom+1 .or. &
&       ipert==natom+3 .or. ipert==natom+4) then
         enl1_k(iband)=enl1_k(iband)-ek1_k(iband)
       end if

!      Compute the 1st-order valence space contribution
!      Note from MT to JZ: I displaced this into nstpaw3

!      Accumulate 1st-order density (only at the last inonsc)
!      Accumulate zero-order potential part of the 2nd-order total energy
       tim_fourwf=5;option=2;if (iscf_mod>0.and.inonsc==nnsclo_now) option=3
       call accrho3(counter,cplex,cwave0,cwave1,cwavef,cwaveprj0,cwaveprj1,&
&       dimcprj,dimffnlk,ffnlk,dimphkxred,eloc0_k(iband),dtfil%filstat,gbound,&
&       gs_hamkq,iband,idir,ipert,isppol,kg_k,kg1_k,kpg_k,kpt,dtset%kptopt,lmnmax,matblk,mgfft,&
&       mpi_enreg,natom,nband_k,ncpgr,nkpg,npw_k,npw1_k,nspinor,ntypat,n4,n5,n6,occ_k,option,&
&       dtset%paral_kgb,pawrhoij1,ph3d,phkxred,prtvol,rhoaug1,tim_fourwf,usecprj,vlocal,tocceig,wtk_k)

!      End of non-zero occupation
     end if

!    Exit loop over inonsc if converged and if non-self-consistent
     if (iscf_mod<0 .and. resid<dtset%tolwfr) exit

!    End loop over inonsc
   end do

!  Write first-order eigenvalues and wavefunctions
   if(mk1mem/=0)then
     cg1(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)=cwave1(:,:)
     if(dim_eig2rf > 0) then
       cg1_active(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)=cwavef(:,:)
       gh1c_set(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)=gh1c_n(:,:)
       gh0c1_set(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)=gh0c1(:,:)
     end if
   else
     call timab(288,1,tsec)
     call WffWriteDataRec(eig1_k,ierr,2*nband_k,wffnew)
     call WffWriteDataRec(cwave1,ierr,2,npw1_k*nspinor,wffnew)
     call timab(288,2,tsec)
   end if

!  PAW: write first-order projected wavefunctions
   if (psps%usepaw==1.and.usecprj==1) then
     call cprj_put(gs_hamkq%atindx,cwaveprj,cprj1,natom,iband,ibg1,ikpt,iorder_cprj1,isppol,&
&     mband,mk1mem,mpi_enreg,natom,1,nband_k,dimcprj,nspinor,nsppol,0,dtfil%unpaw1)
   end if

   if(prtvol>=10)then
     call status(counter,dtfil%filstat,iexit,level,'get residk    ')
   end if

!  End loop over bands
 end do

!Find largest resid over bands at this k point
 residk=maxval(resid_k(:))
 if(prtvol>2 .or. ikpt<=nkpt_max)then
   do ii=0,(nband_k-1)/8
     write(message, '(1p,8e10.2)' ) &
&     (resid_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
     call wrtout(std_out,message,'PERS')
   end do
 end if

 call timab(139,2,tsec)
 call timab(130,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'after loops   ')

 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(cwave1)
 ABI_DEALLOCATE(gh0c1)
 ABI_DEALLOCATE(gvnlc)
 ABI_DEALLOCATE(gvnl1)
 ABI_DEALLOCATE(gh1c_n)
 if (gs_hamkq%usepaw==1) then
   ABI_DEALLOCATE(dcwavef)
   ABI_DEALLOCATE(gsc)
   ABI_DEALLOCATE(gscq)
   call cprj_free(cwaveprj)
   call cprj_free(cwaveprj1)
   ABI_DEALLOCATE(cwaveprj)
   ABI_DEALLOCATE(cwaveprj1)
   if (usecprj==1) then
     call cprj_free(cwaveprj0)
     ABI_DEALLOCATE(cwaveprj0)
   end if
 end if

!###################################################################

!DEBUG
!write(std_out,*)'vtowfk3: iscf_mod, nband_k',iscf_mod, nband_k
!ENDDEBUG

!Write the number of one-way 3D ffts skipped until now (in case of fixed
!occupation numbers
 if(iscf_mod>0 .and. (prtvol>2 .or. ikpt<=nkpt_max))then
   write(message, '(a,i8)' )&
&   ' vtowfk3 : number of one-way 3D ffts skipped in vtowfk until now =',nskip
   call wrtout(std_out,message,'PERS')
 end if

 if(prtvol<=2 .and. ikpt==nkpt_max+1)then
   write(message, '(a,a,a)' ) ch10,&
&   ' vtowfk3 : prtvol=0, 1 or 2, do not print more k-points.',ch10
   call wrtout(std_out,message,'PERS')
 end if

!###################################################################

 if (residk>dtset%tolwfr .and. iscf_mod<=0 .and. iscf_mod/=-3) then
   write(message, '(a,a,a,a,2i5,a,es13.5)' ) ch10,&
&   ' vtowfk3: WARNING -',ch10,&
&   '  Wavefunctions not converged for nnsclo,ikpt=',nnsclo_now,ikpt,&
&   ' max resid=',residk
   call wrtout(std_out,message,'PERS')
 end if

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
   write(message,'(a1,a,a1,a,i2,a)') ch10,&
&   ' vtowfk3 : exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(130,2,tsec)
 call timab(128,2,tsec)

!DEBUG
!write(std_out,*)' vtowfk3 : exit '
!call flush(6)
!if(count==26)stop
!stop
!ENDDEBUG

end subroutine vtowfk3
!!***
