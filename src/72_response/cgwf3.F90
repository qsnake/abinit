!{\src2tex{textfont=tt}}
!!****f* ABINIT/cgwf3
!! NAME
!! cgwf3
!!
!!
!! FUNCTION
!! Update one single wavefunction (cwavef), non self-consistently.
!! Uses a conjugate-gradient algorithm.
!! Try to keep close to the formulas in PRB55, 10337 (1997), for the
!! non-self-consistent case, except that we are computing here
!! the second-derivative of the total energy, and not E(2). There
!! is a factor of 2 between the two quantities ...
!! The wavefunction that is generated is always orthogonal to cgq .
!! It is orthogonal to the active Hilbert space, and will be complemented
!! by contributions from the active space in the calling routine, if needed.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,DRH,XW,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  band=which particular band we are converging.
!!  berryopt=option for Berry phase
!!  cgq(2,mcgq)=wavefunction coefficients for ALL bands at k+Q
!!  cplex=1 if vlocal1 is real, 2 if vlocal1 is complex
!!  cwave0(2,npw*nspinor)=GS wavefunction at k, in reciprocal space
!!  cwaveprj0(natom,nspinor*usecprj)=GS wave function at k projected with nl projectors
!!  dimcprj(natom*usepaw)=array of dimensions of arrays cprj, cprjq (ordered by atom-type)
!!  dimekb=first dimension of ekb (see ekb_typ below)
!!  dime1kb=first dimension of e1kb (see e1kb below)
!!  dimffnlk=second dimension of ffnl (1+number of derivatives)
!!  dimffnl1=second dimension of ffnl1 and ffnlkq (1+number of derivatives)
!!  dimphkxred=second dimension of phkxred
!!  dkinpw(npw)=derivative of the (modified) kinetic energy for each plane wave at k (Hartree)
!!  eig0nk=0-order eigenvalue for the present wavefunction at k
!!  eig0_kq(nband)=GS eigenvalues at k+Q (hartree)
!!  ekb_typ(dimekb,1,nspinor**2)=
!!  ->Norm conserving : (Real) Kleinman-Bylander energies (hartree) for the displaced atom
!!          for number of basis functions (l,n) (lnmax)
!!          dimekb=lnmax
!!    ->PAW : (Real, symmetric) Frozen part of Dij coefficients
!!                               to connect projectors for the displaced atom
!!          for number of basis functions (l,m,n) (lmnmax)
!!          dimekb=lmnmax*(lmnmax+1)/2
!!          These are complex numbers in particular cases (spin-orbit)
!!               ekb_typ(:,:,1) contains Dij^up-up
!!               ekb_typ(:,:,2) contains Dij^dn-dn
!!               ekb_typ(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!               ekb_typ(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  e1kbfr(dime1kb,dimekb2,use1ekb*nspinor**2)=frozen part of 1st der. of ekb
!                     for the current pertubation (not depending on VHxc^(1))
!!  e1kbsc(dime1kb,dimekb2,use1ekb*nspinor**2)=self-consistent part of 1st der. of ekb
!!                                 for the current pertubation (depending on VHxc^(1))
!!  ffnlk(npw,dimffnlk,lmnmax,1+usepaw*(ntypat-1))=nonloc form factors at k, for the displaced atom.
!!  ffnlkq(npw1,dimffnl1,lmnmax,1)=nonloc form fact at k+Q for the displaced atom
!!  ffnl1(npw1,dimffnl1,lmnmax,ntypat)=nonloc form factors at k+Q
!!  filstat=name of the status file
!!  grad_berry(2,mpw1,dtefield%nband_occ) = the gradient of the Berry phase term
!!  gbound(2*mgfft+8,2)=G sphere boundary at k
!!  gscq(2,mgscq)=<g|S|Cnk+q> coefficients for ALL bands (PAW) at k+Q
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+Q
!!  icgq=shift to be applied on the location of data in the array cgq
!!  igscq=shift to be applied on the location of data in the array gscq
!!  idir=direction of the perturbation
!!  indlmn_typ(6,lmnmax,1)=indlmn info for the displaced atom
!!  ipert=type of the perturbation
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere at k.
!!  kg1_k(3,npw1)=coordinates of planewaves in basis sphere at k+Q.
!!  kinpw1(npw1)=(modified) kinetic energy for each plane wave at k+Q (Hartree)
!!  kpg_k(npw,nkpg)= (k+G) components at k (only if useylm=1)
!!  kpg1_k(npw1,nkpg1)= (k+G) components at k+Q (only if useylm=1)
!!  kpt(3)=coordinates of k point.
!!  lmnmax= max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mcgq=second dimension of the cgq array
!!  mgfft=maximum size of 1D FFTs
!!  mgscq=second dimension of gscq
!!  mpi_enreg=informations about MPI parallelization
!!  mpw1=maximum number of planewave for first-order wavefunctions
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband=number of bands.
!!  nbdbuf=number of buffer bands for the minimisation
!!  nkpg,nkpg1=second dimensions of kpg_k and kpg1_k (0 if useylm=0)
!!  nline=number of line minimizations per band.
!!  npw=number of planewaves in basis sphere at given k.
!!  npw1=number of planewaves in basis sphere at k+Q
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell.
!!  n4,n5,n6 dimensions of vlocal and vlocal1
!!  opt_gvnl1=option controlling the use of gvnl1 array:
!!            0: used as an output
!!            1: used as an input:    - used only for ipert=natom+2
!!                 NCPP: contains the ddk 1-st order WF
!!                 PAW: contains frozen part of 1st-order hamiltonian
!!            2: used as input/ouput:    - used only for PAW and ipert=natom+2
!!                 At input: contains the ddk 1-st order WF (times i)
!!                 At output: contains frozen part of 1st-order hamiltonian
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  ph3d(2,npw1,matblk)=3-dim structure factors, for each atom and plane wave.
!!  phkxred(2,dimphkxred)=phase factors exp(2 pi kpoint.xred) at k
!!  prtvol=control print volume and debugging output
!!  quit= if 1, proceeds to smooth ending of the job.
!!  sciss=scissor shift (Ha)
!!  sij_typ(dimekb,gs_hamkq%usepaw)=overlap matrix components (PAW)
!!  tolwfr=tolerance on largest wf residual
!!  usecprj= 1 if cwaveprj0 array is stored in memory
!!  usee1kb=1 if ekb derivatives (e1kbsc, e1kbfr) exist
!!  usedcwavef=flag controlling the use of dcwavef array (PAW only):
!!             0: not used (not allocated)
!!             1: used as input
!!             2: used as output
!!  vlocal(n4,n5,n6)= GS local pot in real space, on the augmented fft grid
!!  vlocal1(cplex*n4,n5,n6)= RF local pot in real space, on the augmented fft grid
!!  wfoptalg=govern the choice of algorithm for wf optimisation (0 or 10, at present)
!!
!! OUTPUT
!!  eig1_k(2*nband**2)=matrix of first-order eigenvalues (hartree)
!!                     eig1(:,ii,jj)=<C0 ii|H1|C0 jj> for norm-conserving psps
!!                     eig1(:,ii,jj)=<C0 ii|H1-(eig0_k+eig0_k+q)/2.S(1)|C0 jj> for PAW
!!  ghc(2,npw1*nspinor)=<G|H0-eig0_k.I|C1 band,k> (NCPP) or <G|H0-eig0_k.S0|C1 band,k> (PAW)
!!  gvnlc(2,npw1*nspinor)=<G|Vnl|C1 band,k>
!!  gvnl1(2,npw1*nspinor)=  part of <G|K1+Vnl1|C0 band,k> not depending on VHxc1           (NCPP)
!!                       or part of <G|K1+Vnl1-eig0k.S1|C0 band,k> not depending on VHxc1 (PAW)
!!  resid=wf residual for current band
!!  gh1c_n= <G|H1|C0 band,k> (NCPP) or <G|H1-eig0k.S1|C0 band,k> (PAW). This vector is not projected
!!     on the subspace orhtogonal to the cg.
!!  === if gs_hamkq%usepaw==1 ===
!!  gsc(2,npw1*nspinor*usepaw)=<G|S0|C1 band,k>
!!
!! SIDE EFFECTS
!!  Input/Output:
!!  cwavef(2,npw1*nspinor)=first-order  wavefunction at k,q, in reciprocal space (updated)
!!  === if gs_hamkq%usepaw==1 ===
!!  cwaveprj(natom,nspinor)= wave functions at k projected with nl projectors
!!  === if usedcwavef>0 ===
!!  dcwavef(2,npw1*nspinor)=change of wavefunction due to change of overlap:
!!         dcwavef is delta_Psi(1)=-1/2.Sum_{j}[<C0_k+q_j|S(1)|C0_k_i>.|C0_k+q_j>]
!!         see PRB 78, 035105 (2008), Eq. (42)
!!         input if usedcwavef=1, output if usedcwavef=2
!!
!! PARENTS
!!      vtowfk3
!!
!! CHILDREN
!!      cprj_alloc,cprj_axpby,cprj_free,cprj_set_zero,dotprod_g,getdc1,getgh1c
!!      getghc,leave_new,precon,projbd,sqnorm_g,status,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cgwf3(band,berryopt,cgq,cplex,cwavef,cwave0,cwaveprj,cwaveprj0,dcwavef,&
& dimcprj,dimekb,dime1kb,dimffnlk,dimffnl1,dimphkxred,dkinpw,eig0nk,eig0_kq,eig1_k,&
& ekb_typ,e1kbfr,e1kbsc,ffnlk,ffnlkq,ffnl1,filstat,gbound,ghc,gh1c_n,grad_berry,gsc,gscq,&
& gs_hamkq,gvnlc,gvnl1,icgq,idir,indlmn_typ,&
& ipert,igscq,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,&
& kpt,lmnmax,matblk,mcgq,mgfft,mgscq,mpi_enreg,mpsang,mpssoang,mpw1,natom,nband,nbdbuf,&
& nkpg,nkpg1,nline,npw,npw1,nspinor,ntypat,n4,n5,n6,opt_gvnl1,ortalg,paral_kgb,ph3d,phkxred,prtvol,&
& quit,resid,sciss,sij_typ,tolwfr,usecprj,usedcwavef,usee1kb,vlocal,vlocal1,wfoptalg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgwf3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_53_spacepar
 use interfaces_66_wfs
 use interfaces_72_response, except_this_one => cgwf3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,berryopt,cplex,dime1kb,dimekb,dimffnl1,dimffnlk
 integer,intent(in) :: dimphkxred,icgq,idir,igscq,ipert,lmnmax,matblk,mcgq
 integer,intent(in) :: mgfft,mgscq,mpsang,mpssoang,mpw1,n4,n5,n6,natom,nband
 integer,intent(in) :: nbdbuf,nkpg,nkpg1,nline,npw,npw1,nspinor,ntypat
 integer,intent(in) :: opt_gvnl1,ortalg,paral_kgb,prtvol,quit,usecprj,usedcwavef,usee1kb
 integer,intent(in) :: wfoptalg
 real(dp),intent(in) :: eig0nk,sciss,tolwfr
 real(dp),intent(out) :: resid
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
!arrays
 integer,intent(in) :: dimcprj(natom*gs_hamkq%usepaw),gbound(2*mgfft+8,2)
 integer,intent(in) :: indlmn_typ(6,lmnmax,1),kg1_k(3,npw1),kg_k(3,npw)
 real(dp),intent(in) :: cgq(2,mcgq),dkinpw(npw)
 real(dp),intent(in) :: e1kbfr(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
 real(dp),intent(in) :: e1kbsc(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
 real(dp),intent(in) :: eig0_kq(nband),ekb_typ(dimekb,1,nspinor**2)
 real(dp),intent(in) :: ffnl1(npw1,dimffnl1,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlk(npw,dimffnlk,lmnmax,1+gs_hamkq%usepaw*(ntypat-1))
 real(dp),intent(in) :: ffnlkq(npw1,dimffnl1,lmnmax,1)
 real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband),gscq(2,mgscq)
 real(dp),intent(in) :: kinpw1(npw1),kpg1_k(npw1,nkpg1),kpg_k(npw,nkpg),kpt(3)
 real(dp),intent(in) :: phkxred(2,dimphkxred)
 real(dp),intent(in) :: sij_typ(dimekb,gs_hamkq%usepaw)
 real(dp),intent(inout) :: cwave0(2,npw*nspinor),cwavef(2,npw1*nspinor)
 real(dp),intent(inout) :: dcwavef(2,npw1*nspinor*((usedcwavef+1)/2))
 real(dp),intent(inout) :: ph3d(2,npw1,matblk),vlocal(n4,n5,n6)
 real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6)
 real(dp),intent(out) :: eig1_k(2*nband**2),gh1c_n(2,npw1*nspinor)
 real(dp),intent(out) :: ghc(2,npw1*nspinor)
 real(dp),intent(out) :: gsc(2,npw1*nspinor*gs_hamkq%usepaw)
 real(dp),intent(out) :: gvnl1(2,npw1*nspinor),gvnlc(2,npw1*nspinor)
 type(cprj_type),intent(inout) :: cwaveprj(natom,nspinor)
 type(cprj_type),intent(inout) :: cwaveprj0(natom,nspinor*usecprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=15,tim_getgh1c=1,tim_getghc=2
 integer,save :: nskip=0
 integer :: cpopt,iband,iexit,igs,ii,iline,indx_cgq,iprint,ipw
 integer :: ipws,ispinor,istwf_k,jband,nvloc,optlocal,optnl,sij_opt
 integer :: tim_projbd,useoverlap,usepaw,usevnl
 real(dp) :: cgwftol,d2edt2,d2te,d2teold,dedt,deltae,deold,dotgg
 real(dp) :: dotgp,doti,dotr,eshift,eshiftkq,gamma,optekin,prod1,prod2
 real(dp) :: theta,u1h0me0u1
 logical :: gen_eigenpb
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: conjgr(:,:),cwaveq(:,:),cwwork(:,:),direc(:,:),dummy(:,:)
 real(dp),allocatable :: gberry(:,:),gh1c(:,:),gh_direc(:,:),gresid(:,:)
 real(dp),allocatable :: gs1c(:,:),gvnl_direc(:,:),pcon(:),sconjgr(:,:)
 real(dp),allocatable :: scprod(:,:),wfraug(:,:,:,:),work(:,:),work1(:,:),work2(:,:)
 type(cprj_type),allocatable :: conjgrprj(:,:),cprj_dummy(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!======================================================================
!========= LOCAL VARIABLES DEFINITIONS AND ALLOCATIONS ================
!======================================================================

 call timab(122,1,tsec)
 if(prtvol<0) call status(0,filstat,iexit,level,'enter         ')
 iprint=0;if(prtvol==-level)iprint=1

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' cgwf3 : enter '
   call wrtout(std_out,msg,'PERS')
 end if

!Tell us what is going on:
 if(prtvol>=10)then
   write(msg, '(a,i6,2x,a,i3,a)' ) &
&   ' --- cgwf3 is called for band',band,'for',nline,' lines'
   call  wrtout(std_out,msg,'PERS')
 end if

!if PAW, one has to solve a generalized eigenproblem
 usepaw=gs_hamkq%usepaw
 gen_eigenpb=(usepaw==1)
 useoverlap=0;if (gen_eigenpb) useoverlap=1

!Use scissor shift on 0-order eigenvalue
 eshift=eig0nk-sciss

!Additional initializations
 nvloc=1 !For the time being, non-collinear potential is not allowed in RF
 istwf_k=gs_hamkq%istwf_k
 optekin=0;if (wfoptalg>=10) optekin=1
 tim_projbd=2

!Memory allocations
 ABI_ALLOCATE(gh1c,(2,npw1*nspinor))
 ABI_ALLOCATE(pcon,(npw1))
 ABI_ALLOCATE(scprod,(2,nband))
 ABI_ALLOCATE(wfraug,(2,n4,n5,n6))
 if (berryopt==4) then
   ABI_ALLOCATE(gberry,(2,npw1*nspinor))
   gberry(:,1:npw1*nspinor)=grad_berry(:,1:npw1*nspinor,band)
 end if

!Several checking statements
 if (prtvol==-level) then
   ABI_ALLOCATE(work,(2,npw1*nspinor))
   ABI_ALLOCATE(work1,(2,npw1*nspinor))
!  ===== Check <Psi_k+q^(0)|S(0)|Psi_k+q^(0)>=delta
   if (.not.gen_eigenpb) work1(:,:)=cgq(:,1+npw1*nspinor*(band-1)+icgq:npw1*nspinor*band+icgq)
   if (     gen_eigenpb) work1(:,:)=gscq(:,1+npw1*nspinor*(band-1)+igscq:npw1*nspinor*band+igscq)
   do iband=1,nband
     work(:,:)=cgq(:,1+npw1*nspinor*(iband-1)+icgq:npw1*nspinor*iband+icgq)
     call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1*nspinor,2,work1,work)
     write(msg,'(a,i3,a,2es22.15)') "<Psi_k+q,i^(0)|S(0)|Psi_k+q,j^(0)> for band j=",iband," is ",dotr,doti
     call wrtout(std_out,msg,'PERS')
   end do
!  ===== Check Pc.Psi_k+q^(0)=0
!  allocate(work(2,npw1*nspinor))
   do iband=1,nband
     work(:,:)=cgq(:,1+npw1*nspinor*(iband-1)+icgq:npw1*nspinor*iband+icgq)
     call projbd(cgq,work,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&     ortalg,iprint,gscq,scprod,0,tim_projbd,useoverlap)
     call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,work)
     write(msg,'(a,i3,a,es22.15)') "Norm of Pc.Psi_k+q_j^(0) for band j=",iband," is ",dotr
     call wrtout(std_out,msg,'PERS')
   end do
!  ===== Check Pc.Psi_k^(0)=0
   work(:,:)=cwave0(:,:)
   call projbd(cgq,work,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&   ortalg,iprint,gscq,scprod,0,tim_projbd,useoverlap)
   call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,work)
   write(msg,'(a,i3,a,es22.15)') "Norm of Pc.Psi_k^(0) for band ",band," is ",dotr
   call wrtout(std_out,msg,'PERS')
!  ===== Check Pc^*.S(0).Psi_k+q^(0)=0
   if (gen_eigenpb) then
     do iband=1,nband
       work(:,:)=gscq(:,1+npw1*nspinor*(iband-1)+igscq:npw1*nspinor*iband+igscq)
       call projbd(gscq,work,-1,igscq,icgq,istwf_k,mgscq,mpi_enreg,mcgq,nband,npw1,nspinor,&
&       ortalg,iprint,cgq,scprod,0,tim_projbd,useoverlap)
       call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,work)
       write(msg,'(a,i3,a,es22.15)') "Norm of Pc^*.S(0).Psi_k+q_j^(0) for band j=",iband," is ",dotr
       call wrtout(std_out,msg,'PERS')
     end do
   end if
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(work1)
 end if

!======================================================================
!========== INITIALISATION OF MINIMIZATION ITERATIONS =================
!======================================================================

!Compute H(1) applied to GS wavefunction Psi(0)
 if (gen_eigenpb) then
   sij_opt=1
   ABI_ALLOCATE(gs1c,(2,npw1*nspinor))
 else
   sij_opt=0
 end if
 usevnl=1 ; optlocal=1 ; optnl=2
 call getgh1c(berryopt,cplex,cwave0,cwaveprj0,dimcprj,dimekb,dime1kb,dimffnlk,dimffnl1,dimffnl1,dimphkxred,dkinpw,&
& ekb_typ,e1kbfr,e1kbsc,ffnlk,ffnlkq,ffnl1,filstat,gbound,gh1c,gberry,gs1c,gs_hamkq,gvnl1,idir,indlmn_typ,&
& ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,eshift,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
& nkpg,nkpg1,npw,npw1,nspinor,ntypat,n4,n5,n6,optlocal,optnl,opt_gvnl1,paral_kgb,ph3d,phkxred,prtvol,&
& sij_opt,sij_typ,tim_getgh1c,usecprj,usee1kb,usevnl,vlocal1,wfraug)

 if (gen_eigenpb) then
   if (ipert/=natom+2) then  ! S^(1) is zero for ipert=natom+2
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(eshift,gh1c,gvnl1,gs1c,npw1,nspinor)
     do ipw=1,npw1*nspinor
       gh1c (1:2,ipw)=gh1c (1:2,ipw)-eshift*gs1c(1:2,ipw)
       gvnl1(1:2,ipw)=gvnl1(1:2,ipw)-eshift*gs1c(1:2,ipw)
     end do
!    $OMP END PARALLEL DO
   end if

!  If generalized eigenPb and dcwavef requested, compute it:
!  dcwavef is delta_Psi(1)=-1/2.Sum_{j}[<C0_k+q_j|S(1)|C0_k_i>.|C0_k+q_j>]
!  see PRB 78, 035105 (2008), Eq. (42)
   if (usedcwavef==2) then
     call getdc1(cgq,cprj_dummy,dcwavef,cprj_dummy,0,icgq,istwf_k,mcgq,0,&
&     mpi_enreg,natom,nband,npw1,nspinor,0,ortalg,gs1c)
   end if
 end if

!Check that Pc^*.(H^(0)-E.S^(0)).delta_Psi^(1) is zero
 if (prtvol==-level.and.usedcwavef==2) then
   ABI_ALLOCATE(cwwork,(2,npw1*nspinor))
   cwwork=dcwavef
!  - Apply H^(0)-E.S^(0) to delta_Psi^(1)
   sij_opt=0;if (gen_eigenpb) sij_opt=-1
   cpopt=-1 ! no use of cprj !!!
   ABI_ALLOCATE(work,(2,npw1*nspinor))
   ABI_ALLOCATE(work1,(2,npw1*nspinor*((sij_opt+1)/2)))
   ABI_ALLOCATE(work2,(2,npw1*nspinor))
   call getghc(cpopt,cwwork,conjgrprj,dimffnl1,ffnl1,filstat,&
&   work,work1,gs_hamkq,work2,kg1_k,kinpw1,eshift,lmnmax,matblk,&
&   mgfft,mpi_enreg,mpsang,mpssoang,natom,1,npw1,nspinor,ntypat,nvloc,n4,n5,n6,&
&   paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
   cwwork=work
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
!  -Apply Pc^*
   call projbd(gscq,cwwork,-1,igscq,icgq,istwf_k,mgscq,mpi_enreg,mcgq,nband,npw1,nspinor,&
&   ortalg,iprint,cgq,scprod,0,tim_projbd,useoverlap)
   call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,cwwork)
   ABI_DEALLOCATE(cwwork)
   write(msg,'(a,i3,a,es22.15)') '|Pc^*.(H^(0)-E.S^(0)).delta_Psi^(1)|^2 (band ',band,')=',dotr
   call wrtout(std_out,msg,'PERS')
 end if

!$OMP PARALLEL DO PRIVATE(ipw) &
!$OMP&SHARED(gh1c,gh1c_n,npw1,nspinor)
 do ipw=1,npw1*nspinor
   gh1c_n(1:2,ipw)=gh1c(1:2,ipw)
 end do
!$OMP END PARALLEL DO

!Projecting out all bands (this could be avoided)
!Note the subtlety:
!-For the generalized eigenPb, S|cgq> is used in place of |cgq>,
!in order to apply P_c+ projector (see PRB 73, 235101 (2006), Eq. (71), (72))
 if(gen_eigenpb)then
   call projbd(gscq,gh1c,-1,igscq,icgq,istwf_k,mgscq,mpi_enreg,mcgq,nband,npw1,nspinor,&
&   ortalg,iprint,cgq,scprod,0,tim_projbd,useoverlap)
 else
   call projbd(cgq,gh1c,-1,icgq,0,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&   ortalg,iprint,dummy,scprod,0,tim_projbd,useoverlap)
 end if

!The array eig1_k contains:
!<u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
!or <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
 jband=(band-1)*2*nband
 if (gen_eigenpb) then
   indx_cgq=icgq
   do iband=1,nband
     eshiftkq=half*(eig0_kq(iband)-eig0nk)
     call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1*nspinor,2,cgq(:,indx_cgq+1:indx_cgq+npw1*nspinor),gs1c)
     eig1_k(2*iband-1+jband)=scprod(1,iband)-eshiftkq*dotr
     eig1_k(2*iband  +jband)=scprod(2,iband)-eshiftkq*doti
     indx_cgq=indx_cgq+npw1*nspinor
   end do
 else
   do iband=1,nband
     eig1_k(2*iband-1+jband)=scprod(1,iband)
     eig1_k(2*iband  +jband)=scprod(2,iband)
   end do
 end if

!No more need of gs1c
 if (gen_eigenpb)  then
   ABI_DEALLOCATE(gs1c)
 end if

!Filter the wavefunctions for large modified kinetic energy (see routine mkkin.f)
 do ispinor=1,nspinor
   ipws=(ispinor-1)*npw1
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cwavef,kinpw1,ipws,npw1)
   do ipw=1+ipws,npw1+ipws
     if(kinpw1(ipw-ipws)>huge(zero)*1.d-11)then
       cwavef(1:2,ipw)=zero
     end if
   end do
!  $OMP END PARALLEL DO
 end do

!Apply the orthogonality condition: <C1 k,q|C0 k+q>=0
!Project out all bands from cwavef, i.e. apply P_c projector on cwavef
!(this is needed when there are some partially or unoccupied states)
 call projbd(cgq,cwavef,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
& ortalg,iprint,gscq,scprod,0,tim_projbd,useoverlap)

!If PAW, the orthogonality condition is
!<C1 k,q|S0|C0 k+q>+1/2<C0 k|S1|C0 k+q>=0
 if (usepaw==1.and.usedcwavef>0) then
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cwavef,dcwavef,npw1,nspinor)
   do ipw=1,npw1*nspinor
     cwavef(1:2,ipw)=cwavef(1:2,ipw)+dcwavef(1:2,ipw)
   end do
!  $OMP END PARALLEL DO
 end if

!Treat the case of buffer bands
 if(band>max(1,nband-nbdbuf))then
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cwavef,ghc,gvnlc,npw1,nspinor)
   do ipw=1,npw1*nspinor
     cwavef(1:2,ipw)=zero
     ghc(1:2,ipw)   =zero
     gvnlc(1:2,ipw) =zero
   end do
!  $OMP END PARALLEL DO
   if (gen_eigenpb) then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(gsc,npw1,nspinor)
     do ipw=1,npw1*nspinor
       gsc(1:2,ipw)=zero
     end do
!    $OMP END PARALLEL DO
   end if
   if (usepaw==1) then
     call cprj_set_zero(cwaveprj)
   end if
!  A small negative residual will be associated with these
   resid=-0.1_dp
!  Number of one-way 3D ffts skipped
   nskip=nskip+nline

 else
!  If not a buffer band, perform the optimisation

   ABI_ALLOCATE(conjgr,(2,npw1*nspinor))
   ABI_ALLOCATE(direc,(2,npw1*nspinor))
   ABI_ALLOCATE(gresid,(2,npw1*nspinor))
   ABI_ALLOCATE(cwaveq,(2,npw1*nspinor))
   if (usepaw==1) then
     ABI_ALLOCATE(conjgrprj,(natom,nspinor))
     call cprj_alloc(conjgrprj,0,dimcprj)
   end if

   cwaveq(:,:)=cgq(:,1+npw1*nspinor*(band-1)+icgq:npw1*nspinor*band+icgq)
   dotgp=one

!  Here apply H(0) at k+q to input orthogonalized 1st-order wfs
   if (prtvol<0) call status(0,filstat,iexit,level,'call getghc(1)')
   sij_opt=0;if (gen_eigenpb) sij_opt=1
   cpopt=-1+usepaw
   call getghc(cpopt,cwavef,cwaveprj,dimffnl1,ffnl1,filstat,&
&   ghc,gsc,gs_hamkq,gvnlc,kg1_k,kinpw1,eshift,lmnmax,matblk,&
&   mgfft,mpi_enreg,mpsang,mpssoang,natom,1,npw1,nspinor,ntypat,nvloc,n4,n5,n6,&
&   paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
!  ghc also includes the eigenvalue shift
   if (gen_eigenpb) then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(eshift,ghc,gsc,npw1,nspinor)
     do ipw=1,npw1*nspinor
       ghc(1:2,ipw)=ghc(1:2,ipw)-eshift*gsc(1:2,ipw)
     end do
!    $OMP END PARALLEL DO
   else
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(cwavef,eshift,ghc,npw1,nspinor)
     do ipw=1,npw1*nspinor
       ghc(1:2,ipw)=ghc(1:2,ipw)-eshift*cwavef(1:2,ipw)
     end do
!    $OMP END PARALLEL DO
   end if
   if (prtvol<0) call status(0,filstat,iexit,level,'after getghc(1)')

!  Initialize resid, in case of nline==0
   resid=zero

!  ======================================================================
!  ====== BEGIN LOOP FOR A GIVEN BAND: MINIMIZATION ITERATIONS ==========
!  ======================================================================

   do iline=1,nline
     if (prtvol<0) call status(iline,filstat,iexit,level,'loop iline    ')

!    ======================================================================
!    ================= COMPUTE THE RESIDUAL ===============================
!    ======================================================================
!    Note that gresid (=steepest-descent vector, Eq.(26) of PRB 55, 10337 (1996))
!    is precomputed to garantee cancellation of errors
!    and allow residuals to reach values as small as 1.0d-24 or better.
     if (berryopt == 4) then
       if (ipert==natom+2) then
         gvnl1=zero
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(ghc,gh1c,npw1,nspinor,gresid,gberry,band)
         do ipw=1,npw1*nspinor
           gresid(1:2,ipw)=-ghc(1:2,ipw)-gh1c(1:2,ipw)
         end do
!        $OMP END PARALLEL DO
       else
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(ghc,gh1c,npw1,nspinor,gresid,gberry,band)
         do ipw=1,npw1*nspinor
           gresid(1,ipw)=-ghc(1,ipw)-gh1c(1,ipw)+gberry(2,ipw)
           gresid(2,ipw)=-ghc(2,ipw)-gh1c(2,ipw)-gberry(1,ipw)
         end do
!        $OMP END PARALLEL DO
       end if
     else
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(ghc,gh1c,npw1,nspinor,gresid)
       do ipw=1,npw1*nspinor
         gresid(1:2,ipw)=-ghc(1:2,ipw)-gh1c(1:2,ipw)
       end do
!      $OMP END PARALLEL DO
     end if

!    ======================================================================
!    =========== PROJECT THE STEEPEST DESCENT DIRECTION ===================
!    ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
!    ======================================================================
!    Project all bands from gresid into direc:
!    The following projection over the subspace orthogonal to occupied bands
!    is not optional in the RF case, unlike the GS case.
!    However, the order of operations could be changed, so that
!    as to make it only applied at the beginning, to H(1) psi(0),
!    so, THIS IS TO BE REEXAMINED
!    Note the subtlety:
!    -For the generalized eigenPb, S|cgq> is used in place of |cgq>,
!    in order to apply P_c+ projector (see PRB 73, 235101 (2006), Eq. (71), (72)
     if (prtvol<0) call status(iline,filstat,iexit,level,'projbd(1)     ')
     if (gen_eigenpb) then
       call projbd(gscq,gresid,-1,igscq,icgq,istwf_k,mgscq,mpi_enreg,mcgq,nband,npw1,nspinor,&
&       ortalg,iprint,cgq,scprod,0,tim_projbd,useoverlap)
     else
       call projbd(cgq,gresid,-1,icgq,0,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&       ortalg,iprint,dummy,scprod,0,tim_projbd,useoverlap)
     end if
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(direc,npw1,nspinor,gresid)
     do ipw=1,npw1*nspinor
       direc(1:2,ipw)=gresid(1:2,ipw)
     end do
!    $OMP END PARALLEL DO

!    ======================================================================
!    ============== CHECK FOR CONVERGENCE CRITERIA ========================
!    ======================================================================

!    Compute second-order derivative of the energy using a variational expression
     call dotprod_g(prod1,doti,istwf_k,mpi_enreg,npw1*nspinor,1,cwavef,gresid)
     call dotprod_g(prod2,doti,istwf_k,mpi_enreg,npw1*nspinor,1,cwavef,gh1c)
     d2te=two*(-prod1+prod2)
!    DEBUG
!    write(std_out,*)' cgwf3: prod1,prod2,d2te=',prod1,prod2,d2te
!    ENDDEBUG

!    Compute <u_m(1)|H(0)-e_m(0)|u_m(1)>
!    (<u_m(1)|H(0)-e_m(0).S|u_m(1)> if gen. eigenPb),
!    that should be positive,
!    except when the eigenvalue eig_mk(0) is higher than
!    the lowest non-treated eig_mk+q(0). For insulators, this
!    has no influence, but for metallic occupations,
!    the conjugate gradient algorithm breaks down. The solution adopted here
!    is very crude, and rely upon the fact that occupancies of such
!    levels should be smaller and smaller with increasing nband, so that
!    a convergence study will give the right result.
!    The same trick is also used later.
     u1h0me0u1=-prod1-prod2
!    Some tolerance is allowed, to account for very small numerical inaccuracies and cancellations.
     if(u1h0me0u1<-tol12)then
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(cwavef,ghc,gvnlc,npw1,nspinor)
       do ipw=1,npw1*nspinor
         cwavef(1:2,ipw)=zero
         ghc(1:2,ipw)   =zero
         gvnlc(1:2,ipw) =zero
       end do
!      $OMP END PARALLEL DO
       if (gen_eigenpb) then
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(gsc,npw1,nspinor)
         do ipw=1,npw1*nspinor
           gsc(1:2,ipw)=zero
         end do
!        $OMP END PARALLEL DO
       end if
       if (usepaw==1) then
         call cprj_set_zero(cwaveprj)
       end if
!      A negative residual will be the signal of this problem ...
       resid=-one
       write(msg, '(a)' )&
&       ' cgwf3: problem of minimisation (likely metallic), set resid to -1'
       call wrtout(std_out,msg,'PERS')
!      Number of one-way 3D ffts skipped
       nskip=nskip+(nline-iline+1)
!      Exit from the loop on iline
       exit
     end if

!    Compute residual (squared) norm
     call sqnorm_g(resid,istwf_k,mpi_enreg,npw1*nspinor,gresid)
     if (prtvol==-level)then
       write(msg,'(a,a,i3,f14.6,a,a,4es12.4)') ch10,&
&       ' cgwf3 : iline,eshift     =',iline,eshift,ch10,&
&       '         resid,prod1,prod2,d2te=',resid,prod1,prod2,d2te
       call wrtout(std_out,msg,'PERS')
     end if

!    If residual sufficiently small stop line minimizations
     if (resid<tolwfr) then
       if(prtvol>=10)then
         write(msg, '(a,i4,a,i2,a,es12.4)' ) &
&         ' cgwf3: band',band,' converged after ',iline,&
&         ' line minimizations : resid =',resid
         call wrtout(std_out,msg,'PERS')
       end if
!      Number of two-way 3D ffts skipped
       nskip=nskip+(nline-iline+1)
!      Exit from the loop on iline
       exit
     end if

!    If user require exiting the job, stop line minimisations
     if (quit==1) then
       write(msg, '(a,i4)' ) &
&       ' cgwf3: user require exiting => skip update of band ',band
       call wrtout(std_out,msg,'PERS')
!      Number of two-way 3D ffts skipped
       nskip=nskip+(nline-iline+1)
!      Exit from the loop on iline
       exit
     end if

!    Check that d2te is decreasing on succeeding lines:
     if (iline/=1) then
       if (d2te>d2teold+1.d-12) then
         write(msg, '(a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&         '  New trial energy at line',iline,' = ',d2te,ch10,&
&         '  is higher than former:',d2teold,ch10
         MSG_WARNING(msg)
       end if
     end if
     d2teold=d2te

!    DEBUG Keep this debugging feature !
!    call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,direc)
!    write(std_out,*)' cgwf3 : before precon, direc**2=',dotr
!    if (gen_eigenpb) then
!    call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1*nspinor,1,cwaveq,&
!    &                 gscq(:,1+npw1*nspinor*(band-1)+igscq:npw1*nspinor*band+igscq))
!    else
!    call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,cwaveq)
!    end if
!    write(std_out,*)' cgwf3 : before precon, cwaveq**2=',dotr
!    ENDDEBUG

!    ======================================================================
!    ======== PRECONDITION THE STEEPEST DESCENT DIRECTION =================
!    ======================================================================

!    If wfoptalg>=10, the precondition matrix is kept constant
!    during iteration ; otherwise it is recomputed
     if (wfoptalg<10.or.iline==1) then
       if (prtvol<0) call status(iline,filstat,iexit,level,'call precon   ')
       call precon(cwaveq,zero,istwf_k,kinpw1,mpi_enreg,npw1,nspinor,0,pcon,direc)
     else
       do ispinor=1,nspinor
         igs=(ispinor-1)*npw1
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(igs,npw,direc,pcon)
         do ipw=1+igs,npw1+igs
           direc(1:2,ipw)=direc(1:2,ipw)*pcon(ipw-igs)
         end do
!        $OMP END PARALLEL DO
       end do
     end if

!    DEBUG Keep this debugging feature !
!    call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,direc)
!    write(std_out,*)' cgwf3 : after precon, direc**2=',dotr
!    ENDDEBUG

!    ======================================================================
!    ======= PROJECT THE PRECOND. STEEPEST DESCENT DIRECTION ==============
!    ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
!    ======================================================================

!    Projecting again out all bands:
!    -For the simple eigenPb, gscq is used as dummy argument
     if (prtvol<0) call status(0,filstat,iexit,level,'prjbd(2)      ')
     call projbd(cgq,direc,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&     ortalg,iprint,gscq,scprod,0,tim_projbd,useoverlap)

!    DEBUG Keep this debugging feature !
!    call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,direc)
!    write(std_out,*)' cgwf3 : after projbd, direc**2=',dotr
!    ENDDEBUG

!    ======================================================================
!    ================= COMPUTE THE CONJUGATE-GRADIENT =====================
!    ======================================================================

!    get dot of direction vector with residual vector
     call dotprod_g(dotgg,doti,istwf_k,mpi_enreg,npw1*nspinor,1,direc,gresid)

!    At first iteration, gamma is set to zero
     if (iline==1) then
       gamma=zero
       dotgp=dotgg
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(conjgr,direc,npw1,nspinor)
       do ipw=1,npw1*nspinor
         conjgr(1:2,ipw)=direc(1:2,ipw)
       end do
!      $OMP END PARALLEL DO
     else
!      At next iterations, h = g + gamma * h
       gamma=dotgg/dotgp
       dotgp=dotgg
       if (prtvol==-level)then
         write(msg,'(a,2es16.6)') 'cgwf3: dotgg,gamma =',dotgg,gamma
         call wrtout(std_out,msg,'PERS')
       end if
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(conjgr,direc,gamma,npw1,nspinor)
       do ipw=1,npw1*nspinor
         conjgr(1:2,ipw)=direc(1:2,ipw)+gamma*conjgr(1:2,ipw)
       end do
!      $OMP END PARALLEL DO
       if (prtvol==-level)then
         write(msg,'(a)') &
&         'cgwf3: conjugate direction has been found'
         call wrtout(std_out,msg,'PERS')
       end if
     end if

!    ======================================================================
!    ===== COMPUTE CONTRIBUTIONS TO 1ST AND 2ND DERIVATIVES OF ENERGY =====
!    ======================================================================
!    ...along the search direction

!    Compute dedt, Eq.(29) of of PRB55, 10337 (1997),
!    with an additional factor of 2 for the difference
!    between E(2) and the 2DTE
     call dotprod_g(dedt,doti,istwf_k,mpi_enreg,npw1*nspinor,1,conjgr,gresid)
     dedt=-two*two*dedt

     ABI_ALLOCATE(gvnl_direc,(2,npw1*nspinor))
     ABI_ALLOCATE(gh_direc,(2,npw1*nspinor))
     if (gen_eigenpb)  then
       ABI_ALLOCATE(sconjgr,(2,npw1*nspinor))
     end if
     if (prtvol<0) call status(iline,filstat,iexit,level,'call getghc(2)')
     sij_opt=0;if (gen_eigenpb) sij_opt=1
     cpopt=-1+usepaw
     call getghc(cpopt,conjgr,conjgrprj,dimffnl1,ffnl1,filstat,&
&     gh_direc,sconjgr,gs_hamkq,gvnl_direc,kg1_k,kinpw1,eshift,lmnmax,matblk,&
&     mgfft,mpi_enreg,mpsang,mpssoang,natom,1,npw1,nspinor,ntypat,nvloc,n4,n5,n6,&
&     paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
!    ghc also includes the eigenvalue shift
     if (gen_eigenpb) then
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(sconjgr,eshift,gh_direc,npw1,nspinor)
       do ipw=1,npw1*nspinor
         gh_direc(1:2,ipw)=gh_direc(1:2,ipw)-eshift*sconjgr(1:2,ipw)
       end do
!      $OMP END PARALLEL DO
     else
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(conjgr,eshift,gh_direc,npw1,nspinor)
       do ipw=1,npw1*nspinor
         gh_direc(1:2,ipw)=gh_direc(1:2,ipw)-eshift*conjgr(1:2,ipw)
       end do
!      $OMP END PARALLEL DO
     end if
     if (prtvol<0) call status(iline,filstat,iexit,level,'after getghc(2)')

!    compute d2edt2, Eq.(30) of of PRB55, 10337 (1997),
!    with an additional factor of 2 for the difference
!    between E(2) and the 2DTE, and neglect of local fields (SC terms)
     call dotprod_g(d2edt2,doti,istwf_k,mpi_enreg,npw1*nspinor,1,conjgr,gh_direc)
     d2edt2=two*two*d2edt2
     if(prtvol==-level)then
       write(msg,'(a,2es14.6)') 'cgwf3: dedt,d2edt2=',dedt,d2edt2
       call wrtout(std_out,msg,'PERS')
     end if

!    ======================================================================
!    ======= COMPUTE MIXING FACTOR - CHECK FOR CONVERGENCE ===============
!    ======================================================================

!    see Eq.(31) of PRB55, 10337 (1997)
!    
!    if(d2edt2<=zero)then
     if(d2edt2<-tol12)then
!      This may happen when the eigenvalue eig_mk(0) is higher than
!      the lowest non-treated eig_mk+q(0). The solution adopted here
!      is very crude, and rely upon the fact that occupancies of such
!      levels should be smaller and smaller with increasing nband, so that
!      a convergence study will give the right result.
       theta=zero
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(cwavef,ghc,gvnlc,npw1,nspinor)
       do ipw=1,npw1*nspinor
         cwavef(1:2,ipw)=zero
         ghc(1:2,ipw)   =zero
         gvnlc(1:2,ipw) =zero
       end do
!      $OMP END PARALLEL DO
       if (gen_eigenpb) then
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(gsc,npw1,nspinor)
         do ipw=1,npw1*nspinor
           gsc(1:2,ipw)=zero
         end do
!        $OMP END PARALLEL DO
       end if
       if (usepaw==1) then
         call cprj_set_zero(cwaveprj)
       end if
!      A negative residual will be the signal of this problem ...
       resid=-two
       write(msg, '(a)' )&
&       ' cgwf3: problem of minimisation (likely metallic), set resid to -2'
       call wrtout(std_out,msg,'PERS')
     else
!      Here, the value of theta that gives the minimum
       theta=-dedt/d2edt2
!      DEBUG
!      write(std_out,*)' cgwf3: dedt,d2edt2=',dedt,d2edt2
!      ENDDEBUG
     end if

!    Check that result is above machine precision
     if (one+theta==one) then
       write(msg, '(a,es16.4)' ) ' cgwf3: converged with theta=',theta
       call wrtout(std_out,msg,'PERS')
!      Number of one-way 3D ffts skipped
       nskip=nskip+2*(nline-iline)
!      Exit from the loop on iline
       exit
     end if

!    ======================================================================
!    ================ GENERATE NEW |wf>, H|wf>, Vnl|Wf ... ================
!    ======================================================================

!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(conjgr,cwavef,ghc,gvnlc,gh_direc,npw1,nspinor,theta)
     do ipw=1,npw1*nspinor
       cwavef(1:2,ipw)=cwavef(1:2,ipw)+conjgr(1:2,ipw)    *theta
       ghc(1:2,ipw)   =ghc(1:2,ipw)   +gh_direc(1:2,ipw)  *theta
       gvnlc(1:2,ipw) =gvnlc(1:2,ipw) +gvnl_direc(1:2,ipw)*theta
     end do
!    $OMP END PARALLEL DO
     if (gen_eigenpb) then
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(conjgr,sconjgr,gsc,npw1,nspinor,theta)
       do ipw=1,npw1*nspinor
         gsc(1:2,ipw)=gsc(1:2,ipw)+sconjgr(1:2,ipw)*theta
       end do
!      $OMP END PARALLEL DO
     end if
     if (usepaw==1) then
       call cprj_axpby(theta,one,conjgrprj,cwaveprj)
     end if
     ABI_DEALLOCATE(gh_direc)
     ABI_DEALLOCATE(gvnl_direc)
     if (gen_eigenpb)  then
       ABI_DEALLOCATE(sconjgr)
     end if

!    ======================================================================
!    =========== CHECK CONVERGENCE AGAINST TRIAL ENERGY ===================
!    ======================================================================

!    Check reduction in trial energy deltae, Eq.(28) of PRB55, 10337 (1997)
     deltae=half*d2edt2*theta**2+theta*dedt

     if (iline==1) then
       deold=deltae
!      Use a lower value for comparison with RESPFN
!      cgwftol=tol14
       cgwftol=0.01_dp   ! Value used until v3.3
       if (prtvol==-level) cgwftol=tol14
     else if (abs(deltae)<cgwftol*abs(deold) .and. iline/=nline ) then
!      else if (abs(deltae)<0.005_dp*abs(deold) .and. iline/=nline ) then
       if(prtvol>=10.or.prtvol==-level)then
         write(msg, '(a,i4,1x,a,1p,e12.4,a,e12.4,a)' ) &
&         ' cgwf3: line',iline,&
&         ' deltae=',deltae,' < cgwftol*',deold,' =>skip lines'
         call wrtout(std_out,msg,'PERS')
       end if
!      Number of one-way 3D ffts skipped
       nskip=nskip+2*(nline-iline)
!      Exit from the loop on iline
       exit
     end if

!    ======================================================================
!    ================== END LOOP FOR GIVEN BAND ===========================
!    ======================================================================

!    Note that there are three "exit" instruction inside the loop.
   end do ! iline

!  Check that final cwavef (Psi^(1)) satisfies the orthogonality condition
   if (prtvol==-level) then
     sij_opt=0 ; usevnl=0 ; optlocal=1 ; optnl=2 ; if (gen_eigenpb)  sij_opt=1
     ABI_ALLOCATE(work,(2,npw1*nspinor))
     ABI_ALLOCATE(work1,(2,npw1*nspinor))
     ABI_ALLOCATE(work2,(2,npw1*nspinor*sij_opt))
     do iband=1,nband
       if (gen_eigenpb) then
         work(:,:)=gscq(:,1+npw1*nspinor*(iband-1)+igscq:npw1*nspinor*iband+igscq)
       else
         work(:,:)=cgq(:,1+npw1*nspinor*(iband-1)+icgq:npw1*nspinor*iband+icgq)
       end if
       call dotprod_g(dotgp,dotgg,istwf_k,mpi_enreg,npw1*nspinor,2,cwavef,work)
       if (gen_eigenpb) then
         call getgh1c(berryopt,cplex,cwave0,cwaveprj0,dimcprj,dimekb,dime1kb,dimffnlk,dimffnl1,dimffnl1,dimphkxred,dkinpw,&
&         ekb_typ,e1kbfr,e1kbsc,ffnlk,ffnlkq,ffnl1,filstat,gbound,work1,gberry,work2,gs_hamkq,dummy,idir,indlmn_typ,&
&         ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,eshift,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&         nkpg,nkpg1,npw,npw1,nspinor,ntypat,n4,n5,n6,optlocal,optnl,opt_gvnl1,paral_kgb,ph3d,phkxred,prtvol,&
&         sij_opt,sij_typ,tim_getgh1c,usecprj,usee1kb,usevnl,vlocal1,wfraug)
         work(:,:)=cgq(:,1+npw1*nspinor*(iband-1)+icgq:npw1*nspinor*iband+icgq)
         call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1*nspinor,2,work2,work)
       else
         dotr=zero;doti=zero
       end if
       dotr=dotgp+half*dotr
       doti=dotgg+half*doti
       if (gen_eigenpb) then
         write(msg,'(2a,i3,a,2es22.15)') '<Psi^(1)_i,k,q|Psi^(0)_j,k+q>',&
&         '+ 1/2<Psi^(0)_i,k|S^(1)|Psi^(0)_j,k+q>, for j= ',iband,' is ',dotr,doti
       else
         write(msg,'(a,i3,a,2es22.15)') '<Psi^(1)_i,k,q|Psi^(0)_j,k+q>, for j= ',iband,' is ',dotr,doti
       end if
       call wrtout(std_out,msg,'PERS')
     end do
     ABI_DEALLOCATE(work)
     ABI_DEALLOCATE(work1)
     ABI_DEALLOCATE(work2)
   end if

!  Check that final cwavef Psi^(1) is Pc.Psi^(1)+delta_Psi^(1)
   if (prtvol==-level)then
     ABI_ALLOCATE(cwwork,(2,npw1*nspinor))
     cwwork=cwavef
!    -Apply Pc to Psi^(1)
     if (gen_eigenpb) then
       call projbd(cgq,cwwork,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&       ortalg,iprint,gscq,scprod,0,tim_projbd,useoverlap)
     else
       call projbd(cgq,cwwork,-1,icgq,0,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&       ortalg,iprint,dummy,scprod,0,tim_projbd,useoverlap)
     end if
!    -Add delta_Psi^(1)
     if (usedcwavef>0) cwwork=cwwork+dcwavef
!    -Compare to Psi^(1)
     cwwork=cwwork-cwavef
     call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,cwwork)
     ABI_DEALLOCATE(cwwork)
     if (gen_eigenpb) then
       write(msg,'(a,i3,a,es22.15)') &
&       '|(Pc.Psi^(1)_i,k,q + delta_Psi^(1)_i,k) - Psi^(1)_i,k,q|^2 (band ',band,')=',dotr
     else
       write(msg,'(a,i3,a,es22.15)') '|Pc.Psi^(1)_i,k,q - Psi^(1)_i,k,q|^2 = ',dotr
     end if
     call wrtout(std_out,msg,'PERS')
   end if

!  Check that final cwavef (Psi^(1)) solves the Sternheimer equation
   if(prtvol==-level)then
     ABI_ALLOCATE(cwwork,(2,npw1*nspinor))
     cwwork=cwavef
!    -Apply Pc to Psi^(1)
     if (gen_eigenpb) then
       call projbd(cgq,cwwork,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&       ortalg,iprint,gscq,scprod,0,tim_projbd,useoverlap)
     else
       call projbd(cgq,cwwork,-1,icgq,0,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&       ortalg,iprint,dummy,scprod,0,tim_projbd,useoverlap)
     end if
!    - Apply H^(0)-E.S^(0)
     sij_opt=0;if (gen_eigenpb) sij_opt=1
     cpopt=-1
     if (usepaw==0)  then
       ABI_ALLOCATE(conjgrprj,(natom,0))
     end if
     ABI_ALLOCATE(work,(2,npw1*nspinor))
     ABI_ALLOCATE(work1,(2,npw1*nspinor*((sij_opt+1)/2)))
     ABI_ALLOCATE(work2,(2,npw1*nspinor))
     call getghc(cpopt,cwwork,conjgrprj,dimffnl1,ffnl1,filstat,&
&     work,work1,gs_hamkq,work2,kg1_k,kinpw1,eshift,lmnmax,matblk,&
&     mgfft,mpi_enreg,mpsang,mpssoang,natom,1,npw1,nspinor,ntypat,nvloc,n4,n5,n6,&
&     paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
     if (gen_eigenpb) then
       cwwork=work-eshift*work1
     else
       cwwork=work-eshift*cwwork
     end if
     ABI_DEALLOCATE(work)
     ABI_DEALLOCATE(work1)
     ABI_DEALLOCATE(work2)
     if (usepaw==0)  then
       ABI_DEALLOCATE(conjgrprj)
     end if
!    -Apply Pc^*
     if (gen_eigenpb) then
       call projbd(gscq,cwwork,-1,igscq,icgq,istwf_k,mgscq,mpi_enreg,mcgq,nband,npw1,nspinor,&
&       ortalg,iprint,cgq,scprod,0,tim_projbd,useoverlap)
     else
       call projbd(cgq,cwwork,-1,icgq,0,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&       ortalg,iprint,dummy,scprod,0,tim_projbd,useoverlap)
     end if
!    - Add Pc^*(H^(1)-E.S^(1)).Psi^(0)
     cwwork=cwwork+gh1c
     call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,cwwork)
     ABI_DEALLOCATE(cwwork)
     write(msg,'(a,i3,a,es22.15,2a)') &
&     'Sternheimer equation test for band ',band,'=',dotr,ch10,&
&     '  (should be zero for large nline)'
     call wrtout(std_out,msg,'PERS')
   end if

   if (allocated(gh_direc))  then
     ABI_DEALLOCATE(gh_direc)
   end if
   if (allocated(gvnl_direc))  then
     ABI_DEALLOCATE(gvnl_direc)
   end if
   ABI_DEALLOCATE(conjgr)
   ABI_DEALLOCATE(cwaveq)
   ABI_DEALLOCATE(direc)
   ABI_DEALLOCATE(gresid)
   if (usepaw==1) then
     call cprj_free(conjgrprj)
     ABI_DEALLOCATE(conjgrprj)
   end if

!  End condition of not being a buffer band
 end if

 if (prtvol<0) call status(0,filstat,iexit,level,'after iline   ')

!At the end of the treatment of a set of bands, write the number
!of one-way 3D ffts skipped
 if(mpi_enreg%paral_compil_kpt==1)then
   if(band==nband .and. prtvol>=10)then
     write(msg, '(a,i8)' )&
&     ' cgwf3: number of one-way 3D ffts skipped in cgwf3 until now =',nskip
     call wrtout(std_out,msg,'PERS')
   end if
 end if

 ABI_DEALLOCATE(gh1c)
 ABI_DEALLOCATE(pcon)
 ABI_DEALLOCATE(scprod)
 ABI_DEALLOCATE(wfraug)
 if (berryopt==4)  then
   ABI_DEALLOCATE(gberry)
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
   write(msg,'(a1,a,a1,a,i2,a)') ch10,&
&   ' cgwf3 : exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(std_out,msg,'PERS')
   call leave_new('PERS')
 end if

 call timab(122,2,tsec)
 if(prtvol<0) call status(0,filstat,iexit,level,'exit          ')

 DBG_EXIT("COLL")

end subroutine cgwf3
!!***
