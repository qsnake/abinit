!{\src2tex{textfont=tt}}
!!****f* ABINIT/getgh1c
!!
!! NAME
!! getgh1c
!!
!! FUNCTION
!! Compute <G|H^(1)|C> (or <G|H^(1)-Eps.S^(1)|C>) for input vector |C> expressed in reciprocal space.
!! (H^(1) is the 1st-order pertubed Hamiltonian, S^(1) is the 1st-order perturbed overlap operator).
!! Result is put in array gh1c.
!! If required, part of <G|K(1)+Vnonlocal^(1)|C> not depending on VHxc^(1) is also returned in gvnl1c.
!! If required, <G|S^(1)|C> is returned in gs1c (S=overlap - PAW only)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  berryopt=option for Berry phase
!!  cplex=1 if vlocal1 is real, 2 if vlocal1 is complex
!!  cwave(2,npw*nspinor)=input wavefunction, in reciprocal space
!!  cwaveprj(natom,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C> (and 1st derivatives)
!!  dimcprj(natom*usepaw)=array of dimensions of array cwaveprj (ordered by atom-type)
!!  dimekb=first dimension of ekb (see ekb_typ below)
!!  dime1kb=first dimension of e1kbsc and e1kbfr (see e1kbsc,e1kbfr below)
!!  dimffnlk =second dimension of ffnl   (1+number of derivatives)
!!  dimffnlkq=second dimension of ffnlkq (1+number of derivatives)
!!  dimffnl1 =second dimension of ffnl1  (1+number of derivatives)
!!  dimphkxred=second dimension of phkxred
!!  dkinpw(npw)=derivative of the (modified) kinetic energy for each plane wave at k (Hartree)
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
!!  ffnlk(npw,dimffnlk,lmnmax,1+usepaw(1-ntypat))=nonloc form factors at k, for the displaced atom.
!!  ffnlkq(npw1,dimffnl1,lmnmax,1)=nonloc form fact at k+q for the displaced atom
!!  ffnl1(npw1,dimffnl1,lmnmax,ntypat)=nonloc form factors at k+q
!!  filstat=name of the status file
!!  gbound(2*mgfft+8,2)=G sphere boundary at k
!!  grad_berry(2,npw1*nspinor)= the gradient of the Berry phase term
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  idir=direction of the perturbation
!!  indlmn_typ(6,lmnmax,1)=indlmn info for the displaced atom
!!  ipert=type of the perturbation
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere at k.
!!  kg1_k(3,npw1)=coordinates of planewaves in basis sphere at k+q.
!!  kinpw1(npw1)=(modified) kinetic energy for each plane wave at k+q (Hartree)
!!  kpg_k(npw,nkpg)= (k+G) components at k (only if useylm=1)
!!  kpg1_k(npw1,nkpg1)= (k+G) components at k+q (only if useylm=1)
!!  kpt(3)=coordinates of k point.
!!  lmnmax= max number of (l,n) comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in unit cell.
!!  nkpg,nkpg1=second dimensions of kpg_k and kpg1_k (0 if useylm=0)
!!  npw=number of planewaves in basis sphere at given k.
!!  npw1=number of planewaves in basis sphere at k+q
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell.
!!  n4,n5,n6 dimensions of vlocal and vlocal1
!!  optlocal=0: local part of H^(1) is not computed in gh1c=<G|H^(1)|C>
!!           1: local part of H^(1) is computed in gh1c=<G|H^(1)|C>
!!  optnl=0: non-local part of H^(1) is not computed in gh1c=<G|H^(1)|C>
!!        1: non-local part of H^(1) depending on VHxc^(1) is not computed in gh1c=<G|H^(1)|C>
!!        2: non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
!!  opt_gvnl1=option controlling the use of gvnl1 array:
!!            0: used as an output
!!            1: used as an input:    - used only for ipert=natom+2
!!                 NCPP: contains the ddk 1-st order WF
!!                 PAW: contains frozen part of 1st-order hamiltonian
!!            2: used as input/ouput:    - used only for PAW and ipert=natom+2
!!                 At input: contains the ddk 1-st order WF (times i)
!!                 At output: contains frozen part of 1st-order hamiltonian
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  ph3d(2,npw1,matblk)=3-dim structure factors, for each atom and plane wave.
!!  phkxred(2,dimphkxred)=phase factors exp(2 pi kpoint.xred) at k
!!  prtvol=control print volume and debugging output
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H^(1)|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S^(1)|C> have to be computed in gs1c in addition to gh1c
!!                       if -1, matrix elements <G|H^(1)-lambda.S^(1)|C> have to be computed in gh1c (gs1c not used)
!!  sij_typ(dimekb,1)=overlap matrix components (only if sij_opt/=0)
!!  tim_getgh1c=timing code of the calling subroutine (can be set to 0 if not attributed)
!!  usecprj=1 if cwaveprj coefficients (and 1st derivatives) are already in memory (PAW only)
!!  usee1kb=1 if ekb derivatives (e1kbsc, e1kbfr) exist
!!  usevnl=1 if gvnl1=(part of <G|K^(1)+Vnl^(1)-lambda.S^(1)|C> not depending on VHxc^(1)) has to be input/output
!!  vlocal1(cplex*n4,n5,n6*optlocal)= 1st-order local pot in real space, on the augmented fft grid
!!
!! OUTPUT
!! gh1c(2,npw*nspinor)= <G|H^(1)|C> or  <G|H^(1)-lambda.S^(1)|C>
!!                     (only kinetic+non-local parts if optlocal=0)
!! if (usevnl==1)
!!  gvnl1(2,npw1*nspinor)=  part of <G|K^(1)+Vnl^(1)|C> not depending on VHxc^(1)              (sij_opt/=-1)
!!                       or part of <G|K^(1)+Vnl^(1)-lambda.S^(1)|C> not depending on VHxc^(1) (sij_opt==-1)
!! if (sij_opt=1)
!!  gs1c(2,npw*nspinor)=<G|S^(1)|C> (S=overlap).
!!
!! SIDE EFFECTS
!! wfraug(2,n4,n5,n6)=is a dummy array (used only if optlocal=1)
!!
!! PARENTS
!!      cgwf3,nstpaw3,nstwf3
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,fourwf,leave_new,nonlop,status,timab
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine getgh1c(berryopt,cplex,cwave,cwaveprj,dimcprj,dimekb,dime1kb,dimffnlk,dimffnlkq,dimffnl1,dimphkxred,dkinpw,&
&          ekb_typ,e1kbfr,e1kbsc,ffnlk,ffnlkq,ffnl1,filstat,gbound,gh1c,grad_berry,gs1c,gs_hamkq,gvnl1,idir,&
&          indlmn_typ,ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&
&          mpssoang,natom,nkpg,nkpg1,npw,npw1,nspinor,ntypat,n4,n5,n6,optlocal,optnl,opt_gvnl1,paral_kgb,ph3d,&
&          phkxred,prtvol,sij_opt,sij_typ,tim_getgh1c,usecprj,usee1kb,usevnl,vlocal1,wfraug)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getgh1c'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_53_ffts
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,cplex,dime1kb,dimekb,dimffnl1,dimffnlk,dimffnlkq
 integer,intent(in) :: dimphkxred,idir,ipert,lmnmax,matblk,mgfft,mpsang
 integer,intent(in) :: mpssoang,n4,n5,n6,natom,nkpg,nkpg1,npw,npw1,nspinor
 integer,intent(in) :: ntypat,optlocal,optnl,opt_gvnl1,paral_kgb,prtvol,sij_opt
 integer,intent(in) :: tim_getgh1c,usecprj,usee1kb,usevnl
 real(dp),intent(in) :: lambda
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
!arrays
 integer,intent(in) :: dimcprj(natom*gs_hamkq%usepaw),gbound(2*mgfft+8,2)
 integer,intent(in) :: indlmn_typ(6,lmnmax,1),kg1_k(3,npw1),kg_k(3,npw)
 real(dp),intent(in) :: dkinpw(npw)
 real(dp),intent(in) :: e1kbfr(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
 real(dp),intent(in) :: e1kbsc(dime1kb,gs_hamkq%dimekb2,usee1kb*nspinor**2)
 real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
 real(dp),intent(in) :: ffnl1(npw1,dimffnl1,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlk(npw,dimffnlk,lmnmax,1+gs_hamkq%usepaw*(ntypat-1))
 real(dp),intent(in) :: ffnlkq(npw1,dimffnlkq,lmnmax,1)
 real(dp),intent(in) :: grad_berry(2,npw1*nspinor*(berryopt/4)),kinpw1(npw1)
 real(dp),intent(in) :: kpg1_k(npw1,nkpg1),kpg_k(npw,nkpg),kpt(3)
 real(dp),intent(in) :: phkxred(2,dimphkxred),sij_typ(dimekb,gs_hamkq%usepaw)
 real(dp),intent(inout) :: cwave(2,npw*nspinor),ph3d(2,npw1,matblk)
 real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6*optlocal)
 real(dp),intent(inout) :: wfraug(2,n4,n5,n6*optlocal)
 real(dp),intent(inout),target :: gvnl1(2,npw1*nspinor*usevnl)
 real(dp),intent(out) :: gh1c(2,npw1*nspinor)
 real(dp),intent(out) :: gs1c(2,npw1*nspinor*((sij_opt+1)/2))
 type(cprj_type),intent(inout),target :: cwaveprj(natom,nspinor*usecprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=16
 integer :: choice,cpopt,dimekb2_der,iatm,iexit,ipw,ipws,ispinor,istr
 integer :: matblk_der,n1,n2,n3,natom_der,nnlout,ntypat_der,paw_opt,shift1
 integer :: shift2,shift3,signs,tim_fourwf,tim_nonlop
 logical :: has_kin,usevnl2
 real(dp) :: weight
 character(len=500) :: message
!arrays
 integer :: atindx1_der(1),atindx_der(1),dimcprj_der(1),nattyp_der(1)
 integer :: nloalg_der(5)
 real(dp) :: phkxredin(2,1),phkxredout(2,1),tsec(2),xred_der(3)
 real(dp) :: sij_dum(1,1),svectout_dum(1,1),vectout_dum(1,1)
 real(dp),allocatable :: cwave_sp(:,:),enlout(:),ffnlk_der(:,:,:,:),gh1c_sp(:,:)
 real(dp),allocatable :: gvnl2(:,:),ph1d_der(:,:),ph3din(:,:,:),ph3dout(:,:,:),work(:,:)
 real(dp),pointer :: gvnl1_(:,:)
 type(cprj_type),allocatable :: cprj_der(:,:)
 type(cprj_type),allocatable,target :: cwaveprj_tmp(:,:)
 type(cprj_type),pointer :: cwaveprj_ptr(:,:)

! *********************************************************************

!Keep track of total time spent in getgh1c
 call timab(196+tim_getgh1c,1,tsec)
 if(prtvol<0)then
   call status(0,filstat,iexit,level,'enter         ')
 end if

!Compatibility tests
 if(gs_hamkq%usepaw==0.and.usecprj==1)then
   write(message,'(4a)') ch10,&
&   ' getgh1c : BUG -',ch10,&
&   '   usecprj==1 not allowed for NC psps !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(gs_hamkq%usepaw==1.and.(ipert>=0.and.ipert<=natom)) then
   if (usee1kb==0)then
     write(message,'(4a)') ch10,&
&     ' getgh1c : BUG -',ch10,&
&     '   ekb derivatives must be allocated for ipert<=natom !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if (usecprj==1) then
     if (cwaveprj(1,1)%ncpgr/=1)then
       write(message,'(4a)') ch10,&
&       ' getgh1c : BUG -',ch10,&
&       '   Projected WFs (cprj) derivatives are not correctly stored !'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end if
 end if
 if(gs_hamkq%usepaw==1.and.(ipert==natom+2)) then
   if (usee1kb==0)then
     write(message,'(4a)') ch10,&
&     ' getgh1c : BUG -',ch10,&
&     '   ekb derivatives must be allocated for ipert=natom+2 !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if (usevnl==0) then
     write(message,'(4a)') ch10,&
&     ' getgh1c : BUG -',ch10,&
&     '   gvnl1 must be allocated for ipert=natom+2 !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
 if(ipert==natom+2.and.opt_gvnl1==0) then
   write(message,'(4a)') ch10,&
&   ' getgh1c : BUG -',ch10,&
&   '   opt_gvnl1=0 not compatible with ipert=natom+2 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (mpi_enreg%paral_spin==1) then
   write(message,'(4a)') ch10,&
&   ' getgh1c : BUG -',ch10,&
&   '   Not compatible with parallelization over spinorial components !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 tim_nonlop=8
 if (tim_getgh1c==1.and.ipert<=natom) tim_nonlop=7
 if (tim_getgh1c==2.and.ipert<=natom) tim_nonlop=5
 if (tim_getgh1c==1.and.ipert>natom)  tim_nonlop=8

!======================================================================
!== Apply the 1st-order local potential to the wavefunction
!======================================================================

!Phonon perturbation
!or Electric field perturbation
!or Strain perturbation
!-------------------------------------------
 if(ipert<=natom+4.and.ipert/=natom+1.and.optlocal>0) then

   weight=one ; tim_fourwf=4
   call fourwf(cplex,vlocal1,cwave,gh1c,wfraug,gbound,gs_hamkq%gbound,&
&   gs_hamkq%istwf_k,kg_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&   npw,npw1,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&   use_gpu_cuda=gs_hamkq%use_gpu_cuda)
   if(nspinor==2)then
     ABI_ALLOCATE(cwave_sp,(2,npw))
     ABI_ALLOCATE(gh1c_sp,(2,npw1))
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(cwave,cwave_sp,npw)
     do ipw=1,npw
       cwave_sp(1,ipw)=cwave(1,ipw+npw)
       cwave_sp(2,ipw)=cwave(2,ipw+npw)
     end do
!    $OMP END PARALLEL DO
     call fourwf(cplex,vlocal1,cwave_sp,gh1c_sp,wfraug,gbound,gs_hamkq%gbound,&
&     gs_hamkq%istwf_k,kg_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     npw,npw1,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&     use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(gh1c,gh1c_sp,npw1)
     do ipw=1,npw1
       gh1c(1,ipw+npw1)=gh1c_sp(1,ipw)
       gh1c(2,ipw+npw1)=gh1c_sp(2,ipw)
     end do
!    $OMP END PARALLEL DO
     ABI_DEALLOCATE(cwave_sp)
     ABI_DEALLOCATE(gh1c_sp)
   end if

!  k-point perturbation or rfmgd (or no local part, i.e. optlocal=0)
!  -------------------------------------------
 else if(ipert==natom+1.or.ipert==natom+5.or.optlocal==0)then

!  In the case of ddk operator, no local contribution
!  (also because no self-consistency)
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(gh1c,npw1,nspinor)
   do ipw=1,npw1*nspinor
     gh1c(:,ipw)=zero
   end do
!  $OMP END PARALLEL DO

 end if

!======================================================================
!== Apply the 1st-order non-local potential to the wavefunction
!======================================================================

!Use of gvnl1 depends on usevnl
 if (usevnl==1) then
   gvnl1_ => gvnl1
 else
   ABI_ALLOCATE(gvnl1_,(2,npw1*nspinor))
 end if
 usevnl2=(gs_hamkq%usepaw==1.and.optnl>=2.and.&
& ((ipert>0.and.ipert<=natom).or.(ipert==natom+2)))

!Phonon perturbation
!-------------------------------------------
 if(ipert<=natom.and.(optnl>0.or.sij_opt/=0)) then

   signs=2 ; nnlout=1 ; natom_der=1 ; nattyp_der(1)=1 ; ntypat_der=1
   dimekb2_der=1 
   ABI_ALLOCATE(enlout,(nnlout))
   matblk_der=1
   xred_der(:)=gs_hamkq%xred(:,ipert)
   atindx_der(1)=1 ; atindx1_der(1)=1
   n1=gs_hamkq%ngfft(1) ; n2=gs_hamkq%ngfft(2) ; n3=gs_hamkq%ngfft(3)
   iatm=gs_hamkq%atindx(ipert)

!  Store at the right place the 1d phases
   ABI_ALLOCATE(ph1d_der,(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
   shift1=(iatm-1)*(2*n1+1)
   ph1d_der(:,1:2*n1+1)=gs_hamkq%ph1d(:,1+shift1:2*n1+1+shift1)
   shift2=(iatm-1)*(2*n2+1)+natom*(2*n1+1)
   ph1d_der(:,1+2*n1+1:2*n2+1+2*n1+1)=gs_hamkq%ph1d(:,1+shift2:2*n2+1+shift2)
   shift3=(iatm-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
   ph1d_der(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=&
&   gs_hamkq%ph1d(:,1+shift3:2*n3+1+shift3)

!  Will compute the 3D phase factors inside nonlop
   ABI_ALLOCATE(ph3din,(2,npw,1))
   ABI_ALLOCATE(ph3dout,(2,npw1,1))
   nloalg_der(:)=gs_hamkq%nloalg(:)
   nloalg_der(1)=-abs(gs_hamkq%nloalg(1))
   nloalg_der(4)=1

!  Retrieve here phkxred for kpt and kpq
   if (gs_hamkq%usepaw==0) then
     phkxredin(:,1)=phkxred(:,1)
   else
     phkxredin(:,1)=phkxred(:,iatm)
   end if
   phkxredout(:,1)=gs_hamkq%phkxred(:,iatm)

!  PAW: retrieve ffnlk and cwaveprj for the displaced atom
   if (gs_hamkq%usepaw==1) then
     ABI_ALLOCATE(ffnlk_der,(npw,dimffnlk,lmnmax,1))
     ffnlk_der(:,:,:,1)=ffnlk(:,:,:,gs_hamkq%typat(ipert))
     if (usecprj==1) then
       dimcprj_der(1)=dimcprj(iatm)
       ABI_ALLOCATE(cprj_der,(1,nspinor))
       call cprj_alloc(cprj_der,1,dimcprj_der)
       do ispinor=1,nspinor
         call cprj_copy(cwaveprj(iatm:iatm,ispinor:ispinor),cprj_der(1:1,ispinor:ispinor))
       end do
     end if
     if (usecprj==1) then
       cwaveprj_ptr => cwaveprj
     else
       ABI_ALLOCATE(cwaveprj_tmp,(natom,nspinor))
       call cprj_alloc(cwaveprj_tmp,1,dimcprj)
       cwaveprj_ptr => cwaveprj_tmp
     end if
   end if

!  Application of 1st-order nl operator

!  PAW:
   if (gs_hamkq%usepaw==1) then

!    1- Compute derivatives due to projectors |p_i>^(1)
!    Only displaced atom contributes
     cpopt=-1+5*usecprj ; choice=2
     paw_opt=1;if (sij_opt/=0) paw_opt=sij_opt+3
     call nonlop(atindx1_der,choice,cpopt,cprj_der,dimekb,dimekb2_der,dimffnlk,dimffnlkq,&
&     ekb_typ,enlout,ffnlk_der,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&     indlmn_typ,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&     lambda,lmnmax,matblk_der,mgfft,mpi_enreg,mpsang,mpssoang,natom_der,nattyp_der,&
&     gs_hamkq%ngfft,nkpg,nkpg1,nloalg_der,nnlout,npw,npw1,nspinor,nspinor,ntypat_der,&
&     0,paw_opt,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,signs,&
&     sij_typ,gs1c,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1_,&
&     use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!    2- Compute derivatives due to frozen part of D_ij^(1) (independent of VHxc^(1))
!    All atoms contribute (unfortunately !)
     if (optnl>=1) then
       ABI_ALLOCATE(gvnl2,(2,npw1*nspinor))
       cpopt=1+3*usecprj ; choice=1 ; paw_opt=1
       call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj_ptr,dime1kb,gs_hamkq%dimekb2,dimffnlk,dimffnl1,&
&       e1kbfr,enlout,ffnlk,ffnl1,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&       gs_hamkq%indlmn,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&       lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,&
&       gs_hamkq%ngfft,nkpg,nkpg1,gs_hamkq%nloalg,nnlout,npw,npw1,nspinor,nspinor,ntypat,&
&       0,paw_opt,phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,signs,&
&       sij_dum,svectout_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl2,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(gvnl1_,gvnl2,npw1)
       do ipw=1,npw1*nspinor
         gvnl1_(:,ipw)=gvnl1_(:,ipw)+gvnl2(:,ipw)
       end do
!      $OMP END PARALLEL DO
     end if
!    3- Compute derivatives due to part of D_ij^(1) depending on VHxc^(1)
!    All atoms contribute (unfortunately !)
     if (optnl>=2) then
       cpopt=4 ; choice=1 ; paw_opt=1
       call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj_ptr,dime1kb,gs_hamkq%dimekb2,dimffnlk,dimffnl1,&
&       e1kbsc,enlout,ffnlk,ffnl1,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&       gs_hamkq%indlmn,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&       lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,&
&       gs_hamkq%ngfft,nkpg,nkpg1,gs_hamkq%nloalg,nnlout,npw,npw1,nspinor,nspinor,ntypat,&
&       0,paw_opt,phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,signs,&
&       sij_dum,svectout_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl2,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     else if (optnl==1) then
       ABI_DEALLOCATE(gvnl2)
     end if

!    Norm-conserving psps:
   else
!    Compute only derivatives due to projectors |p_i>^(1)
     cpopt=-1 ; choice=2 ; paw_opt=0
     call nonlop(atindx1_der,choice,cpopt,cwaveprj,dimekb,dimekb2_der,dimffnlk,dimffnlkq,&
&     ekb_typ,enlout,ffnlk,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&     indlmn_typ,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&     lambda,lmnmax,matblk_der,mgfft,mpi_enreg,mpsang,mpssoang,natom_der,nattyp_der,&
&     gs_hamkq%ngfft,nkpg,nkpg1,nloalg_der,nnlout,npw,npw1,nspinor,nspinor,ntypat_der,&
&     0,paw_opt,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,signs,&
&     sij_dum,svectout_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1_,&
&     use_gpu_cuda=gs_hamkq%use_gpu_cuda)
   end if

   ABI_DEALLOCATE(enlout)
   ABI_DEALLOCATE(ph1d_der)
   ABI_DEALLOCATE(ph3din)
   ABI_DEALLOCATE(ph3dout)
   if (gs_hamkq%usepaw==1) then
     ABI_DEALLOCATE(ffnlk_der)
     if (usecprj==1) then
       call cprj_free(cprj_der)
       ABI_DEALLOCATE(cprj_der)
     end if
     if (usecprj==0) then
       call cprj_free(cwaveprj_tmp)
       ABI_DEALLOCATE(cwaveprj_tmp)
     end if
     nullify(cwaveprj_ptr)
   end if

!  k-point perturbation
!  -------------------------------------------
 else if(ipert==natom+1.and.(optnl>0.or.sij_opt/=0))then

!  Remember, q=0, so can take all RF data...
   nnlout=1
   ABI_ALLOCATE(enlout,(nnlout))
   tim_nonlop=8
   signs=2; choice=5
   if (gs_hamkq%usepaw==1) then
     cpopt=-1+5*usecprj; paw_opt=1; if (sij_opt/=0) paw_opt=sij_opt+3
     call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj,gs_hamkq%dimekb1,gs_hamkq%dimekb2,&
&     dimffnl1,dimffnl1,gs_hamkq%ekb,enlout,ffnl1,ffnl1,gs_hamkq%gmet,&
&     gs_hamkq%gprimd,idir,gs_hamkq%indlmn,gs_hamkq%istwf_k,kg1_k,kg1_k,&
&     kpg1_k,kpg1_k,gs_hamkq%kpoint,gs_hamkq%kpoint,lambda,lmnmax,matblk,&
&     mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&     nkpg1,nkpg1,gs_hamkq%nloalg,nnlout,npw1,npw1,nspinor,nspinor,ntypat,0,paw_opt,&
&     gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,&
&     signs,gs_hamkq%sij,gs1c,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1_,&
&     use_gpu_cuda=gs_hamkq%use_gpu_cuda)
   else
     cpopt=-1 ; paw_opt=0
     call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj,gs_hamkq%dimekb1,gs_hamkq%dimekb2,&
&     dimffnl1,dimffnl1,gs_hamkq%ekb,enlout,ffnl1,ffnl1,gs_hamkq%gmet,&
&     gs_hamkq%gprimd,idir,gs_hamkq%indlmn,gs_hamkq%istwf_k,kg1_k,kg1_k,&
&     kpg1_k,kpg1_k,gs_hamkq%kpoint,gs_hamkq%kpoint,lambda,lmnmax,matblk,&
&     mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&     nkpg1,nkpg1,gs_hamkq%nloalg,nnlout,npw1,npw1,nspinor,nspinor,ntypat,0,paw_opt,&
&     gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,&
&     signs,sij_dum,svectout_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1_,&
&     use_gpu_cuda=gs_hamkq%use_gpu_cuda)
   end if
   ABI_DEALLOCATE(enlout)

!  Electric field perturbation without Berry phase
!  -------------------------------------------
 else if(ipert==natom+2.and.berryopt/=4.and.(optnl>0.or.sij_opt/=0))then
!  gvnl1 was already initialized in the calling routine, by reading a ddk file
!  It contains |i du^(0)/dk_band>

   if (gs_hamkq%usepaw==1) then
     nnlout=1
     ABI_ALLOCATE(enlout,(nnlout))
     if (usecprj==1) then
       cwaveprj_ptr => cwaveprj
     else
       ABI_ALLOCATE(cwaveprj_tmp,(natom,nspinor))
       call cprj_alloc(cwaveprj_tmp,1,dimcprj)
       cwaveprj_ptr => cwaveprj_tmp
     end if
     if (opt_gvnl1==2.and.optnl>=1) then

!      PAW: Compute application of S^(0) to ddk WF
       cpopt=-1 ; choice=1 ; paw_opt=3 ; signs=2
       ABI_ALLOCATE(work,(2,npw1*nspinor))
       call nonlop(gs_hamkq%atindx1,choice,cpopt,cprj_der,gs_hamkq%dimekb1,gs_hamkq%dimekb2,&
&       dimffnlk,dimffnlk,gs_hamkq%ekb,enlout,ffnlk,ffnlk,gs_hamkq%gmet,gs_hamkq%gprimd,0,&
&       gs_hamkq%indlmn,gs_hamkq%istwf_k,kg_k,kg_k,kpg_k,kpg_k,gs_hamkq%kpoint,gs_hamkq%kpoint,&
&       lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,&
&       gs_hamkq%ngfft,nkpg,nkpg,gs_hamkq%nloalg,nnlout,npw,npw,nspinor,nspinor,ntypat,&
&       0,paw_opt,gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,&
&       signs,gs_hamkq%sij,work,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,gvnl1_,vectout_dum,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(gvnl1_,work,npw1,nspinor)
       do ipw=1,npw1*nspinor
         gvnl1_(:,ipw)=work(:,ipw)
       end do
!      $OMP END PARALLEL DO

!      PAW: Compute part of H^(1) due to derivative of S
       cpopt=4*usecprj ; choice=51 ; paw_opt=3 ; signs=2
       call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj_ptr,gs_hamkq%dimekb1,gs_hamkq%dimekb2,&
&       dimffnl1,dimffnl1,gs_hamkq%ekb,enlout,ffnl1,ffnl1,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&       gs_hamkq%indlmn,gs_hamkq%istwf_k,kg1_k,kg1_k,kpg1_k,kpg1_k,gs_hamkq%kpoint,gs_hamkq%kpoint,&
&       lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,&
&       gs_hamkq%ngfft,nkpg1,nkpg1,gs_hamkq%nloalg,nnlout,npw1,npw1,nspinor,nspinor,ntypat,&
&       0,paw_opt,gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,signs,&
&       gs_hamkq%sij,work,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,vectout_dum,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(gvnl1_,work,npw1,nspinor)
!      Note the multiplication by i
       do ipw=1,npw1*nspinor
         gvnl1_(1,ipw)=gvnl1_(1,ipw)-work(2,ipw)
         gvnl1_(2,ipw)=gvnl1_(2,ipw)+work(1,ipw)
       end do
!      $OMP END PARALLEL DO

!      PAW: Compute part of H^(1) due to derivative of electric field part of Dij
       cpopt=2 ; choice=1 ; paw_opt=1 ; signs=2
       call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj_ptr,dime1kb,gs_hamkq%dimekb2,&
&       dimffnlk,dimffnlk,e1kbfr,enlout,ffnlk,ffnlk,gs_hamkq%gmet,gs_hamkq%gprimd,0,&
&       gs_hamkq%indlmn,gs_hamkq%istwf_k,kg_k,kg_k,kpg_k,kpg_k,gs_hamkq%kpoint,gs_hamkq%kpoint,&
&       lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,&
&       gs_hamkq%ngfft,nkpg,nkpg,gs_hamkq%nloalg,nnlout,npw,npw,nspinor,nspinor,ntypat,&
&       0,paw_opt,gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,signs,&
&       sij_dum,svectout_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,work,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(gvnl1_,work,npw1,nspinor)
       do ipw=1,npw1*nspinor
         gvnl1_(:,ipw)=gvnl1_(:,ipw)+work(:,ipw)
       end do
!      $OMP END PARALLEL DO
       ABI_DEALLOCATE(work)

     end if ! opt_gvnl1==2

!    PAW: Compute derivatives due to part of D_ij^(1) depending on VHxc^(1)
     if (optnl>=2) then
       ABI_ALLOCATE(gvnl2,(2,npw1*nspinor))
       cpopt=-1+3*usecprj;if (opt_gvnl1==2) cpopt=2
       choice=1 ; paw_opt=1 ; signs=2
       call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj_ptr,dime1kb,gs_hamkq%dimekb2,dimffnlk,dimffnlk,&
&       e1kbsc,enlout,ffnlk,ffnlk,gs_hamkq%gmet,gs_hamkq%gprimd,0,&
&       gs_hamkq%indlmn,gs_hamkq%istwf_k,kg_k,kg_k,kpg_k,kpg_k,gs_hamkq%kpoint,gs_hamkq%kpoint,&
&       lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,&
&       gs_hamkq%ngfft,nkpg,nkpg,gs_hamkq%nloalg,nnlout,npw,npw,nspinor,nspinor,ntypat,&
&       0,paw_opt,gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,signs,&
&       sij_dum,svectout_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl2,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     end if
     if (sij_opt==1) gs1c=zero
     ABI_DEALLOCATE(enlout)
     if (usecprj==0) then
       call cprj_free(cwaveprj_tmp)
       ABI_DEALLOCATE(cwaveprj_tmp)
     end if
     nullify(cwaveprj_ptr)
   end if  ! PAW

!  Magnetic field perturbation (mimics electric field for testing)
!  -------------------------------------------
 else if(ipert==natom+5.and.(optnl>0.or.sij_opt/=0)) then
!  gvnl1 was already initialized in the calling routine, by reading a ddk file
   if (sij_opt==1) gs1c=zero

!  Electric field perturbation with Berry phase
!  -------------------------------------------
 else if(ipert==natom+2.and.berryopt==4.and.(optnl>0.or.sij_opt/=0))then

   if (optnl>=1) then
     do ipw=1,npw1*nspinor
       gvnl1_(1,ipw)=-grad_berry(2,ipw)
       gvnl1_(2,ipw)= grad_berry(1,ipw)
     end do
   end if
   if (sij_opt==1) gs1c=zero

!  Strain perturbation
!  -------------------------------------------
 else if ((ipert==natom+3.or.ipert==natom+4).and.(optnl>0.or.sij_opt/=0)) then

   if (optnl>=1) then
     if(ipert==natom+3) then
       istr=idir
     else
       istr = idir+3
     end if

!    Remember, q=0, so can take all RF data)
!    Copied from d/dk above; changes may be needed; tim_nonlop may need changing
     signs=2 ; choice=3 ; nnlout=6
     cpopt=-1+4*usecprj 
     paw_opt=0 
     ABI_ALLOCATE(enlout,(nnlout))
     call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj,gs_hamkq%dimekb1,gs_hamkq%dimekb2,&
&     dimffnl1,dimffnl1,gs_hamkq%ekb,enlout,ffnl1,ffnl1,gs_hamkq%gmet,&
&     gs_hamkq%gprimd,istr,gs_hamkq%indlmn,gs_hamkq%istwf_k,kg1_k,kg1_k,&
&     kpg1_k,kpg1_k,gs_hamkq%kpoint,gs_hamkq%kpoint,lambda,lmnmax,matblk,&
&     mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&     nkpg1,nkpg1,gs_hamkq%nloalg,nnlout,npw1,npw1,nspinor,nspinor,ntypat,0,paw_opt,&
&     gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,&
&     signs,sij_dum,svectout_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1_,&
&     use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     ABI_DEALLOCATE(enlout)
   end if

   if (sij_opt==1) gs1c=zero

!  No local part
!  -------------------------------------------
 else if (usevnl>0.or.(sij_opt/=0)) then

   if (optnl>=1) then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(gvnl1_,npw1,nspinor)
     do ipw=1,npw1*nspinor
       gvnl1_(:,ipw)=zero
     end do
!    $OMP END PARALLEL DO
   end if
   if (sij_opt/=0) gs1c=zero

 end if

!======================================================================
!== Apply the 1st-order kinetic operator to the wavefunction
!== (add it to nl contribution)
!======================================================================

!Phonon perturbation
!or Electric field perturbation
!-------------------------------------------
!No kinetic contribution

!k-point perturbation
!or Strain perturbation
!-------------------------------------------
 has_kin=(ipert==natom+1.or.ipert==natom+3.or.ipert==natom+4)
 if (has_kin) then

!  Remember that npw=npw1 for ddk perturbation
   do ispinor=1,nspinor
!    $OMP PARALLEL DO PRIVATE(ipw,ipws) &
!    $OMP&SHARED(cwave,ispinor,gvnl1_,dkinpw,kinpw1,npw,nspinor)
     do ipw=1,npw
       ipws=ipw+npw*(ispinor-1)
       if(kinpw1(ipw)<huge(zero)*1.d-11)then
         gvnl1_(1,ipws)=gvnl1_(1,ipws)+dkinpw(ipw)*cwave(1,ipws)
         gvnl1_(2,ipws)=gvnl1_(2,ipws)+dkinpw(ipw)*cwave(2,ipws)
       else
         gvnl1_(1,ipws)=zero
         gvnl1_(2,ipws)=zero
       end if
     end do
!    $OMP END PARALLEL DO
   end do

 end if

!======================================================================
!== Sum contributions to get the application of H^(1) to the wf
!======================================================================
!Also filter the wavefunctions for large modified kinetic energy

!Add non-local+kinetic to local part
 if (optnl>=1.or.has_kin) then
   do ispinor=1,nspinor
     ipws=(ispinor-1)*npw1
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(gh1c,gvnl1_,kinpw1,ipws,npw1)
     do ipw=1+ipws,npw1+ipws
       if(kinpw1(ipw-ipws)<huge(zero)*1.d-11)then
         gh1c(1,ipw)=gh1c(1,ipw)+gvnl1_(1,ipw)
         gh1c(2,ipw)=gh1c(2,ipw)+gvnl1_(2,ipw)
       else
         gh1c(1,ipw)=zero
         gh1c(2,ipw)=zero
       end if
     end do
!    $OMP END PARALLEL DO
   end do
 end if

!PAW: add non-local part due to first order change of VHxc
 if (usevnl2) then
   do ispinor=1,nspinor
     ipws=(ispinor-1)*npw1
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(gh1c,gvnl2,kinpw1,ipws,npw1)
     do ipw=1+ipws,npw1+ipws
       if(kinpw1(ipw-ipws)<huge(zero)*1.d-11)then
         gh1c(1,ipw)=gh1c(1,ipw)+gvnl2(1,ipw)
         gh1c(2,ipw)=gh1c(2,ipw)+gvnl2(2,ipw)
       end if
     end do
!    $OMP END PARALLEL DO
   end do
   ABI_DEALLOCATE(gvnl2)
 end if

 if (usevnl==1) then
   nullify(gvnl1_)
 else
   ABI_DEALLOCATE(gvnl1_)
 end if

 if(prtvol<0)then
   call status(0,filstat,iexit,level,'exit          ')
 end if

 call timab(196+tim_getgh1c,1,tsec)

end subroutine getgh1c
!!***
