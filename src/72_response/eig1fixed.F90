!{\src2tex{textfont=tt}}
!!****f* ABINIT/eig1fixed
!! NAME
!! eig1fixed
!!
!!
!! FUNCTION
!! Computes the fixed contributions to the diagonal first-order eigenvalues
!! eig1, consisting of the 1st-order nonlocal pseudopotential contribution
!! for phonon and strain perturbations, and the 1st-order kinetic contribution
!! for strain.  This routine is only called for metallic occupation and
!! for Q=0.  Its output is used to compute the first-order Fermi energy
!! which is required under these conditions.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (DRH, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  band=which particular band we are converging.
!!  cwave0(2,npw*nspinor)=GS wavefunction at k, in reciprocal space
!!  dimekb=first dimension of ekb (see ekb_typ)
!!  dimffnlk=second dimension of ffnl (1+number of derivatives)
!!  dimffnl1=second dimension of ffnl1 and ffnlkq (1+number of derivatives)
!!  dkinpw(npw)=derivative of the (modified) kinetic energy for each
!!    plane wave at k (Hartree)
!!  ekb_typ(dimekb,1,nspinor**2)=
!!  ->Norm conserving : (Real) Kleinman-Bylander energies (hartree)
!!                             for the displaced atom
!!          for number of basis functions (l,n) (lnmax)
!!          dimekb=lnmax
!!    ->PAW : (Real, symmetric) Frozen part of Dij coefficients
!!                               to connect projectors
!!                               for the displaced atom
!!          for number of basis functions (l,m,n) (lmnmax)
!!          dimekb=lmnmax*(lmnmax+1)/2
!!  ffnlk(npw,dimffnlk,lmnmax,1)=nonloc form factors at k, for the displaced atom.
!!  ffnlkq(npw1,dimffnl1,lmnmax,1)=nonloc form fact at k+q for the displaced atom
!!  ffnl1(npw1,dimffnl1,lmnmax,ntypat)=nonloc form factors at k+q
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
!!  lmnmax= max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband=number of bands.
!!  nkpg,nkpg1=second dimensions of kpg_k and kpg1_k (0 if useylm=0)
!!  npw=number of planewaves in basis sphere at given k.
!!  npw1=number of planewaves in basis sphere at k+q
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell.
!!  ph3d(2,npw1,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!
!! OUTPUT
!!  eig1_k(2*nband**2)=matrix of first-order eigenvalues (hartree)
!!   ( eig1(:,ii,jj)=<C0 ii|H1(nl+kin only)|C0 jj> ) (diagonal only)
!!  gvnl1(2,npw1*nspinor)=<G|Vnl1|C0 band,k> This is incidental and not
!!   used.  This array could have been allocated locally, but is already
!!   allocated in the calling routine.
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!
!! PARENTS
!!      wfkfermi3
!!
!! CHILDREN
!!      dotprod_g,nonlop,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eig1fixed(band,cwave0,dimekb,dimffnlk,dimffnl1,dkinpw,eig1_k,ekb_typ,&
& ffnlk,ffnlkq,ffnl1,gs_hamkq,gvnl1,idir,indlmn_typ,&
& ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lmnmax,matblk,mgfft,mpi_enreg,&
& mpsang,mpssoang,natom,nband,nkpg,nkpg1,npw,npw1,nspinor,ntypat,ph3d,prtvol)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eig1fixed'
 use interfaces_14_hidewrite
 use interfaces_53_spacepar
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,dimekb,dimffnl1,dimffnlk,idir,ipert,lmnmax,matblk
 integer,intent(in) :: mgfft,mpsang,mpssoang,natom,nband,nkpg,nkpg1,npw,npw1
 integer,intent(in) :: nspinor,ntypat,prtvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
!arrays
 integer,intent(in) :: indlmn_typ(6,lmnmax,1),kg1_k(3,npw1),kg_k(3,npw)
 real(dp),intent(in) :: dkinpw(npw),ekb_typ(dimekb,1,nspinor**2)
 real(dp),intent(in) :: ffnl1(npw1,dimffnl1,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlk(npw,dimffnlk,lmnmax,1)
 real(dp),intent(in) :: ffnlkq(npw1,dimffnl1,lmnmax,1),kinpw1(npw1)
 real(dp),intent(in) :: kpg1_k(npw1,nkpg1),kpg_k(npw,nkpg),kpt(3)
 real(dp),intent(inout) :: cwave0(2,npw*nspinor),ph3d(2,npw1,matblk)
 real(dp),intent(out) :: eig1_k(2*nband**2),gvnl1(2,npw1*nspinor)

!Local variables-------------------------------
!scalars
 integer,parameter :: test=0
 integer :: choice,cpopt,ipw,ipws
 integer :: ispinor,istr,istwf_k,matblk_der,n1,n2,n3,natom_der,nnlout
 integer :: ntypat_der,paw_opt,print,shift1,shift2,shift3,signs,tim_nonlop
 real(dp) :: arg,doti,dotr,dum
 character(len=500) :: message
!arrays
 integer :: atindx1_der(1),atindx_der(1),nattyp_der(1),nloalg_der(5)
 real(dp) :: dum_sij(1,1),dum_svectout(1,1),phkxredin(2,1),phkxredout(2,1),xred_der(3)
 real(dp),allocatable :: enlout(:),ph1d_der(:,:),ph3din(:,:,:),ph3dout(:,:,:)
 type(cprj_type) :: cprj_dum(1,1)

! *********************************************************************

!DEBUG Keep this debugging feature !
!write(std_out,*)' eig1fixed : enter'
!ENDDEBUG

!Tell us what is going on:
 if(prtvol>=10)then
   write(message, '(a,i6)' ) &
&   ' --- eig1fixed is called for band',band
   call  wrtout(std_out,message,'PERS')
 end if

 print=0
 n1=gs_hamkq%ngfft(1) ; n2=gs_hamkq%ngfft(2) ; n3=gs_hamkq%ngfft(3)

 istwf_k=gs_hamkq%istwf_k

!*************************************************************************
!
!Apply H(1) to phi(0)

!Phonon perturbation
 if(ipert <= natom) then

!  Compute dVnonlocal/d(atomic displacement) |C(n,k)>
   signs=2 ; choice=2 ; nnlout=3 ; natom_der=1 ; nattyp_der(1)=1 ; ntypat_der=1
   cpopt=-1
   paw_opt=0 
   ABI_ALLOCATE(enlout,(nnlout))
   matblk_der=1
   xred_der(:)=gs_hamkq%xred(:,ipert)
   atindx_der(1)=1
   atindx1_der(1)=1

!  Store at the right place the 1d phases
   ABI_ALLOCATE(ph1d_der,(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
   shift1=(gs_hamkq%atindx(ipert)-1)*(2*n1+1)
   ph1d_der(:,1:2*n1+1)=gs_hamkq%ph1d(:,1+shift1:2*n1+1+shift1)
   shift2=(gs_hamkq%atindx(ipert)-1)*(2*n2+1)+natom*(2*n1+1)
   ph1d_der(:,1+2*n1+1:2*n2+1+2*n1+1)=gs_hamkq%ph1d(:,1+shift2:2*n2+1+shift2)
   shift3=(gs_hamkq%atindx(ipert)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
   ph1d_der(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=&
&   gs_hamkq%ph1d(:,1+shift3:2*n3+1+shift3)

!  Will compute the 3D phase factors inside nonlop
   ABI_ALLOCATE(ph3din,(2,npw,1))
   ABI_ALLOCATE(ph3dout,(2,npw1,1))
   nloalg_der(:)=gs_hamkq%nloalg(:)
   nloalg_der(1)=-abs(gs_hamkq%nloalg(1))
   nloalg_der(4)=1

!  Compute here phkxred for kpt and kpq
   arg=two_pi*( kpt(1)*gs_hamkq%xred(1,ipert)+ &
&   kpt(2)*gs_hamkq%xred(2,ipert)+ &
&   kpt(3)*gs_hamkq%xred(3,ipert)   )
   phkxredin(1,1)=cos(arg)  ; phkxredin(2,1)=sin(arg)
   arg=two_pi * ( gs_hamkq%kpoint(1) * gs_hamkq%xred(1,ipert) +&
&   gs_hamkq%kpoint(2) * gs_hamkq%xred(2,ipert) +&
&   gs_hamkq%kpoint(3) * gs_hamkq%xred(3,ipert)   )
   phkxredout(1,1)=cos(arg) ; phkxredout(2,1)=sin(arg)

   tim_nonlop=7
   call nonlop(atindx1_der,choice,cpopt,cprj_dum,dimekb,ntypat_der,dimffnlk,dimffnl1,ekb_typ,&
&   enlout,ffnlk,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&   indlmn_typ,istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,dum,lmnmax,&
&   matblk_der,mgfft,mpi_enreg,mpsang,mpssoang,&
&   natom_der,nattyp_der,gs_hamkq%ngfft,nkpg,nkpg1,&
&   nloalg_der,nnlout,npw,npw1,nspinor,nspinor,ntypat_der,0,paw_opt,&
&   phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,&
&   signs,dum_sij,dum_svectout,&
&   tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave0,gvnl1,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
   ABI_DEALLOCATE(enlout)
   ABI_DEALLOCATE(ph1d_der)
   ABI_DEALLOCATE(ph3din)
   ABI_DEALLOCATE(ph3dout)

!  section for strain perturbation
 else if(ipert==natom+3 .or. ipert==natom+4)then
   if(ipert==natom+3) then
     istr=idir
   else
     istr = idir+3
   end if

!  Non-local (remember, q=0, so can take all RF data)
!  copied from d/dk above; changes may be needed; tim_nonlop may need changing
   signs=2 ; choice=3 ; nnlout=6 ; tim_nonlop=8
   cpopt=-1 
   paw_opt=0 
   ABI_ALLOCATE(enlout,(nnlout))
   call nonlop(gs_hamkq%atindx1,choice,cpopt,cprj_dum,&
&   gs_hamkq%dimekb1,gs_hamkq%dimekb2,dimffnl1,dimffnl1,gs_hamkq%ekb,&
&   enlout,ffnl1,ffnl1,gs_hamkq%gmet,gs_hamkq%gprimd,istr,&
&   gs_hamkq%indlmn,istwf_k,kg1_k,kg1_k,kpg1_k,kpg1_k,gs_hamkq%kpoint,gs_hamkq%kpoint,&
&   dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,&
&   gs_hamkq%ngfft,nkpg1,nkpg1,gs_hamkq%nloalg,nnlout,npw1,npw1,&
&   nspinor,nspinor,ntypat,0,paw_opt,gs_hamkq%phkxred,gs_hamkq%phkxred,&
&   gs_hamkq%ph1d,ph3d,ph3d,signs,&
&   dum_sij,dum_svectout,tim_nonlop,&
&   gs_hamkq%ucvol,gs_hamkq%useylm,cwave0,gvnl1,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
   ABI_DEALLOCATE(enlout)

!  Kinetic contribution. Remember that npw=npw1 for strain perturbation

   do ispinor=1,nspinor
!    $OMP PARALLEL DO PRIVATE(ipw,ipws) &
!    $OMP&SHARED(cwave0,ispinor,gvnl1,dkinpw,kinpw1,npw,nspinor)
     do ipw=1,npw
       ipws=ipw+npw*(ispinor-1)
       if(kinpw1(ipw)<huge(0.0_dp)*1.d-11)then
         gvnl1(1,ipws)=gvnl1(1,ipws)+dkinpw(ipw)*cwave0(1,ipws)
         gvnl1(2,ipws)=gvnl1(2,ipws)+dkinpw(ipw)*cwave0(2,ipws)
       else
         gvnl1(1,ipws)=0.0_dp
         gvnl1(2,ipws)=0.0_dp
       end if
     end do
!    $OMP END PARALLEL DO
   end do

!  end section for strain perturbation
 end if

!End of application of H(1) nonlocal + kinetic to phi(0)
!*************************************************************************

 call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw*nspinor,1,cwave0,gvnl1)

 eig1_k(2*band-1 +(band-1)*2*nband)=dotr

!DEBUG
!write(std_out,*)' eig1fixed : exit'
!ENDDEBUG

end subroutine eig1fixed
!!***
