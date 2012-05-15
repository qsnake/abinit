!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkresi
!! NAME
!! mkresi
!!
!! FUNCTION
!! Make residuals from knowledge of wf in G space and application
!! of Hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mcg)=<G|Cnk>=Fourier coefficients of wavefunction
!!  dimffnl=second dimension of ffnl (1+number of derivatives)g
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  filstat=name for the status file
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point
!!  isppol isppol=1 for unpolarized, 2 for spin-polarized
!!  kg_k(3,npw)=planewave reduced coordinates in basis sphere.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mcg=second dimension of the cg array
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in unit cell.
!!  nband=number of bands involved in subspace matrix.
!!  npw=number of planewaves in basis sphere at this k point.
!!  nspinor=number of spinors (on current proc)
!!  ntypat=number of types of atoms.
!!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear)
!!  n4,n5,n6 used for dimensionning of vlocal
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vlocal(n4,n5,n6,nvloc)=local potential in real space, on the augmented fft grid
!!
!! OUTPUT
!!  eig_k(nband)$= \langle C_n \mid H \mid C_n \rangle $ for each band.
!!  resid_k(nband)=residual for each band
!!   $= \langle C_n \mid H H \mid C_n \rangle- \langle C_n \mid H \mid C_n \rangle^2 $.
!!
!! PARENTS
!!      energy
!!
!! CHILDREN
!!      dotprod_g,getghc,sqnorm_g,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkresi(cg,dimffnl,eig_k,ffnl,filstat,gs_hamk,icg,ikpt,isppol,kg_k,kinpw,lmnmax,&
&                 matblk,mcg,mgfft,mpi_enreg,mpsang,mpssoang,natom,nband,npw,nspinor,&
&                 ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,resid_k,usepaw,vlocal)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkresi'
 use interfaces_18_timing
 use interfaces_53_spacepar
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimffnl,icg,ikpt,isppol,lmnmax,matblk,mcg,mgfft,mpsang
 integer,intent(in) :: mpssoang,n4,n5,n6,natom,nband,npw,nspinor,ntypat,nvloc
 integer,intent(in) :: paral_kgb,prtvol,usepaw
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_hamk
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: cg(2,mcg),ffnl(npw,dimffnl,lmnmax,ntypat),kinpw(npw)
 real(dp),intent(inout) :: ph3d(2,npw,matblk),vlocal(n4,n5,n6,nvloc)
 real(dp),intent(out) :: eig_k(nband),resid_k(nband)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_getghc=3
 integer :: cpopt,iband,ipw,istwf_k
 real(dp) :: doti,dotr
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 type(cprj_type) :: cwaveprj(1,1)

! *************************************************************************

!DEBUG
!write(std_out,*)' mkresi : enter '
!if(.true.)stop
!ENDDEBUG

!Keep track of total time spent in mkresi
 call timab(13,1,tsec)

 istwf_k=gs_hamk%istwf_k

 do iband=1,nband

   if(mpi_enreg%paral_compil_kpt==1)then
     if (mpi_enreg%paralbd >1) then
!      Skip this band if not the proper processor
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= mpi_enreg%me) then
         cycle
       end if
     end if
   end if

!  DEBUG
!  write(std_out,*)' mkresi : before getghc '
!  if(.true.)stop
!  ENDDEBUG

   ABI_ALLOCATE(cwavef,(2,npw*nspinor))
   ABI_ALLOCATE(ghc,(2,npw*nspinor))
   ABI_ALLOCATE(gvnlc,(2,npw*nspinor))
   if (usepaw==1)  then
     ABI_ALLOCATE(gsc,(2,npw*nspinor))
   end if

!  cwavef(:,:)=cg(:,1+npw*nspinor*(iband-1)+icg:npw*nspinor*iband+icg)
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cg,cwavef,iband,icg,npw,nspinor)
   do ipw=1,npw*nspinor
     cwavef(1,ipw)=cg(1,ipw+(iband-1)*npw*nspinor+icg)
     cwavef(2,ipw)=cg(2,ipw+(iband-1)*npw*nspinor+icg)
   end do
!  $OMP END PARALLEL DO
   cpopt=-1
   call getghc(cpopt,cwavef,cwaveprj,dimffnl,ffnl,filstat,ghc,gsc,&
&   gs_hamk,gvnlc,kg_k,kinpw,dotr,lmnmax,&
&   matblk,mgfft,mpi_enreg,mpsang,mpssoang,&
&   natom,1,npw,nspinor,ntypat,nvloc,n4,n5,n6,&
&   paral_kgb,ph3d,prtvol,usepaw,tim_getghc,0,vlocal)
   ABI_DEALLOCATE(gvnlc)
!  DEBUG
!  write(std_out,*)' mkresi : after getghc '
!  if(.true.)stop
!  ENDDEBUG

!  Compute the residual, <Cn|(H-<Cn|H|Cn>)**2|Cn>:
!  First get eigenvalue <Cn|H|Cn>:
   call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw*nspinor,1,cwavef,ghc)
   eig_k(iband)=dotr

!  Next need <G|(H-<Cn|H|Cn>)|Cn> (in ghc):
!  ghc(:,:)=ghc(:,:)-eig_k(iband)*cwavef(:,:)
   if (usepaw==0) then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(cwavef,eig_k,ghc,iband,npw,nspinor)
     do ipw=1,npw*nspinor
       ghc(1,ipw)=ghc(1,ipw)-eig_k(iband)*cwavef(1,ipw)
       ghc(2,ipw)=ghc(2,ipw)-eig_k(iband)*cwavef(2,ipw)
     end do
!    $OMP END PARALLEL DO
   else
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(gsc,eig_k,ghc,iband,npw,nspinor)
     do ipw=1,npw*nspinor
       ghc(1,ipw)=ghc(1,ipw)-eig_k(iband)*gsc(1,ipw)
       ghc(2,ipw)=ghc(2,ipw)-eig_k(iband)*gsc(2,ipw)
     end do
!    $OMP END PARALLEL DO
   end if

!  Then simply square the result:
   call sqnorm_g(dotr,istwf_k,mpi_enreg,npw*nspinor,ghc)
   resid_k(iband)=dotr

   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(ghc)
   if (usepaw==1) ABI_DEALLOCATE(gsc)

 end do

 call timab(13,2,tsec)

end subroutine mkresi
!!***
