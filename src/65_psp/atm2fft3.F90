!{\src2tex{textfont=tt}}
!!****f* ABINIT/atm2fft3
!! NAME
!! atm2fft3
!!
!! FUNCTION
!! This routine sums 1st-order atomic functions (density or potential)
!! defined (in rec. space) on a radial grid to get global 1st-order
!! quantities on the fine FFT grid.
!!
!! Possible options:
!!   optn=1: compute a sum of local 1st-order atomic densities
!!   optv=1: compute a sum of local 1st-order atomic potentials
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms ordered by type
!!  cplex: if 1, real space 1-order functions on FFT grid
!!  gauss(2,ntypat)= params for gaussian atm density (optn2=3) for each atom type
!!  gmet(3,3)=reciprocal space metric
!!  gsqcut=cutoff on |G|^2: see setup1 for definition (doubled sphere)
!!  idir=direction of atomic displacement (in case of phonons perturb.)
!!       used only if ndir=1 (see below)
!!  ipert=nindex of perturbation
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in unit cell.
!!  ndir=number of directions of atomic displacement:
!!       can be 1 (idir direction in then used) or 3 (all directions)
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  ntypat=number of types of atoms.
!!  optn,optn2,optv= (see NOTES below)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  typat(natom)=type of each atom
!!  ucvol=unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vspl(mqgrid,2,ntypat)=q^2 v(q) spline of an atomic potential
!!                        (used only if optv=1)
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  ======= if optv==1 =======
!!    atmvloc1(cplex*nfft)=sum of local 1st-order atomic potentials in real space
!!  ======= if optn==1 =======
!!   --- if optatm==1
!!    atmrho1(cplex*nfft)=sum of 1st-order atomic densities in real space
!!
!! NOTES
!! Details on possible options:
!! ============================
!! optv: controls the computation of a local 1st-order potential as sum of atomic potentials
!!          Vloc(r)=Sum_R[V1^AT(r-R)]
!! optn: controls the computation of a 1st-order density as sum of atomic densities
!!          n(r)=Sum_R[n1^AT(r-R)]
!!          n^AT is stored in reciprocal space:
!!          if optn2=1: n^AT is the atomic PAW PS core density stored in array pawtab%tcorespl()
!!                   2: n^AT is the atomic PAW PS valence density stored in array pawtab%tvalespl()
!!                   3: n^AT is a gaussian density: n(g)=gauss(1,ityp)*exp[-(gauss(2,ityp)*G)^2]
!! Note: optv and optn can be activated together
!!
!! Typical uses:
!! =============
!! Computation of:
!!  - 1st-order local potential: optv=1
!!  - 1st-order PS core density: optn=1, optn2=1
!!
!! PARENTS
!!      dyxc13,loper3,nstpaw3,pawgrnl
!!
!! CHILDREN
!!      fourdp,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine atm2fft3(atindx,atmrho1,atmvloc1,cplex,gauss,gmet,gsqcut,idir,ipert,&
&                   mgfft,mpi_enreg,mqgrid,natom,ndir,nfft,ngfft,ntypat,optn,optn2,optv,&
&                   paral_kgb,pawtab,ph1d,qgrid,qphon,typat,ucvol,usepaw,vspl,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'atm2fft3'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mgfft,mqgrid,natom,ndir,nfft,ntypat
 integer,intent(in) :: optn,optn2,optv,paral_kgb,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),ngfft(18),typat(natom)
 real(dp),intent(in) :: gauss(2,ntypat*(optn2/3)),gmet(3,3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid),qphon(3)
 real(dp),intent(in) :: vspl(mqgrid,2,ntypat*optv),xred(3,natom)
 real(dp),intent(out) :: atmrho1(cplex*nfft,ndir*optn)
 real(dp),intent(out) :: atmvloc1(cplex*nfft,ndir*optv)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,iatm,iatom,id,id1,id2,id3
 integer :: ig1,ig1max,ig1min,ig2,ig2max,ig2min,ig3,ig3max,ig3min
 integer :: ii,itypat,jj,me_fft,me_g0,n1,n2,n3,nproc_fft,shift1,shift2,shift3
 logical :: qeq0,qeq05
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,alf2pi2,bb,cc,cutoff,dd,diff,dq,dq2div6,dqdiv6,dqm1,g1,g2,g3
 real(dp) :: gauss1,gauss2,gmag,gq1,gq2,gq3,gsq,gsquar,n_at,ph12i,ph12r,ph1i
 real(dp) :: ph1r,ph2i,ph2r,ph3i,ph3r,phim_igia,phqim,phqre,phre_igia,qxred2pi
 real(dp) :: sfi,sfqi,sfqr,sfr,v_at,xnorm
 character(len=500) :: msg
!arrays
 integer :: jdir(ndir)
 real(dp) :: gq(3)
 real(dp),allocatable :: workn(:,:,:),workv(:,:,:)
!no_abirules
!Define G^2 based on G space metric gmet.
 gsq(g1,g2,g3)=g1*g1*gmet(1,1)+g2*g2*gmet(2,2)+g3*g3*gmet(3,3) &
 &       +two*(g1*g2*gmet(1,2)+g2*g3*gmet(2,3)+g3*g1*gmet(3,1))

! *************************************************************************

 DBG_ENTER("COLL")

 if(ipert==natom+1.or.ipert==natom+2.or.ipert==natom+5.or.ipert==natom+6) then

!  (In case of d/dk or an electric/magnetic field)
   if (optn==1) atmrho1(1:cplex*nfft,1:ndir)=zero
   if (optv==1) atmvloc1(1:cplex*nfft,1:ndir)=zero

 else

!  Useful quantities
   iatom=ipert;iatm=atindx(iatom)
   itypat=typat(iatom)
   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   me_fft=ngfft(11)
   nproc_fft=ngfft(10)
   if (ndir==1) then
     jdir(1)=idir
   else
     do id=1,ndir
       jdir(id)=id
     end do
   end if

   qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)
   qeq05=(abs(abs(qphon(1))-half)<tol12.or. &
&   abs(abs(qphon(2))-half)<tol12.or. &
&   abs(abs(qphon(3))-half)<tol12)
   if (nproc_fft>1.and.qeq05) then
     msg='  not compatible with FFT parallelism'
     MSG_ERROR(msg)
   end if

   if (optn2==3)then
     gauss1=gauss(1,itypat)
     gauss2=gauss(2,itypat)
     alf2pi2=(two_pi*gauss2)**2
   end if

   dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
   dqm1=one/dq
   dqdiv6=dq/six
   dq2div6=dq**2/six
   cutoff=gsqcut*tolfix
   id1=n1/2+2
   id2=n2/2+2
   id3=n3/2+2
   shift1=1+n1+(iatm-1)*(2*n1+1)
   shift2=1+n2+(iatm-1)*(2*n2+1)+natom*(2*n1+1)
   shift3=1+n3+(iatm-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
   ig1max=-1;ig2max=-1;ig3max=-1
   ig1min=n1;ig2min=n2;ig3min=n3

!  Determination of phase qxred*
   qxred2pi=two_pi*(qphon(1)*xred(1,iatom)+ &
&   qphon(2)*xred(2,iatom)+ &
&   qphon(3)*xred(3,iatom) )
   phqre=cos(qxred2pi)
   phqim=sin(qxred2pi)

!  Zero out temporary arrays
   if (optv==1) then
     ABI_ALLOCATE(workv,(2,nfft,ndir))
     workv(:,:,:)=zero
   end if
   if (optn==1) then
     ABI_ALLOCATE(workn,(2,nfft,ndir))
     workn(:,:,:)=zero
   end if

   ii=0
   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     gq3=dble(ig3)+qphon(3)
     gq(3)=gq3

     do i2=1,n2
       if (((i2-1)/(n2/nproc_fft))==me_fft) then
         ig2=i2-(i2/id2)*n2-1
         gq2=dble(ig2)+qphon(2)
         gq(2)=gq2

         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1
           gq1=dble(ig1)+qphon(1)
           gq(1)=gq1

           ii=ii+1
           gsquar=gsq(gq1,gq2,gq3)

!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then

!            Identify min/max indexes (to cancel unbalanced contributions later)
             if (qeq05) then
               ig1max=max(ig1max,ig1);ig1min=min(ig1min,ig1)
               ig2max=max(ig2max,ig2);ig2min=min(ig2min,ig2)
               ig3max=max(ig3max,ig3);ig3min=min(ig3min,ig3)
             end if

             gmag=sqrt(gsquar)
             me_g0=0;if (ig1==0.and.ig2==0.and.ig3==0.and.qeq0) me_g0=1

             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Compute structure factor
             ph1r=ph1d(re,ig1+shift1);ph1i=ph1d(im,ig1+shift1)
             ph2r=ph1d(re,ig2+shift2);ph2i=ph1d(im,ig2+shift2)
             ph3r=ph1d(re,ig3+shift3);ph3i=ph1d(im,ig3+shift3)
             ph12r=ph1r*ph2r-ph1i*ph2i
             ph12i=ph1r*ph2i+ph1i*ph2r
             phre_igia=ph12r*ph3r-ph12i*ph3i
             phim_igia=ph12r*ph3i+ph12i*ph3r

!            Compute V^AT(g+q) and/or n^AT(g+q) for given type of atom
!            Evaluate spline fit: p. 86 Numerical Recipes, Press et al;
!            Note the error in book for sign of "aa" term in derivative.
             if (optv==1.or.optn2/=3) then
               bb = diff*dqm1
               aa = one-bb
               cc = aa*(aa**2-one)*dq2div6
               dd = bb*(bb**2-one)*dq2div6
             end if
             if (optv==1) then
               if (me_g0==1) then
                 v_at=zero
               else
                 v_at=(aa*vspl(jj,1,itypat)+bb*vspl(jj+1,1,itypat)+&
&                 cc*vspl(jj,2,itypat)+dd*vspl(jj+1,2,itypat)) &
&                 /gsquar
               end if
             end if
             if (optn==1) then
               if (optn2==1) then
                 n_at=(aa*pawtab(itypat)%tcorespl(jj,1)+bb*pawtab(itypat)%tcorespl(jj+1,1)+&
&                 cc*pawtab(itypat)%tcorespl(jj,2)+dd*pawtab(itypat)%tcorespl(jj+1,2))
               else if (optn2==2) then
                 n_at=(aa*pawtab(itypat)%tvalespl(jj,1)+bb*pawtab(itypat)%tvalespl(jj+1,1)+&
&                 cc*pawtab(itypat)%tvalespl(jj,2)+dd*pawtab(itypat)%tvalespl(jj+1,2))
               else if (optn2==3) then
                 n_at=gauss1*exp(-gsquar*alf2pi2)
               else
                 n_at=zero
               end if
             end if

             do id=1,ndir


!              Exp(-i.2pi.g.xred)  * -i.2pi.(g+q)
               sfr=-two_pi*gq(jdir(id))*phim_igia
               sfi=-two_pi*gq(jdir(id))*phre_igia

!              Exp(-i.2pi.q.xred)            => -i.2pi.(g+q).Exp(-i.2pi.(g+q).xred)
               sfqr= sfr*phqre+sfi*phqim
               sfqi=-sfr*phqim+sfi*phqre

!              Assemble 1st-order potential/density in g space
               if (optv==1) then
                 workv(re,ii,id)=sfqr*v_at
                 workv(im,ii,id)=sfqi*v_at
               end if
               if (optn==1) then
                 workn(re,ii,id)=sfqr*n_at
                 workn(im,ii,id)=sfqi*n_at
               end if

             end do !id

!            End skip G**2 outside cutoff
           end if

!          End loop on n1, n2, n3
         end do
       end if ! this plane is selected
     end do
   end do

!  Identify unbalanced g-vectors
   if (qeq05) then  !This doesn't work in parallel
     ig1=-1;if (mod(n1,2)==0) ig1=1+n1/2
     ig2=-1;if (mod(n2,2)==0) ig2=1+n2/2
     ig3=-1;if (mod(n3,2)==0) ig3=1+n3/2
     if (abs(abs(qphon(1))-half)<tol12) then
       if (abs(ig1min)<abs(ig1max)) ig1=abs(ig1max)
       if (abs(ig1min)>abs(ig1max)) ig1=n1-abs(ig1min)
     end if
     if (abs(abs(qphon(2))-half)<tol12) then
       if (abs(ig2min)<abs(ig2max)) ig2=abs(ig2max)
       if (abs(ig2min)>abs(ig2max)) ig2=n2-abs(ig2min)
     end if
     if (abs(abs(qphon(3))-half)<tol12) then
       if (abs(ig3min)<abs(ig3max)) ig3=abs(ig3max)
       if (abs(ig3min)>abs(ig3max)) ig3=n3-abs(ig3min)
     end if
   end if

!  Get 1st-order potential/density back to real space
!  Non-symetrized non-zero elements have to be nullified
!  Divide by unit cell volume
   xnorm=one/ucvol
   if (optv==1) then
     do id=1,ndir
!      Eliminate unbalanced g-vectors
       if (qeq0) then       !q=0
         call zerosym(workv(:,:,id),2,mpi_enreg,n1,n2,n3)
       else if (qeq05) then !q=1/2; this doesn't work in parallel
         call zerosym(workv(:,:,id),2,mpi_enreg,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3)
       end if
       call fourdp(cplex,workv(:,:,id),atmvloc1(:,id),1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       atmvloc1(:,id)=atmvloc1(:,id)*xnorm
     end do
     ABI_DEALLOCATE(workv)
   end if
   if (optn==1) then
     do id=1,ndir
!      Eliminate unbalanced g-vectors
       if (qeq0) then       !q=0
         call zerosym(workn(:,:,id),2,mpi_enreg,n1,n2,n3)
       else if (qeq05) then !q=1/2; this doesn't work in parallel
         call zerosym(workn(:,:,id),2,mpi_enreg,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3)
       end if
       call fourdp(cplex,workn(:,:,id),atmrho1(:,id),1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       atmrho1(:,id)=atmrho1(:,id)*xnorm
     end do
     ABI_DEALLOCATE(workn)
   end if

!  End the condition of non-electric-field
 end if

 DBG_EXIT("COLL")

end subroutine atm2fft3
!!***
