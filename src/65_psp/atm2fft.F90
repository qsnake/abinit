!{\src2tex{textfont=tt}}
!!****f* ABINIT/atm2fft
!! NAME
!! atm2fft
!!
!! FUNCTION
!! This routine sums atomic functions (density or potential) defined
!  (in rec. space) on a radial grid to get global quantities on the
!! fine FFT grid. It can also compute contribution to energy derivatives
!! of these atomic functions.
!!
!! Possible options:
!!   optn=1: compute a sum of local atomic densities or contrib. to energy derivatives
!!   optv=1: compute a sum of local atomic potentials or contrib. to energy derivatives
!!
!!   optatm =1: computes sum of atomic potentials/densities
!!   optgr  =1: computes contribution of atomic pot./dens. to forces
!!   optstr =1: computes contribution of atomic pot./dens. to stress tensor
!!   optdyfr=1: computes contribution of atomic pot./dens. to frozen part of dyn. matrix
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  eei=local pseudopotential part of total energy
!!      (used only if optv=1 and optstr=1)
!!  gauss(2,ntypat)= params for gaussian atm density (optn2=3) for each atom type
!!  gmet(3,3)=reciprocal space metric
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on |G|^2: see setup1 for definition (doubled sphere)
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  ntypat=number of types of atoms.
!!  optatm,optdyfr,optgr,optn,optn2,optstr,optv= (see NOTES below)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qprtrb(3)= integer wavevector of possible perturbing potential
!!             in basis of reciprocal lattice translations
!!  rhog(2,nfft)=electron density rho(G) in reciprocal space
!!               (used only if optv=1 and (optgr=1 or optstr=1 or optdyfr=1))
!!  ucvol=unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vspl(mqgrid,2,ntypat)=q^2 v(q) spline of an atomic potential
!!                        (used only if optv=1)
!!  vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!            perturbing potential is added of the form V(G)=(vprtrb(1)+I*vprtrb(2))/2
!!            at the values G=qprtrb and (vprtrb(1)-I*vprtrb(2))/2 at G=-qprtrb
!!  vg(2,nfft)= potential V(G) in reciprocal space
!!              (used only if optn=1 and (optgr=1 or optstr=1 or optdyfr=1))
!!
!! OUTPUT
!!  ======= if optv==1 =======
!!  ============================
!!   --- if optatm==1
!!    atmvloc(nfft)=sum of local atomic potentials in real space
!!   --- if optgr==1
!!    grv(3,natom)=contribution of atomic potentials to forces
!!   --- if optstr==1
!!    strv(6)=contribution of atomic potentials to stress tensor
!!            cart. coordinates, symmetric tensor, 6 comp. in order 11,22,33,32,31,21
!!   --- if optdyfr==1
!!    dyfrv(3,3,natom)=contribution of atomic potentials to frozen part of dyn. matrix
!!
!!  ======= if optn==1 =======
!!  ============================
!!   --- if optatm==1
!!    atmrho(nfft)=sum of atomic densities in real space
!!   --- if optgr==1
!!    grn(3,natom)=contribution of atomic densities to forces
!!   --- if optstr==1
!!    strn(6)=contribution of atomic densities to stress tensor
!!            cart. coordinates, symmetric tensor, 6 comp. in order 11,22,33,32,31,21
!!   --- if optdyfr==1
!!    dyfrn(3,3,natom)=contribution of atomic densities to frozen part of dyn. matrix
!!
!! NOTES
!! Details on possible options:
!! ============================
!! optv: controls the computation of a local potential as sum of atomic potentials
!!          Vloc(r)=Sum_R[V^AT(r-R)]
!!       or its contributions to energy derivatives, i.e. derivatives of Int[Vloc(r).rho(r).dr]
!!          rho(r) is stored in reciprocal space in array rhog()
!!          V^AT is stored in reciprocal space in array vspl (in practice vspl(q)=q^2.V^AT(q))
!!
!! optn: controls the computation of a density as sum of atomic densities
!!          n(r)=Sum_R[n^AT(r-R)]
!!       or its contributions to energy derivatives, i.e. derivatives of Int[n(r).V(r).dr]
!!          V(r) is stored in reciprocal space in array vg()
!!          n^AT is stored in reciprocal space:
!!          if optn2=1: n^AT is the atomic PAW PS core density stored in array pawtab%tcorespl()
!!                   2: n^AT is the atomic PAW PS valence density stored in array pawtab%tvalespl()
!!                   3: n^AT is a gaussian density: n(g)=gauss(1,ityp)*exp[-(gauss(2,ityp)*G)^2]
!! Note: optv and optn can be activated together
!!
!! Options controlling which contrib. to Etot derivatives are computed:
!!   optatm =1: computes Vloc(r) or n(r) as sum of atomic potentials/densities
!!   optgr  =1: computes contribution of atomic Vloc(r) or n(r) to forces
!!   optstr =1: computes contribution of atomic Vloc(r) or n(r) to stress tensor
!!   optdyfr=1: computes contribution of atomic Vloc(r) or n(r) to fr part of dyn. matrix
!! Note: optatm, optgr, optstr and optdyfr can be activated together
!!
!! Typical uses:
!! =============
!! Computation of:
!!  - local potential: optv=1, optatm=1
!!  - contrib. of local potential to Etot derivatives: optv=1, rhog=total valence density
!!                                                     optgr=1 or optstr=1 or optdyfr=1
!!  - PS core density: optn=1, optn2=1, optatm=1
!!  - contrib. of NLCC to Etot derivatives: optn=1, optn2=1, vg=XC potential
!!                                          optgr=1 or optstr=1 or optdyfr=1
!!  - sum of atomic valence densities: optn=1, optn2=2 or 3, optatm=1
!!  - correction of forces due to potential residual: optn=1, optn2=2 or 3, optgr=1
!!                                                    vg=potential residual
!!    etc...
!!
!! PARENTS
!!      dyfro3,extraprho,forces,fresidrsp,prcref,prcref_PMA,respfn,setvtr
!!      stress
!!
!! CHILDREN
!!      fourdp,timab,wrtout,xcomm_init,xsum_mpi,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine atm2fft(atindx1,atmrho,atmvloc,dyfrn,dyfrv,eei,gauss,gmet,gprimd,grn,grv,gsqcut,&
&                  mgfft,mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&
&                  optatm,optdyfr,optgr,optn,optn2,optstr,optv,paral_kgb,&
&                  pawtab,ph1d,qgrid,qprtrb,rhog,strn,strv,ucvol,usepaw,vg,vprtrb,vspl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'atm2fft'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,natom,nfft,ntypat,optatm,optdyfr,optgr,optn
 integer,intent(in) :: optn2,optstr,optv,paral_kgb,usepaw
 real(dp),intent(in) :: eei,gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(18),qprtrb(3)
 real(dp),intent(in) :: gauss(2,ntypat*(optn2/3)),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid)
 real(dp),intent(in) :: rhog(2,nfft*optv*max(optgr,optstr,optdyfr))
 real(dp),intent(in) :: vg(2,nfft*optn*max(optgr,optstr,optdyfr)),vprtrb(2)
 real(dp),intent(in) :: vspl(mqgrid,2,ntypat*optv)
 real(dp),intent(out) :: atmrho(nfft*optn),atmvloc(nfft*optv)
 real(dp),intent(out) :: dyfrn(3,3,natom*optn*optdyfr)
 real(dp),intent(out) :: dyfrv(3,3,natom*optv*optdyfr),grn(3,natom*optn*optgr)
 real(dp),intent(out) :: grv(3,natom*optv*optgr),strn(6*optn*optstr)
 real(dp),intent(out) :: strv(6*optv*optstr)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ierr,ig1,ig1_,ig2,ig2_,ig3,ig3_,ii
 integer :: itypat,jj,me_fft,me_g0,n1,n2,n3,nproc_fft
 integer :: old_paral_level,shift1,shift2,shift3,spaceComm=0
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,alf2pi2,bb,cc,cutoff,dbl_ig1,dbl_ig2,dbl_ig3,dd,diff,dn_at,dq
 real(dp) :: dq2div6,dqdiv6,dqm1,dv_at,ee,ff,gauss1,gauss2,gmag,gsq,gsquar,n_at
 real(dp) :: ph12i,ph12r,ph1i,ph1r,ph2i,ph2r,ph3i,ph3r,sfi,sfr,term,tmpni,tmpnr
 real(dp) :: tmpvi,tmpvr,v_at,xnorm
 character(len=500) :: message
!arrays
 real(dp) :: gcart(3),tsec(2)
 real(dp),allocatable :: dyfrn_indx(:,:,:),dyfrv_indx(:,:,:),grn_indx(:,:)
 real(dp),allocatable :: grv_indx(:,:),phim_igia(:),phre_igia(:),workn(:,:)
 real(dp),allocatable :: workv(:,:)
!no_abirules
!Define G^2 based on G space metric gmet.
 gsq(i1,i2,i3)=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+dble(i3*i3)*gmet(3,3) &
& +two*(dble(i1*i2)*gmet(1,2)+dble(i2*i3)*gmet(2,3)+dble(i3*i1)*gmet(3,1))

! *************************************************************************

!DEBUG
!write(std_out,*)' atm2fft : enter '
!ENDDEBUG

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)

!Zero out arrays to permit accumulation over atom types
 if (optv==1.and.optatm==1) then
   ABI_ALLOCATE(workv,(2,nfft))
   workv(:,:)=zero
 end if
 if (optn==1.and.optatm==1) then
   ABI_ALLOCATE(workn,(2,nfft))
   workn(:,:)=zero
 end if
 if (optv==1.and.optgr==1) then
   ABI_ALLOCATE(grv_indx,(3,natom))
   grv_indx(:,:)=zero
 end if
 if (optn==1.and.optgr==1) then
   ABI_ALLOCATE(grn_indx,(3,natom))
   grn_indx(:,:)=zero
 end if
 if (optv==1.and.optdyfr==1) then
   ABI_ALLOCATE(dyfrv_indx,(3,3,natom))
   dyfrv_indx(:,:,:)=zero
 end if
 if (optn==1.and.optdyfr==1) then
   ABI_ALLOCATE(dyfrn_indx,(3,3,natom))
   dyfrn_indx(:,:,:)=zero
 end if
 if (optv==1.and.optstr==1) strv(:)=zero
 if (optn==1.and.optstr==1) strn(:)=zero
!
 dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
 dqm1=1.0_dp/dq
 dqdiv6=dq/6.0_dp
 dq2div6=dq**2/6.0_dp
 cutoff=gsqcut*tolfix
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 ABI_ALLOCATE(phre_igia,(natom))
 ABI_ALLOCATE(phim_igia,(natom))

 ia1=1
 do itypat=1,ntypat
!  ia1,ia2 sets range of loop over atoms:
   ia2=ia1+nattyp(itypat)-1
   ii=0

   if (optn2==3)then
     gauss1=gauss(1,itypat)
     gauss2=gauss(2,itypat)
     alf2pi2=(two_pi*gauss2)**2
   end if

   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     ig3_=ig3;if (ig3_==(n3/2+1)) ig3_=0
     do i2=1,n2
       ig2=i2-(i2/id2)*n2-1
       ig2_=ig2;if (ig2_==(n2/2+1)) ig2_=0
       if (((i2-1)/(n2/nproc_fft))==me_fft) then
         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1
           ig1_=ig1;if (ig1_==(n1/2+1)) ig1_=0
           ii=ii+1
           gsquar=gsq(ig1,ig2,ig3)

!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then

             gmag=sqrt(gsquar)
             me_g0=0;if (ig1==0.and.ig2==0.and.ig3==0) me_g0=1

             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Compute structure factor for all atoms of given type:
             do ia=ia1,ia2
               shift1=1+n1+(ia-1)*(2*n1+1)
               shift2=1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1)
               shift3=1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
               ph1r=ph1d(1,ig1+shift1);ph1i=ph1d(2,ig1+shift1)
               ph2r=ph1d(1,ig2+shift2);ph2i=ph1d(2,ig2+shift2)
               ph3r=ph1d(1,ig3+shift3);ph3i=ph1d(2,ig3+shift3)
               ph12r=ph1r*ph2r-ph1i*ph2i
               ph12i=ph1r*ph2i+ph1i*ph2r
               phre_igia(ia)=ph12r*ph3r-ph12i*ph3i
               phim_igia(ia)=ph12r*ph3i+ph12i*ph3r
             end do

!            Assemble structure factors for this type of atom= sum[exp(-i.piG.R)]
             if (optatm==1.or.optstr==1) then
               sfr=zero;sfi=zero
               do ia=ia1,ia2
                 sfr=sfr+phre_igia(ia)
                 sfi=sfi-phim_igia(ia)
               end do
             end if

!            Compute V^AT(G) and/or n^AT(G) for given type of atom
!            Evaluate spline fit: p. 86 Numerical Recipes, Press et al;
!            NOTE: error in book for sign of "aa" term in derivative;
!            !           also see splfit routine.
             if (optv==1.or.optn2/=3) then
               bb = diff*dqm1
               aa = 1.0_dp-bb
               cc = aa*(aa**2-1.0_dp)*dq2div6
               dd = bb*(bb**2-1.0_dp)*dq2div6
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

!            Compute sum of local atomic potentials or densities
!            ---------------------------------------------------
             if(optatm==1) then
!              Accumulate V^AT(G)*SF(G) or n^AT(G)*SF(G)
               if (optv==1) then
                 workv(re,ii)=workv(re,ii)+sfr*v_at
                 workv(im,ii)=workv(im,ii)+sfi*v_at
               end if
               if (optn==1) then
                 workn(re,ii)=workn(re,ii)+sfr*n_at
                 workn(im,ii)=workn(im,ii)+sfi*n_at
               end if

!              Compute contrib. to forces and/or frozen part of dyn. matrix
!              ------------------------------------------------------------
             else if (optgr==1.or.optdyfr==1) then
               dbl_ig1=dble(ig1_);dbl_ig2=dble(ig2_);dbl_ig3=dble(ig3_)
!              Compute (2Pi)*V^AT(G)*rho(G) or (2Pi)*n^AT(G)*V(G)
               if (optv==1) then
                 tmpvr=(two_pi*v_at)*rhog(re,ii)
                 tmpvi=(two_pi*v_at)*rhog(im,ii)
               end if
               if (optn==1) then
                 tmpnr=(two_pi*n_at)*vg(re,ii)
                 tmpni=(two_pi*n_at)*vg(im,ii)
               end if
!              === contrib. to forces
               if (optgr==1) then
!                Accumulate -(2Pi.G)*V^AT(G)*rho(G)*SF(G)
!                or -(2Pi)*n^AT(G)*V(G)*SF(G) into forces
                 if (optv==1) then
                   do ia=ia1,ia2
                     term=tmpvi*phre_igia(ia)+tmpvr*phim_igia(ia)
                     grv_indx(1,ia)=grv_indx(1,ia)-dbl_ig1*term
                     grv_indx(2,ia)=grv_indx(2,ia)-dbl_ig2*term
                     grv_indx(3,ia)=grv_indx(3,ia)-dbl_ig3*term
                   end do
                 end if
                 if (optn==1) then
                   do ia=ia1,ia2
                     term=tmpni*phre_igia(ia)+tmpnr*phim_igia(ia)
                     grn_indx(1,ia)=grn_indx(1,ia)-dbl_ig1*term
                     grn_indx(2,ia)=grn_indx(2,ia)-dbl_ig2*term
                     grn_indx(3,ia)=grn_indx(3,ia)-dbl_ig3*term
                   end do
                 end if
               end if
!              === contrib. to frozen part of dyn. matrix
               if (optdyfr==1) then
!                Accumulate -(2Pi^2.Gi.Gj)*V^AT(G)*rho(G)*SF(G)
!                or -(2Pi^2.Gi.Gj)*n^AT(G)*V(G)*SF(G) into dyn. matrix
                 if (optv==1) then
                   do ia=ia1,ia2
                     term=two_pi*(tmpvr*phre_igia(ia)-tmpvi*phim_igia(ia))
                     dyfrv_indx(1,1,ia)=dyfrv_indx(1,1,ia)-dbl_ig1*dbl_ig1*term
                     dyfrv_indx(1,2,ia)=dyfrv_indx(1,2,ia)-dbl_ig1*dbl_ig2*term
                     dyfrv_indx(1,3,ia)=dyfrv_indx(1,3,ia)-dbl_ig1*dbl_ig3*term
                     dyfrv_indx(2,2,ia)=dyfrv_indx(2,2,ia)-dbl_ig2*dbl_ig2*term
                     dyfrv_indx(2,3,ia)=dyfrv_indx(2,3,ia)-dbl_ig2*dbl_ig3*term
                     dyfrv_indx(3,3,ia)=dyfrv_indx(3,3,ia)-dbl_ig3*dbl_ig3*term
                   end do
                 end if
                 if (optn==1) then
                   do ia=ia1,ia2
                     term=two_pi*(tmpnr*phre_igia(ia)-tmpni*phim_igia(ia))
                     dyfrn_indx(1,1,ia)=dyfrn_indx(1,1,ia)-dbl_ig1*dbl_ig1*term
                     dyfrn_indx(1,2,ia)=dyfrn_indx(1,2,ia)-dbl_ig1*dbl_ig2*term
                     dyfrn_indx(1,3,ia)=dyfrn_indx(1,3,ia)-dbl_ig1*dbl_ig3*term
                     dyfrn_indx(2,2,ia)=dyfrn_indx(2,2,ia)-dbl_ig2*dbl_ig2*term
                     dyfrn_indx(2,3,ia)=dyfrn_indx(2,3,ia)-dbl_ig2*dbl_ig3*term
                     dyfrn_indx(3,3,ia)=dyfrn_indx(3,3,ia)-dbl_ig3*dbl_ig3*term
                   end do
                 end if
               end if
             end if

!            Compute contrib. to stress tensor
!            ---------------------------------
             if (optstr==1) then
!              Get (dV^AT(q)/dq)/q and/or (dn^AT(q)/dq)/q:
!              Note: correction of Numerical Recipes sign error before (3._dp*aa**2-1._dp)
!              ee*dqm1 + ff*dqdiv6 is the best estimate of dV(q)/dq from splines
               if (optv==1) then
                 if (me_g0==1) then
                   dv_at=zero
                 else
                   ee=vspl(jj+1,1,itypat)-vspl(jj,1,itypat)
                   ff=(3._dp*bb**2-1._dp)*vspl(jj+1,2,itypat)&
&                   -(3._dp*aa**2-1._dp)*vspl(jj  ,2,itypat)
                   dv_at=((ee*dqm1+ff*dqdiv6)/gmag-2.0_dp*v_at)/gsquar
                 end if
               end if
               if (optn==1) then
                 if (me_g0==1) then
                   if (optn2==1) then
                     dn_at=pawtab(itypat)%dncdq0
                   else if (optn2==2) then
                     dn_at=pawtab(itypat)%dnvdq0
                   else if (optn2==3) then
                     dn_at=-two*gauss1*alf2pi2
                   end if
                 else
                   if (optn2==1) then
                     ee=pawtab(itypat)%tcorespl(jj+1,1)-pawtab(itypat)%tcorespl(jj,1)
                     ff=(3._dp*bb**2-1._dp)*pawtab(itypat)%tcorespl(jj+1,2) &
&                     -(3._dp*aa**2-1._dp)*pawtab(itypat)%tcorespl(jj,2)
                   else if (optn2==2) then
                     ee=pawtab(itypat)%tvalespl(jj+1,1)-pawtab(itypat)%tvalespl(jj,1)
                     ff=(3._dp*bb**2-1._dp)*pawtab(itypat)%tvalespl(jj+1,2) &
&                     -(3._dp*aa**2-1._dp)*pawtab(itypat)%tvalespl(jj,2)
                   else if (optn2==3) then
                     dn_at=-two*gauss1*alf2pi2*exp(-gsquar*alf2pi2)
                   else
                   end if
                   dn_at=(ee*dqm1+ff*dqdiv6)/gmag
                 end if
               end if
!              Compute G in cartesian coordinates
               gcart(1)=gprimd(1,1)*dble(ig1)+gprimd(1,2)*dble(ig2)+&
&               gprimd(1,3)*dble(ig3)
               gcart(2)=gprimd(2,1)*dble(ig1)+gprimd(2,2)*dble(ig2)+&
&               gprimd(2,3)*dble(ig3)
               gcart(3)=gprimd(3,1)*dble(ig1)+gprimd(3,2)*dble(ig2)+&
&               gprimd(3,3)*dble(ig3)
!              Accumulate -dV^AT/dG*rho(G)*SF(G)*Gi.Gj/G
!              or -dn^AT/dG*V(G)*SF(G)*Gi.Gj/G
!              into stress tensor
               if (optv==1) then
                 term=(rhog(re,ii)*sfr+rhog(im,ii)*sfi)*dv_at
                 strv(1)=strv(1)-term*gcart(1)*gcart(1)
                 strv(2)=strv(2)-term*gcart(2)*gcart(2)
                 strv(3)=strv(3)-term*gcart(3)*gcart(3)
                 strv(4)=strv(4)-term*gcart(3)*gcart(2)
                 strv(5)=strv(5)-term*gcart(3)*gcart(1)
                 strv(6)=strv(6)-term*gcart(2)*gcart(1)
               end if
               if (optn==1) then
                 term=(vg(re,ii)*sfr+vg(im,ii)*sfi)*dn_at
                 strn(1)=strn(1)-term*gcart(1)*gcart(1)
                 strn(2)=strn(2)-term*gcart(2)*gcart(2)
                 strn(3)=strn(3)-term*gcart(3)*gcart(3)
                 strn(4)=strn(4)-term*gcart(3)*gcart(2)
                 strn(5)=strn(5)-term*gcart(3)*gcart(1)
                 strn(6)=strn(6)-term*gcart(2)*gcart(1)
               end if
             end if

!            End skip G**2 outside cutoff:
           end if

!          End loop on n1, n2, n3
         end do
       end if ! this plane is for me_fft
     end do
   end do

!  Symmetrize the dynamical matrix with respect to indices
   if (optdyfr==1) then
     if (optv==1) then
       do ia=ia1,ia2
         dyfrv_indx(2,1,ia)=dyfrv_indx(1,2,ia)
         dyfrv_indx(3,1,ia)=dyfrv_indx(1,3,ia)
         dyfrv_indx(3,2,ia)=dyfrv_indx(2,3,ia)
       end do
     end if
     if (optn==1) then
       do ia=ia1,ia2
         dyfrn_indx(2,1,ia)=dyfrn_indx(1,2,ia)
         dyfrn_indx(3,1,ia)=dyfrn_indx(1,3,ia)
         dyfrn_indx(3,2,ia)=dyfrn_indx(2,3,ia)
       end do
     end if
   end if

   ia1=ia2+1

!  End loop on type of atoms
 end do

 ABI_DEALLOCATE(phre_igia)
 ABI_DEALLOCATE(phim_igia)

!Get local potential or density back to real space
 if(optatm==1)then
!  Allow for the addition of a perturbing potential
   if ((optv==1).and.(vprtrb(1)**2+vprtrb(2)**2) > 1.d-30) then
!    Find the linear indices which correspond with the input wavevector qprtrb
!    The double modulus handles both i>=n and i<0, mapping into [0,n-1];
!    then add 1 to get range [1,n] for each
     i3=1+mod(n3+mod(qprtrb(3),n3),n3)
     i2=1+mod(n2+mod(qprtrb(2),n2),n2)
     i1=1+mod(n1+mod(qprtrb(1),n1),n1)
!    Compute the linear index in the 3 dimensional array
     ii=i1+n1*((i2-me_fft*n2/nproc_fft-1)+(n2/nproc_fft)*(i3-1))
!    Add in the perturbation at G=qprtrb
     workv(re,ii)=workv(re,ii)+0.5_dp*vprtrb(1)
     workv(im,ii)=workv(im,ii)+0.5_dp*vprtrb(2)
!    Same thing for G=-qprtrb
     i3=1+mod(n3+mod(-qprtrb(3),n3),n3)
     i2=1+mod(n2+mod(-qprtrb(2),n2),n2)
     i1=1+mod(n1+mod(-qprtrb(1),n1),n1)
!    ii=i1+n1*((i2-1)+n2*(i3-1))
     workv(re,ii)=workv(re,ii)+0.5_dp*vprtrb(1)
     workv(im,ii)=workv(im,ii)-0.5_dp*vprtrb(2)
     write(message, '(a,1p,2e12.4,a,0p,3i4,a)' )&
&     ' atm2fft: perturbation of vprtrb=', vprtrb,&
&     ' and q=',qprtrb,' has been added'
     call wrtout(std_out,message,'COLL')
   end if
!  Non-symetrized non-zero elements have to be nullified
!  Transform back to real space
!  Divide by unit cell volume
   xnorm=1.0_dp/ucvol
   if (optv==1) then
     call zerosym(workv,2,mpi_enreg,n1,n2,n3)
     call fourdp(1,workv,atmvloc,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     atmvloc(:)=atmvloc(:)*xnorm
     ABI_DEALLOCATE(workv)
   end if
   if (optn==1) then
     call zerosym(workn,2,mpi_enreg,n1,n2,n3)
     call fourdp(1,workn,atmrho,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     atmrho(:)=atmrho(:)*xnorm
     ABI_DEALLOCATE(workn)
   end if
 end if

!Additional treatment in case of parallelization
 if ((mpi_enreg%paral_compil_fft==1).and.&
& (optgr==1.or.optstr==1.or.optdyfr==1)) then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
   call timab(48,1,tsec)
   if (optv==1) then
     if (optgr==1)then
       call xsum_mpi(grv_indx,spaceComm,ierr)
     end if
     if (optstr==1)then
       call xsum_mpi(strv,spaceComm,ierr)
     end if
     if (optdyfr==1)then
       call xsum_mpi(dyfrv_indx,spaceComm,ierr)
     end if
   end if
   if (optn==1) then
     if (optgr==1)then
       call xsum_mpi(grn_indx,spaceComm,ierr)
     end if
     if (optstr==1)then
       call xsum_mpi(strn,spaceComm,ierr)
     end if
     if (optdyfr==1)then
       call xsum_mpi(dyfrn_indx,spaceComm,ierr)
     end if
   end if
   call timab(48,2,tsec)
   mpi_enreg%paral_level=old_paral_level
 end if

!Forces: re-order atoms
 if (optgr==1) then
   if (optv==1) then
     do ia=1,natom
       grv(1:3,atindx1(ia))=grv_indx(1:3,ia)
     end do
     ABI_DEALLOCATE(grv_indx)
   end if
   if (optn==1) then
     do ia=1,natom
       grn(1:3,atindx1(ia))=grn_indx(1:3,ia)
     end do
     ABI_DEALLOCATE(grn_indx)
   end if
 end if

!Stress tensor: normalize and add term -eei/ucvol on diagonal
 if (optstr==1) then
   if (optv==1) then
     strv(1)=(strv(1)-eei)/ucvol
     strv(2)=(strv(2)-eei)/ucvol
     strv(3)=(strv(3)-eei)/ucvol
     strv(4)=strv(4)/ucvol
     strv(5)=strv(5)/ucvol
     strv(6)=strv(6)/ucvol
   end if
   if (optn==1) strn(:)=strn(:)/ucvol
 end if

!Dynamical matrix: re-order atoms
 if (optdyfr==1) then
   if (optv==1) then
     do ia=1,natom
       dyfrv(1:3,1:3,atindx1(ia))=dyfrv_indx(1:3,1:3,ia)
     end do
     ABI_DEALLOCATE(dyfrv_indx)
   end if
   if (optn==1) then
     do ia=1,natom
       dyfrn(1:3,1:3,atindx1(ia))=dyfrn_indx(1:3,1:3,ia)
     end do
     ABI_DEALLOCATE(dyfrn_indx)
   end if
 end if

!DEBUG
!write(std_out,*)' atm2fft : exit '
!ENDDEBUG

end subroutine atm2fft
!!***
