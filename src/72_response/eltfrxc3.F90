!{\src2tex{textfont=tt}}
!!****f* ABINIT/eltfrxc3
!! NAME
!! eltfrxc3
!!
!! FUNCTION
!! Compute the 2nd derivatives of exchange-correlation energy
!! with respect to all pairs of strain and strain-atomic displacement
!! for the frozen wavefunction contribution to the elastic
!! and internal strain tensors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  enxc=exchange and correlation energy (hartree)
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=number of fft grid points
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=2nd dimension of kxc
!!  nspden=number of spin components of rhor
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of xccc3d (0 if no core charge, nfft otherwise)
!!  rhor(nfft,nspden)=electron density in r space
!!   (if spin polarized, array contains total density in first half and
!!    spin-up density in second half)
!!   (for non-collinear magnetism, first element: total density,
!!    3 next ones: mx,my,mz)
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  vxc(nfft,nspden)=xc potential (spin up in first half and spin down in
!!   second half if nspden=2)
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  eltfrxc(6+3*natom,6) = xc frozen wavefunction contribution to the
!!   elastic tensor
!!
!! SIDE EFFECTS
!!
!! NOTES
!!      Much of the code in versions of this routine prior to 4.4.5
!!      has been transfered to its child eltxccore.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      dotprod_vn,eltxccore,metric,mkcor3,mkvxcstr3,redgr,timab,xcomm_init
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eltfrxc3(eltfrxc,enxc,kxc,mpi_enreg,natom,&
& nfft,ngfft,nkxc,nspden,ntypat,n1xccc,n3xccc,paral_kgb,rhor,rprimd,&
& typat,vxc,xcccrc,xccc1d,xccc3d,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eltfrxc3'
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
 use interfaces_72_response, except_this_one => eltfrxc3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1xccc,n3xccc,natom,nfft,nkxc,nspden,ntypat,paral_kgb
 real(dp),intent(in) :: enxc
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: kxc(nfft,nkxc),rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: vxc(nfft,nspden),xccc1d(n1xccc,6,ntypat),xccc3d(n3xccc)
 real(dp),intent(in) :: xcccrc(ntypat),xred(3,natom)
 real(dp),intent(out) :: eltfrxc(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: mshift=401
 integer :: cplex,fgga,idir,ierr,ifft,ii,ipert,is1,is2,ispden,ispden_c,jj,ka,kb
 integer :: kd,kg,n1,n2,n3,n3xccc_loc,nfftot,old_paral_level,option,spaceComm
 real(dp) :: d2eacc,d2ecdgs2,d2exdgs2,d2gsds1ds2,d2gstds1ds2,decdgs,dexdgs
 real(dp) :: dgsds10,dgsds20,dgstds10,dgstds20,spnorm,tmp0,tmp0t,ucvol
 real(dp) :: valuei
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp) :: gmet(3,3),gprimd(3,3),qphon(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: d2gm(:,:,:,:),dgm(:,:,:),eltfrxc_tmp(:,:)
 real(dp),allocatable :: rho0_redgr(:,:,:),vxc10(:,:),vxc10_core(:)
 real(dp),allocatable :: vxc1is_core(:),vxc_core(:),work(:),workgr(:,:)
 real(dp),allocatable :: xccc3d1(:)

! *************************************************************************
!Initialize variables
 cplex=1
 qphon(:)=zero
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)

!HACK - should be fixed globally
 if(n1xccc==0) then
   n3xccc_loc=0
 else
   n3xccc_loc=n3xccc
 end if

 if(nkxc==23) then
   fgga=1
 else
   fgga=0
 end if

 ABI_ALLOCATE(eltfrxc_tmp,(6+3*natom,6))
 ABI_ALLOCATE(vxc10,(nfft,nspden))
 ABI_ALLOCATE(xccc3d1,(cplex*nfft))

 if(n1xccc/=0) then
   ABI_ALLOCATE(vxc_core,(nfft))
   ABI_ALLOCATE(vxc10_core,(nfft))
   ABI_ALLOCATE(vxc1is_core,(nfft))

   if(nspden==1) then
     vxc_core(:)=vxc(:,1)
   else
     vxc_core(:)=0.5_dp*(vxc(:,1)+vxc(:,2))
   end if
 end if

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!For GGA case, prepare quantities needed to evaluate contributions
!arising from the strain dependence of the gradient operator itself

 if(fgga==1) then
   ABI_ALLOCATE(rho0_redgr,(3,nfft,nspden))
   ABI_ALLOCATE(work,(nfft))
   ABI_ALLOCATE(workgr,(nfft,3))

!  Set up metric tensor derivatives
   ABI_ALLOCATE(dgm,(3,3,6))
   ABI_ALLOCATE(d2gm,(3,3,6,6))
!  Loop over 2nd strain index
   do is2=1,6
     kg=idx(2*is2-1);kd=idx(2*is2)
     do jj = 1,3
       dgm(:,jj,is2)=-(gprimd(kg,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kg,jj))
     end do

!    Loop over 1st strain index
     do is1=1,6
       ka=idx(2*is1-1);kb=idx(2*is1)
       d2gm(:,:,is1,is2)=0._dp
       do jj = 1,3
         if(ka==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&         +gprimd(kb,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kb,jj)
         if(ka==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&         +gprimd(kb,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(kb,jj)
         if(kb==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&         +gprimd(ka,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(ka,jj)
         if(kb==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&         +gprimd(ka,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(ka,jj)
       end do
       d2gm(:,:,is1,is2)=0.5_dp*d2gm(:,:,is1,is2)
     end do
   end do

!  Compute the reduced gradients of the zero-order charge density.
!  Note that in the spin-polarized case, we are computing the reduced
!  gradients of 2 X the spin-up or spin-down charge.  This simplifies
!  subsequent code for the non-spin-polarized case.
   if(nspden==1) then
     work(:)=rhor(:,1)
   else
     work(:)=2.0_dp*rhor(:,2)
   end if
   if(n1xccc/=0) then
     work(:)=work(:)+xccc3d(:)
   end if
   call redgr (work,workgr,mpi_enreg,nfft,ngfft,paral_kgb)
   do ifft=1,nfft
     rho0_redgr(:,ifft,1)=workgr(ifft,:)
   end do
   if(nspden==2) then
     work(:)=2.0_dp*(rhor(:,1)-rhor(:,2))
     if(n1xccc/=0) then
       work(:)=work(:)+xccc3d(:)
     end if
     call redgr(work,workgr,mpi_enreg,nfft,ngfft,paral_kgb)
     do ifft=1,nfft
       rho0_redgr(:,ifft,2)=workgr(ifft,:)
     end do
   end if
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(workgr)
 end if !GGA


!Null the elastic tensor accumulator
 eltfrxc(:,:)=zero;eltfrxc_tmp(:,:)=zero

!Normalization factor
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 if(nspden==1) then
   spnorm=one
 else
   spnorm=half
 end if

!Big loop over 2nd strain index
 do is2=1,6

!  Translate strain index as needed by mkcor3 below.
   if(is2<=3) then
     ipert=natom+3
     idir=is2
   else
     ipert=natom+4
     idir=is2-3
   end if

!  Generate first-order core charge for is2 strain if core charges are present.
   if(n1xccc/=0)then
     call mkcor3(cplex,idir,ipert,natom,ntypat,n1,n1xccc,&
&     n2,n3,qphon,rprimd,typat,ucvol,&
&     xcccrc,xccc1d,xccc3d1,xred)
   else
     xccc3d1(:)=zero
   end if
!  Compute the first-order potentials.
!  Standard first-order potential for LDA and GGA with core charge
   if(fgga==0 .or. (fgga==1 .and. n1xccc/=0)) then
     option=0
     call mkvxcstr3(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,&
&     nkxc,nspden,n3xccc_loc,option,paral_kgb,qphon,rhor,rhor,rprimd,vxc10,xccc3d1)

     if(n1xccc/=0)then
       if(nspden==1) then
         vxc10_core(:)=vxc10(:,1)
         vxc1is_core(:)=vxc10(:,1)
       else
         vxc10_core(:)=0.5_dp*(vxc10(:,1)+vxc10(:,2))
         vxc1is_core(:)=0.5_dp*(vxc10(:,1)+vxc10(:,2))
       end if
     end if
   end if

!  For GGA, first-order potential with doubled gradient operator strain
!  derivative terms needed for elastic tensor but not internal strain.
   if(fgga==1) then
     option=2
     call mkvxcstr3(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,&
&     nkxc,nspden,n3xccc_loc,option,paral_kgb,qphon,rhor,rhor,rprimd,vxc10,xccc3d1)

     if(n1xccc/=0)then
       if(nspden==1) then
         vxc10_core(:)=vxc10(:,1)
       else
         vxc10_core(:)=0.5_dp*(vxc10(:,1)+vxc10(:,2))
       end if
     end if
   end if

!  Additional term for diagonal strains.
   if(is2<=3) then
     vxc10(:,:)=vxc10(:,:)+vxc(:,:)
     if(n1xccc/=0) then
       vxc10_core(:)=vxc10_core(:)+2.0_dp*vxc_core(:)
       vxc1is_core(:)=vxc1is_core(:)+vxc_core(:)
     end if
   end if

!  For GGA, compute the contributions from the strain derivatives acting
!  on the gradient operators.
   if(fgga==1) then
     do ispden=1,nspden
       ispden_c=nspden-ispden+1

       do ifft=1,nfft

!        Collect the needed derivatives of Exc.  The factors introduced
!        deal with the difference between density as used here and
!        spin density as used with these kxc terms in other contexts.
         dexdgs=0.25_dp*kxc(ifft,2+ispden)
         d2exdgs2=0.03125_dp*kxc(ifft,6+ispden)
         decdgs=0.125_dp*kxc(ifft,12)
         d2ecdgs2=0.015625_dp*kxc(ifft,15)

!        Loop over 1st strain index
         do is1=1,6

!          The notation here is .gs... for the derivatives of the squared-
!          gradient of (2X) each spin density, and .gst... for the total
!          density.  Note the hack that the the total density is given
!          by the same expression for either the non-polarized or spin-
!          polarized case, implemented with the "complementary" index ispden_c
!          in the expression for tmp0t below.
           dgsds10=zero;dgsds20=zero;d2gsds1ds2=zero
           dgstds10=zero;dgstds20=zero;d2gstds1ds2=zero
           do jj=1,3
             do ii=1,3
               tmp0=rho0_redgr(ii,ifft,ispden)*rho0_redgr(jj,ifft,ispden)

               tmp0t=(rho0_redgr(ii,ifft,ispden)+rho0_redgr(ii,ifft,ispden_c))&
&               *(rho0_redgr(jj,ifft,ispden)+rho0_redgr(jj,ifft,ispden_c))

               dgsds10=dgsds10+dgm(ii,jj,is1)*tmp0
               dgsds20=dgsds20+dgm(ii,jj,is2)*tmp0

               dgstds10=dgstds10+dgm(ii,jj,is1)*tmp0t
               dgstds20=dgstds20+dgm(ii,jj,is2)*tmp0t

               d2gsds1ds2=d2gsds1ds2+d2gm(ii,jj,is1,is2)*tmp0

               d2gstds1ds2=d2gstds1ds2+d2gm(ii,jj,is1,is2)*tmp0t
             end do
           end do

!          Volume derivative terms added
           if(is1<=3) then
             d2gsds1ds2=d2gsds1ds2+dgsds20
             d2gstds1ds2=d2gstds1ds2+dgstds20
           end if
           if(is2<=3) then
             d2gsds1ds2=d2gsds1ds2+dgsds10
             d2gstds1ds2=d2gstds1ds2+dgstds10
           end if

!          Add the gradient derivative terms to eltfrxc.

           eltfrxc(is1,is2)=eltfrxc(is1,is2)+spnorm*&
&           (d2exdgs2*(dgsds10*dgsds20)+ dexdgs*d2gsds1ds2&
&           +d2ecdgs2*(dgstds10*dgstds20)+ decdgs*d2gstds1ds2)
         end do !is1
       end do !ifft
     end do !ispden
   end if !GGA

!  Compute valence electron 1st-order charge contributions.  Recall that
!  the diagonal strain derivatives of the valence charge are minus the
!  zero-order density.  The explicit symmetrization avoids the need
!  to store vxc10 for strain indices other than is2.

   call dotprod_vn(1,rhor,d2eacc,valuei,mpi_enreg,nfft,nfftot,nspden,1,&
&   vxc10,ucvol)
   do is1=1,3
     eltfrxc_tmp(is1,is2)=eltfrxc_tmp(is1,is2)-0.5_dp*d2eacc
     eltfrxc_tmp(is2,is1)=eltfrxc_tmp(is2,is1)-0.5_dp*d2eacc
   end do

!  Compute additional core contributions from is1 perturbation
!  Internal strain terms calculated here.
   if(n1xccc/=0) then
     call eltxccore(eltfrxc,is2,natom,nfft,ntypat,&
&     n1,n1xccc,n2,n3,rprimd,typat,ucvol,vxc_core,vxc10_core,vxc1is_core,&
&     xcccrc,xccc1d,xred)
   end if

!  Additional term for diagonal strains
   if(is2<=3) then
     do is1=1,3
       eltfrxc_tmp(is1,is2)=eltfrxc_tmp(is1,is2)+enxc
     end do
   end if

 end do !is2 outermost strain loop

!XG 030920 MPIWF : should accumulate eltfrxc accross processors
!Init mpi_comm
 old_paral_level=mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm)
 call timab(48,1,tsec)
 call xsum_mpi(eltfrxc,spaceComm ,ierr)
 call timab(48,2,tsec)
 mpi_enreg%paral_level=old_paral_level


!Normalize accumulated 2nd derivatives

 eltfrxc(:,:)=eltfrxc_tmp(:,:)+eltfrxc(:,:)*ucvol/dble(nfftot)

 ABI_DEALLOCATE(eltfrxc_tmp)
 ABI_DEALLOCATE(vxc10)

 if(n1xccc/=0) then
   ABI_DEALLOCATE(xccc3d1)
   ABI_DEALLOCATE(vxc_core)
   ABI_DEALLOCATE(vxc10_core)
   ABI_DEALLOCATE(vxc1is_core)
 end if

 if(fgga==1) then
   ABI_DEALLOCATE(rho0_redgr)
   ABI_DEALLOCATE(dgm)
   ABI_DEALLOCATE(d2gm)
 end if

end subroutine eltfrxc3
!!***
