!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyxc13
!! NAME
!! dyxc13
!!
!!
!! FUNCTION
!! Compute 2nd-order non-linear xc core-correction (part1)
!! to the dynamical matrix.
!! In case of derivative with respect to k or
!! electric field perturbation, the 1st-order local potential vanishes.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  kxc(nfft,nkxc)=first-order derivative of the xc potential
!!   if(nkxc=1): kxc(:,1)=dvxc/d$\rho$
!!   if(nkxc=3): kxc(:,1)=dvxc($\uparrow$)/d$\rho(\uparrow)$,
!!               kxc(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!               kxc(:,3)=dvxc($\downarrow$)/d$\rho(\downarrow)$
!!   if(nkxc=23): GGA case, see rhohxc_coll.f
!!  mgfft=maximum size of 1D FFTs
!!  mpert=maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(3)=fft grid dimensions.
!!  nkxc=second dimension of the kxc array
!!   (=1 for non-spin-polarized case, =3 for spin-polarized case)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  rfdir(3)=array that define the directions of perturbations
!!  rfpert(mpert)=array defining the type of perturbations that have to be computed
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xred(3,natom)=fractional coordinates for atoms in unit cell
!!
!! OUTPUT
!!  blkflgfrx1(3,natom,3,natom)=flag to indicate whether an element has been computed or not
!!  dyfrx1(2,3,natom,3,natom)=2nd-order non-linear xc
!!    core-correction (part1) part of the dynamical matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      atm2fft3,dotprod_vn,mkcor3,mkvxc3,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dyxc13(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,kxc,mgfft,mpert,mpi_enreg,mqgrid,&
&          natom,nfft,ngfft,nkxc,nspden,ntypat,n1xccc,paral_kgb,pawtab,&
&          ph1d,qgrid,qphon,rfdir,rfpert,rprimd,timrev,typat,ucvol,usepaw,xcccrc,xccc1d,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dyxc13'
 use interfaces_18_timing
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_65_psp
 use interfaces_72_response, except_this_one => dyxc13
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mpert,mqgrid,n1xccc,natom,nfft,nkxc,nspden,ntypat
 integer,intent(in) :: paral_kgb,timrev,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),ngfft(18),rfdir(3),rfpert(mpert),typat(natom)
 real(dp),intent(in) :: gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid),qphon(3)
 real(dp),intent(in) :: rprimd(3,3),xccc1d(n1xccc,6,ntypat),xcccrc(ntypat)
 real(dp),intent(in) :: xred(3,natom)
 integer,intent(out) :: blkflgfrx1(3,natom,3,natom)
 real(dp),intent(out) :: dyfrx1(2,3,natom,3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom1,iatom2,idir1,idir2,ifft,n1,n2,n3,n3xccc,nfftot
 integer :: option,optn,optn2,optv,upperdir
 real(dp) :: valuei,valuer
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: dum_atmvloc1(:),dum_gauss(:),dum_vspl(:)
 real(dp),allocatable :: rhor1(:,:),vxc10(:,:),xcccwk1(:),xcccwk2(:)

! *********************************************************************

 call timab(182,1,tsec)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 nfftot=n1*n2*n3

!Zero out the output arrays :
 blkflgfrx1(:,:,:,:)=0
 dyfrx1(:,:,:,:,:)=zero

 cplex=2-timrev ; n3xccc=nfft
 ABI_ALLOCATE(vxc10,(cplex*nfft,nspden))

 optv=0;optn=1;optn2=1

!Loop on the perturbation j1
 do iatom1=1,natom
   do idir1=1,3

!    Compute the derivative of the core charge with respect to j1
     ABI_ALLOCATE(xcccwk1,(cplex*n3xccc))

!    PAW: 1st-order core charge in reciprocal space
     if (usepaw==1) then
       call atm2fft3(atindx,xcccwk1,dum_atmvloc1,cplex,dum_gauss,gmet,gsqcut,idir1,iatom1,&
&       mgfft,mpi_enreg,mqgrid,natom,1,nfft,ngfft,&
&       ntypat,optn,optn2,optv,paral_kgb,pawtab,ph1d,qgrid,&
&       qphon,typat,ucvol,usepaw,dum_vspl,xred)

!      Norm-conserving psp: 1st-order core charge in real space
     else
       call mkcor3(cplex,idir1,iatom1,natom,ntypat,n1,n1xccc,&
&       n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xcccwk1,xred)
     end if

!    Compute the corresponding potential
     option=0
     ABI_ALLOCATE(rhor1,(cplex*nfft,nspden))
     call mkvxc3(cplex,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,&
&     n3xccc,option,paral_kgb,qphon,rhor1,rprimd,vxc10,xcccwk1)
     ABI_DEALLOCATE(rhor1)
     ABI_DEALLOCATE(xcccwk1)

!    vxc10 will couple with xcccwk2, that behaves like
!    a total density (ispden=1). Only the spin-up + spin-down
!    average of vxc10 is needed.
     if (nspden/=1)then
       do ifft=1,cplex*nfft
         vxc10(ifft,1)=(vxc10(ifft,1)+vxc10(ifft,2))*half
       end do
     end if


!    Loop on the perturbation j2
     do iatom2=1,iatom1
       upperdir=3
       if(iatom1==iatom2)upperdir=idir1
       do idir2=1,upperdir
         if( (rfpert(iatom1)==1 .and. rfdir(idir1) == 1) .or. &
&         (rfpert(iatom2)==1 .and. rfdir(idir2) == 1)    )then

!          Compute the derivative of the core charge with respect to j2
           ABI_ALLOCATE(xcccwk2,(cplex*n3xccc))

!          PAW: 1st-order core charge in reciprocal space
           if (usepaw==1) then
             call atm2fft3(atindx,xcccwk2,dum_atmvloc1,cplex,dum_gauss,gmet,gsqcut,idir2,iatom2,&
&             mgfft,mpi_enreg,mqgrid,natom,1,nfft,ngfft,&
&             ntypat,optn,optn2,optv,paral_kgb,pawtab,ph1d,qgrid,&
&             qphon,typat,ucvol,usepaw,dum_vspl,xred)

!            Norm-conserving psp: 1st-order core charge in real space
           else
             call mkcor3(cplex,idir2,iatom2,natom,ntypat,n1,n1xccc,&
&             n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xcccwk2,xred)
           end if

!          Get the matrix element j1,j2

           call dotprod_vn(cplex,xcccwk2,valuer,valuei,mpi_enreg,nfft,nfftot,1,2,vxc10,ucvol)

           ABI_DEALLOCATE(xcccwk2)

           dyfrx1(1,idir1,iatom1,idir2,iatom2)= valuer
           dyfrx1(2,idir1,iatom1,idir2,iatom2)= valuei
           dyfrx1(1,idir2,iatom2,idir1,iatom1)= valuer
           dyfrx1(2,idir2,iatom2,idir1,iatom1)=-valuei
           blkflgfrx1(idir1,iatom1,idir2,iatom2)=1
           blkflgfrx1(idir2,iatom2,idir1,iatom1)=1
         end if
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(vxc10)

 call timab(182,2,tsec)

end subroutine dyxc13
!!***
