!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_efg
!! NAME
!! calc_efg
!!
!! FUNCTION
!! calculation and output of electric field gradient tensor at each atomic site
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JZ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3)=matrix relating cartesian coords to crystal coords in reciprocal space
!!  natom=number of atoms in cell.
!!  nfft=number of points on fft grid
!!  ngfft(18)=details of fft
!!  nspden=number of spin densities
!!  ntypat=number of atom types
!!  paral_kgb
!!  ptcharge(ntypat)=user input charges on atoms to make simple point charge calc
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  prtefg=1 to print summary output, 2 for detailed output
!!  quadmom(ntypat)=quadrupole moments in barns of different atomic nuclei
!!  rhor(nfft,nspden)=electron density on grid (strictly $\tilde{n}+\hat{n}$)
!!  rprimd(3,3)=matrix relating cartesian coordinates to crystal coordinates
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume in Bohr^3
!!  usepaw=1 if we are using PAW formalism, 0 else
!!  xred(3,natom)=vectors locating each atom in the unit cell, in crystal coords
!!  zion(ntypat)=net core charge on each type of atom
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      dsyev,gridgcart,leave_new,make_efg_el,make_efg_ion,make_efg_onsite
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine calc_efg(gprimd,natom,nfft,ngfft,nspden,ntypat,paral_kgb,&
&                    paw_an,pawang,pawrad,pawrhoij,pawtab,&
&                    ptcharge,prtefg,quadmom,rhor,rprimd,typat,ucvol,usepaw,xred,zion)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_efg'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => calc_efg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,paral_kgb,prtefg,usepaw
 real(dp),intent(in) :: ucvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),ptcharge(ntypat)
 real(dp),intent(in) :: quadmom(ntypat),rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: zion(ntypat)
 real(dp),intent(inout) :: xred(3,natom)
 type(paw_an_type),intent(in) :: paw_an(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: INFO,LDA,LWORK,N,iatom
 real(dp) :: cq,eta,vxx,vyy,vzz
 character(len=500) :: message
!arrays
 real(dp) :: eigval(3),matr(3,3),work(8)
 real(dp),allocatable :: efg(:,:,:),efg_el(:,:,:),efg_ion(:,:,:),efg_paw(:,:,:)
 real(dp),allocatable :: efg_point_charge(:,:,:),gcart(:,:,:,:)

! ************************************************************************

!DEBUG
!write(std_out,*)' calc_efg : enter'
!ENDDEBUG
!Compatibility tests
 if (usepaw /= 1) then
   write (message,'(4a)')' calc_efg : ERROR- ',ch10,&
&   ' usepaw /= 1 but EFG calculation requires PAW ',ch10
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 ABI_ALLOCATE(gcart,(ngfft(1),ngfft(2),ngfft(3),3))
 call gridgcart(gcart,gprimd,ngfft) ! obtain G vectors in cartesian coords on grid

 ABI_ALLOCATE(efg,(3,3,natom))
 ABI_ALLOCATE(efg_el,(3,3,natom))
 ABI_ALLOCATE(efg_ion,(3,3,natom))
 ABI_ALLOCATE(efg_paw,(3,3,natom))
 ABI_ALLOCATE(efg_point_charge,(3,3,natom))
 efg_el(:,:,:) = zero
 efg_ion(:,:,:) = zero
 efg_paw(:,:,:) = zero
 efg_point_charge(:,:,:) = zero

 call make_efg_el(efg_el,gcart,natom,nfft,ngfft,nspden,paral_kgb,rhor,rprimd,xred)
 call make_efg_ion(efg_ion,natom,ntypat,rprimd,typat,ucvol,xred,zion)
 call make_efg_onsite(efg_paw,natom,ntypat,paw_an,pawang,pawrhoij,pawrad,pawtab,typat)

!calculate efg due to pure point charges, as input in variable ptcharge(ntypat)
!note here all atoms of the same type will have the same valence; in the future this
!could be made more flexible by having ptcharge(natom) but that will require a slightly
!different version than the existing make_efg_ion routine
 if(prtefg > 2) then  
   call make_efg_ion(efg_point_charge,natom,ntypat,rprimd,typat,ucvol,xred,ptcharge)
 end if

 efg(:,:,:) = efg_el(:,:,:) + efg_ion(:,:,:) + efg_paw(:,:,:)

 write(message,'(a,a,a)' ) ch10,' Electric Field Gradient Calculation ',ch10
 call wrtout(ab_out,message,'COLL')

 LDA=3; LWORK=8;N=3 ! these parameters are needed for the LAPACK dsyev routine 
 do iatom = 1, natom
   matr(:,:) = efg(:,:,iatom)
   call dsyev('V','U',N,matr,LDA,eigval,work,LWORK,INFO) ! get eigenvalues and eigenvectors
   if (eigval(3) > abs(eigval(1)) ) then ! In NMR, the convention is that whatever component is
!    largest in magnitude is called Vzz, next comes Vxx, then Vyy
     vzz = eigval(3)
     vxx = eigval(1)
     vyy = eigval(2)
   else
     vzz = eigval(1)
     vxx = eigval(3)
     vyy = eigval(2)
   end if 
   if (abs(quadmom(typat(iatom))) > tol8 ) then ! only relevant when quadmom > 0 for a given atom
!    cq = (eQ)*Vzz/h, where Q is the electric quadrupole moment and Vzz is the largest in magnitude
!    principal component of the EFG tensor. Q is input in quadmom in barns, and Vzz is computed in atomic
!    units. The factor 2349647.81 Ha^{-1}Bohr^2 fm^{-2} sec^-1 converts from atomic units to frequency (see
!    http://www.ismar.org/ISMARpedia/index.php/Nuclear_Quadrupole_Resonance for discussion); we divide by
!    10^6 to convert to MHz from Hz and multiply by 100 to convert from fm^2 to Barns.
     cq = vzz*quadmom(typat(iatom))*2349647.81/1.0E4
     if(abs(cq) > tol6 )then ! if Cq is non-zero, eta is meaningful, otherwise it s numerical noise
       eta = abs(vxx - vyy)/abs(vzz)
     else 
       eta=zero
     end if
   else
     cq =zero
     eta =zero
   end if
!  we always write Cq and eta, these are the NMR observables
   write(message,'(a,i3,a,i3,a,f13.6,a,f13.6)') ' Atom ',iatom,', typat ',typat(iatom),': Cq = ',cq,' MHz     eta = ',eta
   call wrtout(ab_out,message,'COLL')
   if (prtefg > 1) then ! print detailed results on component EFG's
     write(message,'(a,a,f13.6,a,a,3f13.6)')ch10,'      efg eigval : ',eigval(1),ch10,&
&     '-         eigvec : ',matr(1,1),matr(2,1),matr(3,1)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,f13.6,a,a,3f13.6)')'      efg eigval : ',eigval(2),ch10,&
&     '-         eigvec : ',matr(1,2),matr(2,2),matr(3,2)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,f13.6,a,a,3f13.6)')'      efg eigval : ',eigval(3),ch10,&
&     '-         eigvec : ',matr(1,3),matr(2,3),matr(3,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a,3f13.6)')ch10,'      total efg : ',efg(1,1,iatom),efg(1,2,iatom),efg(1,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6)')'      total efg : ',efg(2,1,iatom),efg(2,2,iatom),efg(2,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6,a)')'      total efg : ',efg(3,1,iatom),efg(3,2,iatom),efg(3,3,iatom),ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a,3f13.6)')ch10,'      efg_el : ',efg_el(1,1,iatom),efg_el(1,2,iatom),efg_el(1,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6)')'      efg_el : ',efg_el(2,1,iatom),efg_el(2,2,iatom),efg_el(2,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6,a)')'      efg_el : ',efg_el(3,1,iatom),efg_el(3,2,iatom),efg_el(3,3,iatom),ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6)')'      efg_ion : ',efg_ion(1,1,iatom),efg_ion(1,2,iatom),efg_ion(1,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6)')'      efg_ion : ',efg_ion(2,1,iatom),efg_ion(2,2,iatom),efg_ion(2,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6,a)')'      efg_ion : ',efg_ion(3,1,iatom),efg_ion(3,2,iatom),efg_ion(3,3,iatom),ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6)')'      efg_paw : ',efg_paw(1,1,iatom),efg_paw(1,2,iatom),efg_paw(1,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6)')'      efg_paw : ',efg_paw(2,1,iatom),efg_paw(2,2,iatom),efg_paw(2,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6,a)')'      efg_paw : ',efg_paw(3,1,iatom),efg_paw(3,2,iatom),efg_paw(3,3,iatom),ch10
     call wrtout(ab_out,message,'COLL')
   end if
   if (prtefg > 2) then ! write output of pure pointcharge calculation
     matr(:,:) = efg_point_charge(:,:,iatom)
     call dsyev('V','U',N,matr,LDA,eigval,work,LWORK,INFO) ! get eigenvalues and eigenvectors
     if (eigval(3) > abs(eigval(1)) ) then ! In NMR, the convention is that whatever component is
!      largest in magnitude is called Vzz, next comes Vxx, then Vyy
       vzz = eigval(3)
       vxx = eigval(1)
       vyy = eigval(2)
     else
       vzz = eigval(1)
       vxx = eigval(3)
       vyy = eigval(2)
     end if 
     if (abs(quadmom(typat(iatom))) > tol8 ) then ! only relevant when quadmom > 0 for a given atom
!      cq = e2Qq/h, where Vzz = eq and quadmom = Q; the other factors convert from atomic units to MHz
       cq = vzz*quadmom(typat(iatom))*2349647.81/1.0E4
       if(abs(cq) > tol6 )then ! if Cq is non-zero, eta is meaningful, otherwise it s numerical noise
         eta = abs(vxx - vyy)/abs(vzz)
       else 
         eta=zero
       end if
     else
       cq =zero
       eta =zero
     end if
!    we always write Cq and eta, these are the NMR observables
     write(message,'(a,i3,a,i3,a,f13.6,a,f13.6)') ' Atom ',iatom,', typat ',typat(iatom),&
&     ': Point charge Cq = ',cq,' MHz     eta = ',eta
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a,f13.6,a,a,3f13.6)')ch10,'      point charge efg eigval : ',eigval(1),ch10,&
&     '-         eigvec : ',matr(1,1),matr(2,1),matr(3,1)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,f13.6,a,a,3f13.6)')'      point charge efg eigval : ',eigval(2),ch10,&
&     '-         eigvec : ',matr(1,2),matr(2,2),matr(3,2)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,f13.6,a,a,3f13.6)')'      point charge efg eigval : ',eigval(3),ch10,&
&     '-         eigvec : ',matr(1,3),matr(2,3),matr(3,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a,3f13.6)')ch10,'      point charge efg : ',efg_point_charge(1,1,iatom),&
&     efg_point_charge(1,2,iatom),efg_point_charge(1,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6)')'      point charge efg : ',efg_point_charge(2,1,iatom),&
&     efg_point_charge(2,2,iatom),efg_point_charge(2,3,iatom)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3f13.6,a)')'      point charge efg : ',efg_point_charge(3,1,iatom),&
&     efg_point_charge(3,2,iatom),efg_point_charge(3,3,iatom),ch10
     call wrtout(ab_out,message,'COLL')
   end if
 end do
 write(message,'(3a)')ch10,ch10,ch10
 call wrtout(ab_out,message,'COLL')

 ABI_DEALLOCATE(gcart)
 ABI_DEALLOCATE(efg)
 ABI_DEALLOCATE(efg_el)
 ABI_DEALLOCATE(efg_ion)
 ABI_DEALLOCATE(efg_paw)
 ABI_DEALLOCATE(efg_point_charge)

!DEBUG
!write(std_out,*)' calc_efg : exit '
!stop
!ENDDEBUG

 end subroutine calc_efg
!!***
