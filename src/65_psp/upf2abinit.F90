!{\src2tex{textfont=tt}}
!!****f* ABINIT/upf2abinit
!! NAME
!! upf2abinit
!!
!! FUNCTION
!!  This routine wraps a call to a PWSCF module, which reads in
!!  a UPF (PWSCF / Espresso) format pseudopotential, then transfers
!!  data to abinit internal variables.
!!  "UPF PWSCF format" (pspcod=11)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filpsp = name of file with UPF data
!!  psps = sturcture with global dimension data for pseudopotentials, header
!!     info ...
!!
!! OUTPUT
!!  pspxc = index of xc functional for this pseudo
!!  lmax_ = maximal angular momentum
!!  lloc = local component chosen for pseudopotential
!!  mmax = maximum number of points in real space radial grid
!!  znucl = charge of species nucleus
!!  zion = valence charge
!!  epsatm = integral of local potential - coulomb potential of zion
!!  xcccrc = radius for non linear core correction
!!  ekb(dimekb)= Kleinman Bylander energies, see pspatm.F90
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  nproj= number of projectors for each channel
!!  xccc1d(n1xccc*(1-usepaw),6)=1D core charge function and five derivatives,
!!                              from psp file (used in NC only)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!      cc_derivatives,nderiv,psp11lo,psp11nl,read_pseudo,smooth,spline
!!      symbol2znucl,upfxc2abi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine upf2abinit (filpsp, &
     &  znucl, zion, pspxc, lmax_, lloc, mmax, &
     &  psps, epsatm, xcccrc, indlmn, ekb, ffspl, nproj_l, &
     &  vlspl, xccc1d)

 use m_profiling

  use defs_basis
  use defs_datatypes
  use m_splines
  use pseudo_pwscf ! pwscf module with all data explicit!

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'upf2abinit'
 use interfaces_11_qespresso_ext
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_45_psp_parser
 use interfaces_65_psp, except_this_one => upf2abinit
!End of the abilint section

  implicit none

!Arguments -------------------------------

  character(len=fnlen), intent(in) :: filpsp
  type(pseudopotential_type),intent(in) :: psps
  ! used contents:
  !   psps%lmnmax
  !   psps%mqgrid_ff
  !   psps%mqgrid_vl
  !   psps%dimekb
  !   psps%n1xccc
  !   psps%qgrid_ff
  !   psps%qgrid_vl

  integer, intent(out) :: pspxc, lmax_, lloc, mmax
  real(dp), intent(out) :: znucl, zion
  real(dp), intent(out) :: epsatm, xcccrc
  !arrays
  integer, intent(out)  :: indlmn(6,psps%lmnmax)
  integer, intent(out)  :: nproj_l(psps%mpssoang)
  real(dp), intent(out) :: ekb(psps%dimekb)
  real(dp), intent(out) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  real(dp), intent(out) :: vlspl(psps%mqgrid_vl,2)
  real(dp), intent(out) :: xccc1d(psps%n1xccc,6)

!Local variables -------------------------

  integer :: ir, iproj, ll, iunit, ios
  real(dp) :: yp1, ypn, rcov, amu

  logical, allocatable :: found_l(:)
  real(dp), allocatable :: work_space(:),work_spl(:)
  real(dp), allocatable :: ff(:), ff1(:), ff2(:), rad_cc(:)
  
  ! ######### in module pseudo: ############
  !
  !  only npsx = 1 is used here
  !  grids are allocated for much larger fixed length (ndm=2000)
  !  number of species (6) and projectors (8) as well...
  !
  !  psd(npsx) = specied string
  !  pseudotype = uspp / nc string
  !  dft(npsx) = exchange correlation string (20 chars)
  !  lmax(npsx) = maximum l channel
  !  mesh(npsx) = number of points for local pot
  !  nbeta(npsx) = number of projectors (beta functions for uspp)
  !  nlcc(npsx) = flag for presence of NL core correction
  !  zp(npsx) = valence ionic charge
  !  r(ndm,npsx) = radial mesh
  !  rab(ndm,npsx) = dr / di for radial mesh
  !  rho_atc(ndm,npsx) = NLCC pseudocharge density
  !  vloc0(ndm,npsx) = local pseudopotential
  !  betar(ndm, nbrx, npsx) = projector functions in real space mesh
  !  lll(nbrx,npsx) = angular momentum channel for each projector
  !  ikk2(nbrx,npsx) = maximum index for each projector function
  !  dion(nbrx,nbrx,npsx) = dij or Kleinman Bylander energies
  !  
  !  ########  end description of pseudo module contents ##########  

! *********************************************************************

!call pwscf routine for reading in UPF
 iunit = 1063
 open (unit=iunit, file=filpsp, status='old',form='formatted',iostat=ios)
 if (ios.ne.0) stop 'error opening upf file '

!read in psp data to static data in pseudo module, for ipsx == 1
 call read_pseudo(1,iunit)
 close (iunit)

!convert to Ha units
 vloc0  = half * vloc0
!betar = half * betar ! ???
 dion   = half * dion

!if upf file is a USPP one, stop
 if (pseudotype == 'US') then
   stop ' upf2abinit: Error: USPP UPF files not supported'
 end if

!copy over to abinit internal arrays and vars
 call upfxc2abi(dft(1), pspxc)
 lmax_ = lmax(1)

!Check if the local component is one of the angular momentum channels
!effectively if one of the ll is absent from the NL projectors
 ABI_ALLOCATE(found_l,(0:lmax_))
 found_l = .true.
 do ll = 0, lmax_
   if (any(lll(1:nbeta(1),1) == ll)) then
     found_l(ll) = .false.
   end if
 end do
 if (count(found_l) /= 1) then
   lloc = -1
 else
   do ll = 0, lmax_
     if (found_l(ll)) then
       lloc = ll
       exit
     end if
   end do
 end if
 ABI_DEALLOCATE(found_l)
!FIXME: do something about lloc == -1

 call symbol2znucl(amu,rcov,psd(1),znucl) 
 zion = zp(1)
 mmax = mesh(1)

 call psp11lo(rab(1:mmax,1),epsatm,mmax,psps%mqgrid_vl,psps%qgrid_vl,&
& vlspl(:,1),r(1:mmax,1),vloc0,yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_space,(psps%mqgrid_vl))
 ABI_ALLOCATE(work_spl,(psps%mqgrid_vl))
 call spline (psps%qgrid_vl,vlspl(:,1),psps%mqgrid_vl,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl)
 
!this has to do the FT of the projectors to reciprocal space
 call psp11nl(ffspl, indlmn, mmax, psps%lnmax, psps%lmnmax, psps%mqgrid_ff, &
& nbeta(1), betar(1:mmax,1:nbeta(1),1), lll(1:nbeta(1),1), ikk2(1:nbeta(1),1), &
& psps%qgrid_ff, r(1:mmax,1), rab(1:mmax,1), psps%useylm)

 nproj_l = 0
 do iproj = 1, nbeta(1)
   ll = lll(iproj,1)
   nproj_l(ll+1) = nproj_l(ll+1) + 1
 end do

!shape = dimekb  vs. shape = n_proj
 do ll = 1, nbeta(1)
   ekb(ll) = dion(ll,ll,1)
 end do

 xcccrc = zero
 xccc1d = zero
!if we find a core density, do something about it
!rho_atc contains the nlcc density
!rho_at contains the total density
 if (nlcc(1)) then
   ABI_ALLOCATE(ff,(mmax))
   ABI_ALLOCATE(ff1,(mmax))
   ABI_ALLOCATE(ff2,(mmax))
   ff(1:mmax) = rho_atc(1:mmax,1) ! model core charge without derivative factor

   ff1 = zero
   call nderiv(one,ff,ff1,mmax,1) ! first derivative
   ff1(1:mmax) = ff1(1:mmax) / rab(1:mmax,1)
   call smooth(ff1, mmax, 15) ! run 15 iterations of smoothing

   ff2 = zero
   call nderiv(one, ff1, ff2, mmax, 1) ! second derivative
   ff2(1:mmax) = ff2(1:mmax) / rab(1:mmax,1)
   call smooth(ff2, mmax, 15) ! run 10 iterations of smoothing?

!  determine a good rchrg = xcccrc
   do ir = mmax, 1, -1
     if (abs(ff(ir)) > 1.e-6) then
       xcccrc = r(ir,1)
       exit
     end if
   end do
   ABI_ALLOCATE(rad_cc,(mmax))
   rad_cc = r(1:mmax,1)
   rad_cc(1) = zero ! force this so that the core charge covers whole spline interval.
   call cc_derivatives(rad_cc,ff,ff1,ff2,mmax,psps%n1xccc,xcccrc,xccc1d)
   ABI_DEALLOCATE(rad_cc)
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(ff1)
   ABI_DEALLOCATE(ff2)

 end if !if nlcc present

end subroutine upf2abinit 
!!***
