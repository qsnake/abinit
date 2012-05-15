!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp9in
!! NAME
!! psp9in
!!
!! FUNCTION
!! Initialize pspcod=9 (Pseudopotentials from the XML format):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG, AF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  lloc=angular momentum choice of local pseudopotential
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=maximum number of points in real space grid in the psp file
!!   angular momentum of nonlocal pseudopotential
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=dimension of q (or G) grid for arrays.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  optnlxccc=option for nl XC core correction (input variable)
!!  qgrid(mqgrid)=values of q (or |G|) on grid from 0 to qmax
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r} dr]$ (hartree)
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpsang)=number of projection functions for each angular momentum
!!  qchrg is not used, and could be suppressed later
!!  vlspl(mqgrid,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr)
!!  xccc1d(n1xccc,6)=1D core charge function and five derivatives, from psp file
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!      cc_derivatives,close_xml_t,open_xml_file,parse,psp5lo,psp5nl
!!      smoothvlocal,spline,splint,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psp9in(filpsp,ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&
&                  mmax,mpsang,mqgrid,nproj,n1xccc,qchrg,qgrid,&
&                  useylm,vlspl,xcccrc,xccc1d,zion,znucl)

 use m_profiling

 use defs_basis
 use m_splines
#if defined HAVE_TRIO_FOX
 use fox_sax
 use m_xml_pseudo_types
 use m_xml_pseudo
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp9in'
 use interfaces_14_hidewrite
 use interfaces_65_psp, except_this_one => psp9in
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,lmnmax,lnmax,mpsang,mqgrid,n1xccc
 integer,intent(in) :: useylm
 integer,intent(out) :: mmax
 real(dp),intent(in) :: zion,znucl
 real(dp),intent(out) :: epsatm,qchrg,xcccrc
 character(len=fnlen),intent(in) :: filpsp
!arrays
 integer,intent(out) :: indlmn(6,lmnmax),nproj(mpsang)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: ekb(lnmax),ffspl(mqgrid,2,lnmax),vlspl(mqgrid,2)
 real(dp),intent(out) :: xccc1d(n1xccc,6)

!Local variables-------------------------------
!no_abirules
#if defined HAVE_TRIO_FOX
 integer :: ipsang,ir,jj,jpsang,maxn_pots,mm,npts
 integer :: il,index
 real(dp) :: a1,ainc,al,amesh,cmax,dnorm,fchrg,phi,r2,ratio
 real(dp) :: rmax,rmin,scale,step,yp1,ypn,z
 logical :: vlocsmooth
 character(len=500) :: message
 real(dp),allocatable :: ekb_tmp(:),ffspl_tmp(:,:,:),funrad(:),rad(:)
 real(dp),allocatable :: ratm(:),rfhi(:),vloc(:),vps(:,:),vpspll(:,:)
 real(dp),allocatable ::       ff1(:),ff2(:)
 real(dp),allocatable :: vlocal(:),wfll(:,:),work_space(:),work_spl(:)
!integer :: ii,irad,i1xccc,lhigh,ll,mmax2,pspcod
!real(dp) :: al_announced,der1,dern,factor,ff1fact
!real(dp),allocatable ::       ff3(:),ff4(:)
!real(dp),allocatable :: gg(:),gg1(:),gg2(:),gg3(:),gg4(:),radbis(:)
!----------------------------------------------------------------------------
! real(dp) z    : Atomic number of the element.
! real(dp) rmin : Value to determine the innermost point of the mesh
!                 in the FHI code ( rfhi(1) = rmin/dble(atomic-number) )
! real(dp) rmax : Outermost point of the mesh in the FHI code.
! real(dp) ainc : Default value for the ratio for the geometrical sequence
!                 used to define the grid in the FHI code.
! real(dp) amesh: Actual value for the ratio for the geometrical sequence
!                 used to define the grid in the FHI code.
! integer  mmax : Maximum number of points in the radial grid.
! real(dp) rfhi : Points of the radial grid in the FHI code,
!                 rfhi(ir) = rfhi(ir-1) * amesh .
! real(dp) scale: Scale used to define the points of the logarithmic
!                 radial grid of the ATM3 code.
! real(dp) step : Step used to define the points of the logarithmic
!                 radial grid of the ATM3 code.
! integer  npts : Number of points in the logarithmic radial grid of
!                 the ATM3 code.
! real(dp) ratm : Points of the radial grid in the ATM code,
!                 ratm(ir) = scale *  [ exp( step*ir) - 1 ]
!
! real(dp) rad    : Radial grid (geometrical sequence).
! real(dp) vpspll : Semilocal component of the pseudopotential.
!                   Units in Hartress.
!                   First argument: index in the radial grid.
!                   Second argument: angular momentum.
! real(dp) wfll   : Radial part of the pseudowave functions.
!                   wfll = u_n,l (r) = 1/r R_n,l(r)
!                   Units: electrons/bohr**(3/2).
!                   Normalized to 1.
!                   First argument: index in the radial grid.
!                   Second argument: angular momentum.
!---------------------------------------------------------------------------
! Variables for the non-linear core corrections ------------------------
!integer              :: nptscore
!real(dp),allocatable :: radcore(:)
 real(dp),allocatable :: ratmcore(:)
!----------------------------------------------------------------------------
! real(dp) radcore : Normalized radial grid where the pseudo-core charge density
!                    and its derivatives are known.
!                    It varies from 0 to 1 (from o to xcccrc)
! real(dp) ratmcore: Normalized logarithmic grid for the pseudo-core charge.
!                    It varies from 0 to 1 (from o to xcccrc)
! xcccrc           : XC core correction cutoff radius (bohr)
!                    It is defined as the radius where the pseudo-core
!                    charge density becomes zero
!                    (here we hae set up a tolerance of 1.d-12).
!----------------------------------------------------------------------------
 integer                  :: iostat
 type(xml_t)              :: fxml
 type(pseudo_t), pointer  :: psxml
#endif

! ***************************************************************************

 if(.false.)write(std_out,*)filpsp ! Just to keep filpsp when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)lloc   ! Just to keep lloc when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)lmax   ! Just to keep lmax when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)qgrid  ! Just to keep qgrid when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)useylm ! Just to keep useylm when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)zion   ! Just to keep zion when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)znucl  ! Just to keep znucl when HAVE_TRIO_FOX is false
#if defined HAVE_TRIO_FOX
!DEBUG
!write(std_out,*)' psp9in : enter and stop '
!stop
!ENDDEBUG

 call open_xml_file(fxml,filpsp,iostat)
 if (iostat /=0) stop "Cannot open file"

 call parse(fxml,pcdata_chunk,startElement_handler=begin_element,endElement_handler=end_element)

 call close_xml_t(fxml)
 if (iostat /=0) stop "Cannot close file"

 psxml => pseudo

!Define the mesh used by the FHI pseudopotential code
!The atomic number of the element is read from the header of the XML file
!The default values for rmin, rmax, and  ainc are:
 rmin = 0.00625d0
 rmax = 80.0d0
 ainc = 1.0247d0
 z    = psxml%header%atomicnumber
!---

!Determine the maximum number of points in the grid ---
 amesh    = dble(ainc)
 a1       = log(amesh)
 cmax     = dble(rmax)/rmin
 mmax     = log(cmax*z)/a1
 if(mod(mmax,2) .eq. 0) mmax = mmax + 1
!---

!Allocate the radial grid in the FHI code, and define the points ---
 ABI_ALLOCATE( rfhi,(mmax))
 ABI_ALLOCATE( rad,(mmax))
 rfhi(1) = dble(rmin)/z
 do ir = 2, mmax
   rfhi(ir) = amesh*rfhi(ir-1)
!  write(std_out,'(i5,f20.12)')ir, rfhi(ir)
 end do
 rad(:) = rfhi(:)
!---

!Now, define the mesh used by ATM3 pseudopotential code.
!We will assume at the present stage, that all the radial grids
!for the different magnitudes are the same
 npts  = psxml%pot(1)%V%grid%npts + 1
 scale = psxml%pot(1)%V%grid%scale
 step  = psxml%pot(1)%V%grid%step

!Allocate the radial grid in the FHI code, and define the points ---
 ABI_ALLOCATE( ratm,(npts))
 do ir = 1, npts
   ratm(ir) = scale * (exp(step*(ir-1))-1)
 end do
!---

!Allocate the radial variables (semilocal pseudopotentials and
!pseudoatomic orbitals)
 ABI_ALLOCATE( vpspll,(mmax,mpsang))
 ABI_ALLOCATE(wfll,(mmax,mpsang))
 ABI_ALLOCATE( vps,(npts,0:lmax))

!Interpolate atomic pseudopotential for each l, filling up arrays vpspll
!Note: put each l into vpspll(:,l+1)
 do il = 1, mpsang

   ABI_ALLOCATE( funrad,(npts))
   ABI_ALLOCATE( ff2,(npts))

   funrad(2:npts) = psxml%pot(il)%V%data(1:npts-1)
   r2        = ratm(2) / ( ratm(3) - ratm(2) )
   funrad(1) = funrad(2) - r2 * ( funrad(3) - funrad(2) )
   yp1       = ( funrad(2) - funrad(1) ) / &
&   ( ratm(2) - ratm(1) )
   ypn       = ( funrad(npts) - funrad(npts-1) )/&
&   ( ratm(npts) - ratm(npts-1) )

!  Store the semilocal components of the pseudopotential in the variable vps,
!  required to compute the local part of the pseudopotential.
!  vps is read in r*V format. Here we transform it to V format
   vps(2:npts,il-1) = funrad(2:npts) / ratm(2:npts)
   vps(1,il-1)      = vps(2,il-1)

!  !DEBUG
!  write(std_out,*)
!  write(std_out,*)'# il = ', il-1
!  do ir = 2, npts
!  write(std_out,'(3f20.12)')ratm(ir), &
!  &          psxml%pot(il)%V%data(ir-1)/(ratm(ir)*2.0d0), funrad(ir)
!  end do
!  !ENDDEBUG

!  Be careful, the interpolation is made with funrad,
!  and not with vps.
!  funrad is still in r*V format. That is why we multiply
!  afterwards by (1/r) again.

   ABI_ALLOCATE( work_space,(npts))
   call spline ( ratm, funrad, npts, yp1, ypn, ff2 )
   ABI_DEALLOCATE( work_space)

   call splint( npts, ratm, funrad, ff2, mmax, rad, vpspll(:,il) )

!  Transform the pseudopotential from r*V format to V format,
!  and from Rydbergs to Hartress

   do ir = 1, mmax
     vpspll(ir,il) = vpspll(ir,il) / ( 2.0d0 * rad(ir) )
   end do

   ABI_DEALLOCATE( funrad)
   ABI_DEALLOCATE( ff2)
 end do

!!DEBUG
!do il = 1, mpsang
!write(std_out,*)
!write(std_out,*)'# il = ', il-1
!do ir = 1, mmax
!write(std_out,'(2f20.12)')rad(ir), vpspll(ir,il)
!end do
!end do
!!ENDDEBUG

!Interpolate the radial wave functions ---
 do il = 1, mpsang

   ABI_ALLOCATE( funrad,(npts))
   ABI_ALLOCATE( ff2,(npts))

   funrad(2:npts) = psxml%pswf(il)%V%data(1:npts-1)
   r2        = ratm(2) / ( ratm(3) - ratm(2) )
   funrad(1) = funrad(2) - r2 * ( funrad(3) - funrad(2) )
   yp1       = ( funrad(2) - funrad(1) ) / &
&   ( ratm(2) - ratm(1) )
   ypn       = ( funrad(npts) - funrad(npts-1) )/&
&   ( ratm(npts) - ratm(npts-1) )

!  !DEBUG
!  write(std_out,*)
!  write(std_out,*)'# il = ', il-1
!  do ir = 2, npts
!  write(std_out,'(3f20.12)')ratm(ir), psxml%pswf(il)%V%data(ir-1), funrad(ir)
!  end do
!  !ENDDEBUG

   ABI_ALLOCATE( work_space,(npts))
   call spline ( ratm, funrad, npts, yp1, ypn, ff2 )
   ABI_DEALLOCATE( work_space)

   call splint( npts, ratm, funrad, ff2, mmax, rad, wfll(:,il) )

   ABI_DEALLOCATE( funrad)
   ABI_DEALLOCATE( ff2)

!  Normalize the pseudowave functions to 1.
   dnorm = 0.0d0
   do ir = 2, mmax
     phi   = wfll(ir,il)
     dnorm = dnorm + phi**2 * ( rad(ir) - rad(ir-1) )
   end do
   dnorm = sqrt(dnorm)

   do ir = 1, mmax
     wfll(ir,il) = wfll(ir,il)/dnorm
   end do

 end do

!!DEBUG
!do il = 1, mpsang
!write(std_out,*)
!write(std_out,*)'# il = ', il-1
!do ir = 1, mmax
!write(std_out,'(2f20.12)')rad(ir), wfll(ir,il)
!end do
!end do
!stop
!!ENDDEBUG

!File format of formatted fhi psp input, as adapted for use
!by the ABINIT code (the 3 first lines
!have already been read in calling -pspatm- routine) :

!(1) title (character) line
!(2) znucl,zion,pspdat
!(3) pspcod,pspxc,lmax,lloc,mmax,r2well
!(4) rchrg,fchrg,qchrg
!Note : prior to version 2.2, this 4th line started with  4--  ,
!and no core-correction was available.
!(5)-(18) -empty-
!(19) mmax, amesh ( mesh increment r(m+1)/r(m) )
!Then, for ll=0,lmax :
!for  irad=1,mmax  : irad, r(irad), upsp(irad,ll), vpsp(irad,ll)

!Take care of the non-linear core corrections

 select case(psxml%header%core_corrections)
   case("yes")
!    In Abinit, at least for the Troullier-Martins pseudopotential,
!    the pseudocore charge density and its derivatives (xccc1d)
!    is introduced in a linear grid (radcore).
!    This grid is normalized, so the radial coordinates run between
!    from 0 and 1 (from 0 to xcccrc, where xcccrc is the radius
!    where the pseudo-core becomes zero).

!    Allocate some of the arrays
     ABI_ALLOCATE( ratmcore,(npts))
     ABI_ALLOCATE( funrad,(npts))
     ABI_ALLOCATE( ff1,(npts))
     ABI_ALLOCATE( ff2,(npts))
     ABI_ALLOCATE( work_space,(npts))

!    Store the value of the pseudo-core charge.
!    It is read in a logarithmic grid.
     funrad(2:npts) = psxml%core_charge%data(1:npts-1)
     funrad(2:npts) = funrad(2:npts)/(4.0d0*pi*ratm(2:npts)**2)
     r2        = ratm(2) / ( ratm(3) - ratm(2) )
     funrad(1) = funrad(2) - r2 * ( funrad(3) - funrad(2) )
!    yp1       = ( funrad(2) - funrad(1) ) / &
!    &                  ( ratm(2) - ratm(1) )
     yp1       = zero
     ypn       = ( funrad(npts) - funrad(npts-1) )/&
&     ( ratm(npts) - ratm(npts-1) )

!    determine xcccrc where the pseudocore becomes 0
     do jj=npts,1,-1
       if (funrad(jj) > 1.0d-40) then
         xcccrc=ratm(jj)
         exit
       end if
     end do

!    Find the first derivative of the pseudo-core charge
!    in the logarithmic grid.
     ff1(1)=yp1
     do jj=2,npts-1
       ff1(jj)=(funrad(jj+1)-funrad(jj-1))/(ratm(jj+1)-ratm(jj-1))
     end do
     ff1(npts)=zero ! Could try to set this to 0 as well

!    Find the second derivative of the pseudo-core charge
!    in the logarithmic grid.
!    Be careful, this is very noisy at r->0
     call spline ( ratm, funrad, npts, yp1, ypn, ff2 )

!    call cc_derivatives to get 3rd 4th and 5th derivatives,
!    and everything splined onto regular grid [0:xcccrc]
!    in xccc1d
     write(std_out,*) 'psp9in: about to call cc_derivatives'
     write(std_out,*) npts,n1xccc,xcccrc

     call cc_derivatives(ratm,funrad,ff1,ff2,npts,n1xccc,xcccrc,xccc1d)

!    DEBUG
!    do ii = 2, npts
!    write(std_out,'(3f20.12)')ratmcore(ii),funrad(ii), ff2(ii)
!    enddo
!    
!    write(std_out,*) 'psp9in NLCC data ', n1xccc
!    do ii = 1, n1xccc
!    write(std_out,'(7e20.8)')xcccrc*(ii-1.d0)/(n1xccc-1.d0),xccc1d(ii,1),xccc1d(ii,2),xccc1d(ii,3),xccc1d(ii,4),xccc1d(ii,5),xccc1d(ii,6)
!    enddo
!    ENDDEBUG

     ABI_DEALLOCATE( work_space)
     ABI_DEALLOCATE( ratmcore)
     ABI_DEALLOCATE( funrad)
     ABI_DEALLOCATE( ff1)
     ABI_DEALLOCATE( ff2)
   case("no")
     write(message, '(a)' ) '  No XC core correction.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     xcccrc = 0.0d0 ; fchrg = 0.0d0 ; qchrg = 0.0d0
 end select

!Deallocate the different parts of psxml
 maxn_pots=size(psxml%pswf)

 do il = 1, maxn_pots
   if(associated( psxml%pswf(il)%V%data))   ABI_DEALLOCATE (psxml%pswf(il)%V%data)
   if(associated( psxml%pot(il)%V%data))    ABI_DEALLOCATE (psxml%pot(il)%V%data)
 end do
 if(associated( psxml%core_charge%data))    ABI_DEALLOCATE (psxml%core_charge%data)
 if(associated (psxml%valence_charge%data)) ABI_DEALLOCATE (psxml%valence_charge%data)

!Compute in real(dp) al : the announced amesh is inaccurate.
 ratio=rad(mmax)/rad(1)
 al=log(ratio)/dble(mmax-1)

!DEBUG
!write(std_out,*)' psp6in : al ; al_announced =',al,al_announced
!allocate(radbis(mmax))
!write(std_out,*)' psp6in : lloc  ',lloc
!do ipsang=1,lmax+1
!write(std_out,*)' psp6in : ipsang  ',ipsang
!do irad=1,mmax
!write(std_out,*)irad,rad(irad),wfll(irad,ipsang),vpspll(irad,ipsang)
!end do
!end do
!deallocate(radbis)
!ENDDEBUG

!Define the local component of the pseudopotential
 ABI_ALLOCATE(vloc,(mmax))

 ABI_ALLOCATE( vlocal,(npts))
 vlocsmooth = .true.

 if ( vlocsmooth ) then
   call smoothvlocal( lmax, npts, scale, step, vlocal, vps, psxml%header%zval )
!  vlocal is the soft local pseudopotential used in Siesta.
!  vlocal is computed in Ry. We transform it into Hartrees.
   vlocal(1:npts) = vlocal(1:npts) / 2.0d0
   yp1  = ( vlocal(2) - vlocal(1) ) / &
&   ( ratm(2) - ratm(1) )
   ypn  = ( vlocal(npts) - vlocal(npts-1) )/&
&   ( ratm(npts) - ratm(npts-1) )

   ABI_ALLOCATE( work_space,(npts))
   ABI_ALLOCATE( ff2,(npts))

   call spline ( ratm, vlocal, npts, yp1, ypn, ff2 )

   call splint( npts, ratm, vlocal, ff2, mmax, rad, vloc )

!  do ir = 1, mmax
!  write(std_out,'(2f20.12)') rad(ir), vloc(ir)
!  end do

   ABI_DEALLOCATE( work_space)
   ABI_DEALLOCATE( ff2)
 else
!  vloc(:)=Vlocal(r), lloc=0, 1, or 2 or -1 for avg.
!  Copy appropriate nonlocal psp for use as local one
   vloc( 1:mmax ) = vpspll( 1:mmax , lloc+1 )
 end if

!--------------------------------------------------------------------
!Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform)
!to get q^2 V(q).

 call psp5lo(al,epsatm,mmax,mqgrid,qgrid,&
& vlspl(:,1),rad,vloc,yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_space,(mqgrid))
 ABI_ALLOCATE(work_spl,(mqgrid))
 call spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl)

!--------------------------------------------------------------------
!Take care of non-local part

 ABI_ALLOCATE(ekb_tmp,(mpsang))
 ABI_ALLOCATE(ffspl_tmp,(mqgrid,2,mpsang))

!Zero out all Kleinman-Bylander energies to initialize
 ekb_tmp(:)=0.0d0

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then
   write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 else

!  ----------------------------------------------------------------------
!  Compute KB form factors and fit splines

   call psp5nl(al,ekb_tmp,ffspl_tmp,lmax,mmax,mpsang,mqgrid,qgrid,rad,vloc,&
&   vpspll,wfll)

 end if

!Define the number of projector per angular momentum shell (1 by default)
 do ipsang = 1, lmax + 1
   nproj(ipsang) = 1
 end do

 jj=0;indlmn(:,:)=0
 do ipsang=1,lmax+1
!  nproj had been set at 1, by default
   if(abs(ekb_tmp(ipsang))<tol10)then
     nproj(ipsang)=0
   end if
!  Possible values for nproj in this routine : 0 or 1.
   if(nproj(ipsang)==1)then
     if (useylm==1) then
       do mm=1,2*ipsang-1
         index=jj+mm     ! XG020223 : Marc or Francois, I think this is erroneous :
!        it should be outside of the mm loop, or
!        index=jj+1
         indlmn(1,index)=ipsang-1
         indlmn(2,index)=mm-ipsang
         indlmn(3,index)=1
         indlmn(4,index)=mm+(ipsang-1)*(ipsang-1)
         indlmn(5,index)=jj+1
         indlmn(6,index)=1
       end do
       jj=jj+(2*ipsang-1)
     else
       index=jj+1
       indlmn(1,index)=ipsang-1
       indlmn(2,index)=0
       indlmn(3,index)=1
       indlmn(4,index)=ipsang+(ipsang-1)*(ipsang-1)
       indlmn(5,index)=jj+1
       indlmn(6,index)=1
       jj=jj+1
     end if
   end if
 end do

!Transfer ekb and ffspl to their definitive location
 jpsang=1
 do ipsang=1,lmax+1
   if(nproj(ipsang)/=0)then
     ekb(jpsang)=ekb_tmp(ipsang)
     ffspl(:,:,jpsang)=ffspl_tmp(:,:,ipsang)
     jpsang=jpsang+1
     if(jpsang>lnmax)then
       write(message,'(6a,2i6)') ch10,&
&       ' psp6in : BUG -',ch10,&
&       '  Problem with the dimension of the ekb and ffspl arrays.',ch10,&
&       '  ipsang,lnmax=',ipsang,lnmax
     end if
   end if
 end do

 ABI_DEALLOCATE(ekb_tmp)
 ABI_DEALLOCATE(ffspl_tmp)

!DEBUG
!write(std_out,*)' psp9in : enter '
!write(std_out,*)' psp9in : indlmn(1:6,jj)'
!do jj=1,lmnmax
!write(std_out,*)indlmn(1:6,jj)
!end do
!ENDDEBUG

 ABI_DEALLOCATE( vpspll)
 ABI_DEALLOCATE(rad)
 ABI_DEALLOCATE(vloc)
 ABI_DEALLOCATE(wfll)
 ABI_DEALLOCATE(rfhi)
 ABI_DEALLOCATE(ratm)
 ABI_DEALLOCATE(vps)

 ABI_DEALLOCATE(vlocal)

#else
!Initialize some arguments, for portability at compile time
 indlmn=0 ; mmax=0 ; nproj=0
 ekb=zero ; epsatm=zero ; ffspl=zero ; qchrg=zero ; vlspl=zero ; xcccrc=zero ; xccc1d=zero
#endif

end subroutine psp9in
!!***
