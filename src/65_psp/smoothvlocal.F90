!{\src2tex{textfont=tt}}
!!****f* ABINIT/smoothvlocal
!! NAME
!! smoothvlocal
!!
!! FUNCTION
!!  Constructs the local pseudopotential used by SIESTA. Different from
!!  the individual v_l components, it must be as smooth as possible in
!!  order to be well represented on a 3D real space grid in SIESTA.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      psp9in
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine smoothvlocal(  lmax, npts, scale, step, vlocal, vps, zval)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smoothvlocal'
 use interfaces_65_psp, except_this_one => smoothvlocal
!End of the abilint section

 implicit none

!Arguments --------------------------------
 integer, intent(in)     :: npts, lmax
 real(dp), intent(in)    :: zval, scale, step, vps(npts,0:lmax)
 real(dp), intent(out)   :: vlocal(npts)

!Local variables --------------------------
 real(dp)                :: rpb, ea, rgauss, rgauss2
 real(dp), allocatable   :: rofi(:), drdi(:), s(:), chlocal(:)

 integer :: ir, nrgauss, nchloc

! *************************************************************************

!Allocate radial functions -----
 ABI_ALLOCATE( rofi,(npts))
 ABI_ALLOCATE( drdi,(npts))
 ABI_ALLOCATE( s,(npts))
 ABI_ALLOCATE( chlocal,(npts))

 rpb = scale
 ea  = exp(step)
 do ir = 1, npts
   rofi(ir) = scale * ( exp( step*(ir-1) ) - 1 )
   drdi(ir) = step  * rpb
   s(ir)    = sqrt( step*rpb )
   rpb      = rpb * ea
!  write(std_out,'(i5,3f20.12)')ir, rofi(ir), drdi(ir), s(ir)
 end do

 call radii_ps( vps, rofi, zval, npts, lmax, nrgauss, rgauss, rgauss2)

!Calculate local pseudopotential
 if ( rgauss2 .gt. 1.30d0 * rgauss ) then
!  In this case the atom core is so big that we do not have an asymptotic
!  of 2*Zval/r until Rgauss2 > Rc . To retain the same asymptotic
!  behaviour as in the pseudopotentials we modified the definition
!  of the local potential
!  
   call vlocal2( zval, npts, step, rofi, drdi, s, vps(:,0),  &
&   nrgauss, vlocal, nchloc, chlocal )
 else
!  
!  In this case the pseudopotential reach to an asymptotic
!  behaviour 2*Zval/r for a radius approximately equal to Rc.
!  
   call vlocal1( zval, npts, step, rofi, drdi, s, rgauss, vlocal, &
&   nchloc, chlocal )
 end if

 ABI_DEALLOCATE( rofi)
 ABI_DEALLOCATE( drdi)
 ABI_DEALLOCATE( s)
 ABI_DEALLOCATE( chlocal)

end subroutine smoothvlocal
!!***

!!****f* ABINIT/vlocal2
!! NAME
!! vlocal2
!!
!! FUNCTION
!!     This routine generates the local pseudopotential appropriate
!!     for species with  a large core.
!!
!! NOTES
!!     Written by D. Sanchez-Portal, Aug. 1998
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      smoothvlocal
!!
!! CHILDREN
!!
!! SOURCE
subroutine vlocal2( zval, nrval, a, rofi, drdi, s, vps, nrgauss, &
&                    vlocal,nchloc,chlocal )

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vlocal2'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 real(dp), intent(in)    :: zval, a
 integer,  intent(in)    :: nrval
 integer,  intent(inout) :: nrgauss
 real(dp), intent(in)    :: rofi(:), drdi(:), s(:), vps(:)
 real(dp), intent(out)   :: vlocal(:), chlocal(:)
 integer,  intent(out)   :: nchloc

!Local variables -------------------------
 real(dp) ::  vlc, r, dev, dev2, dev3, var1, var2, var3, v1, v2, v3, v4, &
&             dm11, dm12, dm13, dm21, dm22, dm23, dm31, dm32, dm33,      &
&             g0, g1, g2, g3, g4, d2g, d2u, cons, a2b4, qtot
 integer  ::  ndevfit, ir

 real(dp), parameter  :: eps=1.0d-5

! *********************************************************************

!Continuity up to second derivative***
 ndevfit=2
!Continuity up to third derivative****
!ndevfit=3

 nrgauss = nrgauss + 3        !! For good measure...

 do ir = 1, nrval
   vlocal(ir) = vps(ir) * rofi(ir)
 end do

 ir   = nrgauss
 dev  = ( vlocal(ir+1) - vlocal(ir-1) ) * 0.5d0
 dev2 = ( vlocal(ir+1) + vlocal(ir-1) - 2.0d0 * vlocal(ir) )
 dev3 = ( vlocal(ir+2) - 2.0d0 * vlocal(ir+1)                  &
& + 2.0d0 * vlocal(ir-1) - vlocal(ir-2) ) * 0.5d0
 dev3 = ( dev3 - 3.0d0 * a * dev2 + 2.0d0 * (a**2) * dev ) / ( drdi(ir)**3 )
 dev2 = ( dev2 - a * dev ) / ( drdi(ir)**2 )
 dev  = dev / drdi(ir)

!Local potential is Vloc(r)=v3*exp(v1*r^2+v2*r^3)
!inside Rgauss and equals the
!all-electron atomic potential outside Rgauss
!We impose the continuity up to second derivative

 if( ndevfit .eq. 2 ) then
   vlc = vlocal(nrgauss)
   r   = rofi(nrgauss)

   var1 = dev  / vlc - 1.0d0 / r
   var2 = dev2 / vlc - 2.0d0 * var1 / r  - ( var1**2 )

   dm11 = 2.0d0 * r
   dm12 = 3.0d0 * r * r
   dm21 = 2.0d0
   dm22 = 6.0d0 * r

   v1 = ( dm22 * var1 - dm12 * var2 ) /( 6.0d0 * r * r )
   v2 = ( dm11 * var2 - dm21 * var1 ) /( 6.0d0 * r * r )
   v3 = vlc / ( r * exp( ( v1 + v2*r ) * r * r ) )


!  elseif(ndevfit.eq.3) then
 else

!  We can also construct a local potential
!  Vloc(r)=v4*exp(v1*r^2+v2*r^3+v3*r^4),
!  this new coefficient allows us to impose the continuity
!  of the potential up  to the third derivative.

   vlc  = vlocal( nrgauss )
   r    = rofi( nrgauss )

   var1 = dev  / vlc - 1.d0 / r
   var2 = dev2 / vlc - 2.0d0 * var1 / r - ( var1**2 )
   var3 = dev3 / vlc - 3.0d0 * var1 * var2 - ( var1**3 ) &
&   - 3.0d0 *( var1**2 + var2 ) / r

   dm11 = 2.0d0 * r
   dm12 = 3.0d0 * r * r
   dm13 = 4.0d0 * r * r * r
   dm21 = 2.0d0
   dm22 = 6.0d0 * r
   dm23 = 12.0d0 * r * r
   dm31 = 0.0d0
   dm32 = 6.0d0
   dm33 = 24.0d0 * r

   v1 = ( ( var1 * dm22 * dm33 + var2 * dm13 * dm32 + var3 * dm12 * dm23 ) &
&   -(var3*dm22*dm13+var1*dm32*dm23+var2*dm12*dm33))/(48.0_dp*r*r*r)
   v2 = ( ( var2 * dm11 * dm33 + var3 * dm21 * dm13 + var1 * dm23 * dm31 ) &
&   -(var2*dm31*dm13+var3*dm23*dm11+var1*dm21*dm33))/(48.0_dp*r*r*r)
   v3 = ( ( var3 * dm11 * dm22 + var2 * dm12 * dm31 + var1 * dm32 * dm21 ) &
&   -(var1*dm22*dm31+var3*dm21*dm12+var2*dm11*dm32))/(48.0_dp*r*r*r)
   v4 = vlc / ( r * exp( ( v1 + v2 * r + v3 * r * r ) * r * r ) )

 end if

 do ir = 1, nrval
   r = rofi(ir)
   if( ir .le. nrgauss ) then
!    **   If second derivative fit***
     if( ndevfit .eq. 2 ) then
       vlocal(ir) = v3 * exp( ( v1 + v2*r ) * r * r )
!      **   If third derivative fit****
     else if(ndevfit.eq.3) then
       vlocal(ir) = v4 * exp ( ( v1 + v2 * r + v3 * r * r ) * r * r )
!      ****
     end if
   else
     vlocal(ir) = vps(ir)
   end if
 end do

!Once we have the local potential we define the 'local-pseudopotential
!charge' which help us to calculate the electrostatic interation
!between the ions
!
!Poisson's eq.:
!
!1/r* d2(rV)/dr2 = -8*pi*rho
!
 a2b4 = 0.25d0 * a * a
 qtot = 0.d0
 do ir = 1, nrval-1
   g2 = vlocal(ir) * rofi(ir)
!  
!  To determine the chlocal cutoff, use the reduced_vlocal cutoff
!  
   if( abs ( g2 + 2.0d0 * zval ) .lt. eps ) exit   !exit loop

   if( ir .gt. nrgauss ) then
     if( ( ir .gt. 2 ) .and. ( ir .lt. (nrval-1) ) ) then
       g0 = vlocal(ir-2) * rofi(ir-2) / s(ir-2)
       g1 = vlocal(ir-1) * rofi(ir-1) / s(ir-1)
       g2 = vlocal(ir)   * rofi(ir)   / s(ir)
       g3 = vlocal(ir+1) * rofi(ir+1) / s(ir+1)
       g4 = vlocal(ir+2) * rofi(ir+2) / s(ir+2)
       d2g = ( 16.d0 * ( g1 + g3 ) - ( g0 + g4 ) -30.d0 * g2 ) / 12.d0
     else
       g1 = vlocal(ir-1) * rofi(ir-1) / s(ir-1)
       g2 = vlocal(ir)   * rofi(ir)   / s(ir)
       g3 = vlocal(ir+1) * rofi(ir+1) / s(ir+1)
       d2g = g1 + g3 - 2.0d0 * g2
     end if
     d2u         = d2g - a2b4 * g2
     r           = rofi(ir)
     cons        = 8.0d0 * pi * r * drdi(ir) * s(ir)
     chlocal(ir) = (-d2u) / cons
     qtot        = qtot + 0.5d0 * d2u * r / s(ir)
   else
!    If second derivative fit
     if( ndevfit .eq. 2 )  then
       r  = rofi(ir)
       g0 = v3 * exp( ( v1 + v2 * r ) * r **2 )
       g1 = ( 2.d0 * v1 + 3.0d0 * v2 * r )
       g2 =   2.d0 * v1 + 6.0d0 * v2 * r
       g3 = ( g2 + g1 * g1 * r * r + 2.0d0 * g1 ) * g0
       cons        = 8.0d0 * pi
       chlocal(ir) = (-g3) / cons
       qtot        = qtot  + 0.5d0 * g3 * r * r * drdi(ir)
!      **** If third derivative fit
     else if ( ndevfit .eq. 3 )  then
       r  = rofi(ir)
       g0 = v4 * exp( ( v1 + v2 * r + v3 * r * r ) * r * r )
       g1 = ( 2.0d0 * v1 + 3.0d0 * v2 * r + 4.0d0  * v3 * r * r )
       g2 = ( 2.0d0 * v1 + 6.0d0 * v2 * r + 12.0d0 * v3 * r * r )
       g3 = ( g2 + g1 * g1 * r * r + 2.0d0 * g1 ) * g0

       cons        = 8.0d0 * pi
       chlocal(ir) = -g3 / cons
       qtot        = qtot  + 0.5d0 * g3 * r * r * drdi(ir)
     end if
   end if
 end do
!
!This sets the cutoff point for chlocal in a rather
!arbitrary way, as that in which Vlocal "equals" 2Z/r
!
 nchloc = ir

 do ir = 1, nchloc-1
   chlocal(ir) = zval * chlocal(ir) / qtot
 end do
 do ir = nchloc, nrval
   chlocal(ir) = 0.0_dp
 end do

end subroutine vlocal2
!!***

!!****f* ABINIT/vlocal1
!! NAME
!! vlocal1
!!
!! FUNCTION
!!  This routine generates a smooth local pseudopotential.
!!
!! NOTES
!!  Written by D. Sanchez-Portal, Aug. 1998
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      smoothvlocal
!!
!! CHILDREN
!!
!! SOURCE
subroutine vlocal1( zval, nrval, a, rofi, drdi, s, rgauss, vlocal,     &
&                    nchloc, chlocal)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vlocal1'
 use interfaces_65_psp, except_this_one => vlocal1
!End of the abilint section

 implicit none

 real(dp), intent(in)    :: zval, a
 integer,  intent(in)    :: nrval
 real(dp), intent(in)    :: rofi(:), drdi(:), s(:)
 real(dp), intent(out)   :: vlocal(:)
 real(dp), intent(out)   :: chlocal(:)
 integer,  intent(out)   :: nchloc
 real(dp), intent(inout) :: rgauss      !!???

! *Internal variables*

 real(dp) :: van, factor, alp, cutoff1, cutoff2,       &
&            qtot, eps, chc, r, Rchloc, rhor1, rhor
 integer  :: ir
 character loctype*3

 parameter(eps = 1.0d-4)

! *************************************************************************

!**   Usual local potential
!(generated with an optimum Vandebilt function)**
 loctype = 'new'

!***  The very first local potential used by SIESTA was
!the electrostatic potential generated by a gaussian
!distribution ===> loctype='old'
!loctype='old'
!***

!Local-potential size parameter 'rgauss'
!We choose as a smooth pseudopotential the one generated
!by a 'Vanderbilt-function' charge distribution. We have to select
!the size of this distribution somehow.
!'Vanderbilt-functions' are of the form :
!p(r)=N*exp(-(sinh(van*r)/sinh(van))**2)
!when van---> 0 we will obtain a 'gaussian'
!when van---> Inf. we will obtain a step function
!Some test has revealed that the best election to achieve
!a good convergence in real and reciprocal space is b in the
!range 0.5-1.0 .
!*

!So, the 'gaussian' charge distribution
!must go to zero at a distance 'rgauss'.

 if( loctype .eq. 'new' ) then
!  We take a 'Vanderbilt-function' as local potential
!  van=1.0_dp all the parameter have optimized for this value
   van     = 1.0d0
   cutoff1 = 3.63d0
   cutoff2 = 5.48d0
!  **   99% of charge inside Rgauss**
!  factor=1.627_dp
!  **   99.9% of charge inside Rgauss
   factor = 1.815d0
!  * Scaling factor for local-pseudopot. charge**
   alp = factor / rgauss
!  write(std_out,'(/,a,f10.3,a)')                          &
!  &    'VLOCAL1: 99.0% of the norm of Vloc inside ',   &
!  &    (alp*cutoff1)**2,' Ry'
!  write(std_out,'(a,f10.3,a)')                            &
!  &    'VLOCAL1: 99.9% of the norm of Vloc inside ',   &
!  &    (alp*cutoff2)**2,' Ry'
 else
!  This is just a gaussian !!!!!!!!!!!!!!!!!
   van    = 0.00001d0
   rgauss = 0.80d0
   factor = 2.0d0
!  * Scaling factor for local-pseudopot. charge**
   alp    = factor / rgauss
 end if

 qtot = 0.0d0
 rhor1 = vander( van, alp * rofi(1) )     ! This is 1...
 do ir = 1, nrval
   r           = rofi(ir)
   rhor        = vander( van, alp * r)
   chlocal(ir) = (-4.0d0) * pi * rhor * r * r
   qtot        = qtot + rhor * drdi(ir) * r * r
 end do

 qtot   = 4.0d0 * pi * qtot
 nchloc = 0
 do ir = nrval, 1, -1
   chc         = zval * chlocal(ir) / qtot
   chlocal(ir) = chc
   if( ( abs(chc) .gt. eps ) .and. ( nchloc .eq. 0 ) ) then
     nchloc = ir + 1
   end if
 end do
 Rchloc = rofi(nchloc)
!
!Note that the above cutoff is for 4*pi*r*r*rho_local(r)...
!
 call vhrtre( chlocal, vlocal, rofi, drdi, s, nrval, a )

 do ir = 2, nrval
   r = rofi(ir)
   chlocal(ir) = chlocal(ir) / ( 4.0d0 * pi * r * r )
!  
!  Poor man's cutoff!! Largely irrelevant?
!  
   if ( r .gt. 1.1d0 * Rchloc ) then
     vlocal(ir) = (-2.0d0) * zval / rofi(ir)
   end if
 end do
 chlocal(1) = -rhor1 * zval / qtot

 end subroutine vlocal1
!!***

!!****f* ABINIT/vhrtre
!! NAME
!! vhrtre
!!
!! FUNCTION
!!   Finds the Hartree potential created by a radial electron density,
!!   using Numerov's method to integrate the radial Poisson equation:
!!      d2(r*V)/dr2 = -4*pi*rho*r = -(4*pi*r2*rho)/r
!!
!! NOTES
!!   Imported from SIESTA package
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      smoothvlocal
!!
!! CHILDREN
!!
!! SOURCE
subroutine vhrtre(R2RHO,V,R,DRDI,SRDRDI,NR,A)

 use m_profiling
!***********************************************************************
!Finds the Hartree potential created by a radial electron density,
!using Numerov's method to integrate the radial Poisson equation:
!d2(r*V)/dr2 = -4*pi*rho*r = -(4*pi*r2*rho)/r
!Input:
!real*8  R2RHO(NR) : 4*pi*r**2*rho, with rho the electron density
!real*8  R(NR)     : Logarithmic radial mesh R(i)=B*(exp(A*(i-1)-1)
!real*8  DRDI(NR)  : dr/di at the mesh points
!real*8  SRDRDI(NR): sqrt(dr/di) at the mesh points
!integer NR        : Number of radial mesh points, including r(1)=0
!real*8  A         : The parameter A in r(i)=B*(exp(A*(i-1)-1)
!Output:
!real*8  V(NR)     : Electrostatic potential created by rho, in Ryd
!The constants of integration are fixed so that
!V=finite at the origin and V(NR)=Q/R(NR),
!where Q is the integral of rho up to R(NR)
!Algorithm: see routine NUMOUT
!***********************************************************************

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vhrtre'
!End of the abilint section

 implicit none

 INTEGER          :: NR
 REAL(dp) :: R2RHO(NR),V(NR),R(NR),DRDI(NR),SRDRDI(NR),A

 INTEGER          :: IR
 REAL(dp) :: A2BY4, BETA, DV, DY, DZ, Q, QBYY, QPARTC, QT,  &
& T, V0, Y, YBYQ

!Find some constants
 A2BY4 = A * A / 4.D0
 YBYQ  = 1.D0 - A * A / 48.D0
 QBYY  = 1.D0 / YBYQ

!Use Simpson's rule to find the total charge QT, and the
!potential at the origin V0:
!QT = Int(4*pi*r**2*rho*dr) = Int((4*pi*r**2*rho)*(dr/di)*di)
!V0 = Int(4*pi*r*rho*dr) =  Int((4*pi*r**2*rho)/r*(dr/di)*di)
 V0 = 0.D0
 QT = 0.D0
 do IR = 2, NR-1, 2
   DZ = DRDI(IR) * R2RHO(IR)
   QT = QT + DZ
   V0 = V0 + DZ / R(IR)
 end do
 V0 = V0 + V0
 QT = QT + QT
 do IR = 3, NR-2, 2
   DZ = DRDI(IR) * R2RHO(IR)
   QT = QT + DZ
   V0 = V0 + DZ / R(IR)
 end do
 DZ =DRDI(NR) * R2RHO(NR)
 QT =( QT + QT + DZ ) / 3.D0
 V0 =( V0 + V0 + DZ / R(NR) ) / 3.D0

!Fix V(1) and V(2) to start Numerov integration. To find a
!particular solution of the inhomog. eqn, V(2) is fixed by
!setting rV(2)=0. Notice that V=finite => rV=0 at r=0
 V(1)=2.D0*V0    ! Factor 2 because we use Rydbergs
 T    = SRDRDI(2) / R(2)
 BETA = DRDI(2) * T * R2RHO(2)
 DY   = 0.D0
 Y    = 0.D0
 Q    = ( Y - BETA / 12.D0 ) * QBYY
 V(2) = 2.D0 * T * Q

!Integrate Poisson's equation outwards, using Numerov's method
 do IR = 3,NR
   DY    = DY + A2BY4 * Q - BETA
   Y     = Y + DY
   T     = SRDRDI(IR) / R(IR)
   BETA  = T * DRDI(IR) * R2RHO(IR)
   Q     = ( Y - BETA / 12.D0 ) * QBYY
   V(IR) = 2.D0 * T * Q
 end do 

!Add a solution (finite at r=0) of the homogeneous equation
!d2(r*V)/dr2=0 => rV=const*r => V=const, to ensure that
!V(NR)=Q/R(NR). Notice that V(1) is set independently
 QPARTC = R(NR) * V(NR) / 2.D0
 DZ = QT - QPARTC
 DV = 2.D0 * DZ / R(NR)
 do IR = 2, NR
   V(IR) = V(IR) + DV
 end do

 end subroutine vhrtre
!!***
