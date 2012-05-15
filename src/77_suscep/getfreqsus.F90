!{\src2tex{textfont=tt}}
!!****f* ABINIT/getfreqsus
!! NAME
!! getfreqsus
!!
!! FUNCTION
!! Create an array with frequencies and weights for the imaginary frequency
!! integrations (fluctuation-dissipation theorem) in xcacfd routine.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (YP).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nfreqs= number of frequencies/weights for the mesh
!!   (developers may use nfreqs=1 for tests)
!!  optfreq=
!!   0: use preassigned mesh (see defs_suscep module)
!!    nfreqs= 2: pick-up 2 highest frequencies of H_2 mesh
!!    nfreqs= 8: pick-up 8 frequencies inside Be_2 mesh, depending on freq1
!!    nfreqs= 9: pick-up 9 frequencies inside H_2 mesh, depending on freq1
!!    nfreqs=11: pick-up 11 highest frequencies of Be_2 mesh
!!    nfreqs=16: use full He mesh
!!    nfreqs=18: use full H_2 mesh
!!    nfreqs=20: use full He mesh good up to 8 Ha
!!    nfreqs=24: use full Be_2 mesh
!!   1: create linear mesh and weights for quadrature by Taylor rule
!!    freq1=starting frequency
!!    freq2=frequency increment
!!   2: create mesh and weights using Gauss-Legendre quadrature
!!    A first Gauss-Legendre mesh is built for interval [0,freq1], then
!!    a second one is obtained by transforming the first for the
!!    [freq1,+\infty[ interval. freq2 may be use to compress or expand the
!!    mesh on the second interval (a value of 1.0 is adequate for
!!    most cases). For practical reasons, nfreqs must be even.
!!
!! OUTPUT
!!  freqs(:)  = array of mesh frequencies (hartree)
!!  weights(:)= array of integration weights
!!
!! NOTES
!! This routine has been completely rewritten by YP without any consideration
!! for backward compatibility.
!!
!! The optfreq=0 option should be avoided as much as possible, as
!! it will probably be removed in future versions of ABINIT.
!!
!! Linear mesh generation (optfreq=1) has not been tested.
!!
!! TODO
!! a - enable the reading of the mesh from a file
!! b - document the Gauss-Legendre mesh transform
!! c - create Gauss mesh automatically from maximal and minimal eigenvalues
!!     in the Kohn-Sham spectrum, there pay attention to zero or nearly zero
!!     HOMO-LUMO gap
!!
!! PARENTS
!!      suscep
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine getfreqsus(freqs,weights,nfreqs,optfreq,freq1,freq2)

 use m_profiling

 use defs_basis
  use defs_suscep

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getfreqsus'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfreqs,optfreq
 real(dp),intent(in) :: freq1,freq2
!arrays
 real(dp),intent(out) :: freqs(nfreqs),weights(nfreqs)

!Local variables ------------------------------
!scalars
 integer,parameter :: gl_maxiter=128
 integer :: i,j,mfreqs
 real(dp) :: xl,xm
 character(len=500) :: message
!arrays
 real(dp),allocatable :: p1(:),p2(:),p3(:),pp(:),w(:),x(:),z(:),z1(:)
 logical,allocatable :: unfinished(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' getfreqsus : enter '
!write(std_out,*)' getfreqsus : nfreqs, optfreq, freq1, freq2 =',nfreqs,optfreq,freq1,freq2
!ENDDEBUG

!Check arguments
 if ( nfreqs == 0 ) then
   write(message, '(5a)' ) ch10,&
&   ' getfreqsus: ERROR -',ch10,&
&   '  freqs and weights must be of nonzero size, however nfreqs = 0',ch10
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if ( nfreqs == 1 ) then
   freqs(1)   = freq1
   weights(1) = 1.0_dp
   write(message, '(4a,e12.5,2a)' ) ch10,&
&   ' getfreqsus: INFO -',ch10,&
&   '  Test mode with freqs=(',freq1,'), weights=( 0.10000E+01)',ch10
   call wrtout(std_out,message,'COLL')
   return
 end if

 select case(optfreq)

   case(0)
!    optfreq = 0: use preassigned meshes
!    
!    NOTE: May be removed from future versions.
!    

     select case(nfreqs)

       case(2)

         write(message, '(5a)' ) ch10,&
&         ' getfreqsus: INFO -',ch10,&
&         '  Using preassigned mesh for H2 with 2 frequencies',ch10
         call wrtout(std_out,message,'COLL')

         freqs(:)   = freqs_pa_h2(17:18)
         weights(:) = weights_pa_h2(17:18)

       case(8)

         write(message, '(5a)' ) ch10,&
&         ' getfreqsus: INFO -',ch10,&
&         '  Using preassigned mesh for Be/Be2 with 8 frequencies',ch10
         call wrtout(std_out,message,'COLL')

         if ( freq1 < freqs_pa_be2(8) ) then
           freqs(:)   = freqs_pa_be2(1:8)
           weights(:) = weights_pa_be2(1:8)
         else if ( freq1 < freqs_pa_be2(16) ) then
           freqs(:)   = freqs_pa_be2(9:16)
           weights(:) = weights_pa_be2(9:16)
         else
           freqs(:)   = freqs_pa_be2(17:24)
           weights(:) = weights_pa_be2(17:24)
         end if

       case(9)

         write(message, '(5a)' ) ch10,&
&         ' getfreqsus: INFO -',ch10,&
&         '  Using preassigned mesh for H2 with 9 frequencies',ch10
         call wrtout(std_out,message,'COLL')

         if ( freq1 < freqs_pa_h2(9) ) then
           freqs(:)   = freqs_pa_h2(1:9)
           weights(:) = weights_pa_h2(1:9)
         else
           freqs(:)   = freqs_pa_h2(10:18)
           weights(:) = weights_pa_h2(10:18)
         end if

       case(11)

         write(message, '(5a)' ) ch10,&
&         ' getfreqsus: INFO -',ch10,&
&         '  Using preassigned mesh for Be/Be2 with 11 frequencies',ch10
         call wrtout(std_out,message,'COLL')

         freqs(:)   = freqs_pa_be2(14:24)
         weights(:) = weights_pa_be2(14:24)

       case(16)

         write(message, '(5a)' ) ch10,&
&         ' getfreqsus: INFO -',ch10,&
&         '  Using preassigned mesh for He with 16 frequencies',ch10
         call wrtout(std_out,message,'COLL')

         freqs(:)   = freqs_pa_he(:)
         weights(:) = weights_pa_he(:)

       case(18)

         write(message, '(5a)' ) ch10,&
&         ' getfreqsus: INFO -',ch10,&
&         '  Using preassigned mesh for H2 with 18 frequencies',ch10
         call wrtout(std_out,message,'COLL')

         freqs(:)   = freqs_pa_h2(:)
         weights(:) = weights_pa_h2(:)

       case(20)

         write(message, '(5a)' ) ch10,&
&         ' getfreqsus: INFO -',ch10,&
&         '  Using preassigned mesh for He/8Ha with 20 frequencies',ch10
         call wrtout(std_out,message,'COLL')

         freqs(:)   = freqs_pa_he_8ha(:)
         weights(:) = weights_pa_he_8ha(:)

       case(24)

         write(message, '(5a)' ) ch10,&
&         ' getfreqsus: INFO -',ch10,&
&         '  Using preassigned mesh for Be/Be2 with 24 frequencies',ch10
         call wrtout(std_out,message,'COLL')

         freqs(:)   = freqs_pa_be2(:)
         weights(:) = weights_pa_be2(:)

         case default

         write(message, '(4a,i3,a)' ) ch10,&
&         ' getfreqsus: ERROR -',ch10,&
&         '  No preassigned mesh for nfreqs=',nfreqs,ch10
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')

     end select

   case(1)
!    optfreq = 1: build linear mesh from freq1 with freq2 as increment

     do i=1,nfreqs
       freqs(i) = freq1+(i-1)*freq2
     end do

!    Weights for Taylor integration
     weights(:)      = 1.0_dp
     weights(1)      = 0.5_dp
     weights(nfreqs) = 0.5_dp
     weights(:)      = weights(:)/freq2

   case(2)
!    optfreq = 2: build Gauss-Legendre mesh for [0,freq1] and transform it
!    to [freq1,+\infty[

!    Check arguments
     if ( mod(nfreqs,2) == 1 ) then
       write(message, '(4a,i3,a)' ) ch10,&
&       ' getfreqsus: ERROR -',ch10,&
&       '  freqs and weights must be of even size, however nfreqs =',&
&       nfreqs,ch10
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     mfreqs  = (nfreqs/2+1)/2
     xm = 0.5_dp*(freq1+zero)
     xl = 0.5_dp*(freq1-zero)

     ABI_ALLOCATE(x,(nfreqs/2))
     ABI_ALLOCATE(w,(nfreqs/2))
     ABI_ALLOCATE(unfinished,(mfreqs))
     ABI_ALLOCATE(p1,(mfreqs))
     ABI_ALLOCATE(p2,(mfreqs))
     ABI_ALLOCATE(p3,(mfreqs))
     ABI_ALLOCATE(pp,(mfreqs))
     ABI_ALLOCATE(z,(mfreqs))
     ABI_ALLOCATE(z1,(mfreqs))
     unfinished(:) = .true.
     x(:) = 0.0
     w(:) = 0.0
     do i=1,mfreqs
       z(i) = cos(pi*(i-0.25_dp)/(nfreqs/2+0.5_dp))
     end do

!    First build generic mesh for [-1,1] interval

     do i = 1,gl_maxiter
       where ( unfinished )
       p1 = 1.0
       p2 = 0.0
       end where
       do j=1,nfreqs/2
         where ( unfinished )
         p3=p2
         p2=p1
         p1=((2.0_dp*j - 1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
         end where
       end do

       where ( unfinished )
       pp = nfreqs*(z*p1-p2)/(z*z-1.0_dp)
       z1 = z
       z  = z1-p1/pp
       unfinished = ( abs(z-z1) > tol14 )
       end where

       if ( .not. any(unfinished) ) exit
     end do

     if ( i > gl_maxiter ) then
       write(message, '(4a,i2,2a)' ) ch10,&
&       ' getfreqsus: ERROR -',ch10,&
&       '  Too many iterations (',i,') while building Gauss-Legendre mesh',ch10
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     x(1:mfreqs) = -z(:)
     w(1:mfreqs) = 8.0_dp/((1.0_dp-z(:)**2)*pp**2)
     x(nfreqs/2:nfreqs/2-mfreqs+1:-1) = z
     w(nfreqs/2:nfreqs/2-mfreqs+1:-1) = w(1:mfreqs)

!    Now build the mesh for [0,freq1] ...
     freqs(1:nfreqs/2)   = xl*x(:)+xm
     weights(1:nfreqs/2) = xl*w(:)

!    ... and for [freq1,+\infty]
     freqs(nfreqs/2+1:nfreqs)   = 2.0_dp*freq1/((1-x(:))**freq2)
     weights(nfreqs/2+1:nfreqs) = &
&     freq2*(2.0_dp**freq2)*freq1*w(:)/((1-x(:))**(freq2+1.0_dp))

!    Deallocate temporary variables
     ABI_DEALLOCATE(p1)
     ABI_DEALLOCATE(p2)
     ABI_DEALLOCATE(p3)
     ABI_DEALLOCATE(pp)
     ABI_DEALLOCATE(z)
     ABI_DEALLOCATE(z1)
     ABI_DEALLOCATE(unfinished)
     ABI_DEALLOCATE(x)
     ABI_DEALLOCATE(w)

!    Any other option raises an error
     case default

     write(message, '(4a,i2,a)' ) ch10,&
&     ' getfreqsus: ERROR -',ch10,&
&     '  Unrecognized value for optfreq: ',optfreq,ch10
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')

 end select

!Print out the mesh
 write(std_out,'(1x,a,a)') 'getfreqsus: frequencies and weights',ch10
 write(std_out,'(3x,a5,2(1x,a12))') 'Index','Frequency','Weight'
 write(std_out,'(3x,a5,2(1x,a12))') '-----','--------------','--------------'

 do i=1,nfreqs
   write(std_out,'(3x,i5,2(1x,e12.5))') i,freqs(i),weights(i)
 end do

!Calculate integral of 1/(1+x*x) to check numerical accuracy of the mesh
 ABI_ALLOCATE(x,(nfreqs))
 x(:) = 1.0_dp/(1.0_dp+(freqs(:)/freq1)**2)
 xm = dot_product(x(:),weights(:))
 ABI_DEALLOCATE(x)

 write(message, '(4a,e12.5,a)' ) ch10,&
& ' getfreqsus: INFO -',ch10,&
& '  Residue for current mesh (integral of 1/(1+x*x)): ',xm-pi/2.0_dp,ch10
 call wrtout(std_out,message,'COLL')

!DEBUG
!write(std_out,*)' getfreqsus : exit '
!ENDDEBUG

end subroutine getfreqsus
!!***
