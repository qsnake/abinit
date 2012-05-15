!{\src2tex{textfont=tt}}
!!****f* ABINIT/old_setmesh
!!
!! NAME
!! old_setmesh
!!
!! FUNCTION
!! Calculate the size of the FFT grid for the GW calculation.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, YMN, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! mG0= number of shell that must be added to take into account umklapp processes
!! enforce_sym = flag to enforce a FFT compatible with all the symmetry operations, both rotations and 
!!  fractional translations (should be default)
!! gmet(3,3) = reciprocal space metric.
!! gvec(3,npwvec) = G-vectors array.
!! method = integer flag for FFT grid (see below)
!! npwvec= number of G vectors in the array gvec max(npwwfn,npwsigx)
!! npwsigx = size of the dielectric or self-energy matrix.
!! npwwfn = number of G-vectors in the wavefunctions.
!! nsym = number of symmetry operations 
!! symrel(3,3,nsym) = symmetry operations in real space
!! tnons(3,nsym) = fractional translations
!!
!! OUTPUT
!! nfftot= ngfft(1)*ngfft(2)*ngfft(3)=number of points on the  FFT grid.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! NOTE
!! Four methods have been coded for the calculation of the mesh (see parameter "method" below):
!!
!!  method=0 FFT mesh defined by the user
!!  method=1 roughly takes the FFT box that encloses the larger of the two spheres of radius
!!            aliasing_factor*rwfn and rsigx, where rwfn and rsigx are the radius of the spheres
!!            with npwwfn and npwsigx planewaves respectively. default aliasing_factor is 1
!!  method=2 calculates the optimal FFT grid that allows aliasing only outside the sphere of the
!!            npwsigx planewaves (finer than method=1 with aliasing_factor=1).
!!  method=3 calculates the FFT grid needed to expand the density.
!!            (even finer than method=2, roughly corresponds to method=1 with aliasing_factor=2).
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!      leave_new,size_goed_fft,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine old_setmesh(gmet,gvec,ngfft,npwvec,npwsigx,npwwfn,nfftot,method,mG0,nsym,symrel,tnons,enforce_sym)

 use m_profiling

 use defs_basis
 use m_fft_mesh

 use defs_fftdata,    only : size_goed_fft
 use m_numeric_tools, only : denominator, mincm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'old_setmesh'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enforce_sym,method,npwsigx,npwvec,npwwfn,nsym
 integer,intent(out) :: nfftot
!arrays
 integer,intent(in) :: gvec(3,npwvec),mG0(3),symrel(3,3,nsym)
 integer,intent(inout) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),tnons(3,nsym)

!Local variables ------------------------------
!scalars
 integer :: aliasing_factor,fftalg,ierr,ig,ig1,ig1max,ig2,ig2max,ig3,ig3max,ii
 integer :: index,is,m1,m2,m3,mm1,mm2,mm3,n1,n2,n3,nt
 real(dp) :: ecuteff,ecutsigx,ecutwfn,g1,g2,g3,gsq,gsqmax,reff,rsigx,rwfn
 logical :: fft_ok
 character(len=500) :: msg
!arrays
 integer :: fftnons(3),fftsym(3),mdum(3)

!************************************************************************

 if (ANY(mg0(:)<0)) then 
   write(msg,'(4a,3i4)')ch10,&
&   ' setshell : ERROR - ',ch10,&
&   ' called with wrong value of mG0 = ',mG0(:)
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
 end if 
!
!Calculate the limits of the sphere of npwwfn G-vectors in each direction.
 m1=MAXVAL(ABS(gvec(1,1:npwwfn)))
 m2=MAXVAL(ABS(gvec(2,1:npwwfn)))
 m3=MAXVAL(ABS(gvec(3,1:npwwfn)))
!
!Calculate the limits of the sphere of npsigx G-vectors in each direction.
!Ensure that G+G0 will fit into the FFT grid, where G is any of the npwsigx/npweps vectors 
!and G0 is (i,j,k) [-nG0shell<i,j,k<nG0shell]. This is required when npwsigx>npwwfn since 
!we have to take into account umklapp G0 vectors to evaluate the oscillator matrix elements 
!(see rho_tw_g) or to symmetrize these quantities (see also cigfft).
 mm1=MAXVAL(ABS(gvec(1,1:npwsigx)))
 mm2=MAXVAL(ABS(gvec(2,1:npwsigx)))
 mm3=MAXVAL(ABS(gvec(3,1:npwsigx)))
 mm1=mm1+mG0(1) ; mm2=mm2+mG0(2) ; mm3=mm3+mG0(3)

 write(msg,'(2(2a,i8,a,3i6),2a,3i3)')ch10,&
& ' old_setmesh: npwwfn           = ',npwwfn, '; Max (m1,m2,m3)    = ',m1,m2,m3,ch10,&
& '          npweps/npwsigx   = ',npwsigx,'; Max (mm1,mm2,mm3) = ',mm1,mm2,mm3,ch10,&
& '          mG0 added        = ',mG0(:)
 call wrtout(std_out,msg,'COLL')
!
!Different FFT grids according to method.
 select case (method)

   case (0) 
     n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
     write(msg,'(3a,3(a,i3))')ch10,&
&     ' old_setmesh: COMMENT - ',ch10,&
&     ' mesh size enforced by user = ',n1,'x',n2,'x',n3
     call wrtout(std_out,msg,'COLL') 
     call wrtout(ab_out,msg,'COLL')
     ngfft(1)=n1 
     ngfft(2)=n2 
     ngfft(3)=n3
     ngfft(4)=2*(ngfft(1)/2)+1
     ngfft(5)=2*(ngfft(2)/2)+1
     ngfft(6)=   ngfft(3)
     nfftot=n1*n2*n3
     RETURN

   case (1)
     aliasing_factor=1 
     write(msg,'(2a,i3)')ch10,' using method 1 with aliasing_factor = ',aliasing_factor
     call wrtout(std_out,msg,'COLL')
     m1=m1*aliasing_factor
     m2=m2*aliasing_factor
     m3=m3*aliasing_factor
     
   case (2,3)
!    Calculate the radius of the sphere of npwwfn G-vectors.
     ecutwfn=-one
     do ig=1,npwwfn
       g1=REAL(gvec(1,ig))
       g2=REAL(gvec(2,ig))
       g3=REAL(gvec(3,ig))
       gsq=       gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&       two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
       ecutwfn=MAX(ecutwfn,gsq)
     end do
     rwfn=SQRT(ecutwfn) ; ecutwfn=two*ecutwfn*pi**2
!    Calculate the radius of the sphere of npwsigx/npweps G-vectors.
     ecutsigx=-one
     do ig=1,npwsigx
       g1=REAL(gvec(1,ig))
       g2=REAL(gvec(2,ig))
       g3=REAL(gvec(3,ig))
       gsq=      gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&       two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
       ecutsigx=MAX(ecutsigx,gsq)
     end do
     rsigx=SQRT(ecutsigx) ; ecutsigx=two*ecutsigx*pi**2
     write(msg,'(a,f7.3,3a,f7.3,a)')&
&     ' calculated ecutwfn           = ',ecutwfn,' [Ha] ',ch10,&
&     ' calculated ecutsigx/ecuteps = ',ecutsigx,' [Ha]'
     call wrtout(std_out,msg,'COLL')
!    
!    In the calculation of the GW self-energy or of the RPA dielectric matrix,
!    we make products Rho_12(r)=u_1*(r) u_2(r) of wavefunctions whose Fourier
!    coefficients lie in the sphere of radius rwfn. Such products will have non
!    vanishing Fourier coefficients in the whole sphere of radius 2*rwfn since:
!    Rho(G) = \sum_T u_1*(T) u_2(T+G) 
!    However, we only need the Fourier coefficients of Rho_12 that lie in the sphere 
!    of radius rsigx. we can thus allow aliasing outside that sphere, so that the FFT box 
!    may only enclose a sphere of radius:

!    reff=min(rsigx+rwfn,2.0*rwfn) 
     reff=rsigx+rwfn
!    
!    Extreme case: this yields back the GS FFT grid if full wavefunctions are considered
     if (method==3) reff=two*rwfn
     ecuteff=two*(pi*reff)**2
     write(msg,'(2a,i2,a,f7.3,a)')ch10,&
&     ' using method = ',method,' with ecuteff = ',ecuteff,' [Ha]'
     call wrtout(std_out,msg,'COLL')
!    
!    Search the limits of that sphere in each direction...
!    FIXME this might be wrong in case of augmentation in the FFT
     gsqmax=reff**2
     ig1max=2*m1+1
     ig2max=2*m2+1
     ig3max=2*m3+1
!    this is the correct coding
!    $ ig1max=MAX(2*m1+1,2*mm1+1,mm1+m1+1)
!    $ ig1max=MAX(2*m1+1,2*mm1+1,mm1+m1+1)
!    $ ig1max=MAX(2*m1+1,2*mm1+1,mm1+m1+1)
     m1=-1 ; m2=-1 ; m3=-1
     do ig1=0,ig1max
       do ig2=0,ig2max
         do ig3=0,ig3max
           g1=REAL(ig1)
           g2=REAL(ig2)
           g3=REAL(ig3)
           gsq=     gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&           two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
           if (gsq>gsqmax+tol6) CYCLE ! tol6 to improve portability
           m1=MAX(m1,ig1)
           m2=MAX(m2,ig2)
           m3=MAX(m3,ig3)
         end do
       end do
     end do

     case default
     write(msg,'(4a)')ch10,&
&     ' old_setmesh : BUG- ',ch10,&
&     ' method > 3 or < 0 is not allowed in old_setmesh ' 
     call wrtout(std_out,msg,'COLL') 
     call leave_new('COLL')
 end select 
!write(std_out,*) m1,m2,m3
!write(std_out,*) mm1,mm2,mm3
!
!Warning on low npwwfn
 if (m1<mm1 .or. m2<mm2 .or. m3<mm3) then
   write(msg,'(8a)')ch10,&
&   ' old_setmesh: COMMENT -',ch10,&
&   '  Note that npwwfn is small with respect to npweps or with respect to npwsigx.',ch10,&
&   '  Such a small npwwfn is a waste:',ch10,&
&   '  You could raise npwwfn without loss in cpu time.'
   call wrtout(std_out,msg,'COLL')
 end if
!
!Keep the largest of the m/mm and and find the FFT grid which is compatible 
!with the library and, if required, with the symmetry operations.
 m1=MAX(m1,mm1)
 m2=MAX(m2,mm2)
 m3=MAX(m3,mm3)

 if (enforce_sym==0) then 
!  
!  Determine the best size for the FFT grid *without* considering the symm ops.
!  Ideally n=2*m+1 but this  could not be allowed by the FFT library 
   call size_goed_fft(m1,n1)
   call size_goed_fft(m2,n2)
   call size_goed_fft(m3,n3)
!  
!  Calculate the number of FFT grid points
   nfftot=n1*n2*n3
!  Check if the FFT is compatible, write ONLY a warning if it breaks the symmetry
   fftnons(1)=n1
   fftnons(2)=n2
   fftnons(3)=n3
   fft_ok=.TRUE.
   rd: do ii=1,3 
     do is=1,nsym 
       nt=denominator(tnons(ii,is), ierr)
       if (((fftnons(ii)/nt)*nt) /= fftnons(ii)) then 
         fft_ok=.FALSE.
         exit rd
       end if 
     end do 
   end do rd
   if (.not.fft_ok) then 
     write(msg,'(5a)')ch10,&
&     ' old_setmesh: WARNING -',ch10,&
&     ' FFT mesh is not compatible with non-symmorphic translations',ch10
     call wrtout(std_out,msg,'COLL') 
     call wrtout(ab_out, msg,'COLL')
   end if 
!  Check only rotations
   if (.not.(check_rot_fft(nsym,symrel,n1,n2,n3))) then
     write(msg,'(5a)')ch10,&
&     ' old_setmesh: WARNING -',ch10,&
&     ' FFT mesh is not compatible with rotations',ch10
     call wrtout(std_out, msg,'COLL') 
     call wrtout(ab_out,msg,'COLL')
   end if 
 else 
!  
!  Determine the best size for the FFT grid considering symm ops.
!  Ideally n=2*m+1 but this  could not be allowed by the FFT library 
   write(msg,'(2a)')ch10,&
&   ' finding a FFT mesh compatible with all the symmetries'
   call wrtout(std_out,msg,'COLL')
   fftnons(:)=1
!  Found a FFT mesh compatible with the non-symmorphic operations 
   do ii=1,3
     fftnons(ii)=1
     do is=1,nsym 
       nt=denominator(tnons(ii,is), ierr)
       if (((fftnons(ii)/nt)*nt)/=fftnons(ii)) fftnons(ii)=mincm(fftnons(ii),nt)
     end do 
   end do 
   write(msg,'(a,3i3)')' old_setmesh: divisor mesh',fftnons(:)
   call wrtout(std_out,msg,'COLL')
!  
!  Check if also rotations preserve the grid.
!  Initial guess from previous m values.
   call size_goed_fft(m1,fftsym(1))
   call size_goed_fft(m2,fftsym(2))
   call size_goed_fft(m3,fftsym(3))
   mdum(1)=m1 ; mdum(2)=m2 ; mdum(3)=m3

   index=0
   do ! If a FFT division gets too large the code stops in size_goed_fft
     if ( check_rot_fft(nsym,symrel,fftsym(1),fftsym(2),fftsym(3)) .and.&
&     (MOD(fftsym(1),fftnons(1))==0) .and.&
&     (MOD(fftsym(2),fftnons(2))==0) .and.&
&     (MOD(fftsym(3),fftnons(3))==0)&         
&     ) exit 
     ii=MOD(index,3)+1
     mdum(ii)=mdum(ii)+1
     call size_goed_fft(mdum(ii),fftsym(ii))
     index=index+1
   end do 
!  
!  Got a good FFT grid, Calculate the number of FFT grid points
   n1=fftsym(1) ; n2=fftsym(2) ; n3=fftsym(3)
   nfftot=n1*n2*n3

!  !DEBUG dont uncomment this part for the time being
   if ( .not.( check_rot_fft(nsym,symrel,n1,n2,n3)).or. &
&   ( MOD(fftsym(1),fftnons(1))/=0) .and.      &
&   ( MOD(fftsym(2),fftnons(2))/=0) .and.      &
&   ( MOD(fftsym(3),fftnons(3))/=0) &
&   ) then 
     write(msg,'(a)')' old_setmesh : BUG during the generation of a symmetric FFT'
     call wrtout(std_out,msg,'COLL') 
     call leave_new('COLL')
   end if  
!  !ENDDEBUG
 end if ! enforce_sym

 write(msg,'(3(a,i3),2a,i8,a)')&
& ' old_setmesh: FFT mesh size selected  = ',n1,'x',n2,'x',n3,ch10,&
& '          total number of points  = ',nfftot,ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 ngfft(1)=n1 
 ngfft(2)=n2 
 ngfft(3)=n3
 ngfft(4)=2*(ngfft(1)/2)+1
 ngfft(5)=2*(ngfft(2)/2)+1
 ngfft(6)=   ngfft(3)
!
!Check the value of fftalg ie ngfft(7), 
!presently only SG library is allowed, see size_goed_fft.F90
 fftalg=ngfft(7)/100 
 if (fftalg==2 .or. fftalg==3 .or. fftalg==4) then 
   write(msg,'(8a)')ch10,&
&   ' fftmesh : ERROR - ',ch10,&
&   ' Presently only S. Goedecker routines are allowed in GW calculation',ch10,&
&   ' Action : check the value of fftalg (ngfft(7)) in your input file',ch10,&
&   ' or modify old_setmesh.F90 to be sure that the FFT mesh is compatible with the FFT library '
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
 end if 

end subroutine old_setmesh
!!***
