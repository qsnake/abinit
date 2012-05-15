!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fft_mesh
!! NAME
!!  m_fft_mesh
!!
!! FUNCTION
!!  This module defines a routine used to perform the setup of the FFT mesh 
!!  employed for the oscillator matrix elements used in GW. It also provides
!!  a set of tools to test the grid, rotate the mesh according to the symmetry 
!!  operations of the space group (SG) etc.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG, GMR, VO, LR, RWG, YMN, RS)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_fft_mesh

 use m_profiling

 use defs_basis
 use m_errors

 use defs_fftdata,     only : size_goed_fft
 use m_fftw3,          only : fftw3_gain_wisdom
 use m_numeric_tools,  only : denominator, mincm, iseven, pfactorize
 use m_crystal,        only : crystal_structure

 implicit none

 private 

 public :: setmesh             ! Perform the setup of the FFT mesh for the GW oscillator strengths.
 public :: check_rot_fft       ! Test whether the mesh is compatible with the rotational part of the space group.
 public :: fft_check_rotrans   ! Test whether the mesh is compatible with the symmetries of the space group.
 public :: rotate_FFT_mesh     ! Calculate the FFT index of the rotated mesh.
 public :: cigfft              ! Calculate the FFT index of G-G0. 
 public :: ig2gfft             ! Returns the component of a G in the FFT Box from its sequential index.
 public :: g2ifft              ! Returns the index of the G in the FFT box from its reduced coordinates.
 public :: get_gftt            ! Calculate the G"s in the FFT box from ngfft
 public :: print_ngfft         ! Print the content of ngfft(18) in explicative format.
 public :: fftalg_info         ! Returns strings with info on the FFT library specified by fftalg.
 public :: ceigr               ! e^{iG.r} on the FFT mesh (complex valued).
 public :: eigr                ! e^{iG.r} on the FFT mesh (complex valued).
 public :: ceikr               ! e^{ik.r} on the FFT mesh (complex valued).

#define FFTALGA_SIZE 4
 character(len=*),private,parameter :: fftalga2name(1:FFTALGA_SIZE)= &
& (/"Goedecker     ", &
&   "Vendor FFT    ", &
&   "FFTW3         ", &
&   "Goedecker2002 "/)

#define FFTALGB_SIZE 1
 character(len=*),private,parameter :: fftalgb2name(0:FFTALGB_SIZE)= &
& (/"C2C",&
&   "R2C"/)

#define FFTALGC_SIZE 2
 character(len=*),private,parameter :: fftalgc2name(0:FFTALGC_SIZE)= &
& (/"No pad         ",&
&   "zero-pad       ",&
&   "zero-pad+cache "/)

CONTAINS  !========================================================================================
!!***

!!****f* m_fft_mesh/setmesh
!!
!! NAME
!! setmesh
!!
!! FUNCTION
!! Calculate the size of the FFT grid for the GW calculation.
!!
!! INPUTS
!!  gmet(3,3)=Reciprocal space metric.
!!  gvec(3,npwvec)=G-vectors in reduced coordinates.
!!  npwvec=Number of G vectors in the array gvec max(npwwfn,npwsigx)
!!  npwsigx=Size of the dielectric or self-energy matrix.
!!  npwwfn=Number of G-vectors in the wavefunctions.
!!  method=Integer flag for FFT grid (see below)
!!  mG0=Number of shells that must be added to take into account umklapp processes.
!!  Cryst<Crystal_structure>=Data type gathering information on unit cell and symmetries
!!    %nsym=Number of symmetry operations in the SG.
!!    %symrel(3,3,nsym)=Symmetry operations in real space.
!!    %tnons(3,nsym)=Fractional translations.
!!  enforce_sym=Flag to enforce a FFT which fulfils all symmetry operations, both the 
!!   rotational part and fractional translations.
!!
!! OUTPUT
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see also ~abinit/doc/input_variables/vargs.htm#ngfft
!! nfftot= ngfft(1)*ngfft(2)*ngfft(3)=Total number of points in the FFT grid.
!!
!! NOTES
!! Four methods are implemented for the calculation of the mesh:
!!  method=0 --> FFT mesh defined by the user, useful for debugging.
!!  method=1     Roughly takes the FFT box which encloses the larger of the two spheres of radius
!!               aliasing_factor*rwfn and rsigx, where rwfn and rsigx are the radius of the spheres
!!               with npwwfn and npwsigx planewaves respectively. The default aliasing_factor is 1.
!!  method=2 --> Calculates the optimal FFT grid which allows aliasing only outside the sphere of the
!!               npwsigx planewaves (finer than method=1 with aliasing_factor=1).
!!  method=3 --> Calculates the FFT grid needed to expand the density.
!!               (even finer than method=2, roughly corresponds to method=1 with aliasing_factor=2).
!!
!!  See defs_fftdata for a list of allowable sizes of FFT. 
!!
!! PARENTS
!!      m_shirley,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine setmesh(gmet,gvec,ngfft,npwvec,npwsigx,npwwfn,nfftot,method,mG0,Cryst,enforce_sym,silent)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setmesh'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enforce_sym,method,npwsigx,npwvec,npwwfn
 integer,intent(out) :: nfftot
 logical,optional,intent(in) :: silent
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npwvec),mG0(3)
 integer,intent(inout) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3)

!Local variables ------------------------------
!scalars
 integer :: aliasing_factor,fftalg,fftalga,fftalgc,ig,ig1,ig1max,ig2,ig2max,ig3,ig3max,ii,idx,ierr
 integer :: is,m1,m2,m3,mm1,mm2,mm3,n1,n2,n3,nsym,nt
 real(dp) :: ecuteff,ecutsigx,ecutwfn,g1,g2,g3,gsq,gsqmax,reff,rsigx,rwfn
 logical :: fft_ok
 character(len=500) :: msg
!arrays
 integer :: fftnons(3),fftsym(3),mdum(3)
 integer,allocatable :: pfactors(:),powers(:)
 integer,pointer :: symrel(:,:,:)
 real(dp),pointer :: tnons(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 if (ANY(mg0<0)) then
   write(msg,'(a,3i4)')' called with wrong value of mG0 = ',mG0
   MSG_BUG(msg)
 end if

 nsym   =  Cryst%nsym
 symrel => Cryst%symrel
 tnons  => Cryst%tnons
 !
 ! Calculate the limits of the sphere of npwwfn G-vectors in each direction.
 m1=MAXVAL(ABS(gvec(1,1:npwwfn)))
 m2=MAXVAL(ABS(gvec(2,1:npwwfn)))
 m3=MAXVAL(ABS(gvec(3,1:npwwfn)))
 !
 ! Calculate the limits of the sphere of npsigx G-vectors in each direction.
 ! Ensure that G+G0 will fit into the FFT grid, where G is any of the npwsigx/npweps vectors 
 ! and G0 is (i,j,k) [-nG0shell<i,j,k<nG0shell]. This is required when npwsigx>npwwfn since 
 ! we have to take into account umklapp G0 vectors to evaluate the oscillator matrix elements 
 ! (see rho_tw_g) or to symmetrize these quantities (see also cigfft).
 mm1=MAXVAL(ABS(gvec(1,1:npwsigx)))
 mm2=MAXVAL(ABS(gvec(2,1:npwsigx)))
 mm3=MAXVAL(ABS(gvec(3,1:npwsigx)))

 mm1=mm1+mG0(1) 
 mm2=mm2+mG0(2) 
 mm3=mm3+mG0(3)

 ! To avoid possible wrap-around errors in cigfft, it is safe to start 
 ! with odd divisions so that the FFT box is centered on Gamma 
 ! This holds only if npwsigx > npwwfn.
 !if (iseven(mm1)) mm1=mm1+1
 !if (iseven(mm2)) mm2=mm2+1
 !if (iseven(mm3)) mm3=mm3+1

 write(msg,'(2(2a,i8,a,3i6),2a,3i3)')ch10,&
&  ' setmesh: npwwfn        = ',npwwfn, '; Max (m1,m2,m3)   = ',m1,m2,m3,ch10,&
&  '          npweps/npwsigx= ',npwsigx,'; Max (mm1,mm2,mm3)= ',mm1,mm2,mm3,ch10,&
&  '          mG0 added     = ',mG0(:)
 call wrtout(std_out,msg,'COLL')
 !
 ! === Different FFT grids according to method ==
 select case (method)

 case (0) 
   ! * FFT mesh defined by user, useful for testing.
   !fnam='__fft.in__'
   !inquire(file=fnam,exist=fftfile)
   !if (fftfile) then 
   ! unt=io_unused()
   ! open(file=fnam,unit=unt,form='formatted')
   ! read(unt,*)n1,n2,n3 
   ! close(unt)
   !else 
   ! write(msg,'(5a)')' setmesh : ERROR : ',ch10,' FFT file ',TRIM(fnam),' not found'
   ! call wrtout(std_out,msg,'COLL'); call leave_new('COLL')
   !end if 
   n1=ngfft(1) 
   n2=ngfft(2) 
   n3=ngfft(3)
   write(msg,'(3(a,i3))')' Mesh size enforced by user = ',n1,'x',n2,'x',n3
   MSG_COMMENT(msg)

   ngfft(1)=n1 
   ngfft(2)=n2 
   ngfft(3)=n3
   ngfft(4)=2*(ngfft(1)/2)+1
   ngfft(5)=2*(ngfft(2)/2)+1
   ngfft(6)=   ngfft(3)
   !ngfft(4:6)=ngfft(1:3)
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
  
   ecutwfn=-one  ! Calculate the radius of the sphere of npwwfn G-vectors.
   do ig=1,npwwfn
     g1=REAL(gvec(1,ig))
     g2=REAL(gvec(2,ig))
     g3=REAL(gvec(3,ig))
     gsq=       gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&          two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
     ecutwfn=MAX(ecutwfn,gsq)
   end do
   rwfn=SQRT(ecutwfn); ecutwfn=two*ecutwfn*pi**2

   ! * Calculate the radius of the sphere of (npwsigx|npweps) G-vectors.
   ecutsigx=-one
   do ig=1,npwsigx
     g1=REAL(gvec(1,ig))
     g2=REAL(gvec(2,ig))
     g3=REAL(gvec(3,ig))
     gsq=      gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&         two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
     ecutsigx=MAX(ecutsigx,gsq)
   end do
   rsigx=SQRT(ecutsigx); ecutsigx=two*ecutsigx*pi**2

   write(msg,'(a,f7.3,3a,f7.3,a)')&
&    ' calculated ecutwfn          = ',ecutwfn, ' [Ha] ',ch10,&
&    ' calculated ecutsigx/ecuteps = ',ecutsigx,' [Ha]'
   call wrtout(std_out,msg,'COLL')
   ! 
   ! In the calculation of the GW self-energy or of the RPA dielectric matrix,
   ! we have products $ \rho_{12}(r)=u_1*(r) u_2(r) $ of wavefunctions whose Fourier
   ! coefficients lie in the sphere of radius rwfn. Such products will have non
   ! vanishing Fourier coefficients in the whole sphere of radius 2*rwfn since:
   !  $ rho_{12}(G) = \sum_T u_1*(T) u_2(T+G) $. 
   ! However, we only need the Fourier coefficients of $rho_{12}$ that lie in the sphere 
   ! of radius rsigx. We can thus allow aliasing outside that sphere, so that the FFT box 
   ! will only enclose a sphere of radius reff given by:

   reff=rsigx+rwfn
   if (method==3) reff=two*rwfn ! Yields back the GS FFT grid if full wavefunctions are considered.
   ecuteff=two*(pi*reff)**2
   gsqmax=reff**2

   write(msg,'(a,i2,a,f7.3,a)')' using method = ',method,' with ecuteff = ',ecuteff,' [Ha]'
   call wrtout(std_out,msg,'COLL')
   ! 
   ! === Search the limits of the reff sphere in each direction ===
   !ig1max=2*m1+1
   !ig2max=2*m2+1
   !ig3max=2*m3+1
   if (method==2) then
     ig1max=mm1+m1+1
     ig2max=mm2+m2+1
     ig3max=mm3+m3+1
   else if (method==3) then
     ig1max=MAX(2*m1+1,2*mm1+1,mm1+m1+1)
     ig2max=MAX(2*m2+1,2*mm2+1,mm2+m2+1)
     ig3max=MAX(2*m3+1,2*mm3+1,mm3+m3+1)
   end if

   m1=-1; m2=-1; m3=-1
   do ig1=0,ig1max
     do ig2=0,ig2max
       do ig3=0,ig3max
         g1=REAL(ig1)
         g2=REAL(ig2)
         g3=REAL(ig3)
         gsq=     gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&            two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
         if (gsq>gsqmax+tol6) CYCLE ! tol6 to improve portability
         m1=MAX(m1,ig1)
         m2=MAX(m2,ig2)
         m3=MAX(m3,ig3)
       end do
     end do
   end do

 case default
   write(msg,'(a,i3)')' Method > 3 or < 0 not allowed in setmesh while method= ',method
   MSG_BUG(msg)
 end select 
 !
 ! * Warning if low npwwfn.
 if (m1<mm1 .or. m2<mm2 .or. m3<mm3) then
   write(msg,'(5a)')&
&    ' Note that npwwfn is small with respect to npweps or with respect to npwsigx. ',ch10,&
&    ' Such a small npwwfn is a waste: ',ch10,&
&    ' You could raise npwwfn without loss in cpu time. '
   MSG_COMMENT(msg)
 end if
 !
 ! Keep the largest of the m/mm and and find the FFT grid which is compatible 
 ! with the library and, if required, with the symmetry operations.
 m1=MAX(m1,mm1)
 m2=MAX(m2,mm2)
 m3=MAX(m3,mm3)

 if (enforce_sym==0) then 
   ! === Determine the best size for the FFT grid *without* considering the symm ops ===
   ! * Ideally n=2*m+1 but this  could not be allowed by the FFT library.
   call size_goed_fft(m1,n1)
   call size_goed_fft(m2,n2)
   call size_goed_fft(m3,n3)
   nfftot=n1*n2*n3

   ! * Check if the FFT is compatible, write ONLY a warning if it breaks the symmetry
   fftnons(1)=n1
   fftnons(2)=n2
   fftnons(3)=n3
   fft_ok=.TRUE.
   rd: do ii=1,3 
     do is=1,nsym 
       nt=denominator(tnons(ii,is), ierr)
       if (((fftnons(ii)/nt)*nt) /= fftnons(ii)) then 
         fft_ok=.FALSE.; EXIT rd
       end if 
     end do 
   end do rd
   !
   ! * Warn if not compatibile with tnons or rotational part.
   if (.not.fft_ok) then 
    msg = ' FFT mesh is not compatible with non-symmorphic translations '
    MSG_WARNING(msg)
   end if 
   if (.not.(check_rot_fft(nsym,symrel,n1,n2,n3))) then
     msg = ' FFT mesh is not compatible with rotations '
     MSG_WARNING(msg)
   end if 

 else 
   ! === Determine the best size for the FFT grid considering symm ops ===
   ! * Ideally n=2*m+1 but this could not be allowed by the FFT library (at present only Goedecker)
   write(msg,'(2a)')ch10,' finding a FFT mesh compatible with all the symmetries'
   call wrtout(std_out,msg,'COLL')

   ! 1) Find a FFT mesh compatible with the non-symmorphic operations 
   fftnons(:)=1
   do ii=1,3
     fftnons(ii)=1
     do is=1,nsym 
       nt=denominator(tnons(ii,is), ierr)
       if (((fftnons(ii)/nt)*nt)/=fftnons(ii)) fftnons(ii)=mincm(fftnons(ii),nt)
     end do 
   end do 
   write(msg,'(a,3i3)')' setmesh: divisor mesh',fftnons(:)
   call wrtout(std_out,msg,'COLL')
   ! 
   ! 2) Check if also rotations preserve the grid.
   ! * Use previous m values as Initial guess.
   call size_goed_fft(m1,fftsym(1))
   call size_goed_fft(m2,fftsym(2))
   call size_goed_fft(m3,fftsym(3))
   mdum(1)=m1 
   mdum(2)=m2 
   mdum(3)=m3

   idx=0
   do ! If a FFT division gets too large the code stops in size_goed_fft.
     if ( check_rot_fft(nsym,symrel,fftsym(1),fftsym(2),fftsym(3)) .and. &
&         (MOD(fftsym(1),fftnons(1))==0) .and.                           &
&         (MOD(fftsym(2),fftnons(2))==0) .and.                           &
&         (MOD(fftsym(3),fftnons(3))==0)                                 &         
&       ) EXIT
     ii=MOD(idx,3)+1
     mdum(ii)=mdum(ii)+1
     call size_goed_fft(mdum(ii),fftsym(ii))
     idx=idx+1
   end do 
   ! 
   ! Got a good FFT grid, Calculate the number of FFT grid points
   n1=fftsym(1) 
   n2=fftsym(2) 
   n3=fftsym(3); nfftot=n1*n2*n3

   if (.not.( check_rot_fft(nsym,symrel,n1,n2,n3)) &
&       .or.( MOD(fftsym(1),fftnons(1))/=0) .and.  &
&           ( MOD(fftsym(2),fftnons(2))/=0) .and.  &
&           ( MOD(fftsym(3),fftnons(3))/=0)        &
&     ) then 
     MSG_BUG('Not able to generate a symmetric FFT')
   end if  
 end if ! enforce_sym

 write(msg,'(3(a,i3),2a,i8,a)')&
&  ' setmesh: FFT mesh size selected  = ',n1,'x',n2,'x',n3,ch10,&
&  '          total number of points  = ',nfftot,ch10
 call wrtout(std_out,msg,'COLL')
 if (.not.PRESENT(silent)) call wrtout(ab_out,msg,'COLL')

 ngfft(1)=n1 
 ngfft(2)=n2 
 ngfft(3)=n3
 ngfft(4)=2*(ngfft(1)/2)+1
 ngfft(5)=2*(ngfft(2)/2)+1
 ngfft(6)=   ngfft(3)
 !ngfft(4:6) = ngfft(1:3)
 !
 ! === Check the value of fftalg i.e ngfft(7) ===
 ! * Presently only Goedecker"s library or FFTW3 are allowed, see size_goed_fft.F90
 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=MOD(fftalg,10)

 if ( ALL(fftalga /= (/1,3/)) ) then 
   write(msg,'(6a)')ch10,&
&    " Only Goedecker's routines with fftalg=1xx or FFTW3 routines are allowed in GW calculations. ",ch10,&
&    " Action : check the value of fftalg in your input file, ",ch10,&
&    " or modify setmesh.F90 to make sure the FFT mesh is compatible with the FFT library. "
   MSG_ERROR(msg)
 end if 

 ! TODO Had to change setmesh to avoid bad values for FFTW3 
 if (fftalga==3) then ! check whether mesh is optimal for FFTW3
   ABI_ALLOCATE(pfactors,(5))
   ABI_ALLOCATE(powers,(6))
   pfactors = (/2, 3, 5, 7, 11/)
   do ii=1,3
     call pfactorize(ngfft(ii),5,pfactors,powers)
     if (powers(6)/=1 .or. powers(4)/=0 .or. powers(5)/=0) then 
       write(msg,'(a,i0,a)')&
&        "ngfft(ii) ",ngfft(ii)," contains powers of 7-11 or greater; FFTW3 is not optimal "
       MSG_WARNING(msg)
     end if
   end do

   ABI_DEALLOCATE(pfactors)
   ABI_DEALLOCATE(powers)
 end if

 DBG_EXIT("COLL")

end subroutine setmesh
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/check_rot_fft
!! NAME
!!   check_rot_fft
!!
!! FUNCTION
!!  Return .TRUE. if the given grid in real space is compatible 
!!  with the rotational part of the space group symmetries.
!!
!! INPUTS
!!  nsym=Number of symmetry operations
!!  symrel(3,3,nsym)=Symmetry operations in real space.
!!  nr1,nr2,nr3=FFT divisions.
!!
!! SOURCE

logical function check_rot_fft(nsym,symrel,nr1,nr2,nr3)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_rot_fft'
!End of the abilint section

 implicit none

!Arguments
!Scalar
 integer,intent(in) :: nr1,nr2,nr3,nsym 
!Arrays
 integer,intent(in) :: symrel(3,3,nsym)

!local variables
 integer :: is 

!************************************************************************

 ! The grid is compatible with the symmetries (only rotational part) if 
 ! for each symmetry, each n_i and n_j ==> $n_i*R_{ij}/n_j$ is an integer
 check_rot_fft=.TRUE.
 do is=1,nsym
   if ( MOD(symrel(2,1,is)*nr2, nr1) /=0 .or. &
&       MOD(symrel(3,1,is)*nr3, nr1) /=0 .or. &
&       MOD(symrel(1,2,is)*nr1, nr2) /=0 .or. &
&       MOD(symrel(3,2,is)*nr3, nr2) /=0 .or. &
&       MOD(symrel(1,3,is)*nr1, nr3) /=0 .or. &
&       MOD(symrel(2,3,is)*nr2, nr3) /=0      &
&     ) then
     check_rot_fft=.FALSE.; EXIT
   end if 
 end do

end function check_rot_fft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/fft_check_rotrans
!! NAME
!! fft_check_rotrans
!!
!! FUNCTION
!!  Checks if the real space FFT mesh is compatible both with the rotational 
!!  and the translational part of space group of the crystal.
!!
!! INPUTS
!!  nsym=Number of symmetries.
!!  symrel(3,3,nsym)=Symmetries in real space in reduced coordinates.
!!  tnons(3,nsym)=Fractional translations.
!!  ngfft(18)=Information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  err(3,nsym)=The max error for each symmetry. (given in terms of the FFT vectors)
!!  isok=.FALSE. if the FFT mesh does not fulfil all symmetry properties of the crystal.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!
!! SOURCE

subroutine fft_check_rotrans(nsym,symrel,tnons,ngfft,err,isok)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_check_rotrans'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 logical,intent(out) :: isok
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(out) :: err(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym,ix,iy,iz,jx,jy,jz,ngfft1,ngfft2,ngfft3
 !character(len=500) :: msg
!arrays
 integer :: Rm1(3,3,nsym),r1_FFT(3),red2fft(3,3)
 real(dp) :: Rm1_FFT(3,3,nsym),fft2red(3,3),r2_FFT(3),tnons_FFT(3,nsym)

! *************************************************************************

 ! === Precalculate R^-1 and fractional translations in FFT coordinates ===
 ngfft1=ngfft(1)
 ngfft2=ngfft(2) 
 ngfft3=ngfft(3)

 red2fft=RESHAPE((/ngfft1,0,0,0,ngfft2,0,0,0,ngfft3/),(/3,3/))
 fft2red=RESHAPE((/(one/ngfft1),zero,zero,zero,(one/ngfft2),zero,zero,zero,(one/ngfft3)/),(/3,3/))
 !
 ! === For a fully compatible mesh, each Rm1_FFT should be integer ===
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),Rm1(:,:,isym))
   Rm1(:,:,isym)=TRANSPOSE(Rm1(:,:,isym))
   Rm1_FFT(:,:,isym)=MATMUL(MATMUL(red2fft,Rm1(:,:,isym)),fft2red)
   tnons_FFT(:,isym)=MATMUL(red2fft,tnons(:,isym))
 end do

 err(:,:)=smallest_real 
 do iz=0,ngfft3-1
   R1_FFT(3)=DBLE(iz)
   do iy=0,ngfft2-1
     R1_FFT(2)=DBLE(iy)
     do ix=0,ngfft1-1
       R1_FFT(1)=DBLE(ix)
       do isym=1,nsym  ! Form R^-1 (r-\tau) in the FFT basis ===
         R2_FFT(:)=MATMUL(Rm1_FFT(:,:,isym),R1_FFT(:)-tnons_FFT(:,isym))
         jx=NINT(R2_FFT(1)); err(1,isym)=MAX(err(1,isym),ABS(R2_FFT(1)-jx)/ngfft1)
         jy=NINT(R2_FFT(2)); err(2,isym)=MAX(err(2,isym),ABS(R2_FFT(2)-jy)/ngfft2)
         jz=NINT(R2_FFT(3)); err(3,isym)=MAX(err(3,isym),ABS(R2_FFT(3)-jz)/ngfft3)
       end do 
     end do 
   end do 
 end do

 isok=.TRUE.
 do isym=1,nsym 
   if (ANY(err(:,isym)>tol6)) then
     isok=.FALSE.
     !write(msg,'(a,i3,a,3es14.6)')' symmetry ',isym,') not compatible with FFT grid, error ',err(:,isym)
     !MSG_WARNING(msg)
   end if
 end do

end subroutine fft_check_rotrans
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/rotate_FFT_mesh
!! NAME
!! rotate_FFT_mesh
!!
!! FUNCTION
!!  Find the FFT index of $ R{-1}(r-\tau) $ for each point in the FFT box.
!!  $R$ is a symmetry operation in real space, $\tau$ is the associated
!!  fractional translation.
!!
!! INPUTS
!!  nsym=Number of symmetries.
!!  symrel(3,3,nsym)=Symmetries in real space in reduced coordinates.
!!  tnons(3,nsym)=Fractional translations.
!!  ngfft(18)=Information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  irottb(ngfftot,nsym)=Indeces of $R^{-1}(r-\tau)$ in the FFT box.
!!  iscompatibleFFT=.FALSE. if the FFT mesh does not fulfil all symmetry properties of the crystal.
!!
!! NOTES
!!  The evaluation of the rotated point $R^{-1}(r-\tau)$ is done using real arithmetic.
!!  As a consequence, if the FFT mesh does not fulfil the symmetry properties 
!!  of the crystal, the array irottb will contain the index of the FFT point which 
!!  is the closest one to $R^{-1}(r-\tau)$. This might lead to inaccuracies in the 
!!  final results, in particular in the description of degenerate states.
!!
!! PARENTS
!!      bethe_salpeter,calc_sigc_me,calc_sigx_me,cchi0q0_intraband
!!      classify_bands,cohsex_me,m_wfs,rdm,screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_FFT_mesh(nsym,symrel,tnons,ngfft,irottb,iscompatibleFFT)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_FFT_mesh'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 logical,intent(out) :: iscompatibleFFT
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(in) :: ngfft(18)
 integer,intent(out) :: irottb(ngfft(1)*ngfft(2)*ngfft(3),nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: ir1,isym,ix,iy,iz,jx,jy,jz,ngfft1,ngfft2,ngfft3
 !character(len=500) :: msg
!arrays
 integer :: Rm1(3,3,nsym),r1_FFT(3),red2fft(3,3)
 real(dp) :: Rm1_FFT(3,3,nsym),err(3,nsym),fft2red(3,3),r2_FFT(3)
 real(dp) :: tnons_FFT(3,nsym)

! *************************************************************************

 ! === Precalculate R^-1 and fractional translations in FFT coordinates ===
 ngfft1=ngfft(1)
 ngfft2=ngfft(2) 
 ngfft3=ngfft(3)

 red2fft=RESHAPE((/ngfft1,0,0,0,ngfft2,0,0,0,ngfft3/),(/3,3/))
 fft2red=RESHAPE((/(one/ngfft1),zero,zero,zero,(one/ngfft2),zero,zero,zero,(one/ngfft3)/),(/3,3/))
 !
 ! === For a fully compatible mesh, each Rm1_FFT should be integer ===
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),Rm1(:,:,isym))
   Rm1(:,:,isym)=TRANSPOSE(Rm1(:,:,isym))
   Rm1_FFT(:,:,isym)=MATMUL(MATMUL(red2fft,Rm1(:,:,isym)),fft2red)
   tnons_FFT(:,isym)=MATMUL(red2fft,tnons(:,isym))
 end do

 err(:,:)=smallest_real 
 do iz=0,ngfft3-1
   R1_FFT(3)=DBLE(iz)
   do iy=0,ngfft2-1
     R1_FFT(2)=DBLE(iy)
     do ix=0,ngfft1-1
       R1_FFT(1)=DBLE(ix)
       ir1=1+ix+iy*ngfft1+iz*ngfft1*ngfft2
       do isym=1,nsym
         ! === Form R^-1 (r-\tau) in the FFT basis ===
         R2_FFT(:)=MATMUL(Rm1_FFT(:,:,isym),R1_FFT(:)-tnons_FFT(:,isym))
         jx=NINT(R2_FFT(1)); err(1,isym)=MAX(err(1,isym),ABS(R2_FFT(1)-jx)/ngfft1)
         jy=NINT(R2_FFT(2)); err(2,isym)=MAX(err(2,isym),ABS(R2_FFT(2)-jy)/ngfft2)
         jz=NINT(R2_FFT(3)); err(3,isym)=MAX(err(3,isym),ABS(R2_FFT(3)-jz)/ngfft3)
         jx=MODULO(jx,ngfft1)
         jy=MODULO(jy,ngfft2)
         jz=MODULO(jz,ngfft3)
         irottb(ir1,isym)=1+jx+jy*ngfft1+jz*ngfft1*ngfft2
       end do 
     end do 
   end do 
 end do

 iscompatibleFFT=.TRUE.
 do isym=1,nsym 
   if (ANY(err(:,isym)>tol6)) then
     iscompatibleFFT=.FALSE.
     !write(msg,'(a,i3,a,3es14.6)')' symmetry ',isym,') not compatible with FFT grid, error ',err(:,isym)
     !MSG_WARNING(msg)
   end if
 end do

end subroutine rotate_FFT_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/cigfft
!! NAME
!! cigfft
!!
!! FUNCTION
!! For each of the (2*nG0sh+1)**3 vectors G0 around the origin, 
!! calculate G-G0 and its FFT index number for all the NPWVEC vectors G.
!!
!! INPUTS
!! mG0(3)= For each reduced direction gives the max G0 component to account for umklapp processes.
!! npwvec=Number of plane waves
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! gvec(3,npwvec)=Reduced coordinates of G vectors.
!!
!! OUTPUT
!! igfft(npwvec,2*mG0(1)+1,2*mG0(2)+1,2*mG0(3)+1)=For each G, and each G0 vector, 
!!  it gives the FFT grid index of the G-G0 vector.
!! ierr=Number of G-G0 vectors falling outside the inout FFT box.
!!
!! PARENTS
!!      m_oscillators,rdm
!!
!! CHILDREN
!!
!! SOURCE

subroutine cigfft(mG0,npwvec,ngfft,gvec,igfft,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cigfft'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwvec
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: gvec(3,npwvec)
 integer,intent(in) :: mg0(3),ngfft(18)
 integer,intent(out) :: igfft(npwvec,2*mg0(1)+1,2*mg0(2)+1,2*mg0(3)+1)

!Local variables ------------------------------
!scalars
 integer :: gmg01,gmg02,gmg03,ig,ig01,ig02,ig03,n1,n2,n3 
 character(len=500) :: msg
!arrays
 integer :: gmg0(3)
!************************************************************************

 DBG_ENTER("COLL")

 if (ANY(mg0<0)) then
   write(msg,'(a,3i4)')' Found negative value of mg0= ',mg0
   MSG_BUG(msg)
 end if 

 n1=ngfft(1) 
 n2=ngfft(2) 
 n3=ngfft(3) 
 ierr=0

 do ig=1,npwvec
   do ig01=-mg0(1),mG0(1)
     gmg0(1) = gvec(1,ig)-ig01
     do ig02=-mg0(2),mg0(2)
       gmg0(2) = gvec(2,ig)-ig02
       do ig03=-mg0(3),mg0(3)
         gmg0(3) = gvec(3,ig)-ig03
         ! === Calculate FFT index of G-G0 ===
         ! * Consider possible wrap around errors.
         gmg01=MODULO(gmg0(1),n1)
         gmg02=MODULO(gmg0(2),n2)
         gmg03=MODULO(gmg0(3),n3)
         igfft(ig,ig01+mg0(1)+1,ig02+mg0(2)+1,ig03+mg0(3)+1) = 1+gmg01+gmg02*n1+gmg03*n1*n2
         if ( ANY(gmg0>ngfft(1:3)/2) .or. ANY(gmg0<-(ngfft(1:3)-1)/2) ) then
           igfft(ig,ig01+mg0(1)+1,ig02+mg0(2)+1,ig03+mg0(3)+1) = 0
           ierr=ierr+1
         end if
       end do
     end do
   end do
 end do !ig

 if (ierr/=0) then 
   write(msg,'(a,i4,3a)')&
&    ' Found ',ierr,' G-G0 vectors falling outside the FFT box. ',ch10,&
&    ' igfft will be set to zero for these particular G-G0 '
   MSG_WARNING(msg)
 end if

 DBG_EXIT("COLL")

end subroutine cigfft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/ig2gfft
!! NAME
!!  ig2gfft
!!
!! FUNCTION
!!  Return the reduced component of a G-vector in the FFT mesh starting from is index.
!!
!! INPUTS
!!  ig = The index >=1, <=ng
!!  ng = The number of FFT points along this direction.
!!
!! OUTPUT
!!  gc = The reduced component
!!
!! SOURCE

function ig2gfft(ig,ng) result (gc)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ig2gfft'
!End of the abilint section

 integer,intent(in) :: ig,ng
 integer :: gc

!************************************************************************

 ! Use the following indexing (N means ngfft of the adequate direction)
 ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gc
 ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index ig
 !
 if (ig<=0 .or. ig > ng) then 
   MSG_BUG("ig<=0 .or. ig > ng")
 end if

 if ( ig  > ng/2 + 1) then 
   gc = ig - ng -1
 else 
   gc = ig -1 
 end if

end function ig2gfft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/g2ifft
!! NAME
!!  g2ifft
!!
!! FUNCTION
!! Returns the index of G in the FFT box from its reduced coordinates. 0 if not in the BOX.
!!
!! INPUTS
!!  gg(3)=Reduced coordinated of the G vector.
!!  ngfft(18) = Info on the FFT box.
!!
!! OUTPUT
!!  gidx=Index in the FFT box. 0 if G is outside the box.
!!
!! SOURCE

function g2ifft(gg,ngfft) result (gidx)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'g2ifft'
!End of the abilint section

 integer :: gidx
!arrays
 integer,intent(in) :: gg(3),ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: n1,n2,n3,ig1,ig2,ig3

!************************************************************************

 ! Use the following indexing (N means ngfft of the adequate direction)
 ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gg
 ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index
 !
 if ( ANY(gg>ngfft(1:3)/2) .or. ANY(gg<-(ngfft(1:3)-1)/2) ) then ! out of the box.
   gidx=0 
 else 
   n1=ngfft(1)
   n2=ngfft(2)
   n3=ngfft(3)
                     
   ig1=MODULO(gg(1),n1)
   ig2=MODULO(gg(2),n2)
   ig3=MODULO(gg(3),n3)

   gidx = 1 + ig1 + n1*(ig2+ig3*n2)
 end if

end function g2ifft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/get_gftt
!! NAME
!!  get_gftt
!!
!! FUNCTION
!!  Returns the set of G-vectors in the FFT mesh and the maximal kinetic energy of k+G.
!!
!! INPUTS
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! kpt(3)=input k vector (reduced coordinates --in terms of reciprocal lattice primitive translations)
!! gmet(3,3)=reciprocal space metric (bohr^-2)
!!
!! OUTPUT
!!  gsq_max=Max value of (k+G)^2 for G in the FFT box
!!  gfft(3,nfft_tot) = The reduced components of the G in the FFT mesh (nfft_tot=PRODUCT(ngfft(1:3))
!!
!! PARENTS
!!      bethe_salpeter,calc_sigc_me,cchi0,cchi0q0,check_completeness,cohsex_me
!!      screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_gftt(ngfft,kpt,gmet,gsq_max,gfft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_gftt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: gsq_max
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(out) :: gfft(3,ngfft(1)*ngfft(2)*ngfft(3))
 real(dp),intent(in) :: kpt(3),gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: ifft,ig1,ig2,ig3,g1,g2,g3,i1,i2,i3
 real(dp) :: dsq

!************************************************************************

 dsq(i1,i2,i3)=gmet(1,1)*(kpt(1)+dble(i1))**2&
& +gmet(2,2)*(kpt(2)+dble(i2))**2&
& +gmet(3,3)*(kpt(3)+dble(i3))**2&
& +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
& +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
& +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))

 ifft=0; gsq_max=smallest_real
 do ig3=1,ngfft(3)    ; g3 = ig2gfft(ig3,ngfft(3))
   do ig2=1,ngfft(2)  ; g2 = ig2gfft(ig2,ngfft(2))
     do ig1=1,ngfft(1); g1 = ig2gfft(ig1,ngfft(1))
       ifft=ifft+1
       gfft(1,ifft) = g1
       gfft(2,ifft) = g2
       gfft(3,ifft) = g3
       gsq_max = MAX(dsq(ig1,ig2,ig3),gsq_max)
     end do
   end do
 end do

end subroutine get_gftt
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/print_ngfft
!! NAME
!! print_ngfft
!!
!! FUNCTION
!!  Print the content of ngfft(18) in explicative format.
!!
!! INPUTS
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft.
!!  [unit]=unit number for output (defaults to std_out).
!!  [prtvol]=verbosity level (defaults to 0).
!!  [mode_paral]=either "COLL" or "PERS" ("COLL" is default).
!!
!! OUTPUT
!!  Only writing 
!!
!! PARENTS
!!      bethe_salpeter,getng,m_fft_prof,m_wfs,screening,setup_bse,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_ngfft(ngfft,header,unit,mode_paral,prtvol)

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_ngfft'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=*),intent(in),optional :: header
 character(len=4),intent(in),optional :: mode_paral
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_prtvol=0;       if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='COLL';  if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=ch10//' ==== FFT mesh description (ngfft) ==== '
 if (PRESENT(header)) msg=ch10//' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(2(a,3i5,a),a,i5,2a,i5)')&
&  '  FFT mesh divisions ........................ ',ngfft(1),ngfft(2),ngfft(3),ch10,&
&  '  Augmented FFT divisions ................... ',ngfft(4),ngfft(5),ngfft(6),ch10,&
&  '  FFT algorithm ............................. ',ngfft(7),ch10,&
&  '  FFT cache size ............................ ',ngfft(8)
 call wrtout(my_unt,msg,my_mode)

 if (my_prtvol>0) then
   write(msg,'(6(a,i5,a),a,4i5)')&
&    '  FFT parallelization level ................. ',ngfft(9),ch10,&
&    '  Number of processors in my FFT group ...... ',ngfft(10),ch10,&
&    '  Index of me in my FFT group ............... ',ngfft(11),ch10,&
&    '  No of xy planes in R space treated by me .. ',ngfft(12),ch10,&
&    '  No of xy planes in G space treated by me .. ',ngfft(13),ch10,&
&    '  MPI communicator for FFT .................. ',ngfft(14),ch10,&
&    '  Value of ngfft(15:18) ..................... ',ngfft(15:18)
   call wrtout(my_unt,msg,my_mode)
 end if

end subroutine print_ngfft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/fftalg_info
!! NAME
!! fftalg_info
!!
!! FUNCTION
!!  Returns info on the FFT library specified by fftalg (ngfft(7))
!!
!! INPUTS
!!  fftalg=Input variable.
!!
!! OUTPUT
!!  library=String with the name of FFT library
!!  cplex_mode= String defining whether the FFT library supports real<-->complex transforms.
!!  padding_mode=Padding mode.
!!
!! PARENTS
!!      m_fft_prof
!!
!! CHILDREN
!!
!! SOURCE

subroutine fftalg_info(fftalg,library,cplex_mode,padding_mode)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftalg_info'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg 
 character(len=*),intent(out) :: library,cplex_mode,padding_mode

!Local variables-------------------------------
!scalars
 integer :: fftalga,fftalgb,fftalgc

! *************************************************************************

 library       = "Unknown"
 cplex_mode    = "Unknown"
 padding_mode  = "Unknown"

 fftalga=fftalg/100
 if (fftalga>0 .and. fftalga<=FFTALGA_SIZE) library    = fftalga2name(fftalga)

 fftalgb=mod(fftalg,100)/10
 if (fftalgb>=0 .and. fftalgb<=FFTALGB_SIZE) cplex_mode = fftalgb2name(fftalgb)

 fftalgc=mod(fftalg,10)
 if (fftalgc>=0 .and. fftalgc<=FFTALGC_SIZE) padding_mode = fftalgc2name(fftalgc)

end subroutine fftalg_info
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/ceigr
!! NAME
!! ceigr
!!
!! FUNCTION
!!  Helper function to calculate e^{iG.r} on the FFT mesh.
!!
!! INPUTS
!!  gg(3)=G vector in reduced coordinates.
!!  nfft=Total number of points in the FFT mesh.
!!  ngfft(18)=information about 3D FFT, 
!!
!! OUTPUT
!!  ceigr(nfft)=e^{ik.r} on the FFT mesh.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ceigr(gg,nfft,ngfft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ceigr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
!arrays
 integer,intent(in) :: gg(3)
 integer,intent(in) :: ngfft(18)
 complex(dpc) :: ceigr(nfft)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,fft_idx
 real(dp) :: gdotr

! *************************************************************************

 if (ALL(gg==0)) then
   ceigr=cone; RETURN
 end if 

 fft_idx=0
 do iz=0,ngfft(3)-1
   do iy=0,ngfft(2)-1
     do ix=0,ngfft(1)-1
       gdotr= two_pi*( gg(1)*(ix/DBLE(ngfft(1))) &
&                     +gg(2)*(iy/DBLE(ngfft(2))) &
&                     +gg(3)*(iz/DBLE(ngfft(3))) )
       fft_idx = fft_idx+1
       ceigr(fft_idx)=DCMPLX(DCOS(gdotr),DSIN(gdotr))
     end do
   end do
 end do

end function ceigr
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/eigr
!! NAME
!! eigr
!!
!! FUNCTION
!!  Helper function to calculate e^{iG.r} on the FFT mesh.
!!
!! INPUTS
!!  gg(3)=G vector in reduced coordinates.
!!  nfft=Total number of points in the FFT mesh.
!!  ngfft(18)=information about 3D FFT, 
!!
!! OUTPUT
!!  eigr(2*nfft)=e^{ig.r} on the FFT mesh.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function eigr(gg,nfft,ngfft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
!arrays
 integer,intent(in) :: gg(3)
 integer,intent(in) :: ngfft(18)
 real(dp) :: eigr(2*nfft)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,fft_idx
 real(dp) :: gdotr

! *************************************************************************

 if (ALL(gg==0)) then
   eigr(1:2*nfft:2)=one
   eigr(2:2*nfft:2)=zero
 end if 

 fft_idx=1
 do iz=0,ngfft(3)-1
   do iy=0,ngfft(2)-1
     do ix=0,ngfft(1)-1
       gdotr= two_pi*( gg(1)*(ix/DBLE(ngfft(1))) &
&                     +gg(2)*(iy/DBLE(ngfft(2))) &
&                     +gg(3)*(iz/DBLE(ngfft(3))) )
       eigr(fft_idx  )=DCOS(gdotr)
       eigr(fft_idx+1)=DSIN(gdotr)
       fft_idx = fft_idx+2
     end do
   end do
 end do

end function eigr
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/ceikr
!! NAME
!! ceikr
!!
!! FUNCTION
!!  Helper function to calculate e^{ik.r} on the FFT mesh.
!!
!! INPUTS
!!  kk(3)=k-point in reduced coordinates.
!!  nfft=Total number of points in the FFT mesh.
!!  ngfft(18)=information about 3D FFT, 
!!
!! OUTPUT
!!  ceikr(nfft)=e^{ik.r} on the FFT mesh.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ceikr(kk,nfft,ngfft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ceikr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
!arrays
 real(dp),intent(in) :: kk(3)
 integer,intent(in) :: ngfft(18)
 complex(dpc) :: ceikr(nfft)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,fft_idx
 real(dp) :: kdotr

! *************************************************************************

 !if (ALL(ABS(kk<tol12)) then
 !  ceikr=cone; RETURN
 !end if 

 fft_idx=0
 do iz=0,ngfft(3)-1
   do iy=0,ngfft(2)-1
     do ix=0,ngfft(1)-1
       kdotr= two_pi*( kk(1)*(ix/DBLE(ngfft(1))) &
&                     +kk(2)*(iy/DBLE(ngfft(2))) &
&                     +kk(3)*(iz/DBLE(ngfft(3))) )
       fft_idx = fft_idx+1
       ceikr(fft_idx)=DCMPLX(DCOS(kdotr),DSIN(kdotr))
     end do
   end do
 end do

end function ceikr
!!***

!----------------------------------------------------------------------


END MODULE m_fft_mesh
!!***
