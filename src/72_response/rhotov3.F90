!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhotov3
!!
!! NAME
!! rhotov3
!!
!! FUNCTION
!! This routine is called to compute, from a given 1st-order total density
!! the trial (local) 1st-order potential and the residual potential.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order WF on FFT grid are REAL; if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of
!!    atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=type of the perturbation
!!  kxc(nfft,nkxc)=exchange-correlation kernel
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat1(cplex*nfft,nspden*usepaw)= -PAW only- 1st-order compensation density
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used
!!  optene=option for the computation of additional energies
!!  optres=0: the trial potential residual is computed ; the input potential value is kept
!!         1: the new value of the trial potential is computed in place of the input value
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rhog1(2,nfft)=RF electron density in reciprocal space
!!  rhor(nfft,nspden)=array for GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=RF electron density in real space (electrons/bohr**3).
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp1(cplex*nfft)=first-order derivative of the ionic potential
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!
!! OUTPUT
!!  ==== if optene>=0
!!    ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!    ehart1=1st-order Hartree part of 2nd-order total energy
!!    exc1=1st-order exchange-correlation part of 2nd-order total energy
!!  ==== if optene==0.or.2
!!    elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!!  ==== if optene==1.or.2
!!    To be completed
!!  ==== if optres==0
!!    vresid1(cplex*nfft,nspden)=potential residual
!!    vres2=square of the norm of the residual
!!
!! SIDE EFFECTS
!!  vhartr1(cplex*nfft)=1-order Hartree potential
!!  ==== if optres==1
!!    vtrial1(cplex*nfft,nspden)= new value of 1st-order trial potential
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      dotprod_vn,hartre,hartrestr,leave_new,mkvxc3,mkvxcstr3,sqnorm_v,timab
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


 subroutine rhotov3(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&
&           kxc,mpi_enreg,natom,nfft,ngfft,nhat1,nkxc,nspden,n3xccc,&
&           optene,optres,paral_kgb,qphon,rhog,rhog1,rhor,rhor1,&
&           rprimd,ucvol,usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,xccc3d1)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhotov3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_72_response, except_this_one => rhotov3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,n3xccc,natom,nfft,nkxc,nspden,optene
 integer,intent(in) :: optres,paral_kgb,usepaw,usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: ehart01,ehart1,elpsp1,exc1,vres2
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: nhat1(cplex*nfft,nspden*usepaw),qphon(3),rhog(2,nfft)
 real(dp),intent(in) :: rhog1(2,nfft),rhor(nfft,nspden)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rprimd(3,3),vpsp1(cplex*nfft)
 real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
 real(dp),intent(out) :: vhartr1(cplex*nfft),vresid1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft,ispden,nfftot,option
 real(dp) :: doti,elpsp10
 character(len=500) :: message
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp) :: tsec(20)
 real(dp),allocatable :: rhor1_wk(:,:),vhartr01(:),vxc1(:,:)
 real(dp),allocatable :: vxc10(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' rhotov3 : enter, stop '
!write(std_out,*)rprimd
!ENDDEBUG

!Tests
 if(nspden==4)then
   write(message, '(a,a,a,a)' )ch10,&
&   ' rhotov3 : ERROR -',ch10,&
&   '  Does not work yet for nspden=4 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(optene>0)then
   write(message, '(a,a,a,a)' )ch10,&
&   ' rhotov3 : BUG -',ch10,&
&   '  optene>0 not yet implemented !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 call timab(157,1,tsec)

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!------ Compute 1st-order Hartree potential (and energy) ----------------------

 call hartre(cplex,gmet,gsqcut,0,mpi_enreg,nfft,ngfft,paral_kgb,qphon,rhog1,vhartr1)
 if (optene>=0) then
   call dotprod_vn(cplex,rhor1,ehart1,doti,mpi_enreg,nfft,nfftot,1,1,vhartr1,ucvol)
 end if

 if (optene>=0) ehart01=zero
 if(ipert==natom+3 .or. ipert==natom+4) then
   ABI_ALLOCATE(vhartr01,(cplex*nfft))
   call hartrestr(gmet,gprimd,gsqcut,idir,ipert,mpi_enreg,natom,nfft,ngfft,paral_kgb,rhog,vhartr01)
   if (optene>=0) then
     call dotprod_vn(cplex,rhor1,ehart01,doti,mpi_enreg,nfft,nfftot,1,1,vhartr01,ucvol)
     ehart01=two*ehart01
     ehart1=ehart1+ehart01
   end if
!  Note that there is a factor 2.0_dp difference with the similar GS formula
   vhartr1(:)=vhartr1(:)+vhartr01(:)
   ABI_DEALLOCATE(vhartr01)
 end if

!------ Compute 1st-order XC potential (and energy) ----------------------
!(including the XC core correction)

 option=0;if (optene<0) option=1
 ABI_ALLOCATE(vxc1,(cplex*nfft,nspden))
 ABI_ALLOCATE(vxc10,(cplex*nfft,nspden))

 if(ipert==natom+3 .or. ipert==natom+4) then
   call mkvxcstr3(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,&
&   nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,vxc10,xccc3d1)
 else
   call mkvxc3(cplex,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&
&   paral_kgb,qphon,rhor1,rprimd,vxc10,xccc3d1)
 end if

 if (optene==0.or.optene==2) then
   if (usepaw==0) then
     call dotprod_vn(cplex,rhor1,elpsp10,doti,mpi_enreg,nfft,nfftot,nspden,1,vxc10,ucvol)
     call dotprod_vn(cplex,rhor1,elpsp1 ,doti,mpi_enreg,nfft,nfftot,1     ,1,vpsp1,ucvol)
   else
     if (usexcnhat/=0) then
       ABI_ALLOCATE(rhor1_wk,(cplex*nfft,1))
       rhor1_wk(:,1)=rhor1(:,1)-nhat1(:,1)
       call dotprod_vn(cplex,rhor1   ,elpsp10,doti,mpi_enreg,nfft,nfftot,nspden,1,vxc10,ucvol)
       call dotprod_vn(cplex,rhor1_wk,elpsp1 ,doti,mpi_enreg,nfft,nfftot,1     ,1,vpsp1,ucvol)
     else
       ABI_ALLOCATE(rhor1_wk,(cplex*nfft,nspden))
       rhor1_wk(:,:)=rhor1(:,:)-nhat1(:,:)
       call dotprod_vn(cplex,rhor1_wk,elpsp10,doti,mpi_enreg,nfft,nfftot,nspden,1,vxc10,ucvol)
       call dotprod_vn(cplex,rhor1_wk,elpsp1 ,doti,mpi_enreg,nfft,nfftot,1     ,1,vpsp1,ucvol)
     end if
     ABI_DEALLOCATE(rhor1_wk)
   end if
!  Note that there is a factor 2.0_dp difference with the similar GS formula
   elpsp1=two*(elpsp1+elpsp10)
 end if

!Compute XC contribution exc1 (except the XC core-correction)
 if (optene>=0) then
   option=2
   if (usepaw==1.and.usexcnhat==0) then
     ABI_ALLOCATE(rhor1_wk,(cplex*nfft,nspden))
     rhor1_wk(:,:)=rhor1(:,:)-nhat1(:,:)
     call mkvxc3(cplex,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&
&     paral_kgb,qphon,rhor1_wk,rprimd,vxc1,xccc3d1)
     ABI_DEALLOCATE(rhor1_wk)
   else
     call mkvxc3(cplex,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&
&     paral_kgb,qphon,rhor1,rprimd,vxc1,xccc3d1)
   end if
   call dotprod_vn(cplex,rhor1,exc1,doti,mpi_enreg,nfft,nfftot,nspden,1,vxc1,ucvol)
   if(n3xccc/=0 .or. ipert==natom+3 .or. ipert==natom+4)then
     vxc1(:,:)=vxc1(:,:)+vxc10(:,:)
   end if
 else
   vxc1=vxc10
 end if

 ABI_DEALLOCATE(vxc10)

!DEBUG (do not take away)
!Compute NSC energy ensc1 associated with rhor1 in vtrial1, for debugging purposes
!call dotprod_vn(cplex,rhor1,ensc1,doti,nfft,nfftot,nspden,1,vtrial1,ucvol)
!write(std_out,*)' ek0+eeig0+eloc0=',ek0+eeig0+eloc0
!write(std_out,*)' ensc1=',ensc1
!Compute NSC energy associated with vtrial1, for debugging purposes
!call dotprod_vn(cplex,rhor1,ensc1,doti,mpi_enreg,nfft,nfftot,nspden,1,vtrial1,ucvol)
!ensc1=ensc1+half*enl1
!write(std_out,*)' rhotov3 : check NSC energy, diff=',&
!&  ek0+edocc+eeig0+eloc0+enl0+ensc1
!write(std_out,*)' evarNSC=',ek0+edocc+eeig0+eloc0+enl0
!write(std_out,*)' ensc1,exc1=',ensc1,exc1
!ENDDEBUG

!Here, vhartr1 contains Hartree potential, vpsp1 contains local psp,
!while vxc1 contain xc potential

!------ Produce residual vector and square of norm of it -------------
!(only if requested ; if optres==0)

 if (optres==0) then
   do ispden=1,nspden
     do ifft=1,cplex*nfft
       vresid1(ifft,ispden)=vhartr1(ifft)+vxc1(ifft,ispden)+vpsp1(ifft)-vtrial1(ifft,ispden)
     end do
   end do

!  Compute square norm vres2 of potential residual vresid
   call sqnorm_v(cplex,mpi_enreg,nfft,vres2,nspden,optres,vresid1)

 else

!  ------ Produce new value of trial potential-------------
!  (only if requested ; if optres==1)

   do ispden=1,nspden
     do ifft=1,cplex*nfft
       vtrial1(ifft,ispden)=vhartr1(ifft)+vxc1(ifft,ispden)+vpsp1(ifft)
     end do
   end do

 end if

 ABI_DEALLOCATE(vxc1)

 call timab(157,2,tsec)

!DEBUG
!write(std_out,*)' rhotov3 : exit '
!stop
!ENDDEBUG

end subroutine rhotov3
!!***
