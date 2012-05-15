!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhotov
!! NAME
!! rhotov
!!
!! FUNCTION
!! This routine is called to compute, from a given total density
!! the trial (local) potential and the residual potential.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | fixmom=input variable that governs fixed moment calculation
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!   | ntypat=number of types of atoms in unit cell.
!!   | occopt=option for occupancies
!!   | typat(natom)=type (integer) for each atom
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=option for the computation of additional energies
!!  optres=0: the trial potential residual is computed ; the input potential value is kept
!!         1: the new value of the trial potential is computed in place of the input value
!!  optxc=option to be used for the call to rhohxc
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  ucvol = unit cell volume (Bohr**3)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp(nfft)=array for holding local psp
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  ==== if optres==0
!!    vtrial(nfft,nspden)= old value of trial potential
!!
!! OUTPUT
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_hartree=Hartree part of total energy (hartree units)
!!   | e_xc=exchange-correlation energy (hartree)
!!  ==== if dtset%usewvl==1
!!   | e_vxc=the energy of the exchange-correlation potential
!!  ==== if optene==0.or.2
!!   | e_localpsp=local psp energy (hartree)
!!  ==== if optene==1.or.2
!!   | e_xcdc=exchange-correlation double-counting energy (hartree)
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if optxc==2.
!!  strsxc(6)=xc contribution to stress tensor (hartree/bohr^3)
!!  vxc(nfft,nspden)=Vxc(r) (already computed above; gets recomputed below too)
!!  vxcavg=mean of the vxc potential
!!  ==== if optres==0
!!    vresidnew(nfft,nspden)=potential residual
!!    vnew_mean(nspden)=mean of the potential formed from vpsp, vhartr and vxc, might be spin-dependent
!!    vres_mean(nspden)=mean of the potential residual, might be spin-dependent
!!    vres2=square of the norm of the residual
!!  vxctau(nfftf,dtset%nspden*dtset%usekden,4)=derivative of XC energy density with respect to
!!    kinetic energy density (metaGGA cases) (optional output)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  vhartr(nfft)=array for holding Hartree potential
!!  ==== if optres==1
!!    vtrial(nfft,nspden)= new value of trial potential
!!
!! NOTES
!!  In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfft,ngfft,mgfft) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      prctfvw1,prctfvw2,prctfw3,scfcv
!!
!! CHILDREN
!!      dotprod_vn,mean_fftr,psolver_rhohxc,rhohxc,rhohxcpositron,sqnorm_v
!!      timab,wvl_newvtr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rhotov(dtset,energies,gprimd,gsqcut,kxc,mpi_enreg,nfft,ngfft,&
&  nhat,nhatgr,nhatgrdim,nkxc,vresidnew,n3xccc,optene,optres,optxc,&
&  rhog,rhor,rprimd,strsxc,ucvol,usepaw,usexcnhat,&
&  vhartr,vnew_mean,vpsp,vres_mean,vres2,vtrial,vxcavg,vxc,wvl,xccc3d,&
&  electronpositron,taug,taur,vxctau) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_energies, only : energies_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhotov'
 use interfaces_18_timing
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_62_poisson
 use interfaces_67_common, except_this_one => rhotov
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n3xccc,nfft,nhatgrdim,nkxc,optene,optres,optxc,usepaw
 integer,intent(in) :: usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: vres2,vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim),rhog(2,nfft)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),vhartr(nfft),vpsp(nfft)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc)
 real(dp),intent(out) :: kxc(nfft,nkxc),strsxc(6),vnew_mean(dtset%nspden)
 real(dp),intent(out) :: vres_mean(dtset%nspden),vresidnew(nfft,dtset%nspden)
 real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
 real(dp),intent(in),optional :: taur(nfft,dtset%nspden*dtset%usekden)
 real(dp),intent(out),optional :: vxctau(nfft,dtset%nspden*dtset%usekden,4)

!Local variables-------------------------------
!scalars
 integer :: nk3xc
!integer :: jj,kk
 integer :: ifft,ipositron,ispden,nfft_loc,nfftot,offset
 real(dp) :: doti,e_xcdc_vxctau,ucvol_local
!real(dp) :: sum
 logical :: with_vxctau
!arrays
 real(dp) :: tsec(2),vmean(dtset%nspden)
 real(dp) :: vzeeman(dtset%nspden)
 real(dp),allocatable :: rhowk(:,:),vnew(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)'++++++++++++++++'
!write(std_out,*)'->rhotov : enter'
!sum=0.0
!do jj=1,dtset%nspden
!do kk=1,nfft
!sum=sum+abs(rhor(kk,jj))
!end do
!end do
!write(std_out,*) 'SUM(rhor)',sum
!sum=0.0
!do jj=1,dtset%nspden
!do kk=1,nfft
!sum=sum+abs(vtrial(kk,jj))
!end do
!end do
!write(std_out,*) 'SUM(vtrial)',sum
!sum=0.0
!do jj=1,dtset%nspden
!do kk=1,nfft
!sum=sum+abs(vxc(kk,jj))
!end do
!end do
!write(std_out,*) 'SUM(vxc)',sum
!ENDDEBUG

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = .false.
 if (present(vxctau) .and. present(taur) .and. dtset%usekden /= 0) with_vxctau = .true.

 call timab(940,1,tsec)

!Get size of FFT grid
 if (dtset%usewvl == 0) then
   nfftot=ngfft(1)*ngfft(2)*ngfft(3)
   ucvol_local = ucvol
 else
!  We need to tune the volume when wavelets are used because, not
!  all FFT points are used.
#if defined HAVE_DFT_BIGDFT
   nfftot = wvl%Glr%d%n1i * wvl%Glr%d%n2i * wvl%Glr%d%n3i
#endif
   ucvol_local = product(wvl%h) * real(nfftot, dp) / real(8, dp)
 end if

 ipositron=0
 if (present(electronpositron)) ipositron=electronpositron_calctype(electronpositron)

!------Compute Hartree and xc potentials----------------------------------

 if (ipositron/=1) then
!  Compute xc potential (separate up and down if spin-polarized)
   if (dtset%icoulomb == 0) then
!    Use the periodic solver to compute Hxc.
     nk3xc=1
     call timab(941,1,tsec)
     if (ipositron==0) then
       if(with_vxctau)then
         call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&         rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur,vxctau=vxctau)
       else
         call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&         rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur)
       end if
     else
       if(with_vxctau)then
         call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&         rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,&
&         electronpositron=electronpositron,taug=taug,taur=taur,vxctau=vxctau)
       else
         call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&         rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,&
&         electronpositron=electronpositron,taug=taug,taur=taur)
       end if
     end if
     call timab(941,2,tsec)
     call timab(942,1,tsec)
     call dotprod_vn(1,rhor,energies%e_hartree,doti,mpi_enreg,nfft,nfftot,1,1,vhartr,ucvol)
     energies%e_hartree=half*energies%e_hartree
     call timab(942,2,tsec)
   else
!    Use the free boundary solver.
     call timab(943,1,tsec)
     call PSolver_rhohxc(dtset, energies%e_hartree, energies%e_xc, energies%e_vxc, &
&     mpi_enreg, rhor, rprimd, vhartr, vxc, vxcavg, wvl)
     call timab(943,2,tsec)
   end if
 else
   call timab(944,1,tsec)
   energies%e_hartree=zero;energies%e_xc=zero
   call rhohxcpositron(electronpositron,gprimd,kxc,mpi_enreg,nfft,ngfft,nhat,nkxc,dtset%nspden,n3xccc,&
&   dtset%paral_kgb,rhor,strsxc,ucvol,usexcnhat,usepaw,vhartr,vxc,vxcavg,xccc3d,dtset%xc_denpos)
   call timab(944,2,tsec)
 end if

 call timab(945,1,tsec)
 if (ipositron/=0) then
   call dotprod_vn(1,rhor,electronpositron%e_hartree,doti,mpi_enreg,&
&   nfft,nfftot,1,1,electronpositron%vha_ep,ucvol)
   vhartr=vhartr+electronpositron%vha_ep
 end if

!------Compute parts of total energy depending on potentials--------
 if (optene==0.or.optene==2) then

!  Compute local psp energy energies%e_localpsp
   call dotprod_vn(1,rhor,energies%e_localpsp,doti,mpi_enreg,nfft,nfftot,1,1,vpsp,ucvol_local)

 end if

 if (optene==1.or.optene==2) then

!  Compute double-counting XC energy energies%e_xcdc
   if (ipositron/=1) then
     if (usepaw==0.or.usexcnhat/=0) then
       call dotprod_vn(1,rhor,energies%e_xcdc,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local)
       if(with_vxctau)then
         call dotprod_vn(1,taur,e_xcdc_vxctau,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxctau(:,:,1),ucvol_local)
         energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
       end if
     else
       ABI_ALLOCATE(rhowk,(nfft,dtset%nspden))
       rhowk=rhor-nhat
       call dotprod_vn(1,rhowk,energies%e_xcdc,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local)
       ABI_DEALLOCATE(rhowk)
     end if
     if (ipositron==2) energies%e_xcdc=energies%e_xcdc-electronpositron%e_xcdc
   else
     energies%e_xcdc=zero
   end if

 end if

!------Produce residual vector and square norm of it-------------
!(only if requested ; if optres==0)

!Set up array for Zeeman field
 vzeeman(:) = zero
 if(dtset%nspden==2)vzeeman(2) = dtset%zeemanfield(3) ! For collinear ispden=2 is rho_up only 
 if(dtset%nspden==4)then
   do ispden=2,4
     vzeeman(ispden)=dtset%zeemanfield(ispden-1)
   end do !ispden
 end if

 if (optres==0) then
!  Compute potential residual
   ABI_ALLOCATE(vnew,(nfft,dtset%nspden))
   vmean(:)=zero ; vnew_mean(:)=zero

   if (dtset%usewvl == 0) then
     do ispden=1,min(dtset%nspden,2)
!      $OMP PARALLEL DO PRIVATE(ifft) &
!      $OMP&SHARED(ispden,nfft,dtset%nspden,vhartr,vnew,vpsp,vresidnew,vtrial,vxc)
       do ifft=1,nfft
         vnew(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)+vzeeman(ispden)
         vresidnew(ifft,ispden)=vnew(ifft,ispden)-vtrial(ifft,ispden)
       end do
!      $OMP END PARALLEL DO
     end do
     if(dtset%nspden==4)then
       do ispden=3,4
!        $OMP PARALLEL DO PRIVATE(ifft) &
!        $OMP&SHARED(ispden,nfft,vresidnew,vtrial,vxc)
         do ifft=1,nfft
           vnew(ifft,ispden)=vxc(ifft,ispden)+vzeeman(ispden)
           vresidnew(ifft,ispden)=vxc(ifft,ispden)-vtrial(ifft,ispden)
         end do
!        $OMP END PARALLEL DO
       end do
     end if
     offset   = 0
     nfft_loc = nfft
   else
     call wvl_newvtr(dtset, mpi_enreg, nfft_loc, offset, vhartr, vpsp, vnew, vxc, wvl)
     vresidnew = vnew - vtrial
   end if
!  Compute mean values of potential and residual
   call mean_fftr(vnew(1+offset, 1),vnew_mean,mpi_enreg,nfft_loc,nfftot,dtset%nspden)
   call mean_fftr(vresidnew(1+offset, 1),vmean,mpi_enreg,nfft_loc,nfftot,dtset%nspden)

!  DEBUG
!  write(std_out,*)' rhotov : large values of the potential: ifft, vnew,vhartr,vpsp,vxc,vzeeman '
!  do ifft=1,nfft_loc
!  if(abs(vnew(ifft, 1))>1.0e1)then
!  write(std_out,'(i5,5es16.6)')ifft,vnew(ifft, 1),vhartr(ifft),vpsp(ifft),vxc(ifft,1),vzeeman(1)
!  endif
!  enddo
!  ENDDEBUG

   ABI_DEALLOCATE(vnew)

!  Subtract the mean of the residual
!  Must take into account fixed occupation number in case of spin-polarized
   do ispden=1,dtset%nspden
     if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&     abs(dtset%fixmom+99.99_dp)<1.0d-10)then
       vres_mean(ispden)=(vmean(1)+vmean(2))*half
     else
       vres_mean(ispden)=vmean(ispden)
     end if
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(ispden,nfft,vresidnew,vres_mean)
     do ifft=1,nfft
       vresidnew(ifft,ispden)=vresidnew(ifft,ispden)-vres_mean(ispden)
     end do
!    $OMP END PARALLEL DO
   end do

!  Compute square norm vres2 of potential residual vresid
   call sqnorm_v(1,mpi_enreg,nfft_loc,vres2,dtset%nspden,optres,vresidnew(1+offset, 1))

 else    ! optres==0


!  ------Produce new value of trial potential-------------
!  (only if requested ; if optres==1)

   if (dtset%usewvl == 0) then
     do ispden=1,min(dtset%nspden,2)
!      $OMP PARALLEL DO PRIVATE(ifft) &
!      $OMP&SHARED(ispden,nfft,dtset%nspden,vhartr,vnew,vpsp,vxc)
       do ifft=1,nfft
         vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)+vzeeman(ispden)
       end do
!      $OMP END PARALLEL DO
     end do
     if(dtset%nspden==4) then
       do ifft=1,nfft
         vtrial(ifft,3:4)=vxc(ifft,3:4)+vzeeman(3:4)
       end do 
     end if
   else
!    output offset and nfft_loc are unused here.
     call wvl_newvtr(dtset, mpi_enreg, nfft_loc, offset, vhartr, vpsp, vtrial, vxc, wvl)
   end if

 end if

 call timab(945,2,tsec)
 call timab(940,2,tsec)

!DEBUG
!write(std_out,*)' rhotov : exit '
!stop
!ENDDEBUG

end subroutine rhotov
!!***
