!{\src2tex{textfont=tt}}
!!****f* ABINIT/Psolver_rhohxc
!! NAME
!! Psolver_rhohxc
!!
!! FUNCTION
!! Given rho(r), compute Hartree potential considering the system as
!! an isolated one. This potential is obtained from the convolution
!! of 1/r and rho(r), treated in Fourier space. This method is a wrapper around
!! Psolver() developped for BigDFT.
!! It can compute the xc energy and potential if required. This computation is
!! built on the drivexc() routine of ABINIT but access it directly from real
!! space. The present routine is a real space counter part to rhohxc().
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=MPI-parallelisation information.
!!  rhor(nfft,nspden)=electron density in real space in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  enhartr=returned Hartree energy (hartree).
!!  enxc=returned exchange and correlation energy (hartree).
!!  envxc=returned energy of the Vxc potential (hartree).
!!  vhartr(nfft)=Hartree potential.
!!  vxc(nfft,nspden)=xc potential
!!  vxcavg=<Vxc>=unit cell average of Vxc = (1/ucvol) Int [Vxc(r) d^3 r].
!!
!! NOTE
!!  In PSolver, with nspden == 2, rhor(:,1) = density up and
!!                                rhor(:,2) = density down.
!!  But in ABINIT (dtset%usewvl != 1) rhor(:,1) = total density and
!!                                    rhor(:,2) = density up .
!!  In ABINIT (dtset%usewvl != 1), the same convention is used as in PSolver.
!!
!! PARENTS
!!      energy,rhotov,setvtr
!!
!! CHILDREN
!!      leave_new,mean_fftr,mkdenpos,psolver,psolver_kernel,wrtout,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine Psolver_rhohxc(dtset, enhartr, enxc, envxc, mpi_enreg, rhor, rprimd, &
     & vhartr, vxc, vxcavg, wvl)

 use m_profiling

  use defs_basis
  use defs_abitypes
  use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use poisson_solver
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Psolver_rhohxc'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_62_poisson, except_this_one => Psolver_rhohxc
!End of the abilint section

  implicit none

  !Arguments ------------------------------------
  !scalars
  real(dp), intent(out)         :: enxc, envxc, enhartr, vxcavg
  type(dataset_type),intent(in) :: dtset
  type(MPI_type),intent(inout)  :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  !arrays
  real(dp),intent(in)    :: rprimd(3,3)
  real(dp),intent(inout) :: rhor(dtset%nfft, dtset%nspden)
  real(dp),intent(out)   :: vhartr(dtset%nfft)
  real(dp),intent(out)   :: vxc(dtset%nfft, dtset%nspden)

  !Local variables-------------------------------
  !scalars
  integer :: comm,me,nfft_tot,nproc, i, iwarn=0, opt_mkdenpos=0
  real(dp) :: hgrid(3), tmpDown, tmpUp, tmpPot
  character(len=500) :: message
  character(len = 1) :: datacode, bndcode
  !arrays
  real(dp) :: vxcmean(1)
  real(dp), pointer :: kernel(:)

! *********************************************************************

 if (dtset%icoulomb == 0) then
!  The kernel is built with 'P'eriodic boundary counditions.
   bndcode = 'P'
 else if (dtset%icoulomb == 1) then
!  The kernel is built with 'F'ree boundary counditions.
   bndcode = 'F'
 else if (dtset%icoulomb == 2) then
!  The kernel is built with 'S'urface boundary counditions.
   bndcode = 'S'
 end if

 if(dtset%nspden > 2)then
   write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&   ' PSolver_rhohxc: BUG -',ch10,&
&   '  Only non-spin-polarised or collinear spin is allowed,',ch10,&
&   '  while the argument nspden = ', dtset%nspden
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)

!We get the kernel.
 call PSolver_kernel(dtset, 2, kernel, mpi_enreg, rprimd, wvl)
!If the kernel is not created, we do it now.
 if (.not.associated(kernel)) then
   call PSolver_kernel(dtset, 1, kernel, mpi_enreg, rprimd, wvl)
 end if

!We do the computation.
 write(message, "(A,A,A,3I6)") "Psolver_rhohxc(): compute potentials (Vhartree and Vxc)...", ch10, &
& " | dimension:", dtset%ngfft(1:3)
 call wrtout(std_out, message,'COLL')

!We save total rhor in vhartr
 vhartr(:)  = rhor(:, 1)

#if defined HAVE_DFT_BIGDFT
 if (dtset%usewvl == 0) then
!  In non-wavelet case, we change the rhor values.
   if (dtset%nspden == 2) then
     do i = 1, dtset%nfft, 1
!      We change rhor for PSolver call.
       tmpDown = rhor(i, 1) - rhor(i, 2)
       tmpUp   = rhor(i, 2)
       rhor(i, 1) = tmpUp
       rhor(i, 2) = tmpDown
     end do
   end if


!  Make the density positive everywhere (but do not care about gradients)
   call mkdenpos(iwarn, dtset%nfft, dtset%nspden, opt_mkdenpos, rhor, dtset%xc_denpos)

   hgrid(1) = rprimd(1, 1) / dtset%ngfft(1)
   hgrid(2) = rprimd(2, 2) / dtset%ngfft(2)
   hgrid(3) = rprimd(3, 3) / dtset%ngfft(3)
!  This may not work with MPI in the planewave code...
   call PSolver(bndcode, 'G', me, nproc, dtset%ngfft(1), &
&   dtset%ngfft(2), dtset%ngfft(3), dtset%ixc, hgrid(1), hgrid(2), hgrid(3), &
&   rhor, kernel, vxc, enhartr, enxc, envxc, 0.d0, .false., dtset%nspden)
!  WARNING: This is wrong, but... that's life.
   nfft_tot = product(dtset%ngfft(1:3))
 else
   hgrid(:) = 0.5d0 * wvl%h(:)

!  Data are always distributed when using the wavelets, even if nproc = 1.
!  The size given is the complete size of the box, not the distributed size
!  stored in ngfft.
   if (nproc > 1) then
     datacode = 'D'
   else
     datacode = 'G'
   end if

   call PSolver(bndcode, datacode, me, nproc, &
&   wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i, &
&   dtset%ixc, hgrid(1), hgrid(2), hgrid(3), &
&   rhor, kernel, vxc, enhartr, enxc, envxc, 0.d0, .false., dtset%nspden)
   nfft_tot = wvl%Glr%d%n1i * wvl%Glr%d%n2i * wvl%Glr%d%n3i
 end if

!PSolver work in place, we set back the rhor values.
 do i = 1, dtset%nfft, 1
   tmpPot     = rhor(i, 1)
!  Rhor total was saved in vhartr and current rhor(:,2) is down spin
   rhor(i, 1) = vhartr(i)
   if (dtset%nspden == 2) rhor(i, 2) = rhor(i, 1) - rhor(i, 2)
   vhartr(i)  = tmpPot
 end do

!!!  write(message, "(A,A,3F16.6)") "Psolver_rhohxc(): e_hartr, e_xc, e_vxc", ch10, &
!!!       & enhartr, enxc, envxc
!!!  call wrtout(std_out, message,'COLL')

!Compute vxcavg
 if (associated(mpi_enreg%nscatterarr)) then
   call mean_fftr(vxc, vxcmean, mpi_enreg, wvl%Glr%d%n1i * wvl%Glr%d%n2i * &
&   mpi_enreg%nscatterarr(me, 2), nfft_tot, dtset%nspden)
 else
   call mean_fftr(vxc, vxcmean, mpi_enreg, dtset%nfft, nfft_tot, dtset%nspden)
 end if
 vxcavg = vxcmean(1)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' Psolver_rhohxc: BUG -',ch10,&
& '  BigDFT is not compiled. Use --enable-bigdft during configure.'
 call wrtout(std_out, message, 'COLL')
 call leave_new('COLL')
#endif

end subroutine Psolver_rhohxc
!!***
