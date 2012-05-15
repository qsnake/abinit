!{\src2tex{textfont=tt}}
!!****f* ABINIT/Psolver_hartree
!! NAME
!! Psolver_hartree
!!
!! FUNCTION
!! Given rho(r), compute Hartree potential considering the system as
!! an isolated one. This potential is obtained from the convolution
!! of 1/r and rho(r), treated in Fourier space. This method is a wrapper around
!! Psolver() developped for BigDFT.
!! It does not compute the xc energy nor potential. See psolver_rhohxc() to do it.
!! WARNING : the XC energy and potential computation capability has been
!! for spin-polarized case, as everything is done as if nspden=1
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
!!  vhartr(nfft)=Hartree potential.
!!
!! NOTE
!!  In PSolver, with nspden == 2, rhor(:,1) = density up and
!!                                rhor(:,2) = density down.
!!  But in ABINIT (dtset%usewvl != 1) rhor(:,1) = total density and
!!                                    rhor(:,2) = density up .
!!  In ABINIT (dtset%usewvl != 1), the same convention is used as in PSolver.
!!
!! PARENTS
!!      mklocl_realspace,mklocl_wavelets,nres2vres
!!
!! CHILDREN
!!      leave_new,psolver,psolver_kernel,wrtout,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine Psolver_hartree(dtset, enhartr, mpi_enreg, rhor, rprimd, vhartr, wvl)

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
#define ABI_FUNC 'Psolver_hartree'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_62_poisson, except_this_one => Psolver_hartree
!End of the abilint section

  implicit none

  !Arguments ------------------------------------
  !scalars
  real(dp), intent(out)         :: enhartr
  type(dataset_type),intent(in) :: dtset
  type(MPI_type),intent(inout)  :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  !arrays
  real(dp),intent(in)    :: rprimd(3,3)
  real(dp),intent(in)    :: rhor(dtset%nfft,dtset%nspden)
  real(dp),intent(out)   :: vhartr(dtset%nfft)

  !Local variables-------------------------------
  !scalars
  integer  :: comm,me,nproc,n1,n2,n3
  real(dp) :: hgrid(3)
  real(dp) :: enxc, evxc
  character(len=500) :: message
  character(len = 1) :: datacode, bndcode
  !arrays
  real(dp), pointer :: kernel(:)
  real(dp), dimension(1) :: pot_ion_dummy

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

 if(dtset%nspden > 2 .and. dtset%usewvl/=0 )then
   write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&   ' PSolver_hartree: BUG -',ch10,&
&   '  Only non-spin-polarised or collinear spin is allowed for wavelets,',ch10,&
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
 write(message, "(A,A,A,3I6)") "Psolver_hartree(): compute potential (Vhartree)...", ch10, &
& " | dimension:", dtset%ngfft(1:3)
 call wrtout(std_out, message,'COLL')

#if defined HAVE_DFT_BIGDFT
 if (dtset%usewvl == 0) then
   vhartr(:)  = rhor(:, 1)

   hgrid(1) = rprimd(1, 1) / dtset%ngfft(1)
   hgrid(2) = rprimd(2, 2) / dtset%ngfft(2)
   hgrid(3) = rprimd(3, 3) / dtset%ngfft(3)

   datacode = 'G'
   n1=dtset%ngfft(1)
   n2=dtset%ngfft(2)
   n3=dtset%ngfft(3)
!  This may not work with MPI in the planewave code...
 else
   if(dtset%nspden==1)vhartr(:)  = rhor(:, 1)
   if(dtset%nspden==2)vhartr(:)  = rhor(:, 1) + rhor(:, 2)
   hgrid(:) = 0.5d0 * wvl%h(:)
!  The data are 'D'istributed in the wavelet case or 'G'lobal otherwise.
   if (nproc > 1) then
     datacode = 'D'
   else
     datacode = 'G'
   end if
   n1=wvl%Glr%d%n1i
   n2=wvl%Glr%d%n2i
   n3=wvl%Glr%d%n3i
 end if

!We attack PSolver with the total density contained in vhartr.
!This is also valid for spin-polarized (collinear and non-collinear)
!systems. Thus we enter nspden (last arg of PSolver) as being 1.
!Warning : enxc and evxc are meaningless.
 call PSolver(bndcode, datacode, me, nproc, n1, n2, n3,&
& 0, hgrid(1), hgrid(2), hgrid(3), vhartr, kernel, pot_ion_dummy, &
& enhartr, enxc, evxc, 0.d0, .false., 1)

#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' Psolver_hartree: BUG -',ch10,&
& '  BigDFT is not compiled. Use --enable-bigdft during configure.'
 call wrtout(std_out, message, 'COLL')
 call leave_new('COLL')
#endif

end subroutine Psolver_hartree
!!***
