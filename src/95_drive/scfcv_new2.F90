!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv_new2
!! NAME
!! scfcv_new2
!!
!! FUNCTION
!! WARNING : Temporary wrapper to scfcv
!! Self-consistent-field convergence.
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! ground state and optionally to compute forces and energy.
!! This routine is called to compute forces for given atomic
!! positions or else to do non-SCF band structures.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, GMR, AR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!!  electronpositron <type(electronpositron_type)>=quantities for the
!!      electron-positron annihilation
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wffnew,wffnow=struct info for wf disk files.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic
!!     coordinates at output, current xred is transferred to xred_old
!!
!! NOTES
!! It is worth to explain THE USE OF FFT GRIDS:
!! ============================================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      gstate,mover
!!
!! CHILDREN
!!      scfcv_new
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


! ! THIS IS THE OBJECTIVE
! subroutine scfcv_new(ab_scfcv,dtset,electronpositron,paw_dmft,rprimd,wffnew,&
! &  wffnow,xred,xred_old)
! ! THIS IS THE OBJECTIVE

subroutine scfcv_new2(scfcv_args,electronpositron,&
&  rhog,rhor,&
&  rprimd,&
&  xred,xred_old)

 use m_profiling

 use defs_basis
 use defs_datatypes
! use m_wffile
 use defs_abitypes
! use defs_scftypes
! use defs_wvltypes
! use defs_parameters
! use defs_rectypes
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 !use m_rec
 !use m_io_tools, only : flush_unit
 !use m_paw_dmft, only: paw_dmft_type

 use defs_scfcvargs

#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfcv_new2'
 use interfaces_95_drive, except_this_one => scfcv_new2
!End of the abilint section

 implicit none


!Arguments ------------------------------------
!scalars
 type(ab_scfcvargs),intent(inout) :: scfcv_args

! type(ab_scfcv_args_in),intent(in) :: scfcv_args%ab_scfcv_in
! type(ab_scfcv_args_inout),intent(inout) :: ab_scfcv_inout
! type(dataset_type),intent(inout) :: dtset
 type(electronpositron_type),pointer :: electronpositron
! type(paw_dmft_type), intent(inout) :: paw_dmft
! type(wffile_type),intent(inout) :: wffnew,wffnow
!arrays
 real(dp), intent(inout) :: rprimd(3,3)
 real(dp), intent(inout) :: xred(3,scfcv_args%dtset%natom),xred_old(3,scfcv_args%dtset%natom)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 !real(dp), intent(out) :: resid(scfcv_args%dtset%mband*scfcv_args%dtset%nkpt*scfcv_args%dtset%nsppol)

!Local variables -------------------------
!scalars
! integer :: ii
! logical :: DEBUG=.FALSE.

! *********************************************************************

 call scfcv_new(scfcv_args%ab_scfcv_in,scfcv_args%ab_scfcv_inout,&
& scfcv_args%dtset,electronpositron,&
& scfcv_args%paw_dmft,&
& rhog,rhor,&
& rprimd,&
& scfcv_args%wffnew,scfcv_args%wffnow,&
& xred,xred_old)


!!!!!!!!!!!!!!!!!!!!
!!! COMMENTED
!!!!!!!!!!!!!!!!!!!!


!if (DEBUG) then
!write(std_out,*) 'INTENT(IN) ARGUMENTS ON SCFCV'
!write(std_out,*) 'atindx=',scfcv_args%ab_scfcv_in%atindx
!write(std_out,*) 'atindx1=',scfcv_args%ab_scfcv_in%atindx1
!write(std_out,*) 'cpus=',scfcv_args%ab_scfcv_in%cpus
!write(std_out,*) 'ecore=',scfcv_args%ab_scfcv_in%ecore
!write(std_out,*) 'fatvshift=',scfcv_args%ab_scfcv_in%fatvshift
!write(std_out,*) 'iapp=',scfcv_args%ab_scfcv_in%iapp
!write(std_out,*) 'indsym=',scfcv_args%ab_scfcv_in%indsym
!write(std_out,*) 'kg=',scfcv_args%ab_scfcv_in%kg
!write(std_out,*) 'nattyp=',scfcv_args%ab_scfcv_in%nattyp
!write(std_out,*) 'ndtpawuj=',scfcv_args%ab_scfcv_in%ndtpawuj
!write(std_out,*) 'npwarr=',scfcv_args%ab_scfcv_in%npwarr
!!  write(std_out,*) 'pawang=',scfcv_args%ab_scfcv_in%pawang
!!  write(std_out,*) 'pawrad=',scfcv_args%ab_scfcv_in%pawrad
!!  write(std_out,*) 'pawtab=',scfcv_args%ab_scfcv_in%pawtab
!write(std_out,*) 'phnons=',scfcv_args%ab_scfcv_in%phnons
!!  write(std_out,*) 'psps=',scfcv_args%ab_scfcv_in%psps
!write(std_out,*) 'pwind=',scfcv_args%ab_scfcv_in%pwind
!write(std_out,*) 'pwind_alloc=',scfcv_args%ab_scfcv_in%pwind_alloc
!write(std_out,*) 'pwnsfac=',scfcv_args%ab_scfcv_in%pwnsfac
!write(std_out,*) 'ylm=',scfcv_args%ab_scfcv_in%ylm
!write(std_out,*) 'ylmgr=',scfcv_args%ab_scfcv_in%ylmgr
!end if

!!WVL - reformat the wavefunctions in the case of xred != xred_old
!if (scfcv_args%dtset%usewvl == 1 .and. maxval(xred_old - xred) > zero) then
!!  WVL - Before running scfcv, on non-first geometry step iterations,
!!  we need to reformat the wavefunctions, taking into acount the new
!!  coordinates.
!!  We prepare to change rhog (to be removed) and rhor.
!deallocate(rhog)
!deallocate(rhor)

!call wvl_wfsinp_reformat(scfcv_args%dtset,  scfcv_args%ab_scfcv_inout%mpi_enreg,&
!&   scfcv_args%ab_scfcv_in%psps, rprimd, scfcv_args%ab_scfcv_inout%wvl, xred, xred_old)
!scfcv_args%ab_scfcv_inout%nfftf = scfcv_args%dtset%nfft

!allocate(rhog(2, scfcv_args%dtset%nfft))
!allocate(rhor(2, scfcv_args%dtset%nfft))
!call wvl_mkrho(scfcv_args%dtset, scfcv_args%ab_scfcv_inout%mpi_enreg,&
!&   rhor,scfcv_args%ab_scfcv_inout%wvl%wfs)
!end if

!!Prepare the names of the auxiliary files whose name depend on the
!!itimimage, iimage and itime loops.
!call dtfil_init2(scfcv_args%ab_scfcv_inout%dtfil,scfcv_args%ab_scfcv_in%iapp,scfcv_args%ab_scfcv_inout%mpi_enreg)

!call scfcv(scfcv_args%ab_scfcv_in%atindx,&
!& scfcv_args%ab_scfcv_in%atindx1,&
!& scfcv_args%ab_scfcv_inout%cg,&
!& scfcv_args%ab_scfcv_in%cpus,&
!& scfcv_args%ab_scfcv_inout%dtefield,&
!& scfcv_args%ab_scfcv_inout%dtfil,&
!& scfcv_args%ab_scfcv_inout%dtpawuj,&
!& scfcv_args%dtset,&                     ! scfcv_args%dtset IS NOT in scfcv_args%ab_scfcv_inout
!& scfcv_args%ab_scfcv_in%ecore,&
!& scfcv_args%ab_scfcv_inout%eigen,&
!& scfcv_args%electronpositron,&          ! electronpositron IS NOT in ab_scfcv_inout
!& scfcv_args%ab_scfcv_in%fatvshift,&
!& scfcv_args%ab_scfcv_inout%hdr,&
!& scfcv_args%ab_scfcv_in%iapp,&
!& scfcv_args%ab_scfcv_in%indsym,&
!& scfcv_args%ab_scfcv_inout%initialized,&
!& scfcv_args%ab_scfcv_inout%irrzon,&
!& scfcv_args%ab_scfcv_in%kg,&
!& scfcv_args%ab_scfcv_in%cg,&
!& scfcv_args%ab_scfcv_inout%mpi_enreg,&
!& scfcv_args%ab_scfcv_in%nattyp,&
!& scfcv_args%ab_scfcv_in%ndtpawuj,&
!& scfcv_args%ab_scfcv_inout%nfftf,&
!& scfcv_args%ab_scfcv_in%npwarr,&
!& scfcv_args%ab_scfcv_inout%occ,&
!& scfcv_args%paw_dmft,&                  ! paw_dmft IS NOT in scfcv_args%ab_scfcv_inout
!& scfcv_args%ab_scfcv_in%pawang,&
!& scfcv_args%ab_scfcv_inout%pawfgr,&
!& scfcv_args%ab_scfcv_in%pawrad,&
!& scfcv_args%ab_scfcv_inout%pawrhoij,&
!& scfcv_args%ab_scfcv_in%pawtab,&
!& scfcv_args%ab_scfcv_in%phnons,&
!& scfcv_args%ab_scfcv_in%psps,&
!& scfcv_args%ab_scfcv_in%pwind,&
!& scfcv_args%ab_scfcv_in%pwind_alloc,&
!& scfcv_args%ab_scfcv_in%pwnsfac,&
!& scfcv_args%ab_scfcv_inout%rec_set,&
!& resid,&
!& scfcv_args%ab_scfcv_inout%results_gs,&
!& rhog,&
!& rhor,&
!& rprimd,&
!& scfcv_args%ab_scfcv_inout%scf_history,&
!& scfcv_args%ab_scfcv_inout%symrec,&
!& scfcv_args%ab_scfcv_inout%taug,&
!& scfcv_args%ab_scfcv_inout%taur,&
!& scfcv_args%wffnew,&
!& scfcv_args%wffnow,&
!& scfcv_args%ab_scfcv_inout%wvl,&
!& xred,&
!& xred_old,&
!& scfcv_args%ab_scfcv_in%ylm,&
!& scfcv_args%ab_scfcv_in%ylmgr)

!if(DEBUG)then
!write(std_out,*) 'rprimd='
!do ii=1,3
!write(std_out,*) rprimd(:,ii)
!end do
!write(std_out,*) 'xred='
!do ii=1,scfcv_args%dtset%natom
!write(std_out,*) xred(:,ii)
!end do
!end if

!!!!!!!!!!!!!!!!!!!!
!!! COMMENTED
!!!!!!!!!!!!!!!!!!!!

end subroutine scfcv_new2
!!***
