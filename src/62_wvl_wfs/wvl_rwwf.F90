!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_read
!! NAME
!! wvl_read
!!
!! FUNCTION
!! Simple wrapper around the read disk methods of BigDFT for wavefunctions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  mpi_enreg=informations about MPI parallelization
!!  option= -2 for write with BigDFT format,
!!          -1 for reading wavelets coefficients with BigDFT format,
!!          2 for write,
!!          1 for read.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>=struct info for wavefunction
!!  wfs <type(wvl_wf_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      wvl_wfsinp_disk
!!
!! CHILDREN
!!      etsf_io_basisdata_put,etsf_io_electrons_put,etsf_io_low_error_to_str
!!      etsf_io_main_put,leave_new,writemywaves,writeonewave,wrtout,xcomm_world
!!      xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_abitypes
 use defs_wvltypes
 use m_wffile 

#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only: readonewave, reformatonewave
#endif
#if defined HAVE_TRIO_ETSF_IO
  use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_read'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
!End of the abilint section

  implicit none

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(hdr_type), intent(in)                :: hdr0
  type(hdr_type), intent(in)                :: hdr
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type),intent(in)              :: wff
  type(wvl_wf_type), intent(inout)          :: wfs
  type(wvl_internal_type), intent(in)       :: wvl
  !arrays
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(in)                      :: xred(3, dtset%natom)

!Local variables-------------------------------
  character(len = 500)  :: message
  logical               :: reformat
  integer               :: i, iBand, iCoeff, iCoarse, iFine, bandSize
  integer               :: comm,me
  real(dp), allocatable, target :: wvl_band(:)
  integer               :: coord(3), n_old(3)
  real(dp), allocatable :: xcart(:,:), psifscf(:,:,:), psigold(:,:,:,:,:,:)
  real(dp), allocatable :: xred_(:,:), xcart_old(:,:)
  real(dp)              :: hgrid_old(3), displ
#if defined HAVE_TRIO_ETSF_IO
  character(len=etsf_io_low_error_len)  :: errmess
  type(etsf_basisdata) :: basis_folder
  type(etsf_main)      :: main_folder
  type(etsf_electrons) :: electrons_folder
  logical              :: lstat
  type(etsf_io_low_error) :: error
#endif

! *********************************************************************

 if (abs(option) /= 1) then
   write(message,'(a,a,a,a,a,a,i5,a)') ch10,&
&   ' wvl_read: BUG -', ch10, &
&   '  Option argument is wrong,', ch10, &
&   '  awaited values are -1 or  1 but option = ', option, '.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

#if defined HAVE_DFT_BIGDFT
 call xcomm_world(mpi_enreg,comm,myrank=me)
!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 ABI_ALLOCATE(xcart_old,(3, dtset%natom))
 ABI_ALLOCATE(xred_,(3, dtset%natom))
 xred_ = xred
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred_)

 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_read:  read wavefunctions from file.'
 call wrtout(std_out,message,'COLL')

 if (option > 0) then
   bandSize = wfs%Glr%wfd%nvctr_c + 7 * wfs%Glr%wfd%nvctr_f
!  Read in the ABINIT way.
   if (wff%accesswff == IO_MODE_FORTRAN .or. (wff%accesswff == IO_MODE_FORTRAN_MASTER .and. wff%master==wff%me)) then
     ABI_ALLOCATE(psifscf,(wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i))
     do iBand = 1, dtset%mband * dtset%nsppol, 1
       call readonewave(wff%unwff, .false., iBand, me, &
&       wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3, &
&       wvl%h(1), wvl%h(2), wvl%h(3), wvl%atoms, &
&       wfs%Glr%wfd, xcart_old, xcart, &
&       wfs%psi(bandSize * (iBand - me * wfs%orbs%norbp - 1) + 1: &
&       bandSize * (iBand - me * wfs%orbs%norbp - 1) + bandSize), &
&       wfs%orbs%eval(iBand), psifscf)
     end do
     ABI_DEALLOCATE(psifscf)
#if defined HAVE_TRIO_ETSF_IO
   else if (wff%accesswff == IO_MODE_ETSF) then
!    Read a NetCDF file.
!    coordinates_of_grid_points is used to store the geometric
!    position of coefficients of wavelets i, as integer in
!    wvl%ni(:) dimensions.
!    coefficients_of_wavefunctions is used to store the psi values for
!    each wavelet.

!    We check if we need reformating
     if (abs(hdr%rprimd(1,1) / hdr%ngfft(1) - &
&     hdr0%rprimd(1,1) / hdr0%ngfft(1)) < tol8 .and. &
&     maxval(abs(hdr%nwvlarr - hdr0%nwvlarr)) == 0 .and. & 
&     maxval(abs(hdr%ngfft   - hdr0%ngfft  )) == 0 ) then
       reformat = .false.
       write(message, '(a,a,a,a)' ) ch10,&
&       ' wvl_read:  wavefunctions needs NO reformatting.'
       call wrtout(std_out,message,'COLL')           

     else
       reformat = .true.
       write(message, '(a,a,a,a)' ) ch10,&
&       ' wvl_read:  wavefunctions needs reformatting.'
       call wrtout(std_out,message,'COLL')

!      Values for old system
       xred_ = hdr0%xred
       call xredxcart(dtset%natom, 1, hdr0%rprimd, xcart_old, xred_)
       displ = zero
       do i = 1, dtset%natom, 1
         displ = (xcart_old(1, i) - xcart(1, i)) ** 2 + &
&         (xcart_old(2, i) - xcart(2, i)) ** 2 + &
&         (xcart_old(3, i) - xcart(3, i)) ** 2
       end do
       displ = sqrt(displ)

!      Do the reformating by expanding old wavefunctions on the grid
!      and interpolate it on the new grid using BigDFT reformatonewave().
       n_old = (hdr0%ngfft(1:3) - 31) / 2
       hgrid_old(1) = hdr0%rprimd(1,1) / n_old(1)
       hgrid_old(2) = hdr0%rprimd(2,2) / n_old(2)
       hgrid_old(3) = hdr0%rprimd(3,3) / n_old(3)
       ABI_ALLOCATE(psigold,(0:n_old(1), 2, 0:n_old(2), 2, 0:n_old(3), 2))
       ABI_ALLOCATE(psifscf,(wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i))
       ABI_ALLOCATE(basis_folder%coordinates_of_basis_grid_points%data2D,(3,hdr0%nwvlarr(1)))
     end if

!    We read the basis set definition
     ABI_ALLOCATE(basis_folder%number_of_coefficients_per_grid_point%data1D,(hdr0%nwvlarr(1)))
     call etsf_io_basisdata_get(wff%unwff, basis_folder, lstat, error)
     if (.not. lstat) then
       call etsf_io_low_error_to_str(errmess, error)
       write(message, "(A,A,A,A)") ch10, " wvl_read: ERROR -", ch10, &
&       errmess(1:min(475, len(errmess)))
       call wrtout(std_out, message, 'COLL')
     end if

!    We allocate a temporary array to read the wavefunctions
!    and reorder then into wfs%psi
     ABI_ALLOCATE(wvl_band,(sum(basis_folder%number_of_coefficients_per_grid_point%data1D)))
!    We just read the file from disk to memory, band per band.
     do iBand = 1, dtset%mband, 1
       main_folder%wfs_coeff__state_access = iBand
       main_folder%coefficients_of_wavefunctions%data1D => wvl_band
       call etsf_io_main_get(wff%unwff, main_folder, lstat, error)
       if (.not. lstat) then
         call etsf_io_low_error_to_str(errmess, error)
         write(message, "(A,A,A,A)") ch10, " wvl_read: ERROR -", ch10, &
&         errmess(1:min(475, len(errmess)))
         call wrtout(std_out, message, 'COLL')
       end if

!      Now we reorder
       if (reformat) then
!        Expand the wf on the old grid.
         psigold = zero
         iCoeff = 1
         do i = 1, hdr0%nwvlarr(1), 1
           coord = basis_folder%coordinates_of_basis_grid_points%data2D(:, i) / 2
           psigold(coord(1), 1, coord(2), 1, coord(3), 1) = wvl_band(iCoeff)
           iCoeff = iCoeff + 1
           if (basis_folder%number_of_coefficients_per_grid_point%data1D(i) == 8) then
             psigold(coord(1), 2, coord(2), 1, coord(3), 1) = wvl_band(iCoeff + 0)
             psigold(coord(1), 1, coord(2), 2, coord(3), 1) = wvl_band(iCoeff + 1)
             psigold(coord(1), 2, coord(2), 2, coord(3), 1) = wvl_band(iCoeff + 2)
             psigold(coord(1), 1, coord(2), 1, coord(3), 2) = wvl_band(iCoeff + 3)
             psigold(coord(1), 2, coord(2), 1, coord(3), 2) = wvl_band(iCoeff + 4)
             psigold(coord(1), 1, coord(2), 2, coord(3), 2) = wvl_band(iCoeff + 5)
             psigold(coord(1), 2, coord(2), 2, coord(3), 2) = wvl_band(iCoeff + 6)
             iCoeff = iCoeff + 7
           end if
         end do

         call reformatonewave(me, displ, wfs%Glr%wfd, wvl%atoms, hgrid_old(1), &
&         hgrid_old(2), hgrid_old(3), n_old(1), n_old(2), &
&         n_old(3), xcart_old, psigold, &
&         wvl%h(1), wvl%h(2), wvl%h(3), &
&         wvl%Glr%d%n1, wvl%Glr%d%n2, &
&         wvl%Glr%d%n3, xcart, psifscf, &
&         wfs%psi(bandSize * (iBand - me * wfs%orbs%norbp - 1) + 1: &
&         bandSize * (iBand - me * wfs%orbs%norbp - 1) + bandSize))
       else
!        No reformating, just reorder
         iCoeff = 1
         iFine = bandSize * (iBand - me * wfs%orbs%norbp - 1) + wfs%Glr%wfd%nvctr_c + 1
         do iCoarse = 1, hdr0%nwvlarr(1), 1
           wfs%psi(bandSize * (iBand - me * wfs%orbs%norbp - 1) + iCoarse) = &
&           wvl_band(iCoeff)
           iCoeff  = iCoeff  + 1
           if (basis_folder%number_of_coefficients_per_grid_point%data1D(iCoarse) == 8) then
             wfs%psi(iFine:iFine + 6) = wvl_band(iCoeff:iCoeff + 6)
             iFine  = iFine  + 7
             iCoeff = iCoeff + 7
           end if
         end do
       end if
     end do

!    We read wfs%eval (TODO maybe removed later).
     electrons_folder%eigenvalues%data1D => wfs%orbs%eval
     call etsf_io_electrons_get(wff%unwff, electrons_folder, lstat, error)
     if (.not. lstat) then
       call etsf_io_low_error_to_str(errmess, error)
       write(message, "(A,A,A,A)") ch10, " wvl_read: ERROR -", ch10, &
&       errmess(1:min(475, len(errmess)))
       call wrtout(std_out, message, 'COLL')
     end if
     
!    Deallocate temporary arrays
     ABI_DEALLOCATE(wvl_band)
     ABI_DEALLOCATE(basis_folder%number_of_coefficients_per_grid_point%data1D)
     if (reformat) then
       ABI_DEALLOCATE(basis_folder%coordinates_of_basis_grid_points%data2D)
       ABI_DEALLOCATE(psigold)
       ABI_DEALLOCATE(psifscf)
     end if
#endif

   else
     write(message,'(a,a,a,a,a,a,i5,a)') ch10,&
&     ' wvl_read: BUG -', ch10, &
&     '  wff%accesswff argument is wrong,', ch10, &
&     '  awaited values are -1, 0 (or 3 if netcdf/etsf_io is available) but value = ', wff%accesswff, '.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
 else
   call readmywaves(me, "wavefunctions", wfs%orbs, &
&   wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3, &
&   wvl%h(1), wvl%h(2), wvl%h(3), wvl%atoms, &
&   xcart_old, xcart, wfs%Glr%wfd, wfs%psi)
 end if

 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xcart_old)
 ABI_DEALLOCATE(xred_)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_outwf : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif

end subroutine wvl_read
!!***


!!****f* ABINIT/wvl_write
!! NAME
!! wvl_write
!!
!! FUNCTION
!! Simple wrapper around the write disk methods of BigDFT for wavefunctions.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  mpi_enreg=informations about MPI parallelization
!!  option= -2 for write with BigDFT format,
!!          -1 for reading wavelets coefficients with BigDFT format,
!!          2 for write,
!!          1 for read.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>=struct info for wavefunction
!!  wfs <type(wvl_wf_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      outwf
!!
!! CHILDREN
!!      etsf_io_basisdata_put,etsf_io_electrons_put,etsf_io_low_error_to_str
!!      etsf_io_main_put,leave_new,writemywaves,writeonewave,wrtout,xcomm_world
!!      xredxcart
!!
!! SOURCE
subroutine wvl_write(dtset, eigen, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)

 use m_profiling

  use defs_basis
  use defs_datatypes
 use defs_abitypes
  use defs_wvltypes
 use m_wffile
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : writeonewave
#endif
#if defined HAVE_TRIO_ETSF_IO
  use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_write'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
!End of the abilint section

  implicit none

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type),intent(in)              :: wff
  type(wvl_wf_type), intent(in)             :: wfs
  type(wvl_internal_type), intent(in)       :: wvl
  !arrays
  real(dp), intent(in), target              :: eigen(dtset%mband)
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(inout)                   :: xred(3, dtset%natom)

!Local variables-------------------------------
  character(len = 500)  :: message
  integer               :: comm,me
  integer               :: iGrid, iCoeff, iCoarse, iFine
  integer               :: iorb, jj, j0, j1, i, ii, i0, i1, i2, i3
  integer               :: iseg, nseg, ipsi, npsi
  integer               :: iband
  real(dp), allocatable, target :: wvl_band(:)
  real(dp), allocatable :: xcart(:,:)
#if defined HAVE_TRIO_ETSF_IO
  integer, allocatable  :: coeff_map(:,:,:)
  integer, allocatable, target :: wvl_coord(:,:), wvl_ncoeff(:)
  character(len=etsf_io_low_error_len)  :: errmess
  type(etsf_basisdata) :: basis_folder
  type(etsf_main)      :: main_folder
  type(etsf_electrons) :: electrons_folder
  logical              :: lstat
  type(etsf_io_low_error) :: error
#endif

! *********************************************************************

 if (abs(option) /= 2) then
   write(message,'(a,a,a,a,a,a,i5,a)') ch10,&
&   ' wvl_write: BUG -', ch10, &
&   '  Option argument is wrong,', ch10, &
&   '  awaited values are -2 or  2 but option = ', option, '.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

#if defined HAVE_DFT_BIGDFT
 call xcomm_world(mpi_enreg,comm,myrank=me)
!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_write:  Write wavefunctions to file.'
 call wrtout(std_out,message,'COLL')

 if (option > 0) then
!  Write in the ABINIT way.
   if (wff%accesswff == IO_MODE_FORTRAN .or. (wff%accesswff == IO_MODE_FORTRAN_MASTER .and. wff%master==wff%me)) then
     iseg = wfs%Glr%wfd%nseg_c
     nseg = wfs%Glr%wfd%nseg_c + wfs%Glr%wfd%nseg_f
     ipsi = wfs%Glr%wfd%nvctr_c
     npsi = wfs%Glr%wfd%nvctr_c + 7 * wfs%Glr%wfd%nvctr_f
     do iorb = 1, dtset%mband
       call writeonewave(wff%unwff, .false., iorb, wvl%Glr%d%n1, &
&       wvl%Glr%d%n2, wvl%Glr%d%n3, &
&       wvl%h(1), wvl%h(2), wvl%h(3), dtset%natom, &
&       xcart, wfs%Glr%wfd%nseg_c, wfs%Glr%wfd%nvctr_c, &
&       wfs%Glr%wfd%keyg(:,1:iseg), wfs%Glr%wfd%keyv(1:iseg), wfs%Glr%wfd%nseg_f, &
&       wfs%Glr%wfd%nvctr_f, wfs%Glr%wfd%keyg(:, iseg + 1:nseg), &
&       wfs%Glr%wfd%keyv(iseg + 1:nseg), &
&       wfs%psi(npsi * (iorb - me * wfs%orbs%norbp - 1) + 1: &
&       npsi * (iorb - me * wfs%orbs%norbp - 1) + ipsi), &
&       wfs%psi(npsi * (iorb - me * wfs%orbs%norbp - 1) + ipsi + 1: &
&       npsi * (iorb - me * wfs%orbs%norbp - 1) + npsi), &
&       wfs%orbs%eval(iorb))
     end do
#if defined HAVE_TRIO_ETSF_IO
   else if (wff%accesswff == IO_MODE_ETSF) then
!    Write a NetCDF file.
!    coordinates_of_grid_points is used to store the geometric
!    position of coefficients of wavelets i, as integer in
!    wvl%dpSize(:) dimensions.
!    coefficients_of_wavefunctions is used to store the psi values for
!    each wavelet.

!    We write the basis set.
!    =======================
     ABI_ALLOCATE(wvl_coord,(3, wfs%Glr%wfd%nvctr_c))
     ABI_ALLOCATE(wvl_ncoeff,(wfs%Glr%wfd%nvctr_c))
!    Will store the grid index for a given geometric point
     ABI_ALLOCATE(coeff_map,(0:wvl%Glr%d%n1,0:wvl%Glr%d%n2,0:wvl%Glr%d%n3))
!    coarse part
     iGrid = 0
     coeff_map = 0
     do iseg = 1, wfs%Glr%wfd%nseg_c
       jj = wfs%Glr%wfd%keyv(iseg)
       j0 = wfs%Glr%wfd%keyg(1, iseg)
       j1 = wfs%Glr%wfd%keyg(2, iseg)
       ii = j0 - 1
       i3 = ii / ((wvl%Glr%d%n1 + 1) * &
&       (wvl%Glr%d%n2 + 1))
       ii = ii - i3 * (wvl%Glr%d%n1 + 1) * &
&       (wvl%Glr%d%n2 + 1)
       i2 = ii / (wvl%Glr%d%n1 + 1)
       i0 = ii - i2 * (wvl%Glr%d%n1 + 1)
       i1 = i0 + j1 - j0
       do i = i0, i1
         iGrid = iGrid + 1
         coeff_map(i, i2, i3) = iGrid
         wvl_coord(:, iGrid) = (/ i * 2, i2 * 2, i3 * 2 /)
         wvl_ncoeff(iGrid) = 1
       end do
     end do
!    fine part
     do iseg = 1, wfs%Glr%wfd%nseg_f
       jj = wfs%Glr%wfd%keyv(wfs%Glr%wfd%nseg_c + iseg)
       j0 = wfs%Glr%wfd%keyg(1, wfs%Glr%wfd%nseg_c + iseg)
       j1 = wfs%Glr%wfd%keyg(2, wfs%Glr%wfd%nseg_c + iseg)
       ii = j0 - 1
       i3 = ii / ((wvl%Glr%d%n1 + 1) * &
&       (wvl%Glr%d%n2 + 1))
       ii = ii - i3 * (wvl%Glr%d%n1 + 1) * &
&       (wvl%Glr%d%n2 + 1)
       i2 = ii / (wvl%Glr%d%n1 + 1)
       i0 = ii - i2 * (wvl%Glr%d%n1 + 1)
       i1 = i0 + j1 - j0
       do i = i0, i1
         iGrid = coeff_map(i , i2 , i3)
         wvl_ncoeff(iGrid) = wvl_ncoeff(iGrid) + 7
       end do
     end do
     basis_folder%coordinates_of_basis_grid_points%data2D => wvl_coord
     basis_folder%number_of_coefficients_per_grid_point%data1D => wvl_ncoeff
     call etsf_io_basisdata_put(wff%unwff, basis_folder, lstat, error)
     if (.not. lstat) then
       call etsf_io_low_error_to_str(errmess, error)
       write(message, "(A,A,A,A)") ch10, " wvl_write: ERROR -", ch10, &
&       errmess(1:min(475, len(errmess)))
       call wrtout(std_out, message, 'COLL')
     end if

!    We write the wavefuctions.
!    ==========================
!    We write the wavefunctions band per band.
     ABI_ALLOCATE(wvl_band,(wfs%Glr%wfd%nvctr_c + 7 * wfs%Glr%wfd%nvctr_f))
!    We reorder the wfs%psi(:,iband) into wvl_band
     do iband = 1, dtset%mband, 1
!      iCoeff is the index of the coefficient we are writing in wvl_band
       iCoeff  = 1
!      iCoarse is the index of the coarse part in wfs%psi (=iGrid)
       iCoarse = 1
!      iFine is the index of the fine part in wfs%psi
       iFine   = wfs%Glr%wfd%nvctr_c + 1
       npsi    = wfs%Glr%wfd%nvctr_c + 7 * wfs%Glr%wfd%nvctr_f
       do iGrid = 1, wfs%Glr%wfd%nvctr_c, 1
         wvl_band(iCoeff) = &
&         wfs%psi((iband - me * wfs%orbs%norbp - 1) * npsi + iCoarse)
         iCoeff  = iCoeff + 1
         iCoarse = iCoarse + 1
         if (wvl_ncoeff(iGrid) == 8) then
           wvl_band(iCoeff:iCoeff + 6) = &
&           wfs%psi((iband - me * wfs%orbs%norbp - 1) * npsi + iFine: &
&           (iband - me * wfs%orbs%norbp - 1) * npsi + iFine + 6)
           iCoeff = iCoeff + 7
           iFine  = iFine + 7
         end if
       end do
       main_folder%wfs_coeff__state_access = iband
       main_folder%coefficients_of_wavefunctions%data1D => wvl_band
       call etsf_io_main_put(wff%unwff, main_folder, lstat, error)
       if (.not. lstat) then
         call etsf_io_low_error_to_str(errmess, error)
         write(message, "(A,A,A,A)") ch10, " wvl_write: ERROR -", ch10, &
&         errmess(1:min(475, len(errmess)))
         call wrtout(std_out, message, 'COLL')
       end if
     end do

!    We write the electronic informations.
!    =====================================
     electrons_folder%eigenvalues%data1D => eigen
     electrons_folder%eigenvalues__number_of_states = dtset%mband
     call etsf_io_electrons_put(wff%unwff, electrons_folder, lstat, error)
     if (.not. lstat) then
       call etsf_io_low_error_to_str(errmess, error)
       write(message, "(A,A,A,A)") ch10, " wvl_write: ERROR -", ch10, &
&       errmess(1:min(475, len(errmess)))
       call wrtout(std_out, message, 'COLL')
     end if

!    We deallocate all arrays
     ABI_DEALLOCATE(wvl_band)
     ABI_DEALLOCATE(wvl_coord)
     ABI_DEALLOCATE(wvl_ncoeff)
     ABI_DEALLOCATE(coeff_map)
#endif
   else
     write(message,'(a,a,a,a,a,a,i5,a)') ch10,&
&     ' wvl_write: BUG -', ch10, &
&     '  wff%accesswff argument is wrong,', ch10, &
&     '  awaited values are -1, 0 (or 3 if netcdf/etsf_io is available) but value = ', wff%accesswff, '.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
 else
!  Write in the BigDFT way.
   call  writemywaves(me, "wavefunctions", wfs%orbs, &
&   wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3, &
&   wvl%h(1), wvl%h(2), wvl%h(3),wvl%atoms, &
&   xcart, wfs%Glr%wfd, wfs%psi)
 end if

 ABI_DEALLOCATE(xcart)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_outwf : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif

end subroutine wvl_write
!!***
