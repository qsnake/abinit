!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_io_screening
!! NAME
!!  m_io_screening
!!
!! FUNCTION
!!  This module contains the definition of the header of the 
!!  _SCR and _SUSC file as well as methods used to read/write/echo.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG)
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

MODULE m_io_screening

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
#if defined HAVE_MPI2
 use mpi
#endif
 use m_xmpi
 use m_errors

 use m_gwdefs,          only : epsilonm1_parameters
 use m_copy,            only : deep_copy
 use m_io_tools,        only : get_unit
 use m_numeric_tools,   only : print_arr, remove_copies, imax_loc
 use m_header,          only : hdr_copy, hdr_clean, hdr_nullify, hdr_mpio_skip
 use m_bz_mesh,         only : isequalk

 implicit none

 private
!!***

#if defined HAVE_MPI1
 include 'mpif.h'
#endif


!!****t* m_io_screening/ScrHdr_type
!! NAME
!!  ScrHdr_type
!!
!! FUNCTION
!!  The structure defining the header of the _SCR or _SUSC file.
!!
!! SOURCE

 type,public :: ScrHdr_type

  !Other variables that can be added are, for the moment, commented out. 
  !Most of them are related to the Abinit implementation  and should not compare in the ETSF specs

  integer :: ID           ! Matrix identifier: O if not yet defined, 1 for chi0, 
                          ! 2 for chi, 3 for epsilon, 4 for espilon^{-1}
  integer :: ikxc         ! Kxc kernel used, 0 for None (RPA), >0 for static TDDFT (=ixc), <0 for frequency-dependent TDDFT 
  integer :: inclvkb      ! q-->0 treatment, 0 for None, 1-2 for transversal gauge, 3 for longitudinal
  integer :: headform     ! format of the SCR header
  integer :: fform        ! File format:
  integer :: gwcalctyp    ! Calculation type (G0W0, G0W, GW ...)
  integer :: nI,nJ        ! Number of spin components (rows,columns) in chi|eps^-1. (1,1) if collinear. 
                          !  The internal representation of the matrix is eps(nI*npwe,nJ*npwe) 
  integer :: nqibz        ! Number of q-points in the IBZ.
  integer :: nqlwl        ! Number of points for the treatment of the long wavelength limit.
  integer :: nomega       ! Total number of frequencies.
! New variables related to the pole fits
  integer :: npoles       ! Number of poles for a pole-fitted screening
  integer :: ncoeff       ! Number of coefficients (total) = npoles*3 + 1 for phase
! end of new variables
  integer :: nbnds_used   ! Number of bands used during the screening calculation (only for info)
  integer :: npwe         ! Number of G vectors reported on the file.
  integer :: npwwfn_used  ! Number of G vectors for wavefunctions used during the screening calculation (only for info)
  integer :: spmeth       ! Method used to approximate the delta function in the expression for Im Chi_0
  integer :: test_type    ! 1 for TEST-PARTICLE, 2 for TEST-ELECTRON.
  integer :: tordering    ! 0 if not defined, 1 for Time-Ordered, 2 for Advanced, 3 for Retarded.

  real(dp) :: soenergy    ! Scissor Energy, zero if not used
  real(dp) :: spsmear     ! Smearing of the delta in case of spmeth==2
  real(dp) :: zcut        ! Imaginary shift to avoid the poles along the real axis.

!arrays
  character(len=80) :: title(2)
  ! Title describing the content of the file.

  integer,pointer  :: gvec(:,:)                 
  ! gvec(3,npwe) 
  ! G vectors in r.l.u.

  real(dp),pointer :: qibz(:,:)
  ! qibz(3,nqibz)
  ! q-points in r.l.u.

  real(dp),pointer :: qlwl(:,:)                 
  ! qlwl(3,nqlwl)
  ! q-points for the long wave-length limit treatment (r.l.u)

  complex(dpc),pointer :: lwing(:,:,:)         
  ! lwing(npwe,nomega,nqlwl)
  ! Lower wings for the different q"s -->0 

  complex(dpc),pointer :: omega(:)             
  ! omega(nomega) 
  ! All frequencies calculated both along the real and the imaginary axis.

  complex(dpc),pointer :: uwing(:,:,:)   
  ! uwing(npwe,nomega,nqlwl)
  ! Upper wings for the different q"s -->0 

  type(Hdr_type) :: Hdr   ! The abinit header.

 end type ScrHdr_type
!!***

 public :: scr_hdr_io        ! I/O of the header (read/write/echo).
 public :: print_ScrHdr      ! Print the SCR related part of the header.
 public :: init_ScrHdr       ! Initialize the header.
 public :: nullify_HScr      ! Set all pointers to NULL.
 public :: scrhdr_comm       ! Transmit the header.
 public :: free_scrhdr       ! Free the header.
 public :: copy_scrhdr       ! Deep copy of the SCR|SUSC header.
 public :: merge_ScrHdr      ! Merge two or more SCR headers.
 public :: write_screening   ! Write a q-slice of the matrix in G-space.
 public :: write_pole_screening ! Write a q-slice of the pole coefficients in G-space.
 public :: read_screening    ! Read the content of the (SCR|SUSC) file placed after the header. 
 public :: read_pole_screening ! Read the content of a pole-fit SCR file 
 
! TODO change  the prefix DP => DIEL
 integer,public,parameter :: HSCR_LATEST_HEADFORM=57 ! This is still the old one used in abinit


CONTAINS  !=========================================================================================================================
!!***

!!****f* m_io_screening/scr_hdr_io
!! NAME
!!  scr_hdr_io
!!
!! FUNCTION
!! This subroutine deals with the I/O of the ScrHdr_type structured variables (read/write/echo).
!! According to the value of rdwr, it reads the header of a file, writes it, or echo the value 
!! of the structured variable to a file. Note that, when reading, different records of ScrHdr
!! are allocated here, according to the values of the read variables. Records of ScrHdr should be 
!! deallocated correctly by a call to hdr_clean when ScrHdr is not used anymore.
!!
!! INPUTS
!!  accesswff=Option defining the file format of the external file.
!!  comm=MPI communicator.
!!  master=rank of the master node in comm, usually 0
!!  rdwr= if 1, read the ScrHdr structured variable from the header of the file,
!!        if 2, write the header to unformatted file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the ScrHdr without rewinding (unformatted)
!!        if 6, write the ScrHdr without rewinding (unformatted)
!!  unt=unit number of the file (unformatted if rdwr=1, 2, 5 or 6 formatted if rdwr=3,4)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  ScrHdr <type(ScrHdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!! In all cases, the file is supposed to be open already
!! When reading (rdwr=1) or writing (rdwr=2), rewind the file
!! When echoing (rdwr=3) does not rewind the file.
!! When reading (rdwr=5) or writing (rdwr=6), DOES NOT rewind the file
!!
!! In writing mode, the routine is supposed to called by the master node.
!! no check is done, it is up to the developer.
!!
!! PARENTS
!!      m_io_screening,m_screen,m_screening,mrgscr,screening,setup_bse
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine scr_hdr_io(fform,rdwr,unt,comm,master,accesswff,HScr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scr_hdr_io'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr,unt,accesswff,comm,master
 type(ScrHdr_type),intent(inout) :: Hscr

!Local variables-------------------------------
!scalars
 integer :: my_rank,nprocs,prtvol
 logical :: ltest,I_read
 character(len=500) :: msg
!arrays

! *************************************************************************

 DBG_ENTER("COLL")
 !@ScrHdr_type  

 !ltest=( ANY(accesswff==(/IO_MODE_FORTRAN,IO_MODE_MPI/) )
 ltest=( ANY(accesswff == (/IO_MODE_FORTRAN/)) )
 ABI_CHECK(ltest,'Wrong value for accesswff')

 ABI_UNUSED(master) ! FIXME

 ! === Initialize MPI info for comm ===
 my_rank = xcomm_rank(comm)
 nprocs  = xcomm_size(comm)

 if (rdwr==1.or.rdwr==5) then
   !
   ! === Reading the header of an unformatted file ===
   I_read=.TRUE. !(localrdwf==1.or.my_rank==master)
   if (I_read) then
     ! * Read the abinit header, rewinding of the file (if any) is done here.
     ! TODO write a wrapper, using accesswff as input
     if (accesswff==IO_MODE_FORTRAN) then
       call hdr_io(fform,Hscr%Hdr,rdwr,unt)
     else if (accesswff==IO_MODE_ETSF) then
       call hdr_io_etsf(fform,Hscr%Hdr,rdwr,unt)
     end if

     ! * Reset the variables absent in old versions.
     Hscr%fform=fform
     call reset_HSCR_(Hscr)

     select case (fform)
     case (1002)
       ! Old epsilon^-1 file used in version <5.6
       read(unt,ERR=10)Hscr%title
       read(unt,ERR=10)Hscr%npwe,Hscr%npwwfn_used,Hscr%nbnds_used,Hscr%nqibz,Hscr%nomega

       ABI_ALLOCATE(Hscr%gvec,(3,Hscr%npwe))
       read(unt,ERR=10)Hscr%gvec(:,:)

       ABI_ALLOCATE(Hscr%qibz,(3,Hscr%nqibz))
       read(unt,ERR=10)Hscr%qibz(:,:)

       ABI_ALLOCATE(Hscr%omega,(Hscr%nomega))
       read(unt,ERR=10)Hscr%omega(:)

       ! Quantities not present in the old file format, nqlwl is set to 0
       ABI_ALLOCATE(Hscr%qlwl,(3,Hscr%nqlwl))
       ABI_ALLOCATE(Hscr%lwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))
       ABI_ALLOCATE(Hscr%uwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))

     case (1102)
       ! File format for epsilon^-1, espilon, chi0 
       read(unt,ERR=10)Hscr%title
       read(unt,ERR=10)Hscr%ID,Hscr%nI,Hscr%nJ,Hscr%tordering,Hscr%test_type
       read(unt,ERR=10)Hscr%npwe,Hscr%npwwfn_used,Hscr%nbnds_used,Hscr%nqibz,Hscr%nomega,Hscr%nqlwl

       ABI_ALLOCATE(Hscr%gvec,(3,Hscr%npwe))
       read(unt,ERR=10)Hscr%gvec(:,:)

       ABI_ALLOCATE(Hscr%qibz,(3,Hscr%nqibz))
       read(unt,ERR=10)Hscr%qibz(:,:)

       ABI_ALLOCATE(Hscr%omega,(Hscr%nomega))
       read(unt,ERR=10)Hscr%omega(:)

       ! Read data for q-->0 limit.
       ABI_ALLOCATE(Hscr%qlwl,(3,Hscr%nqlwl))
       ABI_ALLOCATE(Hscr%lwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))
       ABI_ALLOCATE(Hscr%uwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))

       if (Hscr%nqlwl>0) then
         read(unt,ERR=10)Hscr%qlwl(:,:)
         read(unt,ERR=10)Hscr%lwing(:,:,:)
         read(unt,ERR=10)Hscr%uwing(:,:,:)
       end if

     case (2002)
       ! Pole-fit screening
       read(unt,ERR=10)Hscr%title
       read(unt,ERR=10)Hscr%npwe,Hscr%npwwfn_used,Hscr%nbnds_used,Hscr%nqibz,&
&       Hscr%nomega,Hscr%npoles,Hscr%ncoeff

       ABI_ALLOCATE(Hscr%gvec,(3,Hscr%npwe))
       read(unt,ERR=10)Hscr%gvec(:,:)

       ABI_ALLOCATE(Hscr%qibz,(3,Hscr%nqibz))
       read(unt,ERR=10)Hscr%qibz(:,:)

       ABI_ALLOCATE(Hscr%omega,(Hscr%nomega))
       read(unt,ERR=10)Hscr%omega(:)

       ! Quantities not present in the old file format, nqlwl is set to 0
       ABI_ALLOCATE(Hscr%qlwl,(3,Hscr%nqlwl))
       ABI_ALLOCATE(Hscr%lwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))
       ABI_ALLOCATE(Hscr%uwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))

     case default
       write(msg,'(a,i0)')' scr_hdr_io: Wrong fform read = ',fform
       MSG_BUG(msg)
     end select
   end if !I_read

 else if (rdwr==2.or.rdwr==6) then
   ! === Writing the header of an unformatted file ===

   ! Write the abinit header, rewinding of the file (if any) is done here.
   if (accesswff==IO_MODE_FORTRAN) then
     call hdr_io(fform,Hscr%Hdr,rdwr,unt)
   else if (accesswff==IO_MODE_ETSF) then
     call hdr_io_etsf(fform,Hscr%Hdr,rdwr,unt)
   end if

   ! TODO should always use the latest version.
   write(unt)Hscr%title

   select case (fform)
   case (1002) ! Old epsilon^-1 file used in version <5.6
     write(unt)Hscr%npwe,Hscr%npwwfn_used,Hscr%nbnds_used,Hscr%nqibz,Hscr%nomega

   case (1102) ! File format for epsilon^-1, espilon, chi0 
     !TODO add new variables
     write(unt)Hscr%ID,Hscr%nI,Hscr%nJ,Hscr%tordering,Hscr%test_type
     write(unt)Hscr%npwe,Hscr%npwwfn_used,Hscr%nbnds_used,Hscr%nqibz,Hscr%nomega,Hscr%nqlwl 

   case (2002) ! Pole-fit file
     write(unt)Hscr%npwe,Hscr%npwwfn_used,Hscr%nbnds_used,Hscr%nqibz, &
&       Hscr%nomega,Hscr%npoles,Hscr%ncoeff

   case default
     write(msg,'(a,i5)')' scr_hdr_io: Wrong value for fform = ',fform
     MSG_BUG(msg)
   end select

   write(unt)Hscr%gvec(:,:)
   write(unt)Hscr%qibz(:,:)
   write(unt)Hscr%omega(:)

   ! === Add q-points for heads and wings for q-->0  ====
   if ((fform/=1002.and.fform/=2002).and.Hscr%nqlwl>0) then 
     write(unt)Hscr%qlwl(:,:)
     write(unt)Hscr%lwing(:,:,:)
     write(unt)Hscr%uwing(:,:,:)
   end if
  
 else if (rdwr==3.or.rdwr==4) then
   !
   ! === Echo the header to a formatted file ===
   write(unt,'(a)')&
&   ' ==============================================================================='
   if (rdwr==3) then 
     prtvol=0
     write(unt,'(a)') ' ECHO of part of the ABINIT-SCR file header '
   end if
   if (rdwr==4) then 
     prtvol=1
     write(unt,'(a)') ' ECHO of the ABINIT-SCR file header '
   end if

   call print_ScrHdr(Hscr,unit=unt,prtvol=prtvol,mode_paral='COLL')

   if (rdwr==3) write(unt,'(a)')' End the ECHO of part of the ABINIT-SCR file header '
   if (rdwr==4) write(unt,'(a)')' End the ECHO of the ABINIT-SCR file header '
   write(unt,'(a)')' ==============================================================================='

   ! * Echo the abinit header. Rewinding of the file (if any) is done here.
   if (prtvol>0) then
     if (accesswff==IO_MODE_FORTRAN) then
       call hdr_io(fform,Hscr%Hdr,rdwr,unt)
     else if (accesswff==IO_MODE_ETSF) then
       call hdr_io_etsf(fform,Hscr%Hdr,rdwr,unt)
     end if
   end if

 else 
   write(msg,'(a,i5)')' scr_hdr_io: Wrong value for rdwr = ',rdwr
   MSG_BUG(msg)
 end if ! read/write/echo

 DBG_EXIT("COLL")

 RETURN
 !
 ! === Something went wrong while reading! ===
 10 continue
 MSG_ERROR("The header of the (SCR|SUSC) file seems to be corrupted.")

end subroutine scr_hdr_io
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/print_ScrHdr
!! NAME
!! print_ScrHdr
!!
!! FUNCTION
!!  Prints info on the header of the SCR|SUSC file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_io_screening,mrgscr
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine print_ScrHdr(Hscr,header,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_ScrHdr'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(ScrHdr_type),intent(in) :: Hscr

!Local variables-------------------------------
!scalars
 integer :: iomega,iqibz,iqlwl,unt,verbose
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 unt=std_out; if (PRESENT(unit      )) unt    =unit
 verbose=0  ; if (PRESENT(prtvol    )) verbose=prtvol
 mode='COLL'; if (PRESENT(mode_paral)) mode   =mode_paral

 if (PRESENT(header)) then 
   msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
   call wrtout(unt,msg,mode)
 end if

 write(msg,'(1x,a)')TRIM(Hscr%title(1))
 call wrtout(unt,msg,mode)
 write(msg,'(1x,a)')TRIM(Hscr%title(2))
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Identifier                ',Hscr%ID
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Kxc kernel                ',Hscr%ikxc
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Treatment of q-->0 limit  ',Hscr%inclvkb
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' headform                  ',Hscr%headform
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' fform                     ',Hscr%fform
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' gwcalctyp                 ',Hscr%gwcalctyp
 call wrtout(unt,msg,mode)
 write(msg,'(a,2i8)')' Number of components      ',Hscr%nI,Hscr%nJ
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of q-points        ',Hscr%nqibz
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of q-directions    ',Hscr%nqlwl
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of frequencies     ',Hscr%nomega
 call wrtout(unt,msg,mode)
 if (Hscr%fform==2002) then
   write(msg,'(a)')    ' - Generalised pole model -'
   call wrtout(unt,msg,mode)
   write(msg,'(a,i8)') ' - Number of poles           ',Hscr%npoles
   call wrtout(unt,msg,mode)
   write(msg,'(a,i8)') ' - Number of coefficients    ',Hscr%ncoeff
   call wrtout(unt,msg,mode)
   write(msg,'(a,i8)') ' - with one phase value per G,G'''
   call wrtout(unt,msg,mode)
 end if
 write(msg,'(a,i8)') ' Number of bands used      ',Hscr%nbnds_used
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Dimension of matrix       ',Hscr%npwe
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of planewaves used ',Hscr%npwwfn_used
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Spectral method           ',Hscr%spmeth
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Test_type                 ',Hscr%test_type
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Time-ordering             ',Hscr%tordering
 call wrtout(unt,msg,mode)

 write(msg,'(a,es16.6)')' Scissor Energy             ',Hscr%soenergy
 call wrtout(unt,msg,mode)
 write(msg,'(a,es16.6)')' Spectral smearing          ',Hscr%spsmear
 call wrtout(unt,msg,mode)
 write(msg,'(a,es16.6)')' Complex Imaginary Shift    ',Hscr%zcut
 call wrtout(unt,msg,mode)

 if (verbose==0) then
   call wrtout(unt,' The header contains additional records.',mode)
 else
   write(msg,'(2a)')ch10,' q-points [r.l.u.]:'
   call wrtout(unt,msg,mode)
   do iqibz=1,Hscr%nqibz 
     write(msg,'(i5,3f12.6)')iqibz,Hscr%qibz(:,iqibz)
     call wrtout(unt,msg,mode)
   end do

   write(msg,'(2a)')ch10,' Frequencies used [eV]:'
   call wrtout(unt,msg,mode)
   do iomega=1,Hscr%nomega
     write(msg,'(i3,2f7.2)')iomega,REAL(Hscr%omega(iomega))*Ha_eV,AIMAG(Hscr%omega(iomega))*Ha_eV
     call wrtout(unt,msg,mode)
   end do

   if (Hscr%nqlwl>0) then
     write(msg,'(2a)')ch10,' q-points for long-wavelength limit [r.l.u.]:'
     call wrtout(unt,msg,mode)

     do iqlwl=1,Hscr%nqlwl
       write(msg,'(i5,3f12.6)')iqlwl,Hscr%qlwl(:,iqlwl)
       call wrtout(unt,msg,mode)
       if (verbose>0) then
         do iomega=1,Hscr%nomega
           write(msg,'(2x,a,i4,a,2f9.4,a)')&
&            ' Upper and lower wings at the ',iomega,'th omega',Hscr%omega(iomega)*Ha_eV,' [eV]'
           call wrtout(unt,msg,mode)
           call print_arr(Hscr%uwing(:,iomega,iqlwl),max_r=9,unit=ab_out)
           call print_arr(Hscr%lwing(:,iomega,iqlwl),max_r=9,unit=ab_out)
           write(msg,'(a)')ch10 
           call wrtout(unt,msg,mode)
         end do
       end if
     end do
   end if
  ! G-vectors are not printed out.
 end if 

end subroutine print_ScrHdr
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/reset_HSCR_
!! NAME
!!  reset_HSCR_ [PRIVATE]
!!
!! FUNCTION
!!  Initialize variables of Hscr using default values.
!!  in order to maintain backward compatibility.
!!
!! INPUTS
!!
!! OUTPUT
!!  (see side effects)
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine reset_HSCR_(Hscr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'reset_HSCR_'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ScrHdr_type),intent(inout) :: Hscr
! *************************************************************************

 !@ScrHdr_type
 if (Hscr%fform==1002.OR.Hscr%fform==2002) then 
   Hscr%ID        = 4         ! 4 for e^-1 since it was the only format available.
   Hscr%ikxc      = 0         ! 0 for None (RPA),
   Hscr%headform  = 56        ! Oldest one.
   Hscr%inclvkb   = 0         ! q-->0 treatment, 0 for None
   Hscr%gwcalctyp =-1         ! Not present in old fformat
   Hscr%nI        = 1         ! Collinear case.
   Hscr%nJ        = 1     
   Hscr%nqlwl     = 0         ! No. of q-->0 points.
   Hscr%spmeth    = 0         ! Was not released yet.
   Hscr%npoles    = 0
   Hscr%ncoeff    = 0
   Hscr%test_type = 0         ! None
   Hscr%tordering = 1         ! 1 for Time-Ordered, 2 for Advanced, 3 for Retarded.

   Hscr%soenergy  = zero      ! Not available.
   Hscr%spsmear   = zero      ! Not released yet.  
   Hscr%zcut      = 3.67493260d-03  ! Default, 0.1eV.
 end if

end subroutine reset_HSCR_
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/nullify_HScr
!! NAME
!!  nullify_HScr 
!!
!! FUNCTION
!!  Initialize all pointers to NULL.
!!
!! INPUTS
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  All pointers set to null
!!
!! PARENTS
!!      m_io_screening,m_screening
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine nullify_HScr(Hscr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_HScr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ScrHdr_type),intent(inout) :: Hscr

! *************************************************************************
 !@ScrHdr_type

!integer 
 nullify(Hscr%gvec)

!real
 nullify(Hscr%qibz)
 nullify(Hscr%qlwl)

!complex
 nullify(Hscr%lwing)
 nullify(Hscr%omega)
 nullify(Hscr%uwing)

 call hdr_nullify(Hscr%Hdr)

end subroutine nullify_HScr
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/init_ScrHdr
!! NAME
!!  init_ScrHdr
!!
!! FUNCTION
!!  Initialize the Hscr datatype and most of its content from the 
!!  Epsilonm1_parameters data type Ep.
!!
!! INPUTS
!!  ID=Identifier used to define the type of Response function (e^-1, chi0)
!!  ikxc=Integer flag definining the type of XC kernel (0 if None i.e RPA)
!!  test_type=Integer flag defining the type of probing charge (0 for None)
!!  tordering=The time-ordering of the Response function.
!!  gvec(3,Ep%npwe)=The G-vectors used. 
!!  Ep<Epsilonm1_parameters>=Parameters defining the calculation of the screening.
!!  Hdr_abinit<Hdr_type>=The abinit header.
!!
!! OUTPUT
!!  Hscr<type(ScrHdr_type)>=the header, initialized.
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine init_ScrHdr(ID,ikxc,test_type,tordering,title,ngvec,gvec,Ep,Hdr_abinit,Hscr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_ScrHdr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ID,ikxc,test_type,tordering,ngvec
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(Hdr_type),intent(in) :: Hdr_abinit
 type(ScrHdr_type),intent(inout) :: HScr
!arrays
 integer,intent(in) :: gvec(3,ngvec)
 character(len=80),intent(in) :: title(2)

! *************************************************************************

 !@ScrHdr_type
 ABI_CHECK(ngvec==Ep%npwe,'ngvec/=Ep%npwe')

 call nullify_HScr(HScr)

 ! === Deep copy of the abinit header ===
 call hdr_copy(Hdr_abinit,Hscr%Hdr)

 ! === Initialize quantities related to the screening file ===
 Hscr%ID         =ID           
 Hscr%ikxc       =ikxc         
 Hscr%inclvkb    =Ep%inclvkb      
 Hscr%fform      =1002
 !$Hscr%fform      =1102
 Hscr%headform   =HSCR_LATEST_HEADFORM
 Hscr%gwcalctyp  =Ep%gwcalctyp
 Hscr%nI         =Ep%nI
 Hscr%nJ         =Ep%nJ  
 Hscr%nqibz      =Ep%nqcalc  ! Ep%nqcalc==Ep%nqibz unless we splitted the calculation into different runs. 
 Hscr%nqlwl      =Ep%nqlwl        
 Hscr%nomega     =Ep%nomega       
 Hscr%npoles     =Ep%npoles
 Hscr%ncoeff     =Ep%ncoeff
 Hscr%nbnds_used =Ep%nbnds   
 Hscr%npwe       =Ep%npwe 
 Hscr%npwwfn_used=Ep%npwwfn 
 Hscr%spmeth     =Ep%spmeth
 Hscr%test_type  =test_type    
 Hscr%tordering  =tordering    

 Hscr%soenergy   =Ep%soenergy
 Hscr%spsmear    =Ep%spsmear
 Hscr%zcut       =Ep%zcut

 Hscr%title(:)=title(:)

 ABI_ALLOCATE(Hscr%gvec,(3,Ep%npwe))
 Hscr%gvec=gvec(1:3,1:Ep%npwe)

 ABI_ALLOCATE(Hscr%qibz,(3,Ep%nqcalc))
 Hscr%qibz=Ep%qcalc
 
 ABI_ALLOCATE(Hscr%qlwl,(3,Ep%nqlwl))
 Hscr%qlwl=Ep%qlwl

 ABI_ALLOCATE(Hscr%omega,(Ep%nomega))
 Hscr%omega=Ep%omega 

 ABI_ALLOCATE(Hscr%lwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))
 ABI_ALLOCATE(Hscr%uwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))

 !TODO these quantities should be initialized correctly from input arguments
 if (Ep%nqlwl>0) then
   Hscr%lwing=czero
   Hscr%uwing=czero
 end if

end subroutine init_ScrHdr
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/scrhdr_comm
!! NAME
!! scrhdr_comm
!!
!! FUNCTION
!! This subroutine transmit the header structured datatype initialized 
!! on one processor (or a group of processor), to the other processors. 
!! It also allocates the needed part of the header.
!!
!! INPUTS
!!  master=ID of the master node.
!!  my_rank=ID of the node that receives the data.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  (no output)
!!
!! SIDE EFFECTS
!!  Hscr<type(ScrHdr_type)>=the SCR header. For the master, it is already
!!   initialized entirely, while for the other procs, everything has
!!   to be transmitted.
!!
!! NOTES
!! This routine is called only in the case of MPI version of the code.
!!
!! PARENTS
!!      m_io_screening,m_screen,setup_bse
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine scrhdr_comm(Hscr,master,my_rank,comm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scrhdr_comm'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: master,my_rank,comm
 type(ScrHdr_type),intent(inout) :: Hscr

!Local variables-------------------------------
 integer :: ierr

! *************************************************************************

 DBG_ENTER("COLL")
 !@ScrHdr_type  

!scalars
 call xcast_mpi(Hscr%ID,         master,comm,ierr)
 call xcast_mpi(Hscr%ikxc,       master,comm,ierr)
 call xcast_mpi(Hscr%inclvkb,    master,comm,ierr)
 call xcast_mpi(Hscr%headform,   master,comm,ierr)
 call xcast_mpi(Hscr%fform,      master,comm,ierr)
 call xcast_mpi(Hscr%gwcalctyp,  master,comm,ierr)
 call xcast_mpi(Hscr%nI,         master,comm,ierr)
 call xcast_mpi(Hscr%nJ,         master,comm,ierr)
 call xcast_mpi(Hscr%nqibz,      master,comm,ierr)
 call xcast_mpi(Hscr%nqlwl,      master,comm,ierr)
 call xcast_mpi(Hscr%nomega,     master,comm,ierr)
 call xcast_mpi(Hscr%nbnds_used, master,comm,ierr)
 call xcast_mpi(Hscr%npwe,       master,comm,ierr)
 call xcast_mpi(Hscr%npwwfn_used,master,comm,ierr)
 call xcast_mpi(Hscr%spmeth,     master,comm,ierr)
 call xcast_mpi(Hscr%test_type,  master,comm,ierr)
 call xcast_mpi(Hscr%tordering,  master,comm,ierr)

 call xcast_mpi(Hscr%soenergy,   master,comm,ierr)
 call xcast_mpi(Hscr%spsmear,    master,comm,ierr)
 call xcast_mpi(Hscr%zcut,       master,comm,ierr)

 ! Communicate the Abinit header.
 call hdr_comm(Hscr%Hdr,master,my_rank,comm)

!arrays
 call xcast_mpi(Hscr%title, master,comm,ierr)

 if (my_rank/=master) then 
   ABI_ALLOCATE(Hscr%gvec,(3,Hscr%npwe))
   ABI_ALLOCATE(Hscr%qibz,(3,Hscr%nqibz))
   ABI_ALLOCATE(Hscr%qlwl,(3,Hscr%nqlwl))
   ABI_ALLOCATE(Hscr%omega,(Hscr%nomega))
   ABI_ALLOCATE(Hscr%lwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))
   ABI_ALLOCATE(Hscr%uwing,(Hscr%npwe,Hscr%nomega,Hscr%nqlwl))
 end if

 call xcast_mpi(Hscr%gvec, master,comm,ierr)
 call xcast_mpi(Hscr%qibz, master,comm,ierr)
 call xcast_mpi(Hscr%qlwl, master,comm,ierr)
 call xcast_mpi(Hscr%omega,master,comm,ierr)

 if (Hscr%nqlwl>0) then
   call xcast_mpi(Hscr%lwing,master,comm,ierr)
   call xcast_mpi(Hscr%uwing,master,comm,ierr)
 end if

 DBG_EXIT("COLL")

end subroutine scrhdr_comm
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/free_scrhdr
!! NAME
!! free_scrhdr
!!
!! FUNCTION
!! Deallocate the components of the header structured datatype
!!
!! INPUTS
!! hdr <type(hdr_type)>=the header
!!
!! OUTPUT
!!  (only deallocate)
!!
!! PARENTS
!!      m_io_screening,m_screen,m_screening,mrgscr,screening,setup_bse
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine free_scrhdr(Hscr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'free_scrhdr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ScrHdr_type),intent(inout) :: Hscr

! ************************************************************************* 

 !@ScrHdr_type
 DBG_ENTER("COLL")

 if (associated(Hscr%gvec ))   then
   ABI_DEALLOCATE(Hscr%gvec)
 end if
 if (associated(Hscr%qibz ))   then
   ABI_DEALLOCATE(Hscr%qibz)
 end if
 if (associated(Hscr%qlwl ))   then
   ABI_DEALLOCATE(Hscr%qlwl)
 end if
 if (associated(Hscr%omega))   then
   ABI_DEALLOCATE(Hscr%omega)
 end if
                                                   
 if (associated(Hscr%lwing))   then
   ABI_DEALLOCATE(Hscr%lwing)
 end if
 if (associated(Hscr%uwing))   then
   ABI_DEALLOCATE(Hscr%uwing)
 end if

 call hdr_clean(Hscr%Hdr)

 DBG_EXIT("COLL")

end subroutine free_scrhdr
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/copy_ScrHdr
!! NAME
!! copy_ScrHdr
!!
!! FUNCTION
!! Deep copy of the header of the _SCR or _SUSC file.
!!
!! INPUTS
!!
!! PARENTS
!!      m_io_screening,m_screen,m_screening,mrgscr
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine copy_ScrHdr(Hscr_in,Hscr_cp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_ScrHdr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ScrHdr_type),intent(in) :: Hscr_in
 type(ScrHdr_type),intent(inout) :: Hscr_cp

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
! *************************************************************************

 !@ScrHdr_type
 if (Hscr_in%headform>HSCR_LATEST_HEADFORM.or.Hscr_in%headform<0) then 
   write(msg,'(a,i4,2a,i10,2a)')&
&    ' copy_ScrHdr deals with versions of the header only up to ',HSCR_LATEST_HEADFORM,ch10,&
&    ' However headform= ',Hscr_in%headform,ch10,&
&    ' Change the source to add the changes done in the new version. '  
   MSG_BUG(msg)
 end if

 ! === Integer values ===
 Hscr_cp%ID          = Hscr_in%ID           
 Hscr_cp%ikxc        = Hscr_in%ikxc         
 Hscr_cp%inclvkb     = Hscr_in%inclvkb      
 Hscr_cp%headform    = Hscr_in%headform        
 Hscr_cp%fform       = Hscr_in%fform        
 Hscr_cp%gwcalctyp   = Hscr_in%gwcalctyp    
 Hscr_cp%nI          = Hscr_in%nI
 Hscr_cp%nJ          = Hscr_in%nJ
 Hscr_cp%nqibz       = Hscr_in%nqibz        
 Hscr_cp%nqlwl       = Hscr_in%nqlwl        
 Hscr_cp%nomega      = Hscr_in%nomega       
 Hscr_cp%npoles      = Hscr_in%npoles      
 Hscr_cp%ncoeff      = Hscr_in%ncoeff       
 Hscr_cp%nbnds_used  = Hscr_in%nbnds_used   
 Hscr_cp%npwe        = Hscr_in%npwe         
 Hscr_cp%npwwfn_used = Hscr_in%npwwfn_used  
 Hscr_cp%spmeth      = Hscr_in%spmeth       
 Hscr_cp%test_type   = Hscr_in%test_type    
 Hscr_cp%tordering   = Hscr_in%tordering    

 ! === Real variables ====
 Hscr_cp%soenergy = Hscr_in%soenergy    
 Hscr_cp%spsmear  = Hscr_in%spsmear     
 Hscr_cp%zcut     = Hscr_in%zcut        

 ! === Copy the abinit Header ===
 call hdr_copy(Hscr_in%Hdr,Hscr_cp%Hdr)

 Hscr_cp%title(:) = Hscr_in%title(:)

 ! Copy pointers.
 call deep_copy( Hscr_in%gvec , Hscr_cp%gvec  )

 call deep_copy( Hscr_in%qibz , Hscr_cp%qibz  )
 call deep_copy( Hscr_in%qlwl , Hscr_cp%qlwl  )

 call deep_copy( Hscr_in%lwing, Hscr_cp%lwing )
 call deep_copy( Hscr_in%omega, Hscr_cp%omega )
 call deep_copy( Hscr_in%uwing, Hscr_cp%uwing )

end subroutine copy_ScrHdr
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/merge_ScrHdr
!! NAME
!! merge_ScrHdr
!!
!! FUNCTION
!! This subroutine merges diffrent header structured variable (Scrhdr)
!!
!! INPUTS
!!  Hscr_in(:) <ScrHdr_type)>=List of headers to be merged.
!!
!! OUTPUT
!!  Hscr_out<ScrHdr_type>=The output merged header.
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine merge_ScrHdr(Hscr_in,Hscr_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'merge_ScrHdr'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ScrHdr_type),intent(out) :: Hscr_out
 type(ScrHdr_type),intent(in) :: Hscr_in(:)

!Local variables-------------------------------
!scalars
 integer :: nhds,restart,restartpaw,ihd,ii,nqtot,nqneq
 logical :: isok
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: qset(:,:)
! *************************************************************************

 !@ScrHdr_type
 DBG_ENTER("COLL")

 nhds=SIZE(Hscr_in)

 ! === Make an initial copy of the header ===
 ! * If multiple headers, select the header containing q-->0 so that we copy also heads and wings
 ii=imax_loc(Hscr_in(:)%nqlwl)
 call copy_ScrHdr(Hscr_in(ii),Hscr_out)
 !do ihd=1,nhds
 !  if (ihd/=ii.and.(Hscr_in(ihd)%nqlwl/=0).and.Hscr_in(ihd)%headform/=0) then
 !    msg='Only a single header should contain heads and wings'
 !    MSG_ERROR(msg)
 !  end if
 !end do

 if (nhds==1) RETURN

 ! === Check the consistency of the abinit Headers ===
 ! * FFT grid might be q-point dependent so we stop only when restart==0
 isok=.TRUE.
 if (Hscr_in(1)%fform/=1102) then
   do ihd=2,nhds
     call hdr_check(Hscr_in(1)%fform,Hscr_in(ihd)%fform,Hscr_in(1)%Hdr,Hscr_in(ihd)%Hdr,'COLL',restart,restartpaw)
     if (restart==0) then
       isok=.FALSE.
       write(msg,'(a,i3,a)')' Abinit header no.',ihd,' is not consistent with the first header '
       MSG_WARNING(msg)
     end if
   end do
 else
   write(msg,'(a)')' You are merging _SUS files and hdr_check is not performed - proceed at you own peril!!! '
   MSG_WARNING(msg)
 end if

 if (.not.isok) then
   MSG_ERROR('Cannot continue, Check headers')
 end if

 if (Hscr_in(1)%fform/=1102) then

   ! === Now check variables related to polarizability|epsilon^{-1} ===
   ! 1) Tests quantities that must be equal 
   ii = assert_eq(Hscr_in(:)%ID,       'Headers have different Identifiers')
   ii = assert_eq(Hscr_in(:)%ikxc,     'Headers have different ikxc'       )
   ii = assert_eq(Hscr_in(:)%headform, 'Headers have different headform'   )
   ii = assert_eq(Hscr_in(:)%fform,    'Headers have different fform'      )
   ii = assert_eq(Hscr_in(:)%gwcalctyp,'Headers have different gwcalctyp'  )
   ii = assert_eq(Hscr_in(:)%nI,       'Headers have different nI'         )
   ii = assert_eq(Hscr_in(:)%nJ,       'Headers have different nJ'         )
   ii = assert_eq(Hscr_in(:)%nomega,   'Headers have different nomega'     )
   ii = assert_eq(Hscr_in(:)%test_type,'Headers have different test_type'  )
   ii = assert_eq(Hscr_in(:)%tordering,'Headers have different tordering'  )
  
   ! This is not mandatory but makes life easier!
   ii = assert_eq(Hscr_in(:)%npwe,'Headers have different number of G-vectors'  )
  
   do ihd=2,nhds
     if (ANY(ABS(Hscr_in(ihd)%omega-Hscr_in(1)%omega)>tol6)) then
       write(msg,'(a,i3,a)')' Frequencies in the first and the ',ihd,'-th header differ '
       MSG_ERROR(msg)
     end if
     if (ANY(Hscr_in(ihd)%gvec(:,:)-Hscr_in(1)%gvec(:,:)/=0)) then 
       write(msg,'(a,i3,a)')' Incompatible G-vector list found in the ',ihd,'-th header '
       MSG_ERROR(msg)
     end if
   end do !ihd
  
   ! === If error is not fatal, just warn ===
   if (ANY(Hscr_in(:)%npwwfn_used/=Hscr_in(1)%npwwfn_used)) then
     msg = '  Files have been produced with a different number of planewaves for the wavefunctions. '
     MSG_COMMENT(msg)
   end if
  
   if (ANY(Hscr_in(:)%nbnds_used/=Hscr_in(1)%nbnds_used)) then
     msg = '  Files have been produced with a different number of bands. '
     MSG_COMMENT(msg)
   end if

   if (ANY(Hscr_in(:)%spmeth/=Hscr_in(1)%spmeth)) then
     msg = '  Files have been produced with different algorithms. '
     MSG_COMMENT(msg)
   end if
  
   if (ANY(ABS(Hscr_in(:)%soenergy-Hscr_in(1)%soenergy)>tol6)) then
     msg = ' Files have benn produced with different values of soenergy. '
     MSG_COMMENT(msg)
   end if

   if (ANY(ABS(Hscr_in(:)%spsmear-Hscr_in(1)%spsmear)>tol6)) then
     msg = ' Files have been produced with different values of spsmear. '
     MSG_COMMENT(msg)
   end if
  
   if (ANY(ABS(Hscr_in(:)%zcut-Hscr_in(1)%zcut)>tol6)) then
     msg = ' Files have been produced with different values of zcut. '
     MSG_COMMENT(msg)
   end if

 end if ! SUS file check


 ! === Now merge the list of q-points ===
 nqtot=SUM(Hscr_in(:)%nqibz)
 ABI_ALLOCATE(qset,(3,nqtot))
 ii=0
 do ihd=1,nhds
   qset(:,ii+1:ii+Hscr_in(ihd)%nqibz)=Hscr_in(ihd)%qibz(:,:)
   ii=ii+Hscr_in(ihd)%nqibz
 end do

 call remove_copies(nqtot,qset,nqneq,isequalk)

 if (nqneq/=nqtot) then
   write(msg,'(6a,i3,a,i3,a)')ch10,&
&    ' merge_ScrHdr: COMMENT ',ch10,&
&    ' Headers contain duplicated q-points ',ch10,&
&    ' Found ',nqneq,' distinct q-points among the total ',nqtot,' points reported in the headers. ' 
   call wrtout(std_out,msg,'COLL')
 end if

 Hscr_out%nqibz = nqneq
 ABI_DEALLOCATE(Hscr_out%qibz)
 ABI_ALLOCATE(Hscr_out%qibz,(3,nqneq))
 Hscr_out%qibz(:,:)=qset(:,1:nqneq)
 ABI_DEALLOCATE(qset)

 DBG_EXIT("COLL")

end subroutine merge_ScrHdr
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/write_screening
!! NAME
!! write_screening
!!
!! FUNCTION
!! For a single q-point, write either \tilde epsilon^{-1} on the _SCR file 
!! or chi0 on the _SUSC file. The file is supposed to have been open in the calling routine.
!!
!! INPUTS
!!  unt=The unit number of the file to be written (supposed to be already open)
!!  epsm1(npwe,npwe,nomega)=The matrix to be written, for different frequencies, and a single q-point.
!!  nomega=Number of frequencies
!!  npwe=Number of plane waves in epsm1.
!!  accesswff=Integer flag defining the format of the output file. Available options:
!!    1--> Plain Fortran file 
!!    3--> ETSF format (TODO not yet coded)
!!
!! NOTES
!!  On some architecture, the code crashes when trying to write or read a record containing the 
!!  entire (G1,G2) matrix thus we use smaller records containing the columns of the two-point function.
!!
!! OUTPUT
!!  (only writing on file)
!!
!! PARENTS
!!      m_screen,m_screening,mrgscr,screening
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine write_screening(unt,accesswff,npwe,nomega,epsm1)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_screening'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwe,unt,accesswff
!arrays
 complex(gwpc),intent(in) :: epsm1(npwe,npwe,nomega) 

!Local variables-------------------------------
!scalars
 integer :: ipwe,iomega,istat
 logical :: ltest
!arrays
 complex(dpc),allocatable :: epsm1d(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 ltest=(accesswff==IO_MODE_FORTRAN.or.accesswff==IO_MODE_ETSF)
 ABI_CHECK(ltest,'Wrong value for accesswff')
 ABI_CHECK(accesswff==IO_MODE_FORTRAN,'accesswff=3 not coded')
 !
 ! Write a record for each omega, Always use double precision.
 ABI_ALLOCATE(epsm1d,(npwe,1))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,' out of memory in epsm1d')

 do iomega=1,nomega
   do ipwe=1,npwe
     epsm1d(:,1)=epsm1(:,ipwe,iomega) !spc ==> dpc
     write(unt)epsm1d(1:npwe,1)
   end do
 end do

 ABI_DEALLOCATE(epsm1d)

 DBG_EXIT("COLL")

end subroutine write_screening
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/write_pole_screening
!! NAME
!! write_pole_screening
!!
!! FUNCTION
!! For a single q-point, write the pole coefficients of \epsilon^{-1} to file
!! The file is supposed to have been open in the calling routine.
!!
!! INPUTS
!!  unt=The unit number of the file to be written (supposed to be already open)
!!  epsm1(npwe,npwe,nomega)=The matrix to be written, for different frequencies, and a single q-point.
!!  nomega=Number of poles (+1 for phase)
!!  npwe=Number of plane waves in epsm1.
!!  accesswff=Integer flag defining the format of the output file. Available options:
!!    1--> Plain Fortran file 
!!    3--> ETSF format (TODO not yet coded)
!!
!! NOTES
!!  On some architecture, the code crashes when trying to write or read a record containing the 
!!  entire (G1,G2) matrix thus we use smaller records containing the columns of the two-point function.
!!
!! OUTPUT
!!  (only writing on file)
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine write_pole_screening(unt,accesswff,npwe,ncoeff,epsm1)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_pole_screening'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncoeff,npwe,unt,accesswff
!arrays
 real(gwpc),intent(in) :: epsm1(npwe,npwe,ncoeff) 

!Local variables-------------------------------
!scalars
 integer :: ipwe,icoeff,istat
 logical :: ltest
!arrays
 real(dpc),allocatable :: epsm1d(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 ltest=(accesswff==IO_MODE_FORTRAN.or.accesswff==IO_MODE_ETSF)
 ABI_CHECK(ltest,'Wrong value for accesswff')
 ABI_CHECK(accesswff==IO_MODE_FORTRAN,'accesswff=3 not coded')
 !
 ! Write a record for each omega, Always use double precision.
 ABI_ALLOCATE(epsm1d,(npwe,1))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,' out of memory in epsm1d')

 do icoeff=1,ncoeff
   do ipwe=1,npwe
     epsm1d(:,1)=epsm1(:,ipwe,icoeff) !spc ==> dpc
     write(unt)epsm1d(1:npwe,1)
   end do
 end do

 ABI_DEALLOCATE(epsm1d)

 DBG_EXIT("COLL")

end subroutine write_pole_screening
!!***

!----------------------------------------------------------------------


!!****f* m_io_screening/read_screening
!! NAME
!! read_screening
!!
!! FUNCTION
!! Read either a screening (\tilde epsilon^{-1}) file in the SCR format or 
!! the irreducible polarizability (chi0) in the SUSC format.
!!
!! INPUTS
!!  accesswff=Integer flag defining the format of the output file. Available options:
!!    1--> Plain Fortran file 
!!    3--> ETSF format (TODO not yet coded)
!!  iqiA[optional]=Used if only a particular q-point is required. In this case iqiA define the index
!!   of the required q-point in the array qibz(3,Hscr%nqibz)
!!  nqibzA=number of asked q-points (used to dimension the output arrays). 
!!   Equal to Hscr%nqibz if the full matrix is required
!!  comm=MPI communicator.
!!  npweA=number of asked planewaves
!!  nomegaA=number of asked frequencies
!!
!! OUTPUT
!!  epsm1(npweA,npweA,nomegaA,nqibzA) = \tilde\epsilon^{-1}(Ng,Ng,Nw,Nq)
!!
!! NOTES
!!  * If the epsilon matrix read is bigger than npweA x npweA, it will be truncated; 
!!    if it is smaller, an error will occur
!!  * If the number of frequencies asked for is smaller than that reported in the file, the matrix 
!!    will be truncated. If nomegaA > Hscr%nomega an error will occur 
!! 
!! PARENTS
!!      m_screen,m_screening,mrgscr
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine read_screening(fname,npweA,nqibzA,nomegaA,epsm1,accesswff,comm,&
& iqiA) ! Optional

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_screening'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none
    
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,nomegaA,npweA,nqibzA,comm
 integer,optional,intent(in) :: iqiA
 character(len=fnlen),intent(in) :: fname
!arrays
 complex(gwpc),intent(inout) :: epsm1(npweA,npweA,nomegaA,nqibzA)
  
!Local variables-------------------------------
!scalars
 integer,parameter :: master0=0
 integer :: ipwe,fform,iomega,ios,iqibz,istat,unt
 integer :: rdwr,my_rank,nprocs
#if 0 && defined HAVE_FC_STREAM_IO
 integer :: recsize
 integer(kind=16) mypos,startpos,big_npwe
#endif
#ifdef HAVE_MPI_IO
 integer :: ierr,test_fform
 integer :: bsize_frm,mpi_type_frm
 integer :: mpi_fh,buf_dim,mat_ggw,mat_ggwq
 integer(MPI_OFFSET_KIND) :: offset,displ_wq
 !complex(dpc) :: ctmp
#endif
 character(len=500) :: msg
 logical :: read_qslice
 type(ScrHdr_type) :: Hscr
!arrays
#ifdef HAVE_MPI_IO
!integer :: statux(MPI_STATUS_SIZE)
 integer(MPI_OFFSET_KIND),allocatable :: offset_wq(:,:)
#if ! defined HAVE_GW_DPC
 complex(dpc),allocatable :: bufdc(:,:,:)
#endif
#endif
 complex(dpc),allocatable :: bufdc2d(:,:)
!#if ! defined HAVE_GW_DPC
 complex(dpc),allocatable :: bufdc3d(:,:,:)
!#endif

! *************************************************************************

 DBG_ENTER("COLL")

 my_rank = xcomm_rank(comm)
 nprocs  = xcomm_size(comm)
 !
 ! === Read the header of the file ===
 rdwr=1 
 !if (accesswff==IO_MODE_MPI .or. .TRUE.) then
 if (accesswff==IO_MODE_MPI) then
#if defined HAVE_MPI_IO
   bsize_frm    = xmpio_bsize_frm    ! bsize_frm= Byte length of the Fortran record marker.
   mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

   if (my_rank==master0) then ! Master reads the header via Fortran IO then bcast the data.
     unt=get_unit()
     open(unit=unt,file=fname,status='old',form='unformatted',iostat=ios)
     msg = ' Opening file '//TRIM(fname)//' as old unformatted '
     ABI_CHECK(ios==0,msg)
     call scr_hdr_io(fform,rdwr,unt,xmpi_self,master0,IO_MODE_FORTRAN,Hscr)
     close(unt)
   end if

   call scrhdr_comm(HScr,master0,my_rank,comm)

   call xcast_mpi(fform,master0,comm,ierr)
   ! 
   ! Open the file with MPI-IO
   call MPI_FILE_OPEN(comm,fname, MPI_MODE_RDONLY, MPI_INFO_NULL,mpi_fh,ierr)
   ABI_CHECK_MPI(ierr,"MPI_FILE_OPEN") 
   !
   ! Retrieve the offset of the section immediately below the header.
   call scrhdr_mpio_skip(mpi_fh,test_fform,offset) 
   ABI_CHECK(test_fform==fform,"mismatch in fform!")
   !
   ! Offsets of the Fortran markers corresponding to the (w,q) slices.
   ABI_ALLOCATE(offset_wq,(HScr%nomega,HScr%nqibz))
   displ_wq = offset
   do iqibz=1,Hscr%nqibz
     do iomega=1,Hscr%nomega
       offset_wq(iomega,iqibz) = displ_wq 
       displ_wq = displ_wq + Hscr%npwe**2 * xmpi_bsize_dpc + Hscr%npwe * 2 * bsize_frm
     end do
   end do

   !do iqibz=1,Hscr%nqibz
   !  do iomega=1,Hscr%nomega
   !    displ_wq = offset_wq(iqibz,iomega)
   !    call xmpio_read_frm(mpi_fh,displ_wq,xmpio_at_all,bsize_frm,fmarker,ierr,advance=.FALSE.)

   !    call MPI_FILE_READ_AT(mpi_fh,displ_wq,ctmp,1,MPI_DOUBLE_COMPLEX,statux,ierr)
   !    write(std_out,*)"fmarker, ctmp ",fmarker,ctmp
   !  end do
   !end do

#else
   MSG_ERROR("MPI-IO support not enabled at configure-time")
#endif
 else 
   ! Plain Fortran IO, all nodes read.
   unt=get_unit()
   open(unit=unt,file=fname,status='old',form='unformatted',iostat=ios)
   msg = ' Opening file '//TRIM(fname)//' as old unformatted '
   ABI_CHECK(ios==0,msg)
   call scr_hdr_io(fform,rdwr,unt,comm,master0,accesswff,Hscr)
 end if
 !
 ! === Slice or full array? ===
 read_qslice=.FALSE.
 if (PRESENT(iqiA)) then 
   read_qslice=.TRUE.
   write(msg,'(a,i0,2a)')' reading q-slice corresponding to iq = ',iqiA,' from file: ',TRIM(fname)
   call wrtout(std_out,msg,'PERS')
   if (iqiA<=0.or.iqiA>Hscr%nqibz) then
     MSG_BUG('iqiA out of range')
   end if
 end if 
 !
 ! === Do some check ===
 if (Hscr%npwe>npweA) then
   write(msg,'(a,i8,2a,i8)')&
&    ' Total number of G-vectors reported on file = ',Hscr%npwe,ch10,&
&    ' Reading a smaller matrix of dimension      = ',npweA
   MSG_COMMENT(msg)
 end if

 if (npweA>Hscr%npwe) then
   write(msg,'(2(a,i0))')' Dimension of matrix = ',Hscr%npwe," requiring a too big matrix = ",npweA
   MSG_ERROR(msg)
 end if

 ABI_CHECK(nqibzA  <= Hscr%nqibz, ' Requiring too much q-points ')
 ABI_CHECK(nomegaA <= Hscr%nomega,' Requiring too much frequencies ')

 !if (accesswff==IO_MODE_MPI .or. .TRUE.) then
 if (accesswff==IO_MODE_MPI) then
#ifdef HAVE_MPI_IO 
   if (read_qslice) then

     call xmpio_create_fsubarray_3D((/HScr%npwe,HScr%npwe,HScr%nomega/),(/npweA,npweA,nomegaA/),MPI_DOUBLE_COMPLEX,mat_ggw,ierr)
     ABI_CHECK_MPI(ierr,"xmpio_create_fsubarray_3D") 

     offset = offset_wq(1,iqiA) + bsize_frm

     call MPI_FILE_SET_VIEW(mpi_fh,offset,MPI_BYTE,mat_ggw,"native",MPI_INFO_NULL,ierr)
     ABI_CHECK_MPI(ierr,"SET_VIEW") 

     buf_dim = (npweA)**2 * nomegaA 

#ifdef HAVE_GW_DPC
     ! Read in-place.
     call MPI_FILE_READ_ALL(mpi_fh,epsm1,buf_dim,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
     ABI_CHECK_MPI(ierr,"FILE_READ_ALL") 
#else
     ! Have to allocate workspace for dpc data.
     ABI_ALLOCATE(bufdc3d,(npweA,npweA,nomegaA))
     call MPI_FILE_READ_ALL(mpi_fh,bufdc3d,buf_dim,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
     ABI_CHECK_MPI(ierr,"FILE_READ_ALL") 
     epsm1(:,:,:,1) = bufdc3d
     ABI_DEALLOCATE(bufdc3d)
#endif
     call MPI_TYPE_FREE(mat_ggw,ierr)

   else  ! Full matrix (G,G',w,q) is needed.

#ifdef HAVE_GW_DPC
     ! Can read all data at once.
     MSG_COMMENT("Using 4D Access")

     call xmpio_create_fsubarray_4D((/HScr%npwe,HScr%npwe,HScr%nomega,HScr%nqibz/),(/npweA,npweA,nomegaA,HScr%nqibz/),&
&      MPI_DOUBLE_COMPLEX,mat_ggwq,ierr)
     ABI_CHECK_MPI(ierr,"xmpio_create_fsubarray_4D") 

     offset = offset_wq(1,1) + bsize_frm
     call MPI_FILE_SET_VIEW(mpi_fh,offset,MPI_BYTE,mat_ggwq,"native",MPI_INFO_NULL,ierr)
     ABI_CHECK_MPI(ierr,"SET_VIEW") 
                                                                                                       
     buf_dim = (npweA)**2 * nomegaA * HScr%nqibz
     call MPI_FILE_READ_ALL(mpi_fh,epsm1,buf_dim,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
     ABI_CHECK_MPI(ierr,"FILE_READ_ALL") 

     call MPI_TYPE_FREE(mat_ggwq,ierr)

     !do iqibz=1,Hscr%nqibz
     !  do iomega=1,Hscr%nomega
     !    write(std_out,*)epsm1(1:2,1:2,iomega,iqibz)
     !  end do
     !end do
#else
     ! Have to allocate workspace for dpc data.
     ABI_ALLOCATE(bufdc3d,(npweA,npweA,nomegaA))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,'out of memory in bufdc3d, change the code to loop over frequencies!')

     call xmpio_create_fsubarray_3D((/HScr%npwe,HScr%npwe,HScr%nomega/),(/npweA,npweA,nomegaA/),MPI_DOUBLE_COMPLEX,mat_ggw,ierr)
     ABI_CHECK_MPI(ierr,"xmpio_create_fsubarray_3D") 

     do iqibz=1,Hscr%nqibz
       offset = offset_wq(1,iqibz) + bsize_frm
       call MPI_FILE_SET_VIEW(mpi_fh,offset,MPI_BYTE,mat_ggw,"native",MPI_INFO_NULL,ierr)
       ABI_CHECK_MPI(ierr,"SET_VIEW") 

       buf_dim = (2*npweA)**2 * nomegaA
       call MPI_FILE_READ_ALL(mpi_fh,bufdc3d,buf_dim,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
       ABI_CHECK_MPI(ierr,"FILE_READ_ALL")

       epsm1(:,:,:,iqibz) = bufdc3d
       !do iomega=1,Hscr%nomega
       !  write(std_out,*)epsm1(1:2,1:2,iomega,iqibz)
       !end do
     end do

     call MPI_TYPE_FREE(mat_ggw,ierr)
     ABI_DEALLOCATE(bufdc3d)
#endif
   end if

   call MPI_FILE_CLOSE(mpi_fh,ierr)
   ABI_DEALLOCATE(offset_wq)
#endif

 else 
   !
   ! === Now read epsilon^-1 ===
   ! Allocate a single column to save memory.
   ! TODO re-merge the two cases.
   ABI_ALLOCATE(bufdc2d,(Hscr%npwe,1))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,'out of memory in bufdc2d')
   !
   ! === Two coding for different case just to keep it readable ===
   ! * Communication is done inside the loops to avoid problems with 
   !   the size of the MPI packet. Much slower but safer.
   SELECT CASE (read_qslice)

   CASE (.TRUE.) ! Read only a slice of the full array (useful if the entire array is huge).
     !if (dim_wings==1) STOP 'not implemented'
     !TODO this has to be done in a cleaner way.
#if 0 && defined HAVE_FC_STREAM_IO
!    Loop using stream I/O (Fortran 2003) to jump in file
!    First find the byte position after the header in the file by matching the first record
     read(unt,ERR=10) bufdc2d(1:Hscr%npwe,1)
     close(unt)
     open(unit=unt,file=fname,status='old',form='unformatted',access='stream',iostat=ios)     
     recsize = sizeof(epsm1(1,1,1,1)) ! TODO fix this: Nonstandard fortran
     ! Find starting position
     startpos=0
     do
       startpos = startpos + 1
       read(unt,pos=startpos,ERR=10) epsm1(1:npweA,1,1,1)
       if (ALL(epsm1(1:npweA,1,1,1)==bufdc2d(1:npweA,1))) EXIT
     end do
     ! Find stride
     read(unt,pos=startpos,ERR=10) epsm1(1:Hscr%npwe,1,1,1)
     inquire(unt,pos=mypos)
     write(std_out,'(2(a,I0))') ' Stride is: ',mypos-startpos,' should be: ',Hscr%npwe*recsize 
     write(std_out,'(a,I0)') ' Header ends at: ',startpos-1
     big_npwe = Hscr%npwe
     startpos = startpos + recsize*(big_npwe*big_npwe*Hscr%nomega*(iqiA-1)) &
&                        + 8*(big_npwe*Hscr%nomega*(iqiA-1))
     write(std_out,'(2(a,I0))') ' Entering read loop, recsize: ',recsize,' npweA: ',npweA
     write(std_out,'(2(a,I0))') '                     Hscr%npwe: ',Hscr%npwe,' startpos: ',startpos
     do iomega=1,nomegaA
       if (mod(iomega,10)==0) write(std_out,'(2(a,I0))') ' Reading omega: ',iomega
       do ipwe=1,npweA
         ! Fortran record are written with 4B header and 4B footer (64 bit machines)
         ! Stride is the record length*Hscr%npwe + 8 bytes at the end
         mypos = startpos + recsize*(big_npwe*big_npwe*(iomega-1) + big_npwe*(ipwe-1)) &
&                         + 8*big_npwe*(iomega-1) + 8*(ipwe-1)
         read(unt,pos=mypos,ERR=10) epsm1(1:npweA,ipwe,iomega,1)
       end do
     end do
#else 
     qread_loop: &
&    do iqibz=1,Hscr%nqibz
       if (iqibz==iqiA) then 
         do iomega=1,nomegaA
           do ipwe=1,Hscr%npwe
             read(unt,ERR=10) bufdc2d(1:Hscr%npwe,1)
             if (ipwe<=npweA) epsm1(1:npweA,ipwe,iomega,1)=bufdc2d(1:npweA,1)
           end do
         end do
         EXIT qread_loop ! Got data. Do not need to read file till the end.
       else ! Skip other q-points i.e bufdc2d(1:Hscr%npwe,1:Hscr%npwe)
         do iomega=1,Hscr%nomega
           do ipwe=1,Hscr%npwe
            read(unt,ERR=10)
           end do
         end do
       end if ! iqibz==iqiA 
     end do qread_loop ! iqibz
#endif

   CASE (.FALSE.) ! Read the entire array.

     do iqibz=1,Hscr%nqibz
       do iomega=1,nomegaA
         do ipwe=1,Hscr%npwe
          read(unt,ERR=10) bufdc2d(1:Hscr%npwe,1)
          if (ipwe<=npweA) epsm1(1:npweA,ipwe,iomega,iqibz)=bufdc2d(1:npweA,1)
         end do
       end do
       ! Skip other frequencies
       do iomega=nomegaA+1,Hscr%nomega
         do ipwe=1,Hscr%npwe
           read(unt,ERR=10)
         end do
       end do
     end do !iqibz

   END SELECT
   !
   close(unt)                    
 end if ! (MPI-IO | FORTRAN)
 !
 ! === Free memory and close file ===
 if (allocated(bufdc2d))  then
   ABI_DEALLOCATE(bufdc2d)
 end if
 if (allocated(bufdc3d))  then
   ABI_DEALLOCATE(bufdc3d)
 end if

 call free_scrhdr(Hscr)
 DBG_EXIT("COLL")

 RETURN
 !
 ! === Something went wrong in Fortran IO part! ===
10 continue 
 MSG_ERROR("File seems to be corrupted.")

end subroutine read_screening
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/read_pole_screening
!! NAME
!! read_pole_screening
!!
!! FUNCTION
!! Read a pole-fit screening (\tilde epsilon^{-1})
!!
!! INPUTS
!!  accesswff=Integer flag defining the format of the output file. Available options:
!!    1--> Plain Fortran file 
!!    3--> ETSF format (TODO not yet coded)
!!  iqiA[optional]=Used if only a particular q-point is required. In this case iqiA define the index
!!   of the required q-point in the array qibz(3,Hscr%nqibz)
!!  nqibzA=number of asked q-points (used to dimension the output arrays). 
!!   Equal to Hscr%nqibz if the full matrix is required
!!  comm=MPI communicator.
!!  npweA=number of asked planewaves
!!  nomegaA=number of asked frequencies
!!
!! OUTPUT
!!  epsm1(npweA,npweA,ncoeffA,nqibzA) = \tilde\epsilon^{-1}(Ng,Ng,Nw,Nq)
!!
!! NOTES
!!  * If the epsilon matrix read is bigger than npweA x npweA, it will be truncated; 
!!    if it is smaller, an error will occur
!!  * If the number of frequencies asked for is smaller than that reported in the file, the matrix 
!!    will be truncated. If nomegaA > Hscr%nomega an error will occur 
!! 
!! PARENTS
!!      m_screening,mrgscr
!!
!! CHILDREN
!!      hdr_mpio_skip,mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine read_pole_screening(fname,npweA,nqibzA,ncoeffA,epsm1,accesswff,comm,&
& iqiA) ! Optional

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_pole_screening'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none
    
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,ncoeffA,npweA,nqibzA,comm
 integer,optional,intent(in) :: iqiA
 character(len=fnlen),intent(in) :: fname
!arrays
 real(gwpc),intent(inout) :: epsm1(npweA,npweA,ncoeffA,nqibzA)
  
!Local variables-------------------------------
!scalars
 integer,parameter :: master0=0
 integer :: ipwe,fform,icoeff,ios,iqibz,istat,unt
 integer :: rdwr,my_rank,nprocs
#if 0 && defined HAVE_FC_STREAM_IO
 integer :: recsize
 integer(kind=16) mypos,startpos,big_npwe
#endif
#ifdef HAVE_MPI_IO
 integer :: ierr,test_fform
 integer :: bsize_frm,mpi_type_frm
 integer :: mpi_fh,buf_dim,mat_ggw,mat_ggwq
 integer(MPI_OFFSET_KIND) :: offset,displ_wq
 !complex(dpc) :: ctmp
#endif
 character(len=500) :: msg
 logical :: read_qslice
 type(ScrHdr_type) :: Hscr
!arrays
#ifdef HAVE_MPI_IO
!integer :: statux(MPI_STATUS_SIZE)
 integer(MPI_OFFSET_KIND),allocatable :: offset_wq(:,:)
#if ! defined HAVE_GW_DPC
 real(dpc),allocatable :: bufd(:,:,:)
#endif
#endif
 real(dpc),allocatable :: bufd2d(:,:)
!#if ! defined HAVE_GW_DPC
 real(dpc),allocatable :: bufd3d(:,:,:)
!#endif

! *************************************************************************

 DBG_ENTER("COLL")

 my_rank = xcomm_rank(comm)
 nprocs  = xcomm_size(comm)
 !
 ! === Read the header of the file ===
 rdwr=1 
 !if (accesswff==IO_MODE_MPI .or. .TRUE.) then
 if (accesswff==IO_MODE_MPI) then
#if defined HAVE_MPI_IO
   bsize_frm    = xmpio_bsize_frm    ! bsize_frm= Byte length of the Fortran record marker.

mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

   if (my_rank==master0) then ! Master reads the header via Fortran IO then bcast the data.
     unt=get_unit()
     open(unit=unt,file=fname,status='old',form='unformatted',iostat=ios)
     msg = ' Opening file '//TRIM(fname)//' as old unformatted '
     ABI_CHECK(ios==0,msg)
     call scr_hdr_io(fform,rdwr,unt,xmpi_self,master0,IO_MODE_FORTRAN,Hscr)
     close(unt)
   end if

   call scrhdr_comm(HScr,master0,my_rank,comm)

   call xcast_mpi(fform,master0,comm,ierr)
   ! 
   ! Open the file with MPI-IO
   call MPI_FILE_OPEN(comm,fname, MPI_MODE_RDONLY, MPI_INFO_NULL,mpi_fh,ierr)
   ABI_CHECK_MPI(ierr,"MPI_FILE_OPEN") 
   !
   ! Retrieve the offset of the section immediately below the header.
   call scrhdr_mpio_skip(mpi_fh,test_fform,offset) 
   ABI_CHECK(test_fform==fform,"mismatch in fform!")
   !
   ! Offsets of the Fortran markers corresponding to the (w,q) slices.
   ABI_ALLOCATE(offset_wq,(HScr%nomega,HScr%nqibz))
   displ_wq = offset
   do iqibz=1,Hscr%nqibz
     do icoeff=1,Hscr%ncoeff
       offset_wq(icoeff,iqibz) = displ_wq 
       displ_wq = displ_wq + Hscr%npwe**2 * xmpi_bsize_dpc + Hscr%npwe * 2 * bsize_frm
     end do
   end do

   !do iqibz=1,Hscr%nqibz
   !  do iomega=1,Hscr%nomega
   !    displ_wq = offset_wq(iqibz,iomega)
   !    call xmpio_read_frm(mpi_fh,displ_wq,xmpio_at_all,bsize_frm,fmarker,ierr,advance=.FALSE.)

   !    call MPI_FILE_READ_AT(mpi_fh,displ_wq,ctmp,1,MPI_DOUBLE_COMPLEX,statux,ierr)
   !    write(std_out,*)"fmarker, ctmp ",fmarker,ctmp
   !  end do
   !end do

#else
   MSG_ERROR("MPI-IO support not enabled at configure-time")
#endif
 else 
   ! Plain Fortran IO, all nodes read.
   unt=get_unit()
   open(unit=unt,file=fname,status='old',form='unformatted',iostat=ios)
   msg = ' Opening file '//TRIM(fname)//' as old unformatted '
   ABI_CHECK(ios==0,msg)
   call scr_hdr_io(fform,rdwr,unt,comm,master0,accesswff,Hscr)
 end if
 !
 ! === Slice or full array? ===
 read_qslice=.FALSE.
 if (PRESENT(iqiA)) then 
   read_qslice=.TRUE.
   write(msg,'(a,i0,2a)')' reading q-slice corresponding to iq = ',iqiA,' from file: ',TRIM(fname)
   call wrtout(std_out,msg,'PERS')
   if (iqiA<=0.or.iqiA>Hscr%nqibz) then
     MSG_BUG('iqiA out of range')
   end if
 end if 
 !
 ! === Do some check ===
 if (Hscr%npwe>npweA) then
   write(msg,'(a,i8,2a,i8)')&
&    ' Total number of G-vectors reported on file = ',Hscr%npwe,ch10,&
&    ' Reading a smaller matrix of dimension      = ',npweA
   MSG_COMMENT(msg)
 end if

 if (npweA>Hscr%npwe) then
   write(msg,'(2(a,i0))')' Dimension of matrix = ',Hscr%npwe," requiring a too big matrix = ",npweA
   MSG_ERROR(msg)
 end if

 ABI_CHECK(nqibzA  <= Hscr%nqibz, ' Requiring too much q-points ')
 ABI_CHECK(ncoeffA <= Hscr%ncoeff,' Requiring too many coefficients ')

 !if (accesswff==IO_MODE_MPI .or. .TRUE.) then
 if (accesswff==IO_MODE_MPI) then
#ifdef HAVE_MPI_IO 
   if (read_qslice) then

     call xmpio_create_fsubarray_3D((/HScr%npwe,HScr%npwe,HScr%ncoeff/),(/npweA,npweA,ncoeffA/),MPI_DOUBLE_PRECISION,mat_ggw,ierr)
     ABI_CHECK_MPI(ierr,"xmpio_create_fsubarray_3D") 

     offset = offset_wq(1,iqiA) + bsize_frm

     call MPI_FILE_SET_VIEW(mpi_fh,offset,MPI_BYTE,mat_ggw,"native",MPI_INFO_NULL,ierr)
     ABI_CHECK_MPI(ierr,"SET_VIEW") 

     buf_dim = (npweA)**2 * ncoeffA 

#ifdef HAVE_GW_DPC
     ! Read in-place.
     call MPI_FILE_READ_ALL(mpi_fh,epsm1,buf_dim,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
     ABI_CHECK_MPI(ierr,"FILE_READ_ALL") 
#else
     ! Have to allocate workspace for dpc data.
     ABI_ALLOCATE(bufd3d,(npweA,npweA,ncoeffA))
     call MPI_FILE_READ_ALL(mpi_fh,bufd3d,buf_dim,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
     ABI_CHECK_MPI(ierr,"FILE_READ_ALL") 
     epsm1(:,:,:,1) = bufd3d
     ABI_DEALLOCATE(bufd3d)
#endif
     call MPI_TYPE_FREE(mat_ggw,ierr)

   else  ! Full matrix (G,G',w,q) is needed.

#ifdef HAVE_GW_DPC
     ! Can read all data at once.
     MSG_COMMENT("Using 4D Access")

     call xmpio_create_fsubarray_4D((/HScr%npwe,HScr%npwe,HScr%ncoeff,HScr%nqibz/),(/npweA,npweA,ncoeffA,HScr%nqibz/),&
&      MPI_DOUBLE_COMPLEX,mat_ggwq,ierr)
     ABI_CHECK_MPI(ierr,"xmpio_create_fsubarray_4D") 

     offset = offset_wq(1,1) + bsize_frm
     call MPI_FILE_SET_VIEW(mpi_fh,offset,MPI_BYTE,mat_ggwq,"native",MPI_INFO_NULL,ierr)
     ABI_CHECK_MPI(ierr,"SET_VIEW") 
                                                                                                       
     buf_dim = (npweA)**2 * ncoeffA * HScr%nqibz
     call MPI_FILE_READ_ALL(mpi_fh,epsm1,buf_dim,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
     ABI_CHECK_MPI(ierr,"FILE_READ_ALL") 

     call MPI_TYPE_FREE(mat_ggwq,ierr)

     !do iqibz=1,Hscr%nqibz
     !  do iomega=1,Hscr%nomega
     !    write(std_out,*)epsm1(1:2,1:2,iomega,iqibz)
     !  end do
     !end do
#else
     ! Have to allocate workspace for dpc data.
     ABI_ALLOCATE(bufd3d,(npweA,npweA,ncoeffA))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,'out of memory in bufd3d, change the code to loop over frequencies!')

     call xmpio_create_fsubarray_3D((/HScr%npwe,HScr%npwe,HScr%ncoeff/),(/npweA,npweA,ncoeffA/),MPI_DOUBLE_PRECISION,mat_ggw,ierr)
     ABI_CHECK_MPI(ierr,"xmpio_create_fsubarray_3D") 

     do iqibz=1,Hscr%nqibz
       offset = offset_wq(1,iqibz) + bsize_frm
       call MPI_FILE_SET_VIEW(mpi_fh,offset,MPI_BYTE,mat_ggw,"native",MPI_INFO_NULL,ierr)
       ABI_CHECK_MPI(ierr,"SET_VIEW") 

       buf_dim = (2*npweA)**2 * ncoeffA
       call MPI_FILE_READ_ALL(mpi_fh,bufd3d,buf_dim,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
       ABI_CHECK_MPI(ierr,"FILE_READ_ALL")

       epsm1(:,:,:,iqibz) = bufd3d
       !do iomega=1,Hscr%nomega
       !  write(std_out,*)epsm1(1:2,1:2,iomega,iqibz)
       !end do
     end do

     call MPI_TYPE_FREE(mat_ggw,ierr)
     ABI_DEALLOCATE(bufd3d)
#endif
   end if

   call MPI_FILE_CLOSE(mpi_fh,ierr)
   ABI_DEALLOCATE(offset_wq)
#endif

 else 
   !
   ! === Now read epsilon^-1 ===
   ! Allocate a single column to save memory.
   ! TODO re-merge the two cases.
   ABI_ALLOCATE(bufd2d,(Hscr%npwe,1))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,'out of memory in bufd2d')
   !
   ! === Two coding for different case just to keep it readable ===
   ! * Communication is done inside the loops to avoid problems with 
   !   the size of the MPI packet. Much slower but safer.
   SELECT CASE (read_qslice)

   CASE (.TRUE.) ! Read only a slice of the full array (useful if the entire array is huge).
     !if (dim_wings==1) STOP 'not implemented'
     !TODO this has to be done in a cleaner way.
#if 0 && defined HAVE_FC_STREAM_IO
!    Loop using stream I/O (Fortran 2003) to jump in file
!    First find the byte position after the header in the file by matching the first record
     read(unt,ERR=10) bufd2d(1:Hscr%npwe,1)
     close(unt)
     open(unit=unt,file=fname,status='old',form='unformatted',access='stream',iostat=ios)     
     recsize = sizeof(epsm1(1,1,1,1)) ! TODO fix this: Nonstandard fortran
     ! Find starting position
     startpos=0
     do
       startpos = startpos + 1
       read(unt,pos=startpos,ERR=10) epsm1(1:npweA,1,1,1)
       if (ALL(epsm1(1:npweA,1,1,1)==bufd2d(1:npweA,1))) EXIT
     end do
     ! Find stride
     read(unt,pos=startpos,ERR=10) epsm1(1:Hscr%npwe,1,1,1)
     inquire(unt,pos=mypos)
     write(std_out,'(2(a,I0))') ' Stride is: ',mypos-startpos,' should be: ',Hscr%npwe*recsize 
     write(std_out,'(a,I0)') ' Header ends at: ',startpos-1
     big_npwe = Hscr%npwe
     startpos = startpos + recsize*(big_npwe*big_npwe*Hscr%ncoeff*(iqiA-1)) &
&                        + 8*(big_npwe*Hscr%ncoeff*(iqiA-1))
     write(std_out,'(2(a,I0))') ' Entering read loop, recsize: ',recsize,' npweA: ',npweA
     write(std_out,'(2(a,I0))') '                     Hscr%npwe: ',Hscr%npwe,' startpos: ',startpos
     do icoeff=1,ncoeffA
       if (mod(icoeff,10)==0) write(std_out,'(2(a,I0))') ' Reading coeff: ',icoeff
       do ipwe=1,npweA
         ! Fortran record are written with 4B header and 4B footer (64 bit machines)
         ! Stride is the record length*Hscr%npwe + 8 bytes at the end
         mypos = startpos + recsize*(big_npwe*big_npwe*(icoeff-1) + big_npwe*(ipwe-1)) &
&                         + 8*big_npwe*(icoeff-1) + 8*(ipwe-1)
         read(unt,pos=mypos,ERR=10) epsm1(1:npweA,ipwe,icoeff,1)
       end do
     end do
#else 
     qread_loop: &
&    do iqibz=1,Hscr%nqibz
       if (iqibz==iqiA) then 
         do icoeff=1,ncoeffA
           do ipwe=1,Hscr%npwe
             read(unt,ERR=10) bufd2d(1:Hscr%npwe,1)
             if (ipwe<=npweA) epsm1(1:npweA,ipwe,icoeff,1)=bufd2d(1:npweA,1)
           end do
         end do
         EXIT qread_loop ! Got data. Do not need to read file till the end.
       else ! Skip other q-points i.e bufd2d(1:Hscr%npwe,1:Hscr%npwe)
         do icoeff=1,Hscr%ncoeff
           do ipwe=1,Hscr%npwe
            read(unt,ERR=10)
           end do
         end do
       end if ! iqibz==iqiA 
     end do qread_loop ! iqibz
#endif

   CASE (.FALSE.) ! Read the entire array.

     do iqibz=1,Hscr%nqibz
       do icoeff=1,ncoeffA
         do ipwe=1,Hscr%npwe
          read(unt,ERR=10) bufd2d(1:Hscr%npwe,1)
          if (ipwe<=npweA) epsm1(1:npweA,ipwe,icoeff,iqibz)=bufd2d(1:npweA,1)
         end do
       end do
       ! Skip other coefficients
       do icoeff=ncoeffA+1,Hscr%ncoeff
         do ipwe=1,Hscr%npwe
           read(unt,ERR=10)
         end do
       end do
     end do !iqibz

   END SELECT
   !
   close(unt)                    
 end if ! (MPI-IO | FORTRAN)
 !
 ! === Free memory and close file ===
 if (allocated(bufd2d))  then
   ABI_DEALLOCATE(bufd2d)
 end if
 if (allocated(bufd3d))  then
   ABI_DEALLOCATE(bufd3d)
 end if

 call free_scrhdr(Hscr)
 DBG_EXIT("COLL")

 RETURN
 !
 ! === Something went wrong in Fortran IO part! ===
10 continue 
 MSG_ERROR("File seems to be corrupted.")

end subroutine read_pole_screening
!!***

!!****f* m_header/scrhdr_mpio_skip
!! NAME
!!  scrhdr_mio_skip
!!
!! FUNCTION
!!   Skip the header of the (SCR|SUSC) file in MPI-IO mode. This routine uses local MPI-IO calls hence 
!!   it can be safely called by master node only. Note however that in this case the
!!   offset has to be communicated to the other nodes.
!!
!! INPUTS
!!  mpio_fh=MPI-IO file handler
!!  fmarker_bsize   = Byte length of Fortran record marker.
!!  fmarker_mpi_type= MPI type of the Fortran record marker
!!
!! OUTPUT
!!  fform=kind of the array in the file
!!  offset=The offset of the Fortran record located immediately below the Abinit header.
!!
!! SOURCE

subroutine scrhdr_mpio_skip(mpio_fh,fform,offset) 

 use defs_basis

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scrhdr_mpio_skip'
!End of the abilint section

 integer,intent(in) :: mpio_fh
 integer,intent(out) :: fform
 integer(kind=XMPI_OFFSET_KIND),intent(out) :: offset

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,mpi_type_frm
#ifdef HAVE_MPI_IO
 integer :: ierr,isk,nqlwl
 character(len=500) :: msg
!arrays
 integer(kind=MPI_OFFSET_KIND) :: fmarker,positloc
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *************************************************************************

 offset = 0

 bsize_frm    = xmpio_bsize_frm    ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

 call hdr_mpio_skip(mpio_fh,fform,offset) 

#ifdef HAVE_MPI_IO
 select case (fform)
 case (1002)
   ! Old epsilon^-1 file used in version <5.6
   do isk=1,5
     call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)
   end do

 case (1102)
   ! File format for epsilon^-1, espilon, chi0 
   call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)
   call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)

   ! read nqlwl
   !read(unt,ERR=10)Hscr%npwe,Hscr%npwwfn_used,Hscr%nbnds_used,Hscr%nqibz,Hscr%nomega,Hscr%nqlwl
   positloc  = offset + bsize_frm + 5*xmpi_bsize_int
   call MPI_FILE_READ_AT(mpio_fh,positloc,nqlwl,1,MPI_INTEGER,statux,ierr)

   do isk=1,3
     call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)
   end do

   if (nqlwl>0) then ! skip wings
     do isk=1,3
       call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)
     end do
   end if
                                                                                                
 case default
   write(msg,'(a,i0)')' scrhdr_mpio_skip: Wrong fform read = ',fform
   MSG_BUG(msg)
 end select

#else
 MSG_ERROR("scrhdr_mpio_skip cannot be used when MPI-IO is not enabled")
#endif

end subroutine scrhdr_mpio_skip
!!***

!----------------------------------------------------------------------

END MODULE m_io_screening
!!***
