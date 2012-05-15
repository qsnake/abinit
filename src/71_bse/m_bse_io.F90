!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_bse_io
!! NAME
!! m_bse_io
!!
!! FUNCTION
!!  This module provides routines to read the Bethe-Salpeter Hamiltonian from file
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_bse_io

 use m_profiling

 use defs_basis
 use m_xmpi
 use m_errors 

 use defs_abitypes,    only : Hdr_type
 use m_fstrings,       only : toupper
 use m_io_tools,       only : get_unit
 use m_numeric_tools,  only : arth
 use m_special_funcs,  only : gaussian
 use m_header,         only : hdr_clean, hdr_mpio_skip
 use m_bs_defs,        only : excparam

#if defined HAVE_MPI2
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 private
!!***

!----------------------------------------------------------------------
!!***

 public  :: exc_read_rblock_fio      ! Reads the entired resonant sub-block from file using Fortran IO.
 public  :: exc_read_rcblock         ! Reads a distributed sub-block of the excitonic Hamiltonian from file.
 public  :: exc_fullh_from_blocks    ! Initialized the specified sub-blocks of the *full* matrix (reso+anti-reso) from file.
 public  :: rrs_of_glob              ! [+1,-1,0] if (row_glob,col_glob) belongs to the [ resonant, anti-resonant, (anti)coupling block ]
 public  :: ccs_of_glob              ! [+1,-1,0] if (row_glob,col_glob) belongs to the [ coupling, anti-coupling, (anti)resonant block ]
 public  :: offset_in_file           ! Function used to describe the way the Hamiltonian is stored on disk.
 public  :: exc_write_bshdr          ! Writes the Header of the (BSR|BSC) files storing the excitonic Hamiltonian.
 public  :: exc_read_bshdr           ! Reads the Header of the (BSR|BSC) files.
 public  :: exc_skip_bshdr           ! Skip the Header of the (BSR|BSC) files. Fortran version.
 public  :: exc_skip_bshdr_mpio      ! Skip the Header of the (BSR|BSC) files. MPI-IO  version.
 public  :: exc_read_eigen           ! Read selected energies and eigenvectors from the BSEIG file.
 public  :: exc_amplitude            ! Calculate the amplitude function F(w) = \sum_t |<t|exc_vec>|^2 \delta(ww- ene_t) where t is the eh transition.

CONTAINS  !====================================================================


!----------------------------------------------------------------------

!!****f* m_bse_io/exc_write_bshdr
!! NAME
!!  exc_write_bshdr
!!
!! FUNCTION
!!   Writes the header of the (BSR|BSC) files storing the excitonic Hamiltonian.
!!
!! INPUTS
!!  funt=Fortran unit number.
!!  Bsp<excparam>=Structure storing the parameters of the run.
!!  Hdr<hdr_type>=The abinit header.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      exc_build_block
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_write_bshdr(funt,Bsp,Hdr)
                                                                                                      
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_write_bshdr'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none
                                                                                                      
 !Arguments ------------------------------------
 integer,intent(in) :: funt
 type(excparam),intent(in) :: BSp
 type(hdr_type),intent(inout) :: Hdr

!Local variables ------------------------------
!scalars
 integer,parameter :: rdwr2=2 
 integer :: fform_1002
 ! *************************************************************************

 fform_1002=1002 ! TODO: change setup_bse so that Hdr_bse reflects the parameters of the run.
 call hdr_io_int(fform_1002,Hdr,rdwr2,funt)
 write(funt) BSp%nreh,BSp%nkbz

end subroutine exc_write_bshdr
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/exc_read_bshdr
!! NAME
!!  exc_read_bshdr
!!
!! FUNCTION
!!  Reads the header of the (BSR|BSC) files storing the excitonic Hamiltonian.
!!  and performs basilar consistency checks.
!!
!! INPUTS
!!  funt=Unit number.
!!  Bsp<excparam>=Structure storing the parameters of the run.
!!  Hdr<hdr_type>=The abinit header.
!!
!! OUTPUT
!!  fform=Integer defining the file format.
!!  ierr=Status error.
!!
!! PARENTS
!!      exc_diago,m_bse_io
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_read_bshdr(funt,Bsp,fform,ierr)
                                                                                                      
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_read_bshdr'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none
                                                                                                      
 !Arguments ------------------------------------
 integer,intent(in) :: funt
 integer,intent(out) :: fform,ierr
 type(excparam),intent(in) :: BSp

!Local variables ------------------------------
!scalars
 integer,parameter :: rdwr1=1
 integer :: nkbz_read
 type(hdr_type) :: Hdr
!arrays
 integer :: nreh_read(SIZE(BSp%nreh))

 ! *************************************************************************

 ierr=0

 ! Read the header and perform consistency checks.
 call hdr_io_int(fform,Hdr,rdwr1,funt)

 read(funt,ERR=10) nreh_read, nkbz_read

 call hdr_clean(Hdr)

 ABI_CHECK(ANY(nreh_read==BSp%nreh),"Wrong number of e-h transitions")

 RETURN

10 ierr=1

end subroutine exc_read_bshdr
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/exc_skip_bshdr
!! NAME
!!  exc_skip_bshdr
!!
!! FUNCTION
!!   Skip the header of the (BSR|BSC) files storing the excitonic Hamiltonian. Fortran version.
!!
!! INPUTS
!!  funt=Unit number.
!!
!! OUTPUT
!!  ierr=Status error.
!!
!! SIDE EFFECTS
!!  Skip the header.
!!
!! PARENTS
!!      exc_build_block
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_skip_bshdr(funt,ierr)
                                                                                                      
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_skip_bshdr'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none
                                                                                                      
 !Arguments ------------------------------------
 integer,intent(in) :: funt
 integer,intent(out) :: ierr

 ! *************************************************************************

 call hdr_skip_int(funt,ierr)
 if (ierr/=0) RETURN
 read(funt)

end subroutine exc_skip_bshdr
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/exc_skip_bshdr_mpio
!! NAME
!!  exc_skip_bshdr_mpio
!!
!! FUNCTION
!!   Skip the header of the (BSR|BSC) files storing the excitonic Hamiltonian. MPI-IO version.
!!
!! INPUTS
!!  mpi_fh=MPI-IO file handler.
!!  at_option
!!
!! SIDE EFFECTS
!!  ehdr_offset
!!
!! PARENTS
!!      exc_build_block,exc_diago,m_bse_io
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_skip_bshdr_mpio(mpi_fh,at_option,ehdr_offset)
                                                                                                      
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_skip_bshdr_mpio'
!End of the abilint section

 implicit none
                                                                                                      
 !Arguments ------------------------------------
 integer,intent(in) :: mpi_fh,at_option
 integer(XMPI_OFFSET_KIND),intent(inout) :: ehdr_offset

!Local variables ------------------------------
!scalars
 integer :: fform,ierr
#ifdef HAVE_MPI_IO   
 integer(XMPI_OFFSET_KIND) :: fmarker
#endif

 ! *************************************************************************

 call hdr_mpio_skip(mpi_fh,fform,ehdr_offset) 

#ifdef HAVE_MPI_IO
 call xmpio_read_frm(mpi_fh,ehdr_offset,at_option,fmarker,ierr)
 !write(std_out,*)"fmarker last record ",fmarker
#else
 MSG_ERROR("You should not be here")
#endif

end subroutine exc_skip_bshdr_mpio
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/exc_read_eigen
!! NAME
!!  exc_read_eigen
!!
!! FUNCTION
!!  Read selected energies and eigenvectors from the BSEIG file.
!!
!! INPUTS
!!  eig_fname=The name of the file storing the excitonic eigenvectors.
!!  hsize=Size of the Hamiltonian.
!!  nvec=Number of excitonic states to analyze.
!!  vec_idx(nvec)=List with the indeces of the excitonic states sorted in ascending order.
!!  [Bsp]<excparam>=Structure storing the parameters of the run. If present the 
!!    routine will perform additional consistency checks to make sure that 
!!    the content of the file is consistent with the present run.
!!
!! OUTPUT
!!  [ene_list(nvec)]=Excitonic energies
!!  vec_list(hsize,nvec)=Excitonic eigenvectors.
!!
!! PARENTS
!!      exc_plot,m_bse_io
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_read_eigen(eig_fname,hsize,nvec,vec_idx,vec_list,ene_list,Bsp)
                                                                                                      
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_read_eigen'
!End of the abilint section

 implicit none
                                                                                                      
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nvec,hsize
 character(len=*),intent(in) :: eig_fname
 type(excparam),optional,intent(in) :: BSp
! arrays
 integer,intent(in) :: vec_idx(nvec)
 real(dp),optional,intent(out) :: ene_list(nvec)
 complex(dpc),intent(out) :: vec_list(hsize,nvec)

!Local variables ------------------------------
!scalars
 integer :: eig_unt,ios,hsize_read,neig_read,ll,vec
 character(len=500) :: msg
!arrays
 !real(dp),allocatable :: exc_ene(:)
 complex(dpc),allocatable :: exc_ene_cplx(:)

 ! *************************************************************************

 ABI_UNUSED(BSp%nline)

 eig_unt = get_unit() 
 open(unit=eig_unt,file=eig_fname,form='unformatted',status='old',iostat=ios)
 msg="Opening file "//TRIM(eig_fname)//" as old unformatted"
 ABI_CHECK(ios==0,msg) 

 read(eig_unt)hsize_read,neig_read
 !write(std_out,*)hsize_read, neig_read

 if (hsize_read/=hsize) then
   write(msg,'(a,2(1x,i0))')" hsize_read/=hsize: ",hsize_read,hsize
   MSG_ERROR(msg)
 end if
 !
 ! Read eigenvalues, ignore possibly small imaginary part.
 ABI_ALLOCATE(exc_ene_cplx,(neig_read))
 read(eig_unt) exc_ene_cplx 

 if (PRESENT(ene_list)) then
   do vec=1,nvec
     ll = vec_idx(vec)
     ene_list(vec) = DBLE(exc_ene_cplx(ll))
   end do
 end if
 ABI_DEALLOCATE(exc_ene_cplx)

 vec=1
 do ll=1,neig_read ! Read the selected excitons.
   if (ll==vec_idx(vec))  then
     read(eig_unt) vec_list(:,vec)
     if (vec==nvec) EXIT
     vec=vec+1
   else 
     read(eig_unt)
   end if
 end do

 close(eig_unt)

 if (vec/=nvec) then
   write(msg,'(a,2(1x,i0))')" vec_idx is wrong, vec/=nvec ",vec,nvec+1
   MSG_ERROR(msg)
 end if

end subroutine exc_read_eigen
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/exc_read_rcblock
!! NAME
!! exc_read_rcblock
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file 
!!
!! INPUTS
!!  fname=File name.
!!  diago_is_real=.TRUE. if diagonal elements are real (used only if is_resonant==.TRUE.)
!!  nresh(nsppol)=Number of resonant transition for the two spins.
!!  is_resonant=Set to .TRUE. if the block is resonant. 
!!  hsize=Dimension of the block.
!!  nsppol=2 for spin polarized systems. 1 otherwise.
!!  my_t1,my_t2=The first and the last colums of the matrix treated by this node. 
!!  use_mpio=.TRUE. is MPI-IO routines are used.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  hmat(hsize,my_t1:my_t2)=The block read from file fname.
!!
!! PARENTS
!!      exc_iterative_diago,haydock
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_read_rcblock(fname,Bsp,is_resonant,diago_is_real,nsppol,nreh,hsize,my_t1,my_t2,hmat,use_mpio,comm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_read_rcblock'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,hsize,my_t1,my_t2,nsppol
 logical,intent(in) :: is_resonant,use_mpio,diago_is_real
 character(len=*),intent(in) :: fname
 type(excparam),intent(in) :: Bsp
!arrays
 integer,intent(in) :: nreh(nsppol)
 complex(dpc),intent(out) :: hmat(hsize,my_t1:my_t2)

!Local variables ------------------------------
!scalars
 integer :: it,itp,funit,ios,nproc,my_rank,neh1,neh2
 integer :: fform,master,my_nt
 integer :: row,col,block,spad,spin_dim,ierr,size_exp
 character(len=500) :: msg
 !type(Hdr_type) :: bse_Hdr
!arrays
 complex(dpc),allocatable :: buffer_dpc(:)
 logical :: have_row,have_col
#ifdef HAVE_MPI_IO
 integer :: mpi_err,mpi_fh,ham_type,my_nel,old_type,etype,offset_err,amode
 integer :: irec,nrec !,recblock ! ,ncount
 integer(XMPI_OFFSET_KIND) :: ehdr_offset,my_offset,my_offpad,fsize
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecord(:)
 integer :: glob_sizes(2),my_cols(2) !,my_starts(2),my_ends(2)
 integer :: block_sizes(2,3)
 integer :: status(MPI_STATUS_SIZE)
#endif

!************************************************************************

 nproc  = xcomm_size(comm)
 my_rank= xcomm_rank(comm)
 master=0

 neh1 = nreh(1)
 neh2 = neh1; if (nsppol==2) neh2=nreh(2)

 size_exp=neh1; if (nsppol==2) size_exp=SUM(nreh)

 ABI_CHECK(hsize==size_exp,"Wrong hsize")
 if (neh1/=neh2) then
   write(msg,'(a)')" BSE code does not support different number of transitions for the two spin channels"
   MSG_ERROR(msg)
 end if

 my_nt = my_t2-my_t1+1

!BEGINDEBUG
 hmat = HUGE(zero)
!ENDDEBUG

 if (.not.use_mpio) then

   if (is_resonant) then
     call wrtout(std_out,". Reading resonant block from file: "//TRIM(fname)//" using Fortran-IO","COLL") 
   else
     call wrtout(std_out,". Reading coupling block from file: "//TRIM(fname)//" using Fortran-IO","COLL") 
   end if

   funit = get_unit()
   open(funit,file=fname,form='unformatted',status="old",iostat=ios)
   ABI_CHECK(ios==0," Opening file: "//TRIM(fname))
   !
   ! Read the header and perform consistency checks.
   call exc_read_bshdr(funit,Bsp,fform,ierr)
   !
   ! Construct full excitonic Hamiltonian using symmetries.
   if (nsppol==1) then
     ! 
     ABI_ALLOCATE(buffer_dpc,(neh1))
     do itp=1,hsize
       read(funit) buffer_dpc(1:itp)
       !
       ! Fill the upper triangle if I have this column.
       if (itp>=my_t1 .and. itp<=my_t2) then
         do it=1,itp
           hmat(it,itp) = buffer_dpc(it)
         end do
         ! Force the diagonal to be real.
         if (is_resonant .and.diago_is_real) hmat(itp,itp) = DBLE(hmat(itp,itp))
       end if
       !
       ! Reconstruct the rows below the diagonal (diagonal part is not touched here).
       if (is_resonant) then ! Use Hermiticity.
         do it=1,itp-1 
           if (it>=my_t1 .and. it<=my_t2) hmat(itp,it) = CONJG(buffer_dpc(it))
         end do
       else  ! Coupling is symmetric.
         do it=1,itp-1 
           if (it>=my_t1 .and. it<=my_t2) hmat(itp,it) = buffer_dpc(it)
         end do
       end if
       !
     end do ! itp
     ABI_DEALLOCATE(buffer_dpc)
     !
   else 
     !
     ! The file contains 
     ! A) The up-up and the down-down block in packed form 
     ! (only the upper triangle is stored since the blocks are Hermitian)
     ! B) The entire up-down exchange block (no symmetry here)
     !
     ! A) Construct resonant blocks from the upper triangles stored on file.
     ! FIXME this part wont work if we have a different number of e-h pairs
     if (.not.is_resonant) then
       MSG_ERROR("exc_read_rcblock does not support coupling.")
     end if
     ! It should be checked.
     spin_dim=neh1
     do block=1,2
       ABI_ALLOCATE(buffer_dpc,(neh1))
       if (block==1) spad=0
       if (block==2) spad=neh1
       do itp=1,spin_dim
         !
         read(funit) buffer_dpc(1:itp)
         !
         ! Fill the upper triangle if I have this column
         col = itp+spad
         if (col>=my_t1 .and. col<=my_t2) then
           do it=1,itp
             row = it + spad
             hmat(row,col) = buffer_dpc(it)
           end do
           ! Force the diagonal to be real.
           if (is_resonant .and.diago_is_real) hmat(col,col) = DBLE(hmat(col,col)) 
         end if
         !
         ! Reconstruct the rows below the diagonal (diagonal part is not touched here).
         row = itp + spad
         if (is_resonant) then ! Hermitian
           do it=1,itp-1 
             col = it + spad
             if (col>=my_t1 .and. col<=my_t2) hmat(row,col) = CONJG(buffer_dpc(it))
           end do
         else  ! Coupling is symmetric.
           do it=1,itp-1 
             col = it + spad
             if (col>=my_t1 .and. col<=my_t2) hmat(row,col) = buffer_dpc(it)
           end do
         end if
         !
       end do ! itp
       ABI_DEALLOCATE(buffer_dpc)
     end do ! block
     !
     ! B) Kx_{down up} = Kx_{up down}^H.
     ! FIXME this part wont work if we have a different number of e-h pairs
     spad=neh1
     spin_dim=neh1
     ABI_ALLOCATE(buffer_dpc,(neh1))
     do itp=1,spin_dim
       read(funit) buffer_dpc(1:spin_dim)
       have_col = (spad+itp>=my_t1 .and. spad+itp<=my_t2)
       if (have_col) hmat(1:spin_dim,spad+itp) = buffer_dpc(1:spin_dim)
       ! Construct and store the lower block 
       if (is_resonant) then ! Hermitian
         do it=1,spin_dim
           have_row = (it>=my_t1 .and. it<=my_t2)
           if (have_row) hmat(spad+itp,it) = CONJG(buffer_dpc(it))
         end do
       else ! Symmetric
         do it=1,spin_dim
           have_row = (it>=my_t1 .and. it<=my_t2)
           if (have_row) hmat(spad+itp,it) = buffer_dpc(it)
         end do
       end if
     end do
     ABI_DEALLOCATE(buffer_dpc)
   end if

   close(funit)

 else 
#ifdef HAVE_MPI_IO
   if (is_resonant) then
     call wrtout(std_out,". Reading resonant block from file "//TRIM(fname)//" using MPI-IO","COLL") 
   else
     call wrtout(std_out,". Reading coupling block from file "//TRIM(fname)//" using MPI-IO","COLL") 
   end if

   amode=MPI_MODE_RDONLY
   call MPI_FILE_OPEN(comm,fname,amode,MPI_INFO_NULL,mpi_fh,mpi_err)
   msg = " FILE_OPEN "//TRIM(fname)
   ABI_CHECK_MPI(mpi_err,msg)

   call MPI_FILE_GET_SIZE(mpi_fh,fsize,mpi_err)
   write(std_out,*)" file size is ",fsize
   !
   ! Skip the header and find the offset for reading the matrix.
   call exc_skip_bshdr_mpio(mpi_fh,xmpio_at_all,ehdr_offset)
   !
   ! Read my columns from file.
   old_type=MPI_DOUBLE_COMPLEX
   glob_sizes = (/hsize,hsize/); my_cols=(/my_t1,my_t2/)

   if (nsppol==1) then
     call xmpio_create_coldistr_from_fpacked(glob_sizes,my_cols,old_type,ham_type,my_offpad,offset_err)
   else
     MSG_WARNING("nsppol==2 => calling fp3blocks")
     write(std_out,*)"neh, hsize",neh1,neh2,hsize

#if 1
     nrec=neh1+2*neh2
     ABI_ALLOCATE(bsize_frecord,(nrec))
     bsize_frecord(1:neh1)           = (/(irec*xmpi_bsize_dpc, irec=1,neh1)/) 
     bsize_frecord(neh1+1:neh1+neh2) = (/(irec*xmpi_bsize_dpc, irec=1,neh2)/)
     bsize_frecord(neh1+neh2+1:)     = neh1*xmpi_bsize_dpc
     call xmpio_check_frmarkers(mpi_fh,ehdr_offset,xmpio_at_all,nrec,bsize_frecord,ierr)
     ABI_CHECK(ierr==0,"Error in Fortran markers")
     ABI_DEALLOCATE(bsize_frecord)
     MSG_COMMENT("Marker check ok")
     call xbarrier_mpi(comm)
#endif

     block_sizes(:,1) = (/neh1,neh1/)
     block_sizes(:,2) = (/neh2,neh2/)
     block_sizes(:,3) = (/neh1,neh2/)
     MSG_ERROR("fp3blocks is buggy")
     call xmpio_create_coldistr_from_fp3blocks(glob_sizes,block_sizes,my_cols,old_type,ham_type,my_offpad,offset_err)
   end if

   if (offset_err/=0) then 
     write(msg,"(3a)")&
&      " Global position index cannot be stored in a standard Fortran integer ",ch10,&
&      " Excitonic matrix cannot be read with a single MPI-IO call."
     MSG_ERROR(msg)
   end if
   !
   ! The offset used for reading.
   my_offset = ehdr_offset + my_offpad 
   write(std_out,*)"my_offset= ",my_offset

   etype=MPI_BYTE
   call MPI_FILE_SET_VIEW(mpi_fh, my_offset, etype, ham_type, 'native', MPI_INFO_NULL, mpi_err)
   ABI_CHECK_MPI(mpi_err,"SET_VIEW")
   !
   ! Release the MPI filetype.
   call MPI_TYPE_FREE(ham_type,mpi_err)
   ABI_CHECK_MPI(mpi_err,"TYPE_FREE")

   my_nel = my_nt*hsize
   call MPI_FILE_READ_ALL(mpi_fh, hmat, my_nel, MPI_DOUBLE_COMPLEX, status, mpi_err)
   ABI_CHECK_MPI(mpi_err,"READ_ALL")

   !call MPI_GET_COUNT(status, MPI_DOUBLE_COMPLEX, ncount, mpi_err)
   !write(std_out,*)"count, my_nel ",ncount,my_nel
   !ABI_CHECK_MPI(mpi_err,"READ_ALL")
   !
   ! Close the file.
   call MPI_FILE_CLOSE(mpi_fh, mpi_err)               
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")
   !
   ! Use the symmetries of the block to reconstruct the local buffer.
   ! Coupling does not require in-place symmetrization since it is symmetric.
   if (is_resonant) then
     do itp=my_t1,my_t2
       if (itp+1<=hsize) hmat(itp+1:,itp) = DCONJG(hmat(itp+1:,itp)) ! Lower triangle using Hermiticity.
       if (diago_is_real) hmat(itp,itp) = DBLE(hmat(itp,itp))        ! The diagonal is forced to be real when energies are real.
     end do
   end if

   call xbarrier_mpi(comm)
#else 
   MSG_ERROR("MPI-IO support not enabled")
#endif
 end if

!BEGINDEBUG
  if ( ANY(hmat==HUGE(zero)) ) then
    write(std_out,*)"COUNT",COUNT(hmat==HUGE(zero))," hsize= ",hsize
    MSG_ERROR("Something wrong in the reading")
  end if
!ENDDEBUG

 call xbarrier_mpi(comm)

end subroutine exc_read_rcblock
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/exc_fullh_from_blocks
!! NAME
!!  exc_fullh_from_blocks 
!!
!! FUNCTION
!!   Construct the matrix F H
!!
!! INPUTS
!!  funt
!!  nsppol
!!  block_type
!!    "Resonant"
!!    "Coupling"
!!  row_sign
!!   -1 to read ( R   C )
!!              (-C* -R*)
!!
!!   +1 to read ( R   C )
!!              ( C*  R*)
!!  diago_is_real=Used when block_type=resonat to specify wheter the diagonal matrix elements are
!!  real or complex (when QP linewidth are included)
!!  nreh(nsppol)
!!  exc_size=Size of the full excitonic Hamiltonian.
!!
!! SIDE EFFECTS 
!!  exc_ham(exc_size,exc_size)
!!
!! NOTES  
!! 
!! PARENTS
!!      exc_diago
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_fullh_from_blocks(funt,block_type,nsppol,row_sign,diago_is_real,nreh,exc_size,exc_ham)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_fullh_from_blocks'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: funt,exc_size,nsppol,row_sign
 logical,intent(in) :: diago_is_real
 character(len=*),intent(in) :: block_type
!arrays
 integer,intent(in) :: nreh(nsppol)
 complex(dpc),intent(inout) :: exc_ham(exc_size,exc_size)

!Local variables-------------------------------
!scalars
 integer :: it,itp,istat,szbuf,neh,pad_c1,pad_r1,spin_dim,spad_r,spad_c
 integer :: block,spad,row1,col1,row2,col2,spin_stride
 complex(dpc) :: cttp
 character(len=500) :: msg
!arrays
 complex(dpc),allocatable :: cbuff_dpc(:)

! *********************************************************************

 szbuf=exc_size ! FIXME oversized!
 neh = nreh(1)

 if (nsppol==2) then
   MSG_WARNING("nsppol==2 is very experimental")
 end if
 if (ANY(nreh(1)/=nreh)) then
   MSG_ERROR(" different nreh are not supported")
 end if

 ABI_ALLOCATE(cbuff_dpc,(exc_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in cbuff_dpc")
 !
 ! The two cases nsppol==1,2 can be merged but we keep them 
 ! separated to keep to code readable.

 SELECT CASE (toupper(block_type))
 CASE ("RESONANT")

   if (nsppol==1) then 
     !
     ! Construct resonant block from the upper triangle stored on file.
     neh = nreh(1)
     do itp=1,neh
       read(funt) cbuff_dpc(1:itp)
       do it=1,itp-1 ! The diagonal is treated below.
         cttp = cbuff_dpc(it) 
         exc_ham(it     ,itp)     =                cttp   ! R             
         exc_ham(itp    ,it )     =          CONJG(cttp)  ! row_sign R*
         exc_ham(neh+it ,neh+itp) = row_sign*CONJG(cttp)
         exc_ham(neh+itp,neh+it ) = row_sign*cttp
       end do
       if (diago_is_real) then
         exc_ham(itp    ,itp)     =          DBLE(cbuff_dpc(itp))  
         exc_ham(neh+itp,neh+itp) = row_sign*DBLE(cbuff_dpc(itp))
       else
         exc_ham(itp,itp)         =           cbuff_dpc(itp)
         exc_ham(neh+itp,neh+itp) = row_sign*(CONJG(cbuff_dpc(itp)))
       end if
     end do
     !
     !
   else 
     ! FIXME this part wont work if we have a different number of e-h pairs
     ABI_CHECK(ALL(nreh==nreh(1)),"Different number of transitions")
     ! The file contains 
     ! A) The up-up and the down-down block in packed form 
     ! (only the upper triangle is stored since these blocks are Hermitian)
     ! B) The entire up-down exchange block (no symmetry here)
     !
     ! The resonant block is given by:
     !     |  (v'c' up)    | (v'c' dwn) | 
     !     ------------------------------           where v_{-+} = v_{+-}^H when the momentum of the photon is neglected. 
     !     | [T-W+v]++     |      v+-   | (vc up)   Note that v_{+-} is not Hermitian due to the presence of different spins.
     ! R = ------------------------------           Actually it reduces to a Hermitian matrix when the system is not spin polarized.
     !     |     v-+       | [T-W+v]--  | (vc dwn)  [T-W+v] is Hermitian provided the the QP energies are purely real.
     !     ------------------------------   
     !
     ! *) Fill the diagonal blocks.
     !    only the upper triangle is stored on file.
     !    row1,col1 refer to the resonant block.
     !    row2,col2 refer to the anti-resonant block.
     do block=1,2
       if (block==1) then 
         spad=0
         spin_stride=SUM(nreh)
       end if
       if (block==2) then 
         spad=nreh(1)
         spin_stride=2*nreh(1) + nreh(2)
       end if
       do itp=1,nreh(block)
         read(funt) cbuff_dpc(1:itp)
         col1 = itp+spad
         col2 = itp+spin_stride
         do it=1,itp-1
           cttp = cbuff_dpc(it)
           row1 = it + spad
           row2 = it + spin_stride
           exc_ham(row1,col1) =                cttp    ! [T-W+v]
           exc_ham(col1,row1) =          CONJG(cttp)
           exc_ham(row2,col2) = row_sign*CONJG(cttp)   ! row_sign [T-W+v]*
           exc_ham(col2,row2) = row_sign*cttp
         end do
         if (diago_is_real) then
           exc_ham(col1,col1) =          DBLE(cbuff_dpc(itp))
           exc_ham(col2,col2) = row_sign*DBLE(cbuff_dpc(itp))
         else 
           exc_ham(col1,col1) = cbuff_dpc(itp)
           exc_ham(col2,col2) = row_sign*CONJG(cbuff_dpc(itp))
         end if
       end do
     end do
     !
     ! Read v+- and reconstruct resonant and anti-resonat blocks.
     ! TODO recheck this 
     pad_r1=SUM(nreh)
     pad_c1=2*nreh(1) + nreh(2)

     spin_dim=nreh(1)
     do itp=1,spin_dim
       read(funt) cbuff_dpc(1:spin_dim)
       exc_ham(1:spin_dim,nreh(1)+itp) = cbuff_dpc(1:spin_dim)                              ! upper reso
       exc_ham(1+pad_r1:pad_r1+spin_dim,pad_c1+itp) = row_sign*CONJG(cbuff_dpc(1:spin_dim)) ! upper anti-reso
       col1 = itp+nreh(1)
       col2 = itp+(2*nreh(1) + nreh(2))
       do it=1,spin_dim
         cttp = cbuff_dpc(it) 
         row1 = it
         exc_ham(col1,row1) =    CONJG(cttp)  ! lower reso.
         row2 = it + SUM(nreh)
         exc_ham(col2,row2) = row_sign*cttp   ! lower anti-reso.
       end do
     end do
   end if

 CASE ("COUPLING")
   !
   if (nsppol==1) then
     !
     do itp=1,neh
       read(funt) cbuff_dpc(1:itp)
       do it=1,itp
         cttp = cbuff_dpc(it)
         exc_ham(it    ,neh+itp) = cttp
         exc_ham(itp   ,neh+it ) = cttp
         exc_ham(neh+it ,itp   ) = row_sign*CONJG(cttp)
         exc_ham(neh+itp,it    ) = row_sign*CONJG(cttp)
       end do
     end do
     !
   else
     !  The coupling block is given by:
     !      |  (c'v' up)   |    (c'v dwn)     | 
     !      -----------------------------------           where v_{-+} = v_{+-}^t when the momentum of the photon is neglected. 
     !      | [-W+v]++     |      v+-         | (vc up)   The entire matrix v_{+-} is stored on file.
     !  C = -----------------------------------            
     !      |     v-+      |    [-W+v]--      | (vc dwn)
     !      -----------------------------------   
     !
     ! *) Fill blocks that are diagonal in spin coordinates.
     !    row1,col1 refer to the resonat block.
     !    row2,col2 refer to the anti-resonant block.
     do block=1,2
       if (block==1) then 
         spad_r=0
         spad_c=SUM(nreh)
       end if
       if (block==2) then 
         spad_r=nreh(1)
         spad_c=2*nreh(1)+nreh(2)
       end if
       do itp=1,nreh(block)
         read(funt) cbuff_dpc(1:itp)
         col1 = itp+spad_c
         row2 = itp+spad_r
         do it=1,itp-1
           cttp = cbuff_dpc(it)
           row1 = it + spad_r
           col2 = it + spad_c
           exc_ham(row1,col1) =                cttp  ! upper coupling
           exc_ham(row2,col2) =                cttp  ! lower coupling (symmetric)
           exc_ham(col1,row1) = row_sign*CONJG(cttp) ! lower anti-coupling
           exc_ham(col2,row2) = row_sign*CONJG(cttp) ! upper anti-couling
         end do
         ! TODO recheck this 
         row1 = itp+spad_r
         exc_ham(row1,col1) = cbuff_dpc(itp)                  ! Diagonals of the block.
         col2 = itp+spad_c
         exc_ham(col2,row2) = row_sign*CONJG(cbuff_dpc(itp))
       end do
     end do
     !
     ! Read Full v+- and reconstruct resonant and anti-resonat blocks.
     ! TODO recheck this 
     spad=2*nreh(1) + nreh(2)
     pad_r1=SUM(nreh)
     pad_c1=2*nreh(1) + nreh(2)

     spin_dim=nreh(1)
     do itp=1,spin_dim
       read(funt) cbuff_dpc(1:spin_dim)
       exc_ham(1:spin_dim,spad+itp) = cbuff_dpc(1:spin_dim)                                  ! upper block reso
       exc_ham(1+pad_r1:pad_r1+spin_dim,nreh(1)+itp) = row_sign*CONJG(cbuff_dpc(1:spin_dim)) ! upper block anti-reso
       col1 = itp+spad
       row2 = itp+nreh(1)
       do it=1,spin_dim
         cttp = cbuff_dpc(it) 
         row1 = it
         exc_ham(col1,row1) =  row_sign*CONJG(cttp)  ! lower anti-reso.
         col2 = it + SUM(nreh)
         exc_ham(row2,col2) = cttp                   ! lower reso.
       end do
     end do
     !
   end if

 CASE DEFAULT
   msg = " Unknown block_type: "//TRIM(block_type)
   MSG_ERROR(msg)
 END SELECT

 ABI_DEALLOCATE(cbuff_dpc)

end subroutine exc_fullh_from_blocks
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/rrs_of_glob
!! NAME
!!  rrs_of_glob
!!
!! FUNCTION
!!   [+1,-1,0] if (row_glob,col_glob) belongs to the [ resonant, anti-resonant, (anti)coupling block ]
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function rrs_of_glob(row_glob,col_glob,size_glob)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rrs_of_glob'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: rrs_of_glob
 integer,intent(in) :: row_glob,col_glob
 integer,intent(in) :: size_glob(2)

!Local variables ------------------------------
!scalars
 integer :: nreh1,nreh2

! *************************************************************************

 nreh1=size_glob(1)/2 ! Matrix is square and nreh1==nreh2 but oh well.
 nreh2=size_glob(2)/2

 if (row_glob<=nreh1 .and. col_glob<=nreh2) then
   rrs_of_glob=+1  ! Resonant.
 else if (row_glob >nreh1 .and. col_glob >nreh2) then
   rrs_of_glob=-1  ! anti-Resonant.
 else 
   rrs_of_glob=0
 end if

end function rrs_of_glob
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/ccs_of_glob
!! NAME
!!! ccs_of_glob
!!
!! FUNCTION
!!  [+1,-1,0] if (row_glob,col_glob) belongs to the [ coupling, anti-coupling, (anti)resonant block ]
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ccs_of_glob(row_glob,col_glob,size_glob)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ccs_of_glob'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: ccs_of_glob
 integer,intent(in) :: row_glob,col_glob
 integer,intent(in) :: size_glob(2)

!Local variables ------------------------------
!scalars
 integer :: nreh1,nreh2

! *************************************************************************

 nreh1=size_glob(1)/2 ! Matrix is square and nreh1==nreh2 but oh well.
 nreh2=size_glob(2)/2

 if (row_glob<=nreh1 .and. col_glob >nreh2) then      ! Coupling.
   ccs_of_glob = +1
 else if (row_glob >nreh1 .and. col_glob<=nreh2) then ! anti-Coupling
   ccs_of_glob = -1
 else 
   ccs_of_glob = 0
 end if

end function ccs_of_glob
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/offset_in_file
!! NAME
!!  offset_in_file
!!
!! FUNCTION
!!  Return the offset of the matrix element (row_glob,col_glob)
!!  size_glob(2) gives the number of row and column of the global matrix
!!  nsblocks is the number of sublocks, used for nsppol==2 (not used if 1)
!!  sub_block(2,2,nsblocks)= For each subblock the coordinates of the first and last element.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function offset_in_file(row_glob,col_glob,size_glob,nsblocks,sub_block,bsize_elm,bsize_frm)
                                                                                                      
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'offset_in_file'
!End of the abilint section

 implicit none
                                                                                                      
 !Arguments ------------------------------------
 integer(XMPI_OFFSET_KIND) :: offset_in_file
 integer,intent(in) :: row_glob,col_glob,nsblocks,bsize_elm,bsize_frm
 integer,intent(in) :: size_glob(2),sub_block(2,2,nsblocks)

!Local variables ------------------------------
!scalars
 integer :: ii,jj,ijp_glob,swap
 integer(XMPI_OFFSET_KIND) :: my_offset
                                                                                                      
 ! *************************************************************************

 if (nsblocks==1) then
   ii = row_glob
   jj = col_glob
   if (ii>size_glob(1)/2) ii = ii - size_glob(1)/2 ! Wrap the index.
   if (jj>size_glob(2)/2) jj = jj - size_glob(2)/2
   if (jj<ii) then ! Exchange the indeces since the symmetric element is read.
     swap = jj
     jj   = ii
     ii   = swap
   end if
   ijp_glob = ii + jj*(jj-1)/2  ! Index for packed storage mode.
   my_offset = (ijp_glob-1)*bsize_elm + (jj-1)*2*bsize_frm 
 else 
   ABI_UNUSED(sub_block(1,1,1))
   MSG_ERROR("nsppol==2 not coded")
 end if

 offset_in_file = my_offset

end function offset_in_file
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/exc_read_rblock_fio
!! NAME
!!  exc_read_rblock_fio
!!
!! FUNCTION
!!  Reads the resonant block from file using Fortran IO.
!!
!! INPUTS
!!  funt=Fortran unit number.
!!  nsppol=Number of spins
!!  exc_size=Size of the resonant bock.
!!  diago_is_real=.TRUE. if diagonal elements are real.
!!  nreh(nsppol)=Number of resonant transitions for each spin.
!!
!! OUTPUT
!!  ierr=Status error
!!  exc_mat(exc_size,exc_size)=The resonant block.
!!
!! OUTPUT
!!
!! PARENTS
!!      exc_diago
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_read_rblock_fio(funt,diago_is_real,nsppol,nreh,exc_size,exc_mat,ierr)
                                                                                                      
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_read_rblock_fio'
!End of the abilint section

 implicit none
                                                                                                      
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: funt,nsppol,exc_size
 logical,intent(in) :: diago_is_real
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nreh(nsppol)
 complex(dpc),intent(out) :: exc_mat(exc_size,exc_size)

!Local variables ------------------------------
!scalars
 integer :: itp,it,spin_dim,block,col,row,spad
 complex(dpc) :: ctemp
!arrays
 complex(dpc),allocatable :: cbuff_dpc(:)

 ! *************************************************************************

 ierr=0
 !
 ! Construct full resonant block using Hermiticity. File is always in double precision.
 ABI_ALLOCATE(cbuff_dpc,(exc_size))
                                                                                        
 if (nsppol==1) then ! Construct resonant block from the upper triangle stored on file.
   !
   do itp=1,nreh(1)
     read(funt,ERR=10) cbuff_dpc(1:itp)
     do it=1,itp-1 ! Diagonal is treated below.
       ctemp = cbuff_dpc(it)
       exc_mat(it,itp) = ctemp
       exc_mat(itp,it) = CONJG(ctemp)
     end do
     if (diago_is_real) then
       exc_mat(itp,itp) = DBLE(cbuff_dpc(itp))  
     else
       exc_mat(itp,itp) = cbuff_dpc(itp)
     end if
   end do
   !
 else
   ! The file contains 
   ! A) The up-up and the down-down block in packed form 
   ! (only the upper triangle is stored since these blocks are Hermitian)
   ! B) The entire up-down exchange block (no symmetry here)
   !
   ! A) Construct resonant blocks from the upper triangles stored on file.
   ! FIXME this part wont work if we have a different number of e-h pairs
   ABI_CHECK(ALL(nreh==nreh(1)),"Different number of transitions")
                                                                                        
   spin_dim=nreh(1)
   do block=1,2
     if (block==1) spad=0
     if (block==2) spad=nreh(1)
     do itp=1,nreh(block)
       read(funt,ERR=10) cbuff_dpc(1:itp)
       col = itp+spad
       do it=1,itp-1 ! diagonal is treated below.
         ctemp = cbuff_dpc(it)
         row = it + spad
         exc_mat(row,col) = ctemp
         exc_mat(col,row) = CONJG(ctemp)
       end do
       if (diago_is_real) then
         exc_mat(col,col) = DBLE(cbuff_dpc(itp))
       else 
         exc_mat(col,col) = cbuff_dpc(itp)
       end if
     end do
   end do
   !
   spad=nreh(1)
   spin_dim=nreh(1)
   do itp=1,spin_dim
     read(funt,ERR=10) cbuff_dpc(1:spin_dim)
     exc_mat(1:spin_dim,spad+itp) = cbuff_dpc(1:spin_dim)
   end do
 end if

 RETURN
 !
 ! Raise the error.
10 continue 
 ierr=1

end subroutine exc_read_rblock_fio
!!***

!----------------------------------------------------------------------

!!****f* m_bse_io/exc_amplitude
!! NAME
!!  exc_amplitude
!!
!! FUNCTION
!!  Calculate the amplitude function of the excitonic eigenstate |exc_vec\> 
!!    F(w) = \sum_t |<t|exc_vec>|^2 \delta(ww- ene_t) where the sum over t is done 
!!  of the full set of transitions used to construct the BS Hamiltoniam.
!!
!! INPUTS
!!  Bsp<excparam>=Structure storing the parameters of the run.
!!  eig_fname=The name of the file storing the excitonic eigenvectors.
!!  nvec=Number of excitonic states to analyze.
!!  vec_idx(nvec)=List with the indeces of the excitonic states sorted in ascending order.
!!  out_fname=The name of the file where the results are written.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!
!! CHILDREN
!!      exc_read_eigen
!!
!! SOURCE

subroutine exc_amplitude(Bsp,eig_fname,nvec,vec_idx,out_fname)
                                                                                                      
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_amplitude'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nvec
 character(len=*),intent(in) :: eig_fname,out_fname
 type(excparam),intent(in) :: BSp
! arrays
 integer,intent(in) :: vec_idx(nvec)

!Local variables ------------------------------
!scalars
 integer :: ios,istat,vec,art_idx 
 integer :: spin,iw,it,nw,pos_w,neg_w,out_unt,rt_idx,hsize
 real(dp) :: ene_rt,ampl_eh,ampl_he,xx,stdev,w_max,w_min,step
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: wmesh(:),amplitude(:),ene_list(:)
 complex(dpc),allocatable :: vec_list(:,:)

 ! *************************************************************************

 ! Setup of the frequency mesh for F(w).
 w_min=greatest_real; w_max=smallest_real
 do spin=1,BSp%nsppol
   do it=1,BSp%nreh(spin)
     ene_rt  = Bsp%Trans(it,spin)%en
     w_min = MIN(w_min,ene_rt)
     w_max = MAX(w_max,ene_rt)
   end do
 end do

 step = Bsp%domega
 if (Bsp%use_coupling==0) then
   nw = (w_max - w_min)/step + 1
   ABI_ALLOCATE(wmesh,(nw))
   wmesh = arth(w_min,step,nw)
 else
   ! Both positive and negative frequencies are needed.
   pos_w = (w_max - w_min)/step + 1
   neg_w = pos_w; if (ABS(w_min) < tol6) neg_w=neg_w-1 ! zero should not included twice.
   nw = neg_w + pos_w
   ABI_ALLOCATE(wmesh,(nw))
   wmesh(1:neg_w)  = arth(-w_max,step,neg_w)
   wmesh(neg_w+1:) = arth( w_min,step,pos_w)
 end if
 !
 ! Read selected eigenvectors.
 hsize = SUM(Bsp%nreh); if (Bsp%use_coupling>0) hsize=2*hsize

 ABI_ALLOCATE(ene_list,(nvec))
 ABI_ALLOCATE(vec_list,(hsize,nvec))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in vec_list")

 call exc_read_eigen(eig_fname,hsize,nvec,vec_idx,vec_list,ene_list,Bsp=Bsp)

 ABI_DEALLOCATE(ene_list)

 out_unt = get_unit()
 open(unit=out_unt,file=out_fname,form='formatted',iostat=ios)
 msg="Opening file "//TRIM(out_fname)//" as unformatted"
 ABI_CHECK(ios==0,msg) 

 write(out_unt,*)"# Amplitude functions F(w) = \sum_t |<t|exc_vec>|^2 \delta(ww- ene_t), w is given in eV. "
 write(out_unt,*)"# Number of excitonic eigenvectors analyzed: ",nvec

 ABI_ALLOCATE(amplitude,(nw))
 stdev = BSp%broad ! Broadening for the gaussian.

 do vec=1,nvec
   !
   ! amplitude(ww) = \sum_t |<t|exc_vec>|^2 \delta(ww- ene_t)
   amplitude = zero
   do spin=1,BSp%nsppol
     do it=1,BSp%nreh(spin)
      ene_rt  = Bsp%Trans(it,spin)%en
      rt_idx  = it + (spin-1)*Bsp%nreh(1)
      ampl_eh = (ABS(vec_list(rt_idx,vec)))**2
      if (Bsp%use_coupling>0) then ! Need the index and the amplitude of the anti-resonant component.
        art_idx  = it + (spin-1)*Bsp%nreh(1) + SUM(Bsp%nreh)
        ampl_he = (ABS(vec_list(art_idx,vec)))**2
      end if
      !
      do iw=1,nw ! Accumulate
        xx = wmesh(iw) - ene_rt
        amplitude(iw) = amplitude(iw) + ampl_eh*gaussian(xx,stdev)
        if (Bsp%use_coupling>0) then
          xx = wmesh(iw) + ene_rt
          amplitude(iw) = amplitude(iw) + ampl_he*gaussian(xx,stdev)
        end if
      end do
      !
     end do
   end do
   !
   ! Write results
   write(out_unt,*)"# Amplitude function F(w) for exc_vector index ",vec_idx(vec)
   do iw=1,nw 
     write(out_unt,*)wmesh(iw)*Ha_eV,amplitude(iw)
   end do
   write(out_unt,*)"#"
 end do

 close(out_unt)

 ABI_DEALLOCATE(amplitude)
 ABI_DEALLOCATE(wmesh)
 ABI_DEALLOCATE(vec_list)

end subroutine exc_amplitude
!!***

!----------------------------------------------------------------------

END MODULE m_bse_io
!!***
