!{\src2tex{textfont=tt}}
!!****f* ABINIT/exc_diago_resonant
!! NAME
!!  exc_diago_resonant
!!
!! FUNCTION
!!  Calculates eigenvalues and eigenvectors of the Hermitian excitonic Hamiltonian (coupling is neglected).
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT and EXC groups (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  Bsp
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!  comm=MPI communicator.
!!  bseig_fname=The name of the output file
!!  prtvol=Verbosity level.
!!
!! OUTPUT
!!  Eigenvalues and eigenvectors are written on file bseig_fname
!!
!! PARENTS
!!      exc_diago_driver
!!
!! CHILDREN
!!      destruction_matrix_scalapack,end_scalapack,exc_fullh_from_blocks
!!      exc_read_bshdr,exc_skip_bshdr_mpio,flush_unit,hermitianize,idx_glob
!!      init_matrix_scalapack,init_scalapack,mpi_file_close,mpi_file_open
!!      mpi_file_read_all,mpi_file_set_view,mpi_type_free,slk_pzgemm
!!      slk_pzhegvx,slk_single_fview_read_mask,slk_write,slk_zinvert,wrtout
!!      xbarrier_mpi,xgemm,xhdp_invert,xhegv,xhegvx,xmpio_read_frm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine exc_diago_resonant(Bsp,BS_files,Hdr_bse,prtvol,comm)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use m_bse_io
 use defs_scalapack
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_numeric_tools,  only : hermitianize
 use m_io_tools,       only : get_unit, flush_unit
 use m_abilasi,        only : xheev, xheevx
 use defs_abitypes,    only : hdr_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_diago_resonant'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,prtvol
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse

!Local variables ------------------------------
!scalars
 integer :: ii,it,mi,hreso_unt,eig_unt,ios,exc_size,neh1,neh2
 integer :: istat,nsppol,il,iu,mene_found,nstates
 integer :: nprocs,my_rank,master,fform,nene_printed,ierr
 real(dp) :: exc_gap,exc_maxene,abstol
 real(dp) :: vl,vu
 logical :: use_scalapack,do_full_diago,diagonal_is_real
 character(len=500) :: msg
 character(len=fnlen) :: hreso_fname,bseig_fname
!arrays
 real(dp),allocatable :: exc_ene(:)
 complex(dpc),allocatable :: exc_mat(:,:),exc_vec(:,:)
#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
 integer :: amode,mpi_fh,istwf_k,tbloc,tmp_unt
 integer :: itloc,jj,jtloc,itglob,jtglob
 integer(XMPI_OFFSET_KIND) :: ehdr_offset,fmarker
 integer :: block_sizes(2,3),array_of_sizes(2),gsub(2,2)
 logical,parameter :: is_fortran_file=.TRUE.
 real(dp),external :: PDLAMCH
 type(matrix_scalapack)    :: Slk_mat,Slk_vec 
 type(processor_scalapack) :: Slk_processor
#endif

!************************************************************************

 DBG_ENTER("PERS")

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master=0

 if (BSp%have_complex_ene) then ! QP lifetimes are not included
   MSG_ERROR("complex energies not coded yet")
 end if

 if (ANY(Bsp%nreh/=Bsp%nreh(1))) then
   write(std_out,*)" Bsp%nreh: ",Bsp%nreh
   write(msg,'(a)')" BSE code does not support different number of transitions for the two spin channels"
   MSG_ERROR(msg)
 end if

 nsppol   = Hdr_bse%nsppol
 exc_size = SUM(BSp%nreh)
 nstates  = BSp%nstates; do_full_diago=(Bsp%nstates==exc_size)

 neh1 = Bsp%nreh(1); neh2 = neh1
 if (Hdr_bse%nsppol==2) neh2 = Bsp%nreh(2)

 use_scalapack = .FALSE.
#if defined HAVE_LINALG_SCALAPACK
 use_scalapack = (nprocs > 1)
#endif
 !use_scalapack = .FALSE.
 !use_scalapack = .TRUE.

 if (.not.use_scalapack .and. my_rank/=master) GOTO 10 ! Inversion is done by master only.

 nene_printed = MIN(32*nsppol,nstates); if (prtvol>10) nene_printed = nstates

 if (BS_files%in_hreso /= BSE_NOFILE) then
   hreso_fname = BS_files%in_hreso
 else 
   hreso_fname = BS_files%out_hreso
 end if

 bseig_fname = BS_files%out_eig
 if (BS_files%in_eig /= BSE_NOFILE) then
   MSG_ERROR("BS_files%in_eig is defined!")
 end if

 write(msg,'(a,i8)')' Direct diagonalization of the resonant excitonic Hamiltonian, Matrix size= ',exc_size
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")

 ABI_ALLOCATE(exc_ene,(exc_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out of memory: excitonic eigenvalues')

 SELECT CASE (use_scalapack)
 CASE (.FALSE.)

   write(msg,'(a)')". Using LAPACK sequential version. "
   call wrtout(std_out,msg,"PERS")
   call wrtout(ab_out,msg,"COLL")

   write(msg,'(a,f8.1,a)')' Allocating excitonic eigenvalues. Memory required: ', exc_size*dp*b2Mb,' Mb. '
   call wrtout(std_out,msg,"COLL")
                                                                                                        
   write(msg,'(a,f8.1,a)')' Allocating excitonic hamiltonian.  Memory required: ',exc_size**2*dpc*b2Mb,' Mb.'
   call wrtout(std_out,msg,"COLL")
   call flush_unit(std_out)

   ABI_ALLOCATE(exc_mat,(exc_size,exc_size))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,'out of memory: excitonic hamiltonian')
   !
   ! Read data from file.
   hreso_unt = get_unit()
   open(unit=hreso_unt,file=hreso_fname,status='old',form="unformatted",iostat=ios)
   msg=' Opening file '//TRIM(hreso_fname)//' as old-formatted.'
   ABI_CHECK(ios==0,msg)
   !
   ! Read the header and perform consistency checks.
   call exc_read_bshdr(hreso_unt,Bsp,fform,ierr)
   ABI_CHECK(ierr==0,"Fatal error, cannot continue")
   !
   ! Construct full resonant block using Hermiticity. 
   diagonal_is_real = .not.Bsp%have_complex_ene
   call exc_read_rblock_fio(hreso_unt,diagonal_is_real,nsppol,Bsp%nreh,exc_size,exc_mat,ierr)
   ABI_CHECK(ierr==0,"Fatal error, cannot continue")

   close(hreso_unt)

   if (do_full_diago) then 
     call wrtout(std_out," Full diagonalization with XHEEV... ","COLL")
     call xheev("Vectors","Upper",exc_size,exc_mat,exc_ene)
   else
     call wrtout(std_out," Partial diagonalization with XHEEVX... ","COLL")
     abstol=zero; il=1; iu=nstates
     ABI_ALLOCATE(exc_vec,(exc_size,nstates))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,"out of memory in exc_vec")
     call xheevx("Vectors","Index","Upper",exc_size,exc_mat,vl,vu,il,iu,abstol,mene_found,exc_ene,exc_vec,exc_size)
     exc_mat(:,1:nstates) = exc_vec
     ABI_DEALLOCATE(exc_vec)
   end if
   !
   ! ==============================================
   ! === Now exc_mat contains the eigenvectors ====
   ! ==============================================

   ! * Write the final results.
   call wrtout(std_out,' Writing eigenvalues and eigenvectors on file: '//TRIM(bseig_fname),"COLL")

   eig_unt=get_unit()
   open(eig_unt,file=bseig_fname,form='unformatted',iostat=ios)
   ABI_CHECK(ios==0,"Opening file:"//TRIM(bseig_fname))

   !% fform = 1002 ! FIXME
   !% call hdr_io_int(fform,Hdr_bse,2,eig_unt)

   write(eig_unt) exc_size, nstates
   write(eig_unt) CMPLX(exc_ene(1:nstates),kind=dpc)
   do mi=1,nstates
     write(eig_unt) exc_mat(1:exc_size,mi)
   end do

   close(eig_unt)

   ABI_DEALLOCATE(exc_mat)

 CASE (.TRUE.)

#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
   if (nsppol==2) then 
     MSG_WARNING("nsppol==2 + scalapack not coded yet")
   end if

   istwf_k=1; tbloc=50
   write(msg,'(2(a,i0))')". Using scaLAPACK version with nprocs= ",nprocs,"; block size= ",tbloc
   call wrtout(std_out,msg,"PERS")
   call wrtout(ab_out,msg,"COLL")

   write(msg,'(a,f8.1,a)')' Allocating excitonic eigenvalues. Memory required: ',exc_size*dp*b2Mb,' Mb. '
   call wrtout(std_out,msg,"PERS")
   !
   ! Init scaLAPACK environment.     
   call init_scalapack(Slk_processor,comm)
   !
   ! Init scaLAPACK matrices
   call init_matrix_scalapack(Slk_mat,exc_size,exc_size,Slk_processor,istwf_k,tbloc=tbloc)

   call init_matrix_scalapack(Slk_vec,exc_size,exc_size,Slk_processor,istwf_k,tbloc=tbloc) 
   !
   ! Open the file with MPI-IO and skip the record.
   amode=MPI_MODE_RDONLY

   call MPI_FILE_OPEN(comm, hreso_fname, amode, MPI_INFO_NULL, mpi_fh, ierr)
   msg = " MPI_IO error opening file: "//TRIM(hreso_fname)
   ABI_CHECK_MPI(ierr,msg)

   ! Skip the header and find the offset for reading the matrix.
   call exc_skip_bshdr_mpio(mpi_fh,xmpio_at_all,ehdr_offset)
   !
   ! Read scaLAPACK matrix from the file.
   if (nsppol==1) then
     call slk_read(Slk_mat,"Upper","Hermitian",is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset)
   else 
     array_of_sizes = (/exc_size,exc_size/)
     block_sizes(:,1) = (/neh1,neh1/)
     block_sizes(:,2) = (/neh2,neh2/)
     block_sizes(:,3) = (/neh1,neh2/)
     MSG_ERROR("Not tested")
     !call slk_read_from_blocks(Slk_mat,array_of_sizes,block_sizes,is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset)
   end if

   call MPI_FILE_CLOSE(mpi_fh, ierr)               
   ABI_CHECK_MPI(ierr,"FILE_CLOSE")

   if (do_full_diago) then 
     call wrtout(std_out," Performing full diagonalization with scaLAPACK...","COLL")
     
     call slk_pzheev("Vectors","Upper",Slk_mat,Slk_vec,exc_ene)

   else 
     call wrtout(std_out," Performing partial diagonalization with scaLAPACK...","COLL")
     il=1; iu=nstates; abstol=zero !ABSTOL = PDLAMCH(comm,'U')
     call slk_pzheevx("Vectors","Index","Upper",Slk_mat,vl,vu,il,iu,abstol,Slk_vec,mene_found,exc_ene)
   end if

   call destruction_matrix_scalapack(Slk_mat)

   call wrtout(std_out,' Writing eigenvalues/vectors on file: '//TRIM(bseig_fname),"COLL")
   call flush_unit(std_out)

   ! Write distributed matrix on file bseig_fname with Fortran records.
   if (my_rank==master) then ! Write exc eigenvalues. Vectors will be appended in slk_write.
     eig_unt=get_unit()
     open(eig_unt,file=bseig_fname,form='unformatted')
     !% fform = 1002 ! TODO 
     !% call hdr_io_int(fform,Hdr_bse,2,eig_unt)
     write(eig_unt) exc_size, nstates
     write(eig_unt) CMPLX(exc_ene(1:nstates),kind=dpc)
     close(eig_unt)
   end if

   call xbarrier_mpi(comm)
   !
   ! Open the file with MPI-IO and skip the record.
   amode=MPI_MODE_RDWR
                                                                                
   call MPI_FILE_OPEN(comm, bseig_fname, amode, MPI_INFO_NULL, mpi_fh, ierr)
   msg = " MPI_IO error opening file: "//TRIM(hreso_fname)
   ABI_CHECK_MPI(ierr,msg)

   !call MPI_FILE_SYNC(mpi_fh,ierr)

   ehdr_offset = 0
   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,ierr)
   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,ierr)

   write(std_out,*)"Writing nstates ",nstates
   gsub(:,1) = (/1,1/)
   gsub(:,2) = (/exc_size,nstates/)
   call slk_write(Slk_vec,"All",is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset,glob_subarray=gsub)

   call MPI_FILE_CLOSE(mpi_fh, ierr)               
   ABI_CHECK_MPI(ierr,"FILE_CLOSE")

   call destruction_matrix_scalapack(Slk_vec)
   call end_scalapack(Slk_processor)
   call xbarrier_mpi(comm)
#else
   MSG_BUG("You should not be here!")
#endif

 END SELECT

 write(msg,'(a,i4)')' Excitonic eigenvalues in eV up to n= ',nene_printed 
 call wrtout(std_out,msg,"PERS")
 call wrtout(ab_out,msg,"COLL")
 do it=0,(nene_printed-1)/8
   write(msg,'(8f10.5)') ( exc_ene(ii)*Ha_eV, ii=1+it*8,MIN(it*8+8,nene_printed) )
   call wrtout(std_out,msg,"PERS")
   call wrtout(ab_out,msg,"COLL")
 end do
                                                                                   
 exc_gap    = MINVAL(exc_ene(1:nstates))
 exc_maxene = MAXVAL(exc_ene(1:nstates))

 write(msg,'(a,2(a,f7.2,2a),a)')ch10,&                          
&  " First excitonic eigenvalue= ",exc_gap*Ha_eV,   " [eV]",ch10,&
&  " Last  excitonic eigenvalue= ",exc_maxene*Ha_eV," [eV]",ch10,ch10
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")

 ABI_DEALLOCATE(exc_ene)
 call flush_unit(std_out)

10 call xbarrier_mpi(comm)

 DBG_EXIT("PERS")

end subroutine exc_diago_resonant
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/exc_print_eig
!! NAME
!!  exc_print_eig
!!
!! FUNCTION
!!  Print excitonic eigenvalues on std_out and ab_out.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gw_gap=GW direct gap.
!!  bseig_fname=The name of file containing eigenvalues and eigenvectors
!!
!! OUTPUT
!!  exc_gap=Excitonic direct gap.
!!  Additional info on the Excitonic spectrum are reported on standard output.
!!
!! PARENTS
!!      exc_diago_driver
!!
!! CHILDREN
!!      destruction_matrix_scalapack,end_scalapack,exc_fullh_from_blocks
!!      exc_read_bshdr,exc_skip_bshdr_mpio,flush_unit,hermitianize,idx_glob
!!      init_matrix_scalapack,init_scalapack,mpi_file_close,mpi_file_open
!!      mpi_file_read_all,mpi_file_set_view,mpi_type_free,slk_pzgemm
!!      slk_pzhegvx,slk_single_fview_read_mask,slk_write,slk_zinvert,wrtout
!!      xbarrier_mpi,xgemm,xhdp_invert,xhegv,xhegvx,xmpio_read_frm
!!
!! SOURCE

subroutine exc_print_eig(BSp,bseig_fname,gw_gap,exc_gap)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use m_errors

 use defs_abitypes,  only : hdr_type
 use m_io_tools,     only : get_unit
 use m_header,       only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_print_eig'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 complex(dpc),intent(in) :: gw_gap 
 complex(dpc),intent(out) :: exc_gap 
 character(len=*),intent(in) :: bseig_fname
 type(excparam),intent(in) :: BSp

!Local variables ------------------------------
!scalars
 integer :: nstates_read,ii,j,k,eig_unt,ieig,ios,hsize_exp
 integer :: hsize_read !,nstates
 complex(dpc) :: bind_energy,ctemp
 character(len=500) :: msg
 !type(Hdr_type) :: tmp_Hdr
!arrays
 integer,allocatable :: iperm(:)
 real(dp),allocatable :: exc_rene(:)
 complex(dpc),allocatable :: exc_cene(:)

!************************************************************************

 exc_gap = czero

 eig_unt = get_unit()
 open(unit=eig_unt,file=bseig_fname,form='unformatted',status='old',iostat=ios)
 ABI_CHECK(ios==0,"Opening file: "//TRIM(bseig_fname))

 !% call hdr_io_int(fform,tmp_Hdr,1,eig_unt)
 !% call hdr_clean(tmp_Hdr)

 read(eig_unt) hsize_read, nstates_read

 if (BSp%use_coupling==0) hsize_exp =   SUM(Bsp%nreh)
 if (BSp%use_coupling>0)  hsize_exp = 2*SUM(Bsp%nreh)

 if (hsize_exp /= hsize_read) then
   write(msg,'(2(a,i0))')" Wrong dimension: read: ",hsize_read," expected= ",hsize_exp
   MSG_ERROR(msg)
 end if

 ABI_ALLOCATE(exc_cene,(nstates_read))
 read(eig_unt) exc_cene(:)

 ABI_ALLOCATE(exc_rene,(nstates_read))
 exc_rene = DBLE(exc_cene)

 ABI_ALLOCATE(iperm,(nstates_read))
 iperm = (/(ii, ii=1,nstates_read)/) 

 call sort_dp(nstates_read,exc_rene,iperm,tol6)

 ABI_DEALLOCATE(exc_rene)
 ABI_DEALLOCATE(iperm)

 ! put in ascending order
 do ii=nstates_read,2,-1
   do j=1,ii-1
     if (DBLE(exc_cene(j)) > DBLE(exc_cene(j+1))) then
       ctemp = exc_cene(j)
       exc_cene(j) = exc_cene(j+1)
       exc_cene(j+1) = ctemp
     end if
   end do
 end do

 exc_gap = DCMPLX(ABS(DBLE(exc_cene(1))),AIMAG(exc_cene(1)))

 do ii=1,nstates_read
   if (ABS(DBLE(exc_cene(ii))) < DBLE(exc_gap)) then
     exc_gap = CMPLX(ABS(DBLE(exc_cene(ii))),AIMAG(exc_cene(ii)))
   end if
 end do

 bind_energy = gw_gap - exc_gap

 write(msg,"(3(a,2f6.2,2a))")&
&  " GW  direct gap     ",gw_gap*Ha_eV,     " [eV] ",ch10,&
&  " EXC direct gap     ",exc_gap*Ha_eV,    " [eV] ",ch10,&
&  " EXC binding energy ",bind_energy*Ha_eV," [eV] ",ch10
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")
 
 msg=' Excitonic eigenvalues up to the GW energy gap [eV]'
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")

 do ii=1,nstates_read
   if (DBLE(exc_cene(ii)) > zero) EXIT
 end do

 do j=ii,nstates_read
   if (DBLE(exc_cene(j)) > DBLE(gw_gap)) EXIT
 end do
 j=j-1

 do ieig=ii,j
   write(msg,'(i3,a,2f6.2,a)')ieig," (",exc_cene(ieig)*Ha_eV,")"
   call wrtout(std_out,msg,"COLL")
   call wrtout(ab_out,msg,"COLL")
 end do
 
 ii=ii-1
 do j=ii,1,-1
   if (ABS(DBLE(exc_cene(j))) > DBLE(gw_gap)) EXIT
 end do
 j=j+1

 ! This coding is not portable, write to ab_out has been disabled.
 if (ii>0) then 
   do k=ii,j,-1
     write(msg,'(i3,a,2f6.2,a)')k," (",exc_cene(k)*Ha_eV,")"
     call wrtout(std_out,msg,"COLL")
   end do
 end if

 close(eig_unt)

end subroutine exc_print_eig
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/exc_diago_coupling
!! NAME
!!  exc_diago_coupling
!!
!! FUNCTION
!!  Calculate excitonic eigenvalues and eigenvectors by performing a direct diagonalization.
!!  of the non Hermitian excitonic Hamiltoninan (resonant + coupling). 
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  bseig_fname=The name of the output file.
!!  Bsp
!!    neh=Rank of the resonant block of the Hamiltoninan (equal to the rank of the coupling part)
!!  comm=MPI communicator.
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!
!! OUTPUT
!!  Excitonic eigenvectors and eigenvalues are written on file BS_files%out_eig.
!!
!! PARENTS
!!      exc_diago_driver
!!
!! CHILDREN
!!      destruction_matrix_scalapack,end_scalapack,exc_fullh_from_blocks
!!      exc_read_bshdr,exc_skip_bshdr_mpio,flush_unit,hermitianize,idx_glob
!!      init_matrix_scalapack,init_scalapack,mpi_file_close,mpi_file_open
!!      mpi_file_read_all,mpi_file_set_view,mpi_type_free,slk_pzgemm
!!      slk_pzhegvx,slk_single_fview_read_mask,slk_write,slk_zinvert,wrtout
!!      xbarrier_mpi,xgemm,xhdp_invert,xhegv,xhegvx,xmpio_read_frm
!!
!! SOURCE

subroutine exc_diago_coupling(Bsp,BS_files,Hdr_bse,prtvol,comm)

 use m_profiling
      
 use defs_basis
 use defs_abitypes
 use m_bs_defs
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools,       only : get_unit
 use m_numeric_tools,  only : hermitianize
 use m_blas,           only : xdotc, xgemm
 use m_abilasi,        only : xgeev, xginv, xhdp_invert
 use m_bse_io,         only : exc_fullh_from_blocks, offset_in_file, rrs_of_glob, ccs_of_glob, exc_read_bshdr, exc_skip_bshdr_mpio

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_diago_coupling'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,prtvol
 type(excfiles),intent(in) :: BS_files
 type(excparam),intent(in) :: BSp
 type(Hdr_type),intent(in) :: Hdr_bse

!Local variables ------------------------------
!scalars
 integer(i8b) :: bsize_ham
 integer,parameter :: ldvl=1
 integer :: ii,exc_size,hreso_unt,ios,hcoup_unt,eig_unt,nsppol,nstates
 integer :: bsz,block,bs1,bs2,jj
 integer :: fform,row_sign
 integer :: mi,it,istat,nprocs,master,my_rank !itp
 integer :: nene_printed,ierr
 real(dp) :: exc_gap,exc_maxene,temp
 logical :: diago_is_real,do_full_diago
 character(len=500) :: msg
 character(len=fnlen) :: hreso_fname,hcoup_fname,bseig_fname
!arrays
 complex(dpc),allocatable :: exc_ham(:,:),exc_rvect(:,:),exc_ene(:),ovlp(:,:)
 complex(dpc),allocatable :: cbuff(:,:)
 complex(dpc) :: vl_dpc(ldvl,1)

!************************************************************************

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master=0

 nsppol = Hdr_bse%nsppol
 if (nsppol==2) then
   MSG_WARNING("nsppol==2 with coupling is still under development")
 end if

 exc_size = 2*SUM(BSp%nreh)
 nstates  = BSp%nstates
 do_full_diago = (exc_size==nstates)
 ABI_CHECK(do_full_diago,"Partial diago not coded yet")

 bseig_fname = BS_files%out_eig
 if (BS_files%in_eig /= BSE_NOFILE) then
   MSG_ERROR("BS_files%in_eig is defined!")
 end if
 !
 ! Only master performs the diagonalization since ScaLAPACK does not provide the parallel version of ZGEEV. 
 if (my_rank/=master) GOTO 10

 write(msg,'(a,i0)')' Direct diagonalization of the full excitonic Hamiltonian, Matrix size= ',exc_size
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")

 bsize_ham = 2*dpc*exc_size**2
 write(msg,'(a,f9.2,a)')' Allocating full excitonic Hamiltonian. Memory requested: ',bsize_ham*b2Gb,' Gb. '
 call wrtout(std_out,msg,"COLL")

 ABI_ALLOCATE(exc_ham,(exc_size,exc_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out of memory: full excitonic hamiltonian')

 write(msg,'(3a,f8.1,3a,f8.1,a)')&
&  ' Allocating excitonic eigenvalues and eigenvectors. ',ch10,&
&  ' Memory-space requested: ',2*dpc*exc_size*b2Gb,' Gb. ',ch10,&
&  ' Memory-space requested: ',bsize_ham*b2Gb,' Gb. '
 call wrtout(std_out,msg,"COLL")

 ABI_ALLOCATE(exc_ene,(exc_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out of memory: exc_ene')

 if (BS_files%in_hreso /= BSE_NOFILE) then 
   hreso_fname = BS_files%in_hreso
 else 
   hreso_fname = BS_files%out_hreso
 end if

 call wrtout(std_out,' Reading resonant excitonic Hamiltonian from '//TRIM(hreso_fname),"COLL")
 
 hreso_unt = get_unit()
 open(unit=hreso_unt,file=hreso_fname,status='old',form="unformatted",iostat=ios)
 msg=' Opening file '//TRIM(hreso_fname)//' as old-formatted.'
 ABI_CHECK(ios==0,msg)
 !
 ! Read the header and perform consistency checks.
 call exc_read_bshdr(hreso_unt,Bsp,fform,ierr)
 ABI_CHECK(ierr==0,"Wrong header")
 !
 ! Construct resonant and anti-resonant part of the excitonic Hamiltonian using Hermiticity. File is always in double precision.
 ! Fill exc_ham with ( R  0 )
 !                   ( 0 -R*)
!BEGINDEBUG
 exc_ham = HUGE(one)
!ENDDEBUG

 row_sign=-1; diago_is_real=(.not.BSp%have_complex_ene)
 call exc_fullh_from_blocks(hreso_unt,"Resonant",nsppol,row_sign,diago_is_real,BSp%nreh,exc_size,exc_ham)
 close(hreso_unt)

 if (BS_files%in_hcoup /= BSE_NOFILE) then
   hcoup_fname =  BS_files%in_hcoup
 else 
   hcoup_fname =  BS_files%out_hcoup
 end if

 call wrtout(std_out,' Reading coupling excitonic Hamiltonian from '//TRIM(hcoup_fname),"COLL")

 hcoup_unt = get_unit()
 open(unit=hcoup_unt,file=hcoup_fname,status='old',form="unformatted",iostat=ios)
 msg=' Opening file '//TRIM(hcoup_fname)//' as old-formatted.'
 ABI_CHECK(ios==0,msg)
 !
 ! Read the header and perform consistency checks.
 call exc_read_bshdr(hcoup_unt,Bsp,fform,ierr)
 ABI_CHECK(ierr==0,"Wrong header")
 !
 ! Fill exc_ham with ( 0  C) to have ( R   C )
 !                   (-C* 0)         (-C* -R*)
 row_sign=-1; diago_is_real=(.not.BSp%have_complex_ene) ! not used here
 call exc_fullh_from_blocks(hcoup_unt,"Coupling",nsppol,row_sign,diago_is_real,BSp%nreh,exc_size,exc_ham)

!BEGINDEBUG
 if (ANY(exc_ham==HUGE(one))) then
   write(msg,'(a,2(1x,i0))')"There is a bug in exc_fullh_from_blocks",COUNT(exc_ham==HUGE(one)),exc_size**2
   MSG_WARNING(msg)
   bsz = Bsp%nreh(1)
   ABI_ALLOCATE(cbuff,(bsz,bsz))
   block=0
   do jj=1,2*nsppol
     do ii=1,2*nsppol
       block=block+1
       bs1 = (ii-1)*bsz+1
       bs2 = (jj-1)*bsz+1
       cbuff = exc_ham(bs1:bs1+bsz-1,bs2:bs2+bsz-1)
       if (ANY(cbuff==HUGE(one))) then
         write(std_out,*)" for block ",ii,jj," found ",COUNT(cbuff==HUGE(one))," wrong entries"
       end if
     end do
   end do

   ABI_DEALLOCATE(cbuff)
   MSG_ERROR("Cannot continue")
 end if
!ENDDEBUG

 close(hcoup_unt)
 !
 ! ======================================================
 ! ==== Calculate right eigenvectors and eigenvalues ====
 ! ======================================================
 ABI_ALLOCATE(exc_rvect,(exc_size,exc_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory: excitonic eigenvectors")

 if (do_full_diago) then
   call wrtout(std_out,"Complete direct diagonalization with xgeev...","COLL")
   call xgeev("No_left_eigen","Vectors",exc_size,exc_ham,exc_size,exc_ene,vl_dpc,ldvl,exc_rvect,exc_size)
 else 
   MSG_ERROR("Not implemented error")
 end if

 ABI_DEALLOCATE(exc_ham)

 exc_gap    = MINVAL(ABS(DBLE (exc_ene(1:nstates))))
 exc_maxene = MAXVAL(ABS(DBLE (exc_ene(1:nstates))))
 temp       = MAXVAL(ABS(AIMAG(exc_ene(1:nstates))))

 write(msg,'(2(a,f7.2,2a),a,es9.2,2a)')&
&  " First excitonic eigenvalue: ",exc_gap*Ha_eV,   " [eV].",ch10,&
&  " Last  excitonic eigenvalue: ",exc_maxene*Ha_eV," [eV].",ch10,& 
&  " Largest imaginary part:     ",temp*Ha_eV,      " [eV] ",ch10
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")

 nene_printed = MIN(32*nsppol,nstates); if (prtvol>10) nene_printed = nstates

 ! This is not portable as the the eigenvalues calculated by ZGEEV are not sorted.
 ! Even two subsequent calculations with the same input on the same machine
 ! might produce different orderings. Might sort the eigenvalues though, just for printing.

 write(msg,'(a,i4)')' Complex excitonic eigenvalues in eV up to n= ',nene_printed 
 call wrtout(std_out,msg,"PERS")

 do it=0,(nene_printed-1)/4
   write(msg,'(8f10.5)') ( exc_ene(ii)*Ha_eV, ii=1+it*4,MIN(it*4+4,nene_printed) )
   call wrtout(std_out,msg,"PERS")
 end do

 call wrtout(std_out,ch10//" Writing eigenvalues and eigenvectors on file "//TRIM(bseig_fname),"COLL")

 eig_unt = get_unit()
 open(eig_unt,file=bseig_fname,form='unformatted',iostat=ios)
 msg=" Error Opening file: "//TRIM(bseig_fname)
 ABI_CHECK(ios==0,msg)

 !% fform = 1002 ! FIXME
 !% call hdr_io_int(fform,Hdr_bse,2,eig_unt)

 write(eig_unt)exc_size,nstates
 write(eig_unt)CMPLX(exc_ene(1:nstates),kind=dpc)
 do mi=1,nstates
   write(eig_unt) exc_rvect(:,mi)
 end do

 ABI_DEALLOCATE(exc_ene)

 ABI_ALLOCATE(ovlp,(nstates,nstates))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out of memory in ovlp matrix')

 call wrtout(std_out,' Calculating overlap matrix... ',"COLL")

 !do itp=1,nstates
 !  do it=1,nstates
 !    ovlp(it,itp) = xdotc(exc_size,exc_rvect(:,it),1,exc_rvect(:,itp),1)
 !  end do
 !end do
 call xgemm("C","N",exc_size,nstates,nstates,cone,exc_rvect,exc_size,exc_rvect,exc_size,czero,ovlp,nstates)
 ABI_DEALLOCATE(exc_rvect)
      
 call wrtout(std_out," Inverting overlap matrix... ","COLL")

 ! Version for generic complex matrix.
 !call xginv(ovlp,exc_size)

 ! The overlap is Hermitian definite positive.
 call xhdp_invert("Upper",ovlp,nstates)
 call hermitianize(ovlp,"Upper")

 call wrtout(std_out,' Writing overlap matrix S^-1 on file: '//TRIM(bseig_fname),"COLL")

 do it=1,nstates
   write(eig_unt) CMPLX(ovlp(:,it),kind=dpc)
 end do

 ABI_DEALLOCATE(ovlp)

 close(eig_unt)

10 call xbarrier_mpi(comm)
          
end subroutine exc_diago_coupling
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/exc_diago_coupling_hegv
!! NAME
!!  exc_diago_coupling_hegv
!!
!! FUNCTION
!!  Calculate excitonic eigenvalues and eigenvectors by performing a direct diagonalization.
!!  of the non Hermitian excitonic Hamiltonian (resonant + coupling). 
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  bseig_fname=The name of the output file.
!!  Bsp
!!    neh=Rank of the resonant block of the Hamiltoninan (equal to the rank of the coupling part)
!!  comm=MPI communicator.
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!
!! OUTPUT
!!  Excitonic eigenvectors and eigenvalues are written on file BS_files%out_eig.
!!
!! PARENTS
!!      exc_diago_driver
!!
!! CHILDREN
!!      destruction_matrix_scalapack,end_scalapack,exc_fullh_from_blocks
!!      exc_read_bshdr,exc_skip_bshdr_mpio,flush_unit,hermitianize,idx_glob
!!      init_matrix_scalapack,init_scalapack,mpi_file_close,mpi_file_open
!!      mpi_file_read_all,mpi_file_set_view,mpi_type_free,slk_pzgemm
!!      slk_pzhegvx,slk_single_fview_read_mask,slk_write,slk_zinvert,wrtout
!!      xbarrier_mpi,xgemm,xhdp_invert,xhegv,xhegvx,xmpio_read_frm
!!
!! SOURCE

subroutine exc_diago_coupling_hegv(Bsp,BS_files,Hdr_bse,prtvol,comm)

 use m_profiling
      
 use defs_basis
 use defs_abitypes
 use m_bs_defs
 use defs_scalapack
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools,       only : get_unit, flush_unit
 use m_blas,           only : xdotc, xgemm
 use m_abilasi,        only : xhegvx, xginv, xhdp_invert, xhegv
 use m_numeric_tools,  only : print_arr, hermitianize
 use m_bse_io,         only : exc_fullh_from_blocks, offset_in_file, rrs_of_glob, ccs_of_glob, exc_read_bshdr, exc_skip_bshdr_mpio

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_diago_coupling_hegv'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,prtvol
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse

!Local variables ------------------------------
!scalars
 integer(i8b) :: bsize_ham
 integer :: itype,il,iu,spin,row1,row2,pad_r1,pad_r2,neh1,neh2
 integer :: ii,exc_size,hreso_unt,ios,hcoup_unt,eig_unt 
 integer :: fform,neig_found,nstates
 integer :: mi,it,istat,nprocs,master,my_rank
 integer :: nene_printed,nsppol,row_sign,ierr
 real(dp) :: exc_gap,exc_maxene,abstol,vl,vu
 character(len=500) :: msg
 character(len=fnlen) :: reso_fname,coup_fname,bseig_fname
 logical :: use_scalapack,do_full_diago,diago_is_real
!arrays
 real(dp),allocatable :: exc_ene(:) !,test_ene(:)
 complex(dpc),allocatable :: exc_ham(:,:),exc_rvect(:,:),fmat(:,:),ovlp(:,:)
#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
 integer,parameter :: istwfk1=1
 integer :: amode,mpi_fh,tbloc,tmp_unt,mene_found,mpi_err,my_nel,nsblocks
 integer :: iloc,jj,jloc,iglob,jglob,etype,slk_mask_type,offset_err,el,rrs_kind,ccs_kind
 integer :: max_r,max_c
 integer(XMPI_OFFSET_KIND) :: ehdr_offset,fmarker,my_offset
 integer :: gsub(2,2)
 logical,parameter :: is_fortran_file=.TRUE.
 complex(dpc) :: ctmp
 integer,allocatable :: sub_block(:,:,:)
 integer,pointer :: myel2loc(:,:)
 complex(dpc),allocatable :: tmp_cbuffer(:)
 character(50) :: uplo
 real(dp),external :: PDLAMCH
 type(matrix_scalapack)    :: Slk_F,Slk_Hbar,Slk_vec,Slk_ovlp,Slk_tmp
 type(processor_scalapack) :: Slk_processor
#endif

!************************************************************************

!#define DEV_MG_DEBUG_THIS

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master=0

 nsppol = Hdr_bse%nsppol
 if (nsppol==2) then 
   MSG_WARNING("nsppol==2 is still under development!")
 end if

 neh1 = BSp%nreh(1); neh2=neh1 
 if (nsppol==2) neh2 = BSp%nreh(2)
                                                                                                         
 exc_size = 2*SUM(Bsp%nreh)
 nstates  = Bsp%nstates
 do_full_diago=(nstates==exc_size)

 write(msg,'(a,i0)')'. Direct diagonalization of the full excitonic Hamiltonian, Matrix size= ',exc_size
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")

 bseig_fname = BS_files%out_eig
 if (BS_files%in_eig /= BSE_NOFILE) then
   MSG_ERROR("BS_files%in_eig is defined!")
 end if

 if (BS_files%in_hreso /= BSE_NOFILE) then 
   reso_fname = BS_files%in_hreso
 else 
   reso_fname = BS_files%out_hreso
 end if
 call wrtout(std_out,' Reading resonant excitonic Hamiltonian from '//TRIM(reso_fname),"COLL")

 if (BS_files%in_hcoup /= BSE_NOFILE) then
   coup_fname =  BS_files%in_hcoup
 else 
   coup_fname =  BS_files%out_hcoup
 end if
 call wrtout(std_out,' Reading coupling excitonic Hamiltonian from '//TRIM(coup_fname),"COLL")

 use_scalapack = .FALSE.
#ifdef HAVE_LINALG_SCALAPACK
 use_scalapack = (nprocs > 1)
#endif
 !use_scalapack = .FALSE.
 !use_scalapack = .TRUE.

 ABI_ALLOCATE(exc_ene,(exc_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out of memory: exc_ene')

 if (.not.use_scalapack .and. my_rank/=master) GOTO 10

 SELECT CASE (use_scalapack)

 CASE (.FALSE.)
   write(msg,'(a)')". Using LAPACK sequential version to solve FHv = ev with H positive definite. "
   call wrtout(std_out,msg,"PERS")
   call wrtout(ab_out,msg,"COLL")

   bsize_ham = 2*dpc*exc_size**2
   write(msg,'(a,f9.2,a)')' Allocating full excitonic Hamiltonian. Memory requested: ',2*bsize_ham*b2Gb,' Gb. '
   call wrtout(std_out,msg,"COLL")

   ABI_ALLOCATE(exc_ham,(exc_size,exc_size))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,'out of memory: full excitonic hamiltonian')

   ABI_ALLOCATE(fmat,(exc_size,exc_size))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,'out of memory: fmat')

   write(msg,'(3a,f8.1,3a,f8.1,a)')&
&    ' Allocating excitonic eigenvalues and eigenvectors. ',ch10,&
&    ' Memory-space requested: ',2*dpc*exc_size*b2Gb,' Gb. ',ch10,&
&    ' Memory-space requested: ',bsize_ham*b2Gb,' Gb. '
   call wrtout(std_out,msg,"COLL")

   hreso_unt = get_unit()
   open(unit=hreso_unt,file=reso_fname,status='old',form="unformatted",iostat=ios)
   msg=' Opening file '//TRIM(reso_fname)//' as old-formatted.'
   ABI_CHECK(ios==0,msg)
   !
   ! Read the header and perform consistency checks.
   call exc_read_bshdr(hreso_unt,Bsp,fform,ierr)
   ABI_CHECK(ierr==0,"Wrong header")
   !
   ! Construct Hbar = ( R   C )
   !                  ( C*  R*)
   !
   row_sign=+1; diago_is_real=(.not.BSp%have_complex_ene)
   call exc_fullh_from_blocks(hreso_unt,"Resonant",nsppol,row_sign,diago_is_real,Bsp%nreh,exc_size,exc_ham)
   close(hreso_unt)

   hcoup_unt = get_unit()
   open(unit=hcoup_unt,file=coup_fname,status='old',form="unformatted",iostat=ios)
   msg=' Opening file '//TRIM(coup_fname)//' as old-formatted.'
   ABI_CHECK(ios==0,msg)
   !
   ! Read the header and perform consistency checks.
   call exc_read_bshdr(hcoup_unt,Bsp,fform,ierr)
   ABI_CHECK(ierr==0,"Wrong header")

   row_sign=+1; diago_is_real=(.not.BSp%have_complex_ene) ! not used here.
   call exc_fullh_from_blocks(hcoup_unt,"Coupling",nsppol,row_sign,diago_is_real,Bsp%nreh,exc_size,exc_ham)
   close(hcoup_unt)

#ifdef DEV_MG_DEBUG_THIS
write(666)exc_ham
#endif
   !
   ! Fill fmat = (1  0)
   !             (0 -1)
   fmat = czero
   do spin=1,nsppol
     pad_r1 = (spin-1)*Bsp%nreh(1)
     pad_r2 = SUM(Bsp%nreh)
     if (spin==2) pad_r2 = pad_r2 + Bsp%nreh(1)
     do it=1,Bsp%nreh(spin)
       row1 = it + pad_r1
       row2 = it + pad_r2
       fmat(row1,row1) =  cone
       fmat(row2,row2) = -cone
     end do
   end do
   !
   ! ==================================================
   ! ==== Solve generalized EV problem F H u = e u ====
   ! ==================================================
   ! The eigenvectors Z are normalized as follows: if ITYPE = 1 or 2, Z**T*B*Z = I; if ITYPE = 3, Z**T*inv(B)*Z = I.
   !
   itype=2
   if (do_full_diago) then
     call wrtout(std_out," Full diagonalization with XHEGV... ","COLL")
     call xhegv(itype,"Vectors","Upper",exc_size,fmat,exc_ham,exc_ene)
   else 
     call wrtout(std_out," Partial diagonalization with XHEGVX... ","COLL")
     ABI_ALLOCATE(exc_rvect,(exc_size,nstates))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,"out of memory: excitonic eigenvectors")
     il=1; iu=1; abstol=zero
     call xhegvx(itype,"Vectors","All","Upper",exc_size,fmat,exc_ham,vl,vu,il,iu,abstol,neig_found,exc_ene,exc_rvect,exc_size)
   end if

   ABI_DEALLOCATE(exc_ham)

   if (do_full_diago) then
     ABI_ALLOCATE(exc_rvect,(exc_size,nstates))
     istat = ABI_ALLOC_STAT
     exc_rvect = fmat(:,1:nstates)
   end if

   ABI_DEALLOCATE(fmat)

   call wrtout(std_out," Writing eigenvalues and eigenvectors on file: "//TRIM(bseig_fname),"COLL")

   eig_unt = get_unit()
   open(eig_unt,file=bseig_fname,form='unformatted',iostat=ios)
   msg=" Error opening: "//TRIM(bseig_fname)
   ABI_CHECK(ios==0,msg)

   !% fform = 1002 ! FIXME
   !% call hdr_io_int(fform,Hdr_bse,2,eig_unt)

   write(eig_unt) exc_size, nstates
   write(eig_unt) CMPLX(exc_ene(1:nstates),kind=dpc)
   do mi=1,nstates
     write(eig_unt) CMPLX(exc_rvect(:,mi),kind=dpc)
   end do

#ifdef DEV_MG_DEBUG_THIS
   write(888)exc_rvect
   write(888)exc_ene
#endif

   ABI_ALLOCATE(ovlp,(nstates,nstates))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,'out of memory in ovlp matrix')

   call wrtout(std_out,' Calculating overlap matrix...',"COLL")

   call xgemm("C","N",exc_size,nstates,nstates,cone,exc_rvect,exc_size,exc_rvect,exc_size,czero,ovlp,nstates)
   ABI_DEALLOCATE(exc_rvect)

#ifdef DEV_MG_DEBUG_THIS
write(667)ovlp
#endif
        
   call wrtout(std_out," Inverting overlap matrix... ","COLL")
   !
   ! The overlap is Hermitian definite positive.
   call xhdp_invert("Upper",ovlp,nstates)
   call hermitianize(ovlp,"Upper")

   ! Version for generic complex matrix.
   !call xginv(ovlp,nstates)

#ifdef DEV_MG_DEBUG_THIS
write(668,*)ovlp
#endif

   call wrtout(std_out,' Writing overlap matrix O^-1 on file: '//TRIM(bseig_fname),"COLL")
   do it=1,nstates
     write(eig_unt) ovlp(:,it)
   end do

   ABI_DEALLOCATE(ovlp)
   close(eig_unt)

 CASE (.TRUE.)

#if defined HAVE_LINALG_SCALAPACK && defined HAVE_MPI_IO
   !
   ! Init scaLAPACK matrix Hbar = ( R   C )
   !                              ( C*  R*)
   ! Battle plan:
   !   Here the reading is complicated by the fact that R and C are stored on two different files 
   !   and moreover the matrices are in packed storage mode. 
   !   For initializing the local part of the resonant and anti-resonant block we have to allocate 
   !   a temporary buffer. Then we read the buffer from file and the corresponding elements of the 
   !   scaLAPACK matrix are initialized taking into account the symmetries of R (Hermitian)
   !   The same procedure is used to read the coupling and the anti-coupling part (Symmetric).
   !
   tbloc=50
   write(msg,'(2(a,i0))')". Using MPI-IO + scaLAPACK version with nprocs= ",nprocs,"; block size= ",tbloc
   call wrtout(std_out,msg,"PERS")
   call wrtout(ab_out,msg,"COLL")
   call flush_unit(std_out)
   !
   ! Init scaLAPACK environment.     
   call init_scalapack(Slk_processor,comm)
   !
   ! Open the Resonant file with MPI-IO and skip the record.
   amode=MPI_MODE_RDONLY
                                                                                
   call MPI_FILE_OPEN(comm, reso_fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
   msg = " MPI_IO error opening file: "//TRIM(reso_fname)
   ABI_CHECK_MPI(mpi_err,msg)
   !
   ! Skip the header and find the offset for reading the matrix.
   call exc_skip_bshdr_mpio(mpi_fh,xmpio_at_all,ehdr_offset)

   !call hdr_mpio_skip(mpi_fh,fform,ehdr_offset) 
   !call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,mpi_err)
   !write(std_out,*)"fmarker = ",fmarker
   !
   ! Read  = ( R  - )
   !         ( -  R*)
   call init_matrix_scalapack(Slk_Hbar,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc) 

   nullify(myel2loc)
   nsblocks=nsppol
   ABI_ALLOCATE(sub_block,(2,2,nsblocks))
   ABI_CHECK(nsppol==1,"nsppol==2 not coded yet")

   call slk_single_fview_read_mask(Slk_Hbar,rrs_of_glob,offset_in_file,nsblocks,sub_block,my_nel,myel2loc,etype,slk_mask_type,&
&    offset_err,is_fortran_file)

   if (offset_err/=0) then 
     write(msg,"(3a)")&
&      " Global position index cannot be stored in a standard Fortran integer ",ch10,&
&      " Excitonic matrix cannot be read with a single MPI-IO call."
     MSG_ERROR(msg)
   end if

   ! Shift the offset because the view starts at the fist matrix element!
   ! TODO should rationalize the treatment of the offset 
   my_offset = ehdr_offset + xmpio_bsize_frm
   call MPI_FILE_SET_VIEW(mpi_fh, my_offset, etype, slk_mask_type, 'native', MPI_INFO_NULL, mpi_err)
   ABI_CHECK_MPI(mpi_err,"SET_VIEW")

   call MPI_TYPE_FREE(slk_mask_type,mpi_err)
   ABI_CHECK_MPI(mpi_err,"MPI_TYPE_FREE")
   !
   ! Read my portion of the R,-R* sublocks and store the values in a temporary buffer.
   ABI_ALLOCATE(tmp_cbuffer,(my_nel))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0," out of memory tmp_cbuffer")

   call xbarrier_mpi(comm)
  
   call MPI_FILE_READ_ALL(mpi_fh, tmp_cbuffer, my_nel, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpi_err)
   ABI_CHECK_MPI(mpi_err,"READ_ALL")
   !
   ! Symmetrize my Resonant part.
   do el=1,my_nel
     iloc = myel2loc(1,el)
     jloc = myel2loc(2,el)
     call idx_glob(Slk_Hbar,iloc,jloc,iglob,jglob)
     ctmp = tmp_cbuffer(el)
     if (iglob==jglob.and..not.Bsp%have_complex_ene) ctmp = DBLE(ctmp) ! Force the diagonal to be real.
     rrs_kind = rrs_of_glob(iglob,jglob,Slk_Hbar%sizeb_global)
     if (rrs_kind==1.and.jglob<iglob) then ! Lower resonant
       ctmp = DCONJG(ctmp)
     else if (rrs_kind==-1.and.jglob>=iglob) then  ! Lower Anti-resonant (Diagonal is included).
       ctmp = DCONJG(ctmp)
     end if
     Slk_Hbar%buffer_cplx(iloc,jloc) = ctmp
   end do

   ABI_DEALLOCATE(tmp_cbuffer)
   ABI_DEALLOCATE(myel2loc)

   call MPI_FILE_CLOSE(mpi_fh, mpi_err)               
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")
   !
   ! Read  = ( -  C)
   !         (-C* -)
   !
   call MPI_FILE_OPEN(comm, coup_fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
   msg = " MPI_IO error opening file: "//TRIM(coup_fname)
   ABI_CHECK_MPI(mpi_err,msg)
   !
   ! Skip the header and find the offset for reading the matrix.
   call exc_skip_bshdr_mpio(mpi_fh,xmpio_at_all,ehdr_offset)

   !call hdr_mpio_skip(mpi_fh,fform,ehdr_offset) 
   !call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,mpi_err)
   !write(std_out,*)"coup fmarker = ",fmarker

   nullify(myel2loc)
   call slk_single_fview_read_mask(Slk_Hbar,ccs_of_glob,offset_in_file,nsblocks,sub_block,my_nel,myel2loc,etype,slk_mask_type,&
&    offset_err,is_fortran_file)

   ABI_DEALLOCATE(sub_block)

   if (offset_err/=0) then 
     write(msg,"(3a)")&
&      " Global position index cannot be stored in a standard Fortran integer ",ch10,&
&      " Excitonic matrix cannot be read with a single MPI-IO call."
     MSG_ERROR(msg)
   end if
   !
   ! Shift the offset because the view starts at the fist matrix element!
   ! TODO should rationalize the treatment of the offset so that the client code 
   ! will automatically receive my_offset.
   my_offset = ehdr_offset + xmpio_bsize_frm
   call MPI_FILE_SET_VIEW(mpi_fh, my_offset, etype, slk_mask_type, 'native', MPI_INFO_NULL, mpi_err)
   ABI_CHECK_MPI(mpi_err,"SET_VIEW")

   call MPI_TYPE_FREE(slk_mask_type,mpi_err)
   ABI_CHECK_MPI(mpi_err,"MPI_TYPE_FREE")
   !
   ! Read my portion of the C-C* blocks and store the values in a temporary buffer.
   ABI_ALLOCATE(tmp_cbuffer,(my_nel))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0," out of memory tmp_cbuffer")

   call MPI_FILE_READ_ALL(mpi_fh, tmp_cbuffer, my_nel, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpi_err)
   ABI_CHECK_MPI(mpi_err,"READ_ALL")
   !
   ! Symmetrize my coupling part.
   ! Coupling block is symmetric => No symmetrization of the lower triangle.
   do el=1,my_nel
     iloc = myel2loc(1,el)
     jloc = myel2loc(2,el)
     call idx_glob(Slk_Hbar,iloc,jloc,iglob,jglob)
     ccs_kind = ccs_of_glob(iglob,jglob,Slk_Hbar%sizeb_global)
     ctmp = tmp_cbuffer(el)
     if (ccs_kind==-1) ctmp = DCONJG(ctmp) ! Anti-coupling (Diagonal is included).
     Slk_Hbar%buffer_cplx(iloc,jloc) = ctmp
   end do
                                                                                                   
   ABI_DEALLOCATE(tmp_cbuffer)
   ABI_DEALLOCATE(myel2loc)

   !max_r=20; max_c=10
   !call print_arr(Slk_Hbar%buffer_cplx,max_r=max_r,max_c=max_c,unit=std_out)

#ifdef DEV_MG_DEBUG_THIS
   ABI_ALLOCATE(exc_ham,(exc_size,exc_size))
   istat = ABI_ALLOC_STAT
   read(666)exc_ham

   write(std_out,*)"Error Hbar: ",MAXVAL(ABS(exc_ham-Slk_Hbar%buffer_cplx))
   ABI_DEALLOCATE(exc_ham)
#endif
                                                                                                   
   call MPI_FILE_CLOSE(mpi_fh, mpi_err)               
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")
   !
   ! Init scaLAPACK matrix F
   call init_matrix_scalapack(Slk_F,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc)
   !
   ! Global F = (1  0)
   !            (0 -1)
   do jloc=1,Slk_F%sizeb_local(2) 
     do iloc=1,Slk_F%sizeb_local(1) 
       call idx_glob(Slk_F,iloc,jloc,iglob,jglob)
       if (iglob==jglob) then
         if (iglob<=SUM(Bsp%nreh)) then 
           Slk_F%buffer_cplx(iloc,jloc) =  cone
         else 
           Slk_F%buffer_cplx(iloc,jloc) = -cone
         end if
       else
         Slk_F%buffer_cplx(iloc,jloc) =  czero
       end if
     end do
   end do
   !
   ! ===========================================================
   ! ==== Solve generalized EV problem H u = F Hbar u = e u ====
   ! ===========================================================
   call init_matrix_scalapack(Slk_vec,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc) 
   !
   itype=2; vl=1; vu=1; il=1; iu=nstates
   abstol=zero !ABSTOL = PDLAMCH(comm,'U')

#if 1
   if (do_full_diago) then
     call slk_pzhegvx(itype,"Vectors","All","Upper",Slk_F,Slk_Hbar,vl,vu,il,iu,abstol,Slk_vec,mene_found,exc_ene)
   else 
     MSG_WARNING("Partial diago is still under testing")
     call slk_pzhegvx(itype,"Vectors","Index","Upper",Slk_F,Slk_Hbar,vl,vu,il,iu,abstol,Slk_vec,mene_found,exc_ene)
   end if
#else
   call xhegv(itype,"Vectors","Upper",exc_size,Slk_F%buffer_cplx,Slk_Hbar%buffer_cplx,exc_ene)
   Slk_vec%buffer_cplx = Slk_F%buffer_cplx
#endif

#ifdef DEV_MG_DEBUG_THIS
   if (PRODUCT(Slk_Hbar%sizeb_local) /= exc_size**2) then
     MSG_ERROR("Wrong size")
   end if

   ABI_ALLOCATE(exc_ham,(exc_size,exc_size))
   istat = ABI_ALLOC_STAT
   read(888)exc_ham

   write(std_out,*)"Error rvec: ",MAXVAL(ABS(exc_ham-Slk_vec%buffer_cplx))
   ABI_DEALLOCATE(exc_ham)

   ABI_ALLOCATE(test_ene,(exc_size))
   istat = ABI_ALLOC_STAT
   read(888)test_ene
   write(std_out,*)"Error ene: ",MAXVAL(ABS(exc_ene-test_ene))
   ABI_DEALLOCATE(test_ene)
#endif

   call destruction_matrix_scalapack(Slk_F)
   call destruction_matrix_scalapack(Slk_Hbar)

   call wrtout(std_out,ch10//" Writing eigenvalues and eigenvectors on file: "//TRIM(bseig_fname),"COLL")
   !
   ! Open the file with Fortran-IO to write the Header.
   if (my_rank==master) then
     eig_unt = get_unit()
     open(eig_unt,file=bseig_fname,form='unformatted',iostat=ios)
     ABI_CHECK(ios==0," Opening file: "//TRIM(bseig_fname))
     !%fform = 1002 ! TODO
     !%call hdr_io_int(fform,Hdr_bse,2,eig_unt)
                                                               
     write(eig_unt) exc_size,nstates
     write(eig_unt) CMPLX(exc_ene(1:nstates),kind=dpc)
     !do mi=1,exc_size
     !  write(eig_unt) CMPLX(exc_rvect(:,mi),kind=dpc)
     !end do
     close(eig_unt)
   end if
   !
   ! Open the file with MPI-IO and write the distributed eigevectors.
   call xbarrier_mpi(comm)
   amode=MPI_MODE_RDWR
   call MPI_FILE_OPEN(comm, bseig_fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err,"FILE_OPEN: "//TRIM(bseig_fname))
   !
   ! Skip the header and find the offset for writing the matrix.
   ehdr_offset=0
   !call hdr_mpio_skip(mpi_fh,fform,ehdr_offset) 

   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,mpi_err)
   write(std_out,*)" fmarker1 = ",fmarker
   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,mpi_err)
   write(std_out,*)" fmarker2 = ",fmarker

   write(std_out,*)" Writing nstates ",nstates
   gsub(:,1) = (/1,1/)
   gsub(:,2) = (/exc_size,nstates/)
   call slk_write(Slk_vec,"All",is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset,glob_subarray=gsub)

   call wrtout(std_out,' Calculating overlap matrix... ',"COLL")
   if (.not.do_full_diago) then
     MSG_ERROR(" Init of Slk_ovlp is wrong")
   end if

   call init_matrix_scalapack(Slk_ovlp,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc) 

   ! Calculate the overlap matrix.
   ! FIXME 
   ! The ESLL manual says that "matrices matrix1 and matrix2 must have no common elements; otherwise, results are unpredictable."
   ! However the official scaLAPACK documentation does not report this (severe) limitation.

   !call init_matrix_scalapack(Slk_tmp,exc_size,exc_size,Slk_processor,istwfk1,tbloc=tbloc) 
   !Slk_tmp%buffer_cplx = Slk_vec%buffer_cplx
   !call slk_pzgemm("C","N",Slk_tmp,cone,Slk_vec,czero,Slk_ovlp)
   !call destruction_matrix_scalapack(Slk_tmp)

   call slk_pzgemm("C","N",Slk_vec,cone,Slk_vec,czero,Slk_ovlp)

#ifdef DEV_MG_DEBUG_THIS
   ABI_ALLOCATE(exc_ham,(exc_size,exc_size))
   istat = ABI_ALLOC_STAT
   read(667)exc_ham

   write(std_out,*)"Error Ovlp: ",MAXVAL(ABS(exc_ham-Slk_ovlp%buffer_cplx))
   !Slk_ovlp%buffer_cplx = exc_ham
#endif

   !max_r=20; max_c=10
   !call print_arr(Slk_ovlp%buffer_cplx,max_r=max_r,max_c=max_c,unit=std_out)

   call destruction_matrix_scalapack(Slk_vec)

   call wrtout(std_out," Inverting overlap matrix... ","COLL")
   uplo="Upper"

#if 0
!DEBUG
   call xhdp_invert(uplo,Slk_ovlp%buffer_cplx,exc_size)

   !call slk_symmetrize(Slk_ovlp,uplo,"Hermitian")
   call hermitianize(Slk_ovlp%buffer_cplx,uplo)

   exc_ham = MATMUL(exc_ham,Slk_ovlp%buffer_cplx)
   do it=1,exc_size 
     exc_ham(it,it) = exc_ham(it,it) - cone
   end do

   write(std_out,*)"Error Inversion: ",MAXVAL(ABS(exc_ham))
   ABI_DEALLOCATE(exc_ham)
!END DEBUG

#else 
   ! call slk_zdhp_invert(Slk_ovlp,uplo)
   ! call hermitianize(Slk_ovlp%buffer_cplx,uplo)
   ! !call slk_symmetrize(Slk_ovlp,uplo,"Hermitian")

   call slk_zinvert(Slk_ovlp)  ! Version for generic complex matrix.
#endif

   if (allocated(exc_ham))  then
     ABI_DEALLOCATE(exc_ham)
   end if

#ifdef DEV_MG_DEBUG_THIS
   ABI_ALLOCATE(exc_ham,(exc_size,exc_size))
   istat = ABI_ALLOC_STAT
   read(668)exc_ham
   write(std_out,*)"Error in Inv Ovlp: ",MAXVAL(ABS(exc_ham-Slk_ovlp%buffer_cplx))

   !exc_ham = exc_ham-Slk_ovlp%buffer_cplx
   !do it=1,exc_size
   !  if ( MAXVAL(ABS(exc_ham(:,it))) > 0.1 ) write(std_out,*)"it: ",it,exc_ham(:,it)
   !end do

   !Slk_ovlp%buffer_cplx = exc_ham
   ABI_DEALLOCATE(exc_ham)

   !write(std_out,*)"MAX ERR",MAXVAL(ABS(Slk_ovlp%buffer_cplx - TRANSPOSE(DCONJG(Slk_ovlp%buffer_cplx))))
#endif
       
   call wrtout(std_out,' Writing overlap matrix S^-1 on file: '//TRIM(bseig_fname),"COLL")

   call slk_write(Slk_ovlp,"All",is_fortran_file,mpi_fh=mpi_fh,offset=ehdr_offset)
                                                                                                                
   call MPI_FILE_CLOSE(mpi_fh, mpi_err)               
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")

   call destruction_matrix_scalapack(Slk_ovlp)
                                               
   call end_scalapack(Slk_processor)
#else
   MSG_BUG("You should not be here!")
#endif

 END SELECT

 exc_gap    = MINVAL(ABS(exc_ene(1:nstates)))
 exc_maxene = MAXVAL(ABS(exc_ene(1:nstates)))

 write(msg,'(2(a,f7.2,2a))')&
&  " First excitonic eigenvalue: ",exc_gap*Ha_eV,   " [eV].",ch10,&
&  " Last  excitonic eigenvalue: ",exc_maxene*Ha_eV," [eV].",ch10
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")

 nene_printed = MIN(32*nsppol,nstates); if (prtvol>10) nene_printed = nstates
 write(msg,'(a,i4)')' Complex excitonic eigenvalues in eV up to n= ',nene_printed 
 call wrtout(std_out,msg,"PERS")

 do it=0,(nene_printed-1)/4
   write(msg,'(4f10.5)') ( exc_ene(ii)*Ha_eV, ii=1+it*4,MIN(it*4+4,nene_printed) )
   call wrtout(std_out,msg,"COLL")
 end do

 ABI_DEALLOCATE(exc_ene)

10 call xbarrier_mpi(comm)

end subroutine exc_diago_coupling_hegv
!!***
