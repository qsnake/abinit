!{\src2tex{textfont=tt}}
!!****f* ABINIT/wfk_read_ene
!! NAME
!! wfk_read_ene
!!
!! FUNCTION
!!  Read eigenvalues and occupations numbers from the WFK file. 
!!  Only master opens and reads the file in parallel runs. Data are then
!!  MPI broadcasted inside comm.
!!
!! INPUTS
!!  wfk_fname=Name of the WFK file. 
!!  accesswff=Option defining the access mode. 
!!  comm=MPI Communicator.
!!  prtvol=Verbosity level.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Hdr<Hdr_type>=The abinit header.
!!
!! SIDE EFFECTS
!!  anbnds= 
!!    in:  Number of bands to be read. If <=0, all bands stored on file are returned. 
!!         Cannot be greater than the number of bands stored on file.
!!    out: Number of bands read. Gives the first dimension of energies_p.
!!  energies_p(anbnds,Hdr%nkpt,Hdr%nsppol)=
!!    in:  pointer to null() 
!!    out: arrays with energies.
!!
!! PARENTS
!!      setup_screening,setup_sigma
!!
!! CHILDREN
!!      hdr_io,hdr_io_etsf,initmpi_seq,rwwf,wffclose,wffopen,xbarrier_mpi
!!      xcast_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wfk_read_ene(wfk_fname,accesswff,anbnds,energies_p,Hdr,prtvol,comm) 

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors

 use m_io_tools,    only : get_unit
 use m_wffile,      only : wffile_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_read_ene'
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,comm,prtvol
 integer,intent(inout) :: anbnds
 character(len=fnlen),intent(in) :: wfk_fname
 type(Hdr_type),intent(out) :: Hdr
!arrays
 real(dp),pointer :: energies_p(:,:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_rwwf0=0,formeig0=0,headform0=0,icg0=0
 integer :: ierr,ik_ibz,spin,master,mpw,fform,rdwr 
 integer :: my_rank,wfk_unt,nprocs
 integer :: npw_k,nband_disk,nband_k,mband,my_mkmem,mcg,option,optkg
 logical :: master_bcast
 character(len=500) :: msg
 type(Wffile_type) :: Wff
 type(MPI_type) :: MPI_enreg_seq
!arrays
 real(dp),allocatable :: eig_k(:),occ_k(:),cg_k(:,:) 
 integer,allocatable :: kg_k(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (accesswff==IO_MODE_ETSF) then 
   write(msg,'(3a)')&
&   ' when accesswff==3, support for the ETSF I/O library ',ch10,&
&   ' must be compiled. Use --enable-etsf-io when configuring '
#if !defined HAVE_TRIO_ETSF_IO
   MSG_ERROR(msg)
#endif
 end if
 !
 ! === MPI ENVIRONMENT ===
 my_rank = xcomm_rank(comm)
 nprocs  = xcomm_size(comm)
 master=0
 !
 ! Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 !
 ! * Init Wff structure.
 wfk_unt = get_unit()
 call WffOpen(accesswff,comm,wfk_fname,ierr,Wff,master,my_rank,wfk_unt)
 ABI_CHECK(ierr==0,"Error in the initialization of Wff")
 !
 ! * Read the Header. Hdr is local to each node when hdr_io returns.
 rdwr=1
 if (ANY(Wff%accesswff == (/IO_MODE_FORTRAN, IO_MODE_FORTRAN_MASTER, IO_MODE_MPI/) )) then
   call hdr_io(fform,Hdr,rdwr,Wff)
 else if (Wff%accesswff==IO_MODE_ETSF) then
   call hdr_io_etsf(fform,Hdr,rdwr,Wff%unwff)
 end if
 ABI_CHECK(fform/=0,"hdr_io returned fform/=0")
 !
 ! Output the header of the GS wavefunction file.
 if (prtvol>0) call hdr_io(fform,Hdr,4,std_out) 
 !
 ! Allocate the pointer to be returned in output.
 mpw = MAXVAL(Hdr%npwarr); mband=MAXVAL(Hdr%nband)

 if (anbnds<=0) then
   anbnds=mband
 else if (anbnds>mband) then
   ABI_CHECK(anbnds<=MINVAL(Hdr%nband),"anbnds > MINVAL(nband)")
   anbnds=mband
 end if

 ABI_ALLOCATE(energies_p,(anbnds,Hdr%nkpt,Hdr%nsppol))

 master_bcast = (wff%accesswff==IO_MODE_FORTRAN_MASTER) !.or.(wff%accesswff==IO_MODE_MPI)

 if ( wff%accesswff==IO_MODE_FORTRAN .or. &
& (wff%accesswff==IO_MODE_FORTRAN_MASTER .and.wff%master==wff%me).or. &
& (wff%accesswff==IO_MODE_MPI)    ) then
!& (wff%accesswff==IO_MODE_MPI.and.wff%master==wff%me)    ) then

   do spin=1,Hdr%nsppol
     do ik_ibz=1,Hdr%nkpt

       npw_k     = Hdr%npwarr(ik_ibz) 
       nband_disk= Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)
       nband_k   = Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)
       mband     = nband_k  
       my_mkmem  = 0
       mcg       = npw_k*Hdr%nspinor*mband*my_mkmem

       ABI_ALLOCATE(eig_k,((2*mband)**formeig0*mband))
       ABI_ALLOCATE(occ_k,(mband))
       option=3 ! To read eigenvalues only.
       optkg=0  ! G-vectors not needed.

       ABI_ALLOCATE(cg_k,(2,mcg))
       ABI_ALLOCATE(kg_k,(3,optkg*npw_k))
       !
       ! Read the block of bands for this (k,s).
       call rwwf(cg_k,eig_k,formeig0,headform0,icg0,ik_ibz,spin,kg_k,mband,mcg,MPI_enreg_seq,nband_k,&
&         nband_disk,npw_k,Hdr%nspinor,occ_k,option,optkg,tim_rwwf0,Wff)

       energies_p(:,ik_ibz,spin) = eig_k(1:anbnds)

       ABI_DEALLOCATE(eig_k)
       ABI_DEALLOCATE(occ_k)
       ABI_DEALLOCATE(kg_k)
       ABI_DEALLOCATE(cg_k)

     end do !ik_ibz
   end do !spin
 end if
 !
 ! * Broadcast data if needed.
 if (master_bcast .and. nprocs>1) then
   call xcast_mpi(energies_p,master,comm,ierr)
   call xbarrier_mpi(comm)
 end if
 !
 ! * Close the file and do NOT delete it!
 call WffClose(Wff,ierr)
 ABI_CHECK(ierr==0,"Error in WffClose")

 DBG_EXIT("COLL")

end subroutine wfk_read_ene
!!***
