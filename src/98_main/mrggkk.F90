!{\src2tex{textfont=tt}}
!!****p* ABINIT/mrggkk
!! NAME
!! mrggkk
!!
!! FUNCTION
!! This program merges a GS file and several 1WF or GKK files for
!! different q-vectors and perturbations.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!! GKK file structure is composed of header records and eigenvalue arrays,
!! in binary or ascii:
!!   GS header = hdr
!!   GS eigenvalues = eigen
!!   number of perturbations = ntot
!!   for each perturbation
!!      1WF header = hdr1
!!      1st order eigenvalues = eigen1
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,destroy_mpi_enreg,hdr_clean,hdr_io,herald,initmpi_seq
!!      leave_new,mpi_comm_rank,mpi_comm_size,nullify_mpi_enreg,rwwf,wrtout
!!      xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program mrggkk

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
 use m_build_info
#if defined HAVE_MPI2
 use mpi
#endif

 use m_header,          only : hdr_clean
 !use m_io_tools,        only : read_line

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mrggkk'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------

!Local variables-------------------------------
!scalars
 integer,parameter :: unit1wf=22,unitgkk=24,unitgs=21,unitout=23,tim_rwwf=0,optkg=0
 integer :: binascii,fform,formeig,headform,i1wf,icg,igkk,ikpt,ios,isppol,mband
 integer :: mpw,n1wf,nband_disk,ngkk,ntot,ntotgkk,option,rdwr
 integer :: iband, jband
 integer :: rdwrout
 integer :: ierr,ipos
 real(dp) :: tolgkk=tol6
 character(len=1),parameter :: comment="#"
 character(len=24) :: codename
 character(len=500) :: message
 character(len=fnlen) :: file1wf,filegkk,filegs,outfile
 type(MPI_type) :: mpi_enreg,mpi_enreg_seq
 type(hdr_type) :: hdr,hdr1
 type(wffile_type) :: wff_dum
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: cg(:,:),eigen(:),occ(:)

! *************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI : one should write a separate routine -init_mpi_enreg-
!for doing that !!
 call xmpi_init()

!Default for sequential use
 call nullify_mpi_enreg(mpi_enreg)
 mpi_enreg%world_comm=0
 mpi_enreg%world_group=0
 mpi_enreg%me=0
 mpi_enreg%nproc=1
 mpi_enreg%nproc_atom=1
 mpi_enreg%num_group_fft = 0 ! in some cases not initialized but referenced in xdef_comm.F90
 mpi_enreg%paral_compil=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%mode_para="n"
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0
 mpi_enreg%paral_spin = 0
!MG080916 If we want to avoid MPI preprocessing options, %proc_distr should be always allocated and
!set to mpi_enreg%me. In such a way we can safely test its value inside loops parallelized over k-points
!For the time being, do not remove this line since it is needed in outkss.F90.
!nullify(mpi_enreg%proc_distrb)
!nullify(mpi_enreg%bandfft_kpt,mpi_enreg%tab_kpt_distrib)

!Initialize MPI
#if defined HAVE_MPI

 mpi_enreg%world_comm=xmpi_world
 mpi_enreg%world_group=MPI_GROUP_NULL
 call MPI_COMM_RANK(xmpi_world,mpi_enreg%me,ierr)
 call MPI_COMM_SIZE(xmpi_world,mpi_enreg%nproc,ierr)
!write(std_out,*)' mrgscr : nproc,me=',mpi_enreg%nproc,mpi_enreg%me
 mpi_enreg%paral_compil=1
 mpi_enreg%paral_compil_respfn=0
 mpi_enreg%paral_level=2
#endif

!Signal MPI I/O compilation has been activated
#if defined HAVE_MPI_IO
 mpi_enreg%paral_compil_mpio=1
 if(mpi_enreg%paral_compil==0)then
   write(message,'(6a)') ch10,&
&   ' mrgscr : ERROR -',ch10,&
&   '  In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
&   '  Action : recompile your code with different CPP flags.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
#endif

!Initialize paral_compil_kpt, actually always equal to paral_compil
!(paral_compil_kpt should be suppressed after big cleaning)
 mpi_enreg%paral_compil_kpt=0
 if(mpi_enreg%paral_compil==1) mpi_enreg%paral_compil_kpt=1

!Other values of mpi_enreg are dataset dependent, and should NOT be initialized inside mrggkk.F90.

!* Init fake MPI type with values for sequential case.
 call initmpi_seq(MPI_enreg_seq)

!dummy wff datatype to feed to rwwf
 wff_dum%unwff = unitgs
 wff_dum%formwff = 0  !scalar eigenvalues
 wff_dum%accesswff = IO_MODE_FORTRAN
 wff_dum%kgwff = 0
 wff_dum%fname = ""

 codename='MRGGKK'//repeat(' ',18)

!write greating,read the file names, etc.
 call herald(codename,abinit_version,std_out)
!YP: calling dump_config() makes tests fail => commented
!call dump_config()

 write(message,'(17a)')&
& ' Files file format: ',ch10,ch10,&
& '  Name of the output file',ch10,&
& '  Integer flag: 0 --> binary output,   1 --> ascii formatted output',ch10,&
& '  Name of the groud state wavefunction file WF',ch10,&
& '  Number of 1WF, of GKK files, and number of 1WF files in all the GKK files',ch10,&
& '  Names of the 1WF files...',ch10,&
& '  Names of the GKK files...',ch10,ch10,&
& ' Enter name of output file: '
 call wrtout(std_out,message,'COLL')

!get file with filenames and number of 1wf files
 read(*,'(a)') outfile
 ipos=INDEX(outfile,comment)
 if (ipos/=0) outfile=outfile(:ipos-1)

 read(*,*) binascii

 read(*,'(a)') filegs
 ipos=INDEX(filegs,comment)
 if (ipos/=0) filegs=filegs(:ipos-1)

 read(*,*) n1wf,ngkk,ntotgkk

 write(message,'(7a,i4,2a,i4,2a,i4,a)')&
& ' Output                     = ',trim(outfile),ch10,&
& ' Ground State file          = ',trim(filegs),ch10,&
& ' Number of 1WF files        = ',n1wf,ch10,&
& ' Number of GKK files        = ',ngkk,ch10,&
& ' Total Number of 1WF in GKK = ',ntotgkk,ch10
 call wrtout(std_out,message,'COLL')


!output without rewinding the file
 if (binascii == 0) then
!  open output file
   open(unit=unitout,file=outfile,form='unformatted',iostat=ios)
   rdwrout = 6
 else if (binascii == 1) then
!  rdwrout=4 ! use for screen output and change writes of eigen to (*,*)
!  MJV 27/5/2008 removed 'new' constraint on gkk files: presume competent user!
   open(unit=unitout,file=outfile,form='formatted',iostat=ios)
   rdwrout = 4
 else if (binascii == 2) then
!  this is for simple "short" output of the matrices, without headers or imaginary part
   open(unit=unitout,file=outfile,form='formatted',iostat=ios)
   rdwrout = 4
 else
   MSG_ERROR(' binascii must be 0 or 1')
 end if

 if (ios/=0) then
   MSG_ERROR('opening file: '//trim(outfile)//' for output')
 end if

 rewind (unitout)

!-------------------------------------------------------
!now read and write information for GS file
!-------------------------------------------------------

!open GS wf file
 call wrtout(std_out,' normal input for GS file',"COLL")

 open(unit=unitgs,file=filegs,form='unformatted',status='old',iostat=ios)
 if (ios/=0) then
   MSG_ERROR('opening file: '//trim(filegs))
 end if

 rewind (unitgs)
!read header of GS wf file
 rdwr = 5
 call hdr_io(fform,hdr,rdwr,unitgs)
 ABI_CHECK(fform/=0,"fform==0")

!copy header of GS file to output
 if (binascii /= 2) call hdr_io(fform,hdr,rdwrout,unitout)

 call wrtout(std_out,' header echoed to output file',"COLL")

!retrieve GS eigenvalues from GS wf file and echo to output
 mband = maxval(hdr%nband)
 mpw = maxval(hdr%npwarr)
 ABI_ALLOCATE(cg,(2,mpw*hdr%nspinor*mband))
 ABI_ALLOCATE(eigen,(mband))
 ABI_ALLOCATE(kg_k,(3,0))
 ABI_ALLOCATE(occ,(mband))
 option = 1
 formeig = 0
 icg = 0
 headform=hdr%headform
 wff_dum%unwff = unitgs
 wff_dum%formwff = formeig !scalar eigenvalues

 do isppol=1,hdr%nsppol
   do ikpt=1,hdr%nkpt
     call rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,&
&     mband,mpw*hdr%nspinor*mband,mpi_enreg,hdr%nband(ikpt),nband_disk,&
&     hdr%npwarr(ikpt),hdr%nspinor,occ,option,optkg,tim_rwwf,wff_dum)
     if (binascii==0) then
       write(unitout) eigen(1:hdr%nband(ikpt))
     else
       write(unitout,*) eigen(1:hdr%nband(ikpt))
     end if
   end do
 end do

 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(occ)

!close GS wf file
 close (unitgs)
 call hdr_clean(hdr)

 ntot = n1wf + ntotgkk
 if (binascii==0) then
   write (unitout) ntot
 else
   write (unitout,*) ntot
 end if

!-------------------------------------------------------
!now read and write information for 1WF files
!-------------------------------------------------------

 formeig = 1
 wff_dum%unwff = unit1wf
 wff_dum%formwff = formeig !scalar eigenvalues

 do i1wf=1,n1wf
!  for each 1wf file, get name...
   read(*,'(a)') file1wf
   ipos=INDEX(file1wf,comment)
   if (ipos/=0) file1wf=file1wf(:ipos-1)

!  open 1wf file
   call wrtout(std_out,' normal input for 1WF file ',"COLL")

   open(unit=unit1wf,file=file1wf,form='unformatted',status='old',iostat=ios)
   ABI_CHECK(ios==0,'opening: '//trim(file1wf)//' as old')

   rewind (unit1wf)

!  read in header of _WF1 file
   rdwr = 5
   call hdr_io(fform,hdr1,rdwr,unit1wf)
   if (fform == 0) then
     write(message,'(a,i0,a)')' 1WF header number ',i1wf,' was mis-read. fform == 0'
     MSG_ERROR(message)
   end if

!  copy header of 1WF file to output
!  WARNING: cant use normal hdr_io because it rewinds the output file and
!  destroys previous content.
   if (binascii /= 2) then
     call hdr_io(fform,hdr1,rdwrout,unitout)
   else
     write (unitout,'(a,3E20.10)') "qpt ", hdr1%qptn
     write (unitout,'(a,I6)') "pertnum ", hdr1%pertcase
   end if

!  retrieve 1WF <psi_k+q | H | psi_k> from 1wf file and echo to output
   mband = maxval(hdr1%nband)
   mpw = maxval(hdr1%npwarr)
   ABI_ALLOCATE(cg,(2,mpw*hdr1%nspinor*mband))
   ABI_ALLOCATE(eigen,(2*mband*mband))
   ABI_ALLOCATE(kg_k,(3,0))
   ABI_ALLOCATE(occ,(mband))
   option = 1
   headform=hdr1%headform
   do isppol=1,hdr1%nsppol
     do ikpt=1,hdr1%nkpt
!      write(std_out,*) 'isppol,ikpt = ', isppol,ikpt
       call rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,&
&       mband,mpw*hdr1%nspinor*mband,mpi_enreg,hdr1%nband(ikpt),nband_disk,&
&       hdr1%npwarr(ikpt),hdr1%nspinor,occ,option,optkg,tim_rwwf,wff_dum)
       if (binascii==0) then
         write(unitout) eigen(1:2*hdr1%nband(ikpt)**2)
       else if (binascii==1) then
         write(unitout,*) eigen(1:2*hdr1%nband(ikpt)**2)
       else if (binascii==2) then
         do iband=1, hdr1%nband(ikpt)
           do jband=1, hdr1%nband(ikpt)
             if (abs(eigen(2*hdr1%nband(ikpt)*(iband-1)+2*(jband-1)+1))>tolgkk) then
               write(unitout,'(E18.7, 2x)', ADVANCE='NO') eigen(2*hdr1%nband(ikpt)*(iband-1)+2*(jband-1)+1)
             else
               write(unitout,'(I18, 2x)', ADVANCE='NO') 0
             end if 
             if (abs(eigen(2*hdr1%nband(ikpt)*(iband-1)+2*(jband-1)+2))>tolgkk) then
               write(unitout,'(E18.7, 2x)', ADVANCE='NO') eigen(2*hdr1%nband(ikpt)*(iband-1)+2*(jband-1)+2)
             else 
               write(unitout,'(I18, 2x)', ADVANCE='NO') 0
             end if 
           end do
           write(unitout,*)
         end do
         write(unitout,*)
       end if
     end do
     if (binascii==2) write(unitout,'(2a)') ch10, ch10
   end do

   ABI_DEALLOCATE(cg)
   ABI_DEALLOCATE(eigen)
   ABI_DEALLOCATE(kg_k)
   ABI_DEALLOCATE(occ)

   close (unit1wf)
!  clean header to deallocate everything
   call hdr_clean(hdr1)
 end do

!-------------------------------------------------------
!now read and write information for small GKK files
!-------------------------------------------------------
 formeig = 1
 do igkk=1,ngkk
!  for each gkk file, get name...
   read(*,'(a)') filegkk
   ipos=INDEX(filegkk,comment)
   if (ipos/=0) filegkk=filegkk(:ipos-1)

!  open gkk file
   call wrtout(std_out,' normal input for GKK file',"COLL")

   open(unit=unitgkk,file=filegkk,form='unformatted',status='old',iostat=ios)
   if (ios/=0) then
     MSG_ERROR('Opening: '//trim(filegkk)//' as old')
   end if

   rewind (unitgkk)

!  read in header of GS file and eigenvalues
   call hdr_io(fform,hdr,5,unitgkk)
!  
!  could force a comparison of header with global header above for consistency
!  

   ABI_ALLOCATE(eigen,(mband))
   call wrtout(std_out,'mrggkk : try to reread GS eigenvalues','COLL')

   do isppol=1,hdr%nsppol
     do ikpt=1,hdr%nkpt
       read (unitgkk,IOSTAT=ierr) eigen(1:hdr%nband(ikpt))
       if (ierr /= 0) write (std_out,*) 'error reading eigen from gkk file'
     end do
   end do
   read(unitgkk,IOSTAT=ierr) n1wf
   ABI_CHECK(ierr==0,'error reading eigen from gkk file')
   ABI_DEALLOCATE(eigen)

   ABI_ALLOCATE(eigen,(2*mband*mband))
   do i1wf=1,n1wf
!    read in header of 1WF file
     rdwr = 5
     call hdr_io(fform,hdr1,rdwr,unitgkk)
     if (fform == 0) then
       write(message,'(a,i0,a)')' 1WF header number ',i1wf,' was mis-read. fform == 0'
       MSG_ERROR(message)
     end if

!    copy header of 1WF file to output
     if (binascii /= 2) then
       call hdr_io(fform,hdr1,rdwrout,unitout)
     else
       write (unitout,'(a,3E20.10)') "qpt ", hdr1%qptn
       write (unitout,'(a,I6)') "pertnum ", hdr1%pertcase
     end if

!    retrieve 1WF <psi_k+q | H | psi_k> from gkk file and echo to output
     do isppol=1,hdr1%nsppol
       do ikpt=1,hdr1%nkpt
!        write(std_out,*) 'isppol,ikpt = ', isppol,ikpt
         read (unitgkk,IOSTAT=ierr) eigen(1:2*hdr1%nband(ikpt)**2)
         if (ierr /= 0) write (std_out,*) 'error reading eigen2 from gkk file',isppol,ikpt
         if (binascii==0) then
           write (unitout) eigen(1:2*hdr1%nband(ikpt)**2)
         else if (binascii==1) then
           write (unitout,*) eigen(1:2*hdr1%nband(ikpt)**2)
         else if (binascii==2) then
           do iband=1, hdr1%nband(ikpt)
             do jband=1, hdr1%nband(ikpt)
               if (abs(eigen(2*hdr1%nband(ikpt)*(iband-1)+2*(jband-1)+1))>tolgkk) then
                 write(unitout,'(E18.7, 2x)', ADVANCE='NO') eigen(2*hdr1%nband(ikpt)*(iband-1)+2*(jband-1)+1)
               else
                 write(unitout,'(I18, 2x)', ADVANCE='NO') 0
               end if 
               if (abs(eigen(2*hdr1%nband(ikpt)*(iband-1)+2*(jband-1)+2))>tolgkk) then
                 write(unitout,'(E18.7, 2x)', ADVANCE='NO') eigen(2*hdr1%nband(ikpt)*(iband-1)+2*(jband-1)+2)
               else 
                 write(unitout,'(I18, 2x)', ADVANCE='NO') 0
               end if 
             end do
             write(unitout,*)
           end do
           write(unitout,*)
         end if
       end do
       if (binascii==2) write(unitout,'(2a)') ch10, ch10
     end do
     call hdr_clean(hdr1)
   end do
!  end loop over 1wf segments in small gkk file
   ABI_DEALLOCATE(eigen)

   close (unitgkk)
   call hdr_clean(hdr)

 end do
!end loop over small gkk files

 close (unitout)

 write(message,'(2a)')ch10,' Done'
 call wrtout(std_out,message,'COLL')

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program mrggkk
!!***
