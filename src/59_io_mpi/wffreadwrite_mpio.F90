!{src2tex{texfont=tt}}
!!****f* ABINIT/WffReadWrite_mpio
!! NAME
!! WffReadWrite_mpio
!!
!! FUNCTION
!!  This procedure read or write cg in the file _WFK using MPI_IO
!!  when cg are dispatched amoung commcart communicator
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MB,MD,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  wff=struct info for wavefunction
!!  nband_disk=number of bands on disk files to be write
!!  icg=shift to be given to the location of the cg array
!!  mcg=second dimention of cg
!!  mpi_enreg=information about parallelisation
!!  depl_mpi_to_seq=for each proc, index of cg in sequential mode
!!  npwso=npw*nspinor number of plane waves treated by this node.
!!  npwsotot=npwtot*nspinor Total number of planewaves Used to calculate the size of data to be written.
!!  rdwr=1 if reading, 2 if writing
!!
!! OUTPUT
!!  ierr=error status
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions,
!!
!! NOTES
!!  cg is written like the following:
!!    BeginMarker cg ( iband = 1 )  EndMarker
!!    BeginMarker cg ( iband = 2 )  EndMarker
!!    ...
!!    BeginMarker cg( iband = nband_disk ) EndMarker
!!
!!  BeginMarker and EndMarker give the value of the total length of cg for one band
!!
!!  For MPI-IO library the performance is improved by the use a "view" of the file for each proc.

!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!      mpi_file_close,mpi_file_open,mpi_file_read,mpi_file_read_all
!!      mpi_file_set_view,mpi_file_write,mpi_file_write_all,mpi_type_commit
!!      mpi_type_free,mpi_type_struct
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine WffReadWrite_mpio(wff,rdwr,cg,mcg,icg,nband_disk,npwso,npwsotot,depl_mpi_to_seq,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffReadWrite_mpio'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,mcg,nband_disk,npwso,npwsotot,rdwr
 integer,intent(out) :: ierr
 type(wffile_type),intent(inout) :: wff
!arrays
 integer,intent(in) :: depl_mpi_to_seq(npwso)
 real(dp),intent(inout) :: cg(2,mcg)

!Local variables-------------------------------
!scalars
#if defined HAVE_MPI_IO
 integer,parameter :: MAXBAND=500, check_markers=1
 integer :: filetype,iband,ibandmax,ibandmin,iblock,ii,iloc,ipw,jerr,jj,loc_depl_band
 integer :: nb,nband_block,step,totsize1bandByte,totsize1bandcg,wfftempo
 integer(kind=MPI_OFFSET_KIND) :: delim_record,offset
 character(len=500) :: msg
!arrays
 integer,allocatable :: BlockDepl(:),BlockLength(:),BlockType(:)
 integer,allocatable :: map(:),tempo_map(:)
 real(dp),allocatable :: buf(:),tempo_buf(:)
 integer*2,allocatable :: bufdelim2(:)
 integer*4,allocatable :: bufdelim4(:)
 integer*8,allocatable :: bufdelim8(:)
#if defined HAVE_FC_INT_QUAD
 integer*16,allocatable :: bufdelim16(:)
#endif
#endif

! *********************************************************************

 ierr=0

#if defined HAVE_MPI_IO
!----------------------------------------------
!! Prepare WF data
!----------------------------------------------
!Init offset of record
 wff%off_recs = wff%offwff

!Total size to be written (in number of bands and in bytes)
 totsize1bandcg=2*npwsotot
!call xsum_mpi(totsize1bandcg,wff%spaceComm_mpiio,ierr)

 totsize1bandByte=totsize1bandcg*wff%nbOct_dp+2*wff%nbOct_recMarker

!Check file size
 offset=wff%offwff+nband_disk*totsize1bandByte
 if (offset>Huge(offset)) then
   msg='File is too large for MPI-IO specifications !'
   MSG_ERROR(msg)
 end if

!Open file
 call MPI_FILE_OPEN(wff%spaceComm_mpiio,wff%fname,MPI_MODE_RDWR,MPI_INFO_NULL,wfftempo,ierr)

!----------------------------------------------------------
!Loop blocks of bands (to decrease offsets inside the file)
!----------------------------------------------------------
 ibandmax=0;ibandmin=1
 step=min(huge(check_markers)/totsize1bandByte,MAXBAND,nband_disk)
 do iblock=1,nband_disk/step+1
   ibandmax=min(ibandmin+step,nband_disk)
   nband_block=ibandmax-ibandmin+1
   offset=wff%offwff+(ibandmin-1)*totsize1bandByte

!  ----------------------------------------------
!  Read/Write bands
!  ----------------------------------------------

!  Build map; for better performance, map must be in increasing order
   ABI_ALLOCATE(map,(2*npwso*nband_block))
   ABI_ALLOCATE(buf,(2*npwso*nband_block))
   if (rdwr==1) then
!    If reading, only build map
     nb=0;loc_depl_band=0
     ABI_ALLOCATE(tempo_map,(2*npwso))
     do iband=ibandmin,ibandmax
       tempo_map(1:2*npwso)=-1
       jj=1;ipw=(iband-1)*npwso+icg
       do ii=1,npwso
         iloc=loc_depl_band+wff%nbOct_recMarker+2*(depl_mpi_to_seq(ii)-1)*wff%nbOct_dp
         tempo_map(jj  )=iloc              ! Real part
         tempo_map(jj+1)=iloc+wff%nbOct_dp ! Imag part
         jj=jj+2
       end do
       do ii=1,2*npwso ! Now, elimate holes
         if (tempo_map(ii)/=-1) then
           nb=nb+1
           map(nb)=tempo_map(ii)
         end if
       end do
       loc_depl_band=loc_depl_band+totsize1bandByte ! Location in bytes
     end do
   else if (rdwr==2) then
!    If writing, build map and store cg in a buffer
     nb=0;loc_depl_band=0
     ABI_ALLOCATE(tempo_map,(2*npwso))
     ABI_ALLOCATE(tempo_buf,(2*npwso))
     do iband=ibandmin,ibandmax
       tempo_map(1:2*npwso)=-1
       jj=1;ipw=(iband-1)*npwso+icg
       do ii=1,npwso
         iloc=loc_depl_band+wff%nbOct_recMarker+2*(depl_mpi_to_seq(ii)-1)*wff%nbOct_dp
         tempo_map(jj  )=iloc              ! Real part
         tempo_map(jj+1)=iloc+wff%nbOct_dp ! Imag part
         tempo_buf(jj:jj+1)=cg(1:2,ipw+ii)
         jj=jj+2
       end do
       do ii=1,2*npwso ! Now, elimate holes
         if (tempo_map(ii)/=-1) then
           nb=nb+1
           map(nb)=tempo_map(ii)
           buf(nb)=tempo_buf(ii)
         end if
       end do
       loc_depl_band=loc_depl_band+totsize1bandByte ! Location in bytes
     end do
     ABI_DEALLOCATE(tempo_map)
     ABI_DEALLOCATE(tempo_buf)
   end if  ! rdwr

!  Build and commit MPI datatype
   ABI_ALLOCATE(BlockLength,(nb+2))
   ABI_ALLOCATE(BlockDepl,(nb+2))
   ABI_ALLOCATE(BlockType,(nb+2))
   BlockLength(1)=1;BlockDepl(1)=0;BlockType(1)=MPI_LB
   do ii=2,nb+1
     BlockLength(ii)=1
     BlockDepl(ii)=map(ii-1)
     BlockType(ii)=MPI_DOUBLE_PRECISION
   end do
   BlockLength(nb+2)=1;BlockDepl(nb+2)=totsize1bandByte*nband_block;BlockType(nb+2)=MPI_UB
   call MPI_TYPE_STRUCT(nb+2,BlockLength,BlockDepl,BlockType,filetype,ierr)
   call MPI_TYPE_COMMIT(filetype,ierr)
   ABI_DEALLOCATE(BlockLength)
   ABI_DEALLOCATE(BlockDepl)
   ABI_DEALLOCATE(BlockType)

!  Read/Write data on disk
   call MPI_FILE_SET_VIEW(wfftempo,offset,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
   if (rdwr==1) then
     call MPI_FILE_READ_ALL (wfftempo,buf,nb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
   else
     call MPI_FILE_WRITE_ALL(wfftempo,buf,nb,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
   end if

!  In case of reading, retrieve cg
   if (rdwr==1) then
     nb=0;loc_depl_band=0
     ABI_ALLOCATE(tempo_buf,(2*npwso))
     do iband=ibandmin,ibandmax
       do ii=1,2*npwso ! Now, elimate holes
         if (tempo_map(ii)/=-1) then
           nb=nb+1;tempo_buf(ii)=buf(nb)
         end if
       end do
       jj=1;ipw=(iband-1)*npwso+icg
       do ii=1,npwso
         iloc=loc_depl_band+wff%nbOct_recMarker+2*(depl_mpi_to_seq(ii)-1)*wff%nbOct_dp
         cg(1:2,ipw+ii)=tempo_buf(jj:jj+1)
         jj=jj+2
       end do
       loc_depl_band=loc_depl_band+totsize1bandByte ! Location in bytes
     end do
     ABI_DEALLOCATE(tempo_map)
     ABI_DEALLOCATE(tempo_buf)
   end if ! rdwr

!  Free memory
   ABI_DEALLOCATE(map)
   ABI_DEALLOCATE(buf)
   call MPI_TYPE_FREE(filetype,ierr)

!  ----------------------------------------------
!  Check/Write record markers (only master proc)
!  ----------------------------------------------
   if ((rdwr==1.and.check_markers==1).or.(rdwr==2)) then

!    Define view for the file
     nb=2*nband_block
     ABI_ALLOCATE(BlockLength,(nb+2))
     ABI_ALLOCATE(BlockDepl,(nb+2))
     ABI_ALLOCATE(BlockType,(nb+2))
     BlockLength(1)=1;BlockDepl(1)=0;BlockType(1)=MPI_LB
     jj=2
     do ii=1,nband_block
       BlockType(jj:jj+1)  =wff%marker_mpi_type
       BlockLength(jj:jj+1)=1
       BlockDepl(jj  )=(ii-1)*totsize1bandByte
       BlockDepl(jj+1)= ii   *totsize1bandByte-wff%nbOct_recMarker
       jj=jj+2
     end do
     BlockLength(nb+2)=1;BlockDepl(nb+2)=nband_block*totsize1bandByte;BlockType(nb+2)=MPI_UB
     call MPI_TYPE_STRUCT(nb+2,BlockLength,BlockDepl,BlockType,filetype,ierr)
     call MPI_TYPE_COMMIT(filetype,ierr)
     call MPI_FILE_SET_VIEW(wfftempo,offset,MPI_BYTE,filetype,"native",MPI_INFO_NULL,ierr)
     ABI_DEALLOCATE(BlockLength)
     ABI_DEALLOCATE(BlockDepl)
     ABI_DEALLOCATE(BlockType)

!    Read/Write all markers (depend on Fortran marker MPI type)
     if (wff%me_mpiio==0) then
       jerr=0;delim_record=totsize1bandByte-2*wff%nbOct_recMarker
       if (wff%nbOct_recMarker==4) then
         ABI_ALLOCATE(bufdelim4,(nb))
         if (rdwr==2) bufdelim4(:)=delim_record
         if (rdwr==1) then
           call MPI_FILE_READ (wfftempo,bufdelim4,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
           if (any(bufdelim4(:)/=delim_record)) jerr=1
         else
           call MPI_FILE_WRITE(wfftempo,bufdelim4,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
         end if
         ABI_DEALLOCATE(bufdelim4)
       else if (wff%nbOct_recMarker==8) then
         ABI_ALLOCATE(bufdelim8,(nb))
         if (rdwr==2) bufdelim8(:)=delim_record
         if (rdwr==1) then
           call MPI_FILE_READ (wfftempo,bufdelim8,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
           if (any(bufdelim8(:)/=delim_record)) jerr=1
         else
           call MPI_FILE_WRITE(wfftempo,bufdelim8,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
         end if
         ABI_DEALLOCATE(bufdelim8)
#if defined HAVE_FC_INT_QUAD
       else if (wff%nbOct_recMarker==16) then
         ABI_ALLOCATE(bufdelim16,(nb))
         if (rdwr==2) bufdelim16(:)=delim_record
         if (rdwr==1) then
           call MPI_FILE_READ (wfftempo,bufdelim16,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
           if (any(bufdelim16(:)/=delim_record)) jerr=1
         else
           call MPI_FILE_WRITE(wfftempo,bufdelim16,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
         end if
         ABI_DEALLOCATE(bufdelim16)
#endif
       else if (wff%nbOct_recMarker==2) then
         ABI_ALLOCATE(bufdelim2,(nb))
         if (rdwr==2) bufdelim2(:)=delim_record
         if (rdwr==1) then
           call MPI_FILE_READ (wfftempo,bufdelim2,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
           if (any(bufdelim2(:)/=delim_record)) jerr=1
         else
           call MPI_FILE_WRITE(wfftempo,bufdelim2,2*nband_block,wff%marker_mpi_type,MPI_STATUS_IGNORE,ierr)
         end if
         ABI_DEALLOCATE(bufdelim2)
       end if
       if (rdwr==1.and.jerr==1) then
         write(unit=msg,fmt='(2a)') 'Error when reading record markers of file ',trim(wff%fname)
         MSG_PERS_ERROR(msg)
       end if
     end if  ! me_mpiio=0

!    Free memory
     call MPI_TYPE_FREE(filetype,ierr)

   end if ! rdwr

!  -----------------------------------------
!  End loop on blocks of bands
!  -----------------------------------------
   if (ibandmax>=nband_disk) exit
   ibandmin=ibandmax+1
 end do

!-----------------------------------------
!End statements
!-----------------------------------------
!Close file
 call MPI_FILE_CLOSE(wfftempo,ierr)

!Update offset
 wff%offwff=wff%offwff+totsize1bandByte*nband_disk
#endif

#if !defined HAVE_MPI_IO
!Dummy check to avoid warning from compilers.
 ABI_UNUSED((/wff%accesswff,rdwr,size(cg),mcg,icg,nband_disk,npwso,depl_mpi_to_seq(1),npwsotot/))
#endif

end subroutine WffReadWrite_mpio
!!***
