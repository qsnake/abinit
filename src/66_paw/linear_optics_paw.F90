!{\src2tex{textfont=tt}}
!!****p* ABINIT/linear_optics_paw
!! NAME
!! linear_optics_paw
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! linear susceptiblity using matrix elements <-i Nabla> obtained from a
!! PAW ground state calculation. It uses formula 17 from Gadoc et al,
!! Phys. Rev. B 73, 045112 (2006) together with a scissors correction. It uses
!! a Kramers-Kronig transform to compute the real part from the imaginary part, and
!! it will work on all types of unit cells. It outputs all tensor elements of
!! both the real and imaginary parts.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (VR, PGhosh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  filnam: base of file names to read data from
!!  mpi_enreg: mpi set up variable, not used in this code
!!
!! OUTPUT
!!  _real and _imag output files
!!
!! NOTES
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      hdr_clean,hdr_io,initmpi_seq,kramerskronig,matrginv,metric,wffclose
!!      wffopen,wffreadeigk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine linear_optics_paw(filnam,filnam_out,mpi_enreg_seq)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_wffile
 use m_header,  only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'linear_optics_paw'
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=fnlen),intent(in) :: filnam,filnam_out
 type(MPI_type),intent(inout) :: mpi_enreg_seq

!Local variables-------------------------------
!no_abirules
 integer :: accesswff,bantot,bdtot0_index,bdtot_index,fform0,fform1,formeig0,headform
 integer :: iband,ierr,ii,ikpt,iom,iout,isppol,isym,jband,jj,master,me,mband
 integer :: method,mom,nband1,nband_k,nkpt,nspinor,nsppol,nsym,occopt,only_check
 integer :: rdwr,spaceComm,tim_rwwf
 integer,allocatable :: nband(:),symrel(:,:,:)
 real(dp) :: del,dom,fij,gdelta,omin,omax,paijpbij(2),sciss,wij,ucvol
 real(dp) :: diffwp, diffwm
 real(dp) :: e2rot(3,3),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),rprimdinv(3,3),symd(3,3),symdinv(3,3)
 real(dp),allocatable :: e1(:,:,:),e2(:,:,:,:),epsilon_tot(:,:,:,:),eigen0(:),eigtmp(:),eig0_k(:)
 real(dp),allocatable :: eig0tmp(:),kpts(:,:),occ(:),occ_k(:),oml1(:),wtk(:)
 complex,allocatable :: eps_work(:)
!
 character(len=fnlen) :: filnam0,filnam1,filnam_gen
 character(len=6) :: codvsn
!
 type(hdr_type) :: hdr
 type(wffile_type) :: wff0,wff1
 real(dp),allocatable :: psinablapsi(:,:,:,:)

! *********************************************************************************

 DBG_ENTER("COLL")

!* Fake MPI_type for the sequential part.
!This routine should not be parallelized as communicating gbig and other
!tables takes more time than recalculating them in sequential.
 call initmpi_seq(MPI_enreg_seq)
 
!write(std_out,'(a)')' Give the name of the output file ...'
!read(5, '(a)') filnam_out
!write(std_out,'(a)')' The name of the output file is :',filnam_out

!Read data file
 open(15,file=filnam,form='formatted')
 rewind(15)
 read(15,*)
 read(15,'(a)')filnam_gen       ! generic name for the files
 filnam1=trim(filnam_gen)//'_OPT' ! nabla matrix elements file

!Open the Wavefunction and optic files
!These default values are typical of sequential use
 accesswff=IO_MODE_FORTRAN ; spaceComm=abinit_comm_serial ; master=0 ; me=0
 call WffOpen(accesswff,spaceComm,filnam1,ierr,wff1,master,me,11)
 read(11,iostat=ierr)codvsn,headform,fform1
 if ((ierr /=0).or.(fform1/=610)) then
   write(std_out,*)'format prior version 6.1'
   fform1=613
 end if
 call WffClose(wff1,ierr)

!Open the conducti optic files
 call WffOpen(accesswff,spaceComm,filnam1,ierr,wff1,master,me,11)

!Read the header from Ground state file
 rdwr=1
 if (fform1==613) then
   filnam0=trim(filnam_gen)//'_WFK'
   call WffOpen(accesswff,spaceComm,filnam0,ierr,wff0,master,me,10)
   call hdr_io(fform0,hdr,rdwr,wff0)
 else 
   call hdr_io(fform1,hdr,rdwr,wff1)
 end if

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 nkpt=hdr%nkpt
 ABI_ALLOCATE(kpts,(3,nkpt))
 ABI_ALLOCATE(wtk,(nkpt))
 kpts(:,:)=hdr%kptns(:,:)
 wtk(:)=hdr%wtk(:)
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 rprimdinv(:,:) = rprimd(:,:)
 call matrginv(rprimdinv,3,3) ! need the inverse of rprimd to symmetrize the tensors
 ABI_ALLOCATE(nband,(nkpt*nsppol))
 ABI_ALLOCATE(occ,(bantot))
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 nsym=hdr%nsym
 ABI_ALLOCATE(symrel,(3,3,nsym))
 symrel(:,:,:)=hdr%symrel(:,:,:)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

!get ucvol etc.
 iout = -1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimdinv         =',rprimdinv(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimdinv(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimdinv(1:3,3)
 write(std_out,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband

!get eigen0
 if (fform1==613) then
!  Prepare the reading of Wff files
   formeig0=0 ; tim_rwwf=0
   ABI_ALLOCATE(eigtmp,(2*mband*mband))
   ABI_ALLOCATE(eig0tmp,(mband))
!  Read the eigenvalues of ground-state
   ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
   bdtot0_index=0 ; bdtot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband1=nband(ikpt+(isppol-1)*nkpt)
       call WffReadEigK(eig0tmp,formeig0,headform,ikpt,isppol,mband,mpi_enreg_seq,nband1,tim_rwwf,wff0)
       eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)
       bdtot0_index=bdtot0_index+nband1
     end do
   end do
   call WffClose(wff0,ierr)
   ABI_DEALLOCATE(eig0tmp)
 else
   ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
   read(11)(eigen0(iband),iband=1,mband*nkpt*nsppol)
 end if

 read(15,*)sciss
 read(15,*)dom,omin,omax,mom
 close(15)
 ABI_ALLOCATE(oml1,(mom))
 ABI_ALLOCATE(e1,(3,3,mom))
 ABI_ALLOCATE(e2,(2,3,3,mom))
 ABI_ALLOCATE(epsilon_tot,(2,3,3,mom))
 ABI_ALLOCATE(eps_work,(mom))
 del=(omax-omin)/(mom-1)
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do
 write(std_out,'(a,i8,4f10.5,a)')' npts,omin,omax,width,sciss      =',mom,omin,omax,dom,sciss,' Ha'

 ABI_ALLOCATE(psinablapsi,(2,3,mband,mband))

!loop over spin components
 do isppol=1,nsppol
   bdtot_index = 0
!  loop over k points
   do ikpt=1,nkpt
!    
!    number of bands for this k point
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
!    eigenvalues for this k-point
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!    occupation numbers for this k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    values of -i*nabla matrix elements for this k point
     psinablapsi=zero
     read(11)((psinablapsi(1:2,1,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(11)((psinablapsi(1:2,2,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(11)((psinablapsi(1:2,3,iband,jband),iband=1,nband_k),jband=1,nband_k) 

!    occupation numbers for k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    accumulate e2 for this k point, Eq. 17 from PRB 73, 045112 (2006)
     do iband = 1, nband_k
       do jband = 1, nband_k
         fij = occ_k(iband) - occ_k(jband) !occ number difference
         wij = eig0_k(iband) - eig0_k(jband) !energy difference
         if (abs(fij) > zero) then ! only consider states of differing occupation numbers
           do ii = 1, 3
             do jj = 1, 3
               paijpbij(1) = psinablapsi(1,ii,iband,jband)*psinablapsi(1,jj,iband,jband) + &
&               psinablapsi(2,ii,iband,jband)*psinablapsi(2,jj,iband,jband)
               paijpbij(2) = psinablapsi(2,ii,iband,jband)*psinablapsi(1,jj,iband,jband) - &
&               psinablapsi(1,ii,iband,jband)*psinablapsi(2,jj,iband,jband)
               do iom = 1, mom
!                original version
!                diffw = wij + sciss - oml1(iom) ! apply scissors term here
!                gdelta = exp(-diffw*diffw/(4.0*dom*dom))/(2.0*dom*sqrt(pi)) ! delta fnc resolved as Gaussian
!                e2(1,ii,jj,iom) = e2(1,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(1)*gdelta/(oml1(iom)*oml1(iom))
!                e2(2,ii,jj,iom) = e2(2,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(2)*gdelta/(oml1(iom)*oml1(iom))
                 diffwm = wij - sciss + oml1(iom) ! apply scissors term here
                 diffwp = wij + sciss - oml1(iom) ! apply scissors term here
                 gdelta = exp(-diffwp*diffwp/(4.0*dom*dom))/(2.0*dom*sqrt(pi))
                 e2(1,ii,jj,iom) = e2(1,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(1)*gdelta/(wij*wij)
                 e2(2,ii,jj,iom) = e2(2,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(2)*gdelta/(wij*wij)
               end do ! end loop over spectral points
             end do ! end loop over jj = 1, 3
           end do ! end loop over ii = 1, 3
         end if ! end selection on fij /= 0
       end do ! end loop over jband
     end do ! end loop over iband

     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(occ_k)
     bdtot_index=bdtot_index+nband_k
   end do ! end loop over k points
 end do ! end loop over spin polarizations

!here apply nsym symrel transformations to reconstruct full tensor from IBZ part
 epsilon_tot(:,:,:,:) = zero
 do isym = 1, nsym
   symd(:,:)=matmul(rprimd(:,:),matmul(symrel(:,:,isym),rprimdinv(:,:)))
   symdinv(:,:)=symd(:,:)
   call matrginv(symdinv,3,3)
   do iom = 1, mom
     e2rot(:,:)=matmul(symdinv(:,:),matmul(e2(1,:,:,iom),symd(:,:)))
     epsilon_tot(2,:,:,iom) = epsilon_tot(2,:,:,iom)+e2rot(:,:)/nsym
   end do
 end do

!generate e1 from e2 via KK transforma
 method=0 ! use naive integration ( = 1 for simpson)
 only_check=0 ! compute real part of eps in kk routine
 do ii = 1, 3
   do jj = 1, 3
     eps_work(:) = cmplx(0.0,epsilon_tot(2,ii,jj,:))
     call kramerskronig(mom,oml1,eps_work,method,only_check)
     epsilon_tot(1,ii,jj,:) = real(eps_work(:))
     if (ii /= jj) epsilon_tot(1,ii,jj,:) = epsilon_tot(1,ii,jj,:)- 1.0
   end do ! end loop over jj
 end do ! end loop over ii

 open(18,file=trim(filnam_out)//'_imag',form='formatted')
 open(19,file=trim(filnam_out)//'_real',form='formatted')

 write(18,'(a12,6a13)')' # Energy/Ha ','eps_2_xx','eps_2_yy','eps_2_zz',&
& 'eps_2_yz','eps_2_xz','eps_2_xy'
 write(19,'(a12,6a13)')' # Energy/Ha ','eps_1_xx','eps_1_yy','eps_1_zz',&
& 'eps_1_yz','eps_1_xz','eps_1_xy'

 do iom = 1, mom
   write(18,'(ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4)') oml1(iom),' ',&
&   epsilon_tot(2,1,1,iom),' ',epsilon_tot(2,2,2,iom),' ',epsilon_tot(2,3,3,iom),' ',&
&   epsilon_tot(2,2,3,iom),' ',epsilon_tot(2,1,3,iom),' ',epsilon_tot(2,1,2,iom)
   write(19,'(ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4)') oml1(iom),' ',&
&   epsilon_tot(1,1,1,iom),' ',epsilon_tot(1,2,2,iom),' ',epsilon_tot(1,3,3,iom),' ',&
&   epsilon_tot(1,2,3,iom),' ',epsilon_tot(1,1,3,iom),' ',epsilon_tot(1,1,2,iom)
 end do

 close(18)
 close(19)

 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(oml1)
 ABI_DEALLOCATE(e2)
 ABI_DEALLOCATE(e1)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(psinablapsi)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(kpts)

 call hdr_clean(hdr)

 DBG_EXIT("COLL")

 end subroutine linear_optics_paw
!!***
