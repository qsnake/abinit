!{\src2tex{textfont=tt}}
!!****f* ABINIT/outwant
!! NAME
!! outwant
!!
!! FUNCTION
!! This routine creates an output file containing all the
!! information needed to run WanT as a post-processing program
!! The resulting file is 'launch.dat'.
!!
!! The routine writes to the disk (unformatted file unitwnt) the following informations:
!!     alat - lattice parameter
!!     rprim - primitive translation vectors
!!     ntypat - nr of atom types in elementary cell
!!     tpat - nr of types of atoms in the elementary cell
!!     xcart - cartesian coordinates of the atoms in the elem. cell
!!     ecut - energy cut-off
!!     mband - # of bands taken in calculation (same for each K-pt)
!!     nk(3) - # of k-pts for each direction (uniform grid in the WHOLE BZ)
!!     s0(3) - the origin of the K-space
!!     kg_tmp(3,mpw*mkmem ) - reduced planewave coordinates
!!     imax - Maximum index  of a G vector among all k points (see explanation bellow)
!!     nkpt - total no of K-pts
!!     nsppol - nr of spin polarisations (1 or 2)
!!     eig(mband, nkp_tot) - eigenvalues/band/K_point
!!     ngfft(3) - nr of points used for FFT in each direction
!!     wfc(i)- cmplx(cg(1,i),cg(2,i)) - wavefunction
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (CMorari)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig(mband*nkpt*nsppol) = array for holding eigenvalues (Hartree)
!!  cg(2,mcg) = planewave coefficients of wavefunction
!!  kg(3, mpw*mkmem) = reduced planewave coordinates
!!  npwarr(nkpt) = number of planewaves in basis at this k-point
!!  mband = maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  nkpt = number of k - points
!!  nsppol = 1 for unpolarized, 2 for spin polarized
!!  nspinor = number of spinorial components of the wavefunction (on current proc)
!!  mkmem = number of k points which can fit in memory; set to 0 if use disk
!!  mpw = maximum dimensioned size of npw
!!  wff = information about wf disk file (needed for  mkmem=0 case )
!!  prtwant = if set to 1, print 0 in S0 output
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      clsopn,hdr_skip,leave_new,matr3inv,rwwf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outwant(dtfil,dtset,eig,cg,kg,npwarr,mband,mcg,mpi_enreg,nkpt,nsppol,&
                   mkmem,mpw,wff,prtwant)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outwant'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: mband,mcg,mkmem,mpw,nkpt,nsppol,prtwant
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(wffile_type),intent(inout) :: wff
!arrays
 integer :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp) :: cg(2,mcg),eig(mband*nkpt*nsppol)

!Local variables-------------------------------
! the following variables are not used; they are written to 'launch.dat'
! in order to be compatible with the WANT format
!scalars
 integer :: bandtot,i,icount,ierr,ifind,ig,ii,iij,ik,ik_,imax
 integer :: index,index1,ispin_,isppol,iunit,iwf,iwf_k,j,k
 integer :: maxat,mcg_disk,nbd_disk,ngm,ngw_,nk_,nkp,npw,npw_k,nspin_
 integer :: unitwnt
 real(dp) :: alat,scal,scal_,tt
 logical :: twrite=.true.
 character(len=20) :: section_name
 character(len=3) :: nameat
 character(len=500) :: message
 character(len=fnlen) :: filewnt
!arrays
 integer :: ikg(3),nk(3)
 integer,allocatable :: iwfi(:,:),kg_k(:,:),kg_tmp(:,:),tpat(:)
 real(dp) :: drprim(3,3),gmat(3,3),gmod(3),s0(3),t1(3),t2(3)
 real(dp),allocatable :: cg_k(:,:),eig_k(:),occ(:),xcoord(:,:,:)
 complex,allocatable :: wfc(:)

! ***************************************************************************

!WARNING: not tested for nsppol,nspinor >1
!
!Initialisations
 nameat = ' '
 bandtot=mband*nkpt*nsppol*dtset%nspinor
 filewnt='launch.dat'
 unitwnt=47

!Open the file
 open (unitwnt, file=filewnt, form='unformatted', status='unknown')
 write(message,'(3a)')ch10,' Opening file for WanT input: ',filewnt
 call wrtout(std_out,message,'COLL')

!Comments
 if(prtwant>1) then
   write(std_out,*) 'Wrong value for prtwant. Reseting to 1'
   prtwant=1
 elseif(prtwant==1) then
   do i=1,3
     s0(i)=0._dp
   end do
 end if

!Discussion of 'alat' ABINIT/ WanT
 if(dtset%acell_orig(1,1)==dtset%acell_orig(2,1).and.&
 dtset%acell_orig(1,1)==dtset%acell_orig(3,1)) then
   alat=dtset%acell_orig(1,1)
   do i=1,3
     do j=1,3
       drprim( i, j) = dtset%rprim_orig( i, j, 1 )
     end do
   end do
 else
!  Redefining the drprim( i, j)
   alat=dtset%acell_orig(1,1)
   do i=1,3
     do j=1,3
       drprim( i, j) = dtset%rprim_orig( i, j, 1 )*dtset%acell_orig(j, 1)/alat
     end do
   end do
 end if

!Now finding the no of k-pt for each direction PARALEL with the
!generators of the first B.Z.
!First decide if we have the Gamma point in the list; its index in the list is ... index
 nk(:)=1
 ifind=0
 icount=2
 do i=1,nkpt
   index1=0
   do j=1,3
     if(dtset%kptns(j,i)<tol8) index1=index1+1
   end do
   if(index1==3) then
     index=i
     ifind=1
     cycle
   end if
 end do
 if(ifind==0) then
   write(std_out,*) 'GAMMA POINT NOT IN THE LIST OF KPTS?'
   do ii=1,nkpt
     write(std_out,*) (dtset%kptns(j,ii),j=1,3)
   end do
   call leave_new ('COLL')
 end if

 call matr3inv(drprim,gmat)

!Modules for each vector in recip. space; nb: g(index coord, index point)
 do j=1,3
   gmod(j)=0.D0
   do i=1,3
     gmod(j)=gmod(j)+gmat(i,j)**2
   end do
   gmod(j)=sqrt(gmod(j))
 end do
 if(nkpt==2) then
   do j=1,3
     do ii=1,3
       t1(ii)=dtset%kptns(ii,1)-dtset%kptns(ii,2)
     end do
     tt=0._dp
     do iij=1,3
       t2(iij)=0._dp
       do ii=1,3
         t2(iij)=t2(iij)+t1(ii)*gmat(ii,iij)
       end do
       tt=tt + t2(iij)**2
     end do
     tt=sqrt(tt)
     scal=0._dp
     do ii=1,3
       scal=scal+t2(ii)*gmat(j,ii)
     end do
     scal=abs(scal)
!    Compare scal(tt,gmat) with simple product of modules -> paralel or not
     if(abs(scal-tt*gmod(j))<tol8) nk(j)=2
   end do

 elseif(nkpt>2) then

   do i=1,nkpt
     if(i.ne.index) then
       do ii=1,3
         t1(ii)=dtset%kptns(ii,index)-dtset%kptns(ii,i)
       end do
       tt=0._dp
       do iij=1,3
         t2(iij)=0._dp
         do ii=1,3
           t2(iij)=t2(iij)+t1(ii)*gmat(ii,iij)
         end do
         tt=tt + t2(iij)**2
       end do
       tt=sqrt(tt)
!      check for each direction in the BZ
       do j=1,3
         scal=0._dp
         do ii=1,3
           scal=scal+t2(ii)*gmat(j,ii)
         end do
         scal=abs(scal)
!        Compare scal(t1,gmat) with simple product of modules -> paralel or not
         if(abs(scal-tt*gmod(j))<tol8) nk(j)=nk(j)+1
       end do
     end if
   end do
 end if
 index=1
 do i=1,3
   index=index*nk(i)
 end do
 if(index.ne.nkpt) then
   write(std_out,*) 'OutwanT: Wrong assignemt of kpts', index,nkpt
   call leave_new('COLL')
 end if

!End counting/assigning no of kpts/direction
!Reordering the coordinates of all atoms - xcoord array
 ABI_ALLOCATE(tpat,(dtset%ntypat))
 tpat(:)=zero
 do i=1,dtset%natom
   do j=1,dtset%ntypat
     if(dtset%typat(i)==j) tpat(j)=tpat(j)+1
   end do
 end do
 maxat=maxval(tpat(:))
 ABI_ALLOCATE(xcoord,(3,maxat,dtset%ntypat))
 index=1
 do i=1, dtset%ntypat
   do k=1,tpat(i)
     do j=1,3
       xcoord(j,k,i)=dtset%xred_orig(j,index,1)
     end do
     index=index+1
   end do
 end do
!
!Defining the kg_tmp list
!Preparing the output of reduced coords., in a single list (kg_tmp(3,imax))
!We start with kg_tmp(:,i)=kg(:,i=1,npwarr(1)) then the new coordinates are added
!ONLY if they are not allready in the list. An index is associated
!for each kg_tmp which allow us to recover kg(3,mpw*nkpt) from
!the smaller list kg_tmp(3, imax)
 ABI_ALLOCATE(kg_tmp,(3,mpw*nkpt))
 ABI_ALLOCATE(iwfi,(nkpt,mpw))
 if(mkmem==0) then
   mcg_disk=mpw*dtset%nspinor*mband*nsppol
!  What do we need if the WF is stored on disk
   nbd_disk=mband
   ABI_ALLOCATE(kg_k,(3,mpw))
   ABI_ALLOCATE(eig_k,(2*mband))
   ABI_ALLOCATE(cg_k,(2,mcg_disk))
   ABI_ALLOCATE(occ,(mband))
   rewind(unit=dtfil%unkg)
 end if
 kg_tmp(:,:)=zero
 iwfi(:,:)=zero
 imax=npwarr(1)
 index=0

!if(mkmem>0)  write(std_out,*) (kg(1,j),j=npwarr(1)-10, npwarr(1)), 'unu'
!if(mkmem>0)  write(std_out,*) (kg(1,j),j=npwarr(2)-10, npwarr(2)), 'doi'
!if(mkmem>0)  write(std_out,*) (kg(1,j),j=npwarr(3)-10, npwarr(3)), 'trei'
!iunit=90
 do i=1,  nkpt
   if(i>1) then
     index=index+npwarr(i-1)
   end if
   if(mkmem==0) then
!    this will read the kg point(k)-by-point(k)
!    npw_k is npwarr(i)
     read(dtfil%unkg) npw_k
     read(dtfil%unkg)
     read (dtfil%unkg) kg_k(1:3,1:npw_k)
   end if
   do j=1, npwarr(i)
     if(i.eq.1) then
       iwfi(i,j)=j
       if(mkmem>0) kg_tmp(:,j)=kg(:,j)
       if(mkmem==0) kg_tmp(:,j)=kg_k(:,j)
     else
       ifind=0
       if(mkmem>0) ikg(:)=kg(:,index+j)
       if(mkmem==0) ikg(:)=kg_k(:,j)

       do k=1,imax
         if(ikg(1)==kg_tmp(1,k)) then
           if(ikg(2)==kg_tmp(2,k)) then
             if(ikg(3)==kg_tmp(3,k)) then
               ifind=1
               iwfi(i,j)=k
             end if
           end if
         end if
       end do

       if(ifind==0) then
         imax=imax+1
         kg_tmp(:,imax)=ikg(:)
         iwfi(i,j)=imax
       end if
     end if
   end do
 end do
 ngm=imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PART ONE: writing the header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(unitwnt) alat
 write( unitwnt ) ( drprim( i, 1 ), i = 1, 3 )  ! save A1
 write( unitwnt ) ( drprim( i, 2 ), i = 1, 3 )  ! save A2
 write( unitwnt ) ( drprim( i, 3 ), i = 1, 3 )  ! save A3
!write(std_out,* ) ( drprim( i, 1 ), i = 1, 3 )  ! save A1
!write(std_out,* ) ( drprim( i, 2 ), i = 1, 3 )  ! save A2
!write(std_out,* ) ( drprim( i, 3 ), i = 1, 3 )  ! save A3

 write(unitwnt) dtset%ntypat
!write(std_out,*) dtset%ntypat, 'NTYPAT', tpat

 do i = 1, dtset%ntypat
   write(unitwnt) tpat(i), nameat
   write(unitwnt) ((xcoord(j,k,i),j=1,3), k=1, tpat(i))
!  write(std_out,*) tpat(i), nameat
!  write(std_out,*) ((xcoord(j,k,i),j=1,3),k=1,tpat(i)), 'XCART'
 end do
 ABI_DEALLOCATE(xcoord)

!energy cut-off in Rydberg (WANT option)
 write (unitwnt) 2._dp*dtset%ecut, mband
!write(std_out,*)   2._dp*dtset%ecut, mband
 write (unitwnt) ( nk(i), i = 1, 3 ), ( s0(j), j = 1, 3 ),ngm
!write(std_out,*) ( nk(i), i = 1, 3 ), ( s0(j), j = 1, 3 ),imax
 write (unitwnt) ( kg_tmp( 1, i ), kg_tmp( 2, i ), kg_tmp( 3, i ), i = 1, ngm )
 write (unitwnt) mpw, mband, dtset%nkpt/dtset%nsppol
!write(std_out,*) mpw, mband,  dtset%nkpt/dtset%nsppol

 do i=1, nkpt
   write(unitwnt) (iwfi(i,j), j=1,mpw)
 end do
 ABI_DEALLOCATE(kg_tmp)

!Eigenvalues in HARTREE
 write (unitwnt)  ( eig( i ), i = 1, bandtot)
 write (unitwnt) ( npwarr( ik ), ik = 1, nkpt )
 write (unitwnt) ( mband, ik = 1, nkpt )
 write (unitwnt) (dtset%ngfft(i),i=1,3), imax, imax
!write(std_out,*)  ( eig( i ), i = 1, bandtot )
!write(std_out,*) ( npwarr( ik ), ik = 1, nkpt )
!write(std_out,*) ( mband, ik = 1, nkpt )
!write(std_out,*) (dtset%ngfft(i),i=1,3), imax ,imax
!a list with the band structure; usefull for 'windows' and 'disentangle' programs
!from WanT distribution
 iunit=77
 open(iunit, file='band.gpl', status='unknown')
 index=1
 do i=1,mband
   index=1
   do j=1,nkpt
     write(iunit,*) index, Ha_eV*eig(i+(j-1)*mband), eig(i+(j-1)*mband)
     index=index+1
   end do
   write(iunit,*)
 end do
 close(iunit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PART TWO: Writing the wavefunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!not used
 ngw_=0
 ik_=0
 nk_=0
 ispin_=0
 nspin_=0
 scal_=1._dp
!!!!!!!!!!!!!!!!!!!!!!!!!!
 iwf = 1
 iwf_k=1
 ABI_ALLOCATE(wfc,(imax))

 if(mkmem==0) then
   call clsopn(wff)
   call hdr_skip(wff,ierr)
 end if

!Loop over k-pt
 do nkp=1,nkpt
!  Not relevant
   write(unitwnt) twrite, ik_, section_name
!  Only 'mband' is relevant here
   write(unitwnt) ngw_, mband, ik_, nk_, nk_,ispin_, nspin_, scal_
   write(unitwnt) imax
!  Not relevant
   write(unitwnt) twrite
!  Loop over bands
   if(mkmem==0) then
     iwf_k=1
     npw=npwarr(nkp)
     call rwwf(cg_k,eig_k,0,0,0,nkp,isppol,kg_k,mband,mcg_disk,mpi_enreg,nbd_disk,nbd_disk,&
&     npw,dtset%nspinor,occ,-2,0,0,wff)
!    icg=icg+mcg_disk -> not needed since there is no rewind of unit wff
   end if

!  Preparing WF
   do k=1,mband
     if(mkmem >0) then
       wfc(:)=zero
!      From cg to wf:
       do i=iwf, iwf+npwarr(nkp)-1
         index=i-iwf+1
         wfc(iwfi(nkp,index))=cmplx(cg(1,i), cg(2,i), kind(0._dp))
       end do
       iwf=iwf+npwarr(nkp)
     elseif(mkmem==0) then
       wfc(:)=zero
!      From cg to wf:
       do i=iwf_k, iwf_k+npwarr(nkp)-1
         index=i-iwf_k+1
         wfc(iwfi(nkp,index))=cmplx(cg_k(1,i), cg_k(2,i), kind(0._dp))
       end do
       iwf_k=iwf_k+npwarr(nkp)
     else
       write(std_out,*) 'Wrong mkmem in outwant'
       call leave_new ('COLL')
     end if
     write(unitwnt) (wfc(ig), ig=1,imax)
!    End loop over bands
   end do

!  Not relevant
   write(unitwnt) twrite
!  Not relevant
   do i=1,mband
     write(unitwnt) i
   end do

!  End loop over k-pts
 end do
 ABI_DEALLOCATE(wfc)
 close(unit=unitwnt)
 write(message,'(2a)') ' Closing file ',ch10
 call wrtout(std_out,message,'COLL')

!End
end subroutine outwant
!!***
