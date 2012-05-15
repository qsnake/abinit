!{\src2tex{textfont=tt}}
!!****f* ABINIT/eltfrkin3
!! NAME
!! eltfrkin3
!!
!! FUNCTION
!! Compute the frozen-wavefunction kinetic enegy contribution to the
!! elastic tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DRH, DCA, XG, GM, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>
!!            =Fourier coefficients of wavefunction
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha) (NOT NEEDED !)
!!  effmass=effective mass for electrons (1. in common case)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  kptns(3,nkpt)=coordinates of k points in terms of reciprocal space
!!   primitive translations
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimension for number of planewaves
!!  nband(nkpt*nsppol)=number of bands being considered per k point
!!  nkpt=number of k points
!!  ngfft(18)=contain all needed information about 3D FFT, i
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2)
!!    at each k point
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  unkg=unit number for (k+G) sphere data
!!  wfftgs=struct info for disk file containing GS wavefunctions if mkmem==0
!!  wtk(nkpt)=k point weights
!!
!! OUTPUT
!!  eltfrkin(6,6)=non-symmetrized kinetic energy contribution to the
!!                    elastic tensor
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      d2kindstr2,hdr_skip,leave_test,metric,rdnpw,rwwf,sphereboundary,timab
!!      xcomm_world,xdefineoff,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eltfrkin3(cg,eltfrkin,ecut,ecutsm,effmass,&
&  istwfk,kg,kptns,mband,mgfft,mkmem,mpi_enreg,&
&  mpw,nband,nkpt,ngfft,npwarr,nspinor,nsppol,occ,&
&  rprimd,unkg,wfftgs,wtk)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eltfrkin3'
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_72_response, except_this_one => eltfrkin3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mkmem,mpw,nkpt,nspinor,nsppol,unkg
 real(dp),intent(in) :: ecut,ecutsm,effmass
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wfftgs
!arrays
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),rprimd(3,3),wtk(nkpt)
 real(dp),intent(out) :: eltfrkin(6,6)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,formeig,iband,icg,ierr,ii,ikg
 integer :: ikpt,index,ipw,isppol,istwf_k,jj,master,mcg_disk,me,n1,n2
 integer :: n3,nband_k,nkinout,npw_k,nspinor_,spaceComm
 integer :: tim_rwwf
 real(dp) :: ucvol
!arrays
 integer,allocatable :: gbound(:,:),kg_dum(:,:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),cwavef(:,:),eig_dum(:),ekinout(:)
 real(dp),allocatable :: eltfrkink(:,:),occ_dum(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' eltfrkin3 : enter '
!ENDDEBUG

!Default for sequential use
 master=0
 call xme_init(mpi_enreg,me)
!Init mpi_comm
 call xcomm_world(mpi_enreg,spaceComm)

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if (mkmem==0) then

!  Read wavefunction file header
   call hdr_skip(wfftgs,ierr)

!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wfftgs,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)

   mcg_disk=mpw*nspinor*mband
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))

 end if

 eltfrkin(:,:)=0.0_dp
 bdtot_index=0
 icg=0

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(cwavef,(2,mpw*nspinor))
 ABI_ALLOCATE(eltfrkink,(6,6))

!Define k-points distribution

!LOOP OVER SPINS
 do isppol=1,nsppol


!  Rewind kpgsph data file if needed:
   if (mkmem==0) rewind(unkg)

   ikg=0

!  Loop over k points
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)

     if(mpi_enreg%paral_compil_kpt==1)then
!      Skip this k-point if not the proper processor
!      BEGIN TF_CHANGES
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&       - me))/=0) then
!        END TF_CHANGES
         bdtot_index=bdtot_index+nband_k
         cycle
       end if
     end if

     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     kpoint(:)=kptns(:,ikpt)

     kg_k(:,:) = 0
     if (mkmem==0) then

       nspinor_=nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor_,0,unkg)

!      Read k+g data
       read (unkg) ((kg_k(ii,ipw),ii=1,3),ipw=1,npw_k)

       call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

!      Read the wavefunction block for ikpt,isppol
       tim_rwwf=14
       ABI_ALLOCATE(eig_dum,(mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(mband))
       call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&       npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wfftgs)
       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)

     else

!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(ikg,kg,kg_k,npw_k)
       do ipw=1,npw_k
         kg_k(1,ipw)=kg(1,ipw+ikg)
         kg_k(2,ipw)=kg(2,ipw+ikg)
         kg_k(3,ipw)=kg(3,ipw+ikg)
       end do
!      $OMP END PARALLEL DO

       call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

!      End if for choice governed by mkmem
     end if

     index=1+icg

     eltfrkink(:,:)=0.0_dp

     nkinout=6*6
     ABI_ALLOCATE(ekinout,(nkinout))
     ekinout(:)=zero

     do iband=1,nband_k

       if(mpi_enreg%paral_compil_kpt==1)then
!        BEGIN TF_CHANGES
         if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= me) then
!          END TF_CHANGES
           cycle
         end if
       end if

       if(mkmem/=0)then
         cwavef(:,1:npw_k*nspinor)=&
&         cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
       else
         cwavef(:,1:npw_k*nspinor)=&
&         cg_disk(:,1+(iband-1)*npw_k*nspinor:iband*npw_k*nspinor)
       end if

       call d2kindstr2(cwavef,ecut,ecutsm,effmass,ekinout,gmet,gprimd,&
&       istwf_k,kg_k,kpoint,npw_k,nspinor)

       eltfrkink(:,:)=eltfrkink(:,:)+ &
&       occ(iband+bdtot_index)* reshape(ekinout(:), (/6,6/) )

     end do !iband

     ABI_DEALLOCATE(ekinout)

     eltfrkin(:,:)=eltfrkin(:,:)+wtk(ikpt)*eltfrkink(:,:)

     ABI_DEALLOCATE(gbound)

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
!      Handle case in which kg, cg, are kept in core
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     end if

!    End loops on isppol and ikpt
   end do
 end do

!Fill in lower triangle
 do jj=2,6
   do ii=1,jj-1
     eltfrkin(jj,ii)=eltfrkin(ii,jj)
   end do
 end do

 if(mkmem==0)then
   ABI_DEALLOCATE(cg_disk)
 end if

 if(mpi_enreg%paral_compil_kpt==1)then
!  BEGIN TF_CHANGES
   call leave_test()
!  END TF_CHANGES
!  Accumulate eltfrkin on all proc.
   call timab(48,1,tsec)
   call xsum_mpi(eltfrkin,spaceComm,ierr)
   call timab(48,2,tsec)
 end if

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(eltfrkink)
 ABI_DEALLOCATE(kg_k)

!DEBUG
!write(std_out,*)' eltfrkin3 : exit '
!ENDDEBUG

end subroutine eltfrkin3
!!***
