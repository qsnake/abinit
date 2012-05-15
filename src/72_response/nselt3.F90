!{\src2tex{textfont=tt}}
!!****f* ABINIT/nselt3
!! NAME
!! nselt3
!!
!! FUNCTION
!! This routine compute the non-stationary expression for the
!! second derivative of the total energy, wrt strain for a whole row of
!! mixed strain derivatives.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DRH, XG, DCA, GMR, MM, AR, MV, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass=effective mass for electrons (1. in common case)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the
!!     storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points in the reduced BZ
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpert =maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  maximum dimension for q points in grids for nonlocal form factors
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array. If /=0,
!!   the exchange-correlation kernel must be computed.
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!    perturbation
!!  ntypat=number of types of atoms in unit cell.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band
!!   and k in the reduced Brillouin zone (usually =2)
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  prtbbb=if 1, band-by-band decomposition (also dim of d2bbb)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rhor(nfft,nspden)=GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  type(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  unkg=unit number of k+G data file
!!  unkg1=unit number of k+G+q data file
!!  wffnow= struct info for wf disk file
!!  wfftgs=struct info r for ground-state wf disk file
!!  unylm=unit number for disk file containing Ylm(k) if mkmem==0
!!  unylm1=unit number for disk file containing Ylm(k+q) if mk1mem==0
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang)= real spherical harmonics for each
!!    G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang)= real spherical harmonics for each
!!    G and k+q point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical for each
!!    G and k point
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm)= gradients of real spherical for each
!!    G and k+g point
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!!  d2lo(2,3,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,3,mpert,3,mpert)=non-local contributions to the 2DTEs
!!
!! NOTES
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      clsopn,dotprod_vn,hartrestr,hdr_skip,leave_test,mkcor3,mkvxcstr3,nstwf4
!!      rdnpw,stresssym,vlocalstr,wrtout,xcomm_init,xdefineoff,xmaster_init
!!      xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine nselt3(atindx,atindx1,blkflg,cg,cg1,cplex,&
& d2bbb,d2lo,d2nl,ecut,ecutsm,effmass,&
& gmet,gprimd,gsqcut,idir,&
& ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mband,mgfft,&
& mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
& natom,nattyp,nband_rbz,nfft,ngfft,&
& nkpt_rbz,nkxc,nloalg,npwarr,npwar1,nspden,nspinor,nsppol,&
& nsym1,ntypat,occ_rbz,&
& paral_kgb, ph1d,prtbbb,psps,qphon,rhog,&
& rhor,rhor1,rmet,rprimd,symrc1,type,ucvol,&
& unkg,unkg1,&
& wffnow,wfftgs,unylm,unylm1,wtk_rbz,&
& xred,ylm,ylm1,ylmgr,ylmgr1)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nselt3'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_72_response, except_this_one => nselt3
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mband,mgfft,mk1mem
 integer,intent(in) :: mkmem,mpert,mpsang,mpw,mpw1,natom,nfft,nkpt_rbz
 integer,intent(in) :: nkxc,nspden,nspinor,nsppol,nsym1,ntypat
 integer,intent(in) :: paral_kgb,prtbbb,unkg,unkg1,unylm,unylm1
 real(dp),intent(in) :: ecut,ecutsm,effmass,gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow,wfftgs
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom)
 integer,intent(in) :: istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem),nattyp(ntypat)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),ngfft(18)
 integer,intent(in) :: nloalg(5),npwar1(nkpt_rbz),npwarr(nkpt_rbz)
 integer,intent(in) :: symrc1(3,3,nsym1),type(natom)
 integer,intent(out) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpt_rbz(3,nkpt_rbz),kxc(nfft,nkxc)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qphon(3),rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,nspden)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang)
 real(dp),intent(out) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: ban2tot,bantot,bd2tot_index,bdtot_index,formeig
 integer :: icg,icg1,idir1,ierr,ifft,ii,ikg,ikg1,ikpt
 integer :: ilm,ipert1,ispden,isppol,istr1,istwf_k
 integer :: master,mbd2kpsp,mbdkpsp,me,muig,n1,n2,n3,n3xccc,n4,n5,n6
 integer :: nband_k,nfftot,npw1_k,npw_k,nsp,option,spaceComm
 real(dp) :: doti,dotr
 real(dp) :: wtk_k
 character(len=500) :: message
!arrays
 integer :: ikpt_fbz(3)
 integer,allocatable :: kg1_k(:,:),kg_k(:,:)
 real(dp) :: kpoint(3),restr(6)
 real(dp),allocatable :: d2bbb_k(:,:,:,:),d2nl_k(:,:,:)
 real(dp),allocatable :: occ_k(:)
 real(dp),allocatable :: vhartr01(:),vpsp1(:),vxc1(:,:),xccc3d1(:),ylm1_k(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr1_k(:,:,:),ylmgr_k(:,:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' nselt3 : enter '
!stop
!ENDDEBUG

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)
!Init master
 call xmaster_init(mpi_enreg,master)

!Unit numbers

!Zero only portion of nonlocal matrix to be computed here
 d2nl(:,:,natom+3:natom+4,idir,ipert)=0.0_dp
 bdtot_index=0
 bd2tot_index=0
 icg=0
 icg1=0
 mbdkpsp=mband*nkpt_rbz*nsppol
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol

!Update list of computed matrix elements
 if((ipert==natom+3) .or. (ipert==natom+4)) then
!  Eventually expand when strain coupling to other perturbations is
!  implemented
   do ipert1=natom+3,natom+4
     do idir1=1,3
       blkflg(idir1,ipert1,idir,ipert)=1
     end do
   end do
 end if

!allocate(enl1nk(mbdkpsp))
 ABI_ALLOCATE(d2bbb_k,(2,3,mband,mband*prtbbb))
 ABI_ALLOCATE(d2nl_k,(2,3,mpert))


 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))


 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
 nfftot=n1*n2*n3

!Prepare GS k wf file for reading if mkmem==0
 if (mkmem==0) then
   call clsopn(wfftgs)
   call hdr_skip(wfftgs,ierr)

!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wfftgs,mpi_enreg,nband_rbz,npwarr,nspinor,nsppol,nkpt_rbz)

 end if

!Prepare RF wf files for reading and writing if mkmem==0
 if (mk1mem==0) then

   call clsopn(wffnow)

!  Read unwfnow header
   call hdr_skip(wffnow,ierr)

!  Define offsets, in case of MPI I/O
   formeig=1
   call xdefineOff(formeig,wffnow,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)

 end if

 bantot = 0
 ban2tot = 0

!LOOP OVER SPINS
 do isppol=1,nsppol

   if (nsppol/=1) then
     write(message,*)' ****  In nselt3 for isppol=',isppol
     call wrtout(std_out,message,'COLL')
   end if

!  Rewind kpgsph data file if needed:
   if (mkmem==0) rewind(unkg)
   if (mk1mem==0) rewind(unkg1)
   if (mkmem==0.and.psps%useylm==1) rewind(unylm)
   if (mk1mem==0.and.psps%useylm==1) rewind(unylm1)
   ikg=0
   ikg1=0

   ikpt_fbz(1:3)=0

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)

     bantot = bantot + nband_k
     ban2tot = ban2tot + 2*nband_k**2

     if(mpi_enreg%paral_compil_kpt==1)then
!      BEGIN TF_CHANGES
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&       -me))/=0) then
!        END TF_CHANGES
         bdtot_index=bdtot_index+nband_k
         bd2tot_index=bd2tot_index+2*nband_k**2
!        Skip the rest of the k-point loop
         cycle
       end if
     end if

     ABI_ALLOCATE(occ_k,(nband_k))

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang))
     ABI_ALLOCATE(ylm1_k,(npw1_k,mpsang*mpsang))
     if (ipert==natom+3.or.ipert==natom+4) then
       ABI_ALLOCATE(ylmgr_k,(npw_k,3,mpsang*mpsang))
       ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,mpsang*mpsang))
     end if

!    enl1_k(:)=0.0_dp
     d2nl_k(:,:,:)=0.0_dp
     if(prtbbb==1)d2bbb_k(:,:,:,:)=0.0_dp
     kpoint(:)=kpt_rbz(:,ikpt)
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)

     if (mkmem==0) then
!      Read (k+G) basis sphere data (same for each spin)
       nsp=nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,unkg)

!      Read sphere data centered at k in unkg, then k+g data
       read (unkg) ((kg_k(ii,muig),ii=1,3),muig=1,npw_k)

!      Eventually read (k+G) spherical harmonics
       if (psps%useylm==1) then
         read(unylm)
         if (ipert==natom+3.or.ipert==natom+4) then
           read(unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang),&
&           (((ylmgr_k(muig,ii,ilm),muig=1,npw_k),ii=1,3),ilm=1,mpsang*mpsang)
         else
           read(unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
         end if
       end if

     else

       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       if (psps%useylm==1) then
         do ilm=1,mpsang*mpsang
           ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
         end do
         if (ipert==natom+3.or.ipert==natom+4) then
           do ilm=1,mpsang*mpsang
             do ii=1,3
               ylmgr_k(1:npw_k,ii,ilm)=ylmgr(1+ikg:npw_k+ikg,ii,ilm)
             end do
           end do
         end if
       end if

!      End if for choice governed by mkmem
     end if

     wtk_k=wtk_rbz(ikpt)

     kg1_k(:,:) = 0
     if (mk1mem==0) then
!      Read (k+q+G) basis sphere data (same for each spin)
       nsp=nspinor
       call rdnpw(ikpt,isppol,nband_k,npw1_k,nsp,0,unkg1)

!      Read sphere data centered at k in unkg, then k+g data
       read (unkg1) ((kg1_k(ii,muig),ii=1,3),muig=1,npw1_k)

!      Eventually read (k+q+G) spherical harmonics
       if (psps%useylm==1) then
         read(unylm1)
         if (ipert==natom+3.or.ipert==natom+4) then
           read(unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang),&
&           (((ylmgr1_k(muig,ii,ilm),muig=1,npw1_k),ii=1,3),ilm=1,mpsang*mpsang)
         else
           read(unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang)
         end if
       end if

     else

       kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
       if (psps%useylm==1) then
         do ilm=1,mpsang*mpsang
           ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
         end do
         if (ipert==natom+3.or.ipert==natom+4) then
           do ilm=1,mpsang*mpsang
             do ii=1,3
               ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
             end do
           end do
         end if
       end if

!      End if for choice governed by mk1mem
     end if

!    Compute the eigenvalues, wavefunction,
!    contributions to kinetic energy, nonlocal energy, forces,
!    and update of rhor1 to this k-point and this spin polarization.

!    Note that nstwf4 is called with kpoint, while kpt is used inside vtowfk3
     call nstwf4(atindx,atindx1,cg,cg1,d2nl_k,ecut,ecutsm,effmass,&
&     gmet,gprimd,icg,icg1,ikpt,isppol,istwf_k,kg_k,kg1_k,kpoint,&
&     mband,mgfft,mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,natom,&
&     nattyp,nband_k,ngfft,nloalg,npw_k,npw1_k,nspinor,nsppol,ntypat,&
&     occ_k,ph1d,psps,rmet,ucvol,wffnow,wfftgs,wtk_k,xred,ylm_k,ylmgr_k)
     d2nl(:,:,:,idir,ipert)=d2nl(:,:,:,idir,ipert)+d2nl_k(:,:,:)
     if(prtbbb==1)then
       d2bbb(:,:,idir,ipert,:,:) = d2bbb(:,:,idir,ipert,:,:) + &
&       d2bbb_k(:,:,:,:)
     end if

     ABI_DEALLOCATE(occ_k)

!    Keep track of total number of bands (all k points so far, even for
!    k points not treated by me)
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2

!    Shift array memory
     if (mkmem/=0) then
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     end if
     if (mk1mem/=0) then
       icg1=icg1+npw1_k*nspinor*nband_k
       ikg1=ikg1+npw1_k
     end if
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     if (ipert==natom+3.or.ipert==natom+4)  then
       ABI_DEALLOCATE(ylmgr_k)
       ABI_DEALLOCATE(ylmgr1_k)
     end if

!    End big k point loop
   end do

!  End loop over spins
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
!  BEGIN TF_CHANGES
   call leave_test()
!  END TF_CHANGES
   write(message,*) ' nselt3: loop on k-points and spins done in parallel'
   call wrtout(std_out,message,'COLL')
 end if

!Treat now varying occupation numbers
!if(occopt>=3 .and. occopt <=8) then
!SUPPRESSED metallic coding of vtorho

!Treat fixed occupation numbers
!else

!Accumulation over parallel processed now carried out for all terms
!in nstdy3.f

!End of test on varying or fixed occupation numbers
!end if

!The imaginary part of d2nl will be must be set to zero here since
!time-reversal symmetry will always be true for the strain peturbation.
!The symmetry-reduced kpt set will leave a non-zero imaginary part.

 d2nl(2,:,natom+3:natom+4,idir,ipert)=0.0_dp

!Symmetrize the non-local contributions,
!as was needed for the stresses in a ground-state calculation

 if (nsym1>1) then
!  Pack like symmetric-storage cartesian stress tensor
   ii=0
   do ipert1=natom+3,natom+4
     do idir1=1,3
       ii=ii+1
       restr(ii)=d2nl(1,idir1,ipert1,idir,ipert)
     end do
   end do
!  Do the symmetrization using the ground state routine
   call stresssym(gprimd,nsym1,restr,symrc1)
!  Unpack symmetrized stress tensor
   ii=0
   do ipert1=natom+3,natom+4
     do idir1=1,3
       ii=ii+1
       d2nl(1,idir1,ipert1,idir,ipert)=restr(ii)
     end do
   end do
 end if !nsym>1

!----------------------------------------------------------------------------
!Now, treat the local contribution

 ABI_ALLOCATE(vpsp1,(cplex*nfft))
 n3xccc=0
 if(psps%n1xccc/=0)n3xccc=nfft
 ABI_ALLOCATE(xccc3d1,(cplex*n3xccc))
 ABI_ALLOCATE(vxc1,(cplex*nfft,nspden))
 ABI_ALLOCATE(vhartr01,(nfft))

 xccc3d1(:)=0.0_dp

!Double loop over strain perturbations
 do ipert1=natom+3,natom+4
   do idir1=1,3
     if(ipert1==natom+3) then
       istr1=idir1
     else
       istr1=idir1+3
     end if

!    Get first-order local potential.
     call vlocalstr(gmet,gprimd,gsqcut,istr1,mgfft,mpi_enreg,&
&     psps%mqgrid_vl,natom,nattyp,nfft,ngfft,ntypat,paral_kgb,ph1d,psps%qgrid_vl,&
&     ucvol,psps%vlspl,vpsp1)

!    Get first-order hartree potential.
     call hartrestr(gmet,gprimd,gsqcut,idir1,ipert1,mpi_enreg,natom,nfft,ngfft,&
&     paral_kgb,rhog,vhartr01)

!    Get first-order exchange-correlation potential
     if(psps%n1xccc/=0)then
       call mkcor3(cplex,idir1,ipert1,natom,ntypat,n1,psps%n1xccc,&
&       n2,n3,qphon,rprimd,type,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
     end if ! psps%n1xccc/=0

     option=0
     call mkvxcstr3(cplex,idir1,ipert1,kxc,mpi_enreg,natom,nfft,ngfft,&
&     nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,vxc1,xccc3d1)

!    Combines density j2 with local potential j1
     do ispden=1,min(nspden,2)
       do ifft=1,cplex*nfft
         vxc1(ifft,ispden)=vxc1(ifft,ispden)+vpsp1(ifft)+vhartr01(ifft)
       end do
     end do
     call dotprod_vn(cplex,rhor1,dotr,doti,mpi_enreg,nfft,nfftot,nspden,2,vxc1,ucvol)
     write(std_out,*)
     d2lo(1,idir1,ipert1,idir,ipert)=dotr
     d2lo(2,idir1,ipert1,idir,ipert)=doti
   end do ! istr1
 end do ! ipert1

 ABI_DEALLOCATE(vxc1)
 ABI_DEALLOCATE(xccc3d1)
 ABI_DEALLOCATE(vhartr01)

 ABI_DEALLOCATE(d2bbb_k)
 ABI_DEALLOCATE(d2nl_k)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)
 ABI_DEALLOCATE(vpsp1)

!DEBUG
!write(std_out,*)' nselt3: exit '
!ENDDEBUG

end subroutine nselt3
!!***
