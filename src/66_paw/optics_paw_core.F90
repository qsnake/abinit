!{\src2tex{textfont=tt}}
!!****f* ABINIT/optics_paw_core
!! NAME
!! optics_paw_core
!!
!! FUNCTION
!! Compute matrix elements need for X spectr. (in the PAW context) and store them in a file
!!  Matrix elements = <Phi_core|Nabla|Phi_j>
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (SM,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  indlmn(6,lmnmax,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  lmnmax=max. number of (l,m,n) numbers over all types of atom
!!  mband=maximum number of bands
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang =1+maximum angular momentum for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  (only writing in a file)
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      bound_deriv,compmesh,cprj_alloc,cprj_diskinit_r,cprj_free,cprj_get
!!      deducer0,destroy_radmesh,hdr_io,int_ang,leave_new,nderiv_gen,simp_gen
!!      spline,splint,timab,wffclose,wffopen,wrtout,xcomm_init,xcomm_world
!!      xexch_mpi,xme_init,xsum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,eigen0,hdr,indlmn,lmnmax,&
&               mband,mcprj,mkmem,mpi_enreg,mpsang,natom,nkpt,nsppol,pawrad,pawtab)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
 use m_radmesh,   only : ifromr, compmesh, destroy_radmesh, simp_gen, nderiv_gen, deducer0
 use m_splines,   only : spline,splint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'optics_paw_core'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_66_paw, except_this_one => optics_paw_core
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmnmax,mband,mcprj,mkmem,mpsang,natom,nkpt,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom),indlmn(6,lmnmax,dtset%ntypat)
 type(cprj_type) :: cprj(natom,mcprj)
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: basis_size,bdtot_index,cplex,etiq,i1,i2,iatom,ib,ibg
 integer :: ierr,ij_size,ikpt,il,ilm,ilmn,iln
 integer :: iorder_cprj,ios,ispinor,isppol,istwf_k,itypat
 integer :: jb,jbsp,jl,jlm,jlmn,jln,lmn_size,lmncmax
 integer :: me,me_kpt,mesh_size,my_nspinor,nband_k,nmesh,nphicor,npts
 integer :: sender,spaceComm_bandspin,spaceComm_k,spaceComm_w
 integer :: accesswff,fformopt
 logical :: ex,mykpt
 real(dp) :: cpnm1,cpnm2,intg,noccor,r1,r2
 character(len=8) :: dum
 character(len=80) :: fline
 character(len=500) :: message
!arrays
 integer,allocatable :: indlmn_(:,:),indlmn_core(:,:),lcor(:),meshtp(:),meshsz(:),ncor(:)
 real(dp) :: ang_phipphj(mpsang**2,mpsang**2,8)
 real(dp) :: tsec(2)
 real(dp),allocatable :: dphi(:),energy_cor(:),ff(:),int1(:,:),int2(:,:)
 real(dp),allocatable :: logstp(:),phitmp(:),phi_cor(:,:),psinablapsi(:,:,:,:,:)
 real(dp),allocatable :: rad(:),radstp(:),tnm(:,:,:,:,:),work(:)
 type(coeff3_type), allocatable :: phipphj(:)
 type(cprj_type),allocatable :: cprj_k(:,:)
 type(pawrad_type) :: tmpmesh
 type(wffile_type) :: wff1

! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if(dtset%paral_kgb==1.and.mkmem==0)then
   write(message, '(4a)' )ch10,&
&   ' optics_paw_core :  -',ch10,&
&   '  Not compatible with mkmem=0 and band-fft parallelism !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!------------------------------------------------------------------------------------------------
!0-Reading core wavefunctions
!------------------------------------------------------------------------------------------------

 itypat=1

!New format for core WF file
 ex=.false.
 inquire(file='corewf.abinit',iostat=ios,exist=ex)
 if (ios/=0) then
   write(std_out, '(6a,i8,2a)') ch10,&
&   ' inquire : ERROR -',ch10,&
&   '  Checks for existence of file  corewf.abinit',ch10,&
&   '  but INQUIRE statement returns error code',ios,ch10,&
&   '  Action : identify which problem appears with this file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (ex) then
   open(100,file='corewf.abinit',form='formatted')
   read(100,*) ! skip title
   read(100,*) ! skip method,nspinor,nsppol
   read(100,*) ! skip zatom,zcore,pspdat
   read(100,*) ! skip pspcod,pspxc,lmax
   read(100,*) ! skip pspfmt,creatorID
   read(100,*) nphicor
   read(100,*) ! skip orbitals
   read(100,*) nmesh
   ABI_ALLOCATE(meshsz,(nmesh))
   ABI_ALLOCATE(meshtp,(nmesh))
   ABI_ALLOCATE(radstp,(nmesh))
   ABI_ALLOCATE(logstp,(nmesh))
   do iln=1,nmesh
     r2=zero;read(100,'(a80)') fline
     read(unit=fline,fmt=*,err=20,end=20) ib,i1,i2,r1,r2
     20 continue
     if (ib<=nmesh) then
       meshtp(ib)=i1;meshsz(ib)=i2
       radstp(ib)=r1;logstp(ib)=r2
     end if
   end do
   read(100,*) ! skip rmax(core)
   ABI_ALLOCATE(ncor,(nphicor))
   ABI_ALLOCATE(lcor,(nphicor))
   ABI_ALLOCATE(energy_cor,(nphicor))
   ABI_ALLOCATE(phi_cor,(meshsz(1),nphicor))
   do iln=1,nphicor
     read(100,*) ! skip comment
     read(100,*) i1
     read(100,*) ncor(iln),lcor(iln)
     read(100,*) energy_cor(iln)
     ABI_ALLOCATE(phitmp,(meshsz(i1)))
     read(100,*) phitmp
     if ((pawrad(itypat)%mesh_type/=meshtp(i1)) &
&     .or.(pawrad(itypat)%rstep/=radstp(i1)) &
&     .or.(pawrad(itypat)%lstep/=logstp(i1))) then
       tmpmesh%mesh_type=meshtp(i1);tmpmesh%mesh_size=meshsz(i1)
       tmpmesh%rstep=radstp(i1);tmpmesh%lstep=logstp(i1)
       ABI_ALLOCATE(tmpmesh%rad,(meshsz(i1)))
       ABI_ALLOCATE(tmpmesh%radfact,(meshsz(i1)))
       ABI_ALLOCATE(tmpmesh%simfact,(meshsz(i1)))
       call compmesh(tmpmesh,-1._dp)
       npts=pawrad(itypat)%mesh_size
       if (tmpmesh%rmax<pawrad(itypat)%rmax+tol8) npts=ifromr(pawrad(itypat),tmpmesh%rmax)-1
       ABI_ALLOCATE(work,(meshsz(i1)))
       call bound_deriv(phitmp,tmpmesh,meshsz(i1),r1,r2)
       call spline(tmpmesh%rad,phitmp,meshsz(i1),r1,r2,work)
       call splint(meshsz(i1),tmpmesh%rad,phitmp,work,npts,pawrad(itypat)%rad(1:npts),phi_cor(1:npts,iln))
       if (npts<pawrad(itypat)%mesh_size) phi_cor(npts+1:pawrad(itypat)%mesh_size,iln)=zero
       ABI_DEALLOCATE(work)
       call destroy_radmesh(tmpmesh)
     else
       phi_cor(1:meshsz(i1),iln)=phitmp(1:meshsz(i1))
       if (meshsz(i1)<pawrad(itypat)%mesh_size) phi_cor(meshsz(i1)+1:pawrad(itypat)%mesh_size,iln)=zero
     end if
     ABI_DEALLOCATE(phitmp)
   end do
   ABI_DEALLOCATE(meshsz)
   ABI_DEALLOCATE(meshtp)
   ABI_DEALLOCATE(radstp)
   ABI_DEALLOCATE(logstp)
 else

!  Old format for core WF file
   ex=.false.
   inquire(file='corewf.dat',iostat=ios,exist=ex)
   if (ios/=0) then
     write(std_out, '(6a,i8,2a)') ch10,&
&     ' inquire : ERROR -',ch10,&
&     '  Checks for existence of file  corewf.dat',ch10,&
&     '  but INQUIRE statement returns error code',ios,ch10,&
&     '  Action : identify which problem appears with this file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if (ex) then
     open(100,file='corewf.dat',form='formatted')
     do while (dum/='atompaw ')
       read(100,'(a8)') dum
     end do
     read(100,'(2i4)') npts,nphicor
     ABI_ALLOCATE(ncor,(nphicor))
     ABI_ALLOCATE(lcor,(nphicor))
     ABI_ALLOCATE(energy_cor,(nphicor))
     ABI_ALLOCATE(phi_cor,(npts,nphicor))
     ABI_ALLOCATE(rad,(npts))
     do iln=1,nphicor
       read(100,'("# n=",i4," l=",i4," nocc=",f15.7," energy=",f15.7)') &
&       ncor(iln),lcor(iln),noccor,energy_cor(iln)
       do jln=1,npts
         read(100,*) rad(jln),phi_cor(jln,iln)
       end do
       read(100,*)
     end do
     ABI_DEALLOCATE(rad)
     close(100)
   else

!    No core WF file found !
     write(message, '(6a)' ) ch10,&
&     ' inquire : ERROR -',ch10,&
&     '  Checks for existence of files corewf.abinit or corewf.dat (old format)',ch10,&
&     '  but INQUIRE finds file does not exist.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')

   end if
 end if

!set an array 'a la' indlmn
 lmncmax=0
 do ib=1,nphicor
   il=lcor(ib)
   lmncmax=lmncmax+2*il+1
 end do
 ABI_ALLOCATE(indlmn_core,(6,lmncmax))
 ilmn=0;iln=0
 do ib=1,nphicor
   il=lcor(ib)
   iln=iln+1
   do ilm=1,2*il+1
     indlmn_core(1,ilmn+ilm)=il
     indlmn_core(2,ilmn+ilm)=ilm-(il+1)
     indlmn_core(3,ilmn+ilm)=1
     indlmn_core(4,ilmn+ilm)=il*il+ilm
     indlmn_core(5,ilmn+ilm)=iln
     indlmn_core(6,ilmn+ilm)=1
   end do
   ilmn=ilmn+2*il+1
 end do

!----------------------------------------------------------------------------------
!1-Computation of phipphj=<phi_i|nabla|phi_core>
!----------------------------------------------------------------------------------

!1-A Integration of the angular part : all angular integrals have been
!computed outside Abinit and tabulated for each (l,m) value
!----------------------------------------------------------------------------------

 call int_ang(ang_phipphj,mpsang)

 ABI_ALLOCATE(phipphj,(dtset%ntypat))

!We consider the impurity to be the first atom
!loop on atoms type
 do itypat=1,dtset%ntypat

   mesh_size=pawrad(itypat)%mesh_size
   lmn_size=pawtab(itypat)%lmn_size
   basis_size=pawtab(itypat)%basis_size
   ij_size=lmn_size*lmn_size

   ABI_ALLOCATE(indlmn_,(6,lmnmax))
   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(rad,(mesh_size))
   ABI_ALLOCATE(int2,(lmn_size,lmncmax))
   ABI_ALLOCATE(int1,(lmn_size,lmncmax))
   ABI_ALLOCATE(dphi,(mesh_size))
   ABI_ALLOCATE(phipphj(itypat)%value,(3,lmn_size,lmncmax))

   indlmn_(:,:)=indlmn(:,:,itypat)
   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)

!  1-B Computation of int1=\int phi phi_core /r dr
!  ----------------------------------------------------------------------------------
   do jln=1,nphicor
     do iln=1,basis_size
       ff(2:mesh_size)=(pawtab(itypat)%phi(2:mesh_size,iln)*phi_cor(2:mesh_size,jln))/rad(2:mesh_size)
       call deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int1(iln,jln)=intg
     end do
   end do

!  1-C Computation of int2=\int phi/r d/dr(phi_core/r) - phi phi_core/r dr
!  ----------------------------------------------------------------------------------
   do jln=1,nphicor
     ff(1:mesh_size)=phi_cor(1:mesh_size,jln)
     call nderiv_gen(dphi,ff,1,pawrad(itypat))

     do iln=1,basis_size
       ff(2:mesh_size)=pawtab(itypat)%phi(2:mesh_size,iln)*dphi(2:mesh_size) &
&       -pawtab(itypat)%phi(2:mesh_size,iln)*phi_cor(2:mesh_size,jln)/ &
&       rad(2:mesh_size)
       call deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int2(iln,jln)=intg
     end do
   end do

!  1-D Integration of the radial part
!  ----------------------------------------------------------------------------------
   do jlmn=1,lmncmax
     jlm=indlmn_core(4,jlmn)
     jl=indlmn_core(5,jlmn)
     do ilmn=1,lmn_size
       ilm=indlmn_(4,ilmn)
       il=indlmn_(5,ilmn)
       phipphj(itypat)%value(1,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,1)&
&       + int1(il,jl)*(ang_phipphj(ilm,jlm,2)+ang_phipphj(ilm,jlm,3))
       phipphj(itypat)%value(2,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,4)&
&       + int1(il,jl)*(ang_phipphj(ilm,jlm,5)+ang_phipphj(ilm,jlm,6))
       phipphj(itypat)%value(3,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,7)&
&       + int1(il,jl)*ang_phipphj(ilm,jlm,8)
     end do
   end do

   ABI_DEALLOCATE(indlmn_)
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(rad)
   ABI_DEALLOCATE(int2)
   ABI_DEALLOCATE(int1)
   ABI_DEALLOCATE(dphi)

!  end loop on atoms type
 end do

!----------------------------------------------------------------------------------
!2- Computation of <psi_n|p_i>(<phi_i|-i.nabla|phi_core>)
!----------------------------------------------------------------------------------

!Init parallelism
 call xcomm_world(mpi_enreg,spaceComm_w,myrank=me)
 call xcomm_init(mpi_enreg,spaceComm_k,spaceComm_bandfft=mpi_enreg%comm_kpt)
 call xme_init(mpi_enreg,me_kpt)
 if (mpi_enreg%mode_para=='b') then
   spaceComm_bandspin=mpi_enreg%comm_bandspin
 else
   spaceComm_bandspin=0
 end if
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)

!Prepare temporary PAW file if mkmem==0
 iorder_cprj=0
 call cprj_diskinit_r(atindx1,natom,iorder_cprj,mkmem,natom,0,dimcprj,my_nspinor,dtfil%unpaw)


!open _OPT2 file for proc 0
 accesswff=-1
 fformopt=611
 call WffOpen(accesswff,spaceComm_k,dtfil%fnameabo_app_opt2,ierr,wff1,0,me,124)
 call hdr_io(fformopt,hdr,2,wff1)
 if (me==0) then
   write(124)(eigen0(ib),ib=1,mband*nkpt*nsppol)
   write(124) nphicor
   do iln=1,nphicor
     write(124) ncor(iln),lcor(iln),half*energy_cor(iln)
   end do
 end if
 ABI_ALLOCATE(psinablapsi,(2,3,mband,nphicor,natom))

!LOOP OVER SPINS
 ibg=0
 bdtot_index=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     etiq=ikpt+(isppol-1)*nkpt
     psinablapsi=zero

     mykpt=.true.
     if (mpi_enreg%paral_compil_kpt==1) then
       mykpt=(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me_kpt))==0)
     end if
     if (mykpt) then

!      Allocations depending on k-point
       istwf_k=dtset%istwfk(ikpt)
       cplex=2;if (istwf_k>1) cplex=1
       ABI_ALLOCATE(cprj_k,(natom,my_nspinor*nband_k))
       call cprj_alloc(cprj_k,0,dimcprj)

!      Extract cprj for this k-point according to mkmem
       call cprj_get(atindx1,cprj_k,cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&       mband,mkmem,mpi_enreg,natom,nband_k,nband_k,my_nspinor,nsppol,dtfil%unpaw)

       ABI_ALLOCATE(tnm,(2,3,nband_k,nphicor,natom))
       tnm=zero

!      Loops on bands
       do jb=1,nband_k

         if (mpi_enreg%mode_para=='b') then
           if (mod(jb-1,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
         elseif (mpi_enreg%paral_compil_kpt==1) then
           if (abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-me_kpt)/=0) cycle
         end if
         jbsp=(jb-1)*my_nspinor

         if (cplex==1) then
           do ispinor=1,my_nspinor
             jbsp=jbsp+1
             do iatom=1,natom
               itypat=dtset%typat(iatom)
               lmn_size=pawtab(itypat)%lmn_size
               do jlmn=1,lmn_size
                 do ilmn=1,lmncmax
                   ib=indlmn_core(5,ilmn)
                   cpnm1=cprj_k(iatom,jbsp)%cp(1,jlmn)
                   tnm(2,:,jb,ib,iatom)=tnm(2,:,jb,ib,iatom)+cpnm1*phipphj(itypat)%value(:,jlmn,ilmn)
                 end do !ilmn
               end do !jlmn
             end do !iatom
           end do !ispinor
         else
           do ispinor=1,my_nspinor
             jbsp=jbsp+1
             do iatom=1,natom
               itypat=dtset%typat(iatom)
               lmn_size=pawtab(itypat)%lmn_size
               do jlmn=1,lmn_size
                 do ilmn=1,lmncmax
                   ib=indlmn_core(5,ilmn)
                   cpnm1=cprj_k(iatom,jbsp)%cp(1,jlmn)
                   cpnm2=cprj_k(iatom,jbsp)%cp(2,jlmn)
                   tnm(1,:,jb,ib,iatom)=tnm(1,:,jb,ib,iatom)+cpnm1*phipphj(itypat)%value(:,jlmn,ilmn)
                   tnm(2,:,jb,ib,iatom)=tnm(2,:,jb,ib,iatom)+cpnm2*phipphj(itypat)%value(:,jlmn,ilmn)
                 end do !ilmn
               end do !jlmn
             end do !iatom
           end do !ispinor
         end if

!        End loops on bands
       end do ! jb

!      Reduction in case of parallelism
       if (mpi_enreg%mode_para=='b') then
         call timab(48,1,tsec)
         call xsum_master(tnm,0,spaceComm_bandspin,ierr)
         call timab(48,2,tsec)
       end if

       psinablapsi(:,:,:,:,:)=psinablapsi(:,:,:,:,:)+tnm(:,:,:,:,:)
       ABI_DEALLOCATE(tnm)

       if (mkmem/=0) ibg = ibg + my_nspinor*nband_k

       call cprj_free(cprj_k)
       ABI_DEALLOCATE(cprj_k)

       if (me==0) then
         do iatom=1,natom
           write(124) ((psinablapsi(1:2,1,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
           write(124) ((psinablapsi(1:2,2,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
           write(124) ((psinablapsi(1:2,3,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
         end do
       elseif (mpi_enreg%me_band==0.and.mpi_enreg%me_fft==0) then
         call xexch_mpi(psinablapsi,etiq,me_kpt,psinablapsi,0,spaceComm_k,ierr)
       end if

     elseif (me==0) then
       sender=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
       call xexch_mpi(psinablapsi,etiq,sender,psinablapsi,0,spaceComm_k,ierr)
       do iatom=1,natom
         write(124) ((psinablapsi(1:2,1,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
         write(124) ((psinablapsi(1:2,2,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
         write(124) ((psinablapsi(1:2,3,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
       end do
     end if ! mykpt

     bdtot_index=bdtot_index+nband_k

!    End loop on spin,kpt
   end do ! ikpt
 end do !isppol

!Close file
 call WffClose(wff1,ierr)

!Datastructures deallocations
 do itypat=1,dtset%ntypat
   ABI_DEALLOCATE(phipphj(itypat)%value)
 end do
 ABI_DEALLOCATE(phipphj)
 ABI_DEALLOCATE(ncor)
 ABI_DEALLOCATE(lcor)
 ABI_DEALLOCATE(energy_cor)
 ABI_DEALLOCATE(phi_cor)
 ABI_DEALLOCATE(psinablapsi)

 DBG_EXIT("COLL")

 end subroutine optics_paw_core
!!***
