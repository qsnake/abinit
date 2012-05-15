!{\src2tex{textfont=tt}}
!!****f* ABINIT/berryphase_new
!! NAME
!! berryphase_new
!!
!! FUNCTION
!! This routine computes the Berry Phase polarization
!!  and the finite difference expression of the ddk.
!!  [see for example Na Sai et al., PRB 66, 104108 (2002)]
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT  group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! cprj(natom,mcprj*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!! dtfil <type(datafiles_type)>=variables related to files
!! dtset <type(dataset_type)>=all input variables in this dataset
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!! indlmn(6,lmnmax,ntypat)
!!   array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
!!                                  or i=lmn (if useylm=1)
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! lmnmax  If useylm=0, max number of (l,m,n) comp. over all type of psps (lnproj)
!!         If useylm=1, max number of (l,n)   comp. over all type of psps (lmnproj)
!! mband=maximum number of bands
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mkmem=number of k points which can fit in memory; set to 0 if use disk
!! mpi_enreg=informations about MPI parallelization
!! mpw=maximum dimensioned size of npw
!! natom=number of atoms in cell
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! ntypat=number of types of atoms in unit cell
!! nkpt=number of k-points
!! option = 1: compute Berryphase polarization
!!          2: compute finite difference expression of the ddk
!!          3: compute polarization & ddk
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)> atomic occupancies
!! pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3.
!! unit_out= unit for output of the results (usually the .out file of ABINIT)
!!   The option unit_out = 0 is allowed. In this case, no information is written
!!   to the output file but only to the log file.
!! usecprj=1 if cprj datastructure has been allocated
!! usepaw= 1: use paw framework. 0:do not use paw.
!! wffnow=struct info for wf disk file
!! xred(3,natom)=reduced atomic coordinates
!! zion(ntypat)=valence charge of each type of atom
!!
!! OUTPUT
!! pel(3) = reduced coordinates of the electronic polarization (a. u.)
!! pelev(3)= expectation value polarization term (PAW only) in cartesian coordinates
!! pion(3)= reduced coordinates of the ionic polarization (a. u.)
!!
!! SIDE EFFECTS
!! Input/Output
!! dtefield <type(efield_type)> = variables related to Berry phase
!!       and electric field calculations (see initberry.f).
!!       In case berryopt = 4, the overlap matrices computed
!!       in this routine are stored in dtefield%smat in order
!!       to be used in the electric field calculation.
!!
!! TODO
!!  - Use the analytical relation between the overlap matrices
!!    S(k,k+dk) and S(k+dk,k) to avoid to recompute them
!!    when ifor = 2.
!!
!! NOTES
!! - pel and pion do not take into account the factor 1/ucvol
!! - In case of a ddk calculation, the eigenvalues are not computed.
!! - The ddk computed by this routine should not be used to
!!   compute the electronic dielectric tensor.
!!
!! PARENTS
!!      elpolariz,scfcv
!!
!! CHILDREN
!!      appdig,cprj_alloc,cprj_copy,cprj_free,cprj_get,cprj_mpi_allgather
!!      cprj_mpi_recv,cprj_mpi_send,cprj_put,leave_new,outwf,pawpolev,polcart
!!      rhophi,smatrix,smatrix_k_paw,sym_cprj_kn,wrtout,xallgather_mpi
!!      xcomm_init,xme_init,xrecv_mpi,xredxcart,xsend_mpi,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,&
&  gprimd,hdr,indlmn,kg,lmnmax,mband,mcg,mcprj,&
&  mkmem,mpi_enreg,mpw,natom,npwarr,nsppol,ntypat,&
&  nkpt,option,pawrhoij,pawtab,pel,pelev,pion,pwind,&
&  pwind_alloc,pwnsfac,&
&  rprimd,typat,ucvol,unit_out,usecprj,usepaw,wffnow,xred,zion)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_wffile
 use m_efield

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'berryphase_new'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => berryphase_new
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in) :: lmnmax,mband,mcg,mcprj,mkmem,mpw,natom,nkpt
 integer, intent(in) :: nsppol,ntypat,option
 integer, intent(in) :: pwind_alloc,unit_out,usecprj,usepaw
 real(dp), intent(in) :: ucvol
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(efield_type), intent(inout) :: dtefield
 type(hdr_type), intent(inout) :: hdr
 type(wffile_type), intent(inout) :: wffnow
!arrays
 integer, intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kg(3,mpw*mkmem)
 integer, intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: cg(2,mcg),gprimd(3,3)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(in) :: rprimd(3,3),zion(ntypat)
 real(dp), intent(inout) :: xred(3,natom)
 real(dp), intent(out) :: pel(3),pelev(3),pion(3)
 type(pawrhoij_type), intent(in) :: pawrhoij(natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
 type(cprj_type),intent(in) ::  cprj(natom,mcprj*usecprj)

!Local variables -------------------------
 integer :: count,count1,ddkflag,dest
 integer :: iatom,icg,icg1,idir,idum,ikpt1i_sp
 integer :: ierr,ifor,ikg,ikpt,ikpt1,ikpt_loc!,ikpt2,ikpt2i,npw_k2, itrs
 integer :: icp1, icp2,iproc
! integer :: ii ! appears commented out below in a debug section
 integer :: inibz,ikpt1i
 integer :: isppol,istr,itypat,jkpt,jkstr,job,jsppol
 integer :: maxbd,mcg1_k
 integer :: minbd,mxfh,my_nspinor,nband_k,ncpgr,nfor,npw_k1,nstep,ntotcp,nxfh,n2dim,nproc,pertcase
 integer :: response,shiftbd,source,spaceComm,tag
 integer :: jj,jstr,kk,ineigh_str
 integer :: istep,jstep,kpt_mark(dtefield%fnkpt),nkstr,nstr,iunmark,berrystep
 integer :: jkpt2, jkpt2i
 real(dp) :: det_mod,dkinv,dphase,dtm_real,dtm_imag,fac,gmod,phase0
 real(dp) :: pol,polbtot,polion,politot,poltot,rho
 real(dp) :: dphase_new,dphase_init
 character(len=fnlen) :: fiwf1o
 character(len=500) :: message
 type(wvl_wf_type) :: wfs
 type(wvl_internal_type) :: wvl
 integer,allocatable :: dimlmn(:),ikpt1_recv(:), sflag_k(:)!,pwind_k(:)
 integer,allocatable :: ikpt3(:), ikpt3i(:), sflag_k_mult(:,:),npw_k3(:)
 integer,allocatable :: idxkstr_mult(:,:), pwind_k_mult(:,:),itrs_mult(:)
 real(dp) :: det_average(2),dk(3),dtm_k(2),gpard(3),pel_cart(3),pion_cart(3)
 real(dp) :: polb(nsppol),ptot_cart(3),rel_string(2),xcart(3,natom)
 real(dp) :: delta_str(2),dist_,dstr(2)
 real(dp),allocatable :: buffer(:,:),buffer1(:),buffer2(:)
 real(dp),allocatable :: cg1(:,:),cg1_k(:,:),cgq(:,:)
 real(dp),allocatable :: det_string(:,:),eig_dum(:),occ_dum(:)!,dtm(:,:)
 real(dp),allocatable :: polberry(:),resid(:),pwnsfac_k(:,:)
 real(dp),allocatable :: smat_inv(:,:,:),smat_k(:,:,:),smat_k_paw(:,:,:)
 real(dp),allocatable :: xfhist(:,:,:,:)
 real(dp),allocatable :: str_flag(:)
! real(dp),allocatable :: dist_str(:,:),det_string_test(:,:)
 real(dp),allocatable :: dtm_mult(:,:,:), coef(:,:), polb_mult(:,:)
 type(cprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_buf(:,:),cprj_gat(:,:)
 type(cprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

!BEGIN TF_CHANGES
 integer :: me
!END TF_CHANGES

!no_abirules

! ***********************************************************************

!DEBUG
!write(std_out,*)' berryphase_new : enter '
!do ii=1,pwind_alloc
!write(std_out,*)ii,pwnsfac(:,ii)
!end do
!stop
!ENDDEBUG

!Init MPI
 call xcomm_init(mpi_enreg,spaceComm)
 call xme_init(mpi_enreg,me)
 nproc=xcomm_size(spaceComm)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)

!allocate(pwind_k(mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,mpw))
 ABI_ALLOCATE(sflag_k,(dtefield%nband_occ))
!pwind_k(:) = 0
 pwnsfac_k(1,:) = 1.0_dp ! bra real
 pwnsfac_k(2,:) = 0.0_dp ! bra imag
 pwnsfac_k(3,:) = 1.0_dp ! ket real
 pwnsfac_k(4,:) = 0.0_dp ! ket imag

!if (dtset%nspinor == 2) then
!write(message, '(a,a,a,a)' )ch10,&
!&   ' berryphase : ERROR -',ch10,&
!&   '  This routine does not yet work for nspinor = 2'
!call wrtout(std_out,message,'COLL')
!call leave_new('COLL')
!end if

 if (maxval(dtset%istwfk(:)) /= 1) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' berryphase : BUG -',ch10,&
&   '  This routine does not work yet with istwfk /= 1.',ch10,&
&   '  This should have been tested previously ...'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if (usepaw == 1 .and. usecprj /= 1) then
   write(message, '(4a)' )ch10,&
&   ' berryphase : BUG -',ch10,&
&   '  PAW calculation but cprj datastructure has not been allocated !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 mcg1_k = mpw*mband
 shiftbd = 1
 if (option > 1) then
   ABI_ALLOCATE(cg1,(2,mcg))
   ABI_ALLOCATE(eig_dum,(2*mband*mband*nkpt*nsppol))
   ABI_ALLOCATE(occ_dum,(mband*nkpt*nsppol))
   eig_dum(:) = zero
   occ_dum(:) = dtefield%sdeg
 end if

!initialize variable tied to multiple step computation
 berrystep=dtset%berrystep
 ABI_ALLOCATE(ikpt3,(berrystep))
 ABI_ALLOCATE(ikpt3i,(berrystep))
 ABI_ALLOCATE(sflag_k_mult,(dtefield%nband_occ,berrystep))
 ABI_ALLOCATE(npw_k3,(berrystep))
 ABI_ALLOCATE(pwind_k_mult,(mpw,berrystep))
 ABI_ALLOCATE(itrs_mult,(berrystep))
 ABI_ALLOCATE(coef,(berrystep,berrystep))
 ABI_ALLOCATE(polb_mult,(nsppol,berrystep))
!coefficient for berryphase computation
 coef(:,:) = 0.0_dp
 do jstep = 1, berrystep
   coef(jstep,1) = 1.d0/real(jstep*jstep,dp)
   if(jstep/=1)coef(jstep,1)=coef(jstep,1)/real(1-jstep*jstep,dp)
 end do
 do istep = 2, berrystep
   do jstep = 1, berrystep
     coef(jstep, istep) = real(istep*istep,dp)*coef(jstep,istep-1)
     if(jstep /= istep)coef(jstep, istep)=coef(jstep,istep)/real(istep*istep-jstep*jstep,dp)
   end do
 end do
!the berryphase using the strings of steps dk, 2*dk, ..., istep*dk is :
!coef(1,istep)*berryphase(dk) + coef(2,istep)*berryphase(2*dk) + ... + coef(istep,istep)*berryphase(istep*dk)
!DEBUG
!write(std_out,*)'coef, sum coef'
!do istep=1,step
!write(std_out,*)coef(:,istep), sum(coef(1:istep,istep))
!end do
!ENDDEBUG

!allocate(dtm(2,dtefield%fnkpt*nsppol))
 ABI_ALLOCATE(dtm_mult,(2,dtefield%fnkpt*nsppol,berrystep))
 ABI_ALLOCATE(cg1_k,(2,mcg1_k))

 if (usepaw == 1) then ! cprj allocation
   if (dtset%berryopt == 4) then !finite electric field needs gradients for forces
     ncpgr = 3
   else
     ncpgr = 0
   end if
   ABI_ALLOCATE(dimlmn,(natom))
   do iatom=1,natom
     dimlmn(iatom)=pawtab(typat(iatom))%lmn_size
   end do
   ABI_ALLOCATE(cprj_k,(natom,dtefield%nband_occ))
   ABI_ALLOCATE(cprj_kb,(natom,dtefield%nband_occ))
   ABI_ALLOCATE(cprj_gat,(natom,nproc*dtefield%nband_occ))
   call cprj_alloc(cprj_k,ncpgr,dimlmn)
   call cprj_alloc(cprj_kb,ncpgr,dimlmn)
   call cprj_alloc(cprj_gat,ncpgr,dimlmn)
   if (dtset%kptopt /= 3) then
     ABI_ALLOCATE(cprj_ikn,(natom,dtefield%nband_occ))
     ABI_ALLOCATE(cprj_fkn,(natom,dtefield%nband_occ))
     call cprj_alloc(cprj_ikn,ncpgr,dimlmn)
     call cprj_alloc(cprj_fkn,ncpgr,dimlmn)
   end if

   if ( mpi_enreg%paral_compil_kpt==1 ) then
     n2dim = dtefield%nband_occ
     ntotcp = n2dim*SUM(dimlmn(:))
     ABI_ALLOCATE(cprj_buf,(natom,dtefield%nband_occ))
     call cprj_alloc(cprj_buf,ncpgr,dimlmn)
   end if

   if (dtset%berryopt == 4 .or. abs(dtset%berryopt) == 5 ) then
     write(message,'(2a,i5,2a)')ch10,&
&     ' nkpt = ',nkpt,ch10,' copy cprj to dtefield%cprj '
     call wrtout(std_out,message,'COLL')

     do isppol = 1, nsppol

       ikpt_loc = 0
       ikpt1 = 0
       do while (ikpt_loc < mkmem)

         if (ikpt_loc < mkmem) ikpt1 = ikpt1 + 1
         if ((ikpt1 > nkpt).and.(ikpt_loc < mkmem)) exit
         nband_k = dtset%nband(ikpt1) 

         if ( mpi_enreg%paral_compil_kpt==1 ) then
           if ( (minval(abs(mpi_enreg%proc_distrb(ikpt1,1:nband_k,isppol)-me))/=0) .and. &
&           (ikpt_loc <= mkmem) ) cycle
         end if

         ikpt_loc = ikpt_loc + 1

         ABI_ALLOCATE(ikpt1_recv,(nproc))
         call xallgather_mpi(ikpt1,ikpt1_recv,spaceComm,ierr)
         call cprj_get(atindx1,cprj_k,cprj,natom,1,(ikpt_loc-1)*nband_k,ikpt1,0,isppol,mband,&
&         mkmem,mpi_enreg,natom,nband_k,nband_k,my_nspinor,nsppol,0)
         call cprj_mpi_allgather(cprj_k,cprj_gat,natom,nband_k,dimlmn,ncpgr,nproc,spaceComm,ierr)
         do iproc = 1, nproc
           icp2=nband_k*(iproc-1)
           call cprj_get(atindx1,cprj_k,cprj_gat,natom,1,icp2,ikpt1,0,1,nband_k,&
&           nproc,mpi_enreg,natom,nband_k,nband_k,my_nspinor,1,0)
           icp1 = nband_k*(ikpt1_recv(iproc)-1)
           call cprj_put(atindx1,cprj_k,dtefield%cprj,natom,1,icp1,ikpt1,0,isppol,&
&           nband_k,dtefield%fnkpt,mpi_enreg,natom,nband_k,nband_k,dimlmn,my_nspinor,nsppol,spaceComm,0)
         end do
         ABI_DEALLOCATE(ikpt1_recv)

       end do ! close loop over k-points

     end do ! end loop over nsppol

   end if ! end check on berryopt == 4
 end if

 pel(:) = zero     ; pion(:) = zero

 minbd = 1   ;  maxbd = dtefield%nband_occ

 do idir = 1, 3

!  dtm(:,:) = zero
   dtm_mult(:,:,:) = zero

   if (dtset%rfdir(idir) == 1) then

     if (abs(dtefield%efield_dot(idir)) < tol12) dtefield%sflag(:,:,:,idir) = 0

!    Check whether the polarization or the ddk must be computed

!    nfor = 1 : to compute P, I only need the WF at k + dk
!    nfor = 2 : to compute the ddk, I need the WF at k + dk and k - dk
!    dkinv    : +-1/2dk

     if (option > 1) then

       ddkflag = 1
       nfor = 2
       job = 1
       cg1(:,:) = zero
       if (option == 3) job = 11

     else if (option == 1) then

       ddkflag = 0
       nfor = 1
       job = 10

     end if

     dk(:) = dtefield%dkvecs(:,idir)
     gpard(:) = dk(1)*gprimd(:,1) + dk(2)*gprimd(:,2) + dk(3)*gprimd(:,3)
     gmod = sqrt(dot_product(gpard,gpard))
     if (option > 1) dkinv = one/(two*dk(idir))

     write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&     ' Computing the polarization (Berry phase) for reciprocal vector:',ch10,&
&     dk(:),' (in reduced coordinates)',ch10,&
&     gpard(1:3),' (in cartesian coordinates - atomic units)'
     call wrtout(std_out,message,'COLL')
     if (unit_out /= 0) then
       call wrtout(unit_out,message,'COLL')
     end if

     write(message,'(a,i5,a,a,i5)')&
&     ' Number of strings: ',dtefield%nstr(idir),ch10,&
&     ' Number of k points in string:', dtefield%nkstr(idir)
     call wrtout(std_out,message,'COLL')
     if (unit_out /= 0) then
       call wrtout(unit_out,message,'COLL')
     end if

     if ((option == 2).or.(option == 3)) then

       write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&       ' Computing the ddk (Berry phase) for reciprocal vector:',ch10,&
&       dk(:),' (in reduced coordinates)',ch10,&
&       gpard(1:3),' (in cartesian coordinates - atomic units)'
       call wrtout(std_out,message,'COLL')
       if (unit_out /= 0) then
         call wrtout(unit_out,message,'COLL')
       end if

     end if

     do ifor = 1, nfor

       if (ifor == 2) then
         dk(:) = -1_dp*dk(:)
         job = 1   ! only the inverse of the overlap matrix is required
         dkinv = -1_dp*dkinv
       end if


!      Compute the determinant and/or the inverse of the overlap matrix
!      for each pair of k-points < u_nk | u_nk+dk >

       icg = 0 ; icg1 = 0
       ABI_ALLOCATE(smat_k,(2,dtefield%nband_occ,dtefield%nband_occ))
       ABI_ALLOCATE(smat_inv,(2,dtefield%nband_occ,dtefield%nband_occ))
       ABI_ALLOCATE(smat_k_paw,(2,usepaw*dtefield%nband_occ,usepaw*dtefield%nband_occ))

!      Loop on the values of ikpt_loc and ikpt1 :
!      ikpt1 is incremented one by one, and number the k points in the FBZ
!      ikpt1i refer to the k point numbering in the IBZ
!      ikpt_loc differs from ikpt1 only in the parallel case, and gives
!      the index of the k point in the FBZ, in the set treated by the present processor
!      NOTE : in order to allow synchronisation, ikpt_loc contain information about
!      ikpt AND ISPPOL !
!      It means that the following loop is equivalent to a double loop :
!      do isppol = 1, nsppol
!      do ikpt1 =  1, dtefield%fmkmem
!      
       do ikpt_loc = 1, dtefield%fmkmem_max*nsppol

         ikpt1=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,1)
         isppol=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,2)

         if (ikpt1 > 0 .and. isppol > 0) then

           ikpt1i = dtefield%indkk_f2ibz(ikpt1,1)
           nband_k = dtset%nband(ikpt1i + (isppol-1)*dtset%nkpt)

!          DEBUG
!          Please keep this debugging feature
!          write(std_out,'(a,5i4)' )' berryphase_new : ikpt_loc,ikpt1,isppol,idir,ifor=',&
!          &                                  ikpt_loc,ikpt1,isppol,idir,ifor
!          ENDDEBUG

           inibz=0
           if (dtset%kptns(1,ikpt1i) == dtefield%fkptns(1,ikpt1) .and. &
&           dtset%kptns(2,ikpt1i) == dtefield%fkptns(2,ikpt1) .and. &
&           dtset%kptns(3,ikpt1i) == dtefield%fkptns(3,ikpt1)) inibz=1

           ikg = dtefield%fkgindex(ikpt1)
!          ikpt2 = dtefield%ikpt_dk(ikpt1,ifor,idir)
!          ikpt2i = dtefield%indkk_f2ibz(ikpt2,1)

!          ikpt3(istep) : index of kpt1 + istep*dk in the FBZ
!          ikpt3i(istep) : index of kpt1 + istep*dk in the IBZ
           ikpt3(1) = dtefield%ikpt_dk(ikpt1,ifor,idir)
           ikpt3i(1) = dtefield%indkk_f2ibz(ikpt3(1),1)
           do istep = 1, berrystep-1
             ikpt3(istep+1) = dtefield%ikpt_dk(ikpt3(istep),ifor,idir)
             ikpt3i(istep+1) = dtefield%indkk_f2ibz(ikpt3(istep+1),1)
           end do

!          itrs = 0
!          if (dtefield%indkk_f2ibz(ikpt1,6) == 1 ) itrs = itrs + 1
!          if (dtefield%indkk_f2ibz(ikpt2,6) == 1 ) itrs = itrs + 10

           itrs_mult(:)=0
           if (dtefield%indkk_f2ibz(ikpt1,6) == 1 ) itrs_mult(:) = itrs_mult(:) + 1
           do istep=1,berrystep
             if (dtefield%indkk_f2ibz(ikpt3(istep),6) == 1 ) itrs_mult(istep) = itrs_mult(istep) + 10
           end do

           npw_k1 = npwarr(ikpt1i)
!          npw_k2 = npwarr(ikpt2i)

           do istep = 1, berrystep
             npw_k3(istep)=npwarr(ikpt3i(istep))
           end do

!          ji: the loop is over the FBZ, but sflag and smat only apply to the IBZ
           if (dtset%berryopt == 4 .and. inibz == 1) then
             ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
             smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt1i_sp,ifor,idir)
           else
             smat_k(:,:,:) = zero
           end if

!          pwind_k(1:npw_k1) = pwind(ikg+1:ikg+npw_k1,ifor,idir)
           pwnsfac_k(1,1:npw_k1) = pwnsfac(1,ikg+1:ikg+npw_k1)
           pwnsfac_k(2,1:npw_k1) = pwnsfac(2,ikg+1:ikg+npw_k1)

!          the array needed to compute the overlap matrix between k and k+istep*dk (with multiple steps)
!          the 0-case (no corresponding pw in k and k+dk) could be handled better (k+2*dk could have a corresponding pw ?)
           pwind_k_mult(1:npw_k1,1)=pwind(ikg+1:ikg+npw_k1,ifor,idir)
           do istep = 1, berrystep-1
             do jj=1, npw_k1
               if(pwind_k_mult(jj,istep)/=0)then
                 pwind_k_mult(jj,istep+1) = pwind(dtefield%fkgindex(ikpt3(istep))+pwind_k_mult(jj,istep),ifor,idir)
               else
                 pwind_k_mult(jj,istep+1) = 0
               end if
             end do
           end do

!          DEBUG
!          write(std_out,*)' berryphase_new : dtset%berryopt,inibz,ikpt1i,isppol,dtset%nkpt,ifor,idir', &
!          &          dtset%berryopt,inibz,ikpt1i,isppol,dtset%nkpt,ifor,idir
!          write(std_out,'(a,4i4)' )' berryphase_new : sflag_k(:)=',sflag_k(:)
!          ENDDEBUG

           if (dtset%berryopt == 4 .and. inibz == 1) then
             ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
             sflag_k(:) = dtefield%sflag(:,ikpt1i_sp,ifor,idir)
           else
             sflag_k(:) = 0
           end if

           if (usepaw == 1) then
             icp1=dtefield%cprjindex(ikpt1i,isppol)
             call cprj_get(atindx1,cprj_k,cprj,natom,1,icp1,ikpt1i,0,isppol,&
&             mband,mkmem,mpi_enreg,natom,dtefield%nband_occ,dtefield%nband_occ,&
&             my_nspinor,nsppol,0)

             if ( ikpt1i /= ikpt1 ) then
               call cprj_copy(cprj_k,cprj_ikn)
               call sym_cprj_kn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,indlmn,&
&               dtefield%indkk_f2ibz(ikpt1,2),dtefield%indkk_f2ibz(ikpt1,6),&
&               dtefield%fkptns(:,dtefield%i2fbz(ikpt1i)),&
&               dtefield%lmax,dtefield%lmnmax,natom,dtefield%nband_occ,my_nspinor,&
&               dtefield%nsym,ntypat,typat,dtefield%zarot)
               call cprj_copy(cprj_fkn,cprj_k)
             end if

           end if ! end if usepaw

!          DEBUG
!          write(std_out,'(a,4i4)' )' berryphase_new : sflag_k(:)=',sflag_k(:)
!          ENDDEBUG
           
!          DEBUG
!          write(std_out,'(a,7i4)')'me, idir,ifor, ikpt_loc, ikpt1, isppol = ',&
!          & me,idir,ifor,ikpt_loc,ikpt1,isppol
!          write(std_out,'(a,10i3)')'pwind_k(1:10) = ',pwind_k(1:10)
!          ENDDEBUG
           
           do istep=1,berrystep
             sflag_k_mult(:,istep) = sflag_k(:)
           end do
           
         end if ! end check that ikpt1 > 0 and isppol > 0

!        --------------------------------------------------------------------------------
!        Communication
!        --------------------------------------------------------------------------------

         do istep=1,berrystep

!          if(ikpt_loc <= nsppol*dtefield%fmkmem) then
           if (ikpt1 > 0 .and. isppol > 0) then ! I currently have a true kpt to use

             count = npw_k3(istep)*nband_k
             ABI_ALLOCATE(cgq,(2,count))

             if (mpi_enreg%paral_compil_kpt == 1) then
               if(minval(abs(mpi_enreg%proc_distrb(ikpt3i(istep),1:nband_k,isppol) - me)) == 0) then
!                I already have the datas
                 source = me
               else
!                I need the datas from someone else
                 source = mpi_enreg%proc_distrb(ikpt3i(istep),1,isppol)
               end if
             else
               source = me
             end if
           else
             source = -1 ! I do not have a kpt to use
           end if

           do dest = 0, nproc-1

             if ((dest==me) .and. (ikpt1>0) .and. (isppol>0)) then
!              I am destination and I have something to do
!              if (mpi_enreg%paral_compil_kpt == 1) write(std_out,*) &
!              &               'coucou 2, mpi_enreg%proc_distrb(ikpt3i(istep),1:nband_k,isppol) : ', &
!              &               mpi_enreg%proc_distrb(ikpt3i(istep),1:nband_k,isppol)
!              write(std_out,*)'ikpt3i(istep) ', ikpt3i(istep)
!              write(std_out,*)'nband_k ',nband_k
!              write(std_out,*)'isppol ', isppol
!              write(std_out,*)'mpi_enreg%proc_distrb',mpi_enreg%proc_distrb

               if (source == me) then
!                I am destination and source
!                DEBUG
!                write(std_out,*)'copying ... '
!                write(std_out,*)'me: ',me, 'ikpt3i(istep) ', ikpt3i(istep), 'isppol ', isppol
!                ENDDEBUG

!                pwnsfac
                 idum = dtefield%fkgindex(ikpt3(istep))
                 pwnsfac_k(3,1:npw_k3(istep)) = pwnsfac(1,idum+1:idum+npw_k3(istep))
                 pwnsfac_k(4,1:npw_k3(istep)) = pwnsfac(2,idum+1:idum+npw_k3(istep))

!                cgq (and cprj)
                 icg1 = dtefield%cgindex(ikpt3i(istep),isppol)

                 if (usepaw == 1) then
                   icp2=dtefield%cprjindex(ikpt3i(istep),isppol)
                   call cprj_get(atindx1,cprj_kb,cprj,natom,1,icp2,ikpt3i(istep),0,isppol,&
&                   mband,mkmem,mpi_enreg,natom,dtefield%nband_occ,dtefield%nband_occ,&
&                   my_nspinor,nsppol,0)
                 end if

                 cgq(:,1:count)  = cg(:,icg1+1:icg1+count)
!                if (usepaw == 1) call cprj_copy(cprj_buf,cprj_kb)

!                if ((source /= me)) then
               else
!                I am the destination but not the source -> receive
!                DEBUG
!                write(std_out,'(a)')'receiving ...'
!                write(std_out,'(a,i4,a,i4,a,i4,a,i4)')'me: ',me, 'source ', source,'ikpt3i(istep) ', ikpt3i(istep), 'isppol ', isppol
!                ENDDEBUG

!                receive pwnsfac
                 ABI_ALLOCATE(buffer,(2,npw_k3(istep)))
                 tag = ikpt3(istep) + (isppol - 1)*dtefield%fnkpt
                 call xrecv_mpi(buffer,source,tag,spaceComm,ierr)
                 pwnsfac_k(3,1:npw_k3(istep)) = buffer(1,1:npw_k3(istep))
                 pwnsfac_k(4,1:npw_k3(istep)) = buffer(2,1:npw_k3(istep))
                 ABI_DEALLOCATE(buffer)

!                receive cgq (and cprj)
                 tag = ikpt3i(istep) + (isppol - 1)*nkpt
                 call xrecv_mpi(cgq,source,tag,spaceComm,ierr)

                 if (usepaw == 1) then
                   call cprj_mpi_recv(natom,n2dim,dimlmn,ncpgr,cprj_kb,source,spaceComm,ierr)
                 end if

               end if

             else if (dest /= me) then 

!              jkpt is the kpt which is being treated by dest
!              jsppol is his isppol
               jkpt = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,1)
               jsppol = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,2)

               if (jkpt > 0 .and. jsppol > 0) then ! dest is treating a true kpt

                 jkpt2 = dtefield%ikpt_dk(jkpt,ifor,idir)
                 jkpt2i = dtefield%indkk_f2ibz(jkpt2,1)

!                check if I am his source
                 if((mpi_enreg%proc_distrb(jkpt2i,1,jsppol) == me))  then
!                  I know something about jkpt3i and I must send it
!                  DEBUG
!                  write(std_out,'(a)')'sending ...'
!                  write(std_out,'(a,i4,a,i4,a,i4,a,i4)')'dest: ',dest,' me: ',me,&
!                  &                          ' jkpt2i ',jkpt2i,' jsppol: ',jsppol
!                  ENDDEBUG

!                  pwnsfac
                   tag = jkpt2 + (jsppol - 1)*dtefield%fnkpt
                   count1 = npwarr(jkpt2i)
                   ABI_ALLOCATE(buffer,(2,count1))
                   idum = dtefield%fkgindex(jkpt2)
                   buffer(1,1:count1)  = pwnsfac(1,idum+1:idum+count1)
                   buffer(2,1:count1)  = pwnsfac(2,idum+1:idum+count1)
                   call xsend_mpi(buffer,dest,tag,spaceComm,ierr)
                   ABI_DEALLOCATE(buffer)

!                  cgq (and cprj)
                   icg1 = dtefield%cgindex(jkpt2i,jsppol)

                   if (usepaw == 1) then
                     icp2=dtefield%cprjindex(jkpt2i,jsppol)
                     call cprj_get(atindx1,cprj_buf,cprj,natom,1,icp2,jkpt2i,0,jsppol,&
&                     mband,mkmem,mpi_enreg,natom,dtefield%nband_occ,dtefield%nband_occ,&
&                     my_nspinor,nsppol,0)
                   end if

                   tag = jkpt2i + (jsppol - 1)*nkpt
                   count1 = npwarr(jkpt2i)*nband_k
                   ABI_ALLOCATE(buffer,(2,count1))
                   buffer(:,1:count1)  = cg(:,icg1+1:icg1+count1)
                   call xsend_mpi(buffer,dest,tag,spaceComm,ierr)
                   ABI_DEALLOCATE(buffer)

                   if (usepaw == 1 ) then
!                    note that gradients are not sent, may change this later
                     call cprj_mpi_send(natom,n2dim,dimlmn,ncpgr,cprj_buf,dest,spaceComm,ierr)
                   end if

                 end if ! end check that I am his source

               end if ! end check that jkpt > 0 and jsppol > 0

             end if ! end if statements on dest == me or dest /= me

           end do  ! end loop over dest = 0, nproc - 1

           if (ikpt1 > 0 .and. isppol > 0) then ! if I am treating a kpt, compute the smatrix

             if (usepaw == 1) then
               if (ikpt3(istep) /= ikpt3i(istep)) then ! cprj_kb refers to ikpt3i(istep), must compute ikpt3(istep) value
                 call cprj_copy(cprj_kb,cprj_ikn)

                 call sym_cprj_kn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,indlmn,&
&                 dtefield%indkk_f2ibz(ikpt3(istep),2),dtefield%indkk_f2ibz(ikpt3(istep),6),&
&                 dtefield%fkptns(:,dtefield%i2fbz(ikpt3i(istep))),&
&                 dtefield%lmax,dtefield%lmnmax,natom,&
&                 dtefield%nband_occ,my_nspinor,dtefield%nsym,ntypat,typat,&
&                 dtefield%zarot)
                 call cprj_copy(cprj_fkn,cprj_kb)
               end if
               call smatrix_k_paw(cprj_k,cprj_kb,dtefield,idir,ifor,natom,smat_k_paw,typat)
             end if

             icg1 = 0
             icg = dtefield%cgindex(ikpt1i,isppol)
!            DEBUG
!            if(istep<=2)then
!            if(ikpt1==1)then
!            write(std_out,'(a,2i4,3e15.4)')'istep ikpt3, kpt, cgq', istep, ikpt3(istep), dtefield%fkptns(:,ikpt3(istep))
!            write(std_out,*) cgq
!            write(std_out,*)
!            end if
!            end if
!            ENDDEBUG
             call smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs_mult(istep),job,maxbd,&
&             mcg,count,mcg1_k,minbd,&
&             mpw,dtefield%nband_occ,&
&             npw_k1,npw_k3(istep),my_nspinor,pwind_k_mult(:,istep),pwnsfac_k,sflag_k_mult(:,istep),&
&             shiftbd,smat_inv,smat_k,smat_k_paw,usepaw)

             if ((job == 10).or.(job == 11)) then

               if (sqrt(dtm_k(1)*dtm_k(1) + dtm_k(2)*dtm_k(2)) < tol12) then
                 write(message,'(a,a,a,a,i5,a,a,a)')ch10,&
&                 ' berryphase_new : BUG - ',ch10,&
&                 '  For k-point #',ikpt1,',',ch10,&
&                 '  the determinant of the overlap matrix is found to be 0.'
                 call wrtout(std_out,message,'PERS')
                 call leave_new('PERS')
               end if

               dtm_mult(1,ikpt1+(isppol-1)*dtefield%fnkpt,istep) = dtm_k(1)
               dtm_mult(2,ikpt1+(isppol-1)*dtefield%fnkpt,istep) = dtm_k(2)

             end if

             if (dtset%berryopt == 4 .and. inibz == 1 .and. istep == 1) then
               ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
               dtefield%smat(:,:,:,ikpt1i_sp,ifor,idir) = &
&               smat_k(:,:,:)
               dtefield%sflag(:,ikpt1i_sp,ifor,idir) = &
&               sflag_k_mult(:,1)
             end if

             if (((job == 1).or.(job == 11)).and.(inibz == 1) .and. istep == 1) then
               cg1(:,icg + 1: icg + npw_k1*nband_k) = &
               cg1(:,icg + 1:icg + npw_k1*nband_k) + &
               dkinv*cg1_k(:,1:npw_k1*nband_k)
             end if

             ABI_DEALLOCATE(cgq)

           end if ! end if ikpt1 > 0 and isppol > 0

         end do ! end loop over istep

!        if (ikpt_loc <= dtefield%fmkmem) sflag_k(:) = sflag_k_mult(:,1)
         if (ikpt1 > 0) sflag_k(:) = sflag_k_mult(:,1)

       end do ! close loop over ikpt_loc (k-points, isppol)

       ABI_DEALLOCATE(smat_inv)
       ABI_DEALLOCATE(smat_k)
       ABI_DEALLOCATE(smat_k_paw)

     end do   ! close loop over ifor

     if (mpi_enreg%paral_compil_kpt == 1) then
       ABI_ALLOCATE(buffer1,(2*dtefield%fnkpt*nsppol*berrystep))
       ABI_ALLOCATE(buffer2,(2*dtefield%fnkpt*nsppol*berrystep))
       count = 2*dtefield%fnkpt*nsppol*berrystep
       buffer1(:) = reshape(dtm_mult(:,:,:),(/count/))
       call xsum_mpi(buffer1,buffer2,count,spaceComm,ierr)
       dtm_mult(:,:,:) = reshape(buffer2(:),(/2,dtefield%fnkpt*nsppol,berrystep/))
       ABI_DEALLOCATE(buffer1)
       ABI_DEALLOCATE(buffer2)
     end if

!    DEBUG
!    write(std_out,*)
!    write(std_out,*)'istep = 1, nsppol =',nsppol
!    istep=1
!    isppol=1
!    do jkpt = 1, dtefield%fnkpt
!    write(std_out,'(a,i4,3e15.4,2e15.4)')'jkpt, kpt, dtm_mult(:,kpt,1)', jkpt, dtefield%fkptns(:,jkpt),  dtm_mult(:,jkpt+(isppol-1)*dtefield%fnkpt,istep)
!    end do
!    write(std_out,*)
!    write(std_out,*) "istep = 2"
!    if(berrystep>=2)then
!    istep=2
!    isppol=1
!    do jkpt = 1, dtefield%fnkpt
!    write(std_out,'(a,i4,3e15.4,2e15.4)')'jkpt, kpt, dtm_mult(:,kpt,2)', jkpt, dtefield%fkptns(:,jkpt),  dtm_mult(:,jkpt+(isppol-1)*dtefield%fnkpt,istep)
!    end do
!    end if
!    ENDDEBUG

!    ===========================================================================

!    Compute the Berry phase polarization

     if ((option == 1).or.(option == 3)) then

!      Compute the electronic Berry phase

       polb_mult(:,:)=zero
       do istep = 1,berrystep

         if(berrystep==1) then

           write(message,'(a,a)')ch10,&
&           ' Compute the electronic contribution to polarization'
           call wrtout(std_out,message,'COLL')

         else
           write(message,'(a,a,i4,a)')ch10,&
&           ' Compute the electronic contribution to polarization for a step of istep=',&
&           istep,'*dk'
           call wrtout(std_out,message,'COLL')
         end if

         if(istep /= 1) then
!          construct the strings for a step of istep*dk
!          string lenght
           istr=1
           nkstr=1
           ikpt1=1
           do ikpt=1,dtefield%fnkpt
             do jstep = 1,istep
               ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
             end do
             if (ikpt1 == 1) exit
             nkstr = nkstr + 1
           end do
!          Check that the string length is a divisor of nkpt
           if(mod(dtefield%fnkpt,nkstr) /= 0) then
             write(message,'(a,a,a,a,i5,a,i5,a,i7)')ch10,&
&             ' berryphase_new: BUG -',ch10,&
&             '  For istep = ', istep,&
&             '  The string length = ',nkstr,&
&             ', is not a divisor of fnkpt =',dtefield%fnkpt
             call wrtout(std_out,message,'COLL')
             call leave_new('COLL')
           end if
           nstr = dtefield%fnkpt/nkstr

           write(message,'(a,i1,a,i2,a,i3,a,i3)')&
&           '  berryphase_new: for direction ',idir, ' and istep ', istep, ', nkstr = ',nkstr,&
&           ', nstr = ',nstr
           call wrtout(std_out,message,'COLL')
           call wrtout(ab_out,message,'COLL')

           ABI_ALLOCATE(idxkstr_mult,(nkstr,nstr))
           iunmark = 1
           kpt_mark(:)=0
           do istr=1,nstr
             do while(kpt_mark(iunmark) /= 0)
               iunmark = iunmark + 1
             end do
             idxkstr_mult(1,istr) = iunmark
             kpt_mark(iunmark)=1

             ikpt1 = idxkstr_mult(1,istr)
             do jkstr=2, nkstr
               do jstep = 1, istep
                 ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
               end do
               idxkstr_mult(jkstr,istr) = ikpt1
               kpt_mark(ikpt1) = 1
             end do
           end do
         else
           nstr = dtefield%nstr(idir)
           nkstr = dtefield%nkstr(idir)
           ABI_ALLOCATE(idxkstr_mult,(nkstr,nstr))
           idxkstr_mult(:,:) = dtefield%idxkstr(1:nkstr,1:nstr,idir)
         end if
!        DEBUG
!        do istr=1,nstr
!        write(std_out,*)'string ', idxkstr_mult(:,istr)
!        end do
!        ENDBEBUG

         ABI_ALLOCATE(det_string,(2,nstr))
         ABI_ALLOCATE(polberry,(nstr))
         write(message,'(a,10x,a,10x,a)')ch10,&
&         'istr','polberry(istr)'
         call wrtout(std_out,message,'COLL')

         polbtot = zero
         do isppol = 1, nsppol

           det_string(1,:) = one ; det_string(2,:) = zero
           dtm_k(:) = one
           det_average(:) = zero

           do istr = 1, nstr
             do jkstr = 1, nkstr

               ikpt=idxkstr_mult(jkstr,istr)

               dtm_real=dtm_mult(1,ikpt+(isppol-1)*dtefield%fnkpt,istep)
               dtm_imag=dtm_mult(2,ikpt+(isppol-1)*dtefield%fnkpt,istep)

               dtm_k(1) = det_string(1,istr)*dtm_real - &
&               det_string(2,istr)*dtm_imag
               dtm_k(2) = det_string(1,istr)*dtm_imag + &
&               det_string(2,istr)*dtm_real
               det_string(1:2,istr) = dtm_k(1:2)
!              DEBUG
!              write(std_out,'(a,i4,3e15.4,2e15.4)')'ikpt, kpt, dtm', ikpt, dtefield%fkptns(:,ikpt),  dtm_k
!              ENDDEBUG

             end do

             det_average(:) = det_average(:) + &
&             det_string(:,istr)/dble(nstr)

           end do


!          correction to obtain a smouth logarithm of the determinant
           ABI_ALLOCATE(str_flag,(nstr))
!          DEBUG
!          since we don't have any case of non-nul Chern number,
!          we must change the det_string value "by brute force" if we want debug this
!          allocate(det_string_test(2,dtefield%nstr(idir)))
!          det_string_test(:,:)=det_string(:,:)
!          kk=0
!          det_string(1,1)=cos(2._dp*Pi*real(kk,dp)/real(4,dp))
!          det_string(2,1)=sin(2._dp*Pi*real(kk,dp)/real(4,dp))
!          jj=dtefield%str_neigh(1,1,idir)
!          ll=dtefield%str_neigh(2,1,idir)
!          do while (jj/=1)
!          kk=kk+1
!          det_string(1,jj)=cos(2._dp*Pi*real(kk,dp)/real(4,dp))
!          det_string(2,jj)=sin(2._dp*Pi*real(kk,dp)/real(4,dp))
!          det_string(1,ll)=cos(-2._dp*Pi*real(kk,dp)/real(4,dp))
!          det_string(2,ll)=sin(-2._dp*Pi*real(kk,dp)/real(4,dp))
!          jj=dtefield%str_neigh(1,jj,idir)
!          ll=dtefield%str_neigh(2,ll,idir)
!          enddo
!          ENDDEBUG
           if (istep==1) then
             do ineigh_str = 1,2
               str_flag(:)=0
               delta_str(:) = &
&               dtefield%coord_str(:,dtefield%str_neigh(ineigh_str,1,idir),idir) - dtefield%coord_str(:,1,idir)
               dstr(:)= delta_str(:) - nint(delta_str(:)) - real(dtefield%strg_neigh(ineigh_str,1,:,idir),dp)
               dist_=0._dp
               do kk = 1,2
                 do jj = 1,2
                   dist_ = dist_ + dstr(kk)*dtefield%gmet_str(kk,jj,idir)*dstr(jj)
                 end do
               end do
               dist_=sqrt(dist_)
               do istr = 1,dtefield%nstr(idir)
                 if(str_flag(istr)==0)then
!                  write(std_out,*)'new string'
                   str_flag(istr)=1
                   call rhophi(det_string(:,istr),dphase,rho)
!                  write(std_out,'(i4,e15.4,e15.4,e15.4)')istr, det_string(:,istr),dphase
                   dphase_init=dphase
                   jstr = dtefield%str_neigh(ineigh_str,istr,idir)
                   do while (istr/=jstr)
                     str_flag(jstr)=1
                     call rhophi(det_string(:,jstr),dphase_new,rho)
                     jj=nint((dphase_new-dphase)/(2._dp*Pi))
!                    write(std_out,'(i4,e15.4,e15.4,e15.4,e15.4,i4)')jstr, det_string(:,jstr),dphase_new,dphase_new-dphase,jj
                     dphase_new=dphase_new-2._dp*Pi*real(jj,dp)
                     if(jj/=0)then
                       write(message,'(6a)') ch10,&
&                       ' berryphase_new : WARNING -',ch10,&
&                       '  the berry phase has some huge variation in the space of strings of k-points',ch10,&
&                       '  ABINIT is trying to correct the berry phase, but it is highly experimental'
                       call wrtout(std_out,message,'PERS')
                     end if
!                    if(jj/=0)write(std_out,'(i4,e15.4,e15.4,e15.4,e15.4)')jstr, det_string(:,jstr),dphase_new,dphase_new-dphase
                     dphase=dphase_new
                     jstr=dtefield%str_neigh(ineigh_str,jstr,idir)
                   end do
!                  write(std_out,*)dphase_init, dphase, (dphase-dphase_init)/(2._dp*Pi),nint((dphase-dphase_init)/(2._dp*Pi))
                 end if
               end do
             end do
           end if
           ABI_DEALLOCATE(str_flag)
!          DEBUG
!          deallocate(dist_str)
!          det_string(:,:)=det_string_test(:,:)
!          deallocate(det_string_test)
!          ENDDEBUG

!          First berry phase that corresponds to det_average
!          phase0 = atan2(det_average(2),det_average(1))
           call rhophi(det_average,phase0,rho)
           det_mod = det_average(1)**2+det_average(2)**2

!          Then berry phase that corresponds to each string relative to the average
           do istr = 1, nstr

             rel_string(1) = (det_string(1,istr)*det_average(1) + &
             det_string(2,istr)*det_average(2))/det_mod
             rel_string(2) = (det_string(2,istr)*det_average(1) - &
             det_string(1,istr)*det_average(2))/det_mod
!            dphase = atan2(rel_string(2),rel_string(1))
             call rhophi(rel_string,dphase,rho)
             polberry(istr) = dtefield%sdeg*(phase0 + dphase)/two_pi
             polb_mult(isppol,istep) = polb_mult(isppol,istep) + polberry(istr)/(istep*dtefield%nstr(idir))
             polb(isppol) = zero
             do jstep=1, istep
               polb(isppol)=polb(isppol)+coef(jstep,istep)*polb_mult(isppol,jstep)
             end do

             write(message,'(10x,i4,7x,e16.9)')istr,polberry(istr)
             call wrtout(std_out,message,'COLL')

           end do

           if(berrystep>1)then
             write(message,'(9x,a,7x,e16.9,1x,a,i4,a,i4,a)')&
&             'total',polb_mult(isppol,istep),'(isppol=',isppol,', istep=',istep,')'!,ch10
             call wrtout(std_out,message,'COLL')

             write(message,'(3x,a,7x,e16.9,1x,a,i4,a,i4,a,a)')&
&             '+correction',polb(isppol),'(isppol=',isppol,', istep=1..',istep,')',ch10
             call wrtout(std_out,message,'COLL')

           else

             write(message,'(9x,a,7x,e16.9,1x,a,i4,a)')&
&             'total',polb_mult(isppol,istep),'(isppol=',isppol,')'!,ch10
             call wrtout(std_out,message,'COLL')
           end if

           polbtot = polbtot + polb(isppol)

         end do    ! isppol

!        Fold into interval [-1,1]
         polbtot = polbtot - 2_dp*nint(polbtot/2_dp)

         ABI_DEALLOCATE(det_string)
         ABI_DEALLOCATE(polberry)

!        ==========================================================================

!        Compute the ionic Berry phase

         call xredxcart(natom,1,rprimd,xcart,xred)
         politot = zero
         write(message,'(a)')' Compute the ionic contributions'
         call wrtout(std_out,message,'COLL')

         write(message,'(a,2x,a,2x,a,15x,a)')ch10,&
&         'itom', 'itypat', 'polion'
         call wrtout(std_out,message,'COLL')

         do iatom = 1, natom
           itypat = typat(iatom)

!          The ionic phase can be computed much easier
           polion = zion(itypat)*xred(idir,iatom)

!          Fold into interval (-1,1)
           polion = polion - 2_dp*nint(polion/2_dp)
           politot = politot + polion
           write(message,'(2x,i2,5x,i2,10x,e16.9)') iatom,itypat,polion
           call wrtout(std_out,message,'COLL')
         end do

!        Fold into interval [-1,1] again
         politot = politot - 2_dp*nint(politot/2_dp)
         pion(idir) = politot

         write(message,'(9x,a,7x,es19.9)') 'total',politot
         call wrtout(std_out,message,'COLL')


!        ==========================================================================

!        Compute the total polarization

         poltot = politot + polbtot

         if (berrystep==1)then
           write(message,'(a,a)')ch10,&
&           ' Summary of the results'
           call wrtout(std_out,message,'COLL')
           if (unit_out /= 0) then
             call wrtout(unit_out,message,'COLL')
           end if
         else
           write(message,'(a,a,i4)')ch10,&
&           ' Summary of the results for istep =',istep
           call wrtout(std_out,message,'COLL')
           if (unit_out /= 0) then
             call wrtout(unit_out,message,'COLL')
           end if
         end if

         write(message,'(a,es19.9)')&
&         ' Electronic Berry phase ' ,polbtot
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if

         write(message,'(a,es19.9)') &
&         '            Ionic phase ', politot
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if

         write(message,'(a,es19.9)') &
&         '            Total phase ', poltot
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if

         poltot = poltot - 2.0_dp*nint(poltot/2._dp)
         write(message,'(a,es19.9)') &
&         '    Remapping in [-1,1] ', poltot
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if

!        Transform the phase into a polarization
         fac = 1._dp/(gmod*dtefield%nkstr(idir))
         fac = fac/ucvol
         pol = fac*poltot

         write(message,'(a,a,es19.9,a,a,a,es19.9,a,a)')ch10,&
&         '           Polarization ', pol,' (a.u. of charge)/bohr^2',ch10,&
&         '           Polarization ', pol*(e_Cb)/(Bohr_Ang*1d-10)**2,&
&         ' C/m^2',ch10
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if


         ABI_DEALLOCATE(idxkstr_mult)

       end do !istep

       pel(idir) = polbtot

     end if   ! option == 1 or option == 3

!    Write the ddk WF to a file

     if ((option == 2).or.(option == 3)) then

       pertcase = idir + 3*natom
       response = 1
       call appdig(pertcase,dtfil%fnameabo_1wf,fiwf1o)
       mxfh = 0 ; nxfh = 0 ; nstep = 1
       ABI_ALLOCATE(xfhist,(3,natom+4,2,mxfh))
       ABI_ALLOCATE(resid,(mband*nkpt*nsppol))
       xfhist(:,:,:,:) = zero
       resid(:) = zero

       call outwf(cg1,dtset,eig_dum,fiwf1o,hdr,kg,dtset%kptns,&
&       mband,mcg,mkmem,mpi_enreg,mpw,mxfh,natom,dtset%nband,&
&       nkpt,npwarr,nsppol,nstep,&
&       nxfh,occ_dum,resid,response,dtfil%unwff2,wffnow,wfs,wvl,xfhist)
       ABI_DEALLOCATE(xfhist)
       ABI_DEALLOCATE(resid)

     end if  ! option == 2 or option == 3

   end if   ! rfdir(idir) == 1

 end do    ! Close loop over idir

!Compute polarization in cartesian coordinates
 if ((dtset%rfdir(1) == 1).and.(dtset%rfdir(2) == 1).and.&
& (dtset%rfdir(3) == 1)) then

   if(usepaw.ne.1) then
     pelev=zero
   else
     call pawpolev(natom,ntypat,pawrhoij,pawtab,pelev,typat)
   end if

   call polcart(pel,pel_cart,pelev,pion,pion_cart,3,&
&   ptot_cart,rprimd,ucvol,unit_out)
   call polcart(pel,pel_cart,pelev,pion,pion_cart,3,&
&   ptot_cart,rprimd,ucvol,6)

 end if

!deallocate(pwind_k, dtm)
 ABI_DEALLOCATE(pwnsfac_k)
 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(cg1_k)
 if (option > 1)  then
   ABI_DEALLOCATE(cg1)
   ABI_DEALLOCATE(eig_dum)
   ABI_DEALLOCATE(occ_dum)
 end if

 if (usepaw == 1) then
   ABI_DEALLOCATE(dimlmn)
   call cprj_free(cprj_k)
   call cprj_free(cprj_kb)
   call cprj_free(cprj_gat)
   ABI_DEALLOCATE(cprj_k)
   ABI_DEALLOCATE(cprj_kb)
   ABI_DEALLOCATE(cprj_gat)
   if (dtset%kptopt /= 3) then
     call cprj_free(cprj_ikn)
     call cprj_free(cprj_fkn)
     ABI_DEALLOCATE(cprj_ikn)
     ABI_DEALLOCATE(cprj_fkn)
   end if

   if ( mpi_enreg%paral_compil_kpt==1 ) then
     call cprj_free(cprj_buf)
     ABI_DEALLOCATE(cprj_buf)
   end if 

 end if


 ABI_DEALLOCATE(ikpt3)
 ABI_DEALLOCATE(ikpt3i)
 ABI_DEALLOCATE(sflag_k_mult)
 ABI_DEALLOCATE(npw_k3)
 ABI_DEALLOCATE(pwind_k_mult)
 ABI_DEALLOCATE(itrs_mult)
 ABI_DEALLOCATE(coef)
 ABI_DEALLOCATE(polb_mult)
 ABI_DEALLOCATE(dtm_mult)


!DEBUG
!write(std_out,*)'berryphase_new exit'
!END_DEBUG
end subroutine berryphase_new
!!***
