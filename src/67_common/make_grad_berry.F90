!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_grad_berry.F90
!! NAME
!! make_grad_berry
!!
!! FUNCTION
!! compute gradient contribution from berry phase in finite
!! electric field case
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)=input wavefunctions
!!  cgq(2,mcgq) = wavefunctions at neighboring k points
!!  cprj_k(natom,nband_k*usepaw)=cprj at this k point
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dimlmn(natom)=lmn_size for each atom in input order
!!  dimlmn_srt(natom)=lmn_size for each atom sorted by type
!!  direc(2,npw*nspinor)=gradient vector
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  iband=index of band currently being treated
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point currently being treated
!!  isppol=spin polarization currently treated
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  natom=number of atoms in cell.
!!  matblk=dimension of the array ph3d
!!  mband =maximum number of bands
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw
!!  mcg=second dimension of the cg array
!!  mcgq=second dimension of the cgq array
!!  mgfft=maximum size of 1D FFTs
!!  mkgq = second dimension of pwnsfacq
!!  nkpt=number of k points
!!  mpi_enreg=informations about MPI parallelization
!!  nloalg(5) data concerning nonlop application
!!  npw=number of planewaves in basis sphere at given k.
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=number of spin polarizations
!!  ntypat=number of types of atoms in cell.
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwnsfacq(2,mkgq) = phase factors for the nearest neighbours of the
!!                     current k-point (electric field, MPI //)
!!
!! OUTPUT
!! grad_berry(2,npw*nspinor) :: contribution to gradient in finite electric field case
!!
!! SIDE EFFECTS
!!  dtefield <type(efield_type)> = variables related to Berry phase calculations (see initberry.f)
!!
!! NOTES
!!
!! PARENTS
!!      cgwf
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,cprj_get,nonlop,smatrix,smatrix_k_paw
!!      sym_cprj_kn
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine make_grad_berry(cg,cgq,cprj_k,detovc,dimffnl,dimlmn,dimlmn_srt,direc,dtefield,ffnl,grad_berry,&
&                          gs_hamk,iband,icg,ikpt,isppol,kg_k,lmnmax,matblk,mband,&
&                          mcg,mcgq,mgfft,mkgq,mpi_enreg,mpsang,mpssoang,mpw,natom,nkpt,&
&                          nloalg,npw,npwarr,nspinor,nsppol,ntypat,pwind,ph3d,pwind_alloc,&
&                          pwnsfac,pwnsfacq)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_grad_berry'
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => make_grad_berry
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!scalars
 integer,intent(in) :: dimffnl,iband,icg,ikpt,isppol,lmnmax,matblk,mband,mcg,mcgq,mgfft
 integer,intent(in) :: mkgq,mpw,mpsang,mpssoang,natom,nkpt,npw,nspinor,nsppol,ntypat,pwind_alloc
 type(gs_hamiltonian_type),intent(in) :: gs_hamk
 type(efield_type),intent(inout) :: dtefield
 type(MPI_type),intent(inout) :: mpi_enreg

!arrays
 integer,intent(in) :: dimlmn(natom),dimlmn_srt(natom),kg_k(3,npw),nloalg(5)
 integer,intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp),intent(in) :: cg(2,mcg),cgq(2,mcgq)
 real(dp),intent(inout) :: direc(2,npw*nspinor)
 real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
 real(dp),intent(inout) :: ph3d(2,npw,matblk)
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)
 real(dp),intent(out) :: detovc(2,2,3),grad_berry(2,npw*nspinor)
 type(cprj_type),intent(in) :: cprj_k(natom,dtefield%nband_occ*gs_hamk%usepaw)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,ddkflag,dimenlc1,dimenlr1,iatom,icg1,icp2,idum1
 integer :: idir,ifor,ikgf,ikptf,ikpt2,ikpt2f,ipw,i_paw_band,ispinor,itrs,itypat,job
 integer :: klmn,mcg1_k,mcg_q,nbo,nkpg,npw_k2,nspinortot,paw_opt,shiftbd,signs
 real(dp) :: fac
 character(len=500) :: message
!arrays
 integer :: pwind_k(npw),sflag_k(dtefield%nband_occ)
 real(dp) :: cg1_k(2,npw*nspinor),dtm_k(2),pwnsfac_k(4,mpw),sij_dum(2,0)
 real(dp) :: smat_k(2,dtefield%nband_occ,dtefield%nband_occ)
 real(dp) :: smat_inv(2,dtefield%nband_occ,dtefield%nband_occ),svectout_dum(2,0)
 real(dp),allocatable :: cgq_k(:,:),enl_rij(:,:,:),grad_berry_ev(:,:)
 real(dp),allocatable :: kpg(:,:),qijbkk(:,:,:),smat_k_paw(:,:,:)
 type(cprj_type) :: cprj_dum(1,1)
 type(cprj_type),allocatable :: cprj_kb(:,:),cprj_band_srt(:,:)
 type(cprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)


! *********************************************************************

!DBG_ENTER("COLL")

 nbo = dtefield%nband_occ

!allocations

!Electric field: compute the gradient of the Berry phase part of the energy functional.
!See PRL 89, 117602 (2002), grad_berry(:,:) is the second term of Eq. (4)
 grad_berry(:,:) = zero
 job = 11 ; shiftbd = 1
 mcg_q = mpw*mband*nspinor
 mcg1_k = npw*nspinor

 if (gs_hamk%usepaw /= 0) then
   dimenlr1 = lmnmax*(lmnmax+1)/2
   dimenlc1 = 2*dimenlr1
   ABI_ALLOCATE(qijbkk,(dimenlc1,natom,nspinor**2))
   ABI_ALLOCATE(enl_rij,(dimenlr1,natom,nspinor**2))
   ABI_ALLOCATE(smat_k_paw,(2,nbo,nbo))
   ABI_ALLOCATE(grad_berry_ev,(2,npw*nspinor))
   enl_rij = zero
   qijbkk = zero
   smat_k_paw = zero
   ABI_ALLOCATE(cprj_kb,(natom,nbo))
   call cprj_alloc(cprj_kb,0,dimlmn)
   ABI_ALLOCATE(cprj_band_srt,(natom,1))
   call cprj_alloc(cprj_band_srt,0,dimlmn_srt)
   if (nkpt /= dtefield%fnkpt) then
     ABI_ALLOCATE(cprj_fkn,(natom,nbo))
     ABI_ALLOCATE(cprj_ikn,(natom,nbo))
     call cprj_alloc(cprj_fkn,0,dimlmn)
     call cprj_alloc(cprj_ikn,0,dimlmn)
   end if
   nkpg = 0
   ABI_ALLOCATE(kpg,(npw,nkpg))
 end if

 ikptf = dtefield%i2fbz(ikpt)
 ikgf = dtefield%fkgindex(ikptf)  ! this is the shift for pwind

 do idir = 1, 3
!  skip idir values for which efield_dot(idir)=0
   if (abs(dtefield%efield_dot(idir)) < tol12) cycle
!  Implicitly, we use the gradient multiplied by the number of k points in the FBZ
   fac = dtefield%efield_dot(idir)*dble(dtefield%fnkpt)/&
&   (dble(dtefield%nstr(idir))*four_pi)
   do ifor = 1, 2
!    Handle dtefield%i2fbz properly and ask whether t.r.s. is used
     ikpt2f = dtefield%ikpt_dk(ikptf,ifor,idir)
     if (dtefield%indkk_f2ibz(ikpt2f,6) == 1) then
       itrs = 10
     else
       itrs = 0
     end if
     ikpt2 = dtefield%indkk_f2ibz(ikpt2f,1)
     npw_k2 = npwarr(ikpt2)
     ABI_ALLOCATE(cgq_k,(2,nbo*nspinor*npw_k2))
     pwind_k(1:npw) = pwind(ikgf+1:ikgf+npw,ifor,idir)
     pwnsfac_k(1:2,1:npw) = pwnsfac(1:2,ikgf+1:ikgf+npw)
     sflag_k(:) = dtefield%sflag(:,ikpt+(isppol-1)*nkpt,ifor,idir)
     smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt+(isppol-1)*nkpt,ifor,idir)
     if (mpi_enreg%paral_compil_kpt == 1) then
       icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
       cgq_k(:,1:nbo*nspinor*npw_k2) = &
&       cgq(:,icg1+1:icg1+nbo*nspinor*npw_k2)
       idum1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
       pwnsfac_k(3:4,1:npw_k2) = pwnsfacq(1:2,idum1+1:idum1+npw_k2)
     else
       icg1 = dtefield%cgindex(ikpt2,isppol)
       cgq_k(:,1:nbo*nspinor*npw_k2) = &
&       cg(:,icg1+1:icg1+nbo*nspinor*npw_k2)
       idum1 = dtefield%fkgindex(ikpt2f)
       pwnsfac_k(3:4,1:npw_k2) = pwnsfac(1:2,idum1+1:idum1+npw_k2)
     end if
     if (gs_hamk%usepaw == 1) then
       icp2=nbo*(ikpt2-1)
       call cprj_get(gs_hamk%atindx1,cprj_kb,dtefield%cprj,natom,1,icp2,ikpt,0,isppol,&
&       nbo,dtefield%fnkpt,mpi_enreg,natom,nbo,nbo,nspinor,nsppol,0)
       if (ikpt2 /= ikpt2f) then ! construct cprj_kb by symmetry
         call cprj_copy(cprj_kb,cprj_ikn)
         call sym_cprj_kn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,gs_hamk%indlmn,&
&         dtefield%indkk_f2ibz(ikpt2f,2),dtefield%indkk_f2ibz(ikpt2f,6),&
&         dtefield%fkptns(:,dtefield%i2fbz(ikpt2)),&
&         dtefield%lmax,dtefield%lmnmax,natom,nbo,nspinor,&
&         dtefield%nsym,ntypat,gs_hamk%typat,dtefield%zarot)
         call cprj_copy(cprj_fkn,cprj_kb)
       end if
       call smatrix_k_paw(cprj_k,cprj_kb,dtefield,idir,ifor,natom,smat_k_paw,gs_hamk%typat)
     end if

     icg1 = 0 ; ddkflag = 1
     call smatrix(cg,cgq_k,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,&
&     job,iband,mcg,mcg_q,mcg1_k,iband,mpw,nbo,&
&     npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&     shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
     ABI_DEALLOCATE(cgq_k)
     detovc(:,ifor,idir) = dtm_k(:) !store the determinant of the overlap
     if (sqrt(dtm_k(1)*dtm_k(1) + dtm_k(2)*dtm_k(2)) < tol12) then
       write(message,'(3a,i5,a,i3,a,a,a)') &
&       '  (electric field)',ch10,&
&       '  For k-point #',ikpt,' and band # ',iband,',',ch10,&
&       '  the determinant of the overlap matrix is found to be 0.'
       MSG_BUG(message)
     end if

     if (gs_hamk%usepaw == 1) then
!      this loop applies discretized derivative of projectors
!      note that qijb_kk is sorted by input atom order, but nonlop wants it sorted by type
       do iatom = 1, natom
         itypat = gs_hamk%typat(gs_hamk%atindx1(iatom))
         do klmn = 1, dtefield%lmn2_size(itypat)
!          note: D_ij-like terms have 4 spinor components: 11, 22, 12, and 21. Here the qijb is diagonal
!          in spin space so only the first two are nonzero and they are equal
           do ispinor = 1, nspinor
             qijbkk(2*klmn-1,iatom,ispinor) = dtefield%qijb_kk(1,klmn,gs_hamk%atindx1(iatom),idir)
             qijbkk(2*klmn,  iatom,ispinor) = dtefield%qijb_kk(2,klmn,gs_hamk%atindx1(iatom),idir)
             if (ifor > 1) qijbkk(2*klmn,iatom,ispinor) = -qijbkk(2*klmn,iatom,ispinor)
           end do
         end do ! end loop over lmn2_size
       end do ! end loop over natom

       choice = 1
       signs = 2
       paw_opt = 1
       cpopt = 2 ! use cprj_kb in memory
       nspinortot=min(2,nspinor*(1+mpi_enreg%paral_spin))
       do i_paw_band = 1, nbo
         call cprj_get(gs_hamk%atindx,cprj_band_srt,cprj_kb,natom,i_paw_band,0,ikpt,1,isppol,&
&         nbo,1,mpi_enreg,natom,1,nbo,nspinor,nsppol,0)
         call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_band_srt,dimenlc1,natom,&
&         dimffnl,dimffnl,qijbkk,qijbkk,ffnl,ffnl,gs_hamk%gmet,&
&         gs_hamk%gprimd,idir,gs_hamk%indlmn,gs_hamk%istwf_k,&
&         kg_k,kg_k,kpg,kpg,gs_hamk%kpoint,gs_hamk%kpoint,zero,lmnmax,&
&         matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,&
&         gs_hamk%ngfft,nkpg,nkpg,nloalg,0,npw,npw,nspinor,nspinortot,ntypat,&
&         0,paw_opt,gs_hamk%phkxred,gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,&
&         signs,sij_dum,svectout_dum,0,gs_hamk%ucvol,gs_hamk%useylm,direc,grad_berry_ev,&
&         use_gpu_cuda=gs_hamk%use_gpu_cuda)

!        Add i*fac*smat_inv(i_paw_band,iband)*grad_berry_ev to the gradient
         do ipw = 1, npw*nspinor

           grad_berry(1,ipw) = grad_berry(1,ipw) - &
&           fac*(smat_inv(2,i_paw_band,iband)*grad_berry_ev(1,ipw) + &
&           smat_inv(1,i_paw_band,iband)*grad_berry_ev(2,ipw))

           grad_berry(2,ipw) = grad_berry(2,ipw) + &
&           fac*(smat_inv(1,i_paw_band,iband)*grad_berry_ev(1,ipw) - &
&           smat_inv(2,i_paw_band,iband)*grad_berry_ev(2,ipw))

         end do
       end do
     end if ! end if PAW

!    Add i*fac*cg1_k to the gradient
     do ipw = 1, npw*nspinor
       grad_berry(1,ipw) = grad_berry(1,ipw) - fac*cg1_k(2,ipw)
       grad_berry(2,ipw) = grad_berry(2,ipw) + fac*cg1_k(1,ipw)
     end do
     fac = -1._dp*fac
     dtefield%sflag(:,ikpt+(isppol-1)*nkpt,ifor,idir) = sflag_k(:)
     dtefield%sflag(iband,ikpt+(isppol-1)*nkpt,ifor,idir) = 0
     dtefield%smat(:,:,:,ikpt+(isppol-1)*nkpt,ifor,idir) = smat_k(:,:,:)
   end do  ! ifor

   if (gs_hamk%usepaw == 1) then
!    call nonlop to apply on-site dipole <EV> part to direc
!    note that rij is sorted by input atom order, but nonlop wants it sorted by type
     do iatom = 1, natom
       itypat = gs_hamk%typat(gs_hamk%atindx1(iatom))
       do klmn = 1, dtefield%lmn2_size(itypat)
!        note: D_ij-like terms have 4 spinor components: 11, 22, 12, and 21. Here the enl_rij is diagonal
!        in spin space so only the first two are nonzero and they are equal
         do ispinor = 1, nspinor
           enl_rij(klmn,iatom,ispinor) = dtefield%rij(klmn,itypat,idir)
         end do
       end do ! end loop over lmn2_size
     end do ! end loop over natom
     cpopt = -1 ! compute cprj inside nonlop because we do not have them for direc
     call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,dimenlr1,natom,&
&     dimffnl,dimffnl,enl_rij,enl_rij,ffnl,ffnl,gs_hamk%gmet,&
&     gs_hamk%gprimd,idir,gs_hamk%indlmn,gs_hamk%istwf_k,&
&     kg_k,kg_k,kpg,kpg,gs_hamk%kpoint,gs_hamk%kpoint,zero,lmnmax,&
&     matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,&
&     gs_hamk%ngfft,nkpg,nkpg,nloalg,0,npw,npw,nspinor,nspinortot,ntypat,&
&     0,paw_opt,gs_hamk%phkxred,gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,&
&     signs,sij_dum,svectout_dum,0,gs_hamk%ucvol,gs_hamk%useylm,direc,grad_berry_ev,&
&     use_gpu_cuda=gs_hamk%use_gpu_cuda)

     grad_berry(:,:) = grad_berry(:,:) - dtefield%efield_dot(idir)*grad_berry_ev(:,:)/two_pi

   end if

 end do    ! idir

!deallocations
 if(gs_hamk%usepaw /= 0) then
   ABI_DEALLOCATE(grad_berry_ev)
   ABI_DEALLOCATE(qijbkk)
   ABI_DEALLOCATE(enl_rij)
   ABI_DEALLOCATE(smat_k_paw)
   call cprj_free(cprj_kb)
   call cprj_free(cprj_band_srt)
   ABI_DEALLOCATE(cprj_kb)
   ABI_DEALLOCATE(cprj_band_srt)
   if (nkpt /= dtefield%fnkpt) then
     call cprj_free(cprj_fkn)
     call cprj_free(cprj_ikn)
     ABI_DEALLOCATE(cprj_fkn)
     ABI_DEALLOCATE(cprj_ikn)
   end if
   ABI_DEALLOCATE(kpg)
 end if

!DBG_EXIT("COLL")

end subroutine make_grad_berry
!!***
