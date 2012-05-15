!{\src2tex{textfont=tt}}
!!****f* ABINIT/berry_linemin.F90
!! NAME
!! berry_linemin
!!
!! FUNCTION
!! numerical line minimization for electric field case
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
!!  chc = <C|H_0|C> where |C> is the wavefunction of the current band
!!  cprj_k(natom,nband_k*usepaw)=cprj at this k point
!!  detovc(2,2,3) = overlap determinant for conjugate gradient direction
!!  dhc = Re[<D|H_0|C>]
!!  dhd = <D|H_0|D>
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dimlmn(natom)=lmn_size for each atom in input order
!!  dimlmn_srt(natom)=lmn_size for each atom sorted by type
!!  direc(2,npw*nspinor)=gradient vector
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gen_eigenpb :: logical flag concerning generalized (true) eigen problem 
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  iband=index of band currently being treated
!!  ikpt=number of the k-point currently being treated
!!  iline = index of the current line minimization
!!  isppol=spin polarization currently treated
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  natom=number of atoms in cell.
!!  matblk=dimension of the array ph3d
!!  mband =maximum number of bands
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
!!  ntypat=number of types of atoms in cell.
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwnsfacq(2,mkgq) = phase factors for the nearest neighbours of the
!!                     current k-point (electric field, MPI //)
!!  spaceComm_distrb=mpi data
!!  xnorm=normalizing factor from cgwf.F90
!!
!! OUTPUT
!! bcut(ifor,idir) = branch cut of the ellipse associated with (ifor,idir)
!! costh = cos(thetam)
!! hel(ifor,idir) = helicity of the ellipse associated with (ifor,idir)
!! phase_end = total change in Zak phase, must be equal to
!!             dphase_aux1 + n*two_pi
!! sinth = sin(thetam)
!! thetam = optimal angle theta in line_minimization
!!
!! SIDE EFFECTS
!! dphase_aux1 = change in Zak phase accumulated during the loop over iline
!!               (can be used for debugging in cgwf.f)
!! phase_init = initial Zak phase (before doing the first line minimization)
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,cprj_get,cprj_put,getcprj,linemin
!!      smatrix,smatrix_k_paw,sym_cprj_kn
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine berry_linemin(bcut,chc,cg,cgq,costh,cprj_k,detovc,dphase_aux1,dhc,dhd,&
&                        dimffnl,dimlmn,dimlmn_srt,direc,dtefield,ffnl,gen_eigenpb,&
&                        gs_hamk,hel,iband,ikpt,iline,isppol,kg_k,lmnmax,matblk,&
&                        mband,mcg,mcgq,mgfft,mkgq,mpi_enreg,mpw,&
&                        natom,nkpt,nloalg,npw,npwarr,nspinor,ntypat,pwind,ph3d,&
&                        phase_init,phase_end,pwind_alloc,pwnsfac,pwnsfacq,sinth,&
&                        spaceComm_distrb,thetam,xnorm)

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
#define ABI_FUNC 'berry_linemin'
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => berry_linemin
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!scalars
 integer,intent(in) :: dimffnl,iband,ikpt,iline,isppol,lmnmax,matblk,mband,mcg,mcgq,mgfft
 integer,intent(in) :: mkgq,mpw,natom,nkpt,npw,spaceComm_distrb,nspinor,ntypat,pwind_alloc
 real(dp),intent(in) :: chc,dhc,dhd,xnorm
 real(dp),intent(out) :: costh,sinth,thetam
 logical,intent(in) :: gen_eigenpb
 type(gs_hamiltonian_type),intent(in) :: gs_hamk
 type(efield_type),intent(in) :: dtefield
 type(MPI_type),intent(inout) :: mpi_enreg

!arrays
 integer,intent(in) :: dimlmn(natom),dimlmn_srt(natom),kg_k(3,npw),nloalg(5)
 integer,intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 integer,intent(out) :: hel(2,3)
 real(dp),intent(in) :: cg(2,mcg),cgq(2,mcgq),detovc(2,2,3),ffnl(npw,dimffnl,lmnmax,ntypat)
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)
 real(dp),intent(inout) :: direc(2,npw*nspinor),dphase_aux1(3),phase_init(3)
 real(dp),intent(inout) :: ph3d(2,npw,matblk)
 real(dp),intent(out) :: bcut(2,3),phase_end(3)
 type(cprj_type),intent(in) :: cprj_k(natom,dtefield%nband_occ*gs_hamk%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ddkflag,icg1,icp2,idum1,idir,ifor,ikgf,ikptf,ikpt2,ikpt2f
 integer :: itrs,job,mcg1_k,mcg_q,nbo,nkpg,npw_k2,shiftbd
!arrays
 integer :: pwind_k(npw),sflag_k(dtefield%nband_occ)
 real(dp) :: cg1_k(2,npw*nspinor),detovd(2,2,3),direc_tmp(2,npw*nspinor)
 real(dp) :: dtm_k(2),pwnsfac_k(4,mpw)
 real(dp) :: smat_k(2,dtefield%nband_occ,dtefield%nband_occ)
 real(dp) :: smat_inv(2,dtefield%nband_occ,dtefield%nband_occ)
 real(dp),allocatable :: cgq_k(:,:),kpg(:,:)
 real(dp),allocatable :: smat_k_paw(:,:,:)
 type(cprj_type),allocatable :: cprj_direc(:,:),cprj_kb(:,:),cprj_band_srt(:,:)
 type(cprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)


! *********************************************************************

 DBG_ENTER("COLL")

 nbo = dtefield%nband_occ
 ikptf = dtefield%i2fbz(ikpt)
 ikgf = dtefield%fkgindex(ikptf)  ! this is the shift for pwind
 mcg_q = mpw*mband*nspinor
 mcg1_k = npw*nspinor

 if (gen_eigenpb) then
   ABI_ALLOCATE(smat_k_paw,(2,nbo,nbo))
   smat_k_paw = zero
   ABI_ALLOCATE(cprj_kb,(natom,nbo))
   call cprj_alloc(cprj_kb,0,dimlmn)
   ABI_ALLOCATE(cprj_band_srt,(natom,1))
   call cprj_alloc(cprj_band_srt,0,dimlmn_srt)
   ABI_ALLOCATE(cprj_direc,(natom,nbo))
   call cprj_alloc(cprj_direc,0,dimlmn)
   if (nkpt /= dtefield%fnkpt) then
     ABI_ALLOCATE(cprj_fkn,(natom,nbo))
     ABI_ALLOCATE(cprj_ikn,(natom,nbo))
     call cprj_alloc(cprj_fkn,0,dimlmn)
     call cprj_alloc(cprj_ikn,0,dimlmn)
   end if
   nkpg = 0
   ABI_ALLOCATE(kpg,(npw,nkpg))
 end if

!In case the eletric field is on, the line minimization has to be done numerically

!Compute determinant of the overlap matrix where in the band-th line
!the wavefunction is replaced by the search direction
 job = 10 ; shiftbd = 0
 do idir = 1, 3
   if (abs(dtefield%efield_dot(idir)) < tol12) cycle
   do ifor = 1, 2
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
       idum1=dtefield%fkgindex(ikpt2f)
       pwnsfac_k(3:4,1:npw_k2) = pwnsfac(1:2,idum1+1:idum1+npw_k2)
     end if
     icg1 = 0 ; ddkflag = 0
     if (gen_eigenpb) then
       direc_tmp(1:2,1:npw*nspinor) = direc(1:2,1:npw*nspinor)*xnorm
!      need cprj corresponding to direc_tmp in order to make smat_k_paw properly
       call getcprj(1,0,direc_tmp,cprj_band_srt,gs_hamk%dimekb1,gs_hamk%dimekb2,&
&       dimffnl,gs_hamk%ekb,ffnl,0,gs_hamk%indlmn,gs_hamk%istwf_k,kg_k,kpg,gs_hamk%kpoint,&
&       lmnmax,matblk,mgfft,mpi_enreg,natom,gs_hamk%nattyp,gs_hamk%ngfft,0,&
&       nloalg,npw,nspinor,ntypat,gs_hamk%phkxred,gs_hamk%ph1d,ph3d,gs_hamk%ucvol,&
&       gs_hamk%usepaw,gs_hamk%useylm)
       icp2=nbo*(ikpt2-1)
       call cprj_copy(cprj_k,cprj_direc)
       call cprj_put(gs_hamk%atindx,cprj_band_srt,cprj_direc,natom,iband,0,ikpt,1,1,&
&       nbo,1,mpi_enreg,natom,1,nbo,&
&       dimlmn,1,1,spaceComm_distrb,0)
       call cprj_get(gs_hamk%atindx1,cprj_kb,dtefield%cprj,natom,1,icp2,ikpt,0,1,&
&       nbo,dtefield%fnkpt,mpi_enreg,natom,&
&       nbo,nbo,1,1,0)
       if (ikpt2 /= ikpt2f) then ! construct cprj_kb by symmetry
         call cprj_copy(cprj_kb,cprj_ikn)
         call sym_cprj_kn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,gs_hamk%indlmn,&
&         dtefield%indkk_f2ibz(ikpt2f,2),dtefield%indkk_f2ibz(ikpt2f,6),&
&         dtefield%fkptns(:,dtefield%i2fbz(ikpt2)),&
&         dtefield%lmax,dtefield%lmnmax,natom,nbo,nspinor,&
&         dtefield%nsym,ntypat,gs_hamk%typat,dtefield%zarot)
         call cprj_copy(cprj_fkn,cprj_kb)
       end if
       call smatrix_k_paw(cprj_direc,cprj_kb,dtefield,idir,ifor,natom,smat_k_paw,gs_hamk%typat)
       call smatrix(direc_tmp,cgq_k,cg1_k,ddkflag,dtm_k,icg1,icg1,&
&       itrs,job,iband,npw,mcg_q,mpw,iband,mpw,nbo,&
&       npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&       shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
     else ! not a generalized eigenproblem
       call smatrix(direc,cgq_k,cg1_k,ddkflag,dtm_k,icg1,icg1,&
&       itrs,job,iband,npw,mcg_q,mpw,iband,mpw,nbo,&
&       npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&       shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
     end if
     ABI_DEALLOCATE(cgq_k)
     detovd(:,ifor,idir) = dtm_k(:) ! Store the determinant of the overlap
!    matrix (required to compute theta_min)
   end do  ! ifor
 end do    ! idir

 call linemin(bcut,chc,costh,detovc,detovd,dhc,dhd,&
& dphase_aux1,dtefield%efield_dot,iline,&
& dtefield%fnkpt,dtefield%nstr,hel,phase_end,&
& phase_init,dtefield%sdeg,sinth,thetam)

!deallocations
 if(gen_eigenpb) then
   ABI_DEALLOCATE(smat_k_paw)
   call cprj_free(cprj_kb)
   ABI_DEALLOCATE(cprj_kb)
   call cprj_free(cprj_band_srt)
   ABI_DEALLOCATE(cprj_band_srt)
   call cprj_free(cprj_direc)
   ABI_DEALLOCATE(cprj_direc)
   if (nkpt /= dtefield%fnkpt) then
     call cprj_free(cprj_fkn)
     call cprj_free(cprj_ikn)
     ABI_DEALLOCATE(cprj_fkn)
     ABI_DEALLOCATE(cprj_ikn)
   end if
   ABI_DEALLOCATE(kpg)
 end if

 DBG_EXIT("COLL")

end subroutine berry_linemin
!!***
