!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp2in
!! NAME
!! psp2in
!!
!! FUNCTION
!! Initialize pspcod=2 pseudopotentials (GTH format):
!! continue to read the file, then compute the corresponding
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  ipsp=id in the array of the pseudo-potential.
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  zion=nominal valence of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+Zv/r) dr]$ (hartree)
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpsang)=number of projection functions for each angular momentum
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  dvlspl(mqgrid_vl,2)=dVloc(r)/dr and second derivatives from spline fit (only
!!                      allocated if vlspl_recipSpace is false.
!!
!! SIDE EFFECTS
!!  Input/output
!!  lmax : at input =value of lmax mentioned at the second line of the psp file
!!    at output= 1
!!  psps <type(pseudopotential_type)>=at output, values depending on the read
!!                                    pseudo are set.
!!   | lmnmax(IN)=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!   |           =if useylm=0, max number of (l,n)   comp. over all type of psps
!!   | lnmax(IN)=max. number of (l,n) components over all type of psps
!!   |           angular momentum of nonlocal pseudopotential
!!   | mpsang(IN)= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | mqgrid_ff(IN)=dimension of q (or G) grid for nl form factors (array ffspl)
!!   | mqgrid_vl(IN)=dimension of q (or G) grid or r grid (if vlspl_recipSpace = .false.)
!!   | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!   | qgrid_vl(mqgrid_vl)(IN)=values of q on grid from 0 to qmax (bohr^-1) for Vloc
!!   |                         if vlspl_recipSpace is .true. else values of r on grid from
!!   |                         0 to 2pi / qmax * mqgrid_ff (bohr).
!!   | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!   |            1=using Ylm, 0=using Legendre polynomials
!!   | vlspl_recipSpace(IN)=.true. if pseudo are expressed in reciprocal space.
!!   | gth_params(OUT)=store GTH coefficients and parameters.
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!      eleconf,leave_new,psp2lo,psp2nl,spline,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp2in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps,vlspl,dvlspl,zion)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_splines
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only: eleconf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp2in'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_65_psp, except_this_one => psp2in
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ipsp,lmax
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: indlmn(6,psps%lmnmax),nproj(psps%mpsang)
 real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2),ekb(psps%lnmax)
 real(dp),intent(out) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
 real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)

!Local variables-------------------------------
#if defined HAVE_DFT_BIGDFT
 integer :: mxpl,mxchg
 character(len=2) :: symbol
 real(dp) :: rcov, rprb, ehomo, radfine, amu
 real(dp) :: neleconf(6,0:3)
#endif
!scalars
 integer :: ii,iln,index,ios,ipsang,isc,kk,ll,mm
 real(dp) :: cc1,cc2,cc3,cc4,h1p,h1s,h2s,maxrad,rad_core,rad_cov,rad_long,rloc,rrp,rrs
 real(dp) :: yp1,ypn
 character(len=500) :: message
!arrays
 real(dp),allocatable :: work_space(:),work_spl(:)
 real(dp),allocatable :: dvloc(:)

! ***************************************************************************

!Set various terms to 0 in case not defined below
!GTH values
 rloc=0.d0
 cc1=0.d0
 cc2=0.d0
 cc3=0.d0
 cc4=0.d0
 rrs=0.d0
 h1s=0.d0
 h2s=0.d0
 rrp=0.d0
 h1p=0.d0
 nproj(1:psps%mpsang)=0
 rad_cov  = -1.d0
 rad_core = -1.d0
 rad_long = -1.d0

!DEBUG
!write(message,'(a)') ' psp2in : enter '
!call wrtout(std_out,  message,'COLL')
!ENDDEBUG

!Read and write different lines of the pseudopotential file
 read (tmp_unit,*) rloc,cc1,cc2,cc3,cc4
 write(message, '(a,f12.7)' ) ' rloc=',rloc
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
& '  cc1=',cc1,'; cc2=',cc2,'; cc3=',cc3,'; cc4=',cc4
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 read (tmp_unit,*) rrs,h1s,h2s
 write(message, '(a,f12.7,a,f12.7,a,f12.7)' )&
& '  rrs=',rrs,'; h1s=',h1s,'; h2s=',h2s
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 read (tmp_unit,*) rrp,h1p
 write(message, '(a,f12.7,a,f12.7)' )&
& '  rrp=',rrp,'; h1p=',h1p
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Store the coefficients.
 psps%gth_params%set(ipsp)          = .true.
 psps%gth_params%psppar(0, :, ipsp) = (/ rloc, cc1, cc2, cc3, cc4, 0.d0, 0.d0 /)
 psps%gth_params%psppar(1, :, ipsp) = (/ rrs,  h1s, h2s, 0.d0, 0.d0, 0.d0, 0.d0 /)
 psps%gth_params%psppar(2, :, ipsp) = (/ rrp,  h1p, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)

 if (dtset%usewvl == 1) then
#if defined HAVE_DFT_BIGDFT
   call eleconf(int(psps%znuclpsp(ipsp)), int(zion), symbol, rcov, rprb, &
&   ehomo, neleconf, isc, mxpl, mxchg, amu)
#endif
!  We try to read the values from the pseudo.
   read (tmp_unit, *, iostat = ios) rad_long, rad_core, rad_cov
   if (ios /= 0 .or. rad_long == zero .or. rad_core == zero) then
#if defined HAVE_DFT_BIGDFT
!    assigning the radii by calculating physical parameters
     rad_long = one / sqrt(abs(two * ehomo))
     radfine  = 100.d0
     do ii = 0, 4, 1
       if (psps%gth_params%psppar(ii, 0, ipsp) /= zero) then
         radfine = min(radfine, psps%gth_params%psppar(ii, 0, ipsp))
       end if
     end do
     rad_core = radfine
     rad_cov  = rcov
#endif

     write(message, '(a,a,a,a,a,a,a)' ) '-', ch10,&
&     '- psp2in : COMMENT -',ch10,&
&     "-  the pseudo-potential does not include geometric informations,",ch10,&
&     '-  values have been computed.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   write(message, '(a,f12.7,a,f12.7,a,f12.7)' )&
&   '  radii_cf(1)=', rad_long,'; radii_cf(2)=', rad_core, "; rad_cov=", rad_cov
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 psps%gth_params%radii_cf(ipsp, :)  = (/ rad_long, rad_core, rad_core /)
 psps%gth_params%radii_cov(ipsp)    = rad_cov
 if (rad_core > 0.d0 .and. rad_long > 0.d0) then
   psps%gth_params%hasGeometry(ipsp) = .true.
 else
   psps%gth_params%hasGeometry(ipsp) = .false.
 end if
 psps%gth_params%semicore(ipsp)     = isc
!correct the coarse and the fine radius for projectors
 maxrad = zero
 do ii = 0, 2, 1
!  the maximum radii is useful only for projectors
   if (ii == 1) maxrad = zero
   if (psps%gth_params%psppar(ii, 0, ipsp) /= zero) then
     maxrad = max(maxrad, psps%gth_params%psppar(ii, 0, ipsp))
   end if
 end do
 if (maxrad == zero) then
   psps%gth_params%radii_cf(ipsp,3)=zero
 else
   psps%gth_params%radii_cf(ipsp,3)=max( &
&   min(dtset%wvl_crmult*psps%gth_params%radii_cf(ipsp,1), &
&   15._dp*maxrad)/dtset%wvl_frmult, &
   psps%gth_params%radii_cf(ipsp,2))
 end if

 if (abs(h1s)>1.d-08) nproj(1)=1
 if (abs(h2s)>1.d-08) nproj(1)=2

 if (abs(h1p)>1.d-08) then
   if(psps%mpsang<2)then
     write(message, '(a,a,a,a,es12.4,a,a,a,i2,a)' ) ch10,&
&     ' psp2in : BUG -',ch10,&
&     '  With non-zero h1p (=',h1p,&
&     '), mpsang should be at least 2,',ch10,&
&     '  while mpsang=',psps%mpsang,'.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
   nproj(2)=1
   if (lmax<1) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a,a,a)' ) ch10,&
&     ' psp2in : ERROR -',ch10,&
&     '  Input lmax=',lmax,' disagree with input h1p=',h1p,'.',&
&     '  Your pseudopotential is incoherent.',ch10,&
&     '  Action : correct your pseudopotential file.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
 end if

!Initialize array indlmn array giving l,m,n,lm,ln,s for i=lmn
 index=0;iln=0;indlmn(:,:)=0
 do ipsang=1,lmax+1
   if(nproj(ipsang)>0)then
     ll=ipsang-1
     do kk=1,nproj(ipsang)
       iln=iln+1
       do mm=1,2*ll*psps%useylm+1
         index=index+1
         indlmn(1,index)=ll
         indlmn(2,index)=mm-ll*psps%useylm-1
         indlmn(3,index)=kk
         indlmn(4,index)=ll*ll+(1-psps%useylm)*ll+mm
         indlmn(5,index)=iln
         indlmn(6,index)=1
       end do
     end do
   end if
 end do

!First, the local potential --
!compute q^2V(q) or V(r)
!MJV NOTE: psp2lo should never be called with dvspl unallocated, which
!is possible unless .not.psps%vlspl_recipSpace
 ABI_ALLOCATE(dvloc,(psps%mqgrid_vl))
 call psp2lo(cc1,cc2,cc3,cc4,dvloc,epsatm,psps%mqgrid_vl,psps%qgrid_vl,&
& vlspl(:,1),rloc,psps%vlspl_recipSpace,yp1,ypn,zion)

!DEBUG
!write(std_out,*)' psp2in : after psp2lo '
!stop
!ENDDEBUG

!Fit spline to (q^2)V(q) or V(r) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_space,(psps%mqgrid_vl))
 ABI_ALLOCATE(work_spl,(psps%mqgrid_vl))
 call spline (psps%qgrid_vl,vlspl(:,1),psps%mqgrid_vl,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 if (.not.psps%vlspl_recipSpace) then
   dvlspl(:,1) = dvloc
   call spline (psps%qgrid_vl,dvlspl(:,1),psps%mqgrid_vl,yp1,ypn,work_spl)
   dvlspl(:,2)=work_spl(:)
 end if
 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl)
 ABI_DEALLOCATE(dvloc)


!Second, compute KB energies and form factors and fit splines
 ekb(:)=0.0d0
!First check if any nonlocal projectors are being used
 if (maxval(nproj(1:lmax+1))>0) then
   call psp2nl(ekb,ffspl,h1p,h1s,h2s,psps%lnmax,psps%mqgrid_ff,psps%qgrid_ff,rrp,rrs)
 end if
 
end subroutine psp2in
!!***
