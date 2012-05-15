!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp3in
!! NAME
!! psp3in
!!
!! FUNCTION
!! Initialize pspcod=3 pseudopotentials (HGH psps PRB58,3641(1998)):
!! continue to read the file, then compute the corresponding
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, FD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  pspso=spin-orbit characteristics, govern the content of ffspl and ekb
!!   if =0 : this input requires NO spin-orbit characteristics of the psp
!!   if =2 : this input requires HGH characteristics of the psp
!!   if =3 : this input requires HFN characteristics of the psp
!!  ipsp=id in the array of the pseudo-potential.
!!  zion=nominal valence of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!             if any, spin-orbit components begin at l=mpsang+1
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$ (hartree)
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpssoang)=number of projection functions for each angular momentum
!!  vlspl(mqgrid_ff,2)=q^2 Vloc(q) and second derivatives from spline fit
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
!!   | mpssoang(IN)= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!   | mqgrid_ff(IN)=dimension of q (or G) grid for arrays.
!!   | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!   | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!   |            1=using Ylm, 0=using Legendre polynomials
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!      eleconf,psp2lo,psp3nl,spline,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp3in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax, &
     & nproj, psps, pspso, vlspl, zion)

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
#define ABI_FUNC 'psp3in'
 use interfaces_14_hidewrite
 use interfaces_65_psp, except_this_one => psp3in
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ipsp,pspso
 integer,intent(inout) :: lmax
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: indlmn(6,psps%lmnmax),nproj(psps%mpssoang)
 real(dp),intent(out) :: ekb(psps%lnmax),ffspl(psps%mqgrid_ff,2,psps%lnmax)
 real(dp),intent(out) :: vlspl(psps%mqgrid_ff,2)

!Local variables-------------------------------
#if defined HAVE_DFT_BIGDFT
 integer :: mxpl,mxchg
 character(len=2) :: symbol
 real(dp) :: rcov, rprb, ehomo, radfine, amu
 real(dp) :: neleconf(6,0:3)
#endif
!scalars
 integer :: ii,iln,iln0,index,ios,ipsang,jj,kk,ll,mm,mproj,nn,nso, isc
 real(dp) :: cc1,cc2,cc3,cc4,h11d,h11f,h11p,h11s,h22d,h22p,h22s,h33d,h33p,h33s
 real(dp) :: k11d,k11f,k11p,k22d,k22p,k33d,k33p,maxrad,rad_core,rad_cov,rad_long,rloc
 real(dp) :: rrd,rrf,rrp,rrs,yp1,ypn
 character(len=500) :: message
!arrays
 real(dp),allocatable :: dvlspl(:),ekb_so(:,:),ekb_sr(:,:),ffspl_so(:,:,:,:)
 real(dp),allocatable :: ffspl_sr(:,:,:,:),work_space(:),work_spl(:)

! ***************************************************************************

!Set various terms to 0 in case not defined below
!HGH values
 rloc=zero ; rrs=zero  ; h11p=zero ; k33p=zero ; k11d=zero;
 cc1=zero  ; h11s=zero ; h22p=zero ; rrd=zero  ; k22d=zero;
 cc2=zero  ; h22s=zero ; h33p=zero ; h11d=zero ; k33d=zero;
 cc3=zero  ; h33s=zero ; k11p=zero ; h22d=zero ; h11f=zero;
 cc4=zero  ; rrp=zero  ; k22p=zero ; h33d=zero ; k11f=zero;
 rrf=zero
 nproj(1:psps%mpssoang)=0
 rad_cov  = -1.d0
 rad_core = -1.d0
 rad_long = -1.d0
 isc      = 0

!DEBUG
!write(std_out,*) ' psp3in : enter '
!ENDDEBUG

!Read and write different lines of the pseudopotential file

 read (tmp_unit,*) rloc,cc1,cc2,cc3,cc4
 write(message, '(a,f12.7)' ) ' rloc=',rloc
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
& ' cc1 =',cc1,'; cc2 =',cc2,'; cc3 =',cc3,'; cc4 =',cc4
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!For the time being, the s state line must be present and is read,
!even for local pseudopotentials (zero must appear)
 read (tmp_unit,*) rrs,h11s,h22s,h33s
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
& ' rrs =',rrs,'; h11s=',h11s,'; h22s=',h22s,'; h33s=',h33s
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 if (lmax > 0) then

   read (tmp_unit,*) rrp,h11p,h22p,h33p
   write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
&   ' rrp =',rrp,'; h11p=',h11p,'; h22p=',h22p,'; h33p=',h33p
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*) k11p,k22p,k33p
   write(message, '(a,f12.7,a,f12.7,a,f12.7)' )&
&   '                    k11p=',k11p,'; k22p=',k22p,'; k33p=',k33p
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 end if

 if (lmax > 1) then

   read (tmp_unit,*) rrd,h11d,h22d,h33d
   write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
&   ' rrd =',rrd,'; h11d=',h11d,'; h22d=',h22d,'; h33d=',h33d
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*) k11d,k22d,k33d
   write(message, '(a,f12.7,a,f12.7,a,f12.7)' )&
&   '                    k11d=',k11d,'; k22d=',k22d,'; k33d=',k33d
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 end if

 if (lmax > 2) then

   read (tmp_unit,*) rrf,h11f
   write(message, '(a,f12.7,a,f12.7)' )' rrf =',rrf,'; h11f=',h11f
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*) k11f
   write(message, '(a,f12.7)' )'                    k11f=',k11f
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 end if

 if (abs(h11s)>1.d-08) nproj(1)=1
 if (abs(h22s)>1.d-08) nproj(1)=2
 if (abs(h33s)>1.d-08) nproj(1)=3

 if (abs(h11p)>1.d-08) then
   nproj(2)=1
   if (lmax<1) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,'  does not agree with input h11p=',h11p,ch10,&
&     '  setting lmax to 1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=1
   end if
 end if

 if (abs(h22p)>1.d-08) then
   nproj(2)=2
   if (lmax<1) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,' does not agree with input h22p=',h22p,ch10,&
&     '  setting lmax to 1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=1
   end if
 end if


 if (abs(h33p)>1.d-08) then
   nproj(2)=3
   if (lmax<1) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,' does not agree with input h33p=',h33p,ch10,&
&     '  setting lmax to 1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=1
   end if
 end if

 if (abs(h11d)>1.d-08) then
   nproj(3)=1
   if (lmax<2) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,'  does not agree with input h11d=',h11d,ch10,&
&     '  setting lmax to 2'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=2
   end if
 end if

 if (abs(h22d)>1.d-08) then
   nproj(3)=2
   if (lmax<2) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,'  does not agree with input h22d=',h22d,ch10,&
&     '  setting lmax to 2'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=2
   end if
 end if

 if (abs(h33d)>1.d-08) then
   nproj(3)=3
   if (lmax<2) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,' does not agree with input h33d=',h33d,ch10,&
&     '  setting lmax to 2'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=2
   end if
 end if

 if (abs(h11f)>1.d-08) then
   nproj(4)=1
   if (lmax<3) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,' does not agree with input h11f=',h11f,ch10,&
&     '  setting lmax to 3'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=3
   end if
 end if

 if(pspso/=0) then

   if (abs(k11p)>1.d-08) then
     nproj(psps%mpsang+1)=1
     if (lmax<1) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,'  does not agree with input k11p=',k11p,ch10,&
&       '  setting lmax to 1'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=1
     end if
   end if

   if (abs(k22p)>1.d-08) then
     nproj(psps%mpsang+1)=2
     if (lmax<1) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,' does not agree with input k22p=',k22p,ch10,&
&       '  setting lmax to 1'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=1
     end if
   end if


   if (abs(k33p)>1.d-08) then
     nproj(psps%mpsang+1)=3
     if (lmax<1) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,' does not agree with input k33p=',k33p,ch10,&
&       '  setting lmax to 1'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=1
     end if
   end if

   if (abs(k11d)>1.d-08) then
     nproj(psps%mpsang+2)=1
     if (lmax<2) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,'  does not agree with input k11d=',k11d,ch10,&
&       '  setting lmax to 2'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=2
     end if
   end if

   if (abs(k22d)>1.d-08) then
     nproj(psps%mpsang+2)=2
     if (lmax<2) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,'  does not agree with input k22d=',k22d,ch10,&
&       '  setting lmax to 2'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=2
     end if
   end if

   if (abs(k33d)>1.d-08) then
     nproj(psps%mpsang+2)=3
     if (lmax<2) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,' does not agree with input k33d=',k33d,ch10,&
&       '  setting lmax to 2'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=2
     end if
   end if

   if (abs(k11f)>1.d-08) then
     nproj(psps%mpsang+3)=1
     if (lmax<3) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,' does not agree with input k11f=',k11f,ch10,&
&       '  setting lmax to 3'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=3
     end if
   end if

 end if

!Store the coefficients.
 psps%gth_params%set(ipsp)          = .true.
 psps%gth_params%psppar(0, :, ipsp) = (/ rloc, cc1, cc2, cc3, cc4, zero, zero /)
 psps%gth_params%psppar(1, :, ipsp) = (/ rrs,  h11s, h22s, h33s, zero, zero, zero /)
 psps%gth_params%psppar(2, :, ipsp) = (/ rrp,  h11p, h22p, h33p, zero, zero, zero /)
 psps%gth_params%psppar(3, :, ipsp) = (/ rrd,  h11d, h22d, h33d, zero, zero, zero /)
 psps%gth_params%psppar(4, :, ipsp) = (/ rrf,  h11f, zero, zero, zero, zero, zero /)

!Store the k coefficients
 psps%gth_params%psp_k_par(1, :, ipsp) = (/ zero, zero, zero /)
 psps%gth_params%psp_k_par(2, :, ipsp) = (/ k11p, k22p, k33p /)
 psps%gth_params%psp_k_par(3, :, ipsp) = (/ k11d, k22d, k33d /)
 psps%gth_params%psp_k_par(4, :, ipsp) = (/ k11f, zero, zero /)

 if (dtset%usewvl == 1) then
#if defined HAVE_DFT_BIGDFT
   call eleconf(int(psps%znuclpsp(ipsp)), int(zion), symbol, rcov, rprb, &
&   ehomo, neleconf, isc, mxpl, mxchg, amu)
#endif
!  We try to read the values from the pseudo.
   read (tmp_unit, *, iostat = ios) rad_long, rad_core, rad_cov
   if (ios /= 0 .or. rad_long == zero .or. rad_core == zero .or. rad_cov == zero) then
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
&     '- psp3in : COMMENT -',ch10,&
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

 psps%gth_params%radii_cf(ipsp, 1:2)  = (/ rad_long, rad_core /)
 psps%gth_params%radii_cov(ipsp)    = rad_cov
 if (rad_core > tol12 .and. rad_long > tol12) then
   psps%gth_params%hasGeometry(ipsp) = .true.
 else
   psps%gth_params%hasGeometry(ipsp) = .false.
 end if
 psps%gth_params%semicore(ipsp)     = isc
!correct the coarse and the fine radius for projectors
 maxrad = zero
 do ii = 0, 4, 1
!  the maximum radii is useful only for projectors
   if (ii == 1) maxrad = zero
   if (psps%gth_params%psppar(ii, 0, ipsp) /= zero) then
     maxrad = max(maxrad, psps%gth_params%psppar(ii, 0, ipsp))
   end if
 end do
 psps%gth_params%radii_cf(ipsp,3)   = min(dtset%wvl_crmult * &
& psps%gth_params%radii_cf(ipsp, 1), real(15.0, dp) * maxrad) / dtset%wvl_frmult

!Initialize array indlmn array giving l,m,n,ln,lm,s for i=lmn
 nso=1;if(pspso/=0) nso=2
 index=0;iln=0;indlmn(:,:)=0
 do nn=1,nso
   do ipsang=1+(nn-1)*(lmax+1),nn*lmax+1
     if (nproj(ipsang)>0) then
       ll=ipsang-(nn-1)*lmax-1
       do kk=1,nproj(ipsang)
         iln=iln+1
         do mm=1,2*ll*psps%useylm+1
           index=index+1
           indlmn(1,index)=ll
           indlmn(2,index)=mm-ll*psps%useylm-1
           indlmn(3,index)=kk
           indlmn(4,index)=ll*ll+(1-psps%useylm)*ll+mm
           indlmn(5,index)=iln
           indlmn(6,index)=nn
         end do
       end do
     end if
   end do
 end do

 ABI_ALLOCATE(dvlspl,(psps%mqgrid_ff))
!First, the local potential --  compute on q grid and fit spline
 call psp2lo(cc1,cc2,cc3,cc4,dvlspl,epsatm,psps%mqgrid_ff,psps%qgrid_ff,&
& vlspl(:,1),rloc,.true.,yp1,ypn,zion)

!DEBUG
!write(std_out,*)' psp3in : after psp2lo '
!ENDDEBUG

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_space,(psps%mqgrid_ff))
 ABI_ALLOCATE(work_spl,(psps%mqgrid_ff))
 call spline (psps%qgrid_ff,vlspl(:,1),psps%mqgrid_ff,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl)

!Second, compute KB energies and form factors and fit splines
 ekb(:)=zero

!Check if any nonlocal projectors are being used
 mproj=maxval(nproj)
 if (mproj>0) then

   ABI_ALLOCATE(ekb_sr,(psps%mpsang,mproj))
   ABI_ALLOCATE(ffspl_sr,(psps%mqgrid_ff,2,psps%mpsang,mproj))
   ABI_ALLOCATE(ekb_so,(psps%mpsang,mproj))
   ABI_ALLOCATE(ffspl_so,(psps%mqgrid_ff,2,psps%mpsang,mproj))

   call psp3nl(ekb_sr,ffspl_sr,h11s,h22s,h33s,h11p,h22p,h33p,h11d,h22d,&
&   h33d,h11f,mproj,psps%mpsang,psps%mqgrid_ff,psps%qgrid_ff,rrd,rrf,rrp,rrs)
   if(pspso/=0) &
&   call psp3nl(ekb_so,ffspl_so,zero,zero,zero,k11p,k22p,k33p,k11d,&
&   k22d,k33d,k11f,mproj,psps%mpsang,psps%mqgrid_ff,psps%qgrid_ff,rrd,rrf,rrp,rrs)


!  Convert ekb and ffspl
   iln0=0
   do jj=1,psps%lmnmax
     iln=indlmn(5,jj)
     if (iln>iln0) then
       iln0=iln
       if (indlmn(6,jj)<=1) then
         ekb(iln)=ekb_sr(1+indlmn(1,jj),indlmn(3,jj))
         ffspl(:,:,iln)=ffspl_sr(:,:,1+indlmn(1,jj),indlmn(3,jj))
       else
         ekb(iln)=ekb_so(1+indlmn(1,jj),indlmn(3,jj))
         ffspl(:,:,iln)=ffspl_so(:,:,1+indlmn(1,jj),indlmn(3,jj))
       end if
     end if
   end do

   ABI_DEALLOCATE(ekb_sr)
   ABI_DEALLOCATE(ffspl_sr)
   ABI_DEALLOCATE(ekb_so)
   ABI_DEALLOCATE(ffspl_so)

 end if

!DEBUG
!write(std_out,*)' psp3in : exit '
!stop
!ENDDEBUG

end subroutine psp3in
!!***
