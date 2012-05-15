!{\src2tex{textfont=tt}}
!! ===================================================
!! This module contains functions used to manipulate
!! variables of structured datatype cprj_type.
!! cprj_type variables are <p_lmn|Cnk> projected
!! quantities where |p_lmn> are non-local projectors
!!                  |Cnk> are wave functions
!! ===================================================

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!!****f* ABINIT/cprj_alloc
!! NAME
!! cprj_alloc
!!
!! FUNCTION
!! Allocation of a cprj datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ncpgr=number of gradients to be allocated
!!  nlmn(:)=sizes of cprj%cp
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(cprj_type)>= cprj datastructure
!!
!! PARENTS
!!      accrho3,berry_linemin,berryphase_new,calc_optical_mels,calc_sigc_me
!!      calc_sigx_me,calc_vhxc_me,calc_wf_qp,cchi0,cchi0q0,cchi0q0_intraband
!!      cgwf,cgwf3,check_completeness,classify_bands,cohsex_me,cprj_utils_mpi
!!      ctocprj,debug_tools,dyfnl3,energy,exc_build_block,exc_build_ham
!!      exc_plot,extrapwf,getgh1c,getgsc,initberry,ks_ddiago,loper3
!!      m_cprj_bspline,m_electronpositron,m_shirley,m_wfs,mag_loc_k
!!      make_grad_berry,nstpaw3,optics_paw,optics_paw_core,outkss
!!      partial_dos_fractions_paw,paw_symcprj,pawmkaewf,pawmkrhoij,rdm,scfcv
!!      scfcv3,setup_positron,sigma,smatrix_pawinit,suscep_stat,update_mmat
!!      vtorho,vtorho3,vtowfk,vtowfk3,wfd_pawrhoij,wfd_vnlpsi
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_alloc(cprj,ncpgr,nlmn)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_alloc'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncpgr
!arrays
 integer,intent(in) :: nlmn(:)
 type(cprj_type),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim,nn

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2);nn=size(nlmn,dim=1)
 if (nn/=n1dim) then
   write(std_out,*)"Error in cprj_alloc: wrong sizes !",nn,n1dim
   stop
 end if
!write(std_out,*) "cprj_alloc ndim = ", n1dim, n2dim
 do jj=1,n2dim
   do ii=1,n1dim
     nullify (cprj(ii,jj)%cp)
     nullify (cprj(ii,jj)%dcp)

     nn=nlmn(ii)
     cprj(ii,jj)%nlmn=nn
     ABI_ALLOCATE(cprj(ii,jj)%cp,(2,nn))
!    XG 080820 Was needed to get rid of problems with test paral#R with four procs
     cprj(ii,jj)%cp=zero
!    END XG 080820

     cprj(ii,jj)%ncpgr=ncpgr
     if (ncpgr>0) then
       ABI_ALLOCATE(cprj(ii,jj)%dcp,(2,ncpgr,nn))
       cprj(ii,jj)%dcp=zero
     end if
   end do
 end do
end subroutine cprj_alloc
!!***

!!****f* ABINIT/cprj_free
!! NAME
!! cprj_free
!!
!! FUNCTION
!! Deallocation of a cprj datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(cprj_type)>= cprj datastructure
!!
!! PARENTS
!!      accrho3,berry_linemin,berryphase_new,calc_optical_mels,calc_sigc_me
!!      calc_sigx_me,calc_vhxc_me,calc_wf_qp,cchi0,cchi0q0,cchi0q0_intraband
!!      cgwf,cgwf3,check_completeness,classify_bands,cohsex_me,cprj_utils_mpi
!!      ctocprj,debug_tools,dyfnl3,energy,exc_build_block,exc_build_ham
!!      exc_plot,extrapwf,getgh1c,getgsc,ks_ddiago,loper3,m_cprj_bspline
!!      m_efield,m_electronpositron,m_scf_history,m_shirley,m_wfs,mag_loc_k
!!      make_grad_berry,nstpaw3,optics_paw,optics_paw_core,outkss
!!      partial_dos_fractions_paw,paw_symcprj,pawmkaewf,pawmkrhoij,rdm,scfcv
!!      scfcv3,setup_positron,sigma,smatrix_pawinit,suscep_stat,update_mmat
!!      vtorho,vtorho3,vtowfk,vtowfk3,wfd_pawrhoij,wfd_vnlpsi
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_free(cprj)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_free'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(cprj_type),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2)
!write(std_out,*) "cprj_free ndim = ", n1dim, n2dim
 do jj=1,n2dim
   do ii=1,n1dim
     if (associated(cprj(ii,jj)%cp))  then
       ABI_DEALLOCATE(cprj(ii,jj)%cp)
     end if
     if (associated(cprj(ii,jj)%dcp))  then
       ABI_DEALLOCATE(cprj(ii,jj)%dcp)
     end if
   end do
 end do
end subroutine cprj_free
!!***

!!****f* ABINIT/cprj_nullify
!! NAME
!! cprj_nullify
!!
!! FUNCTION
!! Nullify (set to null) a cprj datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(cprj_type)>= cprj datastructure
!!
!! PARENTS
!!      m_scf_history,suscep_stat
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_nullify(cprj)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_nullify'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(cprj_type),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2)
 do jj=1,n2dim
   do ii=1,n1dim
     cprj(ii,jj)%nlmn=0
     cprj(ii,jj)%ncpgr=0
     nullify(cprj(ii,jj)%cp)
     nullify(cprj(ii,jj)%dcp)
   end do
 end do
end subroutine cprj_nullify
!!***

!!****f* ABINIT/cprj_set_zero
!! NAME
!! cprj_set_zero
!!
!! FUNCTION
!! Set to zero all arrays in a cprj datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(cprj_type)>= cprj datastructure
!!
!! PARENTS
!!      cgwf3
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_set_zero(cprj)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_set_zero'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(cprj_type),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2)
 do jj=1,n2dim
   do ii=1,n1dim
     if (cprj(ii,jj)%nlmn>0)  cprj(ii,jj)%cp(:,:)=zero
     if (cprj(ii,jj)%ncpgr>0) cprj(ii,jj)%dcp(:,:,:)=zero
   end do
 end do
end subroutine cprj_set_zero
!!***

!!****f* ABINIT/cprj_copy
!! NAME
!! cprj_copy
!!
!! FUNCTION
!! Copy a cprj datastructure into another
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  icpgr= (optional argument) if present, only component icpgr of
!!         input cprj gradient is copied into output cprj
!!         Not used if cprj(:,:)%ncpgr<icpgr
!!  cprj_in(:,:) <type(cprj_type)>= input cprj datastructure
!!
!! OUTPUT
!!  cprj_out(:,:) <type(cprj_type)>= output cprj datastructure
!!
!! NOTES
!!  MG: What about an option to report a pointer to cprj_in?
!!
!! PARENTS
!!      berry_linemin,berryphase_new,calc_sigc_me,calc_sigx_me,cchi0q0
!!      cchi0q0_intraband,cgwf,classify_bands,cohsex_me,corrmetalwf1
!!      cprj_utils_mpi,dyfnl3,extrapwf,getgh1c,getgsc,loper3,m_electronpositron
!!      m_wfs,make_grad_berry,nstpaw3,outkss,paw_symcprj,setup_positron
!!      update_mmat
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_copy(cprj_in,cprj_out,&
 &                    icpgr) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_copy'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: icpgr
!arrays
 type(cprj_type),intent(in) :: cprj_in(:,:)
 type(cprj_type),intent(inout) :: cprj_out(:,:)

!Local variables-------------------------------
 integer :: ii,jj,kk,n1dim_in,n1dim_out,n2dim_in,n2dim_out,ncpgr_in,ncpgr_out,nlmn
 logical :: has_icpgr
 character(len=500) :: msg

! *************************************************************************

 n1dim_in=size(cprj_in,dim=1); n1dim_out=size(cprj_out,dim=1)
 n2dim_in=size(cprj_in,dim=2); n2dim_out=size(cprj_out,dim=2)
 ncpgr_in=cprj_in(1,1)%ncpgr;  ncpgr_out=cprj_out(1,1)%ncpgr

 if (n1dim_in/=n1dim_out) then
   write(msg,'(a,2(1x,i0))')" Error in cprj_copy: n1 wrong sizes ",n1dim_in,n1dim_out
   MSG_ERROR(msg)
 end if
 if (n2dim_in/=n2dim_out) then
   write(msg,'(a,2(1x,i0))')" Error in cprj_copy: n2 wrong sizes ",n2dim_in,n2dim_out
   MSG_ERROR(msg)
 end if
 if (ncpgr_in<ncpgr_out)  then
   write(msg,'(a,2(1x,i0))')" Error in cprj_copy: ncpgr wrong sizes ",ncpgr_in,ncpgr_out
   MSG_ERROR(msg)
 end if

 do jj=1,n2dim_in
   do ii=1,n1dim_in
     nlmn=cprj_in(ii,jj)%nlmn
     cprj_out(ii,jj)%nlmn =nlmn
     do kk=1,nlmn
       cprj_out(ii,jj)%cp(1:2,kk)=cprj_in(ii,jj)%cp(1:2,kk)
     end do
   end do
 end do

 if (ncpgr_in>0) then
   has_icpgr=present(icpgr)
   if (has_icpgr) has_icpgr=(ncpgr_out>0.and.icpgr>0.or.icpgr<=ncpgr_in)

   if (has_icpgr) then
     do jj=1,n2dim_in
       do ii=1,n1dim_in
         nlmn=cprj_in(ii,jj)%nlmn
         do kk=1,nlmn
           cprj_out(ii,jj)%dcp(1:2,1,kk)=cprj_in(ii,jj)%dcp(1:2,icpgr,kk)
         end do
       end do
     end do
   else
     if (ncpgr_out>=ncpgr_in) then
       do jj=1,n2dim_in
         do ii=1,n1dim_in
           nlmn=cprj_in(ii,jj)%nlmn
           do kk=1,nlmn
             cprj_out(ii,jj)%dcp(1:2,1:ncpgr_in,kk)=cprj_in(ii,jj)%dcp(1:2,1:ncpgr_in,kk)
           end do
         end do
       end do
     end if
   end if
 end if

end subroutine cprj_copy
!!***


!!****f* ABINIT/cprj_axpby
!! NAME
!! cprj_axpby
!!
!! FUNCTION
!! Apply AXPBY (blas-like) operation with 2 cprj datastructures:
!!  cprjy(:,:) <- alpha.cprjx(:,:)+beta.cprjy(:,:)
!!  alpha and beta are REAL scalars
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  alpha,beta= alpha,beta REAL factors
!!  cprjx(:,:) <type(cprj_type)>= input cprjx datastructure
!!
!! SIDE EFFECTS
!!  cprjy(:,:) <type(cprj_type)>= input/output cprjy datastructure
!!
!! PARENTS
!!      cgwf3,getdc1
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_axpby(alpha,beta,cprjx,cprjy)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_axpby'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: alpha,beta
!arrays
 type(cprj_type),intent(in) :: cprjx(:,:)
 type(cprj_type),intent(inout) :: cprjy(:,:)
!Local variables-------------------------------
 integer :: ii,jj,kk,n1dimx,n1dimy,n2dimx,n2dimy,ncpgrx,ncpgry,nlmn

! *************************************************************************

 n1dimy=size(cprjy,dim=1);n2dimy=size(cprjy,dim=2);ncpgry=cprjy(1,1)%ncpgr
 if (abs(alpha)>tol16) then
   n1dimx=size(cprjx,dim=1);n2dimx=size(cprjx,dim=2);ncpgrx=cprjx(1,1)%ncpgr
   if (n1dimx/=n1dimy) stop "Error in cprj_axpby: n1 wrong sizes ! "
   if (n2dimx/=n2dimy) stop "Error in cprj_axpby: n2 wrong sizes ! "
   if (ncpgrx/=ncpgry) stop "Error in cprj_axpby: ncpgr wrong sizes ! "
 end if

 if (abs(alpha)<=tol16) then
   do jj=1,n2dimy
     do ii=1,n1dimy
       nlmn=cprjy(ii,jj)%nlmn
       do kk=1,nlmn
         cprjy(ii,jj)%cp(1:2,kk)=beta*cprjy(ii,jj)%cp(1:2,kk)
       end do
     end do
   end do
   if (ncpgry>0) then
     do jj=1,n2dimy
       do ii=1,n1dimy
         nlmn=cprjy(ii,jj)%nlmn
         do kk=1,nlmn
           cprjy(ii,jj)%dcp(1:2,1:ncpgry,kk)=beta*cprjy(ii,jj)%dcp(1:2,1:ncpgry,kk)
         end do
       end do
     end do
   end if
 else if (abs(beta)<=tol16) then
   do jj=1,n2dimx
     do ii=1,n1dimx
       nlmn=cprjx(ii,jj)%nlmn
       cprjy(ii,jj)%nlmn=nlmn
       do kk=1,nlmn
         cprjy(ii,jj)%cp(1:2,kk)=alpha*cprjx(ii,jj)%cp(1:2,kk)
       end do
     end do
   end do
   if (ncpgrx>0) then
     do jj=1,n2dimx
       do ii=1,n1dimx
         nlmn=cprjx(ii,jj)%nlmn
         do kk=1,nlmn
           cprjy(ii,jj)%dcp(1:2,1:ncpgrx,kk)=alpha*cprjx(ii,jj)%dcp(1:2,1:ncpgrx,kk)
         end do
       end do
     end do
   end if
 else  ! alpha/=0 and beta/=0
   do jj=1,n2dimx
     do ii=1,n1dimx
       nlmn=cprjx(ii,jj)%nlmn
       cprjy(ii,jj)%nlmn=nlmn
       do kk=1,nlmn
         cprjy(ii,jj)%cp(1:2,kk)=alpha*cprjx(ii,jj)%cp(1:2,kk) &
&         +beta *cprjy(ii,jj)%cp(1:2,kk)
       end do
     end do
   end do
   if (ncpgrx>0) then
     do jj=1,n2dimx
       do ii=1,n1dimx
         nlmn=cprjx(ii,jj)%nlmn
         do kk=1,nlmn
           cprjy(ii,jj)%dcp(1:2,1:ncpgrx,kk)=alpha*cprjx(ii,jj)%dcp(1:2,1:ncpgrx,kk) &
&           +beta *cprjy(ii,jj)%dcp(1:2,1:ncpgrx,kk)
         end do
       end do
     end do
   end if
 end if

end subroutine cprj_axpby
!!***


!!****f* ABINIT/cprj_zaxpby
!! NAME
!! cprj_zaxpby
!!
!! FUNCTION
!! Apply ZAXPBY (blas-like) operation with 2 cprj datastructures:
!!  cprjy(:,:) <- alpha.cprjx(:,:)+beta.cprjy(:,:)
!!  alpha and beta are COMPLEX scalars
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  alpha(2),beta(2)= alpha,beta COMPLEX factors
!!  cprjx(:,:) <type(cprj_type)>= input cprjx datastructure
!!
!! SIDE EFFECTS
!!  cprjy(:,:) <type(cprj_type)>= input/output cprjy datastructure
!!
!! PARENTS
!!      corrmetalwf1,extrapwf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_zaxpby(alpha,beta,cprjx,cprjy)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_zaxpby'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: alpha(2),beta(2)
!arrays
 type(cprj_type),intent(in) :: cprjx(:,:)
 type(cprj_type),intent(inout) :: cprjy(:,:)
!Local variables-------------------------------
 integer :: ii,jj,kk,ll,n1dimx,n1dimy,n2dimx,n2dimy,ncpgrx,ncpgry,nlmn
 real(dp) :: cp1,cp2,norma,normb

! *************************************************************************

 norma=alpha(1)**2+alpha(2)**2
 normb=beta(1) **2+beta(2) **2
 n1dimy=size(cprjy,dim=1);n2dimy=size(cprjy,dim=2);ncpgry=cprjy(1,1)%ncpgr
 if (norma>tol16) then
   n1dimx=size(cprjx,dim=1);n2dimx=size(cprjx,dim=2);ncpgrx=cprjx(1,1)%ncpgr
   if (n1dimx/=n1dimy) stop "Error in cprj_zaxpby: n1 wrong sizes ! "
   if (n2dimx/=n2dimy) stop "Error in cprj_zaxpby: n2 wrong sizes ! "
   if (ncpgrx/=ncpgry) stop "Error in cprj_zaxpby: ncpgr wrong sizes ! "
 end if

 if (norma<=tol16) then
   do jj=1,n2dimy
     do ii=1,n1dimy
       nlmn=cprjy(ii,jj)%nlmn
       do kk=1,nlmn
         cp1=beta(1)*cprjy(ii,jj)%cp(1,kk)-beta(2)*cprjy(ii,jj)%cp(2,kk)
         cp2=beta(1)*cprjy(ii,jj)%cp(2,kk)+beta(2)*cprjy(ii,jj)%cp(1,kk)
         cprjy(ii,jj)%cp(1,kk)=cp1
         cprjy(ii,jj)%cp(2,kk)=cp2
       end do
     end do
   end do
   if (ncpgry>0) then
     do jj=1,n2dimy
       do ii=1,n1dimy
         nlmn=cprjy(ii,jj)%nlmn
         do kk=1,nlmn
           do ll=1,ncpgry
             cp1=beta(1)*cprjy(ii,jj)%dcp(1,ll,kk)-beta(2)*cprjy(ii,jj)%dcp(2,ll,kk)
             cp2=beta(1)*cprjy(ii,jj)%dcp(2,ll,kk)+beta(2)*cprjy(ii,jj)%dcp(1,ll,kk)
             cprjy(ii,jj)%dcp(1,ll,kk)=cp1
             cprjy(ii,jj)%dcp(2,ll,kk)=cp2
           end do
         end do
       end do
     end do
   end if
 else if (normb<=tol16) then
   do jj=1,n2dimx
     do ii=1,n1dimx
       nlmn=cprjx(ii,jj)%nlmn
       cprjy(ii,jj)%nlmn=nlmn
       do kk=1,nlmn
         cprjy(ii,jj)%cp(1,kk)=alpha(1)*cprjx(ii,jj)%cp(1,kk)-alpha(2)*cprjx(ii,jj)%cp(2,kk)
         cprjy(ii,jj)%cp(2,kk)=alpha(1)*cprjx(ii,jj)%cp(2,kk)+alpha(2)*cprjx(ii,jj)%cp(1,kk)
       end do
     end do
   end do
   if (ncpgrx>0) then
     do jj=1,n2dimx
       do ii=1,n1dimx
         nlmn=cprjx(ii,jj)%nlmn
         do kk=1,nlmn
           cprjy(ii,jj)%dcp(1,1:ncpgrx,kk)=alpha(1)*cprjx(ii,jj)%dcp(1,1:ncpgrx,kk) &
&           -alpha(2)*cprjx(ii,jj)%dcp(2,1:ncpgrx,kk)
           cprjy(ii,jj)%dcp(2,1:ncpgrx,kk)=alpha(1)*cprjx(ii,jj)%dcp(2,1:ncpgrx,kk) &
&           +alpha(2)*cprjx(ii,jj)%dcp(1,1:ncpgrx,kk)
         end do
       end do
     end do
   end if
 else
   do jj=1,n2dimx
     do ii=1,n1dimx
       nlmn=cprjx(ii,jj)%nlmn
       cprjy(ii,jj)%nlmn =nlmn
       do kk=1,nlmn
         cp1=alpha(1)*cprjx(ii,jj)%cp(1,kk)-alpha(2)*cprjx(ii,jj)%cp(2,kk) &
&         +beta(1) *cprjy(ii,jj)%cp(1,kk)-beta(2) *cprjy(ii,jj)%cp(2,kk)
         cp2=alpha(1)*cprjx(ii,jj)%cp(2,kk)+alpha(2)*cprjx(ii,jj)%cp(1,kk) &
&         +beta(1) *cprjy(ii,jj)%cp(2,kk)+beta(2) *cprjy(ii,jj)%cp(1,kk)
         cprjy(ii,jj)%cp(1,kk)=cp1
         cprjy(ii,jj)%cp(2,kk)=cp2
       end do
     end do
   end do
   if (ncpgrx>0) then
     do jj=1,n2dimx
       do ii=1,n1dimx
         nlmn=cprjx(ii,jj)%nlmn
         do kk=1,nlmn
           do ll=1,ncpgrx
             cp1=alpha(1)*cprjx(ii,jj)%dcp(1,ll,kk)-alpha(2)*cprjx(ii,jj)%dcp(2,ll,kk) &
&             +beta(1) *cprjy(ii,jj)%dcp(1,ll,kk)-beta(2) *cprjy(ii,jj)%dcp(2,ll,kk)
             cp2=alpha(1)*cprjx(ii,jj)%dcp(2,ll,kk)+alpha(2)*cprjx(ii,jj)%dcp(1,ll,kk) &
&             +beta(1) *cprjy(ii,jj)%dcp(2,ll,kk)+beta(2) *cprjy(ii,jj)%dcp(1,ll,kk)
             cprjy(ii,jj)%dcp(1,ll,kk)=cp1
             cprjy(ii,jj)%dcp(2,ll,kk)=cp2
           end do
         end do
       end do
     end do
   end if
 end if

end subroutine cprj_zaxpby
!!***


!!****f* ABINIT/cprj_lincom
!! NAME
!! cprj_lincom
!!
!! FUNCTION
!! Compute a LINear COMbination of cprj datastructure:
!!  cprj_out(:,:) <--- Sum_i [ alpha_i . cprj_i(:,:) ]
!!  alpha_i are COMPLEX scalars
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  alpha(2,nn)= alpha COMPLEX factors
!!  cprj_in(:,:) <type(cprj_type)>= input cprj_in datastructure
!!  nn= number of cprj involved in the linear combination
!!
!! OUTPUT
!!  cprj_out(:,:) <type(cprj_type)>= output cprj_out datastructure
!!
!! NOTES
!!  cprj_in and cprj_out must be dimensionned as cprj_in(n1,n2*nn) and cprj_in(n1,n2)
!!
!! PARENTS
!!      extrapwf,getdc1
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_lincom(alpha,cprj_in,cprj_out,nn)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_lincom'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp),intent(in) :: alpha(2,nn)
!arrays
 type(cprj_type),intent(in) :: cprj_in(:,:)
 type(cprj_type),intent(inout) :: cprj_out(:,:)
!Local variables-------------------------------
 integer :: ii,in,jj,jn,kk,ll,n1in,n1out,n2in,n2out,ncpgrin,ncpgrout,nlmn
 real(dp) :: cp1,cp2

! *************************************************************************

 n1in=size(cprj_in,dim=1);n1out=size(cprj_out,dim=1)
 n2in=size(cprj_in,dim=2);n2out=size(cprj_out,dim=2)
 ncpgrin=cprj_in(1,1)%ncpgr;ncpgrout=cprj_out(1,1)%ncpgr
 if (n1in/=n1out) stop "Bug in cprj_lincom: n1 wrong sizes ! "
 if (n2in/=n2out*nn) stop "Bug in cprj_lincom: n2 wrong sizes ! "
 if (ncpgrin/=ncpgrout) stop "Bug in cprj_lincom: ncpgr wrong sizes ! "

 do jj=1,n2out
   do ii=1,n1out
     nlmn=cprj_in(ii,jj)%nlmn
     cprj_out(ii,jj)%nlmn=nlmn
     cprj_out(ii,jj)%cp(1:2,1:nlmn)=zero
     jn=jj
     do in=1,nn
       do kk=1,nlmn
         cp1=cprj_out(ii,jj)%cp(1,kk) &
&         +alpha(1,in)*cprj_in(ii,jn)%cp(1,kk)-alpha(2,in)*cprj_in(ii,jn)%cp(2,kk)
         cp2=cprj_out(ii,jj)%cp(2,kk) &
&         +alpha(1,in)*cprj_in(ii,jn)%cp(2,kk)+alpha(2,in)*cprj_in(ii,jn)%cp(1,kk)
         cprj_out(ii,jj)%cp(1,kk)=cp1
         cprj_out(ii,jj)%cp(2,kk)=cp2
       end do
       jn=jn+n2out
     end do
   end do
 end do

 if (ncpgrin>0) then
   do jj=1,n2out
     do ii=1,n1out
       nlmn=cprj_in(ii,jj)%nlmn
       cprj_out(ii,jj)%dcp(1:2,1:ncpgrin,1:nlmn)=zero
       jn=jj
       do in=1,nn
         do kk=1,nlmn
           do ll=1,ncpgrin
             cp1=cprj_out(ii,jj)%dcp(1,ll,kk) &
&             +alpha(1,in)*cprj_in(ii,jn)%dcp(1,ll,kk) &
&             -alpha(2,in)*cprj_in(ii,jn)%dcp(2,ll,kk)
             cp2=cprj_out(ii,jj)%dcp(2,ll,kk) &
&             +alpha(1,in)*cprj_in(ii,jn)%dcp(2,ll,kk) &
             +alpha(2,in)*cprj_in(ii,jn)%dcp(1,ll,kk)
             cprj_out(ii,jj)%dcp(1,ll,kk)=cp1
             cprj_out(ii,jj)%dcp(2,ll,kk)=cp2
           end do
         end do
         jn=jn+n2out
       end do
     end do
   end do
 end if

end subroutine cprj_lincom
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/paw_overlap
!! NAME
!! paw_overlap
!!
!! FUNCTION
!!  Helper function returning the onsite contribution to the overlap between two states.
!!
!! INPUTS
!!   spinor_comm= (optional) communicator over spinorial components
!!   typat(:)=The type of each atom.
!!   Pawtab(ntypat)<type(pawtab_type)>=paw tabulated starting data.
!!   cprj1,cprj2<Cprj_type>
!!     Projected wave functions <Proj_i|Cnk> with all NL projectors for the left and the right wavefunction,respectively.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function paw_overlap(cprj1,cprj2,typat,pawtab,spinor_comm) result(onsite)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_overlap'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: spinor_comm
!arrays
 integer,intent(in) :: typat(:)
 real(dp) :: onsite(2)
 type(cprj_type),intent(in) :: cprj1(:,:),cprj2(:,:)
 type(pawtab_type),intent(in) :: pawtab(:)


!Local variables-------------------------------
!scalars
 integer :: iatom,ilmn,itypat,j0lmn,jlmn,klmn,natom,nspinor,isp
 real(dp) :: sij
 character(len=500) :: msg

! *************************************************************************

 natom=SIZE(typat)

 if (SIZE(cprj1,DIM=1)/=SIZE(cprj2,DIM=1) .or. SIZE(cprj1,DIM=1)/=natom) then
   write(msg,'(a,3i4)')' Wrong size in typat, cprj1, cprj2 : ',natom,SIZE(cprj1),SIZE(cprj2)
   MSG_ERROR(msg)
 end if

 nspinor = SIZE(cprj1,DIM=2)

 onsite=zero
 do iatom=1,natom
   itypat=typat(iatom)
   do jlmn=1,pawtab(itypat)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       sij=pawtab(itypat)%sij(klmn); if (jlmn==ilmn) sij=sij*half
       if (ABS(sij)>tol16) then
         do isp=1,nspinor

           onsite(1)=onsite(1) + sij*(                                 &
&           cprj1(iatom,isp)%cp(1,ilmn) * cprj2(iatom,isp)%cp(1,jlmn) &
&           +cprj1(iatom,isp)%cp(2,ilmn) * cprj2(iatom,isp)%cp(2,jlmn) &
&           +cprj1(iatom,isp)%cp(1,jlmn) * cprj2(iatom,isp)%cp(1,ilmn) &
&           +cprj1(iatom,isp)%cp(2,jlmn) * cprj2(iatom,isp)%cp(2,ilmn) &
&           )

           onsite(2)=onsite(2) + sij*(                                  &
&           cprj1(iatom,isp)%cp(1,ilmn) * cprj2(iatom,isp)%cp(2,jlmn) &
&           -cprj1(iatom,isp)%cp(2,ilmn) * cprj2(iatom,isp)%cp(1,jlmn) &
&           +cprj1(iatom,isp)%cp(1,jlmn) * cprj2(iatom,isp)%cp(2,ilmn) &
&           -cprj1(iatom,isp)%cp(2,jlmn) * cprj2(iatom,isp)%cp(1,ilmn) &
&           )
         end do
       end if
     end do
   end do
 end do
 if (present(spinor_comm)) then
   call xsum_mpi(onsite,spinor_comm,isp)
 end if

end function paw_overlap
!!***

!!****f* ABINIT/cprj_output
!! NAME
!! cprj_output
!!
!! FUNCTION
!! Output a cprj. Useful for debugging.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprj(:,:) <type(cprj_type)>= cprj datastructure
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cprj_output(cprj)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_output'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalar
!arrays
 type(cprj_type),intent(in) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,kk,nlmn,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1)
 n2dim=size(cprj,dim=2)

 write(std_out,'(a)')' cprj_output '
 do jj=1,n2dim
   do ii=1,n1dim
     write(std_out,'(a,i4,a,i4)')'atom ',ii,' band*k ',jj
     nlmn=cprj(ii,jj)%nlmn
     do kk=1,nlmn
       write(std_out,'(2f12.8)')cprj(ii,jj)%cp(1,kk),cprj(ii,jj)%cp(2,kk)
     end do
   end do
 end do

end subroutine cprj_output
!!***
