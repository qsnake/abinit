!{\src2tex{textfont=tt}}
!!****f* ABINIT/opernlb_ylm
!! NAME
!! opernlb_ylm
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   from projected scalars to reciprocal space.
!! * Operate with the non-local projectors and the overlap matrix,
!!   from projected scalars to reciprocal space.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  choice=chooses possible output (see below)
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Vnl (NL operator)
!!  dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfacrelated to Sij (overlap)
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))= reduced projected scalars related to Sij (overlap)
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5, 51, 52 and signs=2)
!!                        - strain component (1:6) in the case (choice=2,signs=2) or (choice=6,signs=1)
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  kpg(npw,nkpg)=(k+G) components (if nkpg=3)
!!  matblk=dimension of the array ph3d
!!  ndgxdtfac=second dimension of dgxdtfac
!!  nincat=number of atoms in the subset here treated
!!  nkpg=second dimension of array kpg (0 or 3)
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  ucvol=unit cell volume (bohr^3)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! --if (paw_opt=0, 1 or 4)
!!    vectout(2,npw*nspinor)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1) <G|V_nonlocal|vect_start>
!!      if (choice=2) <G|dV_nonlocal/d(atm coord)|vect_start>
!!      if (choice=3) <G|dV_nonlocal/d(strain)|vect_start>
!!      if (choice=5) <G|dV_nonlocal/dk|vect_start>
!!      if (choice=51) <G|d(right)V_nonlocal/dk|vect_start>
!!      if (choice=52) <G|d(left)V_nonlocal/dk|vect_start>
!!      if (choice=53) <G|d(twist)V_nonlocal/dk|vect_start>
!!  if (paw_opt=2)
!!    vectout(2,npw*nspinor)=final vector in reciprocal space:
!!      if (choice=1) <G|V_nonlocal-lamdba.(I+S)|vect_start>
!!      if (choice=2) <G|d[V_nonlocal-lamdba.(I+S)]/d(atm coord)|vect_start>
!!      if (choice=3) <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_start>
!!      if (choice=5) <G|d[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
!!      if (choice=51) <G|d(right)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
!!      if (choice=52) <G|d(left)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
!!      if (choice=53) <G|d(twist)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
!! --if (paw_opt=3 or 4)
!!    svectout(2,npw*nspinor)=result of the aplication of Sij (overlap matrix)
!!                  to the input vect.:   (S= overlap matrix)
!!      if (choice=1) <G|I+S|vect_start>
!!      if (choice=2) <G|dS/d(atm coord)|vect_start>
!!      if (choice=3) <G|dS/d(strain)|vect_start>
!!      if (choice=5) <G|dS/dk|vect_start>
!!      if (choice=51) <G|d(right)S/dk|vect_start>
!!      if (choice=52) <G|d(left)S/dk|vect_start>
!!      if (choice=53) <G|d(twist)S/dk|vect_start>
!!
!! NOTES
!! Operate for one type of atom, and within this given type of atom,
!! for a subset of at most nincat atoms.
!!
!! PARENTS
!!      nonlop_ylm
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine opernlb_ylm(choice,cplex,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,&
&                      ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,&
&                      nspinor,paw_opt,ph3d,svect,ucvol,vect)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'opernlb_ylm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: choice,cplex,cplex_fac,dimffnl,ia3,idir,matblk,ndgxdtfac,nincat
 integer,intent(in) :: nkpg,nlmn,npw,nspinor,paw_opt
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: indlmn(6,nlmn),nloalg(5)
 real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat*(paw_opt/3),nspinor)
 real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(in) :: kpg(npw,nkpg),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: svect(2,npw*nspinor*(paw_opt/3)),vect(2,npw*nspinor)
!Local variables-------------------------------
!Arrays
!scalars
 integer :: fdb,fdf,ia,iaph3d,ii,il,ilmn,ipwshft,ispinor
 real(dp) :: scale,wt
 logical :: parity
!arrays
 integer :: ffnl_dir_dat(6)=(/3,4,4,2,2,3/)
 real(dp),allocatable :: dgxdtfac_(:,:,:),dgxdtfacs_(:,:,:),gxfac_(:,:),gxfacs_(:,:)
 complex(dpc),allocatable :: ztab(:)

! *************************************************************************

!Nothing to do when choice=4, 6 or 23
 if (choice==4.or.choice==6.or.choice==23) return

!Inits
 wt=four_pi/sqrt(ucvol)
 if (paw_opt/=3) then
   ABI_ALLOCATE(gxfac_,(2,nlmn))
   gxfac_(:,:)=zero
   if (choice>1) then
     ABI_ALLOCATE(dgxdtfac_,(2,ndgxdtfac,nlmn))
     if(ndgxdtfac>0) dgxdtfac_(:,:,:)=zero
   end if
 end if
 if (paw_opt>=3) then
   ABI_ALLOCATE(gxfacs_,(2,nlmn))
   gxfacs_(:,:)=zero
   if (choice>1) then
     ABI_ALLOCATE(dgxdtfacs_,(2,ndgxdtfac,nlmn))
     if (ndgxdtfac>0) dgxdtfacs_(:,:,:)=zero
   end if
 end if

!Loop on spinorial components
 do ispinor=1,nspinor

   ipwshft=(ispinor-1)*npw

!  Loops (blocking)
!  $OMP PARALLEL DEFAULT(PRIVATE) &
!  $OMP&SHARED(choice,nincat,nloalg,ia3,ipwshft,npw,nlmn,nspinor,wt,indlmn)
!  $OMP&SHARED(ffnl,gxfac,gxfac_sij,gxfac_,gxfacs_,dgxdtfac,dgxdtfac_sij,dgxdtfac_,dgxdtfacs_,kpg,ph3d,vect)
!  $OMP DO

!  Loop on atoms
   do ia=1,nincat
     iaph3d=ia;if (nloalg(1)>0) iaph3d=ia+ia3-1

!    Scale gxfac with 4pi/sqr(omega).(-i)^l
     if (paw_opt/=3) then
       do ilmn=1,nlmn
         il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
         scale=wt;if (il>1) scale=-scale
         if (parity) then
           gxfac_(1:cplex_fac,ilmn)=scale*gxfac(1:cplex_fac,ilmn,ia,ispinor)
           if (cplex_fac==1) gxfac_(2,ilmn)=zero
         else
           gxfac_(2,ilmn)=-scale*gxfac(1,ilmn,ia,ispinor)
           if (cplex_fac==2) then
             gxfac_(1,ilmn)=scale*gxfac(2,ilmn,ia,ispinor)
           else
             gxfac_(1,ilmn)=zero
           end if
         end if
       end do
       if (choice>1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn)=scale*dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
             if (cplex_fac==1) dgxdtfac_(2,1:ndgxdtfac,ilmn)=zero
           else
             do ii=1,ndgxdtfac
               dgxdtfac_(2,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
               if (cplex_fac==2) then
                 dgxdtfac_(1,ii,ilmn)=scale*dgxdtfac(2,ii,ilmn,ia,ispinor)
               else
                 dgxdtfac_(1,ii,ilmn)=zero
               end if
             end do
           end if
         end do
       end if
     end if

!    Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
     if (paw_opt>=3) then
       do ilmn=1,nlmn
         il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
         scale=wt;if (il>1) scale=-scale
         if (parity) then
           gxfacs_(1:cplex,ilmn)=scale*gxfac_sij(1:cplex,ilmn,ia,ispinor)
           if (cplex==1) gxfacs_(2,ilmn)=zero
         else
           gxfacs_(2,ilmn)=-scale*gxfac_sij(1,ilmn,ia,ispinor)
           if (cplex==2) then
             gxfacs_(1,ilmn)=scale*gxfac_sij(2,ilmn,ia,ispinor)
           else
             gxfacs_(1,ilmn)=zero
           end if
         end if
       end do
       if (choice>1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn)=scale*dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
             if (cplex==1) dgxdtfacs_(2,1:ndgxdtfac,ilmn)=zero
           else
             do ii=1,ndgxdtfac
               dgxdtfacs_(2,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
               if (cplex==2) then
                 dgxdtfacs_(1,ii,ilmn)=scale*dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
               else
                 dgxdtfacs_(1,ii,ilmn)=zero
               end if
             end do
           end if
         end do
       end if
     end if

     ABI_ALLOCATE(ztab,(npw))

!    Compute <g|Vnl|c> (or derivatives) for each plane wave:

     if (paw_opt/=3) then

       ztab(:)=czero

       if (choice==1) then
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
         end do
       end if

       if (choice==2) then
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
         end do
         ztab(:)=two_pi*kpg(:,idir)*ztab(:)
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
         end do
       end if

       if (choice==3) then
         if (idir<=3) then
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)&
&             *cmplx(dgxdtfac_(1,1,ilmn)-gxfac_(1,ilmn),dgxdtfac_(2,1,ilmn)-gxfac_(2,ilmn),kind=dp)&
&             -ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         else
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)&
&             -ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         end if
       end if

       if (choice==5) then ! full derivative w.r.t. k
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)&
&           +ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
         end do
       end if

       if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
         end do
       end if

       if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
         end do
       end if

       if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi> -
!        <G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
         fdf = ffnl_dir_dat(2*idir-1)
         fdb = ffnl_dir_dat(2*idir)
         do ilmn=1,nlmn
           ztab(:)=ztab(:) + &
&           ffnl(:,fdf,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) - &
&           ffnl(:,fdb,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
         end do
       end if

       ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
       vect(1,1+ipwshft:npw+ipwshft)=vect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
       vect(2,1+ipwshft:npw+ipwshft)=vect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))

     end if

!    Compute <g|S|c> (or derivatives) for each plane wave:

     if (paw_opt>=3) then

       ztab(:)=czero

       if (choice==1) then
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
         end do
       end if

       if (choice==2) then
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
         end do
         ztab(:)=two_pi*kpg(:,idir)*ztab(:)
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
         end do
       end if

       if (choice==3) then
         if (idir<=3) then
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)&
&             *cmplx(dgxdtfacs_(1,1,ilmn)-gxfacs_(1,ilmn),dgxdtfacs_(2,1,ilmn)-gxfacs_(2,ilmn),kind=dp)&
&             -ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         else
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)&
&             -ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if
       end if

       if (choice==5) then ! full derivative w.r.t. k
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)&
&           +ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
         end do
       end if

       if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
         end do
       end if

       if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
         do ilmn=1,nlmn
           ztab(:)=ztab(:)+ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
         end do
       end if

       if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi> -
!        <G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
         fdf = ffnl_dir_dat(2*idir-1)
         fdb = ffnl_dir_dat(2*idir)
         do ilmn=1,nlmn
           ztab(:)=ztab(:) + &
&           ffnl(:,fdf,ilmn)*cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) - &
&           ffnl(:,fdb,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
         end do
       end if

       ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
       svect(1,1+ipwshft:npw+ipwshft)=svect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
       svect(2,1+ipwshft:npw+ipwshft)=svect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
     end if

     ABI_DEALLOCATE(ztab)

!    End loop on atoms
   end do
!  $OMP END DO
!  $OMP END PARALLEL

!  End loop on spinors
 end do

 if (paw_opt/=3) then
   ABI_DEALLOCATE(gxfac_)
   if (choice>1) ABI_DEALLOCATE(dgxdtfac_)
 end if
 if (paw_opt>=3) then
   ABI_DEALLOCATE(gxfacs_)
   if (choice>1) ABI_DEALLOCATE(dgxdtfacs_)
 end if

end subroutine opernlb_ylm

!!***
