!{\src2tex{textfont=tt}}
!!****f* ABINIT/opernlc_ylm
!! NAME
!! opernlc_ylm
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   in order to reduce projected scalars
!! * Operate with the non-local projectors and the overlap matrix,
!!   in order to reduce projected scalars
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms (gives the absolute index of
!!                 an atom from its rank in a block of atoms)
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_enl=1 if enl factors are real, 2 if they are complex
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  dgxdt(cplex,ndgxdt,nlmn,nincat)=grads of projected scalars (only if optder>0)
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  enl(cplex_enl*dimenl1,dimenl2,nspinortot**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1, 2 or 4 ====
!!                      (Real or complex, hermitian) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  gx(cplex,nlmn,nincat*abs(enl_opt))= projected scalars
!!  iatm=absolute rank of first atom of the current block of atoms
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  itypat=type of atoms
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  ndgxdt=second dimension of dgxdt
!!  ndgxdtfac=second dimension of dgxdtfac
!!  nincat=number of atoms in the subset here treated
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nspinor= number of spinorial components of the wavefunctions (on current proc)
!!  nspinortot=total number of spinorial components of the wavefunctions
!!  optder=0=only gxfac is computed, 1=both gxfac and dgxdtfac are computed
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  sij(nlm*(nlmn+1)/2)=overlap matrix components (only if paw_opt=2, 3 or 4)
!!
!! OUTPUT
!!  if (paw_opt=0, 1, 2 or 4)
!!    gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  if (paw_opt=3 or 4)
!!    gxfac_sij(cplex,nlmn,nincat,nspinor)= reduced projected scalars related to Sij (overlap)
!!  if (optder==1.and.paw_opt=0, 1, 2 or 4)
!!    dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Vnl (NL operator)
!!  if (optder==1.and.paw_opt=3 or 4)
!!    dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Sij (overlap)
!!
!! NOTES
!! This routine operates for one type of atom, and within this given type of atom,
!! for a subset of at most nincat atoms.
!!
!! PARENTS
!!      nonlop_ylm
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine opernlc_ylm(atindx1,cplex,cplex_enl,cplex_fac,dgxdt,dgxdtfac,dgxdtfac_sij,&
&                      dimenl1,dimenl2,enl,gx,gxfac,gxfac_sij,iatm,indlmn,itypat,&
&                      lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,nincat,nlmn,&
&                      nspinor,nspinortot,optder,paw_opt,sij)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'opernlc_ylm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,iatm,itypat,natom,ndgxdt,ndgxdtfac
 integer,intent(in) :: nincat,nspinor,nspinortot,optder,paw_opt
 integer,intent(inout) :: nlmn
 real(dp) :: lambda
 type(MPI_type) , intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,nlmn)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
 real(dp),intent(inout) :: gx(cplex,nlmn,nincat,nspinor)
 real(dp),intent(in) :: sij(((paw_opt+1)/3)*nlmn*(nlmn+1)/2)
 real(dp),intent(out) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(out) :: gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))

!Local variables-------------------------------
!Arrays
!scalars
 integer :: ia,ierr,ijlmn,ijspin,ilm,ilmn,i0lmn,iln,index_enl,ispinor,ispinor_index
 integer :: j0lmn,jilmn,jispin,jjlmn,jlm,jlmn,jspinor,jspinor_index,mu,shift
 real(dp) :: sijr
!arrays
 real(dp) :: enl_(2),gxfi(2),gxi(cplex),gxj(cplex)
 real(dp),allocatable :: dgxdtfac_(:,:,:,:,:),gxfac_(:,:,:,:),gxfj(:,:)

! *************************************************************************

!Parallelization over spinors treatment
 shift=0;if (mpi_enreg%paral_spin==1) shift=mpi_enreg%me_spin

!Accumulate gxfac related to non-local operator (Norm-conserving)
!-------------------------------------------------------------------
 if (paw_opt==0) then                  ! Enl is E(Kleinman-Bylander)
   if (cplex_enl==2    ) stop "opernlc_ylm: BUG - invalid cplex_enl=2 !"
   if (cplex_fac/=cplex) stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex) !"
   do ispinor=1,nspinor
     ispinor_index=ispinor+shift
     do ia=1,nincat
       do ilmn=1,nlmn
         iln=indlmn(5,ilmn)
         enl_(1)=enl(iln,itypat,ispinor_index)
         gxfac(1:cplex,ilmn,ia,ispinor)=enl_(1)*gx(1:cplex,ilmn,ia,ispinor)
       end do
     end do
   end do
 end if

!Accumulate gxfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
 if (paw_opt==1.or.paw_opt==2.or.paw_opt==4) then        ! Enl is psp strength Dij
   gxfac(1:cplex_fac,1:nlmn,1:nincat,1:nspinor)=zero      ! or (Dij-lambda.Sij)
!  === Diagonal term(s) (up-up, down-down)
!  1-Enl is real
   if (cplex_enl==1) then
     do ispinor=1,nspinor
       ispinor_index=ispinor+shift
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl(jjlmn,index_enl,ispinor_index)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
           gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1)=enl(ijlmn,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
             gxfac(1:cplex,ilmn,ia,ispinor)=gxfac(1:cplex,ilmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
             gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxi(1:cplex)
           end do
         end do
       end do
     end do
!    2-Enl is complex
   else
     if (cplex_fac/=cplex_enl) then
       stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex_enl !"
     end if
     if (nspinortot==1) then    !===== when nspinor=1, D_ij=D_ji
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,1)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,1)
           gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)+enl_(1)*gxj(1)
           gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(2)*gxj(1)
           if (cplex==2) then
             gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)-enl_(2)*gxj(2)
             gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(1)*gxj(2)
           end if
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,1)
             gxfac(1,ilmn,ia,1)=gxfac(1,ilmn,ia,1)+enl_(1)*gxj(1)
             gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)+enl_(1)*gxi(1)
             gxfac(2,ilmn,ia,1)=gxfac(2,ilmn,ia,1)+enl_(2)*gxj(1)
             gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(2)*gxi(1)
             if (cplex==2) then
               gxfac(1,ilmn,ia,1)=gxfac(1,ilmn,ia,1)-enl_(2)*gxj(2)
               gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)-enl_(2)*gxi(2)
               gxfac(2,ilmn,ia,1)=gxfac(2,ilmn,ia,1)+enl_(1)*gxj(2)
               gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(1)*gxi(2)
             end if
           end do
         end do
       end do
     else                    !===== when nspinor=2, D_ij=D_ji^*
       do ispinor=1,nspinor
         ispinor_index=ispinor+shift
         do ia=1,nincat
           index_enl=atindx1(iatm+ia)
           do jlmn=1,nlmn
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             enl_(1)=enl(2*jjlmn-1,index_enl,ispinor_index)  ! enl_ii is real for up-up and dn-dn
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
             gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
             gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
               gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)+enl_(1)*gxj(1)
               gxfac(1,jlmn,ia,ispinor)=gxfac(1,jlmn,ia,ispinor)+enl_(1)*gxi(1)
               gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+enl_(2)*gxj(1)
               gxfac(2,jlmn,ia,ispinor)=gxfac(2,jlmn,ia,ispinor)-enl_(2)*gxi(1)
               if (cplex==2) then
                 gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)-enl_(2)*gxj(2)
                 gxfac(1,jlmn,ia,ispinor)=gxfac(1,jlmn,ia,ispinor)+enl_(2)*gxi(2)
                 gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+enl_(1)*gxj(2)
                 gxfac(2,jlmn,ia,ispinor)=gxfac(2,jlmn,ia,ispinor)+enl_(1)*gxi(2)
               end if
             end do
           end do
         end do
       end do

     end if !nspinortot

   end if !complex_enl
!  === Off-diagonal term(s) (up-down, down-up)
!  --- No parallelization over spinors ---
   if (nspinortot==2.and.nspinor==nspinortot) then
     if (cplex_enl/=2    ) stop "opernlc_ylm: BUG - invalid cplex_enl=2 !"
     if (cplex_fac/=cplex) stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex) !"
     do ispinor=1,nspinortot
       jspinor=3-ispinor
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor )
           gxi(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,jspinor)
           gxfac(1,jlmn,ia,jspinor)=gxfac(1,jlmn,ia,jspinor)+enl_(1)*gxi(1)
           gxfac(2,jlmn,ia,jspinor)=gxfac(2,jlmn,ia,jspinor)-enl_(2)*gxi(1)
           if (cplex==2) then
             gxfac(1,jlmn,ia,jspinor)=gxfac(1,jlmn,ia,jspinor)+enl_(2)*gxi(2)
             gxfac(2,jlmn,ia,jspinor)=gxfac(2,jlmn,ia,jspinor)+enl_(1)*gxi(2)
           end if
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
             gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)+enl_(1)*gxj(1)
             gxfac(1,jlmn,ia,jspinor)=gxfac(1,jlmn,ia,jspinor)+enl_(1)*gxi(1)
             gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+enl_(2)*gxj(1)
             gxfac(2,jlmn,ia,jspinor)=gxfac(2,jlmn,ia,jspinor)-enl_(2)*gxi(1)
             if (cplex==2) then
               gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)-enl_(2)*gxj(2)
               gxfac(1,jlmn,ia,jspinor)=gxfac(1,jlmn,ia,jspinor)+enl_(2)*gxi(2)
               gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+enl_(1)*gxj(2)
               gxfac(2,jlmn,ia,jspinor)=gxfac(2,jlmn,ia,jspinor)+enl_(1)*gxi(2)
             end if
           end do
         end do
       end do
     end do

!    --- Parallelization over spinors ---
   else if (nspinortot==2.and.nspinor/=nspinortot) then
     if (cplex_enl/=2) stop "opernlc_ylm: BUG - invalid cplex_enl/=2 !"
     if (cplex_fac/=2) stop "opernlc_ylm: BUG - invalid cplex_fac/=2 !"
     ABI_ALLOCATE(gxfac_,(cplex_fac,nlmn,nincat,nspinortot))
     gxfac_(:,:,:,:)=zero
     ispinor_index=mpi_enreg%me_spin+1
     jspinor_index=3-ispinor_index
     if (ispinor_index==1) then
       ijspin=3;jispin=4
     else
       ijspin=4;jispin=3
     end if
     do ia=1,nincat
       index_enl=atindx1(iatm+ia)
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         do ilmn=1,nlmn
           i0lmn=ilmn*(ilmn-1)/2
           if (ilmn<=jlmn) then
             ijlmn=j0lmn+ilmn
             enl_(1)= enl(2*ijlmn-1,index_enl,ijspin)
             enl_(2)=-enl(2*ijlmn  ,index_enl,ijspin)
           else
             jilmn=i0lmn+jlmn
             enl_(1:2)=enl(2*jilmn-1:2*jilmn,index_enl, jispin)
           end if
           gxi(1:cplex)=gx(1:cplex,ilmn,ia,1)
           gxfac_(1,jlmn,ia,jspinor_index)=gxfac_(1,jlmn,ia,jspinor_index)+enl_(1)*gxi(1)
           gxfac_(2,jlmn,ia,jspinor_index)=gxfac_(2,jlmn,ia,jspinor_index)+enl_(2)*gxi(1)
           if (cplex==2) then
             gxfac_(1,jlmn,ia,jspinor_index)=gxfac_(1,jlmn,ia,jspinor_index)-enl_(2)*gxi(2)
             gxfac_(2,jlmn,ia,jspinor_index)=gxfac_(2,jlmn,ia,jspinor_index)+enl_(1)*gxi(2)
           end if
         end do !ilmn
       end do !jlmn
     end do !iat
     call xsum_mpi(gxfac_,mpi_enreg%comm_spin,ierr)
     gxfac(:,:,:,1)=gxfac(:,:,:,1)+gxfac_(:,:,:,ispinor_index)
     ABI_DEALLOCATE(gxfac_)
   end if

 end if !paw_opt

!Accumulate gxfac related to overlap (Sij) (PAW)
!------------------------------------------- ------------------------
 if (paw_opt==3.or.paw_opt==4) then                    ! Use Sij, overlap contribution
   gxfac_sij(1:cplex,1:nlmn,1:nincat,1:nspinor)=zero
   do ispinor=1,nspinor
     do ia=1,nincat
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         jjlmn=j0lmn+jlmn
         jlm=indlmn(4,jlmn)
         sijr=sij(jjlmn);gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
         gxfac_sij(1:cplex,jlmn,ia,ispinor)=gxfac_sij(1:cplex,jlmn,ia,ispinor)+sijr*gxj(1:cplex)
         do ilmn=1,jlmn-1
           ilm=indlmn(4,ilmn)
           if (ilm==jlm) then
             ijlmn=j0lmn+ilmn
             sijr=sij(ijlmn)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
             gxfac_sij(1:cplex,ilmn,ia,ispinor)=gxfac_sij(1:cplex,ilmn,ia,ispinor)+sijr*gxj(1:cplex)
             gxfac_sij(1:cplex,jlmn,ia,ispinor)=gxfac_sij(1:cplex,jlmn,ia,ispinor)+sijr*gxi(1:cplex)
           end if
         end do
       end do
     end do
   end do
 end if

!Accumulate dgxdtfac related to nonlocal operator (Norm-conserving)
!-------------------------------------------------------------------
 if (optder==1.and.paw_opt==0) then    ! Enl is E(Kleinman-Bylander)
   if (cplex_enl/=1    ) stop "opernlc_ylm: BUG - invalid cplex_enl=2 !"
   if (cplex_fac/=cplex) stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex !"
   do ispinor=1,nspinor
     ispinor_index = ispinor + shift
     do ia=1,nincat
       do ilmn=1,nlmn
         iln=indlmn(5,ilmn)
         enl_(1)=enl(iln,itypat,ispinor_index)
         do mu=1,ndgxdtfac
           dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)=enl_(1)*dgxdt(1:cplex,mu,ilmn,ia,ispinor)
         end do
       end do
     end do
   end do
 end if

!Accumulate dgxdtfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
 if (optder==1.and.(paw_opt==1.or.paw_opt==2.or.paw_opt==4)) then  ! Enl is psp strength Dij
   ABI_ALLOCATE(gxfj,(cplex,ndgxdtfac))
   dgxdtfac(1:cplex,1:ndgxdtfac,1:nlmn,1:nincat,1:nspinor)=zero
!  === Diagonal term(s) (up-up, down-down)
!  1-Enl is real
   if (cplex_enl==1) then
     do ispinor=1,nspinor
       ispinor_index=ispinor+shift
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl(jjlmn,index_enl,ispinor_index)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,ndgxdtfac
             gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
             dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1)=enl(ijlmn,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,ndgxdtfac
               gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
               dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)=dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
               dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
             end do
           end do
         end do
       end do
     end do
!    2-Enl is complex
   else
     if (cplex_fac/=cplex_enl) then
       stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex_enl !"
     end if
     if (nspinortot==1) then    !===== when nspinor=1, D_ij=D_ji
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl(2*ijlmn-1:2*jjlmn,index_enl,1)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,ndgxdtfac
             gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,1)
             dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)+enl_(1)*gxfj(1,mu)
             dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(2)*gxfj(1,mu)
             if (cplex==2) then
               dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)-enl_(2)*gxfj(2,mu)
               dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(1)*gxfj(2,mu)
             end if
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,ndgxdtfac
               gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,1)
               dgxdtfac(1,mu,ilmn,ia,1)=dgxdtfac(1,mu,ilmn,ia,1)+enl_(1)*gxfj(1,mu)
               dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
               dgxdtfac(2,mu,ilmn,ia,1)=dgxdtfac(2,mu,ilmn,ia,1)+enl_(2)*gxfj(1,mu)
               dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(2)*gxfi(1)
               if (cplex==2) then
                 dgxdtfac(1,mu,ilmn,ia,1)=dgxdtfac(1,mu,ilmn,ia,1)-enl_(2)*gxfj(2,mu)
                 dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)-enl_(2)*gxfi(2)
                 dgxdtfac(2,mu,ilmn,ia,1)=dgxdtfac(2,mu,ilmn,ia,1)+enl_(1)*gxfj(2,mu)
                 dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
               end if
             end do
           end do
         end do
       end do
     else                    !===== when nspinor=2, D_ij=D_ji^*
       do ispinor=1,nspinor
         ispinor_index = ispinor + shift
         do ia=1,nincat
           index_enl=atindx1(iatm+ia)
           do jlmn=1,nlmn
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             enl_(1)=enl(2*jjlmn-1,index_enl,ispinor_index)  ! enl_ii is real for up-up and dn-dn
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
             do mu=1,ndgxdtfac
               gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
               dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
             end do
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,ndgxdtfac
                 gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
                 dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
                 dgxdtfac(1,mu,jlmn,ia,ispinor)=dgxdtfac(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
                 dgxdtfac(2,mu,jlmn,ia,ispinor)=dgxdtfac(2,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(1)
                 if (cplex==2) then
                   dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                   dgxdtfac(1,mu,jlmn,ia,ispinor)=dgxdtfac(1,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(2)
                   dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
                   dgxdtfac(2,mu,jlmn,ia,ispinor)=dgxdtfac(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
                 end if
               end do
             end do
           end do
         end do
       end do
     end if !nspinortot
   end if !complex
!  === Off-diagonal term(s) (up-down, down-up)
!  --- No parallelization over spinors ---
   if (nspinortot==2.and.nspinor==nspinortot) then
     if (cplex_enl/=2) stop "opernlc_ylm: BUG - invalid cplex_enl/=2 !"
     if (cplex_fac/=2) stop "opernlc_ylm: BUG - invalid cplex_fac/=2 !"
     do ispinor=1,nspinor
       jspinor=3-ispinor
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor)
           do mu=1,ndgxdtfac
             gxfi(1:cplex)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
             gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,jspinor)
             dgxdtfac(1,mu,jlmn,ia,jspinor)=dgxdtfac(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
             dgxdtfac(2,mu,jlmn,ia,jspinor)=dgxdtfac(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
             if (cplex==2) then
               dgxdtfac(2,mu,jlmn,ia,jspinor)=dgxdtfac(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
               dgxdtfac(1,mu,jlmn,ia,jspinor)=dgxdtfac(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
             end if
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
             do mu=1,ndgxdtfac
               gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
               dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
               dgxdtfac(1,mu,jlmn,ia,jspinor)=dgxdtfac(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
               dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
               dgxdtfac(2,mu,jlmn,ia,jspinor)=dgxdtfac(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
               if (cplex==2) then
                 dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                 dgxdtfac(1,mu,jlmn,ia,jspinor)=dgxdtfac(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
                 dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
                 dgxdtfac(2,mu,jlmn,ia,jspinor)=dgxdtfac(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
               end if
             end do !mu
           end do !ilmn
         end do !jmln
       end do !ia
     end do !ispinor

!    --- Parallelization over spinors ---
   else if (nspinortot==2.and.nspinor/=nspinortot) then
     if (cplex_enl/=2) stop "opernlc_ylm: BUG - invalid cplex_enl/=2 !"
     if (cplex_fac/=2) stop "opernlc_ylm: BUG - invalid cplex_fac/=2 !"
     ABI_ALLOCATE(dgxdtfac_,(cplex_fac,ndgxdtfac,nlmn,nincat,nspinortot))
     dgxdtfac_(:,:,:,:,:)=zero
     ispinor_index=mpi_enreg%me_spin+1
     jspinor_index=3-ispinor_index
     if (ispinor_index==1) then
       ijspin=3;jispin=4
     else
       ijspin=4;jispin=3
     end if
     do ia=1,nincat
       index_enl=atindx1(iatm+ia)
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         do ilmn=1,nlmn
           i0lmn=ilmn*(ilmn-1)/2
           if (ilmn<=jlmn) then
             ijlmn=j0lmn+ilmn
             enl_(1)= enl(2*ijlmn-1,index_enl,ijspin)
             enl_(2)=-enl(2*ijlmn  ,index_enl,ijspin)
           else
             jilmn=i0lmn+jlmn
             enl_(1:2)=enl(2*jilmn-1:2*jilmn,index_enl, jispin)
           end if
           do mu=1,ndgxdtfac
             gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,1)
             dgxdtfac_(1,mu,jlmn,ia,jspinor_index)=dgxdtfac_(1,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(1)
             dgxdtfac_(2,mu,jlmn,ia,jspinor_index)=dgxdtfac_(2,mu,jlmn,ia,jspinor_index)+enl_(2)*gxfi(1)
             if (cplex==2) then
               dgxdtfac_(1,mu,jlmn,ia,jspinor_index)=dgxdtfac_(1,mu,jlmn,ia,jspinor_index)-enl_(2)*gxfi(2)
               dgxdtfac_(2,mu,jlmn,ia,jspinor_index)=dgxdtfac_(2,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(2)
             end if
           end do
         end do !ilmn
       end do !jlmn
     end do !iat
     call xsum_mpi(dgxdtfac_, mpi_enreg%comm_spin,ierr)
     dgxdtfac(:,:,:,:,1)=dgxdtfac(:,:,:,:,1)+ dgxdtfac_(:,:,:,:,ispinor_index)
     ABI_DEALLOCATE(dgxdtfac_)
   end if !nspinortot
   ABI_DEALLOCATE(gxfj)

 end if ! pawopt & optder

!Accumulate dgxdtfac related to overlap (Sij) (PAW)
!-------------------------------------------------------------------
 if (optder==1.and.(paw_opt==3.or.paw_opt==4)) then  ! Use Sij, overlap contribution
   ABI_ALLOCATE(gxfj,(cplex,ndgxdtfac))
   dgxdtfac_sij(1:cplex,1:ndgxdtfac,1:nlmn,1:nincat,1:nspinor)=zero
   do ispinor=1,nspinor
     do ia=1,nincat
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         jjlmn=j0lmn+jlmn
         sijr=sij(jjlmn)
         do mu=1,ndgxdtfac
           gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
           dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
         end do
         do ilmn=1,jlmn-1
           ijlmn=j0lmn+ilmn
           sijr=sij(ijlmn)
           do mu=1,ndgxdtfac
             gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
             dgxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
             dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfi(1:cplex)
           end do
         end do
       end do
     end do
   end do
   ABI_DEALLOCATE(gxfj)
 end if

end subroutine opernlc_ylm
!!***
