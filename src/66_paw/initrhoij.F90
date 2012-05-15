!{\src2tex{textfont=tt}}
!!****f* ABINIT/initrhoij
!! NAME
!! initrhoij
!!
!! FUNCTION
!! Initialize PAW rhoij occupancies (in packed storage)
!! from atomic ones
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex=1 if rhoij are REAL, 2 if they are complex
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  lexexch(ntypat)=l on which local exact-exchange is applied for a given type of atom
!!  lmnmax=max number of (l,m,n) comp. over all type of psps
!!  lpawu(ntypat)=l on which U is applied for a given type of atom (PAW+U)
!!  mpi_enreg=informations about MPI parallelization:
!!  natom=number of atoms
!!  natom_paw=size of PAW arrays, # of atoms on current proc.
!!  nspden=number of spin-density components FOR RHOIJ
!!  nspinor=number of spinorial components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of atom types
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!                                     (containing initial rhoij)
!!  spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!!  typat(natom)=type of each atom
!!  === Optional arguments
!!    ngrhoij=number of gradients to be allocated (OPTIONAL, default=0)
!!    nlmnmix=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
!!    use_rhoij_=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
!!    use_rhoijres=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)

!!
!! OUTPUT
!!  pawrhoij(natom) <type(pawrhoij_type)>=rhoij quantities for each atom
!!                                        in packed storage
!!
!! PARENTS
!!      gstate,respfn,setup_positron
!!
!! CHILDREN
!!      leave_new,rhoij_alloc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initrhoij(cplex,indlmn,lexexch,lmnmax,lpawu,mpi_enreg,natom,natom_paw,&
&                    nspden,nspinor,nsppol,ntypat,pawrhoij,pawtab,spinat,typat,&
&                    ngrhoij,nlmnmix,use_rhoij_,use_rhoijres) ! Optional arguments

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initrhoij'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,lmnmax,natom,natom_paw,nspden,nspinor,nsppol,ntypat
 integer,intent(in),optional :: ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
 character(len=500) :: message
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),lexexch(ntypat),lpawu(ntypat)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: spinat(3,natom)
 type(pawrhoij_type),intent(out) :: pawrhoij(natom_paw)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!Arrays
!scalars
 integer :: iatom,iatom_rhoij,ilmn,ispden,itypat,j0lmn,jl,jlmn,jspden,klmn,klmn1,ngrhoij0,nlmnmix0
 integer :: nselect,nselect1,use_rhoij_0,use_rhoijres0
 real(dp) :: ratio,ro,roshift,zratio,zz
 logical :: test_exexch,test_pawu
!arrays
 integer,allocatable :: nlmn(:)

!************************************************************************

 DBG_ENTER("COLL")

 if (mpi_enreg%nproc_atom>1) then
   if (natom_paw/=mpi_enreg%natom) then
     MSG_BUG("natom_paw not equal to mpi_enreg%natom !")
   end if
 end if

!PAW+U and local exact-exchange restriction
 do itypat=1,ntypat
   if (lpawu(itypat)/=lexexch(itypat).and.&
&   lpawu(itypat)/=-1.and.lexexch(itypat)/=-1) then
     write(message, '(4a)' ) ch10,' initrhoij: ERROR - ',&
&     ch10,'  lpawu must be equal to lexexch !'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
 end do

 ratio=one;if (nspden==2) ratio=half

 ABI_ALLOCATE(nlmn,(ntypat))
 do itypat=1,ntypat
   nlmn(itypat)=pawtab(itypat)%lmn_size
 end do
 ngrhoij0=0;if (present(ngrhoij)) ngrhoij0=ngrhoij
 nlmnmix0=0;if (present(nlmnmix)) nlmnmix0=nlmnmix
 use_rhoij_0=0;if (present(use_rhoij_)) use_rhoij_0=use_rhoij_
 use_rhoijres0=0;if (present(use_rhoijres)) use_rhoijres0=use_rhoijres
 call rhoij_alloc(cplex,nlmn,nspden,nspinor,nsppol,pawrhoij,typat,mpi_enreg=mpi_enreg,&
& ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0)
 ABI_DEALLOCATE(nlmn)

 do iatom_rhoij=1,natom_paw
   iatom=iatom_rhoij;if (mpi_enreg%nproc_atom>1) iatom=mpi_enreg%atom_indx(iatom_rhoij)
   itypat=typat(iatom)
   nselect=0

!  Determine Z (trace of rhoij0 or part of it)
   zz=zero
   do jlmn=1,pawtab(itypat)%lmn_size
     jl=indlmn(1,jlmn,itypat)
     j0lmn=jlmn*(jlmn-1)/2
     test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
     test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       if ((ilmn==jlmn).and.test_pawu.and.test_exexch) &
&       zz=zz+pawtab(itypat)%rhoij0(klmn)
     end do
   end do

!  Compute rhoij from tabulated value and magnetization
   do ispden=1,nspden

     zratio=zero
     roshift=one
     ratio=one
     if (nspden==2) then
       ratio=half
       if ((spinat(3,iatom)>zero.and.ispden==1).or.&
&       (spinat(3,iatom)<zero.and.ispden==2)) then
         zratio=two*abs(spinat(3,iatom))/zz
       end if
     else if (nspden==4.and.ispden>=2) then
       roshift=zero
       zratio=spinat(ispden-1,iatom)/zz
     end if

     nselect=0;nselect1=1-cplex
     do jlmn=1,pawtab(itypat)%lmn_size
       jl=indlmn(1,jlmn,itypat)
       j0lmn=jlmn*(jlmn-1)/2
       test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
       test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         ro=pawtab(itypat)%rhoij0(klmn)
         if ((ilmn==jlmn).and.test_pawu.and.test_exexch) then
           ro=ro*ratio*(roshift+zratio)
         else
           ro=ro*ratio*roshift
         end if

         klmn1=cplex*(klmn-1)+1
         if (abs(ro)>tol10) then
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=ro
         else
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=zero
         end if

         if (ispden==nspden) then
           if (any(abs(pawrhoij(iatom_rhoij)%rhoijp(klmn1,:))>tol10)) then
             nselect=nselect+1;nselect1=nselect1+cplex
             pawrhoij(iatom_rhoij)%rhoijselect(nselect)=klmn
             do jspden=1,nspden
               pawrhoij(iatom_rhoij)%rhoijp(nselect1,jspden)=pawrhoij(iatom_rhoij)%rhoijp(klmn1,jspden)
             end do
           end if
         end if

       end do
     end do

   end do
   pawrhoij(iatom_rhoij)%nrhoijsel=nselect

 end do ! iatom_rhoij

 DBG_EXIT("COLL")

end subroutine initrhoij
!!***
