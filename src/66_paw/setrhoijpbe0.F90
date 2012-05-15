!{\src2tex{textfont=tt}}
!!****f* ABINIT/setrhoijpbe0
!! NAME
!! setrhoijpbe0
!!
!! FUNCTION
!! PAW local exact exchange only:
!! Impose value of rhoij for f electrons using an auxiliairy file
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  indlmn(6,i,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  istep=index of the number of steps in the routine scfcv
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all types of psps
!!  natom=number of atoms in cell.
!!  natom=number of atoms in cell
!!  ntypat=number of types of atoms in unit cell.
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type integer for each atom in cell
!!
!! SIDE EFFECTS
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! NOTES
!!  Only valid for f electrons !!!
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      leave_new,print_ij,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setrhoijpbe0(dtset,indlmn,initialized,istep,istep_mix,lmnmax,natom,ntypat,pawrhoij,pawtab,typat)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setrhoijpbe0'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: initialized,istep,lmnmax,natom,ntypat
 integer,intent(inout) :: istep_mix
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),typat(natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: ll=3
 integer :: iatom,ii,ios,irhoij,ispden,itypat,jj,klmn,nselect,nstep1,nstep1_abs,rhoijshft,rhoijsz
 logical :: test0
 character(len=9),parameter :: filnam='rhoijpbe0'
 character(len=9),parameter :: dspin(6)=(/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
 character(len=500) :: strg, message
!arrays
 real(dp),allocatable :: rhoijtmp1(:,:),rhoijtmp(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Test existence of file and open it
 inquire(file=filnam,iostat=ios,exist=test0)
 if(.not.test0) return

!Test if exact-exch. is on f electrons
 test0=.false.
 do itypat=1,ntypat
   if (pawtab(itypat)%useexexch>0.and.pawtab(itypat)%lexexch/=ll) test0=.true.
 end do
 if (test0) then
   write(message, '(4a,i1,a)' )ch10,&
&   ' setrhoijpbe0 : ERROR -',ch10,&
&   '  Local exact exchange: occ. matrix can only be imposed for l=',ll,' !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Open file
 open(77,file=filnam,form='formatted',iostat=ios)

!Read step number and eventually exit
 nstep1=0;test0=.false.
 do while (.not.test0)
   read(77,'(A)') strg
   test0=(strg(1:1)/="#")
   if (test0) read(unit=strg,fmt=*) nstep1
 end do
 nstep1_abs=abs(nstep1)
 if (nstep1_abs==0.or.istep>nstep1_abs.or.(nstep1>0.and.initialized/=0)) then
   close(77)
!  Reinitalize mixing when rhoij is allowed to change; for experimental purpose...
   if (dtset%userib==1234.and.istep==1+nstep1_abs.and.(nstep1<0.or.initialized==0)) istep_mix=1
   return
 end if

!Loop on atoms
 do iatom=1,natom
   itypat=typat(iatom)
   if (pawtab(itypat)%useexexch>0) then

!    Set sizes depending on ll
     rhoijsz=4*ll+2
     rhoijshft=2*ll*ll

!    Uncompress rhoij
     ABI_ALLOCATE(rhoijtmp,(pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
     do ispden=1,pawrhoij(iatom)%nspden
       rhoijtmp=zero
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         rhoijtmp(klmn,ispden)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
       end do
     end do

!    Read rhoij from file
     ABI_ALLOCATE(rhoijtmp1,(rhoijsz,rhoijsz))
     do ispden=1,pawrhoij(iatom)%nspden
       do ii=1,rhoijsz
         test0=.false.
         do while (.not.test0)
           read(77,'(A)') strg
           test0=(strg(1:1)/="#")
           if (test0)  read(unit=strg,fmt=*) (rhoijtmp1(ii,jj), jj=1,rhoijsz)
         end do
       end do

!      Impose rhoij
       do jj=1,rhoijsz
         do ii=1,jj
           rhoijtmp((jj+rhoijshft)*((jj+rhoijshft)-1)/2+ii+rhoijshft,ispden)=rhoijtmp1(ii,jj)
         end do
       end do

     end do
     ABI_DEALLOCATE(rhoijtmp1)

!    Compress rhoij
     nselect=0
     do klmn=1,pawrhoij(iatom)%lmn2_size
       if (any(abs(rhoijtmp(klmn,:))>tol10)) then
         nselect=nselect+1
         do ispden=1,pawrhoij(iatom)%nspden
           pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
         end do
         pawrhoij(iatom)%rhoijselect(nselect)=klmn
       end if
     end do
     pawrhoij(iatom)%nrhoijsel=nselect
     ABI_DEALLOCATE(rhoijtmp)

!    Print new rhoij
     do ispden=1,pawrhoij(iatom)%nspden
       write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&       ' == Imposed occupation matrix'
       if (pawrhoij(iatom)%nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
       if (pawrhoij(iatom)%nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
       if (pawrhoij(iatom)%nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&       trim(dspin(ispden+2*(pawrhoij(iatom)%nspden/4)))," =="
       call wrtout(std_out,message,'COLL')
       call print_ij(pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&       pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,1,ll,&
&       indlmn(1,1:pawtab(itypat)%lmn_size,itypat),&
&       1,-1,pawrhoij(iatom)%rhoijselect(:),-1.d0,1)
     end do

!    End loop on atoms
   end if
 end do

!Close file
 close (77)

 DBG_EXIT("COLL")

end subroutine setrhoijpbe0
!!***
