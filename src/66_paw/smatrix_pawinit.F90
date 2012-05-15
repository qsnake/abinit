!{\src2tex{textfont=tt}}
!!****f* ABINIT/smatrix_pawinit
!! NAME
!! smatrix_pawinit
!!
!! FUNCTION
!! Routine which computes paw part of the overlap used to compute LMWF wannier
!!  functions and berryphase
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (T Rangel,BAmadon,FJollet,PHermet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  g1(3)= reciprocal vector to put k1+b inside the BZ. bb=k2-k1=b-G1
!!  ("b" is the true b, so we have to correct bb with G1).
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  ikpt1(3)=cartesian coordinates of k1
!!  ikpt2(3)=cartesian coordinates of k2
!!  isppol  = spin polarization
!!  mband=maximum number of bands
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  seed_name= seed_name of files containing cg for all k-points to be used with MPI
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  cm2: Inside sphere part of the overlap needed for constructing wannier function
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!  The mpi part will work with mlwfovlp but not for berryphase_new
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,cprj_get,initylmr,leave_new,sbf8,simp_gen,wrtout
!!      xcomm_init,xmaster_init,xme_init,xproc_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine smatrix_pawinit(atindx1,cm2,cprj,ikpt1,ikpt2,isppol,&
& g1,gprimd,kpt,mband,mbandw,mkmem,mpi_enreg,&
& natom,nband,nkpt,nspinor,nsppol,ntypat,pawang,pawrad,pawtab,rprimd,&
& seed_name,typat,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_radmesh,   only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smatrix_pawinit'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_65_psp
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: ikpt1,ikpt2,isppol,mband,mbandw,mkmem,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: ntypat
 character(len=fnlen) ::  seed_name  !seed names of files containing cg info used in case of MPI
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang

!arrays
 integer,intent(in) :: atindx1(natom),g1(3),nband(nsppol*nkpt),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),kpt(3,nkpt),rprimd(3,3),xred(3,natom)
 real(dp),intent(inout) :: cm2(2,mbandw,mbandw)
 type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: dummy
 integer :: iatom,iband1,iband2,icg1,icg2,idx1,idx2,ii
 integer :: ilmn,ios,iunit,ir
 integer :: iorder_cprj,isel,ispinor,itypat,j0lmn,jj,jlmn,klm,klmn,kln,ll,lm0,lmax
 integer :: lmin,lmn_size,max_lmn,mesh_size,mm,nband_k
 integer :: master,nprocs,spaceComm,rank !for mpi
 real(dp) :: arg,bnorm,delta,intg,ppi,ppr,qijbtemp,qijtot,x1
 real(dp) :: x2,xsum,xtemp,xx,yy,zz
 character(len=500) :: message
 character(len=fnlen) ::  cprj_file  !file containing cg info used in case of MPI
 logical::lfile

!arrays
 integer,allocatable :: dimcprj(:)
 real(dp),parameter :: ili(7)=(/zero,-one,zero,one,zero,-one,zero/)
 real(dp),parameter :: ilr(7)=(/one,zero,-one,zero,one,zero,-one/)
 real(dp) :: bb(3),bb1(3),bbn(3),qijb(2),xcart(3,natom)
 real(dp),allocatable :: ff(:),j_bessel(:,:),ylmb(:),ylmrgr_dum(:,:,:)
 type(cprj_type),allocatable :: cprj_k1(:,:)
 type(cprj_type),allocatable :: cprj_k2(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (mkmem==0) then
   write(message,'(a,a,a,a,i3,a,a)')ch10,&
&   ' smatrix_pawinit: ERROR - ',ch10,&
&   '   Not compatible with mkmem=0 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!
!Allocate cprj_k1 and cprj_k2
!
 ABI_ALLOCATE(dimcprj,(natom))
 do iatom=1,natom
   dimcprj(iatom)=pawtab(typat(iatom))%lmn_size
 end do

 nband_k=nband(ikpt1)
 ABI_ALLOCATE(cprj_k1,(natom,nband_k*nspinor))
 call cprj_alloc(cprj_k1,0,dimcprj)

 nband_k=nband(ikpt2)
 ABI_ALLOCATE(cprj_k2,(natom,nband_k*nspinor))
 call cprj_alloc(cprj_k2,0,dimcprj)
 ABI_DEALLOCATE(dimcprj)



!mpi initialization
 call xcomm_init  (MPI_enreg,spaceComm)
 call xmaster_init(MPI_enreg,master   )
 call xme_init    (MPI_enreg,rank     )
 call xproc_init(MPI_enreg,nprocs)
 lfile=.false.
!
!write(std_out,*) "compute PAW overlap for k-points",ikpt1,ikpt2
 do iatom=1,natom
   xcart(:,iatom)=rprimd(:,1)*xred(1,iatom)+&
&   rprimd(:,2)*xred(2,iatom)+&
&   rprimd(:,3)*xred(3,iatom)
 end do

!
!Calculate indices icg1 and icg2
!
 icg1=0
 do ii=1,isppol
   ll=nkpt
   if(ii==isppol) ll=ikpt1-1
   do jj=1,ll
!    MPI: cycle over kpts not treated by this node
     if ( mpi_enreg%paral_compil_kpt==1) then
       if ( ABS(MPI_enreg%proc_distrb(jj,1,ii)-rank)/=0) CYCLE
!      write(std_out,'("kpt loop2: ikpt",i3," rank ",i3)') jj,rank
     end if
     icg1=icg1+nspinor*nband(jj+(ii-1)*nkpt)
   end do
 end do
 icg2=0
 do ii=1,isppol
   ll=nkpt
   if(isppol==ii) ll=ikpt2-1
   do jj=1,ll
!    MPI: cycle over kpts not treated by this node
     if ( mpi_enreg%paral_compil_kpt==1) then
       if (ABS(MPI_enreg%proc_distrb(jj,1,ii)-rank)/=0) CYCLE
!      write(std_out,'("kpt loop2: ikpt",i3," rank ",i3)') jj,rank
     end if
     icg2=icg2+nspinor*nband(jj+(ii-1)*nkpt)
   end do
 end do
!
!MPI: if ikpt2 not found in this processor then
!read info from an unformatted file
!
 if (mpi_enreg%paral_compil_kpt==1 ) then
   if (ABS(MPI_enreg%proc_distrb(ikpt2,1,isppol)-rank)/=0) then
     lfile=.true.
!    
!    get maximum of lmn_size
     max_lmn=0
     do itypat=1,ntypat
       lmn_size=pawtab(itypat)%lmn_size
       if(lmn_size>max_lmn) max_lmn=lmn_size
     end do
!    
!    get file name and open it
!    
     write(cprj_file,'(a,I5.5,".",I1)') trim(seed_name),ikpt2,isppol
     iunit=1000
!    write(std_out,*)'reading file',trim(cprj_file)
     open (unit=iunit, file=cprj_file,form='unformatted',status='old',iostat=ios)
     if(ios /= 0) then
       write(message,*) " smatrix_pawinit: file",trim(cprj_file), "not found"
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
!    
!    start reading
     do ii=1,mband*nspinor
       do iatom=1,natom
         itypat=typat(iatom)
         lmn_size=pawtab(itypat)%lmn_size
         do ilmn=1,lmn_size
           read(iunit)(cprj_k2(iatom,ii)%cp(jj,ilmn),jj=1,2)
         end do !ilmn
       end do
     end do
!    
!    close file
!    
     close (unit=iunit,iostat=ios)
     if(ios /= 0) then
       write(message,*) " smatrix_pawinit: error closing file ",trim(cprj_file)
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
!    
   end if
 end if !mpi

!Extract cprj_k1 and cprj_k2
!these contain the projectors cprj for just one k-point (ikpt1 or ikpt2)

!Extract cprj for k-point 1
 iorder_cprj=0 !do not change the ordering of cprj
 nband_k=nband(ikpt1)
 dummy=1000 !index of file not implemented here, mkmem==0 not implemented
 call cprj_get(atindx1,cprj_k1,cprj,natom,1,icg1,ikpt1,iorder_cprj,isppol,&
& mband,mkmem,mpi_enreg,natom,nband_k,nband_k,nspinor,nsppol,dummy)

!Extract cprj for k-point 2
 if( lfile .eqv. .false. ) then !if it was not already read above
   iorder_cprj=0 !do not change the ordering of cprj
   nband_k=nband(ikpt2)
   dummy=1000 !index of file not implemented here, mkmem==0 not implemented
   call cprj_get(atindx1,cprj_k2,cprj,natom,1,icg2,ikpt2,iorder_cprj,isppol,&
&   mband,mkmem,mpi_enreg,natom,nband_k,nband_k,nspinor,nsppol,dummy)
 end if

!DEBUG
!if(ikpt2==2) then
!if(ikpt1==1) then
!do iband1=1,mbandw
!do iatom=1,natom
!itypat=typat(atindx1(iatom))
!lmn_size=pawtab(itypat)%lmn_size
!do ilmn=1,1!lmn_size
!!write(500,'(a,i3,a,i3,a,i3,a,2f13.7)')'iband ',iband1,' iatom ',iatom,' ilmn ',ilmn,' cprj',cprj(iatom,iband1+icg1)%cp(:,ilmn)
!write(500,'(a,i3,a,i3,a,i3,a,2f13.7)')'iband ',iband1,' iatom ',atindx1(iatom),' ilmn ',ilmn,' cprj',cprj(atindx1(iatom),iband1+icg1)%cp(:,ilmn)
!write(500,'(a,i3,a,i3,a,i3,a,2f13.7)')'iband ',iband1,' iatom ',iatom,' ilmn ',ilmn,' cprj_k1 ',cprj_k1(iatom,iband1)%cp(:,ilmn)
!end do
!end do
!end do
!end if
!end if
!NDDEBUG

!!!!!!!!!!!!!!!!!
!--- Compute intermediate quantities: "b" vector=k2-k1 and its
!normalized value: bbn (and its norm: bnorm)
!compute also Ylm(b).
 ABI_ALLOCATE(ylmb,(pawang%l_size_max*pawang%l_size_max))
 ABI_ALLOCATE(ylmrgr_dum,(1,1,0))
 bb(:)=kpt(:,ikpt2)-kpt(:,ikpt1)+g1(:)
 bb1=bb
 xx=gprimd(1,1)*bb(1)+gprimd(1,2)*bb(2)+gprimd(1,3)*bb(3)
 yy=gprimd(2,1)*bb(1)+gprimd(2,2)*bb(2)+gprimd(2,3)*bb(3)
 zz=gprimd(3,1)*bb(1)+gprimd(3,2)*bb(2)+gprimd(3,3)*bb(3)
 bnorm=two_pi*dsqrt(xx**2+yy**2+zz**2)
 if(bnorm<tol8) then
!  write(std_out,*) "WARNING: bnorm=",bnorm
   bbn(:)=zero
 else
   xx=xx*two_pi
   yy=yy*two_pi
   zz=zz*two_pi
   bb(1)=xx
   bb(2)=yy
   bb(3)=zz
   bbn(1)=xx/bnorm
   bbn(2)=yy/bnorm
   bbn(3)=zz/bnorm
 end if

!debug  bbn=0
!debug  bnorm=0
!bbn has to ne normalized
 call initylmr(pawang%l_size_max,0,1,(/one/),1,bbn(:),ylmb(:),ylmrgr_dum)
!write(std_out,*) "ylmb(:)",ylmb(:)
!write(std_out,*) pawang%l_size_max
!write(std_out,*) "bbn",bbn(:)
!write(std_out,*) "xx,yy,zz",xx,yy,zz
!write(std_out,*) "bnorm",bnorm
 ABI_DEALLOCATE(ylmrgr_dum)

!------- First Compute Qij(b)-
 cm2=zero
 do iatom=1,natom
   itypat=typat(iatom)
   lmn_size=pawtab(itypat)%lmn_size
!  ---  en coordonnnes reelles cartesiennes (espace reel)
!  ---  first radial part(see pawinit)
   mesh_size=pawrad(itypat)%mesh_size
   ABI_ALLOCATE(j_bessel,(mesh_size,pawang%l_size_max))

!  ---  compute bessel function for (br) for all angular momenta necessary
!  ---  and for all value of r.
!  ---  they are needed for radial part
!  ---  of the integration => j_bessel(ir,:)
   do ir=1,mesh_size
     arg=bnorm*pawrad(itypat)%rad(ir)
     call sbf8(pawang%l_size_max,arg,j_bessel(ir,:))
   end do

!  do jlmn=1,pawang%l_size_max
!  write(665,*) "j_bessel",j_bessel(1:mesh_size,jlmn)
!  enddo
!  write(std_out,*) "bessel function computed"
!  ---  Compute \Sum b.R=xsum for future use
   xtemp=zero
   do mm=1,3
     xtemp=xtemp+xred(mm,iatom)*bb1(mm)
   end do
   xtemp=xtemp*two_pi
   xsum=zero
   do mm=1,3
     xsum=xsum+xcart(mm,iatom)*bbn(mm)*bnorm
   end do
!  write(std_out,*)'xsum',xsum,xtemp,lmn_size

!  ---  Loop on jlmn and ilmn
   qijtot=zero
   do jlmn=1,lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     do ilmn=1,jlmn

       klmn=j0lmn+ilmn
       klm=pawtab(itypat)%indklmn(1,klmn);kln=pawtab(itypat)%indklmn(2,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn);lmax=pawtab(itypat)%indklmn(4,klmn)
!      ---  Sum over angular momenta
!      ---  compute radial part integration for each angular momentum => intg
!      ---  (3j) symbols follows the rule: l belongs to abs(li-lj), li+lj.
       qijb=zero

       do ll=lmin,lmax,2
         lm0=ll*ll+ll+1
         ABI_ALLOCATE(ff,(mesh_size))
         ff(1:mesh_size)=(pawtab(itypat)%phiphj(1:mesh_size,kln)&
&         -pawtab(itypat)%tphitphj(1:mesh_size,kln))&
&         *j_bessel(1:mesh_size,ll+1)
         call simp_gen(intg,ff,pawrad(itypat))
         ABI_DEALLOCATE(ff)
         qijbtemp=zero
         do mm=-ll,ll
           isel=pawang%gntselect(lm0+mm,klm)
           if (isel>0) qijbtemp=qijbtemp&
&           +pawang%realgnt(isel)*ylmb(lm0+mm)
         end do ! mm
!        ---     compute angular part with a summation
!        ---     qijb =\sum_{lm} intg(lm)*qijbtemp
         qijb(1)=qijb(1) +intg*qijbtemp*ilr(ll+1)
         qijb(2)=qijb(2) +intg*qijbtemp*ili(ll+1)
!        if(ilmn==jlmn) write(std_out,*) "intg, qij",intg,qijbtemp
       end do ! ll

!      ---  Add exp(-i.b*R) for each atom.
       if(ilmn==jlmn) qijtot=qijtot+qijb(1)
!      if(ilmn==jlmn) write(std_out,*) "qijtot",qijtot
       x1=qijb(1)*dcos(-xsum)-qijb(2)*dsin(-xsum)
       x2=qijb(1)*dsin(-xsum)+qijb(2)*dcos(-xsum)
!      x1 x2 necessary to avoid changing qijb(1) before
!      computing qijb(2)
       qijb(1)=x1
       qijb(2)=x2 !
!      if(ilmn==jlmn) write(std_out,*) "qij",jlmn,ilmn,qijb(1),qijb(2)

       do iband1=1,mbandw ! limite inferieure a preciser
         do iband2=1,mbandw
           ppr=0.d0
           ppi=0.d0
           do ispinor=1,nspinor
             idx1=iband1*nspinor-(nspinor-ispinor)
             idx2=iband2*nspinor-(nspinor-ispinor) !to take into account spinors
!            write(std_out,*) "iband2",iband2
!            product of (a1+ia2)*(b1-ib2) (minus sign because conjugated)
             ppr=ppr+&
!            real part a_1*b_1+a_2*b_2
&             cprj_k1(iatom,idx1)%cp(1,ilmn)*cprj_k2(iatom,idx2)%cp(1,jlmn)+&
&             cprj_k1(iatom,idx1)%cp(2,ilmn)*cprj_k2(iatom,idx2)%cp(2,jlmn)+&
!            &     cprj(iatom,idx1+icg1)%cp(1,ilmn)*cprj(iatom,idx2+icg2)%cp(1,jlmn)+&
!            &     cprj(iatom,idx1+icg1)%cp(2,ilmn)*cprj(iatom,idx2+icg2)%cp(2,jlmn)+&
!            add term on the other triangle  of the matrix
!            qij is the same for this part because phi are real.
&             cprj_k1(iatom,idx1)%cp(1,jlmn)*cprj_k2(iatom,idx2)%cp(1,ilmn)+&
&             cprj_k1(iatom,idx1)%cp(2,jlmn)*cprj_k2(iatom,idx2)%cp(2,ilmn)
!            &     cprj(iatom,idx1+icg1)%cp(1,jlmn)*cprj(iatom,idx2+icg2)%cp(1,ilmn)+&
!            &     cprj(iatom,idx1+icg1)%cp(2,jlmn)*cprj(iatom,idx2+icg2)%cp(2,ilmn)
             ppi=ppi+&
!            imaginary part a_1*b_2-a_2*b_1
&             cprj_k1(iatom,idx1)%cp(1,ilmn)*cprj_k2(iatom,idx2)%cp(2,jlmn)-&
&             cprj_k1(iatom,idx1)%cp(2,ilmn)*cprj_k2(iatom,idx2)%cp(1,jlmn)+&
!            &     cprj(iatom,idx1+icg1)%cp(1,ilmn)*cprj(iatom,idx2+icg2)%cp(2,jlmn)-&
!            &     cprj(iatom,idx1+icg1)%cp(2,ilmn)*cprj(iatom,idx2+icg2)%cp(1,jlmn)+&
!            add term on the other triangle  of the matrix
&             cprj_k1(iatom,idx1)%cp(1,jlmn)*cprj_k2(iatom,idx2)%cp(2,ilmn)-&
&             cprj_k1(iatom,idx1)%cp(2,jlmn)*cprj_k2(iatom,idx2)%cp(1,ilmn)
!            &     cprj(iatom,idx1+icg1)%cp(1,jlmn)*cprj(iatom,idx2+icg2)%cp(2,ilmn)-&
!            &     cprj(iatom,idx1+icg1)%cp(2,jlmn)*cprj(iatom,idx2+icg2)%cp(1,ilmn)
           end do !ispinor
!          
!          delta: diagonal terms are counted twice ! so
!          we need a 0.5 factor for diagonal elements.
           delta=one
!          write(std_out,*) "ppr and ppi computed",ikpt1,ikpt2,iband1,iband2
           if(ilmn==jlmn) delta=half
           cm2(1,iband1,iband2)= cm2(1,iband1,iband2)+ &
&           (qijb(1)*ppr-qijb(2)*ppi)*delta
           cm2(2,iband1,iband2)= cm2(2,iband1,iband2)+ &
&           (qijb(2)*ppr+qijb(1)*ppi)*delta
         end do ! iband2
       end do ! iband1

     end do ! ilmn
   end do ! jlmn
!  write(std_out,*) "final qijtot",qijtot
   ABI_DEALLOCATE(j_bessel)
 end do ! iatom
 ABI_DEALLOCATE(ylmb)
 call cprj_free(cprj_k1)
 call cprj_free(cprj_k2)
 ABI_DEALLOCATE(cprj_k1)
 ABI_DEALLOCATE(cprj_k2)

 DBG_EXIT("COLL")

 end subroutine    smatrix_pawinit
!!***

