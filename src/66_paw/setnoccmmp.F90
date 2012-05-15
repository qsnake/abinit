!{\src2tex{textfont=tt}}
!!****f* ABINIT/setnoccmmp
!! NAME
!! setnoccmmp
!!
!! FUNCTION
!! PAW+U only:
!!   Compute density matrix nocc_{m,m_prime}
!!   or
!!   Impose value of density matrix using dmatpawu input array, then symetrize it.
!!
!! noccmmp^{\sigma}_{m,m'}=\sum_{ni,nj}[\rho^{\sigma}_{ni,nj}*phiphjint_{ni,nj}]
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (BA,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  compute_dmat= flag: if 1, nocc_{m,mp} is computed
!!  dimdmat=first dimension of dmatpawu array
!!  dmatpawu(dimdmat,dimdmat,nsppol*nspinor,natpawu)=input density matrix to be copied into noccmpp
!!  dmatudiag= flag controlling the use of diagonalization:
!!             0: no diagonalization of nocc_{m,mp}
!!             1: diagonalized nocc_{m,mp} matrix is printed
!!             2: dmatpawu matrix is expressed in the basis where nocc_(m,mp} is diagonal
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  natom=number of atoms in cell
!!  natpawu=number of atoms on which PAW+U is applied
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of independant spin components
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of atom types
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  spinat(3,matom)=initial spin of each atom, in unit of hbar/2
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  impose_dmat= flag: if 1, nocc_{m,mp} is replaced by dmatpawu
!!  useexexch=1 if local-exact-exchange is activated
!!  usepawu=1 if PAW+U is activated
!!
!! OUTPUT
!!   paw_ij(natom)%noccmmp(paw_ij(iatom)%cplex_dij,2*pawtab(itypat)%lpawu+1,2*pawtab(itypat)%lpawu+1,nspden)=density matrix
!!
!! NOTES
!! For non-collinear magnetism,
!! - nocc_{m,mp} is computed as:noccmmp(:,:,:,1)=   n{m,mp}
!!                              noccmmp(:,:,:,2)=   m_x{m,mp}
!!                              noccmmp(:,:,:,3)=   m_y{m,mp}
!!                              noccmmp(:,:,:,4)=   m_z{m,mp}
!! - but nocc_{m,mp} is stored as: noccmmp(:,:,:,1)=   n^{up,up}_{m,mp}
!!                                 noccmmp(:,:,:,2)=   n^{dn,dn}_{m,mp}
!!                                 noccmmp(:,:,:,3)=   n^{up,dn}_{m,mp}
!!                                 noccmmp(:,:,:,4)=   n^{dn,up}_{m,mp}
!!   We choose to have noccmmp complex when ndij=4 (ie nspinor=2)
!!    If ndij=4 and pawspnorb=0, one could keep noccmmp real
!!    with the n11, n22, Re(n12), Im(n21) representation, but it would
!!    less clear to change the representation when pawspnorb is activated.
!!   If ndij=4, nocc_{m,mp} is transformed to the Ylm basis
!!    and then to the J, M_J basis (if cplex_dij==2)
!!
!!  Note that n_{m,mp}=<mp|hat(n)|m> because rhoij=<p_j|...|p_i>
!!
!! PARENTS
!!      afterscfloop,pawdenpot,pawprt,scfcv
!!
!! CHILDREN
!!      dgemm,dsyev,leave_new,mat_mlms2jmj,mat_slm2ylm,wrtout,zgemm,zheev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setnoccmmp(compute_dmat,dimdmat,dmatpawu,dmatudiag,impose_dmat,indsym,natom,natpawu,&
&                     nspinor,nsppol,nsym,ntypat,paw_ij,pawang,pawprtvol,pawrhoij,pawtab,&
&                     spinat,symafm,useexexch,usepawu)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setnoccmmp'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: compute_dmat,dimdmat,dmatudiag,impose_dmat,natom,natpawu
 integer,intent(in) :: nspinor,nsppol,nsym,ntypat,useexexch,usepawu
 type(pawang_type),intent(in) :: pawang
 integer,intent(in) :: pawprtvol
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symafm(nsym)
 real(dp),intent(in) :: dmatpawu(dimdmat,dimdmat,nspinor*nsppol,natpawu*impose_dmat)
 real(dp),intent(in) :: spinat(3,natom)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: limp=0 ! could become an input variable
 integer :: at_indx,cplex_dij,dmatudiag_loc,iafm,iatom,iatpawu,icount,ilm,im1,im2,in1,in2,info,iplex,irot,ispden
 integer :: irhoij,itypat,jlm,jrhoij,jspden,klmn,kspden,lcur,ldim,lmax,lmin,lpawu,lwork,ndij,nmat,nspden,nsploop
 real(dp) :: factafm,mnorm,mx,my,mz,ntot,nup,ndn,snorm,sx,sy,szp,szm
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used with non-collinear magnetism
 logical :: antiferro,noccsym_error,use_afm
 character(len=500) :: message
 real(dp),parameter :: invsqrt2=one/sqrt2

! *************************************************************************
!arrays
 integer :: nsym_used(2)
 real(dp) :: sumocc(2)
 real(dp),allocatable :: eig(:),hdp(:,:,:),hdp2(:,:),noccmmptemp(:,:,:,:),noccmmp_tmp(:,:,:,:)
 real(dp),allocatable :: rwork(:),ro(:),noccmmp2(:,:,:,:),nocctot2(:)
 complex(dpc),allocatable :: zhdp(:,:),zhdp2(:,:),znoccmmp_tmp(:,:),zwork(:)
 character(len=9),parameter :: dspin(6)=  (/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
 character(len=9),parameter :: dspinc(6)= (/"up       ","down     ","up-up    ","down-down","up-dn    ","dn-up    "/)
 character(len=9),parameter :: dspinc2(6)=(/"up       ","down     ","dn-dn    ","up-up    ","dn-up    ","up-dn    "/)
 character(len=9),parameter :: dspinm(6)= (/"dn       ","up i     ","n        ","mx       ","my       ","mz       "/)
 complex(dpc),allocatable :: noccmmp_ylm(:,:,:),noccmmp_jmj(:,:),noccmmp_slm(:,:,:)
 type(coeff4_type),allocatable :: tmp_noccmmp(:)

!*********************************************************************

 DBG_ENTER("COLL")

!Tests
 if (nsppol/=paw_ij(1)%nsppol) then
   write(message, '(4a)' ) ch10,&
&   ' setnoccmmp: BUG - ',ch10,&
&   '   inconsistent values for nsppol !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (compute_dmat>0) then
   if (pawrhoij(1)%nspden/=paw_ij(1)%nspden.and.&
&   pawrhoij(1)%nspden/=4.and.paw_ij(1)%nspden/=1) then
     write(message, '(4a)' ) ch10,&
&     ' setnoccmmp: BUG - ',ch10,&
&     '   inconsistent values for nspden !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
 if (usepawu>0.and.useexexch>0) then
   write(message, '(4a)' ) ch10,&
&   ' setnoccmmp: BUG - ',ch10,&
&   '  usepawu>0 and useexexch>0 not allowed !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (impose_dmat/=0.and.dimdmat==0) then
   write(message, '(4a)' ) ch10,&
&   ' setnoccmmp: BUG - ',ch10,&
&   '   dmatpawu must be allocated when impose_dmat/=0 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (usepawu>0.and.compute_dmat/=0.and.impose_dmat/=0.and.pawang%nsym==0) then
   write(message, '(4a)' ) ch10,&
&   ' setnoccmmp: BUG - ',ch10,&
&   '   pawang%zarot must be allocated !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Some inits
 if (usepawu==0.and.useexexch==0) return
 nspden=paw_ij(1)%nspden
 ndij=paw_ij(1)%ndij
 cplex_dij=paw_ij(1)%cplex_dij
 antiferro=(nspden==2.and.nsppol==1)
 use_afm=((antiferro).or.((nspden==4).and.afm_noncoll))
 dmatudiag_loc=dmatudiag
 if (dmatudiag==2.and.(dimdmat==0.or.impose_dmat==0)) dmatudiag_loc=1

!If needed, store dmatpu in suitable format in tmp_noccmmp
 if (usepawu>0.and.impose_dmat/=0) then
   iatpawu=0
   ABI_ALLOCATE(tmp_noccmmp,(natom))
   do iatom=1,natom
     itypat=pawrhoij(iatom)%itypat
     lpawu=pawtab(itypat)%lpawu
     if (lpawu/=-1) then
       iatpawu=iatpawu+1
       if (ndij/=4) then
         ABI_ALLOCATE(tmp_noccmmp(iatom)%value,(cplex_dij,2*lpawu+1,2*lpawu+1,nsppol))
         tmp_noccmmp(iatom)%value(1,1:2*lpawu+1,1:2*lpawu+1,1:nsppol)=&
&         dmatpawu(1:2*lpawu+1,1:2*lpawu+1,1:nsppol,iatpawu)
       else
         ABI_ALLOCATE(tmp_noccmmp(iatom)%value,(cplex_dij,2*lpawu+1,2*lpawu+1,nspden))
         tmp_noccmmp(iatom)%value=zero
         if(limp==0) then ! default reading
           snorm=sqrt(spinat(1,natom)**2+spinat(1,iatom)**2+spinat(3,iatom)**2)
           if (snorm>tol12) then
             sx=half*spinat(1,iatom)/snorm
             sy=half*spinat(2,iatom)/snorm
             szp=half*(one+spinat(3,iatom)/snorm)
             szm=half*(one-spinat(3,iatom)/snorm)
           else
             sx=zero;sy=zero
             szp=one;szm=zero
           end if
           do im2=1,2*lpawu+1
             do im1=1,2*lpawu+1
               nup=dmatpawu(im1,im2,1,iatpawu);ndn=dmatpawu(im1,im2,2,iatpawu)
               tmp_noccmmp(iatom)%value(1,im1,im2,1)=nup*szp+ndn*szm
               tmp_noccmmp(iatom)%value(1,im1,im2,2)=nup*szm+ndn*szp
               tmp_noccmmp(iatom)%value(1,im1,im2,3)=(nup-ndn)*sx
               tmp_noccmmp(iatom)%value(1,im1,im2,4)=(ndn-nup)*sy
             end do
           end do
         else if(limp>=1) then
           ABI_ALLOCATE(noccmmp_ylm,(2*lpawu+1,2*lpawu+1,ndij))
           noccmmp_ylm=czero
           ABI_ALLOCATE(noccmmp_slm,(2*lpawu+1,2*lpawu+1,ndij))
           noccmmp_slm=czero
           ABI_ALLOCATE(noccmmp_jmj,(2*(2*lpawu+1),2*(2*lpawu+1)))
           noccmmp_jmj=czero
           if(limp==1) then ! read input matrix in J,M_J basis (l-1/2, then l+1/2)
             noccmmp_jmj=czero
             do im1=1,2*lpawu+1
               noccmmp_jmj(im1,im1)=cmplx(dmatpawu(im1,im1,1,iatpawu),zero)
               noccmmp_jmj(im1+lpawu,im1+lpawu)=cmplx(dmatpawu(im1+lpawu,im1+lpawu,2,iatpawu),zero)
             end do
             write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&             ' == Imposed occupation matrix (in the J M_J basis: L-1/2 and L+1/2 states)'
             call mat_mlms2jmj(lpawu,noccmmp_ylm,noccmmp_jmj,ndij,2,2,pawprtvol) !  optspin=1: up spin are first
           end if
           if(limp==2) then ! read input matrix in Ylm basis
             noccmmp_ylm=czero
             do im1=1,2*lpawu+1
               noccmmp_ylm(im1,im1,1)=cmplx(dmatpawu(im1,im1,1,iatpawu),zero)
               noccmmp_ylm(im1,im1,2)=cmplx(dmatpawu(im1,im1,2,iatpawu),zero)
             end do
             write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&             ' == Imposed occupation matrix (in the Ylm basis), for dn and up spin'
           end if
           call mat_slm2ylm(lpawu,noccmmp_ylm,noccmmp_slm,ndij,2,2,pawprtvol) ! optspin=1 because up spin are first
!          interchange upup and dndn
           if(limp>=1) then
             tmp_noccmmp(iatom)%value(1,:,:,1)=real(noccmmp_slm(:,:,2))
             tmp_noccmmp(iatom)%value(2,:,:,1)=imag(noccmmp_slm(:,:,2))
             tmp_noccmmp(iatom)%value(1,:,:,2)=real(noccmmp_slm(:,:,1))
             tmp_noccmmp(iatom)%value(2,:,:,2)=imag(noccmmp_slm(:,:,1))
             tmp_noccmmp(iatom)%value(1,:,:,3)=real(noccmmp_slm(:,:,4))
             tmp_noccmmp(iatom)%value(2,:,:,3)=imag(noccmmp_slm(:,:,4))
             tmp_noccmmp(iatom)%value(1,:,:,4)=real(noccmmp_slm(:,:,3))
             tmp_noccmmp(iatom)%value(2,:,:,4)=imag(noccmmp_slm(:,:,3))
           end if
           if(abs(pawprtvol)>2) then
             write(message, '(2a)' ) ch10,&
&             " Check Imposed density matrix in different basis"
             call wrtout(std_out,message,'COLL')
             call mat_slm2ylm(lpawu,noccmmp_slm,noccmmp_ylm,ndij,1,2,pawprtvol) ! optspin=1 because up spin are first
             call mat_mlms2jmj(lpawu,noccmmp_ylm,noccmmp_jmj,ndij,1,2,pawprtvol) !  optspin=1: up spin are first
           end if
           ABI_DEALLOCATE(noccmmp_ylm)
           ABI_DEALLOCATE(noccmmp_jmj)
           ABI_DEALLOCATE(noccmmp_slm)
         end if
       end if
     end if
   end do
 end if  ! impose_dmat/=0

!Print message
 if (usepawu>0.and.impose_dmat/=0) then
   if (dmatudiag_loc/=2) then
     write(message,'(6a)') ch10,'Occupation matrix for correlated orbitals is kept constant',ch10,&
&     'and equal to dmatpawu from input file !',ch10,&
&     '----------------------------------------------------------'
   else
     write(message,'(6a)') ch10,'Occupation matrix for correlated orbitals is imposed',ch10,&
&     'and equal to dmatpawu in the diagonal basis !',ch10,&
&     '----------------------------------------------------------'
   end if
   call wrtout(std_out,message,'COLL')
 end if

 if (usepawu>0.and.dmatudiag_loc/=0) then
   write(message,'(4a)') ch10,'Diagonalized occupation matrix "noccmmp" is printed !',ch10,&
&   '-------------------------------------------------------------'
   call wrtout(std_out,message,'COLL')
 end if

!Loops over atoms
 do iatom=1,natom
   itypat=pawrhoij(iatom)%itypat
   if (useexexch>0) then
     lcur=pawtab(itypat)%lexexch
   else if (usepawu>0) then
     lcur=pawtab(itypat)%lpawu
   end if
   if (lcur/=-1) then

!    ########################################################################################
!    # Compute nocc_mmp
!    ########################################################################################
     if ((usepawu>0.and.compute_dmat/=0).or.useexexch>0) then

       paw_ij(iatom)%noccmmp(:,:,:,:)=zero

!      Loop over spin components
       ABI_ALLOCATE(noccmmptemp,(cplex_dij,2*lcur+1,2*lcur+1,ndij))
       noccmmptemp(:,:,:,:)=zero
       if(ndij==4)  then
         ABI_ALLOCATE(noccmmp2,(cplex_dij,2*lcur+1,2*lcur+1,ndij))
       end if
       if(ndij==4)  then
         ABI_ALLOCATE(nocctot2,(ndij))
       end if
       do ispden=1,ndij
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           im1=pawtab(itypat)%klmntomn(1,klmn)
           im2=pawtab(itypat)%klmntomn(2,klmn)
           in1=pawtab(itypat)%klmntomn(3,klmn)
           in2=pawtab(itypat)%klmntomn(4,klmn)
           lmin=pawtab(itypat)%indklmn(3,klmn)
           lmax=pawtab(itypat)%indklmn(4,klmn)
           ABI_ALLOCATE(ro,(cplex_dij))
           if (ndij==1) then
             ro(1)=half*pawrhoij(iatom)%rhoijp(jrhoij,1)
           else if (ndij==2) then
             ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)
           else  ! ndij==4
!            Non-collinear magnetism: transfer rhoij to ro_c (keep n, m storage because
!            it is easier for the computation of noccmmp from rhoij)
!            cplex_dij has to be used here, because it is the dimension of rhoijp
             ro(1:cplex_dij)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij-1+cplex_dij,ispden)
           end if
           if(lmin==0.and.lmax==2*lcur) then
             icount=in1+(in2*(in2-1))/2
             if(pawtab(itypat)%ij_proj<icount)  then
               write(message, '(4a)' ) ch10,&
&               '  setnoccmmp : BUG -',ch10,&
&               '  PAW+U: Problem in the loop for calculating noccmmp !',ch10
               call wrtout(std_out,message,'COLL')
               call leave_new('COLL')
             end if
             if(in1/=in2) then
               if(im2<=im1) then
                 noccmmptemp(:,im1,im2,ispden)=noccmmptemp(:,im1,im2,ispden)+ro(:)*pawtab(itypat)%phiphjint(icount)
               end if
             end if
             if(im2>=im1) then
               paw_ij(iatom)%noccmmp(:,im1,im2,ispden)=paw_ij(iatom)%noccmmp(:,im1,im2,ispden) &
&               +ro(:)*pawtab(itypat)%phiphjint(icount)
             end if
           end if
           jrhoij=jrhoij+pawrhoij(iatom)%cplex
           ABI_DEALLOCATE(ro)
         end do ! irhoij
         do im2=1,2*lcur+1
           do im1=1,im2
             paw_ij(iatom)%noccmmp(1,im1,im2,ispden)=paw_ij(iatom)%noccmmp(1,im1,im2,ispden) &
&             +noccmmptemp(1,im2,im1,ispden)
             if(cplex_dij==2) paw_ij(iatom)%noccmmp(2,im1,im2,ispden)=paw_ij(iatom)%noccmmp(2,im1,im2,ispden) &
&             -noccmmptemp(2,im2,im1,ispden)
           end do
         end do
         do im1=1,2*lcur+1
           do im2=1,im1
             paw_ij(iatom)%noccmmp(1,im1,im2,ispden)=paw_ij(iatom)%noccmmp(1,im2,im1,ispden)
             if(cplex_dij==2) paw_ij(iatom)%noccmmp(2,im1,im2,ispden)=-paw_ij(iatom)%noccmmp(2,im2,im1,ispden)
           end do
         end do
       end do ! ispden
       ABI_DEALLOCATE(noccmmptemp)
!      Compute noccmmp2, occupation matrix in the spin basis (upup, dndn, updn, dnup)
       if(ndij==4) then
         noccmmp2(:,:,:,:)=zero
         do im1=1,2*lcur+1
           do im2=1,2*lcur+1
             noccmmp2(1,im1,im2,1)=half*(paw_ij(iatom)%noccmmp(1,im1,im2,1)+paw_ij(iatom)%noccmmp(1,im1,im2,4))
             noccmmp2(2,im1,im2,1)=half*(paw_ij(iatom)%noccmmp(2,im1,im2,1)+paw_ij(iatom)%noccmmp(2,im1,im2,4))
             noccmmp2(1,im1,im2,2)=half*(paw_ij(iatom)%noccmmp(1,im1,im2,1)-paw_ij(iatom)%noccmmp(1,im1,im2,4))
             noccmmp2(2,im1,im2,2)=half*(paw_ij(iatom)%noccmmp(2,im1,im2,1)-paw_ij(iatom)%noccmmp(2,im1,im2,4))
             noccmmp2(1,im1,im2,3)=half*(paw_ij(iatom)%noccmmp(1,im1,im2,2)+paw_ij(iatom)%noccmmp(2,im1,im2,3))
             noccmmp2(2,im1,im2,3)=half*(paw_ij(iatom)%noccmmp(2,im1,im2,2)-paw_ij(iatom)%noccmmp(1,im1,im2,3))
             noccmmp2(1,im1,im2,4)=half*(paw_ij(iatom)%noccmmp(1,im1,im2,2)-paw_ij(iatom)%noccmmp(2,im1,im2,3))
             noccmmp2(2,im1,im2,4)=half*(paw_ij(iatom)%noccmmp(2,im1,im2,2)+paw_ij(iatom)%noccmmp(1,im1,im2,3))
           end do
         end do
         if(abs(pawprtvol)>=1) then
           write(message,'(2a)') ch10,"== Calculated occupation matrix for correlated orbitals in the n, m basis :"
           call wrtout(std_out,message,'COLL')
           do ispden=1,ndij
             write(message,'(3a)') ch10,"Calculated occupation matrix for component ",trim(dspinm(ispden+2*(nspden/4)))
             call wrtout(std_out,message,'COLL')
             do im1=1,lcur*2+1  ! ( order of indices in noccmmp is exchanged in order to have the same convention as rhoij: transposition is done after )
               if(cplex_dij==1)&
&               write(message,'(12(1x,9(1x,f10.5)))')&
&               (paw_ij(iatom)%noccmmp(1,im2,im1,ispden),im2=1,lcur*2+1)
               if(cplex_dij==2)&
&               write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))')&
&               (paw_ij(iatom)%noccmmp(:,im2,im1,ispden),im2=1,lcur*2+1)
               call wrtout(std_out,message,'COLL')
             end do
           end do
         end if ! pawprtvol >=1
       end if

!      Compute total number of electrons per spin
       paw_ij(iatom)%nocctot(:)=zero ! contains nmmp in the n m representation
       if(ndij==4) nocctot2(:)=zero ! contains nmmp in the upup dndn updn dnup  representation
       do ispden=1,ndij
         do im1=1,2*lcur+1
           if(ndij==4) then
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)+paw_ij(iatom)%noccmmp(1,im1,im1,ispden)
             nocctot2(ispden)=nocctot2(ispden)+noccmmp2(1,im1,im1,ispden)
           else
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)+paw_ij(iatom)%noccmmp(1,im1,im1,ispden)
           end if
         end do
       end do
!      noccmmp will now be in the up up , dn dn... representation and now n_mmp=<m|n|mp> instead of <mp|n|m> !
       if(ndij==4) then
         do ispden=1,ndij
           do iplex=1,cplex_dij
             do im1=1,2*lcur+1
               do im2=1,2*lcur+1
                 paw_ij(iatom)%noccmmp(iplex,im1,im2,ispden)=noccmmp2(iplex,im2,im1,ispden) ! now, noccmmp is in the upup dndn updn dnup representation
               end do
             end do
           end do
         end do
         ABI_DEALLOCATE(noccmmp2)
       end if
!      Printing of new nocc_mmp
       if (usepawu>0) write(message, '(2a)' )  ch10,'========== LDA+U DATA =================================================== '
       if (useexexch>0) write(message, '(2a)' )ch10,'======= Local ex-exchange (PBE0) DATA =================================== '
       call wrtout(std_out,message,'COLL')
       write(message,'(2a,i5,a,i4,a)') ch10,"====== For Atom", iatom,&
&       ", occupations for correlated orbitals. l =",lcur,ch10
       call wrtout(std_out,message,'COLL')
       if(ndij==2) then
         do ispden=1,2
           write(message,'(a,i4,3a,f10.5)') "Atom", iatom,". Occupations for spin ",&
&           trim(dspin(ispden))," =",paw_ij(iatom)%nocctot(ispden)
           call wrtout(std_out,message,'COLL')
         end do
         write(message,'(a,i4,a,2x,e16.8)') "=> On atom",iatom,", local Mag. is  ",&
&         paw_ij(iatom)%nocctot(2)-paw_ij(iatom)%nocctot(1)
         call wrtout(std_out,message,'COLL')
       end if
       if(ndij==4) then
         ntot=paw_ij(iatom)%nocctot(1)
         mx=paw_ij(iatom)%nocctot(2)
         my=paw_ij(iatom)%nocctot(3)
         mz=paw_ij(iatom)%nocctot(4)
         mnorm=sqrt(mx*mx+my*my+mz*mz)
         nup=nocctot2(1)
         ndn=nocctot2(2)
         write(message,'(a,i4,a,2x,e16.8)') "=> On atom",iatom,", local Mag. x is ",mx
         call wrtout(std_out,message,'COLL')
         write(message,'(14x,a,2x,e16.8)') "  local Mag. y is ",my
         call wrtout(std_out,message,'COLL')
         write(message,'(14x,a,2x,e16.8)') "  local Mag. z is ",mz
         call wrtout(std_out,message,'COLL')
         write(message,'(14x,a,2x,e16.8)') "  norm of Mag. is ",mnorm
         call wrtout(std_out,message,'COLL')
         write(message,'(14x,a,2x,f10.5)') "  occ. of majority spin is ",half*(ntot+mnorm)  ! to be checked versus direct calc from noccmmp
         call wrtout(std_out,message,'COLL')
         if(abs(pawprtvol)>=1) write(message,'(14x,a,2x,f10.5)') "  occ. for spin up (along z) ",nup
         if(abs(pawprtvol)>=1) then
           call wrtout(std_out,message,'COLL')
         end if
         write(message,'(14x,a,2x,f10.5)') "  occ. of minority spin is ",half*(ntot-mnorm)
         call wrtout(std_out,message,'COLL')
         if(abs(pawprtvol)>=1) write(message,'(14x,a,2x,f10.5)') "  occ. for spin dn (along z) ",ndn
         if(abs(pawprtvol)>=1) then
           call wrtout(std_out,message,'COLL')
         end if
         if(ndij==4)  then
           ABI_DEALLOCATE(nocctot2)
         end if
       end if
       write(message,'(2a)') ch10,"== Calculated occupation matrix for correlated orbitals:"
       call wrtout(std_out,message,'COLL')
       do ispden=1,ndij
         write(message,'(3a)') ch10,"Calculated occupation matrix for component ",trim(dspinc(ispden+2*(nspden/4)))
         call wrtout(std_out,message,'COLL')
         do im1=1,lcur*2+1
           if(cplex_dij==1)&
&           write(message,'(12(1x,9(1x,f10.5)))')&
&           (paw_ij(iatom)%noccmmp(1,im1,im2,ispden),im2=1,lcur*2+1)
           if(cplex_dij==2)&
&           write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))')&
&           (paw_ij(iatom)%noccmmp(:,im1,im2,ispden),im2=1,lcur*2+1)
           call wrtout(std_out,message,'COLL')
         end do
       end do

!      Transformation matrices: real->complex spherical harmonics (for test)
       if(ndij==4.and.abs(pawprtvol)>=0) then
         ABI_ALLOCATE(noccmmp_ylm,(2*lcur+1,2*lcur+1,ndij))
         noccmmp_ylm=czero
         ABI_ALLOCATE(noccmmp_slm,(2*lcur+1,2*lcur+1,ndij))
         noccmmp_slm=czero
         ABI_ALLOCATE(noccmmp_jmj,(2*(2*lcur+1),2*(2*lcur+1)))
         noccmmp_jmj=czero
!        go from real notation for complex noccmmp to complex notation in noccmmp_slm
         noccmmp_slm(:,:,:)=cmplx(paw_ij(iatom)%noccmmp(1,:,:,:)&
&         ,paw_ij(iatom)%noccmmp(2,:,:,:))
         call mat_slm2ylm(lcur,noccmmp_slm,noccmmp_ylm,ndij,1,1,pawprtvol) ! optspin=1: up spin are first

         do ispden=1,ndij
           write(message,'(3a)') ch10,"Calculated Ylm occupation matrix for component ",trim(dspinc(ispden+2*(ndij/4)))
           call wrtout(std_out,message,'COLL')
           do im1=1,lcur*2+1
             write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))') (noccmmp_ylm(im1,im2,ispden),im2=1,lcur*2+1)
             call wrtout(std_out,message,'COLL')
           end do
         end do
         call mat_mlms2jmj(lcur,noccmmp_ylm,noccmmp_jmj,ndij,1,1,pawprtvol) !  optspin=1: up spin are first
         ABI_DEALLOCATE(noccmmp_ylm)
         ABI_DEALLOCATE(noccmmp_jmj)
         ABI_DEALLOCATE(noccmmp_slm)
       end if !ndij==4

     end if ! impose_dmat==0

!    ########################################################################################
!    # Diagonalize nocc_mmp
!    ########################################################################################
     if(usepawu>0.and.dmatudiag_loc>0) then

       lpawu=lcur;ldim=2*lpawu+1
       ABI_ALLOCATE(noccmmp_tmp,(1,ldim,ldim,ndij))
       if (ndij==4)  then
         ABI_ALLOCATE(znoccmmp_tmp,(2*ldim,2*ldim))
       end if

!      Select noccmmp for this atom
       do ispden=1,ndij
         noccmmp_tmp(1,:,:,ispden)=paw_ij(iatom)%noccmmp(1,:,:,ispden)
       end do
       if (ndij==4) then
         do im2=1,ldim
           do im1=1,ldim
             znoccmmp_tmp(im1     ,     im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,1)&
&             ,paw_ij(iatom)%noccmmp(2,im1,im2,1))
             znoccmmp_tmp(ldim+im1,ldim+im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,2)&
&             ,paw_ij(iatom)%noccmmp(2,im1,im2,2))
             znoccmmp_tmp(     im1,ldim+im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,3)&
&             ,paw_ij(iatom)%noccmmp(2,im1,im2,3))
             znoccmmp_tmp(ldim+im1,     im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,4)&
&             ,paw_ij(iatom)%noccmmp(2,im1,im2,4))
           end do
         end do
       end if

!      Diagonalize nocc_mmp
       if (ndij/=4) then
         ABI_ALLOCATE(hdp,(ldim,ldim,ndij))
         hdp=zero
         lwork=3*ldim-1
         ABI_ALLOCATE(rwork,(lwork))
         ABI_ALLOCATE(eig,(ldim))
         do ispden=1,ndij
           call dsyev('v','u',ldim,noccmmp_tmp(1,:,:,ispden),ldim,eig,rwork,lwork,info)
           if(info/=0) then
             write(message,'(4a)') ch10,'Error in diagonalization of noccmmp (DSYEV)!'
             call wrtout(std_out,message,'COLL')
             call leave_new('COLL')
           end if
           do ilm=1,ldim
             hdp(ilm,ilm,ispden)=eig(ilm)
           end do
         end do ! ispden
         ABI_DEALLOCATE(rwork)
         ABI_DEALLOCATE(eig)
       else
         ABI_ALLOCATE(hdp,(2*ldim,2*ldim,1))
         hdp=zero
         lwork=4*ldim-1
         ABI_ALLOCATE(rwork,(6*ldim-2))
         ABI_ALLOCATE(zwork,(lwork))
         ABI_ALLOCATE(eig,(2*ldim))
         call zheev('v','u',2*ldim,znoccmmp_tmp,2*ldim,eig,zwork,lwork,rwork,info)
         if(info/=0) then
           write(message,'(4a)') ch10,'Error in diagonalization of znoccmmp_tmp (zheev) !'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
         do ilm=1,2*ldim
           hdp(ilm,ilm,1)=eig(ilm)
         end do
         ABI_DEALLOCATE(rwork)
         ABI_DEALLOCATE(zwork)
         ABI_DEALLOCATE(eig)
       end if

!      Print diagonalized matrix and eigenvectors
       do ispden=1,size(hdp,3)
         write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,' == Diagonalized Occupation matrix'
         if (ndij==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
         if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
         if (ndij==4) write(message,fmt='(2a,i3,a)')trim(message)," =="
         call wrtout(std_out,message,'COLL')
         do ilm=1,size(hdp,1)
           write(message,'(12(1x,9(1x,f10.5)))') (hdp(ilm,jlm,ispden),jlm=1,size(hdp,2))
           call wrtout(std_out,message,'COLL')
         end do
       end do ! ispden
       if(abs(pawprtvol)>=1) then
         if (ndij/=4) then
           do ispden=1,ndij
             write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,' == Eigenvectors'
             if (ndij==1) write(message,fmt='(2a)')     trim(message),' for spin up =='
             if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message),' for spin ',ispden,' =='
             call wrtout(std_out,message,'COLL')
             do ilm=1,ldim
               write(message,'(12(1x,9(1x,f10.5)))') (noccmmp_tmp(1,ilm,jlm,ispden),jlm=1,ldim)
               call wrtout(std_out,message,'COLL')
             end do
           end do
         else
           write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,' == Eigenvectors (spinors) in the real harmonics basis =='
           call wrtout(std_out,message,'COLL')
           do ilm=1,2*ldim
             write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (znoccmmp_tmp(ilm,jlm),jlm=1,2*ldim)
             call wrtout(std_out,message,'COLL')
           end do
         end if
       end if

!      Back rotation of diagonalized matrix and printing
       if(abs(pawprtvol)>=1) then
         if (ndij/=4) then
           ABI_ALLOCATE(hdp2,(ldim,ldim))
           do ispden=1,ndij
             call dgemm('n','t',ldim,ldim,ldim,one,hdp(:,:,ispden),ldim,noccmmp_tmp(1,:,:,ispden),ldim,zero,hdp2,ldim)
             call dgemm('n','n',ldim,ldim,ldim,one,noccmmp_tmp(1,:,:,ispden),ldim,hdp2,ldim,zero,hdp(:,:,ispden),ldim)
             noccmmp_tmp(1,:,:,ispden)=hdp(:,:,ispden)
           end do ! ispden
           ABI_DEALLOCATE(hdp2)
         else
           ABI_ALLOCATE(zhdp,(2*ldim,2*ldim))
           ABI_ALLOCATE(zhdp2,(2*ldim,2*ldim))
           zhdp(:,:)=cmplx(hdp(:,:,1),zero,kind=dp)
           zhdp2(:,:)=cmplx(zero,zero,kind=dp)
           call zgemm('n','c',2*ldim,2*ldim,2*ldim,cone,zhdp,2*ldim,znoccmmp_tmp,2*ldim,czero,zhdp2,2*ldim)
           zhdp(:,:)=cmplx(zero,zero,kind=dp)
           call zgemm('n','n',2*ldim,2*ldim,2*ldim,cone,znoccmmp_tmp,2*ldim,zhdp2,2*ldim,czero,zhdp,2*ldim)
           znoccmmp_tmp=zhdp
           ABI_DEALLOCATE(zhdp)
           ABI_DEALLOCATE(zhdp2)
         end if
         nmat=ndij ; if(ndij==4.and.cplex_dij==2) nmat=1
         do ispden=1,nmat
           write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&           ' == Rotated back diagonalized matrix'
           if (ndij==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
           if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
!          if (nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
!          &      trim(dspin(ispden+2*(nspden/4)))," =="
           if (ndij==4.and.cplex_dij==2) write(message,fmt='(4a)')     trim(message)," for all component "
           call wrtout(std_out,message,'COLL')
           do ilm=1,ldim*cplex_dij
             if(ndij==1.or.ndij==2)&
&             write(message,'(12(1x,9(1x,f10.5)))')&
&             (noccmmp_tmp(1,ilm,jlm,ispden),jlm=1,ldim)
             if(ndij==4.and.cplex_dij==2)&
&             write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))')&
&             (znoccmmp_tmp(ilm,jlm),jlm=1,ldim*cplex_dij)
             call wrtout(std_out,message,'COLL')
           end do
         end do ! ispden
       end if
       ABI_DEALLOCATE(hdp)

     end if ! dmatudiag_loc

!    ########################################################################################
!    # Impose value of nocc_mmp from dmatpu; symetrize it
!    ########################################################################################
     if (usepawu>0.and.impose_dmat/=0) then

       lpawu=lcur
       nsploop=nsppol;if (ndij==4) nsploop=4
       noccsym_error=.false.

!      Loop over spin components
       do ispden=1,nsploop
         if (ndij/=4) then
           jspden=min(3-ispden,paw_ij(iatom)%nsppol)
         else if (ispden<=2) then
           jspden=3-ispden
         else
           jspden=ispden
         end if

!        Loops over components of nocc_mmp
         do jlm=1,2*lpawu+1
           do ilm=1,2*lpawu+1

             if(nsym>1.and.ndij<4) then

               nsym_used(1:2)=0
               sumocc(1:2)=zero

!              Accumulate values of nocc_mmp over symmetries
               do irot=1,nsym
                 if ((symafm(irot)/=1).and.(.not.use_afm)) cycle
                 kspden=ispden;if (symafm(irot)==-1) kspden=jspden
                 factafm=one;if (ispden>3) factafm=dble(symafm(irot))
                 iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2
                 nsym_used(iafm)=nsym_used(iafm)+1
                 at_indx=indsym(4,irot,iatom)
                 do im2=1,2*lpawu+1
                   do im1=1,2*lpawu+1
!                    Be careful: use here R_rel^-1 in term of spherical harmonics
!                    which is tR_rec in term of spherical harmonics
!                    so, use transpose[zarot]
                     sumocc(iafm)=sumocc(iafm)+factafm*tmp_noccmmp(at_indx)%value(1,im1,im2,kspden) &
&                     *pawang%zarot(im1,ilm,lpawu+1,irot)&
&                     *pawang%zarot(im2,jlm,lpawu+1,irot)
!                    sumocc(iafm)=sumocc(iafm)+factafm*tmp_noccmmp(at_indx)%value(im1,im2,kspden) &
!                    &                     *pawang%zarot(ilm,im1,lpawu+1,irot)&
!                    &                     *pawang%zarot(jlm,im2,lpawu+1,irot)
                   end do
                 end do
               end do ! End loop over symmetries

!              Store new values of nocc_mmp
               paw_ij(iatom)%noccmmp(1,ilm,jlm,ispden)=sumocc(1)/nsym_used(1)
               if (.not.noccsym_error)&
&               noccsym_error=(abs(paw_ij(iatom)%noccmmp(1,ilm,jlm,ispden)-tmp_noccmmp(iatom)%value(1,ilm,jlm,ispden))>tol6)

!              Antiferromagnetic case: has to fill up "down" component of nocc_mmp
               if (antiferro.and.nsym_used(2)>0) paw_ij(iatom)%noccmmp(1,ilm,jlm,2)=sumocc(2)/nsym_used(2)

             else  ! nsym=1

!              Case without symetries
               paw_ij(iatom)%noccmmp(:,ilm,jlm,ispden)= tmp_noccmmp(iatom)%value(:,ilm,jlm,ispden)
             end if

           end do !ilm
         end do !jlm
       end do ! ispden

!      Printing of new nocc_mmp
       do ispden=1,ndij
         if(dmatudiag_loc==2) then
           write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&           ' == Imposed occupation matrix (in the basis of diagonalization!!)'
         else
           write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&           ' == Imposed occupation matrix'
         end if
         if (ndij==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
         if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
         if (ndij==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&         trim(dspinc(ispden+2*(ndij/4)))," =="
         call wrtout(std_out,message,'COLL')
         do ilm=1,2*lpawu+1
           if(cplex_dij==1)&
&           write(message,'(12(1x,9(1x,f10.5)))')&
&           (paw_ij(iatom)%noccmmp(1,ilm,jlm,ispden),jlm=1,2*lpawu+1)
           if(cplex_dij==2)&
&           write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))')&
&           (paw_ij(iatom)%noccmmp(:,ilm,jlm,ispden),jlm=1,2*lpawu+1)
           call wrtout(std_out,message,'COLL')
         end do
       end do

!      WARNING if symmetrization changes the matrix
       if (noccsym_error) then
         write(message, '(4a,i4,6a)' ) ch10,&
&         ' setnoccmmp: WARNING - ',ch10,&
&         '   After symmetrization, imposed occupation matrix for atom ',iatom,ch10,&
&         '   is different from dmatpawu value set in input file !',ch10,&
&         '   It is likely that dmatpawu does not match the symmetry operations of the system.',ch10,&
&         '   Action: change dmatpawu in input file.'
         call wrtout(std_out,message,'COLL')
       end if

     end if ! impose_dmat/=0

!    ########################################################################################
!    # Rotate imposed occupation matrix in the non-diagonal basis
!    ########################################################################################
     if (usepawu>0.and.impose_dmat/=0.and.dmatudiag_loc==2) then

       lpawu=lcur;ldim=2*lpawu+1

!      Rotation of imposed nocc_mmp
       if (ndij/=4) then
         ABI_ALLOCATE(hdp2,(ldim,ldim))
         do ispden=1,ndij
           call dgemm('n','t',ldim,ldim,ldim,one,&
&           paw_ij(iatom)%noccmmp(1,:,:,ispden),ldim,noccmmp_tmp(1,:,:,ispden),ldim,zero,hdp2,ldim)
           call dgemm('n','n',ldim,ldim,ldim,one,&
&           noccmmp_tmp(1,:,:,ispden),ldim,hdp2,ldim,zero,paw_ij(iatom)%noccmmp(1,:,:,ispden),ldim)
         end do ! ispden
         ABI_DEALLOCATE(hdp2)
       else
         ABI_ALLOCATE(zhdp,(2*ldim,2*ldim))
         ABI_ALLOCATE(zhdp2,(2*ldim,2*ldim))
         do im2=1,ldim
           do im1=1,ldim
             zhdp(     im1,     im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,1),zero,kind=dp)  ! to be checked
             zhdp(ldim+im1,ldim+im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,2),zero,kind=dp)  ! to be checked
             zhdp(     im1,ldim+im2)=&
&             cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,3),+paw_ij(iatom)%noccmmp(1,im1,im2,4),kind=dp)  ! to be checked
             zhdp(ldim+im1,     im2)=&
&             cmplx(paw_ij(iatom)%noccmmp(1,im2,im1,3),-paw_ij(iatom)%noccmmp(1,im2,im1,4),kind=dp)  ! to be checked
           end do
         end do
         call zgemm('n','c',2*ldim,2*ldim,2*ldim,cone,zhdp,2*ldim,znoccmmp_tmp,2*ldim,czero,zhdp2,2*ldim)
         call zgemm('n','n',2*ldim,2*ldim,2*ldim,cone,znoccmmp_tmp,2*ldim,zhdp2,2*ldim,czero,zhdp,2*ldim)
         do jlm=1,ldim
           do ilm=1,ldim
             paw_ij(iatom)%noccmmp(1,ilm,jlm,1)= real(znoccmmp_tmp(     ilm,     jlm))  ! to be checked
             paw_ij(iatom)%noccmmp(1,ilm,jlm,2)= real(znoccmmp_tmp(ldim+ilm,ldim+jlm))  ! to be checked
             paw_ij(iatom)%noccmmp(1,ilm,jlm,3)= real(znoccmmp_tmp(     ilm,ldim+jlm))  ! to be checked
             paw_ij(iatom)%noccmmp(1,ilm,jlm,4)=aimag(znoccmmp_tmp(     ilm,ldim+jlm))  ! to be checked
           end do
         end do
         ABI_DEALLOCATE(zhdp)
         ABI_DEALLOCATE(zhdp2)
       end if

!      Printing of rotated imposed matrix
       do ispden=1,ndij
         write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&         ' == Imposed density matrix in original basis'
         if (ndij==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
         if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
         if (ndij==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&         trim(dspin(ispden+2*(ndij/4)))," =="
         call wrtout(std_out,message,'COLL')
         do ilm=1,2*lpawu+1
           write(message,'(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(1,ilm,jlm,ispden),jlm=1,2*lpawu+1)  ! to be checked
           call wrtout(std_out,message,'COLL')
         end do
       end do ! ispden

     end if ! dmatudiag_loc==2

     if (usepawu>0.and.dmatudiag_loc>0) then
       ABI_DEALLOCATE(noccmmp_tmp)
       if (ndij==4)  then
         ABI_DEALLOCATE(znoccmmp_tmp)
       end if
     end if

     paw_ij(iatom)%has_pawu_occ=2

   end if ! lcur
 end do ! iatom

!Memory deallocation
 if (usepawu>0.and.impose_dmat/=0) then
   do iatom=1,natom
     lpawu=pawtab(pawrhoij(iatom)%itypat)%lpawu
     if (lpawu/=-1)  then
       ABI_DEALLOCATE(tmp_noccmmp(iatom)%value)
     end if
   end do
   ABI_DEALLOCATE(tmp_noccmmp)
 end if

 DBG_EXIT("COLL")

end subroutine setnoccmmp
!!***
