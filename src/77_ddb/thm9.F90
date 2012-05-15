!{\src2tex{textfont=tt}}
!!****f* ABINIT/thm9
!!
!! NAME
!! thm9
!!
!! FUNCTION
!! This routine to calculate phonon density of states,
!! thermodynamical properties, Debye-Waller factor, and atomic
!! mean square velocity
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (CL,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3) =length scales by which rprim is to be multiplied
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! atmfrc(2,3,natom,3,natom,nrpt) = Interatomic Forces in real space
!! dielt(3,3)    = dielectric tensor
!! dyewq0(3,3,natom)=Ewald part of the dynamical matrix, at q=0.
!! gmet(3,3)= metric tensor in reciprocal space.
!! gprim(3,3)=normalized coordinates in reciprocal space
!! indsym(4,nsym,natom)=label given by subroutine symatm, indicating atom
!!  label which gets rotated into given atom by given symmetry
!!  (first three elements are related primitive translation--
!!  see symatm where this is computed)
!! iout =unit number for output
!! mpert =maximum number of ipert
!! msym = maximum number of symmetries
!! natom=number of atoms in the unit cell
!! nrpt=number of R points in the Big Box
!! nsym=number of symmetries
!! ntypat=number of atom types
!! outfilename_radix=radix of anaddb output file name: append _THERMO for
!!     thermodynamic quantities
!! rmet(3,3)=metric tensor in real space.
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nprt)=canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!! symrec(3,3,nsym)=symmetry operations
!! symrel(3,3,nsym)=symmetry operations
!! tcpui=initial cpu time
!! trans(3,natom)=atomic translations : xred = rcan + trans
!! twalli=initial wall clock time
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume
!! wghatm(natom,natom,nrpt)
!!  =weights associated to a pair of atoms and to a R vector
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!! eigval(3*natom)=contains the eigenvalues of the dynamical matrix
!! eigvec(2*3*natom*3*natom)= at the end, contain
!!  the eigenvectors of the dynamical matrix.
!! displ(2*3*natom*3*natom)= at the end, contain
!!  the displacements of atoms in cartesian coordinates.
!! d2cart(2,3,mpert,3,mpert)=dynamical matrix
!! phfrq(3*natom)=phonon frequencies (square root of the dynamical
!!  matrix eigenvalues, except if these are negative, and in this
!!  case, give minus the square root of the absolute value
!!  of the matrix eigenvalues). Hartree units.
!!
!! NOTES
!! Work space:
!!   dosinc,eigval,eigvec,displ,d2cart,phfrq
!! dosinc=increment between the channels for the phonon DOS in cm-1
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      end_sortph,gtdyn9,leave_new,matr3inv,mkrdim,phfrq3,smpbz,sortph,symkpt
!!      timein,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine thm9(acell,amu,anaddb_dtset,atmfrc,dielt,displ,&
& dyewq0,d2cart,eigval,eigvec,gmet,gprim,indsym,iout,mpert,msym,natom,&
& nrpt,nsym,ntypat,outfilename_radix,phfrq,rmet,rprim,rpt,symrec,symrel,tcpui,&
& trans,twalli,typat,ucvol,wghatm,xred,zeff,themflag, udispl, ufreq)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_sortph

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'thm9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_56_recipspace
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => thm9
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,msym,natom,nrpt,nsym,ntypat
 integer,intent(in),optional :: themflag,udispl,ufreq
 real(dp),intent(in) :: tcpui,twalli,ucvol
 character(len=fnlen) :: outfilename_radix
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: acell(3),amu(ntypat),dielt(3,3),gmet(3,3),gprim(3,3)
 real(dp),intent(in) :: rmet(3,3),rprim(3,3),rpt(3,nrpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt),xred(3,natom)
 real(dp),intent(in) :: zeff(3,3,natom)
 real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt),dyewq0(3,3,natom)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert),displ(2*3*natom*3*natom)
 real(dp),intent(out) :: eigval(3*natom),eigvec(2,3,natom,3*natom)
 real(dp),intent(out) :: phfrq(3*natom)

!Local variables -------------------------
!scalars
 integer :: convth,facbrv,iatom,iavfrq,ichan,icomp,igrid,iii,iiii,ij,ind
 integer :: iqpt2,isym,itemper,iwchan,jjj,mqpt2,nchan,ngrids
 integer :: nqpt2,nspqpt,ntemper,nwchan,option,timrev
 integer :: thermal_unit
 character(len=fnlen) :: thermal_filename
 real(dp) :: change,cothx,diffbb,dosinc,expm2x,factor,factorw,factorv,gerr
 real(dp) :: ggsum,ggsumsum,ggrestsum
 real(dp) :: gijerr,gijsum,gnorm,ln2shx,qphnrm,relchg,tcpu,tmp,twall,wovert
 logical :: part1,part2
 character(len=500) :: message
!arrays
 integer :: igqpt2(3),ii(6),jj(6),qptrlatt(3,3)
 integer,allocatable :: indqpt1(:),nchan2(:)
 real(dp) :: gprimd(3,3),qphon(3),rprimd(3,3),tens1(3,3),tens2(3,3)
 real(dp),allocatable :: bbij(:,:,:),bij(:,:,:),energy(:),energy0(:),entropy(:)
 real(dp),allocatable :: entropy0(:),free(:),free0(:),gdos(:,:),gg(:,:),gg_sum(:,:),gg_rest(:,:)
 real(dp),allocatable :: ggij(:,:,:,:),gij(:,:,:,:),gw(:,:),qpt2(:,:),spheat(:)
 real(dp),allocatable :: spheat0(:),spqpt2(:,:),wme(:),wtq(:),wtq2(:)
 real(dp),allocatable :: wtq_folded(:),vij(:,:,:)
 logical,allocatable :: wgcnv(:),wgijcnv(:)

! *********************************************************************

!DEBUG
!write(std_out,*)' thm9 : enter '
!ENDDEBUG

!If timing is needed at some point ...
 call timein(tcpu,twall)
 write(message, '(a,f11.4,a,f11.4,a)' )&
& ' thm9 : begin at tcpu',tcpu-tcpui,&
& '  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

 thermal_unit = 400 ! file unit for writing out thermal properties
 thermal_filename=trim(outfilename_radix)//"_THERMO"
 open (file=thermal_filename, unit=thermal_unit)
 write (thermal_unit,*) '#'
 write (thermal_unit,*) '#  Thermodynamical quantities calculated by ANADDB'
 write (thermal_unit,*) '#'

 nchan=anaddb_dtset%nchan
 ntemper=anaddb_dtset%ntemper
 nwchan=anaddb_dtset%nwchan
 ngrids=anaddb_dtset%ngrids
 iavfrq=anaddb_dtset%iavfrq

 ABI_ALLOCATE(bbij,(6,natom,ntemper))
 ABI_ALLOCATE(bij,(6,natom,ntemper))
 ABI_ALLOCATE(vij,(6,natom,ntemper))
 ABI_ALLOCATE(energy,(ntemper))
 ABI_ALLOCATE(energy0,(ntemper))
 ABI_ALLOCATE(entropy,(ntemper))
 ABI_ALLOCATE(entropy0,(ntemper))
 ABI_ALLOCATE(free,(ntemper))
 ABI_ALLOCATE(free0,(ntemper))
!Doubling the size (nchan) of gg_sum is needed, because the maximum sum frequency is the double
!the for maximum frequency. However, for the write statement, it is better to double
!also the size of the other arrays ...
 ABI_ALLOCATE(gdos,(2*nchan,nwchan))
 ABI_ALLOCATE(gg,(2*nchan,nwchan))
 ABI_ALLOCATE(gg_sum,(2*nchan,nwchan))
 ABI_ALLOCATE(gg_rest,(2*nchan,nwchan))
 ABI_ALLOCATE(ggij,(6,natom,nchan,nwchan))
 ABI_ALLOCATE(gij,(6,natom,nchan,nwchan))
 ABI_ALLOCATE(nchan2,(nwchan))
 ABI_ALLOCATE(spheat,(ntemper))
 ABI_ALLOCATE(spheat0,(ntemper))
 ABI_ALLOCATE(wgcnv,(nwchan))
 ABI_ALLOCATE(wgijcnv,(nwchan))
 if (iavfrq>0)  then
   ABI_ALLOCATE(wme,(ntemper))
   ABI_ALLOCATE(gw,(nchan,nwchan))
 end if

 call mkrdim(acell,rprim,rprimd)

!initialize ii and jj arrays
 ii(1)=1 ; ii(2)=2 ; ii(3)=3
 ii(4)=1 ; ii(5)=1 ; ii(6)=2
 jj(1)=1 ; jj(2)=2 ; jj(3)=3
 jj(4)=2 ; jj(5)=3 ; jj(6)=3

!Put zeros in the bins of channels of frequencies
 gdos(:,:)=0._dp
 gg_sum(:,:)=0._dp
 gg_rest(:,:)=0._dp
 gij(:,:,:,:)=0._dp

!Neither part1 nor part2 are assumed converged initially
!None of the channel widths are assumed converged initially

 part1=.false.
 part2=.false.

 wgcnv(:)=.false.
 wgijcnv(:)=.false.

!Thermodynamic quantities are put to zero.
!(If exactly zero, initial convergence tests will fail.)

 free0(:)=1.d-05
 energy0(:)=1.d-05
 entropy0(:)=1.d-05
 spheat0(:)=1.d-05
 bij(:,:,:)=1.d-05
 vij(:,:,:)=1.d-05

!Number of actual channels set
 do iwchan=1,nwchan
   nchan2(iwchan)=(nchan-1)/iwchan+1
 end do

!For different types of Bravais lattices
 facbrv=1
 if(anaddb_dtset%brav==2)facbrv=2
 if(anaddb_dtset%brav==3)facbrv=4

 write(message, '(a)' )' thm9 : initialized'
 call wrtout(std_out,message,'COLL')

!Loops on the q point grids
 do igrid=1,ngrids

   igqpt2(:)=(anaddb_dtset%ng2qpt(:)*igrid)/ngrids
   mqpt2=(igqpt2(1)*igqpt2(2)*igqpt2(3))/facbrv
   ABI_ALLOCATE(qpt2,(3,mqpt2))
   ABI_ALLOCATE(spqpt2,(3,mqpt2))

!  DEBUG
!  write(std_out,*)' igrid',igrid
!  write(std_out,*)' iout',iout
!  write(std_out,*)' mqpt2',mqpt2
!  write(std_out,*)' igqpt2',igqpt2
!  write(std_out,*)' nspqpt',nspqpt
!  write(std_out,*)' nsym',nsym
!  ENDDEBUG

   option=1
   qptrlatt(:,:)=0
   qptrlatt(1,1)=igqpt2(1)
   qptrlatt(2,2)=igqpt2(2)
   qptrlatt(3,3)=igqpt2(3)
   call smpbz(anaddb_dtset%brav,iout,qptrlatt,mqpt2,nspqpt,1,&
&   option,anaddb_dtset%q2shft,spqpt2)
   write(message,'(a)')' thm9 : exit smpbz '
   call wrtout(std_out,message,'COLL')

   ABI_ALLOCATE(indqpt1,(nspqpt))
   ABI_ALLOCATE(wtq,(nspqpt))
   ABI_ALLOCATE(wtq_folded,(nspqpt))

!  Reduce the number of such points by symmetrization
   wtq(:)=1.0_dp

   timrev=1 
   call symkpt(0,gmet,indqpt1,ab_out,spqpt2,nspqpt,nqpt2,nsym,&
&   symrec,timrev,wtq,wtq_folded)
   ABI_ALLOCATE(wtq2,(nqpt2))
   do iqpt2=1,nqpt2
     wtq2(iqpt2)=wtq_folded(indqpt1(iqpt2))
     qpt2(:,iqpt2)=spqpt2(:,indqpt1(iqpt2))
!    DEBUG
!    write(std_out,*)' thm9 : iqpt2, wtq2 :',iqpt2,wtq2(iqpt2)
!    ENDDEBUG
   end do
   ABI_DEALLOCATE(wtq_folded)
   call timein(tcpu,twall)
   write(message, '(a,a,a,f11.4,a,f11.4,a)' )&
&   ' thm9 : in the loop over grids, after symkpt',ch10,&
&   '  tcpu',tcpu-tcpui,' twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')

!  Temporary counters are put zero.

   gg(:,:)=zero
   ggij(:,:,:,:)=zero
   if (iavfrq>0) gw(:,:)=zero

!  Sum over the sampled q points

   do iqpt2=1,nqpt2

     qphon(:)=qpt2(:,iqpt2)
     qphnrm=1.0_dp

!    Get d2cart using the interatomic forces and the
!    long-range coulomb interaction through Ewald summation
     call gtdyn9(acell,atmfrc,dielt,anaddb_dtset%dipdip,&
&     dyewq0,d2cart,gmet,gprim,mpert,natom,&
&     nrpt,qphnrm,qphon,rmet,rprim,rpt,&
&     trans,ucvol,wghatm,xred,zeff)

!    Calculation of the eigenvectors and eigenvalues
!    of the dynamical matrix
!    
     call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&     mpert,msym,natom,nsym,ntypat,phfrq,qphnrm,qphon,&
&     rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol)

     if (present(themflag)) then
       if (themflag ==2)then
         call sortph(eigvec,displ,outfilename_radix, natom,phfrq,udispl,ufreq)
       end if
     end if

!    Sum over the phonon modes

     do iii=1,3*natom

!      Slightly negative frequencies are put to zero
!      Imaginary frequencies are also put to zero
       if(phfrq(iii)<0._dp) phfrq(iii)=0._dp

!      Note: frequencies are now in cm-1

!      Sort frequencies into channels of frequencies
!      for each channel width of frequency

       do iwchan=nwchan,1,-1
         if (.not.wgcnv(iwchan))then

           dosinc=dble(iwchan)
           ichan=int(phfrq(iii)*Ha_cmm1/dosinc)+1

           if(ichan>nchan2(iwchan)) then
             write(message, '(a,a,a,es16.6,a,a,a,i7,a,a,a)' )&
&             ' thm9 : ERROR -',ch10,&
&             '  There is a phonon frequency,',phfrq(iii),' larger than the maximum one,',ch10,&
&             '  as defined by the number of channels of width 1 cm-1, nchan=',nchan,'.',ch10,&
&             '  Action : increase nchan (suggestion : double it).'
             call wrtout(std_out,message,'COLL')
             call leave_new('COLL')
           end if

           gg(ichan,iwchan)=gg(ichan,iwchan)+wtq2(iqpt2)

           if (iavfrq>0) gw(ichan,iwchan)=gw(ichan,iwchan)+wtq2(iqpt2)*phfrq(iii)*Ha_cmm1

!          to calculate two phonon DOS for qshift = 0.0
           do iiii=1,3*natom
             
             if(phfrq(iiii)<0.d0) phfrq(iiii)=0.d0

             ichan=int(abs(phfrq(iii)+phfrq(iiii))*Ha_cmm1/dosinc)+1
             gg_sum(ichan,iwchan)=gg_sum(ichan,iwchan)+wtq2(iqpt2)

             ichan=int(abs(phfrq(iii)-phfrq(iiii))*Ha_cmm1/dosinc)+1
             gg_rest(ichan,iwchan)=gg_rest(ichan,iwchan)+wtq2(iqpt2)

           end do ! end loop for iiii

         end if
       end do

       do iwchan=nwchan,1,-1
         if (.not.wgijcnv(iwchan))then

           dosinc=dble(iwchan)
           ichan=int(phfrq(iii)*Ha_cmm1/dosinc)+1

           if(ichan>nchan2(iwchan)) then
             write(message, '(a,a,a,a,a,a,a)' )&
&             ' thm9 : ERROR -',ch10,&
&             '  Phonon frequencies larger than the maximum one,',ch10,&
&             '  as defined by the number of channels of width 1 cm-1.',ch10,&
&             '  Action : increase nchan (suggestion : double it).'
             call wrtout(std_out,message,'COLL')
             call leave_new('COLL')
           end if

           do iatom=1,natom
             do ij=1,6
               ggij(ij,iatom,ichan,iwchan)=ggij(ij,iatom,ichan,iwchan)&
&               +wtq2(iqpt2)*&
&               (eigvec(1,ii(ij),iatom,iii)*eigvec(1,jj(ij),iatom,iii)+&
&               eigvec(2,ii(ij),iatom,iii)*eigvec(2,jj(ij),iatom,iii) )
             end do
           end do

         end if
       end do

!      End of the sum over the phonon modes
     end do

!    End of the sum over q-points in the irreducible zone
   end do

!  deallocate sortph array
   call end_sortph()

!  Symmetrize the gij
   call matr3inv(rprimd,gprimd)
   do ichan=1,nchan
     do iwchan=nwchan,1,-1
       do iatom=1,natom
         do ij=1,6
!          Uses bbij as work space
           bbij(ij,iatom,1)=ggij(ij,iatom,ichan,iwchan)
           ggij(ij,iatom,ichan,iwchan)=0.0_dp
         end do
       end do
       do iatom=1,natom
         do isym=1,nsym
!          Find the atom that will be applied on atom iatom
           ind=indsym(4,isym,iatom)
           do ij=1,6
             tens1(ii(ij),jj(ij))=bbij(ij,ind,1)
           end do
           tens1(2,1)=tens1(1,2)
           tens1(3,1)=tens1(1,3)
           tens1(3,2)=tens1(2,3)
!          Here acomplishes the tensorial operations
!          1) Goto reduced coordinates for both indices
           do iii=1,3
             do jjj=1,3
               tens2(iii,jjj)=tens1(iii,1)*gprimd(1,jjj)&
&               +tens1(iii,2)*gprimd(2,jjj)&
&               +tens1(iii,3)*gprimd(3,jjj)
             end do
           end do
           do jjj=1,3
             do iii=1,3
               tens1(iii,jjj)=tens2(1,jjj)*gprimd(1,iii)&
&               +tens2(2,jjj)*gprimd(2,iii)&
&               +tens2(3,jjj)*gprimd(3,iii)
             end do
           end do
!          2) Apply the symmetry operation on both indices
           do iii=1,3
             do jjj=1,3
               tens2(iii,jjj)=tens1(iii,1)*symrec(1,jjj,isym)&
&               +tens1(iii,2)*symrec(2,jjj,isym)&
&               +tens1(iii,3)*symrec(3,jjj,isym)
             end do
           end do
           do jjj=1,3
             do iii=1,3
               tens1(iii,jjj)=tens2(1,jjj)*symrec(1,iii,isym)&
&               +tens2(2,jjj)*symrec(2,iii,isym)&
&               +tens2(3,jjj)*symrec(3,iii,isym)
             end do
           end do
!          3) Go back to cartesian coordinates
           do iii=1,3
             do jjj=1,3
               tens2(iii,jjj)=tens1(iii,1)*rprimd(jjj,1)&
&               +tens1(iii,2)*rprimd(jjj,2)&
&               +tens1(iii,3)*rprimd(jjj,3)
             end do
           end do
           do jjj=1,3
             do iii=1,3
               tens1(iii,jjj)=tens2(1,jjj)*rprimd(iii,1)&
&               +tens2(2,jjj)*rprimd(iii,2)&
&               +tens2(3,jjj)*rprimd(iii,3)
             end do
           end do

           do ij=1,6
             ggij(ij,iatom,ichan,iwchan)=ggij(ij,iatom,ichan,iwchan)&
&             + tens1(ii(ij),jj(ij))
           end do

         end do
         do ij=1,6
           ggij(ij,iatom,ichan,iwchan)=ggij(ij,iatom,ichan,iwchan)/dble(nsym)
         end do
       end do
     end do
   end do

   write(message,'(a)' )&
&   ' thm9: g(w) and gij(k|w) calculated given a q sampling grid.'
   call wrtout(std_out,message,'COLL')

!  Sum up the counts in the channels to check the normalization
!  and test convergence with respect to q sampling

   gnorm=(igqpt2(1)*igqpt2(2)*igqpt2(3)*3*natom)/facbrv

   if(.not.part1)then
     do iwchan=nwchan,1,-1

       write(message,'(a,i6)' )' thm9 : iwchan=',iwchan
       call wrtout(std_out,message,'COLL')

       if (wgcnv(iwchan)) cycle
       write(message,'(a)' )' thm9 : compute g, f, e, s, c'
       call wrtout(std_out,message,'COLL')

!      Calculate g(w) and F,E,S,C

       ggsum=0.0_dp
       ggsumsum=0.0_dp
       ggrestsum=0.0_dp
       do ichan=1,nchan2(iwchan)
         ggsum=ggsum+gg(ichan,iwchan)
         ggsumsum=ggsumsum+gg_sum(ichan,iwchan)
         ggrestsum=ggrestsum+gg_rest(ichan,iwchan)
       end do

       if(ggsum/=gnorm)then
         write(message, '(a,a,a,es14.6,i6,a,a,es14.6,a)' )&
&         ' thm9 : BUG -',ch10,&
&         '  Frequencies are missing in g(w) : ggsum,iwchan=',ggsum,iwchan,&
&         ch10,'  gnorm=',gnorm,'.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if

!      Check if the density of states changed by more than dostol

       gerr=zero
       if (ngrids>1) then
         do ichan=1,nchan2(iwchan)
           gerr=gerr&
&           +abs(gg(ichan,iwchan)/ggsum-gdos(ichan,iwchan))
         end do
       end if

       if(gerr>anaddb_dtset%dostol.and.ngrids>1) then
         wgcnv(iwchan)=.false.
       else
         wgcnv(iwchan)=.true.
       end if

!      g(w) is updated
       do ichan=1,nchan2(iwchan)
         gdos(ichan,iwchan)=gg(ichan,iwchan)/ggsum
         gg_sum(ichan,iwchan)=gg_sum(ichan,iwchan)/ggsumsum
         gg_rest(ichan,iwchan)=gg_rest(ichan,iwchan)/ggrestsum
       end do
       if (iavfrq>0) then
         do ichan=1,nchan2(iwchan)
           gw(ichan,iwchan)=gw(ichan,iwchan)/ggsum
         end do
       end if

!      Write gerr for each q sampling and w width

       write(message,'(a,a,i3,3i6,f10.1,f10.5)') ch10, &
&       'thm9: iwchan,igqpt(:),norm,error=',&
&       iwchan,igqpt2(1),igqpt2(2),igqpt2(3),ggsum+tol8,gerr+tol10
       call wrtout(std_out,message,'COLL')
!      call wrtout(std_out,message,'COLL')

!      If the DOS with a channel width is newly converged,
!      print it out and calculate the thermodynamic functions.
       convth=0
       if(wgcnv(iwchan)) then

!        if (ngrids==1) then
!        write(message,'(a,i5,a)') ' gij with channel width=  ',iwchan,':'
!        else
!        write(message,'(a,i5,a)') ' gij with channel width=  ',&
!        &           iwchan,' newly converged'
!        end if
!        call wrtout(std_out,message,'COLL')
!        do iatom=0,natom
!        do ichan=1,nchan2(iwchan)
!        write(std_out,'(i6,f12.4,1x,2i3)') ichan,&
!        &             ggij(1,iatom,ichan,iwchan)+&
!        &             ggij(2,iatom,ichan,iwchan)+&
!        &             ggij(3,iatom,ichan,iwchan)+&
!        &             ggij(4,iatom,ichan,iwchan)+&
!        &             ggij(5,iatom,ichan,iwchan)+&
!        &             ggij(6,iatom,ichan,iwchan), iatom,iwchan
!        end do
!        end do
         if (ngrids==1) then
           if (anaddb_dtset%dossum /= 0 ) then
             write(message,'(a65,i5,a16)') ' DOS, SUM AND DIFFERENCE OF TWO PHONON DOS with channel width= ',iwchan,':'
           else
             write(message,'(a25,i5,a16)') ' DOS  with channel width=  ',iwchan,':'
           end if
         else
           if (anaddb_dtset%dossum /= 0 ) then
             write(message,'(a65,i5,a16)') ' DOS, SUM AND DIFFERENCE OF TWO PHONON DOS  with channel width= ', &
&             iwchan,' newly converged'
           else
             write(message,'(a25,i5,a16)') ' DOS  with channel width=  ',&
&             iwchan,' newly converged'
           end if
         end if
         
         call wrtout(std_out,message,'COLL')
         call wrtout(iout,message,'COLL')
         do ichan=1,nchan2(iwchan)
           if (anaddb_dtset%dossum /= 0 ) then
             write(message,'(i8,f11.1,3(f12.5,3x))') ichan,gg(ichan,iwchan)+tol10,&
&             gdos(ichan,iwchan)+tol10, gg_sum(ichan,iwchan)*(3.0*natom*3.0*natom)+tol10, &
&             gg_rest(ichan,iwchan)*(3.0*natom*(3.0*natom-1))+tol10
           else
             write(message,'(i8,f11.1,1x,f12.5)') ichan,gg(ichan,iwchan)+tol10,&
&             gdos(ichan,iwchan)+tol10
           end if
           call wrtout(std_out,message,'COLL')
         end do

         if (ngrids>1) then
           write(message,'(a24,f10.5)')'   with maximal error = ',gerr+tol10
           call wrtout(std_out,message,'COLL')
           call wrtout(iout,message,'COLL')
         end if

         do itemper=1,ntemper

           tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
!          The temperature (tmp) is given in Kelvin

!          Put zeroes for F, E, S, Cv
           free(itemper)=zero
           energy(itemper)=zero
           entropy(itemper)=zero
           spheat(itemper)=zero
           if (iavfrq>0) wme(itemper)=zero

           dosinc=dble(iwchan)

           do ichan=1,nchan2(iwchan)

!            wovert= hbar*w / 2kT dimensionless
             wovert=dosinc*(dble(ichan)-0.5_dp)/Ha_cmm1/(2._dp*kb_HaK*tmp)
             expm2x=exp(-2.0_dp*wovert)
             ln2shx=wovert+log(1.0_dp-expm2x)
             cothx=(1.0_dp+expm2x)/(1.0_dp-expm2x)
             factor=dble(3*natom)*gdos(ichan,iwchan)
             if (iavfrq>0) factorw=3*natom*gw(ichan,iwchan)

!            This matches the equations published in
!            Lee & Gonze, PRB 51, 8610 (1995)
             free(itemper)=free(itemper) +factor*kb_HaK*tmp*ln2shx
             energy(itemper)=energy(itemper) +factor*kb_HaK*tmp*wovert*cothx
             entropy(itemper)=entropy(itemper)&
&             +factor*kb_HaK*(wovert*cothx - ln2shx)

!            The contribution is much lower than 1.0d-16 when wovert<100.0_dp
             if(wovert<100.0_dp)then
               spheat(itemper)=spheat(itemper)&
&               +factor*kb_HaK*wovert**2/sinh(wovert)**2
             end if
             if (iavfrq>0) wme(itemper)=wme(itemper)+factorw*kb_HaK*wovert**2/sinh(wovert)**2

           end do

           if (iavfrq>0.and.abs(spheat(itemper))>tol8) wme(itemper)=wme(itemper)/spheat(itemper)
         end do

!        Check if the thermodynamic functions change within tolerance,

         if (ngrids>1) then
           write(message,'(a,a,a)')&
&           ' thm9 : check if the thermodynamic functions',ch10,&
&           '    change within tolerance.'
           call wrtout(std_out,message,'COLL')
           convth=1

           do itemper=1,ntemper
             change=free(itemper)-free0(itemper)
             relchg=change/free0(itemper)
             if(change>1d-14*dble(mqpt2) .and. relchg>anaddb_dtset%thmtol)then
               write(message,'(a,es14.4,a,a,es14.4)' )&
&               ' thm9 : free energy relative changes ',relchg,ch10,&
&               '        are larger than thmtol ',anaddb_dtset%thmtol
               call wrtout(std_out,message,'COLL')
               convth=0
             end if
             change=energy(itemper)-energy0(itemper)
             relchg=change/energy0(itemper)
             if(change>1d-14*dble(mqpt2) .and. relchg>anaddb_dtset%thmtol)then
               write(message,'(a,es14.4,a,a,es14.4)' )&
&               ' thm9 : energy relative changes ',relchg,ch10,&
&               '        are larger than thmtol ',anaddb_dtset%thmtol
               call wrtout(std_out,message,'COLL')
               convth=0
             end if
             change=entropy(itemper)-entropy0(itemper)
             relchg=change/entropy0(itemper)
             if(change>1d-14*dble(mqpt2) .and. relchg>anaddb_dtset%thmtol)then
               write(message,'(a,es14.4,a,a,es14.4)' )&
&               ' thm9 : entropy relative changes ',relchg,ch10,&
&               '        are larger than thmtol ',anaddb_dtset%thmtol
               call wrtout(std_out,message,'COLL')
               convth=0
             end if
             change=spheat(itemper)-spheat0(itemper)
             relchg=change/spheat0(itemper)
             if(change>1d-14*dble(mqpt2) .and. relchg>anaddb_dtset%thmtol)then
               write(message,'(a,es14.4,a,a,es14.4)' )&
&               ' thm9 : specific heat relative changes ',relchg,ch10,&
&               '        are larger than thmtol ',anaddb_dtset%thmtol
               call wrtout(std_out,message,'COLL')
               convth=0
             end if

             if(convth==0)exit

!            End of check if the thermodynamic functions change within tolerance
           end do

         else
           convth=1
         end if

!        Update F,E,S,C and eventually write them if converged

         if(convth==1)then
           part1=.true.
           if (iavfrq>0) then
             write(message,'(a,a,a)') ch10,&
&             ' # At  T   F(J/mol-c)     E(J/mol-c) ',&
&             '    S(J/(mol-c.K)) C(J/(mol-c.K)) Omega_mean(cm-1)'
           else
             write(message,'(a,a,a)') ch10,&
&             ' # At  T   F(J/mol-c)     E(J/mol-c) ',&
&             '    S(J/(mol-c.K)) C(J/(mol-c.K))'
           end if
           call wrtout(iout,message,'COLL')
           call wrtout(thermal_unit,message,'COLL')
           write(message,'(a)')&
&           ' # (A mol-c is the abbreviation of a mole-cell, that is, the'
           call wrtout(iout,message,'COLL')
           call wrtout(thermal_unit,message,'COLL')
           write(message,'(a)')&
&           ' #  number of Avogadro times the atoms in a unit cell)'
           call wrtout(iout,message,'COLL')
           call wrtout(thermal_unit,message,'COLL')

           write(message, '(a,a,a)' )&
&           ' thm9 : thermodynamic functions have converged',ch10,&
&           '     see main output file ...'
           call wrtout(std_out,message,'COLL')

         end if

         do itemper=1,ntemper

           free0(itemper)=free(itemper)
           energy0(itemper)=energy(itemper)
           entropy0(itemper)=entropy(itemper)
           spheat0(itemper)=spheat(itemper)

           if(convth==1)then
             tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
             if (iavfrq>0) then
               write(message,'(f6.1,5es15.7)') tmp+tol8,&
&               Ha_eV*e_Cb*Avogadro*free(itemper),&
&               Ha_eV*e_Cb*Avogadro*energy(itemper),&
&               Ha_eV*e_Cb*Avogadro*entropy(itemper),&
&               Ha_eV*e_Cb*Avogadro*spheat(itemper),&
&               wme(itemper)
             else
               write(message,'(f6.1,4es15.7)') tmp+tol8,&
&               Ha_eV*e_Cb*Avogadro*free(itemper),&
&               Ha_eV*e_Cb*Avogadro*energy(itemper),&
&               Ha_eV*e_Cb*Avogadro*entropy(itemper),&
&               Ha_eV*e_Cb*Avogadro*spheat(itemper)
             end if
             call wrtout(iout,message,'COLL')
             call wrtout(thermal_unit,message,'COLL')
           end if

         end do

       end if
       if(convth==1)exit

     end do
   end if

   if(.not.part2)then

!    Atomic temperature factor calculation

     do iwchan=nwchan,1,-1
       if (wgijcnv(iwchan))cycle

!      Calculate gij(k|w) and Bij(k)

!      Check if the density of states changed by more than dostol
       gijsum =zero
       wgijcnv(iwchan)=.true.
       if (ngrids>1) then
         do iatom=1,natom
           do ij=1,6
             gijerr=zero
             do ichan=1,nchan2(iwchan)
               gijsum = gijsum + gij(ij,iatom,ichan,iwchan)
               gijerr=gijerr&
&               +abs(ggij(ij,iatom,ichan,iwchan)/gnorm&
&               -gij(ij,iatom,ichan,iwchan))
             end do
             if(gijerr>anaddb_dtset%dostol) then
               wgijcnv(iwchan)=.false.
               exit
             end if
           end do
         end do
       else
         gijerr=0.d0
       end if

!      gij(k|w) is updated

       do ichan=1,nchan2(iwchan)
         do iatom=1,natom
           do ij=1,6
             gij(ij,iatom,ichan,iwchan)=ggij(ij,iatom,ichan,iwchan)/(gnorm/(3*natom))
           end do
         end do
       end do

!      Write gijerr for each q sampling and w width

       write(message,'(a,a,i3,3i6,f10.5,f10.5)') ch10,&
&       ' iwchan,igqpt(i),gijsum, gij error= ',&
&       iwchan,igqpt2(1),igqpt2(2),igqpt2(3),gijsum,gijerr+tol10
       call wrtout(std_out,message,'COLL')
!      call wrtout(iout,message,'COLL')

!      If the generalized DOS with a channel width is newly converged,
!      print it out and calculate Bij(k).

       if(wgijcnv(iwchan)) then


         if (ngrids==1) then
           write(message,'(a,i5,a)') ' gij with channel width=  ',iwchan,':'
         else
           write(message,'(a,i5,a)') ' gij with channel width=  ',&
&           iwchan,' newly converged'
         end if
         call wrtout(iout,message,'COLL')
         write(message,'(a,2i3,3i6,f10.5)') &
&         'iatom,iwchan,igqpt2(i),gij error= ',&
&         iatom,iwchan,igqpt2(1),igqpt2(2),igqpt2(3),gijerr+tol10
         call wrtout(iout,message,'COLL')

         do itemper=1,ntemper

           tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
!          tmp in K

!          Put zeroes for Bij(k)
           do iatom=1,natom
             do ij=1,6
               bbij(ij,iatom,itemper)=0._dp
               vij(ij,iatom,itemper)=0._dp
             end do
           end do

           dosinc=dble(iwchan)
!          
           do ichan=1,nchan2(iwchan)
!            
!            $wovert= \hbar*w / 2kT$, dimensionless
             wovert=dosinc*(dble(ichan)-0.5_dp)/Ha_cmm1/(2.*kb_HaK*tmp)
             expm2x=exp(-2.0_dp*wovert)
             do iatom=1,natom
               factor=Ha_cmm1/(2.*dosinc*(dble(ichan)-0.5))    &
&               *(1.0_dp+expm2x)/(1.0_dp-expm2x) &
&               /amu(typat(iatom))/amu_emass
               factorv=(0.5*dosinc*(float(ichan)-0.5)/Ha_cmm1)    &
&               *(1.0_dp+expm2x)/(1.0_dp-expm2x) &
&               /amu(typat(iatom))/amu_emass

               do ij=1,6
                 bbij(ij,iatom,itemper)=bbij(ij,iatom,itemper)&
&                 +factor*gij(ij,iatom,ichan,iwchan)
                 vij(ij,iatom,itemper)=vij(ij,iatom,itemper)&
&                 +factorv*gij(ij,iatom,ichan,iwchan)
               end do
             end do

           end do

         end do

!        B matrix is now in atomic unit in the Cartesian coordinates.
!        Check if Bij(k) changed within tolerance.
         convth=1
         if (ngrids>1) then
           do itemper=1,ntemper
             do iatom=1,natom
               do ij=1,6
                 diffbb=bbij(ij,iatom,itemper)-bij(ij,iatom,itemper)
                 if( diffbb                    > 1d-10  .and. &
&                 diffbb/bij(ij,iatom,itemper) > anaddb_dtset%thmtol         ) then
                   write(message,'(a)' )' thm9 : Bij changes are larger than thmtol '
                   call wrtout(std_out,message,'COLL')
                   convth=0
                 end if
                 if(convth==0)exit
               end do
               if(convth==0)exit
             end do
             if(convth==0)exit
           end do
         end if

!        Update Bij(k) and write them.
!        B matrix printed in angstrom^2

         if(convth==1)then
           write(message, '(a,a,a)' )&
&           ' B matrix elements as a function of T',ch10,&
&           '    Angstrom^2, cartesian coordinates'
           call wrtout(std_out,message,'COLL')
           call wrtout(iout,message,'COLL')
         end if

         do itemper=1,ntemper

!          tmp in K
           tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
           do iatom=1,natom

             do ij=1,6
               bij(ij,iatom,itemper)=bbij(ij,iatom,itemper)
             end do

             if(convth==1)then
               write(iout,'(2i3,f6.1,6es11.4)')&
&               iwchan,iatom,tmp+tol10,&
&               Bohr_Ang**2*bij(1,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(2,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(3,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(4,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(5,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(6,iatom,itemper)+tol10
             end if

           end do ! end loop over natom
         end do ! end loop over ntemper

!        Mean square velocity matrix printed in angstrom^2/sec^2
         if(convth==1)then
           write(message, '(a,a,a)' )&
&           ' <vel^2> matrix elements as a function of T',ch10,&
&           '    Angstrom^2/sec^2, cartesian coordinates'
           call wrtout(std_out,message,'COLL')
           call wrtout(iout,message,'COLL')
         end if
         do itemper=1,ntemper
!          tmp in K
           tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*float(itemper-1)
           do iatom=1,natom
             if(convth==1)then
               vij(:,iatom,itemper)=Bohr_Ang**2*vij(:,iatom,itemper)/Time_Sec**2

!              The following check zeros out <v^2> if it is very small, in order to 
!              avoid numerical noise being interpreted by the automatic tests as
!              something real. 1.0D12 only looks like a big number, for <v^2> expressed
!              in Angstrom^2/sec^2, this is very small. Note also that we compare it in
!              absolute value, that's because if any of the phonon frequencies are
!              computed as negative, <v^2> can take a negative value.      
               do icomp=1, 6
                 if (abs(vij(icomp,iatom,itemper)) < 1.0D12) vij(icomp,iatom,itemper)=tol10
               end do
               write(iout,'(2i3,f6.1,6es11.4)')&
&               iwchan,iatom,tmp+tol10,&
&               vij(1,iatom,itemper),&
&               vij(2,iatom,itemper),&
&               vij(3,iatom,itemper),&
&               vij(4,iatom,itemper),&
&               vij(5,iatom,itemper),&
&               vij(6,iatom,itemper)
             end if ! end check on convergence
           end do ! end loop over natom
         end do ! end loop over ntemper

         if(convth==1)part2=.true.

!        End of test on wgijcnv
       end if

!      End of loop over iwchan
     end do

!    End of part2
   end if

   if(part1.and.part2)exit

   call timein(tcpu,twall)
   write(message, '(a,f11.4,a,f11.4,a)' ) &
&   ' tcpu',tcpu-tcpui,' twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')

   ABI_DEALLOCATE(indqpt1)
   ABI_DEALLOCATE(qpt2)
   ABI_DEALLOCATE(spqpt2)
   ABI_DEALLOCATE(wtq)
   ABI_DEALLOCATE(wtq2)

!  End of the Loop on the q point grids
 end do

 ABI_DEALLOCATE(bbij)
 ABI_DEALLOCATE(bij)
 ABI_DEALLOCATE(energy)
 ABI_DEALLOCATE(energy0)
 ABI_DEALLOCATE(entropy)
 ABI_DEALLOCATE(entropy0)
 ABI_DEALLOCATE(free)
 ABI_DEALLOCATE(free0)
 ABI_DEALLOCATE(gdos)
 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(gg_sum)
 ABI_DEALLOCATE(gg_rest)
 ABI_DEALLOCATE(ggij)
 ABI_DEALLOCATE(gij)
 ABI_DEALLOCATE(nchan2)
 ABI_DEALLOCATE(spheat)
 ABI_DEALLOCATE(spheat0)
 ABI_DEALLOCATE(wgcnv)
 ABI_DEALLOCATE(wgijcnv)
 if(allocated(wtq))ABI_DEALLOCATE(wtq)
 if(allocated(wtq2))ABI_DEALLOCATE(wtq2)
 if (iavfrq>0)  then
   ABI_DEALLOCATE(gw)
   ABI_DEALLOCATE(wme)
 end if

 if(.not.part1)then
   write(message, '(a,a,a,a,a,a,a,a,a,a,a)' )&
&   ' thm9 : ERROR -',ch10,&
&   '  No thermodynamical function is printed out :',ch10,&
&   '  the tolerance level that was asked ',ch10,&
&   '  has not been match with the grids specified.',ch10,&
&   '  Action : in the input file, increase the resolution',ch10,&
&   '  of grids ng2qpt, or decrease the accuracy requirement thmtol.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(.not.part2)then
   write(message, '(a,a,a,a,a,a,a,a,a,a,a)' )&
&   ' thm9 : WARNING -',ch10,&
&   '  No atomic factor tensor is printed out :',ch10,&
&   '  the tolerance level that was asked ',ch10,&
&   '  has not been match with the grids specified.',ch10,&
&   '  Action : in the input file, increase the resolution',ch10,&
&   '  of grids ng2qpt, or decrease the accuracy requirement thmtol.'
   call wrtout(std_out,message,'COLL')
 end if


 close (thermal_unit)

end subroutine thm9
!!***
