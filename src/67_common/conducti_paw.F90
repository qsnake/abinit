!{\src2tex{textfont=tt}}
!!****f* ABINIT/conducti_paw
!! NAME
!! conducti_paw
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (VRecoules, PGhosh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!!  bantot
!!  doccde(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy.
!!  dom=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  eigen11(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction x
!!  eigen12(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction y
!!  eigen13(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction z
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  fermie= fermi energy (Hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  kin11= Onsager kinetic coeficient=optical conductivity
!!  kin12= Onsager kinetic coeficient
!!  kin21= Onsager kinetic coeficient
!!  kin22= Onsager kinetic coeficient
!!  Kth=thermal conductivity
!!  mom=number of frequency for conductivity computation
!!  mband=maximum number of bands.
!!  natom = number of atoms in the unit cell.
!!  nband(nkpt*nsppol)=number of bands at each RF k point for each spin.
!!  nkpt=number of k points in the IBZ for this perturbation
!!  ngfft(3)=integer fft box dimensions.
!!  nspinor=number of spinorial components of the wavefunctions.
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  ntypat = number of atom types.
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k.
!!  occopt==option for occupancies
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).sigx(mom,nphicor))
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  Sth=thermopower
!!  tsmear=smearing width (or temperature) in Hartree
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  wind=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!  znucl(natom)=atomic number of atoms
!!  np_sum=noziere-pines sumrule
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      hdr_clean,hdr_io,metric,msig,wffclose,wffopen,wffreadeigk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine conducti_paw(filnam,filnam_out,mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile

 use m_header,       only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'conducti_paw'
 use interfaces_42_geometry
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_67_common, except_this_one => conducti_paw
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=fnlen) :: filnam,filnam_out
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: accesswff,bantot,bdtot0_index,bdtot_index,dosdeltae
 integer :: fform0,fform1,formeig0,headform,iband,ierr,ikpt
 integer :: iom,isppol,jband,l1,l2,master,mband,me,mom
 integer :: natom,nband1,nband_k,nkpt,nspinor,nsppol,ntypat
 integer :: occopt,rdwr,spaceComm,tim_rwwf
 real(dp) :: del,deltae,diff_occ,ecut,fermie,maxocc
 real(dp) :: np_sum,np_sum_k1,np_sum_k2,omin,omax,dom,oml,sig,socc,socc_k
 real(dp) :: Tatm,tphysel,tsmear,ucvol
 character(len=fnlen) :: filnam0,filnam1,filnam_gen
 character(len=6) :: codvsn
 type(hdr_type) :: hdr
 type(wffile_type) :: wff0,wff1
!arrays
 integer,allocatable :: nband(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: cond_nd(:,:,:),dhdk2_r(:,:,:,:),dhdk2_g(:,:,:)
 real(dp),allocatable :: doccde(:),doccde_k(:),eig0_k(:),eig0tmp(:),eigen0(:)
 real(dp),allocatable :: eigtmp(:),occ(:),occ_k(:),wtk(:),oml1(:)
 real(dp),allocatable :: kin11(:,:),kin12(:),kin21(:),kin22(:)
 real(dp),allocatable :: kin11_k(:),kin12_k(:),kin21_k(:),kin22_k(:),Kth(:),Stp(:)
 real(dp),allocatable :: psinablapsi(:,:,:,:),sig_abs(:)

! *********************************************************************************
!BEGIN EXECUTABLE SECTION

!write(std_out,'(a)')' The name of the output file is :',trim(filnam_out)
!Read data file
 open(15,file=filnam,form='formatted')
 rewind(15)
 read(15,*)
 read(15,'(a)')filnam_gen       ! generic name for the files
 filnam1=trim(filnam_gen)//'_OPT'
!Read size of the frequency range
 read(15,*) dom,omin,omax,mom
 close(15)
 write(std_out,'(a,i8,3f10.5,a)')' npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'

!Open the Wavefunction and optic files
!These default values are typical of sequential use
 accesswff=IO_MODE_FORTRAN ; spaceComm=abinit_comm_serial ; master=0 ; me=0
 call WffOpen(accesswff,spaceComm,filnam1,ierr,wff1,master,me,11)
 read(11,iostat=ierr)codvsn,headform,fform1
 if ((ierr /=0).or.((fform1/=610).and.(fform1/=612))) then
   write(std_out,*)'format prior version 6.1'
   fform1=613
 else if (fform1==612) then
   write(std_out,*)'Optic for local potential only'
 end if
 call WffClose(wff1,ierr)

!Open the conducti and/or optic files
 call WffOpen(accesswff,spaceComm,filnam1,ierr,wff1,master,me,11)

!Read the header from Ground state file
 rdwr=1
 if (fform1==613) then
   filnam0=trim(filnam_gen)//'_WFK'
   call WffOpen(accesswff,spaceComm,filnam0,ierr,wff0,master,me,10)
   call hdr_io(fform0,hdr,rdwr,wff0)
 else 
   call hdr_io(fform1,hdr,rdwr,wff1)
 end if

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 ecut=hdr%ecut_eff
 natom=hdr%natom
 nkpt=hdr%nkpt
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 ntypat=hdr%ntypat
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 ABI_ALLOCATE(nband,(nkpt*nsppol))
 ABI_ALLOCATE(occ,(bantot))
 ABI_ALLOCATE(wtk,(nkpt))
 fermie=hdr%fermie
 tsmear=hdr%tsmear
 occ(1:bantot)=hdr%occ(1:bantot)
 wtk(1:nkpt)=hdr%wtk(1:nkpt)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
 write(std_out,'(a,i8)')       ' natom             =',natom
 write(std_out,'(a,3i8)')      ' nkpt,mband,nsppol        =',nkpt,mband,nsppol
 write(std_out, '(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'
 Tatm=tsmear*Ha_K
 write(std_out,'(a,f12.5,a,f12.5,a)') ' Temp              =',tsmear,' Ha ',Tatm,' Kelvin'

 if (fform1==613) then
!  Prepare the reading of Wff files
   formeig0=0 ; tim_rwwf=0
   ABI_ALLOCATE(eigtmp,(2*mband*mband))
   ABI_ALLOCATE(eig0tmp,(mband))
!  Read the eigenvalues of ground-state
   ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
   bdtot0_index=0 ; bdtot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband1=nband(ikpt+(isppol-1)*nkpt)
       call WffReadEigK(eig0tmp,formeig0,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff0)
       eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)
       bdtot0_index=bdtot0_index+nband1
     end do
   end do
   call WffClose(wff0,ierr)
   ABI_DEALLOCATE(eig0tmp)
 else
   ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
   read(11)(eigen0(iband),iband=1,mband*nkpt*nsppol)
 end if
!
!
!---------------------------------------------------------------------------------
!gmet inversion to get ucvol
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!---------------------------------------------------------------------------------
!derivative of occupation wrt the energy.
 ABI_ALLOCATE(doccde,(mband*nkpt*nsppol))
 if (occopt==1) then
   write(std_out,'(a,i4)')  ' occopt            =',occopt
   doccde=zero
 else
   tphysel=zero
   maxocc=two/(nsppol*nspinor)
   dosdeltae=zero
 end if
!---------------------------------------------------------------------------------
!size of the frequency range
 del=(omax-omin)/(mom-1)
 ABI_ALLOCATE(oml1,(mom))
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do
 ABI_ALLOCATE(cond_nd,(mom,3,3))
 ABI_ALLOCATE(kin11,(mom,nsppol))
 ABI_ALLOCATE(kin12,(mom))
 ABI_ALLOCATE(kin21,(mom))
 ABI_ALLOCATE(kin22,(mom))
 ABI_ALLOCATE(sig_abs,(mom))
 ABI_ALLOCATE(kin11_k,(mom))
 ABI_ALLOCATE(kin12_k,(mom))
 ABI_ALLOCATE(kin21_k,(mom))
 ABI_ALLOCATE(kin22_k,(mom))
 ABI_ALLOCATE(Kth,(mom))
 ABI_ALLOCATE(Stp,(mom))
 
!---------------------------------------------------------------------------------
!Conductivity -------
!
 ABI_ALLOCATE(psinablapsi,(2,3,mband,mband))
 kin11   = zero
 kin12   = zero
 kin21   = zero
 kin22   = zero
 np_sum  = zero
 socc    = zero
 sig_abs = zero

 bdtot_index = 0

!LOOP OVER SPINS/K
 deltae  = zero
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(doccde_k,(nband_k))
     ABI_ALLOCATE(dhdk2_r,(nband_k,nband_k,3,3))
     ABI_ALLOCATE(dhdk2_g,(natom,nband_k,nband_k))

     cond_nd   = zero
     kin11_k   = zero
     kin12_k   = zero
     kin21_k   = zero
     kin22_k   = zero
     np_sum_k1 = zero
     np_sum_k2 = zero
     socc_k    = zero
     dhdk2_r   = zero
     dhdk2_g   = zero

!    eigenvalue for k-point
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!    first derivative eigenvalues for k-point
     psinablapsi=zero
     read(11)((psinablapsi(1:2,1,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(11)((psinablapsi(1:2,2,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(11)((psinablapsi(1:2,3,iband,jband),iband=1,nband_k),jband=1,nband_k) 
!    DEBUG
!    write(963,*)isppol,ikpt,((psinablapsi(1:2,1,iband,jband),iband=1,nband_k),jband=1,nband_k)
!    write(963,*)isppol,ikpt,((psinablapsi(1:2,2,iband,jband),iband=1,nband_k),jband=1,nband_k)
!    write(963,*)isppol,ikpt,((psinablapsi(1:2,3,iband,jband),iband=1,nband_k),jband=1,nband_k) 
!    ENDDEBUG

!    occupation numbers for k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    derivative of occupation number for k-point
     doccde_k(:)=doccde(1+bdtot_index:nband_k+bdtot_index)

!    LOOP OVER BANDS
     do iband=1,nband_k
       do jband=1,nband_k
         do l1=1,3
           do l2=1,3
             dhdk2_r(iband,jband,l1,l2)=dhdk2_r(iband,jband,l1,l2)+(&
&             psinablapsi(1,l1,iband,jband)*psinablapsi(1,l2,iband,jband)&
&             +psinablapsi(2,l1,iband,jband)*psinablapsi(2,l2,iband,jband))
           end do
         end do

         do l1=1,3
           dhdk2_g(1,iband,jband)=dhdk2_g(1,iband,jband)+( &
&           psinablapsi(1,l1,iband,jband)*psinablapsi(1,l1,iband,jband) &
&           +psinablapsi(2,l1,iband,jband)*psinablapsi(2,l1,iband,jband))
         end do

         diff_occ = occ_k(iband)-occ_k(jband)
         if (dabs(diff_occ)>=tol8) then

!          Conductivity for each omega
!          omin = zero
           do iom=1,mom
             oml=oml1(iom)
             if (jband>iband) then
               sig= dhdk2_g(1,iband,jband)&
&               *(diff_occ)/oml*(dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)&
&               -dexp(-((eig0_k(iband)-eig0_k(jband)-oml)/dom)**2))
               kin11_k(iom)=kin11_k(iom)+sig
               kin12_k(iom)=kin12_k(iom)-sig*(eig0_k(jband)-fermie)
               kin21_k(iom)=kin21_k(iom)-sig*(eig0_k(iband)-fermie)
               kin22_k(iom)=kin22_k(iom) + &
&               sig*(eig0_k(iband)-fermie)*(eig0_k(jband)-fermie)
             end if
             do l1=1,3
               do l2=1,3
                 cond_nd(iom,l1,l2)=cond_nd(iom,l1,l2) +dhdk2_r(iband,jband,l1,l2)&
&                 *(diff_occ)/oml*dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)
               end do
             end do
           end do

!          Sumrule start
           if (dabs(eig0_k(iband)-eig0_k(jband))>=tol10) then
             np_sum_k1=np_sum_k1 -dhdk2_g(1,iband,jband)&
&             *(diff_occ)/(eig0_k(iband)-eig0_k(jband))
           else
             np_sum_k2=np_sum_k2 - doccde_k(iband)*dhdk2_g(1,iband,jband)
           end if

!          end loop over band
         end if
       end do
       socc_k=socc_k+occ_k(iband)
     end do
     do iom=1,mom
       kin11(iom,isppol)=kin11(iom,isppol)+wtk(ikpt)*kin11_k(iom)
       kin12(iom)=kin12(iom)+wtk(ikpt)*kin12_k(iom)
       kin21(iom)=kin21(iom)+wtk(ikpt)*kin21_k(iom)
       kin22(iom)=kin22(iom)+wtk(ikpt)*kin22_k(iom)
     end do
     np_sum=np_sum + wtk(ikpt)*(np_sum_k1+np_sum_k2)
     socc=socc+wtk(ikpt)*socc_k

!    Validity limit
     deltae=deltae+(eig0_k(nband_k)-fermie)

     bdtot_index=bdtot_index+nband_k
     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(doccde_k)
     ABI_DEALLOCATE(dhdk2_r)
     ABI_DEALLOCATE(dhdk2_g)
!    End loop over k
   end do
!  End loop over Spin
 end do

 write(std_out,'(a,3f10.5)')' sumrule           =',np_sum/socc/three/dble(nsppol),socc
 write(std_out,'(a,f10.5,a,f10.5,a)')&
& ' Emax-Efermi       =',deltae/dble(nkpt*nsppol),' Ha',deltae/dble(nkpt*nsppol)*Ha_eV,' eV'


 open(20,file=trim(filnam_out)//'_Lij',form='formatted')
 write(20,'(a)')' # omega(ua) L12 L21 L22 L22'
 open(30,file=trim(filnam_out)//'_sig',form='formatted')
 if (nsppol==1) then
   write(30,'(a)')' # omega(ua) hbar*omega(eV)    cond(ua)             cond(ohm.cm)-1'
 else
   write(30,'(2a)')' # omega(ua) hbar*omega(eV)      cond(ua)            cond(ohm.cm)-1',&
&   '      cond(ohm.cm)-1 UP      cond(ohm.cm)-1 DN'
 end if
 open(41,file=trim(filnam_out)//'_Kth',form='formatted')
 write(41,'(a)')&
& ' #omega(ua) hbar*omega(eV)  thermal cond(ua)   Kth(W/m/K)   thermopower(ua)   Stp(microohm/K)'
 open(45,file=trim(filnam_out)//'.out',form='formatted')
 write(45,'(a)' )' #Conducti output file:'
 write(45,'(a)' )' #Contains all results produced by conducti utility'
 write(45,'(a)' )' '
 write(45,'(a)')' # omega(ua)       cond(ua)             thermal cond(ua)       thermopower(ua)'

!call isfile(filnam_out,'new')

!Compute thermal conductivity and thermopower
 do iom=1,mom
   oml=oml1(iom)
   do isppol=1,nsppol
     kin11(iom,isppol)=kin11(iom,isppol)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
     if (dabs(kin11(iom,isppol))<10.0d-20) kin11(iom,isppol)=zero
     sig_abs(iom)=sig_abs(iom)+kin11(iom,isppol)
   end do
   kin21(iom)=kin21(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin12(iom)=kin12(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin22(iom)=kin22(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   Kth(iom)=kin22(iom)
   Stp(iom)=zero
   if(sig_abs(iom)/=zero)  then
     Kth(iom)=Kth(iom)-(kin12(iom)*kin21(iom)/sig_abs(iom))
     Stp(iom)=kin12(iom)/(sig_abs(iom)*Tatm)
   end if
   if (dabs(Kth(iom))<10.0d-20) Kth(iom)=zero
   if (dabs(Stp(iom))<10.0d-20) Stp(iom)=zero
   write(20,'(f12.5,4es22.12)')oml,kin12(iom),kin21(iom),kin22(iom),kin22(iom)/Tatm*3.4057d9
   if (nsppol==1) then
     write(30,'(2f12.5,2es22.12)') oml,oml*Ha_eV,sig_abs(iom),sig_abs(iom)*Ohmcm
   else
     write(30,'(2f12.5,4es22.12)') oml,oml*Ha_eV,sig_abs(iom),sig_abs(iom)*Ohmcm,&
&     kin11(iom,1)*Ohmcm,kin11(iom,2)*Ohmcm
   end if
   write(41,'(2f12.5,4es22.12)') oml,oml*Ha_eV,Kth(iom),Kth(iom)*3.4057d9/Tatm,&
&   Stp(iom),Stp(iom)*3.6753d-2
   write(45,'(1f12.5,3es22.12)') oml,sig_abs(iom),Kth(iom),Stp(iom)
 end do

!Calculate the imaginary part of the conductivity (principal value)
!+derived optical properties.
 call msig (sig_abs,mom,oml1,filnam_out)

 ABI_DEALLOCATE(psinablapsi)
 ABI_DEALLOCATE(kin11)
 ABI_DEALLOCATE(kin22)
 ABI_DEALLOCATE(kin12)
 ABI_DEALLOCATE(kin21)
 ABI_DEALLOCATE(kin11_k)
 ABI_DEALLOCATE(kin22_k)
 ABI_DEALLOCATE(kin12_k)
 ABI_DEALLOCATE(kin21_k)
 ABI_DEALLOCATE(Stp)
 ABI_DEALLOCATE(Kth)
 ABI_DEALLOCATE(cond_nd)
 ABI_DEALLOCATE(sig_abs)
 close(15);close(20);close(30)
 close(41);close(45)
 call WffClose(wff1,ierr)
 write(std_out,'(2a)')ch10,'OUTPUT'
 write(std_out,'(a)')trim(filnam_out)//'_Lij : Onsager kinetic coefficients'
 write(std_out,'(a)')trim(filnam_out)//'_sig : Optical conductivity'
 write(std_out,'(a)')trim(filnam_out)//'_Kth : Thermal conductivity and thermopower'
 write(std_out,'(a)')trim(filnam_out)//'_eps : Dielectric fonction'
 write(std_out,'(a)')trim(filnam_out)//'_abs : n, k, reflectivity, absorption'
 

 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(oml1)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(wtk)
 call hdr_clean(hdr)

end subroutine conducti_paw
!!***
