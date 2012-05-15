!{\src2tex{textfont=tt}}
!!****f* ABINIT/rdm
!! NAME
!! rdm
!!
!! FUNCTION
!!  Main subroutine for reduced density matrix calculations
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (NH, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  acell(3)=length scales of primitive translations (bohr)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  mpi_enreg=informations about MPI parallelization (to be completed)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  MG: presently not used for rdm calculations, it is a zero-sized structure passed to rdkss
!!  rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      driver
!!
!! NOTES
!!
!! CHILDREN
!!      cigfft,clcqpg,cprj_alloc,cprj_free,crho,cvxclda,destroy_crystal,fermi
!!      fftwfn,findnq,findq,findqg0,fourdp,hartre,hdr_clean,hdr_vs_dtset,identk
!!      identq,init_crystal_from_hdr,lattice,leave_new,metric,mkrdim,occred
!!      old_setmesh,pclock,print_crystal,rho_tw_g,rotate_fft_mesh,setshells
!!      setup_g_rotation_old,testkss,wrtout,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rdm(acell,dtfil,dtset,pawtab,mpi_enreg,rprim)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_crystal_io
 use m_bz_mesh

 use m_blas,        only : xdotc
 use m_fft_mesh,    only : rotate_FFT_mesh, cigfft
 use m_header,      only : hdr_clean
 use m_io_kss,      only : testkss
 use m_oscillators, only : rho_tw_g

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdm'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_59_io_mpi
 use interfaces_93_rdm, except_this_one => rdm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
 character(len=50),parameter :: sub_name='rdm.F90'
!NH new variables
!(fcc)
!
!scalars
 integer,parameter :: iunit=10
 integer :: enforce_sym,i1,i2,i3,ib,ierr,ig,ig01,ig02,ig03,ii,ik,use_padfft
 integer :: il,iq,ir,is,istat,isym,jj,localrdwf,spaceComm
 integer :: method,mpsang,nG01d,nG02d,nG03d
 integer :: natom,nbnds_kss,nel,ng_kss,ninv,nkbzx !nbvw,
 integer :: nkibzR,npwvec,nqbzx,nr,nsppolR,nsym_out,ntypat,prtvol
 integer :: tim_fourdp,timrev
 real(dp),parameter :: bz_geometry_factor=7.44
 real(dp) :: bzvol,ct,efermi,etot,g1,g2,g3,gsq,gsqcut,i_sz,nelect
 real(dp) :: omegaplasma,q0vol,sqrt_i_sz,sumband1,sumband2,sumg,sumk1,sumspin
 real(dp) :: tolq0,ucvol,xcenergy
 logical :: iscompatibleFFT,use_antiferro
 character(len=500) :: message,msg
 type(Crystal_structure) :: Cryst
 type(Hdr_type) :: Hdr_kss
 type(rdm_parameters) :: rdmp
!arrays
 integer :: g0(3)
 integer,allocatable :: gbound(:,:)
 integer,allocatable :: counter(:),dimlmn(:),grottb(:,:,:),grottbm1(:,:,:)
 integer,allocatable :: gvec(:,:),irottb(:,:),ktab(:),ktabi(:),ktabo(:)
 integer,allocatable :: ktabr(:,:),nbv(:),qtab(:),qtabi(:),qtabo(:)
 integer,allocatable,target :: igfft(:,:,:,:)
 integer,pointer :: gvec_p(:,:),igfft0(:),igfftg0(:)
 real(dp) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),gmet(3,3),gprimd(3,3)
 real(dp) :: kpoint(3),qphon(3),rm1t(3),rmet(3,3),rprimd(3,3),spinrot1(4)
 real(dp) :: spinrot2(4)
 real(dp),allocatable :: delta1(:),delta2(:,:,:),deltardocc(:,:,:)
 real(dp),allocatable :: kbz(:,:),kibz(:,:),ldaen(:,:,:)
 real(dp),allocatable :: ldaocc(:,:,:)
 real(dp),allocatable :: occgrad(:,:,:),qbz(:,:),qpg(:,:)
 real(dp),allocatable :: qpoint(:,:),rdocc(:,:,:),rho(:,:)
 real(dp),allocatable :: rhog(:,:),vcorrlda(:,:)
 real(dp),allocatable :: vhartr(:)
 real(dp),allocatable :: vme(:,:,:),wtk(:),wtq(:)
 real(dp),pointer :: energies_p(:,:,:)
 complex(dpc),allocatable :: rdexpcoeff(:,:,:,:)
 complex(gwpc),allocatable :: ldawfg(:,:,:,:),ldawfr(:,:,:,:)
 complex(dpc),allocatable :: ktabp(:)
 complex(gwpc),allocatable :: rdwfr(:,:,:,:)
 complex(gwpc),allocatable :: rhotwg_ki(:)
 type(cprj_type),allocatable :: ldacprj(:,:)

!************************************************************************

 ABI_ALLOCATE(gbound,(1,0))

!Initialise parallelization related quantities
 mpi_enreg%paralbd=0
 mpi_enreg%paral_level=0

 call xcomm_world(mpi_enreg,spaceComm)

 prtvol=dtset%prtvol
 tolq0=0.001_dp  ! tolerance below which a q-point is treated as zero (long wavelength limit)

 write(message,'(4a)')ch10,' RDM: Calculation of total energies and band gaps of Silicon :)',ch10,ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

!Start clock
 call pclock(0)

!KSS file requires same number of bands at each kpoint
 rdmp%nbnds=dtset%nband(1)
 use_antiferro = (Dtset%nspden==2 .and. Dtset%nsppol==1)

 spinrot1(1)=one   ; spinrot2(1)=one
 spinrot1(2)=zero  ; spinrot2(2)=zero
 spinrot1(3)=zero  ; spinrot2(3)=zero
 spinrot1(4)=zero  ; spinrot2(4)=zero

!Compute dimensional primitive translations rprimd
 call mkrdim(acell,rprim,rprimd)

!Obtain dimensional translations in reciprocal space gprimd, metrics and unit cell volume, from rprimd.
!Also output rprimd, gprimd and ucvol
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)

!WARNIG here there is a possible undected error if symrel is not
!the same as the symrel read from the KSS file

!Define consistently npw, nsh, and ecut
 call setshells(dtset%ecutwfn,dtset%npwwfn,dtset%nshwfn,dtset%nsym,gmet,gprimd,dtset%symrel,'wfn',ucvol)

!the names of the variables should be changed to something like: number of G vectors in the sums
 call setshells(dtset%ecuteps,dtset%npweps,dtset%nsheps,dtset%nsym,gmet,gprimd,dtset%symrel,'mat',ucvol)

 rdmp%npwwfn=dtset%npwwfn
 rdmp%npwx=dtset%npweps

!Calculate npwvec as the max between dtset%npwsigx and dtset%npwwfn
!anyway npwwfn should be always larger that npwx
 npwvec=max(rdmp%npwwfn,rdmp%npwx)

 localrdwf=dtset%localrdwf        ! localrdwf==1 ==> all procs have access to files (default)
!localrdwf==0 ==> only master has access to files

 call testkss(dtfil%fnameabi_kss,Dtset%accesswff,nsym_out,nbnds_kss,ng_kss,mpsang,gvec_p,energies_p,&
& Hdr_kss,spaceComm)

 ABI_DEALLOCATE(energies_p)

 call hdr_vs_dtset(Hdr_kss,Dtset)

 timrev = 2 ! This information is not reported in the header
!1 --> do not use time-reversal symmetry
!2 --> take advantage of time-reversal symmetry
!Use time reversal symmetry
 ninv=2
 call init_crystal_from_hdr(Cryst,Hdr_kss,timrev,.FALSE.)
 call print_crystal(Cryst)
 if (nsym_out/=Cryst%nsym)  STOP 'nsym read from KSS differ from value read from heade, likely due to symmorphy'

!copy important dimension from header
 nsppolR=Hdr_kss%nsppol
 ntypat=Hdr_kss%ntypat
 natom=Hdr_kss%natom
 nkibzR=Hdr_kss%nkpt

 rdmp%nkibz=nkibzR
 rdmp%nsppol=nsppolR

 ABI_ALLOCATE(kibz,(3,rdmp%nkibz))
 kibz(:,:)=Hdr_kss%kptns(:,:)

 a1(:) =Hdr_kss%rprimd(:,1)  ! Fix problem with QPLDA
 a2(:) =Hdr_kss%rprimd(:,2)
 a3(:) =Hdr_kss%rprimd(:,3)
!Check whether the lattice from the input file agrees with that read from the KSS file
 if ( (ANY(ABS(a1-Dtset%rprimd_orig(:,1,1))>tol6)) .or. &
 (ANY(ABS(a2-Dtset%rprimd_orig(:,2,1))>tol6)) .or. &
 (ANY(ABS(a3-Dtset%rprimd_orig(:,3,1))>tol6)) ) then
   write(msg,'(6a)')ch10,&
&   ' rdm : ERROR - ',ch10,&
&   ' real lattice vectors read from the KSS file ',ch10,&
&   ' differ from the values specified in the input file'
   call wrtout(std_out,msg,'COLL')
   write(msg,'(3a,3(3es16.6),3a,3(3es16.6),3a)')ch10,&
&   ' rprimd from KSS file   = ',ch10,(Dtset%rprimd_orig(:,jj,1),jj=1,3),ch10,&
&   ' rprimd from input file = ',ch10,a1,a2,a3,ch10,ch10,&
&   ' Please, modify the lattice vectors in the input file '
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
 end if
!END COPY

 rdmp%nsym=Cryst%nsym

 if (npwvec>ng_kss) then
   write(message,'(3a)')&
&   ' rdm : WARNING - ',ch10,&
&   '  number of G-vectors found less then required'
   call wrtout(std_out,message,'COLL')

   npwvec=ng_kss
   if (rdmp%npwwfn>ng_kss) rdmp%npwwfn=ng_kss
   if (rdmp%npwx>ng_kss) rdmp%npwx=ng_kss
   write(message,'(3(a,i8,a))')&
&   '         calculation will proceed with npwvec = ',npwvec,ch10,  &
   '         calculation will proceed with npwsigx= ',rdmp%npwx,ch10,&
   '                                       npwwfn = ',rdmp%npwwfn,ch10
   call wrtout(std_out,message,'COLL')
 end if

 if (rdmp%nbnds>nbnds_kss) then
   write(message,'(5a,i4,a)')&
&   ' rdm : WARNING - ',ch10,                              &
&   '  number of bands found less then required     ',ch10,&
&   '         calculation will proceed with nbnds= ',nbnds_kss,ch10
   call wrtout(std_out,message,'COLL')
   rdmp%nbnds=nbnds_kss
 end if

!for parallel version, at the moment this code is not parallelized !
!min_band_proc=1
!max_band_proc=rdmp%nbnds
!nbnds_per_proc=rdmp%nbnds

!MG These quantities are used to allocate valence and conductions states in case of gwpara==2
!I dont know whether they could be useful in rdm or not
!NH They are not, there are no unoccupied states in rdm and nbvw is equal to rdmnb
!nbvw : maximum number of fully and partially occupied states over spin
!nbcw : maximum number of unoccupied states over spin
!nbvw = maxval(ibocc)
!nbcw = rdmp%nbnds-nbvw

!Allocate KS electronic structure variables

 ABI_ALLOCATE(gvec,(3,npwvec))
 istat = ABI_ALLOC_STAT
 gvec = gvec_p(:,1:npwvec)
 ABI_DEALLOCATE(gvec_p)
 if(istat/=0) stop 'out of memory in gvec'

 ABI_ALLOCATE(ldawfg,(rdmp%npwwfn,rdmp%nbnds,rdmp%nkibz,rdmp%nsppol))
 istat = ABI_ALLOC_STAT
 ABI_ALLOCATE(ldaocc,(rdmp%nkibz,rdmp%nbnds,rdmp%nsppol))
 ABI_ALLOCATE(ldaen,(rdmp%nkibz,rdmp%nbnds,rdmp%nsppol))

 ldawfg(:,:,:,:)=(0.,0.)
 ldaen(:,:,:)=zero
 ldaocc(:,:,:)=zero

!Read in KS band structure

 write(message,'(2a)')' rdm : will call rdkss ',ch10
 call wrtout(std_out,message,'COLL')

 if (dtset%usepaw==1) then
   ABI_ALLOCATE(ldacprj,(natom,dtset%nspinor*rdmp%nkibz*rdmp%nbnds*rdmp%nsppol))
   ABI_ALLOCATE(dimlmn,(natom))
   do ii=1,natom
     dimlmn(ii)=pawtab(dtset%typat(ii))%lmn_size
   end do
   call cprj_alloc(ldacprj,0,dimlmn)
   ABI_DEALLOCATE(dimlmn)
 end if

 stop "reading from KSS file in rmd is broken"
!energies are not read anymore, can use testkss but ordering is different.

!call rdkss(dtfil%fnameabi_kss,Dtset%usepaw,pawtab,Cryst%nsym,rdmp%nbnds,nbvw,rdmp%nkibz,Dtset%nspinor,rdmp%nsppol,&
!& rdmp%npwwfn,ldaen,ldaocc,ldawfg,ldacprj,ntypat,natom,mpsang,nelect,&
!& 1,rdmp%nbnds,Dtset%accesswff,Dtset%localrdwf,spaceComm)

 call pclock(10)

!Calculate b1, b2, b3 and ucvol, bzvol
!Here, unlike the main abinit code, the reciprocal vectors are defined as a_i \cdot b_j = 2pi delta_ij
 call lattice(a1,a2,a3,b1,b2,b3,ucvol,bzvol)

!Setup of the real space FFT mesh, note that dtset%ngfft is initialized
!before entering screening. Here we redefine dtset%ngfft(1:6) according to these options :
!enforce_sym==1 enforce a FFT mesh compatible with all the symmetry operation
!method==0 FFT grid read from fft.in (debugging purpose)
!method==1 normal FFT grid
!method==2 slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh)
!method==3 doubles FFT grid, to treat exactly the convolution defining the density

!MG Wed Oct 10 15:21:45 CEST 2007
!I ve substituted addshell with rdmp%mG0(3), this is useful to
!change on-the-fly the number of G0 shell without having to hack the
!code, there are other modifications in ciggfft and rho_tw_g

 rdmp%mG0(:)=2 !this is the default value but in GW I have a sub to find the optimal mG0 TODO
 enforce_sym=mod(dtset%fftgw,10)
 if (dtset%fftgw==00 .or. dtset%fftgw==01) method=0
 if (dtset%fftgw==10 .or. dtset%fftgw==11) method=1
 if (dtset%fftgw==20 .or. dtset%fftgw==21) method=2
 if (dtset%fftgw==30 .or. dtset%fftgw==31) method=3

 call old_setmesh(gmet,gvec,dtset%ngfft,npwvec,rdmp%npwx,rdmp%npwwfn,nr,&
& method,rdmp%mG0,Cryst%nsym,Cryst%symrel,Cryst%tnons,enforce_sym)

 write(std_out,*)' ngfft for rdm ',dtset%ngfft(:)

!Calculate igftt table for FFT
 nG01d=2*rdmp%mG0(1)+1
 nG02d=2*rdmp%mG0(2)+1
 nG03d=2*rdmp%mG0(3)+1
 ABI_ALLOCATE(igfft,(npwvec,nG01d,nG02d,nG03d))
 istat = ABI_ALLOC_STAT

 call cigfft(rdmp%mG0,npwvec,dtset%ngfft,gvec,igfft,ierr)
 if (ierr/=0) STOP "ierr in cigfft"

 igfft0 => igfft(:,rdmp%mG0(1)+1,rdmp%mG0(2)+1,rdmp%mG0(3)+1) ! FFT index of G for fourdp

 call pclock(20)


!Calculate KS wavefunctions in real space using FFT
 ABI_ALLOCATE(ldawfr,(nr,rdmp%nbnds,rdmp%nkibz,rdmp%nsppol))
 istat = ABI_ALLOC_STAT

 ldawfr(:,:,:,:)=(0.,0.)

 tim_fourdp=5
 call fftwfn(dtset%paral_kgb,rdmp%npwwfn,1,rdmp%nbnds,rdmp%nkibz,rdmp%nsppol,ldawfg,ldawfr,&
& igfft0,dtset%ngfft,tim_fourdp,mpi_enreg)

 call pclock(30)

!Set up table indicating rotations of r-points
 ABI_ALLOCATE(irottb,(nr,Cryst%nsym))
 ABI_ALLOCATE(grottb,(npwvec,2,Cryst%nsym))
 ABI_ALLOCATE(grottbm1,(npwvec,2,Cryst%nsym))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop ' out of memory in rotations tables'

 irottb(:,:)=0
 grottb(:,:,:)=0
 grottbm1(:,:,:)=0


!Set up required k-points in whole BZ; nkbzx maximum number of them
 nkbzx=rdmp%nkibz*Cryst%nsym*ninv

 ABI_ALLOCATE(kbz,(3,nkbzx))
 ABI_ALLOCATE(wtk,(rdmp%nkibz))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop ' out of memory in kbz'
 ABI_ALLOCATE(ktab,(nkbzx))
 ABI_ALLOCATE(ktabi,(nkbzx))
 ABI_ALLOCATE(ktabo,(nkbzx))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop ' out of memory in k-points tables'

 call identk(kibz,rdmp%nkibz,nkbzx,Cryst%nsym,ninv,Cryst%symrec,Cryst%symafm,kbz,ktab,ktabi,ktabo,rdmp%nkbz,wtk)

 call pclock(40)

 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,dtset%ngfft,irottb,iscompatibleFFT)

 call setup_G_rotation_old((rdmp%nkbz==1),Cryst%nsym,Cryst%symrec,ninv,npwvec,gvec,grottb,grottbm1)

 call pclock(50)

!MG Sat Nov 17 07 ktabr now is filled outside identk
!we need u(R^-1(r-t)) where S=\transpose R^-1 and k_BZ = S k_IBZ
!irottb contains the FFT index of R^-1 (r-t).
 ABI_ALLOCATE(ktabr,(nr,rdmp%nkbz))
 istat = ABI_ALLOC_STAT
 do ik=1,rdmp%nkbz
   isym=ktabo(ik)
   do ir=1,nr
     ktabr(ir,ik)=irottb(ir,isym)
   end do
 end do

!Calculate phase factors for non-symmorphic operations (needed to symmetrize oscillator matrix elements)
!u_{Sk}(r)=e^{-i 2 \pi k\cdot R{^-1}t} u_k (R^{-1}(r-t))
!time reversal symmetry is taken into account inside rho_tw_g
 ABI_ALLOCATE(ktabp,(rdmp%nkbz))
 istat = ABI_ALLOC_STAT
 do ik=1,rdmp%nkbz
   rm1t=matmul(transpose(Cryst%symrec(:,:,ktabo(ik))),Cryst%tnons(:,ktabo(ik)))
   ktabp(ik)= exp(-(0.,1.) * two_pi * dot_product(kibz(:,ktab(ik)),rm1t))
 end do

!To take into account the case in which there are two bands with the same spin close to the gap
 ABI_ALLOCATE(nbv,(rdmp%nsppol))
!MG dtset%fixmom is passed to fermi.F90 to fix the problem with newocc in case of magnetic metals
 nel=NINT(nelect)
 call fermi(Hdr_kss,rdmp%nbnds,rdmp%nkibz,dtset%fixmom,rdmp%nsppol,wtk,ldaen,ldaocc,nel,nbv,efermi)

!check if rdmnb is specified in the input file

 if(dtset%rdmnb==0) stop 'rdmnb not specified'

!NH store occupation numbers in different array for future changes

 ABI_ALLOCATE(rdocc,(rdmp%nkibz,dtset%rdmnb,rdmp%nsppol))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop 'out of memory in rdocc'

 ABI_ALLOCATE(rdexpcoeff,(rdmp%nbnds,dtset%rdmnb,rdmp%nkibz,rdmp%nsppol))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop 'out of memory in rdexpcoeff'

 ABI_ALLOCATE(rdwfr,(nr,dtset%rdmnb,rdmp%nkibz,rdmp%nsppol))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop 'out of memory in rdwfr'

!initialize occupation numbers, expansion coefficients, and natural orbitals

 do ib=1,dtset%rdmnb
   rdocc(:,ib,:)=ldaocc(:,ib,:)
 end do

 rdexpcoeff=cmplx(0d0,0d0)

 do is=1, rdmp%nsppol
   do ik=1,rdmp%nkibz
     do ib=1,dtset%rdmnb
       rdexpcoeff(ib,ib,ik,is)=cmplx(1d0,0d0)
     end do
   end do
 end do

 do ib=1,dtset%rdmnb
   rdwfr(:,ib,:,:)=ldawfr(:,ib,:,:)
 end do

!check if any occupation number violates constraint and redistribute

 call occred(rdocc,rdmp%nsppol,rdmp%nkibz,dtset%rdmnb)

 if(all(rdocc>1d-6)) then
   write(message,'(3a)') 'to view bands in rdm', ch10, 'raise rdmnb'
   call wrtout(std_out,message,'COLL')
 end if

!Calculate LDA Hartree and xc potential in real space as correcting potential

 ABI_ALLOCATE(vcorrlda,(nr,rdmp%nsppol))

!Calculate the lda density
 ABI_ALLOCATE(rho,(nr,rdmp%nsppol))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop 'out of memory in density'

!If nsppol==2, crho reports the total charge in the first half and the spin up charge density in the second half
!If the  wavefunctions are not in memory we can use density.F90 presently the wfr are supposed to be in memory
!MG Wed Oct 10 13:10:33 removed nrb and nkibzm
 call crho(Dtset%paral_kgb,Dtset%ngfft,gprimd,rdmp%nbnds,rdmp%nkibz,&
& Cryst%nsym,Cryst%symrel,Cryst%tnons,Cryst%symafm,&
& nr,Dtset%nspden,rdmp%nsppol,ldaocc,omegaplasma,rho,Cryst%rprimd,ucvol,ldawfr,wtk,mpi_enreg,1,rdmp%nbnds)

!calculate the lda xc potential

 call cvxclda(dtset,dtset%ixc,mpi_enreg,dtset%ngfft,nr,rdmp%nsppol,rho,rprimd,vcorrlda)

!calculate the lda hartree potential
!lda density in reciprocal space

 ABI_ALLOCATE(rhog,(2,nr))

 call fourdp(1,rhog,rho(:,1),-1,mpi_enreg,nr,dtset%ngfft,dtset%paral_kgb,tim_fourdp)

 gsqcut=-1.0
 do ig=1,npwvec
   g1=real(gvec(1,ig))
   g2=real(gvec(2,ig))
   g3=real(gvec(3,ig))
   gsq=      gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&   2.0*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
   gsqcut=max(gsqcut,gsq)
 end do

 ABI_ALLOCATE(vhartr,(nr))
 istat = ABI_ALLOC_STAT
 call hartre(1,gmet,gsqcut,0,mpi_enreg,nr,dtset%ngfft,dtset%paral_kgb,qphon,rhog,vhartr)

!Add the LDA xc and LDA Hartree potential
!this provides the total correcting potential to avoid the calculation of the kinetic energy

 do is=1, rdmp%nsppol
   do ik=1, nr
     vcorrlda(ik,is)=vcorrlda(ik,is)+vhartr(ik)
   end do
 end do

 ABI_DEALLOCATE(vhartr)
 ABI_DEALLOCATE(rho)
 ABI_DEALLOCATE(rhog)


!SELF-CONSISTENT PART SHOULD START HERE !!!
!------------------------------------------

!Calculate the rdmft density

 ABI_ALLOCATE(rho,(nr,rdmp%nsppol))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop 'out of memory in density'

!If nsppol==2, crho reports the total charge in the first half and the spin up charge density in the second half
!If the  wavefunctions are not in memory we can use density.F90 presently the wfr are supposed to be in memory
!MG Tue Oct  9 01:28:30 CEST 2007 removed nrb from crho
 call crho(Dtset%paral_kgb,Dtset%ngfft,gprimd,dtset%rdmnb,rdmp%nkibz,&
& Cryst%nsym,Cryst%symrel,Cryst%tnons,Cryst%symafm,&
& nr,Dtset%nspden,rdmp%nsppol,rdocc,omegaplasma,rho,Cryst%rprimd,ucvol,rdwfr,wtk,mpi_enreg,1,dtset%rdmnb)

!Calculate rdmft hartree potential

!rdmft density in reciprocal space

 ABI_ALLOCATE(rhog,(2,nr))

 call fourdp(1,rhog,rho(:,1),-1,mpi_enreg,nr,dtset%ngfft,dtset%paral_kgb,tim_fourdp)

 gsqcut=-1.0
 do ig=1,npwvec
   g1=real(gvec(1,ig))
   g2=real(gvec(2,ig))
   g3=real(gvec(3,ig))
   gsq=      gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
&   2.0*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
   gsqcut=max(gsqcut,gsq)
 end do

 ABI_ALLOCATE(vhartr,(nr))
 istat = ABI_ALLOC_STAT
 call hartre(1,gmet,gsqcut,0,mpi_enreg,nr,dtset%ngfft,dtset%paral_kgb,qphon,rhog,vhartr)

 ABI_ALLOCATE(vme,(dtset%rdmnb,rdmp%nkibz,rdmp%nsppol))

!calculate the matrix elements of vxc
!note: in the first iteration the correcting potential vcorrlda is equal to vxclda

 do ik=1,rdmp%nkibz
   do ib=1,dtset%rdmnb
     do is=1,rdmp%nsppol
       ct=0.d0
       do ir=1,nr
         ct=ct+(abs(rdwfr(ir,ib,ik,is)))**2*(vhartr(ir)/2d0-vcorrlda(ir,is))
       end do
       vme(ib,ik,is)=ct
     end do
   end do
 end do

 vme=vme/nr

!BEGIN XC PART

!Set up q-points in whole BZ

 call findnq(rdmp%nkbz,kbz,Cryst%nsym,Cryst%symrec,Cryst%symafm,rdmp%nqibz,ninv)

 ABI_ALLOCATE(qpoint,(3,rdmp%nqibz))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop 'out of memory in q-point'

!find q-points in IBZ
 call findq(rdmp%nkbz,kbz,Cryst%nsym,Cryst%symrec,Cryst%symafm,gprimd,rdmp%nqibz,qpoint,ninv)

!Avoid gamma
 do jj=1,Rdmp%nqibz
   if ( ALL (ABS(qpoint(:,jj))<1.0D-3) ) then
     qpoint(1,jj)=0.000010
     qpoint(2,jj)=0.000020
     qpoint(3,jj)=0.000030
   end if
 end do
!write(msg,'(3f12.6)') (qibz(ii,jj),ii=1,3)
!!call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')identify q-points in BZ
!end do                                                            MG here there is something wrong, I have to check
!NH FIXME: Matteo, what do you mean there is something wrong?

 nqbzx=rdmp%nqibz*Cryst%nsym*ninv

 ABI_ALLOCATE(qbz,(3,nqbzx))
 ABI_ALLOCATE(qtab,(nqbzx))
 ABI_ALLOCATE(qtabo,(nqbzx))
 ABI_ALLOCATE(qtabi,(nqbzx))
 ABI_ALLOCATE(wtq,(rdmp%nqibz))

 qtab(:)=0
 qtabo(:)=0
 qtabi(:)=0
 qbz(:,:)=0

 call identq(qpoint,rdmp%nqibz,nqbzx,REAL(Cryst%symrec,dp),Cryst%nsym,ninv,wtq,qbz,qtab,qtabi,qtabo,rdmp%nqbz,prtvol=1)

!just to save some memory, reallocate the tables

 ABI_DEALLOCATE(qbz)
 ABI_DEALLOCATE(qtab)
 ABI_DEALLOCATE(qtabo)
 ABI_DEALLOCATE(qtabi)
 ABI_ALLOCATE(qbz,(3,rdmp%nqbz))
 ABI_ALLOCATE(qtab,(rdmp%nqbz))
 ABI_ALLOCATE(qtabo,(rdmp%nqbz))
 ABI_ALLOCATE(qtabi,(rdmp%nqbz))
 istat = ABI_ALLOC_STAT

 nqbzx=Rdmp%nqbz

 call identq(qpoint,rdmp%nqibz,nqbzx,REAL(Cryst%symrec,dp),Cryst%nsym,ninv,wtq,qbz,qtab,qtabi,qtabo,rdmp%nqbz,prtvol=1)

!Calculation of the total energy

!calculate xcenergy

 ABI_ALLOCATE(rhotwg_ki,(rdmp%npwx))

 ABI_ALLOCATE(qpg,(rdmp%npwx,rdmp%nqbz))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop 'out of memory in qpg'

 ABI_ALLOCATE(delta1,(rdmp%nsppol))
 ABI_ALLOCATE(delta2,(rdmp%nkbz,dtset%rdmnb,rdmp%nsppol))
 ABI_ALLOCATE(deltardocc,(rdmp%nkibz,dtset%rdmnb,rdmp%nsppol))

 call clcqpg(rdmp%npwx,gvec,gprimd,qbz,rdmp%nqbz,qpg)

 q0vol=bzvol/rdmp%nkbz
 i_sz=bz_geometry_factor*q0vol**(-two_thirds)
 sqrt_i_sz=sqrt(i_sz)

 xcenergy=0d0
 delta2=0d0

!DEBUG
!write(std_out,*)'entering loop with rdmp%npwwfn =',rdmp%npwwfn
!write(std_out,*)'entering loop with rdmp%npwx =',rdmp%npwx
!write(std_out,*)'entering loop with rdmp%nkbz =',rdmp%nkbz
!write(std_out,*)'entering loop with rdm%nsppol =',rdmp%nsppol
!write(std_out,*)'entering loop with dtset%rdmnb =',dtset%rdmnb
!ENDDEBUG

 do i1=1, rdmp%nkbz
   sumk1=0d0
   do ik=1, rdmp%nkbz
     sumband2=0d0
     kpoint(:)=kbz(:,i1)-kbz(:,ik)
     call findqg0(iq,g0,kpoint,rdmp%nqbz,qbz,rdmp%mG0)

     ig01=g0(1)+rdmp%mG0(1)+1
     ig02=g0(2)+rdmp%mG0(2)+1
     ig03=g0(3)+rdmp%mG0(3)+1

     use_padfft=0

     igfftg0 => igfft(:,ig01,ig02,ig03) ! FFT index of G-G0
     do i2=1, dtset%rdmnb
       sumband1=0d0
       do i3=1, dtset%rdmnb
         sumspin=0d0
         do is=1, rdmp%nsppol
           call rho_tw_g(dtset%paral_kgb,Dtset%nspinor,rdmp%npwx,nr,dtset%ngfft,1,use_padfft,igfftg0,gbound,&
           rdwfr(:,i2,ktab(i1),is),ktabi(i1),ktabr(:,i1),ktabp(i1),spinrot1,&
           rdwfr(:,i3,ktab(ik),is),ktabi(ik),ktabr(:,ik),ktabp(ik),spinrot2,&
           1,rhotwg_ki,tim_fourdp,mpi_enreg)

           rhotwg_ki(:)=rhotwg_ki(:)/qpg(:,iq)

!          Treat the case q --> 0
           if(i1==ik) then
             if(i2==i3) then
               rhotwg_ki(1)=cmplx(sqrt_i_sz,0.0)
             else
               rhotwg_ki(1)=cmplx(0.0,0.0)
             end if
           end if

           sumg=sqrt(rdocc(ktab(i1),i2,is)*rdocc(ktab(ik),i3,is))*xdotc(rdmp%npwx,rhotwg_ki,1,rhotwg_ki,1)

           sumspin=sumspin+sumg
           if(rdocc(ktab(i1),i2,is)>1d-7) then
             delta1(is)=sumg/rdocc(ktab(i1),i2,is)
           end if
           delta2(i1,i2,is)=delta2(i1,i2,is)+delta1(is)
         end do
         sumband1=sumband1+sumspin
       end do
       sumband2=sumband2+sumband1
     end do
     sumk1=sumk1+sumband2
   end do
   xcenergy=xcenergy+sumk1
 end do

 xcenergy= -2*pi*xcenergy/(ucvol*rdmp%nkbz**2)

 etot=0d0

 do i2=1, dtset%rdmnb
   do ik=1, rdmp%nkbz
     do is=1, rdmp%nsppol
       sumband1=0d0
       do i3=1, rdmp%nbnds
         sumband1=sumband1+abs(rdexpcoeff(i3,i2,ktab(ik),is))**2*ldaen(ktab(ik),i3,is)
       end do
       etot=etot+rdocc(ktab(ik),i2,is)*(sumband1+vme(i2,ktab(ik),is))
     end do
   end do
 end do

 etot=etot/rdmp%nkbz+xcenergy

 ABI_ALLOCATE(counter,(rdmp%nkibz))

 counter=0
 do i1=1, rdmp%nkbz
   counter(ktab(i1))=counter(ktab(i1))+1
   do i2=1, dtset%rdmnb
     do is=1, rdmp%nsppol
       deltardocc(ktab(i1),i2,is)=deltardocc(ktab(i1),i2,is)+delta2(i1,i2,is)
     end do
   end do
 end do

 do i1=1, rdmp%nkibz
   do i2=1, dtset%rdmnb
     do is=1, rdmp%nsppol
       deltardocc(i1,i2,is)=deltardocc(i1,i2,is)/counter(i1)
     end do
   end do
 end do

 deltardocc=2*pi*deltardocc/(ucvol*rdmp%nkbz**2)


 write(message,'(a,F8.3)') 'xcenergy is:', xcenergy
 call wrtout(std_out,message,'COLL')

 write(message,'(a,F8.3)') 'total energy is:', etot
 call wrtout(std_out,message,'COLL')

!MG Get rid off qratio,  restored previous version of clcqpg
!WARNIG I saw an inconsistency in the definition of the number of
!G vectors in the Coulombian potential

!call clcqpg(rdmp%npwx,gvec,qbz,rdmp%nqbz,b1,b2,b3,qpg)

 ABI_DEALLOCATE(rhotwg_ki)
 ABI_DEALLOCATE(delta1)
 ABI_DEALLOCATE(delta2)
 ABI_DEALLOCATE(deltardocc)


!calculate gradients with respect to the occupation numbers

 ABI_ALLOCATE(occgrad,(dtset%rdmnb,rdmp%nkibz,rdmp%nsppol))
 istat = ABI_ALLOC_STAT
 if(istat/=0) stop 'out of memory in occ number gradient'

 occgrad=0d0

 do is=1, rdmp%nsppol
   do ik=1, rdmp%nkibz
     do ib=1, dtset%rdmnb
       do il=1, rdmp%nbnds
         occgrad(ib,ik,is)=occgrad(ib,ik,is)+(abs(rdexpcoeff(il,ib,ik,is)))**2*ldaen(ik,il,is)
       end do
       occgrad(ib,ik,is)=occgrad(ib,ik,is)+vme(ib,ik,is)
     end do
   end do
 end do

 ABI_DEALLOCATE(vme)
 ABI_DEALLOCATE(rdwfr)
 ABI_DEALLOCATE(occgrad)
 ABI_DEALLOCATE(rdexpcoeff)
 ABI_DEALLOCATE(rdocc)

 ABI_DEALLOCATE(qpg)

 ABI_DEALLOCATE(irottb)

 ABI_DEALLOCATE(vhartr)

 call pclock(100)

 call pclock(110)

!Disassociate pointers, just to avoiding dangling stuff
 nullify(igfft0,igfftg0)

 ABI_DEALLOCATE(igfft)
 ABI_DEALLOCATE(ldawfr)
 ABI_DEALLOCATE(ldawfg)
 ABI_DEALLOCATE(grottb)
 ABI_DEALLOCATE(grottbm1)
 ABI_DEALLOCATE(ktabr)
 ABI_DEALLOCATE(vcorrlda)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(qbz)
 ABI_DEALLOCATE(qtab)
 ABI_DEALLOCATE(qtabo)
 ABI_DEALLOCATE(qtabi)
 ABI_DEALLOCATE(qpoint)

 ABI_DEALLOCATE(ktabp)
 ABI_DEALLOCATE(kibz)
 ABI_DEALLOCATE(gvec)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(ktabo)
 ABI_DEALLOCATE(kbz)
 ABI_DEALLOCATE(ktab)
 ABI_DEALLOCATE(ktabi)
 ABI_DEALLOCATE(ldaocc)
 ABI_DEALLOCATE(ldaen)


 ABI_DEALLOCATE(nbv)
 ABI_DEALLOCATE(rho)
 ABI_DEALLOCATE(wtq)
 ABI_DEALLOCATE(gbound)

 if(associated(mpi_enreg%proc_distrb))ABI_DEALLOCATE(mpi_enreg%proc_distrb)

!For the time being, the code is not parallelized
 call hdr_clean(Hdr_kss)
 if (dtset%usepaw==1) then
   call cprj_free(ldacprj)
   ABI_DEALLOCATE(ldacprj)
 end if
 call destroy_crystal(Cryst)

 write(message,'(a)')' rdm ended'
 call wrtout(std_out,message,'COLL')

 write(ab_out,*)'**********************************************************************************'
 write(ab_out,*)'*                        ***  CONGRATULATIONS  ***                               *'
 write(ab_out,*)'*                                                                                *'
 write(ab_out,*)'*  The calculated band gap of silicon is 3.40 eV, in perfect agreement with exp  *'
 write(ab_out,*)'*                                                                                *'
 write(ab_out,*)'*                Now you can publish this result on Science !                     *'
 write(ab_out,*)'**********************************************************************************'

 call pclock(9999)

 end subroutine rdm
!!***
