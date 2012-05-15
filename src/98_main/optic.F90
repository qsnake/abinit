!{\src2tex{textfont=tt}}
!!****p* ABINIT/optic
!! NAME
!! optic
!!
!! FUNCTION
!! Driver routine to call linopt and nlinopt, which calculate
!! the linear and non-linear optical responses in the RPA.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (SSharma,MVer,VRecoules)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!  domega=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  eigen11(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 100
!!  eigen12(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 010
!!  eigen13(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 001
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  entropy= entropy associated with the smearing (adimensional)
!!  fermie= fermi energy (Hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gmet_inv(3,3)=inverse of reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  nomega=number of frequency for conductivity computation
!!  mband=maximum number of bands.
!!  natom = number of atoms in the unit cell.
!!  nband(nkpt*nsppol)=number of bands at each RF k point for each spin.
!!  nelect=number of electrons per unit cell
!!  nkpt=number of k points in the IBZ for this perturbation
!!  ngfft(3)=integer fft box dimensions.
!!  nspinor=number of spinorial components of the wavefunctions.
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  ntypat = number of atom types.
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k.
!!  occopt==option for occupancies
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  tsmear=smearing width (or temperature) in Hartree
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  maxomega=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!  znucl(natom)=atomic number of atoms
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,getwtk,hdr_clean,hdr_io,hdr_skip,herald,initmpi_seq
!!      int2char4,linopt,mati3inv,matr3inv,metric,nlinopt,nullify_mpi_enreg
!!      pmat2cart,pmat_renorm,sym2cart,wffclose,wffopen,wffreadeigk,xmpi_end
!!      xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program optic

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile
 use defs_abitypes
 use m_build_info

 use m_header,          only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'optic'
 use interfaces_27_toolbox_oop
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
 integer :: accesswff,bantot,bdtot0_index,bdtot_index!,dosdeltae
 integer :: fform0,formeig,formeig0,headform,ierr,ii,ikpt,isym
 integer :: isppol,mband,nomega,natom,nband1
 integer :: nsym
 integer :: master,me
 integer :: nkpt,nspinor,nsppol,ntypat
 integer :: occopt,rdwr,spaceComm,tim_rwwf
 integer :: linflag(9),mlinflag,lin1,lin2
 integer :: nlinflag(27),mnlinflag,nlin1,nlin2,nlin3
 integer,allocatable :: nband(:)
 integer,allocatable :: symrel(:,:,:)
 integer,allocatable :: symrec(:,:,:)
 real(dp),allocatable :: symcart(:,:,:)
 real(dp) :: domega,ecut,fermie!,maxocc,entropy
!real(dp) :: nelect
 real(dp) :: tsmear,ucvol,maxomega,sc,tol!,tphysel
 real(dp) :: gmet(3,3),gmet_inv(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: kpt(:,:)
 real(dp),allocatable :: cond_kg(:),cond_nd(:),doccde(:)
 real(dp),allocatable :: eig0tmp(:),eigen0(:),eigen11(:)
 real(dp),allocatable :: eigen12(:),eigtmp(:)
 real(dp),allocatable :: eigen13(:),occ(:),wtk(:)
 complex(dpc),allocatable :: pmat(:,:,:,:,:,:)
!
 character(len=fnlen) :: filnam,filnam0,filnam1,filnam2,filnam3,filnam_out
!  for the moment this is imposed by the format in linopt.f and nlinopt.f
 character(len=256) :: fn_radix,tmp_radix
 character(len=4) :: s1,s2,s3
 character(len=24) :: codename
 type(hdr_type) :: hdr
 type(wffile_type) :: wff0,wff1,wff2,wff3
 type(MPI_type) :: mpi_enreg

! *********************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

 call xmpi_init()
 call nullify_mpi_enreg(mpi_enreg)

!Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg)

 codename='OPTIC '//repeat(' ',18)
 call herald(codename,abinit_version,std_out)
!YP: calling dump_config() makes tests fail => commented
!call dump_config()

!Read data file name
 write(std_out,'(a)')' Please, give the name of the data file ...'
 read(5, '(a)')filnam
 write(std_out,'(a)')' The name of the data file is :',filnam
 write(std_out,'(a)')' Please, give the name of the output file ...'
 read(5, '(a)')filnam_out
 write(std_out,'(a)')' The name of the output file is :',filnam_out
 write(std_out,'(a)')' Please, give the root name for the (non)linear optical data output file ...'
 read(5, '(a)')fn_radix
 write(std_out,'(a)')' The root name of the output files is :',trim(fn_radix)

!Read data file
 open(15,file=filnam,form='formatted')
 rewind(15)
 read(15,'(a)')filnam1       ! first ddk file
 read(15,'(a)')filnam2       ! second ddk file
 read(15,'(a)')filnam3       ! third ddk file
 read(15,'(a)')filnam0       ! ground-state data

!Open the Wavefunction files
!These default values are typical of sequential use
 accesswff=IO_MODE_FORTRAN ; spaceComm=abinit_comm_serial ; master=0 ; me=0
 call WffOpen(accesswff,spaceComm,filnam0,ierr,wff0,master,me,10)
 call WffOpen(accesswff,spaceComm,filnam1,ierr,wff1,master,me,11)
 call WffOpen(accesswff,spaceComm,filnam2,ierr,wff2,master,me,12)
 call WffOpen(accesswff,spaceComm,filnam3,ierr,wff3,master,me,13)

!Read the header from the first ddk file (might have been the GS file ?)
 rdwr=1
 call hdr_io(fform0,hdr,rdwr,wff0)

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
 fermie=hdr%fermie
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

 nsym=hdr%nsym
 ABI_ALLOCATE(symrel,(3,3,nsym))
 ABI_ALLOCATE(symrec,(3,3,nsym))
 symrel(:,:,:) = hdr%symrel(:,:,:)
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

 ABI_ALLOCATE(kpt,(3,nkpt))
 kpt(:,:) = hdr%kptns(:,:)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))
 do ii=1,nkpt
   if (nband(ii) /= mband) then
     write(std_out,*) 'optic : Error : nband must be constant across kpts'
     stop
   end if
 end do

 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
 write(std_out,'(a,i8)')       ' natom             =',natom
 write(std_out,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband
 write(std_out,'(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'

!Prepare the reading of ddk Wff files
 formeig0=0 ; formeig=1 ; tim_rwwf=0
 ABI_ALLOCATE(eigtmp,(2*mband*mband))
 ABI_ALLOCATE(eig0tmp,(mband))
 call hdr_skip(wff1,ierr)
 call hdr_skip(wff2,ierr)
 call hdr_skip(wff3,ierr)

!Read the eigenvalues of ground-state and ddk files
 ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen11,(2*mband*mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen12,(2*mband*mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen13,(2*mband*mband*nkpt*nsppol))
 bdtot0_index=0 ; bdtot_index=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband1=nband(ikpt+(isppol-1)*nkpt)
     call WffReadEigK(eig0tmp,0,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff0)
     eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)
     call WffReadEigK(eigtmp,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff1)
     eigen11(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
     call WffReadEigK(eigtmp,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff2)
     eigen12(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
     call WffReadEigK(eigtmp,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband1,tim_rwwf,wff3)
     eigen13(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
     bdtot0_index=bdtot0_index+nband1
     bdtot_index=bdtot_index+2*nband1**2
   end do
 end do
 call WffClose(wff0,ierr)
 call WffClose(wff1,ierr)
 call WffClose(wff2,ierr)
 call WffClose(wff3,ierr)

 ABI_DEALLOCATE(eigtmp)
 ABI_DEALLOCATE(eig0tmp)

!---------------------------------------------------------------------------------
!gmet inversion
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call matr3inv(gmet,gmet_inv)


!---------------------------------------------------------------------------------
!derivative of occupation wrt the energy.
 ABI_ALLOCATE(wtk,(nkpt))
 call getwtk(kpt,nkpt,nsym,symrec,wtk)

 ABI_ALLOCATE(doccde,(mband*nkpt*nsppol))

 read(15,*)tsmear

!if (occopt==1) then
!write(std_out,'(a,i4)')  ' occopt            =',occopt
!doccde=0.0d0
!else
!tphysel=zero
!maxocc=two/(nsppol*nspinor)
!dosdeltae=zero
!call getnel(doccde,dosdeltae,eigen0,entropy,fermie,maxocc,mband,nband,&
!&   nelect,nkpt,nsppol,occ,occopt,1,tphysel,tsmear,11,wtk)
!!DEBUG
!! write(std_out,'(a,f10.5)')' getnel : nelect   =',nelect
!!ENDDEBUG
!end if

!---------------------------------------------------------------------------------
!size of the frequency range
 read(15,*)domega,maxomega
 nomega=int((maxomega+domega*0.001_dp)/domega)
 maxomega = dble(nomega)*domega
 ABI_ALLOCATE(cond_nd,(nomega))
 ABI_ALLOCATE(cond_kg,(nomega))

!Here should read in the scissor shift if wanted
 read(15,*) sc
 write(std_out,'(a,f10.5,a)')' Scissor shift     =', sc, ' Ha'
!tolerance for singularities (small)
 read(15,*) tol
 write(std_out,'(a,f10.5,a)')' Tolerance on closeness to singularities     =', tol, ' Ha'

 read(15,*) mlinflag
 read(15,*) linflag(1:mlinflag)
 write(std_out,'(a)') ' linear coeffs to be calculated : '
!write(std_out,'(a)') ' xx yy zz yz xz xy'
 write(std_out,'(9i3)') linflag(1:mlinflag)
 read(15,*) mnlinflag
 read(15,*) nlinflag(1:mnlinflag)
 write(std_out,'(a)') ' non-linear coeffs to be calculated : '
!write(std_out,'(a)') ' xxx yyy zzz xyz xxz xxy yyz yxz yxy zyz zxz zxy'
 write(std_out,'(27i4)') nlinflag(1:mnlinflag)

 close(15)

 ABI_ALLOCATE(symcart,(3,3,nsym))
 call sym2cart(gprimd,nsym,rprimd,symrel,symcart)

 ABI_ALLOCATE(pmat,(2,mband,mband,nkpt,3,nsppol))
 write(std_out,*) ' optic : Call pmat2cart'
 call pmat2cart(eigen11,eigen12,eigen13,mband,nkpt,nsppol,pmat,rprimd)

 call pmat_renorm(fermie, eigen0, mband, nkpt, nsppol, pmat, sc)

!IN CALLED ROUTINE
!call linopt(nspin,,nkpt,wkpt,nsymcrys,symcrys,nstval,occv,evalv,efermi,pmat, &
!v1,v2,nmesh,de,sc,brod)
!
!v1,v2=desired component of the dielectric function(integer) 1=x,2=y,3=z
!nmesh=desired number of energy mesh points(integer)
!de=desired step in energy(real); nmesh*de=maximum energy
!sc=scissors shift in Ha(real)
!brod=broadening in Ha(real)
!
 write(std_out,*) ' optic : Call linopt'

 do ii=1,mlinflag
   lin1 = int(linflag(ii)/10.0_dp)
   lin2 = mod(linflag(ii),10)
   write(std_out,*) ' linopt ', lin1,lin2
   call int2char4(lin1,s1)
   call int2char4(lin2,s2)
   tmp_radix = trim(fn_radix)//"_"//trim(s1)//"_"//trim(s2)
   call linopt(nsppol,ucvol,nkpt,wtk,nsym,symcart,mband,occ,eigen0,fermie,pmat, &
   lin1,lin2,nomega,domega,sc,tsmear,tmp_radix)
 end do
!IN CALLED ROUTINE
!call nlinopt(nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,evalv,efermi,pmat, &
!v1,v2,v3,emesh,de,sc,brod,tol)
 write(std_out,*) ' optic : Call nlinopt'

 do ii=1,mnlinflag
   nlin1 = int( nlinflag(ii)/100.0_dp)
   nlin2 = int((nlinflag(ii)-nlin1*100.0_dp)/10.0_dp)
   nlin3 = mod( nlinflag(ii),10)
   call int2char4(nlin1,s1)
   call int2char4(nlin2,s2)
   call int2char4(nlin3,s3)
   tmp_radix = trim(fn_radix)//"_"//trim(s1)//"_"//trim(s2)//"_"//trim(s3)
   write(std_out,*) ' nlinopt ', nlinflag(ii),nlin1,nlin2,nlin3
   call nlinopt(nsppol,ucvol,nkpt,wtk,nsym,symcart,mband,eigen0,fermie,pmat, &
&   nlin1,nlin2,nlin3,nomega,domega,sc,tsmear,tol,tmp_radix)
 end do

 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(eigen11)
 ABI_DEALLOCATE(eigen12)
 ABI_DEALLOCATE(eigen13)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(cond_nd)
 ABI_DEALLOCATE(cond_kg)

 ABI_DEALLOCATE(kpt)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(symcart)
 ABI_DEALLOCATE(pmat)

 call hdr_clean(hdr)

 call xmpi_end()

 end program optic
!!***
