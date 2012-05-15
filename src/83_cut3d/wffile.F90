!{\src2tex{textfont=tt}}
!!****f* ABINIT/wffile
!! NAME
!! wffile
!!
!! FUNCTION
!! Part of cut3d that gives the wavefunction for one kpt,one band
!! and one spin polarisation in real space.  The output depends on
!! the chosen option.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (JFB, MCote, MVer,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Needs an unformatted wave function from abinit.
!! ecut= effective ecut (ecut*dilatmx**2)
!! exchn2n3d= if 1, n2 and n3 are exchanged
!! headform= format of the wf file
!! istwfk= input variable indicating the storage option of each k-point
!! natom = number of atoms in the unit cell
!! nband= size of e_kpt
!! nkpt= number of k-points
!! npwarr= array holding npw for each k point
!! nr1,nr2,nr3 = grid size (nr1 x nr2 x nr3 = filrho dimension)
!! nspinor= number of spinorial components of the wavefunctions
!! nsppol= number of spin polarization
!! ntypat = number of atom type
!! paral_kgb= parallization option, it is set to 0 in the parent subroutine
!! rprim = orientation of the unit cell axes
!! tau = cartesian coordinates
!! typat= input variable typat(natom)
!! wff= structure type containing the density/wavefunction information
!! znucl= znucltypat(ntypat) from alchemy
!!
!! OUTPUT
!! Depends on the option chosen.
!! It is the wave function for the k point, band and spin polarisation
!! chosen.  It can be written in different ways. The option are describe
!! with the option list.  It is possible to output a Data Explorer file.
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!      clsopn,date_and_time,dens_in_sph,fourwf,getkpgnorm,getph,handle_ncerr
!!      hdr_skip,init_bess_spl,initylmg,int2char,kpgio,leave_new,metric,ph1d3d
!!      recip_ylm,rwwf,sort_dp,sphereboundary,splint,wffclose,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wffile(ecut,exchn2n3d,headform,istwfk,kpt,natom,nband,nkpt,npwarr,&
     &nr1,nr2,nr3,nspinor,nsppol,ntypat,paral_kgb,rprimd,tau,typat,wff,znucl)

 use m_profiling

 use defs_basis
!no_abirules
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif
 use defs_datatypes
 use defs_abitypes
 use m_splines
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wffile'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_27_toolbox_oop
 use interfaces_28_numeric_noabirule
 use interfaces_42_geometry
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_occeig
 use interfaces_65_nonlocal
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,headform,natom,nkpt,nr1,nr2,nr3,nspinor,nsppol
 integer,intent(in) :: ntypat,paral_kgb
 real(dp),intent(in) :: ecut
 type(wffile_type),intent(inout) :: wff
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt),npwarr(nkpt),typat(natom)
 real(dp),intent(in) :: kpt(3,nkpt),rprimd(3,3),znucl(ntypat)
 real(dp),intent(inout) :: tau(3,natom)

!Local variables-------------------------------
  character(len=*), parameter :: INPUTfile='cut.in'
!scalars
 integer,save :: tim_fourwf=0,tim_rwwf=0
 integer :: cband,cgshift,ckpt,cplex,cspinor,csppol,formeig,gridshift1
 integer :: gridshift2,gridshift3,ia,iatom,iband,ichoice,ierr,ifile
 integer :: ii1,ii2,ii3,ikpt,ilang,ioffkg,ios,iout,iprompt,ipw
 integer :: ir1,ir2,ir3,isppol,ivect,ixfh,ixint,mband,mbess,mcg,mgfft
 integer :: mkmem,mlang,mpw,mstat,n4,n5,n6,nband_disk,nfit,npw_k
 integer :: nradintmax,nxfh,oldcband,oldckpt,oldcspinor,oldcsppol,option
 integer :: prtsphere,select_exit,unkg,unout=12,unylm=0
 integer :: ikpt_qps,nkpt_qps,nband_qps,iscf_qps
 real(dp) :: arg,bessargmax,bessint_delta,kpgmax,ratsph,tmpi,tmpr,ucvol,weight
 real(dp) :: xnow,ynow,znow
 real(dp) :: eigen_qps
 character(len=1) :: outputchar
 character(len=10) :: string
 character(len=4) :: mode_paral
 character(len=500) :: message
 character(len=fnlen) :: kgnam,output,output1
 logical :: filexist
 type(MPI_type) :: mpi_enreg
!arrays
 integer :: atindx(natom),iatsph(natom),ngfft(18),nradint(natom)
 integer,allocatable :: gbound(:,:),iindex(:),kg(:,:),kg_dum(:,:),kg_k(:,:)
 integer,allocatable :: npwarr1(:),npwarrk1(:),npwtot1(:)
 real(dp) :: cmax(natom),gmet(3,3),gprimd(3,3)
 real(dp) :: phkxred(2,natom),ratsph_arr(natom),rmet(3,3),shift_tau(3)
 real(dp) :: tau2(3,natom),xred(3,natom)
 real(dp) :: kpt_qps(3)
 real(dp),allocatable :: bess_fit(:,:,:),bess_spl(:,:),bess_spl_der(:,:)
 real(dp),allocatable :: cg(:,:),cgcband(:,:),denpot(:,:,:),eigen(:)
 real(dp),allocatable :: fofgout(:,:),fofr(:,:,:,:),k1(:,:)
 real(dp),allocatable :: kpgnorm(:),occ1(:),ph1d(:,:),ph3d(:,:,:),rint(:)
 real(dp),allocatable :: sum_1atom_1ll(:,:),sum_1atom_1lm(:,:)
 real(dp),allocatable :: x_bess(:),xfhist(:,:,:,:),xfit(:),yfit(:),ylm_k(:,:)
 real(dp),allocatable :: ylmgr_dum(:,:,:)
 character(len=fnlen) :: fileqps
 character(len=fnlen),allocatable :: filename(:)
 complex(dp),allocatable :: ccoeff(:,:),wfg(:,:),wfg_qps(:)
  !no_abirules
  !For NetCDF********************************************************
#if defined HAVE_TRIO_NETCDF
  integer :: ncid, ncerr, gridsize1DimID, gridsize2DimID, gridsize3DimID
  integer :: latDimID, nbatomDimID, imagwavefunVarID,realwavefunVarID, &
&  latticevecVarID, originVarID,grid1VarID,grid2VarID,grid3VarID
  integer :: atomposiVarID, atomicnumVarID, titlechoice, posDimID,kpointVarID
  character(len=500) :: filetitle
  integer :: dd,mm,yyyy, igrid
  integer :: values(8)
  character(len=5) :: strzone
  character(len=8) :: strdat
  character(len=10) :: strtime
  character(len=3), parameter :: monnam(12)=(/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
  character(len=11) ::stridate
  real :: originatt(3,3), gridwavefun1(3,2), gridwavefun2(3,2), gridwavefun3(3,2)
  real,allocatable :: partwf(:,:,:)
  real :: kptvar(3)
#endif

! ***********************************************************************

!BEGIN EXECUTABLE SECTION

 mpi_enreg%paralbd=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%me=0
 mpi_enreg%me_fft=0
 mpi_enreg%paral_fft=0
 mpi_enreg%paral_compil_kpt=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%paral_spin=0
 mpi_enreg%fft_option_lob=1
 mpi_enreg%mode_para="n"
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0

 formeig=0
 oldckpt=0
 oldcband=0
 oldcsppol=0
 oldcspinor=0

 iout=-1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

!get xred
 call xredxcart(natom,-1,rprimd,tau,xred)

 do iatom=1,natom
   iatsph(iatom) = iatom
   atindx(iatom) = iatom
 end do

!max ang mom + 1
 mlang = 5

 ABI_ALLOCATE(kg_dum,(3,0))

 ABI_ALLOCATE(ph1d,(2,(2*nr1+1+2*nr2+1+2*nr3+1)*natom))
 call getph(atindx,natom,nr1,nr2,nr3,ph1d,xred)

 do
!  Get k-point, band and spin polarisation for the output
   if(nkpt/=1)then
     write(std_out,*)
     write(std_out,'(a,i4,a)') ' For which k-points? (1 to ',nkpt,')'
     read(5,*)ckpt
!    Check if kpt exist
     if(ckpt<1 .or. ckpt>nkpt) then
       write(std_out,*) 'Invalid k-point ',ckpt
       stop
     end if
   else
     ckpt=nkpt
   end if
   write(std_out,*) ' => Your k-point is : ',ckpt
   write(std_out,*)

   if(nband(ckpt)/=1)then
     write(std_out,*)
     write(std_out,'(a,i5,a)') ' For which band ? (1 to ',nband(ckpt),')'
     read(5,*)cband
!    Check if band number exist

     if(cband<1 .or. cband>nband(ckpt)) then
       write(std_out,*) 'Invalid band number',cband
       stop
     end if
   else
     cband=nband(ckpt)
   end if
   write(std_out,*) ' => Your band number is : ',(cband)
   write(std_out,*)

   if(nsppol/=1)then
     write(std_out,*)
     write(std_out,*) ' For which spin polarisation ?'
     read(5,*)csppol
!    Check if spin polarisation exist
     if(csppol<1 .or. csppol>nsppol) then
       write(std_out,*)'Invalid spin polarisation ',csppol
       stop
     end if
   else
     csppol=1
   end if

   write(std_out,*) ' => Your spin polarisation number is : ',(csppol)
   write(std_out,*)

   if(nspinor/=1) then
     write(std_out,*) ' nspinor = ', nspinor
     write(std_out,*)
     write(std_out,*) ' For which spinor component ?'
     read(5,*) cspinor
!    Check if spin polarisation exist
     if(cspinor<1 .or. cspinor>nspinor) then
       write(std_out,*)'Invalid spinor index ',cspinor
       stop
     end if
     write(std_out,*) ' => Your spinor component is : ',(cspinor)
     write(std_out,*)
   else
     cspinor=1
   end if

!  Reading of the data if the value of ckpt and csppol are
!  different from oldckpt and oldcsppol
!  formeig=0 gstate calculation
!  else formeig=1 for response function calculation
   if(csppol/=oldcsppol .or. ckpt/=oldckpt)then
     mband=maxval(nband)
     mpw=maxval(npwarr)
     mcg=mpw*nspinor*mband
     if (allocated (cg))   then
       ABI_DEALLOCATE(cg)
       ABI_DEALLOCATE(eigen)
       ABI_DEALLOCATE(occ1)
     end if
     ABI_ALLOCATE(cg,(2,mcg))
     ABI_ALLOCATE(eigen,((2*mband)**formeig*mband))
     ABI_ALLOCATE(occ1,(mband))
!    Rewind the file and skip the header
     call clsopn(wff)
     call hdr_skip(wff,ierr)
!    Should use xdefineOff if MPI I/O : only with wff%accesswff=1 (MB)
!    Iteration on nsppol and kpt to skip the unwanted datas (option -1)
     do isppol=1,csppol
       do ikpt=1,nkpt
         if(isppol==csppol .and. ikpt==ckpt)then
           option=1
         else
           option=-1
         end if
         call rwwf(cg,eigen,formeig,headform,0,ikpt,isppol,kg_dum,&
&         mband,mcg,mpi_enreg,nband(ikpt),nband_disk,&
&         npwarr(ikpt),nspinor,occ1,option,0,tim_rwwf,wff)

!        In case one read the last wf, the history will be output too
         if(ikpt==nkpt .and. isppol==nsppol)then
           if(cband==nband(nkpt))then
             read(unit=19,iostat=ios)nxfh
             if(ios>0)then
               write(message, '(a,a,a,a,a,a)' )ch10,&
&               ' wffile : BUG -',ch10,&
&               '  An error occurred reading the input wavefunction file,',ch10,&
&               '  history record.'
               call wrtout(std_out,message,'COLL')
               call leave_new('COLL')
             else if(ios==0)then
               write(message, '(a,a,i4,a)' )ch10,&
&               ' wffile : reading',nxfh,&
&               ' (x,f) history pairs from input wf file.'
               call wrtout(std_out,message,'COLL')
             end if
             if(nxfh>=1)then
               ABI_ALLOCATE(xfhist,(3,natom+4,2,nxfh))
               do ixfh=1,nxfh
                 read(19)xfhist(:,:,:,ixfh)
                 write(message, '(a,a,i6,a)')ch10,&
&                 ' History step number ',ixfh,&
&                 ' , atom number, xcart(1:3), fcart(1:3) ='
                 call wrtout(std_out,message,'COLL')
                 do iatom=1,natom
                   write(std_out,'(i5,6es16.6)')iatom,xfhist(:,iatom,:,ixfh)
                 end do
               end do
               ABI_DEALLOCATE(xfhist)
             end if ! nxfh>1
           end if ! cband==nband(nkpt)
         end if ! ikpt==nkpt .and. isppol==nsppol

         if(option==1)exit  ! When the target wf has been read, exit the wf file reading
       end do
       if(option==1)exit
     end do

   end if

   if (csppol/=oldcsppol .or. ckpt/=oldckpt .or. &
&   cband/=oldcband .or. cspinor/=oldcspinor ) then
!    The data of ckpt,cnsspol are in cg
!    Now we have to do the Fourier Transform of the datas

     ngfft(1)=nr1
     ngfft(2)=nr2
     ngfft(3)=nr3
!    ngfft(4) and ngfft(5) can not be even (see getng.f)
     if (mod(nr1,2)==0)then
       ngfft(4)=nr1+1
     else
       ngfft(4)=nr1
     end if
     if (mod(nr2,2)==0)then
       ngfft(5)=nr2+1
     else
       ngfft(5)=nr2
     end if
     ngfft(6)=nr3
!    XG 020829 : 112 does not work yet for all istwfk values
     ngfft(7)=111
     ngfft(8)=256
     ngfft(9)=0
     ngfft(10)=1
     ngfft(11)=0
     ngfft(12)=ngfft(2)
     ngfft(13)=ngfft(3)
     ngfft(14)=0

!    if iout<0, the output of metric will not be print
     mode_paral='PERS'
     mkmem=nkpt
     mgfft=maxval(ngfft(1:3))
     ABI_ALLOCATE(npwarr1,(nkpt))
     ABI_ALLOCATE(kg,(3,mpw*mkmem))
     ABI_ALLOCATE(npwtot1,(nkpt))

!    Create positions index for pw
     call kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kgnam,kpt,mkmem,nband,nkpt,&
&     mode_paral,mpi_enreg,mpw,npwarr1,npwtot1,nsppol,unkg)

     ioffkg=0
     do ikpt=1,ckpt-1
       ioffkg=ioffkg+npwarr1(ikpt)
     end do
     npw_k=npwarr(ckpt)
     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,1:npw_k)=kg(:,1+ioffkg:npw_k+ioffkg)

     ABI_ALLOCATE(ylm_k,(mpw,mlang*mlang))
     mstat = ABI_ALLOC_STAT
     ABI_ALLOCATE(ylmgr_dum,(mpw,3,mlang*mlang))
     mstat = ABI_ALLOC_STAT

!    call for only the kpoint we are interested in !
     ABI_ALLOCATE(k1,(3,1))
     k1(:,1)=kpt(:,ckpt)
     ABI_ALLOCATE(npwarrk1,(1))
     npwarrk1 = (/npw_k/)
     call initylmg(gprimd,kg_k,k1,1,mpi_enreg,&
&     mlang,mpw,nband,1,npwarrk1,nsppol,0,rprimd,unkg,unylm,ylm_k,ylmgr_dum)
     ABI_DEALLOCATE(k1)
     ABI_DEALLOCATE(npwarrk1)

!    Compute the norms of the k+G vectors
     ABI_ALLOCATE(kpgnorm,(npw_k))
     call getkpgnorm(gprimd,kpt(:,ckpt),kg_k,kpgnorm,npw_k)

     call sphereboundary(gbound,istwfk(ckpt),kg_k,mgfft,npw_k)
!    Do the Fourier Transform
     n4=ngfft(4)
     n5=ngfft(5)
     n6=ngfft(6)
!    cplex=0
     cplex=1
!    Complex can be set to 0 with this option(0) of fourwf

!    Read the QPS file if GW wavefunctions are to be analysed
     write(std_out,*) 'Do you want to analyze a GW wavefunction? (1=yes,0=no)'
     read(5,*) ii1
     write(std_out,*) '=> Your choice is :',ii1
     write(std_out,*)

     if(ii1==1) then
       write(std_out,*) 'What is the name of the QPS file?'
       read(5,*) fileqps
!      Checking the existence of data file
       inquire (file=fileqps,exist=filexist)
       if (.NOT. filexist) then
         write(std_out,*) 'Error, missing data file: ',fileqps
         stop
       end if
       open(111,file=trim(fileqps),status='old')
       read(111,*) iscf_qps
       read(111,*) nkpt_qps
       read(111,*) nband_qps
       read(111,*) ikpt_qps
       
       ABI_ALLOCATE(ccoeff,(nband_qps,nband_qps))
       do ikpt=1,ckpt ! nkpt_qps
         read(111,*) kpt_qps(:)
         do iband=1,nband_qps
           read(111,*) eigen_qps
           read(111,*) ccoeff(:,iband)
         end do
       end do
       close(111)

       ABI_ALLOCATE(wfg,(npw_k,nband_qps))
       ABI_ALLOCATE(wfg_qps,(npw_k))
       do iband=1,nband_qps
         cgshift=(iband-1)*npw_k*nspinor + (cspinor-1)*npw_k
         wfg(:,iband) = cmplx( cg(1,cgshift+1:cgshift+npw_k),cg(2,cgshift+1:cgshift+npw_k) )
       end do

       wfg_qps = matmul( wfg(:,:) , ccoeff(:,cband) )

!      write(std_out,*) 'norm',SUM( abs(wfg(:,cband))**2 )
!      write(std_out,*) 'norm',SUM( abs(wfg_qps(:))**2 )
       ABI_DEALLOCATE(ccoeff)
       ABI_DEALLOCATE(wfg)
       ABI_ALLOCATE(cgcband,(2,npw_k))
       cgcband(1,:)= real(wfg_qps(:))
       cgcband(2,:)= aimag(wfg_qps(:))
       ABI_DEALLOCATE(wfg_qps)

     else

!      The shift is to get the good band values
       cgshift=(cband-1)*npw_k*nspinor + (cspinor-1)*npw_k
       ABI_ALLOCATE(cgcband,(2,npw_k))
       cgcband(:,1:npw_k)=cg(:,cgshift+1:cgshift+npw_k)

     end if

!    Fix the phase of cgcband, for portability reasons
!    call fxphas(cgcband,cgcband,0,npw_k,1,npw_k,0)

     ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
     ABI_ALLOCATE(fofgout,(2,npw_k))
     ABI_ALLOCATE(fofr,(2,n4,n5,n6))
     call fourwf(cplex,denpot,cgcband,fofgout,fofr,&
&     gbound,gbound,&
&     istwfk(ckpt),kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,&
&     npw_k,n4,n5,n6,0,paral_kgb,tim_fourwf,weight,weight)

!    Analyse wavefunction inside atomic sphere

     write(std_out,'(a)' ) ' Do you want the atomic analysis for this state : '
     write(std_out,'(a,2i5,a)' ) ' (kpt,band)= (',ckpt,cband,')? '
     write(std_out,'(a)' ) ' If yes, enter the radius of the atomic spheres, in bohr '
     write(std_out,'(a)' ) ' If no, enter 0 '
     read (*,*) ratsph
     write(std_out,'(a,f16.8,a)' ) ' You entered ratsph=',ratsph,' Bohr '

     if (ratsph >= tol10) then

       write(std_out,'(3a)' ) ch10,' Atomic sphere analysis ',ch10

!      Init bessel function integral for recip_ylm
!      max ang mom + 1
       mlang = 5
       bessint_delta = 0.1_dp
       kpgmax = sqrt(ecut)
       bessargmax = ratsph*two_pi*kpgmax
       mbess = int (bessargmax / bessint_delta) + 1
       bessargmax = bessint_delta*mbess

!      Intervals in radial integration
       nradintmax = mbess
       nradint(1:natom)=nradintmax

       write(std_out,'(a,2es16.6,i6)' ) &
&       ' wffile : kpgmax, bessargmax, nradint = ', kpgmax, bessargmax,nradintmax

!      Initialize general Bessel function array on uniform grid
!      x_bess, from 0 to (2 \pi |k+G|_{max} |r_{max}|)
       ABI_ALLOCATE(bess_spl,(mbess,mlang))
       ABI_ALLOCATE(bess_spl_der,(mbess,mlang))
       ABI_ALLOCATE(x_bess,(nradintmax))
       ABI_ALLOCATE(rint,(nradintmax))

!      call init_bess_spl(mbess,bessargmax,bessint_delta,mlang,
       call init_bess_spl(mbess,bessint_delta,mlang,&
&       bess_spl,bess_spl_der,x_bess)
!      DEBUG
!      write(std_out,*) 'wffile : after init_bess_spl :'
!      write(std_out,'(6F12.5)') bess_spl
!      ENDDEBUG

       ABI_ALLOCATE(bess_fit,(mpw,nradintmax,mlang))
       ABI_ALLOCATE(xfit,(npw_k))
       ABI_ALLOCATE(yfit,(npw_k))
       ABI_ALLOCATE(iindex,(npw_k))
       nfit = npw_k

       do ixint=1,nradintmax
         rint(ixint) = (ixint-1)*ratsph / (nradintmax-1)

         do ipw=1,npw_k
           xfit(ipw) = two_pi * kpgnorm(ipw) * rint(ixint)
           iindex(ipw) = ipw
         end do
         call sort_dp (npw_k,xfit,iindex,tol14)
         do ilang=1,mlang
           call splint(mbess,x_bess,bess_spl(:,ilang),bess_spl_der(:,ilang),&
&           nfit,xfit,yfit)
!          Re-order results for different G vectors
           do ipw=1,npw_k
             bess_fit(iindex(ipw),ixint,ilang) = yfit(ipw)
           end do
         end do ! ipw
       end do ! ixint

!      Construct phases ph3d for all G vectors in present sphere
!      make phkred for all atoms

       do ia=1,natom
         iatom=atindx(ia)
         arg=two_pi*( kpt(1,ckpt)*xred(1,ia) &
&         + kpt(2,ckpt)*xred(2,ia) &
&         + kpt(3,ckpt)*xred(3,ia))
         phkxred(1,iatom)=cos(arg)
         phkxred(2,iatom)=sin(arg)
       end do

       ABI_ALLOCATE(ph3d,(2,npw_k,natom))
!      Get full phases for the following
!      write(std_out,*) 'nr1nr2nr3 ',nr1,nr2,nr3
       call ph1d3d(1,natom,kg_k,natom,natom,npw_k,nr1,nr2,nr3,&
&       phkxred,ph1d,ph3d)
!      phases exp (2 pi i (k+G).x_tau) are now in ph3d

       ABI_ALLOCATE(sum_1atom_1ll,(mlang,natom))
       ABI_ALLOCATE(sum_1atom_1lm,(mlang**2,natom))
       prtsphere=1
       ratsph_arr(:)=ratsph
       call recip_ylm (bess_fit,cgcband,iatsph,&
&       istwfk(ckpt),&
&       nradint,nradintmax,mlang,mpi_enreg,mpw,natom,natom,npw_k,&
&       ntypat,ph3d,prtsphere,rint,&
&       ratsph_arr,sum_1atom_1ll,sum_1atom_1lm,typat,ucvol,ylm_k,znucl)

       call dens_in_sph(cmax,cgcband,gmet,istwfk(ckpt),&
&       kg_k,natom,ngfft,mpi_enreg,npw_k,paral_kgb,ph1d,ratsph_arr,ucvol)

       write(std_out,'(a)' )' Charge in the sphere around each atom '
       do iatom=1,natom
         write(std_out,'(a,i4,a,f14.8)' ) ' Atom number ',iatom,' :  charge =',cmax(iatom)
       end do

       ABI_DEALLOCATE(sum_1atom_1ll)
       ABI_DEALLOCATE(sum_1atom_1lm)
       ABI_DEALLOCATE(ph3d)
       ABI_DEALLOCATE(iindex)
       ABI_DEALLOCATE(yfit)
       ABI_DEALLOCATE(xfit)
       ABI_DEALLOCATE(bess_fit)
       ABI_DEALLOCATE(bess_spl)
       ABI_DEALLOCATE(bess_spl_der)
       ABI_DEALLOCATE(x_bess)
       ABI_DEALLOCATE(rint)
     end if ! ratsph < 0     = end if for atomic sphere analysis

     ABI_DEALLOCATE(cgcband)
     ABI_DEALLOCATE(fofgout)
     ABI_DEALLOCATE(denpot)
     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(npwarr1)
     ABI_DEALLOCATE(kg)
     ABI_DEALLOCATE(npwtot1)
     ABI_DEALLOCATE(kpgnorm)
     ABI_DEALLOCATE(ylm_k)

   end if


   write(std_out,*)
   write(std_out,*) ' 3D wave function was read. ',&
&   'Ready for further treatment.'
   write(std_out,*)
   write(std_out,*) '============================',&
&   '==============================='
   write(std_out,*)

!  ------------------------------------------------------------------------

!  At this moment all the input is done
!  The code knows the geometry of the system,
!  and the data file.


   select_exit = 0
   do while (select_exit == 0)
     write(std_out,*) ' What is your choice ? Type:'
     write(std_out,*) '  0 => exit to k-point / band / spin-pol loop'
     write(std_out,*) '  1 => 3D formatted real and imaginary data'
     write(std_out,*) '       (output the bare 3D data - two column,R,I)'
     write(std_out,*) '  2 => 3D formatted real data'
     write(std_out,*) '       (output the bare 3D data - one column)'
     write(std_out,*) '  3 => 3D formatted imaginary data'
     write(std_out,*) '       (output the bare 3D data - one column)'
     write(std_out,*) '  4 => 3D indexed real and imaginary data'
     write(std_out,*) '       (3D data, preceeded by 3D index)'
     write(std_out,*) '  5 => 3D indexed real data'
     write(std_out,*) '       (bare 3D data, preceeded by 3D index)'
     write(std_out,*) '  6 => 3D indexed imaginary data'
     write(std_out,*) '       (bare 3D data, preceeded by 3D index)'
     write(std_out,*) '  7 => 3D Data Explorer formatted data '
     write(std_out,*) '       (Real file and Imaginary file)'
     write(std_out,*) '  8 => 3D Data Explorer formatted data '
     write(std_out,*) '       (Only the Real file)'
     write(std_out,*) '  9 => 3D Data Explorer formatted data '
     write(std_out,*) '       (Only the Imaginary file)'
     write(std_out,*) ' 10 => 3D Data Explorer formatted data and position files'
     write(std_out,*) ' 11 => XCrysden formatted data (norm of wf) and position files'
     write(std_out,*) ' 12 => NetCDF data and position file'
     write(std_out,*) ' 13 => XCrysden/VENUS wavefunction (real part of data)'
     write(std_out,*) ' 14 => Gaussian/cube wavefunction module'
     read(*,*) ichoice
     write(std_out,'(a,a,i2,a)' ) ch10,' Your choice is ',ichoice,char(10)

     if (ichoice>0 .and. ichoice<15)then
       write(std_out,*) ch10,'  Enter the root of an output file:'
       read(*,*) output1
       write(std_out,*) '  The root of your file is : ',trim(output1)
       output=trim(output1)//'_k'
       call int2char(ckpt,string)
       output=trim(output)//trim(string)//'_b'
       call int2char(cband,string)
       output=trim(output)//trim(string)//'_s'
       call int2char(csppol,string)
       output=trim(output)//trim(string)
       write(std_out,*) '  The corresponding filename is : ',trim(output)
     end if

     select case(ichoice)

       case(1)            ! data R,I
         write(std_out,*)
         write(std_out,*) 'Give 1 file of 3D formatted real and imaginary data'
         write(std_out,*) 'The first column is the real data'
         write(std_out,*) 'The second column is the imaginary data'
         write(std_out,*)
         open(unit=unout,file=output,status='replace',form='formatted')
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               write(unout,'(2f20.16)') fofr(:,ir1,ir2,ir3)
             end do
           end do
         end do
         close(unout)
         exit

       case(2)            ! data R
         write(std_out,*)
         write(std_out,*) 'Give 1 file of 3D formatted real data'
         write(std_out,*) 'The only column is the real data'
         write(std_out,*)
         open(unit=unout,file=output,status='replace',form='formatted')
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               write(unout,'(f20.16)') fofr(1,ir1,ir2,ir3)
             end do
           end do
         end do
         close(unout)
         exit

       case(3)            ! data I
         write(std_out,*)
         write(std_out,*) 'Give 1 file of 3D formatted real data'
         write(std_out,*) 'The only column is the imaginary data'
         write(std_out,*)
         open(unit=unout,file=output,status='replace',form='formatted')
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               write(unout,'(f20.16)') fofr(2,ir1,ir2,ir3)
             end do
           end do
         end do
         close(unout)
         exit

       case(4)            ! coord(x,y,z) data R,I
         write(std_out,*)
         write(std_out,*) 'Give 1 file of 3D formatted data'
         write(std_out,*) 'The first three columns are the x,y,z positions(Angstrom)'
         write(std_out,*) 'The fourth column is the real data'
         write(std_out,*) 'The fifth column is the imaginary data'
         write(std_out,*)
         open(unit=unout,file=output,status='replace',form='formatted')
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               xnow = rprimd(1,1)*(ir1-1)/nr1 + rprimd(1,2)*(ir2-1)/nr2 + rprimd(1,3)*(ir3-1)/nr3
               ynow = rprimd(2,1)*(ir1-1)/nr1 + rprimd(2,2)*(ir2-1)/nr2 + rprimd(2,3)*(ir3-1)/nr3
               znow = rprimd(3,1)*(ir1-1)/nr1 + rprimd(3,2)*(ir2-1)/nr2 + rprimd(3,3)*(ir3-1)/nr3
               write(unout,'(3f16.10,2f20.16)') Bohr_Ang*xnow, Bohr_Ang*ynow, Bohr_Ang*znow,fofr(:,ir1,ir2,ir3)
             end do
           end do
         end do
         close(unout)
         exit

       case(5)            ! coord(x,y,z) data R
         write(std_out,*)
         write(std_out,*) 'Give 1 file of 3D formatted data'
         write(std_out,*) 'The first three columns are the x,y,z positions(Angstrom)'
         write(std_out,*) 'The fourth column is the real data'
         write(std_out,*)
         open(unit=unout,file=output,status='replace',form='formatted')
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               xnow = rprimd(1,1)*(ir1-1)/nr1 + rprimd(1,2)*(ir2-1)/nr2 + rprimd(1,3)*(ir3-1)/nr3
               ynow = rprimd(2,1)*(ir1-1)/nr1 + rprimd(2,2)*(ir2-1)/nr2 + rprimd(2,3)*(ir3-1)/nr3
               znow = rprimd(3,1)*(ir1-1)/nr1 + rprimd(3,2)*(ir2-1)/nr2 + rprimd(3,3)*(ir3-1)/nr3
               write(unout,'(3f16.10,2f20.16)') Bohr_Ang*xnow, Bohr_Ang*ynow, Bohr_Ang*znow,fofr(1,ir1,ir2,ir3)
             end do
           end do
         end do
         close(unout)
         exit

       case(6)            ! coord(x,y,z) data I
         write(std_out,*)
         write(std_out,*) 'Give 1 file of 3D formatted data'
         write(std_out,*) 'The first three columns are the x,y,z positions(Angstrom)'
         write(std_out,*) 'The fourth column is the imaginary data'
         write(std_out,*)
         open(unit=unout,file=output,status='replace',form='formatted')
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               xnow = rprimd(1,1)*(ir1-1)/nr1 + rprimd(1,2)*(ir2-1)/nr2 + rprimd(1,3)*(ir3-1)/nr3
               ynow = rprimd(2,1)*(ir1-1)/nr1 + rprimd(2,2)*(ir2-1)/nr2 + rprimd(2,3)*(ir3-1)/nr3
               znow = rprimd(3,1)*(ir1-1)/nr1 + rprimd(3,2)*(ir2-1)/nr2 + rprimd(3,3)*(ir3-1)/nr3
               write(unout,'(3f16.10,2f20.16)') Bohr_Ang*xnow, Bohr_Ang*ynow, Bohr_Ang*znow,fofr(2,ir1,ir2,ir3)
             end do
           end do
         end do
         close(unout)
         exit

       case(7)            !OpenDX format, data R and data I
         write(std_out,*)
         write(std_out,*) 'Give 2 files of 3D formatted data'
         write(std_out,*) 'The file is ready to be use with OpenDX'
         write(std_out,*) 'The eigenvalues and occupations numbers are in comments'
         write(std_out,*)
         ABI_ALLOCATE(filename,(2))
         filename(1)=trim(output)//'Real.dx'
         filename(2)=trim(output)//'Imag.dx'
         write(std_out,*) '  The name of your files is : '
         write(std_out,*) trim(filename(1)),'  for the real part,'
         write(std_out,*) trim(filename(2)),'  for the imaginary part.'
         write(std_out,*)

         do ifile=1,2
           open(unit=unout,file=filename(ifile),status='replace',form='formatted')
           rewind(unout)
           write(unout,*)'# band,  eigenvalues   and   occupations'
           do iband=1,nband(ckpt)
             write(unout,'(a,i4,2f20.16)')'#',iband,eigen(iband),occ1(iband)
           end do
           write(unout,'(a,i10,a)')'object "donnees" class array type float rank 0 items',&
&           nr1*nr2*nr3,' data follows'
           do ir3=1,nr3
             do ir2=1,nr2
               do ir1=1,nr1
                 write(unout,'(f20.16)')fofr(ifile,ir1,ir2,ir3)
               end do
             end do
           end do

           write(unout,'(a)')'# this is the object defining the grid connections'
           write(unout,'(a,3i5)')'object "gridconnections" class gridconnections counts',&
&           nr3,nr2,nr1
           write(unout,*)
           write(unout,*)
           write(unout,'(a)')'# this is the object defining the grid'
           write(unout,'(a,3i5)')'object "positions" class gridpositions counts',&
&           nr3,nr2,nr1

           write(unout,'(a)') 'origin 0 0 0'
           write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,3)/nr3, ii1=1,3)
           write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,2)/nr2, ii1=1,3)
           write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,1)/nr1, ii1=1,3)

           write(unout,'(a)')'# this is the collective object, one for each grid '
           write(unout,'(a)')'object "densite" class field '
           write(unout,'(a)')'component "positions"   value "positions"'
           write(unout,'(a)')'component "connections" value "gridconnections" '
           write(unout,'(a)')'component "data"        value "donnees"'

           close(unit=unout)
         end do
         ABI_DEALLOCATE(filename)
         exit

       case(8)            !OpenDX format, data R and data I
         write(std_out,*)
         write(std_out,*) 'Give 2 files of 3D formatted data'
         write(std_out,*) 'The file is ready to be use with OpenDX'
         write(std_out,*) 'The eigenvalues and occupations numbers are in comments'
         write(std_out,*)
         ABI_ALLOCATE(filename,(1))
         filename(1)=trim(output)//'Real.dx'
         write(std_out,*) '  The name of your file is : '
         write(std_out,*) trim(filename(1)),'  for the real part,'
         write(std_out,*)


         open(unit=unout,file=filename(1),status='replace',form='formatted')
         rewind(unout)
         write(unout,*)'# band,  eigenvalues   and   occupations'
         do iband=1,nband(ckpt)
           write(unout,'(a,i4,2f20.16)')'#',iband,eigen(iband),occ1(iband)
         end do
         write(unout,'(a,i10,a)')'object "donnees" class array type float rank 0 items',&
&         nr1*nr2*nr3,' data follows'
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               write(unout,'(f20.16)')fofr(1,ir1,ir2,ir3)
             end do
           end do
         end do

         write(unout,'(a)')'# this is the object defining the grid connections'
         write(unout,'(a,3i5)')'object "gridconnections" class gridconnections counts',&
&         nr3,nr2,nr1
         write(unout,*)
         write(unout,*)
         write(unout,'(a)')'# this is the object defining the grid'
         write(unout,'(a,3i5)')'object "positions" class gridpositions counts',&
&         nr3,nr2,nr1

         write(unout,'(a)') 'origin 0 0 0'
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,3)/nr3, ii1=1,3)
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,2)/nr2, ii1=1,3)
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,1)/nr1, ii1=1,3)

         write(unout,'(a)')'# this is the collective object, one for each grid '
         write(unout,'(a)')'object "densite" class field '
         write(unout,'(a)')'component "positions"   value "positions"'
         write(unout,'(a)')'component "connections" value "gridconnections" '
         write(unout,'(a)')'component "data"        value "donnees"'

         close(unit=unout)
         ABI_DEALLOCATE(filename)
         exit

       case(9)            !OpenDX format, data R and data I
         write(std_out,*)
         write(std_out,*) 'Give 2 files of 3D formatted data'
         write(std_out,*) 'The file is ready to be use with OpenDX'
         write(std_out,*) 'The eigenvalues and occupations numbers are in comments'
         write(std_out,*)
         ABI_ALLOCATE(filename,(1))
         filename(1)=trim(output)//'Imag.dx'
         write(std_out,*) '  The name of your file is : '
         write(std_out,*) trim(filename(1)),'  for the imaginary part.'
         write(std_out,*)


         open(unit=unout,file=filename(1),status='replace',form='formatted')
         rewind(unout)
         write(unout,*)'# band,  eigenvalues   and   occupations'
         do iband=1,nband(ckpt)
           write(unout,'(a,i4,2f20.16)')'#',iband,eigen(iband),occ1(iband)
         end do
         write(unout,'(a,i10,a)')'object "donnees" class array type float rank 0 items',&
&         nr1*nr2*nr3,' data follows'
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               write(unout,'(f20.16)')fofr(2,ir1,ir2,ir3)
             end do
           end do
         end do

         write(unout,'(a)')'# this is the object defining the grid connections'
         write(unout,'(a,3i5)')'object "gridconnections" class gridconnections counts',&
&         nr3,nr2,nr1
         write(unout,*)
         write(unout,*)
         write(unout,'(a)')'# this is the object defining the grid'
         write(unout,'(a,3i5)')'object "positions" class gridpositions counts',&
&         nr3,nr2,nr1

         write(unout,'(a)') 'origin 0 0 0'
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,3)/nr3, ii1=1,3)
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,2)/nr2, ii1=1,3)
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,1)/nr1, ii1=1,3)

         write(unout,'(a)')'# this is the collective object, one for each grid '
         write(unout,'(a)')'object "densite" class field '
         write(unout,'(a)')'component "positions"   value "positions"'
         write(unout,'(a)')'component "connections" value "gridconnections" '
         write(unout,'(a)')'component "data"        value "donnees"'

         close(unit=unout)
         ABI_DEALLOCATE(filename)
         exit

       case(10)           !OpenDX format, data R and data I, atoms positions, lattice and cell
         write(std_out,*)
         write(std_out,*) 'Give 5 files of formatted data'
         write(std_out,*) 'The files are ready to be use with Data Explorer'
         write(std_out,*) 'The eigenvalues and occupations numbers are in comments'
         write(std_out,*) 'of the two data files'
         write(std_out,*)
         ABI_ALLOCATE(filename,(2))
         filename(1)=trim(output)//'Real.dx'
         filename(2)=trim(output)//'Imag.dx'
         write(std_out,*) '  The name of your data files is : '
         write(std_out,*) trim(filename(1)),'  for the real part,'
         write(std_out,*) trim(filename(2)),'  for the imaginary part.'
         write(std_out,*)

         do ifile=1,2
           open(unit=unout,file=filename(ifile),status='replace',form='formatted')
           rewind(unout)
           do iband=1,nband(ckpt)
             write(unout,'(a,2f20.16)')'#', eigen(iband),occ1(iband)
           end do
           write(unout,'(a,i10,a)')'object "donnees" class array type float rank 0 items',&
&           nr1*nr2*nr3,' data follows'
           do ir3=1,nr3
             do ir2=1,nr2
               do ir1=1,nr1
                 write(unout,'(f20.16)')fofr(ifile,ir1,ir2,ir3)
               end do
             end do
           end do

           write(unout,'(a)')'# this is the object defining the grid connections'
           write(unout,'(a,3i5)')'object "gridconnections" class gridconnections counts',&
&           nr3,nr2,nr1
           write(unout,*)
           write(unout,*)
           write(unout,'(a)')'# this is the object defining the grid'
           write(unout,'(a,3i5)')'object "positions" class gridpositions counts',&
&           nr3,nr2,nr1

           write(unout,'(a)') 'origin 0 0 0'
           write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,3)/nr3, ii1=1,3)
           write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,2)/nr2, ii1=1,3)
           write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,1)/nr1, ii1=1,3)

           write(unout,'(a)')'# this is the collective object, one for each grid '
           write(unout,'(a)')'object "densite" class field '
           write(unout,'(a)')'component "positions"   value "positions"'
           write(unout,'(a)')'component "connections" value "gridconnections" '
           write(unout,'(a)')'component "data"        value "donnees"'

           close(unit=unout)
         end do
         ABI_DEALLOCATE(filename)
!        
!        write LATTICE_VEC.dx file
!        
         ABI_ALLOCATE(filename,(3))
         filename(1)=trim(output1)//'_LATTICE_VEC.dx'
         filename(2)=trim(output1)//'_ATOM_POS.dx'
         filename(3)=trim(output1)//'_UCELL_FRAME.dx'
         write(std_out,*)
         write(std_out,*)'Give the lattice file, ', trim(filename(1))
         open(unit=unout,file=filename(1),status='replace')
         write(unout,'("#",/,"#",/,"#    LATTICE VECTOR INFO:",/,"#",/,"#")')
         write(unout,'(a)') 'object "lattices" class array type float rank 1 shape 3 items 3 data follows'
         do ivect=1,3
           write(unout,'(3f16.10)')  Bohr_Ang*rprimd(1,ivect),Bohr_Ang*rprimd(2,ivect),Bohr_Ang*rprimd(3,ivect)
         end do
         write(unout,'(a,a)') 'object "lattices_location" class array type float ',&
&         'rank 1 shape 3 items 3 data follows'
         do ivect=1,3
           write(unout,'(3f16.10)')  0_dp,0_dp,0_dp
         end do
         write(unout,'("object   3 class field")')
         write(unout,'(a)') 'component "data" value "lattices"'
         write(unout,'(a)') 'component "positions" value "lattices_location"'
         close(unout)


!        
!        write ATOM_POS.dx file
!        
         write(std_out,*)'Give the atoms positions file, ', trim(filename(2))
         open(unit=unout,file=filename(2),status='unknown')
         write(unout,'("#",/,"#",/,"#    BALL AND STICK INFO:",/,"#",/,"#")')
         write(unout,'(a,i5,a)') 'object "atomcoord" array type float rank 1 shape 3 items ',natom,' data follows'
         do iatom=1,natom
           write(unout,'(3f16.10)')  Bohr_Ang*tau(1:3,iatom)
         end do
!        write(unout,'(a,i5,a)') 'object "data" array type string rank 0 shape 2 items ',&
!        &natom,' data follows'
         write(unout,'(a,i5,a)') 'object "colorcode" array type float rank 0 items ',natom,' data follows'
         do iatom=1,natom
           write(unout,'(f10.4)') znucl(typat(iatom))
         end do
         write(unout,'(a)') 'object "molecule" field'
         write(unout,'(a)') 'component "positions" value "atomcoord"'
         write(unout,'(a)') 'component "data" value "colorcode"'
         close(unout)

!        
!        write UCELL_FRAME.dx file
!        
         write(std_out,*)'Give the enveloppe of the cell file, ',trim(filename(3))
         open(unit=unout,file=filename(3),status='unknown')
         write(unout,'("#",/,"#",/,"#    UNIT CELL FRAME INFO:",/,"#",/,"#")')
         write(unout,'(a)')'object 3 class array type int rank 1 shape 2 items 12 data follows'
         write(unout,'(" 0  1",/," 0  2",/," 0  3",/," 1  4",/," 1  5",/," 3  5")')
         write(unout,'(" 3  6",/," 2  6",/," 2  4",/," 7  5",/," 7  6",/," 7  4")')
         write(unout,'(a)') 'attribute "element type" string "lines"'
         write(unout,'("object  4 class array type float rank 1 shape 3 items    8 data follows")')
         write(unout,'("      .00000000      .00000000      .00000000")')
         write(unout,'(3f20.10)') Bohr_Ang*rprimd(:,1)
         write(unout,'(3f20.10)') Bohr_Ang*rprimd(:,2)
         write(unout,'(3f20.10)') Bohr_Ang*rprimd(:,3)
         write(unout,'(3f20.10)') Bohr_Ang*(rprimd(:,1)+rprimd(:,2))
         write(unout,'(3f20.10)') Bohr_Ang*(rprimd(:,1)+rprimd(:,3))
         write(unout,'(3f20.10)') Bohr_Ang*(rprimd(:,2)+rprimd(:,3))
         write(unout,'(3f20.10)') Bohr_Ang*(rprimd(:,1)+rprimd(:,2)+rprimd(:,3))
         write(unout,'("object 5 array type float rank 0 items 12 data follows")')
         do ivect=1,12
           write(unout,'("1.0")')
         end do
         write(unout,'(a)') 'attribute "dep" string "connections"'
         write(unout,'("object 6 class field")')
         write(unout,'(a)') 'component "data" value 5'
         write(unout,'(a)') 'component "positions" value 4'
         write(unout,'(a)') 'component "connections" value 3'
         close(unout)
         ABI_DEALLOCATE(filename)

         write(std_out,*)
         exit

       case(11)
         write(std_out,*)
         write(std_out,*) 'Give 1 files of formatted data'
         write(std_out,*) 'The files are ready to be used with XCrysDen'
         write(std_out,*)
         gridshift1 = 0
         gridshift2 = 0
         gridshift3 = 0
         write(std_out,*) 'Do you want to shift the grid along the x,y or z axis (y/n)?'
         write(std_out,*)
         shift_tau(:) = 0.0
         read (*,*) outputchar
         if (outputchar == 'y' .or. outputchar == 'Y') then
           write(std_out,*) 'Give the three shifts (x,y,z < ',nr1,nr2,nr3,') :'
           write(std_out,*)
           read (*,*) gridshift1, gridshift2, gridshift3
           shift_tau(:) = gridshift1*rprimd(:,1)/(nr1+1) + gridshift2*rprimd(:,2)/(nr2+1) + gridshift3*rprimd(:,3)/(nr3+1)
         end if

         ABI_ALLOCATE(filename,(1))
         filename(1)=trim(output)
         write(std_out,*) '  The name of your data files is : '
         write(std_out,*) trim(filename(1)),'  for the density (norm of the wfk),'
         write(std_out,*)

         open(unit=unout,file=filename(1),status='replace',form='formatted')
         rewind(unout)
         do iband=1,nband(ckpt)
           write(unout,'(a,2f20.16)')'#', eigen(iband),occ1(iband)
         end do

         write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1+1,nr2+1,nr3+1
         write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', (nr1+1)*(nr2+1)*(nr3+1)
         write(std_out,*) ' znucl = ', znucl, ' typat = ', typat, ' ntypat = ', ntypat

         write(unout,'(1X,A)')  'DIM-GROUP'
         write(unout,*) '3  1'
         write(unout,'(1X,A)') 'PRIMVEC'
         do ir1 = 1,3
           write(unout,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(ir2,ir1), ir2=1,3)
         end do
         write(unout,'(1X,A)') 'PRIMCOORD'
         write(unout,*) natom, ' 1'
!        
!        generate translated coordinates to match density shift
!        
         do iatom = 1,natom
           tau2 (:,iatom) = tau(:,iatom) - shift_tau(:)
         end do

         do iatom = 1,natom
           write(unout,'(i9,3(3X,ES17.10))') int(znucl(typat(iatom))), &
&           Bohr_Ang*tau2(1,iatom), &
&           Bohr_Ang*tau2(2,iatom), &
&           Bohr_Ang*tau2(3,iatom)
         end do
         write(unout,'(1X,A)') 'ATOMS'
         do iatom = 1,natom
           write(unout,'(i9,3(3X,ES17.10))') int(znucl(typat(iatom))), &
&           Bohr_Ang*tau2(1,iatom), &
&           Bohr_Ang*tau2(2,iatom), &
&           Bohr_Ang*tau2(3,iatom)
         end do

!        write(unout,'(1X,A)') 'FRAMES'
         write(unout,'(1X,A)') 'BEGIN_BLOCK_DATAGRID3D'
         write(unout,*) 'datagrids'
         write(unout,'(1X,A)') 'DATAGRID_3D_DENSITY'
         write(unout,*) nr1+1,nr2+1,nr3+1
         write(unout,*) '0.0 0.0 0.0 '
         do ir1 = 1,3
           write(unout,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(ir2,ir1), ir2=1,3)
         end do

         do ir3=gridshift3+1,nr3+1
           ii3=mod(ir3-1,nr3) + 1
           do ir2=gridshift2+1,nr2+1
             ii2=mod(ir2-1,nr2) + 1
             do ir1=gridshift1+1,nr1+1
               ii1=mod(ir1-1,nr1) + 1
               tmpr=fofr(1,ii1,ii2,ii3)
               tmpi=fofr(2,ii1,ii2,ii3)
               write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
             end do
             do ir1=1,gridshift1
               ii1=mod(ir1-1,nr1) + 1
               tmpr=fofr(1,ii1,ii2,ii3)
               tmpi=fofr(2,ii1,ii2,ii3)
               write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
             end do
           end do
           do ir2=1,gridshift2
             ii2=mod(ir2-1,nr2) + 1
             do ir1=gridshift1+1,nr1+1
               ii1=mod(ir1-1,nr1) + 1
               tmpr=fofr(1,ii1,ii2,ii3)
               tmpi=fofr(2,ii1,ii2,ii3)
               write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
             end do
             do ir1=1,gridshift1
               ii1=mod(ir1-1,nr1) + 1
               tmpr=fofr(1,ii1,ii2,ii3)
               tmpi=fofr(2,ii1,ii2,ii3)
               write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
             end do
           end do
         end do
         do ir3=1,gridshift3
           ii3=mod(ir3-1,nr3) + 1
           do ir2=gridshift2+1,nr2+1
             ii2=mod(ir2-1,nr2) + 1
             do ir1=gridshift1+1,nr1+1
               ii1=mod(ir1-1,nr1) + 1
               tmpr=fofr(1,ii1,ii2,ii3)
               tmpi=fofr(2,ii1,ii2,ii3)
               write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
             end do
             do ir1=1,gridshift1
               ii1=mod(ir1-1,nr1) + 1
               tmpr=fofr(1,ii1,ii2,ii3)
               tmpi=fofr(2,ii1,ii2,ii3)
               write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
             end do
           end do
           do ir2=1,gridshift2
             ii2=mod(ir2-1,nr2) + 1
             do ir1=gridshift1+1,nr1+1
               ii1=mod(ir1-1,nr1) + 1
               tmpr=fofr(1,ii1,ii2,ii3)
               tmpi=fofr(2,ii1,ii2,ii3)
               write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
             end do
             do ir1=1,gridshift1
               ii1=mod(ir1-1,nr1) + 1
               tmpr=fofr(1,ii1,ii2,ii3)
               tmpi=fofr(2,ii1,ii2,ii3)
               write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
             end do
           end do
         end do


         write(unout,*)
         write(unout,'(1X,A)') 'END_DATAGRID_3D'
         write(unout,'(1X,A)') 'END_BLOCK_DATAGRID3D'
         close(unout)

         ABI_DEALLOCATE(filename)

         write(std_out,*)
         exit


       case(12)
#if defined HAVE_TRIO_NETCDF
         ABI_ALLOCATE(partwf,(nr1,nr2,nr3))
         ABI_ALLOCATE(filename,(1))
         filename(1)=trim(output)//'.nc'
         write(std_out,*) '  The name of your data files is : '
         write(std_out,*) trim(filename(1)),' is your NetCDF file'
         write(std_out,*)

!        Creating NetCDF file
         ncerr = nf90_create(filename(1), nf90_clobber, ncid)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Creating file")

!        Ask for a title

         write(std_out,*) 'Do you want a title in your NetCDF file? (0 = NO, 1 = YES)'

         read(*,*) titlechoice
         do
           if (titlechoice ==0 .or. titlechoice ==1) exit
           write(std_out,*) 'The answer is not correct, you must enter a integer between 0 and 1 (0 = NO, 1 = YES)'
           read(*,*) titlechoice
         end do
         if (titlechoice ==1) then
           write(std_out,*) 'Enter your file''s title'
           read(*,'(A)') filetitle
         else
           write(std_out,*) 'No title will be add in your NetCDF file'
         end if

         kptvar = kpt(1:3,ckpt)
         originatt(1:3,1:3)=0

         do igrid = 1,3
           gridwavefun1(igrid,1)=0
           gridwavefun2(igrid,1)=0
           gridwavefun3(igrid,1)=0
           gridwavefun1(igrid,2)=Bohr_Ang*rprimd(igrid,3)/nr3
           gridwavefun2(igrid,2)=Bohr_Ang*rprimd(igrid,2)/nr2
           gridwavefun3(igrid,2)=Bohr_Ang*rprimd(igrid,1)/nr1
         end do

!        Defining dimensions

         ncerr = nf90_def_dim(ncid,"gridsize1",nr1, gridsize1DimID)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

         ncerr = nf90_def_dim(ncid,"gridsize2",nr2, gridsize2DimID)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

         ncerr = nf90_def_dim(ncid,"gridsize3",nr3, gridsize3DimID)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

         ncerr = nf90_def_dim(ncid, "lat",3, latDimID)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

         ncerr = nf90_def_dim(ncid, "pos",2, posDimID)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

         ncerr = nf90_def_dim(ncid, "nbatom",natom, nbatomDimID)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

!        Defining variables

         ncerr = nf90_def_var(ncid, "kpoint",nf90_float, latDimID, kpointVarID)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "latticevec",nf90_float, (/ latDimID, latDimID /), latticevecVarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "origin",nf90_float, (/ latDimID, latDimID /), originVarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "atomposi",nf90_float, (/ latDimID, nbatomDimID /), atomposiVarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "atomicnum",nf90_float, nbatomDimID, atomicnumVarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "grid1",nf90_float, (/latDimID,posDimID/), grid1VarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "grid2",nf90_float, (/latDimID,posDimID/), grid2VarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "grid3",nf90_float, (/latDimID,posDimID/), grid3VarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "imagwavefunction",nf90_float, &
&         (/ gridsize1DimID, gridsize2DimID, gridsize3DimID /), imagwavefunVarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

         ncerr = nf90_def_var(ncid, "realwavefunction",nf90_float, &
&         (/ gridsize1DimID, gridsize2DimID, gridsize3DimID /), realwavefunVarID)
         if (ncerr /= nf90_noerr) call  handle_ncerr(ncerr, "Defining variables")

!        Defining attributes

         ncerr = nf90_put_att(ncid,latticevecVarID , "field", "latticevec,vector")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,originVarID , "field", "origin,vector")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,atomposiVarID , "field","atomposi,vector")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,latticevecVarID , "positions","origin" )
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,atomicnumVarID , "positions", "atomposi")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,realwavefunVarID , "positions", &
&         "grid1,product,compact;grid2,product,compact;grid3,product,compact")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,realwavefunVarID , "connections", (/nr3,nr2,nr1/))
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,imagwavefunVarID , "positions", &
&         "grid1,product,compact;grid2,product,compact;grid3,product,compact")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,imagwavefunVarID , "connections", (/nr3,nr2,nr1/))
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,kpointVarID , "long_name", "K-point Value")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,latticevecVarID , "long_name", "Lattice Vectors")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,atomposiVarID , "long_name", "Atomic Positions")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid, atomicnumVarID, "long_name", "Atomic Numbers")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid, imagwavefunVarID, "long_name", " Imaginary Part Wave Function")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid, realwavefunVarID, "long_name", " Real Part Wave Function")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid, originVarID, "long_name", "Origin of the Lattice")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,latticevecVarID, "units", "Angstroms")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         ncerr = nf90_put_att(ncid,atomposiVarID, "units", "Angstroms")
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

         if (titlechoice ==1) then
           ncerr = nf90_put_att(ncid, nf90_global, "title", filetitle)
           if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")
         end if

!        Add the creation date
         call date_and_time(strdat,strtime,strzone,values)
         yyyy=values(1)
         mm=values(2)
         dd=values(3)
         write(stridate(1:2),'(I2)') dd
         stridate(3:3)=" "
         stridate(4:7)=monnam(mm)
         write(stridate(8:11),'(I4)') yyyy

         ncerr = nf90_put_att(ncid, nf90_global,"date", stridate)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

!        Ending the define mode and entering data mode

         ncerr = nf90_enddef(ncid)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Entering data mode")

!        Defining the variables

         ncerr = nf90_put_var(ncid,kpointVarID,kptvar)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         ncerr = nf90_put_var(ncid,latticevecVarID,Bohr_Ang*rprimd)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         ncerr = nf90_put_var(ncid,originVarID,originatt)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         ncerr = nf90_put_var(ncid,atomposiVarID, Bohr_Ang*tau)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         ncerr = nf90_put_var(ncid,atomicnumVarID,znucl(typat))
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         ncerr = nf90_put_var(ncid,grid1VarID,gridwavefun1)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         ncerr = nf90_put_var(ncid,grid2VarID,gridwavefun2)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         ncerr = nf90_put_var(ncid,grid3VarID,gridwavefun3)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         partwf = fofr(2,1:nr1,1:nr2,1:nr3)

         ncerr = nf90_put_var(ncid,imagwavefunVarID, partwf)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         partwf = fofr(1,1:nr1,1:nr2,1:nr3)

         ncerr = nf90_put_var(ncid,realwavefunVarID, partwf)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

         ncerr = nf90_close(ncid)
         if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Closing file")


         write(std_out,*) 'The NetCDF file is done'
!        The NetCDF file is done

         ABI_DEALLOCATE(partwf)
         ABI_DEALLOCATE(filename)

         write(std_out,*)
         exit

#else
         write(std_out,*) 'NetCDF is not defined. You must choose another option'
         exit
#endif

!        ************************************************************

       case(13)
         write(std_out,*)
         write(std_out,*) 'Give 1 files of formatted data'
         write(std_out,*) 'The files are ready to be used with XCrysDen'
         write(std_out,*)
         gridshift1 = 0
         gridshift2 = 0
         gridshift3 = 0
         write(std_out,*) 'Do you want to shift the grid along the x,y or z axis (y/n)?'
         write(std_out,*)
         shift_tau(:) = 0.0
         read (*,*) outputchar
         if (outputchar == 'y' .or. outputchar == 'Y') then
           write(std_out,*) 'Give the three shifts (x,y,z < ',nr1,nr2,nr3,') :'
           write(std_out,*)
           read (*,*) gridshift1, gridshift2, gridshift3
           shift_tau(:) = gridshift1*rprimd(:,1)/(nr1+1) + gridshift2*rprimd(:,2)/(nr2+1) + gridshift3*rprimd(:,3)/(nr3+1)
         end if

         ABI_ALLOCATE(filename,(1))
         filename(1)=trim(output)
         write(std_out,*) '  The name of your data files is : '
         write(std_out,*) trim(filename(1)),'  for the density (norm of the wfk),'
         write(std_out,*)

!        open(unit=unout,file=filename(1),status='replace',form='formatted')
         open(unit=unout,file=filename(1),status='unknown',form='formatted')
         rewind(unout)
         do iband=1,nband(ckpt)
           write(unout,'(a,2f20.16)')'#', eigen(iband),occ1(iband)
         end do

         write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1+1,nr2+1,nr3+1
         write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', (nr1+1)*(nr2+1)*(nr3+1)
         write(std_out,*) ' znucl = ', znucl, ' typat = ', typat, ' ntypat = ', ntypat

         write(unout,'(1X,A)')  'DIM-GROUP'
         write(unout,*) '3  1'
         write(unout,'(1X,A)') 'PRIMVEC'
         do ir1 = 1,3
           write(unout,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(ir2,ir1), ir2=1,3)
         end do
         write(unout,'(1X,A)') 'PRIMCOORD'
         write(unout,*) natom, ' 1'
!        
!        generate translated coordinates to match density shift
!        
         do iatom = 1,natom
           tau2 (:,iatom) = tau(:,iatom) - shift_tau(:)
         end do

         do iatom = 1,natom
           write(unout,'(i9,3(3X,ES17.10))') int(znucl(typat(iatom))), &
&           Bohr_Ang*tau2(1,iatom), &
&           Bohr_Ang*tau2(2,iatom), &
&           Bohr_Ang*tau2(3,iatom)
         end do
         write(unout,'(1X,A)') 'ATOMS'
         do iatom = 1,natom
           write(unout,'(i9,3(3X,ES17.10))') int(znucl(typat(iatom))), &
&           Bohr_Ang*tau2(1,iatom), &
&           Bohr_Ang*tau2(2,iatom), &
&           Bohr_Ang*tau2(3,iatom)
         end do

!        write(unout,'(1X,A)') 'FRAMES'
         write(unout,'(1X,A)') 'BEGIN_BLOCK_DATAGRID3D'
         write(unout,*) 'datagrids'
         write(unout,'(1X,A)') 'DATAGRID_3D_DENSITY'
         write(unout,*) nr1,nr2,nr3
         write(unout,*) '0.0 0.0 0.0 '
         do ir1 = 1,3
           write(unout,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(ir2,ir1), ir2=1,3)
         end do

         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               write(unout,'(ES17.10)') fofr(1,ir1,ir2,ir3)
             end do
           end do
         end do
         write(unout,*)
         write(unout,'(1X,A)') 'END_DATAGRID_3D'
         write(unout,'(1X,A)') 'END_BLOCK_DATAGRID3D'
         close(unout)

         ABI_DEALLOCATE(filename)

         write(std_out,*)
         exit

       case(14)            ! CUBE file format from GAUSSIAN

         write(std_out,*)
         write(std_out,*) 'Output a cube file of 3D volumetric data'
         write(std_out,*)

!        EXAMPLE FROM THE WEB
!        CPMD CUBE FILE.
!        OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
!        3    0.000000    0.000000    0.000000
!        40    0.283459    0.000000    0.000000
!        40    0.000000    0.283459    0.000000
!        40    0.000000    0.000000    0.283459
!        8    0.000000    5.570575    5.669178    5.593517
!        1    0.000000    5.562867    5.669178    7.428055
!        1    0.000000    7.340606    5.669178    5.111259
!        -0.25568E-04  0.59213E-05  0.81068E-05  0.10868E-04  0.11313E-04  0.35999E-05


         open(unit=unout,file=output,status='replace',form='formatted')
         write(unout,'(a)') 'ABINIT generated cube file'
         write(unout,'(a)') 'from cut3d tool'

         write(unout,'(i9,3(1x,f12.6))') natom,0.,0.,0.
         write(unout,'(i9,3(1x,f12.6))') nr1,(rprimd(ir2,1)/nr1, ir2=1,3)
         write(unout,'(i9,3(1x,f12.6))') nr2,(rprimd(ir2,2)/nr2, ir2=1,3)
         write(unout,'(i9,3(1x,f12.6))') nr3,(rprimd(ir2,3)/nr3, ir2=1,3)

         do iatom = 1,natom
           write(unout,'(i9,4(3X,ES17.10))') int(znucl(typat(iatom))),0.d0, &
&           tau(1,iatom), &
&           tau(2,iatom), &
&           tau(3,iatom)
         end do

!        C ordering of the indexes 
         do ir1=1,nr1
           do ir2=1,nr2
             do ir3=1,nr3
               write(unout,'(6(f12.6,2x))') sqrt( fofr(1,ir1,ir2,ir3)**2 + fofr(2,ir1,ir2,ir3)**2 )
             end do
           end do
         end do

         close(unout)
         exit

       case(0)
         write(std_out,*)' Exit inner loop'
!        stop
         select_exit = 1


         case default
         write(std_out,*) ' This choice is not valid.'
         write(std_out,*)
         cycle

     end select

   end do

   ckpt=oldckpt
   cband=oldcband
   csppol=oldcsppol
   cspinor=oldcspinor
!  deallocate the datas
   ABI_DEALLOCATE(fofr)

   write(std_out,*) ' Task ',ichoice,' has been done !'
   write(std_out,*)
   write(std_out,*) ' Run interpolation again? (1=default=yes,0=no)'
   read(*,*) iprompt
   if(iprompt==0) then
     exit
   else
     cycle
   end if
 end do
!deallocate the datas
 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(kg_dum)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(occ1)

!Close the WF file
 call WffClose(wff, ierr)

end subroutine wffile
!!***
