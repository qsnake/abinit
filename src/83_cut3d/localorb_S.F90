!{\src2tex{textfont=tt}}
!!****f* ABINIT/localorb_S
!! NAME
!! localorb_S
!!
!! FUNCTION
!! Construction of Wannier type localized orbitals (WLO)
!! Ref:cond-mat/0509571
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JB)
!! this file is distributed under the terms of the
!! gnu general public license, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! Needs an unformatted wave function from abinit.
!! chr_inputfname = name of the input file for WLO construction
!! ecut = effective ecut (ecut*dilatmx**2)
!! exchn2n3d=if 1, n2 and n3 are exchanged
!! headform=format of the wf file
!! istwfk = input variable indicating the storage option of each k-point
!! kpt = input variable kptns
!! natom = number of atoms in the unit cell
!! nband = input variable nband(nkpt*nsppol) !! MC 090827: the definition below in not consistent with what we find in defs_datatypes.F90
!! nkpt = number of k-points
!! npwarr = array holding npw for each k point
!! nr1,nr2,nr3 = grid size (nr1 x nr2 x nr3 = filrho dimension)
!! nspinor = number of spinorial components of the wavefunctions
!! nsppol = number of spin polarization
!! ntypat = number of atom type
!! paral_kgb = parallization option, it is set to 0 in the parent subroutine 
!! rprim = orientation of the unit cell axes
!! tau = cartesian coordinates (real space)
!! typat = input variable typat(natom)
!! znucl = znucltypat(ntypat) from alchemy
!!
!! OUTPUT
!! xcrysden plottable WLO file,
!! corresponding complex data file
!! other analysis and runtime check files

!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!      instrng,intagm,inupper,kptindex,metric,overlap_wf,wfread,xredxcart
!!      zgesvd
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine localorb_s(chr_inputfname,ecut,exchn2n3d,headform,istwfk,kpt,natom,nband,nkpt,&
                       npwarr,nr1,nr2,nr3,nspinor,nsppol,ntypat,paral_kgb,rprimd,tau,typat,znucl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'localorb_s'
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_42_parser
 use interfaces_83_cut3d, except_this_one => localorb_s
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,headform,natom,nkpt,nr1,nr2,nr3,nspinor,nsppol
 integer,intent(in) :: ntypat,paral_kgb
 real(dp),intent(in) :: ecut
 character(len=50),intent(in) :: chr_inputfname
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt),npwarr(nkpt),typat(natom) !! MC 090827: Should it be nband(nkpt*nsppol) instead?
 real(dp),intent(in) :: kpt(3,nkpt),rprimd(3,3),znucl(ntypat)
 real(dp),intent(inout) :: tau(3,natom)

!Local variables-------------------------------
!,unout=12
       character(*), parameter :: inputfile='cut.in'
! for input file reading :
! for wf :
!scalars
 integer :: ckpt,cspinor,csppol,formeig
 integer :: i,i1,i2,i3
 integer :: iatom,iband1,iband2,icellplotx,icellploty,icellplotz
 integer :: icellx,icelly,icellz,icountorb
 integer :: igx,igx1,igy,igy1,igz,igz1,ii1,ii2,ii3,ikpt_x
 integer :: ikpt_y,ikpt_z,info
 integer :: insmet,iorb,iout
 integer :: iwan,iwrite,j
 integer :: j1,k,l,lenstr,m,marr,matom
 integer :: munitx,munity,munitz,n
 integer :: nbands,ncellx,ncellx1,ncelly,ncelly1,ncellz,ncellz1,ndimcsvd
 integer :: ngx,ngy,ngz,njob,nkx,nky,nkz,norb
 integer :: nstartx,nstartx1,nstarty,nstarty1
 integer :: nstartz,nstartz1,nstopx,nstopx1,nstopy,nstopy1,nstopz,nstopz1,nvar
 integer :: nwanmax1,nwanmax2,nwanmax3,nwantot,nxhw,nxmax,nxmin
 integer :: nxtot,nxw,nyhw,nymax,nymin,nytot,nyw,nzhw,nzmax,nzmin,nztot,nzw
 integer :: oldcband,oldckpt,oldcspinor,oldcsppol
 integer :: tread
 real(dp) :: alpha,efermi,envelop,hx,hy
 real(dp) :: hz,normtot,r1,r2,r3,rcellx,rcelly,rcellz
 real(dp) :: rcoordx,rcoordy,rcoordz,rcos1,rcos2,rcos3,rcosphy,rcostheta,re1
 real(dp) :: re2,re3,re4,re5,re6,re7,rk,rkpt1,rkpt2,rkpt3
 real(dp) :: rkx,rkx1,rkx2,rkx3,rky,rky1,rky2,rky3,rkz,rkz1
 real(dp) :: rkz2,rkz3,rlx,rlx1,rlx1_large,rlx2,rlx2_large,rlx3,rlx3_large
 real(dp) :: rlx_large,rly,rly1,rly1_large,rly2,rly2_large,rly3,rly3_large
 real(dp) :: rly_large,rlz,rlz1,rlz1_large,rlz2,rlz2_large,rlz3,rlz3_large
 real(dp) :: rlz_large,rn1,rn2,rn3,rnorm,roccup,rpi,rrcosphy2,rrsinphy2
 real(dp) :: rsin1,rsin2,rsin3,rsinphy,rsintheta,rwanmax,rxcm,rycm
 real(dp) :: rzcm,shiftx,shifty,shiftz,thetax,thetay,thetaz,tpi
 real(dp) :: tsmear,ucvol,x,x1,y,y1,z,z1
 character(len=1) :: crhsvd
 character(len=30) :: token
 character(len=strlen) :: string
 type(mpi_type) :: mpi_enreg
!arrays
 integer,allocatable :: intarr(:),ltypeorb(:,:),numorb(:,:),nwanindex(:)
 real(dp) :: gmet(3,3),gprimd(3,3)
 real(dp) :: rmet(3,3),vec(3)
 real(dp) :: vec1(3),xcart(3,natom),xmat(3,3),xred(3,natom),ymat(3,3),zmat(3,3)
 real(dp),allocatable :: ai(:,:),ar(:,:),dprarr(:),e_kpt(:),fv1(:),fv2(:)
 real(dp),allocatable :: fv3(:),rcoord(:,:),rcoordall(:,:,:),rprojx(:)
 real(dp),allocatable :: rprojy(:),rprojz(:),rtheta(:,:,:),rwork(:),s(:),wi(:)
 real(dp),allocatable :: wr(:),zi(:,:),zr(:,:)
 character(len=15),allocatable :: chrhead(:)
!no_abirules
 complex(dp) :: cmp1,cmp2,csum
 complex(kind =8)::orbital,corbital
 complex(dp),allocatable,dimension(:,:)::      cm1,cm2,cm4,cm5,coverlap
 complex(dp),allocatable,dimension(:,:)::      cmat1,cmat2
 complex(dp),allocatable,dimension(:,:) ::     cu,cvt
 complex(dp),allocatable:: csvd(:)
 complex(dp),allocatable ::cwannier(:,:,:) 
 complex(dp),allocatable ::crotl(:,:,:,:,:)
 complex(dp),allocatable :: cbfns(:,:,:),cbfns1(:,:,:),corb(:,:,:)

! *************************************************************************

 nbands = maxval(nband)
 ngz = nr3
 ngy = nr2
 ngx = nr1
 tpi=two_pi
 ndimcsvd = 3*nbands + 10
 ABI_ALLOCATE(csvd,(ndimcsvd))

!Finding out nkx,nky,nkz from the kpt array:

 nkx = 1
 nky = 1
 nkz = 1

 do i = 1,nkpt

   re1 = kpt(1,i)
   re2 = kpt(2,i)
   re3 = kpt(3,i)

   if(i.gt.1)then
     do j=1,i-1
       if(re1.eq.kpt(1,j))go to 100
     end do
     nkx = nkx + 1
     100 continue
   end if

   if(i.gt.1)then
     do j =1,i-1
       if(re2.eq.kpt(2,j))go to 200
     end do
     nky = nky + 1
     200 continue
   end if

   if(i.gt.1)then
     do j =1,i-1
       if(re3.eq.kpt(3,j))go to 300
     end do
     nkz = nkz + 1
     300 continue
   end if
 end do ! i=1,nkpt

!==================================================================
!reading the input file for wlo :
!default:
 matom=0
 insmet = 1
 iwrite = 0
 alpha = 1.0_dp
 ncellx1 = 1
 ncelly1 = 1
 ncellz1 = 1
 norb = nbands
 rcellx = 0.0_dp
 rcelly = 0.0_dp
 rcellz = 0.0_dp
 njob = 0
 tsmear=zero
 efermi=zero

 nvar = 13 ! maximum number of variables to be read from the input files
 ABI_ALLOCATE(chrhead,(nvar))
 chrhead(1)='ncenter'
 chrhead(2)='norb'
 chrhead(3)='supercell'
 chrhead(4)='alpha'
 chrhead(5)='lofwrite'
 chrhead(6)='shiftk'
 chrhead(7)='insmet'
 chrhead(8)='efermi'
 chrhead(9)='tsmear'
 chrhead(10)='orbdetails'
 chrhead(11)='wfconstrctn'
 chrhead(12)='wfcell'
 chrhead(13)='jobtype'

 call instrng(chr_inputfname,lenstr,1,strlen,string)
!To make case-insensitive, map characters of string to upper case:
 call inupper(string(1:lenstr))

 marr=3
 ABI_ALLOCATE(intarr,(3))
 ABI_ALLOCATE(dprarr,(3))

 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),chrhead(1),tread,'INT')
 if(tread==1) matom=intarr(1)

 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),chrhead(2),tread,'INT')
 if(tread==1) norb=intarr(1)

 if(norb.gt.nbands)then
   write(std_out,*) 'chosen number of bands :',norb
   write(std_out,*) 'nbands in scf :',nbands
   write(std_out,*) 'the former cannot be greater than the latter.'
   stop
 end if
 nbands = norb

 call intagm(dprarr,intarr,0,marr,3,string(1:lenstr),chrhead(3),tread,'INT')
 if(tread==1)then
   ncellx1=intarr(1)
   ncelly1=intarr(2)
   ncellz1=intarr(3)
 end if

 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),chrhead(4),tread,'DPR')
 if(tread==1) alpha=dprarr(1)

 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),chrhead(5),tread,'INT')
 if(tread==1) iwrite=intarr(1)

 call intagm(dprarr,intarr,0,marr,3,string(1:lenstr),chrhead(6),tread,'DPR')
 if(tread==1)then
   shiftx=dprarr(1)
   shifty=dprarr(2)
   shiftz=dprarr(3)
 end if

 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),chrhead(7),tread,'INT')
 if(tread==1) insmet=intarr(1)

 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),chrhead(8),tread,'DPR')
 if(tread==1) efermi=dprarr(1)

 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),chrhead(9),tread,'DPR')
 if(tread==1) tsmear=dprarr(1)

 call intagm(dprarr,intarr,0,marr,3,string(1:lenstr),chrhead(12),tread,'DPR')
 if(tread==1)then
   rcellx=intarr(1)
   rcelly=intarr(2)
   rcellz=intarr(3)
 end if

 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),chrhead(13),tread,'INT')
 if(tread==1) njob=intarr(1)

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 marr=3*nsppol*matom
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 ABI_ALLOCATE(numorb,(nsppol,matom))
 ABI_ALLOCATE(rcoord,(3,matom))
 ABI_ALLOCATE(ltypeorb,(nbands,nsppol))
 ABI_ALLOCATE(e_kpt,(nbands))
 ABI_ALLOCATE(rcoordall,(3,nbands,nsppol))

 token = 'rcoord'
 call intagm(dprarr,intarr,0,marr,3*matom,string(1:lenstr),token,tread,'DPR')
 if(tread==1)rcoord(:,:)=reshape(dprarr(1:3*matom),(/3,matom/))

 token = 'numorb'
 call intagm(dprarr,intarr,0,marr,nsppol*matom,string(1:lenstr),token,tread,'INT')
 if(tread==1)numorb(:,:)=reshape(intarr(1:nsppol*matom),(/nsppol,matom/))

 ABI_ALLOCATE(rtheta,(nbands,3,nsppol))
 ii1=sum(numorb(1:nsppol,1:matom))
 if(nbands/=ii1)then
   write(std_out,'(a,a,a,i4,a,i4,a,a)')' localorb_S : ERROR ',ch10,&
&   ' The sum of the number of orbitals ',ii1,', is not equal to the number of bands =',nbands,ch10,&
&   ' Stop '
   stop
 end if

 do csppol = 1,nsppol
   ii1=0 ; ii2=1
   do i = 1,matom
     j1=numorb(csppol,i) ; ii1 = ii1 + j1
     do k = ii2,ii1
       rcoordall(1:3,k,csppol)=rcoord(1:3,i)
     end do
     ii2 = ii2 +j1
   end do
 end do

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 marr=nbands*3*nsppol
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 token = 'ltypeorb'
 call intagm(dprarr,intarr,0,marr,nbands*nsppol,string(1:lenstr),token,tread,'INT')
 if(tread==1)ltypeorb(:,:)=reshape(intarr(1:nbands*nsppol),(/nbands,nsppol/))

 token = 'rtheta'
 call intagm(dprarr,intarr,0,marr,nbands*3*nsppol,string(1:lenstr),token,tread,'DPR')
 if(tread==1)rtheta(:,:,:)=reshape(dprarr(1:nbands*3*nsppol),(/nbands,3,nsppol/))

 nwantot = nbands
 token = 'nwantot'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1)nwantot=intarr(1)
 ABI_ALLOCATE(nwanindex,(nwantot))

 if(nwantot == nbands)then
   do i1 = 1,nwantot
     nwanindex(i1) = i1
   end do
 else
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   marr=nwantot
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
   token = 'nwanindex'
   call intagm(dprarr,intarr,0,marr,nwantot,string(1:lenstr),token,tread,'INT')
   if(tread==1)nwanindex(1:nwantot)=intarr(1:nwantot)
 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
!The reading of data is finished
!--------------------------------------------------------------------------------------

 if(rcellx+rcelly+rcellz.ne.0.0)then
   if( anint(rcellx).gt.ncellx1)rcellx=real(ncellx1)
   if( anint(rcelly).gt.ncelly1)rcelly=real(ncelly1)
   if( anint(rcellz).gt.ncellz1)rcellz=real(ncellz1)
 end if

 write(std_out,'(3a)')ch10,' Echoing out the read input variables for WLO calculation:',ch10

 write(std_out,'(a,a12,i4)')' ',chrhead(1),matom
 write(std_out,'(a,a12,i4)')' ',chrhead(2),nbands
 write(std_out,'(a,a12,3i4)')' ',chrhead(3),ncellx1,ncelly1,ncellz1
 if(rcellx+rcelly+rcellz.eq.0.0)then
   write(std_out,'(a,a12,3i4)')' ',chrhead(12),ncellx1,ncelly1,ncellz1
 else
   write(std_out,'(a,a12,3i4)')' ',chrhead(12),rcellx,rcelly,rcellz
 end if
 write(std_out,'(a,a12,es16.6)')' ',chrhead(4),alpha
 write(std_out,'(a,a12,i4)')' ',chrhead(5),iwrite
 write(std_out,'(a,a12,3es16.6)')' ',chrhead(6),shiftx,shifty,shiftz
 write(std_out,'(a,a12,i4)')' ',chrhead(7),insmet
 write(std_out,'(a,a12,es16.6)')' ',chrhead(8),efermi
 write(std_out,'(a,a12,es16.6)')' ',chrhead(9),tsmear
 write(std_out,'(a)')' For each atom :        rcoord(1:3)                                 numorb'
 do i=1,matom
   write(std_out,'(8x,i4,3es16.6,8x,2i3)')i,rcoord(:,i),numorb(1:nsppol,i)
 end do
 write(std_out,'(a)')' For each orbital : spin, ltypeorb     and      rtheta(1:3) '
 do csppol=1,nsppol
   do ii2=1,nbands
     write(std_out,'(8x,i4,2i10,4x,3es16.6)')ii2,csppol,ltypeorb(ii2,csppol),rtheta(ii2,:,csppol)
   end do
 end do
 write(std_out,'(1x,a12,i4)')chrhead(11),nwantot
 write(std_out,'(a)')' They are :'
 write(std_out,'(a,20i4)')' ',(nwanindex(i1),i1 =1,nwantot)

!-----------------------------------------------------------------------------
 ABI_ALLOCATE( corb,(ngx,ngy,ngz))
 ABI_ALLOCATE(cbfns,(ngx,ngy,ngz))
 ABI_ALLOCATE( cbfns1,(ngx,ngy,ngz))

!svd :
 ABI_ALLOCATE(coverlap,(nbands,nbands))
 ABI_ALLOCATE(cm1,(nbands,nbands))
 ABI_ALLOCATE(cm2,(nbands,nbands))
 ABI_ALLOCATE(cm4,(nbands,nbands))
 ABI_ALLOCATE(cm5,(nbands,nbands))
 ABI_ALLOCATE(cmat1,(nbands,nbands))
 ABI_ALLOCATE(cmat2,(nbands,nbands))
 ABI_ALLOCATE(s,(nbands))
 ABI_ALLOCATE(cu,(nbands,nbands))
 ABI_ALLOCATE(cvt,(nbands,nbands))
 ABI_ALLOCATE(rwork,(5*nbands))
 ABI_ALLOCATE(ar,(nbands,nbands))
 ABI_ALLOCATE(ai,(nbands,nbands))
 ABI_ALLOCATE(zr,(nbands,nbands))
 ABI_ALLOCATE(zi,(nbands,nbands))
 ABI_ALLOCATE(fv1,(nbands))
 ABI_ALLOCATE(fv2,(nbands))
 ABI_ALLOCATE(fv3,(nbands))
 ABI_ALLOCATE(wr,(nbands))
 ABI_ALLOCATE(wi,(nbands))

 ABI_ALLOCATE( cwannier,(ncellx1*ngx,ncelly1*ngy,ncellz1*ngz))
 ABI_ALLOCATE( crotl,(nkx,nky,nkz,nbands,nbands))
!-------------------------------------------------------------------------------
!begin executable section

 mpi_enreg%paralbd=0

 formeig=0
 oldckpt=0
 oldcband=0
 oldcsppol=0
 oldcspinor=0

 iout=-1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

!get xred
 call xredxcart(natom,-1,rprimd,tau,xred)

 do i=1,natom
   xcart(:,i)=rprimd(:,1)*xred(1,i)+&
&   rprimd(:,2)*xred(2,i)+&
&   rprimd(:,3)*xred(3,i)
 end do

 tpi = 8.0*atan(1.0)
 rpi  = 4.0*atan(1.0)

 ngz = nr3
 ngy = nr2
 ngx = nr1

 if(ncellx1.ge.2)icellplotx = 2
 if(ncellx1.lt.2)icellplotx = 1

 if(ncelly1.ge.2)icellploty = 2
 if(ncelly1.lt.2)icellploty = 1

 if(ncellz1.ge.2)icellplotz = 2
 if(ncellz1.lt.2)icellplotz = 1

 rlx1_large=rprimd(1,1)*real(ncellx1)
 rly1_large=rprimd(2,1)*real(ncellx1)
 rlz1_large=rprimd(3,1)*real(ncellx1)

 rlx2_large=rprimd(1,2)*real(ncelly1)
 rly2_large=rprimd(2,2)*real(ncelly1)
 rlz2_large=rprimd(3,2)*real(ncelly1)

 rlx3_large=rprimd(1,3)*real(ncellz1)
 rly3_large=rprimd(2,3)*real(ncellz1)
 rlz3_large=rprimd(3,3)*real(ncellz1)

 rlx_large = rlx1_large + rlx2_large + rlx3_large
 rly_large = rly1_large + rly2_large + rly3_large
 rlz_large = rlz1_large + rlz2_large + rlz3_large

 rlx1=rprimd(1,1)
 rly1=rprimd(2,1)
 rlz1=rprimd(3,1)

 rlx2=rprimd(1,2)
 rly2=rprimd(2,2)
 rlz2=rprimd(3,2)

 rlx3=rprimd(1,3)
 rly3=rprimd(2,3)
 rlz3=rprimd(3,3)

 rlx = rlx1 + rlx2 + rlx3
 rly = rly1 + rly2 + rly3
 rlz = rlz1 + rlz2 + rlz3

 write(std_out,'(a,3i6)' )' real space grids:',ngx,ngy,ngz
 write(std_out,'(a,3i6)' )' k space grids:',nkx,nky,nkz
 write(std_out,'(a,3es16.6)' )' rlx rly rlz',rlx,rly,rlz

!=============================================================================
!reciprocal space vectors :
!----------------------------
 rkx1 = rly2*rlz3 - rly3*rlz2
 rky1 = rlz2*rlx3 - rlz3*rlx2
 rkz1 = rlx2*rly3 - rlx3*rly2
 rk   = rkx1*rlx1 + rky1*rly1 + rkz1*rlz1
 rkx1 = tpi*rkx1/rk
 rky1 = tpi*rky1/rk
 rkz1 = tpi*rkz1/rk

 rkx2 = rly3*rlz1 - rly1*rlz3
 rky2 = rlz3*rlx1 - rlz1*rlx3
 rkz2 = rlx3*rly1 - rlx1*rly3
 rk   = rkx2*rlx2 + rky2*rly2 + rkz2*rlz2
 rkx2 = tpi*rkx2/rk
 rky2 = tpi*rky2/rk
 rkz2 = tpi*rkz2/rk

 rkx3 = rly1*rlz2 - rly2*rlz1
 rky3 = rlz1*rlx2 - rlz2*rlx1
 rkz3 = rlx1*rly2 - rlx2*rly1
 rk   = rkx3*rlx3 + rky3*rly3 + rkz3*rlz3
 rkx3 = tpi*rkx3/rk
 rky3 = tpi*rky3/rk
 rkz3 = tpi*rkz3/rk
!=============================================================================
 hx=rprimd(1,1)/real(ngx)
 hy=rprimd(2,2)/real(ngy)
 hz=rprimd(3,3)/real(ngz)

 open(26,file="crot.dat",form='unformatted')
 open(29,file="overlap_abs.dat",form='formatted')
 open(30,file="singularvalues.dat",form='formatted')
 open(31,file="lo.status",form='formatted')
 open(6661,file="lo_center.dat",form='formatted')
 open(6662,file="lo_loclength.dat",form='formatted')

 rewind(26)
 rewind(29)
 rewind(30)
 rewind(6661)
 rewind(6662)

 write(29,*)' in the overlap matrix each row correspond to each orbital'
 write(29,*)' described below orbdetails variable in the input file.'
 write(29,*)' '

!diverting for only WLO construction :
 if(njob .ne. 0) go to 8765

!=====================================================
!Box in which the Gaussian has appreciable magnitude
 nxmax = -10000000
 nymax = -10000000
 nzmax = -10000000
 nxmin = 10000000
 nymin = 10000000
 nzmin = 10000000
 do k = 1,ngz
   do j = 1,ngy
     do i = 1,ngx

       r1 = (real(i)-(1.0-mod(real(ngx),2.0)*0.5) )/real(ngx)
       r2 = (real(j)-(1.0-mod(real(ngy),2.0)*0.5) )/real(ngy)
       r3 = (real(k)-(1.0-mod(real(ngz),2.0)*0.5) )/real(ngz)

       x=rlx*(-0.5)+rlx1*r1+rlx2*r2+rlx3*r3
       y=rly*(-0.5)+rly1*r1+rly2*r2+rly3*r3
       z=rlz*(-0.5)+rlz1*r1+rlz2*r2+rlz3*r3
       envelop=exp(-alpha*x*x)*exp(-alpha*y*y)*exp(-alpha*z*z)
       if(envelop.gt.0.000001)then
         if(i.gt.nxmax)nxmax = i
         if(i.lt.nxmin)nxmin = i
         if(j.gt.nymax)nymax = j
         if(j.lt.nymin)nymin = j
         if(k.gt.nzmax)nzmax = k
         if(k.lt.nzmin)nzmin = k
       end if

     end do
   end do
 end do
 nxw = nxmax - nxmin + 1
 nyw = nymax - nymin + 1
 nzw = nzmax - nzmin + 1

 nxhw = (nxw + anint(mod(real(nxw),2.0)) )/2
 nyhw = (nyw + anint(mod(real(nyw),2.0)) )/2
 nzhw = (nzw + anint(mod(real(nzw),2.0)) )/2

 write(std_out,*) 'Box width outside which localized template is set to zero:'
 write(std_out,*) nxw,' x ',nyw,' x ',nzw
!stop
!##########################################################################
 do csppol = 1,nsppol
   cspinor = 1 ! only scalar wave functions.

   do ikpt_z= 1,nkz
     do ikpt_y= 1,nky
       do ikpt_x= 1,nkx

         rewind(31)
!        ---------------------------------------------
!        k index for ascending order of wave vector k
!        ----------------------------------------------
         call kptindex(shiftx,shifty,shiftz,ikpt_x,ikpt_y,ikpt_z,nkx,nky,nkz,&
&         i,j,k)
         ckpt=i + (j - 1)*nkx + (k - 1)*nkx*nky
         write(31, '(a)' )" Doing : Overlap matrix and its singular value decomposition "
         write(31, '(a,3i4)' )" kpoint index: ",ikpt_x,ikpt_y,ikpt_z
         write(31, '(a,3i4)' )' folded to:',i,j,k

         write(std_out,'(a)' )" Doing : Overlap matrix and its singular value decomposition "
         write(std_out,'(a,3i4)' )" kpoint index: ",ikpt_x,ikpt_y,ikpt_z
         write(std_out,'(a,3i4)' )' folded to:',i,j,k

         rkpt1 = kpt(1,ckpt)
         rkpt2 = kpt(2,ckpt)
         rkpt3 = kpt(3,ckpt)

         if(ikpt_x .eq. nkx+1)rkpt1 = rkpt1 + 1.0
         write(31,*)' kpt:',rkpt1,rkpt2,rkpt3

         rkx = rkpt1*rkx1 + rkpt2*rkx2 + rkpt3*rkx3
         rky = rkpt1*rky1 + rkpt2*rky2 + rkpt3*rky3
         rkz = rkpt1*rkz1 + rkpt2*rkz2 + rkpt3*rkz3

         ncellx = 3
         ncelly = 3  ! we fix 3 unitcells only for the inverse f transform.
         ncellz = 3

         if(nkx .eq.1)ncellx = 1
         if(nky .eq.1)ncelly = 1  ! for molecules or low dimensional systems
         if(nkz .eq.1)ncellz = 1

!        ##############################################################
         do iband2 =  1, nbands  ! FIRST BAND INDEX STARTS HERE
!          ---------------------------------------------
!          write(std_out,*) 'iband:',iband2
!          real space grid :

           rnorm = 0.0
           corb = 0.0
           do icellx = 1,ncellx
             do icelly = 1,ncelly
               do icellz = 1,ncellz

                 rn1 = real(icellx-(ncellx + anint(mod(real(ncellx),2.0)) )/2)
                 rn2 = real(icelly-(ncelly + anint(mod(real(ncelly),2.0)) )/2)
                 rn3 = real(icellz-(ncellz + anint(mod(real(ncellz),2.0)) )/2)

                 munitx = (ngx+anint(mod(real(ngx),2.0)))/2 +&
&                 anint(mod(real(ngx)+1.0,2.0))
                 munity = (ngy+anint(mod(real(ngy),2.0)))/2 +&
&                 anint(mod(real(ngy)+1.0,2.0))
                 munitz = (ngz+anint(mod(real(ngz),2.0)))/2 +&
&                 anint(mod(real(ngz)+1.0,2.0))

                 re1 = rcoordall(1,iband2,csppol)
                 re2 = rcoordall(2,iband2,csppol)
                 re3 = rcoordall(3,iband2,csppol)
                 munitx = munitx + anint(re1*real(ngx))
                 munity = munity + anint(re2*real(ngy))
                 munitz = munitz + anint(re3*real(ngz))

                 nstartx = munitx + anint(rn1)*ngx - nxhw
                 nstopx  = munitx + anint(rn1)*ngx + nxhw

                 nstarty = munity + anint(rn2)*ngy - nyhw
                 nstopy  = munity + anint(rn2)*ngy + nyhw

                 nstartz = munitz + anint(rn3)*ngz - nzhw
                 nstopz  = munitz + anint(rn3)*ngz + nzhw

                 if(nstartx.le.1)nstartx = 1
                 if(nstartx.gt.ngx)nstartx = ngx
                 if(nstarty.le.1)nstarty = 1
                 if(nstarty.gt.ngy)nstarty = ngy
                 if(nstartz.le.1)nstartz = 1
                 if(nstartz.gt.ngz)nstartz = ngz

                 if(nstopx.le.1)nstopx = 1
                 if(nstopx.gt.ngx)nstopx = ngx
                 if(nstopy.le.1)nstopy = 1
                 if(nstopy.gt.ngy)nstopy = ngy
                 if(nstopz.le.1)nstopz = 1
                 if(nstopz.gt.ngz)nstopz = ngz

                 do k = 1,ngz
                   do j = 1,ngy
                     do i = 1,ngx

                       if(i.ge.nstartx.and.i.le.nstopx.and.&
&                       j.ge.nstarty.and.j.le.nstopy.and.&
&                       k.ge.nstartz.and.k.le.nstopz)then

                         r1 = (real(i)-(1.0-mod(real(ngx),2.0)*0.5) )/real(ngx)
                         r2 = (real(j)-(1.0-mod(real(ngy),2.0)*0.5) )/real(ngy)
                         r3 = (real(k)-(1.0-mod(real(ngz),2.0)*0.5) )/real(ngz)

                         x1=rlx*(-0.5)+rlx1*r1+rlx2*r2+rlx3*r3
                         y1=rly*(-0.5)+rly1*r1+rly2*r2+rly3*r3
                         z1=rlz*(-0.5)+rlz1*r1+rlz2*r2+rlz3*r3
!                        ---------------------------------------------
!                        for each atom :
                         icountorb = 0
                         do iatom = 1,matom

                           rcoordx=rcoord(1,iatom)*rlx1+rcoord(2,iatom)*rlx2+rcoord(3,iatom)*rlx3
                           rcoordy=rcoord(1,iatom)*rly1+rcoord(2,iatom)*rly2+rcoord(3,iatom)*rly3
                           rcoordz=rcoord(1,iatom)*rlz1+rcoord(2,iatom)*rlz2+rcoord(3,iatom)*rlz3
!                          ---------------------------------------------
!                          for different orbitals at each atom :
                           do iorb = 1,numorb(csppol,iatom)
                             icountorb = icountorb + 1

!                            *************************************************
                             if(iband2 .eq. icountorb ) then
!                              *************************************************

                               i1 = ltypeorb(icountorb,csppol)
                               thetax = rpi*rtheta(icountorb,1,csppol)/180.0
                               thetay = rpi*rtheta(icountorb,2,csppol)/180.0
                               thetaz = rpi*rtheta(icountorb,3,csppol)/180.0
!                              ---------------------------------------------
                               corbital = 0.0

                               x=rlx*(-0.5)+rlx1*r1+rlx2*r2+rlx3*r3-rcoordx-rn1*rlx1-rn2*rlx2-rn3*rlx3
                               y=rly*(-0.5)+rly1*r1+rly2*r2+rly3*r3-rcoordy-rn1*rly1-rn2*rly2-rn3*rly3
                               z=rlz*(-0.5)+rlz1*r1+rlz2*r2+rlz3*r3-rcoordz-rn1*rlz1-rn2*rlz2-rn3*rlz3

!                              -------------------------------------------
!                              rotating the orbitals as desired :
!                              -----------------------------------------
                               rsin1=sin(thetax)
                               rcos1=cos(thetax)

                               rsin2=sin(thetay)
                               rcos2=cos(thetay)

                               rsin3=sin(thetaz)
                               rcos3=cos(thetaz)
!                              -------------------------------------------
!                              xmat  -->  clockwise rotation wrt x
!                              -------------------------------------------
                               if(abs(thetax) .ne. 0.0)then
                                 xmat(1,1)=1.0
                                 xmat(1,2)=0.0
                                 xmat(1,3)=0.0

                                 xmat(2,1)=0.0
                                 xmat(2,2)=rcos1
                                 xmat(2,3)=-rsin1

                                 xmat(3,1)=0.0
                                 xmat(3,2)=rsin1
                                 xmat(3,3)=rcos1

                                 vec(1) = x ; vec(2) = y ; vec(3) = z
                                 vec1 = matmul(xmat,vec)
                                 x = vec1(1) ; y = vec1(2) ; z = vec1(3)
                               end if
!                              -------------------------------------------
!                              ymat  -->  clockwise rotation wrt y
!                              -------------------------------------------
                               if(abs(thetay) .ne. 0.0)then
                                 ymat(1,1)=rcos2
                                 ymat(1,2)=0.0
                                 ymat(1,3)=-rsin2

                                 ymat(2,1)=0.0
                                 ymat(2,2)=1.0
                                 ymat(2,3)=0.0

                                 ymat(3,1)=rsin2
                                 ymat(3,2)=0.0
                                 ymat(3,3)=rcos2

                                 vec(1) = x ; vec(2) = y ; vec(3) = z
                                 vec1 = matmul(ymat,vec)
                                 x = vec1(1) ; y = vec1(2) ; z = vec1(3)
                               end if
!                              -------------------------------------------
!                              zmat  -->  clockwise rotation wrt z
!                              -------------------------------------------
                               if(abs(thetaz) .ne. 0.0)then
                                 zmat(1,1)=rcos3
                                 zmat(1,2)=-rsin3
                                 zmat(1,3)=0.0

                                 zmat(2,1)=rsin3
                                 zmat(2,2)=rcos3
                                 zmat(2,3)=0.0

                                 zmat(3,1)=0.0
                                 zmat(3,2)=0.0
                                 zmat(3,3)=1.0

                                 vec(1) = x ; vec(2) = y ; vec(3) = z
                                 vec1 = matmul(zmat,vec)

                                 x = vec1(1) ; y = vec1(2) ; z = vec1(3)
                               end if
!                              -------------------------------------------
!                              gaussian envelop :
                               envelop=exp(-alpha*x*x)*exp(-alpha*y*y)*exp(-alpha*z*z)
                               if((x*x + y*y + z*z).eq.0.0)then
                                 rsintheta=0.0
                                 rcostheta=0.0
                               else
                                 rsintheta=sqrt(z*z + y*y)/sqrt(x*x + y*y + z*z)
                                 rcostheta=x/sqrt(x*x + y*y + z*z)
                               end if

                               if((y*y + z*z).eq.0.0)then
                                 rsinphy  =0.0
                                 rcosphy  =0.0
                               else
                                 rsinphy  =z/sqrt(z*z + y*y)
                                 rcosphy  =y/sqrt(z*z + y*y)
                               end if

                               rrsinphy2 =2.0*rsinphy*rcosphy
                               rrcosphy2 =rcosphy*rcosphy-rsinphy*rsinphy
!                              ------------------------------------------------------
!                              sherical harmonics :
!                              ------------------------------------------------------
!                              s
                               if(i1.eq.1)orbital= envelop*sqrt(1.0/(4.0*rpi))

                               if(i1.eq.2.or.i1.eq.3)then
                                 cmp1=envelop*rsintheta*(rcosphy+cmplx(0.0,1.0)*rsinphy)*sqrt(3.0/(8.0*rpi))
                                 cmp2=envelop*(-1.0)*rsintheta*(rcosphy-cmplx(0.0,1.0)*rsinphy)*sqrt(3.0/(8.0*rpi))
!                                py
                                 if(i1.eq.2)orbital= (1.0/sqrt(2.0))*(cmp1-cmp2)
!                                pz

                                 if(i1.eq.3)orbital= (-cmplx(0.0,1.0)/sqrt(2.0))*(cmp1+ cmp2)
                               end if
!                              px
                               if(i1.eq.4)then
                                 cmp1=envelop*rcostheta*sqrt(3.0/(4.0*rpi))
                                 orbital=  cmp1
                               end if

                               if(i1.eq.5.or.i1.eq.6)then
                                 cmp1=envelop*(-1.0)*sqrt(15.0/(8.0*rpi))*rsintheta*rcostheta*(rcosphy+cmplx(0.0,1.0)*rsinphy)
                                 cmp2=envelop*sqrt(15.0/(8.0*rpi))*rsintheta*rcostheta*(rcosphy-cmplx(0.0,1.0)*rsinphy)
!                                dyz
                                 if(i1.eq.5)orbital= (cmplx(0.0,1.0)/sqrt(2.0))*(cmp1+cmp2)
!                                dxz
                                 if(i1.eq.6)orbital= (-1.0/sqrt(2.0))*(cmp1-cmp2)
                               end if

                               if(i1.eq.7.or.i1.eq.8)then
                                 cmp1=envelop*sqrt(15.0/(32.0*rpi))*rsintheta*rsintheta*(rrcosphy2+cmplx(0.0,1.0)*rrsinphy2)
                                 cmp2=envelop*sqrt(15.0/(32.0*rpi))*rsintheta*rsintheta*(rrcosphy2-cmplx(0.0,1.0)*rrsinphy2)
!                                dxy
                                 if(i1.eq.7)orbital= (-cmplx(0.0,1.0)/sqrt(2.0))*(cmp1 - cmp2)
!                                dx2-y2
                                 if(i1.eq.8)orbital= (1.0/sqrt(2.0))*(cmp1 + cmp2)
                               end if

!                              dz2
                               if(i1.eq.9)then
                                 cmp1=envelop*sqrt(5.0/(16.0*rpi))*(3.0*rcostheta*rcostheta-1.0)
                                 orbital= cmp1
                               end if
!                              --------------------------------------------------------------------------
                               corbital = corbital  +  orbital*&
&                               exp(cmplx(0.0,1.0)*(rkx*(rn1*rlx1+rn2*rlx2+rn3*rlx3 - x1)+&
&                               rky*(rn1*rly1+rn2*rly2+rn3*rly3 - y1)+&
&                               rkz*(rn1*rlz1+rn2*rlz2+rn3*rlz3 - z1) ))
!                              corb(i,j,k)= corbital
                               corb(i,j,k)= corb(i,j,k) + corbital
!                              rnorm = rnorm + conjg(corbital)*corbital

!                              *************************************************
                             end if
!                            *************************************************

                           end do  ! orbital type
                         end do ! atom
                         if(icountorb.ne.nbands)then
                           write(std_out,*)icountorb,' -vs- ',nbands
                           write(std_out,*)'total number of orbitals and nband do not match'
                           stop
                         end if
!                        -----------new
                       end if
!                      -------

!                      normalization :
                       if(icellx.eq.ncellx.and.icelly.eq.ncelly.and.icellz.eq.ncellz)then
                         rnorm = rnorm + conjg(corb(i,j,k))*corb(i,j,k)
                       end if

                     end do  ! real space grid
                   end do
                 end do

               end do
             end do  ! sum over neighbouring unitcell ends
           end do

!          normalization :
!          write(std_out,*) 'norm:',rnorm
           corb = corb / sqrt(rnorm)

!          #########################################

           do iband1 =  1, nbands
!            write(std_out,*) '   overlap beqween', iband1,' and ',iband2
!            reading wave function from "_wfk" file and calculating the overlap :

             call overlap_wf(corb,e_kpt,exchn2n3d,csppol,iband1,ckpt,&
&             ecut,headform,istwfk,kpt,nband,nbands,nkpt,npwarr,&
&             nr1,nr2,nr3,nspinor,nsppol,paral_kgb,rprimd,cmp2)

             coverlap(iband2,iband1) = cmp2

           end do  ! iband1
!          #########################################
         end do  ! iband2
!        #####################################################

!        singular value decomposition :
!        -----------------------------------
         write(29,*)" kpoint:",ikpt_x,ikpt_y,ikpt_z
         write(30,*)" kpoint:",ikpt_x,ikpt_y,ikpt_z
         write(29,*)"  "
         write(30,*)"  "

         do i = 1,nbands
           write(29,"(20(f12.5))")(abs(coverlap(i,j)),j=1,nbands)
         end do
         crhsvd = 'a'
         call zgesvd(crhsvd,crhsvd,nbands,nbands,coverlap,nbands,s,cu,nbands,&
         cvt,nbands,csvd, ndimcsvd, rwork, info )

         cm5 = 0.0
         do i = 1,nbands
           if(abs(s(i)).gt.0.0000001)cm5(i,i)= 1.0
         end do
         write(30,"(10f8.4)")(s(j),j=1,nbands)
         write(31,*)' singular values:'
         write(31,"(10f8.4)")(s(j),j=1,nbands)
         write(31,*)' note: if you find some of the singular values consitently'
         write(31,*)'       very close to zero (0.00xxx) please stop the'
         write(31,*)'       calculation and look at the overlap_abs.dat file.'
         write(31,*)'       exclude orbitals (from under the orbdetails variable '
         write(31,*)'       in the input file) corresponding to predominantly'
         write(31,*)'       zero rows of overlap matrix. '

         cm1=matmul(cu,cvt)
         cm2=conjg(cm1)

!        crotl(ikpt_x,ikpt_y,ikpt_z,1:nbands,1:nbands) = cm2(1:nbands,1:nbands)

         write(26)((cm2(i1,i2),i1=1,nbands),i2 = 1,nbands)

!        --------------------------------------------------------------------------
       end do  ! ikpt_x
     end do  ! ikpt_y  ! ikpt loop ends here
   end do  ! ikpt_z

 end do  ! csppol
 close(26)

 8765    continue  ! if jobtype is only WLO construction

!##########################################################################
!#                                                                        #
!#                   Wannier function construction                        #
!#                                                                        #
!##########################################################################
 normtot = 0.0
 do  iwan = 1,nwantot  ! BAND INDEX, ie, WANNIER INDEX STARTS HERE

   open(26,file="crot.dat",form='unformatted')
   rewind(26)

   iband1 = nwanindex(iwan)
   cwannier = 0.0
   rwanmax = -1000.0
!  ======================================================
   do csppol = 1,nsppol
     cspinor = 1 ! only scalar wave functions.

     write(std_out,*)' reading'
     do ikpt_z = 1,nkz
       do ikpt_y = 1,nky
         do ikpt_x = 1,nkx
           rewind(31)
           write(31,*)' Doing : WLO construction.'
           write(31,*)' WLO index:',iband1

           write(std_out,*)' Doing : WLO construction.'
           write(std_out,*)' WLO index:',iband1

           read(26)((cm2(i1,i2),i1=1,nbands),i2 = 1,nbands)

!          ---------------------------------------------
!          k index for ascending order of wave vector k
!          ----------------------------------------------
           call kptindex(shiftx,shifty,shiftz,ikpt_x,ikpt_y,ikpt_z,nkx,nky,nkz,i,j,k)
           ckpt=i + (j - 1)*nkx + (k - 1)*nkx*nky
           write(31, '(a,3i4)' )" kpoint index : ",ikpt_x,ikpt_y,ikpt_z
           write(31, '(a,3i4)' )" kpoint index : ",ikpt_x,ikpt_y,ikpt_z
           write(31, '(a,3i4)' )' folded to:',i,j,k

           write(std_out,'(a,3i4)')" kpoint index : ",ikpt_x,ikpt_y,ikpt_z
           rkpt1 = kpt(1,ckpt)
           rkpt2 = kpt(2,ckpt)
           rkpt3 = kpt(3,ckpt)

           write(31,*)'kpt:',rkpt1,rkpt2,rkpt3

           rkx = rkpt1*rkx1 + rkpt2*rkx2 + rkpt3*rkx3
           rky = rkpt1*rky1 + rkpt2*rky2 + rkpt3*rky3
           rkz = rkpt1*rkz1 + rkpt2*rkz2 + rkpt3*rkz3
!          ----------------------------------
!          rotating the energy eigenstates :
!          ----------------------------------
           cbfns1  = cmplx(0.0_dp,0.0_dp)
           do iband2 = 1,nbands

             call wfread(cbfns,e_kpt,exchn2n3d,csppol,iband2,ckpt,&
&             ecut,headform,istwfk,kpt,nband,nbands,nkpt,npwarr,&
&             nr1,nr2,nr3,nspinor,nsppol,paral_kgb,rprimd)

             if(insmet.eq.2)then
               roccup=(e_kpt(iband2)-efermi)/tsmear ! rkt
               roccup = 1.0_dp/(1.0_dp+exp(roccup))
             else
               roccup = 1.0
             end if
             do igz = 1,ngz
               do igy = 1,ngy
                 do igx = 1,ngx

                   cbfns1(igx,igy,igz)=cbfns1(igx,igy,igz)&
&                   +cm2(iband1,iband2)*cbfns(igx,igy,igz)*sqrt(roccup)

                 end do
               end do
             end do
           end do ! iband2

!          ------------------
!          wf construction:
!          ------------------
           if(rcellx+rcelly+rcellz .ne. 0.0) then  ! 0 is default

             munitx = (ngx*ncellx1+anint(mod(real(ngx*ncellx1),2.0)))/2 +&
&             anint(mod(real(ngx*ncellx1)+1.0,2.0))
             munity = (ngy*ncelly1+anint(mod(real(ngy*ncelly1),2.0)))/2 +&
&             anint(mod(real(ngy*ncelly1)+1.0,2.0))
             munitz = (ngz*ncellz1+anint(mod(real(ngz*ncellz1),2.0)))/2 +&
&             anint(mod(real(ngz*ncellz1)+1.0,2.0))

             re1 = rcoordall(1,iband1,csppol)
             re2 = rcoordall(2,iband1,csppol)
             re3 = rcoordall(3,iband1,csppol)
             munitx = munitx + anint(re1*real(ngx))
             munity = munity + anint(re2*real(ngy))
             munitz = munitz + anint(re3*real(ngz))

             nstartx = munitx  - anint(rcellx*0.5*real(ngx))
             nstopx  = munitx  + anint(rcellx*0.5*real(ngx))

             nstarty = munity  - anint(rcelly*0.5*real(ngy))
             nstopy  = munity  + anint(rcelly*0.5*real(ngy))

             nstartz = munitz  - anint(rcellz*0.5*real(ngz))
             nstopz  = munitz  + anint(rcellz*0.5*real(ngz))

             if(nstartx.le.1)nstartx = 1
             if(nstartx.gt.ngx*ncellx1)nstartx = ngx*ncellx1
             if(nstarty.le.1)nstarty = 1
             if(nstarty.gt.ngy*ncelly1)nstarty = ngy*ncelly1
             if(nstartz.le.1)nstartz = 1
             if(nstartz.gt.ngz*ncellz1)nstartz = ngz*ncellz1

             if(nstopx.le.1)nstopx = 1
             if(nstopx.gt.ngx*ncellx1)nstopx = ngx*ncellx1
             if(nstopy.le.1)nstopy = 1
             if(nstopy.gt.ngy*ncelly1)nstopy = ngy*ncelly1
             if(nstopz.le.1)nstopz = 1
             if(nstopz.gt.ngz*ncellz1)nstopz = ngz*ncellz1

           else

             nstartx = 1
             nstopx  = ngx*ncellx1

             nstarty = 1
             nstopy  = ngy*ncelly1

             nstartz = 1
             nstopz  = ngz*ncellz1

           end if

           do icellz = 1,ncellz1
             do icelly = 1,ncelly1
               do icellx = 1,ncellx1


                 do igz = 1,ngz
                   do igy = 1,ngy
                     do igx = 1,ngx

                       l = ngx*(icellx-1) + igx
                       m = ngy*(icelly-1) + igy
                       n = ngz*(icellz-1) + igz

                       if(l.ge.nstartx.and.l.le.nstopx.and.&
&                       m.ge.nstarty.and.m.le.nstopy.and.&
&                       n.ge.nstartz.and.n.le.nstopz)then

                         if(mod(real(ncellx1),2.0).eq.1)then
                           ii1 = igx
                         else
                           if(mod(real(ngx),2.0).eq.1)then
                             ii1 = igx+(ngx-1)/2
                             if(ii1.gt.ngx)ii1 = ii1 - ngx
                           else
                             ii1 = igx+ngx/2
                             if(ii1.gt.ngx)ii1 = ii1 - ngx
                           end if
                         end if

                         if(mod(real(ncelly1),2.0).eq.1)then
                           ii2 = igy
                         else
                           if(mod(real(ngy),2.0).eq.1)then
                             ii2 = igy+(ngy-1)/2
                             if(ii2.gt.ngy)ii2 = ii2 - ngy
                           else
                             ii2 = igy+ngy/2
                             if(ii2.gt.ngy)ii2 = ii2 - ngy
                           end if
                         end if

                         if(mod(real(ncellz1),2.0).eq.1)then
                           ii3 = igz
                         else
                           if(mod(real(ngz),2.0).eq.1)then
                             ii3 = igz+(ngz-1)/2
                             if(ii3.gt.ngz)ii3 = ii3 - ngz
                           else
                             ii3 = igz+ngz/2
                             if(ii3.gt.ngz)ii3 = ii3 - ngz
                           end if
                         end if

                         r1 = (real(l)-(1.0-mod(real(ngx*ncellx1),2.0)*0.5) )/real(ngx)
                         r2 = (real(m)-(1.0-mod(real(ngy*ncelly1),2.0)*0.5) )/real(ngy)
                         r3 = (real(n)-(1.0-mod(real(ngz*ncellz1),2.0)*0.5) )/real(ngz)

                         x = rlx_large*(-0.5)+rlx1*r1+rlx2*r2+rlx3*r3
                         y = rly_large*(-0.5)+rly1*r1+rly2*r2+rly3*r3
                         z = rlz_large*(-0.5)+rlz1*r1+rlz2*r2+rlz3*r3

                         cwannier(l,m,n)=cwannier(l,m,n)+&
&                         cbfns1(ii1,ii2,ii3)*exp(cmplx(0.0,1.0)*(rkx*x+rky*y+rkz*z))

                         if(ikpt_x.eq.nkx.and.ikpt_y.eq.nky.and.ikpt_z.eq.nkz)then
                           if(abs(cwannier(l,m,n)).gt.rwanmax)then
                             rwanmax = abs(cwannier(l,m,n))
                             nwanmax1 = l
                             nwanmax2 = m
                             nwanmax3 = n
                           end if
                         end if
!                        ---------------new
                       end if
!                      ------------
                     end do
                   end do
                 end do
               end do
             end do
           end do

         end do  ! ikpt_x
       end do  ! ikpt_y  ! ikpt loop ends here
     end do  ! ikpt_z

!    --------------------------------------------------------------------------
!    processing of wlo :
!    -------------------
!    normalization  :
     cwannier=cwannier/real(nkx*nky*nkz)

     re1 = 0.0_dp
     re2 = 0.0_dp
     re3 = 0.0_dp
     re4 = 0.0_dp
     re5 = 0.0_dp
     re6 = 0.0_dp
     re7 = 0.0_dp

     do k = nstartz,nstopz
       do j = nstarty,nstopy
         do i = nstartx,nstopx

           r1 = (real(i)-(1.0-mod(real(ngx*ncellx1),2.0)*0.5) )/real(ngx)
           r2 = (real(j)-(1.0-mod(real(ngy*ncelly1),2.0)*0.5) )/real(ngy)
           r3 = (real(k)-(1.0-mod(real(ngz*ncellz1),2.0)*0.5) )/real(ngz)

           x = rlx_large*(-0.5)+rlx1*r1+rlx2*r2+rlx3*r3
           y = rly_large*(-0.5)+rly1*r1+rly2*r2+rly3*r3
           z = rlz_large*(-0.5)+rlz1*r1+rlz2*r2+rlz3*r3

           csum=cwannier(i,j,k)*conjg(cwannier(i,j,k))
           re1 = re1 + x*csum
           re2 = re2 + y*csum
           re3 = re3 + z*csum
           re4 = re4 + x*x*csum
           re5 = re5 + y*y*csum
           re6 = re6 + z*z*csum
           re7 = re7 + csum

         end do
       end do
     end do
     if(re7.gt.0.0000001)then
       re1 = re1/re7  ; rxcm = re1/re7
       re2 = re2/re7  ; rycm = re2/re7
       re3 = re3/re7  ; rzcm = re3/re7
       re4 = re4/re7
       re5 = re5/re7
       re6 = re6/re7
     end if
     normtot = normtot + re7
     write(std_out,*)" ========================"
     write(std_out,*)" WLO : ",iband1,'csppol:',csppol

     write(std_out,*)' normalization:  ',re7
     write(std_out,*)' qo center x',re1
     write(std_out,*)' qo center y',re2
     write(std_out,*)' qo center z',re3
     write(std_out,*)"-----------------------------"
     write(std_out,*)' loclength_x*x',(re4-re1*re1)
     write(std_out,*)' loclength_y*y',(re5-re2*re2)
     write(std_out,*)' loclength_z*z',(re6-re3*re3)
     write(6662,"(i4,4(1x,f8.5),i2)")iband1,re4-re1*re1,re5-re2*re2,re6-re3*re3,re7,csppol
     write(6661,"(i4,4(1x,f8.5),i2)")iband1,re1,re2,re3,re7,csppol
     if(iwan.eq.nwantot)write(6661,*)' total norm:',normtot
     if(iwan.eq.nwantot)write(6662,*)' total norm:',normtot

!    ---------------------------------
!    Plotting :
     munitx = (ngx*ncellx1+anint(mod(real(ngx*ncellx1),2.0)))/2 +&
&     anint(mod(real(ngx*ncellx1)+1.0,2.0))
     munity = (ngy*ncelly1+anint(mod(real(ngy*ncelly1),2.0)))/2 +&
&     anint(mod(real(ngy*ncelly1)+1.0,2.0))
     munitz = (ngz*ncellz1+anint(mod(real(ngz*ncellz1),2.0)))/2 +&
&     anint(mod(real(ngz*ncellz1)+1.0,2.0))

     write(std_out,'(a,3i5)' ) ' mid points :',munitx,munity,munitz
     re1 = rcoordall(1,iband1,csppol)
     re2 = rcoordall(2,iband1,csppol)
     re3 = rcoordall(3,iband1,csppol)
     munitx = munitx + anint(re1*real(ngx))
     munity = munity + anint(re2*real(ngy))
     munitz = munitz + anint(re3*real(ngz))

     write(std_out,'(a,3i5)' ) ' center sites :',munitx,munity,munitz

     if(icellplotx.gt.1)then
       nstartx1 = munitx - ngx + 1
       nstopx1  = icellplotx*ngx + nstartx1 - 1
       if(nstartx1.le.1)then
         nstartx1 = 1
         nstopx1  = icellplotx*ngx
       end if
       if(nstopx1.ge.ncellx1*ngx)then
         nstartx1 = (ncellx1 -2)*ngx + 1
         nstopx1  = ncellx1*ngx
       end if
     else
       nstartx1 = 1
       nstopx1  = icellplotx*ngx
     end if

     if(icellploty.gt.1)then
       nstarty1 = munity - ngy + 1
       nstopy1  = icellploty*ngy + nstarty1 - 1
       if(nstarty1.le.1)then
         nstarty1 = 1
         nstopy1  = icellploty*ngy
       end if
       if(nstopy1.ge.ncelly1*ngy)then
         nstarty1 = (ncelly1 -2)*ngy + 1
         nstopy1  = ncelly1*ngy
       end if
     else
       nstarty1 = 1
       nstopy1 = icellploty*ngy
     end if

     if(icellplotz.gt.1)then
       nstartz1 = munitz - ngz + 1
       nstopz1  = icellplotz*ngz + nstartz1 - 1
       if(nstartz1.le.1)then
         nstartz1 = 1
         nstopz1  = icellplotz*ngz
       end if
       if(nstopz1.ge.ncellz1*ngz)then
         nstartz1 = (ncellz1 -2)*ngz + 1
         nstopz1  = ncellz1*ngz
       end if
     else
       nstartz1 = 1
       nstopz1  = icellplotz*ngz
     end if

     nstartx = max(nstartx,nstartx1)
     nstarty = max(nstarty,nstarty1)
     nstartz = max(nstartz,nstartz1)

     nstopx = min(nstopx1,nstopx)
     nstopy = min(nstopy1,nstopy)
     nstopz = min(nstopz1,nstopz)

     write(std_out,'(a)' )' plotting range of LO :'
     write(std_out,'(a,2i5)' )' plot x :',nstartx,nstopx
     write(std_out,'(a,2i5)' )' plot y :',nstarty,nstopy
     write(std_out,'(a,2i5)' )' plot z :',nstartz,nstopz

     write(std_out,'(a,3i5)' )' ng :',ngx,ngy,ngz

     write(std_out,'(a,3i5)' )'     ',shape(cwannier)

     rewind(10000+iband1*10+csppol)
     if(iwrite.eq.1)rewind(20000+iband1*10+csppol)

!    plotting in xcrysden format :

     write(10000+iband1*10+csppol,*)"CRYSTAL"
     write(10000+iband1*10+csppol,*)"PRIMVEC"
     write(10000+iband1*10+csppol,"(3f10.5)") rlx1  , rly1 ,rlz1
     write(10000+iband1*10+csppol,"(3f10.5)") rlx2,  rly2  , rlz2
     write(10000+iband1*10+csppol,"(3f10.5)") rlx3, rly3,  rlz3
     write(10000+iband1*10+csppol,*)"CONVVEC"
     write(10000+iband1*10+csppol,"(3f10.5)")  rlx1  , rly1 ,rlz1
     write(10000+iband1*10+csppol,"(3f10.5)") rlx2,  rly2  , rlz2
     write(10000+iband1*10+csppol,"(3f10.5)") rlx3, rly3,  rlz3
     write(10000+iband1*10+csppol,*)"PRIMCOORD"
     write(10000+iband1*10+csppol,*)natom, 1

     do i = 1,natom
       write(10000+iband1*10+csppol,"(i3,3f10.5)")int(znucl(typat(i))),(xcart(j,i),j=1,3)
     end do
     write(10000+iband1*10+csppol,*)"BEGIN_BLOCK_DATAGRID_3D"
     write(10000+iband1*10+csppol,*)"WANNIER_CENTER_3D"
     write(10000+iband1*10+csppol,*)"DATAGRID_3D_WANNIER"

     write(10000+iband1*10+csppol,*)nstopx-nstartx+1,nstopy-nstarty+1,nstopz-nstartz+1

     re4 = real(nstopx-nstartx+1-1)/real(ngx)
     re5 = real(nstopy-nstarty+1-1)/real(ngy)
     re6 = real(nstopz-nstartz+1-1)/real(ngz)

     r1 = (real(nstartx)-(1.0-mod(real(ngx*ncellx1),2.0)*0.5) )/real(ngx)
     r2 = (real(nstarty)-(1.0-mod(real(ngy*ncelly1),2.0)*0.5) )/real(ngy)
     r3 = (real(nstartz)-(1.0-mod(real(ngz*ncellz1),2.0)*0.5) )/real(ngz)

     x = rlx_large*(-0.5)+rlx1*r1+rlx2*r2+rlx3*r3
     y = rly_large*(-0.5)+rly1*r1+rly2*r2+rly3*r3
     z = rlz_large*(-0.5)+rlz1*r1+rlz2*r2+rlz3*r3

     write(10000+iband1*10+csppol,"(3f10.5)")x,y,z
     if(iwrite.eq.1)write(20000+iband1*10+csppol,*)x,y,z

     write(10000+iband1*10+csppol,"(3f10.5)")rlx1*re4,rly1*re4,rlz1*re4
     write(10000+iband1*10+csppol,"(3f10.5)")rlx2*re5,rly2*re5,rlz2*re5
     write(10000+iband1*10+csppol,"(3f10.5)")rlx3*re6,rly3*re6,rlz3*re6

     if(iwrite.eq.1)then
       write(20000+iband1*10+csppol,*)nstartx,nstarty,nstartz
       write(20000+iband1*10+csppol,*)nstopx,nstopy,nstopz
       write(20000+iband1*10+csppol,"(3f10.5)")rlx1*re4,rly1*re4,rlz1*re4
       write(20000+iband1*10+csppol,"(3f10.5)")rlx2*re5,rly2*re5,rlz2*re5
       write(20000+iband1*10+csppol,"(3f10.5)")rlx3*re6,rly3*re6,rlz3*re6
     end if

!    Apparently, the determination of nwanmax1, nwanmax2 and nwanmax3 is no portable accross platforms'
     write(std_out,'(a,3i5)' )'-wanmax:',nwanmax1,nwanmax2,nwanmax3
     cmp1 = cwannier(nwanmax1,nwanmax2,nwanmax3)
     write(std_out,*) cmp1
     cmp1 = abs(cmp1)/cmp1

     re1 = 0.0
     do igz1 = nstartz  ,  nstopz
       do igy1 = nstarty  ,  nstopy
         do igx1 = nstartx  ,  nstopx
           write(10000+iband1*10+csppol,"(1f10.5)")real(cwannier(igx1,igy1,igz1)*cmp1)
           if(iwrite.eq.1)write(20000+iband1*10+csppol,*)cwannier(igx1,igy1,igz1)
           re1 = re1 + abs(cwannier(igx1,igy1,igz1))**2.0
         end do
       end do
     end do
     write(std_out,*)'norm enclosed:',re1

     write(10000+iband1*10+csppol,*)"END_DATAGRID_3D"
     write(10000+iband1*10+csppol,*)"END_BLOCK_DATAGRID_3D"

!    ------------------------------------------------------
!    X Y Z projection calculation
!    ------------------------------------------------------
!    goto 6776
     re1 = abs(rlx1)*real(ncellx1)+abs(rlx2)*real(ncelly1)+abs(rlx3)*real(ncellz1)
     re2 = abs(rly1)*real(ncellx1)+abs(rly2)*real(ncelly1)+abs(rly3)*real(ncellz1)
     re3 = abs(rlz1)*real(ncellx1)+abs(rlz2)*real(ncelly1)+abs(rlz3)*real(ncellz1)

     re4 = max(abs(rlx1)/real(ngx),abs(rlx2)/real(ngy),abs(rlx3)/real(ngz))
     re5 = max(abs(rly1)/real(ngx),abs(rly2)/real(ngy),abs(rly3)/real(ngz))
     re6 = max(abs(rlz1)/real(ngx),abs(rlz2)/real(ngy),abs(rlz3)/real(ngz))

     nxtot = anint(re1/re4)
     nytot = anint(re2/re5)
     nztot = anint(re3/re6)
     ABI_ALLOCATE(rprojx,(nxtot))
     ABI_ALLOCATE(rprojy,(nytot))
     ABI_ALLOCATE(rprojz,(nztot))
     rprojx = 0.0
     rprojy = 0.0
     rprojz = 0.0

     re4 = min(rlx1*real(ncellx1),rlx2*real(ncelly1),rlx3*real(ncellz1))
     if(re4.ge.0.0)re4 = 0.0
     re5 = min(rly1*real(ncellx1),rly2*real(ncelly1),rly3*real(ncellz1))
     if(re5.ge.0.0)re5 = 0.0
     re6 = min(rlz1*real(ncellx1),rlz2*real(ncelly1),rlz3*real(ncellz1))
     if(re6.ge.0.0)re6 = 0.0

     do k = nstartz,nstopz
       do j = nstarty,nstopy
         do i = nstartx,nstopx

           r1 = (real(i)-(1.0-mod(real(ngx*ncellx1),2.0)*0.5) )/real(ngx)
           r2 = (real(j)-(1.0-mod(real(ngy*ncelly1),2.0)*0.5) )/real(ngy)
           r3 = (real(k)-(1.0-mod(real(ngz*ncellz1),2.0)*0.5) )/real(ngz)

           x = rlx1*r1+rlx2*r2+rlx3*r3
           y = rly1*r1+rly2*r2+rly3*r3
           z = rlz1*r1+rlz2*r2+rlz3*r3

           i1 = 1+anint((x - re4 )/(abs(re1)/real(nxtot)))
           i2 = 1+anint((y - re5 )/(abs(re2)/real(nytot)))
           i3 = 1+anint((z - re6 )/(abs(re3)/real(nztot)))

           if(i1.gt.nxtot)write(std_out,*) 'i1=',i1
           if(i2.gt.nytot)write(std_out,*) 'i2=',i2
           if(i3.gt.nztot)write(std_out,*) 'i3=',i3

           rprojx(i1)= rprojx(i1) + cwannier(i,j,k)*conjg(cwannier(i,j,k))
           rprojy(i2)= rprojy(i2) + cwannier(i,j,k)*conjg(cwannier(i,j,k))
           rprojz(i3)= rprojz(i3) + cwannier(i,j,k)*conjg(cwannier(i,j,k))

         end do
       end do
     end do

     rewind(510000+iband1*10+csppol)
     rewind(520000+iband1*10+csppol)
     rewind(530000+iband1*10+csppol)

     do i = 1,nxtot
       x =  -0.5*re1 +(abs(re1)/real(nxtot))*(real(i)-(1.0-mod(real(nxtot),2.0)*0.5))
       write(510000+iband1*10+csppol,"(10f12.6)")x,rprojx(i)
     end do
     do i = 1,nytot
       y =  -0.5*re2+(abs(re2)/real(nytot))*(real(i)-(1.0-mod(real(nytot),2.0)*0.5))
       write(520000+iband1*10+csppol,"(10f12.6)")y,rprojy(i)
     end do
     do i = 1,nztot
       z =  -0.5*re3+(abs(re3)/real(nztot))*(real(i)-(1.0-mod(real(nztot),2.0)*0.5))
       write(530000+iband1*10+csppol,"(10f12.6)")z,rprojz(i)
     end do

     close(510000+iband1*10+csppol)
     close(520000+iband1*10+csppol)
     close(530000+iband1*10+csppol)
     ABI_DEALLOCATE(rprojx)
     ABI_DEALLOCATE(rprojy)
     ABI_DEALLOCATE(rprojz)

!    ------------------------------------------------------------
   end do  ! csppol
   close(26)
!  -------------------------------------------------------------
 end do ! iwan

 close(29)
 close(30)
 close(6661)
 close(6662)
 close(unit=31,status='delete')

 close(19)
 return
end subroutine localorb_s
!!***
