!{\src2tex{textfont=tt}}
!!****f* ABINIT/nres2vres
!!
!! NAME
!! nres2vres
!!
!! FUNCTION
!! Convert a density residual into a potential residual
!! using a first order formula:
!!     V^res(r)=dV/dn.n^res(r)
!!             =V_hartree(n^res)(r) + Kxc.n^res(r)
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | natom= number of atoms in cell
!!  | nspden=number of spin-density components
!!  | ntypat=number of atom types
!!  | typat(natom)=type (integer) for each atom
!! gsqcut=cutoff value on G**2 for sphere inside fft box
!! izero=if 1, unbalanced components of Vhartree(g) are set to zero
!! kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!! mpi_enreg=informations about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT
!! nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!! nkxc=second dimension of the array kxc, see rhohxc.F90 for a description
!! nresid(nfft,nspden)= the input density residual
!! n3xccc=dimension of the xccc3d array (0 or nfft).
!! optnc=option for non-collinear magnetism (nspden=4):
!!       1: the whole 2x2 Vres matrix is computed
!!       2: only Vres^{11} and Vres^{22} are computed
!! optxc=0 if LDA part of XC kernel has only to be taken into account (even for GGA)
!!       1 if XC kernel has to be fully taken into
!!      -1 if XC kernel does not have to be taken into account
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! rhor(nfft,nspden)=electron density in real space
!!                   (used only if Kxc was not computed before)
!! rprimd(3,3)=dimensional primitive translation vectors (bohr)
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!! xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!! vresid(nfft,nspden)= the output potential residual
!!
!! PARENTS
!!      etotfor,forstr
!!
!! CHILDREN
!!      fourdp,hartre,leave_new,metric,mkvxc3,pawmknhat,psolver_hartree,rhohxc
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nres2vres(dtset,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,nhat,&
&                 nkxc,nresid,n3xccc,optnc,optxc,pawang,pawfgrtab,pawrhoij,pawtab,&
&                 rhor,rprimd,usepaw,usexcnhat,vresid,wvl,xccc3d,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nres2vres'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_62_poisson
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,n3xccc,nfft,nkxc,optnc,optxc,usepaw,usexcnhat
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: kxc(nfft,nkxc),nresid(nfft,dtset%nspden)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3),xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(out) :: vresid(nfft,dtset%nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%ntypat*usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,ifft,ider,idir,ipert,ispden,nhatgrdim,nkxc_cur,option
 real(dp) :: dum,dvdn,dvdz,energy,fact,m_dot_mres,m_norm_min,ucvol,vxcavg
 character(len=500) :: message
!arrays
 integer :: nk3xc
 real(dp) :: dummy6(6),gmet(3,3),gprimd(3,3),qq(3),rmet(3,3)
 real(dp),allocatable :: dummy(:),kxc_cur(:,:),m_norm(:),nhatgr(:,:,:)
 real(dp),allocatable :: nresg(:,:),nresid_diag(:,:),rhor0(:,:),vhres(:)
 real(dp),allocatable :: vresid_diag(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' nres2vres : enter ',optxc
!ENDDEBUG

!Compatibility tests:
 if(optxc<-1.or.optxc>1)then
   write(message, '(4a)') ch10,&
&   ' nres2vres : BUG -',ch10,&
   '   Wrong value for optxc !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if
 if((optnc/=1.and.optnc/=2).or.(dtset%nspden/=4.and.optnc/=1))then
   write(message, '(4a)') ch10,&
&   ' nres2vres : BUG -',ch10,&
   '   Wrong value for optnc !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if
 if(dtset%icoulomb==1.and.optxc/=-1)then
   write(message, '(4a)') ch10,&
&   ' nres2vres : ERROR -',ch10,&
   '   This routine is not compatible with icoulomb==1 and optxc/=-1 !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if
 if(dtset%nspden==4.and.dtset%xclevel==2.and.optxc==1.and.nkxc/=23)then
   write(message, '(4a)') ch10,&
&   ' nres2vres : ERROR -',ch10,&
   '   Wrong values for optxc and nkxc !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 qq=zero
 nkxc_cur=0
 m_norm_min=EPSILON(0.0_dp)**2
 if (dtset%xclevel==1.or.optxc==0) nkxc_cur=3-2*mod(dtset%nspden,2)
 if (dtset%xclevel==2.and.optxc==1) nkxc_cur=23
 ABI_ALLOCATE(vhres,(nfft))

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Compute density residual in reciprocal space
 if (dtset%icoulomb==0) then
   ABI_ALLOCATE(nresg,(2,nfft))
   ABI_ALLOCATE(dummy,(nfft))
   dummy(:)=nresid(:,1)
   call fourdp(1,nresg,dummy,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   ABI_DEALLOCATE(dummy)
 end if

!First case: Kxc has already been computed
!-----------------------------------------
 if (nkxc==nkxc_cur.or.optxc==-1) then

!  Compute VH(n^res)(r)
   if (dtset%icoulomb == 0) then
     call hartre(1,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qq,nresg,vhres)
   else
     call PSolver_hartree(dtset,energy,mpi_enreg,nresid(:,1),rprimd,vhres,wvl)
   end if

!  Compute Kxc(r).n^res(r)
   if (optxc/=-1) then

!    Collinear magnetism or non-polarized
     if (dtset%nspden/=4) then
       call mkvxc3(1,kxc,mpi_enreg,nfft,ngfft,nkxc,&
&       dtset%nspden,0,2,dtset%paral_kgb,qq,nresid,rprimd,vresid,dummy)

!      Non-collinear magnetism
!      Has to locally "rotate" n^res(r) (according to magnetization),
!      compute V^res(r) and rotate it back
     else
       ABI_ALLOCATE(nresid_diag,(nfft,2))
       ABI_ALLOCATE(vresid_diag,(nfft,2))
       ABI_ALLOCATE(rhor0,(nfft,dtset%nspden))
       ABI_ALLOCATE(m_norm,(nfft))
!      -- Compute "initial" density
       rhor0(:,:)=rhor(:,:)-nresid(:,:)
!      -- Rotate n^res(r)
       do ifft=1,nfft
         nresid_diag(ifft,1)=nresid(ifft,1)
         m_norm(ifft)=sqrt(rhor0(ifft,2)**2+rhor0(ifft,3)**2+rhor0(ifft,4)**2)
         m_dot_mres=rhor0(ifft,2)*nresid(ifft,2)+rhor0(ifft,3)*nresid(ifft,3) &
&         +rhor0(ifft,4)*nresid(ifft,4)
         if(m_norm(ifft)>m_norm_min)then
           nresid_diag(ifft,2)=half*(nresid_diag(ifft,1)+m_dot_mres/m_norm(ifft))
         else
           nresid_diag(ifft,2)=nresid_diag(ifft,1)
         end if
       end do
!      -- Compute Kxc(r).n^res(r)_rotated
       call mkvxc3(1,kxc,mpi_enreg,nfft,ngfft,nkxc,&
&       2,0,2,dtset%paral_kgb,qq,nresid_diag,rprimd,vresid_diag,dummy)
       ABI_DEALLOCATE(nresid_diag)
!      -- Rotate back V^res(r)
       if (optnc==1) then
         do ifft=1,nfft
           dvdn=(vresid_diag(ifft,1)+vresid_diag(ifft,2))*half
           dvdz=(vresid_diag(ifft,1)-vresid_diag(ifft,2))*half
           if(m_norm(ifft)>m_norm_min)then
             fact=dvdz/m_norm(ifft)
             dum=rhor0(ifft,4)*fact
             vresid(ifft,1)=dvdn+dum
             vresid(ifft,2)=dvdn-dum
             vresid(ifft,3)= rhor0(ifft,2)*fact
             vresid(ifft,4)=-rhor0(ifft,3)*fact
           else
             vresid(ifft,1:2)=dvdn
             vresid(ifft,3:4)=zero
           end if
         end do
       else
         do ifft=1,nfft
           dvdn=(vresid_diag(ifft,1)+vresid_diag(ifft,2))*half
           dvdz=(vresid_diag(ifft,1)-vresid_diag(ifft,2))*half
           if(m_norm(ifft)>m_norm_min)then
             dum=dvdz*rhor0(ifft,4)/m_norm(ifft)
             vresid(ifft,1)=dvdn+dum
             vresid(ifft,2)=dvdn-dum
           else
             vresid(ifft,1:2)=dvdn
           end if
         end do
       end if
       ABI_DEALLOCATE(vresid_diag)
       ABI_DEALLOCATE(rhor0)
       ABI_DEALLOCATE(m_norm)
     end if

   else
     vresid=zero
   end if

 end if

!2nd case: Kxc has to be computed
!--------------------------------
 if (nkxc/=nkxc_cur.and.optxc/=-1) then

!  For GGA, has to recompute gradients of nhat
   nhatgrdim=0
   if (usepaw==1.and.usexcnhat>0.and.dtset%xclevel==2.and.dtset%pawnhatxc>0) then
     nhatgrdim=1
     ABI_ALLOCATE(nhatgr,(nfft,dtset%nspden,3))
     ider=1;cplex=1;ipert=0;idir=0
     call pawmknhat(dum,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,dtset%natom,dtset%natom,&
&     nfft,ngfft,nhatgrdim,dtset%nspden,dtset%ntypat,dtset%paral_kgb,pawang,pawfgrtab,&
&     nhatgr,nhat,pawrhoij,pawrhoij,pawtab,qq,rprimd,ucvol,xred)
   end if

!  Has to use the "initial" density to compute Kxc
   ABI_ALLOCATE(rhor0,(nfft,dtset%nspden))
   rhor0(:,:)=rhor(:,:)-nresid(:,:)

!  Compute VH(n^res) and XC kernel (Kxc) together
   ABI_ALLOCATE(kxc_cur,(nfft,nkxc_cur))
   option=2;if (dtset%xclevel==2.and.optxc==0) option=12

!  to be adjusted for the call rhohxc
   nk3xc=1
   call rhohxc(dtset,energy,gsqcut,izero,kxc_cur,mpi_enreg,nfft,ngfft,&
&   nhat,usepaw,nhatgr,nhatgrdim,nkxc_cur,nk3xc,dtset%nspden,n3xccc,option,nresg,&
&   rhor0,rprimd,dummy6,usexcnhat,vhres,vresid,vxcavg,xccc3d)  !vresid=work space
   if (dtset%nspden/=4)  then
     ABI_DEALLOCATE(rhor0)
   end if
   if (nhatgrdim>0)  then
     ABI_DEALLOCATE(nhatgr)
   end if

!  Compute Kxc(r).n^res(r)

!  Collinear magnetism or non-polarized
   if (dtset%nspden/=4) then

     call mkvxc3(1,kxc_cur,mpi_enreg,nfft,ngfft,nkxc_cur,&
&     dtset%nspden,0,2,dtset%paral_kgb,qq,nresid,rprimd,vresid,dummy)

!    Non-collinear magnetism
!    Has to locally "rotate" n^res(r) (accroding to magnetization),
!    compute V^res(r) and rotate it back
   else
     ABI_ALLOCATE(nresid_diag,(nfft,2))
     ABI_ALLOCATE(vresid_diag,(nfft,2))
     ABI_ALLOCATE(m_norm,(nfft))
!    -- Rotate n^res(r)
     do ifft=1,nfft
       nresid_diag(ifft,1)=nresid(ifft,1)
       m_norm(ifft)=sqrt(rhor0(ifft,2)**2+rhor0(ifft,3)**2+rhor0(ifft,4)**2)
       m_dot_mres=rhor0(ifft,2)*nresid(ifft,2)+rhor0(ifft,3)*nresid(ifft,3) &
&       +rhor0(ifft,4)*nresid(ifft,4)
       nresid_diag(ifft,2)=half*(nresid_diag(ifft,1)+m_dot_mres/m_norm(ifft))
     end do
!    -- Compute Kxc(r).n^res(r)_rotated
     call mkvxc3(1,kxc_cur,mpi_enreg,nfft,ngfft,nkxc_cur,&
&     2,0,2,dtset%paral_kgb,qq,nresid_diag,rprimd,vresid_diag,dummy)
     ABI_DEALLOCATE(nresid_diag)
!    -- Rotate back V^res(r)
     if (optnc==1) then
       do ifft=1,nfft
         dvdn=(vresid_diag(ifft,1)+vresid_diag(ifft,2))*half
         dvdz=(vresid_diag(ifft,1)-vresid_diag(ifft,2))*half
         if(m_norm(ifft)>m_norm_min)then
           fact=dvdz/m_norm(ifft)
           dum=rhor0(ifft,4)*fact
           vresid(ifft,1)=dvdn+dum
           vresid(ifft,2)=dvdn-dum
           vresid(ifft,3)= rhor0(ifft,2)*fact
           vresid(ifft,4)=-rhor0(ifft,3)*fact
         else
           vresid(ifft,1:2)=dvdn
           vresid(ifft,3:4)=zero
         end if
       end do
     else
       do ifft=1,nfft
         dvdn=(vresid_diag(ifft,1)+vresid_diag(ifft,2))*half
         dvdz=(vresid_diag(ifft,1)-vresid_diag(ifft,2))*half
         if(m_norm(ifft)>m_norm_min)then
           dum=dvdz*rhor0(ifft,4)/m_norm(ifft)
           vresid(ifft,1)=dvdn+dum
           vresid(ifft,2)=dvdn-dum
         else
           vresid(ifft,1:2)=dvdn
         end if
       end do
     end if
     ABI_DEALLOCATE(vresid_diag)
     ABI_DEALLOCATE(m_norm)
     ABI_DEALLOCATE(rhor0)
   end if

   ABI_DEALLOCATE(kxc_cur)
 end if

!Assemble potential residual: V^res(r)=VH(n^res)(r) + Kxc(r).n^res(r)
!--------------------------------------------------------------------
 do ispden=1,dtset%nspden/optnc
   vresid(:,ispden)=vresid(:,ispden)+vhres(:)
 end do

 if (dtset%icoulomb==0)  then
   ABI_DEALLOCATE(nresg)
 end if
 ABI_DEALLOCATE(vhres)

!DEBUG
!write(std_out,*)' nres2vres : exit '
!ENDDEBUG

end subroutine nres2vres

!!***
