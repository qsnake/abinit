
!----------------------------------------------------------------------

!!****f* ABINIT/irred_perts
!! NAME
!!  irred_perts
!!
!! FUNCTION
!!  This routine finds the basic perturbations for a given q-points and 
!!  calculated the number of k-points in the extended Brillouin zone. 
!!  Results are written in out_unt using a XML-like format that can
!!  be easily parsed with an external script.
!!
!! INPUTS
!!  Dtset<dataset_type>=All input variables for this dataset.
!!  out_unt=Unit number for output.
!!
!! OUTPUT
!!  Output is written on unit out_unt.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      insy3,mati3inv,metric,mkrdim,symatm,symkpt,symq3,syper3,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine irred_perts(Dtset,out_unt)

 use m_profiling

 use defs_basis
 use m_errors

 use defs_abitypes, only : dataset_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'irred_perts'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_56_recipspace
 use interfaces_72_response, except_this_one => irred_perts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: out_unt
 type(dataset_type),intent(in) :: Dtset
!arrays

!Local variables -------------------------
!scalars
 integer,parameter :: rfmeth2=2,syuse0=0
 integer :: idx,idir,ipert,ikpt,nkpt_rbz,nsym,nsym1
 integer :: natom,mpert,sym,timrev
 integer :: rfddk,rfelfd,rfphon,rfstrs,rfuser,nkpt
 integer :: timrev_pert,nirred
 real(dp),parameter :: tolsym8=tol8
 real(dp) :: ucvol
 logical :: qeq0
 character(len=500) :: msg
!arrays
 integer :: rfdir(3)
 integer,allocatable :: symaf1(:),symaf1_tmp(:),symrc1(:,:,:),symrl1(:,:,:),symrl1_tmp(:,:,:)
 integer,allocatable :: pertsy(:,:),rfpert(:),symq(:,:,:)
 integer,allocatable :: indsym(:,:,:),symrec(:,:,:)
 integer,allocatable :: indkpt1(:),indkpt1_tmp(:),nkrbz(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: tnons1(:,:),tnons1_tmp(:,:)
 real(dp),allocatable :: wtk_folded(:),wtk_rbz(:)
 real(dp),allocatable :: kpt_rbz(:,:)

! *********************************************************************

!Obtain dimensional translations in reciprocal space gprimd, metrics and unit cell volume, from rprimd.
 call mkrdim(Dtset%acell_orig(:,1),Dtset%rprim_orig(:,:,1),rprimd)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
!
!Define the set of admitted perturbations.
 nsym  = Dtset%nsym
 natom = Dtset%natom
 mpert = natom+6
!
!Respfn input variables
 rfdir(:) = Dtset%rfdir(1:3)
 rfddk    = Dtset%rfddk 
 rfelfd   = Dtset%rfelfd
 rfphon   = Dtset%rfphon
 rfstrs   = Dtset%rfstrs
 rfuser   = Dtset%rfuser
!
!Initialize the list of perturbations rfpert.
 ABI_ALLOCATE(rfpert,(mpert))
 rfpert=0

 if (rfphon==1) rfpert(Dtset%rfatpol(1):Dtset%rfatpol(2))=1
 if (rfddk ==1) rfpert(natom+1)=1
 if (rfddk ==2) rfpert(natom+6)=1

 if (rfelfd==1.or.rfelfd==2) rfpert(natom+1)=1
 if (rfelfd==1.or.rfelfd==3) rfpert(natom+2)=1
 if (rfstrs==1.or.rfstrs==3) rfpert(natom+3)=1
 if (rfstrs==2.or.rfstrs==3) rfpert(natom+4)=1
 if (rfuser==1.or.rfuser==3) rfpert(natom+5)=1
 if (rfuser==2.or.rfuser==3) rfpert(natom+6)=1
 
 qeq0=(Dtset%qptn(1)**2+Dtset%qptn(2)**2+Dtset%qptn(3)**2<1.d-14)
!
!========================================
!Determine the symmetrical perturbations
!========================================
!
!* Get symmetries in reciprocal space.
 ABI_ALLOCATE(symrec,(3,3,nsym))
 do sym=1,nsym
   call mati3inv(Dtset%symrel(:,:,sym),symrec(:,:,sym))
 end do
!
!Symmetry tables for the atoms. 
 ABI_ALLOCATE(indsym,(4,nsym,natom))

 call symatm(indsym,natom,nsym,symrec,Dtset%tnons,tolsym8,Dtset%typat,Dtset%xred_orig)
!
!Examine the symmetries of the q wavevector
 ABI_ALLOCATE(symq,(4,2,nsym))
 call symq3(nsym,Dtset%qptn,symq,symrec,timrev)
!
!Find the set of base perturbations.
 ABI_ALLOCATE(pertsy,(3,mpert))

 call syper3(indsym,mpert,natom,nsym,pertsy,rfdir,rfpert,symq,symrec,Dtset%symrel)

 call wrtout(out_unt,' The list of irreducible perturbations for this q vector is:','COLL')
 idx=1
 do ipert=1,mpert
   do idir=1,3
!    if (rfpert(ipert)==1.and.rfdir(idir)==1) then
     if ( pertsy(idir,ipert)==1 ) then ! basis perturbation.
       write(msg, '(i5,a,i2,a,i4)' )idx,')    idir=',idir,'    ipert=',ipert
       call wrtout(out_unt,msg,'COLL')
       idx=idx+1
     end if
!    end if
   end do
 end do
!
!Determine the extended IBZ associated to the perturbation.
 nkpt = Dtset%nkpt
 ABI_ALLOCATE(nkrbz,(3,mpert))
 nkrbz=0

 do ipert=1,mpert
   do idir=1,3
     if ( pertsy(idir,ipert)==1 ) then ! basis perturbation.
!      
!      Determine the subset of symmetry operations (nsym1 operations) that leaves 
!      the perturbation invariant, and initialize corresponding arrays symaf1, symrl1, tnons1.
       ABI_ALLOCATE(symaf1_tmp,(nsym))
       ABI_ALLOCATE(symrl1_tmp,(3,3,nsym))
       ABI_ALLOCATE(tnons1_tmp,(3,nsym))

!      MJV TODO: check whether prepgkk should be used here
       if (Dtset%prepanl /= 1 .and. Dtset%berryopt /=4 ) then
         call insy3(gprimd,idir,indsym,ipert,natom,nsym,nsym1,rfmeth2,&
&         Dtset%symafm,symaf1_tmp,symq,symrec,&
&         Dtset%symrel,symrl1_tmp,syuse0,Dtset%tnons,tnons1_tmp)
       else
         nsym1             = 1
         symaf1_tmp(1)     = 1
         symrl1_tmp(:,:,1) = Dtset%symrel(:,:,1)
         tnons1_tmp(:,1)   = zero
       end if

       ABI_ALLOCATE(symrc1,(3,3,nsym1))
       ABI_ALLOCATE(symaf1,(nsym1))
       ABI_ALLOCATE(symrl1,(3,3,nsym1))
       ABI_ALLOCATE(tnons1,(3,nsym1))

       symaf1(1:nsym1)     = symaf1_tmp(1:nsym1)
       symrl1(:,:,1:nsym1) = symrl1_tmp(:,:,1:nsym1)
       tnons1(:,1:nsym1)   = tnons1_tmp(:,1:nsym1)
!      
!      Get the symmetry matrices in terms of reciprocal basis.
       do sym=1,nsym1
         call mati3inv(symrl1(:,:,sym),symrc1(:,:,sym))
       end do

       ABI_DEALLOCATE(symaf1_tmp)
       ABI_DEALLOCATE(symrl1_tmp)
       ABI_DEALLOCATE(tnons1_tmp)
       ABI_DEALLOCATE(symaf1)
       ABI_DEALLOCATE(symrl1)
       ABI_DEALLOCATE(tnons1)
!      
!      Determine the subset of k-points needed in the "reduced Brillouin zone", and initialize other quantities
       ABI_ALLOCATE(indkpt1_tmp,(nkpt))
       ABI_ALLOCATE(wtk_folded,(nkpt))
       indkpt1_tmp(:)=0
       timrev_pert=timrev

       if (Dtset%ieig2rf>0) then
         call symkpt(0,gmet,indkpt1_tmp,dev_null,Dtset%kptns,nkpt,nkpt_rbz,&
&         1,symrc1,0,Dtset%wtk,wtk_folded)
       else
!        For the time being, the time reversal symmetry is not used for ddk, elfd, mgfd perturbations.
         if (ipert==Dtset%natom+1 .or. ipert==Dtset%natom+2 .or. &
&         ipert==Dtset%natom+5 .or. Dtset%berryopt==4 ) timrev_pert=0

         call symkpt(0,gmet,indkpt1_tmp,dev_null,Dtset%kptns,nkpt,nkpt_rbz,&
         nsym1,symrc1,timrev_pert,Dtset%wtk,wtk_folded)
       end if

       ABI_DEALLOCATE(symrc1)

       ABI_ALLOCATE(indkpt1,(nkpt_rbz))
       ABI_ALLOCATE(kpt_rbz,(3,nkpt_rbz))
       ABI_ALLOCATE(wtk_rbz,(nkpt_rbz))

       indkpt1(:)=indkpt1_tmp(1:nkpt_rbz)
       do ikpt=1,nkpt_rbz
         kpt_rbz(:,ikpt) = Dtset%kptns(:,indkpt1(ikpt))
         wtk_rbz(ikpt)   = wtk_folded(indkpt1(ikpt))
       end do

       ABI_DEALLOCATE(indkpt1_tmp)
       ABI_DEALLOCATE(wtk_folded)
       ABI_DEALLOCATE(indkpt1)
       ABI_DEALLOCATE(kpt_rbz)
       ABI_DEALLOCATE(wtk_rbz)
!      
!      Save the number of points in the RBZ for this perturbation.
       nkrbz(idir,ipert) = nkpt_rbz
!      
     end if
   end do ! idir
 end do ! ipert

 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(symq)
 ABI_DEALLOCATE(rfpert)
!
!Printout of the irreducible perturbations.
 call wrtout(out_unt,'<BEGIN IRREDUCIBLE_PERTURBATIONS>','COLL')

 write(msg,'(a,3es16.8)')' qpt= ',Dtset%qptn(:)
 call wrtout(out_unt,msg,'COLL')

 nirred = COUNT(pertsy==1) 
 write(msg,'(a,i0)')' Number of irreducible perturbations= ',nirred
 call wrtout(out_unt,msg,'COLL')

 idx=1
 do ipert=1,mpert
   do idir=1,3
!    if (rfpert(ipert)==1.and.rfdir(idir)==1) then
     if ( pertsy(idir,ipert)==1 ) then ! basis perturbation.
       write(msg,'(i0,3(a,i0))')idx,') idir= ',idir,' ipert= ',ipert,' nkrbz= ',nkrbz(idir,ipert)
       call wrtout(out_unt,msg,'COLL')
       idx=idx+1
     end if
!    end if
   end do
 end do

 call wrtout(out_unt,'<END IRREDUCIBLE_PERTURBATIONS>','COLL')

 ABI_DEALLOCATE(nkrbz)
 ABI_DEALLOCATE(pertsy)

end subroutine irred_perts
!!***
