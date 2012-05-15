!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhoij_utils
!! This module contains functions used to manipulate
!! variables of structured datatype pawrhoij_type.
!! pawrhoij_type variables are rhoij occupancies
!! matrixes used within PAW formalism
!!
!!***

!!****f* ABINIT/rhoij_alloc
!! NAME
!! rhoij_alloc
!!
!! FUNCTION
!! Initialize and allocate a pawrhoij datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! cplex=1 if rhoij is REAL,2 if COMPLEX
!! nlmn(:)=array of (l,m,n) sizes for rhoij for each type of atom.
!! nspden=number of spin-components for rhoij
!! nsppol=number of spinorial components for rhoij
!! nsppol=number of independant spin-components for rhoij
!! typat(:)=types of atoms
!! mpi_enreg=informations about MPI parallelization (OPTIONAL)
!! ngrhoij=number of gradients to be allocated (OPTIONAL, default=0)
!! nlmnmix=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
!! use_rhoij_=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
!! use_rhoijres=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! PARENTS
!!      bethe_salpeter,extraprho,hdr_comm,hdr_io_netcdf,initrhoij,loper3
!!      m_electronpositron,m_header,m_qparticles,nstpaw3,paw_qpscgw,rhoij_utils
!!      screening,setup_bse,setup_positron,setup_screening,setup_sigma,sigma
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rhoij_alloc(cplex,nlmn,nspden,nspinor,nsppol,pawrhoij,typat,&          ! Mandatory arguments
&                      mpi_enreg,ngrhoij,nlmnmix,use_rhoij_,use_rhoijres) ! Optional arguments

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_alloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden,nspinor,nsppol
 integer,intent(in),optional :: ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
 type(MPI_type),intent(inout),optional :: mpi_enreg
!arrays
 integer,intent(in) :: nlmn(:),typat(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,itypat,lmn2_size,nn1,nn2,nrhoij

! *************************************************************************

 nrhoij=size(pawrhoij);nn1=size(nlmn);nn2=size(typat)
 if (nrhoij>nn2.or.maxval(typat)>nn1) stop "Error in rhoij_alloc: wrong sizes ! "

 do irhoij=1,nrhoij

   itypat=typat(irhoij)
   if (present(mpi_enreg)) then
     if (mpi_enreg%nproc_atom>1) itypat=typat(mpi_enreg%atom_indx(irhoij))
   end if

   lmn2_size=nlmn(itypat)*(nlmn(itypat)+1)/2

!  Scalars initializations
   pawrhoij(irhoij)%cplex=cplex
   pawrhoij(irhoij)%itypat=itypat
   pawrhoij(irhoij)%lmn_size=nlmn(itypat)
   pawrhoij(irhoij)%lmn2_size=lmn2_size
   pawrhoij(irhoij)%nspden=nspden
   pawrhoij(irhoij)%nspinor=nspinor
   pawrhoij(irhoij)%nsppol=nsppol
   pawrhoij(irhoij)%nrhoijsel=0
   pawrhoij(irhoij)%lmnmix_sz=0
   pawrhoij(irhoij)%ngrhoij=0
   pawrhoij(irhoij)%use_rhoij_=0
   pawrhoij(irhoij)%use_rhoijres=0

!  Mandatory pointers allocations
   ABI_ALLOCATE(pawrhoij(irhoij)%rhoijselect,(lmn2_size))
   ABI_ALLOCATE(pawrhoij(irhoij)%rhoijp,(cplex*lmn2_size,nspden))
   pawrhoij(irhoij)%rhoijselect(:)=0
   pawrhoij(irhoij)%rhoijp(:,:)=zero

!  Optional pointers allocations
   if (present(ngrhoij)) then
     if (ngrhoij>0) then
       pawrhoij(irhoij)%ngrhoij=ngrhoij
       ABI_ALLOCATE(pawrhoij(irhoij)%grhoij,(ngrhoij,cplex*lmn2_size,nspden))
       pawrhoij(irhoij)%grhoij=zero
     else
       nullify(pawrhoij(irhoij)%grhoij)
     end if
   else
     nullify(pawrhoij(irhoij)%grhoij)
   end if
   if (present(nlmnmix)) then
     if (nlmnmix>0) then
       pawrhoij(irhoij)%lmnmix_sz=nlmnmix
       ABI_ALLOCATE(pawrhoij(irhoij)%kpawmix,(nlmnmix))
       pawrhoij(irhoij)%kpawmix=0
     else
       nullify(pawrhoij(irhoij)%kpawmix)
     end if
   else
     nullify(pawrhoij(irhoij)%kpawmix)
   end if
   if (present(use_rhoij_)) then
     if (use_rhoij_>0) then
       pawrhoij(irhoij)%use_rhoij_=use_rhoij_
       ABI_ALLOCATE(pawrhoij(irhoij)%rhoij_,(cplex*lmn2_size,nspden))
       pawrhoij(irhoij)%rhoij_=zero
     else
       nullify(pawrhoij(irhoij)%rhoij_)
     end if
   else
     nullify(pawrhoij(irhoij)%rhoij_)
   end if
   if (present(use_rhoijres)) then
     if (use_rhoijres>0) then
       pawrhoij(irhoij)%use_rhoijres=use_rhoijres
       ABI_ALLOCATE(pawrhoij(irhoij)%rhoijres,(cplex*lmn2_size,nspden))
       pawrhoij(irhoij)%rhoijres=zero
     else
       nullify(pawrhoij(irhoij)%rhoijres)
     end if
   else
     nullify(pawrhoij(irhoij)%rhoijres)
   end if

 end do

end subroutine rhoij_alloc
!!***

!!****f* ABINIT/rhoij_free
!! NAME
!! rhoij_free
!!
!! FUNCTION
!! Destroy a pawrhoij datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! PARENTS
!!      bethe_salpeter,gstate,gw_tools,loper3,m_electronpositron,m_header
!!      m_scf_history,nstpaw3,respfn,screening,setup_bse,setup_positron
!!      setup_screening,setup_sigma,sigma
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

subroutine rhoij_free(pawrhoij)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,nrhoij

! *************************************************************************

 nrhoij=size(pawrhoij)

 if (nrhoij>0) then
   do irhoij=1,nrhoij
     pawrhoij(irhoij)%nrhoijsel=0
     pawrhoij(irhoij)%ngrhoij=0
     pawrhoij(irhoij)%lmnmix_sz=0
     pawrhoij(irhoij)%use_rhoij_=0
     pawrhoij(irhoij)%use_rhoijres=0
     if (associated(pawrhoij(irhoij)%rhoijp))       then
       ABI_DEALLOCATE(pawrhoij(irhoij)%rhoijp)
     end if
     if (associated(pawrhoij(irhoij)%rhoijselect))  then
       ABI_DEALLOCATE(pawrhoij(irhoij)%rhoijselect)
     end if
     if (associated(pawrhoij(irhoij)%grhoij))       then
       ABI_DEALLOCATE(pawrhoij(irhoij)%grhoij)
     end if
     if (associated(pawrhoij(irhoij)%kpawmix))      then
       ABI_DEALLOCATE(pawrhoij(irhoij)%kpawmix)
     end if
     if (associated(pawrhoij(irhoij)%rhoij_))       then
       ABI_DEALLOCATE(pawrhoij(irhoij)%rhoij_)
     end if
     if (associated(pawrhoij(irhoij)%rhoijres))     then
       ABI_DEALLOCATE(pawrhoij(irhoij)%rhoijres)
     end if
   end do
 end if

end subroutine rhoij_free
!!***

!!****f* ABINIT/rhoij_nullify
!! NAME
!! rhoij_nullify
!!
!! FUNCTION
!! Nullify (initialize to null) a pawrhoij datastructure
!!
!! COPYRIGHT
!! Copyright (C) 20011-2011 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! PARENTS
!!      m_scf_history
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

subroutine rhoij_nullify(pawrhoij)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,nrhoij

! *************************************************************************

 nrhoij=size(pawrhoij)

 if (nrhoij>0) then
   do irhoij=1,nrhoij
     pawrhoij(irhoij)%nrhoijsel=0
     pawrhoij(irhoij)%ngrhoij=0
     pawrhoij(irhoij)%lmnmix_sz=0
     pawrhoij(irhoij)%use_rhoij_=0
     pawrhoij(irhoij)%use_rhoijres=0
     nullify(pawrhoij(irhoij)%rhoijp)
     nullify(pawrhoij(irhoij)%rhoijselect)
     nullify(pawrhoij(irhoij)%grhoij)
     nullify(pawrhoij(irhoij)%kpawmix)
     nullify(pawrhoij(irhoij)%rhoij_)
     nullify(pawrhoij(irhoij)%rhoijres)
   end do
 end if

end subroutine rhoij_nullify
!!***

!!****f* ABINIT/rhoij_copy
!! NAME
!! rhoij_copy
!!
!! FUNCTION
!! Copy one pawrhoij datastructure into another
!! Can take into accound changes of dimensions
!! Can copy a shared pawrhoij into distributed ones (when parallelism is activated)
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! keep_cplex= optional argument (logical, default=.TRUE.)
!!             if .TRUE. pawrhoij_out(:)%cplex is NOT MODIFIED,
!!             even if different from pawrhoij_in(:)%cplex
!! keep_nspden= optional argument (logical, default=.TRUE.)
!!              if .TRUE. pawrhoij_out(:)%nspden is NOT MODIFIED,
!!              even if different from pawrhoij_in(:)%nspden
!! mpi_enreg= optional argument: informations about MPI parallelization:
!! pawrhoij_in(:)<type(pawrhoij_type)>= input rhoij datastructure
!!
!! SIDE EFFECTS
!! pawrhoij_out(:)<type(pawrhoij_type)>= output rhoij datastructure
!!
!! NOTES
!! When mpi_enreg is present, the input pawrhoij datastructure is distributed
!! on every procs. (output pawrhoij should be part of input one).
!!
!! PARENTS
!!      bethe_salpeter,gstate,hdr_update,inwffil,ioarr,loper3
!!      m_electronpositron,m_header,respfn,rhoij_utils,screening,setup_bse
!!      setup_positron,setup_screening,setup_sigma,sigma
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

subroutine rhoij_copy(pawrhoij_in,pawrhoij_out, &
&                     keep_cplex,keep_nspden,mpi_enreg) ! optional arguments

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in),optional :: keep_cplex,keep_nspden
 type(MPI_type),intent(inout),optional :: mpi_enreg
!arrays
 type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij_out(:)

!Local variables-------------------------------
!scalars
 integer :: cplex_in,cplex_out,dplex_in,dplex_out,i_in,i_out,ilmn,irhoij,ispden,jrhoij
 integer :: lmn2_size_out,lmnmix,ngrhoij,nrhoij_in,nrhoij_out,nselect,nspden_in
 integer :: nspden_out,use_rhoij_,use_rhoijres
 logical :: change_dim,keep_cplex_,keep_nspden_,use_index
 character(len=500) :: msg

! *************************************************************************

!Retrieve sizes
 nrhoij_in=size(pawrhoij_in);nrhoij_out=size(pawrhoij_out)

!Init flags
 keep_cplex_=.true.
 if (present(keep_cplex)) keep_cplex_=keep_cplex
 keep_nspden_=.true.
 if (present(keep_nspden)) keep_nspden_=keep_nspden
 use_index=.false.
 if (present(mpi_enreg)) use_index=((mpi_enreg%nproc_atom>1).and.(nrhoij_in>nrhoij_out))

!Test on sizes
 if (nrhoij_in<nrhoij_out.and.(.not.use_index)) stop "Error in rhoij_copy: wrong sizes ! "

!Loop on rhoij components
 do irhoij=1,nrhoij_out
   jrhoij=irhoij;if (use_index) jrhoij=mpi_enreg%atom_indx(irhoij)

   lmn2_size_out=pawrhoij_in(jrhoij)%lmn2_size
   cplex_in=pawrhoij_in(jrhoij)%cplex
   cplex_out=cplex_in;if(keep_cplex_)cplex_out=pawrhoij_out(irhoij)%cplex
   nspden_in=pawrhoij_in(jrhoij)%nspden
   nspden_out=nspden_in;if(keep_nspden_)nspden_out=pawrhoij_out(irhoij)%nspden

   change_dim=(pawrhoij_out(irhoij)%cplex/=cplex_out.or. &
&   pawrhoij_out(irhoij)%lmn2_size/=lmn2_size_out.or.&
&   pawrhoij_out(irhoij)%nspden/=nspden_out)
   dplex_in=cplex_in-1;dplex_out=cplex_out-1

!  Scalars
   nselect=pawrhoij_in(irhoij)%nrhoijsel
   pawrhoij_out(irhoij)%cplex=cplex_out+0
   pawrhoij_out(irhoij)%nspden=nspden_out+0
   pawrhoij_out(irhoij)%lmn2_size=lmn2_size_out+0
   pawrhoij_out(irhoij)%lmn_size=pawrhoij_in(jrhoij)%lmn_size+0
   if(.not.keep_nspden_) pawrhoij_out(irhoij)%nsppol =pawrhoij_in(jrhoij)%nsppol+0
   if(.not.keep_nspden_) pawrhoij_out(irhoij)%nspinor=pawrhoij_in(jrhoij)%nspinor+0
   pawrhoij_out(irhoij)%nrhoijsel=nselect+0

!  pawrhoij_out(irhoij)%itypat=pawrhoij_in(jrhoij)%itypat+0
   if (pawrhoij_out(irhoij)%itypat/=pawrhoij_in(jrhoij)%itypat) then
     write(unit=msg,fmt='(a,i3,a)') 'Type of atom ',jrhoij,' is different (dont copy it) !'
     MSG_COMMENT(msg)
   end if

!  Mandatory pointer: indexes for non-zero elements selection
   if (change_dim) then
     ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijselect)
     ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijselect,(lmn2_size_out))
   end if
   pawrhoij_out(irhoij)%rhoijselect(1:nselect)=pawrhoij_in(jrhoij)%rhoijselect(1:nselect)+0
   if (nselect<lmn2_size_out) pawrhoij_out(irhoij)%rhoijselect(nselect+1:lmn2_size_out)=0

!  Mandatory pointer: non-zero elements of rhoij
   if (change_dim) then
     ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijp)
     ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijp,(cplex_out*lmn2_size_out,nspden_out))
   end if
   if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
     do ispden=1,nspden_out
       pawrhoij_out(irhoij)%rhoijp(1:cplex_out*nselect,ispden)=pawrhoij_in(jrhoij)%rhoijp(1:cplex_out*nselect,ispden)+zero
       if (nselect<lmn2_size_out) pawrhoij_out(irhoij)%rhoijp(cplex_out*nselect+1:cplex_out*lmn2_size_out,ispden)=zero
     end do
   else
     pawrhoij_out(irhoij)%rhoijp(:,:)=zero
     if (nspden_out==1) then
       if (nspden_in==2) then
         do ilmn=1,nselect
           i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
           pawrhoij_out(irhoij)%rhoijp(i_out,1)=pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&           +pawrhoij_in(jrhoij)%rhoijp(i_in,2)+zero
         end do
       else ! nspden_in==1 or nspden_in=4
         do ilmn=1,nselect
           i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
           pawrhoij_out(irhoij)%rhoijp(i_out,1)=pawrhoij_in(jrhoij)%rhoijp(i_in,1)+zero
         end do
       end if
     else if (nspden_out==2) then
       if (nspden_in==1) then
         do ilmn=1,nselect
           i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
           pawrhoij_out(irhoij)%rhoijp(i_out,1)=half*pawrhoij_in(jrhoij)%rhoijp(i_in,1)+zero
           pawrhoij_out(irhoij)%rhoijp(i_out,2)=pawrhoij_in(jrhoij)%rhoijp(i_out,1)
         end do
       else if (nspden_in==2) then
         do ilmn=1,nselect
           i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
           pawrhoij_out(irhoij)%rhoijp(i_out,1:2)=pawrhoij_in(jrhoij)%rhoijp(i_in,1:2)+zero
         end do
       else ! nspden_in==4
         do ilmn=1,nselect
           i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
           pawrhoij_out(irhoij)%rhoijp(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&           +pawrhoij_in(jrhoij)%rhoijp(i_in,4))+zero
           pawrhoij_out(irhoij)%rhoijp(i_out,2)=half*(pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&           -pawrhoij_in(jrhoij)%rhoijp(i_in,4))+zero
         end do
       end if
     else if (nspden_out==4) then
       if (nspden_in==1) then
         do ilmn=1,nselect
           i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
           pawrhoij_out(irhoij)%rhoijp(i_out,1)=pawrhoij_in(jrhoij)%rhoijp(i_in,1)+zero
         end do
       else if (nspden_in==2) then
         do ilmn=1,nselect
           i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
           pawrhoij_out(irhoij)%rhoijp(i_out,1)=pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&           +pawrhoij_in(jrhoij)%rhoijp(i_in,2)+zero
           pawrhoij_out(irhoij)%rhoijp(i_out,4)=pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&           -pawrhoij_in(jrhoij)%rhoijp(i_in,2)+zero
         end do
       else ! nspden_in==4
         do ilmn=1,nselect
           i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
           pawrhoij_out(irhoij)%rhoijp(i_out,1:4)=pawrhoij_in(jrhoij)%rhoijp(i_in,1:4)+zero
         end do
       end if
     end if
   end if

!  Optional pointer: indexes of rhoij to be mixed
   lmnmix=pawrhoij_in(jrhoij)%lmnmix_sz
   if (pawrhoij_out(irhoij)%lmnmix_sz/=lmnmix) then
     if (pawrhoij_out(irhoij)%lmnmix_sz>0)  then
       ABI_DEALLOCATE(pawrhoij_out(irhoij)%kpawmix)
     end if
     if (lmnmix>0)  then
       ABI_ALLOCATE(pawrhoij_out(irhoij)%kpawmix,(lmnmix))
     end if
     pawrhoij_out(irhoij)%lmnmix_sz=lmnmix
   end if
   if (lmnmix>0) pawrhoij_out(irhoij)%kpawmix(1:lmnmix)=pawrhoij_in(jrhoij)%kpawmix(1:lmnmix)

!  Optional pointer: gradients of rhoij
   ngrhoij=pawrhoij_in(jrhoij)%ngrhoij
   if (pawrhoij_out(irhoij)%ngrhoij/=ngrhoij) then
     if (pawrhoij_out(irhoij)%ngrhoij>0)  then
       ABI_DEALLOCATE(pawrhoij_out(irhoij)%grhoij)
     end if
     if (ngrhoij>0)  then
       ABI_ALLOCATE(pawrhoij_out(irhoij)%grhoij,(ngrhoij,cplex_out*lmn2_size_out,nspden_out))
     end if
     pawrhoij_out(irhoij)%ngrhoij=ngrhoij
   end if
   if (ngrhoij>0) then
     if (change_dim) then
       ABI_DEALLOCATE(pawrhoij_out(irhoij)%grhoij)
       ABI_ALLOCATE(pawrhoij_out(irhoij)%grhoij,(ngrhoij,cplex_out*lmn2_size_out,nspden_out))
     end if
     if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
       do ispden=1,nspden_out
         do ilmn=1,cplex_out*lmn2_size_out
           pawrhoij_out(irhoij)%grhoij(1:ngrhoij,ilmn,ispden)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,ilmn,ispden)
         end do
       end do
     else
       pawrhoij_out(irhoij)%grhoij(:,:,:)=zero
       if (nspden_out==1) then
         if (nspden_in==2) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&             +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,2)
           end do
         else ! nspden_in==4
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1)
           end do
         end if
       else if (nspden_out==2) then
         if (nspden_in==1) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&             +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,2))
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,2)=pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)
           end do
         else ! nspden_in==4
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&             +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,4))
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,2)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&             -pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,4))
           end do
         end if
       else if (nspden_out==4) then
         if (nspden_in==1) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1)
           end do
         else ! nspden_in==2
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&             +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,2))
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,4)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&             -pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,2))
           end do
         end if
       end if
     end if
   end if

!  Optional pointer: residuals of rhoij
   use_rhoijres=pawrhoij_in(jrhoij)%use_rhoijres
   if (pawrhoij_out(irhoij)%use_rhoijres/=use_rhoijres) then
     if (pawrhoij_out(irhoij)%use_rhoijres>0)  then
       ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijres)
     end if
     if (use_rhoijres>0)  then
       ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijres,(cplex_out*lmn2_size_out,nspden_out))
     end if
     pawrhoij_out(irhoij)%use_rhoijres=use_rhoijres
   end if
   if (use_rhoijres>0) then
     if (change_dim) then
       ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijres)
       ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijres,(cplex_out*lmn2_size_out,nspden_out))
     end if
     if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
       do ispden=1,nspden_out
         do ilmn=1,cplex_out*lmn2_size_out
           pawrhoij_out(irhoij)%rhoijres(ilmn,ispden)=pawrhoij_in(jrhoij)%rhoijres(ilmn,ispden)
         end do
       end do
     else
       pawrhoij_out(irhoij)%rhoijres(:,:)=zero
       if (nspden_out==1) then
         if (nspden_in==2) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoijres(i_out,1)=pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&             +pawrhoij_in(jrhoij)%rhoijres(i_in,2)
           end do
         else ! nspden_in==4
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoijres(i_out,1)=pawrhoij_in(jrhoij)%rhoijres(i_in,1)
           end do
         end if
       else if (nspden_out==2) then
         if (nspden_in==1) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoijres(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&             +pawrhoij_in(jrhoij)%rhoijres(i_in,2))
             pawrhoij_out(irhoij)%rhoijres(i_out,2)=pawrhoij_out(irhoij)%rhoijres(i_out,1)
           end do
         else ! nspden_in==4
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoijres(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&             +pawrhoij_in(jrhoij)%rhoijres(i_in,4))
             pawrhoij_out(irhoij)%rhoijres(i_out,2)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&             -pawrhoij_in(jrhoij)%rhoijres(i_in,4))
           end do
         end if
       else if (nspden_out==4) then
         if (nspden_in==1) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoijres(i_out,1)=pawrhoij_in(jrhoij)%rhoijres(i_in,1)
           end do
         else ! nspden_in==2
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoijres(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&             +pawrhoij_in(jrhoij)%rhoijres(i_in,2))
             pawrhoij_out(irhoij)%rhoijres(i_out,4)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&             -pawrhoij_in(jrhoij)%rhoijres(i_in,2))
           end do
         end if
       end if
     end if
   end if

!  Optional pointer: non-symmetrized rhoij
   use_rhoij_=pawrhoij_in(jrhoij)%use_rhoij_
   if (pawrhoij_out(irhoij)%use_rhoij_/=use_rhoij_) then
     if (pawrhoij_out(irhoij)%use_rhoij_>0)  then
       ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoij_)
     end if
     if (use_rhoij_>0)  then
       ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoij_,(cplex_out*lmn2_size_out,nspden_out))
     end if
     pawrhoij_out(irhoij)%use_rhoij_=use_rhoij_
   end if
   if (use_rhoij_>0) then
     if (change_dim) then
       ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoij_)
       ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoij_,(cplex_out*lmn2_size_out,nspden_out))
     end if
     if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
       do ispden=1,nspden_out
         do ilmn=1,cplex_out*lmn2_size_out
           pawrhoij_out(irhoij)%rhoij_(ilmn,ispden)=pawrhoij_in(jrhoij)%rhoij_(ilmn,ispden)
         end do
       end do
     else
       pawrhoij_out(irhoij)%rhoij_(:,:)=zero
       if (nspden_out==1) then
         if (nspden_in==2) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoij_(i_out,1)=pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&             +pawrhoij_in(jrhoij)%rhoij_(i_in,2)
           end do
         else ! nspden_in==4
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoij_(i_out,1)=pawrhoij_in(jrhoij)%rhoij_(i_in,1)
           end do
         end if
       else if (nspden_out==2) then
         if (nspden_in==1) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoij_(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&             +pawrhoij_in(jrhoij)%rhoij_(i_in,2))
             pawrhoij_out(irhoij)%rhoij_(i_out,2)=pawrhoij_out(irhoij)%rhoij_(i_out,1)
           end do
         else ! nspden_in==4
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoij_(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&             +pawrhoij_in(jrhoij)%rhoij_(i_in,4))
             pawrhoij_out(irhoij)%rhoij_(i_out,2)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&             -pawrhoij_in(jrhoij)%rhoij_(i_in,4))
           end do
         end if
       else if (nspden_out==4) then
         if (nspden_in==1) then
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoij_(i_out,1)=pawrhoij_in(jrhoij)%rhoij_(i_in,1)
           end do
         else ! nspden_in==2
           do ilmn=1,lmn2_size_out
             i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
             pawrhoij_out(irhoij)%rhoij_(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&             +pawrhoij_in(jrhoij)%rhoij_(i_in,2))
             pawrhoij_out(irhoij)%rhoij_(i_out,4)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&             -pawrhoij_in(jrhoij)%rhoij_(i_in,2))
           end do
         end if
       end if
     end if
   end if

 end do ! irhoij

end subroutine rhoij_copy
!!***

!!****f* ABINIT/rhoij_allgather
!! NAME
!! rhoij_allgather
!!
!! FUNCTION
!! Gather pawrhoij datastructures from every process to every process
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! mpi_enreg=informations about MPI parallelization:
!! pawrhoij_in(:)<type(pawrhoij_type)>= input rhoij datastructures on every process
!!
!! OUTPUT
!! pawrhoij_gathered(:)<type(pawrhoij_type)>= output rhoij datastructure
!!
!! PARENTS
!!      hdr_update
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

subroutine rhoij_allgather(mpi_enreg,pawrhoij_in,pawrhoij_gathered)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_allgather'
 use interfaces_44_abitypes_defs, except_this_one => rhoij_allgather
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij_gathered(:)

!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_dp_size_all,buf_int_size,buf_int_size_all
 integer :: cplex,ierr,ii,indx_dp,indx_int,irhoij,isp,jrhoij,lmn2_size,lmnmix,lmnmix_out,natom
 integer :: ngrhoij,ngrhoij_out,nrhoij_in,nrhoij_in_sum,nrhoij_out,nselect,nspden
 integer :: use_rhoijres,use_rhoijres_out,use_rhoij_,use_rhoij_out_
 integer,allocatable :: atm_indx_all(:)
 integer,allocatable :: buf_int(:),buf_int_all(:)
 integer,allocatable :: count_dp(:),count_int(:)
 integer,allocatable :: disp_dp(:),disp_int(:)
 real(dp),allocatable :: buf_dp(:),buf_dp_all(:)

! *************************************************************************

!Without parallelism, just copy input to output
 if (mpi_enreg%nproc_atom==1) then
   call rhoij_copy(pawrhoij_in,pawrhoij_gathered)
   return
 end if

!Test on sizes
 nrhoij_in=size(pawrhoij_in);nrhoij_out=size(pawrhoij_gathered)
 nrhoij_in_sum=nrhoij_in
 call xsum_mpi(nrhoij_in_sum,mpi_enreg%comm_atom,ierr)
 if (nrhoij_in_sum/=nrhoij_out) stop "Error (1) in rhoij_allgather: wrong sizes ! "

!Tests on scalars
 lmnmix_out=pawrhoij_gathered(1)%lmnmix_sz
 ngrhoij_out=pawrhoij_gathered(1)%ngrhoij
 use_rhoijres_out=pawrhoij_gathered(1)%use_rhoijres
 use_rhoij_out_=pawrhoij_gathered(1)%use_rhoij_
 if (pawrhoij_in(1)%cplex/=pawrhoij_gathered(1)%cplex) stop "Error (2) in rhoij_allgather: wrong cplex values ! "
 if (pawrhoij_in(1)%lmn2_size/=pawrhoij_gathered(1)%lmn2_size) stop "Error (3) in rhoij_allgather: wrong lmn2_size values ! "
 if (pawrhoij_in(1)%nspden/=pawrhoij_gathered(1)%nspden) stop "Error (4) in rhoij_allgather: wrong nspden values ! "

!Compute sizes of buffers
 buf_int_size=12
 buf_dp_size =0
 do irhoij=1,nrhoij_in
   buf_int_size=buf_int_size &
&   +pawrhoij_in(irhoij)%nrhoijsel &
&   +pawrhoij_in(irhoij)%lmnmix_sz
   buf_dp_size=buf_dp_size &
&   +pawrhoij_in(irhoij)%cplex*pawrhoij_in(irhoij)%nrhoijsel*pawrhoij_in(irhoij)%nspden &
&   +pawrhoij_in(irhoij)%cplex*pawrhoij_in(irhoij)%lmn2_size*pawrhoij_in(irhoij)%nspden &
   *(1+pawrhoij_in(irhoij)%ngrhoij+pawrhoij_in(irhoij)%use_rhoijres+pawrhoij_in(irhoij)%use_rhoij_)
 end do

!Fill input buffers
 ABI_ALLOCATE(buf_int,(buf_int_size))
 indx_int=1
 ABI_ALLOCATE(buf_dp ,(buf_dp_size))
 indx_dp =1
 do irhoij=1,nrhoij_in
   cplex       =pawrhoij_in(irhoij)%cplex
   nselect     =pawrhoij_in(irhoij)%nrhoijsel
   nspden      =pawrhoij_in(irhoij)%nspden
   lmn2_size   =pawrhoij_in(irhoij)%lmn2_size
   lmnmix      =pawrhoij_in(irhoij)%lmnmix_sz;if (lmnmix/=lmnmix_out) lmnmix=0
   ngrhoij     =pawrhoij_in(irhoij)%ngrhoij;if (ngrhoij/=ngrhoij_out) ngrhoij=0
   use_rhoijres=pawrhoij_in(irhoij)%use_rhoijres;if (use_rhoijres/=use_rhoijres_out) use_rhoijres=0
   use_rhoij_  =pawrhoij_in(irhoij)%use_rhoij_;if (use_rhoij_/=use_rhoij_out_) use_rhoij_=0
   buf_int(indx_int)=cplex                           ;indx_int=indx_int+1
   buf_int(indx_int)=nselect                         ;indx_int=indx_int+1
   buf_int(indx_int)=nspden                          ;indx_int=indx_int+1
   buf_int(indx_int)=lmn2_size                       ;indx_int=indx_int+1
   buf_int(indx_int)=lmnmix                          ;indx_int=indx_int+1
   buf_int(indx_int)=ngrhoij                         ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijres                    ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoij_                      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%itypat      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%lmn_size    ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%nsppol      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%nspinor     ;indx_int=indx_int+1
   buf_int(indx_int:indx_int+nselect-1)=pawrhoij_in(irhoij)%rhoijselect(1:nselect)
   indx_int=indx_int+nselect
   do isp=1,nspden
     buf_dp(indx_dp:indx_dp+cplex*nselect-1)=pawrhoij_in(irhoij)%rhoijp(1:cplex*nselect,isp)
     indx_dp=indx_dp+cplex*nselect
   end do
   if (lmnmix>0) then
     buf_int(indx_int:indx_int+lmnmix-1)=pawrhoij_in(irhoij)%kpawmix(1:lmnmix)
     indx_int=indx_int+lmnmix
   end if
   if (ngrhoij>0) then
     do isp=1,nspden
       do ii=1,cplex*lmn2_size
         buf_dp(indx_dp:indx_dp+ngrhoij-1)=pawrhoij_in(irhoij)%grhoij(1:ngrhoij,ii,isp)
         indx_dp=indx_dp+ngrhoij
       end do
     end do
   end if
   if (use_rhoijres>0) then
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(irhoij)%rhoijres(1:cplex*lmn2_size,isp)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
   if (use_rhoij_>0) then
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(irhoij)%rhoij_(1:cplex*lmn2_size,isp)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
 end do
 if (indx_int/=1+buf_int_size) stop "Error (5) in rhoij_allgather: wrong buffer sizes ! "
 if (indx_dp /=1+buf_dp_size)  stop "Error (6) in rhoij_allgather: wrong buffer sizes ! "

!Prepare communications
 ABI_ALLOCATE(count_int,(mpi_enreg%nproc_atom))
 ABI_ALLOCATE(disp_int,(mpi_enreg%nproc_atom))
 ABI_ALLOCATE(count_dp,(mpi_enreg%nproc_atom))
 ABI_ALLOCATE(disp_dp,(mpi_enreg%nproc_atom))
 call xallgather_mpi(buf_int_size,count_int,mpi_enreg%comm_atom,ierr)
 call xallgather_mpi(buf_dp_size ,count_dp ,mpi_enreg%comm_atom,ierr)
 disp_int(1)=0;disp_dp(1)=0
 do ii=2,mpi_enreg%nproc_atom
   disp_int(ii)=disp_int(ii-1)+count_int(ii-1)
   disp_dp (ii)=disp_dp (ii-1)+count_dp (ii-1)
 end do
 buf_int_size_all=sum(count_int)
 buf_dp_size_all =sum(count_dp)
 ABI_ALLOCATE(buf_int_all,(buf_int_size_all))
 ABI_ALLOCATE(buf_dp_all ,(buf_dp_size_all))
 call xsum_mpi(mpi_enreg%natom,natom,mpi_enreg%comm_atom,ierr)
 ABI_ALLOCATE(atm_indx_all,(natom))
 call xallgather_mpi(mpi_enreg%atom_indx,mpi_enreg%natom,atm_indx_all,mpi_enreg%comm_atom,ierr)

!Communicate !!!
 call xallgatherv_mpi(buf_int,buf_int_size, buf_int_all,count_int,disp_int, mpi_enreg%comm_atom,ierr)
 call xallgatherv_mpi(buf_dp ,buf_dp_size , buf_dp_all ,count_dp ,disp_dp , mpi_enreg%comm_atom,ierr)

!Empty output buffer
 indx_int=1;indx_dp=1
 do irhoij=1,nrhoij_out
   jrhoij=atm_indx_all(irhoij)
   cplex       =buf_int(indx_int)    ;indx_int=indx_int+1
   nselect     =buf_int(indx_int)    ;indx_int=indx_int+1
   nspden      =buf_int(indx_int)    ;indx_int=indx_int+1
   lmn2_size   =buf_int(indx_int)    ;indx_int=indx_int+1
   lmnmix      =buf_int(indx_int)    ;indx_int=indx_int+1
   ngrhoij     =buf_int(indx_int)    ;indx_int=indx_int+1
   use_rhoijres=buf_int(indx_int)    ;indx_int=indx_int+1
   use_rhoij_  =buf_int(indx_int)    ;indx_int=indx_int+1
   pawrhoij_gathered(jrhoij)%cplex=cplex
   pawrhoij_gathered(jrhoij)%nrhoijsel=nselect
   pawrhoij_gathered(jrhoij)%nspden=nspden
   pawrhoij_gathered(jrhoij)%lmn2_size=lmn2_size
   pawrhoij_gathered(jrhoij)%lmnmix_sz=lmnmix
   pawrhoij_gathered(jrhoij)%ngrhoij=ngrhoij
   pawrhoij_gathered(jrhoij)%use_rhoijres=use_rhoijres
   pawrhoij_gathered(jrhoij)%use_rhoij_=use_rhoij_
   pawrhoij_gathered(jrhoij)%itypat=buf_int(indx_int)      ;indx_int=indx_int+1
   pawrhoij_gathered(jrhoij)%lmn_size=buf_int(indx_int)    ;indx_int=indx_int+1
   pawrhoij_gathered(jrhoij)%nsppol=buf_int(indx_int)      ;indx_int=indx_int+1
   pawrhoij_gathered(jrhoij)%nspinor=buf_int(indx_int)     ;indx_int=indx_int+1
   pawrhoij_gathered(jrhoij)%rhoijselect(1:nselect)=buf_int(indx_int:indx_int+nselect-1)
   indx_int=indx_int+nselect
   do isp=1,nspden
     pawrhoij_gathered(jrhoij)%rhoijp(1:cplex*nselect,isp)=buf_dp(indx_dp:indx_dp+cplex*nselect-1)
     indx_dp=indx_dp+cplex*nselect
   end do
   if (lmnmix>0) then
     pawrhoij_gathered(jrhoij)%kpawmix(1:lmnmix)=buf_int(indx_int:indx_int+lmnmix-1)
     indx_int=indx_int+lmnmix
   end if
   if (ngrhoij>0) then
     do isp=1,nspden
       do ii=1,cplex*lmn2_size
         pawrhoij_gathered(jrhoij)%grhoij(1:ngrhoij,ii,isp)=buf_dp(indx_dp:indx_dp+ngrhoij-1)
         indx_dp=indx_dp+ngrhoij
       end do
     end do
   end if
   if (use_rhoijres>0) then
     do isp=1,nspden
       pawrhoij_gathered(jrhoij)%rhoijres(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
   if (use_rhoij_>0) then
     do isp=1,nspden
       pawrhoij_gathered(jrhoij)%rhoij_(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
 end do
 if (indx_int/=1+buf_int_size) stop "Error (7) in rhoij_allgather: wrong buffer sizes ! "
 if (indx_dp /=1+buf_dp_size)  stop "Error (8) in rhoij_allgather: wrong buffer sizes ! "

!Free memory
 ABI_DEALLOCATE(atm_indx_all)
 ABI_DEALLOCATE(buf_int)
 ABI_DEALLOCATE(buf_int_all)
 ABI_DEALLOCATE(count_int)
 ABI_DEALLOCATE(disp_int)
 ABI_DEALLOCATE(buf_dp)
 ABI_DEALLOCATE(buf_dp_all)
 ABI_DEALLOCATE(count_dp)
 ABI_DEALLOCATE(disp_dp)

end subroutine rhoij_allgather
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/rhoij_io
!! NAME
!! rhoij_io
!!
!! FUNCTION
!! IO method for pawrhoij datastructures.
!!
!! INPUTS
!!  unitfi=Unit number for IO file (already opened in the caller).
!!  nsppol_in=Number of independent spin polarizations. Only used for reading.
!!  nspinor_in=Number of spinorial components. Only used for reading.
!!  nspden_in=Number of spin-density components. only used for reading.
!!  nlmn_type(ntypat)= Number of (l,m,n) elements for the paw basis for each type of atom. Only used for reading.
!!  typat(natom) =Type of each atom.
!!  headform=Format of the abinit header (only used for reading as we need to know how to read
!!    the data. Writing is always done using the latest headform.
!!  rdwr_mode(len=*)=String defining the IO mode. Possible values (not case sensitive):
!!    "W"= For writing to unitfi
!!    "R"= For reading from unitfi
!!    "E"= For echoing.
!!  [form(len=*)]= String defining the file format. Defaults to Fortran binary mode i.e., "unformatted"
!!  Other possible values are (case insensitive):
!!    "formatted"=For IO on a file open in formatted mode.
!!  [natinc]=Defines the increment in the loop over natom used for echoing the pawrhoij(natom) datastructures.
!!    If not specified, only the first and the last atom will be printed.
!!
!! SIDE EFFECTS
!!  pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure.
!!   if rdwr_mode="W", it will be written on unit unitfi using the file format given by form.
!!   if rdwr_mode="R", pawrhoij will be read and initialized from unit unitfi that has been
!!      opened with form=form.
!!   if rdwr_mode="E", the routines only echoes the content of the structure.
!!
!! PARENTS
!!      hdr_io,m_qparticles
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

subroutine rhoij_io(pawrhoij,unitfi,nsppol_in,nspinor_in,nspden_in,nlmn_type,typat,headform,rdwr_mode,form,natinc)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_fstrings,   only : toupper

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_io'
 use interfaces_32_util
 use interfaces_44_abitypes_defs, except_this_one => rhoij_io
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unitfi,headform,nspden_in,nspinor_in,nsppol_in
 integer,optional,intent(in) :: natinc
 character(len=*),intent(in) :: rdwr_mode
 character(len=*),optional,intent(in) :: form
!arrays
 integer,intent(in) :: typat(:),nlmn_type(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: iatom,natom,ispden,bsize,ii,jj,nselect,my_cplex,my_nspden,my_natinc
 logical :: isbinary
!arrays
 integer,allocatable :: ibuffer(:),nsel44(:,:),nsel56(:)
 real(dp), allocatable :: buffer(:)

! *************************************************************************

 natom = SIZE(pawrhoij)

 isbinary=.TRUE.
 if (PRESENT(form)) then
   if (toupper(form)=="FORMATTED") isbinary=.FALSE.
 end if

 select case (rdwr_mode(1:1))

   case ("R","r") ! Reading the Rhoij tab.

     if ((headform>=44).and.(headform<56)) then

       ABI_ALLOCATE(nsel44,(nspden_in,natom))

       if (isbinary) then
         read(unitfi  ) ((nsel44(ispden,iatom),ispden=1,nspden_in),iatom=1,natom)
       else
         read(unitfi,*) ((nsel44(ispden,iatom),ispden=1,nspden_in),iatom=1,natom)
       end if

       call rhoij_alloc(1,nlmn_type,nspden_in,nspinor_in,nsppol_in,pawrhoij,typat)
       do iatom=1,natom
         pawrhoij(iatom)%nrhoijsel=nsel44(1,iatom)
       end do

       bsize=sum(nsel44)
       ABI_ALLOCATE(ibuffer,(bsize))
       ABI_ALLOCATE(buffer,(bsize))

       if (isbinary) then
         read(unitfi  ) ibuffer(:),buffer(:)
       else
         read(unitfi,*) ibuffer(:),buffer(:)
       end if

       ii=0
       do iatom=1,natom
         nselect=nsel44(1,iatom)
         pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
         do ispden=1,nspden_in
           pawrhoij(iatom)%rhoijp(1:nselect,ispden)=buffer(ii+1:ii+nselect)
           ii=ii+nselect
         end do
       end do
       ABI_DEALLOCATE(ibuffer)
       ABI_DEALLOCATE(buffer)
       ABI_DEALLOCATE(nsel44)

     else if (headform>=56) then
       ABI_ALLOCATE(nsel56,(natom))

       if (headform==56) then
         if (isbinary) then
           read(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex
         else
           read(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex
         end if
         my_nspden=nspden_in
       else
         if (isbinary) then
           read(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
         else
           read(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
         end if
       end if

       call rhoij_alloc(my_cplex,nlmn_type,my_nspden,nspinor_in,nsppol_in,pawrhoij,typat)
       do iatom=1,natom
         pawrhoij(iatom)%nrhoijsel=nsel56(iatom)
       end do
       bsize=sum(nsel56)
       ABI_ALLOCATE(ibuffer,(bsize))
       ABI_ALLOCATE(buffer,(bsize*nspden_in*my_cplex))

       if (isbinary) then
         read(unitfi  ) ibuffer(:),buffer(:)
       else
         read(unitfi,*) ibuffer(:),buffer(:)
       end if

       ii=0;jj=0
       do iatom=1,natom
         nselect=nsel56(iatom)
         pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
         ii=ii+nselect
         do ispden=1,nspden_in
           pawrhoij(iatom)%rhoijp(1:my_cplex*nselect,ispden)=buffer(jj+1:jj+my_cplex*nselect)
           jj=jj+my_cplex*nselect
         end do
       end do
       ABI_DEALLOCATE(ibuffer)
       ABI_DEALLOCATE(buffer)
       ABI_DEALLOCATE(nsel56)
     end if

   case ("W","w") ! Writing the Rhoij tab. Latest format is used.

     ABI_ALLOCATE(nsel56,(natom))
     my_cplex =pawrhoij(1)%cplex
     my_nspden=pawrhoij(1)%nspden
     do iatom=1,natom
       nsel56(iatom)=pawrhoij(iatom)%nrhoijsel
     end do

     if (isbinary) then
       write(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
     else
       write(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
     end if

     bsize=sum(nsel56)
     ABI_ALLOCATE(ibuffer,(bsize))
     ABI_ALLOCATE(buffer,(bsize*my_nspden*my_cplex))
     ii=0;jj=0
     do iatom=1,natom
       nselect=nsel56(iatom)
       ibuffer(ii+1:ii+nselect)=pawrhoij(iatom)%rhoijselect(1:nselect)
       ii=ii+nselect
       do ispden=1,my_nspden
         buffer(jj+1:jj+my_cplex*nselect)=pawrhoij(iatom)%rhoijp(1:my_cplex*nselect,ispden)
         jj=jj+my_cplex*nselect
       end do
     end do

     if (isbinary) then
       write(unitfi  ) ibuffer(:),buffer(:)
     else
       write(unitfi,*) ibuffer(:),buffer(:)
     end if

     ABI_DEALLOCATE(ibuffer)
     ABI_DEALLOCATE(buffer)
     ABI_DEALLOCATE(nsel56)

   case ("E","e") ! Echoing
     my_natinc=1; if(natom>1) my_natinc=natom-1
     if (PRESENT(natinc)) my_natinc = natinc ! user-defined increment.
     ABI_ALLOCATE(ibuffer,(0))
     do iatom=1,natom,my_natinc
       do ispden=1,pawrhoij(iatom)%nspden
         write(unitfi, '(a,i4,a,i1,a)' ) ' rhoij(',iatom,',',ispden,')=  (max 12 non-zero components will be written)'
         call print_ij(pawrhoij(iatom)%rhoijp(:,ispden),&
&         pawrhoij(iatom)%nrhoijsel,&
&         pawrhoij(iatom)%cplex,&
&         pawrhoij(iatom)%lmn_size,2,-1,ibuffer,1,0,&  !TODO should pass unitfi to print_ij. This writes on std_out and ab_out
&        pawrhoij(iatom)%rhoijselect(:),-1.d0,1)
!        &         pawrhoij(iatom)%lmn_size,1,-1,ibuffer,1,0,& ! This write to std_out only.
       end do
     end do
     ABI_DEALLOCATE(ibuffer)

     case default
     MSG_ERROR("Wrong rdwr_mode"//TRIM(rdwr_mode))
 end select

end subroutine rhoij_io
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/rhoij_unpack
!! NAME
!! rhoij_unpack
!!
!! FUNCTION
!!  Unpack the values store in rhoijp copying them to the rhoij_ array.
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is filled with the values stored in the packed array rhoijp.
!!   * If use_rhoij_/=1, rhoij_ is allocated and the corresponding flag is set to 1.
!!
!! PARENTS
!!      paw_qpscgw
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

subroutine rhoij_unpack(rhoij)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_unpack'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: natom,iat,lmn2_size,isel,klmn,nspden,cplex

! *************************************************************************

 natom  = SIZE(rhoij)
 nspden = rhoij(1)%nspden    ! MT jan-2010: this should not be nspden but nsppol or 4 if nspden=4
 cplex  = rhoij(1)%cplex

 do iat=1,natom

   lmn2_size =rhoij(iat)%lmn2_size

   if (rhoij(iat)%use_rhoij_/=1) then ! Have to allocate rhoij_ABI_ALLOCATE(iat,)
     ABI_ALLOCATE(rhoij(iat)%rhoij_,(cplex*lmn2_size,nspden))
     rhoij(iat)%use_rhoij_=1
   end if
   rhoij(iat)%rhoij_ = zero

   do isel=1,rhoij(iat)%nrhoijsel ! Looping over non-zero ij elements.
     klmn = rhoij(iat)%rhoijselect(isel)
     rhoij(iat)%rhoij_(klmn,:) = rhoij(iat)%rhoijp(isel,:)
   end do

 end do ! natom

end subroutine rhoij_unpack
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/rhoij_init_unpacked
!! NAME
!! rhoij_init_unpacked
!!
!! FUNCTION
!!  Initialize field of rhoij datastructure for unpacked values (pawrhoij%rhoij_ array)
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is allocated
!!
!! PARENTS
!!      vtorho,vtorho3
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

subroutine rhoij_init_unpacked(rhoij)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_init_unpacked'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: iat,nrhoij,nsp2

! *************************************************************************

 nrhoij  = SIZE(rhoij)
 nsp2=rhoij(1)%nsppol;if (rhoij(1)%nspden==4) nsp2=4

 do iat=1,nrhoij

   if (associated(rhoij(iat)%rhoij_))  then
     ABI_DEALLOCATE(rhoij(iat)%rhoij_)
   end if
   ABI_ALLOCATE(rhoij(iat)%rhoij_,(rhoij(iat)%cplex*rhoij(iat)%lmn2_size,nsp2))
   rhoij(iat)%use_rhoij_=1
   rhoij(iat)%rhoij_=zero

 end do

end subroutine rhoij_init_unpacked
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/rhoij_destroy_unpacked
!! NAME
!! rhoij_destroy_unpacked
!!
!! FUNCTION
!!  Destroy field of rhoij datastructure for unpacked values (pawrhoij%rhoij_ array)
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is deallocated
!!
!! PARENTS
!!      pawmkrho
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

subroutine rhoij_destroy_unpacked(rhoij)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_destroy_unpacked'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: iat,nrhoij

! *************************************************************************

 nrhoij  = SIZE(rhoij)

 do iat=1,nrhoij

   if (associated(rhoij(iat)%rhoij_))  then
     ABI_DEALLOCATE(rhoij(iat)%rhoij_)
   end if
   nullify(rhoij(iat)%rhoij_)
   rhoij(iat)%use_rhoij_=0

 end do

end subroutine rhoij_destroy_unpacked
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/rhoij_mpi_sum
!! NAME
!! rhoij_mpi_sum
!!
!! FUNCTION
!! Build the MPI sum of the unsymmetrized PAW rhoij_ (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! INPUTS
!!  comm1=MPI communicator. Data will be MPI summed inside comm1
!!  [comm2]=second MPI communicator. If present, rhoij_ will be MPI summed inside comm2 after the collective sum in comm1.
!!
!! SIDE EFFECTS
!!  pawrhoij(:) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  Input: the data calculateed by this processor.
!!  Otput: the final MPI sum over comm1 and comm2.
!!
!! PARENTS
!!      wfd_pawrhoij
!!
!! CHILDREN
!!      timab,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rhoij_mpi_sum(pawrhoij,comm1,comm2)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_mpi_sum'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: comm1
 integer,optional,intent(in) :: comm2
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables ---------------------------------------
!scalars
 integer :: bufdim,iatom,ierr,isppol,jdim,nsp2,natom
 integer :: nproc1,nproc2
 !character(len=500) :: msg
!arrays
 integer,allocatable :: dimlmn(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:)

!************************************************************************

 DBG_ENTER("COLL")

 nproc1 = xcomm_size(comm1)
 nproc2=1; if (PRESENT(comm2)) nproc2 = xcomm_size(comm2)
 if (nproc1==1.and.nproc2==1) RETURN

 call timab(48,1,tsec)

!Fill the MPI buffer from the local rhoij_
 natom = SIZE(pawrhoij)
 ABI_ALLOCATE(dimlmn,(natom))
 dimlmn(1:natom)=pawrhoij(1:natom)%cplex*pawrhoij(1:natom)%lmn2_size
 nsp2=pawrhoij(1)%nsppol; if (pawrhoij(1)%nspden==4) nsp2=4
 bufdim=sum(dimlmn)*nsp2
 ABI_ALLOCATE(buffer1,(bufdim))
 ABI_ALLOCATE(buffer2,(bufdim))
 jdim=0
 do iatom=1,natom
   do isppol=1,nsp2
     buffer1(jdim+1:jdim+dimlmn(iatom))=pawrhoij(iatom)%rhoij_(:,isppol)
     jdim=jdim+dimlmn(iatom)
   end do
 end do
!
!Build sum of  pawrhoij%rhoij_
 call xsum_mpi(buffer1,buffer2,bufdim,comm1,ierr)      ! Sum over the first communicator.
 if (PRESENT(comm2)) call xsum_mpi(buffer2,comm2,ierr) ! Sum over the second communicator.
!
!Unpack the MPI packet filling rhoij_.
 jdim=0
 do iatom=1,natom
   do isppol=1,nsp2
     pawrhoij(iatom)%rhoij_(:,isppol)=buffer2(jdim+1:jdim+dimlmn(iatom))
     jdim=jdim+dimlmn(iatom)
   end do
 end do

 ABI_DEALLOCATE(buffer1)
 ABI_DEALLOCATE(buffer2)
 ABI_DEALLOCATE(dimlmn)

 call timab(48,2,tsec)

 DBG_EXIT("COLL")

end subroutine rhoij_mpi_sum
!!***

