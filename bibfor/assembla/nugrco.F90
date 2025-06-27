! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine nugrco(nu, base)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_jev.h"
#include "asterfort/asmpi_comm_point.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    character(len=14) :: nu
    character(len=2) :: base
! person_in_charge: nicolas.sellenet at edf.fr
! ----------------------------------------------------------------------
!  NUME_DDL - CREATION DU GRAPH DE COMMUNICATION
!  --                     --       --
! ----------------------------------------------------------------------
!
!   ON CREE LE GRAPH DE COMMUNICATION QUI PERMETTRA DANS LES
!    ROUTINES PETSC DE SAVOIR QUEL PROCESSEUR DOIT COMMUNIQUER
!    AVEC QUEL AUTRE
!
! IN  :
!   NU      K14  NOM DU NUME_DDL
!   BASE    K2   BASE(1:1) : BASE POUR CREER LE NUME_DDL
!                    (SAUF LE NUME_EQUA)
!                BASE(2:2) : BASE POUR CREER LE NUME_EQUA
!
    integer(kind=8) :: rang, nbproc, jcomm1, iddl, neql, jgraco
    integer(kind=8) :: iproc, nbedge, iaux, jtmp, nmatch, iproc1
    integer(kind=8) :: iproc2, posit, jordjo, num, neqg, nulodd
    integer(kind=8) :: nbddlj, jjoint, curpos, numpro, jjoin2
    integer(kind=8) :: iddlg, iddll
!
    character(len=4) :: chnbjo
    character(len=24) :: nojoin, nogrco
    integer(kind=8), pointer :: masque(:) => null()
    integer(kind=8), pointer :: posproc(:) => null()
    integer(kind=8), pointer :: nulg(:) => null()
    integer(kind=8), pointer :: nugl(:) => null()
    integer(kind=8), pointer :: nequg(:) => null()
    integer(kind=8), pointer :: nequl(:) => null()
    integer(kind=8), pointer :: pddl(:) => null()
    mpi_int :: mrank, msize
    parameter(nogrco='&&NUGRCO.GRAPH_COMM')
!
    call jemarq()
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    call jeveuo(nu//'.NUML.NUGL', 'L', vi=nugl)
    call jeveuo(nu//'.NUML.NULG', 'L', vi=nulg)
    call jeveuo(nu//'.NUML.PDDL', 'L', vi=pddl)
    call jeveuo(nu//'.NUML.NEQU', 'L', vi=nequl)
    call jeveuo(nu//'.NUME.NEQU', 'L', vi=nequg)
    neql = nequl(1)
    neqg = nequg(1)
    call wkvect('&&NUGRCO.COMM1', 'V V I', nbproc, jcomm1)
!
!---- DETERMINATION DE QUI COMMUNIQUE AVEC QUI
    do iddl = 0, neql-1
        numpro = pddl(iddl+1)
        ASSERT(numpro .lt. nbproc)
        zi(jcomm1+numpro) = zi(jcomm1+numpro)+1
    end do
    zi(jcomm1+rang) = 0
!
    call wkvect(nogrco, 'V V I', nbproc*nbproc, jgraco)
    do iproc = 0, nbproc-1
        if (zi(jcomm1+iproc) .ne. 0) then
            zi(jgraco+iproc+rang*nbproc) = 1
            zi(jgraco+rang+iproc*nbproc) = 1
        end if
    end do
    call asmpi_comm_jev('MPI_SUM', nogrco)
!
!---- RECHERCHE DES COUPLAGES DANS LE GRAPH
    nbedge = 0
    do iaux = 1, nbproc*nbproc
        if (zi(jgraco+iaux-1) .eq. 1) nbedge = nbedge+1
    end do
    nbedge = nbedge/2
!
!---- RECHERCHE DES COUPLAGES MAXIMAUX
    AS_ALLOCATE(vi=masque, size=nbproc*nbproc)
    call wkvect('&&NUGRCO.TMP', 'V V I', nbproc, jtmp)
    nmatch = 1
60  continue
    do iproc1 = 0, nbproc-1
        do iproc2 = 0, nbproc-1
            posit = iproc1*nbproc+iproc2
            if (zi(jgraco+posit) .eq. 1 .and. zi(jtmp+iproc1) .eq. 0 .and. zi(jtmp+iproc2) &
                .eq. 0) then
                zi(jgraco+posit) = 0
                masque(posit+1) = nmatch
                posit = iproc2*nbproc+iproc1
                zi(jgraco+posit) = 0
                masque(posit+1) = nmatch
                nbedge = nbedge-1
                zi(jtmp+iproc1) = 1
                zi(jtmp+iproc2) = 1
            end if
        end do
    end do
    nmatch = nmatch+1
    do iaux = 0, nbproc-1
        zi(jtmp+iaux) = 0
    end do
    if (nbedge .gt. 0) goto 60
    call jedetr('&&NUGRCO.TMP')
!
!---- CREATION DU GRAPH
    nmatch = nmatch-1
    call wkvect(nu//'.NUML.JOIN', base(1:1)//' V I', nmatch, jordjo)
    AS_ALLOCATE(vi=posproc, size=2*nbproc)
    do iaux = 0, nmatch-1
        zi(jordjo+iaux) = -1
    end do
    do iaux = 0, nbproc-1
        num = masque(1+rang*nbproc+iaux)
        ASSERT(num .le. nmatch)
        if (num .ne. 0) then
            zi(jordjo+num-1) = iaux
!
            call codent(num, 'G', chnbjo)
            nojoin = nu//'.NUML.'//chnbjo
            nbddlj = zi(jcomm1+iaux)
            if (nbddlj .ne. 0) then
                call wkvect(nojoin, base(1:1)//' V I', nbddlj, jjoint)
                posproc(1+2*iaux) = jjoint
                posproc(1+2*iaux+1) = 0
            else
                posproc(1+2*iaux) = -1
                posproc(1+2*iaux+1) = -1
            end if
        else
            posproc(1+2*iaux) = -1
            posproc(1+2*iaux+1) = -1
        end if
    end do
!
    do iddl = 0, neqg-1
        nulodd = nugl(iddl+1)
        if (nulodd .ne. 0) then
            numpro = pddl(nulodd)
            if (numpro .ne. rang) then
                jjoint = posproc(1+2*numpro)
                ASSERT(jjoint .ne. -1)
!
                curpos = posproc(1+2*numpro+1)
                zi(jjoint+curpos) = nulodd
                posproc(1+2*numpro+1) = posproc(1+2*numpro+1)+1
            end if
        end if
    end do
!
    do iaux = 0, nmatch-1
        numpro = zi(jordjo+iaux)
        if (numpro .eq. -1) goto 110
!
        if (rang .gt. numpro) then
            nbddlj = posproc(1+2*numpro+1)
            call asmpi_comm_point('MPI_SEND', 'I', numpro, iaux, sci=nbddlj)
!
            jjoint = posproc(1+2*numpro)
            call wkvect('&&NUGRCO.TMP', 'V V I', nbddlj, jjoin2)
            do iddl = 0, nbddlj-1
                iddlg = nulg(1+zi(jjoint+iddl)-1)
                ASSERT(iddlg .ne. 0)
                zi(jjoin2+iddl) = iddlg
            end do
!
            call asmpi_comm_point('MPI_SEND', 'I', numpro, iaux, nbval=nbddlj, &
                                  vi=zi(jjoin2))
            call jedetr('&&NUGRCO.TMP')
        else if (rang .lt. numpro) then
            call asmpi_comm_point('MPI_RECV', 'I', numpro, iaux, sci=nbddlj)
            call wkvect('&&NUGRCO.TMP', 'V V I', nbddlj, jjoin2)
!
            num = iaux+1
            call codent(num, 'G', chnbjo)
            nojoin = nu//'.NUML.'//chnbjo
            call wkvect(nojoin, base(1:1)//' V I', nbddlj, jjoint)
!
            call asmpi_comm_point('MPI_RECV', 'I', numpro, iaux, nbval=nbddlj, &
                                  vi=zi(jjoin2))
            do iddl = 0, nbddlj-1
                iddll = nugl(1+zi(jjoin2+iddl)-1)
                ASSERT(iddll .ne. 0)
                zi(jjoint+iddl) = iddll
            end do
            call jedetr('&&NUGRCO.TMP')
        else
            ASSERT(.false.)
        end if
110     continue
    end do
!
    call jedetr('&&NUGRCO.COMM1')
    call jedetr(nogrco)
    AS_DEALLOCATE(vi=masque)
    AS_DEALLOCATE(vi=posproc)
!
    call jedema()
!
end subroutine
