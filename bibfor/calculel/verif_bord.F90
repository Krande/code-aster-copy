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

subroutine verif_bord(modele, ligrel)
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/liglma.h"
#include "asterfort/dismoi.h"
#include "asterfort/assert.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/utmess.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*), intent(in) :: modele
    character(len=*), intent(in) :: ligrel
!
!-----------------------------------------------------------------------
!   But :
!     Emettre des alarmes si le ligrel ne contient pas toutes les mailles
!     de bord necessaires.
!
!   Entrees:
!     modele     :  sd_modele
!     ligrel     :  sous-ligrel du modele
!
!
!    Algorithme :
!      On parcourt toutes les mailles du modele : ima
!         Si ima n'appartient pas au ligrel
!            Si tous les noeuds de ima sont des noeuds du ligrel => Alarme
!
!-----------------------------------------------------------------------
    character(len=8) :: modele_, noma
    character(len=19) :: ligrel_, ligrmo
    character(len=24) :: valk(4)
    integer(kind=8) :: nbmamo, nbmalg, numa, kma, nbmat
    integer(kind=8) :: iconx1, iconx2, nno, nuno, kno, nbnot

    character(len=24), parameter :: linumamo = '&&VERIF_BORD.NUMAMO'
    character(len=24), parameter :: linutemo = '&&VERIF_BORD.NUTEMO'
    character(len=24), parameter :: linumalg = '&&VERIF_BORD.NUMALG'
    character(len=24), parameter :: linutelg = '&&VERIF_BORD.NUTELG'

    integer(kind=8), pointer :: numamo(:) => null()
    integer(kind=8), pointer :: numalg(:) => null()
    integer(kind=8), pointer :: eximalg(:) => null()
    integer(kind=8), pointer :: exinolg(:) => null()

#define nbno(imail) zi(iconx2+imail) - zi(iconx2+imail-1)
#define connex(imail,j) zi(iconx1-1+zi(iconx2+imail-1)+j-1)

!-----------------------------------------------------------------------
!
    call jemarq()
    modele_ = modele
    ligrmo = modele_//'.MODELE'
    ligrel_ = ligrel

    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbmat)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnot)
    call jeveuo(noma//'.CONNEX', 'L', iconx1)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', iconx2)

    call liglma(ligrmo, nbmamo, linumamo, linutemo)
    call liglma(ligrel_, nbmalg, linumalg, linutelg)
    call jeveuo(linumamo, 'L', vi=numamo)
    call jeveuo(linumalg, 'L', vi=numalg)

!   -- 1. Calcul de eximalg et exinolg :
!      eximalg(numa) = 1 : la maille numa existe dans ligrel
!      exinolg(nuno) = 1 : le noeud numo existe dans ligrel
!   ----------------------------------------------------------
    AS_ALLOCATE(vi=eximalg, size=nbmat)
    AS_ALLOCATE(vi=exinolg, size=nbnot)
    eximalg = 0
    exinolg = 0
    do kma = 1, nbmalg
        numa = numalg(kma)
        eximalg(numa) = 1
        nno = nbno(numa)
        do kno = 1, nno
            nuno = connex(numa, kno)
            ASSERT(nuno .gt. 0 .and. nuno .le. nbnot)
            exinolg(nuno) = 1
        end do
    end do

!   -- 2. boucle sur les mailles de modele :
!   ----------------------------------------
    B1: do kma = 1, nbmamo
        numa = numamo(kma)
        if (eximalg(numa) .eq. 1) cycle B1

        nno = nbno(numa)
        do kno = 1, nno
            nuno = connex(numa, kno)
            if (exinolg(nuno) .eq. 0) cycle B1
        end do
        valk(1) = modele
        valk(2) = int_to_char8(numa)
        call utmess('A', 'CALCULEL4_74', nk=2, valk=valk)
    end do B1

!   -- menage :
!   -----------
    call jedetr(linumamo)
    call jedetr(linutemo)
    call jedetr(linumalg)
    call jedetr(linutelg)
    AS_DEALLOCATE(vi=eximalg)
    AS_DEALLOCATE(vi=exinolg)

    call jedema()
end subroutine
