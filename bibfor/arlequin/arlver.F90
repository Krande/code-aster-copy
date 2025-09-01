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
subroutine arlver(modele, lgma, nbgma, nomsd, model, &
                  cine)
!
!
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/arlelt.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/tri.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: modele
    character(len=8) :: lgma(*)
    integer(kind=8) :: nbgma
    character(len=10) :: nomsd
    character(len=8) :: model
    character(len=8) :: cine
!
! ----------------------------------------------------------------------
!
! ROUTINE ARLEQUIN
!
! FILTRE, REGROUPEMENT, VERIFICATION GROUPES DE MAILLES
!
! ----------------------------------------------------------------------
!
!
! IN  MODELE : NOM DU MODELE
! IN  LGMA   : NOMS DES GROUPES DE MAILLES
! I/O NBGMA  : EN ENTREE -> NOMBRE DE GROUPES DE MAILLES
!              EN SORTIE -> NOMBRE DE MAILLES DANS .GROUPEMA
! I/O NOMSD  : SD DOMAINE
! OUT MODEL  : MODELISATION ASSOCIEE AUX MAILLES ('3D')
! OUT CINE   : CINEMATIQUE ASSOCIEE AUX MAILLES
!              'SOLIDE' OU 'POUTRE'
!               UN POUR CHAQUE GROUPE
!
! SD DE SORTIE (NOMSD) :
! ======================
! .GROUPEMA : LISTE TRIEE DES MAILLES DU DOMAINE NOMSD (MA1,MA2,...)
!
!
    character(len=8) :: noma, k8bid
    character(len=16) :: nomte, modte, cinte
    character(len=19) :: ligrel
    integer(kind=8) :: nbma, ntot, nbligr, numa, ninit
    integer(kind=8) :: icompt, igma, ima, iligr, iret
    integer(kind=8) :: jrepe, jcompt, jte, jgma, jgroup, jtyel
    aster_logical :: eltok
    integer(kind=8) :: liste(nbgma)
!
    character(len=6) :: nompro
    parameter(nompro='ARLVER')
!
! ---------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    ninit = 0
    cinte = ' '
    modte = ' '
!
! --- ACCES MODELE
!
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrel)
    call jeveuo(ligrel//'.REPE', 'L', jrepe)
    call jelira(ligrel//'.LIEL', 'NMAXOC', nbligr, k8bid)
    call jeveuo(ligrel//'.TYFE', 'L', jtyel)
!
! --- NOM DU MAILLAGE
!
    call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=noma)
!
! --- ALLOCATION OBJETS TEMPORAIRES
!
    call wkvect('&&'//nompro//'.COMPTEUR', 'V V I', nbligr, jcompt)
    call wkvect('&&'//nompro//'.TE', 'V V I', nbligr, jte)
!
! --- COMPTER LE NOMBRE DE MAILLES AFFECTEES A CHAQUE TYPE D'ELEMENT
! --- DONNER LE NUMERO ABSOLU DU NUMERO D'ELEMENT AFFECTE A CHAQUE TYPE
!
    do igma = 1, nbgma
        call jeexin(jexnom(noma(1:8)//'.GROUPEMA', lgma(igma)), iret)
        if (iret == 0) then
            ASSERT(.false.)
        end if
        call jeveuo(jexnom(noma(1:8)//'.GROUPEMA', lgma(igma)), 'L', jgma)
        call jelira(jexnom(noma(1:8)//'.GROUPEMA', lgma(igma)), 'LONMAX', nbma, k8bid)
        ninit = ninit+nbma
        do ima = 1, nbma
            numa = zi(jgma-1+ima)
            iligr = zi(jrepe+2*(numa-1))
            if (iligr /= 0) then
                if (zi(jcompt-1+iligr) == 0) then
                    zi(jte-1+iligr) = zi(jtyel-1+numa)
                end if
                zi(jcompt-1+iligr) = zi(jcompt-1+iligr)+1
            end if
        end do
    end do
!
! --- VERIFICATION COHERENCE MODELISATION / CINEMATIQUE
!
    ntot = 0
    do iligr = 1, nbligr
        if (zi(jcompt-1+iligr) /= 0) then
            call jenuno(jexnum('&CATA.TE.NOMTE', zi(jte-1+iligr)), nomte)
            eltok = arlelt(nomte, modte, cinte)
            if (.not. eltok) then
                zi(jcompt-1+iligr) = 0
                goto 30
            else
                if (ntot == 0) then
                    model = modte(1:8)
                    cine = cinte(1:8)
                else
                    if (modte /= model) then
                        ASSERT(.false.)
                    end if
                    if (cinte /= cine) then
                        ASSERT(.false.)
                    end if
                end if
                ntot = ntot+zi(jcompt-1+iligr)
            end if
        end if
30      continue
    end do
!
! --- AUCUNE MAILLE DU GROUPE N'EST UTILISABLE
!
    if (ntot == 0) then
        ASSERT(.false.)
    end if
!
! --- ALLOCATION .GROUPEMA
!
    call wkvect(nomsd//'.GROUPEMA', 'V V I', ntot, jgroup)
!
! --- COPIE DES MAILLES VALIDES
!
    icompt = 0
    do igma = 1, nbgma
        call jeveuo(jexnom(noma//'.GROUPEMA', lgma(igma)), 'L', jgma)
        call jelira(jexnom(noma//'.GROUPEMA', lgma(igma)), 'LONMAX', nbma, k8bid)
        do ima = 1, nbma
            numa = zi(jgma-1+ima)
            iligr = zi(jrepe+2*(numa-1))
            if ((iligr /= 0) .and. (zi(jcompt-1+iligr) /= 0)) then
                zi(jgroup+icompt) = numa
                icompt = icompt+1
            end if
        end do
    end do
!
    if (icompt /= ntot) then
        ASSERT(.false.)
    end if
    nbgma = ntot
    call tri(zi(jgroup), liste, 0, nbgma)
!
! --- MENAGE
!
    call jedetr('&&'//nompro//'.COMPTEUR')
    call jedetr('&&'//nompro//'.TE')
!
    call jedema()
!
end subroutine
