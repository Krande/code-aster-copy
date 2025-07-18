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

subroutine stkgrp(ifl, icl, iv, rv, cv, &
                  cnl, mcl, nbm, numn, numm, &
                  grn, grm, irteti)
    implicit none
!       SECONDE LECTURE DES DONNEES POUR UN MOT CLE DE TYPE GROUPE
!       ----------------------------------------------------------------
!       IN      IFL,ICL,IV,RV,CV,CNL = VOIR LIRITM
!               MCL             = MOTS CLES TYPE GROUPE
!               NBM             = NB DE MOTS CLES TYPE GROUPE
!               GRN             = NOMU.GROUPNOV
!               GRM             = NOMU.GROUPMAV
!               NUMN            = NUMERO DU NOEUD COURANT DANS GRN
!               NUMM            = NUMERO DE MAILLE COURANT DANS GRM
!       OUT     (RETURN)        = MOT CLE SUIVANT (MOT CLE NON RECONNU)
!               (RETURN 1)      = EXIT            (MOT CLE FIN TROUVE)
!               (RETURN 2)      = LIGNE SUIVANTE  (MOT CLE FINSF TROUVE
!                                                  OU ERREUR DETECTE)
!       ----------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/liritm.h"
#include "asterfort/lirtet.h"
#include "asterfort/tesfin.h"
#include "asterfort/tesmcl.h"
#include "asterfort/utmess.h"
!
    real(kind=8) :: rv
    integer(kind=8) :: deblig, nbm
    character(len=8) :: mcl(nbm), nom, b8
    character(len=14) :: cnl
    character(len=*) :: cv
    character(len=24) :: grn, grm, grp, nomg
    save b8
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadg, icl, ifl, ifn, iret
    integer(kind=8) :: irtet, irteti, iv, nbitem, num
    integer(kind=8) :: numm, numn
!-----------------------------------------------------------------------
    data b8/'        '/
!
    call jemarq()
    irteti = 0
!
!
! - ITEM = MOT CLE TYPE GROUPE ?
!
    do i = 1, nbm
        call tesmcl(icl, iv, cv, mcl(i), irtet)
        if (irtet .eq. 1) goto 4
        if (i .eq. 1) then
            grp = grn
            num = numn
        else
            grp = grm
            num = numm
        end if
        goto 10
4       continue
    end do
    goto 3
!
10  continue
    call jeveuo(grp, 'E', iadg)
!
! ----- LIRE ITEM SUIVANT =  NOM DU GROUPE ?
    deblig = 0
    call lirtet(ifl, 2, 1, cnl, nomg, &
                icl, iv, rv, cv, deblig)
!
! ----- LECTURE NOM DU GROUPE SI IL N Y A PAS D'ENTETE
    if (nomg .eq. 'INDEFINI') then
        nomg = ' '
        nomg(1:iv) = cv(1:iv)
        call tesfin(icl, iv, cv, irtet)
        ASSERT(irtet .eq. 0)
        if (irtet .eq. 1) then
            goto 7
        else if (irtet .eq. 2) then
            goto 8
        end if
    else
!
! -----   STOCKAGE PREMIER NOM DE NOEUD / MAILLE OU FIN APRES L'ENTETE
        call tesfin(icl, iv, cv, irtet)
        if (irtet .eq. 1) then
            goto 7
        else if (irtet .eq. 2) then
            goto 8
        end if
        nom = b8
        nom(1:iv) = cv(1:iv)
        zk8(iadg+num) = nom
        num = num+1
    end if
!
! ----- STOCKAGE DES NOMS DES NOEUDS OU MAILLES DU GROUPE
6   continue
!
    call liritm(ifl, icl, iv, rv, cv, &
                cnl, deblig, 2)
!
!
! ----- ITEM = MOT  CLE FIN  OU FINSF ?
    if (deblig .eq. 1) then
        call tesfin(icl, iv, cv, irtet)
        if (irtet .eq. 1) then
            goto 7
        else if (irtet .eq. 2) then
            goto 8
        end if
    end if
!
! ----- STOCKAGE DES NOEUDS OU MAILLES DU GROUPE
    nom = b8
    nom(1:iv) = cv(1:iv)
    zk8(iadg+num) = nom
!
! ----- INCREMENTATION DU NB D'ITEM LUS
    num = num+1
!
! ----- NOEUD OU MAILLE SUIVANT
    goto 6
!
! ----- SORTIE EN FIN OU FINSF
7   continue
    ifn = 0
    goto 9
8   continue
    ifn = 1
!
!
! ----- CREATION ET DIMENSIONNEMENT DE L OBJET GRP.NOM_DU_GROUPE
9   continue
!
    if (i .eq. 1) then
        nbitem = num-numn
        numn = num
!         -- POUR UN GROUPE VIDE, LONMAX=1
        if (nbitem .eq. 0) numn = numn+1
    else
        nbitem = num-numm
        numm = num
!         -- POUR UN GROUPE VIDE, LONMAX=1
        if (nbitem .eq. 0) numm = numm+1
    end if
!
    ASSERT(nbitem .ge. 0)
    call jeexin(jexnom(grp, nomg), iret)
    if (iret .eq. 0) then
        call jecroc(jexnom(grp, nomg))
        call jeecra(jexnom(grp, nomg), 'LONMAX', max(nbitem, 1))
        call jeecra(jexnom(grp, nomg), 'LONUTI', nbitem)
    else
        call utmess('F', 'MODELISA7_11', sk=nomg)
    end if
!
    if (ifn .eq. 0) goto 1
    if (ifn .eq. 1) goto 2
!
1   continue
    irteti = 1
    goto 999
!
!       FINSF
2   continue
    irteti = 2
    goto 999
!
3   continue
    irteti = 0
    goto 999
!
999 continue
    call jedema()
end subroutine
