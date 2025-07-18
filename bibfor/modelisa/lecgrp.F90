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

subroutine lecgrp(ifl, icl, iv, rv, cv, &
                  cnl, mcl, nbm, nbg, dim, &
                  nbt, ier, irteti)
    implicit none
!       PREMIERE LECTURE DES DONNEES POUR UN MOT CLE DE TYPE GROUPE
!       ----------------------------------------------------------------
!       IN      IFL,ICL,IV,RV,CV,CNL = VOIR LIRITM
!               MCL             = MOTS CLE TYPE GROUPE
!               NBG             = NIVEAU DEBUG
!               NBM             = NB DE MOTS CLES TYPE GROUPE
!       OUT     IER             = 0 > LECTURE CORRECTE
!                               = 1 > ERREUR EN LECTURE
!               DIM             = DIMENSIONS DES OBJETS GROUPE..(NB GR.)
!               NBT             = NB TOTAL DE NOEUDS/MAILLES LUS
!               (RETURN)        = MOT CLE SUIVANT (MOT CLE NON RECONNU)
!               (RETURN 1)      = EXIT            (MOT CLE FIN TROUVE)
!               (RETURN 2)      = LIGNE SUIVANTE  (MOT CLE FINSF TROUVE
!                                                  OU ERREUR DETECTE)
!       ----------------------------------------------------------------
!
#include "asterfort/iunifi.h"
#include "asterfort/liritm.h"
#include "asterfort/lirtet.h"
#include "asterfort/tesfin.h"
#include "asterfort/tesmcl.h"
#include "asterfort/utmess.h"
#include "asterfort/verdbl.h"
#include "asterfort/vermot.h"
    real(kind=8) :: rv
    integer(kind=8) :: nbm
    character(len=8) :: mcl(nbm)
    integer(kind=8) :: dim(nbm), nbt(nbm), deblig
    character(len=14) :: cnl
    character(len=24) :: valk(2), nom
    character(len=*) :: cv
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icl, ier, ifl, ifm, ifn, ilec
    integer(kind=8) :: inom, irtet, irteti, iv, nbg
    integer(kind=8) :: nbtav, numtcl
!-----------------------------------------------------------------------
    irteti = 0
!
    ifm = iunifi('MESSAGE')
!
! ----- ITEM = MOT CLE TYPE  GROUPE ?
!
    do i = 1, nbm
        call tesmcl(icl, iv, cv, mcl(i), irtet)
        if (irtet .eq. 1) goto 4
        numtcl = i
        nbtav = nbt(numtcl)
        goto 5
4       continue
    end do
    goto 3
!
! ----- LIRE ITEM SUIVANT ( = NOM DU GROUPE ?)
!
5   continue
    inom = 1
    ilec = 1
    deblig = 0
    call lirtet(ifl, ilec, inom, cnl, nom, &
                icl, iv, rv, cv, deblig)
!
    if (nom .eq. 'INDEFINI') then
!
! -----    LECTURE NOM DU GROUPE SI IL N Y A PAS D'ENTETE
!
        if (nbg .ge. 1) write (ifm, *) '       LIRITM : ICL = ', icl, ' IV = ', iv, ' RV = ', rv, &
            ' CV(1:24) = ', cv(1:24), ' DEBLIG =', deblig
        call verdbl(deblig, cnl, ier, irtet)
        if (irtet .eq. 1) goto 2
        call tesfin(icl, iv, cv, irtet)
        if (irtet .eq. 1) then
            goto 7
        else if (irtet .eq. 2) then
            goto 8
        end if
        call vermot(icl, iv, cv, cnl, ier, &
                    irtet)
        if (irtet .eq. 1) goto 2
        nom = cv(1:iv)
    else
!
! -----    LECTURE PREMIER NOM DE NOEUD / MAILLE OU FIN APRES L'ENTETE
!
        call tesfin(icl, iv, cv, irtet)
        if (irtet .eq. 1) then
            goto 7
        else if (irtet .eq. 2) then
            goto 8
        end if
        call vermot(icl, iv, cv, cnl, ier, &
                    irtet)
        if (irtet .eq. 1) goto 2
        nbt(numtcl) = nbt(numtcl)+1
    end if
!
! ----- LECTURE DES NOMS DE NOEUDS OU MAILLES DU GROUPE
!
6   continue
    call liritm(ifl, icl, iv, rv, cv, &
                cnl, deblig, 1)
    if (nbg .ge. 1) write (ifm, *) '       LIRITM : ICL = ', icl, ' IV = ', iv, ' RV = ', rv, &
        ' CV(1:24) = ', cv(1:24), ' DEBLIG =', deblig
!
    if (deblig .eq. 1) call tesfin(icl, iv, cv, irtet)
    if (irtet .eq. 1) then
        goto 7
    else if (irtet .eq. 2) then
        goto 8
    end if
!
    call vermot(icl, iv, cv, cnl, ier, &
                irtet)
    if (irtet .eq. 1) goto 2
!
    nbt(numtcl) = nbt(numtcl)+1
    goto 6
!
7   continue
    ifn = 0
    goto 9
8   continue
    ifn = 1
!
9   continue
    if ((nbtav-nbt(numtcl)) .eq. 0) then
        valk(1) = cnl
        valk(2) = nom
        call utmess('A', 'MODELISA4_80', nk=2, valk=valk)
!         -- ON VA CREER UN GROUPE VIDE DE LONGUEUR 1 :
        dim(numtcl) = dim(numtcl)+1
        nbt(numtcl) = nbt(numtcl)+1
        goto 2
    end if
!
    dim(numtcl) = dim(numtcl)+1
    if (ifn .eq. 0) goto 1
    if (ifn .eq. 1) goto 2
!
1   continue
    irteti = 1
    goto 999
2   continue
    irteti = 2
    goto 999
3   continue
    irteti = 0
    goto 999
!
999 continue
end subroutine
