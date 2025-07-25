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

subroutine lecmai(ifl, icl, iv, rv, cv, &
                  cnl, mcl, nbm, nbg, fmt, &
                  dim, nbt, ier, irteti)
    implicit none
!       PREMIERE LECTURE DES DONNEES POUR UN MOT CLE DE TYPE MAILLE
!       ----------------------------------------------------------------
!       IN      IFL,ICL,IV,RV,CV,CNL = VOIR LIRITM
!               MCL             = MOTS CLE TYPE MAILLE
!               NBM             = NB DE MOTS CLES TYPE MAILLE
!               NBG             = NIVEAU DEBUG
!               FMT             = NB NOEUDS A LIRE PAR MAILLE
!       OUT     IER             = 0 > LECTURE CORRECTE
!                               = 1 > ERREUR EN LECTURE
!               DIM             = DIMENSIONS DE L OBJET CONNEX (NB MAIL)
!               NBT             = NB TOTAL DE NOEUDS LUS
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
#include "asterfort/verdbl.h"
#include "asterfort/vermot.h"
    integer(kind=8) :: nbm
    real(kind=8) :: rv
    character(len=8) :: mcl(nbm)
    integer(kind=8) :: dim(nbm), nbt(nbm), deblig
    character(len=14) :: cnl
    character(len=24) :: nom
    character(len=*) :: cv
    integer(kind=8) :: fmt(nbm)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icl, ier, ifl, ifm, ilec, inom
    integer(kind=8) :: irtet, irteti, iv, nbg, numtcl
!
!-----------------------------------------------------------------------
    irteti = 0
!
    ifm = iunifi('MESSAGE')
!
! ----- ITEM = MOT CLE TYPE  MAILLE ?
!
    do i = 1, nbm
        call tesmcl(icl, iv, cv, mcl(i), irtet)
        if (irtet .eq. 1) goto 4
        numtcl = i
        goto 7
4       continue
    end do
    goto 3
!
! - LIRE ITEM SUIVANT
!
7   continue
    if (nbg .ge. 1) write (ifm, *) ' ----- LECMAI'
!
! - LECTURE DE L'ENTETE
!
    inom = 0
    ilec = 1
    deblig = 0
    call lirtet(ifl, ilec, inom, cnl, nom, &
                icl, iv, rv, cv, deblig)
    goto 9
5   continue
    call liritm(ifl, icl, iv, rv, cv, &
                cnl, deblig, 1)
    if (nbg .ge. 1) write (ifm, *) '       LIRITM : ICL = ', icl, ' IV = ', iv, ' RV = ', rv, &
        ' CV(1:8) = ', cv(1:8), ' DEBLIG =', deblig
9   continue
!
! - ITEM EN DEBUT DE LIGNE ?
!
    call verdbl(deblig, cnl, ier, irtet)
    if (irtet .eq. 1) goto 2
!
! - ITEM = MOT  CLE FIN  OU FINSF ?
!
    call tesfin(icl, iv, cv, irtet)
    if (irtet .eq. 1) then
        goto 1
    else if (irtet .eq. 2) then
        goto 2
    end if
!
! - ITEM = MOT (# FIN OU FINSF) = NOM DE LA MAILLE ?
!
    call vermot(icl, iv, cv, cnl, ier, &
                irtet)
    if (irtet .eq. 1) goto 2
!
! ---- LECTURE DES NOMS DES NOEUDS DE LA MAILLE
!
    do i = 1, fmt(numtcl)
!
! - LIRE ITEM SUIVANT
!
        call liritm(ifl, icl, iv, rv, cv, &
                    cnl, deblig, 1)
        if (nbg .ge. 1) write (ifm, *) '       LIRITM : ICL = ', icl, ' IV = ', iv, ' RV = ', rv, &
            ' CV(1:8) = ', cv(1:8), ' DEBLIG =', deblig
!
! - ITEM = MOT (# FIN OU FINSF) ?
!
        call vermot(icl, iv, cv, cnl, ier, &
                    irtet)
        if (irtet .eq. 1) goto 2
!
! - INCREMENTATION DU NB TOTAL DE NOEUDS LUS
!
        nbt(numtcl) = nbt(numtcl)+1
    end do
!
    dim(numtcl) = dim(numtcl)+1
    goto 5
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
