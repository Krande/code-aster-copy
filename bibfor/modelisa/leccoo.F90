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

subroutine leccoo(ifl, icl, iv, rv, cv, &
                  cnl, mcl, nbm, nbg, dim, &
                  nbt, ier, irteti)
    implicit none
!       PREMIERE LECTURE DES DONNEES POUR UN MOT CLE DE TYPE COORDONNEE
!       ----------------------------------------------------------------
!       IN      IFL,ICL,IV,RV,CV,CNL = VOIR LIRITM
!               MCL             = MOTS CLE TYPE COORDONNEE
!               NBG             = NIVEAU DEBUG
!               NBM             = NB DE MOTS CLES TYPE COORDONNEE
!       OUT     IER             = 0 > LECTURE CORRECTE
!                               = 1 > ERREUR EN LECTURE
!               DIM             = DIMENSION DE L OBJET COORDO (NB NOEU.)
!               NBT             = NB TOTAL DE COORDONNEES LUES
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
#include "asterfort/vernmb.h"
    integer(kind=8) :: nbm
    real(kind=8) :: rv
    character(len=8) :: mcl(nbm)
    integer(kind=8) :: dim(nbm), nbt(nbm), deblig
    character(len=14) :: cnl
    character(len=*) :: cv
    character(len=24) :: valk(3), nom
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icl, ier, ifl, ifm, ilec, inom
    integer(kind=8) :: irtet, irteti, iv, j, nbg
    integer(kind=8) :: numtcl
!-----------------------------------------------------------------------
    irteti = 0
!
    ifm = iunifi('MESSAGE')
!
! ----- ITEM = MOT CLE  TYPE COORDONNEES ?
!
    do i = 1, nbm
        call tesmcl(icl, iv, cv, mcl(i), irtet)
        if (irtet .eq. 1) goto 4
        numtcl = i
        goto 5
4       continue
    end do
    goto 3
!
! - VERIFICATION DE COMPATIBILITE DES DECLARATIONS DE DIMENSIONS
!
5   continue
    if (nbg .ge. 1) write (ifm, *) ' ----- LECCOO'
    do j = 1, nbm
        if (dim(j) .ne. 0 .and. j .ne. i) then
            valk(1) = cnl
            valk(2) = mcl(i)
            valk(3) = mcl(j)
            call utmess('E', 'MODELISA4_77', nk=3, valk=valk)
            ier = 1
            goto 2
        end if
    end do
!
! - LECTURE DE L'ENTETE
!
    inom = 0
    ilec = 1
    deblig = 0
    call lirtet(ifl, ilec, inom, cnl, nom, &
                icl, iv, rv, cv, deblig)
    goto 9
!
! - LIRE ITEM SUIVANT
!
7   continue
    call liritm(ifl, icl, iv, rv, cv, &
                cnl, deblig, 1)
    if (nbg .ge. 1) then
        write (ifm, *) '       LIRITM : ICL = ', icl, ' IV = ', iv, ' RV = ', rv, &
            ' CV(1:8) = ', cv(1:8), ' DEBLIG =', deblig
        flush (ifm)
    end if
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
! - ITEM = MOT ? > LECTURE NOM DU NOEUD
!
    call vermot(icl, iv, cv, cnl, ier, &
                irtet)
    if (irtet .eq. 1) goto 2
!
! - LECTURE DES COORDONNEES
!
    do i = 1, numtcl
!
! - LIRE ITEM SUIVANT
!
        call liritm(ifl, icl, iv, rv, cv, &
                    cnl, deblig, 1)
        if (nbg .ge. 1) then
            write (ifm, *) '       LIRITM : ICL = ', icl, ' IV = ', iv, ' RV = ', rv, &
                ' CV(1:8) = ', cv(1:8), ' DEBLIG =', deblig
            flush (ifm)
        end if
!
! - ITEM = NOMBRE ?
!
        call vernmb(icl, cnl, ier, irtet)
        if (irtet .eq. 1) goto 2
!
! - INCREMENTATION DU NB TOTAL DE COORDONNEES LUES
!
        nbt(numtcl) = nbt(numtcl)+1
!
    end do
    dim(numtcl) = dim(numtcl)+1
!
    goto 7
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
