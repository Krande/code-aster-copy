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

subroutine stkcoo(ifl, icl, iv, rv, cv, &
                  cnl, mcl, nbm, num, coo, &
                  nno, irteti)
    implicit none
!       SECONDE LECTURE DES DONNEES POUR UN MOT CLE DE TYPE COORDONNEE
!       ----------------------------------------------------------------
!       IN      IFL,ICL,IV,RV,CV,CNL = VOIR LIRITM
!               MCL             = MOTS CLES TYPE COORDONNEE
!               NBM             = NB DE MOTS CLES TYPE COORDONNEE
!               COO             = NOMU.COORDO.VALE
!               NNO             = NOMU.NOMNOE
!               NUM             = NUMERO DU NOEUD COURANT
!       OUT     (RETURN)        = MOT CLE SUIVANT (MOT CLE NON RECONNU)
!               (RETURN 1)      = EXIT            (MOT CLE FIN TROUVE)
!               (RETURN 2)      = LIGNE SUIVANTE  (MOT CLE FINSF TROUVE
!                                                  OU ERREUR DETECTE)
!       ----------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
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
    integer(kind=8) :: deblig
    integer(kind=8) :: nbm
    real(kind=8) :: rv
    character(len=8) :: mcl(nbm), nomn
    character(len=14) :: cnl
    character(len=*) :: cv
    character(len=24) :: coo, nno, nom
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, icl, idec, ifl, iret
    integer(kind=8) :: irtet, irteti, iv, num, numtcl
!
!-----------------------------------------------------------------------
    call jemarq()
    irteti = 0
!
!
! - ITEM = MOTS CLES  TYPE COORDONNEES ?
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
5   continue
    call jeveuo(coo, 'E', iad)
!
! - LECTURE DE L'ENTETE
!
    deblig = 0
    call lirtet(ifl, 2, 0, cnl, nom, &
                icl, iv, rv, cv, deblig)
    goto 9
!
! - LIRE ITEM SUIVANT =  NOM DU NOEUD ?
!
7   continue
    call liritm(ifl, icl, iv, rv, cv, &
                cnl, deblig, 2)
9   continue
!
! - ITEM = MOT  CLE FIN  OU FINSF ?
!
    call tesfin(icl, iv, cv, irtet)
    if (irtet .eq. 1) then
        goto 1
    else if (irtet .eq. 2) then
        goto 2
    end if
! - CREATION DE NOM_DU_NOEUD DANS LE REPERTOIRE NOMNOE
!
    nomn = '        '
    nomn(1:iv) = cv(1:iv)
    call jeexin(jexnom(nno, nomn), iret)
    if (iret .eq. 0) then
        call jecroc(jexnom(nno, nomn))
    else
        call utmess('F', 'MODELISA7_10', sk=nomn)
    end if
!
!
! - INCREMENTATION NUMERO DU NOEUD
!
    num = num+1
    idec = iad+(num-1)*3
!       IDEC = IAD + (NUM-1) * NUMTCL
!
! - STOCKAGE DES  COORDONNEES DU NOEUD
!
    do i = 1, 3
        zr(idec+i-1) = 0.d0
    end do
!
    do i = 1, numtcl
        call liritm(ifl, icl, iv, rv, cv, &
                    cnl, deblig, 2)
        if (icl .eq. 1) rv = iv
        zr(idec+i-1) = rv
    end do
!
! - NOEUD SUIVANT
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
    call jedema()
end subroutine
