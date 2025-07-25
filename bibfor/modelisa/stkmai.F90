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

subroutine stkmai(ifl, icl, iv, rv, cv, &
                  cnl, mcl, nbm, nume, numn, &
                  cnx, typ, fmt, irteti, nommai)
    implicit none
!       SECONDE LECTURE DES DONNEES POUR UN MOT CLE DE TYPE MAILLE
!       ----------------------------------------------------------------
!       IN      IFL,ICL,IV,RV,CV,CNL = VOIR LIRITM
!               MCL             = MOTS CLES TYPE MAILLE
!               NBM             = NB DE MOTS CLES TYPE MAILLE
!               FMT             = NB NOEUDS A LIRE / MAILLE
!               CNX             = NOMU.CONXV
!               TYP             = NOMU.TYPMAIL
!               NUME            = NUMERO DE L ELEMENT COURANT
!               NUMN            = NUMERO DU NOEUD COURANT DANS CNX
!       OUT     (RETURN)        = MOT CLE SUIVANT (MOT CLE NON RECONNU)
!               (RETURN 1)      = EXIT            (MOT CLE FIN TROUVE)
!               (RETURN 2)      = LIGNE SUIVANTE  (MOT CLE FINSF TROUVE
!                                                  OU ERREUR DETECTE)
!       ----------------------------------------------------------------
!
#include "jeveux.h"
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
    integer(kind=8) :: nbm
    character(len=8) :: mcl(nbm), noma, b8
    integer(kind=8) :: deblig, fmt(nbm)
    character(len=14) :: cnl
    character(len=*) :: cv
    character(len=24) :: cnx, typ, nom, nommai
    save b8
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadc, iadt, icl, ifl, iret
    integer(kind=8) :: irtet, irteti, iv, nume, numn
    integer(kind=8) :: numtcl
!-----------------------------------------------------------------------
    data b8/'        '/
    call jemarq()
    irteti = 0
!
!
! - ITEM = MOT CLE  TYPE MAILLE ?
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
!
5   continue
    call jeveuo(cnx, 'E', iadc)
    call jeveuo(typ, 'E', iadt)
!
! - LECTURE DE L'ENTETE
!
    deblig = 0
    call lirtet(ifl, 2, 0, cnl, nom, &
                icl, iv, rv, cv, deblig)
    goto 9
!
! - LIRE ITEM SUIVANT = NOM DE MAILLE ?
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
!
! - CREATION DE CONXV.NOM_DE_MAILLE ET TYPMAIL.NOM_DE_MAILLE
!
    noma = b8
    noma(1:iv) = cv(1:iv)
    call jeexin(jexnom(nommai, noma), iret)
    if (iret .eq. 0) then
        call jecroc(jexnom(nommai, noma))
        call jecroc(jexnom(cnx, noma))
        call jeecra(jexnom(cnx, noma), 'LONMAX', fmt(numtcl))
    else
        call utmess('F', 'MODELISA7_10', sk=noma)
    end if
!
! - STOCKAGE DES NOMS DES NOEUDS DE LA MAILLE ET DU TYPE DE MAILLE
!
    zi(iadt+nume) = numtcl
!
    do i = 1, fmt(numtcl)
        call liritm(ifl, icl, iv, rv, cv, &
                    cnl, deblig, 2)
        nom = b8
        nom(1:iv) = cv(1:iv)
!
        zk8(iadc+numn) = nom(1:8)
!
! - INCREMENTATION DU NB DE NOEUDS LUS
!
        numn = numn+1
    end do
!
! - INCREMENTATION DU NB D ELEMENTS LUS
!
    nume = nume+1
!
! - MAILLE SUIVANTE
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
