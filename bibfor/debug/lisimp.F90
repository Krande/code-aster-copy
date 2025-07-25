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

subroutine lisimp(lischa, ifm)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/liscpp.h"
#include "asterfort/lisdef.h"
#include "asterfort/lislch.h"
#include "asterfort/lislco.h"
#include "asterfort/lisllc.h"
#include "asterfort/lislnf.h"
#include "asterfort/lislta.h"
#include "asterfort/lisltc.h"
#include "asterfort/lisltf.h"
#include "asterfort/lisnnb.h"
    character(len=19) :: lischa
    integer(kind=8) :: ifm
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (LISTE_CHARGES)
!
! IMPRESSION DU CONTENU DE LA SD LISTE_CHARGES
!
! ----------------------------------------------------------------------
!
!
! IN  LISCHA : NOM DE LA SD LISTE_CHARGES
! IN  IFM    : NUMERO UNITE LOGIQUE IMPRESSION
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ichar, nbchar, ibid
    character(len=8) :: charge, typech, nomfct, k8bid
    character(len=16) :: typapp, typfct
    integer(kind=8) :: genrec(1), tabcod(30)
    character(len=24) :: lisgen, nomlis, gencha
    integer(kind=8) :: jlisg, nbgenr(2), igenr, iposit(2)
    character(len=13) :: prefob
    real(kind=8) :: phase
    integer(kind=8) :: npuis
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    call lisnnb(lischa, nbchar)
    if (nbchar .eq. 0) then
        write (ifm, *) '<LISCHA> PAS DE CHARGES'
        goto 99
    else
        write (ifm, *) '<LISCHA> NOMBRE DE CHARGES: ', nbchar
    end if
!
! --- LISTE DES GENRES DISPONIBLES
!
    lisgen = '&&LISIMP.LISGEN'
    call lisdef('LISG', lisgen, ibid, k8bid, nbgenr)
    call jeveuo(lisgen, 'L', jlisg)
!
! --- AFFICHAGE
!
    do ichar = 1, nbchar
!
! ----- INFORMATIONS GENERALES
!
        call lislch(lischa, ichar, charge)
        call lisltc(lischa, ichar, typech)
        call lislta(lischa, ichar, typapp)
        call lislco(lischa, ichar, genrec(1))
        call lisltf(lischa, ichar, typfct)
        call lisllc(lischa, ichar, prefob)
        write (6, *) 'CHARGE NUMERO : ', ichar
        write (6, *) '  * NOM DE LA CHARGE                  : ', charge
        write (6, *) '  * TYPE DE LA CHARGE                 : ', typech
        write (6, *) '  * TYPE D APPLICATION                : ', typapp
        write (6, *) '  * CODE DE LA CHARGE                 : ', genrec(1)
        write (6, *) '  * PREFIXE DE L''OBJET DE LA CHARGE   : ', prefob
        write (6, *) '  * FONCTION MULTIPLICATRICE:'
        write (6, *) '  ** TYPE                : ', typfct
        if (typfct(1:5) .eq. 'FONCT') then
            call lislnf(lischa, ichar, nomfct)
            write (6, *) '  ** NOM (DEFI_FONCTION) : ', nomfct
        else
            call lislnf(lischa, ichar, nomfct)
            write (6, *) '  ** NOM (INTERNE)       : ', nomfct
        end if
!
        if (typfct(7:10) .eq. 'COMP') then
            call liscpp(lischa, ichar, phase, npuis)
            write (6, *) '  ** PHASE               : ', phase
            write (6, *) '  ** PUISSANCE           : ', npuis
        end if
!
! ----- BOUCLE SUR LES GENRES
!
        write (6, *) '  * GENRES DE LA CHARGE:'
        do igenr = 1, nbgenr(1)
            gencha = zk24(jlisg-1+igenr)
            nomlis = '&&LISIMP.NOMLIS'
!
! ------- POSITION ENTIER CODE POUR CE GENRE
!
            call lisdef('POSG', gencha, ibid, k8bid, iposit)
            call isdeco([genrec], tabcod, 30)
            if (tabcod(iposit(1)) .eq. 1) then
!
! --------- GENRE PRESENT DANS CETTE CHARGE
!
                write (6, *) '  ** GENRE        : ', gencha
            end if
        end do
    end do
!
    call jedetr(lisgen)
!
99  continue
!
    call jedema()
end subroutine
