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

subroutine vtdef1(chpout, chpin, base, typc)
!     DEFINITION DE LA STRUCTURE D'UN CHAM_NO OU CHAM_ELEM "CHPOUT"
!                    QUI S'APPUIE SUR LA MEME NUMEROTATION QUE "CHPIN",
!     LE CHAM_... "CHPOUT" EST CREEE SUR LA BASE "BASE".
!     LE CHAM_... "CHPOUT" EST A COEFFICIENTS "TYPE".
!     ------------------------------------------------------------------
! IN : CHPOUT : NOM DU CHAM_NO OU CHAM_ELEM A CREER
! IN : CHPIN  : NOM DU CHAM_NO OU CHAM_ELEM MODELE
! IN : BASE   : NOM DE LA BASE SUR LAQUELLE LE CHAM_... DOIT ETRE CREER
! IN : TYPC   : TYPE DES VALEURS DU CHAM_... A CREER
!                    'R'  ==> COEFFICIENTS REELS
!                    'C'  ==> COEFFICIENTS COMPLEXES
!                    ' '  ==> COEFFICIENTS DU TYPE DU CHAM_... CHPIN
!     ------------------------------------------------------------------
!     PRECAUTIONS D'EMPLOI :
!       1) LE CHAM_... "CHPOUT" NE DOIT PAS EXISTER
!       2) LES COEFFICIENTS DU CHAM_... "CHPOUT" NE SONT PAS AFFECTES
!     -----------------------------------------------------------------
!     ASTER INFORMATIONS:
!       16/01/04 (OB): CREATION.
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sdchgd.h"
#include "asterfort/utmess.h"
    character(len=*) :: chpout, chpin, base, typc
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: nbval, ival, lchpou, lchpin, lchp, nbval1
    character(len=1) :: classe, type
    character(len=4) :: tych, docu
    character(len=19) :: ch19
    character(len=24) :: vale, refe, desc, celk, tamp
!
!     ------------------------------------------------------------------
    data refe/'                   .REFE'/
    data celk/'                   .CELK'/
!     ------------------------------------------------------------------
!
    call jemarq()
    classe = base(1:1)
    ch19 = chpin
!
    call dismoi('TYPE_CHAMP', ch19, 'CHAMP', repk=tych)
!
    if (tych .eq. 'NOEU') then
        docu = 'CHNO'
        tamp = refe
        desc(20:24) = '.REFE'
        vale(20:24) = '.VALE'
    else if (tych(1:2) .eq. 'EL') then
        docu = 'CHML'
        desc(20:24) = '.CELD'
        vale(20:24) = '.CELV'
        tamp = celk
    else
        call utmess('F', 'UTILITAI_21')
    end if
!
!     --------------------------- CELK --------------------------------
!     --- RECUPERATION DES INFORMATIONS DE CHPIN ---
    tamp(1:19) = chpin
    call jelira(tamp, 'LONMAX', nbval)
    call jeveuo(tamp, 'L', lchpin)
!
!     --- AFFECTATION DES INFORMATIONS A CHPOUT ---
    tamp(1:19) = chpout
    call jecreo(tamp, classe//' V K24')
    call jeecra(tamp, 'LONMAX', nbval)
    call jeecra(tamp, 'LONUTI', nbval)
    call jeveuo(tamp, 'E', lchpou)
    nbval1 = nbval-1
    do ival = 0, nbval1
        zk24(lchpou+ival) = zk24(lchpin+ival)
    end do
!
    tamp(1:19) = chpin
    tamp(1:19) = chpout
!
!     --------------------------- DESC --------------------------------
!     --- RECUPERATION DES INFORMATIONS DU DESCRIPTEUR CHPIN ---
    if (tych(1:2) .eq. 'EL') then
        desc(1:19) = chpin
        call jelira(desc, 'LONMAX', nbval)
        call jeveuo(desc, 'L', lchpin)
!
!     --- AFFECTATION DES INFORMATIONS DE DESCRIPTEUR CHPOUT ---
        desc(1:19) = chpout
        call jecreo(desc, classe//' V I')
        call jeecra(desc, 'LONMAX', nbval)
        call jeecra(desc, 'LONUTI', nbval)
        nbval1 = nbval-1
!
!
        call jeveuo(desc, 'E', lchpou)
        do ival = 0, nbval1
            zi(lchpou+ival) = zi(lchpin+ival)
        end do
!
    end if
    desc(1:19) = chpout
    call jeecra(desc, 'DOCU', cval=docu)
!
!     --------------------------- VALE --------------------------------
    vale(1:19) = chpin
    type = typc(1:1)
    if (type .eq. ' ') call jelira(vale, 'TYPE', cval=type)
    if (tych(1:2) .eq. 'EL') then
        ASSERT(type .ne. 'C')
    end if
    call jelira(vale, 'LONMAX', nbval)
    vale(1:19) = chpout
    call jecreo(vale, classe//' V '//type)
    call jeecra(vale, 'LONMAX', nbval)
    call jeecra(vale, 'LONUTI', nbval)
    call jeveuo(vale, 'E', lchp)
!
!
!     --- CHANGER LA GRANDEUR ---
    call sdchgd(chpout, type)
    call jedema()
end subroutine
