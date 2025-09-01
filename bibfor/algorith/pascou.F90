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

subroutine pascou(mate, mateco, carele, sddyna, sddisc)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
!
    character(len=24) :: mate, mateco, carele
    character(len=19) :: sddyna, sddisc
!
! ----------------------------------------------------------------------
!
! ROUTINE DYNA_NON_LINE (UTILITAIRE)
!
! EVALUATION DU PAS DE TEMPS DE COURANT POUR LE MODELE
!
! ----------------------------------------------------------------------
!
!
! IN  MATE   : CHAMP MATERIAU
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! IN  SDDYNA : SD DEDIEE A LA DYNAMIQUE (CF NDLECT)
! IN  SDDISC : SD DISCRETISATION
!
!
!
!
    integer(kind=8), parameter :: nbFieldInMax = 4, nbFieldOutMax = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    integer(kind=8) :: ibid, jcesd, jcesl, n1, i, numins, nbFieldIn, nbFieldOut
    integer(kind=8) :: nbma, ima, iad, nbinst, nbmcfl
    real(kind=8) :: dtcou, valeur, phi, instin
    aster_logical :: booneg, boopos
    character(len=2) :: codret
    character(len=6) :: nompro
    character(len=8) :: mo, stocfl, maicfl, mail
    character(len=19) :: chams, chvarc
    character(len=24) :: chgeom, ligrel, chcara(18)
    real(kind=8), pointer :: ditr(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
!
! ---------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    nompro = 'OP0070'
    chvarc = '&&PASCOU.CH_VARC_R'
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '
!
    call getvid(' ', 'MODELE', scal=mo, nbret=ibid)
!
    call dismoi('NOM_LIGREL', mo, 'MODELE', repk=ligrel)
!
    lpain(1) = 'PMATERC'
    lchin(1) = mateco
!
! --- RECUPERATION DU CHAMP GEOMETRIQUE
    call megeom(mo, chgeom)
!
    lpain(2) = 'PGEOMER'
    lchin(2) = chgeom
!
! --- CHAMP DES VARIABLES DE COMMANDE
    numins = 0
    instin = diinst(sddisc, numins)
    call vrcins(mo, mate, carele, instin, chvarc, &
                codret)

    lpain(3) = 'PVARCPR'
    lchin(3) = chvarc(1:19)
!
! --- CHAMP DE CARACTERISTIQUES ELEMENTAIRES
    call mecara(carele(1:8), chcara)
!
    nbFieldIn = 3
    if (carele(1:8) .ne. ' ') then
        lpain(4) = 'PCACOQU'
        lchin(4) = chcara(7)
        nbFieldIn = nbFieldIn+1
    end if
!
    lpaout(1) = 'PCOURAN'
    lchout(1) = '&&'//nompro//'.PAS_COURANT'
    nbFieldOut = 1
!
    ASSERT(nbFieldIn .le. nbFieldInMax)
    ASSERT(nbFieldOut .le. nbFieldOutMax)
    call calcul('S', 'PAS_COURANT', ligrel, nbFieldIn, lchin, &
                lpain, nbFieldOut, lchout, lpaout, 'V', &
                'OUI')
!
!     PASSAGE D'UN CHAM_ELEM EN UN CHAM_ELEM_S
    chams = '&&'//nompro//'.CHAMS'
!
    call celces(lchout(1), 'V', chams)
!
    call jeveuo(chams//'.CESD', 'L', jcesd)
!
    call dismoi('NOM_LIGREL', mo, 'MODELE', repk=ligrel)
    call dismoi('NB_MA_MAILLA', ligrel, 'LIGREL', repi=nbma)
    call jeveuo(chams//'.CESL', 'L', jcesl)
    call jeveuo(chams//'.CESV', 'L', vr=cesv)
!
!     INITIALISATION DE DTCOU
!
    dtcou = -1.d0
!
! A L'ISSUE DE LA BOUCLE :
! BOONEG=TRUE SI L'ON N'A PAS PU CALCULER DTCOU POUR AU MOINS UN ELMNT
! BOOPOS=TRUE SI L'ON A CALCULE DTCOU POUR AU MOINS UN ELEMENT
    booneg = .false.
    boopos = .false.
    nbmcfl = 1
    do ima = 1, nbma
        call cesexi('C', jcesd, jcesl, ima, 1, &
                    1, 1, iad)
        if (iad .gt. 0) then
            valeur = cesv(iad)
        else if (iad .eq. 0) then
            goto 10
        end if
        if (valeur .lt. 0) then
            booneg = .true.
        else
            boopos = .true.
            if (dtcou .gt. 0) then
                if (valeur .le. dtcou) then
                    dtcou = valeur
                    nbmcfl = ima
                end if
            else
                dtcou = valeur
            end if
        end if
10      continue
    end do
!
    call getvtx('SCHEMA_TEMPS', 'STOP_CFL', iocc=1, scal=stocfl, nbret=n1)
!
! BOOPOS=TRUE SI L'ON A CALCULE DTCOU POUR AU MOINS UN ELEMENT
    if (boopos) then
        if (booneg) then
            call utmess('A', 'DYNAMIQUE_3')
        end if
!
!       VERIFICATION DE LA CONFORMITE DE LA LISTE D'INSTANTS
        call utdidt('L', sddisc, 'LIST', 'NBINST', &
                    vali_=nbinst)
        call jeveuo(sddisc//'.DITR', 'L', vr=ditr)
!
        call dismoi('NOM_MAILLA', mo, 'MODELE', repk=mail)
        maicfl = int_to_char8(nbmcfl)
!
!
        if (ndynlo(sddyna, 'DIFF_CENT')) then
            dtcou = dtcou/(2.d0)
            call utmess('I', 'DYNAMIQUE_5', sk=maicfl, sr=dtcou)
        else
            if (ndynlo(sddyna, 'TCHAMWA')) then
                phi = ndynre(sddyna, 'PHI')
                dtcou = dtcou/(phi*2.d0)
                call utmess('I', 'DYNAMIQUE_6', sk=maicfl, sr=dtcou)
            else
                call utmess('F', 'DYNAMIQUE_1')
            end if
        end if
!
        do i = 1, nbinst-1
            if (ditr(i+1)-ditr(i) .gt. dtcou) then
                if (stocfl(1:3) .eq. 'OUI') then
                    call utmess('F', 'DYNAMIQUE_2')
                else
                    call utmess('A', 'DYNAMIQUE_2')
                end if
            end if
        end do
!
    else if (stocfl(1:3) .eq. 'OUI') then
        call utmess('F', 'DYNAMIQUE_4')
    else if (stocfl(1:3) .eq. 'NON') then
        call utmess('A', 'DYNAMIQUE_4')
    end if
!
    call jedema()
!
end subroutine
