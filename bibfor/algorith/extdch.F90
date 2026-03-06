! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine extdch(typeExtr, hvalIncr, fieldType, cmpName, vale, lst_loca)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/barych.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsred.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmchex.h"
!
    character(len=8), intent(in)  :: typeExtr
    character(len=19), intent(in) :: hvalIncr(*)
    character(len=16), intent(in) :: fieldType
    character(len=16), intent(in) :: cmpName
    real(kind=8), intent(out) :: vale
    integer(kind=8), optional, intent(in) :: lst_loca(:)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
!    CALCUL D'UN EXTREMUM (MIN, MAX, EN VALEUR ABSOLUE OU NON)
!    DE L'INCREMENT D'UN CHAMP ('DEPL', 'SIEL_ELGA', OU 'VARI_ELGA')
!
! --------------------------------------------------------------------------------------------------
!
! IN  TYPEXT : TYPE D'EXTREMUM : MIN(), MAX(), MIN(ABS()), MAX(ABS())
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  NOCHAM : NOM DU CHAMP
! IN  NOCMP  : NOM DE LA COMPOSANTE
! OUT DVAL   : EXTREMUM
! IN  LST_LOCA: Liste des mailles/noeuds a prendre en compte (si absent: tout)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jcnsl
    integer(kind=8) :: nbno, nodeNume
    integer(kind=8) :: nbma, cellNume, ipt, isp, icmp, nbpt, nbsp, nbcmp, idx, nb_idx
    integer(kind=8) :: jmoid, jmoil, jmoiv, imoiad
    integer(kind=8) :: jplud, jplul, jpluv, ipluad
    real(kind=8) :: valeIncr, valePrev, valeCurr
    character(len=16) :: fieldDisc
    character(len=19) :: fieldCurr, fieldPrev
    character(len=19), parameter :: dch = '&&EXTDCH.DELTACH', dchs = '&&EXTDCH.DELTACHS'
    character(len=19), parameter :: fieldPrevS = '&&EXTDCH.CHMOIS', fieldCurrS = '&&EXTDCH.CHPLUS'
    aster_logical :: typeValid, filtered
    real(kind=8), pointer :: cnsv(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    typeValid = typeExtr .eq. 'MIN' .or. typeExtr .eq. 'MAX' .or. &
        typeExtr .eq. 'MIN_ABS' .or. typeExtr .eq. 'MAX_ABS' .or. &
        typeExtr .eq. 'MIN_VAR'
    ASSERT(typeValid)
!
    ASSERT(fieldType .eq. 'VARI_ELGA' .or. fieldType .eq. 'SIEF_ELGA' .or. fieldType .eq. 'DEPL')

! - Extract field
    if (fieldType .eq. 'VARI_ELGA') then
        call nmchex(hvalIncr, 'VALINC', 'VARMOI', fieldPrev)
        call nmchex(hvalIncr, 'VALINC', 'VARPLU', fieldCurr)
        fieldDisc = 'CHAM_ELGA'
    else if (fieldType .eq. 'SIEF_ELGA') then
        call nmchex(hvalIncr, 'VALINC', 'SIGMOI', fieldPrev)
        call nmchex(hvalIncr, 'VALINC', 'SIGPLU', fieldCurr)
        fieldDisc = 'CHAM_ELGA'
    else if (fieldType .eq. 'DEPL') then
        call nmchex(hvalIncr, 'VALINC', 'DEPMOI', fieldPrev)
        call nmchex(hvalIncr, 'VALINC', 'DEPPLU', fieldCurr)
        fieldDisc = 'CHAM_NO'
    else
        ASSERT(ASTER_FALSE)
    end if

    filtered = present(lst_loca)

! - Initialization of bound
    if (typeExtr .eq. 'MIN' .or. typeExtr .eq. 'MIN_ABS' .or. typeExtr .eq. 'MIN_VAR') then
        vale = r8maem()
    end if
    if (typeExtr .eq. 'MAX') then
        vale = r8miem()
    end if
    if (typeExtr .eq. 'MAX_ABS') then
        vale = 0.d0
    end if
!

    if (fieldDisc(1:7) .eq. 'CHAM_EL') then
        call celces(fieldPrev, 'V', fieldPrevS)
        call cesred(fieldPrevS, 0, [0], 1, cmpName, 'V', fieldPrevS)
        call celces(fieldCurr, 'V', fieldCurrS)
        call cesred(fieldCurrS, 0, [0], 1, cmpName, 'V', fieldCurrS)
!
        call jeveuo(fieldPrevS//'.CESD', 'L', jmoid)
        call jeveuo(fieldCurrS//'.CESD', 'L', jplud)
        call jeveuo(fieldPrevS//'.CESL', 'L', jmoil)
        call jeveuo(fieldCurrS//'.CESL', 'L', jplul)
        call jeveuo(fieldPrevS//'.CESV', 'L', jmoiv)
        call jeveuo(fieldCurrS//'.CESV', 'L', jpluv)
!
        nbma = zi(jmoid-1+1)
        ASSERT(zi(jplud-1+1) .eq. nbma)
!
        if (filtered) then
            nb_idx = size(lst_loca)
        else
            nb_idx = nbma
        end if

        do idx = 1, nb_idx
            if (filtered) then
                cellNume = lst_loca(idx)
            else
                cellNume = idx
            end if

            nbpt = zi(jmoid-1+5+4*(cellNume-1)+1)
            nbsp = zi(jmoid-1+5+4*(cellNume-1)+2)
            nbcmp = zi(jmoid-1+5+4*(cellNume-1)+3)
!
            ASSERT(zi(jplud-1+5+4*(cellNume-1)+1) .eq. nbpt)
            ASSERT(zi(jplud-1+5+4*(cellNume-1)+2) .eq. nbsp)
            ASSERT(zi(jplud-1+5+4*(cellNume-1)+3) .eq. nbcmp)
!
            do ipt = 1, nbpt
                do isp = 1, nbsp
                    do icmp = 1, nbcmp
                        call cesexi('C', jmoid, jmoil, cellNume, ipt, isp, 1, imoiad)
                        call cesexi('C', jplud, jplul, cellNume, ipt, isp, 1, ipluad)
                        if (imoiad .gt. 0 .or. ipluad .gt. 0) then
                            ASSERT(imoiad .gt. 0 .and. ipluad .gt. 0)
                            valePrev = zr(jmoiv-1+imoiad)
                            valeCurr = zr(jpluv-1+ipluad)
                            valeIncr = valeCurr-valePrev
                            if (typeExtr(5:7) .eq. 'ABS') then
                                valeIncr = abs(valeIncr)
                            end if
                            if (typeExtr(5:7) .eq. 'VAR') then
                                if (abs(valeIncr) .gt. r8prem()) then
                                    valeIncr = 1.d-3/abs(valeIncr)
                                else
                                    valeIncr = r8maem()
                                end if
                            end if

                            if (typeExtr(1:3) .eq. 'MIN') then
                                vale = min(vale, valeIncr)

                            else if (typeExtr(1:3) .eq. 'MAX') then
                                vale = max(vale, valeIncr)

                            end if
                        end if
                    end do
                end do
            end do
        end do
!
    else if (fieldDisc .eq. 'CHAM_NO') then
        call barych(fieldCurr, fieldPrev, 1.d0, -1.d0, dch, 'V')
        call cnocns(dch, 'V', dchs)
        call cnsred(dchs, 0, [0], 1, cmpName, 'V', dchs)
        call jeveuo(dchs//'.CNSV', 'L', vr=cnsv)
        call jeveuo(dchs//'.CNSL', 'L', jcnsl)
        call jeveuo(dchs//'.CNSD', 'L', vi=cnsd)
        nbno = cnsd(1)

        if (filtered) then
            nb_idx = size(lst_loca)
        else
            nb_idx = nbno
        end if

        do idx = 1, nb_idx
            if (filtered) then
                nodeNume = lst_loca(idx)
            else
                nodeNume = idx
            end if

            if (zl(jcnsl-1+nodeNume)) then
                valeIncr = abs(cnsv(nodeNume))
                if (typeExtr(5:7) .eq. 'ABS') then
                    valeIncr = abs(valeIncr)
                end if
                if (typeExtr(1:3) .eq. 'MIN') then
                    vale = min(vale, valeIncr)
                else if (typeExtr(1:3) .eq. 'MAX') then
                    vale = max(vale, valeIncr)
                end if
            end if
        end do
!
    end if
!

    call jedema()
end subroutine
