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
!
subroutine rcver1(phenoz, tablz, tably)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbexip.h"
#include "asterfort/tbexv1.h"
#include "asterfort/tbexve.h"
#include "asterfort/utmess.h"
    character(len=*) :: tablz, tably, phenoz
!      OPERATEUR POST_RCCM
!
!     IN  TABLZ: TABLE DE REFERENCE
!     IN  TABLY: TABLE A COMPARER A TABLZ
!     IN  PHENOZ : PHENOMENE (THERMIQUE OU MECANIQUE)
!
!      VERIFICATIONS SUR LES 2 TABLES:
!      - MEMES COORDONNEES  (CAS "EVOLUTION" ET "UNITAIRE"):
!          -> SI 'MECANIQUE': MEME NOMBRE DE NOEUDS ET MEMES COORDONNEES
!          -> SI 'THERMIQUE': MEMES NOEUDS EXTREMITES
!      - MEME NOMBRE DE LIGAMENTS, MEMES INSTANTS (CAS "EVOLUTION")
!     ------------------------------------------------------------------
    integer(kind=8) :: n1, nbno1, nbno2, jordo1, jordo2, nbint1, nbint2, i, j
    integer(kind=8) :: vali(2), nbins1, nbins2, jinst1, jinst2
    real(kind=8) :: valr(2), eps, v1, v2
    parameter(eps=1.d-6)
    character(len=8) :: k8b, valk(3), tyva
    character(len=16) :: tabref, tabcom, typmec, valek(5), phenom
    character(len=24) :: ordo1, ordo2, intit1, intit2, inst1, inst2
    aster_logical :: exi1, exi2, exi3, exist
!
    call jemarq()
!
    valek(1) = 'COOR_X          '
    valek(2) = 'COOR_Y          '
    valek(3) = 'COOR_Z          '
    valek(4) = 'INTITULE        '
    valek(5) = 'INST            '
!
    ordo1 = '&&RCVER1_TABREF_ORDO'
    ordo2 = '&&RCVER1_TABCOM_ORDO'
    intit1 = '&&RCVER1_TABREF_INTITU'
    intit2 = '&&RCVER1_TABCOM_INTITU'
    inst1 = '&&RCVER1_TABREF_INST'
    inst2 = '&&RCVER1_TABCOM_INST'
!
    tabref = tablz
    tabcom = tably
    phenom = phenoz
!
! --- VERIFICATION DES VALEURS DES COORDONNEES
!     ----------------------------------------
!
!     VERIFICATION DE LA PRESENCE DES COORDONNEES DANS LA TABLE 'TABCOM'
    call tbexip(tabcom, valek(1), exi1, k8b)
    call tbexip(tabcom, valek(2), exi2, k8b)
    call tbexip(tabcom, valek(3), exi3, k8b)
    if (.not. exi1 .and. .not. exi2 .and. .not. exi3) then
        call utmess('I', 'POSTRCCM_39', sk=tabcom)
        goto 999
    end if
!
    do j = 1, 3
        call tbexve(tabref, valek(j), ordo1, 'V', nbno1, &
                    k8b)
        call jeveuo(ordo1, 'L', jordo1)
        call tbexve(tabcom, valek(j), ordo2, 'V', nbno2, &
                    k8b)
        call jeveuo(ordo2, 'L', jordo2)
!       ON COMPARE LES COORDONNEES DE TOUS LES NOEUDS
        if (phenom .eq. 'MECANIQUE') then
            ASSERT(nbno1 .eq. nbno2)
            do i = 1, nbno1
                if (abs(zr(jordo1+i-1)-zr(jordo2+i-1)) .gt. eps) then
                    valk(1) = tabref(1:8)
                    valk(2) = tabcom(1:8)
                    valk(3) = valek(j) (1:8)
                    valr(1) = zr(jordo1+i-1)
                    valr(2) = zr(jordo2+i-1)
                    call utmess('F', 'POSTRCCM_41', nk=3, valk=valk, nr=2, &
                                valr=valr)
                end if
            end do
!       ON COMPARE LES COORDONNEES DES NOEUDS EXTREMITES
!       (CAR ON N'A PAS FORCEMENT NBNO1 = NBNO2)
        else if (phenom .eq. 'THERMIQUE') then
            do i = 1, 2
                v1 = zr(jordo1+(nbno1-1)*(i-1))
                v2 = zr(jordo2+(nbno2-1)*(i-1))
                if (abs(v1-v2) .gt. eps) then
                    valk(1) = tabref(1:8)
                    valk(2) = tabcom(1:8)
                    valk(3) = valek(j) (1:8)
                    valr(1) = v1
                    valr(2) = v2
                    call utmess('F', 'POSTRCCM_41', nk=3, valk=valk, nr=2, &
                                valr=valr)
                end if
            end do
        end if
        call jedetr(ordo1)
        call jedetr(ordo2)
    end do
!
! --- VERIFICATION DU NOMBRE DE LIGAMENTS
!     -----------------------------------
    call getvtx(' ', 'TYPE_RESU_MECA', scal=typmec, nbret=n1)
    call tbexip(tabcom, valek(4), exist, k8b)
    ASSERT(exist)
    call tbexv1(tabcom, valek(4), intit2, 'V', nbint2, &
                tyva)
!     CAS UNITAIRE: 1 SEUL LIGAMENT
    if (typmec .eq. 'UNITAIRE' .and. nbint2 .ne. 1) then
        call utmess('F', 'POSTRCCM_40', sk=tabcom, si=nbint2)
!     CAS EVOLUTION: MEME NOMBRE DE LIGAMENTS
    else if (typmec .eq. 'EVOLUTION') then
        call tbexv1(tabref, valek(4), intit1, 'V', nbint1, &
                    tyva)
        if (nbint1 .ne. nbint2) then
            valk(1) = tabref(1:8)
            valk(2) = tabcom(1:8)
            vali(1) = nbint1
            vali(2) = nbint2
            call utmess('F', 'POSTRCCM_42', nk=2, valk=valk, ni=2, &
                        vali=vali)
        end if
        call jedetr(intit1)
    end if
    call jedetr(intit2)
!
! --- VERIFICATION DE LA VALEUR DES INSTANTS (CAS 'EVOLUTION')
!     --------------------------------------------------------
    if (typmec .eq. 'EVOLUTION') then
        call tbexip(tabcom, valek(5), exist, k8b)
        ASSERT(exist)
        call tbexve(tabref, valek(5), inst1, 'V', nbins1, &
                    k8b)
        call jeveuo(inst1, 'L', jinst1)
        call tbexve(tabcom, valek(5), inst2, 'V', nbins2, &
                    k8b)
        call jeveuo(inst2, 'L', jinst2)
        ASSERT(nbins1 .eq. nbins2)
        do i = 1, nbins1
            if (abs(zr(jinst1+i-1)-zr(jinst2+i-1)) .gt. eps) then
                valk(1) = tabref(1:8)
                valk(2) = tabcom(1:8)
                valr(1) = zr(jinst1+i-1)
                valr(2) = zr(jinst2+i-1)
                call utmess('F', 'POSTRCCM_43', nk=2, valk=valk, nr=2, &
                            valr=valr)
            end if
        end do
        call jedetr(inst1)
        call jedetr(inst2)
    end if
!
!
999 continue
!
    call jedema()
!
end subroutine
