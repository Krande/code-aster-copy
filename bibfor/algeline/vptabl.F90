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
subroutine vptabl(tabmod, typevp, fmin, fmax, precdc, &
                  nfreq, effmin, effmax)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/tbexip.h"
#include "asterfort/tbextb.h"
#include "asterfort/tbexve.h"
#include "asterfort/utmess.h"
    real(kind=8) :: fmin, fmax, precdc, effmin, effmax
    integer(kind=8) :: nfreq
    character(len=9) :: typevp
    character(len=19) :: tabmod
! DANS LA TABLE TABMOD RECUPERE LE NBRE DE FREQUENCES CORRESPONDANT A
! L'INTERVALLE FREMIN/FREMAX. ON DECALE LES BORNES SUIVANT LA REGLE
! PRATIQUEE DS LES OPERATEURS MODAUX VIA LA PARAMETRE
! VERIF_MODE/PREC_SHIFT
! IN TABMOD    : K19 : TABLE ISSUE DE INFO_MODE
! IN TYPEVP    : K9  : TYPE DE PROBLEME, FREQ OU CHAR_CRIT
! IN FMIN/FMAX : R8  : BORNES DE LA BANDE FREQUENTIELLE DEMANDEE
! IN PRECDC    : R8  : DECALAGE EXTERNE POUR ENGLOBER LA BANDE
! OUT NFREQ    : IN  : NBRE DE MODES DANS LA BANDE
! OUT EFFMIN/EFFMAX : R8  : BORNES DE LA BANDE FREQUENTIELLE EFFECTIVE
!-----------------------------------------------------------------------
! person_in_charge: albert.alarcon at edf.fr
!
! VARIABLES LOCALES
    integer(kind=8) :: ier, nbval1, nbval2, nbval3, jobj1, jobj2, jobj3, jobj4, jobj5, i
    integer(kind=8) :: ifm, niv, ibid, ll, nbval4, nbval5
    real(kind=8) :: vr(2), lprec(2), rdenom, eps
    complex(kind=8) :: cbid
    character(len=4) :: typval
    character(len=8) :: k8bid, lcrit(2)
    character(len=10) :: lcrpa(2)
    character(len=16) :: lipacr(2), para
    character(len=19) :: tabmof
    character(len=24) :: valk(2), nomob1, nomob2, nomob3, nomob4, nomob5
    aster_logical :: lexist
!
!
! --- INITS.
    cbid = (0.d0, 0.d0)
    ibid = 0
    call jemarq()
    call infniv(ifm, niv)
    valk(1) = tabmod
    tabmof = '&&VPTABL.TABMOF'
    nomob1 = '&&VPTABL.NB_MODE'
    nomob2 = '&&VPTABL.BORNE_MIN'
    nomob3 = '&&VPTABL.BORNE_MAX'
    nomob4 = '&&VPTABL.EFFEC_MIN'
    nomob5 = '&&VPTABL.EFFEC_MAX'
    nfreq = -9999
    ll = lxlgut(typevp)
    eps = 100.d0*r8prem()
!
! --- TESTS DE VALIDITE DE LA CARTE
    call tbexip(tabmod, typevp(1:ll)//'_MIN', lexist, k8bid)
    if ((lexist .and. (k8bid(1:1) .ne. 'R')) .or. (.not. lexist)) then
        call utmess('F', 'ALGELINE2_23', sk=valk(1))
    end if
!
    call tbexip(tabmod, typevp(1:ll)//'_MAX', lexist, k8bid)
    if ((lexist .and. (k8bid(1:1) .ne. 'R')) .or. (.not. lexist)) then
        call utmess('F', 'ALGELINE2_23', sk=valk(1))
    end if
!
    call tbexip(tabmod, 'NB_MODE', lexist, k8bid)
    if ((lexist .and. (k8bid(1:1) .ne. 'I')) .or. (.not. lexist)) then
        call utmess('F', 'ALGELINE2_23', sk=valk(1))
    end if
!
    call tbexip(tabmod, 'BORNE_MIN_EFFECT', lexist, k8bid)
    if ((lexist .and. (k8bid(1:1) .ne. 'R')) .or. (.not. lexist)) then
        call utmess('F', 'ALGELINE2_23', sk=valk(1))
    end if
!
    call tbexip(tabmod, 'BORNE_MAX_EFFECT', lexist, k8bid)
    if ((lexist .and. (k8bid(1:1) .ne. 'R')) .or. (.not. lexist)) then
        call utmess('F', 'ALGELINE2_23', sk=valk(1))
    end if
!
! --- EXTRACTION DES LIGNES REPONDANT AU CRITERE DS LA TABLE
! --- INTERMEDIAIRE: TABMOF
    lipacr(1) = typevp(1:ll)//'_MIN'
    lipacr(2) = typevp(1:ll)//'_MAX'
    lcrpa(1) = 'GE'
    lcrpa(2) = 'LE'
    vr(1) = fmin-eps
    vr(2) = fmax+eps
    lprec(1) = precdc
    lprec(2) = precdc
    lcrit(1) = 'RELA'
    lcrit(2) = 'RELA'
!
    call tbextb(tabmod, 'V', tabmof, 2, lipacr, &
                lcrpa, [ibid], vr, [cbid], k8bid, &
                lprec, lcrit, ier)
!
! --- PB EXTRACTION: PAS DE LIGNE CORRESPONDANT AUX CRITERES
    if ((ier .eq. 1) .or. (ier .eq. 2)) then
!
        call utmess('F', 'ALGELINE2_24', sk=valk(1), nr=2, valr=vr)
!
    else
!
! --- EXTRACTION DES COLONNES 'MIN','MAX' ET 'NB_MODE' + VERIFICATIONS
        para = 'NB_MODE'
        call tbexve(tabmof, para, nomob1, 'V', nbval1, &
                    typval)
        if ((typval .ne. 'I') .or. (nbval1 .eq. 0)) then
            call utmess('F', 'ALGELINE2_23', sk=valk(1))
        end if
        para = lipacr(1)
        call tbexve(tabmof, para, nomob2, 'V', nbval2, &
                    typval)
        if ((typval .ne. 'R') .or. (nbval2 .eq. 0)) then
            call utmess('F', 'ALGELINE2_23', sk=valk(1))
        end if
        para = lipacr(2)
        call tbexve(tabmof, para, nomob3, 'V', nbval3, &
                    typval)
        if ((typval .ne. 'R') .or. (nbval3 .eq. 0)) then
            call utmess('F', 'ALGELINE2_23', sk=valk(1))
        end if
        para = 'BORNE_MIN_EFFECT'
        call tbexve(tabmof, para, nomob4, 'V', nbval4, &
                    typval)
        if ((typval .ne. 'R') .or. (nbval4 .eq. 0)) then
            call utmess('F', 'ALGELINE2_23', sk=valk(1))
        end if
        para = 'BORNE_MAX_EFFECT'
        call tbexve(tabmof, para, nomob5, 'V', nbval5, &
                    typval)
        if ((typval .ne. 'R') .or. (nbval5 .eq. 0)) then
            call utmess('F', 'ALGELINE2_23', sk=valk(1))
        end if
!
! --- VERIF NBRE DE LIGNES EGAUX
        if ((nbval1 .ne. nbval2) .or. (nbval1 .ne. nbval3) .or. (nbval2 .ne. nbval3)) then
            call utmess('F', 'ALGELINE2_23', sk=valk(1))
        end if
        if ((nbval4 .ne. nbval5) .or. (nbval1 .ne. nbval4) .or. (nbval1 .ne. nbval5)) then
            call utmess('F', 'ALGELINE2_23', sk=valk(1))
        end if
!
! --- VERIF PAS DE TROU DS LES OBJETS
        call jeveuo(nomob2, 'L', jobj2)
        call jeveuo(nomob3, 'L', jobj3)
        call jeveuo(nomob4, 'L', jobj4)
        call jeveuo(nomob5, 'L', jobj5)
!
! POUR DEBUGGAGE: LISTE DES BORNES DES INTERVALLES SELECTIONNES
!       DO I=1,NBVAL1
!         WRITE(6,*)I,ZR(JOBJ2+I-1),ZR(JOBJ3+I-1),
!     &              ZR(JOBJ4+I-1),ZR(JOBJ5+I-1)
!       ENDDO
!
! --- VERIF DES BORNES EXTREMES
        if (abs(fmin) .gt. eps) then
            rdenom = abs(fmin)
        else
            rdenom = 1.d0
        end if
        if (abs(zr(jobj2)-fmin)/rdenom .gt. precdc) then
            call utmess('F', 'ALGELINE2_25', sk=valk(1))
        end if
!
        if (abs(fmax) .gt. eps) then
            rdenom = abs(fmax)
        else
            rdenom = 1.d0
        end if
        if (abs(zr(jobj3+nbval3-1)-fmax)/rdenom .gt. precdc) then
            call utmess('F', 'ALGELINE2_25', sk=valk(1))
        end if
!
! --- VERIF DES BORNES INTERNES (INITIALES ET EFFECTIVES)
        do i = 1, nbval3-1
            if (abs(zr(jobj3+i-1)) .gt. eps) then
                rdenom = abs(zr(jobj3+i-1))
            else
                rdenom = 1.d0
            end if
            if (abs(zr(jobj3+i-1)-zr(jobj2+i))/rdenom .gt. precdc) then
                call utmess('F', 'ALGELINE2_25', sk=valk(1))
            end if
        end do
        do i = 1, nbval5-1
            if (abs(zr(jobj5+i-1)) .gt. eps) then
                rdenom = abs(zr(jobj5+i-1))
            else
                rdenom = 1.d0
            end if
            if (abs(zr(jobj5+i-1)-zr(jobj4+i))/rdenom .gt. precdc) then
                call utmess('F', 'ALGELINE2_25', sk=valk(1))
            end if
        end do
!
! --- SOMME DES NB_MODES DES BANDES SELECTIONNEES
        call jeveuo(nomob1, 'L', jobj1)
        nfreq = 0
        do i = 1, nbval1
            nfreq = nfreq+zi(jobj1+i-1)
        end do
!
! --- RECUPERATION DES BORNES EFFECTIVES
        effmin = zr(jobj4)
        effmax = zr(jobj5+nbval5-1)
!
        call jedetr(nomob1)
        call jedetr(nomob2)
        call jedetr(nomob3)
        call jedetr(nomob4)
        call jedetr(nomob5)
    end if
    call jedetr(tabmof)
!
    call jedema()
end subroutine
