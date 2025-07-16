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
subroutine extdch(typext, valinc, nocham, nocmp, dval, lst_loca)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
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
    character(len=8), intent(in)  :: typext
    character(len=19), intent(in) :: valinc(*)
    character(len=16), intent(in) :: nocham
    character(len=16), intent(in) :: nocmp
    real(kind=8), intent(out) :: dval
    integer(kind=8), optional, intent(in) :: lst_loca(:)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
!    CALCUL D'UN EXTREMUM (MIN, MAX, EN VALEUR ABSOLUE OU NON)
!    DE L'INCREMENT D'UN CHAMP ('DEPL', 'SIEL_ELGA', OU 'VARI_ELGA')
!
! ----------------------------------------------------------------------
!
!
! IN  TYPEXT : TYPE D'EXTREMUM : MIN(), MAX(), MIN(ABS()), MAX(ABS())
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  NOCHAM : NOM DU CHAMP
! IN  NOCMP  : NOM DE LA COMPOSANTE
! OUT DVAL   : EXTREMUM
! IN  LST_LOCA: Liste des mailles/noeuds a prendre en compte (si absent: tout)
!
!
!
!
!
    integer(kind=8) :: jcnsl
    integer(kind=8) :: nbno, ino
    integer(kind=8) :: nbma, ima, ipt, isp, icmp, nbpt, nbsp, nbcmp, idx, nb_idx
    integer(kind=8) :: jmoid, jmoil, jmoiv, imoiad
    integer(kind=8) :: jplud, jplul, jpluv, ipluad
    real(kind=8) :: valeur, vmoi, vplu
    character(len=6) :: nompro
    character(len=16) :: typch
    character(len=19) :: dch, dchs, chplu, chmoi, chmois, chplus
    parameter(nompro='EXTDCH')
    aster_logical :: bool, filtered
    real(kind=8), pointer :: cnsv(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
!
!      REAL*8  TMP
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    bool = typext .eq. 'MIN' .or. typext .eq. 'MAX' .or. typext .eq. 'MIN_ABS' .or. typext .eq. &
           'MAX_ABS' .or. typext .eq. 'MIN_VAR'
    ASSERT(bool)
!
    ASSERT(nocham .eq. 'VARI_ELGA' .or. nocham .eq. 'SIEF_ELGA' .or. nocham .eq. 'DEPL')
!
!     DECOMPACTION DES VARIABLES CHAPEAUX
    if (nocham .eq. 'VARI_ELGA') then
        call nmchex(valinc, 'VALINC', 'VARMOI', chmoi)
        call nmchex(valinc, 'VALINC', 'VARPLU', chplu)
        typch = 'CHAM_ELGA'
    else if (nocham .eq. 'SIEF_ELGA') then
        call nmchex(valinc, 'VALINC', 'SIGMOI', chmoi)
        call nmchex(valinc, 'VALINC', 'SIGPLU', chplu)
        typch = 'CHAM_ELGA'
    else if (nocham .eq. 'DEPL') then
        call nmchex(valinc, 'VALINC', 'DEPMOI', chmoi)
        call nmchex(valinc, 'VALINC', 'DEPPLU', chplu)
        typch = 'CHAM_NO'
    end if

    filtered = present(lst_loca)

!
!     INITIALISATION DE L'EXTREMUM
    if (typext .eq. 'MIN' .or. typext .eq. 'MIN_ABS' .or. typext .eq. 'MIN_VAR') dval = &
        r8maem()
!
    if (typext .eq. 'MAX') dval = r8miem()
    if (typext .eq. 'MAX_ABS') dval = 0.d0
!
!
!
!     ON APPELLERA MEMAX QUAND CETTE ROUTINE SERA MIEUX PROGRAMMEE
    if (typch(1:7) .eq. 'CHAM_EL') then
!
        chmois = '&&'//nompro//'.CHMOIS'
        chplus = '&&'//nompro//'.CHPLUS'
!
        call celces(chmoi, 'V', chmois)
        call cesred(chmois, 0, [0], 1, nocmp, &
                    'V', chmois)
        call celces(chplu, 'V', chplus)
        call cesred(chplus, 0, [0], 1, nocmp, &
                    'V', chplus)
!
        call jeveuo(chmois//'.CESD', 'L', jmoid)
        call jeveuo(chplus//'.CESD', 'L', jplud)
        call jeveuo(chmois//'.CESL', 'L', jmoil)
        call jeveuo(chplus//'.CESL', 'L', jplul)
        call jeveuo(chmois//'.CESV', 'L', jmoiv)
        call jeveuo(chplus//'.CESV', 'L', jpluv)
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
                ima = lst_loca(idx)
            else
                ima = idx
            end if

            nbpt = zi(jmoid-1+5+4*(ima-1)+1)
            nbsp = zi(jmoid-1+5+4*(ima-1)+2)
            nbcmp = zi(jmoid-1+5+4*(ima-1)+3)
!
            ASSERT(zi(jplud-1+5+4*(ima-1)+1) .eq. nbpt)
            ASSERT(zi(jplud-1+5+4*(ima-1)+2) .eq. nbsp)
            ASSERT(zi(jplud-1+5+4*(ima-1)+3) .eq. nbcmp)
!
            do ipt = 1, nbpt
                do isp = 1, nbsp
                    do icmp = 1, nbcmp
                        call cesexi('C', jmoid, jmoil, ima, ipt, &
                                    isp, 1, imoiad)
                        call cesexi('C', jplud, jplul, ima, ipt, &
                                    isp, 1, ipluad)
!
                        if (imoiad .gt. 0 .or. ipluad .gt. 0) then
!
                            ASSERT(imoiad .gt. 0 .and. ipluad .gt. 0)
!
                            vmoi = zr(jmoiv-1+imoiad)
                            vplu = zr(jpluv-1+ipluad)
                            valeur = vplu-vmoi
!
                            if (typext(5:7) .eq. 'ABS') valeur = abs(valeur)
                            if (typext(5:7) .eq. 'VAR') then
                                if (abs(valeur) .gt. r8prem()) then
                                    valeur = 1.d-3/abs(valeur)
!                      DVAL = MIN(DVAL,TMP)
                                else
                                    valeur = r8maem()
                                end if
!
                            end if
                            if (typext(1:3) .eq. 'MIN') then
                                dval = min(dval, valeur)
!
                            else if (typext(1:3) .eq. 'MAX') then
                                dval = max(dval, valeur)
!
                            end if
!
                        end if
                    end do
                end do
            end do
        end do
!
    else if (typch .eq. 'CHAM_NO') then
!
!       CALCUL DE L'INCREMENT DU CHAMP
!       DCH = CHPLU - CHMOI
        dch = '&&'//nompro//'.DELTACH   '
        dchs = '&&'//nompro//'.DELTACHS  '
        call barych(chplu, chmoi, 1.d0, -1.d0, dch, &
                    'V')
!
        call cnocns(dch, 'V', dchs)
        call cnsred(dchs, 0, [0], 1, nocmp, &
                    'V', dchs)
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
                ino = lst_loca(idx)
            else
                ino = idx
            end if

            if (zl(jcnsl-1+ino)) then
                valeur = abs(cnsv(ino))
                if (typext(5:7) .eq. 'ABS') valeur = abs(valeur)
                if (typext(1:3) .eq. 'MIN') then
                    dval = min(dval, valeur)
                else if (typext(1:3) .eq. 'MAX') then
                    dval = max(dval, valeur)
                end if
            end if
        end do
!
    end if
!
    call jedema()
end subroutine
