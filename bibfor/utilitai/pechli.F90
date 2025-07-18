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

subroutine pechli(resu, modele, mateco)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/memaxm.h"
#include "asterfort/mesomm.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsexpa.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
    character(len=8) :: modele
    character(len=24) :: mateco
    character(len=*) :: resu
!
!     OPERATEUR   POST_ELEM
!
!     TRAITEMENT DU MOT CLE-FACTEUR "CHAR_LIMITE"
!
! ----------------------------------------------------------------------
!
    aster_logical :: chrcst
    integer(kind=8) :: i, iret, jinst, jpilo
    integer(kind=8) :: nbord, jord, numord
    real(kind=8) :: chlim(3), chmax(3), inst, eta, prec, valer(4), f0u, m
    complex(kind=8) :: c16b
!
    character(len=8) :: crit, result, lpain(4), lpaout(1)
    character(len=8) :: k8b, typarr(5), chli(3)
    character(len=16) :: option, noparr(5)
    character(len=24) :: ligrmo, chgeom, depla, chtime
    character(len=24) :: lchin(4), lchout(1), lisord
    character(len=72) :: rep
!
! ----------------------------------------------------------------------
!
!
    call jemarq()
    c16b = (0.d0, 0.d0)
    lisord = '&&PECHLI.VECTORDR'
    chtime = '&&PECHLI.CH_INST_R'
    f0u = 0
!
!
! -- VERIFICATIONS INITIALES
!
    call getvid(' ', 'RESULTAT', scal=result, nbret=iret)
    call rsexpa(result, 2, 'ETA_PILOTAGE', iret)
!
    if (iret .eq. 0) then
        call utmess('F', 'POSTELEM_3', sk=result)
    end if
!
!
! -- EXISTENCE D'UN CHARGEMENT CONSTANT
    call getvtx('CHAR_LIMITE', 'CHAR_CSTE', iocc=1, scal=rep, nbret=iret)
    chrcst = rep .eq. 'OUI'
!
!
!
!
! -- ECRITURE DE LA TABLE RESULTAT
!
    typarr(1) = 'I'
    typarr(2) = 'R'
    typarr(3) = 'R'
    typarr(4) = 'R'
    typarr(5) = 'R'
!
    noparr(1) = 'NUME_ORDRE'
    noparr(2) = 'INST'
    noparr(3) = 'M'
    noparr(4) = 'CHAR_LIMI_SUP'
!
    if (chrcst) then
        noparr(5) = 'PUIS_CHAR_CSTE'
    else
        noparr(5) = 'CHAR_LIMI_ESTIM'
    end if
!
    call tbcrsd(resu, 'G')
    call tbajpa(resu, 5, noparr, typarr)
!
!
!
    call megeom(modele, chgeom)
    ligrmo = modele//'.MODELE'
!
!
! -- EXTRACTION DES NUMEROS D'ORDRE DU CALCUL
!
    call getvr8(' ', 'PRECISION', scal=prec, nbret=iret)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=iret)
    call rsutnu(result, ' ', 0, lisord, nbord, &
                prec, crit, iret)
    if (iret .ne. 0) then
        call utmess('F', 'POSTELEM_1', sk=result)
    end if
    call jeveuo(lisord, 'L', jord)
!
!
! -- CALCUL DES CHARGES LIMITES AUX DIFFERENTS INSTANTS
!
    do i = 1, nbord
        call jemarq()
        call jerecu('V')
!
!      EXTRACTION DU CHAMP DE DEPLACEMENT
        numord = zi(jord-1+i)
        call rsexch('F', result, 'DEPL', numord, depla, &
                    iret)
!
!
!      CREACTION DE LA CARTE DE L INSTANT DE CALCUL
        call rsadpa(result, 'L', 1, 'INST', numord, &
                    0, sjv=jinst, styp=k8b)
        inst = zr(jinst)
        call mecact('V', chtime, 'MODELE', ligrmo, 'INST_R', &
                    ncmp=1, nomcmp='INST', sr=inst)

!     CALCUL DES VALEUR DE M EN FONCTION DE L'INSTANT
        m = 1+10**(1-inst)
!
!
!      CALCUL DES TERMES ELEMENTAIRES
        lpaout(1) = 'PECHLI'
        lchout(1) = '&&PECHLI'
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PDEPLAR'
        lchin(2) = depla
        lpain(3) = 'PMATERC'
        lchin(3) = mateco
        lpain(4) = 'PINSTR'
        lchin(4) = chtime
        option = 'CHAR_LIMITE'
        call calcul('S', option, ligrmo, 4, lchin, &
                    lpain, 1, lchout, lpaout, 'V', &
                    'OUI')
!
!
!      SOMMATION DE TOUS LES TERMES ELEMENTAIRES
!       CHLIM(1) : CHAR_LIMI_SUP
!       CHLIM(2) : CHAR_LIMI_ESTIM
!       CHLIM(3) : MAX UTILE AU CALCUL DE CHAR_LIMI_ESTIM
!
        call mesomm(lchout(1), 3, vr=chlim)
!
        chli(1) = 'CHLI1'
        chli(2) = 'CHLI2'
        chli(3) = 'CHLI3'
        call memaxm('MAX', lchout(1), 'CHLI3', 3, chli, &
                    chmax, 0, [0])
!
!      CALCUL DU CHARGEMENT PERMANENT SI NECESSAIRE
        if (chrcst) then
            call rsadpa(result, 'L', 1, 'ETA_PILOTAGE', numord, &
                        0, sjv=jpilo, styp=k8b)
            eta = zr(jpilo)
            f0u = m*chlim(2)-eta
            chlim(1) = chlim(1)-f0u
        else
            if (chmax(3) .le. r8miem()) then
                chlim(2) = 0.0d0
            else
                chlim(2) = chlim(2)/chmax(3)
            end if
        end if
!
!
!      ECRITURE DANS LA TABLE RESU DE LA CHARGE LIMITE
        valer(1) = inst
        valer(2) = m
        valer(3) = chlim(1)
        if (chrcst) then
            valer(4) = f0u
        else
            valer(4) = chlim(2)
        end if
        call tbajli(resu, 5, noparr, [numord], valer, &
                    [c16b], k8b, 0)
!
        call jedema()
    end do
!
! --- MENAGE
    call jedetr('&&PECHLI.VECTORDR')
    call detrsd('CARTE', '&&PECHLI.CH_INST_R')
!
    call jedema()
end subroutine
