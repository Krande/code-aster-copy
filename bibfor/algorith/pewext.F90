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
subroutine pewext(resu)
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsdot.h"
#include "asterfort/cnsreddepl.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: resu
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, iret, jinst
    integer(kind=8) :: nbord, jord, numord
    real(kind=8) :: inst, prec, f0u0, f1u0, f0u1, f1u1, w, valer(3)
    complex(kind=8) :: c16b
    character(len=8) :: crit, result
    character(len=8) :: k8b, typarr(4)
    character(len=16) :: noparr(4)
    character(len=19) :: depla1, force1
    character(len=19) :: depls0, depls1, forcs0, forcs1
    character(len=24) :: lisord
    integer(kind=8) :: ier
!
!-----------------------------------------------------------------------
!
!
    call jemarq()
    c16b = (0.d0, 0.d0)
    lisord = '&&PEWEXT.VECTORDR'
    call getvid(' ', 'RESULTAT', scal=result, nbret=iret)
!
!
! -- INITIALISATION DE LA TABLE RESULTAT
!
    typarr(1) = 'I'
    typarr(2) = 'R'
    typarr(3) = 'R'
    typarr(4) = 'R'
!
    noparr(1) = 'NUME_ORDRE'
    noparr(2) = 'INST'
    noparr(3) = 'TRAV_ELAS'
    noparr(4) = 'TRAV_REEL'
!
    call tbcrsd(resu, 'G')
    call tbajpa(resu, 4, noparr, typarr)
!
!
!
! -- EXTRACTION DES NUMEROS D'ORDRE DU CALCUL
!
    call getvr8(' ', 'PRECISION', scal=prec, nbret=iret)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=iret)
    call rsutnu(result, ' ', 0, lisord, nbord, &
                prec, crit, iret)
    if (iret .ne. 0) then
        call utmess('F', 'POSTELEM_11', sk=result)
    end if
    call jeveuo(lisord, 'L', jord)
!
!
! -- CALCUL DU TRAVAIL DES FORCES EXTERIEURES AUX DIFFERENTS INSTANTS
!
    depls0 = '&&PEWEXT.DEPLS0'
    depls1 = '&&PEWEXT.DEPLS1'
    forcs0 = '&&PEWEXT.FORCS0'
    forcs1 = '&&PEWEXT.FORCS1'
!
    do i = 1, nbord
        call jemarq()
        call jerecu('V')
        numord = zi(jord-1+i)
!
!       EXTRACTION DE L'INSTANT DE CALCUL
        call rsadpa(result, 'L', 1, 'INST', numord, &
                    0, sjv=jinst, styp=k8b)
        inst = zr(jinst)
!
!       EXTRACTION DU CHAMP DE DEPLCAMENT
        call rsexch('F', result, 'DEPL', numord, depla1, &
                    iret)
!
!       -- TOUS LES CHAMPS DE LA SD_RESULTAT N'ONT PAS FORCEMENT
!          LA MEME NUMEROTATION, C'EST POURQUOI ON PASSE PAR DES
!          CHAMPS SIMPLES :
        call cnocns(depla1, 'V', depls1)
!       suppression des composantes non concernées
        call cnsreddepl(depls1)
!
!       EXTRACTION DU CHAMP DE FORCE NODALE
        call rsexch('F', result, 'FORC_NODA', numord, force1, &
                    iret)
        call cnocns(force1, 'V', forcs1)
!       suppression des composantes non concernées
        call cnsreddepl(forcs1)
!
!       CALCUL DU PRODUIT SCALAIRE F.U
        call cnsdot(depls1, forcs1, f1u1, ier)
        ASSERT(ier .eq. 0)
!
!       CALCUL DE L'INTEGRALE I(F.DU)
        if (i .ge. 2) then
            call cnsdot(depls0, forcs1, f1u0, ier)
            ASSERT(ier .eq. 0)
            call cnsdot(depls1, forcs0, f0u1, ier)
            ASSERT(ier .eq. 0)
            w = w+0.5d0*(f0u1-f1u0+f1u1-f0u0)
        else
            w = 0
        end if
!
        valer(1) = inst
        valer(2) = f1u1/2
        valer(3) = w
        call tbajli(resu, 4, noparr, [numord], valer, &
                    [c16b], k8b, 0)
!
        call copisd('CHAM_NO_S', 'V', depls1, depls0)
        call copisd('CHAM_NO_S', 'V', forcs1, forcs0)
        f0u0 = f1u1
!
        call jedema()
    end do
!
    call detrsd('CHAM_NO_S', depls1)
    call detrsd('CHAM_NO_S', depls0)
    call detrsd('CHAM_NO_S', forcs1)
    call detrsd('CHAM_NO_S', forcs0)
!
    call jedema()
end subroutine
