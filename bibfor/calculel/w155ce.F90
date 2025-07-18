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

subroutine w155ce(nomres, resu, nbordr, liordr)
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
! ======================================================================
!     COMMANDE :  POST_CHAMP / COQU_EXCENT
! ----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/exlima.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/rsexch.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
    character(len=8) :: nomres, resu
    integer(kind=8) :: nbordr, liordr(nbordr)
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iret, i, nuordr, ibid, nocc, iocc
    character(len=8) :: modele, carele, mplan
    character(len=8) :: modeav, lpain(2), lpaout(1)
    character(len=4) :: tsca
    character(len=16) :: motfac, nomsym
    character(len=19) :: chin, chextr, ligrel, resu19, lchin(2), lchout(1)
    integer(kind=8) :: vali(2), iexi
    aster_logical :: lreel, lnoeu, ldetli, lvide
!     ------------------------------------------------------------------
!
    call jemarq()
!
!
    call infmaj()
    call infniv(ifm, niv)
    resu19 = resu
!
!
!
!     -- 1. : Y-A-T-IL QUELQUE CHOSE A FAIRE ?
!     ----------------------------------------
    call getfac('COQU_EXCENT', nocc)
    if (nocc .eq. 0) goto 30
    ASSERT(nocc .lt. 10)
!
!
    modeav = ' '
    ldetli = .false.
    lvide = .true.
    do iocc = 1, nocc
!
!     -- 2.  : NOMSYM, MPLAN :
!     --------------------------------------------------
        motfac = 'COQU_EXCENT'
        call getvtx(motfac, 'NOM_CHAM', iocc=iocc, scal=nomsym, nbret=ibid)
        ASSERT(nomsym .eq. 'EFGE_ELNO' .or. nomsym .eq. 'EFGE_ELGA')
        call getvtx(motfac, 'MODI_PLAN', iocc=iocc, scal=mplan, nbret=ibid)
        ASSERT(mplan .eq. 'OUI')
        lnoeu = nomsym .eq. 'EFGE_ELNO'
!
!
!     -- 3. : BOUCLE SUR LES NUMERO D ORDRE
!     --------------------------------------------------
        do i = 1, nbordr
            nuordr = liordr(i)
            call rsexch(' ', resu19, nomsym, nuordr, chin, &
                        iret)
            if (iret .eq. 0) then
!
!         -- 3.1 : MODELE, CARELE, LIGREL :
                call rslesd(resu, nuordr, model_=modele, cara_elem_=carele)
                if (modele .ne. modeav) then
                    if (ldetli) call detrsd('LIGREL', ligrel)
                    call exlima(' ', 1, 'G', modele, ligrel)
                    modeav = modele
!             -- SI ON CREE UN LIGREL, IL FAUT VERIFIER QUE L'ON S'EN
!                SERT VRAIMENT. SINON, IL FAUT LE DETRUIRE:
                    ldetli = .false.
                    if (ligrel(1:8) .ne. modele) ldetli = .true.
                end if
!
                call rsexch(' ', nomres, nomsym, nuordr, chextr, &
                            iret)
                ASSERT(iret .eq. 100)
!
                call jelira(chin//'.CELV', 'TYPE', cval=tsca)
                if (tsca .eq. 'R') then
                    lreel = .true.
                else if (tsca .eq. 'C') then
                    lreel = .false.
                else
                    ASSERT(.false.)
                end if
!
                if (lnoeu) then
                    if (lreel) then
                        lpain(1) = 'PEFFONR'
                        lpaout(1) = 'PEFFOENR'
                    else
                        lpain(1) = 'PEFFONC'
                        lpaout(1) = 'PEFFOENC'
                    end if
                else
                    if (lreel) then
                        lpain(1) = 'PEFFOGR'
                        lpaout(1) = 'PEFFOEGR'
                    else
                        lpain(1) = 'PEFFOGC'
                        lpaout(1) = 'PEFFOEGC'
                    end if
                end if
!
                lchin(1) = chin
                lchout(1) = chextr
!
                lpain(2) = 'PCACOQU'
                lchin(2) = carele//'.CARCOQUE'
!
                call calcul('C', 'EFGE_EXCENT', ligrel, 2, lchin, &
                            lpain, 1, lchout, lpaout, 'G', &
                            'OUI')
!
                call jeexin(lchout(1)//'.CELV', iexi)
                if (iexi .eq. 0) then
                    vali(1) = iocc
                    vali(2) = nuordr
                    call utmess('A', 'CALCULEL2_19', ni=2, vali=vali)
                else
                    ldetli = .false.
                    lvide = .false.
                    call rsnoch(nomres, nomsym, nuordr)
                end if
            end if
        end do
    end do
!
    if (ldetli) call detrsd('LIGREL', ligrel)
    if (lvide) then
        call utmess('F', 'CALCULEL2_20')
    end if
!
!
30  continue
    call jedema()
end subroutine
