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
subroutine vpnor2(nomcon, nbmode, numord, coef)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/peenc2.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: nbmode, numord(*)
    real(kind=8) :: coef(*)
    character(len=*) :: nomcon
!     NORMALISATION DE TOUS LES CHAMPS D'UN MODE_MECA
!
! IN  NOMCON : NOM DU CONCEPT RESULTAT DE TYPE MODE_MECA
! IN  NBMODE : NOMBRE DE MODES
! IN  NUMORD : NUMERO D'ORDRE
! IN  COEF  : COEFFICIENT REEL A APLLIQUER AUX CHAMPS
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) :: ibid, nbnosy, isy, im, iordr, iret, lvale, neq, ieq
    real(kind=8) :: rcoef
    character(len=8) :: typmod
    character(len=16) :: nomsym
    character(len=19) :: nomd2
    character(len=24) :: vale
!     ------------------------------------------------------------------
    data vale/'                   .VALE'/
!     ------------------------------------------------------------------
!
    call jemarq()
!
    nomd2 = nomcon
!
    call jelira(nomd2//'.DESC', 'NOMMAX', nbnosy)
    if (nbnosy .eq. 0) goto 999
!
    do isy = 1, nbnosy
        call jenuno(jexnum(nomd2//'.DESC', isy), nomsym)
        do im = 1, nbmode
            iordr = numord(im)
            call rsexch(' ', nomcon, nomsym, iordr, vale(1:19), &
                        iret)
            if (iret .eq. 0) then
                call jeexin(vale(1:19)//'.VALE', ibid)
                if (ibid .gt. 0) then
                    vale = vale(1:19)//'.VALE'
                else
                    vale = vale(1:19)//'.CELV'
                end if
!
                call jelira(vale, 'TYPE', cval=typmod)
                if (nomsym(1:4) .eq. 'EFGE' .or. nomsym(1:4) .eq. 'SIGM' .or. nomsym(1:4) &
                    .eq. 'EPSI' .or. nomsym(1:4) .eq. 'SIEF' .or. nomsym(1:4) .eq. 'FORC' &
                    .or. nomsym(1:4) .eq. 'REAC' .or. nomsym(1:4) .eq. 'DEGE') then
                    rcoef = coef(im)
                else if (nomsym(1:4) .eq. 'EQUI') then
                    call utmess('A', 'UTILITAI5_88', sk=nomsym)
                    goto 12
                elseif (nomsym(1:4) .eq. 'EPOT' .or. nomsym(1:4) &
                        .eq. 'ECIN') then
                    rcoef = coef(im)*coef(im)
                    if (typmod(1:1) .eq. 'R') then
                        call peenc2(vale(1:19), rcoef)
                    else
                        call utmess('F', 'UTILITAI5_89')
                    end if
                    goto 12
                else
                    goto 12
                end if
                call jeveuo(vale, 'E', lvale)
                call jelira(vale, 'LONMAX', neq)
                if (typmod(1:1) .eq. 'R') then
                    do ieq = 0, neq-1
                        zr(lvale+ieq) = zr(lvale+ieq)*rcoef
                    end do
                else
                    call utmess('F', 'UTILITAI5_89')
                end if
            end if
12          continue
        end do
    end do
!
!
999 continue
    call jedema()
end subroutine
