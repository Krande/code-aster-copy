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
subroutine te0113(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/get_value_mode_local.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"

    character(len=16) :: option, nomte
!
!     BUT:
!       OPTION : 'ROCH_ELNO'
!       DISPOSER DE PARAMETRES MATERIAU ET DE CARACTERISTIQUES DE POUTRE
!       DANS UN CHAMP ELNO
!       ELEMENTS POU_D_T
! ---------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbparaPR = 6
    integer(kind=8) :: icodre(nbparaPR), iret
    integer(kind=8) :: imate, nno, nbcmp, ino, ival, ndim, ino2

!
    real(kind=8) :: valres(nbparaPR), young, k, nexpo, rp02_min, rm_min
    real(kind=8) :: rp02_moy, coef
    real(kind=8) :: caragene(4), carageo(4)
    character(len=8)  :: valp(4)
    character(len=16) :: nomres(nbparaPR)
!
! ----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno)
    ASSERT(nno .eq. 2)

    call jevech('PMATERC', 'L', imate)
    call jevech('PROCHRR', 'E', ival)

    nbcmp = 1+nbparaPR+7

!   parametres elastique
    nomres(1) = 'E'
    call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                ' ', 'ELAS', 0, '', [0.d0], &
                1, nomres, valres, icodre, 1)
    young = valres(1)

!   parametres de materiau POST_ROCHE
    nomres(1) = 'RAMB_OSGO_FACT'
    nomres(2) = 'RAMB_OSGO_EXPO'
    nomres(3) = 'RP02_MIN'
    nomres(4) = 'RM_MIN'
    nomres(5) = 'RP02_MOY'
    nomres(6) = 'COEF'
    call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                ' ', 'POST_ROCHE', 0, '', [0.d0], &
                nbparaPR, nomres, valres, icodre, 0)
    if (icodre(1) .ne. 0) call utmess('F', 'POSTROCHE_16')
    k = valres(1)
    if (k .lt. 0.d0) call utmess('F', 'POSTROCHE_19', sk=nomres(1), sr=k)
    nexpo = valres(2)
    if (nexpo .lt. 0.d0) call utmess('F', 'POSTROCHE_19', sk=nomres(2), sr=nexpo)

!   pour les paramètres facultatifs, on les mets à une valeur négative si absent
!   afin d'émettre les messages d'erreur dans POST_ROCHE s'ils étaient nécessaires
!   voir issue30703
    if (icodre(3) .ne. 0) then
        rp02_min = -1d0
    else
        if (valres(3) .lt. 0.d0) call utmess('F', 'POSTROCHE_19', sk=nomres(3), sr=valres(3))
        rp02_min = valres(3)
    end if
    if (icodre(4) .ne. 0) then
        rm_min = -1d0
    else
        if (valres(4) .lt. 0.d0) call utmess('F', 'POSTROCHE_19', sk=nomres(4), sr=valres(4))
        rm_min = valres(4)
    end if

    if (icodre(3) .eq. 0 .and. icodre(5) .ne. 0) then
        rp02_moy = 1.25d0*rp02_min
    elseif (icodre(5) .ne. 0) then
        rp02_moy = -1d0
    else
        if (valres(5) .lt. 0.d0) call utmess('F', 'POSTROCHE_19', sk=nomres(5), sr=valres(5))
        rp02_moy = valres(5)
    end if

    coef = valres(6)

!   caractéristiques de poutre
    valp = ['A1 ', 'IY1', 'A2 ', 'IY2']
    call get_value_mode_local('PCAGNPO', valp, caragene, iret)
    valp = ['R1 ', 'EP1', 'R2 ', 'EP2']
    call get_value_mode_local('PCAGEPO', valp, carageo, iret)

    do ino = 1, nno
        if (ino .eq. 1) then
            ino2 = 2
        else
            ino2 = 1
        end if
        zr(ival+(ino-1)*nbcmp-1+1) = young
        zr(ival+(ino-1)*nbcmp-1+2) = k
        zr(ival+(ino-1)*nbcmp-1+3) = nexpo
        zr(ival+(ino-1)*nbcmp-1+4) = rp02_min
        zr(ival+(ino-1)*nbcmp-1+5) = rm_min
        zr(ival+(ino-1)*nbcmp-1+6) = rp02_moy
        zr(ival+(ino-1)*nbcmp-1+7) = coef
!       A
        zr(ival+(ino-1)*nbcmp-1+nbparaPR+2) = caragene(2*(ino-1)+1)
!       I
        zr(ival+(ino-1)*nbcmp-1+nbparaPR+3) = caragene(2*(ino-1)+2)
!       R
        zr(ival+(ino-1)*nbcmp-1+nbparaPR+4) = carageo(2*(ino-1)+1)
!       EP
        zr(ival+(ino-1)*nbcmp-1+nbparaPR+5) = carageo(2*(ino-1)+2)
!       valeur de l'autre noeud pour réduction
!       I2
        zr(ival+(ino-1)*nbcmp-1+nbparaPR+6) = caragene(2*(ino2-1)+2)
!       R2
        zr(ival+(ino-1)*nbcmp-1+nbparaPR+7) = carageo(2*(ino2-1)+1)
!       EP
        zr(ival+(ino-1)*nbcmp-1+nbparaPR+8) = carageo(2*(ino2-1)+2)
    end do
! ----------------------------------------------------------------------
end subroutine
