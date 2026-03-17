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
subroutine nmjalo(sddisc, timeCurr, prec, jalon, lOutOFTime)
!
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/compr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc
    real(kind=8), intent(in) :: timeCurr, prec
    real(kind=8), intent(out) :: jalon
    aster_logical, intent(out) :: lOutOFTime
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (UTILITAIRE)
!
! PROCHAIN INSTANT DE PASSAGE DANS LA LISTE DES JALONS
!
! --------------------------------------------------------------------------------------------------
!
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  INST   : INSTANT RECHERCHE
! IN  PREC   : PRECISION
! OUT JALON  : VALEUR DE L'INSTANT JALON TROUVE
!              VAUT R8VIDE SI L'INSTANT EST AU DELA DE LA LSITE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: tpsipo
    integer(kind=8) :: jipo
    integer(kind=8) :: ipo, nipo
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lOutOFTime = ASTER_TRUE
    jalon = r8vide()

! - Get access
    tpsipo = sddisc(1:19)//'.LIPO'
    call jelira(tpsipo, 'LONMAX', ival=nipo)
    call jeveuo(tpsipo, 'L', jipo)

! - RECHERCHE PROCHAIN JALON
    do ipo = 1, nipo
        if (compr8(zr(jipo-1+ipo), 'GT', timeCurr, prec, 1)) then
            jalon = zr(jipo-1+ipo)
            lOutOFTime = ASTER_FALSE
            goto 20
        end if
    end do
20  continue
!
    call jedema()
end subroutine
