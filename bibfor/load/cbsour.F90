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
subroutine cbsour(load, mesh, model, geomDime, valeType)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/casour.h"
#include "asterfort/copisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/utmess.h"
#include "asterfort/dismoi.h"
#include "asterfort/convertFieldNodeToFieldElga.h"

!
    character(len=8), intent(in) :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
    character(len=4), intent(in) :: valeType
!
    integer(kind=8) :: nbcalc, icalc, nbfac, isour, iocc
    real(kind=8) :: r8bid
    character(len=4) :: tych
    character(len=8) :: scalc
    character(len=16) :: motfac
    character(len=24) :: carte, chsour, chame
!     ------------------------------------------------------------------
!
    motfac = 'SOURCE'
    call getfac(motfac, nbfac)
!
    chame = load//'.CHTH.SOURC'
!
    nbcalc = 0
    if (valeType .eq. 'REEL') then
        do iocc = 1, nbfac
            call getvid(motfac, 'SOUR_CALCULEE', iocc=iocc, scal=chsour, nbret=icalc)
            nbcalc = nbcalc+icalc
        end do
        if (nbcalc .gt. 1) then
            call utmess('F', 'MODELISA3_64')
        else if (nbcalc .eq. 1) then
            call dismoi('TYPE_CHAMP', chsour(1:19), 'CHAMP', repk=tych, arret='F')
            if (tych == "NOEU") then
                call convertFieldNodeToFieldElga(model, chsour(1:19), chame(1:19))
            elseif (tych == "ELGA") then
                call copisd('CHAMP_GD', 'G', chsour(1:19), chame(1:19))
            else
                print *, tych
                ASSERT(ASTER_FALSE)
            end if
        end if
    end if
!
    carte = load//'.CHTH.SOURE'
!
    if (valeType .eq. 'REEL') then
        call getvr8(motfac, 'SOUR', iocc=1, scal=r8bid, nbret=isour)
    else
        call getvid(motfac, 'SOUR', iocc=1, scal=scalc, nbret=isour)
    end if
    if (isour .eq. 1) then
        call casour(load, mesh, model, geomDime, valeType)
    end if
!
end subroutine
