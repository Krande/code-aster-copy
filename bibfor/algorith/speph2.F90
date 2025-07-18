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
subroutine speph2(movrep, napexc, nbmode, nbpf, intmod, &
                  table, specmr, specmi)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
    integer(kind=8) :: napexc, nbmode, nbpf
    real(kind=8) :: specmr(nbpf, *), specmi(nbpf, *)
    aster_logical :: intmod
    character(len=8) :: table
    character(len=16) :: movrep
!
    integer(kind=8) :: ival(2), ideb1, ifin1, i, j, imi, imj, ideb, isj, ifon, if1
    integer(kind=8) :: mxval, lnumi, lnumj, i1
!
    character(len=24) :: chnumi, chnumj, chvale
!
!     ------------------------------------------------------------------
!
    if (movrep .eq. 'ABSOLU') then
        ideb1 = 1
        ifin1 = napexc+nbmode
    else if (movrep .eq. 'RELATIF') then
        ideb1 = napexc+1
        ifin1 = napexc+nbmode
    else if (movrep .eq. 'DIFFERENTIEL') then
        ideb1 = 1
        ifin1 = napexc
    end if
!
    chnumi = table//'.NUMI'
    chnumj = table//'.NUMJ'
    chvale = table//'.VALE'
    call jeveuo(chnumi, 'L', lnumi)
    call jeveuo(chnumj, 'L', lnumj)
    call jelira(chnumi, 'LONMAX', mxval)
!
    j = 0
    do imj = ideb1, ifin1
        j = j+1
!
        ival(2) = imj
!
        ideb = imj
        if (intmod) ideb = ideb1
!
        i = 0
        do imi = ideb, imj
            i = i+1
!
            ival(1) = imi
!
            do i1 = 1, mxval
                if ((zi(lnumi-1+i1) .eq. ival(1)) .and. (zi(lnumj-1+i1) .eq. ival(2))) then
                    call jeveuo(jexnum(chvale, i1), 'L', ifon)
                end if
            end do
!
            isj = j*(j-1)/2+i
!
            do if1 = 1, nbpf
                if (ival(1) .eq. ival(2)) then
                    specmr(if1, isj) = zr(ifon-1+if1)
                    specmi(if1, isj) = 0.d0
                else
                    specmr(if1, isj) = zr(ifon+(if1-1)*2)
                    specmi(if1, isj) = zr(ifon+(if1-1)*2+1)
                end if
            end do
        end do
!
    end do
end subroutine
