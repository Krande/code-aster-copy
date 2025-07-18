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

subroutine tresu_print(refer, legend, llab, nbref, rela, &
                       tole, ssigne, refr, valr, refi, &
                       vali, refc, valc, ignore, compare)
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/bool_to_int.h"
#include "asterc/testresu_print.h"
!
!
! person_in_charge: mathieu.courtois at edf.fr
!
    character(len=16), intent(in) :: refer
    character(len=16), intent(in) :: legend
    aster_logical, intent(in) :: llab
    integer(kind=8), intent(in) :: nbref
    character(len=*), intent(in) :: rela
    real(kind=8), intent(in) :: tole
    character(len=*), intent(in), optional :: ssigne
    real(kind=8), intent(in), optional :: refr(nbref)
    real(kind=8), intent(in), optional :: valr
    integer(kind=8), intent(in), optional :: refi(nbref)
    integer(kind=8), intent(in), optional :: vali
    complex(kind=8), intent(in), optional :: refc(nbref)
    complex(kind=8), intent(in), optional :: valc
    aster_logical, intent(in), optional :: ignore
    real(kind=8), intent(in), optional :: compare
!
!   Interface d'appel à la fonction d'impression en C/Python pour les TEST_RESU
!   Quand plusieurs valeurs de référence sont fournis, on conserve la plus proche
!   de la valeur calculée.
!
    real(kind=8) :: arefr
    real(kind=8) :: avalr, minvr, tmpr, minvc, tmpc
    integer(kind=8) :: arefi
    integer(kind=8) :: avali, minvi, tmpi
    integer(kind=8) :: i, imin
    complex(kind=8) :: arefc
    complex(kind=8) :: avalc
    real(kind=8) :: arg_cmp
    aster_logical :: skip, isrela, valabs
    integer(kind=8) :: typ
!
    valabs = .false.
    if (present(ssigne)) then
        valabs = ssigne .eq. 'OUI'
    end if
!
    typ = 0
    ASSERT(UN_PARMI3(refr, refi, refc))
    ASSERT(UN_PARMI3(valr, vali, valc))
    ASSERT(ENSEMBLE2(refr, valr))
    ASSERT(ENSEMBLE2(refi, vali))
    ASSERT(ENSEMBLE2(refc, valc))
!
    arefr = 0.d0
    avalr = 0.d0
    if (present(refr)) then
        typ = 1
        avalr = valr
        arefr = refr(1)
        if (valabs) then
            avalr = abs(avalr)
            arefr = abs(arefr)
        end if
        minvr = abs(avalr-arefr)
        imin = 1
        do i = 1, nbref-1
            arefr = refr(i+1)
            if (valabs) then
                arefr = abs(arefr)
            end if
            tmpr = abs(avalr-arefr)
            if (tmpr .lt. minvr) then
                tmpr = minvr
                imin = i+1
            end if
        end do
        arefr = refr(imin)
        if (valabs) then
            arefr = abs(arefr)
        end if
    end if
!
    arefi = 0
    avali = 0
    if (present(refi)) then
        typ = 2
        avali = vali
        arefi = refi(1)
        if (valabs) then
            avali = abs(avali)
            arefi = abs(arefi)
        end if
        minvi = abs(avali-arefi)
        imin = 1
        do i = 1, nbref-1
            arefi = refi(i+1)
            if (valabs) then
                arefi = abs(arefi)
            end if
            tmpi = abs(avali-arefi)
            if (tmpi .lt. minvi) then
                tmpi = minvi
                imin = i+1
            end if
        end do
        arefi = refi(imin)
        if (valabs) then
            arefi = abs(arefi)
        end if
    end if
!
    arefc = dcmplx(0.d0, 0.d0)
    avalc = dcmplx(0.d0, 0.d0)
    if (present(refc)) then
        typ = 3
        avalc = valc
        arefc = refc(1)
        if (valabs) then
            avalc = abs(avalc)
            arefc = abs(arefc)
        end if
        minvc = abs(avalc-arefc)
        imin = 1
        do i = 1, nbref-1
            arefc = refc(i+1)
            if (valabs) then
                arefc = abs(arefc)
            end if
            tmpc = abs(avalc-arefc)
            if (tmpc .lt. minvc) then
                tmpc = minvc
                imin = i+1
            end if
        end do
        arefc = refc(imin)
        if (valabs) then
            arefc = abs(arefc)
        end if
    end if
!
    ASSERT(typ .ge. 1 .and. typ .le. 3)
!
    isrela = rela(1:4) .eq. 'RELA'
    skip = .false.
    if (present(ignore)) then
        skip = ignore
    end if
!
    arg_cmp = 1.d0
    if (present(compare)) then
        arg_cmp = compare
    end if
!
    call testresu_print(refer, legend, bool_to_int(llab), bool_to_int(skip), &
                        bool_to_int(isrela), &
                        tole, typ, arefr, avalr, arefi, &
                        avali, arefc, avalc, arg_cmp)
!
end subroutine tresu_print
