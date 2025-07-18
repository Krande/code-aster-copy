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

subroutine afdi2d(irep, eta, car, val, jdc, &
                  jdv, ivr, iv, kma, ncmp, &
                  ntp, jdcinf, jdvinf, isym)
!
!
! --------------------------------------------------------------------------------------------------
!
!           Affectation des valeurs des matrices à tous les éléments
!           demandés par l'utilisateur dans les cartes correspondantes
!           les éléments concernés sont les éléments discrets 2D
!
! --------------------------------------------------------------------------------------------------
    implicit none
    integer(kind=8) :: irep, jdv(3), jdc(3), ivr(*), iv, ncmp, ntp, ifm
    integer(kind=8) :: isym, jdcinf, jdvinf
    real(kind=8) :: eta, val(*)
    character(len=1) :: kma(3)
    character(len=*) :: car
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/impmv.h"
#include "asterfort/infdis.h"
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: jj, ll, dimmat, ibid
    real(kind=8) :: r8bid
    aster_logical :: nonsym
    character(len=7) :: ki
    character(len=8) :: k8bid
    character(len=13) :: carbid
! --------------------------------------------------------------------------------------------------
!
    call infdis('DMXM', dimmat, r8bid, k8bid)
    ifm = ivr(4)
!
!   Boucle sur les types de matrice (K,M,A)
    do jj = 1, 3
        if (car(1:1) .ne. kma(jj)) cycle
!       Repère GLOBAL ou repère LOCAL
        zk8(jdcinf+jj-1) = 'REP'//kma(jj)//'    '
        zr(jdvinf+jj-1) = irep
!       Matrice symétrique ou pas
        zk8(jdcinf+jj+2) = 'SYM'//kma(jj)//'    '
        zr(jdvinf+jj+2) = isym
!       Affectation de la matrice
        zk8(jdcinf+jj+5) = 'DIS'//kma(jj)//'    '
        zr(jdvinf+jj+5) = 1.0d0
!       Coefficient amortissement hystérétique
        if (car(1:1) .eq. 'K') then
            zk8(jdcinf+9) = 'ETAK    '
            zr(jdvinf+9) = eta
        end if
        zk8(jdcinf+10) = 'TYDI    '
!
        nonsym = (isym .eq. 2)
        ntp = jj
!
!       NCMP : nombre de composantes de la matrice non-symétrique
        do ll = 1, dimmat
            call codent(ll, 'G', ki)
            zk8(jdc(jj)+ll-1) = kma(jj)//ki
            zr(jdv(jj)+ll-1) = 0.d0
        end do
! ------------------------------------------------------ [K|M|A]_T_D_N
        if (car(3:7) .eq. 'T_D_N') then
            carbid = '2D_DIS_'//car(3:7)
            ASSERT(isym .eq. 1)
            ncmp = 4
            if (car(1:1) .eq. 'M') then
                zr(jdv(jj)) = val(iv)
                zr(jdv(jj)+2) = val(iv)
                iv = iv+1
            else
                zr(jdv(jj)) = val(iv)
                zr(jdv(jj)+2) = val(iv+1)
                iv = iv+2
            end if
! ------------------------------------------------------ [K|M|A]_T_D_L
        else if (car(3:7) .eq. 'T_D_L') then
            carbid = '2D_DIS_'//car(3:7)
            ASSERT(isym .eq. 1)
            ncmp = 16
            if (car(1:1) .eq. 'M') then
                zr(jdv(jj)) = val(iv)
                zr(jdv(jj)+2) = val(iv)
                zr(jdv(jj)+5) = val(iv)
                zr(jdv(jj)+9) = val(iv)
                iv = iv+1
            else
                zr(jdv(jj)) = val(iv)
                zr(jdv(jj)+2) = val(iv+1)
                zr(jdv(jj)+3) = -val(iv)
                zr(jdv(jj)+5) = val(iv)
                zr(jdv(jj)+7) = -val(iv+1)
                zr(jdv(jj)+9) = val(iv+1)
                iv = iv+2
            end if
! ------------------------------------------------------ [K|M|A]_TR_D_N
        else if (car(3:8) .eq. 'TR_D_N') then
            carbid = '2D_DIS_'//car(3:8)
            ASSERT(isym .eq. 1)
            ncmp = 9
            if (car(1:1) .eq. 'M') then
                zr(jdv(jj)) = val(iv)
                zr(jdv(jj)+2) = val(iv)
                zr(jdv(jj)+3) = val(iv)*val(iv+3)
                zr(jdv(jj)+4) = -val(iv)*val(iv+2)
                zr(jdv(jj)+5) = val(iv+1)+val(iv)*(val(iv+2)*val(iv+ &
                                                                 2)+val(iv+3)*val(iv+3))
                iv = iv+4
            else
                zr(jdv(jj)) = val(iv)
                zr(jdv(jj)+2) = val(iv+1)
                zr(jdv(jj)+5) = val(iv+2)
                iv = iv+3
            end if
! ------------------------------------------------------ [K|M|A]_TR_D_L
        else if (car(3:8) .eq. 'TR_D_L') then
            carbid = '2D_DIS_'//car(3:8)
            ASSERT(isym .eq. 1)
            ncmp = 36
            if (car(1:1) .eq. 'M') then
                zr(jdv(jj)) = val(iv)
                zr(jdv(jj)+2) = val(iv)
                zr(jdv(jj)+5) = val(iv+1)
                zr(jdv(jj)+9) = val(iv)
                zr(jdv(jj)+14) = val(iv)
                zr(jdv(jj)+20) = val(iv+1)
                iv = iv+2
            else
                zr(jdv(jj)) = val(iv)
                zr(jdv(jj)+2) = val(iv+1)
                zr(jdv(jj)+5) = val(iv+2)
                zr(jdv(jj)+6) = -val(iv)
                zr(jdv(jj)+9) = val(iv)
                zr(jdv(jj)+11) = -val(iv+1)
                zr(jdv(jj)+14) = val(iv+1)
                zr(jdv(jj)+17) = -val(iv+2)
                zr(jdv(jj)+20) = val(iv+2)
                iv = iv+3
            end if
! ------------------------------------------------------ [K|M|A]_T_N
        else if (car(3:5) .eq. 'T_N') then
            carbid = '2D_DIS_'//car(3:5)
            if (nonsym) then
                ncmp = 4
            else
                ncmp = 3
            end if
            do ll = 1, ncmp
                call codent(ll, 'G', ki)
                zk8(jdc(jj)+ll-1) = kma(jj)//ki
                zr(jdv(jj)+ll-1) = val(iv)
                iv = iv+1
            end do
            ncmp = 4
! ------------------------------------------------------ [K|M|A]_T_L
        else if (car(3:5) .eq. 'T_L') then
            carbid = '2D_DIS_'//car(3:5)
            if (nonsym) then
                ncmp = 16
            else
                ncmp = 10
            end if
            do ll = 1, ncmp
                call codent(ll, 'G', ki)
                zk8(jdc(jj)+ll-1) = kma(jj)//ki
                zr(jdv(jj)+ll-1) = val(iv)
                iv = iv+1
            end do
            ncmp = 16
! ------------------------------------------------------ [K|M|A]_TR_N
        else if (car(3:6) .eq. 'TR_N') then
            carbid = '2D_DIS_'//car(3:6)
            if (nonsym) then
                ncmp = 9
            else
                ncmp = 6
            end if
            do ll = 1, ncmp
                call codent(ll, 'G', ki)
                zk8(jdc(jj)+ll-1) = kma(jj)//ki
                zr(jdv(jj)+ll-1) = val(iv)
                iv = iv+1
            end do
            ncmp = 9
! ------------------------------------------------------ [K|M|A]_TR_L
        else if (car(3:6) .eq. 'TR_L') then
            carbid = '2D_DIS_'//car(3:6)
            if (nonsym) then
                ncmp = 36
            else
                ncmp = 21
            end if
            do ll = 1, ncmp
                call codent(ll, 'G', ki)
                zk8(jdc(jj)+ll-1) = kma(jj)//ki
                zr(jdv(jj)+ll-1) = val(iv)
                iv = iv+1
            end do
            ncmp = 36
        else
            ASSERT(.false.)
        end if
        call infdis('CODE', ibid, zr(jdvinf+10), carbid)
!
        if (ivr(3) .eq. 2) call impmv(ifm, car(1:8), zr(jdv(jj)), ncmp, isym)
    end do
!
end subroutine
