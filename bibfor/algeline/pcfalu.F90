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
subroutine pcfalu(n, inc, ipc, inx, ipx, &
                  index, imp)
    implicit none
    integer(kind=8) :: n
    integer(kind=8) :: inc(n), index(*)
    integer(kind=4) :: ipc(*), ipx(*)
    integer(kind=8) :: inx(0:n)
!-----------------------------------------------------------------------
! FONCTION  INC;IPC ------------>INX;IPX
!           SYMETRIQUE           PLEINE FORMAT LU
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ii, ii1, ii2, imp, j, kkk
!-----------------------------------------------------------------------
!
!
!  CALCUL DU NOMBRE DE COEFS PAR LIGNE COMPLETE
!        SANS LA DIAGONALE
! ----------------------------------------------
!
!  DEMI-LIGNE
    index(1) = 0
    do i = 2, n
        index(i) = inc(i)-inc(i-1)-1
    end do
!  SYMETRIQUES
    do i = 2, n
        ii1 = inc(i-1)+1
        ii2 = inc(i)
        do ii = ii1, ii2-1
            j = ipc(ii)
            index(j) = index(j)+1
        end do
    end do
!
    if (imp .eq. 1) then
        ii1 = index(1)
        ii2 = index(1)
        do i = 2, n
            ii1 = min(ii1, index(i))
            ii2 = max(ii2, index(i))
        end do
    end if
!
! CALCUL NOUVEAU INX POINTEUR DEBUT DE LIGNE
!                             QUI SE DECALE ENSUITE
! --------------------------------------------------
!
    inx(0) = 0
    inx(1) = 1
    do i = 2, n
        inx(i) = inx(i-1)+index(i-1)
    end do
!
!
! CALCUL NOUVEAUX IPX
!    INX(I)= DEBUT LIGNE I DE L
!    INDEX(I)=POINTEUR DERNIER COEF LIGNE I DE U
! ------------------------------------------------
!
    index(1) = 0
    do i = 2, n
        ii1 = inc(i-1)+1
        ii2 = inc(i)
        kkk = inx(i)-1
        do ii = ii1, ii2-1
            kkk = kkk+1
            ipx(kkk) = ipc(ii)
            j = ipc(ii)
            index(j) = index(j)+1
            ipx(index(j)) = int(i, 4)
        end do
        index(i) = kkk
    end do
!
!
!  TRANSFORMATION DE INX DE DEBUT EN FIN DE LIGNE DE LU
!  ----------------------------------------------------
    inx(0) = 0
    do i = 1, n-1
        inx(i) = inx(i+1)-1
    end do
    inx(n) = index(n)
!
    if (imp .eq. 1) then
        do i = 1, 5
            ii1 = inx(i-1)+1
            ii2 = inx(i)
        end do
    end if
end subroutine
