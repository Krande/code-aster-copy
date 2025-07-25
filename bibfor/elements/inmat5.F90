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
subroutine inmat5(elrefa, nno, nnos, npg, mganos, mgano2)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
!
    character(len=8), intent(in) :: elrefa
    integer(kind=8), intent(in) :: nnos, npg, nno
    real(kind=8), intent(in) :: mganos(MT_NBPGMX, MT_NNOMAX)
    real(kind=8), intent(out) :: mgano2(MT_NBPGMX, MT_NNOMAX)
!
! --------------------------------------------------------------------------------------------------
!
! Finite elements management
!
! Compute the Gauss passage matrix at all nodes
!
! --------------------------------------------------------------------------------------------------
!
! In  elrefa           : name of geometric support for finite element
! In  nno              : number of nodes
! In  nnos             : number of middle nodes
! In  npg              : number of Gauss points
! In  mganos           : Gauss passage matrix at vertex nodes
! Out mgano2           : Gauss passage matrix at all nodes
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: kpg, kno, knos, k
    real(kind=8) :: nosom(MT_NNOMAX, MT_NNOMAX)
    real(kind=8), parameter :: demi = 0.5d0, quart = 0.25d0
!
! --------------------------------------------------------------------------------------------------
!
!
!     -- SI NNO=NNOS, IL N'Y A QU'A COPIER :
!
    if (nnos .eq. nno) then
        do kpg = 1, npg
            do kno = 1, nno
                mgano2(kpg, kno) = mganos(kpg, kno)
            end do
        end do
        goto 100
    end if
!
!
!     1) CALCUL DE NOSOM :
!     ---------------------
    do kno = 1, nno
        do knos = 1, nnos
            nosom(kno, knos) = 0.d0
        end do
    end do
!
!     1.1) LES NOEUDS SOMMETS SONT TOUJOURS LES 1ERS :
    do knos = 1, nnos
        nosom(knos, knos) = 1.d0
    end do
!
!
!     1.2) LES NOEUDS MILIEUX SE DEDUISENT DES SOMMETS :
    if ((elrefa .eq. 'H20') .or. (elrefa .eq. 'H27')) then
        ASSERT(nnos .eq. 8)
        nosom(9, 1) = demi
        nosom(9, 2) = demi
        nosom(10, 2) = demi
        nosom(10, 3) = demi
        nosom(11, 3) = demi
        nosom(11, 4) = demi
        nosom(12, 1) = demi
        nosom(12, 4) = demi
        nosom(13, 1) = demi
        nosom(13, 5) = demi
        nosom(14, 2) = demi
        nosom(14, 6) = demi
        nosom(15, 3) = demi
        nosom(15, 7) = demi
        nosom(16, 4) = demi
        nosom(16, 8) = demi
        nosom(17, 5) = demi
        nosom(17, 6) = demi
        nosom(18, 6) = demi
        nosom(18, 7) = demi
        nosom(19, 7) = demi
        nosom(19, 8) = demi
        nosom(20, 5) = demi
        nosom(20, 8) = demi
!
        if (elrefa .eq. 'H27') then
            nosom(21, 1) = quart
            nosom(21, 2) = quart
            nosom(21, 3) = quart
            nosom(21, 4) = quart
!
            nosom(22, 1) = quart
            nosom(22, 2) = quart
            nosom(22, 5) = quart
            nosom(22, 6) = quart
!
            nosom(23, 2) = quart
            nosom(23, 3) = quart
            nosom(23, 6) = quart
            nosom(23, 7) = quart
!
            nosom(24, 3) = quart
            nosom(24, 4) = quart
            nosom(24, 7) = quart
            nosom(24, 8) = quart
!
            nosom(25, 1) = quart
            nosom(25, 4) = quart
            nosom(25, 5) = quart
            nosom(25, 8) = quart
!
            nosom(26, 5) = quart
            nosom(26, 6) = quart
            nosom(26, 7) = quart
            nosom(26, 8) = quart
!
            do k = 1, 8
                nosom(MT_NNOMAX, k) = demi/4.d0
            end do
        end if
!
!
    else if ((elrefa .eq. 'P15') .or. (elrefa .eq. 'P18') .or. (elrefa .eq. 'P21')) then
        ASSERT(nnos .eq. 6)
        nosom(7, 1) = demi
        nosom(7, 2) = demi
        nosom(8, 2) = demi
        nosom(8, 3) = demi
        nosom(9, 1) = demi
        nosom(9, 3) = demi
        nosom(10, 1) = demi
        nosom(10, 4) = demi
        nosom(11, 2) = demi
        nosom(11, 5) = demi
        nosom(12, 3) = demi
        nosom(12, 6) = demi
        nosom(13, 4) = demi
        nosom(13, 5) = demi
        nosom(14, 5) = demi
        nosom(14, 6) = demi
        nosom(15, 4) = demi
        nosom(15, 6) = demi
!
        if (elrefa .eq. 'P18' .or. (elrefa .eq. 'P21')) then
!
            nosom(16, 2) = quart
            nosom(16, 1) = quart
            nosom(16, 4) = quart
            nosom(16, 5) = quart
!
            nosom(17, 2) = quart
            nosom(17, 5) = quart
            nosom(17, 6) = quart
            nosom(17, 3) = quart
!
            nosom(18, 1) = quart
            nosom(18, 3) = quart
            nosom(18, 6) = quart
            nosom(18, 4) = quart
!
        end if
!
    else if (elrefa .eq. 'T10' .or. elrefa .eq. 'T15') then
        ASSERT(nnos .eq. 4)
        nosom(5, 1) = demi
        nosom(5, 2) = demi
        nosom(6, 2) = demi
        nosom(6, 3) = demi
        nosom(7, 1) = demi
        nosom(7, 3) = demi
        nosom(8, 1) = demi
        nosom(8, 4) = demi
        nosom(9, 2) = demi
        nosom(9, 4) = demi
        nosom(10, 3) = demi
        nosom(10, 4) = demi
!
!
    else if (elrefa .eq. 'P13' .or. elrefa .eq. 'P19') then
        ASSERT(nnos .eq. 5)
        nosom(6, 1) = demi
        nosom(6, 2) = demi
        nosom(7, 2) = demi
        nosom(7, 3) = demi
        nosom(8, 3) = demi
        nosom(8, 4) = demi
        nosom(9, 1) = demi
        nosom(9, 4) = demi
        nosom(10, 1) = demi
        nosom(10, 5) = demi
        nosom(11, 2) = demi
        nosom(11, 5) = demi
        nosom(12, 3) = demi
        nosom(12, 5) = demi
        nosom(13, 4) = demi
        nosom(13, 5) = demi
!
!
    else if (elrefa .eq. 'TR6') then
        ASSERT(nnos .eq. 3)
        nosom(4, 1) = demi
        nosom(4, 2) = demi
        nosom(5, 2) = demi
        nosom(5, 3) = demi
        nosom(6, 3) = demi
        nosom(6, 1) = demi
!
!
    else if (elrefa .eq. 'TR7') then
        ASSERT(nnos .eq. 3)
        nosom(4, 1) = demi
        nosom(4, 2) = demi
        nosom(5, 2) = demi
        nosom(5, 3) = demi
        nosom(6, 3) = demi
        nosom(6, 1) = demi
        nosom(7, 1) = quart
        nosom(7, 2) = quart
        nosom(7, 3) = quart
!
!
    else if (elrefa .eq. 'QU8') then
        ASSERT(nnos .eq. 4)
        nosom(5, 1) = demi
        nosom(5, 2) = demi
        nosom(6, 2) = demi
        nosom(6, 3) = demi
        nosom(7, 3) = demi
        nosom(7, 4) = demi
        nosom(8, 4) = demi
        nosom(8, 1) = demi
!
!
    else if (elrefa .eq. 'QU9') then
        ASSERT(nnos .eq. 4)
        nosom(5, 1) = demi
        nosom(5, 2) = demi
        nosom(6, 2) = demi
        nosom(6, 3) = demi
        nosom(7, 3) = demi
        nosom(7, 4) = demi
        nosom(8, 4) = demi
        nosom(8, 1) = demi
        nosom(9, 1) = quart
        nosom(9, 2) = quart
        nosom(9, 3) = quart
        nosom(9, 4) = quart
!
!
    else if (elrefa .eq. 'SE3') then
        ASSERT(nnos .eq. 2)
        nosom(3, 1) = demi
        nosom(3, 2) = demi
!
!
    else if (elrefa .eq. 'SE4') then
        ASSERT(nnos .eq. 2)
        nosom(3, 1) = 2.d0/3.d0
        nosom(3, 2) = 1.d0/3.d0
        nosom(4, 1) = 1.d0/3.d0
        nosom(4, 2) = 2.d0/3.d0
!
    elseif (elrefa .eq. 'HE9') then
!    Pas besoin de faire un passage 9eme noeud --> Noeud SOMMET
!    Le neuvieme noeud sert uniquement au calcul du pincement

    elseif (elrefa .eq. 'PE7') then
!    Pas besoin de faire un passage 7eme noeud --> Noeud SOMMET
!    Le septième noeud sert uniquement au calcul du pincement

    else
        ASSERT(ASTER_FALSE)
    end if
!
!
!     2) ON MULTIPLIE : MGANO2=MGANOS * NOSOM :
!     ----------------------------------
    do kno = 1, nno
        do kpg = 1, npg
            do knos = 1, nnos
                mgano2(kpg, kno) = mgano2(kpg, kno)+mganos(kpg, knos)*nosom(kno, knos)
            end do
        end do
    end do
!
100 continue
!
end subroutine
