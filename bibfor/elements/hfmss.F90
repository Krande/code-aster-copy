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

subroutine hfmss(ind, vectt, hsfm, hss)
    implicit none
    integer(kind=8) :: ind
    real(kind=8) :: vectt(3, 3)
    real(kind=8) :: hfm(3, 6), hs(2, 6), hsfm(3, 9), hss(2, 9)
!
!     CONSTRUCTION DE  HFM  :  (3,6) AUX X PTS D'INTEGRATION
!                                        X= REDUIT OU NORMAL
!     (DANS L'EXPRESSION DE HFM, PAS DE COMPOSANTES VECTT(3,K))
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    hfm(1, 1) = vectt(1, 1)*vectt(1, 1)
    hfm(1, 2) = vectt(1, 2)*vectt(1, 2)
    hfm(1, 3) = vectt(1, 3)*vectt(1, 3)
    hfm(1, 4) = vectt(1, 1)*vectt(1, 2)
    hfm(1, 5) = vectt(1, 1)*vectt(1, 3)
    hfm(1, 6) = vectt(1, 2)*vectt(1, 3)
!
    hfm(2, 1) = vectt(2, 1)*vectt(2, 1)
    hfm(2, 2) = vectt(2, 2)*vectt(2, 2)
    hfm(2, 3) = vectt(2, 3)*vectt(2, 3)
    hfm(2, 4) = vectt(2, 1)*vectt(2, 2)
    hfm(2, 5) = vectt(2, 1)*vectt(2, 3)
    hfm(2, 6) = vectt(2, 2)*vectt(2, 3)
!
    hfm(3, 1) = 2*vectt(1, 1)*vectt(2, 1)
    hfm(3, 2) = 2*vectt(1, 2)*vectt(2, 2)
    hfm(3, 3) = 2*vectt(1, 3)*vectt(2, 3)
    hfm(3, 4) = vectt(2, 1)*vectt(1, 2)+vectt(1, 1)*vectt(2, 2)
    hfm(3, 5) = vectt(2, 1)*vectt(1, 3)+vectt(1, 1)*vectt(2, 3)
    hfm(3, 6) = vectt(2, 2)*vectt(1, 3)+vectt(1, 2)*vectt(2, 3)
!
!     CONSTRUCTION DE LA MATRICE HFM * S : (3,9) AUX X PTS DE GAUSS
!
    hsfm(1, 1) = hfm(1, 1)
    hsfm(1, 2) = hfm(1, 4)
    hsfm(1, 3) = hfm(1, 5)
    hsfm(1, 4) = hfm(1, 4)
    hsfm(1, 5) = hfm(1, 2)
    hsfm(1, 6) = hfm(1, 6)
    hsfm(1, 7) = hfm(1, 5)
    hsfm(1, 8) = hfm(1, 6)
    hsfm(1, 9) = hfm(1, 3)
!
    hsfm(2, 1) = hfm(2, 1)
    hsfm(2, 2) = hfm(2, 4)
    hsfm(2, 3) = hfm(2, 5)
    hsfm(2, 4) = hfm(2, 4)
    hsfm(2, 5) = hfm(2, 2)
    hsfm(2, 6) = hfm(2, 6)
    hsfm(2, 7) = hfm(2, 5)
    hsfm(2, 8) = hfm(2, 6)
    hsfm(2, 9) = hfm(2, 3)
!
    hsfm(3, 1) = hfm(3, 1)
    hsfm(3, 2) = hfm(3, 4)
    hsfm(3, 3) = hfm(3, 5)
    hsfm(3, 4) = hfm(3, 4)
    hsfm(3, 5) = hfm(3, 2)
    hsfm(3, 6) = hfm(3, 6)
    hsfm(3, 7) = hfm(3, 5)
    hsfm(3, 8) = hfm(3, 6)
    hsfm(3, 9) = hfm(3, 3)
!
    if (ind .eq. 1) then
!
!     CONSTRUCTION DE HS  :  (2,6) AUX PTS D'INTEGRATION REDUITS
!
        hs(1, 1) = 2*vectt(1, 1)*vectt(3, 1)
        hs(1, 2) = 2*vectt(1, 2)*vectt(3, 2)
        hs(1, 3) = 2*vectt(1, 3)*vectt(3, 3)
        hs(1, 4) = vectt(3, 2)*vectt(1, 1)+vectt(3, 1)*vectt(1, 2)
        hs(1, 5) = vectt(1, 1)*vectt(3, 3)+vectt(3, 1)*vectt(1, 3)
        hs(1, 6) = vectt(3, 3)*vectt(1, 2)+vectt(1, 3)*vectt(3, 2)
!
        hs(2, 1) = 2*vectt(2, 1)*vectt(3, 1)
        hs(2, 2) = 2*vectt(2, 2)*vectt(3, 2)
        hs(2, 3) = 2*vectt(3, 3)*vectt(2, 3)
        hs(2, 4) = vectt(2, 1)*vectt(3, 2)+vectt(3, 1)*vectt(2, 2)
        hs(2, 5) = vectt(2, 1)*vectt(3, 3)+vectt(3, 1)*vectt(2, 3)
        hs(2, 6) = vectt(2, 2)*vectt(3, 3)+vectt(2, 3)*vectt(3, 2)
!
!     CONSTRUCTION DE LA MATRICE HS * S : (2,9) AUX PTS
!     D'INTEGRATION REDUITS
        hss(1, 1) = hs(1, 1)
        hss(1, 2) = hs(1, 4)
        hss(1, 3) = hs(1, 5)
        hss(1, 4) = hs(1, 4)
        hss(1, 5) = hs(1, 2)
        hss(1, 6) = hs(1, 6)
        hss(1, 7) = hs(1, 5)
        hss(1, 8) = hs(1, 6)
        hss(1, 9) = hs(1, 3)
!
        hss(2, 1) = hs(2, 1)
        hss(2, 2) = hs(2, 4)
        hss(2, 3) = hs(2, 5)
        hss(2, 4) = hs(2, 4)
        hss(2, 5) = hs(2, 2)
        hss(2, 6) = hs(2, 6)
        hss(2, 7) = hs(2, 5)
        hss(2, 8) = hs(2, 6)
        hss(2, 9) = hs(2, 3)
!
    end if
end subroutine
