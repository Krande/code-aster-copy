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
subroutine ecla2d(nomte, elrefa, fapg, npg, npoini, &
                  nterm1, nsomm1, csomm1, tyma, nbno2, &
                  connx, mxnbn2, mxnbpi, mxnbte, mxnbse, &
                  nbsel, corsel)
    implicit none
#include "jeveux.h"
#include "asterfort/eclac1.h"
#include "asterfort/eclaco.h"
#include "asterfort/eclan1.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: mxnbn2, mxnbpi, mxnbte, mxnbse
    integer(kind=8) :: npg, connx(mxnbn2, mxnbse), nsomm1(mxnbpi, mxnbte)
    integer(kind=8) :: nterm1(mxnbpi), nbno2(mxnbse), npoini, tyma(mxnbse)
    integer(kind=8) :: nbsel, corsel(mxnbse)
    real(kind=8) :: csomm1(mxnbpi, mxnbte)
    character(len=16) :: nomte
    character(len=8) :: elrefa, fapg
! BUT : DECOMPOSER LES TYPE_ELEM 2D EN AUTANT DE SOUS-ELEMENTS QUE
!       DE POINTS DE GAUSS.
!
!     DECOUPAGE DU TR3, TR6, TR7 :
!       NPG= 1 :  3 NOEUDS, 1 TRIA3
!       NPG= 3 :  7 NOEUDS, 3 QUAD4
!       NPG= 6 : 10 NOEUDS, 3 TRIA3 ET 3 QUAD4
!     DECOUPAGE DU QU4, QU8, QU9 :
!       NPG= 1 :  4 NOEUDS, 1 QUAD4
!       NPG= 2 :  4 NOEUDS, 2 QUAD4
!       NPG= 4 :  9 NOEUDS, 4 QUAD4
!       NPG= 9 : 16 NOEUDS, 9 QUAD4
!
! ---------------------------------------------------------------------
! DESCRIPTION DES POINTS INTERMEDIAIRES (POINT_I) :
! ------------------------------------------------
! UN POINT_I EST DEFINI COMME UNE COMBINAISON LINEAIRE DES NOEUDS
! DE LA MAILLE SOUS-JACENTE AU TYPE_ELEM :
! POINT_I = SOMME COEF(K)*NOEUD(K)  (1<= K <=NTERMES)
!           NTERMES <= 27 (HEXA27)
!
! ---------------------------------------------------------------------
    integer(kind=8) :: k, itria3, iquad4
    character(len=24) :: valk(3)
! ---------------------------------------------------------------------
    call jemarq()
!
    call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA3'), itria3)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD4'), iquad4)
!
!     -----------------------------------------------------------------
!     ELEMENT TRIANGLE
!     -----------------------------------------------------------------
    if (elrefa .eq. 'TR3' .or. elrefa .eq. 'TR6' .or. elrefa .eq. 'TR7') then
!
        if (fapg .eq. 'FPG1') then
!           -----------------
            npoini = 3
            tyma(1) = itria3
            nbno2(1) = 3
!
!        -- DEFINITION DES POINT_I :
            nterm1(1) = 1
            call eclan1(1, mxnbpi, nsomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(1, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(2) = 1
            call eclan1(2, mxnbpi, nsomm1, nterm1, 2, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(2, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(3) = 1
            call eclan1(3, mxnbpi, nsomm1, nterm1, 3, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(3, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
!        -- CONNECTIVITE DES SOUS-ELEMENTS :
            call eclaco(1, mxnbn2, connx, nbno2, 1, &
                        2, 3, 0, 0, 0, &
                        0, 0)
!
!
        else if (fapg .eq. 'FPG3') then
!               -----------------
            npoini = 7
            do k = 1, npg
                tyma(k) = iquad4
                nbno2(k) = 4
            end do
!
!
!        -- DEFINITION DES POINT_I :
            nterm1(1) = 1
            call eclan1(1, mxnbpi, nsomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(1, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(2) = 1
            call eclan1(2, mxnbpi, nsomm1, nterm1, 2, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(2, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(3) = 1
            call eclan1(3, mxnbpi, nsomm1, nterm1, 3, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(3, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(4) = 3
            call eclan1(4, mxnbpi, nsomm1, nterm1, 1, &
                        2, 3, 0, 0, 0, &
                        0, 0)
            call eclac1(4, mxnbpi, csomm1, nterm1, 1, &
                        1, 1, 0, 0, 0, &
                        0, 0)
!
            nterm1(5) = 2
            call eclan1(5, mxnbpi, nsomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(5, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(6) = 2
            call eclan1(6, mxnbpi, nsomm1, nterm1, 2, &
                        3, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(6, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(7) = 2
            call eclan1(7, mxnbpi, nsomm1, nterm1, 3, &
                        1, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(7, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
!        -- CONNECTIVITE DES SOUS-ELEMENTS :
            call eclaco(1, mxnbn2, connx, nbno2, 1, &
                        5, 4, 7, 0, 0, &
                        0, 0)
            call eclaco(2, mxnbn2, connx, nbno2, 2, &
                        6, 4, 5, 0, 0, &
                        0, 0)
            call eclaco(3, mxnbn2, connx, nbno2, 6, &
                        3, 7, 4, 0, 0, &
                        0, 0)
!
        else if (fapg .eq. 'FPG6') then
!               -----------------
            npoini = 10
            tyma(1) = itria3
            tyma(2) = itria3
            tyma(3) = itria3
            tyma(4) = iquad4
            tyma(5) = iquad4
            tyma(6) = iquad4
            nbno2(1) = 3
            nbno2(2) = 3
            nbno2(3) = 3
            nbno2(4) = 4
            nbno2(5) = 4
            nbno2(6) = 4
!
!        -- DEFINITION DES POINT_I :
            nterm1(1) = 1
            call eclan1(1, mxnbpi, nsomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(1, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(2) = 1
            call eclan1(2, mxnbpi, nsomm1, nterm1, 2, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(2, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(3) = 1
            call eclan1(3, mxnbpi, nsomm1, nterm1, 3, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(3, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(4) = 2
            call eclan1(4, mxnbpi, nsomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(4, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(5) = 2
            call eclan1(5, mxnbpi, nsomm1, nterm1, 2, &
                        3, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(5, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(6) = 2
            call eclan1(6, mxnbpi, nsomm1, nterm1, 3, &
                        1, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(6, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(7) = 3
            call eclan1(7, mxnbpi, nsomm1, nterm1, 1, &
                        2, 3, 0, 0, 0, &
                        0, 0)
            call eclac1(7, mxnbpi, csomm1, nterm1, 1, &
                        2, 1, 0, 0, 0, &
                        0, 0)
!
            nterm1(8) = 3
            call eclan1(8, mxnbpi, nsomm1, nterm1, 1, &
                        2, 3, 0, 0, 0, &
                        0, 0)
            call eclac1(8, mxnbpi, csomm1, nterm1, 1, &
                        1, 2, 0, 0, 0, &
                        0, 0)
!
            nterm1(9) = 3
            call eclan1(9, mxnbpi, nsomm1, nterm1, 1, &
                        2, 3, 0, 0, 0, &
                        0, 0)
            call eclac1(9, mxnbpi, csomm1, nterm1, 2, &
                        1, 1, 0, 0, 0, &
                        0, 0)
!
            nterm1(10) = 3
            call eclan1(10, mxnbpi, nsomm1, nterm1, 1, &
                        2, 3, 0, 0, 0, &
                        0, 0)
            call eclac1(10, mxnbpi, csomm1, nterm1, 1, &
                        1, 1, 0, 0, 0, &
                        0, 0)
!
!        -- CONNECTIVITE DES SOUS-ELEMENTS :
            call eclaco(1, mxnbn2, connx, nbno2, 1, &
                        4, 6, 0, 0, 0, &
                        0, 0)
            call eclaco(2, mxnbn2, connx, nbno2, 2, &
                        5, 4, 0, 0, 0, &
                        0, 0)
            call eclaco(3, mxnbn2, connx, nbno2, 5, &
                        3, 6, 0, 0, 0, &
                        0, 0)
!
            call eclaco(4, mxnbn2, connx, nbno2, 4, &
                        7, 10, 9, 0, 0, &
                        0, 0)
            call eclaco(5, mxnbn2, connx, nbno2, 7, &
                        5, 8, 10, 0, 0, &
                        0, 0)
            call eclaco(6, mxnbn2, connx, nbno2, 6, &
                        9, 10, 8, 0, 0, &
                        0, 0)
!
        else
            valk(1) = nomte
            valk(2) = elrefa
            valk(3) = fapg
            call utmess('F', 'CALCULEL5_76', nk=3, valk=valk)
!
        end if
!
!
!     -----------------------------------------------------------------
!     ELEMENT QUADRANGLE
!     -----------------------------------------------------------------
    elseif (elrefa .eq. 'QU4' .or. elrefa .eq. 'QU8' .or. elrefa .eq. 'QU9') &
        then
!
        if (fapg .eq. 'FPG1') then
!           -----------------
            npoini = 4
            tyma(1) = iquad4
            nbno2(1) = 4
!
!        -- DEFINITION DES POINT_I :
            nterm1(1) = 1
            call eclan1(1, mxnbpi, nsomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(1, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(2) = 1
            call eclan1(2, mxnbpi, nsomm1, nterm1, 2, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(2, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(3) = 1
            call eclan1(3, mxnbpi, nsomm1, nterm1, 3, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(3, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(4) = 1
            call eclan1(4, mxnbpi, nsomm1, nterm1, 4, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(4, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
!        -- CONNECTIVITE DES SOUS-ELEMENTS :
            call eclaco(1, mxnbn2, connx, nbno2, 1, &
                        2, 3, 4, 0, 0, &
                        0, 0)
!
!
        else if (fapg .eq. 'FIS2') then
!              -----------------------
            npoini = 6
            tyma(1) = iquad4
            nbno2(1) = 4
            tyma(2) = iquad4
            nbno2(2) = 4
!
!        -- DEFINITION DES POINT_I :
            nterm1(1) = 1
            call eclan1(1, mxnbpi, nsomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(1, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(2) = 1
            call eclan1(2, mxnbpi, nsomm1, nterm1, 2, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(2, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(3) = 1
            call eclan1(3, mxnbpi, nsomm1, nterm1, 3, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(3, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(4) = 1
            call eclan1(4, mxnbpi, nsomm1, nterm1, 4, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(4, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(5) = 2
            call eclan1(5, mxnbpi, nsomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(5, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(6) = 2
            call eclan1(6, mxnbpi, nsomm1, nterm1, 3, &
                        4, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(6, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
!
!        -- CONNECTIVITE DES SOUS-ELEMENTS :
            call eclaco(1, mxnbn2, connx, nbno2, 1, &
                        5, 6, 4, 0, 0, &
                        0, 0)
            call eclaco(2, mxnbn2, connx, nbno2, 5, &
                        2, 3, 6, 0, 0, &
                        0, 0)
!
!
        else if (fapg .eq. 'FPG4') then
!               -----------------
            npoini = 9
            do k = 1, npg
                tyma(k) = iquad4
                nbno2(k) = 4
            end do
!
!        -- DEFINITION DES POINT_I :
            nterm1(1) = 1
            call eclan1(1, mxnbpi, nsomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(1, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(2) = 1
            call eclan1(2, mxnbpi, nsomm1, nterm1, 2, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(2, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(3) = 1
            call eclan1(3, mxnbpi, nsomm1, nterm1, 3, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(3, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(4) = 1
            call eclan1(4, mxnbpi, nsomm1, nterm1, 4, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(4, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(5) = 2
            call eclan1(5, mxnbpi, nsomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(5, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(6) = 2
            call eclan1(6, mxnbpi, nsomm1, nterm1, 2, &
                        3, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(6, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(7) = 2
            call eclan1(7, mxnbpi, nsomm1, nterm1, 3, &
                        4, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(7, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(8) = 2
            call eclan1(8, mxnbpi, nsomm1, nterm1, 4, &
                        1, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(8, mxnbpi, csomm1, nterm1, 1, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(9) = 4
            call eclan1(9, mxnbpi, nsomm1, nterm1, 1, &
                        2, 3, 4, 0, 0, &
                        0, 0)
            call eclac1(9, mxnbpi, csomm1, nterm1, 1, &
                        1, 1, 1, 0, 0, &
                        0, 0)
!
!        -- CONNECTIVITE DES SOUS-ELEMENTS :
            call eclaco(1, mxnbn2, connx, nbno2, 8, &
                        1, 5, 9, 0, 0, &
                        0, 0)
            call eclaco(2, mxnbn2, connx, nbno2, 5, &
                        2, 6, 9, 0, 0, &
                        0, 0)
            call eclaco(3, mxnbn2, connx, nbno2, 6, &
                        3, 7, 9, 0, 0, &
                        0, 0)
            call eclaco(4, mxnbn2, connx, nbno2, 7, &
                        4, 8, 9, 0, 0, &
                        0, 0)
!
!
        else if (fapg .eq. 'FPG9') then
!               -----------------
            npoini = 16
            do k = 1, npg
                tyma(k) = iquad4
                nbno2(k) = 4
            end do
!
!        -- DEFINITION DES POINT_I :
            nterm1(1) = 1
            call eclan1(1, mxnbpi, nsomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(1, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(2) = 1
            call eclan1(2, mxnbpi, nsomm1, nterm1, 2, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(2, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(3) = 1
            call eclan1(3, mxnbpi, nsomm1, nterm1, 3, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(3, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(4) = 1
            call eclan1(4, mxnbpi, nsomm1, nterm1, 4, &
                        0, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(4, mxnbpi, csomm1, nterm1, 1, &
                        0, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(5) = 2
            call eclan1(5, mxnbpi, nsomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(5, mxnbpi, csomm1, nterm1, 2, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(6) = 2
            call eclan1(6, mxnbpi, nsomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(6, mxnbpi, csomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(7) = 2
            call eclan1(7, mxnbpi, nsomm1, nterm1, 2, &
                        3, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(7, mxnbpi, csomm1, nterm1, 2, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(8) = 2
            call eclan1(8, mxnbpi, nsomm1, nterm1, 2, &
                        3, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(8, mxnbpi, csomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(9) = 2
            call eclan1(9, mxnbpi, nsomm1, nterm1, 3, &
                        4, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(9, mxnbpi, csomm1, nterm1, 2, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
!
            nterm1(10) = 2
            call eclan1(10, mxnbpi, nsomm1, nterm1, 3, &
                        4, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(10, mxnbpi, csomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
!
!
            nterm1(11) = 2
            call eclan1(11, mxnbpi, nsomm1, nterm1, 4, &
                        1, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(11, mxnbpi, csomm1, nterm1, 2, &
                        1, 0, 0, 0, 0, &
                        0, 0)
!
!
            nterm1(12) = 2
            call eclan1(12, mxnbpi, nsomm1, nterm1, 4, &
                        1, 0, 0, 0, 0, &
                        0, 0)
            call eclac1(12, mxnbpi, csomm1, nterm1, 1, &
                        2, 0, 0, 0, 0, &
                        0, 0)
!
            nterm1(13) = 3
            call eclan1(13, mxnbpi, nsomm1, nterm1, 1, &
                        2, 4, 0, 0, 0, &
                        0, 0)
            call eclac1(13, mxnbpi, csomm1, nterm1, 1, &
                        1, 1, 0, 0, 0, &
                        0, 0)
!
            nterm1(14) = 3
            call eclan1(14, mxnbpi, nsomm1, nterm1, 1, &
                        2, 3, 0, 0, 0, &
                        0, 0)
            call eclac1(14, mxnbpi, csomm1, nterm1, 1, &
                        1, 1, 0, 0, 0, &
                        0, 0)
!
            nterm1(15) = 3
            call eclan1(15, mxnbpi, nsomm1, nterm1, 2, &
                        3, 4, 0, 0, 0, &
                        0, 0)
            call eclac1(15, mxnbpi, csomm1, nterm1, 1, &
                        1, 1, 0, 0, 0, &
                        0, 0)
!
            nterm1(16) = 3
            call eclan1(16, mxnbpi, nsomm1, nterm1, 1, &
                        3, 4, 0, 0, 0, &
                        0, 0)
            call eclac1(16, mxnbpi, csomm1, nterm1, 1, &
                        1, 1, 0, 0, 0, &
                        0, 0)
!
!        -- CONNECTIVITE DES SOUS-ELEMENTS :
            call eclaco(1, mxnbn2, connx, nbno2, 1, &
                        5, 13, 12, 0, 0, &
                        0, 0)
            call eclaco(2, mxnbn2, connx, nbno2, 2, &
                        7, 14, 6, 0, 0, &
                        0, 0)
            call eclaco(3, mxnbn2, connx, nbno2, 8, &
                        3, 9, 15, 0, 0, &
                        0, 0)
            call eclaco(4, mxnbn2, connx, nbno2, 16, &
                        10, 4, 11, 0, 0, &
                        0, 0)
            call eclaco(5, mxnbn2, connx, nbno2, 6, &
                        14, 13, 5, 0, 0, &
                        0, 0)
            call eclaco(6, mxnbn2, connx, nbno2, 7, &
                        8, 15, 14, 0, 0, &
                        0, 0)
            call eclaco(7, mxnbn2, connx, nbno2, 15, &
                        9, 10, 16, 0, 0, &
                        0, 0)
            call eclaco(8, mxnbn2, connx, nbno2, 16, &
                        11, 12, 13, 0, 0, &
                        0, 0)
            call eclaco(9, mxnbn2, connx, nbno2, 13, &
                        14, 15, 16, 0, 0, &
                        0, 0)
!
        else
            valk(1) = nomte
            valk(2) = elrefa
            valk(3) = fapg
            call utmess('F', 'CALCULEL5_76', nk=3, valk=valk)
        end if
!
    else
        valk(1) = nomte
        valk(2) = elrefa
        call utmess('F', 'CALCULEL5_78', nk=2, valk=valk)
    end if
!
!     -- POUR TOUS LES SCHEMAS 2D, IL Y A IDENTITE : KSE -> KPG :
    nbsel = npg
    do k = 1, npg
        corsel(k) = k
    end do
!
    call jedema()
!
end subroutine
