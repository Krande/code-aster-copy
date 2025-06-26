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
subroutine porea1(nno, nc, deplm, deplp, geom, &
                  gamma, vecteu, pgl, xl1, angp)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8rddg.h"
#include "asterfort/angvx.h"
#include "asterfort/angvxy.h"
#include "asterfort/assert.h"
#include "asterfort/matrot.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
    integer(kind=8) :: nno, nc
    real(kind=8) :: deplm(nno*nc), deplp(nno*nc), geom(3, nno), gamma
!
    real(kind=8) :: pgl(3, 3), xl1, angp(3)
    aster_logical :: vecteu
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DE LA MATRICE DE PASSAGE GLOBALE/LOCALE EN TENANT COMPTE
!     DE LA GEOMETRIE REACTUALISEE POUR LES POUTRES AINSI QUE LA
!     LONGUEUR DE LA POUTRE
!     POUR LES OPTIONS FULL_MECA RAPH_MECA ET RIGI_MECA_TANG
!
! --------------------------------------------------------------------------------------------------
!
! IN  NNO    : NOMBRE DE NOEUDS
! IN  NC     : NOMBRE DE COMPOSANTE DU CHAMP DE DEPLACEMENTS
! IN  DEPLM  : DEPLACEMENT AU TEMPS -
! IN  DEPLP  : INCREMENT DE DEPLACEMENT AU TEMPS +
! IN  GEOM   : COORDONNEES DES NOEUDS
! IN  GAMMA  : ANGLE DE VRILLE AU TEMPS -
! IN  VECTEU : TRUE SI FULL_MECA OU RAPH_MECA
! OUT PGL    : MATRICE DE PASSAGE GLOBAL/LOCAL
! OUT XL     : LONGUEUR DE L'ELEMENT
! OUT ANGP   : ANGLES NAUTIQUES ACTUALISEE
!              ATTENTION ANGP(3) EST DIFFERENT DE ANG1(3) QUI A SERVIT
!              POUR CALCUL PGL
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i
    real(kind=8) :: utg(14), xug(6), xd0(3), xd1(3), xdn0(3), xdn1(3), cosangle
    real(kind=8) :: tet1, tet2, xl0, pgl2(3, 3)
    real(kind=8) :: ang1(3), angm(3), ytrm(3), ytrp(3), vect(3), sinangle, angle
!
    integer(kind=8) :: iadzi, iazk24
    character(len=24) :: valkm
    real(kind=8) :: valrm
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
    ASSERT(nno .eq. 2)
!   Calcul du vecteur xlocal au temps t-
    do i = 1, 3
        xug(i) = geom(i, 1)+deplm(i)
        xug(i+3) = geom(i, 2)+deplm(i+nc)
        xd0(i) = xug(i+3)-xug(i)
    end do
!   Déplacement total a t+
    do i = 1, nno*nc
        utg(i) = deplm(i)+deplp(i)
    end do
!   Calcul du vecteur xlocal au temps t+
    do i = 1, 3
        xug(i) = geom(i, 1)+utg(i)
        xug(i+3) = geom(i, 2)+utg(i+nc)
        xd1(i) = xug(i+3)-xug(i)
    end do
!   Angle entre xd0 et xd1
    xdn0(:) = xd0(:)
    xdn1(:) = xd1(:)
    call normev(xdn0, xl0)
    call normev(xdn1, xl1)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    cosangle = ddot(b_n, xdn0, b_incx, xdn1, b_incy)
    call provec(xdn0, xdn1, vect)
    sinangle = (ddot(b_n, vect, b_incx, vect, b_incy))**0.5
    angle = atan2(sinangle, cosangle)*r8rddg()
!   si angle > pi/4
    if (abs(angle) .gt. 22.5d0) then
        call tecael(iadzi, iazk24)
        valkm = zk24(iazk24+3-1)
        valrm = angle
        call utmess('A', 'ELEMENTS_38', sk=valkm, sr=valrm)
    end if
!
    if (vecteu) then
!       mise a jour du 3eme angle nautique au temps t+

!       ON CALCULE LES VECTEURS LOCAUX À L'INSTANT T-
        call angvx(xd0, angm(1), angm(2))
        angm(3) = gamma
        call matrot(angm, pgl)

!       VY À L'INSTANT T-
        ytrm(1) = pgl(2, 1)
        ytrm(2) = pgl(2, 2)
        ytrm(3) = pgl(2, 3)

!       CALCUL ET SAUVEGARDE DES ANGLES NAUTIQUES À L'INSTANT T+
!       ON S'EN SERT DE VY À L'INSTANT T- COMMME VECT_Y À L'INSTANT T+
!       ztrp = xdp ^ ytrm
!       ytrp = ztrp ^ xdp

        call angvxy(xd1, ytrm, angp)

    else
        call angvx(xd1, angp(1), angp(2))
        angp(3) = gamma
    end if
!
!   CALCUL DES TORSIONS AUX NOEUDS EXPRIMEE DABNS LE REPERE GLOBAL À T+
!
    tet1 = ddot(b_n, utg(4), b_incx, xd1, b_incy)
    tet2 = ddot(b_n, utg(nc+4), b_incx, xd1, b_incy)
    tet1 = tet1/xl1
    tet2 = tet2/xl1
!   Matrice de passage global -> local
    ang1(1) = angp(1)
    ang1(2) = angp(2)
    ang1(3) = angp(3)+0.5d0*(tet1+tet2)
    call matrot(ang1, pgl)

!   VERIF ANGLE ENTRE YLOCAL T- ET T+
    if (vecteu) then
!       0.5d0*(tet1+tet2) n'est pas pris en compte dans angm, il faut
!       donc comparer à angp et non à ang1
        call matrot(angp, pgl2)
        ytrp(1) = pgl2(2, 1)
        ytrp(2) = pgl2(2, 2)
        ytrp(3) = pgl2(2, 3)
        cosangle = ddot(b_n, ytrm, b_incx, ytrp, b_incy)
        call provec(ytrm, ytrp, vect)
        sinangle = (ddot(b_n, vect, b_incx, vect, b_incy))**0.5
        angle = atan2(sinangle, cosangle)*r8rddg()
        if (abs(angle) .gt. 22.5d0) then
            call tecael(iadzi, iazk24)
            valkm = zk24(iazk24+3-1)
            valrm = angle
            call utmess('A', 'ELEMENTS_38', sk=valkm, sr=valrm)
        end if
    end if

end subroutine
