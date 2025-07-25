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

subroutine mudirx(nbsom, geom, idim, al1, al2, &
                  axe, ang)
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
    integer(kind=8) :: idim, nbsom
    real(kind=8) :: geom(idim, nbsom), axe(3, 3), ang(2), al1, al2
!
!.......................................................................
! .                                                                    .
! .  - FONCTION REALISEE : CALCULE LES COSINUS DIRECTEURS DE LA MATRICE.
! .                        DE PASSAGE DU REPERE DE L'ELEMENT AU REPERE .
! .                        DE REFERENCE AINSI QUE LES 3 DIRECTIONS     .
! .                        NORMEES DU REPERE DE L'ELEMENT.             .
! .  - ARGUMENTS :                                                     .
! .                                                                    .
! .      ENTREE :  NBSOM  --> NB DE SOMMETS DE L'ELEMENT (3 OU 4)      .
! .                GEOM   --> TABLEAU DES COORDONNEES DES SOMMETS      .
! .                IDIM   --> DIMENSION DE GEOM (2 OU 3)               .
! .                AL1    --> ANGLE 1 DU REPERE DE REFERENCE           .
! .                AL2    --> ANGLE 2 DU REPERE DE REFERENCE           .
! .      SORTIE :                                                      .
! .                AXE    <-- 3 DIRECTIONS NORMEES DU REPERE DE        .
! .                           L'ELEMENT                                .
! .                ANG    <-- COSINUS ET SINUS DE L'ANGLE (REF,ELEM)   .
! .                                                                    .
! .  - ROUTINES APPELEES:                                              .
! .                                                                    .
! .    JEFINI                                                          .
!.......................................................................
!
    real(kind=8) :: xi1, xi2, xi3, x12, y12, z12, x13, y13, z13, x24, y24, z24
    real(kind=8) :: s
    real(kind=8) :: pjxi1, pjxi2, pjxi3, s1, s2, s3, psxin, coepi
!
!-----------------------------------------------------------------------
    coepi = r8dgrd()
    if (nbsom .ne. 3 .and. nbsom .ne. 4) then
        ASSERT(.false.)
    end if
    if (idim .ne. 2 .and. idim .ne. 3) then
        ASSERT(.false.)
    end if
    xi1 = cos(coepi*al2)*cos(coepi*al1)
    xi2 = cos(coepi*al2)*sin(coepi*al1)
    xi3 = sin(coepi*al2)
    s = (xi1**2+xi2**2+xi3**2)**0.5d0
    xi1 = xi1/s
    xi2 = xi2/s
    xi3 = xi3/s
    x12 = geom(1, 2)-geom(1, 1)
    y12 = geom(2, 2)-geom(2, 1)
    z12 = 0.d0
    if (idim .eq. 3) z12 = geom(3, 2)-geom(3, 1)
    s = (x12**2+y12**2+z12**2)**0.5d0
    axe(1, 1) = x12/s
    axe(2, 1) = y12/s
    axe(3, 1) = z12/s
    x13 = geom(1, 3)-geom(1, 1)
    y13 = geom(2, 3)-geom(2, 1)
    z13 = 0.d0
    if (idim .eq. 3) z13 = geom(3, 3)-geom(3, 1)
    if (nbsom .eq. 3) then
        s1 = y12*z13-z12*y13
        s2 = z12*x13-x12*z13
        s3 = x12*y13-y12*x13
        s = (s1**2+s2**2+s3**2)**0.5d0
        axe(1, 3) = s1/s
        axe(2, 3) = s2/s
        axe(3, 3) = s3/s
    end if
    if (nbsom .eq. 4) then
        x24 = geom(1, 4)-geom(1, 2)
        y24 = geom(2, 4)-geom(2, 2)
        z24 = 0.d0
        if (idim .eq. 3) z24 = geom(3, 4)-geom(3, 2)
        s1 = y13*z24-z13*y24
        s2 = z13*x24-x13*z24
        s3 = x13*y24-y13*x24
        s = (s1**2+s2**2+s3**2)**0.5d0
        axe(1, 3) = s1/s
        axe(2, 3) = s2/s
        axe(3, 3) = s3/s
    end if
    axe(1, 2) = axe(2, 3)*axe(3, 1)-axe(3, 3)*axe(2, 1)
    axe(2, 2) = axe(3, 3)*axe(1, 1)-axe(1, 3)*axe(3, 1)
    axe(3, 2) = axe(1, 3)*axe(2, 1)-axe(2, 3)*axe(1, 1)
    psxin = xi1*axe(1, 3)+xi2*axe(2, 3)+xi3*axe(3, 3)
    pjxi1 = xi1-psxin*axe(1, 3)
    pjxi2 = xi2-psxin*axe(2, 3)
    pjxi3 = xi3-psxin*axe(3, 3)
    s = (pjxi1**2+pjxi2**2+pjxi3**2)**0.5d0
    pjxi1 = pjxi1/s
    pjxi2 = pjxi2/s
    pjxi3 = pjxi3/s
    ang(1) = pjxi1*axe(1, 1)+pjxi2*axe(2, 1)+pjxi3*axe(3, 1)
    ang(2) = pjxi1*axe(1, 2)+pjxi2*axe(2, 2)+pjxi3*axe(3, 2)
end subroutine
