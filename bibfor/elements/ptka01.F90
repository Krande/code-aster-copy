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
subroutine ptka01(sk, e, a, xl, xiy, &
                  xiz, xjx, g, alfay, alfaz, &
                  ey, ez, ist)
    implicit none
#include "asterc/r8gaem.h"
#include "asterfort/utmess.h"
    real(kind=8) :: sk(*)
    real(kind=8) :: e, a, xl, xiy, xiz, xjx, g, alfay, alfaz, ey, ez
    integer(kind=8) :: ist
!    -------------------------------------------------------------------
!    * CE SOUS PROGRAMME CALCULE LA MATRICE DE RAIDEUR DE L'ELEMENT DE
!    POUTRE DROITE A SECTION CONSTANTE.
!
!    * DESCRIPTION DE L'ELEMENT:
!      C'EST UN ELEMENT A DEUX NOEUDS ET A SIX DEGRES DE LIBERTES PAR
!      NOEUDS (3 DEPLACEMENTS ET 3 ROTATIONS).
!      IL PEUT PRENDRE EN COMPTE L'EFFORT TRANCHANT ET L'EXCENTRICITE
!      DU CENTRE DE ROTATION (CENTRE DE TORSION) PAR RAPPORT A LA FIBRE
!      NEUTRE (LIEU DES CENTRES DE GRAVITES).
!
!    * REMARQUE :
!      LA MATRICE EST STOCKEE TRIANGULAIRE INFERIEURE DANS UN TABLEAU
!      UNICOLONNE
!    -------------------------------------------------------------------
!  DONNEES NON MODIFIEES
!
! IN TYPE ! NOM    ! TABLEAU !             SIGNIFICATION
! IN -------------------------------------------------------------------
! IN R*8  ! E      !     -   ! MODULE D'ELASTICITE DU MATERIAU
! IN R*8  ! A      !     -   ! AIRE DE LA SECTION DROITE DE L'ELEMENT
! IN R*8  ! XL     !     -   ! LONGUEUR DE L ELEMENT
! IN R*8  ! XIY    !     -   ! MOMENT D INERTIE / Y PRINCIPAL
! IN R*8  ! XIZ    !     -   ! MOMENT D INERTIE / Z PRINCIPAL
! IN R*8  ! XJX    !     -   ! CONSTANTE DE TORSION
! IN R*8  ! G      !     -   ! MODULE DE CISAILLEMENT DU MATERIAU
! IN R*8  ! ALFAY  !     -   ! COEFFICIENT DE CISAILLEMENT AXE Y (+)
! IN R*8  ! ALFAZ  !     -   ! COEFFICIENT DE CISAILLEMENT AXE Z (+)
! IN R*8  ! EY     !     -   ! COMPOSANTE GT SUR Y PRINCIPAL
! IN R*8  ! EZ     !     -   ! COMPOSANTE GT SUR Z PRINCIPAL
! IN  I   ! IST    !    -    ! TYPE DE STRUCTURE DE LA POUTRE
! IN
! IN (+) REMARQUES :
! IN  -  LE COEFFICIENT DE CISAILLEMENT EST L'INVERSE DU COEFFICIENT DE
! IN     FORME ( IL EST DONC SUPERIEUR A 1)
! IN  -  SI ALFAY OU ALFAZ EST NUL ALORS ON CONSIDERE L'ELEMENT DE TYPE
! IN     EULER-BERNOULLI (I.E.  SANS EFFORT TRANCHANT)
!
! OUT TYPE ! NOM   ! TABLEAU !             SIGNIFICATION
! OUT ------------------------------------------------------------------
! OUT R*8 !   SK   ! (78)    ! MATRICE ELEMENTAIRE UNICOLONNE
!
!
! LOC TYPE !  NOM  ! TABLEAU !              SIGNIFICATION
! LOC ------------------------------------------------------------------
! LOC I   ! IP     !   12    ! POINTEUR SUR L'ELEMENT DIAGONAL PRECEDENT
!     ------------------------------------------------------------------
    integer(kind=8) :: ip(12), i
    real(kind=8) :: zero, xl2, xl3, phiy, phiz, eiy, eiz
!--- -------------------------------------------------------------------
    parameter(zero=0.d0)
    data ip/0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66/
!--- -------------------------------------------------------------------
    do i = 1, 78
        sk(i) = zero
    end do
!
! --- SI G ET E SONT NULS : K=0
    if (abs(g) .lt. 1.d0/r8gaem()) then
        if (abs(e) .lt. 1.d0/r8gaem()) goto 999
        call utmess('F', 'ELEMENTS2_54')
    end if
!
!     1/ TRACTION - COMPRESSION
    sk(1) = e*a/xl
    sk(ip(7)+1) = -sk(1)
    sk(ip(7)+7) = sk(1)
!
    if ((ist .eq. 2) .or. (ist .eq. 5)) goto 999
!
!     2/ FLEXION
!     2.1) CALCUL DES CONSTANTES
    xl2 = xl*xl
    xl3 = xl*xl2
    eiy = e*xiy
    eiz = e*xiz
    phiy = (12.d0*eiz*alfay)/(g*a*xl2)
    phiz = (12.d0*eiy*alfaz)/(g*a*xl2)
!
!     2.1) REMPLISSAGE DE LA MATRICE
!     FLEXION DANS LE PLAN XOY
    sk(ip(2)+2) = 12.d0*eiz/((1.d0+phiy)*xl3)
    sk(ip(6)+2) = 6.d0*eiz/((1.d0+phiy)*xl2)
    sk(ip(8)+2) = -sk(ip(2)+2)
    sk(ip(12)+2) = sk(ip(6)+2)
    sk(ip(6)+6) = (4.d0+phiy)*eiz/((1.d0+phiy)*xl)
    sk(ip(8)+6) = -sk(ip(6)+2)
    sk(ip(12)+6) = (2.d0-phiy)*eiz/((1.d0+phiy)*xl)
    sk(ip(8)+8) = sk(ip(2)+2)
    sk(ip(12)+8) = -sk(ip(6)+2)
    sk(ip(12)+12) = sk(ip(6)+6)
!
    if ((ist .eq. 3) .or. (ist .eq. 6)) goto 999
!
!     3/ FLEXION DANS LE PLAN XOZ
    sk(ip(3)+3) = 12.d0*eiy/((1.d0+phiz)*xl3)
    sk(ip(5)+3) = -6.d0*eiy/((1.d0+phiz)*xl2)
    sk(ip(9)+3) = -sk(ip(3)+3)
    sk(ip(11)+3) = sk(ip(5)+3)
    sk(ip(5)+5) = (4.d0+phiz)*eiy/((1.d0+phiz)*xl)
    sk(ip(9)+5) = -sk(ip(5)+3)
    sk(ip(11)+5) = (2.d0-phiz)*eiy/((1.d0+phiz)*xl)
    sk(ip(9)+9) = sk(ip(3)+3)
    sk(ip(11)+9) = -sk(ip(5)+3)
    sk(ip(11)+11) = sk(ip(5)+5)
!
!     4/ TORSION
    sk(ip(4)+4) = g*xjx/xl
    sk(ip(10)+4) = -sk(ip(4)+4)
    sk(ip(10)+10) = sk(ip(4)+4)
!
    if ((ez .eq. zero) .and. (ey .eq. zero)) goto 999
!
!     5/ AVEC EXCENTREMENT
!     RECTIFICATION POUR LA TORSION
    sk(ip(4)+4) = sk(ip(4)+4)+ez*ez*sk(ip(2)+2)+ey*ey*sk(ip(3)+3)
    sk(ip(10)+4) = -sk(ip(4)+4)
    sk(ip(10)+10) = sk(ip(4)+4)
!     TERME INDUIT PAR L'EXCENTREMENT
    sk(ip(4)+2) = -ez*sk(ip(2)+2)
    sk(ip(10)+2) = -sk(ip(4)+2)
    sk(ip(4)+3) = ey*sk(ip(3)+3)
    sk(ip(10)+3) = -sk(ip(4)+3)
    sk(ip(5)+4) = ey*sk(ip(5)+3)
    sk(ip(6)+4) = -ez*sk(ip(6)+2)
    sk(ip(8)+4) = sk(ip(10)+2)
    sk(ip(9)+4) = sk(ip(10)+3)
    sk(ip(11)+4) = sk(ip(5)+4)
    sk(ip(12)+4) = sk(ip(6)+4)
    sk(ip(10)+5) = -sk(ip(5)+4)
    sk(ip(10)+6) = -sk(ip(6)+4)
    sk(ip(10)+8) = sk(ip(4)+2)
    sk(ip(10)+9) = sk(ip(4)+3)
    sk(ip(11)+10) = sk(ip(10)+5)
    sk(ip(12)+10) = sk(ip(10)+6)
999 continue
end subroutine
