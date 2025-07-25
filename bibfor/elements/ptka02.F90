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
subroutine ptka02(id, sk, e, a1, a2, &
                  xl, xiy, xiy2, xiz, xiz2, &
                  xjx, xjx2, g, alfay1, alfay2, &
                  alfaz1, alfaz2, ey, ez, ist)
    implicit none
#include "asterc/r8gaem.h"
#include "asterfort/fun1.h"
#include "asterfort/fun2.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: id, ist
    real(kind=8) :: sk(*), e, a1, a2, xl, xiy, xiy2
    real(kind=8) :: xiz, xiz2, xjx, xjx2, g, alfay1, alfay2
    real(kind=8) :: alfaz1, alfaz2, ey, ez
!    * CE SOUS PROGRAMME CALCULE LA MATRICE DE RAIDEUR DE L'ELEMENT DE
!      POUTRE DROITE A SECTION VARIABLE.
!
!    * DESCRIPTION DE L'ELEMENT:
!      C'EST UN ELEMENT A DEUX NOEUDS ET A SIX DEGRES DE LIBERTES PAR
!      NOEUDS (3 DEPLACEMENTS ET 3 ROTATIONS)
!      IL PEUT PRENDRE EN COMPTE L'EFFORT TRANCHANT ET L'EXCENTRICITE DU
!      CENTRE DE ROTATION (CENTRE DE TORSION) PAR RAPPORT A LA FIBRE
!      NEUTRE (LIEU GEOMETRIQUE DES CENTRES DE GRAVITES).
!
!    * REMARQUE FONDAMENTALE SUR LA VARIATION DE LA SECTION DROITE :
!      DANS LE CAS DE SECTION AFFINE ON SUPPOSE QUE LA COMPOSANTE EN Y
!      VARIE ALORS QUE CELLE SELON Z EST FIXE
!
!    * REMARQUE : LA MATRICE EST STOCKEE TRIANGULAIRE INFERIEURE DANS UN
!      TABLEAU UNICOLONNE
!     ------------------------------------------------------------------
!
!  DONNEES NON MODIFIEES
!
! IN TYPE ! NOM    ! TABLEAU !             SIGNIFICATION
! IN -------------------------------------------------------------------
! IN  I   ! ID     !   -     ! TYPE DE VARIATION SECTION DROITE
! IN      !        !         !   ID = 1 : SECTIONS AFFINES
! IN      !        !         !   ID = 2 : SECTIONS HOMOTHETIQUES
! IN R*8  ! E      !   -     ! MODULE D'ELASTICITE DU MATERIAU
! IN R*8  ! A1     !   -     ! AIRE DE LA SECTION DROITE INITIALE
! IN R*8  ! A2     !   -     ! AIRE DE LA SECTION DROITE FINALE
! IN R*8  ! XL     !   -     ! LONGUEUR DE L ELEMENT
! IN R*8  ! XIY    !   -     ! MOMENT D INERTIE / Y PRINCIPAL SECTION
! IN      !        !         !   INITIALE
! IN R*8  ! XIY2   !   -     ! MOMENT D INERTIE / Y PRINCIPAL SECTION
! IN      !        !         !   FINALE
! IN R*8  ! XIZ    !   -     ! MOMENT D INERTIE / Z PRINCIPAL SECTION
! IN      !        !         !   INITIALE
! IN R*8  ! XIZ2   !   -     ! MOMENT D INERTIE / Z PRINCIPAL SECTION
! IN      !        !         !   FINALE
! IN R*8  ! XJX    !   -     ! CONSTANTE DE TORSION SECTION INITIALE
! IN R*8  ! XJX2   !   -     ! CONSTANTE DE TORSION SECTION FINALE
! IN R*8  ! G      !   -     ! MODULE DE CISAILLEMENT DU MATERIAU
! IN R*8  ! ALFAY1 !   -     ! COEFFICIENT DE CISAILLEMENT AXE Y
! IN      !        !         !   SECTION INITIALE (+)
! IN R*8  ! ALFAY2 !   -     ! COEFFICIENT DE CISAILLEMENT AXE Y
! IN      !        !         !   SECTION FINALE   (+)
! IN R*8  ! ALFAZ1 !   -     ! COEFFICIENT DE CISAILLEMENT AXE Z
! IN      !        !         !   SECTION INITIALE (+)
! IN R*8  ! ALFAZ2 !   -     ! COEFFICIENT DE CISAILLEMENT AXE Z
! IN      !        !         !   SECTION FINALE   (+)
! IN R*8  ! EY     !   -     ! COMPOSANTE TG SUR Y PRINCIPAL
! IN R*8  ! EZ     !   -     ! COMPOSANTE TG SUR Z PRINCIPAL
! IN  I   ! IST    !   -     ! TYPE DE STRUCTURE DE LA POUTRE
! IN
! IN (+) REMARQUES :
! IN  -  LE COEFFICIENT DE CISAILLEMENT EST L'INVERSE DU COEFFICIENT DE
! IN     FORME ( IL EST DONC SUPERIEUR A 1)
! IN  -  SI ALFAY OU ALFAZ EST NUL ALORS ON CONSIDERE L'ELEMENT DE TYPE
! IN     EULER-BERNOULLI (I.E.  SANS EFFORT TRANCHANT)
!
! OUT R*8 ! SK     ! (78)    ! MATRICE ELEMENTAIRE UNICOLONNE
!
! LOC R*8 ! ASY    !   -     ! AIRE REDUITE CISAILLEE SUIVANT Y
! LOC R*8 ! ASZ    !   -     ! AIRE REDUITE CISAILLEE SUIVANT Z
! LOC R*8 ! PHIY   !   -     ! COEFFICIENT DU A L'EFFORT TRANCHANT SUR Z
! LOC R*8 ! PHIZ   !   -     ! COEFFICIENT DU A L'EFFORT TRANCHANT SUR Y
! LOC I   ! IP     !   -     ! POINTEUR SUR L'ELEMENT DIAGONAL PRECEDENT
!     ------------------------------------------------------------------
!
!     SOUS - PROGRAMMES UTILISES
! BBL  FUN1     - AIRES ET CONSTANTE DE TORSION EQUIVALENTES
! BBL  FUN2     - MOMENTS D INERTIE EQUIVALENTS
!     ------------------------------------------------------------------
    integer(kind=8) :: ip(12), i, k
    real(kind=8) :: exl, xl2, xl3, phiy, phiz, asy, asz
    real(kind=8) :: aa, aas1, aas2, as1, as2, tk, vt, q, xkk, zero
!-----------------------------------------------------------------------
    parameter(zero=0.d0)
    data ip/0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66/
! ---------------------------------------------------------------------
    do i = 1, 78
        sk(i) = zero
    end do
!
! --- SI G  ET E SONT NULS : K=0
    if (abs(g) .lt. 1.d0/r8gaem()) then
        if (abs(e) .lt. 1.d0/r8gaem()) goto 999
        call utmess('F', 'ELEMENTS2_54')
    end if
!
!     1/ TRACTION-COMPRESSION
    call fun1(aa, a1, a2, id)
!
    sk(ip(1)+1) = e*aa/xl
    sk(ip(7)+1) = -sk(ip(1)+1)
    sk(ip(7)+7) = sk(ip(1)+1)
!
    if (ist .eq. 2 .or. ist .eq. 5) goto 999
!
!     2/ FLEXION
!     2.1) CALCUL DES CONSTANTES
    xl2 = xl*xl
    xl3 = xl*xl2
    exl = e/xl
    if (alfaz1 .ne. zero) then
        as1 = a1/alfaz1
    else
        as1 = a1
    end if
    if (alfaz2 .ne. zero) then
        as2 = a2/alfaz2
    else
        as2 = a2
    end if
    if (alfaz1 .eq. zero .and. alfaz2 .eq. zero) then
        phiz = zero
    else
        call fun1(asy, as1, as2, id)
        phiz = e/(g*asy*xl2)
    end if
!
    if (alfay1 .ne. zero) then
        aas1 = a1/alfay1
    else
        aas1 = a1
    end if
    if (alfay2 .ne. zero) then
        aas2 = a2/alfay2
    else
        aas2 = a2
    end if
    if (alfay1 .eq. zero .and. alfay2 .eq. zero) then
        phiy = zero
    else
        call fun1(asz, aas1, aas2, id)
        phiy = e/(g*asz*xl2)
    end if
!
!     2.2) REMPLISSAGE DE LA MATRICE
!
!     2.2.1/  FLEXION DANS LE PLAN X0Y
    k = id+2
    call fun2(xiz, xiz2, phiy, xkk, q, &
              vt, k)
    sk(ip(2)+2) = e*xkk/xl3
    sk(ip(6)+2) = xl*q*sk(ip(2)+2)
    sk(ip(8)+2) = -sk(ip(2)+2)
    sk(ip(12)+2) = (1.d0/q-1.d0)*sk(ip(6)+2)
    sk(ip(6)+6) = exl*(vt+q*q*xkk)
    sk(ip(8)+6) = -sk(ip(6)+2)
    sk(ip(12)+6) = exl*(xkk*q*(1.d0-q)-vt)
!
    call fun2(xiz2, xiz, phiy, xkk, q, &
              vt, k)
    sk(ip(8)+8) = sk(ip(2)+2)
    sk(ip(12)+8) = -sk(ip(12)+2)
    sk(ip(12)+12) = exl*(vt+q*q*xkk)
!
    if (ist .eq. 3 .or. ist .eq. 6) goto 999
!
!     2.2.2/  FLEXION DANS LE PLAN X0Z
    if (id .eq. 2) then
        k = 4
    else
        k = 1
    end if
!
    call fun2(xiy, xiy2, phiz, xkk, q, &
              vt, k)
    sk(ip(3)+3) = e*xkk/xl3
    sk(ip(5)+3) = -xl*q*sk(ip(3)+3)
    sk(ip(9)+3) = -sk(ip(3)+3)
    sk(ip(11)+3) = (1.d0/q-1.d0)*sk(ip(5)+3)
    sk(ip(5)+5) = exl*(vt+q*q*xkk)
    sk(ip(9)+5) = -sk(ip(5)+3)
    sk(ip(11)+5) = exl*(xkk*q*(1.d0-q)-vt)
!
    call fun2(xiy2, xiy, phiz, xkk, q, &
              vt, k)
    sk(ip(9)+9) = sk(ip(3)+3)
    sk(ip(11)+9) = -sk(ip(11)+3)
    sk(ip(11)+11) = exl*(vt+q*q*xkk)
!     2.2.3/  TORSION
    call fun1(tk, xjx, xjx2, id+2)
    sk(ip(4)+4) = g*tk/xl
    sk(ip(10)+4) = -sk(ip(4)+4)
    sk(ip(10)+10) = sk(ip(4)+4)
!
    if (ez .eq. zero .and. ey .eq. zero) goto 999
!
!     2.2.4/ AVEC EXCENTREMENT, TORSION DANS CE CAS
    sk(ip(4)+4) = sk(ip(4)+4)+ez*ez*sk(ip(2)+2)+ey*ey*sk(ip(3)+3)
    sk(ip(10)+4) = -sk(ip(4)+4)
    sk(ip(10)+10) = sk(ip(4)+4)
!
    sk(ip(4)+2) = -ez*sk(ip(2)+2)
    sk(ip(10)+2) = -sk(ip(4)+2)
    sk(ip(4)+3) = ey*sk(ip(3)+3)
    sk(ip(10)+3) = -sk(ip(4)+3)
    sk(ip(5)+4) = ey*sk(ip(5)+3)
    sk(ip(6)+4) = -ez*sk(ip(6)+2)
    sk(ip(8)+4) = sk(ip(10)+2)
    sk(ip(9)+4) = sk(ip(10)+3)
    sk(ip(11)+4) = ey*sk(ip(11)+3)
    sk(ip(12)+4) = -ez*sk(ip(12)+2)
    sk(ip(10)+5) = -sk(ip(5)+4)
    sk(ip(10)+6) = -sk(ip(6)+4)
    sk(ip(10)+8) = sk(ip(4)+2)
    sk(ip(10)+9) = sk(ip(4)+3)
    sk(ip(11)+10) = -sk(ip(11)+4)
    sk(ip(12)+10) = -sk(ip(12)+4)
!     SORTIE
999 continue
end subroutine
