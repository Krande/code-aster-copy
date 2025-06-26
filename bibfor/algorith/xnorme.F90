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
subroutine xnorme(indipt, iptbor, vectn, nbfacb, nunoa, &
                  nunob, nunoc, jcoor, coorg)
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/provec.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: iptbor(2), nbfacb, jcoor
    integer(kind=8) :: nunoa, nunob, nunoc, indipt
    real(kind=8) :: vectn(12), coorg(3)
!
!
!             CALCUL ET STOCKAGE DE LA NORMALE D'UNE FACE
!
!     ENTREE
!       INDIPT   : INDICE DU POINT DE BORD
!       VECTN    : VECTEUR CONTENANT LES VECTEURS NORMAUX DES FACES DE
!                  BORD DE LA MAILLE
!       NBFACB   : NOMBRE DE FACES DE BORD DANS LA MAILLE
!       NUNOA    : NOEUD DE LA FACE COURANTE
!       NUNOB    : NOEUD DE LA FACE COURANTE
!       NUNOC    : NOEUD DE LA FACE COURANTE
!       JCOOR    : ADRESSE DES COORDONNEES DES NOEUDS DU MAILLAGE
!       COORG    : COORDONNEES DU CENTRE DE GRAVITE DE LA MAILLE
!
!     SORTIE
!       IPTBOR   : VECTEUR CONTENANT LES INDICES DES POINTS DE BORD DE
!                  LA MAILLE
!     ------------------------------------------------------------------
!
    integer(kind=8) :: k
    real(kind=8) :: ab(3), ac(3), ag(3), normal(3), proj
    blas_int :: b_incx, b_incy, b_n
! ----------------------------------------------------------------------
    call jemarq()
!
!     CALCUL DE LA NORMALE
    do k = 1, 3
!       A,B ET C SONT DES NOEUDS DE LA FACE
        ab(k) = zr(jcoor-1+3*(nunob-1)+k)-zr(jcoor-1+3*(nunoa-1)+k)
        ac(k) = zr(jcoor-1+3*(nunoc-1)+k)-zr(jcoor-1+3*(nunoa-1)+k)
!       G EST LE CENTRE DE GRAVITE DE LA MAILLE
        ag(k) = coorg(k)-zr(jcoor-1+3*(nunoa-1)+k)
    end do
!
    call provec(ab, ac, normal)
!
!     ORIENTATION DE LA NORMALE VERS L'EXTERIEUR
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    proj = ddot(b_n, normal, b_incx, ag, b_incy)
!
    if (proj .gt. 0) then
        normal(1) = -normal(1)
        normal(2) = -normal(2)
        normal(3) = -normal(3)
    end if
!
!     NOMBRE DE FACES DE BORD DANS LA MAILLE
    nbfacb = nbfacb+1
!
!     STOCKAGE DE LA NORMALE
    vectn(1+3*(nbfacb-1)) = normal(1)
    vectn(2+3*(nbfacb-1)) = normal(2)
    vectn(3+3*(nbfacb-1)) = normal(3)
!
!     IPTBOR(2) SERT SI LA MAILLE POSSEDE 2 POINTS
!     DE FOND QUI SONT SUR UNE FACE DE BORD
!     (CAS TRES SPECIAL OU LE FOND DE FISSURE N'EST QUE DANS UNE SEULE
!      MAILLE)
    if (iptbor(1) .eq. 0) then
        iptbor(1) = indipt
    else
        iptbor(2) = indipt
    end if
!
    call jedema()
end subroutine
