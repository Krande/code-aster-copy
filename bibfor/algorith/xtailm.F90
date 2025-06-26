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
subroutine xtailm(ndim, vecdir, numa, typma, jcoor, &
                  jconx1, jconx2, ipt, jtail)
    implicit none
!
#include "jeveux.h"
#include "asterfort/conare.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "blas/ddot.h"
    integer(kind=8) :: ndim, numa, jcoor, jconx1, jconx2, ipt, jtail
    real(kind=8) :: vecdir(ndim)
    character(len=8) :: typma
!
! person_in_charge: samuel.geniaut at edf.fr
!       ----------------------------------------------------------------
!       DETERMINATION DE LA TAILLE MAXIMALE DE LA MAILLE CONNECTEE AU
!       POINT DU FOND IPT SUIVANT LA DIRECTION DU VECTEUR VECDIR
!       ----------------------------------------------------------------
!    ENTREES
!       NDIM   : DIMENSION DU MODELE
!       VECDIR : VECTEUR TANGENT
!       NUMA   : NUMERO DE LA MAILLE COURANTE
!       TYPMA  : TYPE DE LA MAILLE COURANTE
!       JCOOR  : ADRESSE DU VECTEUR DES COORDONNEES DES NOEUDS
!       JCONX1 : ADRESSE DE LA CONNECTIVITE DU MAILLAGE
!       JCONX2 : LONGUEUR CUMULEE DE LA CONNECTIVITE DU MAILLAGE
!       IPT    : INDICE DU POINT DU FOND
!    SORTIE
!       JTAIL  : ADRESSE DU VECTEUR DES TAILLES MAXIMALES DES MAILLES
!                CONNECTEES AUX NOEUDS DU FOND DE FISSURE
!
!
    integer(kind=8) :: ar(12, 3)
    integer(kind=8) :: iar, ino1, ino2, k, nbar, nno1, nno2
    real(kind=8) :: arete(3), p
    blas_int :: b_incx, b_incy, b_n
!     -----------------------------------------------------------------
!
    call jemarq()
!
    call conare(typma, ar, nbar)
!
!     BOUCLE SUR LE NOMBRE D'ARETES DE LA MAILLE NUMA
    do iar = 1, nbar
!
        ino1 = ar(iar, 1)
        nno1 = zi(jconx1-1+zi(jconx2+numa-1)+ino1-1)
        ino2 = ar(iar, 2)
        nno2 = zi(jconx1-1+zi(jconx2+numa-1)+ino2-1)
!
!       VECTEUR REPRESENTANT L'ARETE IAR
        do k = 1, ndim
            arete(k) = zr(jcoor-1+(nno1-1)*3+k)-zr(jcoor-1+(nno2-1)*3+k)
        end do
!
!       PROJECTION DE L'ARETE IAR SUR LE VECTEUR TANGENT
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        p = ddot(b_n, arete, b_incx, vecdir, b_incy)
        p = abs(p)
!
        if (p .gt. zr(jtail-1+ipt)) zr(jtail-1+ipt) = p
!
    end do
!
    call jedema()
end subroutine
