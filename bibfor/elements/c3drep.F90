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
subroutine c3drep(nomte, epais, alpha, beta, coord, &
                  numnoe, pgl)
    implicit none
#include "jeveux.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vdrep2.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
    integer(kind=8) :: numnoe
    character(len=16) :: nomte
    real(kind=8) :: epais, alpha, beta, coord(3, 9), pgl(3, 3)
! person_in_charge: nicolas.sellenet at edf.fr
!     ------------------------------------------------------------------
!
!         CETTE ROUTINE REALISE LA MEME TACHE QUE COQREP MAIS POUR LES
!         COQUES 3D
!         CALCUL DE LA MATRICE DE PASSAGE DU REPERE DE L'ELEMENT A
!         LA VARIETE (LE REPERE DE LA VARIETE EST OBTENU PAR LA MATRICE
!         DE PASSAGE GLOBAL -> LOCAL) AINSI QUE SON INVERSE
!
!     ------------------------------------------------------------------
    integer(kind=8) :: nb1, nb2, npgsr, i, j, k, ind, intsr
!
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectg(2, 3), vectt(3, 3)
    real(kind=8) :: zero, vectpt(9, 2, 3), vectmp(3, 3), pgltmp(3, 3)
    real(kind=8) :: matevn(2, 2, 10), matevg(2, 2, 10), v
    real(kind=8), pointer :: desr(:) => null()
    integer(kind=8), pointer :: desi(:) => null()
!
    zero = 0.d0
    call jeveuo('&INEL.'//nomte(1:8)//'.DESI', 'L', vi=desi)
    call jeveuo('&INEL.'//nomte(1:8)//'.DESR', 'L', vr=desr)
    nb1 = desi(1)
    nb2 = desi(2)
    npgsr = desi(3)
!
!     -- POUR REMPLIR LZR+1090+...  ET CALCULER VECTN :
    call vectan(nb1, nb2, coord, desr, vecta, &
                vectn, vectpt)
!
!     -- POUR REMPLIR LZR+2000+... :
!     -- QUELLE VALEUR POUR IND ? FICHE ???
! ind=0 => calcul aux points d'intégration réduite
! ind=1 => calcul aux points d'intégration normale
    ind = 0
    k = 0
    do intsr = 1, npgsr
        call vectgt(ind, nb1, coord, zero, intsr, &
                    desr, epais, vectn, vectg, vectt)
        do j = 1, 3
            do i = 1, 3
                k = k+1
                desr(1+2000+k-1) = vectt(i, j)
            end do
        end do
    end do
!
    call vdrep2(alpha, beta, desi, desr, matevn, &
                matevg)
!
    vectmp(1, 1) = matevn(1, 1, numnoe)
    vectmp(1, 2) = matevn(1, 2, numnoe)
    vectmp(2, 1) = matevn(2, 1, numnoe)
    vectmp(2, 2) = matevn(2, 2, numnoe)
    vectmp(1, 3) = 0.d0
    vectmp(2, 3) = 0.d0
    vectmp(3, 3) = 1.d0
    vectmp(3, 1) = 0.d0
    vectmp(3, 2) = 0.d0
!
    k = 0
    do j = 1, 3
        do i = 1, 3
            k = k+1
            pgltmp(i, j) = desr(1+1090+(numnoe-1)*9+k-1)
        end do
    end do
    do i = 1, 3
        do j = 1, 3
            v = 0.d0
            do k = 1, 3
                v = v+vectmp(i, k)*pgltmp(k, j)
            end do
            pgl(i, j) = v
        end do
    end do
!
end subroutine
