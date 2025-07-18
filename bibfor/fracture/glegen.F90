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
subroutine glegen(nbre, lobj2, xl, absgam, legen)
    implicit none
!
!     ------------------------------------------------------------------
!
! FONCTION REALISEE:
!
!      CALCUL DES VALEURS DES POLYNOMES DE LEGENDRE
!
!     ------------------------------------------------------------------
! ENTREE:
!        NBRE       : DEGRE DES POLYNOMES DE LEGENDRE
!        LOBJ2      : NOMBRE DE NOEUDS SUR GAMMA0
!        XL         : LONGUEUR DE LA FISSURE
!        ABSGAM     : ABSCISSES CURVILIGNES
!
! SORTIE:
!        LEGEN      : VALEURS DES POLYNOMES DE LEGENDRE
!     ------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: lobj2, nbre, iadabs, iadpo
!
    real(kind=8) :: xl, s1, coef
!
    character(len=24) :: absgam
    real(kind=8) :: legen(1)
    integer(kind=8) :: i, j
    real(kind=8) :: cof1, cof2
!
!   POLYNOMES DE LEGENDRE
!
#define pleg2(x) (3.d0*(x)*(x)-1.d0)/2.d0
#define pleg3(x) (x)*(5.d0*(x)*(x)-3.d0)/2.d0
#define pleg4(x) (35.d0*(x)**4-30.d0*(x)*(x)+3.d0)/8.d0
#define pleg5(x) (x)*(63.d0*(x)**4-70.d0*(x)*(x)+15.d0)/8.d0
#define pleg6(x) (231.d0*(x)**6-315.d0*(x)**4+105.d0*(x)*(x)-5.d0)/16.d0
#define pleg7(x) (x)*(429.d0*(x)**6-693.d0*(x)**4+315.d0*(x)*(x)-35.d0)/16.d0
!
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(absgam, 'L', iadabs)
!
! CALCUL DES MODULES DES CHAMPS THETA = POLYNOMES DE LEGENDRE
!
    call wkvect('&&LEGEND.VALPOL', 'V V R8', (nbre+1)*lobj2, iadpo)
!
! COEFFICIENTS POUR LES 2 PREMIERS POLYNOMES
!
    cof1 = sqrt(1.d0/xl)
    cof2 = sqrt(3.d0/xl)
!
    do i = 1, lobj2
        zr(iadpo+i-1) = 1.d0
        legen(i) = cof1*zr(iadpo+i-1)
        if (nbre .ne. 0) then
            zr(iadpo+lobj2+i-1) = -1+2*zr(iadabs+i-1)/xl
            legen(lobj2+i) = cof2*zr(iadpo+lobj2+i-1)
        end if
    end do
!
    if (nbre .ge. 2) then
        coef = sqrt(5.d0/xl)
        do j = 1, lobj2
            s1 = -1+2*zr(iadabs+j-1)/xl
            legen(2*lobj2+j) = coef*pleg2(s1)
        end do
    end if
!
    if (nbre .ge. 3) then
        coef = sqrt(7.d0/xl)
        do j = 1, lobj2
            s1 = -1+2*zr(iadabs+j-1)/xl
            legen(3*lobj2+j) = coef*pleg3(s1)
        end do
    end if
!
    if (nbre .ge. 4) then
        coef = sqrt(9.d0/xl)
        do j = 1, lobj2
            s1 = -1+2*zr(iadabs+j-1)/xl
            legen(4*lobj2+j) = coef*pleg4(s1)
        end do
    end if
!
    if (nbre .ge. 5) then
        coef = sqrt(11.d0/xl)
        do j = 1, lobj2
            s1 = -1+2*zr(iadabs+j-1)/xl
            legen(5*lobj2+j) = coef*pleg5(s1)
        end do
    end if
!
    if (nbre .ge. 6) then
        coef = sqrt(13.d0/xl)
        do j = 1, lobj2
            s1 = -1+2*zr(iadabs+j-1)/xl
            legen(6*lobj2+j) = coef*pleg6(s1)
        end do
    end if
!
    if (nbre .ge. 7) then
        coef = sqrt(15.d0/xl)
        do j = 1, lobj2
            s1 = -1+2*zr(iadabs+j-1)/xl
            legen(7*lobj2+j) = coef*pleg7(s1)
        end do
    end if
!
    call jedetr('&&LEGEND.VALPOL')
    call jedema()
end subroutine
