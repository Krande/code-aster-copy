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
subroutine nmdcrg(depart, iterat, vresi, xa0, xa1, &
                  xdet)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: depart, iterat
    real(kind=8) :: vresi(*)
    real(kind=8) :: xa0, xa1, xdet
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (EVENEMENTS)
!
! CALCUL DE L'EXTRAPOLATION SUR LES RESIDUS
!
! ----------------------------------------------------------------------
!
!
! EXTRAPOLATION LINEAIRE (XA0 + ITER*XA1) / XDET
!
! ON DONNE UN POIDS DOUBLE AU 3 DERNIERS POINTS CELA REVIENT
! A AJOUTER DES POINTS => MEILLEURE EXTRAPOLATION
!
! IN  DEPART : DEBUT DE L'EXTRAPOLATION
! IN  ITERAT : FIN DE L'EXTRAPOLATION
! IN  VRESI  : LISTE DES VALEURS DU RESIDU ACTUEL [0,ITERAT]
! OUT XA0    : VALEUR DE L'EXTRAPOLATION
! OUT XA1    : VALEUR DE L'EXTRAPOLATION
! OUT XDET   : VALEUR DE L'EXTRAPOLATION
!
!
!
!
    real(kind=8) :: zero, un, deux
    parameter(zero=0.0d0, un=1.0d0, deux=2.0d0)
!
    integer(kind=8) :: i
    real(kind=8) :: sx, sy, sx2, syx
    real(kind=8) :: xn, xx
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    sx = zero
    sy = zero
    sx2 = zero
    syx = zero
    xn = zero
    xa0 = zero
    xa1 = zero
    xdet = zero
!
! --- CALCUL DE L'EXTRAPOLATION
!
    do i = depart, iterat
        xx = log(vresi(i+1))
        if (i .gt. (iterat-3)) then
            xn = xn+deux
            sx = sx+deux*xx
            sy = sy+deux*i
            sx2 = sx2+deux*(xx**2)
            syx = syx+deux*xx*i
        else
            xn = xn+un
            sx = sx+xx
            sy = sy+i
            sx2 = sx2+xx**2
            syx = syx+xx*i
        end if
    end do
    xdet = -sx**2+sx2*xn
    xa0 = sx2*sy-sx*syx
    xa1 = -(sx*sy)+syx*xn
!
    call jedema()
end subroutine
