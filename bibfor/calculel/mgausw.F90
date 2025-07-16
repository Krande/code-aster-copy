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
subroutine mgausw(a, b, dim, nordre, nb, &
                  det, iret)
!
    implicit none
#include "asterf_types.h"
!
    integer(kind=8) :: dim, nb, nordre
    real(kind=8) :: a(dim, dim), b(dim, nb), det
    aster_logical :: iret
!
!
! ----------------------------------------------------------------------
!     RESOLUTION PAR LA METHODE DE GAUSS D'UN SYSTEME LINEAIRE
! ----------------------------------------------------------------------
!     VARIABLES D'ENTREE
!     REAL*8      A(DIM, DIM)     : MATRICE CARREE PLEINE
!     REAL*8      B(DIM, NB)      : SECONDS MEMBRES
!     INTEGER     DIM             : DIMENSION DE A
!     INTEGER     NORDRE          : RANG DE LA MATRICE
!     INTEGER     NB              : NOMBRE DE SECONDS MEMBRES
!     REAL*8      DET             : 0 = PAS DE CALCUL DU DETERMINANT
!
!     VARIABLES DE SORTIE
!     REAL*8      B(DIM, NB)      : A-1 * B
!     REAL*8      DET             : DETERMINANT DE A (SI DEMANDE)
!     LOGICAL     IRET            : .FALSE. SI A SINGULIERE
!
! ----------------------------------------------------------------------
!     ATTENTION : LA MATRICE A EST MODIFIEE
! ----------------------------------------------------------------------
!
!     PARAMETRE
    real(kind=8) :: condmx
    parameter(condmx=1.d16)
!
    integer(kind=8) :: i, j, k
    real(kind=8) :: c, d, cmin, cmax
    aster_logical :: flag, ldet
!
    iret = .true.
!
    if (det .eq. 0.d0) then
        ldet = .false.
    else
        ldet = .true.
        det = 1.d0
    end if
!
    do i = 1, nordre
!
! ----- RECHERCHE DU MEILLEUR PIVOT
!
        j = i
        c = a(i, i)
        flag = .false.
!
        do k = i+1, nordre
            d = a(k, i)
            if (abs(c) .lt. abs(d)) then
                c = d
                j = k
                flag = .true.
            end if
        end do
!
! ----- DETERMINANT
!
        if (ldet) det = det*c
!
! ----- ESTIMATION GROSSIERE DU CONDITIONNEMENT
!
        if (i .eq. 1) then
            cmin = abs(c)
            cmax = cmin
        else
            if (abs(c) .lt. cmin) then
                cmin = abs(c)
                if (cmax .gt. condmx*cmin) then
                    iret = .false.
                    goto 100
                end if
                goto 30
            end if
            if (abs(c) .gt. cmax) then
                cmax = abs(c)
                if (cmax .gt. condmx*cmin) then
                    iret = .false.
                    goto 100
                end if
            end if
        end if
!
30      continue
!
! ----- PERMUTATION
!
        if (flag) then
!
            do k = i, nordre
                d = a(i, k)
                a(i, k) = a(j, k)
                a(j, k) = d
            end do
!
            do k = 1, nb
                d = b(i, k)
                b(i, k) = b(j, k)
                b(j, k) = d
            end do
!
            det = (-1.d0)*det
!
        end if
!
! ----- ELIMINATION
!
        do j = i+1, nordre
!
            if (a(j, i) .ne. 0.d0) then
!
                d = a(j, i)/c
!
                do k = 1, nb
                    b(j, k) = b(j, k)-d*b(i, k)
                end do
!
                do k = i+1, nordre
                    a(j, k) = a(j, k)-d*a(i, k)
                end do
!
            end if
!
        end do
    end do
!
! --- RESOLUTION
!
    do k = 1, nb
        b(nordre, k) = b(nordre, k)/c
!
        do i = nordre-1, 1, -1
            d = 0.d0
            do j = i+1, nordre
                d = d+a(i, j)*b(j, k)
            end do
!
            b(i, k) = (b(i, k)-d)/a(i, i)
!
        end do
    end do
!
100 continue
!
end subroutine
