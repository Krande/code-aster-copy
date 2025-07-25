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
subroutine mexthr(n, a, lda)
    implicit none
#include "asterc/r8prem.h"
    integer(kind=8) :: n, lda
    complex(kind=8) :: a(lda, *)
! IN  : N    : DIMENSION DE LA MATRICE.
!     : LDA  : DIMENSION DE A.
! I/O : A    : MATRICE HERMITIENNE COMPLEXE D'ORDRE N.
!         IN : PARTIE TRIANGULAIRE SUPERIEURE.
!        OUT   PARTIE TRIANGULAIRE INFERIEURE DE LA MATRICE A DEFINIE
!              COMME UNE MATRICE COMPLEXE HERMITIENNE.
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
    real(kind=8) :: eps
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    if (n .lt. 1) then
        write (6, *) 'THE ARGUMENT N = %(I1).  IT MUST BE AT '//&
     &               'LEAST 1.'
        goto 9000
    end if
    if (lda .lt. n) then
        write (6, *) 'THE ARGUMENT LDA = %(I1).  IT MUST BE AT '//&
     &               'LEAST AS LARGE AS N = %(I2).'
        goto 9000
    end if
    eps = 10.0d0*r8prem()
    do i = 1, n
        if (abs(dimag(a(i, i))) .ne. 0.0d0) then
            if (abs(dimag(a(i, i))) .gt. eps*abs(dble(a(i, i)))) then
                write (6, *) a(i, i)
                write (6, *) 'THE MATRIX ELEMENT A(%(I1),%(I1)) '//&
     &                     '= %(Z1).  THE DIAGONAL OF A HERMITIAN '//&
     &                     'MATRIX MUST BE REAL.'
                goto 9000
            else
                write (6, *) a(i, i)
                write (6, *) 'THE MATRIX ELEMENT A(%(I1),%(I1)) '//&
     &                     '= %(Z1).  THE DIAGONAL OF A HERMITIAN '//&
     &                     'MATRIX MUST BE REAL: ITS IMAGINARY PART '//&
     &                     'IS SET TO ZERO.'
                a(i, i) = dble(a(i, i))
            end if
        end if
    end do
!
    do j = 1, n-1
        do i = j+1, n
            a(i, j) = dconjg(a(j, i))
        end do
    end do
!
9000 continue
end subroutine
