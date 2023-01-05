! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine deflg3(gn, feta, xi, me, t, &
                  tl)
    implicit none
!     CALCUL DES DEFORMATIONS LOGARITHMIQUES ET DES TERMES NECESSAIRES
!     AU POST TRAITEMENT DES CONTRAINTES ET A LA RIGIDITE TANGENTE
!     SUIVANT ARTICLE MIEHE APEL LAMBRECHT CMAME 2002
! ----------------------------------------------------------------------
!     IN GN    directions propres du tenseur F
!     IN FETA  utilitaires issus de DEFLG2  f_i=-2/lambda_i**2 puis eta
!     IN XI    utilitaires issus de DEFLG2  xi_ij
!     IN ME    utilitaires issus de DEFLG2  tenseur M d'ordre 4
!     IN T     tenseur des contraintesissu de NMCOMP (avec sqrt(2))
!     OUT TL   tenseur d'ordre 4 T:L
! ----------------------------------------------------------------------
#include "asterfort/r8inir.h"
#include "asterfort/tnsvec.h"
    real(kind=8) :: gn(3, 3), t(6), tl(3, 3, 3, 3)
    real(kind=8) :: dzeta(3, 3), t33(3, 3), me(3, 3, 3, 3), xi(3, 3), feta(4)
    integer :: i, j, k, a, b, c, d
! ----------------------------------------------------------------------
!
!     CALCUL DU TERME T.L
!
    call r8inir(81, 0.d0, tl, 1)
    call r8inir(9, 0.d0, dzeta, 1)
    call tnsvec(6, 3, t33, t, 1.d0/sqrt(2.d0))
!
!     A,B sont les composantes, J,I sont les modes propres
    do i = 1, 3
        do j = 1, 3
            do a = 1, 3
                do b = 1, 3
                    dzeta(i, j) = dzeta(i, j)+t33(a, b)*gn(a, i)*gn(b, j)
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do a = 1, 3
            do b = 1, 3
                do c = 1, 3
                    do d = 1, 3
                        tl(a, b, c, d) = tl(a, b, c, d)+0.25d0*feta(i)*dzeta( &
                                         i, i)*me(a, b, i, i)*me(c, d, i, i)
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do a = 1, 3
                    do b = 1, 3
                        do c = 1, 3
                            do d = 1, 3
                                if ((j .ne. i) .and. (j .ne. k) .and. (k .ne. i)) then
                                    tl(a, b, c, d) = tl(a, b, c, d)+2.d0* &
                                                     feta(4)*dzeta(i, j)*me(a, b, i, k)*me( &
                                                     c, d, j, k)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            do a = 1, 3
                do b = 1, 3
                    do c = 1, 3
                        do d = 1, 3
                            if (j .ne. i) then
                                tl(a, b, c, d) = tl(a, b, c, d)+2.d0*xi(i, j)* &
                                                 dzeta(i, j)*me(a, b, i, j)*me(c, d, j, j)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            do a = 1, 3
                do b = 1, 3
                    do c = 1, 3
                        do d = 1, 3
                            if (j .ne. i) then
                                tl(a, b, c, d) = tl(a, b, c, d)+2.d0*xi(i, j)* &
                                                 dzeta(i, j)*me(a, b, j, j)*me(c, d, i, j)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            do a = 1, 3
                do b = 1, 3
                    do c = 1, 3
                        do d = 1, 3
                            if (j .ne. i) then
                                tl(a, b, c, d) = tl(a, b, c, d)+2.d0*xi(i, j)* &
                                                 dzeta(j, j)*me(a, b, i, j)*me(c, d, i, j)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
!
end subroutine
