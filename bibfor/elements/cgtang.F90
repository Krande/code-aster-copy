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
subroutine cgtang(ndim, nno, npg, geom, dffr, &
                  t)
    implicit none
    integer(kind=8) :: nno, ndim, npg
    real(kind=8) :: geom(ndim, nno), dffr(nno, npg), t(3, 3)
!
!--------------------------------------------------
!  DEFINITION DES TANGENTES POUR UNE VARIETE 1D PLONGEE DANS DU 3D
!
!  IN  : NDIM : DIMENSION DE L'ESPACE
!        NNO  : NOMBRE DE NOEUDS
!        NPG  : NOMBRE DE POINTS DE GAUSS
!        GEOM : COORDONNEES DES NOEUDS
!        DFFR : DERIVEES DES FONCTIONS DE FORME
!  OUT : T    : TANGENTE AUX NOEUDS
!--------------------------------------------------
    integer(kind=8) :: ig, i, n
    real(kind=8) :: nor
!
!
    do ig = 1, npg
        do i = 1, ndim
            t(i, ig) = 0.d0
            do n = 1, nno
                t(i, ig) = t(i, ig)+geom(i, n)*dffr(n, ig)
            end do
        end do
        nor = 0.d0
        do i = 1, ndim
            nor = nor+t(i, ig)**2
        end do
        nor = sqrt(nor)
        do i = 1, ndim
            t(i, ig) = t(i, ig)/nor
        end do
    end do
!
!
!
end subroutine
