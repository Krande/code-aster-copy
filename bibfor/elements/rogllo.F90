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
subroutine rogllo(nb1, nb2, vrg, blam, ctor, &
                  knn)
!
    implicit none
!
! ......................................................................
!     FONCTION  :  ROTATION DES BLOCS 3 3 DE LA MATRICE DE RIGIDITE
!                  DU REPERE GLOBAL AU REPERE LOCAL
!                  COQUE_3D
!
!           (                                  )         (    U    )
!           ( ( M_TRANSLATION ) (  COUPLAGE )  )         (    V    )
!           (                                  )         (    W    )
! ( VRG ) = (                                  )   ( U )=(         )
!       I   (                                  )         ( THETA_X )
!           ( (   COUPLAGE    ) ( M_ROTATION ) )         ( THETA_Y )
!           (                                  )         ( THETA_Z )
!                                                I                   I
!   ON TOURNE  ( M_ROTATION ) SEULEMENT
!
! ......................................................................
!
#include "asterc/r8prem.h"
#include "asterfort/btkb.h"
    integer(kind=8) :: in
    integer(kind=8) :: ii, jj
    integer(kind=8) :: i, j
    integer(kind=8) :: nb1, nb2
    integer(kind=8) :: irig
!
!---- DECLARATIONS RIGIDITE GEOMETRIQUE
!
    real(kind=8) :: vrg(2601)
    real(kind=8) :: blam(9, 3, 3)
    real(kind=8) :: rigrl(3, 3)
    real(kind=8) :: rigrg(3, 3)
    real(kind=8) :: bid33(3, 3)
    real(kind=8) :: ctor, knn, xmin
!
    real(kind=8) :: barl(3, 3)
!
! DEB
!
!---- A CHAQUE ITERATION
!
    xmin = 1.d0/r8prem()
!
!---- EN CHAQUE NOEUD
!
    do in = 1, nb2
!
!------- ON RECUPERE BARLAMBDA
!
        do jj = 1, 3
            do ii = 1, 3
!
                barl(ii, jj) = blam(in, ii, jj)
!
            end do
        end do
!
!-------    ON CONSTRUIT RIGRG
!
        if (in .le. nb1) then
!
!--------------    NOEUDS DE SERENDIP
            do jj = 1, 3
                do ii = 1, 3
                    j = 6*(in-1)+jj+3
                    i = 6*(in-1)+ii+3
                    irig = (6*nb1+3)*(j-1)+i
                    rigrg(ii, jj) = vrg(irig)
                end do
            end do
!
        else
!
!--------------    SUPERNOEUD
            do jj = 1, 3
                do ii = 1, 3
                    j = 6*nb1+jj
                    i = 6*nb1+ii
                    irig = (6*nb1+3)*(j-1)+i
                    rigrg(ii, jj) = vrg(irig)
                end do
            end do
!
        end if
!
!-------    ROTATION DE RIGRG : LOCALES --> GLOBALES
!
!           RIGRL =  ( LAMBDA0 )   * SIGMT * ( LAMBDA0 ) T
!
!
        call btkb(3, 3, 3, rigrg, barl, &
                  bid33, rigrl)
!
!-------    ON COMPARE LES DEUX PREMIERS TERMES DIAGONAUX DE RIGRL
!
        if (abs(rigrl(1, 1)) .lt. xmin) then
            xmin = abs(rigrl(1, 1))
        end if
        if (abs(rigrl(2, 2)) .lt. xmin) then
            xmin = abs(rigrl(2, 2))
        end if
!
    end do
!
    knn = ctor*xmin
!
! FIN
!
end subroutine
