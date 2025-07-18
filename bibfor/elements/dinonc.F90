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

subroutine dinonc(nomte, icodre, valre, klv, raide, &
                  nbpar, param, okdire)
    implicit none
#include "asterf_types.h"
#include "asterfort/utmess.h"
    character(len=16) :: nomte
    integer(kind=8) :: icodre(*)
    integer(kind=8) :: nbpar
    real(kind=8) :: valre(*), klv(*), raide(*), param(6, nbpar)
    aster_logical :: okdire(6)
!
! person_in_charge: jean-luc.flejou at edf.fr
! --------------------------------------------------------------------------------------------------
!           AFFECTATION DES VALEURS ISSUES DU COMPORTEMENT
!
!   si on essaye d'affecter un comportement sur un ddl non autorise on sort en 'f'
!
!   pour que cela fonctionne correctement il faut que les paramètres soient ranges dans le
!   data 'nomre' de la façon suivante
!
!   PARAMETRES SUIVANT X        PARAMETRES SUIVANT Y        ETC
!   P_1_DX  P_2_DX  P_3_DX ...  P_1_DY  P_2_DY  P_3_DY ...
!
!     NOMRE  /'FLIM_X','PUIS_DX',
!             'FLIM_Y','PUIS_DY',
!             'FLIM_Z','PUIS_DZ',
!             'MLIM_X','PUIS_RX',
!             'MLIM_Y','PUIS_RY',
!             'MLIM_Z','PUIS_RZ'/
! --------------------------------------------------------------------------------------------------
!
!  IN
!     nomte : nom de l'élément
!     icodre : 0 si le coeff est présent sinon 1
!     valre : valeur des coefficients
!     klv   : raideur élastique du discret
!     nbpar : nombre de paramètre de la loi par ddl
!
!  OUT
!     raide  : raideur au comportement
!     param  : paramètres de la loi
!     okdire : vrai si la direction est affectée par le comportement
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: ii, jj
!
    do ii = 1, 6
        okdire(ii) = .false.
    end do
!
    if ((nomte .eq. 'MECA_DIS_TR_N') .or. (nomte .eq. 'MECA_DIS_TR_L')) then
        do ii = 0, 5
            do jj = 1, nbpar
                if (icodre(nbpar*ii+jj) .eq. 0) then
                    param(ii+1, jj) = valre(nbpar*ii+jj)
                    okdire(ii+1) = .true.
                end if
            end do
        end do
        raide(1) = klv(1)
        raide(2) = klv(3)
        raide(3) = klv(6)
        raide(4) = klv(10)
        raide(5) = klv(15)
        raide(6) = klv(21)
    end if
    if ((nomte .eq. 'MECA_DIS_T_N') .or. (nomte .eq. 'MECA_DIS_T_L')) then
        do ii = 0, 2
            do jj = 1, nbpar
                if (icodre(nbpar*ii+jj) .eq. 0) then
                    param(ii+1, jj) = valre(nbpar*ii+jj)
                    okdire(ii+1) = .true.
                end if
            end do
        end do
        do ii = 3, 5
            do jj = 1, nbpar
                if (icodre(nbpar*ii+jj) .eq. 0) then
                    call utmess('F', 'DISCRETS_1', sk=nomte)
                end if
            end do
        end do
        raide(1) = klv(1)
        raide(2) = klv(3)
        raide(3) = klv(6)
    end if
    if ((nomte .eq. 'MECA_2D_DIS_TR_N') .or. (nomte .eq. 'MECA_2D_DIS_TR_L')) then
        do ii = 0, 1
            do jj = 1, nbpar
                if (icodre(nbpar*ii+jj) .eq. 0) then
                    param(ii+1, jj) = valre(nbpar*ii+jj)
                    okdire(ii+1) = .true.
                end if
            end do
        end do
        ii = 5
        do jj = 1, nbpar
            if (icodre(nbpar*ii+jj) .eq. 0) then
                param(3, jj) = valre(nbpar*ii+jj)
                okdire(3) = .true.
            end if
        end do
        do ii = 2, 4
            do jj = 1, nbpar
                if (icodre(nbpar*ii+jj) .eq. 0) then
                    call utmess('F', 'DISCRETS_2', sk=nomte)
                end if
            end do
        end do
        raide(1) = klv(1)
        raide(2) = klv(3)
        raide(3) = klv(6)
    end if
    if ((nomte .eq. 'MECA_2D_DIS_T_N') .or. (nomte .eq. 'MECA_2D_DIS_T_L')) then
        do ii = 0, 1
            do jj = 1, nbpar
                if (icodre(nbpar*ii+jj) .eq. 0) then
                    param(ii+1, jj) = valre(nbpar*ii+jj)
                    okdire(ii+1) = .true.
                end if
            end do
        end do
        do ii = 2, 5
            do jj = 1, nbpar
                if (icodre(nbpar*ii+jj) .eq. 0) then
                    call utmess('F', 'DISCRETS_3', sk=nomte)
                end if
            end do
        end do
        raide(1) = klv(1)
        raide(2) = klv(3)
    end if
!
end subroutine
