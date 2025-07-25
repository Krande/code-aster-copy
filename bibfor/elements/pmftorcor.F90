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

subroutine pmftorcor(tygrfi, nbpout, gxjx, gxjxpout, deplm, deplp, xl, fl)

    implicit none
!            CORRECTION DES EFFORTS GENERALISES POUR TORSION
!
! -----------------------------------------------------------
! --- IN :
!       typgrfi    : type des fibres : 1 2 ou 3
!       nbpou      : nombre de sous-poutres
!       gxjx       : module de torsion multifibre
!       gxjxpou(*) : Module de torsion pour multipoutres
!       deplm      : champs de deplacement au temps -
!       deplp      : champs de deplacement au temps +
!       xl         : longueur de l'element

! --- OUT
!       fl         : efforts nodaux corriges pour torsion

! -----------------------------------------------------------
    integer(kind=8) :: tygrfi, nbpout, ii
    real(kind=8) :: gxjx, gxjxpout(*), deplm(*), deplp(*), xl, fl(*)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!    parameter (zero=0.0d+0)

    if (tygrfi .eq. 1) then
        fl(10) = gxjx*(deplm(10)+deplp(10)-deplm(4)-deplp(4))/xl
        fl(4) = -fl(10)
    elseif (tygrfi .eq. 2) then
        do ii = 1, nbpout
            fl(10) = fl(10)+gxjxpout(ii)*(deplm(10)+deplp(10)-deplm(4)-deplp(4))/xl
            fl(4) = fl(4)-gxjxpout(ii)*(deplm(10)+deplp(10)-deplm(4)-deplp(4))/xl
        end do
    elseif (tygrfi .eq. 3) then
        do ii = 1, nbpout
            fl(13) = fl(13)+gxjxpout(ii)*(deplm(13)+deplp(13)-deplm(4)-deplp(4))/xl
            fl(4) = fl(4)-gxjxpout(ii)*(deplm(13)+deplp(13)-deplm(4)-deplp(4))/xl
        end do
    end if

end subroutine
