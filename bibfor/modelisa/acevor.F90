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

subroutine acevor(nbocc, nlg, ier)
!
!
    implicit none
    integer(kind=8) :: nbocc, nlg, ier
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!       Toutes les vérifications sont faites dans le catalogue de la commande
!
! --------------------------------------------------------------------------------------------------
!
! IN  : nbocc  : nombre d'occurence
! OUT : nlg    : nombre total de groupe de maille
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
#include "asterfort/getvtx.h"
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: ng, ioc
! --------------------------------------------------------------------------------------------------
!
    nlg = 0
!
    do ioc = 1, nbocc
        call getvtx('ORIENTATION', 'GROUP_MA', iocc=ioc, nbval=0, nbret=ng)
        nlg = max(nlg, -ng)
    end do
!
end subroutine
