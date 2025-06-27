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

module proj_champ_module
!
    implicit none
!
#include "asterf_types.h"
!
! --------------------------------------------------------------------------------------------------
!   Définition des :
!
!   "user_type" utilisés pour les projections
!       prolongation    :
!       distance        :
!
!   "subroutines" utilisés pour les projections
!       prolongation_get
!       prolongation_init
!
! --------------------------------------------------------------------------------------------------
!
    type :: prolongation
        character(len=8)    :: prol_vale_r
        real(kind=8)        :: vale_r
    end type prolongation
    !
!
! ! !     private :: prolongation_write
!
contains
!
    subroutine prolongation_init(prol)
        implicit none
        !
        type(prolongation), intent(inout) :: prol
        ! Initialisation des valeurs par défaut
        prol%prol_vale_r = 'NON'
        prol%vale_r = 0.0
    end subroutine prolongation_init
!
    subroutine prolongation_get(prol)
        implicit none
        !
#include "asterfort/getvtx.h"
#include "asterfort/getvr8.h"
        !
        type(prolongation), intent(inout) :: prol
        !
        integer(kind=8)          :: ispresent
        character(len=8) :: prol0
        real(kind=8)     :: vscal
        !
        ! Initialisation des valeurs par défaut
        call prolongation_init(prol)
        !
        ! Quel est le type de la prolongation
        call getvtx(' ', 'PROL_ZERO', scal=prol0, nbret=ispresent)
        if (ispresent .ne. 0) then
            prol%prol_vale_r = prol0
        end if
        call getvr8(' ', 'PROL_VALE', scal=vscal, nbret=ispresent)
        if (ispresent .ne. 0) then
            prol%prol_vale_r = 'OUI'
            prol%vale_r = vscal
        end if
        ! Pour debug
! ! !         call prolongation_write(prol)
    end subroutine prolongation_get
!
! ! !     subroutine prolongation_write(prol)
! ! !         implicit none
! ! !         !
! ! !         type(prolongation), intent(in) :: prol
! ! !         !
! ! !         write (*, *) 'Prolongation :'
! ! !         write (*, *) '    prol_vale_r : ', prol%prol_vale_r
! ! !         write (*, *) '    vale_r      : ', prol%vale_r
! ! !     end subroutine prolongation_write
!
end module
