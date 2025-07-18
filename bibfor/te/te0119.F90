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

subroutine te0119(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/teattr.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!                                   OPTION VERI_CARA_ELEM
!
!
!   Vérification du contenu des cartes sur les éléments
!
!
!   Remarque : Si la vérification est faite içi, il est inutile de la faire dans "veri_affe_carte"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: j1, ibid, iadzi, iazk24
    real(kind=8)        :: excent
    character(len=3)    :: cmod
    character(len=8)    :: alias8
    character(len=24)   :: valk(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Récupération du code de la modélisation
    call teattr('S', 'ALIAS8', alias8, ibid)
    cmod = alias8(3:5)
!
!   Vérification que l'excentrement est nul pour COQUE_3D
    if (cmod .eq. 'CQ3') then
        call jevech('PCACOQU', 'L', j1)
        excent = zr(j1-1+6)
        if (nint(excent) .ne. 0) then
            call tecael(iadzi, iazk24)
            valk(1) = zk24(iazk24-1+3) (1:8)
            call utmess('F', 'CALCULEL2_31', sk=valk(1))
        end if
    end if
!
end subroutine
