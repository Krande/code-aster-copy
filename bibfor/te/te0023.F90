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
subroutine te0023(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
! IN OPTION    : K16 :  OPTION DE CALCUL
!     'INI_STRX'
! IN NOMTE     : K16 : NOM DU TYPE ELEMENT
!     POUTRE
!        'MECA_POU_D_E'  'MECA_POU_D_T'  'MECA_POU_D_TG'
!        'MECA_POU_D_EM' 'MECA_POU_D_TGM'
!
! INITIALISATION DU CHAMP STRX_ELGA
!
    integer(kind=8) :: iorien, istrx, i, kpg, npg, ncomp
!     ------------------------------------------------------------------
!
    ASSERT(option .eq. 'INI_STRX')
!
    call jevech('PCAORIE', 'L', iorien)
    call jevech('PSTRX_R', 'E', istrx)
!
    if (nomte .eq. 'MECA_POU_D_TGM' .or. nomte .eq. 'MECA_POU_D_EM' .or. nomte .eq. &
        'MECA_POU_D_SQUE') then
        if (nomte .eq. 'MECA_POU_D_EM' .or. nomte .eq. 'MECA_POU_D_SQUE') then
            npg = 2
        else
            npg = 3
        end if
        if (nomte .eq. 'MECA_POU_D_SQUE') then
            ncomp = 21
        else
            ncomp = 18
        end if
        do kpg = 1, npg
            do i = 1, 15
                zr(istrx-1+ncomp*(kpg-1)+i) = 0.d0
            end do
            do i = 1, 3
                zr(istrx-1+ncomp*(kpg-1)+15+i) = zr(iorien-1+i)
            end do
        end do
    else
        do i = 1, 3
            zr(istrx-1+i) = zr(iorien-1+i)
        end do
    end if
end subroutine
