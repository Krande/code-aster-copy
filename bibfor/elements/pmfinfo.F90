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

subroutine pmfinfo(nbfibr, nbgrfi, tygrfi, nbcarm, nug, jacf, nbassfi)
!
!
! --------------------------------------------------------------------------------------------------
!
!           Informations sur les PMF
!
! person_in_charge: jean-luc.flejou at edf.fr
! --------------------------------------------------------------------------------------------------
!
!   OUT
!       nbfibr      : nombre total de fibre
!       nbgrfi      : nombre de groupe de fibres
!       tygrfi      : type des groupes de fibres
!       nbcarm      : nombre de composantes dans la carte
!       nug         : numéro des groupes de fibres nug(1:nbgrfi)
!       jacf        : pointeur sur les caractéristiques de fibres
!       nbassfi     : nombre d'assemblage de fibre
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(out) :: nbfibr, nbgrfi, tygrfi, nbcarm, nug(*)
    integer(kind=8), intent(out), optional :: jacf, nbassfi
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jnbspi, ii, numgr, jacf_loc
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PNBSP_I', 'L', jnbspi)
    nbfibr = zi(jnbspi)
    nbgrfi = zi(jnbspi+1)
    tygrfi = zi(jnbspi+2)
    nbcarm = zi(jnbspi+3)
    nug(1:nbgrfi) = zi(jnbspi+3+1:jnbspi+3+nbgrfi)
!
    if (tygrfi .eq. 1) then
        if (present(jacf)) then
            call jevech('PFIBRES', 'L', jacf)
        end if
        if (present(nbassfi)) then
            nbassfi = 1
        end if
    else if (tygrfi .eq. 2) then
        if (present(jacf)) then
            call jevech('PFIBRES', 'L', jacf_loc)
            jacf = jacf_loc
        end if
        if (present(nbassfi)) then
            if (absent(jacf)) then
                call jevech('PFIBRES', 'L', jacf_loc)
            end if
            nbassfi = 0
            do ii = 1, nbfibr
                numgr = nint(zr(jacf_loc-1+ii*nbcarm))
                nbassfi = max(nbassfi, numgr)
            end do
            ASSERT(nbassfi .ne. 0)
        end if
    else
        call utmess('F', 'ELEMENTS2_40', si=tygrfi)
    end if
end subroutine
