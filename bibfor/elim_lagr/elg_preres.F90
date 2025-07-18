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

subroutine elg_preres(solve1, base, iret, matpre, matas1, &
                      npvneg, istop)
!
    use elg_data_module
    use matrasse_module
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/elg_calc_matk_red.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/prere1.h"
!-----------------------------------------------------------------------
! But : faire "preres" si ELIM_MAGR='OUI'
!-----------------------------------------------------------------------
!
    character(len=19) :: matas1, solve1
    character(len=8) :: metres
    character(len=*) :: base, matpre
    integer(kind=8) :: istop, iret
    character(len=19) :: matas2, solve2
    integer(kind=8) ::   npvneg, iexi
    character(len=24), pointer :: slvk(:) => null()
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8) :: nlag1, nlag2, nphys
    aster_logical :: elg_is_ok
!
!
    call jemarq()
!
!   L'élimination des Lagrange est elle licite ?
!   Y a t il des Lagrange dans la matrice ?
!   Y a t il moins de Lagrange que de degrés de liberté "physiques" ?
!   A-t-on bien des double Lagrange ?
    nlag1 = get_num_of_dofs(lagrange1_dof, matas1)
    nlag2 = get_num_of_dofs(lagrange2_dof, matas1)
    nphys = get_num_of_dofs(physical_dof, matas1)
    elg_is_ok = (nlag1 == nlag2) .and. (nlag1 > 0) .and. (nlag1 < nphys)
!
    if (elg_is_ok) then
!
!   -- ON CREE LA MATRICE (REDUITE) MATAS2
!   Si elle existe déjà, on la détruit
        matas2 = "ELG_"//matas1(5:19)
        call jeexin(matas2//'.REFA', iexi)
        if (iexi .gt. 0) then
            call dismoi('METH_RESO', matas2, 'MATR_ASSE', repk=metres)
        end if
        call detrsd("MATR_ASSE", matas2)
        call elg_gest_data('NOTE', matas1, matas2, ' ')
        call elg_calc_matk_red(matas1, solve1, matas2, 'V')
!
!   -- ON DUPLIQUE SOLVE1 EN CHANGEANT ELIM_LAGR: OUI -> NON
        solve2 = "ELG_"//solve1(5:19)
        call copisd('SOLVEUR', 'V', solve1, solve2)
        call jeveuo(solve2//'.SLVK', 'E', vk24=slvk)
        slvk(13) = 'NON'
        call jeveuo(matas2//'.REFA', 'E', vk24=refa)
        refa(7) = solve2
!
!   --  ON APPELLE PRERE1 AVEC MATAS2 ET SOLVE2 :
        call prere1(' ', base, iret, matpre, matas2, &
                    npvneg, istop)
!
    else
!    -- ON APPELLE PRERE1 AVEC MATAS1 ET SOLVE1
!       EN AYANT REMIS ELIM_LAGR A 'NON'
        call jeveuo(solve1//'.SLVK', 'E', vk24=slvk)
        slvk(13) = 'NON'
        call prere1(solve1, base, iret, matpre, matas1, &
                    npvneg, istop)
    end if
!
    call jedema()
end
