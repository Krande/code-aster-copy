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

subroutine delete_matrix(matas, typsol)
!
    use elg_data_module
!
#include "asterf_types.h"
!
! person_in_charge: nicolas.sellenet at edf.fr
!
    implicit none
    character(len=19) :: matas
    character(len=5) :: typsol
!
#include "asterfort/amumph.h"
#include "asterfort/apetsc.h"
#include "asterfort/detlsp.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/detrsd.h"
#include "asterfort/jeveuo.h"
#include "jeveux.h"

    character(len=19) :: matas2, solveu2
    integer(kind=8) :: iexi, ibid, iret, iexi2
    complex(kind=8) :: cbid
    aster_logical :: lbid
    character(len=24), pointer :: refa(:) => null()
!
    cbid = dcmplx(0.d0, 0.d0)

!
!----------------------------------------------------------------
!
!  DELETE MATRIX FOR MUMPS AND PETSC
!
!----------------------------------------------------------------
!
    call jemarq()

    call jeexin(matas//'.REFA', iexi)
    if (iexi .gt. 0) then
!       -- DESTRUCTION DE L'EVENTUELLE MATRICE RéDUITE (ELIM_LAGR) :
        call jeveuo(matas//'.REFA', 'L', vk24=refa)
        if (refa(19) .ne. ' ') then
            matas2 = refa(19) (1:19)
            call jeexin(matas2//".REFA", iexi2)
            if (iexi2 .gt. 0) then
                call dismoi('SOLVEUR', matas2, 'MATR_ASSE', repk=solveu2, arret='C', &
                            ier=iret)
                if (iret .eq. 0) then
                    call detlsp(matas2, solveu2)
                end if
                call detrsd('MATR_ASSE', matas2)
                call elg_gest_data('EFFACE', ' ', matas2, ' ')
            end if
        end if
! --- Destruction du solveur
        if (typsol .eq. 'MUMPS') then
            call amumph('DETR_MAT', ' ', matas, [0.d0], [cbid], &
                        ' ', 0, ibid, lbid)
        else if (typsol .eq. 'PETSC') then
            call apetsc('DETR_MAT', ' ', matas, [0.d0], ' ', &
                        0, ibid, iret)
        end if

    end if
    call jedema()

!
end subroutine
