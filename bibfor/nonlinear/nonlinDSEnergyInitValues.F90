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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nonlinDSEnergyInitValues(ds_energy, stin_evol, ds_inout)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/exisd.h"
#include "asterfort/IncrEnergy.h"
#include "asterfort/ltnotb.h"
#include "asterfort/tbGetListPara.h"
#include "asterfort/tbliva.h"
!
    type(NL_DS_Energy), intent(inout) :: ds_energy
    character(len=8), intent(in) :: stin_evol
    type(NL_DS_InOut), intent(inout) :: ds_inout
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Energy management
!
! Initializations for energy values from initial state data structure
!
! --------------------------------------------------------------------------------------------------
!

!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_para, iret, nb_line, i, ibid
    real(kind=8) ::  valr, precision, init_time
    complex(kind=8) :: cbid
    character(len=8) :: k8b, ctype, criterion
    character(len=19) :: nomtab
    character(len=24) :: para
    character(len=24), pointer :: paraName(:) => null()
    character(len=24), pointer :: paraType(:) => null()
!
! --------------------------------------------------------------------------------------------------

    nomtab = ' '
    call ltnotb(stin_evol, 'PARA_CALC', nomtab)
    call exisd('TABLE', nomtab, iret)
    ASSERT(iret .ne. 0)

    init_time = ds_inout%init_time
    precision = ds_inout%precision
    criterion = ds_inout%criterion

    call tbGetListPara(nomtab, nb_para, paraName, paraType, nb_line)

    do i = 1, nb_para
        if (paraName(i) .eq. 'NUME_REUSE') cycle
        if (paraName(i) .eq. 'INST') cycle
        valr = 0.d0
        para = paraName(i)
        call tbliva(nomtab, 1, ['INST'], [ibid], [init_time], &
                    [cbid], [k8b], [criterion], [precision], para, &
                    ctype, ibid, valr, cbid, k8b, &
                    iret)
        call IncrEnergy(ds_energy, para, valr)
    end do

    AS_DEALLOCATE(vk24=paraType)
    AS_DEALLOCATE(vk24=paraName)

end subroutine
