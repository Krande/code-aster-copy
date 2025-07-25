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
subroutine nxdomt(ds_algopara, ds_algorom)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/romAlgoNLRead.h"
!
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
!
! --------------------------------------------------------------------------------------------------
!
! THER_NON_LINE - Algorithm management
!
! Read parameters for algorithm management
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_algopara      : datastructure for algorithm parameters
! IO  ds_algorom       : datastructure for ROM parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: reac_iter, iter_line_maxi
    real(kind=8) :: resi_line_rela
    character(len=16) :: keywf, algo_meth
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<THERNONLINE> . Read parameters for algorithm parameters'
    end if
!
! - Initializations
!
    algo_meth = ' '
    reac_iter = 0
    iter_line_maxi = 0
    resi_line_rela = 1.d-3
!
! - Get method
!
    call getvtx(' ', 'METHODE', scal=algo_meth)
    ds_algopara%method = algo_meth
!
! - Get parameters of method
!
    if ((algo_meth .eq. 'NEWTON') .or. (algo_meth .eq. 'NEWTON_KRYLOV')) then
        keywf = 'NEWTON'
        call getvis(keywf, 'REAC_ITER', iocc=1, scal=reac_iter)
        ASSERT(reac_iter .ge. 0)
        ds_algopara%reac_iter = reac_iter
        call getvr8(keywf, 'RESI_LINE_RELA', iocc=1, scal=resi_line_rela)
        call getvis(keywf, 'ITER_LINE_MAXI', iocc=1, scal=iter_line_maxi)
        ds_algopara%line_search%resi_rela = resi_line_rela
        ds_algopara%line_search%iter_maxi = iter_line_maxi
    else if (algo_meth .eq. 'MODELE_REDUIT') then
        keywf = 'MODELE_REDUIT'
        call getvis(keywf, 'REAC_ITER', iocc=1, scal=reac_iter)
        ASSERT(reac_iter .ge. 0)
        ds_algopara%reac_iter = reac_iter
        call romAlgoNLRead(ds_algorom)
    else
        ASSERT(.false.)
    end if
!
end subroutine
