! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine ntdoth(model, mater, mateco, cara_elem, list_load, result, &
                  nume_store, matcst_, coecst_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterc/getexm.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ntdoch.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rslesd.h"
!
    character(len=24), intent(out) :: model
    character(len=24), intent(out) :: cara_elem
    character(len=24), intent(out) :: mateco, mater
    character(len=19), intent(inout) :: list_load
    character(len=8), optional, intent(in) :: result
    integer, optional, intent(in) :: nume_store
    aster_logical, optional, intent(out) :: matcst_
    aster_logical, optional, intent(out) :: coecst_
!
! --------------------------------------------------------------------------------------------------
!
! Thermics - Read parameters
!
! --------------------------------------------------------------------------------------------------
!
! Out model            : name of model
! Out list_load        : list of loads
! Out mater            : name of material characteristics (field)
! Out mateco           : name of coded material
! Out cara_elem        : name of datastructure for elementary parameters (CARTE)
! In  result           : name of datastructure for results
! In  nume_store       : index to store in results
! Out matcst           : .true. if constant material parameters
! Out coecst           : .true. if constant rigidity matrix
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: k16bid, nomcmd
    character(len=19) :: list_load_resu
    integer :: i_load, nb_load, iexcit, iret
    character(len=8) :: k8bid, repk, materi
    aster_logical :: l_load_user
    aster_logical :: matcst, coecst
    character(len=24) :: lload_info
    integer, pointer :: v_load_info(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    materi = ' '
    mateco = ' '
    cara_elem = ' '
    if (list_load .eq. ' ') then
        list_load = '&&NTDOTH.LISCHA'
    end if
    model = ' '
    l_load_user = .true.
    coecst = .true.
    matcst = .true.
!
    call getres(k8bid, k16bid, nomcmd)
!
    if (nomcmd .eq. 'CALC_CHAMP') then
        call rslesd(result, nume_store, model, materi, cara_elem, &
                    list_load_resu, iexcit)
        l_load_user = iexcit .eq. 1
    else
        call getvid(' ', 'MODELE', scal=model)
        call getvid(' ', 'CHAM_MATER', scal=materi)
        cara_elem = ' '
        if (getexm(' ', 'CARA_ELEM') .eq. 1) then
            call getvid(' ', 'CARA_ELEM', scal=cara_elem, nbret=iret)
            if (iret .eq. 0) then
                cara_elem = ' '
            end if
        end if
    end if
!
! - Coding material parameters
!
    if (materi .ne. ' ') then
        call rcmfmc(materi, mateco, l_ther_=ASTER_TRUE)
    end if
    mater = materi
!
! - Get loads information and create datastructure
!
    call ntdoch(list_load, l_load_user, list_load_resu, 'V')
!
! - Detect non-constant rigidity matrix
!
    lload_info = list_load(1:19)//'.INFC'
    call jeveuo(lload_info, 'L', vi=v_load_info)
    nb_load = v_load_info(1)
    coecst = .true.
    do i_load = 1, nb_load
        if (v_load_info(nb_load+i_load+1) .eq. 3) then
            coecst = .false.
        end if
    end do
!
! - Detect non-constant material parameters
!
    call dismoi('THER_F_INST', mater, 'CHAM_MATER', repk=repk)
    matcst = repk .eq. 'NON'
!
    if (present(coecst_)) then
        coecst_ = coecst
    end if
    if (present(matcst_)) then
        matcst_ = matcst
    end if
!
end subroutine
