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
!
subroutine verins(sddisc, ds_posttimestep)
    use NonLin_Datastructure_type
    implicit none
#include "asterfort/utmess.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utacli.h"
!
    type(NL_DS_PostTimeStep) :: ds_posttimestep
    character(len=19) :: sddisc
    real(kind=8) :: vibr_inst, inst_init, valr(1)
    type(NL_DS_SelectList) :: selectList
    integer(kind=8) :: iinst, nb_nl_inst, nb_found, alarm(2)
    real(kind=8), pointer :: nl_inst_val(:) => null()
    real(kind=8), pointer :: nl_inst_info(:) => null()
    character(len=15) :: mess_alarm(2)
!
! Vérifier les instants demandés par MODE_VIBR sont dans la liste d'instant
! du calcul DYNA_NON_LINE
!
    call jemarq()
    if (ds_posttimestep%l_mode_vibr .or. ds_posttimestep%l_crit_stab) then
! Obtenir les instants obligatoire de DYNA_NON_LINE
        call jeveuo(sddisc(1:19)//'.LIPO', 'L', vr=nl_inst_val)
        call jeveuo(sddisc(1:19)//'.LINF', 'L', vr=nl_inst_info)
        inst_init = nl_inst_val(1)
        nb_nl_inst = nint(nl_inst_info(8))
    end if

    alarm = 0
! Vérifier les instants demandés par CRIT_STAB
    if (ds_posttimestep%l_crit_stab) then
        selectList = ds_posttimestep%crit_stab%selector
        do iinst = 1, selectList%nb_value
            vibr_inst = selectList%list_value(iinst)
            call utacli(vibr_inst, nl_inst_val, nb_nl_inst, &
                        vibr_inst*selectList%precision, nb_found)
!!!! Si un instant est antérieur à l'instant initial
!!!! ou s'il n'est pas dans la lsite de DNL
            if (vibr_inst <= inst_init .or. nb_found == -1) then
                alarm(1) = alarm(1)+1
            end if
        end do
    end if

! Vérifier les instants demandés par MODE_VIBR
    if (ds_posttimestep%l_mode_vibr) then
        selectList = ds_posttimestep%mode_vibr%selector
        do iinst = 1, selectList%nb_value
            vibr_inst = selectList%list_value(iinst)
            call utacli(vibr_inst, nl_inst_val, nb_nl_inst, &
                        vibr_inst*selectList%precision, nb_found)
!!!! Si un instant est antérieur à l'instant initial
!!!! ou s'il n'est pas dans la lsite de DNL
            if (vibr_inst <= inst_init .or. nb_found == -1) then
                alarm(2) = alarm(2)+1
            end if
        end do
    end if

! Emettre le message d'alarme
    mess_alarm = (/' ', ' '/)
    valr(1) = inst_init
    if (alarm(1) .ne. 0) mess_alarm(1) = 'CRIT_STAB'
    if (alarm(2) .ne. 0) mess_alarm(2) = 'MODE_VIBR'
    if (alarm(1) .ne. 0 .or. alarm(2) .ne. 0) then
        call utmess('A', 'UTILITAI8_75', nk=2, valk=mess_alarm, nr=1, valr=valr)
    end if
!
    call jedema()
end subroutine
