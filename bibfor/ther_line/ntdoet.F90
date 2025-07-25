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

subroutine ntdoet(model, nume_dof, l_stat, ds_inout)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmdoin.h"
#include "asterfort/nmetl1.h"
#include "asterfort/nmetl2.h"
#include "asterfort/nmetnc.h"
#include "asterfort/ntetl3.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_InOut), intent(inout) :: ds_inout
    aster_logical, intent(out) :: l_stat
!
! --------------------------------------------------------------------------------------------------
!
! Thermics - Init
!
! Read initial state
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  nume_dof         : name of nume_dof object (numbering equation)
! IO  ds_inout         : datastructure for input/output management
! Out l_stat           : .true. if stationnary
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=8) :: calcri
    character(len=24) :: stin_evol, field_type, algo_name, field_algo
    integer(kind=8) :: nb_field, i_field, neq, init_nume
    real(kind=8) :: temp_init, init_time
    aster_logical :: l_stin_evol, l_state_init, l_field_read, l_init_stat, l_init_vale, l_reuse
    real(kind=8), pointer :: vale(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
! - Initializations
!
    call dismoi('NB_EQUA', nume_dof, 'NUME_DDL', repi=neq)
!
! - Model can compute rigidity ?
!
    call dismoi('CALC_RIGI', model, 'MODELE', repk=calcri)
    if (calcri .ne. 'OUI') then
        call utmess('F', 'CALCULEL2_65', sk=model)
    end if
!
! - Parameters from input/output datastructure
!
    nb_field = ds_inout%nb_field
    l_reuse = ds_inout%l_reuse
    l_init_stat = ds_inout%l_init_stat
    l_init_vale = ds_inout%l_init_vale
    temp_init = ds_inout%temp_init
    l_state_init = ds_inout%l_state_init
!
! - PAS D'ETAT INITIAL EN PRESENCE D'UN CONCEPT REENTRANT
!
    if (l_state_init) then
        if (niv .ge. 2) then
            write (ifm, *) '<THER> LECTURE ETAT INITIAL'
        end if
    else
        if (l_reuse) then
            call utmess('A', 'ETATINIT_1')
        else
            call utmess('I', 'ETATINIT_20')
        end if
    end if
!
! - Get name of result datastructure in ETAT_INIT
!
    l_stin_evol = ds_inout%l_stin_evol
    stin_evol = ds_inout%stin_evol
!
! - Initial storing index and time
!
    call nmdoin(ds_inout)
    init_time = ds_inout%init_time
    init_nume = ds_inout%init_nume
!
! - Print
!
    if (niv .ge. 2) then
        write (ifm, *) '<THER> ... INSTANT INITIAL'
        if (init_nume .eq. -1) then
            write (ifm, *) '<THER> ...... NON DEFINI PAR ETAT_INIT'
        else
            write (ifm, *) '<THER> ...... VALEUR    : ', init_time
            write (ifm, *) '<THER> ...... NUME_ORDRE: ', init_nume
        end if
    end if
!
! - No initial state => stationnary
!
    if (.not. l_state_init) then
        l_stat = .true.
        goto 99
    end if
!
! - Loop on fields
!
    do i_field = 1, nb_field
!
! ----- Read field for ETAT_INIT
!
        if (l_stin_evol) then
!
! --------- Initial temperature: from results datastructure
!
            call nmetl1(i_field, ds_inout)
        else
!
! --------- Name of field (type) in results datastructure
!
            field_type = ds_inout%field(i_field)%type
!
! --------- Name of field in algorithm
!
            algo_name = ds_inout%field(i_field)%algo_name
            call nmetnc(algo_name, field_algo)
!
! --------- Informations about initial state
!
            l_field_read = ds_inout%l_field_read(i_field)
            if (field_type .eq. 'TEMP') then
                call jeveuo(field_algo(1:19)//'.VALE', 'E', vr=vale)
!
! ------------- Initial temperature: from field
!
                if (l_field_read) then
                    call nmetl2(model, i_field, ds_inout)
                end if
!
! ------------- Initial temperature: stationnary computation
!
                if (l_init_stat) then
                    l_stat = .true.
                    ds_inout%field(i_field)%init_type = 'STAT'
                end if
!
! ------------- Initial temperature: constant
!
                if (l_init_vale) then
                    vale(1:neq) = temp_init
                    ds_inout%field(i_field)%init_type = 'VALE'
                end if
            else
                call nmetl2(model, i_field, ds_inout)
            end if
        end if
!
! ----- Read field for ETAT_INIT - Some checks
!
        call ntetl3(i_field, ds_inout, temp_init)
    end do
!
99  continue
!
    call jedema()
!
end subroutine
