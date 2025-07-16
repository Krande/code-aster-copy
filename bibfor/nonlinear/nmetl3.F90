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
subroutine nmetl3(model, compor, i_field, ds_inout, verbose)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/chpver.h"
#include "asterfort/dismoi.h"
#include "asterfort/nmetnc.h"
#include "asterfort/nmsigi.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcomp.h"
#include "asterfort/vrcom2.h"
#include "asterfort/stressChck.h"
!
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: compor
    type(NL_DS_InOut), intent(in) :: ds_inout
    integer(kind=8), intent(in) :: i_field
    aster_logical, intent(in) :: verbose
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Read field for ETAT_INIT - Some checks
!
! --------------------------------------------------------------------------------------------------
!
! In  compor           : name of <CARTE> COMPOR
! In  model            : name of model
! In  i_field          : field index
! In  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=24) :: fieldType, algo_name, fieldRefe
    character(len=4) :: fieldDisc
    character(len=8) :: gran_name
    character(len=24) :: field_algo
    character(len=24) :: modelLigrel, comporPrev
    character(len=4) :: init_type
    aster_logical :: l_state_init, l_stin_evol, l_acti, lModiVari
    integer(kind=8) :: init_nume
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! - Get parameters
    l_stin_evol = ds_inout%l_stin_evol
    init_nume = ds_inout%init_nume
    l_state_init = ds_inout%l_state_init
!
! - Field to read ?
!
    if (ds_inout%field(i_field)%l_read_init) then
! ----- Name of field (type) in results datastructure
        fieldType = ds_inout%field(i_field)%type

! ----- Name of field for initial state
        fieldRefe = ds_inout%field(i_field)%init_name

! ----- Name of field in algorithm
        algo_name = ds_inout%field(i_field)%algo_name
        call nmetnc(algo_name, field_algo)

! ----- Spatial discretization of field
        fieldDisc = ds_inout%field(i_field)%disc_type

! ----- Actual state of field
        init_type = ds_inout%field(i_field)%init_type

! ----- Type of GRANDEUR of field
        gran_name = ds_inout%field(i_field)%gran_name

! ----- Is field should been active ?
        l_acti = ds_inout%l_field_acti(i_field)

! ----- Informations about field
        if (l_acti) then
            if (init_type .eq. ' ') then
                call utmess('F', 'ETATINIT_30', sk=fieldType)
            else
                if (init_type .eq. 'ZERO') then
                    call utmess('I', 'ETATINIT_31', sk=fieldType)
                else if (init_type .eq. 'RESU') then
                    call utmess('I', 'ETATINIT_32', sk=fieldType)
                else if (init_type .eq. 'READ') then
                    call utmess('I', 'ETATINIT_33', sk=fieldType)
                else
                    ASSERT(.false.)
                end if
            end if

! --------- Check GRANDEUR and discretization
            if (gran_name .ne. ' ') then
                call chpver('F', field_algo, fieldDisc, gran_name, iret)
            end if

!  -------- Check stresses
            if (fieldType .eq. 'SIEF_ELGA') then
                if (l_state_init) then
                    call stressChck(fieldRefe, field_algo, ASTER_TRUE, iret)
                    if (iret .eq. 1) then
                        call utmess('A', 'MECANONLINE5_81')
                    end if
                end if
            end if

! --------- For pre-stressed load
            if (fieldType .eq. 'SIEF_ELGA') then
                call nmsigi(modelLigrel, compor, field_algo(1:19))
            end if

!  -------- Check internal variables
            if (fieldType .eq. 'VARI_ELGA') then
                if (l_state_init) then
                    comporPrev = ' '
                    if (l_stin_evol) then
                        call rsexch(' ', ds_inout%stin_evol, 'COMPORTEMENT', init_nume, &
                                    comporPrev, iret)
                        if (iret .ne. 0) then
                            comporPrev = ' '
                        end if
                    end if
                    if (comporPrev .eq. ' ') then
                        call vrcomp(compor, field_algo, modelLigrel, iret, verbose_=verbose, &
                                    lModiVari_=lModiVari)
                    else
                        call vrcomp(compor, field_algo, modelLigrel, iret, &
                                    comporPrevz_=comporPrev, verbose_=verbose, &
                                    lModiVari_=lModiVari)
                    end if
                    if (iret .eq. 1) then
                        call utmess('F', 'MECANONLINE5_82')
                    end if
                    if (lModiVari) then
                        call vrcom2(compor, field_algo, modelLigrel, ASTER_FALSE)
                    end if
                end if
            end if
        else
            if (init_type .eq. 'RESU' .or. init_type .eq. 'READ') then
                call utmess('I', 'ETATINIT_36', sk=fieldType)
            end if
        end if
    end if
!
end subroutine
