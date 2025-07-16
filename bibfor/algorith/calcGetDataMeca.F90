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
subroutine calcGetDataMeca(listLoad, model, materField, mateco, caraElem, &
                           disp_prev, disp_cumu_inst, vari_prev, sigm_prev, &
                           ds_constitutive, l_elem_nonl, nume_harm)
!
    use NonLin_Datastructure_type
    use listLoad_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/chamnoIsSame.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/isOptionPossible.h"
#include "asterfort/nmdoch.h"
#include "asterfort/nmdorc.h"
#include "asterfort/nmlect.h"
#include "asterfort/nonlinDSConstitutiveInit.h"
#include "asterfort/utmess.h"
!
    character(len=24), intent(out) :: listLoad, mateco
    character(len=8), intent(out) :: model, materField, caraElem
    character(len=19), intent(out) :: disp_prev
    character(len=19), intent(out) :: disp_cumu_inst
    character(len=19), intent(out) :: vari_prev
    character(len=19), intent(out) :: sigm_prev
    type(NL_DS_Constitutive), intent(out) :: ds_constitutive
    aster_logical, intent(out) :: l_elem_nonl
    integer(kind=8), intent(out) :: nume_harm
!
! --------------------------------------------------------------------------------------------------
!
! Command CALCUL
!
! Get data for mechanics
!
! --------------------------------------------------------------------------------------------------
!
! Out listLoad         : name of datastructure for list of loads
! Out model            : name of model
! Out mateco           : name of coded material
! Out caraElem         : name of elementary characteristics (field)
! Out disp_prev        : displacement at beginning of step
! Out disp_cumu_inst   : displacement increment from beginning of step
! Out vari_prev        : internal variables at beginning of step
! Out sigm_prev        : stress at beginning of step
! Out ds_constitutive  : datastructure for constitutive laws management
! Out l_elem_nonl      : .true. if all elements can compute non-linear options
! Out nume_harm        : Fourier harmonic number
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: result
    aster_logical :: l_etat_init, verbose
    integer(kind=8) :: nocc, ier
    character(len=19) :: ligrmo
    type(ListLoad_Prep) :: listLoadPrep
!
! --------------------------------------------------------------------------------------------------
!
    listLoad = '&&OP0026.LISCHA'
    caraElem = ' '
    model = ' '
    mateco = ' '
    vari_prev = ' '
    sigm_prev = ' '
    disp_prev = ' '
    disp_cumu_inst = ' '
    l_elem_nonl = ASTER_FALSE

! - Get parameters from command file
    call nmlect(result, model, materField, mateco, caraElem)

! - Get loads/BC and create list of loads datastructure
    listLoad = '&&OP0026.LIST_LOAD'
    listLoadPrep%model = model
    listLoadPrep%lHasPilo = ASTER_FALSE
    listLoadPrep%funcIsCplx = ASTER_FALSE
    listLoadPrep%staticOperator = ASTER_TRUE
    call nmdoch(listLoadPrep, listLoad, "V")

! - Can have internal variables ?
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)
    call isOptionPossible(ligrmo, 'TOU_INI_ELGA', 'PVARI_R', l_some_=l_elem_nonl)

! - Does option FULL_MECA exist ?
    if (l_elem_nonl) then
        call isOptionPossible(ligrmo, 'FULL_MECA', 'PDEPLPR', l_some_=l_elem_nonl)
    end if
!
! - Get displacements
!
    call getvid(' ', 'DEPL', scal=disp_prev, nbret=nocc)
    if (nocc .eq. 0) then
        disp_prev = ' '
    end if
    call getvid(' ', 'INCR_DEPL', scal=disp_cumu_inst, nbret=nocc)
    if (nocc .eq. 0) then
        disp_cumu_inst = ' '
    end if

    if (disp_prev .ne. " " .and. disp_cumu_inst .ne. " ") then
        call chamnoIsSame(disp_prev, disp_cumu_inst, ier)
        if (ier .gt. 0) then
            call utmess('F', 'CALCUL1_14')
        end if
    end if
!
! - Get stresses
!
    call getvid(' ', 'SIGM', scal=sigm_prev, nbret=nocc)
    l_etat_init = nocc .ne. 0
    if (nocc .eq. 0) then
        sigm_prev = ' '
    end if
!
! - Get internal variables
!
    call getvid(' ', 'VARI', scal=vari_prev, nbret=nocc)
    if (nocc .eq. 0) then
        vari_prev = ' '
    end if
    if (vari_prev .ne. ' ' .and. .not. l_elem_nonl) then
        call utmess('I', 'CALCUL1_7')
    end if
!
! - Get Fourier Mode
!
    call getvis(' ', 'MODE_FOURIER', scal=nume_harm, nbret=nocc)
    if (nocc .eq. 0) then
        nume_harm = 0
    end if
!
! - Prepare constitutive laws management datastructure
!
    if (l_elem_nonl) then
        call getvid('COMPORTEMENT', 'CARCRI', iocc=1, nbret=nocc)
        if (nocc > 0) then
            call getvid('COMPORTEMENT', 'CARCRI', iocc=1, nbret=nocc, &
                        scal=ds_constitutive%carcri)
            ASSERT(nocc == 1)
            call getvid('COMPORTEMENT', 'COMPOR', iocc=1, nbret=nocc, &
                        scal=ds_constitutive%compor)
            ASSERT(nocc == 1)
            call getvid('COMPORTEMENT', 'MULT_COMP', iocc=1, nbret=nocc, &
                        scal=ds_constitutive%mult_comp)
            ASSERT(nocc == 1)
            verbose = ASTER_FALSE
        else
            call nmdorc(model, materField, l_etat_init, &
                        ds_constitutive%compor, ds_constitutive%carcri, ds_constitutive%mult_comp)
            verbose = ASTER_TRUE
        end if
        call nonlinDSConstitutiveInit(model, caraElem, ds_constitutive, verbose)
    end if
!
end subroutine
