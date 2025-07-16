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
subroutine chckVari(comporPrevZ, comporCurrZ, variZ, ligrelZ)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/adaptVari.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcomp_chck_cmp.h"
#include "asterfort/vrcomp_chck_rela.h"
#include "asterfort/vrcomp_prep.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: comporPrevZ, comporCurrZ, variZ, ligrelZ
!
! --------------------------------------------------------------------------------------------------
!
! Check consistency of VARI_ELGA field and modify it if possible
!
! --------------------------------------------------------------------------------------------------
!
! In  vari         : internal variable field
! In  comporCurr   : map to describe previous behaviours
! In  ligrel       : finite element descriptor
! In  comporPrev   : map to describe current behaviours
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: typeStop = 'I'
    aster_logical, parameter :: verbose = ASTER_TRUE
    character(len=19) :: variRedu, vari
    character(len=19) :: comporCurr, comporPrev
    character(len=19) :: comporCurrRedu, comporPrevRedu
! - Beahviours can been mixed with each other
    character(len=48), parameter :: comp_comb_1 = &
                                    'LEMAITRE        VMIS_ISOT_LINE  VMIS_ISOT_TRAC'
! - Beahviours can been mixed with all other ones
    character(len=48), parameter :: comp_comb_2 = &
                                    'ELAS            SANS            KIT_CG'
    aster_logical :: newBehaviourOnCell, nbSpgDifferent, lModiVari, error
    aster_logical :: inconsistentBehaviour, nbVariDifferent, hasPreviousBehaviour
    character(len=19) :: ligrel, ligrelPrev
    character(len=8) :: meshCompor, meshField, mesh
    integer(kind=8) :: nbCell
    character(len=8), pointer :: cesk(:) => null()
    integer(kind=8), pointer :: cesd(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    comporPrev = comporPrevZ
    comporCurr = comporCurrZ
    ligrel = ligrelZ
    vari = variZ
    hasPreviousBehaviour = comporPrev .ne. " "
    error = ASTER_FALSE

! - Acces to reduced CARTE DCEL_I (see CESVAR) on current comportement
    call jeveuo(comporCurr//'.CESD', 'L', vi=cesd)
    call jeveuo(comporCurr//'.CESK', 'L', vk8=cesk)
    meshCompor = cesk(1)
    nbCell = cesd(1)

! - Get LIGREL from VARI_ELGA field
    call dismoi('NOM_LIGREL', vari, 'CHAM_ELEM', repk=ligrelPrev)

! - Check meshes
    call dismoi('NOM_MAILLA', vari, 'CHAMP', repk=meshField)
    if (meshCompor .ne. meshField) then
        call utmess('F', 'COMPOR6_1')
    end if
    mesh = meshCompor

! - Prepare fields
    if (hasPreviousBehaviour) then
        call vrcomp_prep(vari, variRedu, &
                         comporCurr, comporCurrRedu, &
                         comporPrev, comporPrevRedu)
    else
        call vrcomp_prep(vari, variRedu, &
                         comporCurr, comporCurrRedu)
        comporPrevRedu = " "
    end if

! - Check if comportments are the same (or compatible)
    newBehaviourOnCell = ASTER_FALSE
    inconsistentBehaviour = ASTER_FALSE
    if (hasPreviousBehaviour) then
        call vrcomp_chck_rela(mesh, nbCell, &
                              comporCurrRedu, comporPrevRedu, &
                              ligrel, ligrelPrev, &
                              comp_comb_1, comp_comb_2, verbose, &
                              newBehaviourOnCell, inconsistentBehaviour, &
                              lModiVari)
    end if
    if (newBehaviourOnCell) then
        error = ASTER_TRUE
        call utmess(typeStop, 'COMPOR6_2')
    end if
    if (inconsistentBehaviour) then
        error = ASTER_TRUE
        call utmess(typeStop, 'COMPOR6_3')
    end if

! - Check if elements have the same number of internal variables and Gauss-subpoints
    nbSpgDifferent = ASTER_FALSE
    nbVariDifferent = ASTER_FALSE
    call vrcomp_chck_cmp(mesh, nbCell, &
                         comporCurr, comporCurrRedu, comporPrevRedu, &
                         variRedu, comp_comb_2, &
                         ligrel, ligrelPrev, &
                         verbose, &
                         nbSpgDifferent, nbVariDifferent, lModiVari)
    if (nbSpgDifferent) then
        error = ASTER_TRUE
        call utmess(typeStop, 'COMPOR6_4')
    end if
    if (nbVariDifferent) then
        error = ASTER_TRUE
        call utmess(typeStop, 'COMPOR6_5')
    end if

! - Change internal variable field
    if (error) then
        call utmess('F', 'MECANONLINE5_2')
    end if
    if (lModiVari) then
        call adaptVari(comporCurr, vari, ligrel)
    end if

! - Clean
    call detrsd('CHAM_ELEM_S', variRedu)
    if (hasPreviousBehaviour) then
        call detrsd('CHAM_ELEM_S', comporPrevRedu)
    end if
    call detrsd('CHAM_ELEM_S', comporCurrRedu)
!
    call jedema()
!
end subroutine
