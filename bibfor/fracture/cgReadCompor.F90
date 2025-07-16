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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine cgReadCompor(result_in, compor, iord0, l_incr)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_init.h"
#include "asterfort/comp_meca_elas.h"
#include "asterfort/cgCreateCompIncr.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in)  :: result_in
    character(len=19), intent(inout) :: compor
    integer(kind=8), intent(in) :: iord0
    aster_logical, intent(out) :: l_incr
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Read or Create CARTE Comportement
!
! In result_in : name of SD RESULTAT
! In compor    : name of <CARTE> COMPORTEMENT
! In iord0     : which NUME_ORDRE
! Out l_incr   : Incremental behavior or not ?
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, nbetin, ier
    integer(kind=8) :: nb_vale, nb_zone, nb_cmp_max, i_zone
    character(len=8) :: mesh, repk, model, discretization
    character(len=16) :: defo_comp
    aster_logical :: l_etat_init, l_impel, l_disc
    integer(kind=8), pointer :: v_compor_desc(:) => null()
    character(len=16), pointer :: v_compor_vale(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    l_etat_init = ASTER_FALSE
    l_impel = ASTER_FALSE
!
    call dismoi('MODELE', result_in, 'RESULTAT', repk=model)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
!
    call getfac('ETAT_INIT', nbetin)
    if (nbetin .ne. 0) l_etat_init = ASTER_TRUE
!
! --- Read COMPOR <CARTE> in RESULT
!
    call rsexch(' ', result_in, 'COMPORTEMENT', iord0, compor, iret)
!
! --- No COMPOR <CARTE> in RESULT: create ELAS COMPOR <CARTE>
!
    if (iret .ne. 0) then
        l_impel = ASTER_TRUE
        call comp_init(mesh, compor, 'V')
        call comp_meca_elas(compor, l_etat_init)
    end if
!
! --- Check COMPORTEMENT / RELATION in result for incremental comportement
!
    if (.not. l_impel) then
        call cgCreateCompIncr(compor, l_etat_init)
    end if
!
! --- Incremental behavior or not ?
!
    call dismoi('ELAS_INCR', compor, 'CARTE_COMPOR', repk=repk)
    if (repk == 'ELAS') then
        l_incr = ASTER_FALSE
    else if (repk == 'INCR' .or. repk == 'MIXTE') then
        l_incr = ASTER_TRUE
    else
        ASSERT(ASTER_FALSE)
    end if
!
! --- Only PETIT is accepted
!
    call jeveuo(compor//'.DESC', 'L', vi=v_compor_desc)
    call jeveuo(compor//'.VALE', 'L', vk16=v_compor_vale)
    call jelira(compor//'.VALE', 'LONMAX', nb_vale)
    nb_zone = v_compor_desc(3)
    nb_cmp_max = nb_vale/v_compor_desc(2)
!
    do i_zone = 1, nb_zone
        defo_comp = v_compor_vale(nb_cmp_max*(i_zone-1)+DEFO)
        if ((defo_comp .ne. "PETIT") .and. (defo_comp .ne. "GREEN_LAGRANGE")) then
            call utmess("F", "RUPTURE3_9", sk=defo_comp)
        end if
        if (defo_comp .eq. "GREEN_LAGRANGE") then
            call utmess("I", "RUPTURE3_12")

            call getvtx('THETA', 'DISCRETISATION', iocc=1, scal=discretization, nbret=ier)
            l_disc = (discretization == "LINEAIRE") .or. (discretization == "LEGENDRE")
            ASSERT(l_disc)

            if (discretization .eq. "LEGENDRE") then
                call utmess("F", "RUPTURE3_13")
            end if

        end if
    end do
!
end subroutine
