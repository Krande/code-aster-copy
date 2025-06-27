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
! person_in_charge: samuel.geniaut at edf.fr
!
subroutine cgleco(resu, modele, mate, iord0, compor, &
                  incr)
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/cgvein.h"

#include "asterfort/comp_init.h"
#include "asterfort/comp_info.h"
#include "asterfort/comp_meca_elas.h"
#include "asterfort/dismoi.h"
#include "asterfort/gverlc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmdocc.h"
#include "asterfort/rsexch.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: iord0
    character(len=8), intent(in) :: resu
    character(len=8), intent(in) :: modele
    character(len=8), intent(in) :: mate
    character(len=19), intent(out) :: compor
    aster_logical, intent(out) :: incr
!
! --------------------------------------------------------------------------------------------------
!
! CALC_G
!
! Comportment selection
!
! --------------------------------------------------------------------------------------------------
!
! In  resu   : name of result
! In  model  : name of model
! In  mate   : name of material field
! In  iord0  : first NUME_ORDRE in result
! Out compor : name of COMPOR <CARTE>
! Out incr   : if incrental comportment
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbFactorKeyword, iret, nbetin
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    character(len=24) :: repk
    character(len=8) :: mesh
    aster_logical :: limpel, l_etat_init
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)

! - Initializations
    limpel = .false.
    incr = .false.
    nbFactorKeyword = 0

    compor = '&&CGLECO.COMPOR'
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
    call getfac('ETAT_INIT', nbetin)
    l_etat_init = nbetin .ne. 0
!
! - How many COMPORTEMENT in CALC_G ?
!
    call getfac(factorKeyword, nbFactorKeyword)
!
! - Get or create COMPOR <CARTE>
!
    if (nbFactorKeyword .eq. 0) then
!
! ----- No COMPORTEMENT: get from RESULT
!
        call rsexch(' ', resu, 'COMPORTEMENT', iord0, compor, iret)
!
! ----- No COMPOR <CARTE> in RESULT: create ELAS COMPOR <CARTE>
!
        if (iret .ne. 0) then
            limpel = .true.
            call comp_init(mesh, compor, 'V')
            call comp_meca_elas(compor, l_etat_init)
        end if
    else
! ----- Get COMPORTEMENT from command file
        call nmdocc(modele, mate, l_etat_init, compor, 'V')
        if (niv .ge. 2) then
            call comp_info(modele, compor)
        end if
    end if
!
! - Incremental comportement or not ?
!
    if (limpel) then
        incr = .false.
    else
        call dismoi('ELAS_INCR', compor, 'CARTE_COMPOR', repk=repk)
        if (repk .eq. 'ELAS') then
            incr = .false.
        else if (repk .eq. 'INCR' .or. repk .eq. 'MIXTE') then
            incr = .true.
        else
            ASSERT(.false.)
        end if
    end if
!
!  -Si Comportement dans CALC_G(alors RELATION est renseigne) ---> emission d'un message d'alarme !
!    normalement le comportement est recupere dans Meca_stat ou stat_non_line.
!
    if (nbFactorKeyword .gt. 0) then
        call utmess('A', 'RUPTURE1_41')
    end if
!
! - Check is CALG_G COMPOR <CARTE> is coherent with result COMPOR <CARTE>
!
    call gverlc(resu, compor, iord0)
!
! - Check COMPORTEMENT / RELATION in result for incremental comportement
!
    if (incr) then
        call cgvein(compor)
    end if
!
    call jedema()
!
end subroutine
