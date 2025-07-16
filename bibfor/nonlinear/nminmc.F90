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
subroutine nminmc(listFuncActi, &
                  model, listLoad, numfix, &
                  meelem, measse)
!
    use NonLin_Datastructure_type
    use NonLinearElem_module, only: elemSuper, asseSuper, elemDiri
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/utmess.h"
#include "asterfort/nmchex.h"
!
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=24), intent(in) :: model
    character(len=19), intent(in) :: listLoad
    character(len=24), intent(in) :: numfix
    character(len=19), intent(in) :: meelem(*), measse(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL)
!
! PRE-CALCUL DES MATRICES ELEMENTAIRES CONSTANTES AU COURS DU CALCUL
!
! ----------------------------------------------------------------------
!
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IN  SDDYNA : SD DYNAMIQUE
! IN  MODELE : NOM DU MODELE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! IN  LISCHA : LISTE DES CHARGEMENTS
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR DEPLACEMENTS
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IO  ds_measure       : datastructure for measure and statistics management
! OUT MEELEM : MATRICES ELEMENTAIRES
! In  ds_system        : datastructure for non-linear system management
! OUT MEASSE : MATRICES ASSEMBLEES
!
! ----------------------------------------------------------------------
!
    aster_logical :: lSuperElement
    integer(kind=8) :: ifm, niv
    character(len=24) :: diriElem
    character(len=24) :: superElem, superAsse
!
! ----------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_18')
    end if

! - Active functionnalities
    lSuperElement = isfonc(listFuncActi, 'MACR_ELEM_STAT')

! - Compute elementary matrices for Dirichet (B matrix for Lagrange multipliers)
    call nmchex(meelem, 'MEELEM', 'MEDIRI', diriElem)
    call elemDiri(model, listLoad, diriElem)

! - Compute elementary matrices for super-elements and assemble them
    if (lSuperElement) then
        call nmchex(meelem, 'MEELEM', 'MESSTR', superElem)
        call nmchex(measse, 'MEASSE', 'MESSTR', superAsse)
        call elemSuper(model, superElem)
        call asseSuper(numfix, listLoad, superElem, superAsse)
    end if
!
end subroutine
