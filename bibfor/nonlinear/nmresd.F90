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
subroutine nmresd(fonact, sddyna, ds_measure, solveu, &
                  numedd, instan, maprec, matass, cndonn, &
                  cnpilo, cncine, solalg, rescvg, &
                  ds_algorom_)
!
    use NonLin_Datastructure_type
    use ROM_Datastructure_type
    use ldlt_xp_data_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmresg.h"
#include "asterfort/nmreso.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/utmess.h"
#include "asterfort/romAlgoNLSystemSolve.h"
#include "asterfort/romAlgoNLCorrEFMatrixModify.h"
#include "asterfort/romAlgoNLCorrEFResiduModify.h"
!
    integer(kind=8) :: fonact(*)
    character(len=19) :: solalg(*)
    character(len=19) :: maprec, matass
    type(ROM_DS_AlgoPara), optional, intent(in) :: ds_algorom_
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19) :: solveu, sddyna
    character(len=19) :: cncine, cndonn, cnpilo
    character(len=24) :: numedd
    real(kind=8) :: instan
    integer(kind=8) :: rescvg
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - CALCUL)
!
! RESOLUTION DU SYSTEME LINEAIRE K.dU = F
!
! --------------------------------------------------------------------------------------------------
!
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IN  SDDYNA : SD DYNAMIQUE
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_algorom       : datastructure for ROM parameters
! IN  SOLVEU : SOLVEUR
! IN  NUMEDD : NUME_DDL
! IN  INSTAN : INSTANT COURANT
! IN  MAPREC : MATRICE DE PRECONDITIONNEMENT (GCPC)
! IN  MATASS : MATRICE ASSEMBLEE
! IN  CNDONN : CHAM_NO DE CHARGE DONNEE
! IN  CNPILO : CHAM_NO DE CHARGE PILOTEE
! IN  CNCINE : CHAM_NO DE CHARGE CINEMATIQUE
! OUT DEPSOL : SOLUTION DU DU SYSTEME K.DU = DF
!                      DEPSOL(1) EN L'ABSENCE DE PILOTAGE
!                      DEPSOL(1) ET DEPSOL(2) AVEC PILOTAGE
! OUT RESCVG : CODE RETOUR RESOLUTION SYSTEME LINEAIRE
!                -1 : PAS DE RESOLUTION
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : NOMBRE MAXIMUM D'ITERATIONS ATTEINT
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lprmo, l_rom
    character(len=19) :: depso1, depso2
    character(len=24) :: mata24, vect24, cncine24
    integer(kind=8) :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_77')
    end if
!
! --- EXTRACTION VARIABLES CHAPEAUX
!
    call nmchex(solalg, 'SOLALG', 'DEPSO1', depso1)
    call nmchex(solalg, 'SOLALG', 'DEPSO2', depso2)
!
! --- FONCTIONNALITES ACTIVEES
!
    lprmo = ndynlo(sddyna, 'PROJ_MODAL')
    l_rom = isfonc(fonact, 'ROM')
!
! --- RESOLUTION GENERALISEE OU PHYSIQUE
!
    call nmtime(ds_measure, 'Init', 'Solve')
    call nmtime(ds_measure, 'Launch', 'Solve')
!
    if (lprmo) then
        call nmresg(numedd, sddyna, instan, cndonn, depso1)
    elseif (l_rom) then
        if (ds_algorom_%phase .eq. 'HROM') then
            rescvg = 0
            mata24 = matass
            vect24 = cndonn
            cncine24 = cncine
            call romAlgoNLSystemSolve(mata24, vect24, cncine24, ds_algorom_, depso1)
        elseif (ds_algorom_%phase .eq. 'CORR_EF') then
            vect24 = cndonn
            call romAlgoNLCorrEFMatrixModify(numedd, matass, ds_algorom_)
            call romAlgoNLCorrEFResiduModify(vect24, ds_algorom_)
            call nmreso(fonact, vect24, cnpilo, cncine, solveu, &
                        maprec, matass, depso1, depso2, rescvg)
        else
            ASSERT(.false.)
        end if
    else
        call nmreso(fonact, cndonn, cnpilo, cncine, solveu, &
                    maprec, matass, depso1, depso2, rescvg)
    end if
!
    call nmtime(ds_measure, 'Stop', 'Solve')
    call nmrinc(ds_measure, 'Solve')
!   la stratégie ap2foi a causé une factorisation : on incrémente
    if (ap2foi_called) then
        call nmrinc(ds_measure, 'Factor')
        ap2foi_called = ASTER_FALSE
    end if

!
end subroutine
