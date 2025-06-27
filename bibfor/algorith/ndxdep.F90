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
subroutine ndxdep(numeDof, listFuncActi, numeTime, sddisc, &
                  sddyna, nlDynaDamping, &
                  sdnume, valinc, solalg, veasse)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nmfext.h"
#include "asterfort/nmmajc.h"
#include "asterfort/nmsolu.h"
!
    integer(kind=8), intent(in) :: listFuncActi(*), numeTime
    character(len=19), intent(in) :: sddisc, sdnume
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=24), intent(in) :: numeDof
    character(len=19), intent(in) :: veasse(*), solalg(*), valinc(*)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! CALCUL DE L'INCREMENT DE DEPLACEMENT A PARTIR DE(S) DIRECTION(S)
! DE DESCENTE
!
! CAS EXPLICITE
!
! --------------------------------------------------------------------------------------------------
!
! IN  NUMEDD : NUME_DDL
! In  listFuncActi     : list of active functionnalities
! IN  NUMINS : NUMERO D'INSTANT
! IN  SDNUME : SD NUMEROTATION
! IN  SDDISC : SD DISCRETISATION
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: timePrev, timeCurr, timeIncr
    character(len=19) :: cnfext
    character(len=19) :: ddepla, deppr1
    character(len=19) :: dvitla, vitpr1
    character(len=19) :: daccla, accpr1
    integer(kind=8) :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> CORRECTION INCR. DEPL.'
    end if

! - Get hat-variables
    call nmchex(solalg, 'SOLALG', 'DDEPLA', ddepla)
    call nmchex(solalg, 'SOLALG', 'DEPPR1', deppr1)
    call nmchex(solalg, 'SOLALG', 'DVITLA', dvitla)
    call nmchex(solalg, 'SOLALG', 'VITPR1', vitpr1)
    call nmchex(solalg, 'SOLALG', 'DACCLA', daccla)
    call nmchex(solalg, 'SOLALG', 'ACCPR1', accpr1)

! - Time
    timePrev = diinst(sddisc, numeTime-1)
    timeCurr = diinst(sddisc, numeTime)
    timeIncr = timeCurr-timePrev

! - Compute external forces
    call nmchex(veasse, 'VEASSE', 'CNFEXT', cnfext)
    call nmfext(0.d0, listFuncActi, veasse, cnfext, &
                sddyna_=sddyna, nlDynaDamping_=nlDynaDamping)
!
! --- CONVERSION RESULTAT dU VENANT DE K.dU = F SUIVANT SCHEMAS
!
    call nmsolu(sddyna, solalg)
!
! --- AJUSTEMENT DE LA DIRECTION DE DESCENTE (AVEC ETA, RHO ET OFFSET)
!
    call copisd('CHAMP_GD', 'V', deppr1, ddepla)
    call copisd('CHAMP_GD', 'V', vitpr1, dvitla)
    call copisd('CHAMP_GD', 'V', accpr1, daccla)
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... DEPL. PRED. (1) : '
        call nmdebg('VECT', deppr1, ifm)
        write (ifm, *) '<MECANONLINE> ... DEPL. SOLU.     : '
        call nmdebg('VECT', ddepla, ifm)
        write (ifm, *) '<MECANONLINE> ... VITE. PRED. (1) : '
        call nmdebg('VECT', vitpr1, ifm)
        write (ifm, *) '<MECANONLINE> ... VITE. SOLU.     : '
        call nmdebg('VECT', dvitla, ifm)
        write (ifm, *) '<MECANONLINE> ... ACCE. PRED. (1) : '
        call nmdebg('VECT', accpr1, ifm)
        write (ifm, *) '<MECANONLINE> ... ACCE. SOLU.     : '
        call nmdebg('VECT', daccla, ifm)
    end if
!
! --- ACTUALISATION DES CHAMPS SOLUTIONS
!
    call nmmajc(listFuncActi, sddyna, sdnume, timeIncr, numeDof, &
                valinc, solalg)
!
    call jedema()
end subroutine
