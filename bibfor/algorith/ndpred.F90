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
subroutine ndpred(sddyna, valinc, solalg)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ndynin.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/vtaxpy.h"
#include "asterfort/vtzero.h"
!
    character(len=19) :: sddyna
    character(len=19) :: solalg(*), valinc(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (DYNAMIQUE)
!
! CALCUL DES PREDICTEURS
!
! ----------------------------------------------------------------------
!
!
! IN  SDDYNA : SD DYNAMIQUE
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
!
!
!
!
    aster_logical :: ldepl, lvite, lacce
    integer(kind=8) :: n
    real(kind=8) :: coefd(3), coefv(3), coefa(3)
    character(len=24) :: vect(3)
    character(len=19) :: depdel, vitdel, accdel
    character(len=19) :: depkm1, vitkm1, acckm1, romkm1, romk
    character(len=19) :: depplu, vitplu, accplu
    character(len=19) :: depmoi, vitmoi, accmoi
    integer(kind=8) :: ifm, niv
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... CALCUL PREDICTEURS'
    end if
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(valinc, 'VALINC', 'DEPKM1', depkm1)
    call nmchex(valinc, 'VALINC', 'VITKM1', vitkm1)
    call nmchex(valinc, 'VALINC', 'ACCKM1', acckm1)
    call nmchex(valinc, 'VALINC', 'ROMKM1', romkm1)
    call nmchex(valinc, 'VALINC', 'ROMK  ', romk)
    call nmchex(solalg, 'SOLALG', 'DEPDEL', depdel)
    call nmchex(solalg, 'SOLALG', 'VITDEL', vitdel)
    call nmchex(solalg, 'SOLALG', 'ACCDEL', accdel)
!
! --- TYPE DE FORMULATION SCHEMA DYNAMIQUE GENERAL
!
    ldepl = ndynin(sddyna, 'FORMUL_DYNAMIQUE') .eq. 1
    lvite = ndynin(sddyna, 'FORMUL_DYNAMIQUE') .eq. 2
    lacce = ndynin(sddyna, 'FORMUL_DYNAMIQUE') .eq. 3
!
! --- COEFFICIENTS POUR PREDICTEURS
!
    coefd(1) = ndynre(sddyna, 'COEF_DEPL_DEPL')
    coefd(2) = ndynre(sddyna, 'COEF_DEPL_VITE')
    coefd(3) = ndynre(sddyna, 'COEF_DEPL_ACCE')
    coefv(1) = ndynre(sddyna, 'COEF_VITE_DEPL')
    coefv(2) = ndynre(sddyna, 'COEF_VITE_VITE')
    coefv(3) = ndynre(sddyna, 'COEF_VITE_ACCE')
    coefa(1) = ndynre(sddyna, 'COEF_ACCE_DEPL')
    coefa(2) = ndynre(sddyna, 'COEF_ACCE_VITE')
    coefa(3) = ndynre(sddyna, 'COEF_ACCE_ACCE')
!
! --- MISE A JOUR CHAMPS GRANDES ROTATIONS
!
    call vtzero(romkm1)
    call vtzero(romk)
!
! --- CALCUL DES PREDICTEURS
!
    call nmchex(valinc, 'VALINC', 'DEPPLU', depplu)
    call nmchex(valinc, 'VALINC', 'VITPLU', vitplu)
    call nmchex(valinc, 'VALINC', 'ACCPLU', accplu)
    call nmchex(valinc, 'VALINC', 'DEPMOI', depmoi)
    call nmchex(valinc, 'VALINC', 'VITMOI', vitmoi)
    call nmchex(valinc, 'VALINC', 'ACCMOI', accmoi)
    call copisd('CHAMP_GD', 'V', depmoi, depkm1)
    call copisd('CHAMP_GD', 'V', vitmoi, vitkm1)
    call copisd('CHAMP_GD', 'V', accmoi, acckm1)
    vect(1) = depkm1
    vect(2) = vitkm1
    vect(3) = acckm1
    if (ldepl) then
        call vtzero(vitplu)
        call vtzero(accplu)
        do n = 1, 3
            call vtaxpy(coefv(n), vect(n), vitplu)
            call vtaxpy(coefa(n), vect(n), accplu)
        end do
    else if (lvite) then
        call vtzero(depplu)
        call vtzero(accplu)
        call copisd('CHAMP_GD', 'V', vitkm1, vitplu)
        do n = 1, 3
            call vtaxpy(coefd(n), vect(n), depplu)
            call vtaxpy(coefa(n), vect(n), accplu)
        end do
    else if (lacce) then
        call vtzero(vitplu)
        call vtzero(depplu)
        do n = 1, 3
            call vtaxpy(coefv(n), vect(n), vitplu)
            call vtaxpy(coefd(n), vect(n), depplu)
        end do
        call copisd('CHAMP_GD', 'V', acckm1, accplu)
    else
        ASSERT(.false.)
    end if
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ...... PRED. DEPL. '
        call nmdebg('VECT', depplu, ifm)
        write (ifm, *) '<MECANONLINE> ...... PRED. VITE. '
        call nmdebg('VECT', vitplu, ifm)
        write (ifm, *) '<MECANONLINE> ...... PRED. ACCE. '
        call nmdebg('VECT', accplu, ifm)
    end if
!
! --- INITIALISATION DE L'INCREMENT DE DEPLACEMENT DEPDEL
!
    call vtzero(depdel)
    if (lacce) then
        do n = 1, 3
            call vtaxpy(coefd(n), vect(n), depdel)
        end do
        call vtaxpy(-1.d0, depkm1, depdel)
    end if
!
! --- INITIALISATION DE L'INCREMENT DE VITESSE/ACCELERATION
!
    call vtzero(vitdel)
    call vtzero(accdel)
!
    call jedema()
!
end subroutine
