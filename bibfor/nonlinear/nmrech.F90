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

subroutine nmrech(fm, f, fopt, fcvg, rhomin, &
                  rhomax, rhoexm, rhoexp, rhom, rho, &
                  rhoopt, ldcopt, ldccvg, opt, act, &
                  stite)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/zbinte.h"
    real(kind=8) :: rhomin, rhomax, rhoexm, rhoexp
    real(kind=8) :: rhom, rho, rhoopt
    real(kind=8) :: fm, f, fopt, fcvg
    aster_logical :: stite
    integer(kind=8) :: ldcopt, ldccvg
    integer(kind=8) :: opt, act
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (RECH. LINE. - METHODE CORDE)
!
! RECHERCHE LINEAIRE AVEC LA METHODE CORDE
!
! ----------------------------------------------------------------------
!
!
! I/O FM     : VALEUR PRECEDENTE DE LA FONCTIONNELLE
! IN  F      : VALEUR COURANTE DE LA FONCTIONNELLE
! I/O FOPT   : VALEUR OPTIMALE DE LA FONCTIONNELLE
! IN  FCVG   : VALEUR DONNANT LA VALEUR DE LA FONCTIONNELLE POUR QUE
!              L'ALGO CONVERGE
! I/O RHOM   : VALEUR PRECEDENTE DU COEF. RECH. LINE.
! IN  RHO    : VALEUR COURANTE DU COEF. RECH. LINE.
! I/O RHOOPT : VALEUR OPTIMALE DU COEF. RECH. LINE.
! IN  LDCCVG : CODE RETOUR INTEGRATION COMPORTEMENT
! OUT LDCOPT : CODE RETOUR INTEGRATION COMPORTEMENT QUAND COEF. RECH.
!              LINE. OPTIMAL
! I/O ACT    : INDICE DE LA SOLUTION (DEUX QUAND PILOTAGE)
! OUT OPT    : INDICE DE LA SOLUTION QUAND COEF. RECH.
!              LINE. OPTIMAL
! OUT STITE  : .TRUE. SI ALGO. A CONVERGE
! IN  RHOMIN : BORNE INFERIEURE DE RECHERCHE
! IN  RHOMAX : BORNE SUPERIEURE DE RECHERCHE
! IN  RHOEXM : INTERVALLE [RHOEXM,RHOEXP] POUR EXCLUSION
! IN  RHOEXP : INTERVALLE [RHOEXM,RHOEXP] POUR EXCLUSION
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: rhotmp
!
! ----------------------------------------------------------------------
!
    stite = .false.
!
! --- PRISE EN COMPTE D'UN RESIDU OPTIMAL SI NECESSAIRE
!
    if (abs(f) .lt. fopt) then
        rhoopt = rho
        ldcopt = ldccvg
        fopt = abs(f)
        opt = act
        act = 3-act
        if (abs(f) .lt. fcvg) then
            stite = .true.
            goto 100
        end if
    end if
!
! --- CALCUL DE RHO(N+1) PAR METHODE DE SECANTE AVEC BORNES
!
    rhotmp = rho
    if (abs(f-fm) .gt. r8prem()) then
        rho = (f*rhom-fm*rho)/(f-fm)
        call zbinte(rho, rhomin, rhomax, rhoexm, rhoexp)
    else if (f*(rho-rhom)*(f-fm) .le. 0.d0) then
        rho = rhomax
    else
        rho = rhomin
    end if
    rhom = rhotmp
    fm = f
!
100 continue
!
end subroutine
