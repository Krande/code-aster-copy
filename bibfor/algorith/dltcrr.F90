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

subroutine dltcrr(result, neq, nbordr, iarchi, texte, &
                  t0, lcrea, typres, masse, rigid, &
                  amort, dep0, vit0, acc0, fexte, &
                  famor, fliai, numedd, nume, nbtyar, &
                  typear)
!
!
!       DYNAMIQUE LINEAIRE TRANSITOIRE - CREATION DES RESULTATS
!       -         -        -             --           -
!
! ----------------------------------------------------------------------
!  IN  : NEQ       : NOMBRE D'EQUATIONS
!  IN  : IARCHI    : PILOTAGE DE L'ARCHIVAGE DES RESULTATS
!  IN  : TEXTE     : COMMENTAIRE A IMPRIMER
!  IN  : T0        : INSTANT DE CALCUL INITIAL
!  IN  : LCREA     : LOGIQUE INDIQUANT SI IL Y A REPRISE
!  IN  : TYPRES    : TYPE DE RESULTAT
!  IN  : MASSE     : MATRICE DE MASSE
!  IN  : RIGID     : MATRICE DE RIGIDITE
!  IN  : AMORT     : MATRICE D'AMORTISSEMENT
!  VAR : DEP0      : TABLEAU DES DEPLACEMENTS A L'INSTANT N
!  VAR : VIT0      : TABLEAU DES VITESSES A L'INSTANT N
!  VAR : ACC0      : TABLEAU DES ACCELERATIONS A L'INSTANT N
!  IN  : NUMEDD    : NUME_DDL DE LA MATR_ASSE RIGID
!  IN  : NUME      : NUMERO D'ORDRE DE REPRISE
!
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dlarch.h"
#include "asterfort/refdaj.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rscrsd.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"

    integer(kind=8) :: neq, nbordr, iarchi, ir
    integer(kind=8) :: nume, nbtyar
!
    real(kind=8) :: dep0(neq), vit0(neq), acc0(neq), t0
    real(kind=8) :: fexte(2*neq), famor(2*neq), fliai(2*neq)
!
    character(len=8) :: masse, rigid, amort
    character(len=8) :: result
    character(len=16) :: typres
    character(len=16) :: typear(nbtyar)
    character(len=24) :: numedd, matric(3)
    character(len=*) :: texte
!
    aster_logical :: lcrea
!
!
!
    integer(kind=8) :: istoc
!
!
!
!====
! 2. CREATION DE LA STRUCTURE DE DONNEE RESULTAT
!====
!
    if (lcrea) then
!
! 2.1. ==> CREATION DE LA STRUCTURE DE DONNEE RESULTAT
!
        call rscrsd('G', result, typres, nbordr)
        matric(1) = rigid
        matric(2) = masse
        matric(3) = amort
        call refdaj('F', result, nbordr, numedd, 'DYNAMIQUE', &
                    matric, ir)
!
! 2.2. ==> ARCHIVAGE INITIAL
!
        iarchi = -1
        istoc = 0
!
        call dlarch(result, neq, istoc, iarchi, texte, &
                    1, t0, nbtyar, typear, masse, &
                    dep0, vit0, acc0, fexte(neq+1), famor(neq+1), &
                    fliai(neq+1))
!
        call utmess('I', 'PROGRESS_1', ni=2, vali=[0, 0], &
                    nr=2, valr=[t0, t0])

!
        iarchi = 0
!
!====
! 3. RECUPERATION
!====
    else
        nbordr = nbordr+nume
        call rsagsd(result, nbordr)
    end if
!
!====
! 4. TITRE
!====
!
    call titre()
!
end subroutine
