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
subroutine dlarch(result, neq, istoc, iarchi, texte, &
                  alarm, temps, nbtyar, typear, masse, &
                  depl, vite, acce, fexte, famor, &
                  fliai)
!
!  ARCHIVAGE DANS L'OBJET CHAMNO DU CHAMP DE DEPLACEMENT,DE VITESSE
!  ET/OU D'ACCELERATION ISSU D'UN CALCUL TRANSITOIRE DIRECT
!
! ---------------------------------------------------------------------
!  IN  : NEQ       : NOMBRE D'EQUATIONS
!  IN  : ISTOC     : PILOTAGE DU STOCKAGE DES RESULTATS
!  IN  : IARCHI    : PILOTAGE DE L'ARCHIVAGE DES RESULTATS
!  IN  : TEXTE     : COMMENTAIRE A IMPRIMER
!  IN  : ALARM     : EMISSION D'ALARME SI >0
!  IN  : TEMPS     : INSTANT D'ARCHIVAGE
!  IN  : NBTYAR    : TAILLE DE TYPEAR
!  IN  : TYPEAR    : TABLEAU INDIQUANT SI ON ARCHIVE LES DIFFERENTS
!                    CHAMPS (DEPL, VIT ET ACC) (NBTYAR)
!  IN  : MASSE     : NOM DE LA MATRICE DE MASSE
!  IN  : DEPL      : TABLEAU DES DEPLACEMENTS
!  IN  : VITE      : TABLEAU DES VITESSES
!  IN  : ACCE      : TABLEAU DES ACCELERATIONS
!
!     ------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
! DECLARATION PARAMETRES D'APPELS
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jelibe.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrem.h"
    integer(kind=8) :: neq, istoc, iarchi, alarm
    integer(kind=8) :: nbtyar
!
    real(kind=8) :: depl(neq), vite(neq), acce(neq)
    real(kind=8) :: fexte(neq), famor(neq), fliai(neq)
    real(kind=8) :: temps
!
    character(len=8) :: masse
    character(len=8) :: result
    character(len=16) :: typear(nbtyar)
    character(len=*) :: texte
!
!
!
    integer(kind=8) :: iaux, jaux, itype
    integer(kind=8) :: lgcomm
    character(len=8) :: k8b
!
    character(len=24) :: chamno
!
!====
! 1. PREALABLES
!====
!
!
! 1.2. ==> INSTANT
!
    if (istoc .eq. 0) then
        iarchi = iarchi+1
        call rsadpa(result, 'E', 1, 'INST', iarchi, &
                    0, sjv=iaux, styp=k8b)
        zr(iaux) = temps
    else
        call rsadpa(result, 'L', 1, 'INST', iarchi, &
                    0, sjv=iaux, styp=k8b)
        temps = zr(iaux)
    end if
!
! 1.3. ==> COMMENTAIRE
!
    lgcomm = lxlgut(texte)
!
!====
! 2. STOCKAGE DES CHAMPS DESIGNES
!====
!
    do 21, itype = 1, nbtyar
!
        if (typear(itype) .ne. '    ') then
!
            call rsexch(' ', result, typear(itype), iarchi, chamno, &
                        iaux)
            if (iaux .eq. 0) then
                if (alarm .gt. 0) then
                    call utmess('A', 'ALGORITH2_64', sk=chamno)
                end if
                goto 21
            else if (iaux .eq. 100) then
                call vtcrem(chamno, masse, 'G', 'R')
            else
                ASSERT(.false.)
            end if
!
            chamno(20:24) = '.VALE'
            call jeveuo(chamno, 'E', jaux)
!
            if (typear(itype) .eq. 'DEPL') then
                do 211, iaux = 1, neq
                    zr(jaux+iaux-1) = depl(iaux)
211                 continue
                    else if (typear(itype) .eq. 'VITE') then
                    do 212, iaux = 1, neq
                        zr(jaux+iaux-1) = vite(iaux)
212                     continue
                        else if (typear(itype) .eq. 'ACCE') then
                        do 213, iaux = 1, neq
                            zr(jaux+iaux-1) = acce(iaux)
213                         continue
                            else if (typear(itype) .eq. 'FORC_EXTE') then
                            do 214, iaux = 1, neq
                                zr(jaux+iaux-1) = fexte(iaux)
214                             continue
                                else if (typear(itype) .eq. 'FORC_AMOR') then
                                do 215, iaux = 1, neq
                                    zr(jaux+iaux-1) = famor(iaux)
215                                 continue
                                    else if (typear(itype) .eq. 'FORC_LIAI') then
                                    do 216, iaux = 1, neq
                                        zr(jaux+iaux-1) = fliai(iaux)
216                                     continue
                                        end if
!
                                        call jelibe(chamno)
                                        call rsnoch(result, typear(itype), iarchi)
!
                                        end if
!
21                                      continue
!
                                        istoc = 1
!
!
                                        end subroutine
