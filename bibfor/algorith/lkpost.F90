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

subroutine lkpost(imate, sigf, nvi, vip)
!
    implicit none
#include "asterfort/cos3t.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lkcrit.h"
#include "asterfort/rcvala.h"
#include "asterfort/get_varc.h"
    integer(kind=8) :: imate, nvi
    real(kind=8) :: sigf(6), vip(nvi)
! =================================================================
! IN  : IMATE  : ADRESSE DU MATERIAU CODE -------------------------
! --- : NVI    : NOMBRE DE VARIABLES INTERNES ---------------------
! OUT : VIP    : MISE A JOUR DES VARIABLES INTERNES DE POST -------
! =================================================================
    integer(kind=8) :: dimpar
    parameter(dimpar=12)
    integer(kind=8) :: cerr(dimpar)
    real(kind=8) :: mater(dimpar), i1, sii, devsig(6), lgleps, rcos3t
    real(kind=8) :: crit0, crite, tempd, tempf, tref
    parameter(lgleps=1.0d-8)
    character(len=16) :: nomc(dimpar)
    integer(kind=8) :: ndi, ndt
    common/tdim/ndt, ndi
!
! - Get temperatures
!
    call get_varc('RIGI', 1, 1, 'T', &
                  tempd, tempf, tref)
!
! =================================================================
! --- RECUPERATION DES PROPRIETES MATERIAUX -----------------------
! =================================================================
    nomc(1) = 'XI_PIC   '
    nomc(2) = 'XI_E     '
    nomc(3) = 'XI_ULT   '
    nomc(4) = 'A_0      '
    nomc(5) = 'M_0      '
    nomc(6) = 'S_0      '
    nomc(7) = 'XIV_MAX  '
    nomc(8) = 'MV_MAX   '
    nomc(9) = 'SIGMA_C  '
    nomc(10) = 'GAMMA_CJS'
    nomc(11) = 'PA       '
    nomc(12) = 'H0_EXT   '
!
    call rcvala(imate, ' ', 'LETK', 1, 'TEMP', &
                [tempd], dimpar, nomc(1), mater(1), cerr(1), 0)
!
! =================================================================
! --- DEFINITION DU NIVEAU DE DEGRADATION DE LA ROCHE SUIVANT LE DOMAINE
! --- DOMAINE = 0 : LE COMPORTEMENT RESTE ELASTIQUE
! --- DOMAINE = 1 : LA ROCHE EST FISSUREE PRE-PIC
! --- DOMAINE = 2 : LA ROCHE EST FISSUREE POST-PIC
! --- DOMAINE = 3 : LA ROCHE EST FRACTUREE
! --- DOMAINE = 4 : LA ROCHE EST DANS SON ETAT RESIDUEL
! =================================================================
    if (vip(1) .eq. 0.0d0) then
        vip(8) = 0
    else if (vip(1) .lt. mater(1)) then
        vip(8) = 1
    else if (vip(1) .lt. mater(2)) then
        vip(8) = 2
    else if (vip(1) .lt. mater(3)) then
        vip(8) = 3
    else
        vip(8) = 4
    end if
!
! =================================================================
! --- VARIABLE DE POST-TRAITEMENT POUR SUIVRE L'EVOLUTION DE
! --- L'ETAT DE CONTRAINTE PAR RAPPORT AUX DIFFERENTS SEUILS
! --- INDIC = 1 : LA POSITION DE L'ETAT DE CONTRAINTE EST EN-DESSOUS
! ---           : DU SEUIL D'ENDOMMAGEMENT INITIAL ET AU-DESSUS DU
! ---           : SEUIL DE VISCOSITE MAXIMAL (ATTENTION NOTION
! ---           : DIFFERENTE DE LA DEFINITION INITIALE)
! --- INDIC = 2 : LA POSITION DE L'ETAT DE CONTRAINTE EST EN-DESSOUS
! ---           : DU SEUIL D'ENDOMMAGEMENT INITIAL ET EN-DESSOUS DU
! ---           : SEUIL DE VISCOSITE MAXIMAL
! --- INDIC = 3 : LA POSITION DE L'ETAT DE CONTRAINTE EST AU-DESSUS
! ---           : DU SEUIL D'ENDOMMAGEMENT INITIAL ET EN-DESSOUS DU
! ---           : SEUIL DE VISCOSITE MAXIMAL
! --- INDIC = 4 : LA POSITION DE L'ETAT DE CONTRAINTE EST AU-DESSUS
! ---           : DU SEUIL D'ENDOMMAGEMENT INITIAL ET AU-DESSUS DU
! ---           : SEUIL DE VISCOSITE MAXIMAL
! =================================================================
    call lcdevi(sigf, devsig)
    sii = norm2(devsig(1:ndt))
    i1 = -sigf(1)-sigf(2)-sigf(3)
    rcos3t = -cos3t(devsig, mater(11), lgleps)
!
    crit0 = lkcrit(mater(4), mater(5), mater(6), mater(10), mater(9), mater(12), rcos3t, i1, sii)
    crite = lkcrit(1.0d0, mater(8), mater(6), mater(10), mater(9), mater(12), rcos3t, i1, sii)
!
    if (crit0 .lt. 0.0d0 .and. crite .gt. 0.0d0) then
        vip(9) = 1
    else if (crit0 .lt. 0.0d0 .and. crite .lt. 0.0d0) then
        vip(9) = 2
    else if (crit0 .gt. 0.0d0 .and. crite .lt. 0.0d0) then
        vip(9) = 3
    else
        vip(9) = 4
    end if
end subroutine
