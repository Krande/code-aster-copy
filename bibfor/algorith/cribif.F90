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
subroutine cribif(mod, dsidep, vbifur, nbrac4, racine)
    implicit none
#include "asterc/r8nnem.h"
#include "asterc/r8prem.h"
#include "asterc/r8rddg.h"
#include "asterfort/fbifur.h"
#include "asterfort/utmess.h"
#include "asterfort/zerop3.h"
#include "asterfort/zeropn.h"
    integer(kind=8) :: nbrac4
    real(kind=8) :: dsidep(6, 6), racine(4), vbifur
    character(len=8) :: mod
! =====================================================================
! --- RECHERCHE DE ZONES DE LOCALISATION PAR LE CRITERE DE RICE -------
! =====================================================================
    integer(kind=8) :: ii, degre, compt, nbrac3, ier
    real(kind=8) :: zero, un, deux, trois, quatre
    real(kind=8) :: a0, a1, a2, a3, a4, lamba, lambb, lambc
    real(kind=8) :: valeur
    real(kind=8) :: ai(4), rac4(8), signe, rac3(3)
! =====================================================================
    parameter(zero=0.0d0)
    parameter(un=1.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
    parameter(quatre=4.0d0)
! =====================================================================
! --- INITIALISATIONS ET COHERENCES -----------------------------------
! =====================================================================
    vbifur = zero
    nbrac4 = 0
    signe = 1.0d0
    valeur = r8nnem()
    racine(1) = 0.0d0
    racine(2) = 0.0d0
    racine(3) = 0.0d0
    racine(4) = 0.0d0
    if ((mod(1:6) .ne. 'D_PLAN') .and. (mod(1:6) .ne. 'C_PLAN') .and. (mod(1:4) .ne. 'AXIS')) then
        call utmess('F', 'ALGORITH2_43')
    end if
! =====================================================================
! --- AFFECTATION DES VARIABLES ---------------------------------------
! =====================================================================
    a0 = dsidep(1, 1)*dsidep(4, 4)-dsidep(1, 4)*dsidep(4, 1)
   a1 = dsidep(1, 1)*(dsidep(4, 2)+dsidep(2, 4))-dsidep(1, 4)*dsidep(2, 1)-dsidep(1, 2)*dsidep(4, 1&
           &)
    a2 = dsidep(1, 1)*dsidep(2, 2)+dsidep(1, 4)*dsidep(4, 2)+dsidep(4, 1)*dsidep(2, 4)-dsidep(1, 2&
         &)*dsidep(4, 4)-dsidep(1, 2)*dsidep(2, 1)-dsidep(4, 4)*dsidep(2, 1)
    a3 = dsidep(2, 2)*(dsidep(1, 4)+dsidep(4, 1))-dsidep(1, 2)*dsidep(2, 4)-dsidep(4, 2)*dsidep(2&
         &, 1)
    a4 = dsidep(4, 4)*dsidep(2, 2)-dsidep(4, 2)*dsidep(2, 4)
    if (a4 .lt. -r8prem()) then
        signe = -1.0d0
    end if
! =====================================================================
! --- CAS OU A4 = 0 ---------------------------------------------------
! =====================================================================
    if (abs(a4) .lt. r8prem()) then
! =====================================================================
! --- ON LOCALISE POUR A4 = 0 -----------------------------------------
! =====================================================================
        vbifur = un
        nbrac4 = 1
        racine(1) = 90.0d0
        racine(2) = 0.0d0
        racine(3) = 0.0d0
        racine(4) = 0.0d0
        goto 9998
    end if
    lamba = trois*a3/quatre/a4
    lambb = a2/deux/a4
    lambc = a1/quatre/a4
! =====================================================================
! --- RESOLUTION DU POLYNOME DE DEGRE 3 -------------------------------
! =====================================================================
    call zerop3(lamba, lambb, lambc, rac3, nbrac3)
    do ii = 1, nbrac3
        valeur = signe*fbifur(a0, a1, a2, a3, a4, rac3(ii))
        if (valeur .lt. -r8prem()) vbifur = un
    end do
! =====================================================================
! --- RECHERCHE DES RACINES DU POLYNOME -------------------------------
! =====================================================================
    if (vbifur .eq. un) then
        degre = 4
        ai(1) = a0/a4
        ai(2) = a1/a4
        ai(3) = a2/a4
        ai(4) = a3/a4
        call zeropn('F', degre, ai(1), rac4, ier)
! =====================================================================
! --- ON RECUPERE LES RACINES REELLES ---------------------------------
! =====================================================================
        compt = 0
        do ii = 1, 4
            if (abs(rac4((ii-1)*2+2)) .lt. r8prem()) then
                compt = compt+1
                racine(compt) = atan2(rac4((ii-1)*2+1), 1.0d0)
                racine(compt) = racine(compt)*r8rddg()
            end if
        end do
        nbrac4 = compt
    end if
! =====================================================================
9998 continue
! =====================================================================
end subroutine
