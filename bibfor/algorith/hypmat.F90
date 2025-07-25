! --------------------------------------------------------------------
! Copyright (C) 2005 UCBL LYON1 - T. BARANGER     WWW.CODE-ASTER.ORG
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
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

subroutine hypmat(fami, kpg, ksp, poum, imate, &
                  c10, c01, c20, k)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
    character(len=*) :: fami, poum
    integer(kind=8) :: kpg, ksp, imate
    real(kind=8) :: c10, c01, c20
    real(kind=8) :: k
!
! ----------------------------------------------------------------------
!
! LOI DE COMPORTEMENT HYPERELASTIQUE DE SIGNORINI
!
! C10 (I1-3) + C01 (I2-3)+ C20 (I1-3)^2 + K/2(J-1)²
!
! RECUPERE DONNEES MATERIAUX
!
! ----------------------------------------------------------------------
!
!
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  FAMI    : FAMILLE DE POINTS DE GAUSS
! IN  KPG     : NUMERO DU POINT DE GAUSS
! IN  KSP     : NUMERO DU SOUS-POINT DE GAUSS
! IN  POUM    : '-' POUR VARIABLES DE COMMANDE
! OUT C10     : PARAMETRE LOI DE COMPORTEMENT
! OUT C02     : PARAMETRE LOI DE COMPORTEMENT
! OUT C20     : PARAMETRE LOI DE COMPORTEMENT
! OUT K       : MODULE DE COMPRESSIBILITE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbres
    parameter(nbres=3)
    character(len=16) :: nomres(nbres)
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    real(kind=8) :: nu
    real(kind=8) :: denom
    aster_logical :: cmpk
!
! ----------------------------------------------------------------------
!
! --- INITIALISATIONS
!
!
! --- SOIT ON PREND LE K DONNE PAR DEFI_MATERIAU, SOIT ON LE CALCULE
! --- A PARTIR DU NU
!
    nomres(1) = 'K'
    call rcvalb(fami, kpg, ksp, poum, imate, &
                ' ', 'ELAS_HYPER', 0, ' ', [0.d0], &
                1, nomres, valres, codres, 0)
!
    if (codres(1) .eq. 0) then
        k = valres(1)
        cmpk = .false.
    else
        nomres(1) = 'NU'
        call rcvalb(fami, kpg, ksp, poum, imate, &
                    ' ', 'ELAS_HYPER', 0, ' ', [0.d0], &
                    1, nomres, valres, codres, 2)
        nu = valres(1)
        denom = 3.d0*(1.d0-2.d0*nu)
        if (denom .le. r8prem()) then
            call utmess('F', 'ELASHYPER_98')
        end if
        cmpk = .true.
    end if
!
! --- RECUPERATION C10, C01 ET C20
!
    nomres(1) = 'C10'
    nomres(2) = 'C01'
    nomres(3) = 'C20'
    call rcvalb(fami, kpg, ksp, poum, imate, &
                ' ', 'ELAS_HYPER', 0, ' ', [0.d0], &
                nbres, nomres, valres, codres, 2)
    c10 = valres(1)
    c01 = valres(2)
    c20 = valres(3)
!
! --- CALCUL DU MODULE DE COMPRESSIBILITE
!
    if (cmpk) then
        k = 4.d0*(c10+c01)*(1.d0+nu)/denom
    end if
!
end subroutine
