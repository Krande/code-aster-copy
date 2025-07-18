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
subroutine rkdhay(mod, nvi, vini, coeft, nmat, &
                  sigi, dvin, iret)
    implicit none
!      MODELE VISCOPLASTIQUE A ECROUISSAGE ISOTROPE COUPLE A DE
!      L ENDOMMAGEMENT ISOTROPE
! ==================================================================
! INTEGRATION DE LA LOI (HAYHURST) PAR UNE METHODE DE RUNGE KUTTA
!
!     CETTE ROUTINE FOURNIT LA DERIVEE DE L ENSEMBLE DES VARIABLES
!     INTERNES DU MODELE
!       IN  MOD     :  TYPE DE MODELISATION
!           NVI     :  NOMBRE DE VARIABLES INTERNES
!           VINI    :  VARIABLES INTERNES A T
!           COEFT   :  COEFFICIENTS MATERIAU INELASTIQUE A T
!           NMAT    :  DIMENSION MAXI DE COEFT
!           SIGI    :  CONTRAINTES A L'INSTANT COURANT, AVEC SQRT(2)
!     OUT:
!           DVIN    :  DERIVEES DES VARIABLES INTERNES A T
!           IRET    :  CODE RETOUR =0 SI OK, =1 SI PB
!     ----------------------------------------------------------------
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/fgequi.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    character(len=8) :: mod
!
    integer(kind=8) :: iret, itens, ndi, nmat, nvi, ndt, ndim
!
    real(kind=8) :: coeft(nmat), vini(nvi), dvin(nvi), smx(6), sigi(*)
    real(kind=8) :: ecrou(2), h, dmg, dmgmi
    real(kind=8) :: devi(6), devcum, decrou(2), ddmg, ddmgmi
    real(kind=8) :: ze, td, sinn, grj0
    real(kind=8) :: eps0, pk, h1, h2, delta1, delta2, h1st, h2st, pkc
    real(kind=8) :: sig0, biga, alphad
    real(kind=8) :: trsig, grj2v, grj1, epsi, terme1, shmax, sequi
    real(kind=8) :: equi(17), rmin, sequid
    blas_int :: b_incx, b_incy, b_n
!     ----------------------------------------------------------------
    parameter(ze=0.0d0)
    parameter(td=1.5d0)
!
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
!
! -- TEMPERATURE
!
!      IF (MOD(1:2).EQ.'3D')THEN
    if (ndt .eq. 6) then
        ndim = 3
    else if (ndt .eq. 4) then
        ndim = 2
        smx(5) = 0.d0
        smx(6) = 0.d0
        devi(5) = 0.d0
        devi(6) = 0.d0
    end if
    iret = 0
    rmin = r8miem()
    shmax = 50.d0
!
! --    COEFFICIENTS MATERIAU INELASTIQUES
!
    eps0 = coeft(1)
    pk = coeft(2)
    h1 = coeft(3)
    h2 = coeft(4)
    delta1 = coeft(5)
    delta2 = coeft(6)
    h1st = coeft(7)
    h2st = coeft(8)
    biga = coeft(9)
    sig0 = coeft(10)
    pkc = coeft(11)
    alphad = coeft(12)
    sequid = coeft(13)
!
    epsi = r8prem()*pk
!
! --  VARIABLES INTERNES
!
    do itens = 1, 2
        ecrou(itens) = vini(itens+7)
    end do
!
    dmg = vini(11)
    dmgmi = vini(10)
    h = ecrou(1)+ecrou(2)
!
!----------------------------------------------------------------
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sigi, b_incx, smx, b_incy)
!
!------------CALCUL DES INVARIANTS DE CONTRAINTE  -------
!     attention FGEQUI ne prend pas en compte les SQRT(2)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, 1.d0/sqrt(2.d0), smx(4), b_incx)
    call fgequi(smx, 'SIGM_DIR', ndim, equi)
!     on retablit le tenseur
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, sqrt(2.d0), smx(4), b_incx)
    trsig = equi(16)
    grj0 = max(equi(3), equi(4))
    grj0 = max(grj0, equi(5))
    grj1 = trsig
    grj2v = equi(1)
    if (sequid .eq. 0.d0) then
        sequi = grj0
    else
        sequi = grj1
    end if
!------------ CALCUL DU TENSEUR DEVIATORIQUE DES CONTRAINTES ---
    do itens = 1, ndt
        if (itens .le. 3) smx(itens) = smx(itens)-grj1/3.d0
    end do
!
!----- EQUATION DONNANT LA DERIVEE DU MICROENDOMMAG
    ddmgmi = (pkc/3)*((1-dmgmi)**4)
!
!----- EQUATION DONNANT LA DERIVEE DE LA DEF VISCO PLAST
!----- CUMULEE
!
    terme1 = (grj2v*(1-h))/(pk*(1-dmgmi)*(1-dmg))
    if (grj2v .le. epsi) then
        devcum = ze
    else if (abs(terme1) .lt. shmax) then
        devcum = eps0*(sinh(terme1))
    else
        iret = 1
        goto 999
    end if
!
!----- EQUATION DONNANT LA DERIVEE DE H
!
    if (grj2v .le. epsi) then
!       DIVISION PAR ZERO EVITEE
        decrou(1) = ze
        decrou(2) = ze
    else
        if (ecrou(1) .le. (h1st-rmin)) then
            decrou(1) = (h1/grj2v)*(h1st-(delta1*ecrou(1)))*devcum
            decrou(2) = (h2/grj2v)*(h2st-(delta2*ecrou(2)))*devcum
        else
            iret = 1
            goto 999
        end if
    end if
!
!
!----- EQUATION DONNANT LA DERIVEE DE L ENDOMMAGEMENT
!
    if (sequi .ge. ze) then
        sinn = alphad*sequi+((1.d0-alphad)*grj2v)
    else
        sinn = (1.d0-alphad)*grj2v
    end if
    if ((sinn/sig0) .lt. shmax) then
        ddmg = biga*sinh(sinn/sig0)
    else
        iret = 1
        goto 999
    end if
!
!------ EQUATION DONNANT LA DERIVEE DE LA DEF VISCO PLAST
!
    if (grj2v .le. epsi) then
        do itens = 1, ndt
            devi(itens) = ze
        end do
    else
        do itens = 1, ndt
            devi(itens) = td*devcum*smx(itens)/grj2v
        end do
    end if
!
! --    DERIVEES DES VARIABLES INTERNES
!
    do itens = 1, 6
        dvin(itens) = devi(itens)
    end do
    dvin(7) = devcum
    dvin(8) = decrou(1)
    dvin(9) = decrou(2)
    dvin(10) = ddmgmi
    dvin(11) = ddmg
! VARIABLE INTERNE INUTILE CAR ON Y STOCKE L'INDICATEUR DE PLASTICITE
! DANS LCDPEC
    dvin(12) = ze
!
999 continue
end subroutine
