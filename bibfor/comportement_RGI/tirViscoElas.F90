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

subroutine tirViscoElas(fl3d, var0, xmat, inputR, inputVR6, ngf, &
                        deltam, avean, A, B, X, ipzero, &
                        epsk16, epsm16, epse16, &
                        sig16, sigke16, raideur66, we1)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!   tir visco-élastique
!-----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "rgi_module.h"
#include "asterfort/dflufin3d.h"
#include "asterfort/conso3d.h"
#include "asterfort/matfluag3d.h"
#include "asterfort/couplagf3d.h"
#include "asterfort/gauss3d.h"
#include "asterfort/utmess.h"
#include "asterfort/getValVect.h"
#include "asterfort/getR6Mat6.h"
    aster_logical, intent(in) :: fl3d
    real(kind=8), intent(in) :: inputR(*), inputVR6(6, *), raideur66(6, 6)
    integer(kind=8), intent(in) :: ngf
    integer(kind=8), intent(inout) :: ipzero(ngf)
    real(kind=8), intent(in) :: var0(*), xmat(*)
    real(kind=8), intent(out) :: deltam, avean, epsk16(6), epsm16(6)
    real(kind=8), intent(out) :: epse16(6), sig16(6), sigke16(6), we1
    real(kind=8), intent(inout) :: A(ngf, ngf+1), B(ngf), X(ngf)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, errgauss
    real(kind=8) :: bw0, pw0, bg0, pg0, CMp0, dfin0, xflu, dfmx, psik
    real(kind=8) :: ccmin0, ccmax0, cc03(3), vcc33(3, 3), vcc33t(3, 3)
    real(kind=8) :: kveve66(6, 6), kvem66(6, 6), kmve66(6, 6), kmm66(6, 6)
    real(kind=8) :: bve6(6), bm6(6), iden33(3, 3), depsk6(6), depsm6(6)
    real(kind=8) :: deps6r3(6), depse6(6)
    real(kind=8) :: delta, rc, epsm11, epser, CWp, CthP, Cthv
    real(kind=8) :: dt1, theta1, tauk1, taum1
    real(kind=8) :: dsw06(6), epsm06(6), sigke06(6), deps6r2(6)
    real(kind=8) :: epsk06(6), epse06(6), sig06(6)
!-----------------------------------------------------------------------
!
    call getValVect(var0, bw0, pw0, bg0, pg0, vectInd=[BIOW, PSHR, BIOG, PRGI])
    call getValVect(xmat, xflu, dfmx, psik, vectInd=[XFLU, DFMX, YKSY])
    call getValVect(inputR, delta, rc, epsm11, epser, CWp, CthP, Cthv, &
                    dt1, theta1, tauk1, taum1, ind1=1)
    call getR6Mat6(inputVR6, dsw06, epsm06, sigke06, deps6r2, &
                   epsk06, epse06, sig06)
!
!   initialisation de matrice identité
    iden33(:, :) = 0.d0
    do i = 1, 3
        iden33(i, i) = 1.d0
    end do

    if (fl3d) then
!
!       effet du chargement sur le potentiel de fluage
        call dflufin3d(sig06, bw0, pw0, bg0, pg0, &
                       dsw06, delta, rc, xflu, dfin0, &
                       CMp0, dfmx)
!
!       actualisation coeffs de consolidation
        call conso3d(epsm11, epser, ccmin0, ccmax0, epsm06, &
                     epse06, cc03, vcc33, vcc33t, CWp, &
                     CMp0, CthP, Cthv)
!
!       prise en compte de la deformation thermique transitoire
        do i = 1, 6
! E. Cheignon : je remplace ett61(i) par 0 car il restait nul dans fluendo3d
!            varf(96+i)=ett61(i)
!            varf(96+i)=0.d0
! E. Cheignon : depstt6 restait nul dans fluendo3d
!            deps6r3(i)=deps6r2(i)-depstt6(i)
            deps6r3(i) = deps6r2(i)
        end do
!
!       construction des matrices de couplage ds base fixe
        call matfluag3d(epse06, epsk06, sig06, psik, tauk1, &
                        taum1, deps6r2, dt1, theta1, kveve66, &
                        kvem66, kmve66, kmm66, bve6, bm6, &
                        deltam, avean, cc03, vcc33, vcc33t, &
                        iden33, iden33)
!
!       assemblage des matrices de couplage de fluage ds base fixe
        call couplagf3d(A, B, ngf, kveve66, kmm66, &
                        kmve66, kvem66, bve6, bm6)
!
!       resolution du tir visco elastique
        call gauss3d(12, A, X, B, ngf, &
                     errgauss, ipzero)
        if (errgauss .eq. 1) call utmess('E', 'COMPOR3_26')
    else
!       actualisation pas de deformation pour deps therm transitoire
!       prise en compte de la deformation thermique transitoire
!       E.Cheignon : pas d'intéret puisque non utilisé après
!        do i = 1, 3
!            cc03(i)=1.d0
!            do j = 1, 3
!                if (i .eq. j) then
!                    vcc33(i,j)=1.d0
!                    vcc33t(i,j)=1.d0
!                else
!                    vcc33(i,j)=0.d0
!                    vcc33t(i,j)=0.d0
!                end if
!            end do
!        end do
!            prise en compte de la deformation thermique transitoire
        do i = 1, 6
!            varf(96+i)=ett61(i)
!            varf(96+i)=0.d0
!            deps6r3(i)=deps6r2(i)-depstt6(i)
            deps6r3(i) = deps6r2(i)
        end do
!
!           tir elastique
        x(:) = 0.d0
    end if
!
!        recuperation des increments et affectation
    do i = 1, 6
!           increment deformation kelvin
        depsk6(i) = x(i)
        epsk16(i) = epsk06(i)+depsk6(i)
!           ipla deformation maxwell
        depsm6(i) = x(i+6)
        epsm16(i) = epsm06(i)+depsm6(i)
!           increment deformation elastique
        depse6(i) = deps6r3(i)-depsm6(i)-depsk6(i)
        epse16(i) = epse06(i)+depse6(i)
    end do
!
!        etat du materiau apres tir visco elastique
!    E.Cheignon : phi1 non utilisé
!    phi1=phi0
    we1 = 0.d0
    do i = 1, 6
        sig16(i) = sig06(i)
        sigke16(i) = sigke06(i)
        do j = 1, 6
            sig16(i) = sig16(i)+raideur66(i, j)*depse6(j)
            sigke16(i) = sigke16(i)+raideur66(i, j)*depsk6(j)*psik
        end do
!            actualisation de la dissipation visqueuse => EC : non utilisé
!        dphi=0.5d0*(sig06(i)+sig16(i))*(epsm16(i)-epsm06(i))
!        logic1=.false.
!        if (dphi .lt. 0.d0) logic1=.true.
!        phi1=phi1+dphi
!            actualisation du potentiel elastique
        we1 = we1+0.5d0*(sig16(i)*epse16(i))
    end do

end subroutine
