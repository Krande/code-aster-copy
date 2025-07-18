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

subroutine calsig(fami, kpg, ksp, ein, mod, &
                  rela_comp, vini, x, dtime, epsd, &
                  detot, nmat, coel, sigi)
!
    implicit none
!
#include "asterfort/lcopli.h"
#include "asterfort/rcvarc.h"
!
!     INTEGRATION DE LOIS DE COMPORTEMENT ELASTO-VISCOPLASTIQUE
!     PAR UNE METHODE DE RUNGE KUTTA
!     A MODIFER SI ELASTCITE ORTHOTROPE
!     CALCUL DES CONTRAINTES A PARTIR DES CHAMPS DE DEFORMATION
!     ----------------------------------------------------------------
!     IN  FAMI    :  FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!         KPG,KSP :  NUMERO DU (SOUS)POINT DE GAUSS
!         EIN     :  DEFORMATION INELASTIQUE
!         MOD     :  TYPE DE MODELISATION
!         COMP    :  COMPORTEMENT
!         VINI    :  VARIABLES INTERNES
!         X       :  INSTANT COURANT
!         DTIME   :  INTERVALLE DE TEMPS
!         EPSD    :  DEFORMATION TOTALE A T
!         DETOT   :  INCREMENT DE DEFORMATION TOTALE
!         NMAT    :  NOMBRE MAXI DE COEFFICIENTS MATERIAU
!         COEL    :  COEFFICENT DE L'OPERATEUR D'ELASTICITE ORTHOTROPE
!     OUT SIGI    :  CONTRAINTES A L'INSTANT COURANT
!     ----------------------------------------------------------------

    character(len=*) :: fami
    character(len=8) :: mod
    character(len=16) :: rela_comp
    integer(kind=8) :: kpg, ksp, icp, nmat, iret, iret1, iret2, iret3
    real(kind=8) :: nu, coel(nmat), hook(6, 6), alphal, alphat, alphan, ethl
    real(kind=8) :: etht
    real(kind=8) :: ein(6), xsdt, x, dtime, eth, alpha, tperd, dtper, tperef
    real(kind=8) :: ethn, dmg
    real(kind=8) :: eel(6), sigi(6), epsd(6), detot(6), demu, e, treel, tf, e0
    real(kind=8) :: vini(*)
!     ----------------------------------------------------------------
!
    call rcvarc(' ', 'TEMP', '-', fami, kpg, &
                ksp, tperd, iret1)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                ksp, tf, iret2)
    call rcvarc(' ', 'TEMP', 'REF', fami, kpg, &
                ksp, tperef, iret3)
!
    dtper = 0.d0
    iret = iret1+iret2+iret3
!
    xsdt = x/dtime
!
    if (nint(coel(nmat)) .eq. 0) then
!
        e = coel(1)
        e0 = e
        dmg = 0.d0
!        ENDOMMAGEMNT EVENTUEL
        if (rela_comp(1:9) .eq. 'VENDOCHAB') then
            dmg = vini(9)
        else if (rela_comp(1:8) .eq. 'HAYHURST') then
            dmg = vini(11)
        end if
        e = e0*(1.d0-dmg)
        nu = coel(2)
        alpha = coel(3)
        if (iret .eq. 0) then
            dtper = tf-tperd
            eth = alpha*(tperd+xsdt*dtper-tperef)
        else
            eth = 0.d0
        end if
        do icp = 1, 6
            eel(icp) = epsd(icp)+detot(icp)*xsdt-ein(icp)-eth
            if (icp .eq. 3) eth = 0.0d0
        end do
!
! --     CAS DES CONTRAINTES PLANES
!
        if (mod(1:6) .eq. 'C_PLAN') then
            eel(3) = -nu*(eel(1)+eel(2))/(1.0d0-nu)
        end if
!
        demu = e/(1.0d0+nu)
        treel = (eel(1)+eel(2)+eel(3))
        treel = nu*demu*treel/(1.0d0-nu-nu)
        do icp = 1, 6
            sigi(icp) = demu*eel(icp)+treel
            if (icp .eq. 3) treel = 0.0d0
        end do
    else if (nint(coel(nmat)) .eq. 1) then
        call lcopli('ORTHOTRO', mod, coel, hook)
        alphal = coel(73)
        alphat = coel(74)
        alphan = coel(75)
        if (iret .eq. 0) then
            ethl = alphal*(tperd+xsdt*dtper-tperef)
            ethn = alphat*(tperd+xsdt*dtper-tperef)
            etht = alphan*(tperd+xsdt*dtper-tperef)
        else
            ethl = 0.d0
            ethn = 0.d0
            etht = 0.d0
        end if
        eel(1) = epsd(1)+detot(1)*xsdt-ein(1)-ethl
        eel(2) = epsd(2)+detot(2)*xsdt-ein(2)-ethn
        eel(3) = epsd(3)+detot(3)*xsdt-ein(3)-etht
        eel(4) = epsd(4)+detot(4)*xsdt-ein(4)
        eel(5) = epsd(5)+detot(5)*xsdt-ein(5)
        eel(6) = epsd(6)+detot(6)*xsdt-ein(6)
        sigi = matmul(hook, eel)
    end if
end subroutine
