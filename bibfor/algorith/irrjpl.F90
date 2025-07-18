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

subroutine irrjpl(model, nmat, mater, sigf, vind, &
                  vinf, dsde)
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/irrfss.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcnrts.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcprte.h"
#include "asterfort/mgauss.h"
    character(len=8) :: model
    integer(kind=8) :: nmat
    real(kind=8) :: mater(nmat, 2), dsde(6, 6), sigf(6)
    real(kind=8) :: vind(*), vinf(*)
!
! person_in_charge: jean-luc.flejou at edf.fr
!       ----------------------------------------------------------------
!       IRRAD3M   :  MATRICE SYMETRIQUE DE COMPORTEMENT TANGENT
!                     ELASTO_PLASTIQUE EN VITESSE A T OU T+DT
!       ----------------------------------------------------------------
!       IN  FAMI   :  FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!           KPG,KSP:  NUMERO DU (SOUS)POINT DE GAUSS
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATER  :  COEFFICIENTS MATERIAU
!           NR     :  TAILLE DE LA MATRICE JACOBIENNE
!           SIGF   :  CONTRAINTES A T+DT
!           VIND   :  VARIABLES INTERNES A T
!           VINF   :  VARIABLES INTERNES A T+DT
!       OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT = DSIG/DEPS
!       ----------------------------------------------------------------
    integer(kind=8) :: ndt, ndi
!     ------------------------------------------------------------------
    common/tdim/ndt, ndi
!
    real(kind=8) :: det, mat(6, 6), i4(6, 6)
    real(kind=8) :: fkooh(6, 6), dphi, pf, etaif, dp, detai, dpi
    real(kind=8) :: ai0, etais, k, n, p0, kappa, r02, zetaf, penpe, pk, pe, spe
    real(kind=8) :: sr
    real(kind=8) :: ddfdds(6, 6), drsds(6, 6), sequiv, dev(6), dfds(6), drids
    real(kind=8) :: drpdp
!
    integer(kind=8) :: iret
    aster_logical :: ldrpdp
!     ------------------------------------------------------------------
    data i4/1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
     &             0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,&
     &             0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,&
     &             0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0,&
     &             0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0,&
     &             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0/
!     ------------------------------------------------------------------
!
!     CALCUL DE LA MATRICE JACOBIENNE ==> methode b (plus sure)
!     a) Faire appel a IRRJAC
!        certaines des équations sont normées   ==> précautions
!        si le jacobien est changé              ==> répercutions
!        CALL IRRJAC (FAMI,KPG,KSP,MOD,NMAT,MATER,YF,DY,NR,DRDY)
!     b) Calcul qui ressemble a IRRJAC
!        independant de IRRJAC
!        adaptation du code de IRRJAC
!
    call lcopil('ISOTROPE', model, mater(1, 1), fkooh)
!     RECUPERATION DES VARIABLES INTERNES A t+
    pf = vinf(1)
    etaif = vinf(2)
!     RECUPERATION DES INCREMENTS DES VARIABLES INTERNES
    dp = vinf(1)-vind(1)
    detai = vinf(2)-vind(2)
    dpi = vinf(3)-vind(3)
!
!     CARACTERISTIQUES MATERIAUX
    ai0 = mater(4, 2)
    etais = mater(5, 2)
    k = mater(7, 2)
    n = mater(8, 2)
    p0 = mater(9, 2)
    kappa = mater(10, 2)
    r02 = mater(11, 2)
    zetaf = mater(12, 2)
    penpe = mater(13, 2)
    pk = mater(14, 2)
    pe = mater(15, 2)
    spe = mater(16, 2)
!     INCREMENT D'IRRADIATION
    dphi = mater(23, 2)
!
!     Calcul de DRSDS : (6,6)
    call irrfss(sigf, ddfdds)
    drsds(1:ndt, 1:ndt) = (dp+dpi)*ddfdds(1:ndt, 1:ndt)
    drsds(1:ndt, 1:ndt) = fkooh(1:ndt, 1:ndt)+drsds(1:ndt, 1:ndt)
!
!     Calcul de DRPDP : scalaire
!     loi de comportement : Calcul du seuil
    if (pf .lt. pk) then
        sr = kappa*r02
    else if (pf .lt. pe) then
        sr = penpe*(pf-pe)+spe
    else
        sr = k*((pf+p0)**n)
    end if
!     Calcul de Sigma equivalent
    call lcdevi(sigf, dev)
    sequiv = lcnrts(dev)
    ldrpdp = .true.
    if (((sequiv .ge. sr) .and. (dp .ge. 0.d0)) .or. (dp .gt. r8prem())) then
        if (pf .lt. pk) then
            ldrpdp = .false.
        else if (pf .lt. pe) then
            drpdp = penpe
        else
            drpdp = n*k*((pf+p0)**(n-1.d0))
        end if
    else
        ldrpdp = .false.
    end if
!     Calcul de DRIDS : scalaire
    if ((etaif-detai) .gt. etais) then
        drids = ai0*dphi*zetaf
    else if (etaif .le. etais) then
        drids = 0.0d0
    else if (detai .gt. 0.0d0) then
        drids = ai0*dphi*zetaf*(etaif-etais)/detai
    else
        drids = 0.0d0
    end if
    if (sequiv .eq. 0.0d0) then
        ddfdds(:, :) = 0.0d0
    else
        dfds(1:ndt) = (1.5d0/sequiv)*dev(1:ndt)
        call lcprte(dfds, dfds, ddfdds)
    end if
!
    if (ldrpdp) then
        ddfdds = (1.0d0/drpdp+drids)*ddfdds
    else
        ddfdds = drids*ddfdds
    end if
!
!     Assemblage de DRSDS et DDFDDS : (6,6)
    mat(1:ndt, 1:ndt) = drsds(1:ndt, 1:ndt)+ddfdds(1:ndt, 1:ndt)
!
!     Inversion de MAT : DSDE(6,6)
    dsde(1:ndt, 1:ndt) = i4(1:ndt, 1:ndt)
    call mgauss('NFVP', mat, dsde, 6, ndt, &
                ndt, det, iret)
!
end subroutine
