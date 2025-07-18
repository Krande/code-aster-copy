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
subroutine cvmmat(fami, kpg, ksp, mod, imat, &
                  nmat, materd, materf, matcst, typma, &
                  ndt, ndi, nr, crit, vim, &
                  nvi, sigd)
    implicit none
!    VISCOCHABOCHE : RECUPERATION DU MATERIAU A T ET T+DT
!                    NB DE CMP DIRECTES/CISAILLEMENT , NB VAR. INTERNES
!                    MATER(*,1) = E , NU , ALPHA
!                    MATER(*,2) = K_0  , A_K , A_R , K    , N   , ALP
!                                 B    , M_R , G_R , MU   , Q_M , Q_0
!                                 QR_0 , ETA , A_I
!                                 M_1  , D1  , G_X1 ,G1_0 , C1
!                                 M_2  , D2  , G_X2 ,G2_0 , C2
!                    VARIABLES INTERNES : X1 , X2 , P , R , Q , XXI , E
!       ----------------------------------------------------------------
!       IN  IMAT   :  ADRESSE DU MATERIAU CODE
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION  DE MATER
!       OUT MATERD :  COEFFICIENTS MATERIAU A T
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!                     MATER(*,1) = CARACTERISTIQUES   ELASTIQUES
!                     MATER(*,2) = CARACTERISTIQUES   PLASTIQUES
!           MATCST :  'OUI' SI  MATERIAU A T = MATERIAU A T+DT
!                     'NON' SINON
!           NDT    :  NB TOTAL DE COMPOSANTES TENSEURS
!           NDI    :  NB DE COMPOSANTES DIRECTES  TENSEURS
!           NR     :  NB DE COMPOSANTES SYSTEME NL
!           NVI    :  NB DE VARIABLES INTERNES
!       ----------------------------------------------------------------
#include "asterfort/rcvalb.h"
#include "asterfort/rcvalt.h"
#include "asterfort/rupmat.h"
!       ----------------------------------------------------------------
    integer(kind=8) :: kpg, ksp, nmat, ndt, ndi, nr, nvi, lgpg
    integer(kind=8) :: ioptio, idnr, i, j, imat
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: epsi, c1d, c2d, vim(*), sigd(6), crit(*)
    character(len=*) :: fami
    character(len=8) :: mod, nomc(28), typma
    integer(kind=8) :: cerr(28)
    character(len=3) :: matcst
    character(len=11) :: meting
!       ----------------------------------------------------------------
    common/opti/ioptio, idnr
    common/meti/meting
    common/coed/c1d, c2d
!       ----------------------------------------------------------------
    data epsi/1.d-15/
!
! -     NB DE COMPOSANTES / VARIABLES INTERNES -------------------------
!
!JMP    APPAREMMENT, NVI EST LE NOMBRE MAXI DE VARIABLES INTERNES
!       X1(6), X2(6), P, R, Q + XXI(6) SI ETA.NE.1.D0 + INDICATEUR
!       CE QUI CONDUIT A NDT+NDT+3(+NDT)+1
!       DU COUP NR DOIT VALOIR :
!          NDT (SIGMA) + 2*NDT+3 SI ETA.EQ.1.D0
!       ET NDT (SIGMA) + 3*NDT+3 SI ETA.NE.1.D0
! -     POUR INTEGRATION PAR METHODE EXPLICITE
! -     ON RESTE DIMENSIONNE EN 3D
!
    if (meting(1:11) .eq. 'RUNGE_KUTTA') then
        ndt = 6
        ndi = 3
        nr = 3*ndt+3
        nvi = 4*ndt+4
! - 3D
    else if (mod(1:2) .eq. '3D') then
        ndt = 6
        ndi = 3
        nr = 3*ndt+3
        nvi = 3*ndt+4
! - D_PLAN AXIS
    else if (mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        ndt = 4
        ndi = 3
        nr = 3*ndt+3
        nvi = 3*ndt+4
! - C_PLAN
    else if (mod(1:6) .eq. 'C_PLAN') then
        ndt = 4
        ndi = 3
        nr = 3*ndt+4
        nvi = 3*ndt+4
    end if
!
!
! -     VISCO-PLASTICITE --->  CALCUL DE LA MATRICE DE COMPORTEMENT
! -     TANGENT  'COHERENT'
!
    typma = 'COHERENT'
!
! -     RECUPERATION MATERIAU -----------------------------------------
!
    nomc(1) = 'E'
    nomc(2) = 'NU'
    nomc(3) = 'ALPHA'
    nomc(4) = 'K_0'
    nomc(5) = 'A_K'
    nomc(6) = 'A_R'
    nomc(7) = 'K'
    nomc(8) = 'N'
    nomc(9) = 'ALP'
    nomc(10) = 'B'
    nomc(11) = 'M_R'
    nomc(12) = 'G_R'
    nomc(13) = 'MU'
    nomc(14) = 'Q_M'
    nomc(15) = 'Q_0'
    nomc(16) = 'QR_0'
    nomc(17) = 'ETA'
    nomc(18) = 'C1'
    nomc(19) = 'M_1'
    nomc(20) = 'D1'
    nomc(21) = 'G_X1'
    nomc(22) = 'G1_0'
    nomc(23) = 'C2'
    nomc(24) = 'M_2'
    nomc(25) = 'D2'
    nomc(26) = 'G_X2'
    nomc(27) = 'G2_0'
    nomc(28) = 'A_I'
!
    do j = 1, 2
        do i = 1, nmat
            materd(i, j) = 0.d0
            materf(i, j) = 0.d0
        end do
    end do
!
! -     RECUPERATION MATERIAU A (T)
!
    call rcvalb(fami, kpg, ksp, '-', imat, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                3, nomc(1), materd(1, 1), cerr(1), 0)
!
!
    if (crit(13) .gt. 0.d0) then
        lgpg = 34
        call rupmat(fami, kpg, ksp, imat, vim, &
                    lgpg, materd(1, 1), sigd)
    end if
!
!
    if (cerr(3) .ne. 0) materd(3, 1) = 0.d0
!
    call rcvalt(fami, kpg, ksp, '-', imat, &
                ' ', 'VISCOCHAB', 0, [' '], [0.d0], &
                5, materd(1, 2), cerr(4), 2)
!
!
!
! -     MISE A JOUR DU COMMUN COED POUR TRAITER LE CAS ANISOTHERME
!
    c1d = materd(15, 2)
    c2d = materd(20, 2)
!
! -     RECUPERATION MATERIAU A (T+DT)
!
    call rcvalb(fami, kpg, ksp, '+', imat, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                3, nomc(1), materf(1, 1), cerr(1), 0)
!
!
    if (crit(13) .gt. 0.d0) then
        lgpg = 34
        call rupmat(fami, kpg, ksp, imat, vim, &
                    lgpg, materf(1, 1), sigd)
    end if
!
    if (cerr(3) .ne. 0) materf(3, 1) = 0.d0
!
    call rcvalt(fami, kpg, ksp, '+', imat, &
                ' ', 'VISCOCHAB', 0, [' '], [0.d0], &
                25, materf(1, 2), cerr(4), 2)
!
! -     PARAMETRES DES LOIS DE COMPORTEMENT A 2 SEUILS
!
!       SI ETA=1, PAS DE MEMOIRE D'ECROUISSAGE
    if (materd(14, 2) .eq. 1.d0) then
        ioptio = 0
        idnr = 0
    else
        ioptio = 2
        idnr = ndt
        nr = nr+ndt
    end if
!
! -     MATERIAU CONSTANT ?
!
    matcst = 'OUI'
    do i = 1, 2
        if (abs(materd(i, 1)-materf(i, 1)) .gt. epsi) then
            matcst = 'NON'
            goto 999
        end if
    end do
    do i = 1, 25
        if (abs(materd(i, 2)-materf(i, 2)) .gt. epsi) then
            matcst = 'NON'
            goto 999
        end if
    end do
!
999 continue
!     NOMBRE DE COEF MATERIAU
    materf(nmat, 2) = 25
    materd(nmat, 2) = 25
!
end subroutine
