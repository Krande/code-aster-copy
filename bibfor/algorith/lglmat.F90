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
subroutine lglmat(mod, imat, nbmat, tempd, materd, &
                  materf, matcst, ndt, ndi, nr, &
                  nvi)
!
    implicit none
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lglnvi.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, nr, nvi, imat, nbmat
    real(kind=8) :: materd(nbmat, 2), materf(nbmat, 2), tempd
    character(len=3) :: matcst
    character(len=8) :: mod
! |---------------------------------------------------------------|
! |-- BUT : RECUPERATION DES DONNEES MATERIAU POUR LA LOI DE -----|
! |------ : COMPORTEMENT DE LAIGLE (MECANIQUE DES ROCHES) --------|
! |---------------------------------------------------------------|
! |----- NB DE CMP DIRECTES/CISAILLEMENT , NB VAR. INTERNES ------|
! |----- MATER(*,1) = E, NU, MU, K -------------------------------|
! |----- MATER(*,2) = GAMMA_ULT, GAMMA_E, M_ULT, M_E, A_E, -------|
! |---------------- : M_PIC, A_PIC, ETA, SIGMA_C, GAMMA, ---------|
! |---------------- : KSI, GAMMA_CJS, SIGMA_P1, SIGMA_P2, PA -----|
! |----- VARIABLE INTERNE : GAMMA_P, EPS_P-  ---------------------|
! |---------------------------------------------------------------|
! =================================================================
! IN  : MOD    : TYPE DE MODELISATION -----------------------------
! --- : IMAT   : ADRESSE DU MATERIAU CODE -------------------------
! --- : NBMAT  : NOMBRE DE PARAMETRES MATERIAU --------------------
! --- : TEMPD  : TEMPERATURE BIDON --------------------------------
! OUT : MATERD : COEFFICIENTS MATERIAU A T ------------------------
! --- : MATERF : COEFFICIENTS MATERIAU A T+DT ---------------------
! ------------ : MATER(*,1) = CARACTERISTIQUES ELASTIQUES ---------
! ------------ : MATER(*,2) = CARACTERISTIQUES PLASTIQUES ---------
! --- : MATCST : 'OUI' --------------------------------------------
! --- : NDT    : NOMBRE TOTAL DE COMPOSANTES DU TENSEUR -----------
! --- : NDI    : NOMBRE DE COMPOSANTES DIRECTES DU TENSEUR --------
! --- : NR     : NOMBRE D'EQUATION DU SYSTEME NL ------------------
! --- : NVI    : NB DE VARIABLES INTERNES -------------------------
! =================================================================
    integer(kind=8) :: ii
    real(kind=8) :: e, nu, mu, k, gamma, ksi, sigc, mult, me, ae
    real(kind=8) :: un, deux, trois, sigmp2, sigmp1, apic, mpic
    real(kind=8) :: cohere
    integer(kind=8) :: cerr(17)
    character(len=16) :: nomc(17)
! =================================================================
! --- INITIALISATION DE PARAMETRES --------------------------------
! =================================================================
    parameter(un=1.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
! =================================================================
    call jemarq()
! =================================================================
! --- NB DE COMPOSANTES / VARIABLES INTERNES ----------------------
! =================================================================
    call lglnvi(mod, ndt, ndi, nvi)
! =================================================================
! - NOMBRE DE CONDITIONS NON-LINEAIRES ----------------------------
! =================================================================
    nr = ndt+4
! =================================================================
! --- DEFINITION DES CHAMPS ---------------------------------------
! =================================================================
    nomc(1) = 'E        '
    nomc(2) = 'NU       '
    nomc(3) = 'ALPHA    '
    nomc(4) = 'GAMMA_ULT'
    nomc(5) = 'GAMMA_E  '
    nomc(6) = 'M_ULT    '
    nomc(7) = 'M_E      '
    nomc(8) = 'A_E      '
    nomc(9) = 'M_PIC    '
    nomc(10) = 'A_PIC    '
    nomc(11) = 'ETA      '
    nomc(12) = 'SIGMA_C  '
    nomc(13) = 'GAMMA    '
    nomc(14) = 'KSI      '
    nomc(15) = 'GAMMA_CJS'
    nomc(16) = 'SIGMA_P1 '
!      NOMC(17) = 'SIGMA_P2 '
    nomc(17) = 'PA       '
! =================================================================
! --- RECUPERATION DES PARAMETRES MATERIAU ------------------------
! =================================================================
    call rcvala(imat, ' ', 'ELAS', 1, 'TEMP', &
                [tempd], 3, nomc(1), materd(1, 1), cerr(1), &
                0)
    call rcvala(imat, ' ', 'LAIGLE', 1, 'TEMP', &
                [tempd], 14, nomc(4), materd(1, 2), cerr(4), &
                0)
! =================================================================
! - CALCUL DES MODULES DE CISAILLEMENT ET DE DEFORMATION VOLUMIQUE-
! =================================================================
    e = materd(1, 1)
    nu = materd(2, 1)
    mu = e/(deux*(un+nu))
    k = e/(trois*(un-deux*nu))
! =================================================================
! - VERIFICATIONS -------------------------------------------------
! =================================================================
    gamma = materd(10, 2)
    ksi = materd(11, 2)
    if ((gamma/ksi) .gt. un) then
        call utmess('F', 'ALGORITH5_11')
    end if
! =================================================================
! --- STOCKAGE DES MODULES CALCULES COMME PARAMETRES MATERIAU -----
! =================================================================
    materd(4, 1) = mu
    materd(5, 1) = k
! =================================================================
! --- CALCUL DE SIGMA_P2 ET DECALAGE DE PA ------------------------
! =================================================================
    materd(15, 2) = materd(14, 2)
    sigc = materd(9, 2)
    mult = materd(3, 2)
    me = materd(4, 2)
    ae = materd(5, 2)
    sigmp2 = sigc*((mult/me**ae)**(un/(ae-un)))
    materd(14, 2) = sigmp2
! =================================================================
! --- VERIFICATION DE LA COHERENCE DES PARAMETRES : ---------------
! --- SIGMA_C, SIGMA_P1, M_PIC, A_PIC, A_E ET M_E -----------------
! =================================================================
    mpic = materd(6, 2)
    apic = materd(7, 2)
    sigmp1 = materd(13, 2)
    cohere =&
     &        abs(sigc/sigmp1*((mpic*sigmp1/sigc+1)**(apic/ae))-me)
    if (cohere .gt. 1.0d-2) then
        call utmess('F', 'ALGORITH5_12')
    end if
! =================================================================
! --- DEFINITION D'UN MATERIAU FINAL ------------------------------
! =================================================================
    do ii = 1, nbmat
        materf(ii, 1) = materd(ii, 1)
        materf(ii, 2) = materd(ii, 2)
    end do
    matcst = 'OUI'
! =================================================================
    call jedema()
! =================================================================
end subroutine
