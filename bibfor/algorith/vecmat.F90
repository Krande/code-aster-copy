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
subroutine vecmat(fami, kpg, ksp, mod, loi, &
                  jmat, nmat, materd, materf, matcst, &
                  typma, ndt, ndi, nr, nvi)
    implicit none
!
!    ---------------------------ATTENTION----------------------------
!    ROUTINE ASSOCIEE A LA LOI VENDOCHAB : VISCOPLASTICITE COUPLE
!    A DE L ENDOMMAGEMENT ISOTROPE
!    ---------------------------ATTENTION----------------------------
!
!    VENDOCHAB : RECUPERATION DU MATERIAU A T(TEMPD) ET T+DT(TEMPF)
!                    NB DE CMP DIRECTES/CISAILLEMENT , NB VAR. INTERNES
!                    MATER(*,1) = E , NU , ALPHA
!                    MATER(*,2) = S , ALPHA_D, BETA_D
!                                 N , UN_SUR_M , UN_SUR_K
!                                 R_D  , A_D   , K_D
!
!                    VARIABLES INTERNES : EVI,  P , R , DMG, ETAT
!    ----------------------------------------------------------------
!           FAMI   :  FAMILLE DE POINT DE GAUSS
!           KPG    :  NUMERO DE POINT DE GAUSS
!           KSP    :  NUMERO DU SOUS-POINT DE GAUSS
!       IN  JMAT   :  ADRESSE DU MATERIAU CODE
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION  DE MATER
!       OUT MATERD :  COEFFICIENTS MATERIAU A T
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!                     MATER(*,1) = CARACTERISTIQUES   ELASTIQUES
!                     MATER(*,2) = CARACTERISTIQUES   PLASTIQUES
!           MATCST :  'NAP' SI  MATERIAU EST UNE S_D NAPPE
!                     'OUI' SI  MATERIAU EST CONSTANT
!                     'NON' SINON
!           NDT    :  NB TOTAL DE COMPOSANTES TENSEURS
!           NDI    :  NB DE COMPOSANTES DIRECTES  TENSEURS
!           NR     :  NB DE COMPOSANTES SYSTEME NL
!           NVI    :  NB DE VARIABLES INTERNES
!    ----------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
    integer(kind=8) :: kpg, ksp, nmat, ndt, ndi, nr, nvi
    integer(kind=8) :: ioptio, idnr
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: epsi
    character(len=8) :: mod, nomc(12), typma
    integer(kind=8) :: cerr(12)
    character(len=3) :: matcst
    character(len=11) :: meting
    character(len=16) :: loi
    character(len=*) :: fami
!   ----------------------------------------------------------------
    common/opti/ioptio, idnr
    common/meti/meting
!   ----------------------------------------------------------------
    integer(kind=8) :: lmat, lfct, ipi, ipif, ik, imat, ivalk, i, j, jpro, il, jmat
    integer(kind=8) :: nbmat
    parameter(lmat=9, lfct=10)
    data epsi/1.d-15/
!   ----------------------------------------------------------------
!
! -     NB DE COMPOSANTES / VARIABLES INTERNES -------------------------
!
! - POUR INTEGRATION PAR METHODE EXPLICITE ON RESTE DIMENSIONNE EN 3D
    if (meting(1:11) .eq. 'RUNGE_KUTTA') then
        ndt = 6
        ndi = 3
        nr = ndt+3
        nvi = ndt+4
! - POUR INTEGRATION PAR METHODE IMPLICITE ON RESTE DIMENSIONNE EN 3D
    else if (meting(1:9) .eq. 'NEWTON') then
        ndt = 6
        ndi = 3
        nr = ndt+2
        nvi = ndt+3
! - 3D
    else if (mod(1:2) .eq. '3D') then
        ndt = 6
        ndi = 3
        nr = ndt+3
        nvi = ndt+4
! - D_PLAN AXIS
    else if (mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        ndt = 4
        ndi = 3
        nr = ndt+3
        nvi = ndt+4
! - C_PLAN
    else if (mod(1:6) .eq. 'C_PLAN') then
        ndt = 4
        ndi = 3
        nr = ndt+4
        nvi = ndt+4
    end if
!
!
! -     VISCO-PLASTICITE-AVEC-ENDOMMAGEMENT
!
    typma = 'COHERENT'
!
! --- VERIFICATION: LE MATERIAU CODE EST IL UNE NAPPE OU PAS?
!     ON REPREND LA PROGRAMMATION DE FOINT1 UTILISE PAR FOINTA
!     UTILISE PAR RCVALB
! --- CALCUL DE JPRO: DANS FOINT1 JPRO=ZI(IPIF+1) OU IPIF EST LE
!     POINTEUR DANS LE MATERIAU CODE UTILISE PAR FOINTA(IPIF,ETC.)
!
!     IPIF EST PASSE DIRECTEMENT EN ARGUMENT DE FOINTA (APPELEE
!     PAR RCVALB) A FOINT1
!
! --- CALCUL DE IPIF: DANS RCVALB ON A
!            IFON = IPI+LMAT-1+LFCT*(IK-1)
!            CALL FOINTA (IFON,NBPAR,NOMPAR,VALPAR,VALRES(IRES))
!
!            ET IPI=ZI(IMAT+2+ZI(IMAT+1)-1)
!            IK=(1..NBF) OU NBF=ZI(IPI+2)
!            (NBF EST LE NOMBRE DE PARAMETRES MAIS ON PREND IK=NBF CAR
!             C'EST LE COEFF NBMATK_D (LE DERNIER) QUI EST SUCEPTIBLE
!            D'ETRE UNE NAPPE
!            ON BOUCLE SUR LES 9 (ZI(IPI+2)) FONCTIONS DU MATERIAU
!            VENDOCHAB, QUAND ON A K_D, ON TESTE SI C EST UNE NAPPE
!
! --- DONC:
!
!
!
!
    matcst = 'OUI'
!
    nbmat = zi(jmat)
!     UTILISABLE SEULEMENT AVEC UN MATERIAU PAR MAILLE
    ASSERT(nbmat .eq. 1)
!
    do j = 1, 2
        do i = 1, nmat
            materd(i, j) = 0.d0
            materf(i, j) = 0.d0
        end do
    end do
!
    if (loi .eq. 'VENDOCHAB') then
!
        imat = jmat+zi(jmat+nbmat+1)
!
        do ik = 1, zi(imat+1)
            if (zk32(zi(imat)+ik-1) .eq. 'VENDOCHAB') then
                ipi = zi(imat+ik+2-1)
                do il = 1, zi(ipi+2)
                    ivalk = zi(ipi+3)
                    if (zk16(ivalk+il-1) (1:3) .eq. 'K_D') then
                        ipif = ipi+lmat-1+lfct*(il-1)
                        jpro = zi(ipif+1)
                        if (zk24(jpro) .eq. 'NAPPE') then
                            matcst = 'NAP'
                        end if
                    end if
                end do
            end if
        end do
!
! -     RECUPERATION MATERIAU -----------------------------------------
!
! --- ATTENTION :
! SI MATCST='NAP' LES VALEURS INITIALES ET FINALES DE K_D (MATERD(9,2)
! ET MATERF(9,2)) SONT PRISES NULLES (0.0D0) PAR DEFFAUT. LA LECTURE
! DES VALEURS DE K_D QUAND CE COEFFICIENT EST RENTRE SOUS FORME DE
! NAPPE EST FAITE  DIRECTEMENT PAR UN APPEL A LA SOUS ROUTINE
! RCVALB DANS LA SOUS ROUTINE DONNANT LES VALEURS DES DERIVEES DES
! VARIABLES INTERNES DANS LA SOUS ROUTINE D'INTEGRATION LOCALE PAR
! UNE METHODE DE RUNGE-KUTTA NMVPRK.F
!
!
        nomc(1) = 'E       '
        nomc(2) = 'NU      '
        nomc(3) = 'ALPHA   '
!
        nomc(4) = 'N       '
        nomc(5) = 'UN_SUR_M'
        nomc(6) = 'UN_SUR_K'
!
        nomc(7) = 'SY      '
        nomc(8) = 'ALPHA_D '
        nomc(9) = 'BETA_D  '
        nomc(10) = 'R_D     '
        nomc(11) = 'A_D     '
        nomc(12) = 'K_D     '
!
!
! -      RECUPERATION MATERIAU A TEMPD (T)
!
        call rcvalb(fami, kpg, ksp, '-', jmat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    3, nomc(1), materd(1, 1), cerr(1), 0)
        if (cerr(3) .ne. 0) then
            materd(3, 1) = 0.d0
        end if
        call rcvalb(fami, kpg, ksp, '-', jmat, &
                    ' ', 'LEMAITRE', 0, ' ', [0.d0], &
                    3, nomc(4), materd(1, 2), cerr(4), 2)
        call rcvalb(fami, kpg, ksp, '-', jmat, &
                    ' ', 'VENDOCHAB', 0, ' ', [0.d0], &
                    5, nomc(7), materd(4, 2), cerr(7), 2)
        if (matcst .eq. 'NAP') then
            materd(9, 2) = 0.0d0
        else
            call rcvalb(fami, kpg, ksp, '-', jmat, &
                        ' ', 'VENDOCHAB', 0, ' ', [0.d0], &
                        1, nomc(12), materd(9, 2), cerr(12), 2)
        end if
!
! -      RECUPERATION MATERIAU A TEMPF (T+DT)
!
        call rcvalb(fami, kpg, ksp, '+', jmat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    3, nomc(1), materf(1, 1), cerr(1), 0)
        if (cerr(3) .ne. 0) then
            materf(3, 1) = 0.d0
        end if
        call rcvalb(fami, kpg, ksp, '+', jmat, &
                    ' ', 'LEMAITRE', 0, ' ', [0.d0], &
                    3, nomc(4), materf(1, 2), cerr(4), 2)
        call rcvalb(fami, kpg, ksp, '+', jmat, &
                    ' ', 'VENDOCHAB', 0, ' ', [0.d0], &
                    5, nomc(7), materf(4, 2), cerr(7), 2)
        if (matcst .eq. 'NAP') then
            materf(9, 2) = 0.0d0
        else
            call rcvalb(fami, kpg, ksp, '+', jmat, &
                        ' ', 'VENDOCHAB', 0, ' ', [0.d0], &
                        1, nomc(12), materf(9, 2), cerr(12), 2)
        end if
!
    else if (loi .eq. 'VISC_ENDO_LEMA') then
!
        nomc(1) = 'E       '
        nomc(2) = 'NU      '
        nomc(3) = 'ALPHA   '
        nomc(4) = 'N       '
        nomc(5) = 'UN_SUR_M'
        nomc(6) = 'UN_SUR_K'
        nomc(7) = 'SY      '
        nomc(8) = 'R_D     '
        nomc(9) = 'A_D     '
!
! -     RECUPERATION MATERIAU A TEMPD (T)
!
        call rcvalb(fami, kpg, ksp, '-', jmat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    3, nomc(1), materd(1, 1), cerr(1), 0)
        if (cerr(3) .ne. 0) then
            materd(3, 1) = 0.d0
        end if
        call rcvalb(fami, kpg, ksp, '-', jmat, &
                    ' ', 'LEMAITRE', 0, ' ', [0.d0], &
                    3, nomc(4), materd(1, 2), cerr(4), 2)
        call rcvalb(fami, kpg, ksp, '-', jmat, &
                    ' ', 'VISC_ENDO', 0, ' ', [0.d0], &
                    3, nomc(7), materd(4, 2), cerr(7), 2)
!
! -     RECUPERATION MATERIAU A TEMPF (T+DT)
!
        call rcvalb(fami, kpg, ksp, '+', jmat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    3, nomc(1), materf(1, 1), cerr(1), 0)
        if (cerr(3) .ne. 0) then
            materf(3, 1) = 0.d0
        end if
        call rcvalb(fami, kpg, ksp, '+', jmat, &
                    ' ', 'LEMAITRE', 0, ' ', [0.d0], &
                    3, nomc(4), materf(1, 2), cerr(4), 2)
        call rcvalb(fami, kpg, ksp, '+', jmat, &
                    ' ', 'VISC_ENDO', 0, ' ', [0.d0], &
                    3, nomc(7), materf(4, 2), cerr(7), 2)
!
    end if
!     NOMBRE DE COEF MATERIAU
    materd(nmat, 2) = 9
    materf(nmat, 2) = 9
!
!
!
! - MATERIAU CONSTANT ? REMARQUE : NON UTILISE POUR L'INSTANT
!
    if (matcst .eq. 'OUI') then
        do i = 1, 2
            if (abs(materd(i, 1)-materf(i, 1)) .gt. epsi) then
                matcst = 'NON'
                goto 999
            end if
        end do
        do i = 1, 9
            if (abs(materd(i, 2)-materf(i, 2)) .gt. epsi) then
                matcst = 'NON'
                goto 999
            end if
        end do
    end if
!
999 continue
end subroutine
