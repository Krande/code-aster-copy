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

subroutine matrc2(nbpar, nompar, valpar, kcis, matc, &
                  vectt)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/coqrep.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbpar
    real(kind=8) :: valpar(*), kcis, matc(5, 5), vectt(3, 3)
    character(len=8) :: nompar(*)
!
!
    real(kind=8) :: valres(26)
    integer(kind=8) :: icodre(26)
    character(len=16) :: nomres(26)
    character(len=32) :: phenom
    real(kind=8) :: young, nu, nult, nutl, alpha, beta
    real(kind=8) :: passag(3, 3), pas2(2, 2), dorth(3, 3), work(3, 3), d(3, 3)
    real(kind=8) :: dcis(2, 2), c, s, d2(2, 2), el, et, glt, gtn, delta
    real(kind=8) :: r8bid4(4)
    integer(kind=8) :: i, j, jmate, nbv, jcoqu
!
    do i = 1, 5
        do j = 1, 5
            matc(i, j) = 0.d0
        end do
    end do
!
    call jevech('PMATERC', 'L', jmate)
!
    call rccoma(zi(jmate), 'ELAS', 1, phenom, icodre(1))
!
    if (phenom .eq. 'ELAS') then
        nbv = 2
        nomres(1) = 'E'
        nomres(2) = 'NU'
    else if (phenom .eq. 'ELAS_ORTH') then
        nomres(1) = 'E_L'
        nomres(2) = 'E_T'
        nomres(3) = 'NU_LT'
        nomres(4) = 'G_LT'
        nomres(5) = 'G_TN'
        nbv = 5
    else
        call utmess('F', 'ELEMENTS_45', sk=phenom)
    end if
!
    if (phenom .eq. 'ELAS') then
!
        call rcvalb('RIGI', 1, 1, '+', zi(jmate), &
                    ' ', phenom, nbpar, nompar, valpar, &
                    nbv, nomres, valres, icodre, 1)
!
!     MATERIAU ISOTROPE
!
        young = valres(1)
        nu = valres(2)
!
!     CONSTRUCTION DE LA MATRICE DE COMPORTEMENT MATC : (5,5)
!
        matc(1, 1) = young/(1.d0-nu*nu)
        matc(1, 2) = matc(1, 1)*nu
        matc(2, 1) = matc(1, 2)
        matc(2, 2) = matc(1, 1)
        matc(3, 3) = young/2.d0/(1.d0+nu)
        matc(4, 4) = matc(3, 3)*kcis
        matc(5, 5) = matc(4, 4)
!
    else if (phenom .eq. 'ELAS_ORTH') then
!
! ----   INTERPOLATION DES COEFFICIENTS EN FONCTION DE LA TEMPERATURE
! ----   ET DU TEMPS
!        -----------
        call rcvalb('RIGI', 1, 1, '+', zi(jmate), &
                    ' ', phenom, nbpar, nompar, valpar, &
                    nbv, nomres, valres, icodre, 1)
!
        el = valres(1)
        et = valres(2)
        nult = valres(3)
        glt = valres(4)
        gtn = valres(5)
        nutl = et*nult/el
        delta = 1.d0-nult*nutl
!
        dorth(1, 1) = el/delta
        dorth(1, 2) = nult*et/delta
        dorth(2, 2) = et/delta
        dorth(2, 1) = dorth(1, 2)
        dorth(3, 3) = glt
        dorth(1, 3) = 0.d0
        dorth(2, 3) = 0.d0
        dorth(3, 1) = 0.d0
        dorth(3, 2) = 0.d0
!
! --- RECUPERATION DES ANGLES DETERMINANT LE REPERE UTILISATEUR
! --- PAR RAPPORT AU REPERE GLOBAL :
!     ============================
        call jevech('PCACOQU', 'L', jcoqu)
!
        alpha = zr(jcoqu+1)*r8dgrd()
        beta = zr(jcoqu+2)*r8dgrd()
!
!     CALCUL DU COSINUS ET DU SINUS DE L'ANGLE ENTRE LE REPERE
!     INTRINSEQUE ET LE REPERE UTILISATEUR
        call coqrep(vectt, alpha, beta, r8bid4, r8bid4, &
                    c, s)
!
!
! ----   TENSEUR D'ELASTICITE DANS LE REPERE INTRINSEQUE :
! ----   D_GLOB = PASSAG_T * D_ORTH * PASSAG
!
        do i = 1, 3
            do j = 1, 3
                passag(i, j) = 0.d0
            end do
        end do
        passag(1, 1) = c*c
        passag(2, 2) = c*c
        passag(1, 2) = s*s
        passag(2, 1) = s*s
        passag(1, 3) = c*s
        passag(3, 1) = -2.d0*c*s
        passag(2, 3) = -c*s
        passag(3, 2) = 2.d0*c*s
        passag(3, 3) = c*c-s*s
!
        call utbtab('ZERO', 3, 3, dorth, passag, &
                    work, d)
!
        do i = 1, 3
            do j = 1, 3
                matc(i, j) = d(i, j)
            end do
        end do
!
        dcis(1, 1) = glt
        dcis(1, 2) = 0.d0
        dcis(2, 1) = 0.d0
        dcis(2, 2) = gtn
        pas2(1, 1) = c
        pas2(2, 2) = c
        pas2(1, 2) = s
        pas2(2, 1) = -s
        call utbtab('ZERO', 2, 2, dcis, pas2, &
                    work, d2)
        do i = 1, 2
            do j = 1, 2
                matc(3+i, 3+j) = d2(i, j)
            end do
        end do

!
    end if
!
end subroutine
