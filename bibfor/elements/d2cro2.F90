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
subroutine d2cro2(zimat, nmnbn, nmplas, nmdpla, nmddpl, &
                  nmprox, cnbn, cplas, rpara, cief, &
                  cdeps, cdtg, cier, cdepsp, dc1, &
                  dc2)
    implicit none
!
!     CALCUL DU MULTIPLICATEUR PLASTIQUE
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
!     DANS LE CAS OU 1 CRITERE PLASTIQUE EST ACTIVE
!     METHODE EXPLICITE AVEC UNE CONDITION DE SECOND ORDRE
!
! IN  ZIMAT : ADRESSE DE LA LISTE DE MATERIAU CODE
! IN  NMNBN : FORCE - BACKFORCE
! IN  NMPLAS : MOMENTS LIMITES DE PLASTICITE
! IN  NMDPLA : DERIVEES DES MOMENTS LIMITES DE PLASTICITE
! IN  NMDDPL : DERIVEES SECONDES DES MOMENTS LIMITES DE PLASTICITE
! IN  NMPROX : NMPROX > 0 : NBN DANS ZONE DE CRITIQUE
! IN  CDTG : MATRICE TANGENTE
! IN  DC1 : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
! IN  DC2 : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
!
! IN/OUT RPARA : LISTES DE PARAMETRES DE TYPE ENTIER
!
! OUT CNBN : NOUVELLE FORCE - BACKFORCE
! OUT CPLAS : NOUVEAUX MOMENTS LIMITES DE PLASTICITE
! OUT CIEF : NOUVEAU CIEF > 0 : NBN HORS DE LA ZONE DE DEFINITION DE MP
! OUT CDEPS : NOUVEL INCREMENT DE DEFORMATION DANS LE REPERE ORTHO
! OUT CIER : NOUVEAU CODE ERREUR
! OUT CDEPSP : NOUVEL INCREMENT DE DEF PLASTIQUE DANS LE REPERE ORTHO
!
#include "asterfort/dclass.h"
#include "asterfort/dfplgl.h"
#include "asterfort/dfuuss.h"
#include "asterfort/dracs2.h"
#include "asterfort/fplass.h"
#include "asterfort/hplass.h"
#include "asterfort/nmnet2.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: i, j, zimat, nmprox(2)
    integer(kind=8) :: nbroot, iran(8), cief, cier
!
    real(kind=8) :: nmnbn(6), nmplas(2, 3), nmdpla(2, 2)
    real(kind=8) :: nmddpl(2, 2), cnbn(6), cplas(2, 3)
    real(kind=8) :: czef, czeg, cdeps(6), cdtg(6, 6), cdepsp(6)
    real(kind=8) :: dc1(6, 6), dc2(6, 6), normm, tdeph(2, 6), h1(6, 6), h2(6, 6)
    real(kind=8) :: df1(6), df2(6), df(6, 2), tdf(2, 6), dfu1(6), dfu2(6)
    real(kind=8) :: dfu(6, 2)
    real(kind=8) :: ddeps(6), tddeps(1, 6), tdeph1(1, 6), tdeph2(1, 6), dcfu1(6)
    real(kind=8) :: dcfu2(6), tdcfu1(1, 6), tdcfu2(1, 6), auxd1(1), auxe1(1)
    real(kind=8) :: auxe2(1), auxf1(1), auxf2(1), a(2), b(2), c(2), d(2), e(2)
    real(kind=8) :: f(2)
    real(kind=8) :: aa(2), bb(2), cc(2), dd(2), ee(2), ff(2), lambda(2, 2), x(8)
    real(kind=8) :: y(8)
    real(kind=8) :: normxy(8), depsp2(6, 2)
    real(kind=8) :: cp(6), cp2(2, 6), rpara(3)
!
!
    czef = rpara(1)
    czeg = rpara(2)
    normm = rpara(3)
!
!     CALCUL LE GRADIENT DES CRITERES DE PLASICITE
    call dfplgl(nmnbn, nmplas, nmdpla, 1, df1)
    call dfplgl(nmnbn, nmplas, nmdpla, 2, df2)
!
    do j = 1, 6
        df(j, 1) = df1(j)
        df(j, 2) = df2(j)
        tdf(1, j) = df(j, 1)
        tdf(2, j) = df(j, 2)
    end do
!
!     CALUL DES DIRECTIONS DE L ECOULEMENT DES DEFORMATIONS PLASTIQUES
    call dfuuss(nmnbn, nmplas, nmdpla, nmprox, 1, &
                dfu1)
    call dfuuss(nmnbn, nmplas, nmdpla, nmprox, 2, &
                dfu2)
!
    do j = 1, 6
        dfu(j, 1) = dfu1(j)
        dfu(j, 2) = dfu2(j)
    end do
!
!     CALCUL LA MATRICE HESSIENNE DES CRITERES DE PLASTICITE
    call hplass(nmnbn, nmplas, nmdpla, nmddpl, 1, &
                h1)
    call hplass(nmnbn, nmplas, nmdpla, nmddpl, 2, &
                h2)
    ddeps = matmul(cdtg, cdeps)
!
    do j = 1, 6
        tddeps(1, j) = ddeps(j)
    end do
!
    tdeph1 = matmul(tddeps, h1)
    tdeph2 = matmul(tddeps, h2)
!
    do j = 1, 6
        tdeph(1, j) = tdeph1(1, j)
        tdeph(2, j) = tdeph2(1, j)
    end do
!
    dcfu1 = matmul(dc1, dfu1)
    dcfu2 = matmul(dc2, dfu2)
!
    do j = 1, 6
        tdcfu1(1, j) = dcfu1(j)
        tdcfu2(1, j) = dcfu2(j)
    end do
!
    cp = matmul(h1, dcfu2)
    auxd1(1) = dot_product(dcfu1, cp)
    cp = matmul(h1, dcfu1)
    auxe1(1) = dot_product(dcfu1, cp)
    cp = matmul(h2, dcfu1)
    auxe2(1) = dot_product(dcfu1, cp)
    cp = matmul(h1, dcfu2)
    auxf1(1) = dot_product(dcfu2, cp)
    cp = matmul(h2, dcfu2)
    auxf2(1) = dot_product(dcfu2, cp)
!
    do j = 1, 6
        cp2(1, j) = tdf(1, j)+0.5d0*tdeph(1, j)
        cp2(2, j) = tdf(2, j)+0.5d0*tdeph(2, j)
    end do
!
    a = matmul(cp2, ddeps)
    a(1) = a(1)+fplass(nmnbn, nmplas, 1)
    a(2) = a(2)+fplass(nmnbn, nmplas, 2)
!
    do j = 1, 6
        cp2(1, j) = -tdf(1, j)-tdeph(1, j)
        cp2(2, j) = -tdf(2, j)-tdeph(2, j)
    end do
!
    b = matmul(cp2, dcfu1)
    c = matmul(cp2, dcfu2)
!
    d(1) = auxd1(1)
    d(2) = d(1)
    e(1) = 0.5d0*auxe1(1)
    e(2) = 0.5d0*auxe2(1)
    f(1) = 0.5d0*auxf1(1)
    f(2) = 0.5d0*auxf2(1)
!
    aa(1) = a(1)-a(2)
    aa(2) = a(1)+a(2)
    bb(1) = b(1)-b(2)
    bb(2) = b(1)+b(2)
    cc(1) = c(1)-c(2)
    cc(2) = c(1)+c(2)
    dd(1) = d(1)-d(2)
    dd(2) = d(1)+d(2)
    ee(1) = e(1)-e(2)
    ee(2) = e(1)+e(2)
    ff(1) = f(1)-f(2)
    ff(2) = f(1)+f(2)
!
!     EVALUE LES RACINES DU POLYNOME
!     NOUVELLE ROUTINE
    call dracs2(aa, bb, cc, dd, ee, &
                ff, nbroot, x, y)
!  si trop de racines on les affiche
!
    if (nbroot .eq. 0) then
        cier = 4
        call r8inir(6, 0.0d0, cdepsp, 1)
        goto 100
    end if
!
    if (nbroot .ge. 1) then
        do i = 1, nbroot
            normxy(i) = abs(x(i))+abs(y(i))
        end do
!
!     EVALUE L ORDRE DES RACINES DANS IRAN
        call dclass(nbroot, normxy, iran)
    end if
!
    call r8inir(2*2, 0.0d0, lambda, 1)
!
    do i = 1, nbroot
        if ((x(iran(i)) .ge. 0.d0) .and. (y(iran(i)) .ge. 0.d0)) then
            lambda(1, 1) = x(iran(i))
            lambda(2, 2) = y(iran(i))
!
            if (abs(lambda(1, 1)-lambda(2, 2)) .le. 1.d-12*(lambda(1, 1)+lambda(2, 2))) then
                lambda(2, 2) = lambda(1, 1)
            end if
!
            depsp2 = matmul(dfu, lambda)
!
            do j = 1, 6
                cdepsp(j) = depsp2(j, 1)+depsp2(j, 2)
            end do
!
!     CALCUL DE CNBN ET CDEPSP QUAND DEUX CRITERES PLAST SONT ACTIVES
            call nmnet2(zimat, nmnbn, cnbn, cplas, czef, &
                        czeg, cief, cdeps, cdtg, cier, &
                        dc1, dc2, depsp2, normm)
!
            if (cier .eq. 0) goto 100
        end if
    end do
!
    cier = 3
    call r8inir(6, 0.0d0, cdepsp, 1)
!
100 continue
!
    rpara(1) = czef
    rpara(2) = czeg
    rpara(3) = normm
!
end subroutine
