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
subroutine d1cro2(zimat, nmnbn, nmplas, nmdpla, nmddpl, &
                  nmprox, cnbn, cplas, rpara, cief, &
                  cdeps, cdtg, cier, cdepsp, dc, &
                  bend)
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
! IN  DC : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
! IN  BEND : FLEXION POSITIVE (1) OU NEGATIVE (-1)
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
#include "asterfort/dfplgl.h"
#include "asterfort/dfuuss.h"
#include "asterfort/draac2.h"
#include "asterfort/fplass.h"
#include "asterfort/hplass.h"
#include "asterfort/nmnet1.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: bend, nbxx, i, j
    integer(kind=8) :: nmprox(2), cief, cier, zimat
!
    real(kind=8) :: nmnbn(6), nmplas(2, 3), nmdpla(2, 2), nmddpl(2, 2)
    real(kind=8) :: cnbn(6), cplas(2, 3), czef, czeg
    real(kind=8) :: cdeps(6), cdtg(6, 6), cdepsp(6)
    real(kind=8) :: dc(6, 6), normm, cp1(1, 6), cp0(1), tdcu(1, 6), h(6, 6)
    real(kind=8) :: df(6)
    real(kind=8) :: lambda, u(6), a0(1), a1(1), a2(1), xx(2), dcu(6), hdcu(6)
    real(kind=8) :: tdf(1, 6), ddeps(6), tddeps(1, 6), rpara(3)
!
    czef = rpara(1)
    czeg = rpara(2)
    normm = rpara(3)
!
!     CALCUL LA MATRICE HESSIENNE DU CRITERE DE PLASTICITE
    call hplass(nmnbn, nmplas, nmdpla, nmddpl, bend, &
                h)
!     CALUL DES DIRECTIONS DE L ECOULEMENT DES DEFORMATIONS PLASTIQUES
    call dfuuss(nmnbn, nmplas, nmdpla, nmprox, bend, &
                u)
!     CALCUL LE GRADIENT DU CRITERE DE PLASICITE
    call dfplgl(nmnbn, nmplas, nmdpla, bend, df)
!
    do j = 1, 6
        tdf(1, j) = df(j)
    end do
!
    ddeps = matmul(cdtg, cdeps)
!
    do j = 1, 6
        tddeps(1, j) = ddeps(j)
    end do
!
    dcu = matmul(dc, u)
!
    do j = 1, 6
        tdcu(1, j) = dcu(j)
    end do
!
    hdcu = matmul(h, dcu)
    cp1 = matmul(tddeps, h)
!
    do j = 1, 6
        cp1(1, j) = tdf(1, j)+0.5d0*cp1(1, j)
    end do
!
    cp0 = matmul(cp1, ddeps)
!
    a0(1) = fplass(nmnbn, nmplas, bend)+cp0(1)
!
    a1 = matmul(tdf, dcu)
    cp0 = matmul(tddeps, hdcu)
    a1(1) = -a1(1)-cp0(1)
    cp0 = matmul(tdcu, hdcu)
    a2(1) = 0.5d0*cp0(1)
!
!     RESOLUTION DE L EQUATION DU SECOND DEGRE
    call draac2(a2(1), a1(1), a0(1), xx(1), xx(2), &
                nbxx)
!
    do i = 1, nbxx
        if (xx(i) .ge. 0.d0) then
            lambda = xx(i)
!
            do j = 1, 6
                cdepsp(j) = lambda*u(j)
            end do
!
!     CALCUL DE CNBN ET CDEPSP QUAND UN CRITERE PLASTIQUE EST ACTIVE
            call nmnet1(zimat, nmnbn, cnbn, cplas, czef, &
                        czeg, cief, cdeps, cdtg, cier, &
                        cdepsp, dc, normm)
!
            if (cier .eq. 0) goto 60
        end if
    end do
!
    cier = 3
!
    call r8inir(6, 0.0d0, cdepsp, 1)
!
60  continue
!
    rpara(1) = czef
    rpara(2) = czeg
    rpara(3) = normm
end subroutine
