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
subroutine dndiss(ipara, nmnbn, nmplas, nmdpla, nmddpl, &
                  nmprox, deps, newnbn, newpla, newdpl, &
                  newddp, newzfg, despit, ddisit, dc1, &
                  dc2, dtg, normm, normn)
!
    implicit none
!
!     CALCUL DES EFFORTS
!     CALCUL DES COURBURES PLASTIQUES
!     CALCUL DE LA DISSIPATION PLASTIQUE
!
! IN  NMNBN : FORCE - BACKFORCE
! IN  NMPLAS : MOMENTS LIMITES DE PLASTICITE
! IN  NMDPLA : DERIVEES DES MOMENTS LIMITES DE PLASTICITE
! IN  NMDDPL : DERIVEES SECONDES DES MOMENTS LIMITES DE PLASTICITE
! IN  NMPROX : NMPROX > 0 : NBN DANS ZONE DE CRITIQUE
! IN  DC1 : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
! IN  DC2 : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
! IN  DTG : MATRICE TANGENTE
! IN  NORMM : NORME SUR LA FONCTION MP = F(N)
! IN  NORMN : NORME SUR LA FONCTION MP = F(N)
!
! IN/OUT IPARA : LISTES DE PARAMETRES DE TYPE ENTIER
!                (ZIMAT NCRIT NEWIEF IER)
!
! OUT NEWNBN : FORCE - BACKFORCE A T+1
! OUT NEWPLA : MOMENTS LIMITES DE PLASTICITE A T+1
! OUT NEWDPL : DERIVEES DES MOMENTS LIMITES DE PLASTICITE A T+1
! OUT NEWDDP : DERIVEES SECONDES DES MOMENTS LIMITES DE PLASTICITE A T+1
! OUT NEWZFG : ZEROS ADIMENSIONNELS POUR LE CRITERE F ET G
! OUT NEWPRO : NEWPRO > 0 : NBN DANS ZONE DE CRITIQUE
! OUT DESPIT : INCREMENT DE DEFORMATION PLASTIQUE
! OUT DDISIT : INCREMENT DE DISSIPATION PLASTIQUE
!
#include "asterfort/assert.h"
#include "asterfort/d0mpfn.h"
#include "asterfort/d1crit.h"
#include "asterfort/d1cro2.h"
#include "asterfort/d2crit.h"
#include "asterfort/d2cro2.h"
#include "asterfort/ddmpfn.h"
#include "asterfort/mppffn.h"
#include "asterfort/r8inir.h"
#include "asterfort/restzo.h"
    integer(kind=8) :: newief, nmprox(*)
    integer(kind=8) :: ncrit, zimat, ier, cief
    integer(kind=8) :: cier, i, j, ipara(4)
!
    real(kind=8) :: nmnbn(*), newnbn(*), nmplas(2, *), newpla(2, *)
    real(kind=8) :: nmdpla(2, *), newdpl(2, *), nmddpl(2, *), newddp(2, *)
    real(kind=8) :: newzef, newzeg, newzfg(2)
    real(kind=8) :: deps(*), dc1(6, 6), dc2(6, 6), dtg(6, 6)
    real(kind=8) :: despit(*), ddisit, tdespi(1, 6), aux(1)
    real(kind=8) :: cnbn(6), cplas(2, 3), czef, czeg, cdeps(6)
    real(kind=8) :: cdtg(6, 6), cdepsp(6)
    real(kind=8) :: nmp(6), normm, normn, rpara(3)
!
    newzef = newzfg(1)
    newzeg = newzfg(2)
    czef = 0.d0
    czeg = 0.d0
!
    zimat = ipara(1)
    ncrit = ipara(2)
    newief = ipara(3)
    ier = ipara(4)
!
    do j = 1, 6
        cdeps(j) = deps(j)
        do i = 1, 6
            cdtg(i, j) = dtg(i, j)
        end do
    end do
!
    if (ncrit .eq. -1) then
!     PREDICTEUR ELASTIQUE HORS CHAMP
        ier = 3
!
    else if (ncrit .eq. 0) then
!     ELASTICITE = AUCUN CRITERE ACTIVE
        call r8inir(6, 0.0d0, despit, 1)
        nmp = matmul(dtg, cdeps)
!
        do j = 1, 6
            newnbn(j) = nmnbn(j)+nmp(j)
        end do
!
!     CALCUL DES MOMENTS LIMITES DE PLASTICITE
!     ET DES ZEROS DES CRITERES
        call mppffn(zimat, newnbn, newpla, newzef, newzeg, &
                    newief, normm)
!
        if (newief .gt. 0) then
            ier = 2
        else
            ier = 0
        end if
    else if (ncrit .eq. 1) then
!     CRITERE 1 ACTIVE
!
!     EST-ON PROCHE DU SOMMET DU CONE
        nmprox(1) = restzo(zimat, nmnbn, 1, normm, normn)
!
        if (nmprox(1) .gt. 0) then
            rpara(1) = czef
            rpara(2) = czeg
            rpara(3) = normm
!
!     CALCUL DU MULTIPLICATEUR PLASTIQUE
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
            call d1cro2(zimat, nmnbn, nmplas, nmdpla, nmddpl, &
                        nmprox, cnbn, cplas, rpara, cief, &
                        cdeps, cdtg, cier, cdepsp, dc1, &
                        1)
!
            czef = rpara(1)
            czeg = rpara(2)
            normm = rpara(3)
!
            if (cier .eq. 3) then
                rpara(1) = czef
                rpara(2) = czeg
                rpara(3) = normm
!
!     CALCUL DU MULTIPLICATEUR PLASTIQUE
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
                call d1crit(zimat, nmnbn, nmplas, nmdpla, nmprox, &
                            cnbn, cplas, rpara, cief, cdeps, &
                            cdtg, cier, cdepsp, dc1, 1)
                czef = rpara(1)
                czeg = rpara(2)
                normm = rpara(3)
            end if
        else
            rpara(1) = czef
            rpara(2) = czeg
            rpara(3) = normm
!
!     CALCUL DU MULTIPLICATEUR PLASTIQUE
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
            call d1crit(zimat, nmnbn, nmplas, nmdpla, nmprox, &
                        cnbn, cplas, rpara, cief, cdeps, &
                        cdtg, cier, cdepsp, dc1, 1)
!
            czef = rpara(1)
            czeg = rpara(2)
            normm = rpara(3)
        end if
!
        do j = 1, 6
            newnbn(j) = cnbn(j)
        end do
!
        do j = 1, 3
            do i = 1, 2
                newpla(i, j) = cplas(i, j)
            end do
        end do
!
        newzef = czef
        newzeg = czeg
        newief = cief
!
        do j = 1, 6
            despit(j) = cdepsp(j)
        end do
!
        ier = cier
!
    else if (ncrit .eq. 2) then
!     CRITERE 2 ACTIVE
        nmprox(2) = restzo(zimat, nmnbn, 2, normm, normn)
!
        if (nmprox(2) .gt. 0) then
            rpara(1) = czef
            rpara(2) = czeg
            rpara(3) = normm
!
!     CALCUL DU MULTIPLICATEUR PLASTIQUE
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
            call d1cro2(zimat, nmnbn, nmplas, nmdpla, nmddpl, &
                        nmprox, cnbn, cplas, rpara, cief, &
                        cdeps, cdtg, cier, cdepsp, dc2, &
                        2)
            czef = rpara(1)
            czeg = rpara(2)
            normm = rpara(3)
!
            if (cier .eq. 3) then
                rpara(1) = czef
                rpara(2) = czeg
                rpara(3) = normm
!
!     CALCUL DU MULTIPLICATEUR PLASTIQUE
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
                call d1crit(zimat, nmnbn, nmplas, nmdpla, nmprox, &
                            cnbn, cplas, rpara, cief, cdeps, &
                            cdtg, cier, cdepsp, dc2, 2)
                czef = rpara(1)
                czeg = rpara(2)
                normm = rpara(3)
            end if
        else
!
            rpara(1) = czef
            rpara(2) = czeg
            rpara(3) = normm
!
!     CALCUL DU MULTIPLICATEUR PLASTIQUE
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
            call d1crit(zimat, nmnbn, nmplas, nmdpla, nmprox, &
                        cnbn, cplas, rpara, cief, cdeps, &
                        cdtg, cier, cdepsp, dc2, 2)
!
            czef = rpara(1)
            czeg = rpara(2)
            normm = rpara(3)
        end if
!
        do j = 1, 6
            newnbn(j) = cnbn(j)
        end do
!
        do j = 1, 3
            do i = 1, 2
                newpla(i, j) = cplas(i, j)
            end do
        end do
!
        newzef = czef
        newzeg = czeg
        newief = cief
!
        do j = 1, 6
            despit(j) = cdepsp(j)
        end do
!
        ier = cier
!
    else if (ncrit .eq. 12) then
!     CRITERES 1 ET 2 ACTIVES
        nmprox(1) = restzo(zimat, nmnbn, 1, normm, normn)
        nmprox(2) = restzo(zimat, nmnbn, 2, normm, normn)
!
        if ((nmprox(1) .gt. 0) .and. (nmprox(2) .gt. 0)) then
            rpara(1) = czef
            rpara(2) = czeg
            rpara(3) = normm
!
!     CALCUL DU MULTIPLICATEUR PLASTIQUE
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
            call d2cro2(zimat, nmnbn, nmplas, nmdpla, nmddpl, &
                        nmprox, cnbn, cplas, rpara, cief, &
                        cdeps, cdtg, cier, cdepsp, dc1, &
                        dc2)
!
            czef = rpara(1)
            czeg = rpara(2)
            normm = rpara(3)
!
            if (cier .gt. 2) then
                rpara(1) = czef
                rpara(2) = czeg
                rpara(3) = normm
!
!     CALCUL DES MULTIPLICATEURS PLASTIQUES
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
                call d2crit(zimat, nmnbn, nmplas, nmdpla, nmprox, &
                            cnbn, cplas, rpara, cief, cdeps, &
                            cdtg, cier, cdepsp, dc1, dc2)
!
                czef = rpara(1)
                czeg = rpara(2)
                normm = rpara(3)
            end if
        else
!
            rpara(1) = czef
            rpara(2) = czeg
            rpara(3) = normm
!
!     CALCUL DES MULTIPLICATEURS PLASTIQUES
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
            call d2crit(zimat, nmnbn, nmplas, nmdpla, nmprox, &
                        cnbn, cplas, rpara, cief, cdeps, &
                        cdtg, cier, cdepsp, dc1, dc2)
!
            czef = rpara(1)
            czeg = rpara(2)
            normm = rpara(3)
!
            if ((cier .eq. 1) .or. (cier .gt. 2)) then
                if ((nmprox(1) .gt. 0) .or. (nmprox(2) .gt. 0)) then
!
                    rpara(1) = czef
                    rpara(2) = czeg
                    rpara(3) = normm
!
!     CALCUL DES MULTIPLICATEURS PLASTIQUES
!     ET DE L INCREMENT DE COURBURE PLASTIQUE
                    call d2cro2(zimat, nmnbn, nmplas, nmdpla, nmddpl, &
                                nmprox, cnbn, cplas, rpara, cief, &
                                cdeps, cdtg, cier, cdepsp, dc1, &
                                dc2)
!
                    czef = rpara(1)
                    czeg = rpara(2)
                    normm = rpara(3)
                end if
            end if
        end if
!
        do j = 1, 6
            newnbn(j) = cnbn(j)
        end do
!
        do j = 1, 3
            do i = 1, 2
                newpla(i, j) = cplas(i, j)
            end do
        end do
!
        newzef = czef
        newzeg = czeg
        newief = cief
!
        do j = 1, 6
            despit(j) = cdepsp(j)
        end do
!
        ier = cier
!
    else
        ASSERT(.false.)
!
    end if
!
    if (ier .eq. 0) then
!     CALCUL DES DERIVEES DES MOMENTS LIMITES DE PLASTICITE
!     CALCUL DES DERIVEES SECONDES DES MOMENTS LIMITES DE PLASTICITE
        call d0mpfn(zimat, newnbn, newdpl)
        call ddmpfn(zimat, newnbn, newddp)
    end if
!
!     CALCUL DE LA DISSIPATION
    do j = 1, 6
        tdespi(1, j) = despit(j)
    end do
!
    do j = 1, 6
        nmp(j) = 0.5d0*(nmnbn(j)+newnbn(j))
    end do
!
    aux = matmul(tdespi, nmp)
!
    ddisit = aux(1)
!
    newzfg(1) = newzef
    newzfg(2) = newzeg
!
    ipara(1) = zimat
    ipara(2) = ncrit
    ipara(3) = newief
    ipara(4) = ier
end subroutine
