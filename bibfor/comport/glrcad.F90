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
subroutine glrcad(zimat, mp1, mp2, delas, rpara, &
                  dmax1, dmax2, dam1, dam2, curvcu, &
                  c1, c2, nbackn, deps, depsp, &
                  df, ddiss, dsidep, normm, normn, &
                  crit, codret)
!
! aslint: disable=W1504
    implicit none
!
!     INTEGRE LA LOI DE COMPORTEMENT GLRC_DAMAGE POUR UN INCREMENT
!     DE DEFORMATION DEPS DEFINIT DANS LE REPERE D ORTHOTROPIE
!
! IN  ZIMAT : ADRESSE DE LA LISTE DE MATERIAU CODE
! IN  DELAS : MATRICE ELASTIQUE EN MEMBRANE, FLEXION ET COUPLAGE
! IN  MP1 : MOMENTS LIMITES ELASTIQUES EN FLEXION POSITIVE
! IN  MP2 : MOMENTS LIMITES ELASTIQUES EN FLEXION NEGATIVE
! IN  RPARA : PARAMETRE MATERIAU DE LA LOI D ENDOMMAGEMENT
! IN  DMAX1 : ENDOMMAGEMENT MAX EN FLEXION +
! IN  DMAX2 : ENDOMMAGEMENT MAX EN FLEXION -
! IN  CURVCU : TENSEUR DES COURBURES ELASTIQUES DANS LE REPERE ORTHO
! IN  C1 : TENSEUR D ECROUISSAGE CINEMATIQUE DE PRAGER EN MEMBRANE
! IN  C2 : TENSEUR D ECROUISSAGE CINEMATIQUE DE PRAGER EN FLEXION
! IN  DEPS : INCREMENT DE DEFORMATION DANS LE REPERE ORTHO
! IN  NORMM : NORME SUR LA FONCTION MP = F(N)
! IN  NORMN : NORME SUR LA FONCTION MP = F(N)
!
! IN/OUT DAM1 : ENDOMMAGEMENT EN FLEXION +
! IN/OUT DAM2 : ENDOMMAGEMENT EN FLEXION -
! IN/OUT NBACKN : EFFORT - EFFORT DE RAPPEL
!
! OUT DEPSP : INCREMENT DE DEFORMATION PLASTIQUE DANS LE REPERE ORTHO
! OUT DF : INCREMENT D EFFORT DANS LE REPERE ORTHO
! OUT DDISS : INCREMENT DE DISSIPATION
! OUT DSIDEP : MATRICE TANGENTE
! person_in_charge: sebastien.fayolle at edf.fr
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/brbagl.h"
#include "asterfort/critnu.h"
#include "asterfort/d0mpfn.h"
#include "asterfort/ddmpfn.h"
#include "asterfort/dndiss.h"
#include "asterfort/dxktan.h"
#include "asterfort/fplass.h"
#include "asterfort/mppffn.h"
#include "asterfort/norrm6.h"
#include "asterfort/r8inir.h"
#include "asterfort/tanmat.h"
#include "blas/dcopy.h"
    aster_logical :: bbok
!
    integer(kind=8) :: ncrit, ncrit2, ier, zimat
    integer(kind=8) :: i, j, kk, kkk, ipara(4), codret, kmax
!
    real(kind=8) :: delas(6, 6), alpha, beta, gamma, k1, k2, dmax1
    real(kind=8) :: dmax2, curvcu(3), c1(6, 6), c2(6, 6), deps(6), mp1(*)
    real(kind=8) :: mp2(*)
    real(kind=8) :: dam1, dam2, nbackn(6), normm, normn
    real(kind=8) :: depsp(6), ddiss, df(6), dsidep(6, 6)
    real(kind=8) :: dc1(6, 6), dc2(6, 6), reps(6), depste(6)
    real(kind=8) :: ddisst, depspt(6), depst2(6), zerode
    real(kind=8) :: dtg(6, 6), curcup(3), dcc1(3, 3), dcc2(3, 3), rpara(5)
    real(kind=8) :: zero, dfp(6), dfp2(6)
    real(kind=8) :: dff(3, 3), crit(*)
!
!---------------------------------------------
    real(kind=8) :: nmnbn(6), newnbn(6)
!         = FORCE - BACKFORCE
    real(kind=8) :: nmplas(2, 3), newpla(2, 3)
!         = PLASMOM(BENDING,_X _Y _XY)
    real(kind=8) :: nmdpla(2, 2), newdpl(2, 2)
!         = DPLASMOM(BENDING,_X _Y)
    real(kind=8) :: nmddpl(2, 2), newddp(2, 2)
!         = DDPLASMOM(BENDING,_X _Y)
    real(kind=8) :: nmzef, newzef
!         ZERO ADIMENSIONNEL POUR LE CRITERE F
    real(kind=8) :: nmzeg, newzeg, newzfg(2)
!         ZERO ADIMENSIONNEL POUR LE CRITERE G
    integer(kind=8) :: nmief, newief
!         NMIEF > 0 : NBN HORS DE LA ZONE DE DEFINITION DE MP
    integer(kind=8) :: nmprox(2), newpro(2)
    blas_int :: b_incx, b_incy, b_n
!         NMPROX > 0 : NBN DANS ZONE DE CRITIQUE
!---------------------------------------------
!
    zero = crit(3)
    kmax = nint(crit(1))
    ncrit = 0
!
    alpha = rpara(1)
    beta = rpara(2)
    gamma = rpara(3)
    k1 = rpara(4)
    k2 = rpara(5)
!
    zerode = zero*norrm6(deps)
!
    do i = 1, 6
!     COPIE DU TENSEUR DES EFFORT - EFFORT DE RAPPEL
        nmnbn(i) = nbackn(i)
    end do
!
!     CALCUL DES MOMENTS LIMITES DE PLASTICITE
!     ET DES ZEROS DES CRITERES
    call mppffn(zimat, nmnbn, nmplas, nmzef, nmzeg, &
                nmief, normm)
!
!     CALCUL DES DERIVEES DES MOMENTS LIMITES DE PLASTICITE
    call d0mpfn(zimat, nmnbn, nmdpla)
!
!     CALCUL DES DERIVEES SECONDES DES MOMENTS LIMITES DE PLASTICITE
    call ddmpfn(zimat, nmnbn, nmddpl)
!
    newzef = nmzef
    newzeg = nmzeg
!
    ASSERT(nmief .le. 0)
!
    do j = 1, 6
        do i = 1, 6
!     DTG : MATRICE TANGENTE
            dtg(i, j) = delas(i, j)
        end do
    end do
!
    ddiss = 0.d0
    call r8inir(6, 0.0d0, df, 1)
    call r8inir(6, 0.0d0, depsp, 1)
!
    do i = 1, 3
        curcup(i) = curvcu(i)
    end do
!
!     METHODE MIXTE POUR S ASSURER
!     D AVOIR f(m,backm)<= 0 A CHAQUE PAS DE TEMPS
!
    do i = 1, 6
!     REPS EST LE RESIDU DE L INCREMENT DE DEFORMATION
        reps(i) = deps(i)
    end do
!
    do j = 1, 6
        do i = 1, 6
!     DC1 : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
!     DC2 : MATRICE ELASTIQUE + CONSTANTES DE PRAGER
            dc1(i, j) = dtg(i, j)+c1(i, j)
            dc2(i, j) = dtg(i, j)+c2(i, j)
        end do
    end do
!
    do kk = 1, kmax
        if (norrm6(reps) .le. zerode) then
!     TEST DE CV DE L ALGO D INTEGRATION
            goto 230
        end if
!
        do i = 1, 6
!     AFFECTATION DE L INCREMENT DE DEFORMATION TEST
            depste(i) = reps(i)
        end do
!
!     CALCUL DE L ENDOMMAGEMENT ET DE LA MATRICE TANGENTE
        call tanmat(alpha, beta, gamma, k1, k2, &
                    dmax1, dmax2, dam1, dam2, curcup, &
                    depste(4), dff)
!
        do j = 1, 3
            do i = 1, 3
                dtg(i+3, j+3) = dff(i, j)
            end do
        end do
!
        do j = 1, 6
            do i = 1, 6
                dc1(i, j) = dtg(i, j)+c1(i, j)
                dc2(i, j) = dtg(i, j)+c2(i, j)
            end do
        end do
!
!     CALCUL DU PREDICTEUR ELASTIQUE ET DU NOMBRE DE CRITERE ATTEINT
        ncrit = critnu(zimat, nmnbn, depste, dtg, normm)
!
        do kkk = 1, kmax
            do j = 1, 6
                depst2(j) = 0.5d0*depste(j)
            end do
!
!     CALCUL DU PREDICTEUR ELASTIQUE ET DU NOMBRE DE CRITERE ATTEINT
            ncrit2 = critnu(zimat, nmnbn, depst2, dtg, normm)
!
            if (ncrit2 .ne. ncrit) then
                do j = 1, 6
                    depste(j) = depst2(j)
                end do
                ncrit = ncrit2
            else
                newzfg(1) = newzef
                newzfg(2) = newzeg
!
                ipara(1) = zimat
                ipara(2) = ncrit
!
!     CALCUL DU NOUVEAU MOMENT
!     DE L INCREMENT DE COURBURE PLASTIQUE ET DE LA DISSIPATION
                call dndiss(ipara, nmnbn, nmplas, nmdpla, nmddpl, &
                            nmprox, depste, newnbn, newpla, newdpl, &
                            newddp, newzfg, depspt, ddisst, dc1, &
                            dc2, dtg, normm, normn)
!
                zimat = ipara(1)
                ncrit = ipara(2)
                newief = ipara(3)
                ier = ipara(4)
!
                newzef = newzfg(1)
                newzeg = newzfg(2)
!
                if (ier .gt. 0) then
                    do j = 1, 6
                        depste(j) = depst2(j)
                    end do
                    ncrit = ncrit2
                else
!     LE POINT EST DANS LA ZONE G < 0
!     SI LE NOUVEAU MOMENT EST A L EXTERIEUR DE LA SURFACE DE PLAST
!     LA METHODE BRINGBACK EST UTILISEE
!     POUR ESSAYER DE LE RAMENER SUR LA SURFACE
!
                    if (fplass(newnbn, newpla, 1) .gt. newzef .or. &
                        fplass(newnbn, newpla, 2) .gt. newzef) then
!
!     PROCEDURE BRINGBACK
                        call brbagl(zimat, newnbn, newpla, newdpl, newddp, &
                                    newzef, newzeg, newief, newpro, depspt, &
                                    ddisst, dc1, dc2, dtg, bbok, &
                                    normm, normn)
!
!     BRINGBACK OK : INCREMENT VALIDE
!
                        if (bbok) goto 123
!
!     BRINGBACK NOK : DICHOTOMIE
                        do j = 1, 6
                            depste(j) = depst2(j)
                        end do
!
                        ncrit = ncrit2
                    else
                        goto 123
                    end if
                end if
            end if
        end do
!
!     NON CONVERGENCE DE L ALGO DE DICHOTOMIE
        codret = 1
!
123     continue
!
!     L INCREMENT EST VALIDE : MISE A JOUR DES VARIABLES
!
        do j = 1, 6
            nmnbn(j) = newnbn(j)
        end do
!
        do j = 1, 3
            do i = 1, 2
                nmplas(i, j) = newpla(i, j)
            end do
        end do
!
        do j = 1, 2
            do i = 1, 2
                nmdpla(i, j) = newdpl(i, j)
                nmddpl(i, j) = newddp(i, j)
            end do
        end do
!
        nmzef = newzef
        nmzeg = newzeg
        nmief = newief
!
        do j = 1, 2
            nmprox(j) = newpro(j)
        end do
!
        do j = 1, 6
            depsp(j) = depsp(j)+depspt(j)
        end do
!
        ddiss = ddiss+ddisst
!
        do j = 1, 3
            curcup(j) = curcup(j)+depste(j+3)-depspt(j+3)
        end do
!
        do j = 1, 6
            dfp2(j) = depste(j)-depspt(j)
        end do
!
        dfp = matmul(dtg, dfp2)
!
        do j = 1, 6
            df(j) = df(j)+dfp(j)
        end do
!
        do j = 1, 6
            reps(j) = reps(j)-depste(j)
        end do
    end do
!
!     NON CONVERGENCE DE L ALGO D INTEGRATION
    codret = 1
!
230 continue
!
    do j = 1, 6
        nbackn(j) = nmnbn(j)
    end do
!
    do i = 1, 3
        do j = 1, 3
            dcc1(j, i) = dc1(3+j, 3+i)
            dcc2(j, i) = dc2(3+j, 3+i)
        end do
    end do
!
    b_n = to_blas_int(36)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, delas, b_incx, dsidep, b_incy)
!
!     REALISE LE CALCUL DE LA MATRICE TANGENTE
    call dxktan(dtg, mp1, mp2, nbackn, ncrit, &
                dcc1, dcc2, dsidep)
!
end subroutine
