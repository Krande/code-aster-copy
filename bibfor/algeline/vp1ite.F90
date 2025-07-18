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

subroutine vp1ite(lmasse, lraide, ldynam, x, imode, &
                  valp, neq, mxiter, tol, iter, &
                  x0, mx, err, iexcl, place, &
                  iquoti, solveu)
    implicit none
!
#include "jeveux.h"
#include "asterfort/ggubs.h"
#include "asterfort/mrmult.h"
#include "asterfort/resoud.h"
#include "asterfort/vpmort.h"
#include "asterfort/vpstur.h"
!
    integer(kind=8) :: neq
    real(kind=8) :: x(neq, 1), mx(neq, *), err, x0(neq)
    real(kind=8) :: valp
    integer(kind=8) :: place, iexcl(*), imode, mxiter, iter
    integer(kind=8) :: lmasse, lraide, ldynam
    character(len=19) :: solveu
!     CALCUL D'UN COUPLE VECTEUR ET VALEUR PROPRE PAR ITERATION INVERSE
!     ------------------------------------------------------------------
! VAR X      :    :  VECTEUR(S) PROPRE(S)
! VAR VALP : R8 :    VALEUR PROPRE, LA VALEUR INITIALE EST CORRIGEE
!     X0     :    :  VECTEUR PROPRE OBTENU A L'ITERATION PRECEDENTE
! IN  MXITER : IS : NOMBRE MAXIMUM D'ITERATION
! OUT ITER   : IS : NOMBRE D'ITERATION EFFECTUEE
!                   EN CAS DE NON CONVERGENCE ITER = -MXITER
! IN  TOL    : R8 : TOLERENCE (CRITERE DE CONVERGENCE SUR LE MODE)
! OUT ERR    : R8 : ERREUR SUR LE DERNIER ITERE
! IN  SOLVEU : K19 : SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
!     ------------------------------------------------------------------
!     QUOTIENT DE RAYLEIGH GENERALISE
!                Y.A.Y / Y.X  = L + Y.( A.X - L.X) / Y.X
!
!     REFERENCE: F.L. BAUER - J.H. WILKINSON - C. REINSCH
!        HANDBOOK FOR AUTOMATIC COMPUTATION - LINEAR ALGEBRA - VOL.2
!        PAGE 73
!     ------------------------------------------------------------------
!
!
    real(kind=8) :: xmx, x0mx, xxx, det0, pvalp
    real(kind=8) :: coef, coeft, rmg
    complex(kind=8) :: cbid
    character(len=1) :: kbid
    character(len=19) :: k19bid, matass, chcine, criter
!     ------------------------------------------------------------------
!
!     INIT. OBJETS ASTER
!-----------------------------------------------------------------------
    integer(kind=8) :: idet0, ieq, ier, iquoti, jter
    real(kind=8) :: dseed, tol
    integer(kind=8) :: iret
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    matass = zk24(zi(ldynam+1))
    chcine = ' '
    criter = ' '
    k19bid = ' '
!
!     --- VECTEUR INITIAL ALEATOIRE ---
    dseed = 526815.0d0
    call ggubs(dseed, neq, x0)
    call vpmort(neq, x0, x, mx, imode)
    call mrmult('ZERO', lmasse, x0, mx(1, imode), 1, &
                .false._1)
    do ieq = 1, neq
        mx(ieq, imode) = mx(ieq, imode)*iexcl(ieq)
    end do
!
    x0mx = 0.d0
    do ieq = 1, neq
        x0mx = x0mx+x0(ieq)*mx(ieq, imode)
    end do
!
    coef = 1.d0/sqrt(abs(x0mx))
    coeft = sign(1.d0, x0mx)*coef
    do ieq = 1, neq
        x0(ieq) = coef*x0(ieq)
        mx(ieq, imode) = coeft*mx(ieq, imode)
    end do
!
    do jter = 1, mxiter
        iter = jter
!
!        --- ELIMINATION DES DDL EXTERNES ---
        do ieq = 1, neq
            x(ieq, imode) = mx(ieq, imode)*iexcl(ieq)
        end do
!
!        --- RESOLUTION DE (K-W.M) X = (M).X ---
        call resoud(matass, k19bid, solveu, chcine, 1, &
                    k19bid, k19bid, kbid, x(1, imode), [cbid], &
                    criter, .false._1, 0, iret)
!
!        --- ORTHOGONALISATION EN CAS DE MODES MULTIPLES  ---
        call vpmort(neq, x(1, imode), x, mx, imode)
!
!        --- CALCUL DE M.XN ---
        call mrmult('ZERO', lmasse, x(1, imode), mx(1, imode), 1, &
                    .false._1)
        do ieq = 1, neq
            mx(ieq, imode) = mx(ieq, imode)*iexcl(ieq)
        end do
!
!        --- CALCUL DE XN.M.XN ---
        xmx = 0.d0
        do ieq = 1, neq
            xmx = xmx+x(ieq, imode)*mx(ieq, imode)
        end do
!
!        --- NORMALISATION DE XN ---
        coef = 1.d0/sqrt(abs(xmx))
        coeft = sign(1.d0, xmx)*coef
        do ieq = 1, neq
            x0(ieq) = coef*x0(ieq)
            mx(ieq, imode) = coeft*mx(ieq, imode)
        end do
!
!        --- CALCUL DE LA NORME DE XN-1.M.XN ---
        xxx = 0.d0
        do ieq = 1, neq
            xxx = xxx+x0(ieq)*mx(ieq, imode)
        end do
!
!        --- CALCUL DE L'ERREUR ---
        coef = xxx/xmx/coef
        err = abs(abs(xxx)-1.d0)
        if (err .lt. tol) goto 900
!
!        --- SAUVEGARDE DE XN DANS XN-1 ---
        do ieq = 1, neq
            x0(ieq) = x(ieq, imode)
        end do
!
!        --- SHIFT ---
        if (iquoti .gt. 0) then
            pvalp = valp+coef
!
            if (pvalp .gt. valp*0.9d0 .and. pvalp .lt. valp*1.1d0) then
! --- POUR OPTIMISER ON NE CALCULE PAS LE DET
                valp = pvalp
                call vpstur(lraide, valp, lmasse, ldynam, det0, &
                            idet0, place, ier, solveu, .false._1, &
                            .true._1)
            end if
!
        end if
!
    end do
!
!     --- SORTIE SANS CONVERGENCE ---
    iter = -mxiter
900 continue
!
!     --- FREQUENCE CORRIGEE ---
    valp = valp+coef
!
!     --- NORMALISATION DU VECTEUR ---
    rmg = 0.d0
    do ieq = 1, neq
        rmg = max(abs(x(ieq, imode)), rmg)
    end do
    do ieq = 1, neq
        x(ieq, imode) = x(ieq, imode)/rmg
    end do
!
end subroutine
