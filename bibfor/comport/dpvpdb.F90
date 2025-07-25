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
subroutine dpvpdb(nbmat, mater, crit, dt, vinm, &
                  vinp, nvi, seqe, i1e, seqm, &
                  i1m, dp, nbre, retcom)
    implicit none
#include "asterfort/dpvpdf.h"
#include "asterfort/dpvpdv.h"
#include "asterfort/dpvpeq.h"
#include "asterfort/dpvpva.h"
    integer(kind=8) :: nbmat, nvi, nbre, retcom
    real(kind=8) :: mater(nbmat, 2)
    real(kind=8) :: crit(3), dt
    real(kind=8) :: vinm(nvi), vinp(nvi)
    real(kind=8) :: seqe, i1e, seqm, i1m, dp
! =====================================================================
! --- IN --- : NBMAT   NOMBRE DE PARAMETRES DU MODELE -----------------
! ---------- : MATER   COEFFICIENTS MATERIAU --------------------------
! ---------- : CRIT    TABLEAU DES PARAMETRES DE CONVERGENCE ----------
! ---------- : DT      PAS DE TEMPS -----------------------------------
! ---------- : VINM    VARIABLES INTERNES AU TEMPS MOINS --------------
! ---------- : VINP    VARIABLES INTERNES AU TEMPS PLUS ---------------
! ---------- : NVI     NOMBRE DE VI -----------------------------------
! ---------- : SEQE   CONTRAINTE EQUIVALENTE DE LA PREDICTION ELASTIQUE
! ---------- : I1E    TRACE DE LA PREDICTION ELASTIQUE ----------------
! ---------- : SEQM   CONTRAINTE EQUIVALENTE A l INSTANT MOINS --------
! ---------- : I1M    TRACE DE LA CONTRAINTE A L INSTANT MOINS---------
! ----OUT -- : DP     INCONNUE - DEFORMATION VISCOPLASTIQUE CUMULEE ---
! ---------  : NBRE   NOMBRE D ITERATIONS POUR LA CONVERGENCE LOCALE --
! ---------  : RETCOM  CODE RETOUR 0 OU 1 SI REDECOUPAGE NECESSAIRE  --
! =====================================================================
! =====================================================================
! --- LOI DE COMPORTEMENT DE TYPE DRUCKER PRAGER VISCOPLASTIQUE -------
! --- VISC_DRUC_PRAG --------------------------------------------------
! --- RESOLUTION NUMERIQUE DE L EQ NON LINEAIRE AVEC BRACKETING ET ----
! --------------LA METHODE DES CORDES
! =====================================================================
    integer(kind=8) :: niter, i
    integer(kind=8) :: signf, signfi
    real(kind=8) :: mu, k
    real(kind=8) :: trois, neuf, zero
    real(kind=8) :: pref, a, n, const
    real(kind=8) :: fonc1, fonc2, fonc3, fonc4
    real(kind=8) :: f, fp, seuil, xinf, xsup, finf, fsup
    real(kind=8) :: fonecp(3), fonecm(3), fonder(3)
    real(kind=8) :: alpham, rm, betam
    real(kind=8) :: dalpdp, drdp, dbetdp
    real(kind=8) :: dp0
    real(kind=8) :: fi
! =====================================================================
    parameter(trois=3.0d0)
    parameter(neuf=9.0d0)
    parameter(zero=0.0d0)
! =====================================================================
! --- AFFECTATION DES VARIABLES ---------------------------------------
! =====================================================================
    mu = mater(4, 1)
    k = mater(5, 1)
    pref = mater(1, 2)
    a = mater(2, 2)
    n = mater(3, 2)
! =====================================================================
!
    const = a*dt/(pref)**n
    retcom = 0
!
! =====================================================================
! --- CALCUL DE DP ----------------------------------------------------
! =====================================================================
    call dpvpva(vinm, nbmat, mater, fonecm)
    call dpvpva(vinp, nbmat, mater, fonecp)
    call dpvpdv(vinp, nbmat, mater, fonder)
!
!
    alpham = fonecm(1)
    rm = fonecm(2)
    betam = fonecm(3)
!
    dalpdp = fonder(1)
    drdp = fonder(2)
    dbetdp = fonder(3)
!
    fonc1 = seqe+alpham*i1e-rm
!
    fonc2 = trois*mu+drdp-dalpdp*i1e+neuf*k*alpham*betam
!
    fonc3 = neuf*k*(alpham*dbetdp+betam*dalpdp)
!
    fonc4 = neuf*k*dalpdp*dbetdp
!
!
    if (fonc1 .gt. zero) then
        fonc1 = fonc1
    else
        fonc1 = zero
    end if
!
    xinf = zero
!
    xsup = a*(abs(fonc1)/pref)**n*dt
!
    finf = dpvpeq(xinf, n, const, fonc1, fonc2, fonc3, fonc4)
!
    fsup = dpvpeq(xsup, n, const, fonc1, fonc2, fonc3, fonc4)
!
!
    niter = int(crit(1))
!
!
    dp0 = xinf
!
!
    f = dpvpeq(dp0, n, const, fonc1, fonc2, fonc3, fonc4)
    fp = dpvpdf(dp0, n, const, fonc1, fonc2, fonc3, fonc4)
!
    seuil = dpvpeq(xinf, n, const, fonc1, fonc2, fonc3, fonc4)
!
    if (abs(finf/seuil) .le. crit(3)) then
        dp0 = xinf
        nbre = 1
        goto 50
    else if (abs(fsup/seuil) .le. crit(3)) then
        dp0 = xsup
        nbre = 1
        goto 50
    end if
!
    do i = 1, niter
!
        if ((abs(f/seuil)) .lt. crit(3)) then
            nbre = i
            goto 50
        end if
!
        dp0 = dp0-f/fp
!
        if (dp0 .ge. xsup .or. dp0 .le. xinf) dp0 = (xinf+xsup)/2
!
        f = dpvpeq(dp0, n, const, fonc1, fonc2, fonc3, fonc4)
        fp = dpvpdf(dp0, n, const, fonc1, fonc2, fonc3, fonc4)
!
!
        if (f .gt. zero) then
            signf = 1
        else
            signf = -1
        end if
!
        fi = dpvpeq(xinf, n, const, fonc1, fonc2, fonc3, fonc4)
        if (fi .gt. zero) then
            signfi = 1
        else
            signfi = -1
        end if
!
        if ((signf*signfi) .lt. zero) xsup = dp0
        if ((signf*signfi) .gt. zero) xinf = dp0
!
        if (abs(finf/seuil) .le. crit(3)) then
            dp0 = xinf
            nbre = 1
            goto 50
        else if (abs(fsup/seuil) .le. crit(3)) then
            dp0 = xsup
            nbre = 1
            goto 50
        end if
!
    end do
    retcom = 1
    goto 30
! =====================================================================
50  continue
    dp = dp0
30  continue
! =====================================================================
end subroutine
