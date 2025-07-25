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
subroutine lchobr(toler, itmax, mod, nbmat, materf, &
                  nr, nvi, depsm, sigm, vim, &
                  seuil, vp, vecp, icomp, sigp, &
                  vip, irtet)
    implicit none
#include "asterc/r8dgrd.h"
#include "asterfort/calcvh.h"
#include "asterfort/hbcalc.h"
#include "asterfort/hbcrel.h"
#include "asterfort/hbderi.h"
#include "asterfort/hbmajs.h"
#include "asterfort/hbvaec.h"
#include "asterfort/lcdevi.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
    integer(kind=8) :: itmax, nbmat, nr, nvi, icomp, irtet
    real(kind=8) :: toler, materf(nbmat, 2), depsm(6), sigm(6)
    real(kind=8) :: vim(*), sigp(6), vip(*), seuil, vp(3), vecp(3, 3)
    character(len=8) :: mod
! ======================================================================
! --- LOI DE COMPORTEMENT DE TYPE HOEK BROWN MODIFIE -------------------
! --- *ELASTICITE ISOTROPE ---------------------------------------------
! --- *CRITERE DE PLASTICITE DE HEOK BROWN -----------------------------
! --- *ECOULEMENT PLASTIQUE DE DRUCKER PRAGER --------------------------
! --- CALCUL DU TENSEUR DES CONTRAINTES DANS LE CAS PLASTIQUE ----------
! ======================================================================
! IN  : TOLER  : VALEUR DE LA TOLERANCE DE CONVERGENCE -----------------
! --- :        : (RESI_INTE) --------------------------------------
! --- : ITMAX  : NOMBRE D'ITERATIONS MAXIMUM A CONVERGENCE -------------
! --- :        : (ITER_INTE_MAXI) --------------------------------------
! --- : MOD    : TYPE DE MODELISATION ----------------------------------
! --- : NBMAT  : NOMBRE DE DONNEES MATERIAU ----------------------------
! --- : MATERF : DONNEES MATERIAU --------------------------------------
! --- : NR     : NOMBRE DE RELATIONS NON LINEAIRES ---------------------
! --- : NVI    : NOMBRE DE VARIABLES INTERNES --------------------------
! --- : DEPSM  : ACCROISSEMENTS DE DEFORMATIONS A L'ITERATION COURANTE -
! --- : SIGM   : CONTRAINTES A L'INSTANT PRECEDENT ---------------------
! --- : VIM    : VARIABLES INTERNES A L'INSTANT PRECEDENT --------------
! --- : SEUIL  : VARIABLE SEUIL ELASTIQUE ------------------------------
! --- : VP     : VALEURS PROPRES DU DEVIATEUR ELASTIQUE ----------------
! --- : VECP   : VECTEURS PROPRES DU DEVIATEUR ELASTIQUE ---------------
! --- : ICOMP  : COMPTEUR POUR LE REDECOUPAGE DU PAS DE TEMPS ----------
! OUT : SIGP   : CONTRAINTES A L'INSTANT COURANT -----------------------
! --- : VIP    : VARIABLES INTERNES A L'INSTANT COURANT ----------------
! --- : IRTET  : CONTROLE DU REDECOUPAGE DU PAS DE TEMPS ---------------
! ======================================================================
!     INFO   MATERF   COEFFICIENTS MATERIAUX A T+DT
!                     MATERF(*,1) = CARACTERISTIQUES ELASTIQUES
!                     MATERF(*,2) = CARACTERISTIQUES PLASTIQUES
!            NDT      NOMBRE DE COMPOSANTES TOTALES DES TENSEURS
!            NDI      NOMBRE DE COMPOSANTES DIRECTES DES TENSEURS
!            NVI      NOMBRE DE VARIABLES INTERNES
! ======================================================================
    integer(kind=8) :: ndt, ndi, iteri, iter
    real(kind=8) :: gm, etam, etap, aux
    real(kind=8) :: se(6)
    real(kind=8) :: seq, i1e, seuil2, plas, dg, sigeqe
    real(kind=8) :: deux, trois
    real(kind=8) :: incrg, gnp, dgnp, etanp, vh, vg
    real(kind=8) :: parame(4), derive(5), pi, fmoins
    blas_int :: b_incx, b_incy, b_n
! ======================================================================
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
! --- INITIALISATION DES PARAMETRES DE CONVERGENCE ---------------------
! ======================================================================
    pi = r8dgrd()
! ======================================================================
! --- INITIALISATION ---------------------------------------------------
! ======================================================================
    gm = vim(1)
    if (gm .lt. 0.0d0) then
        call utmess('F', 'ALGORITH3_88')
    end if
    irtet = 0
    iteri = 0
! =====================================================================
! --- CALCUL DES PARAMETRES D ECROUISSAGE -----------------------------
! =====================================================================
    call hbvaec(gm, nbmat, materf, parame)
    etam = deux*sin(parame(4)*pi)/(trois+sin(parame(4)*pi))
! =====================================================================
    call lcdevi(sigp, se)
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    seq = ddot(b_n, se, b_incx, se, b_incy)
    sigeqe = sqrt(trois*seq/deux)
    i1e = trace(ndi, sigp)
! ======================================================================
    dg = 0.0d0
! ======================================================================
! --- CALCUL DE DELTA GAMMA --------------------------------------------
! ======================================================================
    plas = 1.0d0
    call hbcrel(vp, gm, dg, nbmat, materf, &
                sigeqe, i1e, etam, parame, seuil2)
    fmoins = seuil2
! ======================================================================
! --------- CALCUL DE L INCREMENT DE GAMMA PAR METHODE DE NEWTON -------
! ======================================================================
! --------- INITIALISATION DES VARIABLES -------------------------------
! ======================================================================
    iter = 0
    incrg = 0.d0
    dgnp = dg
    gnp = gm
    etanp = etam
    call calcvh(nbmat, materf, etanp, vp, sigeqe, &
                vh, vg)
    call hbderi(gnp, nbmat, materf, vg, etanp, &
                parame, derive)
! ======================================================================
! --------- PREMIERE ITERATION -----------------------------------------
! ======================================================================
    call hbcalc(seuil2, gnp, dgnp, nbmat, materf, &
                i1e, sigeqe, vp, etanp, vh, &
                vg, parame, derive, incrg)
2   continue
    gnp = gnp+incrg
    dgnp = dgnp+incrg
! ======================================================================
! -- ON OBTIENT DGAMMA_P NEGATIF : ON ESSAIE DE DECOUPER LE PAS DE TEMPS
! ======================================================================
    if (dgnp .lt. 0.d0) then
        if ((icomp .eq. 0) .or. (icomp .eq. 1)) then
            call utmess('I', 'ALGORITH4_57')
            iteri = 1
            goto 100
        else
            call utmess('I', 'ALGORITH4_60')
            goto 100
        end if
    end if
    call hbvaec(gnp, nbmat, materf, parame)
    etanp = deux*sin(parame(4)*pi)/(trois+sin(parame(4)*pi))
    call hbcrel(vp, gnp, dgnp, nbmat, materf, &
                sigeqe, i1e, etanp, parame, seuil2)
! ======================================================================
! ---------- IL Y A CONVERGENCE ----------------------------------------
! ======================================================================
    if ((abs(seuil2) .lt. toler) .or. (abs(seuil2/fmoins) .lt. toler)) then
! ======================================================================
! --------- ON DETECTE LES SOLUTIONS NON ADMISSIBLES -------------------
! ======================================================================
        aux = sigeqe*(etanp+1.0d0)/(trois*materf(4, 1))
        if (dgnp .gt. aux) then
            call utmess('I', 'ALGORITH4_58')
            iteri = 1
            goto 100
        end if
        dg = dgnp
        iteri = 0
! ======================================================================
! --------- LE NOMBRE MAX D ITERATIONS N A PAS ETE ATTEINT -------------
! ======================================================================
    else if (iter .lt. itmax) then
        iter = iter+1
        iteri = 0
        call calcvh(nbmat, materf, etanp, vp, sigeqe, &
                    vh, vg)
        call hbderi(gnp, nbmat, materf, vg, etanp, &
                    parame, derive)
        call hbcalc(seuil2, gnp, dgnp, nbmat, materf, &
                    i1e, sigeqe, vp, etanp, vh, &
                    vg, parame, derive, incrg)
        goto 2
! ======================================================================
! --------- LE NOMBRE MAX D ITERATIONS A ETE ATTEINT -------------------
! ======================================================================
    else
! ======================================================================
! --------- ON ESSAIE DE DECOUPER LE PAS DE TEMPS ----------------------
! ======================================================================
        if ((icomp .eq. 0) .or. (icomp .eq. 1)) then
            call utmess('I', 'ALGORITH4_59')
            iteri = 1
            goto 100
        else
            call utmess('F', 'ALGORITH4_61')
        end if
    end if
100 continue
    if (iteri .eq. 1) goto 1
! ======================================================================
    etap = etanp
    call hbmajs(dg, nbmat, materf, se, i1e, &
                sigeqe, etap, sigp)
    vip(1) = vim(1)+dg
    vip(2) = vim(2)+trois*etap*dg/(etap+1.0d0)
    vip(3) = plas
! ======================================================================
    irtet = 0
    goto 999
1   continue
    irtet = 1
999 continue
! ======================================================================
end subroutine
