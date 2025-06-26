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
subroutine lgldom(nbmat, mater, yf, fiter)
!
    implicit none
#include "jeveux.h"
#include "asterfort/cos3t.h"
#include "asterfort/domrev.h"
#include "asterfort/gdev.h"
#include "asterfort/hlode.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/ucritp.h"
#include "asterfort/varecr.h"
#include "asterfort/wkvect.h"
#include "blas/ddot.h"
    integer(kind=8) :: nbmat
    real(kind=8) :: mater(nbmat, 2), yf(10), fiter
! --- BUT : VALEUR DE F POUR LE CONVEXE ELASTO-PLASTIQUE ----------
! =================================================================
! IN  : NBMAT : NOMBRE DE PARAMETRES MATERIAU ---------------------
! --- : MATER : PARAMETRES MATERIAU -------------------------------
! --- : NR    : NOMBRE DE CONDITIONS NON LINEAIRES ----------------
! --- : YF    : INCREMENTS A L'INSTANT COURANT --------------------
! OUT : FITER : VALEUR DE F(S) A L'INSTANT COURANT ----------------
! =================================================================
! =================================================================
    integer(kind=8) :: ndt, ndi, jpara
    real(kind=8) :: sn(6), i1n, gampn, snii, lgleps, gamcjs, pref
    real(kind=8) :: rcos3t, rhlode, rgdev, sigc
    real(kind=8) :: rucpla
    character(len=16) :: parecr
    blas_int :: b_incx, b_incy, b_n
! =================================================================
! --- INITIALISATION DE PARAMETRES --------------------------------
! =================================================================
    parameter(lgleps=1.0d-8)
! =================================================================
    common/tdim/ndt, ndi
! =================================================================
    call jemarq()
! =================================================================
! --- DEFINITIONS -------------------------------------------------
! =================================================================
    parecr = '&&LGLDOM.PARECR'
    call wkvect(parecr, 'V V R', 5, jpara)
! =================================================================
! --- RECUPERATION DE DONNEES -------------------------------------
! =================================================================
    sigc = mater(9, 2)
    gamcjs = mater(12, 2)
    pref = mater(15, 2)
    sn(1:ndt) = yf(1:ndt)
    i1n = yf(ndt+1)
    gampn = yf(ndt+2)
! =================================================================
! --- CALCUL DE G(S) ----------------------------------------------
! =================================================================
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    snii = ddot(b_n, sn, b_incx, sn, b_incy)
    snii = sqrt(snii)
    rcos3t = cos3t(sn, pref, lgleps)
    rhlode = hlode(gamcjs, rcos3t)
    rgdev = gdev(snii, rhlode)
! =================================================================
! --- CALCUL DE U(SIG, GAMP) --------------------------------------
! =================================================================
    call varecr(gampn, nbmat, mater, zr(jpara))
! =================================================================
! --- SI LE CRITERE PLASTIQUE EST NEGATIF ON REDECOUPE ------------
! =================================================================
    rucpla = ucritp(nbmat, mater, zr(jpara), rgdev, i1n)
    fiter = domrev(gamcjs, sigc, zr(jpara), rgdev, rucpla)
! =================================================================
! --- DESTRUCTION DES VECTEURS INUTILES ---------------------------
! =================================================================
    call jedetr(parecr)
! =================================================================
    call jedema()
! =================================================================
end subroutine
