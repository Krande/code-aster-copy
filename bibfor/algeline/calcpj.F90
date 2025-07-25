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
subroutine calcpj(nbmat, mater, gamp, evp, sigd, &
                  sige, epssig, invare, gamps, evps, &
                  invars, b)
!
    implicit none
#include "jeveux.h"
#include "asterfort/bprime.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/lcdevi.h"
#include "asterfort/trace.h"
#include "asterfort/varecr.h"
#include "asterfort/wkvect.h"
#include "blas/ddot.h"
    integer(kind=8) :: nbmat
    real(kind=8) :: mater(nbmat, 2), gamp, evp, sigd(6), sige(6), epssig
    real(kind=8) :: invare, gamps, invars, evps, b
! --- BUT : CALCUL DE LA PROJECTION AU SOMMET --------------------------
! ======================================================================
! IN  : NDT    : NOMBRE DE COMPOSANTES TOTALES DU TENSEUR --------------
! --- : NDI    : NOMBRE DE COMPOSANTES DIAGONALES DU TENSEUR -----------
! --- : NBMAT  : NOMBRE DE PARAMETRES MATERIAU -------------------------
! --- : MATER  : PARAMETRES MATERIAU -----------------------------------
! --- : GAMP   : DEFORMATION DEVIATOIRE PLASTIQUE CUMULEE --------------
! --- : EVP    : DEFORMATION VOLUMIQUE PLASTIQUE CUMULEE ---------------
! --- : SIIE   : NORME DU TENSEUR --------------------------------------
! --- : EPSSIG : EPSILON -----------------------------------------------
! --- : INVARE : PREMIER INVARIANT DU TENSEUR DES CONTRAINTES ELASTIQUE-
! OUT : GAMPS  : DEFORMATION DEVIATOIRE PLASTIQUE CUMULEE AU SOMMET ----
! --- : INVARS : PREMIER INVARIANT DU TENSEUR DES CONTRAINTES AU SOMMET-
! --- : EVPS   : DEFORMATION VOLUMIQUE PLASTIQUE CUMULEE AU SOMMET -----
! --- : B      : PARAMETRE CONTROLANT LE COMPORTEMENT VOLUMIQUE --------
! ------------ : DU MATERIAU -------------------------------------------
! ======================================================================
! ======================================================================
    integer(kind=8) :: jpara, jpara2, ndt, ndi
    real(kind=8) :: mu, sigc, sig(6), sd(6), sgamp, mgamp
    real(kind=8) :: zero, deux, trois, se(6)
    real(kind=8) :: sigii, siie, invar, gamult, k
    character(len=16) :: parecr, parec2
    blas_int :: b_incx, b_incy, b_n
! ======================================================================
! --- INITIALISATION DE PARAMETRE --------------------------------------
! ======================================================================
    parameter(zero=0.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
    call jemarq()
! ======================================================================
! --- RECUPERATION DE PARAMETRES DU MODELE -----------------------------
! ======================================================================
    mu = mater(4, 1)
    k = mater(5, 1)
    gamult = mater(1, 2)
    sigc = mater(9, 2)
! ======================================================================
! --- INITIALISATION DE PARAMETRES -------------------------------------
! ======================================================================
    parecr = '&&CALCPJ.PARECR'
    parec2 = '&&CALCPJ.PAREC2'
    call wkvect(parecr, 'V V R', 5, jpara)
    call wkvect(parec2, 'V V R', 5, jpara2)
! ======================================================================
! --- CALCUL DES PROJECTIONS AU SOMMET ---------------------------------
! ======================================================================
    call lcdevi(sige, se)
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    siie = ddot(b_n, se, b_incx, se, b_incy)
    siie = sqrt(siie)
    gamps = gamp+sqrt(deux/trois)*siie/(deux*mu)
    call varecr(gamps, nbmat, mater, zr(jpara))
    sgamp = zr(jpara-1+1)
    mgamp = zr(jpara-1+4)
    invars = trois*sigc*sgamp/mgamp
    evps = evp+(invare-invars)/(trois*k)
! ======================================================================
! --- CALCUL DU PARAMETRE DE DILATANCE B A L'INSTANT MOINS -------------
! ======================================================================
! --- CAS OU GAMP > GAMULT(1-EPS) --------------------------------------
! ======================================================================
    if (gamp .gt. gamult) then
        b = zero
    else
! ======================================================================
! --- CAS OU GAMP <= GAMULT(1-EPS) -------------------------------------
! ======================================================================
        b_n = to_blas_int(ndt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        sigii = ddot(b_n, sigd, b_incx, sigd, b_incy)
        if (sigii .lt. epssig) then
            sig(1:ndt) = sige(1:ndt)
        else
            sig(1:ndt) = sigd(1:ndt)
        end if
        call lcdevi(sig, sd)
        invar = trace(ndi, sig)
        call varecr(gamp, nbmat, mater, zr(jpara2))
        b = bprime(nbmat, mater, zr(jpara2), invar, sd, epssig)
    end if
! ======================================================================
! --- DESTRUCTION DES VECTEURS INUTILES --------------------------------
! ======================================================================
    call jedetr(parecr)
    call jedetr(parec2)
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
