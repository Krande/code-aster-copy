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
subroutine dktrig(nomte, xyzl, option, pgl, rig, &
                  ener, multic)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/bsthpl.h"
#include "asterfort/dktbf.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxtbm.h"
#include "asterfort/dxtloc.h"
#include "asterfort/dxtloe.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/r8inir.h"
#include "asterfort/utbtab.h"
#include "asterfort/utctab.h"
#include "asterfort/utpvgl.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    real(kind=8) :: xyzl(3, *), pgl(*), rig(*), ener(*)
    character(len=16) :: option, nomte
!
!     MATRICE DE RIGIDITE DE L'ELEMENT DE PLAQUE DKT
!     ------------------------------------------------------------------
!     IN  XYZL   : COORDONNEES LOCALES DES TROIS NOEUDS
!     IN  OPTION : OPTION RIGI_MECA OU EPOT_ELEM
!     IN  PGL    : MATRICE DE PASSAGE GLOBAL/LOCAL
!     OUT RIG    : MATRICE DE RIGIDITE
!     OUT ENER   : TERMES POUR ENER_POT (EPOT_ELEM)
!     ------------------------------------------------------------------
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: multic, i, jcoqu, jdepg
    real(kind=8) :: wgt, aire
    real(kind=8) :: dm(9), df(9), dmf(9), df2(9), dmf2(9), dc(4), dci(4)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bf(3, 9), bm(3, 6)
    real(kind=8) :: xab1(3, 9), depl(18)
    real(kind=8) :: flex(81), memb(36), mefl(54)
    real(kind=8) :: bsigth(24), enerth, ctor
    aster_logical :: coupmf, indith
    real(kind=8) :: qsi, eta, carat3(21), t2iu(4), t2ui(4), t1ve(9)
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
    enerth = 0.0d0
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
    call jevech('PCACOQU', 'L', jcoqu)
    ctor = zr(jcoqu+3)
!
!     ------ MISE A ZERO DES MATRICES : FLEX ET MEFL -------------------
    call r8inir(81, 0.d0, flex, 1)
    call r8inir(36, 0.d0, memb, 1)
    call r8inir(54, 0.d0, mefl, 1)
!
!     ----- CALCUL DES GRANDEURS GEOMETRIQUES SUR LE TRIANGLE ----------
    call gtria3(xyzl, carat3)
!
!     CALCUL DES MATRICES DE RIGIDITE DU MATERIAU EN FLEXION
!     MEMBRANE ET CISAILLEMENT INVERSEE
    call dxmate('RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, nno, pgl, &
                multic, coupmf, t2iu, t2ui, t1ve)
!     ------------------------------------------------------------------
!     CALCUL DE LA MATRICE DE RIGIDITE DE L'ELEMENT EN MEMBRANE
!     ------------------------------------------------------------------
!
!     ------ CALCUL DE LA MATRICE BM -----------------------------------
    call dxtbm(carat3(9), bm)
    aire = carat3(8)
!
!     ------ CALCUL DU PRODUIT BMT.DM.BM -------------------------------
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, dm, b_incx, dmf2, b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    call dscal(b_n, aire, dmf2, b_incx)
    call utbtab('ZERO', 3, 6, dmf2, bm, &
                xab1, memb)
!
!     ------------------------------------------------------------------
!     CALCUL DES MATRICES DE RIGIDITE DE L'ELEMENT EN FLEXION ET
!     COUPLAGE MEMBRANE/FLEXION
!     ------------------------------------------------------------------
    do i = 1, npg
        qsi = zr(icoopg-1+ndim*(i-1)+1)
        eta = zr(icoopg-1+ndim*(i-1)+2)
        wgt = zr(ipoids+i-1)*carat3(7)
!        ----- CALCUL DE LA MATRICE BF AU POINT QSI ETA ------------
        call dktbf(qsi, eta, carat3, bf)
!        ----- CALCUL DU PRODUIT BFT.DF.BF -------------------------
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, df, b_incx, df2, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, wgt, df2, b_incx)
        call utbtab('CUMU', 3, 9, df2, bf, &
                    xab1, flex)
        if (coupmf) then
!        ----- CALCUL DU PRODUIT BMT.DMF.BF ------------------------
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, dmf, b_incx, dmf2, b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            call dscal(b_n, wgt, dmf2, b_incx)
            call utctab('CUMU', 3, 9, 6, dmf2, &
                        bf, bm, xab1, mefl)
        end if
!
    end do
!
    if (option .eq. 'RIGI_MECA') then
        call dxtloc(flex, memb, mefl, ctor, rig)
!
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PDEPLAR', 'L', jdepg)
        call utpvgl(3, 6, pgl, zr(jdepg), depl)
        call dxtloe(flex, memb, mefl, ctor, coupmf, &
                    depl, ener)
        call bsthpl(nomte, bsigth, indith)
        if (indith) then
            do i = 1, 18
                enerth = enerth+depl(i)*bsigth(i)
            end do
            ener(1) = ener(1)-enerth
        end if
    end if
!
end subroutine
