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
subroutine mnlru(imat, xcdl, parcho, adime, xvect, &
                 ninc, nd, nchoc, h, hf, &
                 xru)
    implicit none
!
!
!     MODE_NON_LINE CALCUL DE R(U)
!     -    -   -              - -
! ----------------------------------------------------------------------
!
! CALCUL R(U) = L(U) + Q(U,U)
! ----------------------------------------------------------------------
! IN  IMAT   : I(2) : DESCRIPTEUR DES MATRICES :
!                       - IMAT(1) => MATRICE DE RAIDEUR
!                       - IMAT(2) => MATRICE DE MASSE
! IN  XCDL   : K14  : INDICE DES CONDITIONS AUX LIMITES
! IN  PARCHO : K14  : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN  ADIME  : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN  XVECT  : K14  : NOM DU VECTEUR SOLUTION
! IN  NINC   : I    : NOMBRE D INCONNUES DU SYSTEME
! IN  ND     : I    : NOMBRE DE DEGRES DE LIBERTE
! IN  NCHOC  : I    : NOMBRE DE CONTACTEURS
! IN  H      : I    : NOMBRE D'HARMONIQUES POUR LE DEPLACEMENT
! IN  HF     : I    : NOMBRE D'HARMONIQUES POUR LA FORCE
! OUT XRU    : K14  : NOM DU VECTEUR DE SORTIE, R(XVECT)=XRU
! ----------------------------------------------------------------------
!
!
#include "jeveux.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnlcst.h"
#include "asterfort/mnline.h"
#include "asterfort/mnlqnl.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: imat(2), ninc, nd, nchoc, h, hf
    character(len=14) :: xcdl, parcho, adime, xvect, xru
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    integer(kind=8) :: iru, ivint
    character(len=14) :: xvint
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
!    call jxveri(' ', ' ')
! ----------------------------------------------------------------------
! --- RECUPERATION DU POINTEUR DE R(XVECT)
! ----------------------------------------------------------------------
    call jeveuo(xru, 'E', iru)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(iru), b_incx)
! ----------------------------------------------------------------------
! --- CREATION D'UN VECTEUR INTERMEDIAIRE
! ----------------------------------------------------------------------
    xvint = '&&MNLRU.INT'
    call wkvect(xvint, 'V V R', ninc-1, ivint)
! ----------------------------------------------------------------------
! --- CALCUL DE R(XVECT)
! ----------------------------------------------------------------------
! --- CALCUL DE L0
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(ivint), b_incx)
    call mnlcst(parcho, adime, ninc, nd, nchoc, &
                h, hf, xvint)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(ivint), b_incx, zr(iru), b_incy)
! --- CALCUL DE L(XVECT)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(ivint), b_incx)
    call mnline(imat, xcdl, parcho, adime, xvect, &
                ninc, nd, nchoc, h, hf, &
                xvint)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, zr(ivint), b_incx, zr(iru), &
               b_incy)
! --- CALCUL DE Q(XVECT,XVECT)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(ivint), b_incx)
    call mnlqnl(imat, xcdl, parcho, adime, xvect, &
                xvect, ninc, nd, nchoc, h, &
                hf, xvint)
! --- R(XVECT) = L0 + L(XVECT) + Q(XVECT,XVECT)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, zr(ivint), b_incx, zr(iru), &
               b_incy)
! ----------------------------------------------------------------------
! --- DESTRUCTION DU VECTEUR INTERMEDIAIRE
! ----------------------------------------------------------------------
    call jedetr(xvint)
!
    call jedema()
!
end subroutine
