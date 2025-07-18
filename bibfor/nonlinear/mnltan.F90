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
subroutine mnltan(lcal, imat, numdrv, matdrv, xcdl, &
                  parcho, adime, xvect, ninc, nd, &
                  nchoc, h, hf, xtang)
    implicit none
!
!
!     MODE_NON_LINE CALCUL D'UN VECTEUR TANGENT
!     -    -   -                        ---
! ----------------------------------------------------------------------
!
! CALCUL LE VECTEUR TANGENT AU POINT XVECT
! ----------------------------------------------------------------------
! IN  LCAL   : L     : SI .TRUE. ALORS ON RECALCUL LA MATRICE
!                                SINON ON REMPLACE SEULEMENT LA
!                                           DERNIERE LIGNE DE LA MATRICE
! IN  IMAT   : I(2) : DESCRIPTEUR DES MATRICES :
!                       - IMAT(1) => MATRICE DE RAIDEUR
!                       - IMAT(2) => MATRICE DE MASSE
! IN  NUMDRV : K14  : NUME_DDL_GENE DE LA MATRICE JACOBIENNE
! IN  MATDRV : K19  : NOM DE  LA MATRICE JACOBIENNE
! IN  XCDL   : K14  : INDICE DES CONDITIONS AUX LIMITES
! IN  PARCHO : K14  : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN  ADIME  : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN  XVECT  : K14  : NOM DU VECTEUR SOLUTION
! IN  NINC   : I    : NOMBRE D INCONNUES DU SYSTEME
! IN  ND     : I    : NOMBRE DE DEGRES DE LIBERTE ACTIFS
! IN  NCHOC  : I    : NOMBRE DE CONTACTEURS
! IN  H      : I    : NOMBRE D'HARMONIQUES POUR LE DEPLACEMENT
! IN  HF     : I    : NOMBRE D'HARMONIQUES POUR LA FORCE
! OUT XTANG  : K14  : NOM DU VECTEUR TANGENT AU POINT SOLUTION
! ----------------------------------------------------------------------
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getran.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnldrv.h"
#include "asterfort/resoud.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    aster_logical :: lcal
    integer(kind=8) :: imat(2), ninc, nd, nchoc, h, hf
    character(len=14) :: numdrv, xcdl, parcho, adime, xvect, xtang
    character(len=19) :: matdrv, solveu
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    integer(kind=8) :: i, itang, iret, ib, ivplu
    real(kind=8) :: norme
    complex(kind=8) :: cbid
    blas_int :: b_incx, b_incy, b_n
    cbid = dcmplx(0.d0, 0.d0)
! ----------------------------------------------------------------------
!
    call jemarq()
    solveu = '&&OP0061.SOLVEUR'
! ----------------------------------------------------------------------
! --- CREATION VECTEURS TEMPORAIRES
! ----------------------------------------------------------------------
    call wkvect('&&mnltan.b', 'V V R', ninc, ib)
    call wkvect('&&mnltan.vecplu', 'V V R', ninc, ivplu)
! ----------------------------------------------------------------------
! --- CREATION D'UN VECTEUR ALEATOIRE (A AJOUTER A LA DERNIERE LIGNE
! ---                                          DE LA MATRICE JACOBIENNE)
! ----------------------------------------------------------------------
    do i = 1, ninc
        call getran(zr(ivplu-1+i))
    end do
! ----------------------------------------------------------------------
! --- RECUPERATION DU VECTEUR TANGENT
! ----------------------------------------------------------------------
    call jeveuo(xtang, 'E', itang)
! ----------------------------------------------------------------------
! --- CALCUL (OU RECUPERATION) DE LA MATRICE JACOBIENNE
! ----------------------------------------------------------------------
    call mnldrv(lcal, imat, numdrv, matdrv, xcdl, &
                parcho, adime, xvect, zr(ivplu), ninc, &
                nd, nchoc, h, hf)
! ----------------------------------------------------------------------
! --- ON CREE UN VECTEUR [0 ... 0 1]
! ----------------------------------------------------------------------
    zr(ib-1+ninc) = 1.d0
! ----------------------------------------------------------------------
! --- ON RESOUD TANGENTE=DRDV\[0 ... 0 1]
! ----------------------------------------------------------------------
    call resoud(matdrv, ' ', solveu, ' ', 1, &
                ' ', ' ', 'v', zr(ib), [cbid], &
                ' ', .false._1, 0, iret)
    b_n = to_blas_int(ninc)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(ib), b_incx, zr(itang), b_incy)
! ----------------------------------------------------------------------
! --- ON NORMALISE LE VECTEUR TANGENT
! ----------------------------------------------------------------------
    b_n = to_blas_int(ninc)
    b_incx = to_blas_int(1)
    norme = dnrm2(b_n, zr(itang), b_incx)
    b_n = to_blas_int(ninc)
    b_incx = to_blas_int(1)
    call dscal(b_n, -1.d0/norme, zr(itang), b_incx)
! ----------------------------------------------------------------------
! --- ON DETRUIT LE VECTEUR TEMPORAIRE
! ----------------------------------------------------------------------
    call jedetr('&&mnltan.b')
    call jedetr('&&mnltan.vecplu')
!
    call jedema()
!
end subroutine
