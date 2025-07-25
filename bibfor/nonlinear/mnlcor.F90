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
subroutine mnlcor(imat, numdrv, matdrv, xcdl, parcho, &
                  adime, ninc, nd, nchoc, h, &
                  hf, itemax, epscor, xvect, cor, &
                  info)
    implicit none
!
!
!     MODE_NON_LINE CORRECTION D'UN POINT SOLUTION
!     -    -   -                        ---
! ----------------------------------------------------------------------
!
! CORRIGE LE POINT DONNE EN ARGUMENT A L'AIDE D'UN ALGORITHME DE NEWTON
! ----------------------------------------------------------------------
! IN      IMAT       : I(2) : DESCRIPTEUR DES MATRICES :
!                              - IMAT(1) => MATRICE DE RAIDEUR
!                              - IMAT(2) => MATRICE DE MASSE
! IN      NUMDRV : K14  : NUME_DDL_GENE DE LA MATRICE JACOBIENNE
! IN      MATDRV : K19  : NOM DE  LA MATRICE JACOBIENNE
! IN      XCDL   : K14  : INDICE DES CONDITIONS AUX LIMITES
! IN      PARCHO : K14  : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN      ADIME  : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN      NINC   : I    : NOMBRE D INCONNUES DU SYSTEME
! IN      ND     : I    : NOMBRE DE DEGRES DE LIBERTE ACTIFS
! IN      NCHOC  : I    : NOMBRE DE CONTACTEURS
! IN      H      : I    : NOMBRE D'HARMONIQUES POUR LE DEPLACEMENT
! IN      HF     : I    : NOMBRE D'HARMONIQUES POUR LA FORCE
! IN      ITEMAX : I    : NOMBRE MAXIMAL D'ITERATIONS DE NEWTON
! IN      EPSCOR : R8   : PRECISION POUR L'ALGORITHME DE NEWTON
! IN/OUT  XVECT  : K14  : NOM DU VECTEUR D'ENTREE/VECTEUR CORRIGE
! OUT     COR    : L    : TRUE SI LA CORRECTION A REUSSI
!                         FALSE SINON
! IN      INFO   : I    : /1 --> PAS D'AFFICHAGE DE LA NORME D'ERREUR
!                         /2 -->       AFFICHAGE DE LA NORME D'ERREUR
! ----------------------------------------------------------------------
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnldrv.h"
#include "asterfort/mnlru.h"
#include "asterfort/mnltan.h"
#include "asterfort/resoud.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dnrm2.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
!
    integer(kind=8) :: imat(2), ninc, nd, nchoc, h, hf, itemax, info
    character(len=14) :: numdrv, xcdl, parcho, adime, xvect
    character(len=19) :: matdrv, solveu
    real(kind=8) :: epscor
    aster_logical :: cor
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    character(len=14) :: xru, xtang, xtemp
    integer(kind=8) :: iru, itang, ivect, cptr, iret, itemp, inddl, ifres
    real(kind=8) :: eps, normr, normc
    complex(kind=8) :: cbid
    blas_int :: b_incx, b_incy, b_n
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
!
    ifres = iunifi('MESSAGE')
    call jemarq()
    solveu = '&&OP0061.SOLVEUR'
! ----------------------------------------------------------------------
! --- CREATION DE VECTEURS UTILES
! ----------------------------------------------------------------------
    xru = '&&MNLCOR.RU'
    xtang = '&&MNLCOR.TANG'
    xtemp = '&&MNLCOR.TEMP'
    call wkvect(xru, 'V V R', ninc, iru)
    call wkvect(xtang, 'V V R', ninc, itang)
    call wkvect(xtemp, 'V V R', ninc, itemp)
    eps = epscor
! ----------------------------------------------------------------------
! --- RECUPERATION DU VECTEUR A CORRIGER
! ----------------------------------------------------------------------
    call jeveuo(xvect, 'E', ivect)
! ----------------------------------------------------------------------
! --- HEURISTIQUE SUR LE CHOIX DE L'EPS
! ----------------------------------------------------------------------
    call jeveuo(parcho//'.NDDL', 'L', inddl)
! ----------------------------------------------------------------------
! --- INITIALISATION DE L'ALGORITHME DE NEWTON
! ----------------------------------------------------------------------
    cptr = 0
    call mnlru(imat, xcdl, parcho, adime, xvect, &
               ninc, nd, nchoc, h, hf, &
               xru)
    zr(iru-1+ninc) = 0.d0
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    normr = dnrm2(b_n, zr(iru), b_incx)
    b_n = to_blas_int(ninc)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(ivect), b_incx, zr(itemp), b_incy)
    normc = normr
900 format(' Norme erreur iteration Newton numero : ', i2, ' : ', 1pe12.5)
    if (info .eq. 2) then
        write (ifres, 900) cptr, normc
    end if
! ----------------------------------------------------------------------
! --- ALGORITHME DE NEWTON
! ----------------------------------------------------------------------
120 continue
    if (normc .ge. eps .and. cptr .lt. itemax) then
! ---   INCREMENTATION DU COMPTEUR
        cptr = cptr+1
! ---   CALCUL DU VECTEUR TANGENT
        if (cptr .eq. 1) then
            call mnltan(.true._1, imat, numdrv, matdrv, xcdl, &
                        parcho, adime, xtemp, ninc, nd, &
                        nchoc, h, hf, xtang)
        else
            call mnltan(.false._1, imat, numdrv, matdrv, xcdl, &
                        parcho, adime, xtemp, ninc, nd, &
                        nchoc, h, hf, xtang)
        end if
! ---   RECALCUL DE LA MATRICE JACOBIENNE
        call mnldrv(.true._1, imat, numdrv, matdrv, xcdl, &
                    parcho, adime, xtemp, zr(itang), ninc, &
                    nd, nchoc, h, hf)
! ---   ON RESOUD LE SYSTEME LINEAIRE (DRDV\XTANG)
        call resoud(matdrv, ' ', solveu, ' ', 1, &
                    ' ', ' ', 'V', zr(iru), [cbid], &
                    ' ', .false._1, 0, iret)
! ---   ON AJOUTE AU VECTEUR SOLUTION
        b_n = to_blas_int(ninc)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, zr(iru), b_incx, zr(itemp), &
                   b_incy)
! ---   ON CALCUL R(NOUVEAU VECTEUR SOLUTION)
        call mnlru(imat, xcdl, parcho, adime, xtemp, &
                   ninc, nd, nchoc, h, hf, &
                   xru)
        zr(iru-1+ninc) = 0.d0
! ---   ON CALCUL LA NORME DE R(NOUVEAU VECTEUR SOLUTION)
        b_n = to_blas_int(ninc-1)
        b_incx = to_blas_int(1)
        normc = dnrm2(b_n, zr(iru), b_incx)
        normr = normc
        b_n = to_blas_int(ninc)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(itemp), b_incx, zr(ivect), b_incy)
        if (info .eq. 2) then
            write (ifres, 900) cptr, normc
        end if
        goto 120
    end if
    if (normr .ge. eps) then
        cor = .false.
    else
        cor = .true.
    end if
! ----------------------------------------------------------------------
! --- DESTRUCTION DES VECTEURS UTILES
! ----------------------------------------------------------------------
    call jedetr(xru)
    call jedetr(xtang)
    call jedetr(xtemp)
!
    call jedema()
!
end subroutine
