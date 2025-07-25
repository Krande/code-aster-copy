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
subroutine vppgen(lmasse, lamor, lraide, masseg, amorg, &
                  raideg, vect, neq, nbvect, iddl)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mrmult.h"
#include "asterfort/wkvect.h"
#include "blas/ddot.h"
    integer(kind=8) :: lmasse, lamor, lraide, neq, nbvect, iddl(*)
    real(kind=8) :: masseg(*), amorg(*), raideg(*), vect(neq, *)
!     CALCUL DES PARAMETRES MODAUX :
!            MASSE, AMORTISSEMENT ET RAIDEUR GENERALISES
!     ------------------------------------------------------------------
! IN  LMASSE : IS : DESCRIPTEUR NORMALISE DE LA MATRICE DE MASSE
!            = 0  ON NE CALCULE PAS LA MASSE GENERALISEE
! IN  LAMOR  : IS : DESCRIPTEUR NORMALISE DE LA MATRICE D'AMORTISSEMENT
!            = 0  ON NE CALCULE PAS L'AMORTISSEMENT GENERALISE
! IN  LRAIDE : IS : DESCRIPTEUR NORMALISE DE LA MATRICE DE RIGIDITE
!            = 0  ON NE CALCULE PAS LA RIGIDITE GENERALISEE
!     ------------------------------------------------------------------
!     REMARQUE : ON FAIT LES CALCULS VECTEURS APRES VECTEURS
!              : C'EST PLUS LONG MAIS PAS DE PB DE TAILLE MEMOIRE
!     ------------------------------------------------------------------
!
!
    real(kind=8) :: rzero
    character(len=24) :: vecaux, vecau1
    integer(kind=8) :: ieq, ivect, laux, laux1
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    data vecaux/'&&VPPGEN.VECTEUR.AUX0'/
    data vecau1/'&&VPPGEN.VECTEUR.AUX1'/
!     ------------------------------------------------------------------
    call jemarq()
    call wkvect(vecaux, 'V V R', neq, laux)
    call wkvect(vecau1, 'V V R', neq, laux1)
    laux = laux-1
    laux1 = laux1-1
!
    rzero = 0.d0
!     ------------------------------------------------------------------
!     ----------------- CALCUL DE LA MASSE GENERALISEE -----------------
!     ------------------------------------------------------------------
    if (lmasse .ne. 0) then
        do ivect = 1, nbvect
            call mrmult('ZERO', lmasse, vect(1, ivect), zr(laux+1), 1, &
                        .false._1)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            masseg(ivect) = ddot(b_n, vect(1, ivect), b_incx, zr(laux+1), b_incy)
        end do
    end if
!     ------------------------------------------------------------------
!     --------------- CALCUL DE L'AMORTISSEMENT GENERALISE -------------
!     ------------------------------------------------------------------
    if (lamor .ne. 0) then
        do ivect = 1, nbvect
            call mrmult('ZERO', lamor, vect(1, ivect), zr(laux+1), 1, &
                        .false._1)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            amorg(ivect) = ddot(b_n, vect(1, ivect), b_incx, zr(laux+1), b_incy)
        end do
    else
        amorg(1:nbvect) = rzero
    end if
!     ------------------------------------------------------------------
!     ---------------- CALCUL DE LA RAIDEUR GENERALISEE ----------------
!     ------------------------------------------------------------------
    if (lraide .ne. 0) then
        do ivect = 1, nbvect
            do ieq = 1, neq
                zr(laux1+ieq) = vect(ieq, ivect)*iddl(ieq)
            end do
            call mrmult('ZERO', lraide, zr(laux1+1), zr(laux+1), 1, &
                        .false._1)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            raideg(ivect) = ddot(b_n, zr(laux+1), b_incx, zr(laux1+1), b_incy)
        end do
    end if
!     ------------------------------------------------------------------
    call jedetr(vecaux)
    call jedetr(vecau1)
!     ------------------------------------------------------------------
    call jedema()
end subroutine
