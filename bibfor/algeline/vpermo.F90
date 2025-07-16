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
subroutine vpermo(lmasse, lraide, nbprop, vecp, valp, &
                  excl, omecor, ernorm)
    implicit none
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mrmult.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dscal.h"
    integer(kind=8) :: lmasse, lraide, nbprop, excl(*)
    real(kind=8) :: vecp(*), valp(*), omecor, ernorm(*)
!     CALCUL DE LA NORME D'ERREUR MODALE
!     ( IE NORME D'ERREUR SUR LES VALEURS ET VECTEURS PROPRES.)
!     ------------------------------------------------------------------
!     PROBLEME GENERALISE :  (LRAIDE)*VECP = VALP *(LMASSE)*VECP
!
!                   !! LRAIDE * VECP  - VALP * LMASSE * VECP !!
!       ERNORM   =     -------------------------------------
!                           !! LRAIDE * VECP !!
!     ------------------------------------------------------------------
!     REFERENCE: BATHE ET WILSON
!     ------------------------------------------------------------------
! IN  LMASSE : IS : DESCRIPTEUR MATRICE DE "MASSE"
! IN  LRAIDE : IS : DESCRIPTEUR MATRICE DE "RAIDEUR"
! IN  NBPROP : IS : NOMBRE DE VALEURS ET DE VECTEURS PROPRES
! IN  VECP   : R8 : TABLEAU DES VECTEURS PROPRES
! IN  VALP   : R8 : TABLEAU DES VALEURS PROPRES
! IN  EXCL   : IS : TABLEAU DES NON-EXCLUS
! IN  FCORIG : R8 : FREQUENCE DE CORPS RIGIDE
! OUT ERNORM : R8 : TABLEAU DES NORMES D'ERREUR
!     ------------------------------------------------------------------
    real(kind=8) :: anorm1, anorm2
!
!
!     --- SEUIL EN PULSATION POUR LES MODES DE CORPS RIGIDE ---
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iaux1, iaux2, j, neq, ivec
    real(kind=8) :: xseuil, rmin, raux
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    call jemarq()
    xseuil = omecor
!
!     ------------------------------------------------------------------
!     ---------------------- DONNEES SUR LES MATRICES ------------------
!     ------------------------------------------------------------------
    neq = zi(lmasse+2)
!     ------------------------------------------------------------------
!     -------------- ALLOCATION DES ZONES DE TRAVAIL -------------------
!     ------------------------------------------------------------------
    call wkvect('&&VPERMO.TAMPON.PROV_1', 'V V R', neq, iaux1)
    call wkvect('&&VPERMO.TAMPON.PROV_2', 'V V R', neq, iaux2)
!     ------------------------------------------------------------------
!     ---------------------- CALCUL DES NORMES D'ERREUR ----------------
!     ------------------------------------------------------------------
    rmin = 100.d0*r8miem()
!
!        --- NON PRISE EN COMPTE DES DDLS EXCLUS
    do i = 1, neq
        raux = excl(i)
        b_n = to_blas_int(nbprop)
        b_incx = to_blas_int(neq)
        call dscal(b_n, raux, vecp(i), b_incx)
    end do
!
    do i = 1, nbprop
        ivec = (i-1)*neq+1
        call mrmult('ZERO', lraide, vecp(ivec), zr(iaux1), 1, &
                    .false._1)
        call mrmult('ZERO', lmasse, vecp(ivec), zr(iaux2), 1, &
                    .false._1)
        anorm1 = 0.d0
        do j = 1, neq
            raux = zr(iaux1+j-1)
            anorm1 = anorm1+raux*raux*excl(j)
        end do
        raux = -valp(i)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, raux, zr(iaux2), b_incx, zr(iaux1), &
                   b_incy)
        anorm2 = 0.d0
        do j = 1, neq
            raux = zr(iaux1+j-1)
            anorm2 = anorm2+raux*raux*excl(j)
        end do
!
!
        if (abs(valp(i)) .gt. xseuil) then
            if (anorm1 .ge. rmin) then
                ernorm(i) = sqrt(anorm2/anorm1)
            else
                ernorm(i) = 1.d+70
            end if
        else
            ernorm(i) = abs(valp(i))*sqrt(anorm2)
        end if
    end do
!
!     ----------------------------------------------------------------
!     -------------- DESALLOCATION DES ZONES DE TRAVAIL --------------
!     ----------------------------------------------------------------
!
    call jedetr('&&VPERMO.TAMPON.PROV_1')
    call jedetr('&&VPERMO.TAMPON.PROV_2')
!
    call jedema()
end subroutine
