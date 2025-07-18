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
subroutine chrgd(nbcmp, jcesd, jcesl, jcesv, imai, &
                 ipt, isp, type_gd, rc, p, &
                 permvec)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/tpsivp.h"
#include "asterfort/cesexi.h"
#include "asterfort/utmess.h"
#include "blas/dgemv.h"
!
    integer(kind=8), intent(in) :: nbcmp, jcesd, jcesl, jcesv, imai, ipt, isp
    real(kind=8), dimension(:, :), intent(in) :: p
    character(len=*), intent(in) :: type_gd
    character, intent(in) :: rc
!
    integer(kind=8), dimension(:), intent(in), optional :: permvec
! ----------------------------------------------------------------------
!
!     BUT : CHANGEMENT DE REPERE D'UNE GRANDEUR REELLE EN UN SOUS-POINT
!           D'UN CHAM_ELEM SIMPLE
! ----------------------------------------------------------------------
!     ARGUMENTS :
!     JCESD, JCESL IN I : ADRESSES DES OBJETS .CESD ET .CESL
!     IMAI     IN  I    : NUMERO DE LA MAILLE
!     IPT, ISP IN  I    : NUMERO DU POINT ET DU SOUS-POINT
!     NBCMP    IN  I    : NOMBRE DE COMPOSANTES A TRAITER
!     TYPE_GD  IN  K16  : TYPE DU CHAMP :'TENS_3D' 'TENS2D' 'VECT_3D' 'VECT_2D'
!     RC       IN  K1   : REEL OU COMPLEXE
!     P (3,3)  IN  R    : MATRICE DE PASSAGE
!     PERMVEC  IN  I    : VECTEUR DE PERMUTATION (POUR DIM 2 & CHANGEMENT CYLINDRIQUE)
! ---------------------------------------------------------------------
!
    integer(kind=8) :: ii, iad, kk
    real(kind=8), dimension(6) :: val1, val1r, val1i
    real(kind=8), dimension(6) :: val, valr, vali
    integer(kind=8), dimension(6) :: permvec_loc
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
!
!   Lecture des composantes du champ (vecteur ou tenseur) au sous-point courant
!       Si tenseur réel ou vecteur réel : lecture dans val1
!       Si tenseur complexe : lecture dans val1r et val1i
!   On utilise un vecteur de longueur fixe (6).
!       Les composantes non lues sont mises à 0.
!       On ne distingue pas dim=2 et dim=3
!
    val(:) = 0.d0
    valr(:) = 0.d0
    vali(:) = 0.d0
    val1(:) = 0.d0
    val1r(:) = 0.d0
    val1i(:) = 0.d0
!
!   Vecteur (optionnel) permettant de permuter les composantes (du tenseur ou du vecteur)
!   après application du changement de repère.  Cette possibilité est utilisée pour le
!   changement de repère vers un repère cylindrique dans les éléments de milieu continu.
    if (present(permvec)) then
        permvec_loc(:) = permvec(:)
    else
        permvec_loc(:) = (/(ii, ii=1, 6)/)
    end if
!
    do ii = 1, nbcmp
        call cesexi('C', jcesd, jcesl, imai, ipt, &
                    isp, ii, iad)
        if (iad .gt. 0) then
            select case (rc)
            case ('R')
                val1(ii) = zr(jcesv-1+iad)
            case ('C')
                val1r(ii) = dreal(zc(jcesv-1+iad))
                val1i(ii) = dimag(zc(jcesv-1+iad))
            case default
                ASSERT(.false.)
            end select
        end if
    end do
!   Changement de base
    select case (type_gd(1:7))
    case ('TENS_3D', 'TENS_2D')
!           Sigma <- P^T Sigma P
        select case (rc)
        case ('R')
            val(:) = val1(:)
            call tpsivp(p, val)
        case ('C')
            valr(:) = val1r(:)
            vali(:) = val1i(:)
            call tpsivp(p, valr)
            call tpsivp(p, vali)
        end select
    case ('VECT_3D', 'VECT_2D')
!           val = P^T val1
        select case (rc)
        case ('R')
            b_lda = to_blas_int(3)
            b_m = to_blas_int(3)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dgemv(trans='T', m=b_m, n=b_n, alpha=1.d0, a=p, &
                       lda=b_lda, x=val1, incx=b_incx, beta=0.d0, y=val, &
                       incy=b_incy)
        case ('C')
            b_lda = to_blas_int(3)
            b_m = to_blas_int(3)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dgemv(trans='T', m=b_m, n=b_n, alpha=1.d0, a=p, &
                       lda=b_lda, x=val1r, incx=b_incx, beta=0.d0, y=valr, &
                       incy=b_incy)
            b_lda = to_blas_int(3)
            b_m = to_blas_int(3)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dgemv(trans='T', m=b_m, n=b_n, alpha=1.d0, a=p, &
                       lda=b_lda, x=val1i, incx=b_incx, beta=0.d0, y=vali, &
                       incy=b_incy)
        end select
    case ('1D_GENE')
!           val = P val1 sur 2 blocs de 3 composantes  X,Y,Z  RX,RY,RZ
        select case (rc)
        case ('R')
            do ii = 1, 3
                do kk = 1, 3
                    val(ii) = val(ii)+p(ii, kk)*val1(kk)
                    val(ii+3) = val(ii+3)+p(ii, kk)*val1(kk+3)
                end do
            end do
        case ('C')
            call utmess('F', 'ALGORITH2_31')
        end select
    case default
        ASSERT(.false.)
    end select
!   Copie des composantes modifiées dans le champ
    do ii = 1, nbcmp
        call cesexi('C', jcesd, jcesl, imai, ipt, &
                    isp, ii, iad)
        if (iad .gt. 0) then
            select case (rc)
            case ('R')
                zr(jcesv-1+iad) = val(permvec_loc(ii))
            case ('C')
                zc(jcesv-1+iad) = dcmplx(valr(permvec_loc(ii)), vali(permvec_loc(ii)))
            end select
        end if
    end do
!
end subroutine chrgd
