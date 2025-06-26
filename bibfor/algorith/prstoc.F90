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
subroutine prstoc(vecsol, vestoc, j, k, iad, &
                  nbvale, nbrefe)
    implicit none
!
!
!         ROUTINE STOCKANT LE VECTEUR PRESSION
!         ISSUE D' UNE RESOLUTION DE LAPLACE
! IN : VECSOL : VECTEUR SOLUTION K19
! IN : J : INDICE DE BOUCLE
! IN : IAD : ADRESSE DU VECTEUR DES NOMS DES CHAMNOS STOCKES
! IN : NBVALE,NBREFE,NBDESC : DIMENSIONS DE VECTEURS POUR UN CHAMNO
!---------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
    integer(kind=8) :: ivalp, irefp, j, k
    integer(kind=8) :: nbrefe, nbvale, iad, nbvec
    character(len=19) :: vecsol, vestoc
    character(len=24) :: chaine
!
! -------------------------------------------------------------------
!----------------CREATION DU VECTEUR PRESSION -----------------------
!
! -----------CREATION DU TABLEAU DE VECTEURS CONTENANT---------------
!--------------------------LA PRESSION-------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: kb
    real(kind=8), pointer :: vale(:) => null()
    character(len=24), pointer :: refe(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    call jemarq()
    chaine = 'CBIDON'
!
    call codent(j, 'D0', chaine(1:5))
    zk24(iad+k-1) = vestoc(1:14)//chaine(1:5)
!
    call wkvect(zk24(iad+k-1) (1:19)//'.VALE', 'V V R', nbvale, ivalp)
    call wkvect(zk24(iad+k-1) (1:19)//'.REFE', 'V V K24', nbrefe, irefp)
!
    call jeveuo(vecsol//'.VALE', 'L', vr=vale)
    call jelira(vecsol//'.VALE', 'LONMAX', nbvec)
    call jeveuo(vecsol//'.REFE', 'L', vk24=refe)
!
!-------------STOCKAGE DANS LE VECTEUR CREE -------------------------
!
    b_n = to_blas_int(nbvec)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vale, b_incx, zr(ivalp), b_incy)
!
    do kb = 1, nbrefe
        zk24(irefp+kb-1) = refe(kb)
    end do
!
!
    call detrsd('CHAM_NO', vecsol)
!
    call jedema()
end subroutine
