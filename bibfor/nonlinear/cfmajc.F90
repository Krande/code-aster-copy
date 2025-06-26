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
subroutine cfmajc(resoco, neq, nbliac)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "blas/daxpy.h"
!
    character(len=24) :: resoco
    integer(kind=8) :: neq
    integer(kind=8) :: nbliac
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - RESOLUTION)
!
! MISE A JOUR DU VECTEUR SOLUTION ITERATION DE CONTACT
!
! ----------------------------------------------------------------------
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NEQ    : NOMBRE D'EQUATIONS
! IN  NBLIAC : NOMBRE DE LIAISONS ACTIVES
!
!
!
!
    integer(kind=8) :: iliai, iliac
    character(len=19) :: mu, cm1a, liac
    integer(kind=8) :: jmu, jcm1a, jliac
    character(len=19) :: ddelt
    real(kind=8), pointer :: vale(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- ACCES STRUCTURES DE DONNEES DE CONTACT
!
    cm1a = resoco(1:14)//'.CM1A'
    mu = resoco(1:14)//'.MU'
    liac = resoco(1:14)//'.LIAC'
    call jeveuo(mu, 'L', jmu)
    call jeveuo(liac, 'L', jliac)
!
! --- ACCES AUX CHAMPS DE TRAVAIL
!
    ddelt = resoco(1:14)//'.DDEL'
    call jeveuo(ddelt(1:19)//'.VALE', 'E', vr=vale)
!
! --- AJOUT DES LIAISONS DE CONTACT ACTIVES
!
    do iliac = 1, nbliac
        iliai = zi(jliac-1+iliac)
        call jeveuo(jexnum(cm1a, iliai), 'L', jcm1a)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -zr(jmu-1+iliac), zr(jcm1a), b_incx, vale, &
                   b_incy)
        call jelibe(jexnum(cm1a, iliai))
    end do
!
    call jedema()
!
end subroutine
