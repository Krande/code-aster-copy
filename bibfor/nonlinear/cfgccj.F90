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
subroutine cfgccj(resoco, nbliai, conjug)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "blas/daxpy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    character(len=24) :: resoco
    integer(kind=8) :: nbliai
    aster_logical :: conjug
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - GCP)
!
! CONJUGAISON
!
! ----------------------------------------------------------------------
!
! ON FAIT LA CONJUGAISON DE POLAK-RIBIERE :
!     - SI L'ETAT DE CONTACT EST LE MEME D'UNE ITERATION SUR L'AUTRE
!     - TOUTES LES FREQ ITERATIONS
!     - SI DIRECT EST UNE DIRECTION DE DESCENTE I.E.
!                                                  (DIRECT' SGRAD+)>0
!  NB1 :  FORMULE DE CONJUGAISON EN PRESENCE DE PRECONDITIONNEUR :
!         GAMMA = (SGRADP' (SGRPRP - SGRPRM)) / (SGRADM' SGRPRM))
!  NB2 : LA CONJUGAISON DE FLETCHER-REEVES EST : GAMMA = NUMER/DENOM
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT
! I/O CONJUG : DIRECTIONS CONJUGUEES
!
!
!
!
!
    integer(kind=8) :: ifm, niv
    real(kind=8) :: numer, numer2, denom, gamma
    character(len=19) :: sgradm, sgradp, sgrprm, sgrprp, direct
    integer(kind=8) :: jsgram, jsgrap, jsgprm, jsgprp, jdirec
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
!
! --- INITIALISATIONS
!
    gamma = 0.d0
!
! --- ACCES VECTEURS DE TRAVAIL
!
    sgradm = resoco(1:14)//'.SGDM'
    sgradp = resoco(1:14)//'.SGDP'
    sgrprm = resoco(1:14)//'.SGPM'
    sgrprp = resoco(1:14)//'.SGPP'
    direct = resoco(1:14)//'.DIRE'
    call jeveuo(sgradm, 'L', jsgram)
    call jeveuo(sgradp, 'L', jsgrap)
    call jeveuo(sgrprm, 'L', jsgprm)
    call jeveuo(sgrprp, 'L', jsgprp)
    call jeveuo(direct, 'E', jdirec)
!
! --- COEFFICIENT DE DIRECTION
!
    if (conjug) then
        b_n = to_blas_int(nbliai)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        numer = ddot(b_n, zr(jsgrap), b_incx, zr(jsgprp), b_incy)
        b_n = to_blas_int(nbliai)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        numer2 = ddot(b_n, zr(jsgrap), b_incx, zr(jsgprm), b_incy)
        b_n = to_blas_int(nbliai)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        denom = ddot(b_n, zr(jsgram), b_incx, zr(jsgprm), b_incy)
        gamma = (numer-numer2)/denom
    end if
!
! --- MISE A JOUR DIRECTION
!
    b_n = to_blas_int(nbliai)
    b_incx = to_blas_int(1)
    call dscal(b_n, gamma, zr(jdirec), b_incx)
    b_n = to_blas_int(nbliai)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, zr(jsgprp), b_incx, zr(jdirec), &
               b_incy)
!
! --- AFFICHAGE
!
    if (conjug) then
        if (niv .eq. 2) then
            write (ifm, *) '<CONTACT><CALC> CONJUGAISON DES DIRECTIONS '//&
     &      'DE RECHERCHE, GAMMA=', gamma
        end if
    end if
!
    call jedema()
!
end subroutine
