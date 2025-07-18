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
subroutine cffpm1(resoco, nbliai, ndim, nesmax)
!
!
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
!
    character(len=24) :: resoco
    integer(kind=8) :: nbliai, ndim, nesmax
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - PENALISATION)
!
! CALCUL DE LA MATRICE FRO1 = E_T*AaT
!
! ----------------------------------------------------------------------
!
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT POSSIBLES
! IN  NDIM   : DIMENSION DU PROBLEME
! IN  NESMAX : NOMBRE MAX DE NOEUDS ESCLAVES
!
!
!
!
    integer(kind=8) :: ndlmax
    parameter(ndlmax=30)
    integer(kind=8) :: jdecal, nbddl
    real(kind=8) :: xmu, jeuini
    integer(kind=8) :: iliai
    character(len=19) :: mu
    integer(kind=8) :: jmu
    character(len=24) :: appoin
    integer(kind=8) :: japptr
    character(len=24) :: apcofr
    integer(kind=8) :: japcof
    character(len=24) :: jeux
    integer(kind=8) :: jjeux
    character(len=19) :: fro1
    integer(kind=8) :: jfro11, jfro12
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    appoin = resoco(1:14)//'.APPOIN'
    jeux = resoco(1:14)//'.JEUX'
    mu = resoco(1:14)//'.MU'
    apcofr = resoco(1:14)//'.APCOFR'
    fro1 = resoco(1:14)//'.FRO1'
    call jeveuo(jeux, 'L', jjeux)
    call jeveuo(mu, 'L', jmu)
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apcofr, 'L', japcof)
!
! --- CALCUL DE LA MATRICE E_T*AaT
!
    do iliai = 1, nbliai
!
! ----- INITIALISATION DES COLONNES
!
        call jeveuo(jexnum(fro1, iliai), 'E', jfro11)
        call r8inir(ndlmax, 0.d0, zr(jfro11), 1)
        if (ndim .eq. 3) then
            call jeveuo(jexnum(fro1, iliai+nbliai), 'E', jfro12)
            call r8inir(ndlmax, 0.d0, zr(jfro12), 1)
        end if
!
! ----- LA LIAISON EST-ELLE ACTIVE ?
!
        jeuini = zr(jjeux+3*(iliai-1)+1-1)
!
! ----- CALCUL
!
        if (jeuini .lt. r8prem()) then
            jdecal = zi(japptr+iliai-1)
            nbddl = zi(japptr+iliai)-zi(japptr+iliai-1)
            xmu = zr(jmu+3*nbliai+iliai-1)
            b_n = to_blas_int(nbddl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, xmu, zr(japcof+jdecal), b_incx, zr(jfro11), &
                       b_incy)
            if (ndim .eq. 3) then
                b_n = to_blas_int(nbddl)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, xmu, zr(japcof+jdecal+ndlmax*nesmax), b_incx, zr(jfro12), &
                           b_incy)
            end if
        end if
!
        call jelibe(jexnum(fro1, iliai))
        if (ndim .eq. 3) then
            call jelibe(jexnum(fro1, iliai+nbliai))
        end if
    end do
!
    call jedema()
!
end subroutine
