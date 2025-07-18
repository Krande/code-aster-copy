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
subroutine te0166(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL FORCE DE PESANTEUR POUR MEPOULI
!                          OPTION : 'CHAR_MECA_PESA_R'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: icodre(1)
    real(kind=8) :: rho(1), a, w(9), l1(3), l2(3), l10(3), l20(3)
    real(kind=8) :: norml1, norml2, norl10, norl20, l0, norm1p, norm2p
    real(kind=8) :: poids(3)
    character(len=8) :: fami, poum
    integer(kind=8) :: i, neu, neum1, kc, ic, ivectu, ipesa, kpg, spt
    integer(kind=8) :: igeom, imate, lsect, idepla, ideplp, iret
!
!
!-----------------------------------------------------------------------
    real(kind=8) :: r8b
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    r8b = 0.d0
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'ELAS', 0, ' ', [r8b], &
                1, 'RHO', rho, icodre, 1)
    call jevech('PCACABL', 'L', lsect)
    a = zr(lsect)
!
    call tecach('ONO', 'PDEPLMR', 'L', iret, iad=idepla)
    if (iret .ne. 0) then
        call utmess('F', 'CALCULEL6_78')
    end if
    call jevech('PDEPLPR', 'L', ideplp)
!
    do i = 1, 9
        w(i) = zr(idepla-1+i)+zr(ideplp-1+i)
    end do
!
    do kc = 1, 3
        l1(kc) = w(kc)+zr(igeom-1+kc)-w(6+kc)-zr(igeom+5+kc)
        l10(kc) = zr(igeom-1+kc)-zr(igeom+5+kc)
    end do
    do kc = 1, 3
        l2(kc) = w(3+kc)+zr(igeom+2+kc)-w(6+kc)-zr(igeom+5+kc)
        l20(kc) = zr(igeom+2+kc)-zr(igeom+5+kc)
    end do
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norml1 = ddot(b_n, l1, b_incx, l1, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norml2 = ddot(b_n, l2, b_incx, l2, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norl10 = ddot(b_n, l10, b_incx, l10, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norl20 = ddot(b_n, l20, b_incx, l20, b_incy)
    norml1 = sqrt(norml1)
    norml2 = sqrt(norml2)
    norl10 = sqrt(norl10)
    norl20 = sqrt(norl20)
    l0 = norl10+norl20
!
    call jevech('PPESANR', 'L', ipesa)
    call jevech('PVECTUR', 'E', ivectu)
!
    norm1p = norml1*l0/(norml1+norml2)
    norm2p = norml2*l0/(norml1+norml2)
    poids(1) = rho(1)*a*norm1p*zr(ipesa)/2.d0
    poids(2) = rho(1)*a*norm2p*zr(ipesa)/2.d0
    poids(3) = poids(1)+poids(2)
!
!
    do neu = 1, 3
        neum1 = neu-1
        do ic = 1, 3
            zr(ivectu+3*neum1+ic-1) = poids(neu)*zr(ipesa+ic)
        end do
    end do
!
end subroutine
