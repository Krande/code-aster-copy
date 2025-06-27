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
subroutine te0315(option, nomte)
    implicit none
!
!.......................................................................
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES DE FLUX FLUIDE EN MECANIQUE
!          ELEMENTS ISOPARAMETRIQUES 1D
!
!          OPTION : 'CHAR_THER_ACCE_R 'OU 'CHAR_THER_ACCE_X'
!                    OU 'CHAR_THER_ACCE_Y'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/rcvalb.h"
#include "asterfort/vff2dn.h"
!
    integer(kind=8) :: icodre(1)
    character(len=8) :: fami, poum
    character(len=16) :: nomte, option
    real(kind=8) :: poids, nx, ny, norm(2), acloc(2, 3)
    real(kind=8) :: acc(2, 4), flufn(4)
    integer(kind=8) :: ipoids, ivf, idfde, igeom
    integer(kind=8) :: nno, kp, npg, ivectt, imate
    integer(kind=8) :: ldec, spt
    aster_logical :: laxi
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacce, idim, itemp, jgano, k, ndim
    integer(kind=8) :: nnos
    real(kind=8) :: r, rho(1)
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PVECTTR', 'E', ivectt)
    fami = 'RIGI'
    spt = 1
    poum = '+'
!
    if (option(16:16) .eq. 'R') then
        call jevech('PACCELR', 'L', iacce)
    else
        if ((option(16:16) .eq. 'X') .or. (option(16:16) .eq. 'Y')) then
            call jevech('PTEMPER', 'L', itemp)
        end if
    end if
!
    k = 0
    do i = 1, nno
        if (option(16:16) .eq. 'R') then
            do idim = 1, 2
                k = k+1
                acloc(idim, i) = zr(iacce+k-1)
            end do
        else if ((option(16:16) .eq. 'X')) then
            k = k+1
            acloc(1, i) = zr(itemp+k-1)
            acloc(2, i) = 0.d0
        else if (option(16:16) .eq. 'Y') then
            k = k+1
            acloc(1, i) = 0.d0
            acloc(2, i) = zr(itemp+k-1)
        end if
    end do
!
    do i = 1, nno
        zr(ivectt+i-1) = 0.d0
    end do
!
!     BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg
        ldec = (kp-1)*nno
!
        call rcvalb(fami, kp, spt, poum, zi(imate), &
                    ' ', 'THER', 0, ' ', [0.d0], &
                    1, 'RHO_CP', rho, icodre, 1)
!
        nx = 0.d0
        ny = 0.d0
!        --- ON CALCULE L ACCEL AU POINT DE GAUSS
        acc(1, kp) = 0.d0
        acc(2, kp) = 0.d0
        do i = 1, nno
            acc(1, kp) = acc(1, kp)+acloc(1, i)*zr(ivf+ldec+i-1)
            acc(2, kp) = acc(2, kp)+acloc(2, i)*zr(ivf+ldec+i-1)
        end do
!
        call vff2dn(ndim, nno, kp, ipoids, idfde, &
                    zr(igeom), nx, ny, poids)
        norm(1) = nx
        norm(2) = ny
        flufn(kp) = 0.d0
!
! CALCUL DU FLUX FLUIDE NORMAL AU POINT DE GAUSS
!
        flufn(kp) = acc(1, kp)*norm(1)+acc(2, kp)*norm(2)
!
! CAS AXISYMETRIQUE
!
        if (laxi) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(igeom+2*(i-1))*zr(ivf+ldec+i-1)
            end do
            poids = poids*r
        end if
!
        do i = 1, nno
            zr(ivectt+i-1) = zr(ivectt+i-1)+poids*flufn(kp)*rho(1)*zr(ivf+ldec+i-1)
        end do
    end do
!
end subroutine
