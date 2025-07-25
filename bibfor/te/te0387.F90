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
subroutine te0387(option, nomte)
! ......................................................................
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'MTAN_THER_PARO_F'
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icode, igeom2, ihechp, jgano, nbres, ndim, nnos
!
!-----------------------------------------------------------------------
    parameter(nbres=3)
    character(len=8) :: nompar(nbres)
    real(kind=8) :: valpar(nbres), poids, poids1, poids2, r, r1, r2
    real(kind=8) :: z, z1, z2, hechp, nx, ny, mat(6)
    integer(kind=8) :: nno, kp, npg, ipoids, ivf, idfde, igeom
    integer(kind=8) :: itemps, imatt, k, i, j, l, li, lj
    aster_logical :: laxi
!-----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PHECHPF', 'L', ihechp)
    call jevech('PMATTTR', 'E', imatt)
    if (nomte(5:8) .eq. 'SE22') then
        igeom2 = igeom+4
    else if (nomte(5:8) .eq. 'SE33') then
        igeom2 = igeom+6
    end if
!
    do kp = 1, npg
        call vff2dn(ndim, nno, kp, ipoids, idfde, &
                    zr(igeom), nx, ny, poids1)
        call vff2dn(ndim, nno, kp, ipoids, idfde, &
                    zr(igeom2), nx, ny, poids2)
        r = 0.d0
        r1 = 0.d0
        r2 = 0.d0
        z = 0.d0
        z1 = 0.d0
        z2 = 0.d0
        do i = 1, nno
            l = (kp-1)*nno+i
            r1 = r1+zr(igeom+2*i-2)*zr(ivf+l-1)
            r2 = r2+zr(igeom2+2*i-2)*zr(ivf+l-1)
            z1 = z1+zr(igeom+2*i-1)*zr(ivf+l-1)
            z2 = z2+zr(igeom2+2*i-1)*zr(ivf+l-1)
        end do
        poids = (poids1+poids2)/2.0d0
        r = (r1+r2)/2.d0
        z = (z1+z2)/2.d0
        if (laxi) poids = poids*r
        valpar(1) = r
        nompar(1) = 'X'
        valpar(2) = z
        nompar(2) = 'Y'
        valpar(3) = zr(itemps)
        nompar(3) = 'INST'
        call fointe('A', zk8(ihechp), 3, nompar, valpar, &
                    hechp, icode)
        ASSERT(icode .eq. 0)
        k = 0
        do i = 1, nno
            li = ivf+(kp-1)*nno+i-1
            do j = 1, i
                lj = ivf+(kp-1)*nno+j-1
                k = k+1
                mat(k) = poids*zr(li)*zr(lj)*hechp
            end do
        end do
        if (nomte(5:8) .eq. 'SE22') then
            zr(imatt-1+1) = zr(imatt-1+1)+mat(1)
            zr(imatt-1+2) = zr(imatt-1+2)+mat(2)
            zr(imatt-1+3) = zr(imatt-1+3)+mat(3)
            zr(imatt-1+4) = zr(imatt-1+4)-mat(1)
            zr(imatt-1+5) = zr(imatt-1+5)-mat(2)
            zr(imatt-1+6) = zr(imatt-1+6)+mat(1)
            zr(imatt-1+7) = zr(imatt-1+7)-mat(2)
            zr(imatt-1+8) = zr(imatt-1+8)-mat(3)
            zr(imatt-1+9) = zr(imatt-1+9)+mat(2)
            zr(imatt-1+10) = zr(imatt-1+10)+mat(3)
        else if (nomte(5:8) .eq. 'SE33') then
            zr(imatt-1+1) = zr(imatt-1+1)+mat(1)
            zr(imatt-1+2) = zr(imatt-1+2)+mat(2)
            zr(imatt-1+3) = zr(imatt-1+3)+mat(3)
            zr(imatt-1+4) = zr(imatt-1+4)+mat(4)
            zr(imatt-1+5) = zr(imatt-1+5)+mat(5)
            zr(imatt-1+6) = zr(imatt-1+6)+mat(6)
            zr(imatt-1+7) = zr(imatt-1+7)-mat(1)
            zr(imatt-1+8) = zr(imatt-1+8)-mat(2)
            zr(imatt-1+9) = zr(imatt-1+9)-mat(4)
            zr(imatt-1+10) = zr(imatt-1+10)+mat(1)
            zr(imatt-1+11) = zr(imatt-1+11)-mat(2)
            zr(imatt-1+12) = zr(imatt-1+12)-mat(3)
            zr(imatt-1+13) = zr(imatt-1+13)-mat(5)
            zr(imatt-1+14) = zr(imatt-1+14)+mat(2)
            zr(imatt-1+15) = zr(imatt-1+15)+mat(3)
            zr(imatt-1+16) = zr(imatt-1+16)-mat(4)
            zr(imatt-1+17) = zr(imatt-1+17)-mat(5)
            zr(imatt-1+18) = zr(imatt-1+18)-mat(6)
            zr(imatt-1+19) = zr(imatt-1+19)+mat(4)
            zr(imatt-1+20) = zr(imatt-1+20)+mat(5)
            zr(imatt-1+21) = zr(imatt-1+21)+mat(6)
        end if
    end do
end subroutine
