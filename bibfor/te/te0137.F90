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
subroutine te0137(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/connec.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/vff2dn.h"
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'RESI_THER_COEF_F'
!                          OPTION : 'RESI_THER_RAYO_F'
!                          ELEMENTS 2D
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: nbres
    parameter(nbres=3)
    character(len=8) :: nompar(nbres), elrefe, alias8
    real(kind=8) :: valpar(nbres), poids, r, z, nx, ny, tpg
    real(kind=8) :: coenp1, sigma, epsil, tz0
    real(kind=8) :: coorse(18), vectt(9)
    integer(kind=8) :: nno, nnos, ndim, kp, npg, ipoids, ivf, idfde, jgano, igeom
    integer(kind=8) :: itemps, iveres, i, j, l, li, iech, iray, itemp, icode, ier
    integer(kind=8) :: nnop2, c(6, 9), ise, nse, ibid
    aster_logical :: laxi
!
!
    call elref1(elrefe)
!
    if (lteatt('LUMPE', 'OUI')) then
        call teattr('S', 'ALIAS8', alias8, ibid)
        if (alias8(6:8) .eq. 'SE3') elrefe = 'SE2'
    end if
!
    call elrefe_info(elrefe=elrefe, fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    tz0 = r8t0()
!
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    if (option(11:14) .eq. 'COEF') then
        call jevech('PCOEFHF', 'L', iech)
    else if (option(11:14) .eq. 'RAYO') then
        call jevech('PRAYONF', 'L', iray)
    end if
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPEI', 'L', itemp)
    call jevech('PRESIDU', 'E', iveres)
!
    call connec(nomte, nse, nnop2, c)
!
    do i = 1, nnop2
        vectt(i) = 0.d0
    end do
!
! BOUCLE SUR LES SOUS-ELEMENTS
!
    do ise = 1, nse
!
        do i = 1, nno
            do j = 1, 2
                coorse(2*(i-1)+j) = zr(igeom-1+2*(c(ise, i)-1)+j)
            end do
        end do
!
        do kp = 1, npg
            call vff2dn(ndim, nno, kp, ipoids, idfde, &
                        coorse, nx, ny, poids)
            r = 0.d0
            z = 0.d0
            tpg = 0.d0
            do i = 1, nno
                l = (kp-1)*nno+i
                r = r+coorse(2*(i-1)+1)*zr(ivf+l-1)
                z = z+coorse(2*(i-1)+2)*zr(ivf+l-1)
                tpg = tpg+zr(itemp-1+c(ise, i))*zr(ivf+l-1)
            end do
            if (laxi) poids = poids*r
            valpar(1) = r
            nompar(1) = 'X'
            valpar(2) = z
            nompar(2) = 'Y'
            nompar(3) = 'INST'
            valpar(3) = zr(itemps)
            if (option(11:14) .eq. 'COEF') then
                call fointe('A', zk8(iech), 3, nompar, valpar, &
                            coenp1, icode)
                ASSERT(icode .eq. 0)
                do i = 1, nno
                    li = ivf+(kp-1)*nno+i-1
                    vectt(c(ise, i)) = vectt(c(ise, i))+poids*zr(li)*coenp1*tpg
                end do
            else if (option(11:14) .eq. 'RAYO') then
                call fointe('A', zk8(iray), 4, nompar, valpar, &
                            sigma, ier)
                ASSERT(ier .eq. 0)
                call fointe('A', zk8(iray+1), 4, nompar, valpar, &
                            epsil, ier)
                ASSERT(ier .eq. 0)
                do i = 1, nno
                    li = ivf+(kp-1)*nno+i-1
                    vectt(c(ise, i)) = vectt(c(ise, i))+poids*zr(li)*sigma*epsil*(tpg+tz0 &
                                                                                  )**4
                end do
            end if
!
        end do
    end do
!
    do i = 1, nnop2
        zr(iveres-1+i) = vectt(i)
    end do
!
end subroutine
