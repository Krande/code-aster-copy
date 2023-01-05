! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine te0070(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/connec.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'RIGI_THER_COEF_R'
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    character(len=8) :: elrefe, alias8
    real(kind=8) :: poids, r, nx, ny, theta
    real(kind=8) :: mrigt(9, 9), coorse(18)
    integer :: nno, kp, npg, ipoids, ivf, idfde, igeom, jgano
    integer :: imattt, i, j, ij, l, li, lj, icoefh, ndim, nnos, itemps
    integer :: c(6, 9), ise, nse, nnop2, ibid
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
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCOEFHR', 'L', icoefh)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PMATTTR', 'E', imattt)
!
    theta = zr(itemps+2)
!
    call connec(nomte, nse, nnop2, c)
!
    do i = 1, nnop2
        do j = 1, nnop2
            mrigt(i, j) = 0.d0
        end do
    end do
!
! --- CALCUL ISO-P2 : BOUCLE SUR LES SOUS-ELEMENTS -------
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
            if (laxi) then
                r = 0.d0
                do i = 1, nno
                    l = (kp-1)*nno+i
                    r = r+coorse(2*(i-1)+1)*zr(ivf+l-1)
                end do
                poids = poids*r
            end if
            ij = imattt-1
            do i = 1, nno
                li = ivf+(kp-1)*nno+i-1
                do j = 1, i
                    lj = ivf+(kp-1)*nno+j-1
                    ij = ij+1
                    mrigt(c(ise, i), c(ise, j)) = mrigt( &
                                               c(ise, i), &
                                               c(ise, j))+poids*theta*zr(li)*zr(lj)*zr(icoef&
                                               &h &
                                               )
                end do
            end do
        end do
    end do
!
! MISE SOUS FORME DE VECTEUR
!
    ij = imattt-1
    do i = 1, nnop2
        do j = 1, i
            ij = ij+1
            zr(ij) = mrigt(i, j)
        end do
    end do
end subroutine
