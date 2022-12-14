! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine te0251(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/connec.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES TANGENTES ELEMENTAIRES
!                          OPTION : 'MTAN_THER_FLUXNL'
!                          ELEMENTS DE FACE 2D
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! ......................................................................
!
    real(kind=8) :: poids, r, nx, ny, theta, alphap, rbid, tpg
    real(kind=8) :: mrigt(9, 9), coorse(18)
    integer :: nno, nnos, ndim, kp, npg, ipoids, ivf, idfde, jgano
    integer :: imattt, i, j, ij, l, li, lj, ibid
    integer :: c(6, 9), ise, nse, nnop2
    integer :: igeom, iflux, itempi, itemps
    aster_logical :: laxi
    character(len=8) :: elrefe, alias8
!
!
    call elref1(elrefe)
!
    if (lteatt('LUMPE','OUI')) then
        call teattr('S', 'ALIAS8', alias8, ibid)
        if (alias8(6:8) .eq. 'SE3') elrefe='SE2'
    endif
!
    call elrefe_info(elrefe=elrefe, fami='RIGI', ndim=ndim, nno=nno, nnos=nnos,&
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    laxi = .false.
    if (lteatt('AXIS','OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PFLUXNL', 'L', iflux)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PMATTTR', 'E', imattt)
!
    if (zk8(iflux) (1:7) .eq. '&FOZERO') goto 120
    theta = zr(itemps+2)
!
    call connec(nomte, nse, nnop2, c)
!
    do i = 1, nnop2
        do j = 1, nnop2
            mrigt(i,j) = 0.d0
        end do
    end do
!
! --- CALCUL ISO-P2 : BOUCLE SUR LES SOUS-ELEMENTS -------
!
    do ise = 1, nse
!
        do i = 1, nno
            do j = 1, 2
                coorse(2* (i-1)+j) = zr(igeom-1+2* (c(ise,i)-1)+j)
            end do
        end do
!
        do kp = 1, npg
            call vff2dn(ndim, nno, kp, ipoids, idfde,&
                        coorse, nx, ny, poids)
            r = 0.d0
            tpg = 0.d0
            do i = 1, nno
                l = (kp-1)*nno + i
                r = r + coorse(2* (i-1)+1)*zr(ivf+l-1)
                tpg = tpg + zr(itempi-1+c(ise,i))*zr(ivf+l-1)
            end do
            if (laxi) poids = poids*r
            call foderi(zk8(iflux), tpg, rbid, alphap)
            ij = imattt - 1
            do i = 1, nno
                li = ivf + (kp-1)*nno + i - 1
                do j = 1, i
                    lj = ivf + (kp-1)*nno + j - 1
                    ij = ij + 1
                    mrigt(c(ise,i),c(ise,j)) = mrigt(&
                                               c(ise, i),&
                                               c(ise, j) ) - poids*theta*alphap*zr(li)* zr(lj&
                                               )
                end do
            end do
        end do
    end do
!
! MISE SOUS FORME DE VECTEUR
!
    ij = imattt - 1
    do i = 1, nnop2
        do j = 1, i
            ij = ij + 1
            zr(ij) = mrigt(i,j)
        end do
    end do
120 continue
end subroutine
