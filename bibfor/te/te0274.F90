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
subroutine te0274(option, nomte)
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN THERMIQUE
!          CORRESPONDANT AU FLUX NON-LINEAIRE
!          ELEMENTS DE FACE 2D
!
!         OPTION : 'CHAR_THER_FLUNL'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
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
    real(kind=8) :: poids, r, nx, ny, theta, alpha, rbid, tpg, coorse(18)
    real(kind=8) :: vectt(9)
    integer :: nno, nnos, jgano, ndim, kp, npg, ipoids, ivf, idfde, igeom, i, j
    integer :: l, li, iflux, ivectt, nnop2, c(6, 9), ise, nse, itempr, itemps
    integer :: ibid
    character(len=8) :: coef, elrefe, alias8
    aster_logical :: laxi
!
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
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
    if (lteatt('AXIS', 'OUI')) then
        laxi = .true.
    else
        laxi = .false.
    end if
!
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
!
! RECUPERATION DE T-
    call jevech('PTEMPER', 'L', itempr)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PTEMPSR', 'L', itemps)
! FLUX NON-LIN
    call jevech('PFLUXNL', 'L', iflux)
    call jevech('PVECTTR', 'E', ivectt)
!
!====
! 1.3 PREALABLES LIES AUX CALCULS
!====
!
    theta = zr(itemps+2)
    coef = zk8(iflux)
    if (coef(1:7) .eq. '&FOZERO') goto 130
!
    call connec(nomte, nse, nnop2, c)
    do i = 1, nnop2
        vectt(i) = 0.d0
    end do
!
!====
! 2. CALCULS TERMES
!====
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
! BOUCLE SUR LES POINTS DE GAUSS
        do kp = 1, npg
            call vff2dn(ndim, nno, kp, ipoids, idfde, &
                        coorse, nx, ny, poids)
            tpg = 0.d0
            do i = 1, nno
! CALCUL DE T-
                l = (kp-1)*nno+i
                tpg = tpg+zr(itempr-1+c(ise, i))*zr(ivf+l-1)
            end do
!
! CALCUL DU JACOBIEN EN AXI
            if (laxi) then
                r = 0.d0
                do i = 1, nno
                    l = (kp-1)*nno+i
                    r = r+coorse(2*(i-1)+1)*zr(ivf+l-1)
                end do
                poids = poids*r
            end if
!
            call foderi(coef, tpg, alpha, rbid)
!
            if (theta < -0.5d0) then
                do i = 1, nno
                    li = ivf+(kp-1)*nno+i-1
                    vectt(c(ise, i)) = vectt(c(ise, i))+poids*alpha*zr(li)
                end do
            else
                do i = 1, nno
                    li = ivf+(kp-1)*nno+i-1
                    vectt(c(ise, i)) = vectt(c(ise, i))+poids*(1.d0-theta)*alpha*zr(li)
                end do
            end if
! FIN BOUCLE SUR LES PTS DE GAUSS
        end do
! FIN BOUCLE SUR LES SOUS-ELEMENTS
    end do
!
    do i = 1, nnop2
        zr(ivectt-1+i) = vectt(i)
    end do
!
130 continue
! FIN ------------------------------------------------------------------
end subroutine
