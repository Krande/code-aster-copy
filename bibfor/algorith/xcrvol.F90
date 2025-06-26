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
subroutine xcrvol(nse, ndim, jcnse, nnose, jpint, &
                  igeom, elrefp, inoloc, nbnoma, jcesd3, &
                  jcesl3, jcesv3, numa2, iheav, nfiss, &
                  vhea, jcesd8, jcesl8, jcesv8, lfiss, &
                  vtot)
!
! person_in_charge: samuel.geniaut at edf.fr
!
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/iselli.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_code.h"
    integer(kind=8) :: nse, ndim, jcnse, nnose, jpint, igeom, inoloc, nfiss, iheav
    character(len=8) :: elrefp
    integer(kind=8) :: nbnoma, jcesd3, jcesl3, jcesv3, numa2, jcesd8, jcesl8, jcesv8
    real(kind=8) :: vhea, vtot
    aster_logical :: lfiss
!
!  BUT: ESTIMATION CRITERE DE RIGIDITE
!
    real(kind=8) :: co(ndim+1, ndim), mat(ndim, ndim), vse, bary(ndim)
    real(kind=8) :: point(ndim), he(nfiss)
    real(kind=8) :: ff(nbnoma), dfdi(nbnoma, ndim), xe(ndim), deriv
    integer(kind=8) :: ise, ino2, i, j, iad, k, hea_se, hea_no
!
! ----------------------------------------------------------------------
!
!     BOUCLE SUR LES SOUS ELEMENTS
    do ise = 1, nse
!       RECUPERATION DES COORDONNEES DES NOEUDS DU SOUS ELEMENT
        do i = 1, ndim+1
            ino2 = zi(jcnse-1+nnose*(ise-1)+i)
! on ne recupere pas les noeuds milieux
            ASSERT(ino2 .le. 2000)
            if (ino2 .gt. 1000) then
                do j = 1, ndim
                    co(i, j) = zr(jpint-1+ndim*(ino2-1000-1)+j)
                end do
            else
                do j = 1, ndim
                    co(i, j) = zr(igeom-1+ndim*(ino2-1)+j)
                end do
            end if
        end do
        do i = 1, ndim
            do j = 1, ndim
                mat(i, j) = co(1, j)-co(i+1, j)
            end do
        end do
!
!       CALCUL DU VOLUME DU SOUS ELEMENTS (DÉTERMINANT)
!
        vse = 0.d0
        if (ndim .eq. 2) then
            vse = abs(mat(1, 1)*mat(2, 2)-mat(2, 1)*mat(1, 2))/2
        else if (ndim .eq. 3) then
            vse = abs( &
                 mat(1, 1)*mat(2, 2)*mat(3, 3)+mat(2, 1)*mat(3, 2)*mat(1, 3)+mat(3, 1)*mat(1, 2)*ma&
                  &t(2, 3)-mat(3, 1)*mat(2, 2)*mat(1, 3)-mat(2, 1)*mat(1, 2)*mat(3, 3)-mat(1, 1)*m&
                  &at(3, 2)*mat(2, 3) &
                  )/6
        end if
!
!       CALCUL DU BARYCENTRE
!
        bary(:) = 0.d0
        do j = 1, ndim
            do i = 1, ndim+1
                bary(j) = bary(j)+co(i, j)
            end do
            bary(j) = bary(j)/(ndim+1)
        end do
!
!        CALCUL DES DERIVEES DES FONCTIONS DE FORME
!
        call reeref(elrefp, nbnoma, zr(igeom), bary, ndim, &
                    xe, ff, dfdi=dfdi)
        deriv = 0.d0
        do i = 1, ndim
            deriv = max(abs(dfdi(inoloc, i)), deriv)
        end do
!       EN QUADRATIQUE : AUGMENTATION DU NOMBRE DE POINTS
        if (.not. iselli(elrefp)) then
            do k = 1, ndim+1
                point(:) = 0.d0
                do j = 1, ndim
                    point(j) = (bary(j)+co(k, j))/2
                end do
                call reeref(elrefp, nbnoma, zr(igeom), point, ndim, &
                            xe, ff, dfdi=dfdi)
                do i = 1, ndim
                    deriv = max(abs(dfdi(inoloc, i)), deriv)
                end do
            end do
        end if
        vse = vse*deriv**2
!
!  EN QUADRATIQUE :: MULTIPLICATION PAR UN TERME CORRECTIF CAR L INTEGRATION EST IMPRECISE
!    ASYMPTOTIQUEMENT DFDI EST PROCHE DE EPS=VSE**1/NDIM
!    L INTEGRALE DE DFDI**2 VARIE EN EPS**3
        if (.not. iselli(elrefp) .and. lfiss) vse = vse*vse**(3/ndim)
!       DETERMINATION DU SIGNE DU SOUS ELEMENT
        do i = 1, nfiss
            call cesexi('S', jcesd3, jcesl3, numa2, 1, &
                        i, ise, iad)
            he(i) = zi(jcesv3-1+iad)
        end do
!       CALCUL DU CODE DU SOUS ELEMENT
        hea_se = xcalc_code(nfiss, he_real=[he])
!       CALCUL DU CODE DU DDL HEAVISIDE
        call cesexi('C', jcesd8, jcesl8, numa2, inoloc, &
                    1, iheav, iad)
        hea_no = zi(jcesv8-1+iad)
        if (hea_se .eq. hea_no) then
            vhea = vhea+vse
        end if
        vtot = vtot+vse
    end do
!
end subroutine
