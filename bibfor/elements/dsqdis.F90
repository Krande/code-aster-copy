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

subroutine dsqdis(xyzl, caraq4, df, dci, an)
    implicit none
#include "jeveux.h"
#include "asterfort/dsxhft.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jquad4.h"
#include "asterfort/mgauss.h"
    real(kind=8) :: xyzl(3, *), df(3, 3), dci(2, 2), an(4, 12), caraq4(*)
!     MATRICE AN(4,12) DU CISAILLEMENT POUR LE DSQ
!     -------------------------------------------------------
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: nc, k, ic, int, j, i, iret
    real(kind=8) :: qsi, eta, peta, meta, pqsi, mqsi, det, jacob(5)
    real(kind=8) :: l(4), c(4), s(4), x(4), y(4)
    real(kind=8) :: hft2(2, 6), tb(6, 12), ta(6, 4), dt(2, 6)
    real(kind=8) :: dib(2, 12), dia(2, 4), aw(4, 12), aa(4, 4), aai(4, 4)
!     ------------------------------------------------------------------
!
    call elrefe_info(elrefe='SE2', fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, &
                     jdfd2=idfd2, jgano=jgano)
    nc = 4
!
    c(1) = caraq4(13)
    c(2) = caraq4(14)
    c(3) = caraq4(15)
    c(4) = caraq4(16)
    s(1) = caraq4(17)
    s(2) = caraq4(18)
    s(3) = caraq4(19)
    s(4) = caraq4(20)
!
    l(1) = caraq4(9)
    l(2) = caraq4(10)
    l(3) = caraq4(11)
    l(4) = caraq4(12)
!
    x(1) = caraq4(1)
    x(2) = caraq4(2)
    x(3) = caraq4(3)
    x(4) = caraq4(4)
    y(1) = caraq4(5)
    y(2) = caraq4(6)
    y(3) = caraq4(7)
    y(4) = caraq4(8)
!
    tb(:, :) = 0.d0
    tb(3, 2) = 0.25d0
    tb(3, 5) = -0.25d0
    tb(3, 8) = 0.25d0
    tb(3, 11) = -0.25d0
    tb(6, 3) = 0.25d0
    tb(6, 6) = -0.25d0
    tb(6, 9) = 0.25d0
    tb(6, 12) = -0.25d0
!
    do ic = 1, nc
!
        dib(:, :) = 0.d0
        dia(:, :) = 0.d0
        do int = 1, 2
!
            if (ic .eq. 1) then
                qsi = -zr(icoopg-1+ndim*(int-1)+1)
                eta = -zr(ipoids-1+int)
            else if (ic .eq. 2) then
                qsi = zr(ipoids-1+int)
                eta = -zr(icoopg-1+ndim*(int-1)+1)
            else if (ic .eq. 3) then
                qsi = zr(icoopg-1+ndim*(int-1)+1)
                eta = zr(ipoids-1+int)
            else if (ic .eq. 4) then
                qsi = -zr(ipoids-1+int)
                eta = zr(icoopg-1+ndim*(int-1)+1)
            end if
!
            call jquad4(xyzl, qsi, eta, jacob)
!
            peta = 1.d0+eta
            meta = 1.d0-eta
            pqsi = 1.d0+qsi
            mqsi = 1.d0-qsi
!
            ta(:, :) = 0.d0
            ta(1, 1) = -meta*c(1)
            ta(1, 3) = -peta*c(3)
            ta(2, 2) = -pqsi*c(2)
            ta(2, 4) = -mqsi*c(4)
            ta(3, 1) = qsi*c(1)
            ta(3, 2) = -eta*c(2)
            ta(3, 3) = -qsi*c(3)
            ta(3, 4) = eta*c(4)
            ta(4, 1) = -meta*s(1)
            ta(4, 3) = -peta*s(3)
            ta(5, 2) = -pqsi*s(2)
            ta(5, 4) = -mqsi*s(4)
            ta(6, 1) = qsi*s(1)
            ta(6, 2) = -eta*s(2)
            ta(6, 3) = -qsi*s(3)
            ta(6, 4) = eta*s(4)
!
            call dsxhft(df, jacob(2), hft2)
!
!           -------- PRODUIT DCI.HFT2 ----------------------------------
            dt(:, :) = 0.d0
            do j = 1, 6
                dt(1, j) = dt(1, j)+dci(1, 1)*hft2(1, j)+dci(1, 2)*hft2(2, j)
                dt(2, j) = dt(2, j)+dci(2, 1)*hft2(1, j)+dci(2, 2)*hft2(2, j)
            end do
!           -------- PRODUIT DT.TB -------------------------------------
            do j = 1, 12
                do k = 1, 6
                    dib(1, j) = dib(1, j)+dt(1, k)*tb(k, j)
                    dib(2, j) = dib(2, j)+dt(2, k)*tb(k, j)
                end do
            end do
!           -------- PRODUIT DT.TA -------------------------------------
            do j = 1, 4
                do k = 1, 6
                    dia(1, j) = dia(1, j)+dt(1, k)*ta(k, j)
                    dia(2, j) = dia(2, j)+dt(2, k)*ta(k, j)
                end do
            end do
        end do
        do j = 1, 12
            aw(ic, j) = (c(ic)*dib(1, j)+s(ic)*dib(2, j))*l(ic)/2.d0
        end do
        do j = 1, 4
            aa(ic, j) = -(c(ic)*dia(1, j)+s(ic)*dia(2, j))*l(ic)/2.d0
        end do
        aa(ic, ic) = aa(ic, ic)+2.d0/3.d0*l(ic)
    end do
!     -------------- INVERSION DE AA -----------------------------------
    aai(:, :) = 0.d0
    do i = 1, 4
        aai(i, i) = 1.d0
    end do
    call mgauss('NFVP', aa, aai, 4, 4, &
                4, det, iret)
!
    aw(1, 1) = aw(1, 1)+1.d0
    aw(1, 2) = aw(1, 2)-x(1)/2.d0
    aw(1, 3) = aw(1, 3)-y(1)/2.d0
    aw(1, 4) = aw(1, 4)-1.d0
    aw(1, 5) = aw(1, 5)-x(1)/2.d0
    aw(1, 6) = aw(1, 6)-y(1)/2.d0
    aw(2, 4) = aw(2, 4)+1.d0
    aw(2, 5) = aw(2, 5)-x(2)/2.d0
    aw(2, 6) = aw(2, 6)-y(2)/2.d0
    aw(2, 7) = aw(2, 7)-1.d0
    aw(2, 8) = aw(2, 8)-x(2)/2.d0
    aw(2, 9) = aw(2, 9)-y(2)/2.d0
    aw(3, 7) = aw(3, 7)+1.d0
    aw(3, 8) = aw(3, 8)-x(3)/2.d0
    aw(3, 9) = aw(3, 9)-y(3)/2.d0
    aw(3, 10) = aw(3, 10)-1.d0
    aw(3, 11) = aw(3, 11)-x(3)/2.d0
    aw(3, 12) = aw(3, 12)-y(3)/2.d0
    aw(4, 1) = aw(4, 1)-1.d0
    aw(4, 2) = aw(4, 2)-x(4)/2.d0
    aw(4, 3) = aw(4, 3)-y(4)/2.d0
    aw(4, 10) = aw(4, 10)+1.d0
    aw(4, 11) = aw(4, 11)-x(4)/2.d0
    aw(4, 12) = aw(4, 12)-y(4)/2.d0
!
!     -------------- AN = AAI.AW ---------------------------------------
    an(:, :) = 0.d0
    do i = 1, 4
        do k = 1, 4
            do j = 1, 12
                an(i, j) = an(i, j)+aai(i, k)*aw(k, j)
            end do
        end do
    end do
!
end subroutine
