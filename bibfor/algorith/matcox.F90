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
subroutine matcox(ndim, pp, ddt1, ddt2, ddt3, &
                  ddt4, p, nno, ddlh, ddls, &
                  jac, ffp, singu, fk, mmat)
    implicit none
#include "asterfort/assert.h"
    real(kind=8) :: ddt1(3, 3), ddt2(3, 3), ddt3(3, 3), ddt4(3, 3), pp(3, 3)
    real(kind=8) :: p(3, 3), mmat(216, 216)
    real(kind=8) :: jac, ffp(27)
    real(kind=8) :: fk(27, 3, 3)
    integer :: ndim, ddlh, ddls, nno, singu
!.......................................................................
!
!         CALCUL DES MATRICES DE CONTACT FROTTEMENT
!                 LOI COHESIVE - POUR X-FEM
!                     (METHODE CONTINUE)
!
!
!
!  ENTREES  --->  PP,P,JAC,FFP,FK, NDIM,DDLH,DDLS,NNO
!           --->  SINGU,DDT1,DDT2,DDT3,DDT4
!  SORTIES  --->  MMAT
!
!.......................................................................
    real(kind=8) :: ddt11(3, 3), ddt21(3, 3), ddt31(3, 3), ddt41(3, 3)
    real(kind=8) :: ddt111(3, 3), ddt211(3, 3), ddt311(3, 3), ddt411(3, 3)
    integer :: i, j, k, l, alpi, alpj
!
!.......................................................................
!
!   D'APRES LE DIMENSIONNMENT DE DDT* LE MULTI-HEAVISIDE EST EXCLUS
!   PAR PRECAUTION ON MET UN ASSERT
    ASSERT(ddlh .eq. ndim)
!
    ddt11(:, :) = 0.d0
    ddt21(:, :) = 0.d0
    ddt31(:, :) = 0.d0
    ddt41(:, :) = 0.d0
!
    ddt111(:, :) = 0.d0
    ddt211(:, :) = 0.d0
    ddt311(:, :) = 0.d0
    ddt411(:, :) = 0.d0
!
    do i = 1, ndim
        do j = 1, ndim
            do l = 1, ndim
                ddt11(i, j) = ddt11(i, j)+pp(i, l)*ddt1(l, j)
                ddt21(i, j) = ddt21(i, j)+pp(i, l)*ddt2(l, j)
                ddt31(i, j) = ddt31(i, j)+p(i, l)*ddt3(l, j)
                ddt41(i, j) = ddt41(i, j)+p(i, l)*ddt4(l, j)
            end do
        end do
    end do
!
    do i = 1, ndim
        do j = 1, ndim
            do l = 1, ndim
                ddt111(i, j) = ddt111(i, j)+ddt11(i, l)*pp(l, j)
                ddt211(i, j) = ddt211(i, j)+ddt21(i, l)*p(l, j)
                ddt311(i, j) = ddt311(i, j)+ddt31(i, l)*pp(l, j)
                ddt411(i, j) = ddt411(i, j)+ddt41(i, l)*p(l, j)
            end do
        end do
    end do
!
    do i = 1, nno
        do j = 1, nno
            do k = 1, ddlh
                do l = 1, ddlh
!
                    mmat(ddls*(i-1)+ndim+k, ddls*(j-1)+ndim+l) = &
                        mmat(ddls*(i-1)+ndim+k, ddls*(j-1)+ndim+l)+ &
                        4.d0*ffp(i)*ddt111(k, l)*ffp(j)*jac+4.d0*ffp(i)* &
                        ddt211(k, l)*ffp(j)*jac+4.d0*ffp(i)*ddt311(k, l)* &
                        ffp(j)*jac+4.d0*ffp(i)*ddt411(k, l)*ffp(j)*jac
!
                end do
!
                do alpj = 1, singu*ndim
                    do l = 1, ndim
                        mmat(ddls*(i-1)+ndim+k, ddls*(j-1)+ndim+ddlh+alpj) = &
                            mmat(ddls*(i-1)+ndim+k, ddls*(j-1)+ndim+ddlh+alpj)+ &
                            4.d0*ffp(i)*ddt111(k, l)*jac*fk(j, alpj, l)+4.d0*ffp(i) &
                            *ddt211(k, l)*jac**fk(j, alpj, l)+4.d0*ffp(i)*ddt311(k, &
                                                   l)*jac*fk(j, alpj, l)+4.d0*ffp(i)*ddt411(k, l)* &
                            jac*fk(j, alpj, l)
                    end do
                end do
!
            end do
            do alpi = 1, singu*ndim
                do k = 1, ndim
                    do l = 1, ddlh
                        mmat(ddls*(i-1)+ndim+ddlh+alpi, ddls*(j-1)+ndim+l) = &
                            mmat(ddls*(i-1)+ndim+ddlh+alpi, ddls*(j-1)+ndim+l)+ &
                            4.d0*ddt111(k, l)*ffp(j)*jac*fk(i, alpi, k)+4.d0 &
                            *ddt211(k, l)*ffp(j)*jac*fk(i, alpi, k)+4.d0*ddt311(k, &
                                            l)*ffp(j)*jac*fk(i, alpi, k)+4.d0*ddt411(k, l)*ffp(j)* &
                            jac*fk(i, alpi, k)
                    end do
!
                    do alpj = 1, singu*ndim
                        do l = 1, ndim
                            mmat(ddls*(i-1)+ndim+ddlh+alpi, ddls*(j-1)+ndim+ddlh+ &
                                 alpj) = mmat(ddls*(i-1)+ndim+ddlh+alpi, ddls*(j-1)+ndim+ &
                                    ddlh+alpj)+4.d0*ddt111(k, l)*jac*fk(i, alpi, k)*fk(j, alpj, l) &
                                    +4.d0*ddt211(k, l)*jac*fk(i, alpi, k)*fk(j, alpj, l)+4.d0* &
                                    ddt311(k, l)*jac*fk(i, alpi, k)*fk(j, alpj, l)+4.d0* &
                                    ddt411(k, l)*jac*fk(i, alpi, k)*fk(j, alpj, l)
                        end do
                    end do
                end do
            end do
!
!
        end do
    end do
!
end subroutine
