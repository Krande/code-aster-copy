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
subroutine q4gniw(qsi, eta, caraq4, dba, wsq)
    implicit none
    real(kind=8) :: qsi, eta, caraq4(*), dba(2, 12), wsq(12)
!     INTERPOLATION DE LA FLECHE AU POINTS D'INTEGRATION POUR LE Q4G
!     ------------------------------------------------------------------
    integer(kind=8) :: j
    real(kind=8) :: pqsi, mqsi, peta, meta, qsic, etac
    real(kind=8) :: x5, x6, x7, x8, y5, y6, y7, y8
    real(kind=8) :: n(12)
!     ------------------------------------------------------------------
    x5 = caraq4(1)
    x6 = caraq4(2)
    x7 = caraq4(3)
    x8 = caraq4(4)
    y5 = caraq4(5)
    y6 = caraq4(6)
    y7 = caraq4(7)
    y8 = caraq4(8)
!
    peta = 1.d0+eta
    meta = 1.d0-eta
    pqsi = 1.d0+qsi
    mqsi = 1.d0-qsi
    etac = 1.d0-eta*eta
    qsic = 1.d0-qsi*qsi
!
!     --- FONCTIONS D'INTERPOLATION DU DKQ DANS LE REPERE REDUIT -------
    n(1) = mqsi*meta/8.d0*(qsic+etac-qsi-eta)
    n(2) = mqsi*meta/8.d0*qsic
    n(3) = mqsi*meta/8.d0*etac
    n(4) = pqsi*meta/8.d0*(qsic+etac+qsi-eta)
    n(5) = -pqsi*meta/8.d0*qsic
    n(6) = pqsi*meta/8.d0*etac
    n(7) = pqsi*peta/8.d0*(qsic+etac+qsi+eta)
    n(8) = -pqsi*peta/8.d0*qsic
    n(9) = -pqsi*peta/8.d0*etac
    n(10) = mqsi*peta/8.d0*(qsic+etac-qsi+eta)
    n(11) = mqsi*peta/8.d0*qsic
    n(12) = -mqsi*peta/8.d0*etac
!
!     --- RAPPEL DE LA DISTORSION ( DBA ) ---------------
!
!     --- FONCTIONS D'INTERPOLATION DU DSQ DANS LE REPERE LOCAL --------
    do j = 1, 12
        wsq(j) = ( &
                 -dba(1, j)*x5+dba(2, j)*x8)*n(2)+(-dba(1, j)*y5+dba(2, j)*y8)*n(2)+(-&
                 & dba(1, j)*x5-dba(2, j)*x6)*n(5)+(-dba(1, j)*y5-dba(2, j)*y6)*n(5)+( &
                 &dba(1, j)*x7-dba(2, j)*x6)*n(8)+(dba(1, j)*y7-dba(2, j)*y6)*n(8)+(dba&
                 &(1, j)*x7+dba(2, j)*x8)*n(11)+(dba(1, j)*y7+dba(2, j)*y8)*n(11 &
                 )
    end do
!
! ----- FONCTIONS D'INTERPOLATION DANS LE REPERE LOCAL -------------
    wsq(1) = wsq(1)+n(1)
    wsq(2) = wsq(2)+(-x5*n(2)+x8*n(3))/2.d0
    wsq(3) = wsq(3)+(-y5*n(2)+y8*n(3))/2.d0
    wsq(4) = wsq(4)+n(4)
    wsq(5) = wsq(5)+(-x5*n(5)-x6*n(6))/2.d0
    wsq(6) = wsq(6)+(-y5*n(5)-y6*n(6))/2.d0
    wsq(7) = wsq(7)+n(7)
    wsq(8) = wsq(8)+(x7*n(8)-x6*n(9))/2.d0
    wsq(9) = wsq(9)+(y7*n(8)-y6*n(9))/2.d0
    wsq(10) = wsq(10)+n(10)
    wsq(11) = wsq(11)+(x7*n(11)+x8*n(12))/2.d0
    wsq(12) = wsq(12)+(y7*n(11)+y8*n(12))/2.d0
end subroutine
