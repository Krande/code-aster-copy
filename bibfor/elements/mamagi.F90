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
subroutine mamagi(nomte, xr, yr)
    implicit none
#include "asterfort/mgauss.h"
    character(len=8) :: nomte
    real(kind=8) :: xr(*), yr(*), det
    real(kind=8) :: psi3(3), bt(12, 27), btb(12, 12)
    real(kind=8) :: atild(3), btild(3), ctild(3)
    real(kind=8) :: coefa(4), coefb(4), coefc(4), coefd(4)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i2, ifon, ig, intef, integ
    integer(kind=8) :: intsn, iret, j, jf, k, kp, l
    integer(kind=8) :: npge, npgsn, nso
    real(kind=8) :: bij, xi1, xi2, xi3
!-----------------------------------------------------------------------
    npge = 3
!
!--- ON PREND LES 3 POINTS D'INTEGRATION DANS L'EPAISSEUR. METHODE DE
!--- NEWTON-COTES
!
    psi3(1) = -1.d0
    psi3(2) = 0.d0
    psi3(3) = 1.d0
!
!-- COEFFICIENTS DES POLYNOMES D'ORDRE 2 D'INTERPOLATION DU TYPE
!     AXI3*XI3 + B*XI3 + C
!
    atild(1) = 0.5d0
    btild(1) = -0.5d0
    ctild(1) = 0.d0
!
    atild(2) = -1.d0
    btild(2) = 0.d0
    ctild(2) = 1.d0
!
    atild(3) = 0.5d0
    btild(3) = 0.5d0
    ctild(3) = 0.d0
!
!     LES 4 OU 3 FONCTIONS LINEAIRES AUX NOEUDS SOMMETS SONT DU TYPE
!     A + B*XI1 + C*XI2 + D*XI1*XI2
!
    if (nomte(1:8) .eq. 'MEC3QU9H' .or. nomte(1:8) .eq. 'MEGRC3Q9') then
!
        nso = 4
        npgsn = 9
!
        coefa(1) = 0.25d0
        coefb(1) = -0.25d0
        coefc(1) = -0.25d0
        coefd(1) = 0.25d0
!
        coefa(2) = 0.25d0
        coefb(2) = 0.25d0
        coefc(2) = -0.25d0
        coefd(2) = -0.25d0
!
        coefa(3) = 0.25d0
        coefb(3) = 0.25d0
        coefc(3) = 0.25d0
        coefd(3) = 0.25d0
!
        coefa(4) = 0.25d0
        coefb(4) = -0.25d0
        coefc(4) = 0.25d0
        coefd(4) = -0.25d0
!
    elseif (nomte(1:8) .eq. 'MEC3TR7H' .or. nomte(1:8) .eq. &
            'MEGRC3T7') then
!
        nso = 3
        npgsn = 7
!
        coefa(1) = 1.d0
        coefb(1) = -1.d0
        coefc(1) = -1.d0
        coefd(1) = 0.d0
!
        coefa(2) = 0.d0
        coefb(2) = 1.d0
        coefc(2) = 0.d0
        coefd(2) = 0.d0
!
        coefa(3) = 0.d0
        coefb(3) = 0.d0
        coefc(3) = 1.d0
        coefd(3) = 0.d0
!
    end if
!
!     CREATION DE LA MATRICE BT(NPGE*NSO,NPGE*NPGSN)
!
    ig = 0
    do integ = 1, npge
        xi3 = psi3(integ)
        do intsn = 1, npgsn
            i1 = 108+intsn
            i2 = 108+9+intsn
            xi1 = xr(i1)
            xi2 = xr(i2)
            ig = ig+1
            do intef = 1, npge
                do ifon = 1, nso
                    jf = nso*(intef-1)+ifon
                    bij = (atild(intef)*xi3**2+btild(intef)*xi3+ctild( &
                           intef))*(coefa(ifon)+coefb(ifon)*xi1+coefc(ifon)* &
                                    xi2+coefd(ifon)*xi1*xi2)
                    bt(jf, ig) = bij
                end do
            end do
        end do
    end do
!
    do i = 1, npge*nso
        do j = 1, npge*nso
            btb(i, j) = 0.d0
            do k = 1, npge*npgsn
                btb(i, j) = btb(i, j)+bt(i, k)*bt(j, k)
            end do
        end do
    end do
!
!     MATRICE DE PASSAGE = (BT*B*B)-1*BT
!
    call mgauss('NFVP', btb, bt, 12, npge*nso, &
                npge*npgsn, det, iret)
!
    do i = 1, npge*nso
        l = npge*npgsn*(i-1)
        do kp = 1, npge*npgsn
            yr(l+kp) = bt(i, kp)
        end do
    end do
!
end subroutine
