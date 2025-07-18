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

subroutine ejinit(nomte, iu, ip)
!
! person_in_charge: jerome.laverne at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
    character(len=16) :: nomte
    integer(kind=8) :: iu(3, 16), ip(8)
! ----------------------------------------------------------------------
!            DECALAGE D'INDICE POUR LES ELEMENTS DE JOINT HM
! ----------------------------------------------------------------------
! IN  NOMTE  NOM DE L'ELEMENT FINI
! OUT IU     DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT
! OUT IP     DECALAGE D'INDICE POUR ACCEDER AUX DDL DE PRESSION
! ----------------------------------------------------------------------
!
! EXEMPLE POUR QUAD8 (6 NOEUDS DEPL + 2 NOEUDS PRESS)
!     NUMEROTATION DE DDL :
!     U1_X, U1_Y   - IU(1,1) IU(2,1)
!     U2_X, U2_Y   - IU(1,2) IU(2,2)
!     U3_X, U3_Y   - IU(1,5) IU(2,5)
!     U4_X, U4_Y   - IU(1,4) IU(2,4)
!     U5_X, U5_Y   - IU(1,3) IU(2,3)
!     P6           - IP(2)
!     U7_X, U7_Y   - IU(1,6) IU(2,6)
!     P8           - IP(1)
!
!     RACCOURCIS
!          IU(1,1:3) => DEPL_X JOINT + (NOEUDS 1,2 et 5)
!          IU(2,1:3) => DEPL_Y JOINT + (NOEUDS 1,2 et 5)
!          IU(1,4:6) => DEPL_X JOINT - (NOEUDS 4,3 et 7)
!          IU(2,4:6) => DEPL_Y JOINT - (NOEUDS 4,3 et 7)
!          IP(1:2)   => PRESS (NOEUDS 8 et 6)
! ----------------------------------------------------------------------
!
    integer(kind=8) :: n, offset
    integer(kind=8) :: uh20(16), ph20(8)
    integer(kind=8) :: up15(12), pp15(6)
    integer(kind=8) :: uq8(6), pq8(3)
    aster_logical :: ifqu8, ifh20, ifp15
! ----------------------------------------------------------------------
    data uh20/1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 17, 18, 19, 20/
    data ph20/13, 14, 15, 16, 17, 18, 19, 20/
!
    data up15/1, 2, 3, 7, 8, 9, 4, 5, 6, 13, 14, 15/
    data pp15/10, 11, 12, 13, 14, 15/
!
    data uq8/1, 2, 5, 4, 3, 7/
    data pq8/8, 6, 7/
! ----------------------------------------------------------------------
!     INDICATEURS DE TYPE DE MAILLE : QUAD8, PENTA15 ET HEXA20
!
    ifqu8 = (nomte .eq. 'EJHYME_PLQU8') .or. (nomte .eq. 'EJHYME_AXQU8') .or. (nomte .eq. 'MFPLQU8')
    ifp15 = (nomte .eq. 'EJHYME_PENTA15') .or. (nomte .eq. 'MEFI_PENTA15')
    ifh20 = (nomte .eq. 'EJHYME_HEXA20') .or. (nomte .eq. 'MEFI_HEXA20')
!
!CCCCCCCCCCCC HEXA20 CCCCCCCCCCCCCCCC
    if (ifh20) then
!       DDL U_XYZ POUR NOEUDS 1,2,3,4,9,10,11,12,5,6,7,8
        do n = 1, 12
            iu(1, n) = 1+(uh20(n)-1)*3
            iu(2, n) = 2+(uh20(n)-1)*3
            iu(3, n) = 3+(uh20(n)-1)*3
        end do
!       DDL PRESS POUR NOEUDS 16,13,14,15
        offset = 12*3
        do n = 1, 4
            ip(n) = 1+ph20(n)-13+offset
        end do
!       DDL U_XYZ POUR NOEUDS 17,18,19,20
        offset = 12*3+4*1
        do n = 13, 16
            iu(1, n) = 1+(uh20(n)-17)*4+offset
            iu(2, n) = 2+(uh20(n)-17)*4+offset
            iu(3, n) = 3+(uh20(n)-17)*4+offset
        end do
!       DDL PRESS POUR NOEUDS 17,18,19,20
        offset = 12*3+4*1
        do n = 5, 8
            ip(n) = 4+(ph20(n)-17)*4+offset
        end do
!
!CCCCCCCCCCCC PENTA15 CCCCCCCCCCCCCCCC
    else if (ifp15) then
!       DDL U_XYZ POUR NOEUDS 1,2,3,7,8,9,4,5,6
        do n = 1, 9
            iu(1, n) = 1+(up15(n)-1)*3
            iu(2, n) = 2+(up15(n)-1)*3
            iu(3, n) = 3+(up15(n)-1)*3
        end do
!       DDL PRESS POUR NOEUDS 12,10,11
        offset = 9*3
        do n = 1, 3
            ip(n) = 1+pp15(n)-10+offset
        end do
!       DDL U_XYZ POUR NOEUDS 13,14,15
        offset = 9*3+3*1
        do n = 10, 12
            iu(1, n) = 1+(up15(n)-13)*4+offset
            iu(2, n) = 2+(up15(n)-13)*4+offset
            iu(3, n) = 3+(up15(n)-13)*4+offset
        end do
!       DDL U_XYZ POUR NOEUDS 13,14,15
        offset = 9*3+3*1
        do n = 4, 6
            ip(n) = 4+(pp15(n)-13)*4+offset
        end do
!
!CCCCCCCCCCCC QUAD8 CCCCCCCCCCCCCCCC
    else if (ifqu8) then
!ccccccccccccccccccccccccccccccccccc
!
        do n = 1, 5
            iu(1, n) = 1+(uq8(n)-1)*2
            iu(2, n) = 2+(uq8(n)-1)*2
        end do
        iu(1, 6) = 1+(uq8(6)-1)*2-1
        iu(2, 6) = 2+(uq8(6)-1)*2-1
!
! 15
        ip(1) = 1+(pq8(1)-1)*2
! 11
        ip(2) = 1+(pq8(2)-1)*2
! 14
        ip(3) = 2+(pq8(3)-1)*2
!ccccccccccccccccccccccccccccccccccccc
    else
!     NOM D'ELEMENT ILLICITE
        ASSERT(ifqu8 .or. ifp15 .or. ifh20)
    end if
!
end subroutine
