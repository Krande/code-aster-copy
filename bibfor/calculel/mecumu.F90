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
subroutine mecumu(scal, ncmp, iad1, iad2, nec, &
                  dg1, dg2)
    implicit none
!
! ----------------------------------------------------------------------
!     "CUMULE" LES CMPS D'1 GRANDEUR AU SENS DE:
!     '   ' + ' 2' = ' 2 '
!     ' 3 ' + ' 2' = ' 3 '
!
! ----------------------------------------------------------------------
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/exisdg.h"
#include "asterfort/utmess.h"
    character(len=8) :: scal
    integer(kind=8) :: ncmp, iad1, iad2, nec, dg1(nec), dg2(nec)
! ----------------------------------------------------------------------
!     ENTREES:
!       SCAL : TYPE SCALAIRE : 'R  ', 'C  ', 'I  ', 'K8 ' ,'K16 '
!       NCMP : NOMBRE DE COMPOSANTES A ACCUMULER.
!       IAD1 : ADRESSE DANS ZI,ZR,... DU SEGMENT A CUMULER
!       IAD2 : ADRESSE DANS ZI,ZR,... DU SEGMENT OU ON CUMULE
!
!       NEC  : NOMBRE D'ENTIERS CODES
!       DG1  : DESCRIPTEUR DE LA GRANDEUR A CUMULER.
!       DG2  : DESCRIPTEUR  DE LA GRANDEUR OU ON CUMULE.
!
!     SORTIES:
!       LES OBJETS SONT MODIFIES.
! ----------------------------------------------------------------------
!
!     FONCTIONS EXTERNES:
!     -------------------
!
!     VARIABLES LOCALES:
!     ------------------
!
! DEB-------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ico
!-----------------------------------------------------------------------
    ico = 0
    if (scal(1:1) .eq. 'I') then
!
!        -- CAS D'1 SEGMENT ENTIER:
!        --------------------------
        do i = 1, ncmp
            if (exisdg(dg1, i)) then
                ico = ico+1
                zi(iad2-1+i) = zi(iad1-1+ico)
            end if
        end do
    else if (scal(1:1) .eq. 'R') then
!
!        -- CAS D'1 SEGMENT REEL  :
!        --------------------------
        do i = 1, ncmp
            if (exisdg(dg1, i)) then
                ico = ico+1
                zr(iad2-1+i) = zr(iad1-1+ico)
            end if
        end do
!
!        -- CAS D'1 SEGMENT COMPLEX:
!        --------------------------
    else if (scal(1:1) .eq. 'C') then
        do i = 1, ncmp
            if (exisdg(dg1, i)) then
                ico = ico+1
                zc(iad2-1+i) = zc(iad1-1+ico)
            end if
        end do
!
!        -- CAS D'1 SEGMENT DE CARACTERES (K8):
!        ---------------------------------
    else if (scal(1:3) .eq. 'K8 ') then
        do i = 1, ncmp
            if (exisdg(dg1, i)) then
                ico = ico+1
                zk8(iad2-1+i) = zk8(iad1-1+ico)
            end if
        end do
!
!        -- CAS D'1 SEGMENT DE CARACTERES (K16):
!        ---------------------------------
    else if (scal(1:3) .eq. 'K16') then
        do i = 1, ncmp
            if (exisdg(dg1, i)) then
                ico = ico+1
                zk16(iad2-1+i) = zk16(iad1-1+ico)
            end if
        end do
!
!        -- CAS D'1 SEGMENT DE CARACTERES (K24):
!        ---------------------------------
    else if (scal(1:3) .eq. 'K24') then
        do i = 1, ncmp
            if (exisdg(dg1, i)) then
                ico = ico+1
                zk24(iad2-1+i) = zk24(iad1-1+ico)
            end if
        end do
    else
        call utmess('F', 'CALCULEL3_38', sk=scal(1:4))
    end if
!
    do i = 1, nec
        dg2(i) = ior(dg2(i), dg1(i))
    end do
!
end subroutine
