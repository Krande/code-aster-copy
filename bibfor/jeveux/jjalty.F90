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

subroutine jjalty(typei, ltypi, cel, inatb, jctab)
! aslint: disable=C1002,W0405
    implicit none
#include "jeveux.h"
#include "asterfort/jxveuo.h"
    integer(kind=8) :: ltypi, inatb, jctab
    character(len=*) :: typei, cel
!-----------------------------------------------------------------------
! ALLOUE LE SEGMENT DE VALEURS EN MEMOIRE ET LE POSITIONNE EN
! FONCTION DU TYPE ASSOCIE
!
! IN   TYPEI  : TYPE DE L'OBJET
! IN   LTYPI  : LONGUEUR DU TYPE
! IN   CEL    : 'L' OU 'E'
! IN   INATB  : TYPE D'OBJET 1:OS 2:CO 3:OC
! OUT  JCTAB  : ADRESSE PAR RAPPORT AU COMMUN DE REFERENCE
!
!-----------------------------------------------------------------------
    integer(kind=8) :: izr(1), izc(1), izl(1), izk8(1), izk16(1), izk24(1)
    integer(kind=8) :: izk32(1), izk80(1), izi4(1)
    equivalence(izr, zr), (izc, zc), (izl, zl), (izk8, zk8), (izk16, zk16), &
        (izk24, zk24), (izk32, zk32), (izk80, zk80), (izi4, zi4)
! DEB ------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    jctab = 0
    if (typei .eq. 'I') then
        call jxveuo(cel, zi, inatb, jctab)
    else if (typei .eq. 'S') then
        call jxveuo(cel, izi4, inatb, jctab)
    else if (typei .eq. 'R') then
        call jxveuo(cel, izr, inatb, jctab)
    else if (typei .eq. 'C') then
        call jxveuo(cel, izc, inatb, jctab)
    else if (typei .eq. 'K') then
        if (ltypi .eq. 8) then
            call jxveuo(cel, izk8, inatb, jctab)
        else if (ltypi .eq. 16) then
            call jxveuo(cel, izk16, inatb, jctab)
        else if (ltypi .eq. 24) then
            call jxveuo(cel, izk24, inatb, jctab)
        else if (ltypi .eq. 32) then
            call jxveuo(cel, izk32, inatb, jctab)
        else if (ltypi .eq. 80) then
            call jxveuo(cel, izk80, inatb, jctab)
        else
            call jxveuo(cel, izk8, inatb, jctab)
        end if
    else if (typei .eq. 'L') then
        call jxveuo(cel, izl, inatb, jctab)
    end if
! FIN ------------------------------------------------------------------
end subroutine
