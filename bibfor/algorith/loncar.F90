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
subroutine loncar(ndim, typma, coord, l)
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
    integer(kind=8) :: ndim
    real(kind=8) :: coord(*), l
    character(len=8) :: typma
!
!
!                      LONGUEUR CARACTÉRISTIQUE D'UNE MAILLE
!
!     ENTREE
!       NDIM    : DIMENSION TOPOLOGIQUE DU MAILLAGE
!       TYPMA   : TYPE DE MAILLE (TYPE_MAILLE)
!       COORD   : COORDONNÉES DES NOEUDS (X Y Z SI NDIM = 3
!                                         X Y   SI NDIM = 2)
!
!     SORTIE
!       L      : LONGUEUR CARACTÉRISTIQUE
!     ------------------------------------------------------------------
!
    integer(kind=8) :: i
    real(kind=8) :: ar(3), m(3)
! ----------------------------------------------------------------------
!
!
!     ATTENTION,
!     NDIM EST LA DIMENSION DU MAILLAGE
!     POUR LES MAILLES DE BORD, CE N'EST PAS LA DIMENSION DE LA MAILLE
!
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
!
    if (typma(1:4) .eq. 'HEXA') then
!
!       LA LONGUEUR CARACTÉRISTIQUE EST LA GRANDE DIAGONALE N1-N7
        l = sqrt((coord(1)-coord(19))**2+(coord(2)-coord(20))**2 &
                 +(coord(3)-coord(21))**2)
!
    else if (typma(1:5) .eq. 'PENTA') then
!
!       LA LONGUEUR CARACTÉRISTIQUE EST ((N3-N1)*(N3-N2)*(N3-N6))^(1/3)
        ar(1) = sqrt((coord(7)-coord(1))**2+(coord(8)-coord(2))**2 &
                     +(coord(9)-coord(3))**2)
        ar(2) = sqrt((coord(7)-coord(4))**2+(coord(8)-coord(5))**2 &
                     +(coord(9)-coord(6))**2)
        ar(3) = sqrt((coord(7)-coord(16))**2+(coord(8)-coord(17))**2 &
                     +(coord(9)-coord(18))**2)
        l = (ar(1)*ar(2)*ar(3))**(1.d0/3.d0)
!
    else if (typma(1:5) .eq. 'PYRAM') then
!
!       M : MILIEU DE LA FACE QUADRANGLE
        do i = 1, 3
            m(i) = (coord(3*(1-1)+i)+coord(3*(2-1)+i)+coord(3*(3-1)+i)+coord(3*(4-1)+i) &
                    )/4.d0
        end do
!
!       LA LONGUEUR CARACTÉRISTIQUE EST ((N1-N3)*(M-N5))^(1/2)
        ar(1) = sqrt((coord(3*(3-1)+1)-coord(3*(1-1)+1))**2+(coord(3* &
                                       (3-1)+2)-coord(3*(1-1)+2))**2+(coord(3*(3-1)+3)-coord(3*(1- &
                                                                                          1)+3))**2)
        ar(2) = sqrt((m(1)-coord(3*(5-1)+1))**2+(m(2)-coord(3*(5-1) &
                                                            +2))**2+(m(3)-coord(3*(5-1)+3))**2)
        l = sqrt(ar(1)*ar(2))
!
    else if (typma(1:5) .eq. 'TETRA') then
!
!       LA LONGUEUR CARACTÉRISTIQUE EST ((N1-N2)*(N1-N3)*(N1-N4))^(1/3)
        do i = 1, 3
            ar(i) = sqrt((coord(1)-coord(3*i+1))**2+(coord(2)-coord(3* &
                                                                i+2))**2+(coord(3)-coord(3*i+3))**2)
        end do
        l = (ar(1)*ar(2)*ar(3))**(1.d0/3.d0)
!
    else if (typma(1:4) .eq. 'QUAD') then
!
!     LA LONGUEUR CARACTÉRISTIQUE EST ((N1-N2)*(N1-N3))^(1/2)
        do i = 1, 2
            ar(i) = (coord(1)-coord(ndim*i+1))**2+(coord(2)-coord(ndim*i+2))**2
            if (ndim .eq. 3) ar(i) = ar(i)+(coord(3)-coord(ndim*i+3))**2
        end do
        l = (sqrt(ar(1)*ar(2)))**(1.d0/2.d0)
!
    else if (typma(1:4) .eq. 'TRIA') then
!
!     LA LONGUEUR CARACTÉRISTIQUE EST ((N1-N2)*(N1-N3))^(1/2)
        do i = 1, 2
            ar(i) = (coord(1)-coord(ndim*i+1))**2+(coord(2)-coord(ndim*i+2))**2
            if (ndim .eq. 3) ar(i) = ar(i)+(coord(3)-coord(ndim*i+3))**2
        end do
        l = (sqrt(ar(1)*ar(2)))**(1.d0/2.d0)
!
    else if (typma(1:3) .eq. 'SEG') then
!
!       LA LONGUEUR CARACTÉRISTIQUE EST (N1-N2)^(1/2)
        l = (coord(1)-coord(ndim+1))**2+(coord(2)-coord(ndim+2))**2
        if (ndim .eq. 3) l = l+(coord(3)-coord(ndim+3))**2
        l = sqrt(l)
!
    else
!
!       TYPE D'ELEMENT FINI PAS TRAITE
        ASSERT(.false.)
!
    end if
!
end subroutine
