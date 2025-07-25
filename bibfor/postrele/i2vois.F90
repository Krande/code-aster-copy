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
subroutine i2vois(conec, type, maille, n, v1, &
                  v2)
    implicit none
!
!
!**********************************************************************
!
!     RECHERCHE DES VOISINS (AU PLUS 2 PAR HYPOTHESE) DES MAILLES
!     D' UN ENSEMBLE DE MAILLES
!
!       CONEC  (IN)  : NOM DE L' OBJET CONNECTIVITE DU MAILLAGE
!
!       TYPE   (IN)  : NOM DE L' OBJET CONTENANT LES TYPES DES MAILLES
!
!       MAILLE (IN)  : TABLEAU DES NUMERO DE MAILLES
!
!       N      (IN)  : NOMBRE DE MAILLES DE L' ENSEMBLE TRAITE
!
!       V1     (OUT) : TABLEAU D' ACCES AUX VOISINS NUMERO 1
!
!       V2     (OUT) : TABLEAU D' ACCES AUX VOISINS NUMERO 2
!
!
!     EXPLICATION DE LA STRUCTURE DE DONNEES
!
!        MAILLE (V1(I)) <---VOISIN NUMERO 1 DE MAILLE(I)
!
!        MAILLE (V2(I)) <---VOISIN NUMERO 2 DE MAILLE(I)
!
!**********************************************************************
!
#include "asterf_types.h"
#include "asterfort/i2extf.h"
    integer(kind=8) :: n, v1(*), v2(*), maille(*)
    character(len=24) :: conec, type
!
    integer(kind=8) :: i, j, mi, mj
    integer(kind=8) :: nig, njg
    integer(kind=8) :: nid, njd
    aster_logical :: nonv1, nonv2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    i = 0
    j = 0
!
    mi = 0
    mj = 0
!
    nig = 0
    nig = 0
    njg = 0
    njg = 0
    nid = 0
    nid = 0
    njd = 0
    njd = 0
!
    nonv1 = .true.
    nonv2 = .true.
!
    do i = 1, n, 1
!
        mi = maille(i)
!
        call i2extf(mi, 1, conec(1:15), type(1:16), nig, &
                    nid)
!
        nonv1 = .true.
        nonv2 = .true.
!
        j = 0
!
20      continue
        if ((nonv1 .or. nonv2) .and. (j .lt. n)) then
!
            j = j+1
!
            if (i .ne. j) then
!
                mj = maille(j)
!
                call i2extf(mj, 1, conec(1:15), type(1:16), njg, &
                            njd)
!
                if ((nig .eq. njg) .or. (nid .eq. njg) .or. (nig .eq. njd) .or. &
                    (nid .eq. njd)) then
!
                    if (nonv1) then
!
                        nonv1 = .false.
!
                        v1(i) = j
!
                    else if (nonv2) then
!
                        nonv2 = .false.
!
                        v2(i) = j
!
                    else
!
                    end if
!
                end if
!
            end if
!
            goto 20
!
        end if
!
        if (nonv1) then
!
            v1(i) = 0
!
        end if
!
        if (nonv2) then
!
            v2(i) = 0
!
        end if
!
    end do
!
end subroutine
