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

subroutine pjeflo(elrefa, ndim, ipb, xr2, disprj)
!
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "jeveux.h"
!
#include "asterfort/assert.h"
!
    integer(kind=8) :: ipb, ndim
    real(kind=8) :: xr2(ndim), disprj
    character(len=*) :: elrefa
! ----------------------------------------------------------------------
! BUT :
!   * calculer disprj : distance de projection d'un point
!                       (normee par la "taille" de la maille)
!  disprj =   0. => le point est interieur a la maille
!  disprj = 999. => la routine reereg n'a pas converge
!  disprj = a>0  => le point est exterieur a la maille.
! ----------------------------------------------------------------------
! in  elrefa   : elrefa de l'element
! in  ndim     : dimension de l'espace
! in  xr2      : coordonnees du point dans l'element de reference
!                (calcule par reereg)
! in  ipbd     : code retour de reereg
! out disprj   : distance de projection (en relatif)
! ----------------------------------------------------------------------
!
    real(kind=8) :: x, y, z, diam
! --------------------------------------------------------------------------------------------------
    disprj = 0.0d0
!   si reereg n'a pas converge, on n'a pas confiance dans xr2 :
    if (ipb .ne. 0) then
        disprj = dble(999)
        goto 80
    end if
!
    if (ndim .ge. 1) x = xr2(1)
    if (ndim .ge. 2) y = xr2(2)
    if (ndim .ge. 3) z = xr2(3)
!
! --------------------------------------------------------------------------------------------------
!   POUR LES HEXA : KSI,ETA,DZETA SONT DANS [-1,1]
    if (elrefa .eq. 'HE8' .or. elrefa .eq. 'H20' .or. elrefa .eq. 'H27' .or. elrefa .eq. 'HE9') then
        ASSERT(ndim .eq. 3)
        if (abs(x) .gt. 1.d0) goto 10
        if (abs(y) .gt. 1.d0) goto 10
        if (abs(z) .gt. 1.d0) goto 10
!       ON EST INTERIEUR
        goto 80
!
10      continue
!       ON EST EXTERIEUR. EST-ON LOIN ?
        disprj = 0.d0
        disprj = max(disprj, abs(x)-1.d0)
        disprj = max(disprj, abs(y)-1.d0)
        disprj = max(disprj, abs(z)-1.d0)
!       -- diam : "dimension" de l'elrefe :
        diam = 2.
        disprj = disprj/diam
! --------------------------------------------------------------------------------------------------
!   POUR LES TETRA :
    elseif (elrefa .eq. 'TE4' .or. elrefa .eq. 'T10') then
        ASSERT(ndim .eq. 3)
        if (x .lt. 0.d0) goto 20
        if (y .lt. 0.d0) goto 20
        if (z .lt. 0.d0) goto 20
        if (x+y+z .gt. 1.d0) goto 20
!
!       ON EST INTERIEUR
        goto 80
!
20      continue
!       ON EST EXTERIEUR. EST-ON LOIN ?
        disprj = 0.d0
        disprj = max(disprj, -x)
        disprj = max(disprj, -y)
        disprj = max(disprj, -z)
        disprj = max(disprj, x+y+z-1.d0)
!       -- diam : "dimension" de l'elrefe :
        diam = 1.
        disprj = disprj/diam
!
! --------------------------------------------------------------------------------------------------
!   POUR LES PYRAM :
    elseif (elrefa .eq. 'PY5' .or. elrefa .eq. 'P13') then
        ASSERT(ndim .eq. 3)
        if (z .lt. 0.d0) goto 30
        if (x+y+z .gt. 1.d0) goto 30
        if (x-y+z .gt. 1.d0) goto 30
        if (-x+y+z .gt. 1.d0) goto 30
        if (-x-y+z .gt. 1.d0) goto 30
!
!       ON EST INTERIEUR
        goto 80
!
30      continue
!       ON EST EXTERIEUR. EST-ON LOIN ?
        disprj = 0.d0
        disprj = max(disprj, -z)
        disprj = max(disprj, x+y+z-1.d0)
        disprj = max(disprj, x-y+z-1.d0)
        disprj = max(disprj, -x+y+z-1.d0)
        disprj = max(disprj, -x-y+z-1.d0)
!       -- diam : "dimension" de l'elrefe :
        diam = 2.
        disprj = disprj/diam
!
! --------------------------------------------------------------------------------------------------
!   POUR LES PENTA :
    elseif (elrefa .eq. 'PE6' .or. elrefa .eq. 'P15' .or. elrefa .eq. 'P18') then
        ASSERT(ndim .eq. 3)
        if (x .lt. -1.d0) goto 40
        if (x .gt. +1.d0) goto 40
        if (y .lt. 0.d0) goto 40
        if (z .lt. 0.d0) goto 40
        if (y+z .gt. 1.d0) goto 40
!
!       ON EST INTERIEUR
        goto 80
!
40      continue
!       ON EST EXTERIEUR. EST-ON LOIN ?
        disprj = 0.d0
        disprj = max(disprj, abs(x)-1.d0)
        disprj = max(disprj, -y)
        disprj = max(disprj, -z)
        disprj = max(disprj, +y+z-1.d0)
!       -- diam : "dimension" de l'elrefe :
        diam = 2.
        disprj = disprj/diam
!
! --------------------------------------------------------------------------------------------------
!   POUR LES TRIA :
    elseif (elrefa .eq. 'TR3' .or. elrefa .eq. 'TR6' .or. elrefa .eq. 'TR7') then
        ASSERT(ndim .eq. 2)
        if (x .lt. 0.d0) goto 50
        if (y .lt. 0.d0) goto 50
        if (x+y .gt. 1.d0) goto 50
!
!       ON EST INTERIEUR
        goto 80
!
50      continue
!       ON EST EXTERIEUR. EST-ON LOIN ?
        disprj = 0.d0
        disprj = max(disprj, -x)
        disprj = max(disprj, -y)
        disprj = max(disprj, +x+y-1.d0)
!       -- diam : "dimension" de l'elrefe :
        diam = 1.
        disprj = disprj/diam
!
! --------------------------------------------------------------------------------------------------
!   POUR LES QUAD :
    elseif (elrefa .eq. 'QU4' .or. elrefa .eq. 'QU8' .or. elrefa .eq. 'QU9') then
        ASSERT(ndim .eq. 2)
        if (x .lt. -1.d0) goto 60
        if (y .lt. -1.d0) goto 60
        if (x .gt. +1.d0) goto 60
        if (y .gt. +1.d0) goto 60
!
!       ON EST INTERIEUR
        goto 80
!
60      continue
!       ON EST EXTERIEUR. EST-ON LOIN ?
        disprj = 0.d0
        disprj = max(disprj, -1.d0-x)
        disprj = max(disprj, -1.d0-y)
        disprj = max(disprj, x-1.d0)
        disprj = max(disprj, y-1.d0)
!       -- diam : "dimension" de l'elrefe :
        diam = 2.
        disprj = disprj/diam
    else
        ASSERT(.false.)
    end if
!
80  continue
end subroutine
