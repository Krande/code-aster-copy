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
subroutine mefgri(ntypg, nbgtot, zg, hg, itypg, &
                  zmin, zmax)
    implicit none
#include "asterf_types.h"
#include "asterfort/utmess.h"
!-----------------------------------------------------------------------
!     APPELANT : FLUST3
!     VERIFICATION DE LA REPARTITION GEOMETRIQUE DES GRILLES
!-----------------------------------------------------------------------
!  IN   : NTYPG  : NOMBRE DE TYPES DE GRILLES
!  IN   : NBGTOT : NOMBRE TOTAL DE GRILLES
!  IN   : ZG     : COORDONNEES 'Z' DES POSITIONS DES GRILLES DANS LE
!                  REPERE AXIAL
!  IN   : HG     : VECTEUR DES HAUTEURS DE GRILLE
!  IN   : ITYPG  : VECTEUR DES TYPES DE GRILLES
!  IN   : ZMIN   : COTE MIN DU FAISCEAU DE TUBES
!  IN   : ZMAX   : COTE MAX DU FAISCEAU DE TUBES
!-----------------------------------------------------------------------
    integer(kind=8) :: ntypg, nbgtot, itypg(nbgtot)
    real(kind=8) :: zg(nbgtot), hg(ntypg), zmin, zmax
!
    character(len=3) :: k3ig, k3jg
    character(len=24) :: valk(2)
    aster_logical :: intnul
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ig, jg
    real(kind=8) :: z1, z1ig, z1jg, z2, z2ig, z2jg
!-----------------------------------------------------------------------
    do ig = 1, nbgtot
        z1 = zg(ig)-hg(itypg(ig))/2.0d0
        z2 = zg(ig)+hg(itypg(ig))/2.0d0
        if ((z1 .lt. zmin) .or. (z2 .gt. zmax)) then
            write (k3ig, '(I3.3)') ig
            call utmess('F', 'ALGELINE_83', sk=k3ig)
        end if
    end do
!
    if (nbgtot .gt. 1) then
        do ig = 1, nbgtot-1
            z1ig = zg(ig)-hg(itypg(ig))/2.0d0
            z2ig = zg(ig)+hg(itypg(ig))/2.0d0
            do jg = ig+1, nbgtot
                z1jg = zg(jg)-hg(itypg(jg))/2.0d0
                z2jg = zg(jg)+hg(itypg(jg))/2.0d0
                intnul = ((z2ig .lt. z1jg) .or. (z2jg .lt. z1ig))
                if (.not. intnul) then
                    write (k3ig, '(I3.3)') ig
                    write (k3jg, '(I3.3)') jg
                    valk(1) = k3ig
                    valk(2) = k3jg
                    call utmess('F', 'ALGELINE_84', nk=2, valk=valk)
                end if
            end do
        end do
    end if
!
! --- FIN DE MEFGRI.
end subroutine
