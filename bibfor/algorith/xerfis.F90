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

subroutine xerfis(ndime, ninter, npts, nptm)
    implicit none
!
#include "asterfort/utmess.h"
    integer(kind=8) :: ndime, ninter, npts, nptm
!
!
!                 AFFICHER DES MESSAGES D'ERREUR LORSQUE LES
!                CONFIGURATIONS DE FISSURE SUR L'ELEMENT SONT :
!       (1) DOUTEUSES (AURAIT DU ETRE RETIRE)/ IMPROBABLES (FAUX)
!       (2) INTERDITES
!     ENTREE
!       NDIME   : DIMENSION DE L'ELEMENT FINI
!       NINTER  : NOMBRE DE POINTS D'INTERSECTION
!       NPTS    : NB DE POINTS INTERSEPTES EN UN DES 3 NOEUDS SOMMETS
!       NPTM    : NB DE POINT INTERSEPTES EN UN DES 3 NOEUDS MILIEUX
!......................................................................
!
!
! --- POUR LES TETRA10
!
    if (ndime .eq. 3) then
!
        if (ninter .eq. 3 .and. npts .eq. 2) then
            call utmess('F', 'XFEM_64')
        end if
!
! --- POUR LES TRIA6
!
    else if (ndime .eq. 2) then
!
! PLUTOT QUE DE VERIFIER LES SEULES CONFIGURATIONS QUE L'ON RETIENT
! ON EXCLUT CELLES QUE L'ON NE VEUT PAS POUR AVOIR UNE MEILLEURE
! VISIBILITE DES DIFFERENTES CONFIGURATIONS QUI SE PRESENTENT
!
!       NBRE DE POINT D'INTERSECTION INCORRECT (1) OU (2)
        if (ninter .le. 1) then
            call utmess('F', 'XFEM_64')
!
!       NBRE DE POINT D'INTERSECTION INCORRECT (1) OU (2)
        else if (ninter .gt. 3) then
            call utmess('F', 'XFEM_64')
!
!       NBRE PT INTER SOMMET > NBRE PT INTER TOTAL (1)
        else if (npts .gt. ninter) then
            call utmess('F', 'XFEM_64')
!
!       LA FISSURE INTERCEPTE DEUX NOEUDS SOMMETS UNIQUEMENT (1) OU (2)
        else if (ninter .eq. 2 .and. npts .eq. 2) then
            call utmess('F', 'XFEM_64')
!
!       LA FISSURE INTERCEPTE LES 3 ARETES STRICTEMENT (2)
        else if (ninter .eq. 3 .and. npts .eq. 0) then
            call utmess('F', 'XFEM_64')
!
!       (2)
        else if (ninter .eq. 3 .and. npts .eq. 1 .and. nptm .ne. 1) then
            call utmess('F', 'XFEM_64')
!
!       LA FISSURE JOUXTE UN BORD DE L'ELEMENT (1)
        else if (ninter .eq. 3 .and. npts .eq. 2 .and. nptm .eq. 1) then
            call utmess('F', 'XFEM_64')
!
!       (1) OU (2)
        else if (ninter .eq. 3 .and. npts .eq. 2 .and. nptm .ne. 1) then
            call utmess('F', 'XFEM_64')
!
!       (1) OU (2)
        else if (ninter .eq. 3 .and. npts .eq. 3) then
            call utmess('F', 'XFEM_64')
!
        end if
!
! --- POUR LES SEG3
!
    else if (ndime .eq. 1) then
!
!       NBRE DE POINT D'INTERSECTION INCORRECT (1) OU (2)
        if (ninter .ne. 1 .and. npts .ne. 0) then
            call utmess('F', 'XFEM_64')
        end if
!
    end if
!
end subroutine
