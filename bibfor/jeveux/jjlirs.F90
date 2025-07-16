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

subroutine jjlirs(jadm, iclas, idos, ius, ist)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "jeveux_private.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: jadm, iclas, ius, ist
! ----------------------------------------------------------------------
!     RELIT LES ENTIERS ENCADRANT UN SEGMENT DE VALEURS
!
! IN  JADM   : ADRESSE DU PREMIER MOT DU SEGMENT DE VALEUR
! IN  ICLAS  : CLASSE DE L'OBJET JEVEUX
! IN  IDOS   : IDENTIFICATEUR D'OBJET SIMPLE OU D'OBJET DE COLLECTION
! OUT IUS    : USAGE DU SEGMENT DE VALEUR
! OUT IST    : STATUT DU SEGMENT DE VALEUR
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
    integer(kind=8) :: istat
    common/istaje/istat(4)
! DEB ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: icla2, idatoc, idos, is, ista1, ista2
!
!-----------------------------------------------------------------------
    ista1 = iszon(jiszon+jadm-1)
    idatoc = iszon(jiszon+jadm-2)
    if (idatoc .ne. idos) then
        call utmess('F', 'JEVEUX1_54', si=jadm)
    end if
!
    if (ista1 .ne. istat(1) .and. ista1 .ne. istat(2)) then
        call utmess('F', 'JEVEUX1_54', si=jadm)
    end if
!
    is = jiszon+iszon(jiszon+jadm-4)
    ista2 = iszon(is-4)
    icla2 = iszon(is-2)
    if (icla2 .ne. iclas) then
        call utmess('F', 'JEVEUX1_55', si=jadm)
    end if
!
    if (ista2 .ne. istat(3) .and. ista2 .ne. istat(4)) then
        call utmess('F', 'JEVEUX1_55', si=jadm)
    end if
!
    ius = ista1
    ist = ista2
! FIN ------------------------------------------------------------------
end subroutine
