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

subroutine gnoms3(noojb, k1, k2, test)
    implicit none
!
!
#include "asterfort/codlet.h"
#include "asterfort/jeexin.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: noojb, test
    integer(kind=8), intent(in) :: k1, k2
! BUT :
!  TROUVER UN NOM POSSIBLE POUR UN OBJET JEVEUX QUI RESPECTE :
!     - CE NOM VAUT NOOJB (DONNE EN ENTREE) SAUF POUR LA SOUS-CHAINE
!           NOOJB(K1:K2)
!     - LE NOM DE L'OBJET N'EXISTE PAS ENCORE DANS LES BASES OUVERTES
!     - LE NOM (K1:K2) EST UN NUMERO ('0001','0002', ...)
!
! VAR : NOOJB : NOM D'UN OBJET JEVEUX  (K*)
! IN  : K1,K2 : INDICES DANS NOOJB DE LA SOUS-CHAINE "NUMERO"
!     -----------------------------------------------------------------
    integer(kind=8) :: entier, nb_essai, i, iret
!----------------------------------------------------------------------

    entier = k2-k1+1
    nb_essai = 36**entier
!   name of partition
    do i = 0, nb_essai
        call codlet(i, "D0", noojb(k1:k2), "F")
        call jeexin(noojb//test, iret)
        if (iret .eq. 0) goto 20
    end do
    call utmess('F', 'MODELISA4_69', si=nb_essai)
20  continue

end subroutine
