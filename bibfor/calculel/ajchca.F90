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
subroutine ajchca(para, cham, lpara, lcham, nbent, &
                  maxent, surch)
    implicit none
!
#include "asterfort/assert.h"
    integer(kind=8) :: maxent, nbent, indice, i
    character(len=*) :: para, cham, lpara(maxent), lcham(maxent), surch
!
!     CETTE ROUTINE PERMET D'AJOUTER UN COUPLE (PARAMETRE,NOM_CHAMP)
!     A UNE LISTE LPAIN/LCHIN (OU LPAOUT/LCHOUT)
!
! ----------------------------------------------------------------------
!     IN  : PARA    : NOM DE PARAMETRE CORRESPONDANT A LPARA
!     IN  : CHAM    : NOM DE CHAMP     CORRESPONDANT A LCHAM
!     I/O : LPARA  : TABLEAU DES PARAMETRES
!     I/O : LCHAM  : TABLEAU DES NOMS DES CHAMPS
!     I/O : NBENT  : NOMBRE D'ENTREES
!     IN  : MAXENT : NOMBRE D'ENTREES MAXIMUM
!     IN  : SURCH :  'O' OU 'N'
!
!     SURCH SERT A DETERMINER CE QUE L'ON FAIT SI LE PARAMETRE PARA
!     APPARTIENT DEJA A LPARA.
!
!     ANCIEN       AJOUT    SURCH      NOUVEAU
!     ' '          CH2      O/N        CH2
!     CH1          ' '      O/N        CH1
!     CH1          CH1      O/N        CH1
!     CH1          CH2      N          CH1
!     CH1          CH2      O          CH2
! ----------------------------------------------------------------------
!
!
!     1. RECHERCHE SI LE PARAMETRE EXISTE DEJA
    indice = 0
    do i = 1, nbent
        if (lpara(i) .eq. para) then
            indice = i
            goto 20
        end if
    end do
20  continue
!
!
!     2. IL S'AGIT D'UN NOUVEAU PARAMETRE ON L'AJOUTE :
!     -------------------------------------------------
    if (indice .eq. 0) then
        ASSERT(nbent .lt. maxent)
        nbent = nbent+1
        lpara(nbent) = para
        lcham(nbent) = cham
        goto 999
    end if
!
!
!     3. LE NOM DU PARAMETRE EST DEJA DANS LPARA :
!     -------------------------------------------------
!     ANCIEN       AJOUT    SURCH      NOUVEAU
!     ' '          CH2      O/N        CH2
!     CH1          ' '      O/N        CH1
!     CH1          CH1      O/N        CH1
!     CH1          CH2      N          CH1
!     CH1          CH2      O          CH2
!
!     3.1 L'ANCIEN CHAMP ETAIT "BLANC" : ON STOCKE LE NOUVEAU:
    if (lcham(indice) .eq. ' ') then
        lcham(indice) = cham
        goto 999
    end if
!
!
!     3.2 LE NOUVEAU CHAMP EST "BLANC" : ON NE LE STOCKE PAS
    if (cham .eq. ' ') then
        goto 999
    end if
!
!
!     3.3 LE NOUVEAU NOM EST LE MEME QUE L'ANCIEN :
    if (lcham(indice) .eq. cham) then
        goto 999
    end if
!
!
!     3.4 L'ANCIEN NOM DU CHAMP N'ETAIT PAS "BLANC"
!         LE NOUVEAU NOM DU CHAMP NON PLUS ET IL EST DIFFERENT :
    if (surch .eq. 'O') then
        lcham(indice) = cham
    else if (surch .eq. 'N') then
    else
        ASSERT(.false.)
    end if
!
!
999 continue
end subroutine
