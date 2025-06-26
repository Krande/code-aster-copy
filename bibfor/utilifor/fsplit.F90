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

subroutine fsplit(fieldType, fieldSupport, fieldQuantity)
    implicit none
#include "asterfort/assert.h"
    character(len=16) :: fieldType
    character(len=4) :: fieldSupport
    character(len=8) :: fieldQuantity
!-----------------------------------------------------------------------
!     ROUTINE FORTRAN PERMETTANT DE DECOUPER LE TYPE DE CHAMP
!     POUR RECUPERER LE SUPPORT DU CHAMP (EL**, NOEU, CART)
!     ET LE NOM DE LA GRANDEUR PHYSIQUE AVEC LE TYPE FORTRAN (R, C, I, F)
!     (ATTENTION A LA CONFUSION, EN ASTER CE QU'ON APPELLE TYPE DE CHAMP
!     EST EN REALITE LE NOM DU CHAMP ET PAS LE TYPE FORTRAN)
!
!     EN PRINCIPE LE TYPE DE CHAMP EST FORMATE EN XXXX_YYYY_Z
!     AVEC XXXX LE SUPPORT DU CHAMP SUR 4 CARACTERES
!          YYYY LE NOM DE LA GRANDEUR PHYSIQUE SUR 4 CARACTERES
!       ET Z    LE TYPE FORTRAN
!     MAIS IL EXISTE DES EXCEPTIONS DU GENRE XXXX_YYYYY..._Z OU LE NOM DE
!     LA GRANDEUR PHYSIQUE A PLUS DE 4 CARACTERES (MAX 8) ET DU GENRE
!     XXXX_YYYYY... OU DE PLUS LE TYPE FORTRAN EST ABSENT
!
!     SI LE TYPE FORTRAN EST ABSENT, ON LE MET A R PAR DEFAUT
!-----------------------------------------------------------------------
! IN   : fieldType   : TYPE DE CHAMP en principe au format XXXX_YYYY_Z
! OUT  : fieldSupport   : SUPPORT DU CHAMP XXXX
! OUT  : fieldQuantity  : NOM DE LA GRANDEUR PHYSIQUE AVEC SON TYPE YYYY..._Z
! ----------------------------------------------------------------------
    integer(kind=8) :: I, J, K, L
!-----------------------------------------------------------------------
!
    fieldSupport = ''
    fieldQuantity = ''
!
    L = len_trim(fieldType)
!
    I = index(fieldType, '_')
    if (I .eq. 0) then
! LE TYPE DE CHAMP POSSEDE AU MOINS UN '_'
        ASSERT(.FALSE.)
    end if
! LE SUPPORT DOIT TENIR SUR 4 CARACTERES
    ASSERT(I .eq. 5)
    fieldSupport = fieldType(1:I-1)
!
    J = index(fieldType(I+1:L), '_')
    if (J .eq. 0) then
! ON EST DANS LE CAS XXXX_YYYYY...
        fieldQuantity = fieldType(I+1:L)//'_R'
    else
! ON EST DANS LE CAS XXXX_YYYY_Z OU XXXX_YYYYY..._Z
        fieldQuantity = fieldType(I+1:L)
! ON VERIFIE QUAND MEME QUE L'ON N'A PAS D'AUTRE '_'
        K = index(fieldType(I+J+1:L), '_')
        ASSERT(K .eq. 0)
    end if
!
end subroutine
