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

subroutine cffact(ldscon, isto, nbliac, &
                  indfac, lechec)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/tldlg3.h"
!
!
    integer(kind=8) :: nbliac, indfac
    integer(kind=8) :: ldscon
    integer(kind=8) :: isto
    aster_logical :: lechec
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - RESOLUTION)
!
! FACTORISATION LDLT DE [-A.C-1.AT]
!
! ----------------------------------------------------------------------
!
! ATTENTION : SI ON RAJOUTE DES LIAISONS ON NE FACTORISE QUE
! LA PARTIE RAJOUTEE (LE RESTE EST ENCORE VALABLE, CF. PROPRIETES
! MAGIQUES DES FACTORISATIONS).
! SI ON ENLEVE LA DERNIERE LIAISON (IDEBUT > NBLIAC),PAS BESOIN DE
! REFACTORISER : L'INSTRUCTION ZI(LDSCON+2) = NBLIAC ECRITE PLUS
! LOIN FERA QUE RLDLGG PRENDRA LA BONNE TAILLE DE MATRICE, QUI
! EST DEJA FACTORISEE (SI ON REFACTORISAIT A PARTIR DE 1, ON
! FACTORISERAIT LA FACTORISEE, CE QUI EST GENANT, CAR
! FACTORISATION EN PLACE)
!
! IN  LDSCON : DESCRIPTEUR DE LA MATRICE DE CONTACT
! IN  ISTO   : INDICATEUR D'ARRET EN CAS DE PIVOT NUL
! IN  NBLIAC : NOMBRE DE LIAISONS ACTIVES
! I/O INDFAC : INDICE DE DEBUT DE LA FACTORISATION
! OUT LECHEC : .TRUE. SI LA FACTORISATION A ECHOUE (MATRICE SINGULIERE)
!
!
!
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: ilideb, ilifin, ier
    integer(kind=8) :: ndeci, isingu, npvneg
!
! ----------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    lechec = .false.
    if (indfac .le. nbliac) then
        if (niv .ge. 2) then
            write (ifm, *) '<CONTACT><CALC> FACTORISATION MATRICE CONTACT '
        end if
        ilideb = indfac
        ilifin = nbliac
        call tldlg3('LDLT', ' ', 2, ldscon, ilideb, ilifin, 0, &
                    ndeci, isingu, npvneg, ier, ' ')
        indfac = ilifin+1
        if (ier .gt. isto) then
            lechec = .true.
        end if
    end if
!
end subroutine
