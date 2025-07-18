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

subroutine elmddl(raide, option, neq, ddl, nddle, &
                  nbddl, vecddl)
!
!
! aslint: disable=W1306
    implicit none
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/pteddl.h"
    character(len=19) :: raide
    character(len=14) :: option
    integer(kind=8) :: neq, nbddl, vecddl(neq), nddle
    character(len=8) :: ddl(nddle)
!
! ----------------------------------------------------------------------
!
! CONSTRUCTION D'UN TABLEAU D'ENTIERS REPERANT LA POSITION DES DDL
! EXCLUS DE LA RECHERCHE DE VALEURS PROPRES
!
! ----------------------------------------------------------------------
!
! IN  RAIDEUR : NOM DE LA MATRICE DE "RAIDEUR"
! IN  OPTION  : TYPE DE DDL A TROUVER
! IN  NEQ     : NPMBRE DE DDL
! IN  DDL     : NOM DU DDL A ELIMINER
! IN  NDDLE   : NOMBRE DE TYPES DE DDL EXCLUS
! OUT NBDDL   : NOMBRE DE DDL A ELIMINER
! OUT VECDDL  : POSITION DES DDL A ELIMINER
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ieq, ifm, niv, i, inter(neq)
    character(len=14) :: nume
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
! --- INITIALISATIONS
!
    nbddl = 0
    call dismoi('NOM_NUME_DDL', raide, 'MATR_ASSE', repk=nume)
    do ieq = 1, neq
        vecddl(ieq) = 1
    end do
!
! --- CALCUL DU NOMBRE DE DDL A ELIMINER
!
    if (nddle .gt. 0) then
        do i = 1, nddle
!
! ------- RECUPERATION DES POSITIONS DES DDL
!
            call pteddl('NUME_DDL', nume, 1, ddl(i), neq, &
                        list_equa=inter)
!
! ------- CALCUL DU NOMBRE DE 'DDL': NBDDL
!
            do ieq = 1, neq
                nbddl = nbddl+inter(ieq)
            end do
!
! ------- STOP SI ON CHERCHE A ELIM UN DDL ABSENT DE LA MODELISATION
!
            if (nbddl .eq. 0) then
                ASSERT(.false.)
            end if
!
! ------- INVERSION : INTER = 0 SI DDL TROUVE ET 1 SINON
!
            do ieq = 1, neq
                inter(ieq) = abs(inter(ieq)-1)
            end do
!
            do ieq = 1, neq
                vecddl(ieq) = vecddl(ieq)*inter(ieq)
            end do
!
        end do
    end if
!
! --- IMPRESSION DES DDL
!
    if (niv .ge. 1) then
        if (nbddl .gt. 0) then
            write (ifm, *) option
            do i = 1, nddle
                write (ifm, 910) ddl(i)
            end do
            write (ifm, 950) nbddl
            write (ifm, 960)
        else
            write (ifm, 901)
        end if
    end if
!
! ----------------------------------------------------------------------
!
    call jedema()
!
901 format(1x, 'PAS DE DDL_TROUVE')
910 format(13x, a8,/)
950 format(1x, 'NOMBRE DE DDL', 10x, i7,/)
960 format(72('-'))
!
end subroutine
