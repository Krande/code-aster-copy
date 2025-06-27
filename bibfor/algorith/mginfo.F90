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

subroutine mginfo(modeMecaZ, numeDof_, nbmode_, nbEqua_, occ_)
!
    implicit none
!
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jelira.h"
!
    character(len=*), intent(in) :: modeMecaZ
    integer(kind=8), optional, intent(out) :: nbmode_, nbEqua_
    character(len=14), optional, intent(out) :: numeDof_
    integer(kind=8), optional, intent(in) :: occ_
!
!
! ----------------------------------------------------------------------
!
! UTILITAIRE
!
! INFORMATIONS SUR MATRICE MODES MECANIQUES
!
! ----------------------------------------------------------------------
!
! IN  MODMEC : NOM DE LA MATRICE DES MODES MECANIQUES
! OUT NUMDDL : NOM DU DDL
! OUT NBMODE : NOMBRE DE MODES
! OUT NEQ    : NOMBRE D'EQUATIONS

    character(len=24) :: matrix
    character(len=4) :: indik4
    character(len=8) :: indik8, modeMeca
    integer(kind=8) :: nbmode, nbEqua, occ, ier, ier2
    character(len=14) :: numeDof
!
! ----------------------------------------------------------------------
!
    occ = 1
    if (present(occ_)) then
        occ = occ_
    end if
    call codent(occ, 'D0', indik8)
    call codent(occ, 'D0', indik4)

    modeMeca = modeMecaZ
    nbEqua = 0
    nbMode = 0
    numeDof = ' '
    call dismoi('NUME_CHAM_'//indik8, modeMeca, 'RESU_DYNA', repk=numeDof, &
                arret='C', ier=ier)
    if (ier /= 0) then
        call dismoi('NUME_DDL', modeMeca, 'RESU_DYNA', repk=numeDof)
        call dismoi('NB_EQUA', numeDof, 'NUME_DDL', repi=nbEqua)
    else if (.not. present(occ_)) then
        call dismoi('NB_EQUA', numeDof, 'NUME_DDL', repi=nbEqua)
    else
        call dismoi('REF_RIGI_'//indik4, modeMeca, 'RESU_DYNA', repk=matrix, &
                    arret='C', ier=ier2)
        if (ier2 /= 0) call dismoi('REF_RIGI_PREM', modeMeca, &
                                   'RESU_DYNA', repk=matrix)
        call dismoi('NOM_NUME_DDL', matrix, 'MATR_ASSE', repk=numeDof)
        call dismoi('NB_EQUA', matrix, 'MATR_ASSE', repi=nbEqua)
    end if
    call jelira(modeMeca(1:8)//'           .ORDR', 'LONMAX', nbmode)

    if (present(nbEqua_)) then
        nbEqua_ = nbEqua
    end if
    if (present(nbMode_)) then
        nbMode_ = nbMode
    end if
    if (present(numeDof_)) then
        numeDof_ = numeDof
    end if
!
end subroutine
