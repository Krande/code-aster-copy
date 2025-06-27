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
function iposdg(dg, cmp)
! aslint: disable=
    implicit none
    integer(kind=8) :: iposdg
!
    integer(kind=8) :: dg(*), cmp
!
!
!***********************************************************************
!
!     REND LA POSITION D'1 CMP DANS UN DESCRIPTEUR-GRANDEUR DG
!
!       DG   (IN) : TABLE DES ENTIERS CODES
!
!       CMP  (IN) : NUMERO DE LA COMPOSANTE
!
!     EXEMPLE
!
!      POUR LE GRANDEUR DEPLA_R :
!
!         L' ENVELOPPE DU DESCRIPTEUR EST
!
!            DX DY DZ DRX DRY DRZ LAGR
!
!         SUPPOSONS QUE LA DESCRIPTION LOCALE SOIT : DX DZ DRY
!
!         ALORS IPOSDG(DG,NUM('DZ') ---> 2
!
!         LA COMPOSAMTE DZ APPARAIT EN POSITION 2 DANS LA DESCRIPTION
!
!***********************************************************************
!
    integer(kind=8) :: paquet, valec, nbec, reste, code, cmpt, i, lshift
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    nbec = (cmp-1)/30+1
    cmpt = 0
    reste = cmp-30*(nbec-1)
    code = lshift(1, reste)
!
    if (iand(dg(nbec), code) .eq. code) then
!
        do paquet = 1, nbec-1, 1
!
            valec = dg(paquet)
!
            do i = 1, 30, 1
!
                code = lshift(1, i)
!
                if (iand(valec, code) .eq. code) then
!
                    cmpt = cmpt+1
!
                end if
!
            end do
!
        end do
!
        valec = dg(nbec)
!
        do i = 1, reste, 1
!
            code = lshift(1, i)
!
            if (iand(valec, code) .eq. code) then
!
                cmpt = cmpt+1
!
            end if
!
        end do
!
    end if
!
    iposdg = cmpt
!
end function
