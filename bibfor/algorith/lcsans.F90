! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine lcsans(option, sigp, dsidep)
!---&s---1---------2---------3---------4---------5---------6---------7--
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG, FULL_MECA ,RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE
!
!_______________________________________________________________________
!
!
! ROUTINE CALCULANT UN COMPORTEMENT VIDE:
!
!    - CONTRAINTES FINALES NULLES         : SIGP (NSTRS)
!
!_______________________________________________________________________
!
implicit none
!
character(len=16) :: option
real(kind=8) :: sigp(6), dsidep(6, 6)
!
    if ((option(1:9).eq.'FULL_MECA') .or. (option(1:9).eq.'RAPH_MECA')) then
        sigp = 0.d0
    endif
!
!_______________________________________________________________________
!
! CONSTRUCTION DE LA MATRICE TANGENTE
!_______________________________________________________________________
!
    if ((option(1:9).eq.'FULL_MECA') .or. (option(1:9).eq.'RIGI_MECA')) then
        dsidep = 0.d0
    endif
!
end subroutine
