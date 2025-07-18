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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmch4p(veelem)
!
    implicit none
!
#include "asterfort/nmcha0.h"
!
    character(len=19) :: veelem(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (INITIALISATION)
!
! CREATION DES VARIABLES CHAPEAUX - VEELEM
!
! ----------------------------------------------------------------------
!
! OUT VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
!
! ----------------------------------------------------------------------
!
    character(len=19) :: vefsdo, vebudi, vedido, vesstf
    character(len=19) :: vedipi, vefedo, vefepi, veondp
    character(len=19) :: vedidi, vediri
    character(len=19) :: verefe
    character(len=19) :: veimpe
!
    data vefedo, vefsdo/'&&NMCH4P.VEFEDO', '&&NMCH4P.VEFSDO'/
    data vedido, vefepi/'&&NMCH4P.VEDIDO', '&&NMCH4P.VEFEPI'/
    data vedipi/'&&NMCH4P.VEDIPI'/
    data vebudi, vedidi/'&&NMCH4P.VEBUDI', '&&NMCH4P.VEDIDI'/
    data veondp/'&&NMCH4P.VEONDP'/
    data vediri/'&&NMCH4P.VEDIRI'/
    data vesstf/'&&NMCH4P.VESSTF'/
    data verefe/'&&NMCH4P.VEREFE'/
    data veimpe/'&&NMCH4P.VEIMPE'/
!
! ----------------------------------------------------------------------
!
    call nmcha0('VEELEM', 'ALLINI', ' ', veelem)
    call nmcha0('VEELEM', 'CNDIRI', vediri, veelem)
    call nmcha0('VEELEM', 'CNBUDI', vebudi, veelem)
    call nmcha0('VEELEM', 'CNDIDO', vedido, veelem)
    call nmcha0('VEELEM', 'CNDIPI', vedipi, veelem)
    call nmcha0('VEELEM', 'CNFEDO', vefedo, veelem)
    call nmcha0('VEELEM', 'CNFEPI', vefepi, veelem)
    call nmcha0('VEELEM', 'CNONDP', veondp, veelem)
    call nmcha0('VEELEM', 'CNFSDO', vefsdo, veelem)
    call nmcha0('VEELEM', 'CNIMPE', veimpe, veelem)
    call nmcha0('VEELEM', 'CNDIDI', vedidi, veelem)
    call nmcha0('VEELEM', 'CNSSTF', vesstf, veelem)
    call nmcha0('VEELEM', 'CNREFE', verefe, veelem)
!
end subroutine
