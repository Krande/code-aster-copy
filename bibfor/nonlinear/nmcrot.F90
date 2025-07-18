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
subroutine nmcrot(result, sd_obsv)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/nonlinDSTableIOSetPara.h"
#include "asterfort/nonlinDSTableIOGetName.h"
#include "asterfort/nonlinDSTableIOCreate.h"
#include "asterfort/nonlinDSTableIOClean.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: result
    character(len=19), intent(in) :: sd_obsv
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear operators - Observation
!
! Create table
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of results datastructure
! In  sd_obsv          : datastructure for observation parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPara = 17
    type(NL_DS_TableIO) :: tableio
    character(len=24), parameter :: paraName(nbPara) = (/'NOM_OBSERVATION  ', 'TYPE_OBJET       ', &
                                                         'NUME_REUSE       ', 'NUME_OBSE        ', &
                                                         'INST             ', 'NOM_CHAM         ', &
                                                         'EVAL_CHAM        ', 'NOM_CMP          ', &
                                                         'NOM_VARI         ', 'EVAL_CMP         ', &
                                                         'NOEUD            ', 'MAILLE           ', &
                                                         'EVAL_ELGA        ', 'POINT            ', &
                                                         'SOUS_POINT       ', 'VALE             ', &
                                                         'NUME_GLOBAL_NOEUD'/)
    character(len=8), parameter :: paraType(nbPara) = (/'K16', 'K16', &
                                                        'I  ', 'I  ', &
                                                        'R  ', 'K16', &
                                                        'K8 ', 'K8 ', &
                                                        'K16', 'K8 ', &
                                                        'K8 ', 'K8 ', &
                                                        'K8 ', 'I  ', &
                                                        'I  ', 'R  ', &
                                                        'I  '/)
    character(len=24) :: obsv_tabl
    character(len=24), pointer :: v_obsv_tabl(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    tableio%tablSymbName = 'OBSERVATION'
    tableio%resultName = result
!
! - Create vector for table name
!
    obsv_tabl = sd_obsv(1:14)//'     .TABL'
    call wkvect(obsv_tabl, 'V V K24', 1, vk24=v_obsv_tabl)
!
! - Create list of parameters
!
    call nonlinDSTableIOSetPara(tableio_=tableio, &
                                nbPara_=nbPara, &
                                paraName_=paraName, &
                                paraType_=paraType)
!
! - Get name of table in results datastructure
!
    call nonlinDSTableIOGetName(tableio)
!
! - Create table in results datastructure
!
    call nonlinDSTableIOCreate(tableio)
!
! - Save name of table
!
    v_obsv_tabl(1) = tableio%tablName
!
! - Clean table in results datastructure
!
    call nonlinDSTableIOClean(tableio)
!
end subroutine
