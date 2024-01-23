! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine te0067(option, nomte)
!
    use Metallurgy_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/nzcomp_prep.h"
#include "asterfort/nzcomp.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: THERMIQUE - AXIS*, PLAN*
! Option: META_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: metaType
    real(kind=8) :: dt10, dt21, inst2
    real(kind=8) :: tno1, tno0, tno2
    integer :: nbNode, iNode, itempe, itempa, jvTime
    integer :: jvMater, nbVari, numeComp, nbPhase
    integer :: itempi
    integer :: jvPhaseIn, jvPhaseOut, jvComporMeta
    integer :: jvMaterCode
    type(META_MaterialParameters) :: metaPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nbNode)

! - Input/Output fields
    call jevech('PMATERC', 'L', jvMater)
    jvMaterCode = zi(jvMater)
    call jevech('PTEMPAR', 'L', itempa)
    call jevech('PTEMPER', 'L', itempe)
    call jevech('PTEMPIR', 'L', itempi)
    call jevech('PTIMMTR', 'L', jvTime)
    call jevech('PPHASIN', 'L', jvPhaseIn)
    call jevech('PCOMPOR', 'L', jvComporMeta)
    call jevech('PPHASNOU', 'E', jvPhaseOut)

! - Parameters from map
    metaType = zk16(jvComporMeta-1+ZMETATYPE)
    read (zk16(jvComporMeta-1+ZNUMECOMP), '(I16)') numeComp
    read (zk16(jvComporMeta-1+ZNBPHASE), '(I16)') nbPhase
    read (zk16(jvComporMeta-1+ZNBVARI), '(I16)') nbVari

! - Preparation
    call nzcomp_prep(jvMaterCode, metaType, metaPara)

! - Time parameters: 0 - 1 - 2
    dt10 = zr(jvTime+1)
    dt21 = zr(jvTime+2)
    inst2 = zr(jvTime)+dt21

! - Loop on nodes
    do iNode = 1, nbNode
! ----- Temperatures: 0 - 1 - 2
        tno1 = zr(itempe+iNode-1)
        tno0 = zr(itempa+iNode-1)
        tno2 = zr(itempi+iNode-1)

! ----- General switch
        call nzcomp(jvMaterCode, metaPara, &
                    numeComp, nbPhase, nbVari, &
                    dt10, dt21, inst2, &
                    tno0, tno1, tno2, &
                    zr(jvPhaseIn+nbVari*(iNode-1)), &
                    zr(jvPhaseOut+nbVari*(iNode-1)))
    end do
!
end subroutine
