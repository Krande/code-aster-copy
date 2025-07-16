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
subroutine te0064(option, nomte)
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
#include "asterfort/nzcompTemper.h"
#include "asterfort/tecach.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: THERMIQUE - 3D
! Option: META_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: deltaTime01, deltaTime12, time2
    real(kind=8) :: temp1, tempInit, temp2
    integer(kind=8) :: nbNode, iNode, itempe, jvTempInit, jvTime
    integer(kind=8) :: jvMater
    character(len=16) :: metaType
    integer(kind=8) :: numeComp
    integer(kind=8) :: nbVari, nbPhase
    integer(kind=8) :: nbVariTemper, nbPhaseTemper, nbVariPrev
    integer(kind=8) :: itempi
    integer(kind=8) :: jvPhaseIn, jvPhaseOut, jvPhasePrev
    integer(kind=8) :: jvMaterCode
    integer(kind=8) :: itab(7), iret, jvComporMetaTemper
    character(len=16), pointer :: comporMeta(:) => null()
    aster_logical :: hasTemper, prevMetaIsTemper
    type(META_MaterialParameters) :: metaPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nbNode)
    ASSERT(nbNode .le. MT_NNOMAX3D)

! - Input/Output fields
    call jevech('PMATERC', 'L', jvMater)
    jvMaterCode = zi(jvMater)
    call jevech('PTEMPER', 'L', itempe)
    call jevech('PTEMPIR', 'L', itempi)
    call jevech('PTIMMTR', 'L', jvTime)
    call jevech('PPHASIN', 'L', jvPhaseIn)
    call jevech('PCOMPME', 'L', vk16=comporMeta)
    call jevech('PPHASOUT', 'E', jvPhaseOut)

! - Get parameters from metallurgy behaviour (without tempering)
    metaType = comporMeta(ZMETATYPE)
    read (comporMeta(ZNUMECOMP), '(I16)') numeComp
    read (comporMeta(ZNBPHASE), '(I16)') nbPhase
    read (comporMeta(ZNBVARI), '(I16)') nbVari

! - Specific input/output fields
    hasTemper = ASTER_FALSE
    call tecach("ONN", "PCOMPMT", 'L', iret, nval=7, itab=itab)
    if (iret .eq. 0) then
        jvComporMetaTemper = itab(1)
        hasTemper = ASTER_TRUE
        if (metaType .ne. 'ACIER') then
            ASSERT(ASTER_FALSE)
        end if
        call jevech('PPHASEP', 'L', jvPhasePrev)
        call tecach("ONN", "PPHASEP", 'L', iret, nval=7, itab=itab)
        nbVariPrev = itab(6)
        if (nbVariPrev .eq. NBVARIJMA) then
            prevMetaIsTemper = ASTER_TRUE
        elseif (nbVariPrev .eq. NBVARIWAECKEL) then
            prevMetaIsTemper = ASTER_FALSE
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        call jevech('PTEMPAR', 'L', jvTempInit)
    end if

! - Get parameters from metallurgy behaviour (with tempering)
    if (hasTemper) then
        metaType = zk16(jvComporMetaTemper-1+ZMETATYPE)
        read (zk16(jvComporMetaTemper-1+ZNUMECOMP), '(I16)') numeComp
        read (zk16(jvComporMetaTemper-1+ZNBPHASE), '(I16)') nbPhaseTemper
        read (zk16(jvComporMetaTemper-1+ZNBVARI), '(I16)') nbVariTemper
        if (metaType .ne. 'ACIER_REVENU') then
            ASSERT(ASTER_FALSE)
        end if
    end if

! - Get material and TRC parameters
    call nzcomp_prep(jvMaterCode, metaType, metaPara)

! - Time parameters: 1 - 2
    deltaTime12 = zr(jvTime+2)
    time2 = zr(jvTime)+deltaTime12

! - Loop on nodes
    do iNode = 1, nbNode
! ----- Temperatures: 1 - 2
        temp1 = zr(itempe+iNode-1)
        temp2 = zr(itempi+iNode-1)

! ----- General switches
        if (hasTemper) then
            call nzcompTemper(metaPara, numeComp, &
                              nbVari, nbVariTemper, nbVariPrev, &
                              deltaTime12, &
                              temp1, temp2, &
                              prevMetaIsTemper, &
                              zr(jvPhasePrev+nbVariPrev*(iNode-1)), &
                              zr(jvPhaseIn+nbVari*(iNode-1)), &
                              zr(jvPhaseOut+nbVariTemper*(iNode-1)))

        else
            deltaTime01 = zr(jvTime+1)
            tempInit = zr(jvTempInit+iNode-1)
            call nzcomp(jvMaterCode, metaPara, numeComp, &
                        nbPhase, nbVari, &
                        deltaTime01, deltaTime12, time2, &
                        tempInit, temp1, temp2, &
                        zr(jvPhaseIn+nbVari*(iNode-1)), &
                        zr(jvPhaseOut+nbVari*(iNode-1)))
        end if
    end do
!
end subroutine
