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
subroutine dyGetKineLoad(matrRigiz, matrMassz, matrDampz, lDamp, listLoadz, kineLoad, &
                         integScheme_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/idenob.h"
#include "asterfort/ischar_iden.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisico.h"
#include "asterfort/lislch.h"
#include "asterfort/lislco.h"
#include "asterfort/lisnch.h"
#include "asterfort/lisnnb.h"
#include "asterfort/utmess.h"
#include "asterfort/isParallelMatrix.h"
#include "asterfort/asmpi_any.h"
!
    character(len=*), intent(in) :: matrRigiz, matrMassz, matrDampz, listLoadz
    aster_logical, intent(in) :: lDamp
    character(len=24), intent(out) :: kineLoad
    integer(kind=8), optional, intent(in) :: integScheme_
!
! --------------------------------------------------------------------------------------------------
!
! Get kinematic load
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: matrRigi, matrMass, matrDamp, listload
    character(len=24) :: listLoadInfoName, listLoadName
    character(len=8) :: answer, kineLoadName
    character(len=8) :: matrRigiMesh, matrMassMesh, matrDampMesh, kineLoadMesh
    aster_logical :: lKineLoadInRigi, lKineLoadInMass, lKineLoadInDamp, lKineLoad, lTransient
    aster_logical :: l_parallel_matrix, lKineLoadGlobal
    integer(kind=8) :: matrRigiNbEqua, matrMassNbEqua, matrDampNbEqua
    integer(kind=8) :: iLoad, nbLoad, iLoadKine, genrec
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    character(len=24), pointer :: vlistLoadName(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    matrRigi = matrRigiz
    matrMass = matrMassz
    matrDamp = matrDampz
    listLoad = listLoadz
    kineLoad = ' '
    lTransient = ASTER_FALSE
    if (present(integScheme_)) then
        lTransient = ASTER_TRUE
    end if

! - Get CHAM_CINE from list of loads
    if (lTransient) then
        call lisnch(listLoad, nbLoad)
        listLoadName = listLoad(1:19)//'.LCHA'
        listLoadInfoName = listLoad(1:19)//'.INFC'
        if (nbLoad .ne. 0) then
            call jeveuo(listLoadInfoName, 'L', vi=listLoadInfo)
            call jeveuo(listLoadName, 'L', vk24=vlistLoadName)
        end if
        iLoadKine = 0
        kineLoadName = ' '
        do iLoad = 1, nbLoad
            if (ischar_iden(listLoadInfo, iLoad, nbLoad, 'DIRI', 'ELIM', vlistLoadName(iLoad))) then
                if (iLoadKine .ne. 0) then
                    call utmess('F', 'DYNALINE2_13')
                end if
                iLoadKine = iLoad
            end if
        end do
        lKineLoad = iLoadKine .ne. 0
        if (lKineLoad) then
            kineLoadName = vlistLoadName(iLoadKine) (1:8)
        end if
    else
        call lisnnb(listLoad, nbLoad)
        iLoadKine = 0
        kineLoadName = ' '
        do iLoad = 1, nbLoad
            call lislco(listLoad, iLoad, genrec)
            if (lisico('DIRI_ELIM', genrec)) then
                if (iLoadKine .ne. 0) then
                    call utmess('F', 'DYNALINE2_13')
                end if
                iLoadKine = iLoad
            end if
        end do
        lKineLoad = iLoadKine .ne. 0
        if (lKineLoad) then
            call lislch(listLoad, iLoadKine, kineLoadName)
        end if
    end if
! - One subdomain (at least) has a kinematic load
    lKineLoadGlobal = asmpi_any(lKineLoad, ASTER_TRUE)

! - Some parameters from matrices
    l_parallel_matrix = isParallelMatrix(matrRigi)
    ASSERT(l_parallel_matrix .eqv. isParallelMatrix(matrMass))
    if (lDamp) then
        ASSERT(l_parallel_matrix .eqv. isParallelMatrix(matrDamp))
    end if
!
    lKineLoadInDamp = ASTER_FALSE
    matrDampMesh = ' '
    matrDampNbEqua = 0
    call dismoi('EXIS_CINE', matrRigi, 'MATR_ASSE', repk=answer)
    lKineLoadInRigi = answer .eq. 'OUI'
    call dismoi('NOM_MAILLA', matrRigi, 'MATR_ASSE', repk=matrRigiMesh)
    call dismoi('NB_EQUA', matrRigi, 'MATR_ASSE', repi=matrRigiNbEqua)
    call dismoi('EXIS_CINE', matrMass, 'MATR_ASSE', repk=answer)
    lKineLoadInMass = answer .eq. 'OUI'
    call dismoi('NOM_MAILLA', matrMass, 'MATR_ASSE', repk=matrMassMesh)
    call dismoi('NB_EQUA', matrMass, 'MATR_ASSE', repi=matrMassNbEqua)
    if (lDamp) then
        call dismoi('EXIS_CINE', matrDamp, 'MATR_ASSE', repk=answer)
        lKineLoadInDamp = answer .eq. 'OUI'
        call dismoi('NOM_MAILLA', matrDamp, 'MATR_ASSE', repk=matrDampMesh)
        call dismoi('NB_EQUA', matrDamp, 'MATR_ASSE', repi=matrDampNbEqua)
    end if

! - Some parameters from kineLoad
    kineLoadMesh = ' '
    if (lKineLoad) then
        call dismoi('NOM_MAILLA', kineLoadName, 'CHARGE', repk=kineLoadMesh)
    end if

! - Only for Newmark
    if (present(integScheme_)) then
        if (lKineLoadInRigi .and. &
            (integScheme_ .ne. 1 .and. integScheme_ .ne. 2)) then
            call utmess('F', 'DYNALINE2_1')
        end if
    end if

! - Check consistency between matrices
    if (matrRigiMesh .ne. matrMassMesh) then
        call utmess('F', 'DYNALINE2_2')
    end if
    if (matrRigiNbEqua .ne. matrMassNbEqua) then
        call utmess('F', 'DYNALINE2_3')
    end if
    if (lDamp) then
        if (matrRigiMesh .ne. matrDampMesh) then
            call utmess('F', 'DYNALINE2_2')
        end if
        if (matrRigiNbEqua .ne. matrDampNbEqua) then
            call utmess('F', 'DYNALINE2_3')
        end if
    end if

! - Check consistency between matrices for kinematic loads
    if (lKineLoadInRigi) then
        if (lKineLoadInMass) then
            if (.not. idenob(matrRigi//'.CCID', matrMass//'.CCID')) then
                call utmess('F', 'DYNALINE2_4')
            end if
        else
            call utmess('F', 'DYNALINE2_5')
        end if
        if (lDamp) then
            if (lKineLoadInDamp) then
                if (.not. idenob(matrRigi//'.CCID', matrDamp//'.CCID')) then
                    call utmess('F', 'DYNALINE2_4')
                end if
            else
                call utmess('F', 'DYNALINE2_5')
            end if
        end if
    else
        if (lKineLoadInMass) then
            call utmess('F', 'DYNALINE2_6')
        end if
        if (lDamp) then
            if (lKineLoadInDamp) then
                call utmess('F', 'DYNALINE2_6')
            end if
        end if
    end if

! - Check consistency between kinematic load and kinematic load in matrices
! - Consistency between load from user and load in matrix
    if (lKineLoadInRigi) then
        ! There is at least one subdomain with a kinematic load
        if (lKineLoadGlobal) then
            ! Check the supporting mesh on the subdomain with a kinematic load
            if (lKineLoad .and. kineLoadMesh .ne. matrRigiMesh) then
                call utmess('F', 'DYNALINE2_7')
            end if
        else
            call utmess('F', 'DYNALINE2_9')
        end if
    else
        if (lKineLoad) then
            call utmess('F', 'DYNALINE2_10')
        end if
    end if

! - Convert in nodal field
    if (lKineLoad) then
        kineLoad = '&&COMDLT.KINELOAD'
    end if

!
end subroutine
