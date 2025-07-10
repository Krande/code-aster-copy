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
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.a
! --------------------------------------------------------------------
!
subroutine compareFieldShape(fieldModelZ, fieldZ, &
                             projectOnLigrel, paraNameZ, &
                             iret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/cestas.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: fieldModelZ, fieldZ
    aster_logical, intent(in) :: projectOnLigrel
    character(len=*), intent(in) :: paraNameZ
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! Compare shape of field with another one
!
! --------------------------------------------------------------------------------------------------
!
! In  fieldModel       : field for model
! In  field            : field to compare
! In  projectOnLigrel  : project field on FED from model
! In  paraName         : name of parameter to project field
! Out iret             : error code
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: fieldS = '&&SGCOMP.FIELD'
    character(len=19), parameter :: fieldModelS = '&&SGCOMP.FIELDS'
    character(len=24) :: modelLigrel
    character(len=8) :: meshRefe, mesh, discType, paraName
    character(len=1) :: scalarType
    character(len=1), parameter :: jeveuxBase = "G"
    aster_logical :: noSameNpg, noSameCmp, noSameNspg
    integer(kind=8) :: iadRefe, iad, iretZero
    integer(kind=8) :: nncp, iCell, nbCell
    integer(kind=8) :: npgRefe, nspgRefe, ncmpRefe
    integer(kind=8) :: npg, nspg, ncmp
    integer(kind=8) :: jvCesdRefe, jvCeslRefe, jvCesvRefe
    integer(kind=8) :: jvCesd, jvCesl, jvCesv
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    iret = 0
    paraName = paraNameZ

! - Checks
    call dismoi("TYPE_SCA", fieldZ, "CHAMP", repk=scalarType)
    if (scalarType .ne. "R") then
        call utmess("F", "FIELD0_31")
    end if
    call dismoi('TYPE_CHAMP', fieldZ, 'CHAMP', repk=discType)
    if (discType .ne. "ELGA") then
        call utmess("F", "FIELD0_32")
    end if
    call dismoi('NOM_MAILLA', fieldModelZ, 'CHAMP', repk=meshRefe)
    call dismoi('NOM_MAILLA', fieldZ, 'CHAMP', repk=mesh)
    if (meshRefe .ne. mesh) then
        call utmess('F', 'FIELD0_33')
    end if

! - Create reduced fields
    call celces(fieldModelZ, 'V', fieldModelS)
    call cestas(fieldModelS)
    call celces(fieldZ, 'V', fieldS)
    call cestas(fieldS)

! - Access to reduced fields
    call jeveuo(fieldModelS//'.CESD', 'L', jvCesdRefe)
    call jeveuo(fieldModelS//'.CESL', 'L', jvCeslRefe)
    call jeveuo(fieldModelS//'.CESV', 'L', jvCesvRefe)
    call jeveuo(fieldS//'.CESD', 'L', jvCesd)
    call jeveuo(fieldS//'.CESL', 'L', jvCesl)
    call jeveuo(fieldS//'.CESV', 'E', jvCesv)
    nbCell = zi(jvCesd-1+1)

! - Check
    noSameNpg = ASTER_FALSE
    noSameCmp = ASTER_FALSE
    noSameNspg = ASTER_FALSE

    do iCell = 1, nbCell
        call cesexi('C', jvCesdRefe, jvCeslRefe, iCell, 1, 1, 1, iadRefe)
        if (iadRefe .gt. 0) then
            npgRefe = zi(jvCesdRefe-1+5+4*(iCell-1)+1)
            nspgRefe = zi(jvCesdRefe-1+5+4*(iCell-1)+2)
            ncmpRefe = zi(jvCesdRefe-1+5+4*(iCell-1)+3)

            call cesexi('C', jvCesd, jvCesl, iCell, 1, 1, 1, iad)
            if (iad .gt. 0) then
                npg = zi(jvCesd-1+5+4*(iCell-1)+1)
                nspg = zi(jvCesd-1+5+4*(iCell-1)+2)
                ncmp = zi(jvCesd-1+5+4*(iCell-1)+3)

! ------------- Check number of Gauss points
                if (npgRefe .ne. 0 .and. npg .ne. 0) then
                    if (npgRefe .ne. npg) then
                        noSameNpg = ASTER_TRUE
                        cycle
                    end if
                end if

! ------------- Check number of Gauss sub-points
                if (nspgRefe .ne. 0 .and. nspg .ne. 0) then
                    if (nspgRefe .ne. nspg) then
                        noSameNspg = ASTER_TRUE
                        cycle
                    end if
                end if

! ------------- Check number of components
                if (ncmpRefe .ne. 0 .and. ncmp .ne. 0) then
                    if (ncmpRefe .ne. ncmp) then
                        noSameCmp = ASTER_TRUE
                        cycle
                    end if
                end if
            end if
        end if
    end do

! - Conclusion
    if (noSameCmp) then
        iret = 1
        call utmess("I", 'FIELD0_40')
    end if
    if (noSameNpg) then
        iret = 1
        call utmess("I", 'FIELD0_41')
    end if
    if (noSameNspg) then
        iret = 1
        call utmess("I", 'FIELD0_42')
    end if

! - Project field
    if (projectOnLigrel) then
        call dismoi('NOM_LIGREL', fieldModelZ, 'CHAMP', repk=modelLigrel)
        if (paraName .eq. " ") then
            call utmess("F", "FIELD0_34")
        end if
        call cescel(fieldS, modelLigrel, 'TOU_INI_ELGA', paraName, "OUI", &
                    nncp, jeveuxBase, fieldZ, "F", iretZero)
        if (iretZero .ne. 0) then
            call utmess("F", "FIELD0_35")
        end if
    end if

! - Clean
    call detrsd('CHAM_ELEM_S', fieldS)
    call detrsd('CHAM_ELEM_S', fieldModelS)
!
    call jedema()
!
end subroutine
