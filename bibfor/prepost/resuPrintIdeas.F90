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
subroutine resuPrintIdeas(fileUnit, dsName, lResu, &
                          storeNb, storeListIndx, &
                          fieldListNb, fieldListType, &
                          title_, titleKeywf_, &
                          titleKeywfIocc_, realFormat_, &
                          cmpUserNb_, cmpUserName_, &
                          nodeUserNb_, nodeUserNume_, &
                          cellUserNb_, cellUserNume_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/resuPrintIdeasElem.h"
#include "asterfort/resuPrintIdeasNode.h"
#include "asterfort/rsexch.h"
#include "asterfort/titre2.h"
#include "asterfort/jedetr.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: fileUnit
    character(len=*), intent(in) :: dsName
    aster_logical, intent(in) :: lResu
    integer(kind=8), intent(in) :: storeNb, storeListIndx(:)
    integer(kind=8), intent(in) :: fieldListNb
    character(len=*), intent(in) :: fieldListType(*)
    character(len=*), optional, intent(in) :: title_, titleKeywf_
    integer(kind=8), optional, intent(in) :: titleKeywfIocc_
    character(len=*), optional, intent(in) :: realFormat_
    integer(kind=8), optional, intent(in) :: cmpUserNb_
    character(len=8), optional, pointer :: cmpUserName_(:)
    integer(kind=8), optional, intent(in) :: cellUserNb_
    integer(kind=8), optional, pointer :: cellUserNume_(:)
    integer(kind=8), optional, intent(in) :: nodeUserNb_
    integer(kind=8), optional, pointer :: nodeUserNume_(:)
!
! --------------------------------------------------------------------------------------------------
!
! Print result or field in a file (IMPR_RESU)
!
! Print field or result for 'IDEAS'
!
! --------------------------------------------------------------------------------------------------
!
! In  fileUnit         : index of file (logical unit)
! In  dsName           : name of datastructure (result or field)
! In  lResu            : flag if datastructure is a result
! In  title            : title of result
! In  titleKeywf       : keyword for sub-title
! In  titleKeywfIocc   : index of keyword for sub-title
! In  storeNb          : number of storing slots
! In  storeListIndx    : index of stroring slots
! In  fieldListNb      : length of list of fields to save
! In  fieldListType    : list of fields type to save
! In  realFormat       : format of real numbers
! In  cmpUserNb        : number of components to select
! Ptr cmpUserName      : pointer to the names of components to select
! In  cellUserNb       : number of cells require by user
! Ptr cellUserNume     : pointer to the list of index of cells require by user
! In  nodeUserNb       : number of nodes require by user
! Ptr nodeUserNume     : pointer to the list of index of nodes require by user
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: fieldType, fieldName
    character(len=16) :: resultType
    character(len=4) :: fieldSupport
    character(len=24) :: subtitleJvName
    integer(kind=8) :: storeIndx, iret
    integer(kind=8) :: iStore, iField
    integer(kind=8) :: cmpUserNb, nodeUserNb, cellUserNb
    character(len=8), pointer :: cmpUserName(:) => null()
    integer(kind=8), pointer :: nodeUserNume(:) => null()
    integer(kind=8), pointer :: cellUserNume(:) => null()
    character(len=80) :: title
    character(len=8) :: realFormat
    character(len=16) :: titleKeywf
    integer(kind=8) :: titleKeywfIocc
!
! --------------------------------------------------------------------------------------------------
!
    subtitleJvName = '&&IRECRI.SOUS_TITRE.TITR'
    cmpUserNb = 0
    if (present(cmpUserNb_)) then
        cmpUserNb = cmpUserNb_
        cmpUserName => cmpUserName_
    end if
    cellUserNb = 0
    if (present(cellUserNb_)) then
        cellUserNb = cellUserNb_
        cellUserNume => cellUserNume_
    end if
    nodeUserNb = 0
    if (present(nodeUserNb_)) then
        nodeUserNb = nodeUserNb_
        nodeUserNume => nodeUserNume_
    end if
    realFormat = '1PE12.5'
    if (present(realFormat_)) then
        realFormat = realFormat_
    end if
    title = ' '
    if (present(title_)) then
        title = title_
    end if
    titleKeywfIocc = 0
    if (present(titleKeywfIocc_)) then
        titleKeywfIocc = titleKeywfIocc_
    end if
    titleKeywf = ' '
    if (present(titleKeywf_)) then
        titleKeywf = titleKeywf_
    end if
!
! - Loop on storing slots
!
    do iStore = 1, storeNb
        call jemarq()
        call jerecu('V')
        storeIndx = storeListIndx(iStore)
! ----- Loop on fields
        if (fieldListNb .ne. 0) then
            do iField = 1, fieldListNb
                fieldType = fieldListType(iField)
! ------------- Extract field from results if necessary
                fieldName = ' '
                if (lResu) then
                    call rsexch(' ', dsName, fieldType, storeIndx, fieldName, iret)
                    if (iret .ne. 0) then
                        cycle
                    end if
                    call dismoi('TYPE_RESU', dsName, 'RESULTAT', repk=resultType)
                else
                    fieldName = dsName
                end if
! ------------- Check support
                call dismoi('TYPE_CHAMP', fieldName, 'CHAMP', repk=fieldSupport)
                if ((fieldSupport .eq. 'NOEU') .or. (fieldSupport(1:2) .eq. 'EL')) then
                else if (fieldSupport .eq. 'CART') then
                    if (.not. lResu) then
                        call utmess('A', 'RESULT3_3')
                        cycle
                    end if
                else
                    if (resultType .eq. 'MODE_GENE' .or. resultType .eq. 'HARM_GENE') then
                        call utmess('A+', 'RESULT3_2', sk=fieldType)
                    else
                        call utmess('A+', 'RESULT3_1', sk=fieldType)
                    end if
                end if
! ------------- Create sub-title
                if (lResu) then
                    call titre2(dsName, fieldName, subtitleJvName, titleKeywf, titleKeywfIocc, &
                                realFormat, fieldType, storeIndx)
                else
                    call titre2(dsName, fieldName, subtitleJvName, titleKeywf, titleKeywfIocc, &
                                realFormat)
                end if
! ------------- Print field
                if (fieldSupport .eq. 'NOEU') then
                    call resuPrintIdeasNode(fileUnit, dsName, &
                                            title, storeIndx, &
                                            fieldType, fieldName, &
                                            cmpUserNb, cmpUserName, &
                                            nodeUserNb, nodeUserNume)
                else
                    call resuPrintIdeasElem(fileUnit, dsName, &
                                            title, storeIndx, &
                                            fieldType, fieldName, &
                                            cmpUserNb, cmpUserName, &
                                            cellUserNb, cellUserNume)
                end if
            end do
        end if
        call jedema()
    end do
!
! - Clean
!
    call jedetr(subtitleJvName)
end subroutine
