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
! ==================================================================================================
!
! Module to read External State Variables
!
! ==================================================================================================
!
module ExternalStateVariableRead_module
! ==================================================================================================
    use ExternalStateVariable_type
    use ExternalStateVariable_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: readDataFromUser, freeDataFromUser
    private :: getTotalNumberOfCmp, getIndxInCata, getRefeValue, getAffeType
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterc/getfac.h"
#include "asterc/r8vide.h"
#include "asterfort/ExternalStateVariable_type.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/dismoi.h"
#include "asterfort/utmess.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! readDataFromUser
!
! Read external state variables from user
!
! Ptr exteVariCata     : catalog of external state variables
! out exteVariAffe     : external state variables from user
!
! --------------------------------------------------------------------------------------------------
    subroutine readDataFromUser(exteVariCata, exteVariAffe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
        type(EXTE_VARI_AFFE), intent(out) :: exteVariAffe
! ----- Locals
        character(len=16), parameter :: factorKeyword = "AFFE_VARC"
        integer(kind=8) :: nbFactorKeyword, nbCmpTotal, nocc
        integer(kind=8) :: iFactorKeyword
        character(len=8) :: exteVariName
        type(EXTE_VARI_DESC) :: exteVariDesc
!   ------------------------------------------------------------------------------------------------
!
! ----- Number of assigned external state variables by user
        call getfac(factorKeyword, nbFactorKeyword)
        exteVariAffe%nbAffe = nbFactorKeyword

! ----- Count total number of components
        call getTotalNumberOfCmp(exteVariCata, factorKeyword, nbFactorKeyword, nbCmpTotal)
        exteVariAffe%nbCmpTotal = nbCmpTotal

! ----- Get properties
        if (nbFactorKeyword .ne. 0) then
            allocate (exteVariAffe%exteVariList(nbFactorKeyword))
            do iFactorKeyword = 1, nbFactorKeyword
! ------------- Type of external state variable
                call getvtx(factorKeyword, 'NOM_VARC', iocc=iFactorKeyword, &
                            scal=exteVariName, nbret=nocc)
                ASSERT(nocc .eq. 1)

! ------------- Get index in catalog
                call getIndxInCata(exteVariCata, exteVariName, exteVariDesc)

! ------------- Get reference value
                call getRefeValue(factorKeyword, iFactorKeyword, exteVariDesc)

! ------------- Get type of assignation
                call getAffeType(factorKeyword, iFactorKeyword, exteVariCata, exteVariDesc)

! ------------- Assign
                exteVariAffe%exteVariList(iFactorKeyword) = exteVariDesc

            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getTotalNumberOfCmp
!
! Get total number of components
!
! Ptr exteVariCata     : cataolog of external state variables
! In  factorKeword     : factor keyword to read
! In  nbFactorKeyword  : number of factor keywords
! Out nbCmpTotal       : total number of components for external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine getTotalNumberOfCmp(exteVariCata, factorKeyword, nbFactorKeyword, nbCmpTotal)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
        character(len=16), intent(in) :: factorKeyword
        integer(kind=8), intent(in) :: nbFactorKeyword
        integer(kind=8), intent(out) :: nbCmpTotal
! ----- Locals
        integer(kind=8) :: iExteVariCata, iFactorKeyword
        integer(kind=8) :: nbCmp, nocc
        aster_logical :: l_found
        character(len=8) :: exteVariName
!   ------------------------------------------------------------------------------------------------
!
        nbCmpTotal = 0
        do iExteVariCata = 1, EXTEVARI_NB_MAXI
            nbCmp = exteVariCata(iExteVariCata)%nbCmp
            l_found = ASTER_FALSE
            do iFactorKeyword = 1, nbFactorKeyword
                call getvtx(factorKeyword, 'NOM_VARC', iocc=iFactorKeyword, &
                            scal=exteVariName, nbret=nocc)
                ASSERT(nocc .eq. 1)
                if (exteVariName .eq. exteVariCata(iExteVariCata)%name) then
                    l_found = ASTER_TRUE
                    exit
                end if
            end do
            if (l_found) then
                nbCmpTotal = nbCmpTotal+nbCmp
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getIndxInCata
!
! Get index in catalog of external state variables
!
! Ptr exteVariCata     : cataolog of external state variables
! In  exteVariName     : name of external state variable
! IO  exteVariDesc     : descriptor of external state variable
!
! --------------------------------------------------------------------------------------------------
    subroutine getIndxInCata(exteVariCata, exteVariName, exteVariDesc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
        character(len=8), intent(in) :: exteVariName
        type(EXTE_VARI_DESC), intent(inout) :: exteVariDesc
! ----- Locals
        integer(kind=8) :: iExteVariCata, cataIndx
!   ------------------------------------------------------------------------------------------------
!
        cataIndx = 0
        do iExteVariCata = 1, EXTEVARI_NB_MAXI
            if (exteVariName .eq. exteVariCata(iExteVariCata)%name) then
                cataIndx = iExteVariCata
                exit
            end if
        end do
        ASSERT(cataIndx .gt. 0)
        exteVariDesc%cataIndx = cataIndx
        exteVariDesc%exteVariName = exteVariName
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getRefeValue
!
! Get reference value
!
! In  factorKeword     : factor keyword to read
! In  iFactorKeyword   : index of factor keyword to read
! IO  exteVariDesc     : descriptor of external state variable
!
! --------------------------------------------------------------------------------------------------
    subroutine getRefeValue(factorKeyword, iFactorKeyword, exteVariDesc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: factorKeyword
        integer(kind=8), intent(in) :: iFactorKeyword
        type(EXTE_VARI_DESC), intent(inout) :: exteVariDesc
! ----- Locals
        integer(kind=8) :: nocc
        real(kind=8):: valeRefe
!   ------------------------------------------------------------------------------------------------
!
        call getvr8(factorKeyword, 'VALE_REF', iocc=iFactorKeyword, &
                    scal=valeRefe, nbret=nocc)
        if (nocc .eq. 0) then
            valeRefe = r8vide()
        end if
        exteVariDesc%valeRefe = valeRefe
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getAffeType
!
! Get type of assignment
!
! In  factorKeword     : factor keyword to read
! In  iFactorKeyword   : index of factor keyword to read
! Ptr exteVariCata     : catalog of external state variables
! IO  exteVariDesc     : descriptor of external state variable
!
! --------------------------------------------------------------------------------------------------
    subroutine getAffeType(factorKeyword, iFactorKeyword, exteVariCata, exteVariDesc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: factorKeyword
        integer(kind=8), intent(in) :: iFactorKeyword
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
        type(EXTE_VARI_DESC), intent(inout) :: exteVariDesc
! ----- Locals
        integer(kind=8) :: n1, n2, cataIndx
        character(len=8) :: affeType, exteVariName
        character(len=8) :: funcResult, resultName, fieldName, dsName
        character(len=16) :: funcExtrLeft, funcExtrRight, fieldType
        character(len=8) :: physQuantity, physQuantityUser
!   ------------------------------------------------------------------------------------------------
!
        cataIndx = exteVariDesc%cataIndx
        ASSERT(cataIndx .gt. 0)
        dsName = " "
        resultName = " "
        fieldName = " "
        affeType = " "
        call getvid(factorKeyword, 'CHAM_GD', iocc=iFactorKeyword, scal=fieldName, nbret=n1)
        call getvid(factorKeyword, 'EVOL', iocc=iFactorKeyword, scal=resultName, nbret=n2)
        ASSERT(n1+n2 .le. 1)
        if (n1 .eq. 1) then
            affeType = 'CHAMP'
            dsName = fieldName
        else if (n2 .eq. 1) then
            affeType = 'EVOL'
            dsName = resultName
        else
            ASSERT(ASTER_FALSE)
        end if
        exteVariDesc%affeType = affeType
        exteVariDesc%dsName = dsName

! ----- Get parameters
        fieldType = " "
        funcExtrLeft = " "
        funcExtrRight = " "
        funcResult = " "
        if (affeType .eq. 'CHAMP') then
            fieldType = exteVariCata(cataIndx)%fieldType
            physQuantity = exteVariCata(cataIndx)%physQuantity
            call dismoi('NOM_GD', fieldName, 'CHAMP', repk=physQuantityUser)
            if (physQuantity .ne. physQuantityUser) then
                exteVariName = exteVariDesc%exteVariName
                call utmess('F', 'MATERIAL2_50', &
                            nk=3, valk=[exteVariName, physQuantity, physQuantityUser])
            end if
        elseif (affeType .eq. 'EVOL') then
            call getvtx(factorKeyword, 'PROL_GAUCHE', iocc=iFactorKeyword, &
                        scal=funcExtrLeft, nbret=n1)
            call getvtx(factorKeyword, 'PROL_DROITE', iocc=iFactorKeyword, &
                        scal=funcExtrRight, nbret=n1)
            call getvid(factorKeyword, 'FONC_INST', iocc=iFactorKeyword, &
                        scal=funcResult, nbret=n1)
            if (n1 .eq. 0) then
                funcResult = ' '
            end if
            exteVariDesc%funcExtrLeft = funcExtrLeft
            exteVariDesc%funcExtrRight = funcExtrRight
            exteVariDesc%funcResult = funcResult
            call getvtx(factorKeyword, 'NOM_CHAM', iocc=iFactorKeyword, &
                        scal=fieldType, nbret=n1)
            if (n1 .eq. 0) then
                fieldType = exteVariCata(cataIndx)%fieldType
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
        exteVariDesc%fieldType = fieldType
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! freeDataFromUser
!
! Free data about external state variables from user
!
! IO exteVariAffe     : external state variables from user
!
! --------------------------------------------------------------------------------------------------
    subroutine freeDataFromUser(exteVariAffe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_AFFE), intent(inout) :: exteVariAffe
! ----- Locals
!   ------------------------------------------------------------------------------------------------
!
        if (exteVariAffe%nbAffe .ne. 0) then
            deallocate (exteVariAffe%exteVariList)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module ExternalStateVariableRead_module
