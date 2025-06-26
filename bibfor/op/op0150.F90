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
subroutine op0150()
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lect58.h"
#include "asterfort/lrcomm.h"
#include "asterfort/lridea.h"
#include "asterfort/rsmode.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/resuReadCheckFields.h"
#include "asterfort/resuReadMed.h"
#include "asterfort/resuReadStorageAccess.h"
#include "asterfort/resuReadDebug.h"
#include "asterfort/resuReadCreateREFD.h"
#include "asterfort/resuReadGetParameters.h"
#include "asterfort/resuSaveParameters.h"
#include "asterfort/resuGetLoads.h"
#include "asterfort/resuGetEmpiricParameters.h"
#include "asterfort/resuReadPrepareDatastructure.h"
!
! --------------------------------------------------------------------------------------------------
!
! LIRE_RESU
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: storeIndxNb, storeCreaNb, storeTimeNb
    character(len=19) :: storeTime, storeIndx
    character(len=8) :: storeCrit
    real(kind=8) :: storeEpsi
    character(len=4) :: storePara
    character(len=10) :: storeAccess
    integer(kind=8) :: empiNumePlan, empiSnapNb
    character(len=24) :: empiFieldType
    integer(kind=8) :: fieldStoreNb(100)
    character(len=16) :: fieldList(100)
    integer(kind=8) :: nbOcc, iField, fieldNb
    character(len=8) :: resultName, resultNameReuse
    character(len=8) :: meshAst, model, caraElem, materField
    character(len=8) :: matrRigi, matrMass
    character(len=16) :: nomcmd, resultType2, resultType
    integer(kind=8) :: fileUnit, nbret
    character(len=16) :: fileFormat
    character(len=24) :: listLoad
    aster_logical:: lReuse, lLireResu, lVeriVari
    character(len=8) :: answer
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
! - Initializations
!
    lLireResu = ASTER_TRUE
!
! - Get results datastructure
!
    call getres(resultName, resultType2, nomcmd)
    call getvtx(' ', 'TYPE_RESU', scal=resultType, nbret=nbOcc)
    ASSERT(resultType .eq. resultType2)
!
! - Reuse mode
!
    lReuse = ASTER_FALSE
    if (getexm(' ', 'RESULTAT') .eq. 1) then
        call getvid(' ', 'RESULTAT', scal=resultNameReuse, nbret=nbOcc)
        if (nbOcc .ne. 0) then
            lReuse = ASTER_TRUE
            if (resultName .ne. resultNameReuse) then
                call utmess('F', 'SUPERVIS2_79', sk='RESULTAT')
            end if
        end if
    end if
!
! - Output format
!
    call getvtx(' ', 'FORMAT', scal=fileFormat, nbret=nbOcc)
    ASSERT(nbOcc .eq. 1)
    call getvis(' ', 'UNITE', scal=fileUnit, nbret=nbOcc)
    ASSERT(nbOcc .eq. 1)
!
! - Get list of fields to read
!
    call getfac('FORMAT_MED', nbOcc)
    if (nbOcc .gt. 0) then
        fieldNb = nbOcc
        if (fieldNb .gt. 100) then
            call utmess('F', 'UTILITAI2_86')
        else
            do iField = 1, fieldNb
                call getvtx('FORMAT_MED', 'NOM_CHAM', iocc=iField, scal=fieldList(iField), &
                            nbret=nbOcc)
                ASSERT(nbOcc .eq. 1)
            end do
        end if
    else
        call getvtx(' ', 'NOM_CHAM', nbval=100, vect=fieldList, nbret=fieldNb)
        if (fieldNb .lt. 0) then
            call utmess('F', 'UTILITAI2_86')
        end if
    end if

! - Get standard parameters
    call resuReadGetParameters(meshAst, model, caraElem, materField)

! - Get loads/BC and create list of loads
    call resuGetLoads(model, resultType, listLoad)
!
! - Get parameters for MODE_EMPI
!
    call resuGetEmpiricParameters(resultType, fieldNb, fieldList, &
                                  empiNumePlan, empiSnapNb, empiFieldType)
!
! - Get storage access from command file
!
    call resuReadStorageAccess(storeAccess, &
                               storeIndxNb, storeIndx, &
                               storeTimeNb, storeTime, &
                               storeEpsi, storeCrit)
!
! - Prepare datastructure
!
    call resuReadPrepareDatastructure(resultName, resultType, lReuse, &
                                      storeIndxNb, storeTimeNb, &
                                      storeIndx, storeTime, &
                                      storeCreaNb, storePara)
!
! - Create .REFD object and save matrices (dynamic results)
!
    matrRigi = ' '
    matrMass = ' '
    if (resultType(1:9) .eq. 'MODE_MECA') then
        call getvid(' ', 'MATR_RIGI', scal=matrRigi, nbret=nbOcc)
        if (nbOcc .eq. 0) then
            matrRigi = ' '
        end if
        call getvid(' ', 'MATR_MASS', scal=matrMass, nbret=nbOcc)
        if (nbOcc .eq. 0) then
            matrMass = ' '
        end if
    end if
    call resuReadCreateREFD(resultName, resultType, matrRigi, matrMass)
!
! - Check if fields are allowed for the result
!
    call resuReadCheckFields(resultName, resultType, fieldNb, fieldList)
!
! - Read
!
    if (fileFormat .eq. 'IDEAS') then
        call lridea(fileUnit, &
                    resultName, resultType, &
                    model, meshAst, &
                    fieldNb, fieldList, &
                    storeAccess, &
                    storeIndxNb, storeTimeNb, &
                    storeIndx, storeTime, &
                    storeCrit, storeEpsi, &
                    storePara)
        fieldStoreNb = storeCreaNb
    else if (fileFormat .eq. 'IDEAS_DS58') then
        call lect58(fileUnit, &
                    resultName, resultType, meshAst, &
                    fieldNb, fieldList, &
                    storeAccess, &
                    storeIndxNb, storeTimeNb, &
                    storeIndx, storeTime, &
                    storeCrit, storeEpsi)
        fieldStoreNb = storeCreaNb
    else if (fileFormat .eq. 'MED') then
        call resuReadMed(fileUnit, &
                         resultName, &
                         model, meshAst, &
                         fieldNb, fieldList, &
                         storeAccess, storeCreaNb, &
                         storeIndxNb, storeIndx, &
                         storeTimeNb, storeTime, &
                         storeEpsi, storeCrit, &
                         storePara, fieldStoreNb)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Save standard parameters in results datastructure
!
    call resuSaveParameters(resultName, resultType, &
                            model, caraElem, materField, listLoad, &
                            empiNumePlan, empiSnapNb, empiFieldType)
!
! - Non-linear behaviour management
!
    if (resultType .eq. 'EVOL_NOLI') then
        call getvtx(' ', 'VERI_VARI', scal=answer, nbret=nbret)
        lVeriVari = answer .eq. 'OUI'
        call lrcomm(lReuse, resultName, model, caraElem, materField, lLireResu, lVeriVari)
    end if
!
! - Debug
!
    if (niv .ge. 2) then
        call resuReadDebug(resultName, &
                           fieldNb, fieldList, fieldStoreNb, &
                           storePara, storeEpsi, storeCrit)
    end if
!
! - Save title
!
    call titre()
!
! - Set same numbering for matrices and nodal field
!
    call rsmode(resultName)
!
    call jedema()
end subroutine
