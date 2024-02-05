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
module MetallurgyOperator_module
! ==================================================================================================
    use Metallurgy_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: metaGetInitialState, metaPrepTRCWorkingField
    public :: metaComputeInitialField, metaPrepareInitialState, metaCompMetaElno
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/chpver.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/mecact.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/rcadme.h"
#include "asterfort/rs_getnume.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! metaGetInitialState
!
! Get field for initial phases
!
! In  resultName       : name of result datastructure
! Out metaInit         : name of field for initial phases
! Out numeFieldInit    : storing index of initial field in result datastructure
!                        (0 if already exists)
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetInitialState(resultName, metaInit, numeFieldInit)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: resultName
        character(len=24), intent(out) :: metaInit
        integer, intent(out) :: numeFieldInit
! ----- Local
        character(len=8) :: resultInit
        character(len=24) :: selectCrit, metaFieldUser
        character(len=16), parameter :: factorKeyword = "ETAT_INIT"
        integer :: nbRet, iret, nbtrou, numeStoreInit
        real(kind=8) :: timeInit, selectTole
!   ------------------------------------------------------------------------------------------------
!
        metaInit = "Unknown"
        numeFieldInit = 0
        call getvid(factorKeyword, 'META_INIT_ELNO', iocc=1, scal=metaFieldUser, nbret=nbRet)
        if (nbRet .gt. 0) then
! --------- Check type
            call chpver('F', metaFieldUser, 'CART', 'VAR2_R', iret)
            metaInit = "&&SMEVOL_ZINIT"
            call copisd('CHAMP_GD', 'V', metaFieldUser, metaInit)
        else
            call getvid(factorKeyword, 'EVOL_THER', iocc=1, scal=resultInit, nbret=nbRet)
            if (resultInit .ne. resultName) then
                call utmess('F', 'META1_2')
            end if
            call getvis(factorKeyword, 'NUME_INIT', iocc=1, scal=numeStoreInit, nbret=nbRet)
            if (nbRet .eq. 0) then
                call getvr8(factorKeyword, 'INST_INIT', iocc=1, scal=timeInit, nbret=nbRet)
                call getvr8(factorKeyword, 'PRECISION', iocc=1, scal=selectTole, nbret=nbRet)
                call getvtx(factorKeyword, 'CRITERE', iocc=1, scal=selectCrit, nbret=nbRet)
                call rs_getnume(resultName, timeInit, selectCrit, selectTole, numeStoreInit, nbtrou)
                if (nbtrou .eq. 0 .or. nbtrou .gt. 1) then
                    call utmess('F', 'META1_51', sr=timeInit)
                end if
            end if
            call rsexch('F', resultName, 'META_ELNO', numeStoreInit, metaInit, iret)
            if (iret .ne. 0) then
                call utmess('F', 'META1_52')
            end if
            numeFieldInit = numeStoreInit
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaPrepTRCWorkingField
!
! Prepare field to manage TRC curves in elementary computation
!
! In  model            : model
! In  materialField    : field for material parameters
! In  chftrc           : name of field to manage TRC curves
!
! --------------------------------------------------------------------------------------------------
    subroutine metaPrepTRCWorkingField(model, materialField, chftrc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: model, materialField
        character(len=24), intent(out) :: chftrc
! ----- Local
        integer :: iMaterVale, materValeLen
        character(len=8), pointer :: materVale(:) => null()
        integer :: nbhist, icodre, jvDummy
        character(len=8) :: materPara
        aster_logical :: hasTRC
        character(len=24) :: modelLigrel
        integer, parameter :: nbPara = 2
        character(len=8), parameter :: paraName(nbPara) = (/'I1  ', 'I2  '/)
        integer :: iadtrc(nbPara), adrsJv(nbPara)
!   ------------------------------------------------------------------------------------------------
!
        chftrc = " "
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! ----- Get access to material field
        call jeveuo(materialField(1:8)//'.CHAMP_MAT .VALE', 'E', vk8=materVale)
        call jelira(materialField(1:8)//'.CHAMP_MAT .VALE', 'LONMAX', materValeLen)

! ----- Get number of histories in TRC diagram
        nbhist = 0
        hasTRC = ASTER_FALSE
        iadtrc = 0
        do iMaterVale = 1, materValeLen
            materPara = materVale(iMaterVale)
            if (materPara .ne. '        ') then
                call rcadme(materPara, 'META_ACIER', 'TRC', iadtrc, icodre, 0)
                if (icodre .eq. 0) hasTRC = ASTER_TRUE
                nbhist = max(nbhist, iadtrc(1))
            end if
        end do

! ----- Create working field for TRC
        if (hasTRC) then
            call wkvect('&&SMEVOL_FTRC', 'V V R', 9*nbhist, jvDummy)
            call wkvect('&&SMEVOL_TRC', 'V V R', 15*nbhist, jvDummy)
            call jeveut('&&SMEVOL_FTRC', "E", adrsJv(1))
            call jeveut('&&SMEVOL_TRC', "E", adrsJv(2))
            chftrc = '&&SMEVOL.ADRESSES'
            call mecact('V', chftrc, 'LIGREL', modelLigrel, 'ADRSJEVN', &
                        ncmp=nbPara, lnomcmp=paraName, vi=adrsJv)
        else
            chftrc = " "
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaComputeInitialField
!
! Compute initial field of metallurgy
!
! From the map describing the initial phases provided by the user, we create the complete map
! with the other necessary parameters.
!
! In  model            : model
! In  materialCoding   : field for coding material parameters
! In  comporMeta       : map of behaviour for metallurgy
! In  temp             : temperature
! In  metaInitUser     : map describing the initial phases provided by the user
! In  metaInit         : name of map for initial metallurgy computation
!
! --------------------------------------------------------------------------------------------------
    subroutine metaComputeInitialField(model, materialCoding, comporMeta, temp, &
                                       metaInitUser, metaInit)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: model
        character(len=24), intent(in) :: materialCoding, comporMeta, temp, metaInitUser
        character(len=24), intent(in) :: metaInit
! ----- Local
        integer, parameter :: nbout = 1, nbin = 4
        character(len=8) :: lpaout(nbout), lpain(nbin)
        character(len=24) :: lchin(nbin)
        character(len=24) :: modelLigrel
        character(len=1), parameter :: base = "V"
        character(len=16), parameter :: option = "META_INIT_ELNO"
!   ------------------------------------------------------------------------------------------------
!
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
        lpain(1) = 'PMATERC'
        lchin(1) = materialCoding
        lpain(2) = 'PCOMPOR'
        lchin(2) = comporMeta
        lpain(3) = 'PTEMPER'
        lchin(3) = temp
        lpain(4) = 'PPHASIN'
        lchin(4) = metaInitUser
        lpaout(1) = 'PPHASNOU'
        call calcul('S', option, modelLigrel, nbin, lchin, &
                    lpain, nbout, metaInit, lpaout, base, &
                    'OUI')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaPrepareInitialState
!
! Prepare initial state
!
! In  resultName       : name of result datastructure
! In  numeStore        : current index of stroing to proceed
! In  model            : model
! In  materialCoding   : field for coding material parameters
! In  comporMeta       : map of behaviour for metallurgy
! In  metaInitUser     : map describing the initial phases provided by the user
!
! --------------------------------------------------------------------------------------------------
    subroutine metaPrepareInitialState(resultName, numeStore, &
                                       model, materialCoding, comporMeta, metaInitUser)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: resultName
        integer, intent(in) :: numeStore
        character(len=8), intent(in) :: model
        character(len=24), intent(in) :: materialCoding, comporMeta, metaInitUser
! ----- Local
        integer :: iret, jvPara
        character(len=24) :: metaInit, temp, resultField
        real(kind=8) :: time
!   ------------------------------------------------------------------------------------------------
!

! ----- Prepare metallurgical field (dynamic field)
        metaInit = '&&SMEVOL.PHAS_META1'
        call copisd('CHAM_ELEM_S', 'V', comporMeta, metaInit)

! ----- Get temperature field
        call rsexch('F', resultName, 'TEMP', numeStore, temp, iret)

! ----- Get current time
        call rsadpa(resultName, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
        time = zr(jvPara)

! ----- Compute initial field of metallurgy
        call metaComputeInitialField(model, materialCoding, comporMeta, temp, &
                                     metaInitUser, metaInit)

! ----- Save map for initial metallurgy computation in result datastructure
        call rsexch(' ', resultName, 'META_ELNO', numeStore, resultField, iret)
        call copisd('CHAMP_GD', 'G', metaInit, resultField)
        call rsnoch(resultName, 'META_ELNO', numeStore)
        call utmess('I', 'ARCHIVAGE_6', sk='META_ELNO', si=numeStore, sr=time)

! ----- Save Save COMPORMETA in result datastructure
        call rsexch(' ', resultName, 'COMPORMETA', numeStore, resultField, iret)
        call copisd('CHAMP_GD', 'G', comporMeta, resultField)
        call rsnoch(resultName, 'COMPORMETA', numeStore)
        call utmess('I', 'ARCHIVAGE_6', sk='COMPORMETA', si=numeStore, sr=time)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaCompMetaElno
!
! Compute option META_ELNO
!
! In  resultName       : name of result datastructure
! In  model            : model
! In  materialCoding   : field for coding material parameters
! In  comporMeta       : map of behaviour for metallurgy
! In  chftrc           : name of field to manage TRC curves
! In  metaIn           : input field of metallurgy
! Out metaOut          : output field of metallurgy
!
! --------------------------------------------------------------------------------------------------
    subroutine metaCompMetaElno(model, materialCoding, &
                                comporMeta, chftrc, &
                                time_0, time_1, time_2, &
                                temp_0, temp_1, temp_2, &
                                metaIn, metaOut)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: model
        character(len=24), intent(in) :: materialCoding, comporMeta, chftrc
        real(kind=8), intent(in) :: time_0, time_1, time_2
        character(len=24), intent(in) :: temp_0, temp_1, temp_2
        character(len=24), intent(in) :: metaIn
        character(len=24), intent(out) :: metaOut
! ----- Local
        integer, parameter :: nbout = 1, nbin = 8
        character(len=8) :: lpaout(nbout), lpain(nbin)
        character(len=24) :: lchin(nbin), lchout(nbout)
        character(len=24) :: modelLigrel
        character(len=1), parameter :: base = "V"
        character(len=16), parameter :: option = "META_ELNO"
        character(len=24), parameter :: chtime = '&&SMEVOL.CH_INST_R'
        integer, parameter :: nbTimePara = 3
        character(len=8), parameter :: timeParaName(nbTimePara) = &
                                       (/'INST    ', 'DELTA01 ', 'DELTA12 '/)
        real(kind=8) :: timeParaVale(nbTimePara)
!   ------------------------------------------------------------------------------------------------
!
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
        metaOut = "&&SMEVOL.PHAS_META3"

! ----- Create map for time
        timeParaVale(1) = time_1
        timeParaVale(2) = time_1-time_0
        timeParaVale(3) = time_2-time_1
        call mecact('V', chtime, 'MODELE', modelLigrel, 'INST_R  ', &
                    ncmp=nbTimePara, lnomcmp=timeParaName, vr=timeParaVale)

! ----- Input fields
        lpain(1) = 'PMATERC'
        lchin(1) = materialCoding
        lpain(2) = 'PCOMPOR'
        lchin(2) = comporMeta
        lpain(3) = 'PTEMPAR'
        lchin(3) = temp_0
        lpain(4) = 'PTEMPER'
        lchin(4) = temp_1
        lpain(5) = 'PTEMPIR'
        lchin(5) = temp_2
        lpain(6) = 'PTIMMTR'
        lchin(6) = chtime
        lpain(7) = 'PPHASIN'
        lchin(7) = metaIn
        lpain(8) = 'PFTRC'
        lchin(8) = chftrc

! ----- Output field
        lpaout(1) = 'PPHASNOU'
        lchout(1) = metaOut

! ----- Field with dynamical size (VARI_R) => allocate from compor map
        call copisd('CHAM_ELEM_S', 'V', comporMeta, metaOut)

! ----- Compute
        call calcul('S', option, modelLigrel, &
                    nbin, lchin, lpain, &
                    nbout, lchout, lpaout, &
                    base, 'OUI')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module MetallurgyOperator_module
