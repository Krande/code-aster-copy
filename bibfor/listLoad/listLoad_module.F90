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
! Module for the management of list of loads
!
! ==================================================================================================
!
module listLoad_module
! ==================================================================================================
    use listLoad_type
    use loadAcouCompute_module
    use loadAcouCompute_type
    use loadMecaCompute_module
    use loadMecaCompute_type
    use loadTherCompute_module
    use loadTherCompute_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: getLoadKeyword, getLoadFunc, getLoadName, getLoadApply, getLoadIndexKeyword
    public :: getNbLoadsFromUser, getNbLoadsFromList
    public :: getLoadKine, getLoadGround, getLoadEvolChar
    public :: getLoadDiri, createUnitFunc
    public :: getLoadPlaneWave, detectMecaNeumLoad
    private :: getLoadIdenNeum
    private :: getLoadFuncUser, getLoadIdenFromList
    public :: creaListLoad, creaListLoadFSIOne, creaListLoadFromList
    public :: getNbLoadType, listLoadDebug, resizeListLoad
    public :: setLoadToList, addLoadToList, nameListLoad
    public :: getMecaNeum, getTherNeum, getAcouNeum
    public :: checkConsistency, getLoadParameters, addLoadMeca, addLoadTher, addLoadAcou
    public :: getListLoadLigrel
    private :: getListLoadAccess, copyListLoad
! ==================================================================================================
    private
#include "asterc/getexm.h"
#include "asterc/getfac.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/focstc.h"
#include "asterfort/focste.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
#include "LoadTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getLoadKine
!
! Get identifier for standard Neumann loads
!
! In  loadApply         : how to apply load
! In  loadIsFunc        : flag if load is a function
! In  paraIsTime        : flag if function of time
! Out loadIden          : identifier for load in list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadIdenNeum(loadApply, loadIsFunc, paraIsTime, &
                               loadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: loadApply
        aster_logical, intent(in) :: loadIsFunc, paraIsTime
        character(len=24), intent(out) :: loadIden
!   ------------------------------------------------------------------------------------------------
!
        loadIden = "None"
        if (loadApply .eq. 'FIXE_PILO') then
            if (loadIsFunc) then
                loadIden = 'NEUM_PILO_F'
                if (paraIsTime) then
                    call utmess('F', 'CHARGES9_1')
                end if
            else
                loadIden = 'NEUM_PILO'
            end if
        else if (loadApply .eq. 'SUIV_PILO') then
            if (loadIsFunc) then
                loadIden = 'NEUM_SUIP_F'
                if (paraIsTime) then
                    call utmess('F', 'CHARGES9_1')
                end if
            else
                loadIden = 'NEUM_SUIP'
            end if
        else if (loadApply .eq. 'SUIV') then
            loadIden = 'NEUM_SUIV'
        else if (loadApply .eq. 'FIXE_CSTE') then
            if (loadIsFunc) then
                if (paraIsTime) then
                    loadIden = 'NEUM_FT'
                else
                    loadIden = 'NEUM_FO'
                end if
            else
                loadIden = 'NEUM_CSTE'
            end if
        else if (loadApply .eq. 'DIDI') then
            if (loadIsFunc) then
                if (paraIsTime) then
                    loadIden = 'NEUM_FT'
                else
                    loadIden = 'NEUM_FO'
                end if
            else
                loadIden = 'NEUM_CSTE'
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMecaNeum
!
! Get identifier for mechanical Neumann loads
!
! In  staticOperator    : flag for static operator
! In  loadPreObject     : base JEVEUX name for object
! In  loadApply         : how to apply load
! In  loadIsFunc        : flag if load is a function
! IO  nbLoadIden        : number of identifier of loads in listLoadIden
! IO  listLoadIden      : list of loads's identifier
! In  neumExcl          : flag to exclude AFFE_CHAR_MECA (Neumann/by type)
!
! --------------------------------------------------------------------------------------------------
    subroutine getMecaNeum(staticOperator, loadPreObjectZ, &
                           loadApply, loadIsFunc, &
                           nbLoadIden, listLoadIden, &
                           neumExcl_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: staticOperator
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=16), intent(in) :: loadApply
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(inout) :: nbLoadIden
        character(len=24), intent(inout) :: listLoadIden(LOAD_NBIDEN_MAXI)
        aster_logical, optional, intent(in) :: neumExcl_(LOAD_NEUM_NBTYPE)
! ----- Locals
        integer(kind=8) :: iTypeNeum
        character(len=8) :: answer
        character(len=24) :: loadIden, loadField
        aster_logical :: paraIsTime, paraIsVite, paraIsAcce
        aster_logical :: lSuivLoad, lPiloLoad, loadExist
        aster_logical :: neumExcl(LOAD_NEUM_NBTYPE)
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

        neumExcl = ASTER_FALSE
        if (present(neumExcl_)) then
            neumExcl = neumExcl_
        end if

        do iTypeNeum = 1, LOAD_NEUM_NBTYPE
            loadIden = 'None'

! --------- Detect this kind of load
            call isMecaLoadExist(iTypeNeum, loadPreObjectZ, &
                                 loadExist, loadField)

! --------- Manage exclusion
            if (loadExist) then
                if (neumExcl(iTypeNeum)) then
                    call utmess('A', 'CHARGES9_63', sk=mecaLoadKeyword(iTypeNeum))
                    loadExist = ASTER_FALSE
                end if
            end if

            if (loadExist) then
! ------------- Parameter dependence
                paraIsTime = ASTER_FALSE
                paraIsVite = ASTER_FALSE
                paraIsAcce = ASTER_FALSE

                if (iTypeNeum .ne. LOAD_NEUM_VECT_ASSE) then
                    answer = "NON"
                    call dismoi('PARA_INST', loadField, 'CARTE', repk=answer)
                    paraIsTime = answer .eq. 'OUI'
                    answer = "NON"
                    call dismoi('PARA_VITE', loadField, 'CARTE', repk=answer)
                    paraIsVite = answer .eq. 'OUI'
                    answer = "NON"
                    call dismoi('PARA_ACCE', loadField, 'CARTE', repk=answer)
                    paraIsAcce = answer .eq. 'OUI'
                end if

! ------------- Get identifier of load
                if (iTypeNeum .eq. LOAD_NEUM_PRE_SIGM) then
                    loadIden = "NEUM_SIGM_INT"
                else
                    call getLoadIdenNeum(loadApply, loadIsFunc, paraIsTime, &
                                         loadIden)
                end if

! ------------- Add load to list
                ASSERT(loadIden .ne. 'None')
                nbLoadIden = nbLoadIden+1
                ASSERT(nbLoadIden .lt. LOAD_NBIDEN_MAXI)
                listLoadIden(nbLoadIden) = loadIden

! ------------- Automatic checks
                call getApplyTypeForce(iTypeNeum, lSuivLoad, lPiloLoad)
                if (loadApply == 'FIXE_PILO' .or. loadApply == 'SUIV_PILO') then
                    if (.not. lPiloLoad) then
                        call utmess('F', 'CHARGES9_2')
                    end if
                end if
                if (loadApply == 'SUIV' .or. loadApply == 'SUIV_PILO') then
                    if (.not. lSuivLoad) then
                        call utmess('F', 'CHARGES9_3')
                    end if
                end if

! ------------- Some loads have to be undead loads
                if (.not. (loadApply .eq. 'SUIV')) then
                    if (iTypeNeum .eq. LOAD_NEUM_THM_ECHA .or. &
                        iTypeNeum .eq. LOAD_NEUM_THM_ECHAH) then
                        call utmess('F', 'CHARGES9_4')
                    end if
                    if (paraIsVite .or. paraIsAcce) then
                        call utmess('F', 'CHARGES9_5')
                    end if
                end if

! ------------- Some loads only for dynamic
                if (paraIsVite .or. paraIsAcce) then
                    if (staticOperator) then
                        call utmess('F', 'CHARGES9_6')
                    end if
                end if
            end if
        end do
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadKine
!
! Get identifier for kinematic boundary conditions (AFFE_CHAR_CINE)
!
! In  loadApply         : how to apply load
! In  loadCommand       : command defining load (AFFE_CHAR_*)
! In  loadIsFunc        : flag if load is a function
! IO  nbLoadIden        : number of identifier of loads in listLoadIden
! IO  listLoadIden      : list of loads's identifier
! In  kineExcl          : flag to exclude AFFE_CHAR_CINE
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadKine(loadApply, loadCommand, loadIsFunc, &
                           nbLoadIden, listLoadIden, &
                           kineExcl_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: loadApply
        character(len=8), intent(in) :: loadCommand
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(inout) :: nbLoadIden
        character(len=24), intent(inout) :: listLoadIden(LOAD_NBIDEN_MAXI)
        aster_logical, optional, intent(in) :: kineExcl_
! ----- Locals
        character(len=24) :: loadIden
        aster_logical :: kineExcl
!   ------------------------------------------------------------------------------------------------
!
        kineExcl = ASTER_FALSE
        if (present(kineExcl_)) then
            kineExcl = kineExcl_
        end if
        loadIden = "None"
        if (loadCommand(1:5) .eq. 'CIME_' .or. &
            loadCommand(1:5) .eq. 'CITH_' .or. &
            loadCommand(1:5) .eq. 'CIAC_') then
            if (loadApply .eq. 'SUIV') then
                call utmess('F', 'CHARGES9_3')
            else if (loadApply .eq. 'FIXE_PILO') then
                call utmess('F', 'CHARGES9_7')
            else if (loadApply .eq. 'DIDI') then
                if (loadIsFunc) then
                    if (loadCommand(5:7) .eq. "_FT") then
                        loadIden = 'CINE_FT_DIDI'
                    else
                        loadIden = 'CINE_FO_DIDI'
                    end if
                else
                    loadIden = 'CINE_CSTE_DIDI'
                end if
            else if (loadApply .eq. 'FIXE_CSTE') then
                if (loadIsFunc) then
                    if (loadCommand(5:7) .eq. "_FT") then
                        loadIden = 'CINE_FT'
                    else
                        loadIden = 'CINE_FO'
                    end if
                else
                    loadIden = 'CINE_CSTE'
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if (loadIden .ne. "None") then
            if (kineExcl) then
                call utmess('A', 'CHARGES9_61')
            else
                nbLoadIden = nbLoadIden+1
                ASSERT(nbLoadIden .lt. LOAD_NBIDEN_MAXI)
                listLoadIden(nbLoadIden) = loadIden
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadDiri
!
! Get identifier for kinematic boundary conditions with Lagrange multiplers (AFFE_CHAR_MECA)
!
! In  loadPreObject     : base JEVEUX name for object
! In  loadApply         : how to apply load
! In  loadIsFunc        : flag if load is a function
! IO  nbLoadIden        : number of identifier of loads in listLoadIden
! IO  listLoadIden      : list of loads's identifier
! In  diriExcl          : flag to exclude AFFE_CHAR_MECA/DDL_IMPO
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadDiri(loadPreObjectZ, loadApply, loadIsFunc, &
                           nbLoadIden, listLoadIden, &
                           diriExcl_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=16), intent(in) :: loadApply
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(inout) :: nbLoadIden
        character(len=24), intent(inout) :: listLoadIden(LOAD_NBIDEN_MAXI)
        aster_logical, optional, intent(in) :: diriExcl_
! ----- Locals
        character(len=24) :: loadField, loadIden
        integer(kind=8) :: iret
        character(len=8) :: answer
        aster_logical :: paraIsTime, diriExcl
!   ------------------------------------------------------------------------------------------------
!
        diriExcl = ASTER_FALSE
        if (present(diriExcl_)) then
            diriExcl = diriExcl_
        end if
        loadIden = 'None'
        loadField = loadPreObjectZ(1:13)//'.CIMPO.DESC'
        call exisd('CHAMP_GD', loadField, iret)
        if (iret .eq. 1) then
            call dismoi('PARA_INST', loadField, 'CARTE', repk=answer)
            paraIsTime = answer .eq. "OUI"
            if (loadApply .eq. 'SUIV') then
                loadIden = 'DIRI_SUIV'
            else if (loadApply .eq. 'FIXE_PILO') then
                if (loadIsFunc) then
                    loadIden = 'DIRI_PILO_F'
                    if (paraIsTime) then
                        call utmess('F', 'CHARGES9_1')
                    end if
                else
                    loadIden = 'DIRI_PILO'
                end if
            else if (loadApply .eq. 'DIDI') then
                if (loadIsFunc) then
                    loadIden = 'DIRI_FO_DIDI'
                    if (paraIsTime) then
                        loadIden = 'DIRI_FT_DIDI'
                    end if
                else
                    loadIden = 'DIRI_CSTE_DIDI'
                end if
            else if (loadApply .eq. 'FIXE_CSTE') then
                if (loadIsFunc) then
                    loadIden = 'DIRI_FO'
                    if (paraIsTime) then
                        loadIden = 'DIRI_FT'
                    end if
                else
                    loadIden = 'DIRI_CSTE'
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if (loadIden .ne. "None") then
            if (diriExcl) then
                call utmess('A', 'CHARGES9_62')
            else
                nbLoadIden = nbLoadIden+1
                ASSERT(nbLoadIden .lt. LOAD_NBIDEN_MAXI)
                listLoadIden(nbLoadIden) = loadIden
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadGround
!
! Get identifier for EXCIT_SOL load
!
! In  staticOperator    : flag for static operator
! In  loadPreObject     : base JEVEUX name for object
! In  loadApply         : how to apply load
! In  loadIsFunc        : flag if load is a function
! IO  nbLoadIden        : number of identifier of loads in listLoadIden
! IO  listLoadIden      : list of loads's identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadGround(staticOperator, loadPreObjectZ, loadApply, loadIsFunc, &
                             nbLoadIden, listLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: staticOperator
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=16), intent(in) :: loadApply
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(inout) :: nbLoadIden
        character(len=24), intent(inout) :: listLoadIden(LOAD_NBIDEN_MAXI)
! ----- Locals
        character(len=24) :: loadIden, loadField
        integer(kind=8) :: iret
!   ------------------------------------------------------------------------------------------------
!
        loadIden = 'None'
        loadField = loadPreObjectZ(1:13)//'.VEISS'
        call jeexin(loadField, iret)
        if (iret .ne. 0) then
            if (loadApply .eq. 'SUIV') then
                call utmess('F', 'CHARGES9_11')
            else if (loadApply .eq. 'DIDI') then
                call utmess('F', 'CHARGES9_12')
            else if (loadApply .eq. 'FIXE_PILO') then
                call utmess('F', 'CHARGES9_13')
            else if (loadApply .eq. 'FIXE_CSTE') then
                ASSERT(.not. loadIsFunc)
                loadIden = 'EXCIT_SOL'
            else if (loadApply .eq. 'SUIV_PILO') then
                call utmess('F', 'CHARGES9_14')
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if (loadIden .ne. "None") then
            if (staticOperator) then
                call utmess('F', 'CHARGES9_15')
            end if
            nbLoadIden = nbLoadIden+1
            ASSERT(nbLoadIden .lt. LOAD_NBIDEN_MAXI)
            listLoadIden(nbLoadIden) = loadIden
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadEvolChar
!
! Get identifier for EVOL_CHAR load
!
! In  loadPreObject     : base JEVEUX name for object
! In  loadApply         : how to apply load
! In  loadIsFunc        : flag if load is a function
! IO  nbLoadIden        : number of identifier of loads in listLoadIden
! IO  listLoadIden      : list of loads's identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadEvolChar(loadPreObjectZ, loadApply, loadIsFunc, &
                               nbLoadIden, listLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=16), intent(in) :: loadApply
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(inout) :: nbLoadIden
        character(len=24), intent(inout) :: listLoadIden(LOAD_NBIDEN_MAXI)
! ----- Locals
        character(len=24) :: loadIden, loadField
        integer(kind=8) :: iret
!   ------------------------------------------------------------------------------------------------
!
        loadIden = 'None'
        loadField = loadPreObjectZ(1:13)//'.EVOL.CHAR '
        call jeexin(loadField, iret)

        if (iret .ne. 0) then
            if (loadApply .eq. 'SUIV') then
                loadIden = 'NEUM_SUIV'
            else if (loadApply .eq. 'DIDI') then
                call utmess('F', 'CHARGES9_9')
            else if (loadApply .eq. 'FIXE_PILO') then
                call utmess('F', 'CHARGES9_21')
            else if (loadApply .eq. 'FIXE_CSTE') then
                ASSERT(.not. loadIsFunc)
                loadIden = 'NEUM_CSTE'
            else if (loadApply .eq. 'SUIV_PILO') then
                call utmess('F', 'CHARGES9_22')
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if (loadIden .ne. "None") then
            nbLoadIden = nbLoadIden+1
            ASSERT(nbLoadIden .lt. LOAD_NBIDEN_MAXI)
            listLoadIden(nbLoadIden) = loadIden
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadPlaneWave
!
! Get identifier for ONDE_PLANE load
!
! In  loadPreObject     : base JEVEUX name for object
! In  loadApply         : how to apply load
! IO  nbLoadIden        : number of identifier of loads in listLoadIden
! IO  listLoadIden      : list of loads's identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadPlaneWave(loadPreObjectZ, loadApply, &
                                nbLoadIden, listLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=16), intent(in) :: loadApply
        integer(kind=8), intent(inout) :: nbLoadIden
        character(len=24), intent(inout) :: listLoadIden(LOAD_NBIDEN_MAXI)
! ----- Locals
        character(len=24) :: loadIden
        aster_logical :: loadExist
!   ------------------------------------------------------------------------------------------------
!
        loadIden = 'None'
        call isMecaLoadExist(LOAD_NEUM_PWAVE, loadPreObjectZ, &
                             loadExist)
        if (loadExist) then
            if (loadApply .eq. 'SUIV') then
                call utmess('F', 'CHARGES9_32')
            else if (loadApply .eq. 'DIDI') then
                call utmess('F', 'CHARGES9_9')
            else if (loadApply .eq. 'FIXE_PILO') then
                call utmess('F', 'CHARGES9_31')
            else if (loadApply .eq. 'FIXE_CSTE') then
                loadIden = 'NEUM_ONDE'
            else if (loadApply .eq. 'SUIV_PILO') then
                call utmess('F', 'CHARGES9_33')
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if (loadIden .ne. "None") then
            nbLoadIden = nbLoadIden+1
            ASSERT(nbLoadIden .lt. LOAD_NBIDEN_MAXI)
            listLoadIden(nbLoadIden) = loadIden
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNbLoadsFromUser
!
! Get number of loads from User
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  loadKeyword       : keyword to read loads
! Out nbLoadList        : number of loads in list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine getNbLoadsFromUser(phenom, loadKeyword, nbLoadList)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        character(len=16), intent(in) :: loadKeyword
        integer(kind=8), intent(out) :: nbLoadList
! ----- Locals
        integer(kind=8) :: iret_cable_cine, iret_cable, nocc
        integer(kind=8) :: nbFactorKeyword, iKeyword
        character(len=8) :: loadName
        character(len=24) :: preStressName, kineName
        character(len=13) :: loadPreObject
!   ------------------------------------------------------------------------------------------------
!
        nbLoadList = 0

! ----- Count number of loads
        if (loadKeyword .eq. '  ') then
            call getvid(loadKeyword, 'CHARGE', iocc=0, nbret=nocc)
            nbLoadList = abs(nocc)
        elseif (loadKeyword .eq. 'EXCIT') then
            call getfac(loadKeyword, nbFactorKeyword)
            do iKeyword = 1, nbFactorKeyword
                call getvid(loadKeyword, 'CHARGE', iocc=iKeyword, &
                            scal=loadName, nbret=nocc)
                if (phenom .eq. 'MECA') then
                    loadPreObject = loadName(1:8)//'.CHME'
                    if (nocc .eq. 1) then
                        preStressName = loadPreObject(1:13)//'.SIGIN'
                        call exisd("CHAMP", preStressName, iret_cable)
                        if (iret_cable .eq. 0) then
! --------------------- All loads
                            nbLoadList = nbLoadList+1
                        else
! --------------------- For DEFI_CABLE_BP: count load only if kinematic
! --------------------- (because Neumann is not load but initial stress !)
                            kineName = loadPreObject(1:13)//'.CIMPO'
                            call exisd("CHAMP", kineName, iret_cable_cine)
                            if (iret_cable_cine .ne. 0) then
                                nbLoadList = nbLoadList+1
                            end if
                        end if
                    end if
                elseif (phenom .eq. 'THER') then
                    if (nocc .eq. 1) then
                        nbLoadList = nbLoadList+1
                    end if
                elseif (phenom .eq. 'ACOU') then
                    if (nocc .eq. 1) then
                        nbLoadList = nbLoadList+1
                    end if
                else
                    ASSERT(ASTER_FALSE)
                end if
            end do
        else
            nbLoadList = 0
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadKeyword
!
! Get factor keyword to read loads
!
! Out loadKeyword     : keyword to read loads
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadKeyword(loadKeyword)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(out) :: loadKeyword
!   ------------------------------------------------------------------------------------------------
!
        loadKeyword = "None"
        if (getexm('EXCIT', 'CHARGE') .eq. 1) then
            loadKeyword = "EXCIT"
        end if
        if (getexm(' ', 'CHARGE') .eq. 1) then
            loadKeyword = " "
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadFunc
!
! Get multiplier function for current load
!
! In  listLoadPrep      : parameters to construct list of loads
! In  jvBase            : JEVEUX jvBase where to create objects
! In  funcCste          : constant function
! In  loadKeyword       : factor keyword to read loads
! In  iKeyword          : index of factor keyword to read loads
! Out loadFunc          : multiplier function for current load
! Out hasMultFunc       : flag when a multiplier function has been used from user
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadFunc(listLoadPrep, &
                           jvBase, funcCste, &
                           loadKeyword, iKeyword, &
                           loadFunc, hasMultFunc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(ListLoad_Prep), intent(in) :: listLoadPrep
        character(len=1), intent(in) :: jvBase
        character(len=8), intent(in) :: funcCste
        character(len=16), intent(in) :: loadKeyword
        integer(kind=8), intent(in) :: iKeyword
        character(len=8), intent(out) :: loadFunc
        aster_logical, intent(out) :: hasMultFunc
!   ------------------------------------------------------------------------------------------------
!
        loadFunc = " "
        call getLoadFuncUser(listLoadPrep%funcIsCplx, jvBase, funcCste, &
                             loadKeyword, iKeyword, &
                             loadFunc)
        hasMultFunc = loadFunc(1:2) .ne. '&&'
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadFuncUser
!
! Get multiplier function for current load (from user)
!
! In  funcIsCplx        : get/create complex function
! In  jvBase            : JEVEUX jvBase where to create objects
! In  funcCste          : constant function
! In  loadKeyword       : factor keyword to read loads
! In  iKeyword          : index of factor keyword to read loads
! Out loadFunc          : multiplier function for current load
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadFuncUser(funcIsCplx, jvBase, funcCste, &
                               loadKeyword, iKeyword, &
                               loadFunc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: funcIsCplx
        character(len=1), intent(in) :: jvBase
        character(len=16), intent(in) :: loadKeyword
        integer(kind=8), intent(in) :: iKeyword
        character(len=8), intent(in) :: funcCste
        character(len=8), intent(out) :: loadFunc
! ----- Locals
        integer(kind=8) :: nbFuncCplx, nbFuncReal, nbCoefCplx, nbCoefReal, nbCoefAcce
        character(len=4) :: knum
        complex(kind=8) :: coefCplx
        character(len=19) :: funcCste19
        real(kind=8) :: coefReal, coefImag
!   ------------------------------------------------------------------------------------------------
!
        loadFunc = " "
        funcCste19 = funcCste
        if (funcIsCplx) then
! --------- Get function from user
            call getvid(loadKeyword, 'FONC_MULT_C', iocc=iKeyword, &
                        scal=loadFunc, nbret=nbFuncCplx)
            call getvid(loadKeyword, 'FONC_MULT', iocc=iKeyword, &
                        scal=loadFunc, nbret=nbFuncReal)

! --------- Generate function from user's coefficients
            if ((nbFuncCplx .eq. 0) .and. (nbFuncReal .eq. 0)) then
                call codent(iKeyword, 'D0', knum)
                loadFunc = '&&NC'//knum
                call getvc8(loadKeyword, 'COEF_MULT_C', iocc=iKeyword, &
                            scal=coefCplx, nbret=nbCoefCplx)
                if (nbCoefCplx .eq. 0) then
                    call getvr8(loadKeyword, 'COEF_MULT', iocc=iKeyword, &
                                scal=coefReal, nbret=nbCoefReal)
                    ASSERT(nbCoefReal .eq. 0)
                    call focste(loadFunc, 'TOUTRESU', coefReal, jvBase)
                else
                    coefReal = dble(coefCplx)
                    coefImag = dimag(coefCplx)
                    call focstc(loadFunc, 'TOUTRESU', coefReal, coefImag, jvBase)
                end if
            end if
        else
! --------- Detect real coefficient for function
            nbFuncReal = 0
            call getvid(loadKeyword, 'FONC_MULT', iocc=iKeyword, nbval=0, &
                        nbret=nbFuncReal)

! --------- Detect 'ACCE' coefficient for function
            nbCoefAcce = 0
            call getvid(loadKeyword, 'ACCE', iocc=iKeyword, nbval=0, &
                        nbret=nbCoefAcce)

! --------- Detect COEF_MULT
            nbCoefReal = 0
            if (nbFuncReal .eq. 0 .and. nbCoefAcce .eq. 0) then
                call getvr8(loadKeyword, 'COEF_MULT', iocc=iKeyword, nbval=0, &
                            nbret=nbCoefReal)
            end if

            if (nbFuncReal .eq. 0 .and. nbCoefAcce .eq. 0 .and. nbCoefReal .eq. 0) then
                call createUnitFunc(funcCste, "G", loadFunc)
            else
                if (nbFuncReal .ne. 0) then
                    call getvid(loadKeyword, 'FONC_MULT', iocc=iKeyword, scal=loadFunc)
                end if
                if (nbCoefAcce .ne. 0) then
                    call getvid(loadKeyword, 'ACCE', iocc=iKeyword, scal=loadFunc)
                end if
                if (nbCoefReal .ne. 0) then
                    call codent(iKeyword, 'D0', knum)
                    loadFunc = '&&NC'//knum
                    call getvr8(loadKeyword, 'COEF_MULT', iocc=iKeyword, &
                                scal=coefReal, nbret=nbCoefReal)
                    call focste(loadFunc, 'TOUTRESU', coefReal, jvBase)
                end if
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! createUnitFunc
!
! Create unitary function
!
! In  funcCste          : constant function
! In  jvBase            : JEVEUX jvBase where to create objects
! Out loadFunc          : multiplier function for current load
!
! --------------------------------------------------------------------------------------------------
    subroutine createUnitFunc(funcCste, jvBase, loadFunc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: funcCste
        character(len=1), intent(in) :: jvBase
        character(len=8), intent(out) :: loadFunc
! ----- Locals
        integer(kind=8) :: iret
        real(kind=8), parameter :: coefReal = 1.d0
        character(len=19) :: funcCste19
!   ------------------------------------------------------------------------------------------------
!
        loadFunc = " "
        funcCste19 = funcCste
        call jeexin(funcCste19//'.PROL', iret)
        if (iret .eq. 0) then
            call focste(funcCste, 'TOUTRESU', coefReal, jvBase)
        end if
        loadFunc = funcCste
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadName
!
! Get name of load
!
! In  loadKeyword       : keyword to read loads
! In  iKeyword          : index of keyword to read loads
! In  iLoadList         : index of load in list of loads
! In  nbLoadList        : number of loads in list of loads
! Ptr loadDble          : pointer to list to detect double loads
! Out loadName          : name of load
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadName(loadKeyword, iKeyword, &
                           iLoadList, nbLoadList, loadDble, &
                           loadName)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: loadKeyword
        integer(kind=8), intent(in) :: iKeyword, iLoadList, nbLoadList
        character(len=8), pointer :: loadDble(:)
        character(len=8), intent(out) :: loadName
! ----- Locals
        character(len=8), pointer :: listLoadUser(:) => null()
        integer(kind=8) :: iLoadDble
!   ------------------------------------------------------------------------------------------------
!
        loadName = " "
        if (loadKeyword .eq. ' ') then
            AS_ALLOCATE(vk8=listLoadUser, size=nbLoadList)
            call getvid(loadKeyword, 'CHARGE', nbval=nbLoadList, vect=listLoadUser)
            loadName = listLoadUser(iKeyword)
            AS_DEALLOCATE(vk8=listLoadUser)
        elseif (loadKeyword .eq. 'EXCIT') then
            call getvid(loadKeyword, 'CHARGE', iocc=iKeyword, scal=loadName)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Only one load in the list
        do iLoadDble = 1, nbLoadList
            if (loadName .eq. loadDble(iLoadDble)) then
                call utmess('F', 'CHARGES9_10')
            end if
        end do
        loadDble(iLoadList) = loadName
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadIndexKeyword
!
! Get index of keyword
!
! In  loadKeyword       : keyword to read loads
! IO  iKeyword          : index of keyword to read loads
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadIndexKeyword(loadKeyword, iKeyword)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: loadKeyword
        integer(kind=8), intent(inout) :: iKeyword
! ----- Locals
        integer(kind=8) :: nocc
!   ------------------------------------------------------------------------------------------------
!
        iKeyword = iKeyword+1
30      continue
        call getvid(loadKeyword, 'CHARGE', iocc=iKeyword, nbval=0, nbret=nocc)
        if (nocc .eq. 0) then
            iKeyword = iKeyword+1
            goto 30
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadApply
!
! Get how to apply load
!
! In  loadKeyword       : keyword to read loads
! In  iKeyword          : index of keyword to read loads
! In  iLoadList         : index of load in list of loads
! In  nbLoadList        : number of loads in list of loads
! Out loadApply         : how to apply load
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadApply(loadKeyword, iKeyword, &
                            loadApply)
!   ---------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: loadKeyword
        integer(kind=8), intent(in) :: iKeyword
        character(len=16), intent(out) :: loadApply
!   ------------------------------------------------------------------------------------------------
!
        loadApply = " "
        if (getexm(loadKeyword, 'TYPE_CHARGE') .eq. 1) then
            call getvtx(loadKeyword, 'TYPE_CHARGE', iocc=iKeyword, scal=loadApply)
        else
            loadApply = 'FIXE_CSTE'
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! creaListLoad
!
! Create list of loads
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  jvBase            : JEVEUX jvBase where to create objects
! In  nbLoad            : number of loads
! In  listLoad          : list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine creaListLoad(phenom, jvBase, nbLoad, listLoadZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(in) :: listLoadZ
        integer(kind=8), intent(in) :: nbLoad
! ----- Locals
        integer(kind=8) :: nbLoadEff
        character(len=24) :: loadNameJv, loadInfoJv, loadFuncJv
        integer(kind=8), pointer :: listLoadInfo(:) => null()
        character(len=24), pointer :: listLoadName(:) => null(), listLoadFunc(:) => null()
!   ------------------------------------------------------------------------------------------------
!

! ----- Datastructure access
        loadNameJv = listLoadZ(1:19)//'.LCHA'
        loadInfoJv = listLoadZ(1:19)//'.INFC'
        loadFuncJv = listLoadZ(1:19)//'.FCHA'

! ----- No loads datastructure
        if (nbLoad .eq. 0) then
            nbLoadEff = 1
        else
            nbLoadEff = nbLoad
        end if
        call detrsd('LISTE_CHARGES', listLoadZ)

        if (phenom .eq. 'MECA') then
            call wkvect(loadNameJv, jvBase//' V K24', nbLoadEff, vk24=listLoadName)
            call wkvect(loadInfoJv, jvBase//' V IS', 4*nbLoadEff+7, vi=listLoadInfo)
            call wkvect(loadFuncJv, jvBase//' V K24', nbLoadEff, vk24=listLoadFunc)
        elseif (phenom .eq. 'THER') then
            call wkvect(loadNameJv, jvBase//' V K24', nbLoadEff, vk24=listLoadName)
            call wkvect(loadInfoJv, jvBase//' V IS', 2*nbLoadEff+1, vi=listLoadInfo)
            call wkvect(loadFuncJv, jvBase//' V K24', nbLoadEff, vk24=listLoadFunc)
        elseif (phenom .eq. 'ACOU') then
            call wkvect(loadNameJv, jvBase//' V K24', nbLoadEff, vk24=listLoadName)
            call wkvect(loadInfoJv, jvBase//' V IS', 2*nbLoadEff+1, vi=listLoadInfo)
            call wkvect(loadFuncJv, jvBase//' V K24', nbLoadEff, vk24=listLoadFunc)
        else
            ASSERT(ASTER_FALSE)
        end if
        listLoadInfo(1) = nbLoad
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNbLoadsFromList
!
! Get number of loads from list of loads
!
! In  listLoad          : list of loads
! Out nbLoad            : number of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine getNbLoadsFromList(listLoadZ, nbLoad)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoadZ
        integer(kind=8), intent(out) :: nbLoad
! ----- Locals
        integer(kind=8) :: iret
        character(len=19) :: listLoad
        character(len=24) :: loadInfoJv
        integer(kind=8), pointer :: listLoadInfo(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        nbLoad = 0
        listLoad = listLoadZ
        if (listLoad .ne. " ") then
            loadInfoJv = listLoadZ(1:19)//'.INFC'
            call jeexin(loadInfoJv, iret)
            if (iret .ne. 0) then
                call jeveuo(loadInfoJv, 'L', vi=listLoadInfo)
                nbLoad = listLoadInfo(1)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! addLoadToList
!
! Add load to list of loads
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  listLoad          : list of loads
! In  jvBase            : JEVEUX jvBase where to create objects
! In  loadName          : name of load
! In  loadFunc          : name of function
! In  nbLoadIden        : number of identifier of loads in listLoadIden
! In  listLoadIden      : list of loads's identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine addLoadToList(phenom, listLoadZ, jvBase, &
                             loadNameZ, loadFuncZ, &
                             nbLoadIden, listLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        character(len=*), intent(in) :: listLoadZ
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(in) :: loadNameZ, loadFuncZ
        integer(kind=8), intent(in) :: nbLoadIden
        character(len=24), intent(in) :: listLoadIden(LOAD_NBIDEN_MAXI)
! ----- Locals
        integer(kind=8), parameter :: nbLoadToAdd = 1
        integer(kind=8) :: nbLoad, indxLoad
!   ------------------------------------------------------------------------------------------------
!

! ----- Resize list of loads
        call resizeListLoad(phenom, jvBase, nbLoadToAdd, listLoadZ)
        call getNbLoadsFromList(listLoadZ, nbLoad)
        indxLoad = nbLoad

! ----- Set load in list
        call setLoadToList(phenom, listLoadZ, &
                           indxLoad, loadNameZ, loadFuncZ, &
                           nbLoadIden, listLoadIden)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! resizeListLoad
!
! Resize load to list of loads
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  jvBase            : JEVEUX jvBase where to create objects
! In  nbLoadToAdd       : number of loads to add
! In  listLoad          : list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine resizeListLoad(phenom, jvBase, nbLoadToAdd, listLoadZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        character(len=1), intent(in) :: jvBase
        integer(kind=8), intent(in) :: nbLoadToAdd
        character(len=*), intent(in) :: listLoadZ
! ----- Locals
        integer(kind=8) :: nbLoadOld, nbLoadNew
        character(len=24), parameter :: listLoadNew = "&&LISTCHAR"
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(phenom .eq. 'MECA')
        ASSERT(nbLoadToAdd .gt. 0)

! ----- Create new list of loads
        call getNbLoadsFromList(listLoadZ, nbLoadOld)
        nbLoadNew = nbLoadOld+nbLoadToAdd
        call creaListLoad(phenom, jvBase, nbLoadNew, listLoadNew)

! ----- Copy
        call copyListLoad(phenom, listLoadZ, listLoadNew)
        call copisd('LISTE_CHARGES', jvBase, listLoadNew, listLoadZ)
        call detrsd('LISTE_CHARGES', listLoadNew)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNbLoadType
!
! Get number of loads with predefined type
!
! In  listLoad          : list of loads
! In  LoadIden          : identifier of loads to seek
! In  nbLoadIden        : number of loads for the identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine getNbLoadType(listLoadZ, loadIdenZ, nbLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoadZ, loadIdenZ
        integer(kind=8), intent(out) :: nbLoadIden
! ----- Locals
        integer(kind=8) :: nbLoad, iLoad
        character(len=24) :: loadIden, loadInfoJv
        integer(kind=8), pointer :: listLoadInfo(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        nbLoadIden = 0
        loadIden = loadIdenZ

! ----- Datastructure access
        loadInfoJv = listLoadZ(1:19)//'.INFC'
        call jeveuo(loadInfoJv, 'L', vi=listLoadInfo)

! ----- Datastructure informations
        call getNbLoadsFromList(listLoadZ, nbLoad)

        do iLoad = 1, nbLoad
            if (loadIden .eq. "PILO") then
                if (listLoadInfo(iLoad+1) .eq. 5 .or. &
                    listLoadInfo(iLoad+1) .eq. 6) then
                    nbLoadIden = nbLoadIden+1
                end if
                if (listLoadInfo(nbLoad+iLoad+1) .eq. 5 .or. &
                    listLoadInfo(nbLoad+iLoad+1) .eq. 8 .or. &
                    listLoadInfo(nbLoad+iLoad+1) .eq. 9 .or. &
                    listLoadInfo(nbLoad+iLoad+1) .eq. 11) then
                    nbLoadIden = nbLoadIden+1
                end if
            else if (loadIden .eq. "DIRI_SUIV") then
                if (listLoadInfo(iLoad+1) .eq. 4) then
                    nbLoadIden = nbLoadIden+1
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! listLoadDebug
!
! Advanced debug of list of loads datastructure
!
! In  listLoad          : list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine listLoadDebug(listLoadZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoadZ
! ----- Locals
        integer(kind=8) :: nbLoad, iLoad
        character(len=8) :: loadName, loadFunc
        integer(kind=8) :: loadNumeKine, loadNumeNeum
        character(len=24) :: listLoad
        character(len=24) :: loadInfoJv, loadNameJv, loadFuncJv
        integer(kind=8), pointer :: listLoadInfo(:) => null()
        character(len=24), pointer :: listLoadName(:) => null(), listLoadFunc(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        listLoad = listLoadZ
        WRITE (6, *) " ===================== "
        WRITE (6, *) " DEBUG - LIST OF LOADS "
        WRITE (6, *) " ===================== "

! ----- Datastructure access
        loadNameJv = listLoadZ(1:19)//'.LCHA'
        loadInfoJv = listLoadZ(1:19)//'.INFC'
        loadFuncJv = listLoadZ(1:19)//'.FCHA'
        call jeveuo(loadInfoJv, 'L', vi=listLoadInfo)
        call jeveuo(loadNameJv, 'L', vk24=listLoadName)
        call jeveuo(loadFuncJv, 'L', vk24=listLoadFunc)

! ----- Datastructure informations
        nbLoad = listLoadInfo(1)
        WRITE (6, *) " Number of loads: ", nbLoad

        do iLoad = 1, nbLoad
            loadName = listLoadName(iLoad) (1:8)
            loadFunc = listLoadFunc(iLoad) (1:8)
            loadNumeKine = listLoadInfo(1+iLoad)
            loadNumeNeum = listLoadInfo(1+iLoad+nbLoad)
            WRITE (6, *) " Load <", iLoad, ">"
            WRITE (6, *) "   Name of datastructure from AFFE_CHAR_* : <", loadName, ">"
            WRITE (6, *) "   Name of function : <", loadFunc, ">"
            WRITE (6, *) "   loadNumeKine : <", loadNumeKine, ">"
            WRITE (6, *) "   loadNumeNeum : <", loadNumeNeum, ">"
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadIdenFromList
!
! Get load identifiers from list of loads
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  listLoad          : list of loads
! In  iLoad             : index of load
! Out nbLoadIden        : number of identifier of loads in listLoadIden
! Out listLoadIden      : list of loads's identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadIdenFromList(phenom, listLoadZ, iLoad, &
                                   nbLoadIden, listLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        character(len=*), intent(in) :: listLoadZ
        integer(kind=8), intent(in) :: iLoad
        integer(kind=8), intent(out) :: nbLoadIden
        character(len=24), intent(out) :: listLoadIden(LOAD_NBIDEN_MAXI)
! ----- Locals
        integer(kind=8) :: nbLoad, iLoadIden
        character(len=24) :: listLoad, loadInfoJv
        integer(kind=8), pointer :: listLoadInfo(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        listLoad = listLoadZ
        nbLoadIden = 0
        listLoadIden = " "
        ASSERT(phenom .eq. 'MECA')

! ----- Datastructure access
        call getNbLoadsFromList(listLoad, nbLoad)
        loadInfoJv = listLoad(1:19)//'.INFC'
        call jeveuo(loadInfoJv, 'L', vi=listLoadInfo)

        iLoadIden = 0
        if (listLoadInfo(iLoad+1) .eq. 4) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'DIRI_SUIV'
        end if
        if (listLoadInfo(iLoad+1) .eq. -1) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'CINE_CSTE'
        end if
        if (listLoadInfo(iLoad+1) .eq. -2) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'CINE_FO'
        end if
        if (listLoadInfo(iLoad+1) .eq. -3) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'CINE_FT'
        end if
        if (listLoadInfo(iLoad+1) .eq. 5) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'DIRI_PILO'
        end if
        if (listLoadInfo(iLoad+1) .eq. 6) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'DIRI_PILO_F'
        end if
        if (listLoadInfo(iLoad+1) .eq. 1) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'DIRI_CSTE'
            if (listLoadInfo(3*nbLoad+2+iLoad+1) .eq. 1) then
                listLoadIden(iLoadIden) = 'DIRI_CSTE_DIDI'
            end if
        end if
        if (listLoadInfo(iLoad+1) .eq. 2) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'DIRI_FO'
            if (listLoadInfo(3*nbLoad+2+iLoad+1) .eq. 1) then
                listLoadIden(iLoadIden) = 'DIRI_FO_DIDI'
            end if
        end if
        if (listLoadInfo(iLoad+1) .eq. 3) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'DIRI_FT'
            if (listLoadInfo(3*nbLoad+2+iLoad+1) .eq. 1) then
                listLoadIden(iLoadIden) = 'DIRI_FT_DIDI'
            end if
        end if
        if (listLoadInfo(nbLoad+iLoad+1) .eq. 6) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'NEUM_ONDE'
        end if
        if ((listLoadInfo(nbLoad+iLoad+1) .eq. 55) .and. (listLoadInfo(4*nbLoad+5) .eq. 99)) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'NEUM_SIGM_INT'
        end if
        if (listLoadInfo(nbLoad+iLoad+1) .eq. 5) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'NEUM_PILO'
        end if
        if (listLoadInfo(nbLoad+iLoad+1) .eq. 4) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'NEUM_SUIV'
        end if
        if (listLoadInfo(nbLoad+iLoad+1) .eq. 2) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'NEUM_FO'
        end if
        if (listLoadInfo(nbLoad+iLoad+1) .eq. 3) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'NEUM_FT'
        end if
        if (listLoadInfo(nbLoad+iLoad+1) .eq. 1) then
            iLoadIden = iLoadIden+1
            ASSERT(iLoadIden .le. LOAD_NBIDEN_MAXI)
            listLoadIden(iLoadIden) = 'NEUM_CSTE'
        end if
        nbLoadIden = iLoadIden
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setLoadToList
!
! Set load to list of loads
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  listLoad          : list of loads
! In  jvBase            : JEVEUX jvBase where to create objects
! In  indxLoadInList    : index of load in list of loads
! In  loadName          : name of load
! In  loadFunc          : name of function
! In  nbLoadIden        : number of identifier of loads in listLoadIden
! In  listLoadIden      : list of loads's identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine setLoadToList(phenom, listLoadZ, &
                             indxLoadInList, loadNameZ, loadFuncZ, &
                             nbLoadIden, listLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        character(len=*), intent(in) :: listLoadZ
        integer(kind=8), intent(in) :: indxLoadInList
        character(len=*), intent(in) :: loadNameZ, loadFuncZ
        integer(kind=8), intent(in) :: nbLoadIden
        character(len=24), intent(in) :: listLoadIden(LOAD_NBIDEN_MAXI)
! ----- Locals
        integer(kind=8) :: nbLoad, iLoadIden
        character(len=24) :: loadIden
        character(len=24) :: loadNameJv, loadInfoJv, loadFuncJv
        integer(kind=8), pointer :: listLoadInfo(:) => null()
        character(len=24), pointer :: listLoadName(:) => null(), listLoadFunc(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Datastructure informations
        call getNbLoadsFromList(listLoadZ, nbLoad)
        ASSERT(indxLoadInList .gt. 0)
        ASSERT(indxLoadInList .le. nbLoad)
        loadNameJv = listLoadZ(1:19)//'.LCHA'
        loadInfoJv = listLoadZ(1:19)//'.INFC'
        loadFuncJv = listLoadZ(1:19)//'.FCHA'
        call jeveuo(loadNameJv, 'E', vk24=listLoadName)
        call jeveuo(loadInfoJv, 'E', vi=listLoadInfo)
        call jeveuo(loadFuncJv, 'E', vk24=listLoadFunc)

! ----- Set basic properties of load
        listLoadName(indxLoadInList) = loadNameZ
        listLoadFunc(indxLoadInList) = loadFuncZ

! ----- Set type of load
        if (phenom .eq. 'MECA') then
            do iLoadIden = 1, nbLoadIden
                loadIden = listLoadIden(iLoadIden)
                if (loadIden .eq. 'CINE_CSTE') then
                    listLoadInfo(indxLoadInList+1) = -1
                else if (loadIden .eq. 'CINE_FO') then
                    listLoadInfo(indxLoadInList+1) = -2
                else if (loadIden .eq. 'CINE_FT') then
                    listLoadInfo(indxLoadInList+1) = -3
                else if (loadIden .eq. 'CINE_CSTE_DIDI') then
                    listLoadInfo(indxLoadInList+1) = -1
                    listLoadInfo(3*nbLoad+2+indxLoadInList+1) = 1
                else if (loadIden .eq. 'CINE_FO_DIDI') then
                    listLoadInfo(indxLoadInList+1) = -2
                    listLoadInfo(3*nbLoad+2+indxLoadInList+1) = 1
                else if (loadIden .eq. 'CINE_FT_DIDI') then
                    listLoadInfo(indxLoadInList+1) = -3
                    listLoadInfo(3*nbLoad+2+indxLoadInList+1) = 1
! --------------------------------------------------------------------------------------------------
                else if (loadIden .eq. 'DIRI_CSTE') then
                    listLoadInfo(indxLoadInList+1) = 1
                else if (loadIden .eq. 'DIRI_CSTE_DIDI') then
                    listLoadInfo(indxLoadInList+1) = 1
                    listLoadInfo(3*nbLoad+2+indxLoadInList+1) = 1
                else if (loadIden .eq. 'DIRI_FO') then
                    listLoadInfo(indxLoadInList+1) = 2
                else if (loadIden .eq. 'DIRI_FO_DIDI') then
                    listLoadInfo(indxLoadInList+1) = 2
                    listLoadInfo(3*nbLoad+2+indxLoadInList+1) = 1
                else if (loadIden .eq. 'DIRI_FT') then
                    listLoadInfo(indxLoadInList+1) = 3
                else if (loadIden .eq. 'DIRI_FT_DIDI') then
                    listLoadInfo(indxLoadInList+1) = 3
                    listLoadInfo(3*nbLoad+2+indxLoadInList+1) = 1
                else if (loadIden .eq. 'DIRI_SUIV') then
                    listLoadInfo(indxLoadInList+1) = 4
                else if (loadIden .eq. 'DIRI_PILO ') then
                    listLoadInfo(indxLoadInList+1) = 5
                else if (loadIden .eq. 'DIRI_PILO_F') then
                    listLoadInfo(indxLoadInList+1) = 6
                else if (loadIden .eq. 'NEUM_CSTE') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 1
                else if (loadIden .eq. 'NEUM_FO') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 2
                else if (loadIden .eq. 'NEUM_FT') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 3
                else if (loadIden .eq. 'NEUM_SUIV') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 4
                else if (loadIden .eq. 'NEUM_PILO') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 5
                else if (loadIden .eq. 'NEUM_ONDE') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 6
                else if (loadIden .eq. 'NEUM_PILO_F') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 8
                else if (loadIden .eq. 'NEUM_SUIP') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 9
                else if (loadIden .eq. 'NEUM_SUIP_F') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 11
                else if (loadIden .eq. 'EXCIT_SOL') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 20
                else if (loadIden .eq. 'NEUM_SIGM_INT') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 55
                    listLoadInfo(4*nbLoad+5) = 99
                else
                    write (6, *) 'LISCAD: ', loadIden
                    ASSERT(ASTER_FALSE)
                end if
            end do
        elseif (phenom .eq. 'THER') then
            do iLoadIden = 1, nbLoadIden
                loadIden = listLoadIden(iLoadIden)
                if (loadIden .eq. 'CINE_CSTE') then
                    listLoadInfo(indxLoadInList+1) = -1
                else if (loadIden .eq. 'CINE_FO') then
                    listLoadInfo(indxLoadInList+1) = -2
                else if (loadIden .eq. 'CINE_FT') then
                    listLoadInfo(indxLoadInList+1) = -3
                else if (loadIden(1:9) .eq. 'DIRI_CSTE') then
                    listLoadInfo(indxLoadInList+1) = 1
                else if (loadIden(1:9) .eq. 'DIRI_FO') then
                    listLoadInfo(indxLoadInList+1) = 2
                else if (loadIden(1:9) .eq. 'DIRI_FT') then
                    listLoadInfo(indxLoadInList+1) = 3
                else if (loadIden .eq. 'NEUM_CSTE') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 1
                else if (loadIden .eq. 'NEUM_FO') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 2
                else if (loadIden .eq. 'NEUM_FT') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 3
                else
                    write (6, *) 'LISCAD: ', loadIden
                    ASSERT(ASTER_FALSE)
                end if
            end do
        elseif (phenom .eq. 'ACOU') then
            do iLoadIden = 1, nbLoadIden
                loadIden = listLoadIden(iLoadIden)
                if (loadIden .eq. 'CINE_CSTE') then
                    listLoadInfo(indxLoadInList+1) = -1
                else if (loadIden .eq. 'CINE_FO') then
                    listLoadInfo(indxLoadInList+1) = -2
                else if (loadIden .eq. 'CINE_FT') then
                    listLoadInfo(indxLoadInList+1) = -3
                else if (loadIden(1:9) .eq. 'DIRI_CSTE') then
                    listLoadInfo(indxLoadInList+1) = 1
                else if (loadIden .eq. 'NEUM_CSTE') then
                    listLoadInfo(nbLoad+indxLoadInList+1) = 1
                else
                    write (6, *) 'LISCAD: ', loadIden
                    ASSERT(ASTER_FALSE)
                end if
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! copyListLoad
!
! Copy list of loads
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  listLoadIn        : list of loads (input)
! In  listLoadOut       : list of loads (output)
!
! --------------------------------------------------------------------------------------------------
    subroutine copyListLoad(phenom, listLoadInZ, listLoadOutZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        character(len=*), intent(in) :: listLoadInZ, listLoadOutZ
! ----- Locals
        integer(kind=8) :: nbLoadIn, nbLoadOut, iLoadIn
        character(len=24) :: loadNameInJv, loadFuncInJv
        character(len=24), pointer :: listLoadNameIn(:) => null(), listLoadFuncIn(:) => null()
        character(len=8) :: loadName, loadFunc
        integer(kind=8) :: nbLoadIden
        character(len=24) :: listLoadIden(LOAD_NBIDEN_MAXI)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(phenom .eq. 'MECA')

! ----- Datastructures access
        loadNameInJv = listLoadInZ(1:19)//'.LCHA'
        loadFuncInJv = listLoadInZ(1:19)//'.FCHA'
        call jeveuo(loadNameInJv, 'E', vk24=listLoadNameIn)
        call jeveuo(loadFuncInJv, 'E', vk24=listLoadFuncIn)
        call getNbLoadsFromList(listLoadInZ, nbLoadIn)
        call getNbLoadsFromList(listLoadOutZ, nbLoadOut)
        ASSERT(nbLoadOut .ge. nbLoadIn)

! ----- Copy
        do iLoadIn = 1, nbLoadIn
            loadName = listLoadNameIn(iLoadIn) (1:8)
            loadFunc = listLoadFuncIn(iLoadIn) (1:8)
            call getLoadIdenFromList(phenom, listLoadInZ, iLoadIn, &
                                     nbLoadIden, listLoadIden)
            call setLoadToList(phenom, listLoadOutZ, &
                               iLoadIn, loadName, loadFunc, &
                               nbLoadIden, listLoadIden)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! nameListLoad
!
! Generate name of datastructure for list of load
!
! Out listLoad          : name of datastructure for list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine nameListLoad(listLoad)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=24), intent(out) :: listLoad
! ----- Locals
        character(len=24) :: noobj
!   ------------------------------------------------------------------------------------------------
!
        listLoad = " "
        noobj = '12345678'//'.1234'//'.EXCIT'//".LCHA"
        call gnomsd(' ', noobj, 10, 13)
        listLoad = noobj(1:19)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getTherNeum
!
! Get identifier for thermal Neumann loads
!
! In  linearOperator    : flag for linear operator
! In  implicitTheta     : flag for implicit theta scheme
! In  loadPreObject     : base JEVEUX name for object
! In  loadApply         : how to apply load
! In  loadIsFunc        : flag if load is a function
! IO  nbLoadIden        : number of identifier of loads in listLoadIden
! IO  listLoadIden      : list of loads's identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine getTherNeum(linearOperator, implicitTheta, &
                           loadPreObjectZ, &
                           loadApply, loadIsFunc, hasMultFunc, &
                           nbLoadIden, listLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: linearOperator, implicitTheta
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=16), intent(in) :: loadApply
        aster_logical, intent(in) :: loadIsFunc, hasMultFunc
        integer(kind=8), intent(inout) :: nbLoadIden
        character(len=24), intent(inout) :: listLoadIden(LOAD_NBIDEN_MAXI)
! ----- Locals
        integer(kind=8) :: iTypeNeut
        character(len=8) :: answer
        character(len=24) :: loadIden, loadField
        aster_logical :: paraIsTime, paraIsTemp, loadExist
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

        do iTypeNeut = 1, LOAD_NEUT_NBTYPE
            loadIden = 'None'

! --------- Detect this kind of load
            call isTherLoadExist(iTypeNeut, loadPreObjectZ, &
                                 loadExist, loadField)

            if (loadExist) then
! ------------- Parameter dependence
                paraIsTime = ASTER_FALSE
                paraIsTemp = ASTER_FALSE
                if (loadIsFunc) then
                    call dismoi('PARA_INST', loadField, 'CARTE', repk=answer)
                    paraIsTime = answer .eq. 'OUI'
                    call dismoi('PARA_TEMP', loadField, 'CARTE', repk=answer)
                    paraIsTemp = answer .eq. 'OUI'
                end if

! ------------- Get identifier of load
                call getLoadIdenNeum(loadApply, loadIsFunc, paraIsTime, &
                                     loadIden)

! ------------- Automatic checks
                if (hasMultFunc) then
                    ! Débranchement provisoire (issue34322)
                    if (therLoadVectF(iTypeNeut) .ne. "NoVector" .and. paraIsTemp) then
                        call utmess('F', 'CHARGES9_43')
                    end if
                    if (iTypeNeut .eq. LOAD_NEUT_ECHANGE) then
                        call utmess('F', 'CHARGES9_39')
                    end if
                    if (.not. linearOperator .and. (.not. implicitTheta)) then
                        call utmess('F', 'CHARGES9_40', sk=therLoadKeyword(iTypeNeut))
                    end if
                end if

                if (iTypeNeut .eq. LOAD_NEUT_EVOL_CHAR) then
                    if (.not. linearOperator .and. (.not. implicitTheta)) then
                        call utmess('F', 'CHARGES9_41')
                    else
                        ASSERT(.not. loadIsFunc)
                        loadIden = 'NEUM_CSTE'
                    end if
                end if
                if (iTypeNeut .eq. LOAD_NEUT_FLUX_NL .or. &
                    iTypeNeut .eq. LOAD_NEUT_SOUR_NL .or. &
                    iTypeNeut .eq. LOAD_NEUT_RAYO) then
                    if (linearOperator) then
                        call utmess('F', 'CHARGES9_42', sk=therLoadKeyword(iTypeNeut))
                    end if
                end if

! ------------- Add load to list
                ASSERT(loadIden .ne. 'None')
                nbLoadIden = nbLoadIden+1
                ASSERT(nbLoadIden .lt. LOAD_NBIDEN_MAXI)
                listLoadIden(nbLoadIden) = loadIden
            end if
        end do
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkConsistency
!
! Get consistency of two list of loads
!
! In  listLoad1         : name of datastructure for first list of loads
! In  listLoad2         : name of datastructure for second list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine checkConsistency(listLoad1Z, listload2Z, lConsistent)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoad1Z, listload2Z
        aster_logical, intent(out) :: lConsistent
! ----- Locals
        integer(kind=8) :: nbLoad1, nbLoad2, iLoad1, iLoad2
        character(len=8) :: loadFunc
        integer(kind=8), pointer :: listLoadInfo1(:) => null()
        character(len=24), pointer :: listLoadName1(:) => null(), listLoadFunc1(:) => null()
        integer(kind=8), pointer :: listLoadInfo2(:) => null()
        character(len=24), pointer :: listLoadName2(:) => null(), listLoadFunc2(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        lConsistent = ASTER_TRUE
        call getListLoadAccess(listLoad1Z, &
                               nbLoad1, listLoadName1, listLoadInfo1, listLoadFunc1)
        call getListLoadAccess(listLoad2Z, &
                               nbLoad2, listLoadName2, listLoadInfo2, listLoadFunc2)

        if ((nbLoad1 .ne. 0) .and. (nbLoad2 .ne. 0)) then
            if (nbLoad1 .ne. nbLoad2) then
                call utmess('I', 'CHARGES9_49', ni=2, vali=[nbLoad1, nbLoad2])
            end if
! --------- Same loads ?
            do iLoad1 = 1, nbLoad1
                do iLoad2 = 1, nbLoad2
                    if (listLoadName1(iLoad1) .eq. listLoadName2(iLoad2) (1:8)) then
                        goto 30
                    end if
                end do
                call utmess('I', 'CHARGES9_50')
                lConsistent = ASTER_FALSE
30              continue
            end do

! --------- Same functions ?
            do iLoad1 = 1, nbLoad1
                do iLoad2 = 1, nbLoad2
                    loadFunc = listLoadFunc2(iLoad2) (1:8)
                    if (loadFunc(1:2) .eq. '&&') then
                        loadFunc = ' '
                    end if
                    if (listLoadFunc1(iLoad1) .eq. loadFunc) then
                        goto 60
                    end if
                    if (loadFunc .eq. ' ') then
                        goto 60
                    end if
                end do
                call utmess('I', 'CHARGES9_51')
                lConsistent = ASTER_FALSE
60              continue
            end do

! --------- Same (loads,functions) ?
            do iLoad1 = 1, nbLoad1
                do iLoad2 = 1, nbLoad2
                    if (listLoadName1(iLoad1) .eq. listLoadName2(iLoad2) (1:8)) then
                        if (listLoadFunc1(iLoad1) .eq. listLoadFunc2(iLoad2)) then
                            goto 95
                        end if
                        call utmess('I', 'CHARGES9_52')
                        lConsistent = ASTER_FALSE
                    end if
                end do
95              continue
            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListLoadAccess
!
! Get access to list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine getListLoadAccess(listLoadZ, &
                                 nbLoad, listLoadName, listLoadInfo, listLoadFunc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoadZ
        integer(kind=8), intent(out) :: nbLoad
        integer(kind=8), pointer :: listLoadInfo(:)
        character(len=24), pointer :: listLoadName(:), listLoadFunc(:)
! ----- Locals
        character(len=24) :: loadNameJv, loadInfoJv, loadFuncJv
!   ------------------------------------------------------------------------------------------------
!
        call getNbLoadsFromList(listLoadZ, nbLoad)
        loadNameJv = listLoadZ(1:19)//'.LCHA'
        loadInfoJv = listLoadZ(1:19)//'.INFC'
        loadFuncJv = listLoadZ(1:19)//'.FCHA'
        call jeveuo(loadNameJv, "E", vk24=listLoadName)
        call jeveuo(loadInfoJv, "E", vi=listLoadInfo)
        call jeveuo(loadFuncJv, "E", vk24=listLoadFunc)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLoadParameters
!
! Get other load parameters
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  model             : model
! In  loadName          : name of load
! Out loadPreObject     : base JEVEUX name for object
! Out loadCommand       : command defining load (AFFE_CHAR_*)
! Out loadIsFunc        : flag if load is a function
!
! --------------------------------------------------------------------------------------------------
    subroutine getLoadParameters(phenom, modelZ, loadNameZ, &
                                 loadPreObjectZ, loadCommandZ, loadIsFunc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        character(len=*), intent(in) :: loadNameZ, modelZ
        character(len=*), intent(out) :: loadPreObjectZ, loadCommandZ
        aster_logical, intent(out) :: loadIsFunc
! ----- Locals
        character(len=16) :: loadPhenomenon
        character(len=8) :: loadModel
!   ------------------------------------------------------------------------------------------------
!
        loadPreObjectZ = " "
        loadCommandZ = " "
        loadIsFunc = ASTER_FALSE

        if (phenom .eq. "MECA") then
            loadPreObjectZ = loadNameZ(1:8)//'.CHME'
        elseif (phenom .eq. "THER") then
            loadPreObjectZ = loadNameZ(1:8)//'.CHTH'
        elseif (phenom .eq. "ACOU") then
            loadPreObjectZ = loadNameZ(1:8)//'.CHAC'
        else
            ASSERT(ASTER_FALSE)
        end if
        call dismoi('TYPE_CHARGE', loadNameZ, 'CHARGE', repk=loadCommandZ)
        loadIsFunc = loadCommandZ(5:6) .eq. '_F'
        call dismoi('PHENOMENE', loadNameZ, 'CHARGE', repk=loadPhenomenon)

! ----- Checks
        if (phenom .eq. "MECA") then
            if (loadPhenomenon .ne. 'MECANIQUE') then
                call utmess('F', 'CHARGES9_17')
            end if
        elseif (phenom .eq. "THER") then
            if (loadPhenomenon .ne. 'THERMIQUE') then
                call utmess('F', 'CHARGES9_37')
            end if
        elseif (phenom .eq. "ACOU") then
            if (loadPhenomenon .ne. 'ACOUSTIQUE') then
                call utmess('F', 'CHARGES9_38')
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
        if (modelZ .ne. ' ') then
            call dismoi('NOM_MODELE', loadNameZ, 'CHARGE', repk=loadModel)
            if (loadModel .ne. modelZ) then
                call utmess('F', 'CHARGES9_16')
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! addLoadMeca
!
! Get mechanical loads and add in list
!
! In  staticOperator    : flag for static operator-
! In  listLoad          : list of loads
! In  loadName          : name of load
! In  loadFunc          : name of function
! In  loadApply         : how to apply load
! In  loadCommand       : command defining load (AFFE_CHAR_*)
! In  loadPreObject     : base JEVEUX name for object
! In  loadIsFunc        : flag if load is a function
! In  indxLoadInList    : index of load in list of loads
! In  kineExcl          : flag to exclude AFFE_CHAR_CINE
! In  diriExcl          : flag to exclude AFFE_CHAR_MECA/DDL_IMPO
! In  neumExcl          : flag to exclude AFFE_CHAR_MECA (Neumann/by type)
!
! --------------------------------------------------------------------------------------------------
    subroutine addLoadMeca(staticOperator, listLoadZ, &
                           loadNameZ, loadFuncZ, &
                           loadApplyZ, loadCommandZ, loadPreObjectZ, &
                           loadIsFunc, &
                           indxLoadInList, &
                           kineExcl_, diriExcl_, neumExcl_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: staticOperator
        character(len=*), intent(in) :: listLoadZ
        character(len=*), intent(in) :: loadNameZ, loadFuncZ
        character(len=*), intent(in) :: loadApplyZ, loadCommandZ, loadPreObjectZ
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(inout) :: indxLoadInList
        aster_logical, optional, intent(in) :: kineExcl_, diriExcl_, neumExcl_(LOAD_NEUM_NBTYPE)
! ----- Locals
        character(len=4), parameter :: phenom = "MECA"
        character(len=24) :: loadIden, listLoadIden(LOAD_NBIDEN_MAXI)
        integer(kind=8) :: nbLoadIden
        aster_logical :: kineExcl, diriExcl, neumExcl(LOAD_NEUM_NBTYPE)
!   ------------------------------------------------------------------------------------------------
!
        nbLoadIden = 0
        listLoadIden = " "
        loadIden = "None"

! ----- For exclusion
        kineExcl = ASTER_FALSE
        diriExcl = ASTER_FALSE
        neumExcl = ASTER_FALSE
        if (present(kineExcl_)) then
            kineExcl = kineExcl_
        end if
        if (present(diriExcl_)) then
            diriExcl = diriExcl_
        end if
        if (present(neumExcl_)) then
            neumExcl = neumExcl_
        end if

! ----- Dirichlet loads (AFFE_CHAR_CINE)
        call getLoadKine(loadApplyZ, loadCommandZ, loadIsFunc, &
                         nbLoadIden, listLoadIden, &
                         kineExcl)

! ----- Dirichlet loads (AFFE_CHAR_MECA)
        call getLoadDiri(loadPreObjectZ, loadApplyZ, loadIsFunc, &
                         nbLoadIden, listLoadIden, &
                         diriExcl)

! ----- Standard Neuman loads
        call getMecaNeum(staticOperator, loadPreObjectZ, &
                         loadApplyZ, loadIsFunc, &
                         nbLoadIden, listLoadIden, &
                         neumExcl)

! ----- Non-standard load: EXCIT_SOL
        call getLoadGround(staticOperator, loadPreObjectZ, loadApplyZ, loadIsFunc, &
                           nbLoadIden, listLoadIden)

! ----- Non-standard load: EVOL_CHAR
        call getLoadEvolChar(loadPreObjectZ, loadApplyZ, loadIsFunc, &
                             nbLoadIden, listLoadIden)

! ----- Non-standard load: ONDE_PLANE
        call getLoadPlaneWave(loadPreObjectZ, loadApplyZ, &
                              nbLoadIden, listLoadIden)

! ----- Add new load(s) in list of loads
        if (nbLoadIden .gt. 0) then
            indxLoadInList = indxLoadInList+1
            call setLoadToList(phenom, listLoadZ, &
                               indxLoadInList, loadNameZ, loadFuncZ, &
                               nbLoadIden, listLoadIden)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! addLoadTher
!
! Get thermal loads and add in list
!
! In  linearOperator    : flag for linear operator
! In  implicitTheta     : flag for implicit theta scheme
! In  listLoad          : list of loads
! In  loadName          : name of load
! In  loadFunc          : name of function
! In  loadApply         : how to apply load
! In  loadCommand       : command defining load (AFFE_CHAR_*)
! In  loadPreObject     : base JEVEUX name for object
! In  loadIsFunc        : flag if load is a function
! In  hasMultFunc       : flag when a multiplier function has been used from user
! In  indxLoadInList    : index of load in list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine addLoadTher(linearOperator, implicitTheta, listLoadZ, &
                           loadNameZ, loadFuncZ, &
                           loadApplyZ, loadCommandZ, loadPreObjectZ, &
                           loadIsFunc, hasMultFunc, &
                           indxLoadInList)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: linearOperator, implicitTheta
        character(len=*), intent(in) :: listLoadZ
        character(len=*), intent(in) :: loadNameZ, loadFuncZ
        character(len=*), intent(in) :: loadApplyZ, loadCommandZ, loadPreObjectZ
        aster_logical, intent(in) :: loadIsFunc, hasMultFunc
        integer(kind=8), intent(inout) :: indxLoadInList
! ----- Locals
        character(len=4), parameter :: phenom = "THER"
        character(len=24) :: loadIden, listLoadIden(LOAD_NBIDEN_MAXI)
        integer(kind=8) :: nbLoadIden
!   ------------------------------------------------------------------------------------------------
!
        nbLoadIden = 0
        listLoadIden = " "
        loadIden = "None"
        ASSERT(loadApplyZ .eq. 'FIXE_CSTE')

! ----- Dirichlet loads (AFFE_CHAR_CINE)
        call getLoadKine(loadApplyZ, loadCommandZ, loadIsFunc, &
                         nbLoadIden, listLoadIden)

! ----- Dirichlet loads (AFFE_CHAR_THER)
        call getLoadDiri(loadPreObjectZ, loadApplyZ, loadIsFunc, &
                         nbLoadIden, listLoadIden)

! ----- Standard Neuman loads
        call getTherNeum(linearOperator, implicitTheta, &
                         loadPreObjectZ, &
                         loadApplyZ, loadIsFunc, hasMultFunc, &
                         nbLoadIden, listLoadIden)

! ----- Non-standard load: EVOL_CHAR
        call getLoadEvolChar(loadPreObjectZ, loadApplyZ, loadIsFunc, &
                             nbLoadIden, listLoadIden)

! ----- Add new load(s) in list of loads
        if (nbLoadIden .gt. 0) then
            indxLoadInList = indxLoadInList+1
            call setLoadToList(phenom, listLoadZ, &
                               indxLoadInList, loadNameZ, loadFuncZ, &
                               nbLoadIden, listLoadIden)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! creaListLoadFSIOne
!
! Create list of pseudo-thermic load (for FSI) - One load
!
! In  model             : model
! In  nbLoad            : number of loads (only one !)
! In  loadName          : name of load
! In  listLoad          : list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine creaListLoadFSIOne(modelZ, nbLoad, loadNameZ, listLoadZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbLoad
        character(len=*), intent(in) :: modelZ, loadNameZ, listLoadZ
! ----- Locals
        character(len=4), parameter :: phenom = "THER"
        character(len=16), parameter :: loadApply = 'FIXE_CSTE'
        character(len=8), parameter :: funcCste = '&&NTDOCH'
        aster_logical, parameter :: hasMultFunc = ASTER_FALSE
        integer(kind=8) :: indxLoadInList
        character(len=13) :: loadPreObject
        character(len=16) :: loadCommand
        character(len=8) :: loadFunc
        aster_logical :: loadIsFunc, linearOperator, implicitTheta
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nbLoad .le. 1)
        linearOperator = ASTER_TRUE
        implicitTheta = ASTER_TRUE

! ----- Create list of loads datastructure
        call creaListLoad(phenom, "V", nbLoad, listLoadZ)

! ----- Get other load parameters
        call getLoadParameters(phenom, modelZ, loadNameZ, &
                               loadPreObject, loadCommand, loadIsFunc)

! ----- Create constant function
        call createUnitFunc(funcCste, "V", loadFunc)

! ----- Add thermal load
        indxLoadInList = 0
        call addLoadTher(linearOperator, implicitTheta, listLoadZ, &
                         loadNameZ, loadFunc, &
                         loadApply, loadCommand, loadPreObject, &
                         loadIsFunc, hasMultFunc, &
                         indxLoadInList)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! addLoadAcou
!
! Get acoustic loads and add in list
!
! In  listLoad          : list of loads
! In  loadName          : name of load
! In  loadFunc          : name of function
! In  loadApply         : how to apply load
! In  loadCommand       : command defining load (AFFE_CHAR_*)
! In  loadPreObject     : base JEVEUX name for object
! In  loadIsFunc        : flag if load is a function
! In  indxLoadInList    : index of load in list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine addLoadAcou(listLoadZ, &
                           loadNameZ, loadFuncZ, &
                           loadApplyZ, loadCommandZ, loadPreObjectZ, &
                           loadIsFunc, &
                           indxLoadInList)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoadZ
        character(len=*), intent(in) :: loadNameZ, loadFuncZ
        character(len=*), intent(in) :: loadApplyZ, loadCommandZ, loadPreObjectZ
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(inout) :: indxLoadInList
! ----- Locals
        character(len=4), parameter :: phenom = "THER"
        character(len=24) :: loadIden, listLoadIden(LOAD_NBIDEN_MAXI)
        integer(kind=8) :: nbLoadIden
!   ------------------------------------------------------------------------------------------------
!
        nbLoadIden = 0
        listLoadIden = " "
        loadIden = "None"
        ASSERT(loadApplyZ .eq. 'FIXE_CSTE')

! ----- Dirichlet loads (AFFE_CHAR_CINE)
        call getLoadKine(loadApplyZ, loadCommandZ, loadIsFunc, &
                         nbLoadIden, listLoadIden)

! ----- Dirichlet loads (AFFE_CHAR_ACOU)
        call getLoadDiri(loadPreObjectZ, loadApplyZ, loadIsFunc, &
                         nbLoadIden, listLoadIden)

! ----- Standard Neuman loads
        call getAcouNeum(loadPreObjectZ, &
                         loadApplyZ, loadIsFunc, &
                         nbLoadIden, listLoadIden)

! ----- Add new load(s) in list of loads
        if (nbLoadIden .gt. 0) then
            indxLoadInList = indxLoadInList+1
            call setLoadToList(phenom, listLoadZ, &
                               indxLoadInList, loadNameZ, loadFuncZ, &
                               nbLoadIden, listLoadIden)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getAcouNeum
!
! Get identifier for acoustic Neumann loads
!
! In  loadPreObject     : base JEVEUX name for object
! In  loadApply         : how to apply load
! In  loadIsFunc        : flag if load is a function
! IO  nbLoadIden        : number of identifier of loads in listLoadIden
! IO  listLoadIden      : list of loads's identifier
!
! --------------------------------------------------------------------------------------------------
    subroutine getAcouNeum(loadPreObjectZ, &
                           loadApply, loadIsFunc, &
                           nbLoadIden, listLoadIden)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=16), intent(in) :: loadApply
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(inout) :: nbLoadIden
        character(len=24), intent(inout) :: listLoadIden(LOAD_NBIDEN_MAXI)
! ----- Locals
        integer(kind=8) :: iTypeNeua
        character(len=8) :: answer
        character(len=24) :: loadIden, loadField
        aster_logical :: paraIsTime, loadExist
!   ------------------------------------------------------------------------------------------------
!
        do iTypeNeua = 1, LOAD_NEUA_NBTYPE
            loadIden = 'None'

! --------- Detect this kind of load:
            call isAcouLoadExist(iTypeNeua, loadPreObjectZ, &
                                 loadExist, loadField)

            if (loadExist) then
! ------------- Parameter dependence
                answer = "NON"
                if (loadIsFunc) then
                    call dismoi('PARA_INST', loadField, 'CARTE', repk=answer)
                end if
                paraIsTime = answer .eq. 'OUI'

! ------------- Get identifier of load
                call getLoadIdenNeum(loadApply, loadIsFunc, paraIsTime, &
                                     loadIden)

! ------------- Add load to list
                ASSERT(loadIden .ne. 'None')
                nbLoadIden = nbLoadIden+1
                ASSERT(nbLoadIden .lt. LOAD_NBIDEN_MAXI)
                listLoadIden(nbLoadIden) = loadIden
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! creaListLoadFromList
!
! Create list of loads from list of loads
!
! In  phenom            : phenomenon (MECA/THER/ACOU)
! In  listLoadPrep      : parameters to construct list of loads
! In  listLoad          : name of datastructure for list of loads
! In  jvBase            : JEVEUX base where to create objects
! In  nbLoad            : number of loads
! Ptr loadList          : pointer to list of loads
! In  kineExcl          : flag to exclude AFFE_CHAR_CINE
! In  diriExcl          : flag to exclude AFFE_CHAR_MECA/DDL_IMPO
!
! --------------------------------------------------------------------------------------------------
    subroutine creaListLoadFromList(phenom, listLoadPrep, &
                                    listLoadZ, jvBase, &
                                    nbLoadList, loadList, &
                                    kineExcl, diriExcl)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: phenom
        type(ListLoad_Prep), intent(in) :: listLoadPrep
        character(len=*), intent(in) :: listLoadZ
        character(len=1), intent(in) :: jvBase
        integer(kind=8), intent(in) :: nbLoadList
        character(len=8), pointer :: loadList(:)
        aster_logical, intent(in) :: kineExcl, diriExcl
! ----- Locals
        integer(kind=8) :: iLoadList, indxLoadInList, iLoadDble
        character(len=24) :: listLoad
        character(len=8), parameter :: funcCste = '&&SSCHGE'
        character(len=16), parameter :: loadApply = 'FIXE_CSTE'
        character(len=16) ::  loadCommand
        character(len=8) :: loadName, loadFunc
        character(len=13) :: loadPreObject
        aster_logical :: loadIsFunc
        character(len=8), pointer :: loadDble(:) => null()
        aster_logical, parameter :: neumExcl(LOAD_NEUM_NBTYPE) = ASTER_FALSE
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()
        ASSERT(phenom .eq. 'MECA')

! ----- Initializations
        listLoad = listLoadZ
        indxLoadInList = 0

! ----- Create list of loads datastructure
        call creaListLoad(phenom, jvBase, nbLoadList, listLoad)

! ----- List of loads to avoid same loads
        if (nbLoadList .ne. 0) then
            AS_ALLOCATE(vk8=loadDble, size=nbLoadList)
        end if

! ----- Add loads
        do iLoadList = 1, nbLoadList

! --------- Get current load
            loadName = loadList(iLoadList)

! --------- Only one load in the list
            do iLoadDble = 1, nbLoadList
                if (loadName .eq. loadDble(iLoadDble)) then
                    call utmess('F', 'CHARGES9_10')
                end if
            end do

! --------- Create constant function
            call createUnitFunc(funcCste, jvBase, loadFunc)

! --------- Get how to apply load
! --------- Only FIXE_CSTE

! --------- Get other load parameters
            call getLoadParameters(phenom, listLoadPrep%model, loadName, &
                                   loadPreObject, loadCommand, loadIsFunc)

! --------- Add mechanical loads
            call addLoadMeca(listLoadPrep%staticOperator, listLoad, &
                             loadName, loadFunc, &
                             loadApply, loadCommand, loadPreObject, &
                             loadIsFunc, &
                             indxLoadInList, &
                             kineExcl, diriExcl, neumExcl)

        end do

! ----- Debug
        call listLoadDebug(listLoad)

! ----- Clean
        AS_DEALLOCATE(vk8=loadDble)
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! detectMecaNeumLoad
!
! Detect mechanical Neumann loads
!
! In  listLoad          : name of datastructure for list of loads
! In  indxNeumType      : index of the type
! Out nbLoadDetect      : number of loads detect
! Ptr listLoadDetect    : pointer to list of loads
!
! --------------------------------------------------------------------------------------------------
    subroutine detectMecaNeumLoad(listLoadZ, indxNeumType, &
                                  nbLoadDetect, listLoadDetect)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoadZ
        integer(kind=8), intent(in) :: indxNeumType
        integer(kind=8), intent(out) :: nbLoadDetect
        character(len=8), pointer :: listLoadDetect(:)
! ----- Locals
        character(len=24) :: listLoad
        integer(kind=8) :: iLoad, nbLoad
        character(len=8) :: loadName
        character(len=24) :: loadNameJv
        character(len=24), pointer :: listLoadName(:) => null()
        character(len=13) :: loadPreObject
        aster_logical :: loadExist
!   ------------------------------------------------------------------------------------------------
!

! ----- Initializations
        listLoad = listLoadZ
        nbLoadDetect = 0

! ----- Access to datastructure
        call getNbLoadsFromList(listLoad, nbLoad)
        loadNameJv = listLoad(1:19)//'.LCHA'
        call jeveuo(loadNameJv, 'L', vk24=listLoadName)

! ----- Output
        if (nbLoad .gt. 0) then
            AS_ALLOCATE(vk8=listLoadDetect, size=nbLoad)
        end if

! ----- List of loads
        nbLoadDetect = 0
        do iLoad = 1, nbLoad
! --------- Current load
            loadName = listLoadName(iLoad) (1:8)
            loadPreObject = loadName(1:8)//'.CHME'

! --------- Detect this load
            call isMecaLoadExist(indxNeumType, loadPreObject, &
                                 loadExist)

! --------- Add load
            if (loadExist) then
                nbLoadDetect = nbLoadDetect+1
                ASSERT(nbLoadDetect .le. nbLoad)
                listLoadDetect(nbLoadDetect) = loadName
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListLoadLigrel
!
! Get all LIGREL from list of loads
!
! In  listLoad          : name of datastructure for list of loads
! IO  nbLigr            : number of LIGREL in list
! Ptr listLigr          : list of LIGREL
!
! --------------------------------------------------------------------------------------------------
    subroutine getListLoadLigrel(listLoadZ, nbLigr, listLigr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoadZ
        integer(kind=8), intent(inout) :: nbLigr
        character(len=24), pointer :: listLigr(:)
! ----- Locals
        character(len=8) :: loadCommand, loadName
        character(len=24) :: listLoad, loadLigrel, loadNameJv
        integer(kind=8) :: iLoad, nbLoad, nbLigrNew, nbLigrLoad, iret
        character(len=24), pointer :: listLigrSave(:) => null(), listLoadName(:) => null()
        character(len=13) :: loadPreObject
!   ------------------------------------------------------------------------------------------------
!
        listLoad = listLoadZ

! ----- Get number of loads
        call getNbLoadsFromList(listLoad, nbLoad)

! ----- Count number of LIGREL in loads
        nbLigrLoad = 0
        if (nbLoad .gt. 0) then
            loadNameJv = listLoad(1:19)//'.LCHA'
            call jeveuo(loadNameJv, 'L', vk24=listLoadName)
            do iLoad = 1, nbLoad
                loadName = listLoadName(iLoad) (1:8)
                if (loadName .ne. " ") then
                    call dismoi('TYPE_CHARGE', loadName, 'CHARGE', repk=loadCommand)
                    loadPreObject = loadName(1:8)//'.CH'//loadCommand(1:2)
                    loadLigrel = loadPreObject(1:13)//'.LIGRE'
                    call jeexin(loadLigrel(1:19)//'.LIEL', iret)
                    if (iret .gt. 0) then
                        nbLigrLoad = nbLigrLoad+1
                    end if
                end if
            end do
        end if

        if (nbLigr .eq. 0) then
! --------- Create list
            AS_ALLOCATE(vk24=listLigr, size=nbLigrLoad)
        else
! --------- Increase length og list
            nbLigrNew = nbLigr+nbLigrLoad
            if (nbLigrNew .gt. nbLigr) then
                AS_ALLOCATE(vk24=listLigrSave, size=nbLigr)
                listLigrSave(1:nbLigr) = listLigr(1:nbLigr)
                AS_DEALLOCATE(vk24=listLigr)
                AS_ALLOCATE(vk24=listLigr, size=nbLigrNew)
                listLigr(1:nbLigr) = listLigrSave(1:nbLigr)
                AS_DEALLOCATE(vk24=listLigrSave)
            end if
        end if

! ----- Add LIGREL
        if (nbLigrLoad .gt. 0) then
            loadNameJv = listLoad(1:19)//'.LCHA'
            call jeveuo(loadNameJv, 'L', vk24=listLoadName)
            do iLoad = 1, nbLoad
                loadName = listLoadName(iLoad) (1:8)
                if (loadName .ne. " ") then
                    call dismoi('TYPE_CHARGE', loadName, 'CHARGE', repk=loadCommand)
                    loadPreObject = loadName(1:8)//'.CH'//loadCommand(1:2)
                    loadLigrel = loadPreObject(1:13)//'.LIGRE'
                    call jeexin(loadLigrel(1:19)//'.LIEL', iret)
                    if (iret .gt. 0) then
                        nbLigr = nbLigr+1
                        listLigr(nbLigr) = loadLigrel
                    end if
                end if
            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module listLoad_module
