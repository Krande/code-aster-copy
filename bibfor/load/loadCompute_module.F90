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
! ==================================================================================================
!
! Module for the management of computation of loads
!
! ==================================================================================================
!
module loadCompute_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: compEvolChar, getApplyTypeForce, getRHSOption
    private :: getFieldFromEvol, loadWind
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/barych.h"
#include "asterfort/calcul.h"
#include "asterfort/chpnua.h"
#include "asterfort/cnocre.h"
#include "asterfort/copisd.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/gcnco2.h"
#include "asterfort/gcncon.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/load_neum_comp.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/nuachp.h"
#include "asterfort/pronua.h"
#include "asterfort/reajre.h"
#include "asterfort/rsinch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtgpld.h"
#include "jeveux.h"
#include "LoadTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! compEvolChar
!
! Compute loads from EVOL_CHAR keyword
!
! In  model             : name of model
! In  caraElem          : name of elementary characteristics (field)
! In  time              : time
! In  iLoad             : index of load in loads datastructure
! In  loadName          : name of current load
! In  loadApply         : type of application for load
!                        'Dead' - Dead loads (not dependent on displacements)
!                        'Pilo' - Loads for continuation (not dependent on displacements)
!                        'Suiv' - Undead loads (dependent on displacements)
! In  ligrel_calc       : LIGREL to compute
! In  nb_in_prep        : number of input fields (before specific ones)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! IO  resuElem          : number of elementary results
! In  vectElem          : name of elementary vectors
! In  disp_prev         : displacement at beginning of current time
! In  disp_cumu_inst    : displacement increment from beginning of current time
! In  strx_prev         : fibers information at beginning of current time
! In  vite_curr         : speed at current time
!
! --------------------------------------------------------------------------------------------------
    subroutine compEvolChar(model, caraElem, time, jvBase, &
                            iLoad, loadName, loadApply, ligrel_calc, &
                            nb_in_prep, lpain, lchin, &
                            resuElem, vectElem, &
                            disp_prev_, disp_cumu_inst_, strx_prev_, vite_curr_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=24), intent(in) :: model, caraElem
        real(kind=8), intent(in) :: time
        character(len=1), intent(in) :: jvBase
        integer, intent(in) :: iLoad
        character(len=8), intent(in) :: loadName
        character(len=4), intent(in) :: loadApply
        character(len=24), intent(in)  :: ligrel_calc
        integer, intent(in) :: nb_in_prep
        character(len=*), intent(inout) :: lpain(LOAD_NEUM_NBMAXIN)
        character(len=*), intent(inout) :: lchin(LOAD_NEUM_NBMAXIN)
        character(len=19), intent(inout) :: resuElem
        character(len=19), intent(in) :: vectElem
        character(len=19), optional, intent(in) :: disp_prev_, disp_cumu_inst_
        character(len=19), optional, intent(in) :: strx_prev_, vite_curr_
! ----- Local
        character(len=19), parameter :: inputFieldName = '&&NMDEPR'
        character(len=1), parameter :: stop = 'S'
! ----- Loads in EVOl_CHAR: no function
        aster_logical, parameter :: loadIsFunc = ASTER_FALSE
        integer, parameter :: loadNume = 1
        integer :: ier, nbField, iexist
        character(len=8) :: evol_char, newnom
        character(len=16) :: type_sd, option
        character(len=19) :: mapName
        character(len=24) :: loadObjectJv
        character(len=8), pointer :: loadObject(:) => null()
        integer :: indxNeumType
        character(len=19) :: disp_prev, disp_cumu_inst, strx_prev, vite_curr
        aster_logical :: existBody2D, existBody3D, existSurf, existLine, existPressure
        aster_logical :: existNodalForce, existWind
        aster_logical :: lSuivLoad
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Initializations
        disp_prev = " "
        if (present(disp_prev_)) then
            disp_prev = disp_prev_
        end if
        disp_cumu_inst = " "
        if (present(disp_cumu_inst_)) then
            disp_cumu_inst = disp_cumu_inst_
        end if
        strx_prev = " "
        if (present(strx_prev_)) then
            strx_prev = strx_prev_
        end if
        vite_curr = " "
        if (present(vite_curr_)) then
            vite_curr = vite_curr_
        end if

! ----- Get object from AFFE_CHAR_MECA
        loadObjectJv = loadName//'.CHME.EVOL.CHAR'
        call jeexin(loadObjectJv, ier)
        if (ier .eq. 0) then
            goto 99
        end if
        call jeveuo(loadObjectJv, 'L', vk8=loadObject)
        evol_char = loadObject(1)

        mapName = " "
        indxNeumType = -1

! ----- Some checks
        call dismoi('NB_CHAMP_UTI', evol_char, 'RESULTAT', repi=nbField)
        ASSERT(nbField .gt. 0)
        call gettco(evol_char, type_sd)
        ASSERT(type_sd .eq. 'EVOL_CHAR')

! ----- Get body force (3D)
        call getFieldFromEvol(evol_char, time, &
                              "FVOL_3D", inputFieldName, &
                              existBody3D)
        if (existBody3D) then
            mapName = ".F3D3D"
            indxNeumType = 2
        end if

! ----- Get body force (2d)
        call getFieldFromEvol(evol_char, time, &
                              "FVOL_2D", inputFieldName, &
                              existBody2D)
        if (existBody2D) then
            mapName = ".F2D2D"
            indxNeumType = 5
        end if

        if (existBody2D .and. existBody3D) then
            call utmess('F', 'CHARGES8_13')
        end if

! ----- Compute body forces (CHAR_MECA_FR2D2D / CHAR_MECA_FR3D3D)
        if (existBody3D .or. existBody2D) then
            call getApplyTypeForce(indxNeumType, lSuivLoad)
            if (loadApply .eq. "Suiv" .and. .not. lSuivLoad) then
                call utmess('F', 'CHARGES8_14')
            end if
            call getRHSOption(indxNeumType, loadApply, loadIsFunc, option)
            call load_neum_comp(stop, iLoad, loadName, loadNume, loadApply, &
                                ligrel_calc, LOAD_NEUM_NBMAXIN, nb_in_prep, lpain, lchin, &
                                jvBase, resuElem, vectElem, iden_direct=mapName, &
                                name_inputz=inputFieldName)
        end if

! ----- Get surfacic force
        call getFieldFromEvol(evol_char, time, &
                              "FSUR_3D", inputFieldName, &
                              existSurf)
        if (existSurf) then
            mapName = ".F2D3D"
            indxNeumType = 3
        end if

! ----- Get lineic force
        call getFieldFromEvol(evol_char, time, &
                              "FSUR_2D", inputFieldName, &
                              existLine)
        if (existLine) then
            mapName = ".F1D2D"
            indxNeumType = 6
        end if

        if (existSurf .and. existLine) then
            call utmess('F', 'CHARGES8_13')
        end if

! ----- Compute surfacic forces (CHAR_MECA_FR2D3D / CHAR_MECA_FR1D2D)
        if (existSurf .or. existLine) then
            call getApplyTypeForce(indxNeumType, lSuivLoad)
            if (loadApply .eq. "Suiv" .and. .not. lSuivLoad) then
                call utmess('F', 'CHARGES8_14')
            end if
            call getRHSOption(indxNeumType, loadApply, loadIsFunc, option)
            call load_neum_comp(stop, iLoad, loadName, loadNume, loadApply, &
                                ligrel_calc, LOAD_NEUM_NBMAXIN, nb_in_prep, lpain, lchin, &
                                jvBase, resuElem, vectElem, iden_direct=mapName, &
                                name_inputz=inputFieldName)
        end if

! ----- Get pressure
        call getFieldFromEvol(evol_char, time, &
                              "PRES", inputFieldName, &
                              existPressure)
        if (existPressure) then
            mapName = ".PRESS"
            indxNeumType = 10
        end if

! ----- Compute pressure (CHAR_MECA_PRES_R)
        if (existPressure) then
            call getApplyTypeForce(indxNeumType, lSuivLoad)
            if (loadApply .eq. "Suiv" .and. .not. lSuivLoad) then
                call utmess('F', 'CHARGES8_14')
            end if
            call getRHSOption(indxNeumType, loadApply, loadIsFunc, option)
            call load_neum_comp(stop, iLoad, loadName, loadNume, loadApply, &
                                ligrel_calc, LOAD_NEUM_NBMAXIN, nb_in_prep, lpain, lchin, &
                                jvBase, resuElem, vectElem, iden_direct=mapName, &
                                name_inputz=inputFieldName)
        end if

! ----- Get nodal force (VECT_ASSE)
        call getFieldFromEvol(evol_char, time, &
                              "FORC_NODA", inputFieldName, &
                              existNodalForce)
        if (existNodalForce) then
            option = "Copy_Load"
            mapName = "None"
        end if

! ----- Compute nodal force (VECT_ASSE)
        if (existNodalForce) then
            if (loadApply .eq. "Suiv") then
                call utmess('F', 'CHARGES8_14')
            end if
            newnom = resuElem(10:16)
            call gcnco2(newnom)
            resuElem(10:16) = newnom(2:8)
            call corich('E', resuElem, ichin_=iLoad)
            call copisd('CHAMP_GD', jvBase, inputFieldName, resuElem)
            call exisd('CHAMP_GD', resuElem, iexist)
            ASSERT((iexist .gt. 0) .or. (stop .eq. 'C'))
            call reajre(vectElem, resuElem, jvBase)
        end if

! ----- Get lineic forces for wind (CHAR_MECA_SR1D1D)
        call getFieldFromEvol(evol_char, time, &
                              "VITE_VENT", inputFieldName, &
                              existWind)

! ----- Compute pressure (CHAR_MECA_SR1D1D)
        if (loadApply .eq. "Suiv" .and. existWind) then
            call loadWind(model, caraElem, &
                          disp_prev, disp_cumu_inst, strx_prev, vite_curr, &
                          iLoad, inputFieldName, ligrel_calc, &
                          resuElem, vectElem)
        end if
!
99      continue
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getFieldFromEvol
!
! Get field from EVOL_CHAR datastructure
!
! In  ndim              : space dimension (2 or 3)
! In  evol_char         : name of datastructure from AFFE_CHAR_MECA for EVOL_CHAR
! In  time              : time
! In  nameInEvol        : type of load in EVOL_CHAR
! In  inputFieldName    : name of datastructure for defining load
! Out exist             : flag if nameInEvol is in evol_char
!
! --------------------------------------------------------------------------------------------------
    subroutine getFieldFromEvol(evol_char, time, &
                                nameInEvolZ, inputFieldName, &
                                exist)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: evol_char
        real(kind=8), intent(in) :: time
        character(len=*), intent(in) :: nameInEvolZ
        character(len=19), intent(in) :: inputFieldName
        aster_logical, intent(out):: exist
! ----- Local
        integer :: ier
        character(len=16) :: nameInEvol
        real(kind=8), parameter :: prec = 1.0d-10
        character(len=8), parameter :: crit = 'ABSOLU'
!   ------------------------------------------------------------------------------------------------
!
        exist = ASTER_FALSE
        nameInEvol = nameInEvolZ
        call rsinch(evol_char, nameInEvol, 'INST', time, inputFieldName, &
                    'EXCLU', 'EXCLU', 0, 'V', prec, crit, ier)
        if (ier .le. 2) then
            exist = ASTER_TRUE
        else if (ier .eq. 11 .or. ier .eq. 12 .or. ier .eq. 20) then
            call utmess('F', 'CHARGES8_2', sr=time, sk=nameInEvol)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getApplyTypeForce
!
! Get status of application of this kind of force
!
! In  indxNeumType      : index of the type
! Out lSuivLoad         : Flag if load can been undead load type
! Out lPiloLoad         : Flag if load can been used for continuation methods
!
! --------------------------------------------------------------------------------------------------
    subroutine getApplyTypeForce(indxNeumType, lSuivLoad_, lPiloLoad_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: indxNeumType
        aster_logical, optional, intent(out) :: lSuivLoad_, lPiloLoad_
! ----- Local
        aster_logical :: lSuivLoad, lPiloLoad
        aster_logical, parameter :: l_suiv(LOAD_NEUM_NBTYPE) = (/ &
                                                               .false._1, .false._1, &
                                                               .false._1, .false._1, &
                                                               .false._1, .false._1, &
                                                               .true._1, .true._1, &
                                                               .true._1, .true._1, &
                                                               .false._1, .true._1, &
                                                               .false._1, .false._1, &
                                                               .false._1, .false._1, &
                                                               .false._1, .true._1, &
                                                               .true._1, .true._1/)
        aster_logical, parameter :: l_pilo(LOAD_NEUM_NBTYPE) = (/ &
                                                               .true._1, .true._1, &
                                                               .true._1, .true._1, &
                                                               .true._1, .true._1, &
                                                               .true._1, .true._1, &
                                                               .false._1, .true._1, &
                                                               .false._1, .true._1, &
                                                               .true._1, .false._1, &
                                                               .true._1, .true._1, &
                                                               .false._1, .true._1, &
                                                               .false._1, .false._1/)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1 .and. indxNeumType .le. LOAD_NEUM_NBTYPE)
        lSuivLoad = l_suiv(indxNeumType)
        lPiloLoad = l_pilo(indxNeumType)
        if (present(lSuivLoad_)) then
            lSuivLoad_ = lSuivLoad
        end if
        if (present(lPiloLoad_)) then
            lPiloLoad_ = lPiloLoad
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getRHSOption
!
! Get name of option for vector
!
! In  indxNeumType      : index of the type
! In  loadApply         : type of application for load
!                        'Dead' - Dead loads (not dependent on displacements)
!                        'Pilo' - Loads for continuation (not dependent on displacements)
!                        'Suiv' - Undead loads (dependent on displacements)
! In  loadIsFunc        : flag if load is a function
! Out option            : option for RHS
!
! --------------------------------------------------------------------------------------------------
    subroutine getRHSOption(indxNeumType, loadApply, loadIsFunc, option)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: indxNeumType
        character(len=4), intent(in) :: loadApply
        aster_logical, intent(in) :: loadIsFunc
        character(len=16), intent(out) :: option
! ----- Local
        character(len=16), parameter :: optionRHSFunc(LOAD_NEUM_NBTYPE) = (/ &
                                        'CHAR_MECA_FORC_F', 'CHAR_MECA_FF3D3D', &
                                        'CHAR_MECA_FF2D3D', 'CHAR_MECA_FF1D3D', &
                                        'CHAR_MECA_FF2D2D', 'CHAR_MECA_FF1D2D', &
                                        'CHAR_MECA_FF1D1D', 'CHAR_MECA_PESA_R', &
                                        'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRES_F', &
                                        'CHAR_MECA_FRELEC', 'CHAR_MECA_FFCO3D', &
                                        'CHAR_MECA_FFCO2D', 'CHAR_MECA_EPSI_F', &
                                        'CHAR_MECA_FLUX_F', 'Copy_Load       ', &
                                        'FORC_NODA       ', 'CHAR_MECA_EFON_F', &
                                        'CHAR_ECHA_THM_F ', 'CHAR_ECHA_HR_F  '/)
        character(len=16), parameter :: optionRHSReal(LOAD_NEUM_NBTYPE) = (/ &
                                        'CHAR_MECA_FORC_R', 'CHAR_MECA_FR3D3D', &
                                        'CHAR_MECA_FR2D3D', 'CHAR_MECA_FR1D3D', &
                                        'CHAR_MECA_FR2D2D', 'CHAR_MECA_FR1D2D', &
                                        'CHAR_MECA_FR1D1D', 'CHAR_MECA_PESA_R', &
                                        'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRES_R', &
                                        'CHAR_MECA_FRELEC', 'CHAR_MECA_FRCO3D', &
                                        'CHAR_MECA_FRCO2D', 'CHAR_MECA_EPSI_R', &
                                        'CHAR_MECA_FLUX_R', 'Copy_Load       ', &
                                        'FORC_NODA       ', 'CHAR_MECA_EFON_R', &
                                        'CHAR_ECHA_THM_R ', 'CHAR_ECHA_HR_R  '/)
        character(len=16), parameter :: optionSRHSFunc(LOAD_NEUM_NBTYPE) = (/ &
                                        'No_Load         ', 'No_Load         ', &
                                        'No_Load         ', 'No_Load         ', &
                                        'No_Load         ', 'No_Load         ', &
                                        'CHAR_MECA_SF1D1D', 'No_Load         ', &
                                        'No_Load         ', 'CHAR_MECA_PRSU_F', &
                                        'No_Load         ', 'CHAR_MECA_SFCO3D', &
                                        'No_Load         ', 'No_Load         ', &
                                        'No_Load         ', 'No_Load         ', &
                                        'No_Load         ', 'CHAR_MECA_EFSU_F', &
                                        'CHAR_ECHA_THM_F ', 'CHAR_ECHA_HR_F  '/)
        character(len=16), parameter :: optionSRHSReal(LOAD_NEUM_NBTYPE) = (/ &
                                        'No_Load         ', 'No_Load         ', &
                                        'No_Load         ', 'No_Load         ', &
                                        'No_Load         ', 'No_Load         ', &
                                        'CHAR_MECA_SR1D1D', 'CHAR_MECA_PESA_R', &
                                        'CHAR_MECA_ROTA_R', 'CHAR_MECA_PRSU_R', &
                                        'No_Load         ', 'CHAR_MECA_SRCO3D', &
                                        'No_Load         ', 'No_Load         ', &
                                        'No_Load         ', 'No_Load         ', &
                                        'No_Load         ', 'CHAR_MECA_EFSU_R', &
                                        'CHAR_ECHA_THM_R ', 'CHAR_ECHA_HR_R  '/)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1 .and. indxNeumType .le. LOAD_NEUM_NBTYPE)
        if (loadIsFunc) then
            if (loadApply .eq. "Suiv") then
                option = optionSRHSFunc(indxNeumType)
            else if (loadApply .eq. "Dead") then
                option = optionRHSFunc(indxNeumType)
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            if (loadApply .eq. "Suiv") then
                option = optionSRHSReal(indxNeumType)
            else if (loadApply .eq. "Dead") then
                option = optionRHSReal(indxNeumType)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! loadWind
!
! Compute load for wind
!
! In  model             : name of model
! In  caraElem          : name of elementary characteristics (field)
! In  disp_prev         : displacement at beginning of current time
! In  disp_cumu_inst    : displacement increment from beginning of current time
! In  strx_prev         : fibers information at beginning of current time
! In  vite_curr         : speed at current time
! In  iLoad             : index of load in loads datastructure
! In  inputFieldName    : name of datastructure for defining load
! In  ligrel_calc       : LIGREL to compute
! IO  resuElem       : number of elementary results
! In  vectElem       : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine loadWind(model, caraElem, &
                        disp_prev, disp_cumu_inst, strx_prev, vite_curr, &
                        iLoad, inputFieldName, ligrel_calc, &
                        resuElem, vectElem)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=24), intent(in) :: model, caraElem
        character(len=19), intent(in) :: disp_prev, disp_cumu_inst, strx_prev, vite_curr
        integer, intent(in) :: iLoad
        character(len=19), intent(in) :: inputFieldName
        character(len=24), intent(in) :: ligrel_calc
        character(len=19), intent(inout) :: resuElem
        character(len=19), intent(in) :: vectElem
! ----- Local
        integer, parameter :: nbCmp = 3
        character(len=8), parameter :: cmpName(3) = (/'DX', 'DY', 'DZ'/)
        character(len=1), parameter :: jvBaseTemporary = "V"
        character(len=16), parameter :: option = 'CHAR_MECA_SR1D1D'
        character(len=19), parameter :: field_no_refe = '&&MNVGME.RESU_PROJE'
        character(len=19), parameter :: nuage1 = '&&NUAGE1', nuage2 = '&&NUAGE2'
        character(len=19), parameter :: method = 'NUAGE_DEG_1'
        integer, parameter :: nbFieldIn = 8, nbFieldOut = 1
        character(len=8) :: lpain(nbFieldIn), lpaout(nbFieldOut)
        character(len=24) :: lchin(nbFieldIn)
        integer :: nbEqua, ndim, nbno, dime, ibid, ier
        character(len=8) :: mesh_1, mesh_2, mesh_defo, answer, newnom
        integer, pointer :: mesh1Dime(:) => null()
        character(len=19) :: nume_equa, field_no_refe1
        character(len=24) :: chgeom, chcara(18)
        character(len=24), pointer :: fieldRefe(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        mesh_defo = '.0000000'
        field_no_refe1 = '.0000000'
        newnom = '.0000000'

! ----- Mesh information
        call jelira(inputFieldName//'.VALE', 'LONMAX', ival=nbequa)
        call dismoi('NOM_MAILLA', inputFieldName, 'CHAMP', repk=mesh_1)
        call dismoi('Z_CST', mesh_1, 'MAILLAGE', repk=answer)
        ndim = 3
        if (answer .eq. 'OUI') then
            ndim = 2
        end if

! ----- Check
        call jeveuo(mesh_1//'.DIME', 'E', vi=mesh1Dime)
        nbno = mesh1Dime(1)
        dime = mesh1Dime(6)
        if (nbno*dime .ne. nbequa) then
            call utmess('F', 'CHARGES8_10')
        end if

! ----- Mesh deformation
        call gcncon('.', mesh_defo)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh_2)
        call copisd('MAILLAGE', jvBaseTemporary, mesh_2, mesh_defo)
        call vtgpld('CUMU', 1.d0, mesh_2//'.COORDO', disp_prev, jvBaseTemporary, &
                    mesh_defo//'.COORDO1')
        call vtgpld('CUMU', 1.d0, mesh_defo//'.COORDO1', disp_cumu_inst, jvBaseTemporary, &
                    mesh_defo//'.COORDO')
        call detrsd('CHAMP_GD', mesh_defo//'.COORDO1')

! ----- Create reference field
        ibid = 0
        call cnocre(mesh_defo, 'DEPL_R', 0, [ibid], nbCmp, &
                    cmpName, [ibid], jvBaseTemporary, ' ', field_no_refe)

! ----- Create SD NUAGE
        call chpnua(ndim, inputFieldName, ' ', nuage1)
        call chpnua(ndim, field_no_refe, ' ', nuage2)

! ----- Projection on deformed mesh
        call pronua(method, nuage1, nuage2)
        call nuachp(nuage2, ' ', field_no_refe)

! ----- Set right mesh
        call dismoi("NUME_EQUA", field_no_refe, "CHAM_NO", repk=nume_equa)
        call jeveuo(nume_equa//'.REFN', 'E', vk24=fieldRefe)
        fieldRefe(1) = mesh_2

! ----- Generate relative speed field
        call jeexin(vite_curr(1:19)//'.VALE', ier)
        if (ier .gt. 0) then
            call gcncon('.', field_no_refe1)
            call copisd('CHAMP_GD', jvBaseTemporary, field_no_refe, field_no_refe1)
            call barych(field_no_refe1, vite_curr(1:19), 1.0d0, -1.0d0, &
                        field_no_refe, jvBaseTemporary)
        end if

! ----- Input fields
        call megeom(model, chgeom)
        call mecara(caraElem, chcara)
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PVITER'
        lchin(2) = field_no_refe
        lpain(3) = 'PVENTCX'
        lchin(3) = chcara(14)
        lpain(4) = 'PDEPLMR'
        lchin(4) = disp_prev
        lpain(5) = 'PDEPLPR'
        lchin(5) = disp_cumu_inst
        lpain(6) = 'PCAGNPO'
        lchin(6) = chcara(6)
        lpain(7) = 'PCAORIE'
        lchin(7) = chcara(1)
        lpain(8) = 'PSTRXMR'
        lchin(8) = strx_prev

! ----- Output fields
        lpaout(1) = 'PVECTUR'

! ----- Generate new RESU_ELEM name
        newnom = resuElem(10:16)
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        call corich('E', resuElem, ichin_=iLoad)

! ----- Compute
        call calcul('S', option, ligrel_calc, nbFieldIn, lchin, &
                    lpain, nbFieldOut, resuElem, lpaout, jvBaseTemporary, &
                    'OUI')
        call reajre(vectElem, resuElem, jvBaseTemporary)

! ----- Clean
        call detrsd('NUAGE', nuage1)
        call detrsd('NUAGE', nuage2)
        call jedetc(jvBaseTemporary, mesh_defo, 1)
        call detrsd('CHAMP_GD', field_no_refe1)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module loadCompute_module
