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
module NonLinearDyna_module
! ==================================================================================================
    use NonLinearDyna_type
    use NonLin_Datastructure_type
    use Damping_type
    use Damping_module, only: dampModalPrintParameters, dampModalGetParameters, &
                              dampModalPreparation, dampComputeMatrix
    use NonLinearElem_module, only: elemMass, elemElas, asseMass
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: compMatrInit, compDampMatrix, asseMassMatrix
    public :: compAcceForce, compViteForce, compResiForce
    public :: dampGetParameters, dampPrintParameters
    public :: isDampMatrUpdate, isMassMatrAssemble
    public :: shiftMassMatrix
    private :: needElasMatrix, massGetType, compMassMatrix
! ==================================================================================================
    private
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/asmatr.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmchex.h"
#include "asterfort/utmess.h"
#include "asterfort/vtzero.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! dampGetParameters
!
! Get damping parameters from command file
!
! In  model            : model
! In  materialField    : field for material parameters
! In  caraElem         : field for elementary characteristics
! Out nlDynaDamping    : damping parameters
!
! --------------------------------------------------------------------------------------------------
    subroutine dampGetParameters(model, materialField, caraElem, &
                                 nlDynaDamping)
        ! - Parameters
        character(len=24), intent(in) :: model, materialField, caraElem
        type(NLDYNA_DAMPING), intent(out) :: nlDynaDamping
        ! - Local
        character(len=16), parameter :: factorKeyword = 'AMOR_MODAL'
        character(len=16) :: answer
        integer(kind=8) :: iret, nbOcc
        aster_logical :: hasMatrDamp, hasVectDamp, hasDamp
        aster_logical :: lDampRayleigh, lDampRayleighTang
        aster_logical :: lDampContact, lDampFEModel, lDampDiscret
        aster_logical :: lElemDampFromUser
        character(len=8) :: dampFromUser, dampElem
        aster_logical :: lElemDampToCompute
        aster_logical :: lDampModal, lDampModalReacVite
        type(MODAL_DAMPING) :: modalDamping
!   ------------------------------------------------------------------------------------------------
        !

        ! - Damping - Rayleigh
        lDampRayleigh = ASTER_FALSE
        call dismoi('EXI_AMOR_ALPHA', materialField, 'CHAM_MATER', repk=answer)
        if (answer .eq. "OUI") then
            lDampRayleigh = ASTER_TRUE
        end if
        call dismoi('EXI_AMOR_BETA', materialField, 'CHAM_MATER', repk=answer)
        if (answer .eq. "OUI") then
            lDampRayleigh = ASTER_TRUE
        end if
        if (lDampRayleigh) then
            call utmess('I', 'MECANONLINE5_7')
        end if

        ! - Damping - From contact/friction (DIS_CONTACT, JOINT)
        lDampContact = ASTER_FALSE
        call dismoi('EXI_AMOR_NOR', materialField, 'CHAM_MATER', repk=answer)
        if (answer .eq. "OUI") then
            lDampContact = ASTER_TRUE
        end if
        call dismoi('EXI_AMOR_TAN', materialField, 'CHAM_MATER', repk=answer)
        if (answer .eq. "OUI") then
            lDampContact = ASTER_TRUE
        end if

        ! ----- Update Damping matrix
        lDampRayleighTang = ASTER_FALSE
        if (lDampRayleigh) then
            call getvtx(' ', 'AMOR_RAYL_RIGI', scal=answer, nbret=iret)
            if (answer .eq. 'TANGENTE') then
                lDampRayleighTang = ASTER_TRUE
            end if
        end if

        ! - Damping - From super-elements or FLUI_ABSO Elements
        lDampFEModel = ASTER_FALSE
        call dismoi('EXI_AMOR', model, 'MODELE', repk=answer)
        if (answer .eq. "OUI") then
            lDampFEModel = ASTER_TRUE
        end if

        ! - Damping - From DIS_T elements
        lDampDiscret = ASTER_FALSE
        call dismoi('EXI_AMOR', caraElem, 'CARA_ELEM', repk=answer)
        if (answer .eq. "OUI") then
            lDampDiscret = ASTER_TRUE
        end if

        ! - Damping - From user
        call getvid(' ', 'MATR_ELEM_AMOR', scal=dampElem, nbret=iret)
        lElemDampFromUser = ASTER_FALSE
        if (iret .eq. 1) then
            lElemDampFromUser = ASTER_TRUE
            dampFromUser = dampElem
        end if

        ! - Modal damping
        lDampModalReacVite = ASTER_FALSE
        call getfac(factorKeyword, nbOcc)
        lDampModal = nbOcc .gt. 0
        if (lDampModal) then
            call dampModalGetParameters(factorKeyword, modalDamping)
            call dampModalPreparation(modalDamping)
        end if

        ! - Which cases damping is matrix or vector ?
        lElemDampToCompute = lDampRayleigh .or. lDampContact .or. lDampFEModel .or. lDampDiscret
        hasMatrDamp = lElemDampToCompute .or. lElemDampFromUser
        hasVectDamp = lDampModal
        hasDamp = hasMatrDamp .or. hasVectDamp

        ! - Save parameters
        nlDynaDamping%hasDamp = hasDamp
        nlDynaDamping%hasMatrDamp = hasMatrDamp
        nlDynaDamping%hasVectDamp = hasVectDamp
        nlDynaDamping%lElemDampToCompute = lElemDampToCompute
        nlDynaDamping%lDampRayleigh = lDampRayleigh
        nlDynaDamping%lDampRayleighTang = lDampRayleighTang
        nlDynaDamping%lDampContact = lDampContact
        nlDynaDamping%lDampFEModel = lDampFEModel
        nlDynaDamping%lDampDiscret = lDampDiscret
        nlDynaDamping%lDampModal = lDampModal
        nlDynaDamping%modalDamping = modalDamping
        nlDynaDamping%lElemDampFromUser = lElemDampFromUser
        nlDynaDamping%dampFromUser = dampFromUser
        nlDynaDamping%dampAsse = "&&NMCH6P.AMORT"
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isDampMatrUpdate
!
! Do the damping matrices need to be updated ?
!
! In  nlDynaDamping    : damping parameters
! In  l_renumber       : flag to renumbering
! Out lDampMatrUpdate  : flag if damp elementary matrices have to be updated
!
! --------------------------------------------------------------------------------------------------
    subroutine isDampMatrUpdate(nlDynaDamping, l_renumber, lDampMatrUpdate)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        aster_logical, intent(in) :: l_renumber
        aster_logical, intent(out) :: lDampMatrUpdate
        ! - Local
        aster_logical :: lDampRayleighTang, lDampMatrix
!   ------------------------------------------------------------------------------------------------
        !
        lDampMatrUpdate = ASTER_FALSE
        lDampMatrix = nlDynaDamping%hasMatrDamp
        lDampRayleighTang = nlDynaDamping%lDampRayleighTang
        if (lDampMatrix) then
            if (l_renumber .or. lDampRayleighTang .or. nlDynaDamping%lDampContact) then
                lDampMatrUpdate = ASTER_TRUE
            end if
        end if
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isMassMatrAssemble
!
! Do the mass matrices have to be assemble ?
!
! In  listFuncActi     : list of active functionnalities
! In  l_update_matr    : flag to update matrix
! Out lMassAssemble    : flag if mass elementary matrices have to be assembled
!
! --------------------------------------------------------------------------------------------------
    subroutine isMassMatrAssemble(listFuncActi, l_update_matr, lMassAssemble)
        ! - Parameters
        integer(kind=8), intent(in) :: listFuncActi(*)
        aster_logical, intent(in) :: l_update_matr
        aster_logical, intent(out) :: lMassAssemble
        ! - Local
        aster_logical :: lDyna
!   ------------------------------------------------------------------------------------------------
        !
        lMassAssemble = ASTER_FALSE
        lDyna = isfonc(listFuncActi, 'DYNAMIQUE')
        if (lDyna .and. l_update_matr) then
            lMassAssemble = ASTER_TRUE
        end if
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! needElasMatrix
!
! Do you need elastic matrix ?
!
! In  nlDynaDamping    : damping parameters
! Out lElas            : flag if elastic elementary matrices have to be calculated
!
! --------------------------------------------------------------------------------------------------
    subroutine needElasMatrix(nlDynaDamping, lElas)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        aster_logical, intent(out) :: lElas
        ! - Local
        aster_logical :: lDampRayleigh, lDampRayleighTang
!   ------------------------------------------------------------------------------------------------
        !
        lElas = ASTER_FALSE
        lDampRayleigh = nlDynaDamping%lDampRayleigh
        lDampRayleighTang = nlDynaDamping%lDampRayleighTang
        if (lDampRayleigh .and. .not. lDampRayleighTang) then
            lElas = ASTER_TRUE
        end if
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! massGetType
!
! Get type of mass matrix
!
! In  sddyna           : name of dynamic parameters datastructure
! Out massOption       : name of option for mass matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine massGetType(sddyna, massOption)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        character(len=19), intent(in) :: sddyna
        character(len=16), intent(out) :: massOption
        ! - Local
        aster_logical :: lexpl
!   ------------------------------------------------------------------------------------------------
        !
        massOption = ' '
        lexpl = ndynlo(sddyna, 'EXPLICITE')
        if (lexpl) then
            if (ndynlo(sddyna, 'MASS_DIAG')) then
                massOption = 'MASS_MECA_EXPLI'
            else
                massOption = 'MASS_MECA'
            end if
        else
            massOption = 'MASS_MECA'
        end if
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compMatrInit
!
! Prepare matrices at initial state
!
! In  listFuncActi     : list of active functionnalities
! In  sddyna           : name of dynamic parameters datastructure
! In  nlDynaDamping    : damping parameters
! In  sddisc           : datastructure for time discretization
! In  listLoad         : name of datastructure for list of loads
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_system        : datastructure for non-linear system management
! In  numeDof          : name of numbering (NUME_DDL)
! In  numeDofFix       : name of numbering (NUME_DDL)
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_meelem      : hat-variable for elementary matrix
! In  hval_measse      : hat-variable for matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compMatrInit(listFuncActi, &
                            sddyna, nlDynaDamping, &
                            sddisc, listLoad, &
                            model, caraElem, &
                            ds_material, ds_constitutive, &
                            ds_measure, ds_system, &
                            numeDof, numeDofFix, &
                            hval_incr, hval_algo, &
                            hval_meelem, hval_measse)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        integer(kind=8), intent(in) :: listFuncActi(*)
        character(len=19), intent(in) :: sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        character(len=19), intent(in) :: sddisc, listLoad
        character(len=24), intent(in) :: model, caraElem
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_Measure), intent(inout) :: ds_measure
        type(NL_DS_System), intent(in) :: ds_system
        character(len=14), intent(in) :: numeDof, numeDofFix
        character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
        character(len=19), intent(in) :: hval_meelem(*), hval_measse(*)
        ! - Local
        integer(kind=8) :: ifm, niv
        aster_logical :: lElas, lVarc, lExpl, lWithDirichlet
        aster_logical :: lDampMatrix
        integer(kind=8), parameter :: numeInstInit = 0
        real(kind=8) :: timeInit
        character(len=16) :: massOption
!   ------------------------------------------------------------------------------------------------
        !
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_22')
        end if

        ! - Initial tome
        timeInit = diinst(sddisc, numeInstInit)

        ! - Active functionnalities
        lVarc = isfonc(listFuncActi, 'EXI_VARC')
        lExpl = ndynlo(sddyna, 'EXPLICITE')
        lDampMatrix = nlDynaDamping%hasMatrDamp

        ! - Compute elementary matrices for elasticity ?
        call needElasMatrix(nlDynaDamping, lElas)

        ! - Compute elementary matrices for elasticity
        if (lElas) then
            if (lVarc) then
                call utmess('F', 'MECANONLINE3_2')
            end if
            call elemElas(listFuncActi, &
                          model, caraElem, &
                          ds_material, ds_constitutive, &
                          ds_measure, ds_system, &
                          sddyna, &
                          hval_incr, hval_algo)
        end if

        ! - Compute mass matrix
        call massGetType(sddyna, massOption)
        lWithDirichlet = lExpl
        call compMassMatrix(model, caraElem, &
                            ds_material, ds_constitutive, &
                            timeInit, listLoad, &
                            numeDof, numeDofFix, &
                            hval_meelem, hval_measse, &
                            massOption, lWithDirichlet)

        ! - Compute damping matrix
        if (lDampMatrix) then
            call compDampMatrix(model, caraElem, &
                                ds_material, ds_constitutive, &
                                timeInit, listLoad, numeDof, nlDynaDamping, &
                                ds_system, hval_incr, hval_meelem, sddyna)
        end if
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compDampMatrix
!
! Compute damping matrix
!
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  time             : value of time
! In  listLoad         : name of datastructure for list of loads
! In  numeDof          : name of numbering (NUME_DDL)
! In  nlDynaDamping    : damping parameters
! In  ds_system        : datastructure for non-linear system management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_meelem      : hat-variable for elementary matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compDampMatrix(model, caraElem, &
                              ds_material, ds_constitutive, &
                              time, listLoad, numeDof, nlDynaDamping, &
                              ds_system, hval_incr, hval_meelem, sddyna)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        character(len=24), intent(in) :: model, caraElem
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        real(kind=8), intent(in) :: time
        character(len=19), intent(in) :: listLoad, sddyna
        character(len=14), intent(in) :: numeDof
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        type(NL_DS_System), intent(in) :: ds_system
        character(len=19), intent(in) :: hval_incr(*), hval_meelem(*)
        ! - Local
        character(len=1), parameter :: jvBase = "V"
        character(len=24) :: dampAsse, variPrev
        character(len=24) :: massElem, rigiElem, dampElem
!   ------------------------------------------------------------------------------------------------
        !
        dampAsse = nlDynaDamping%dampAsse
        call detrsd('MATR_ASSE', dampAsse)

        ! - Compute damping matrix from elementary matrices
        if (nlDynaDamping%lElemDampToCompute) then
            call nmchex(hval_meelem, 'MEELEM', 'MEMASS', massElem)
            rigiElem = ds_system%merigi
            call nmchex(hval_incr, 'VALINC', 'VARMOI', variPrev)
            call dampComputeMatrix(model, caraElem, &
                                   ds_material%mater, ds_material%mateco, &
                                   ds_constitutive%compor, &
                                   variPrev, time, listLoad, numeDof, &
                                   rigiElem, massElem, &
                                   dampAsse, sddyna)
        end if

        ! - Get damping matrix from user
        if (nlDynaDamping%lElemDampFromUser) then
            dampElem = nlDynaDamping%dampFromUser
            if (nlDynaDamping%lElemDampToCompute) then
                call asmatr(1, dampElem, ' ', numeDof, &
                            listLoad, 'CUMU', jvBase, 1, dampAsse)
            else
                call asmatr(1, dampElem, ' ', numeDof, &
                            listLoad, 'ZERO', jvBase, 1, dampAsse)
            end if
        end if
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compMassMatrix
!
! Compute mass matrix
!
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  time             : value of time
! In  listLoad         : name of datastructure for list of loads
! In  numeDof          : name of numbering (NUME_DDL)
! In  numeDofFix       : name of numbering (NUME_DDL)
! In  hval_meelem      : hat-variable for elementary matrix
! In  hval_measse      : hat-variable for matrix
! In  massOption       : type of mass to compute
! In  lWithDirichlet   : flag to add Dirichlet matrix to mass (for explicit schemes)
!
! --------------------------------------------------------------------------------------------------
    subroutine compMassMatrix(model, caraElem, &
                              ds_material, ds_constitutive, &
                              time, listLoad, &
                              numeDof, numeDofFix, &
                              hval_meelem, hval_measse, &
                              massOption_, lWithDirichlet_)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        character(len=24), intent(in) :: model, caraElem
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        real(kind=8), intent(in) :: time
        character(len=19), intent(in) :: listLoad
        character(len=14), intent(in) :: numeDof, numeDofFix
        character(len=19), intent(in) :: hval_meelem(*), hval_measse(*)
        character(len=16), optional, intent(in) :: massOption_
        aster_logical, optional, intent(in) :: lWithDirichlet_
        ! - Local
        character(len=16) :: massOption
        aster_logical:: lWithDirichlet
        character(len=24) :: massAsse
        character(len=24) :: massElem, diriElem
!   ------------------------------------------------------------------------------------------------
        !
        massOption = 'MASS_MECA'
        lWithDirichlet = ASTER_FALSE
        if (present(massOption_)) then
            ASSERT(present(lWithDirichlet_))
            massOption = massOption_
            lWithDirichlet = lWithDirichlet_
        end if

        ! - Compute elementary matrices for mass
        call nmchex(hval_meelem, 'MEELEM', 'MEMASS', massElem)
        call elemMass(massOption, &
                      model, caraElem, &
                      ds_material%mater, ds_material%mateco, &
                      ds_constitutive%compor, &
                      time, massElem)

        ! - Assemble elementary matrices for mass
        diriElem = " "
        if (lWithDirichlet) then
            call nmchex(hval_meelem, 'MEELEM', 'MEDIRI', diriElem)
        end if
        call nmchex(hval_measse, 'MEASSE', 'MEMASS', massAsse)
        call asseMass(lWithDirichlet, listLoad, &
                      numeDof, numeDofFix, &
                      diriElem, massElem, massAsse)
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! asseMassMatrix
!
! Assemble mass matrix
!
! In  listLoad         : name of datastructure for list of loads
! In  numeDof          : name of numbering (NUME_DDL)
! In  numeDofFix       : name of numbering (NUME_DDL)
! In  hval_meelem      : hat-variable for elementary matrix
! In  hval_measse      : hat-variable for matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine asseMassMatrix(listLoad, &
                              numeDof, numeDofFix, &
                              hval_meelem, hval_measse)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        character(len=19), intent(in) :: listLoad
        character(len=14), intent(in) :: numeDof, numeDofFix
        character(len=19), intent(in) :: hval_meelem(*), hval_measse(*)
        ! - Local
        aster_logical, parameter :: lWithDirichlet = ASTER_FALSE
        character(len=24), parameter :: diriElem = " "
        character(len=24) :: massAsse, massElem
!   ------------------------------------------------------------------------------------------------
        !
        call nmchex(hval_meelem, 'MEELEM', 'MEMASS', massElem)
        call nmchex(hval_measse, 'MEASSE', 'MEMASS', massAsse)
        call asseMass(lWithDirichlet, listLoad, &
                      numeDof, numeDofFix, &
                      diriElem, massElem, massAsse)
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dampPrintParameters
!
! Print damping parameters
!
! In  nlDynaDamping    : damping parameters
!
! --------------------------------------------------------------------------------------------------
    subroutine dampPrintParameters(nlDynaDamping)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
!   ------------------------------------------------------------------------------------------------
        !
        if (nlDynaDamping%hasDamp) then
            call utmess('I', 'MECANONLINE15_30')
            if (nlDynaDamping%hasMatrDamp) then
                call utmess('I', 'MECANONLINE15_31')
            end if
            if (nlDynaDamping%hasVectDamp) then
                call utmess('I', 'MECANONLINE15_32')
            end if
            if (nlDynaDamping%lDampRayleigh) then
                if (nlDynaDamping%lDampRayleighTang) then
                    call utmess('I', 'MECANONLINE15_34')
                else
                    call utmess('I', 'MECANONLINE15_33')
                end if
            end if
            if (nlDynaDamping%lDampContact) then
                call utmess('I', 'MECANONLINE15_35')
            end if
            if (nlDynaDamping%lDampFEModel) then
                call utmess('I', 'MECANONLINE15_36')
            end if
            if (nlDynaDamping%lDampDiscret) then
                call utmess('I', 'MECANONLINE15_37')
            end if
            if (nlDynaDamping%lDampModal) then
                call utmess('I', 'MECANONLINE15_38')
                call dampModalPrintParameters(nlDynaDamping%modalDamping)
            end if
        else
            call utmess('I', 'MECANONLINE15_29')
        end if
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! shiftMassMatrix
!
! Shift mass matrix
!
! In  sddyna           : name of dynamic parameters datastructure
! In  hval_measse      : hat-variable for matrix
! In  numeTime         : index of current time step
!
! --------------------------------------------------------------------------------------------------
    subroutine shiftMassMatrix(sddyna, numeTime, hval_measse)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        integer(kind=8), intent(in) :: numeTime
        character(len=19), intent(in) :: sddyna
        character(len=19), intent(in) :: hval_measse(*)
        ! - Local
        real(kind=8) :: coefVale(3)
        character(len=24) :: matrName(3)
        character(len=4) :: coefType(3)
        aster_logical :: lExpl, lShiftMass, lFirstStep
        real(kind=8) :: coefShiftMass
        character(len=19) :: rigiAsse, massAsse
!   ------------------------------------------------------------------------------------------------
        !
        lFirstStep = numeTime .le. 1
        coefShiftMass = ndynre(sddyna, 'COEF_MASS_SHIFT')
        lExpl = ndynlo(sddyna, 'EXPLICITE')
        lShiftMass = ndynlo(sddyna, 'COEF_MASS_SHIFT')

        ! - Get name of matrices
        call nmchex(hval_measse, 'MEASSE', 'MERIGI', rigiAsse)
        call nmchex(hval_measse, 'MEASSE', 'MEMASS', massAsse)

        ! - To combine
        coefType(1) = 'R'
        coefType(2) = 'R'
        coefVale(1) = 1.d0
        coefVale(2) = coefShiftMass
        matrName(1) = massAsse
        matrName(2) = rigiAsse

        ! - Combination
        if (lShiftMass .and. lFirstStep) then
            if (lExpl) then
                call mtcmbl(2, coefType, coefVale, matrName, massAsse, ' ', ' ', 'ELIM=')
            else
                call mtcmbl(2, coefType, coefVale, matrName, massAsse, 'LAGR', ' ', 'ELIM=')
            end if
        end if
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compAcceForce
!
! Compute force (RHS) for acceleration (inertial effect)
!
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_measse      : hat-variable for matrix
! In  acceForce        : name of field for force (inertial effect)
!
! --------------------------------------------------------------------------------------------------
    subroutine compAcceForce(hval_incr, hval_measse, acceForce)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        character(len=19), intent(in) :: hval_incr(*), hval_measse(*)
        character(len=19), intent(in) :: acceForce
        ! - Local
        character(len=19) :: massAsse, acceCurr
        integer(kind=8) :: jvMass
        real(kind=8), pointer :: acce(:) => null()
        real(kind=8), pointer :: force(:) => null()
!   ------------------------------------------------------------------------------------------------
        !
        call vtzero(acceForce)

        ! - Get name of matrices and vectors
        call nmchex(hval_measse, 'MEASSE', 'MEMASS', massAsse)
        call nmchex(hval_incr, 'VALINC', 'ACCPLU', acceCurr)

        ! - Get access
        call jeveuo(acceCurr(1:19)//'.VALE', 'L', vr=acce)
        call jeveuo(massAsse(1:19)//'.&INT', 'L', jvMass)
        call jeveuo(acceForce(1:19)//'.VALE', 'E', vr=force)

        ! - Compute
        call mrmult('ZERO', jvMass, acce, force, 1, ASTER_TRUE)
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compViteForce
!
! Compute force (RHS) for speed (damping effect)
!
! In  nlDynaDamping    : damping parameters
! In  hval_incr        : hat-variable for incremental values fields
! In  name_vite        : name of speed hat-variable ('VITMOI', 'VITPLU')
! In  viteForce        : name of field for force (damping effect)
!
! --------------------------------------------------------------------------------------------------
    subroutine compViteForce(nlDynaDamping, hval_incr, name_vite, viteForce)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        character(len=19), intent(in) :: hval_incr(*)
        character(len=6), intent(in) :: name_vite
        character(len=19), intent(in) :: viteForce
        ! - Local
        character(len=19) :: dampAsse, viteCurr
        integer(kind=8) :: jvDamp
        real(kind=8), pointer :: vite(:) => null()
        real(kind=8), pointer :: force(:) => null()
!   ------------------------------------------------------------------------------------------------
        !
        call vtzero(viteForce)

        ! - Get name of matrices and vectors
        dampAsse = nlDynaDamping%dampAsse
        call nmchex(hval_incr, 'VALINC', name_vite, viteCurr)

        ! - Get access
        call jeveuo(viteCurr(1:19)//'.VALE', 'L', vr=vite)
        call jeveuo(dampAsse(1:19)//'.&INT', 'L', jvDamp)
        call jeveuo(viteForce(1:19)//'.VALE', 'E', vr=force)

        ! - Compute
        call mrmult('ZERO', jvDamp, vite, force, 1, ASTER_TRUE)
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compResiForce
!
! Compute force for residual computation
!
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_measse      : hat-variable for matrix
! In  resiForce        : name of field for force (inertial effect)
!
! --------------------------------------------------------------------------------------------------
    subroutine compResiForce(hval_incr, hval_measse, resiForce)
!   ------------------------------------------------------------------------------------------------
        ! - Parameters
        character(len=19), intent(in) :: hval_incr(*), hval_measse(*)
        character(len=19), intent(in) :: resiForce
        ! - Local
        character(len=19) :: massAsse, viteCurr
        integer(kind=8) :: jvMass
        real(kind=8), pointer :: vite(:) => null()
        real(kind=8), pointer :: force(:) => null()
!   ------------------------------------------------------------------------------------------------
        !
        call vtzero(resiForce)

        ! - Get name of matrices and vectors
        call nmchex(hval_measse, 'MEASSE', 'MEMASS', massAsse)
        call nmchex(hval_incr, 'VALINC', 'VITPLU', viteCurr)

        ! - Get access
        call jeveuo(viteCurr(1:19)//'.VALE', 'L', vr=vite)
        call jeveuo(massAsse(1:19)//'.&INT', 'L', jvMass)
        call jeveuo(resiForce(1:19)//'.VALE', 'E', vr=force)

        ! - Compute
        call mrmult('ZERO', jvMass, vite, force, 1, ASTER_TRUE)
        !
!   -----------------------------------------------------------------------------------------------
    end subroutine
end module NonLinearDyna_module
