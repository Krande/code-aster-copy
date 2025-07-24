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
! Module to compute elementary components in non-linear
!
! ==================================================================================================
!
module NonLinearElem_module
! ==================================================================================================
    use NonLin_Datastructure_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: elemSuper, asseSuper
    public :: elemMass, asseMass
    public :: elemElas
    public :: elemDamp, asseDamp
    public :: elemDiri, elemNeum
    public :: elemGeom
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/asmatr.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/meamme.h"
#include "asterfort/mecgme.h"
#include "asterfort/medime.h"
#include "asterfort/memame.h"
#include "asterfort/memare.h"
#include "asterfort/merige.h"
#include "asterfort/mtdscr.h"
#include "asterfort/nmrigi.h"
#include "asterfort/utmess.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! elemSuper
!
! Compute elementary matrices for super-elements
!
! In  model            : model
! In  materialField    : field for material parameters
! In  caraElem         : field for elementary characteristics
! In  superElem        : name of elementary matrices for super elements
!
! --------------------------------------------------------------------------------------------------
    subroutine elemSuper(model, superElem)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=24), intent(in) :: model
        character(len=24), intent(in) :: superElem
! - Local
        character(len=1), parameter :: jvBase = "V"
        character(len=16), parameter :: option = "RIGI_MECA"
        integer(kind=8) :: ifm, niv
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_85')
        end if
        call memare(jvBase, superElem, model, option, ASTER_TRUE)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! asseSuper
!
! Assemble elementary matrices for super-elements
!
! In  numeDof          : name of numbering (NUME_DDL)
! In  listLoad         : name of datastructure for list of loads
! In  superElem        : name of elementary matrices for super elements
! In  superAsse        : name of assembled matrice for super elements
!
! --------------------------------------------------------------------------------------------------
    subroutine asseSuper(numeDof, listLoad, superElemZ, superAsseZ)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=14), intent(in) :: numeDof
        character(len=19), intent(in) :: listLoad
        character(len=*), intent(in) :: superElemZ, superAsseZ
! - Local
        character(len=1), parameter :: jvBase = "V"
        integer(kind=8) :: ifm, niv
        character(len=19) :: superElem
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_73')
        end if
        superElem = superElemZ(1:19)
        call asmatr(1, superElem, ' ', numeDof, &
                    listLoad, 'ZERO', jvBase, 1, superAsseZ)
        call mtdscr(superAsseZ)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! elemMass
!
! Compute elementary matrices for mass
!
! In  massOption       : type of mass to compute
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  materialField    : field for material parameters
! In  materialCoding   : field for coding material parameters
! In  behaviourField   : field for behaviour parameters
! In  sddyna           : name of dynamic parameters datastructure
! In  time             : value of time
! In  massElem         : name of elementary matrices for mass
!
! --------------------------------------------------------------------------------------------------
    subroutine elemMass(massOption, &
                        model, caraElem, &
                        materialField, materialCoding, &
                        behaviourField, &
                        time, massElem)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: massOption
        character(len=24), intent(in) :: model, caraElem
        character(len=24), intent(in) :: materialField, materialCoding
        character(len=24), intent(in) :: behaviourField
        real(kind=8), intent(in) :: time
        character(len=24), intent(in) :: massElem
! - Local
        character(len=1), parameter :: jvBase = "V"
        integer(kind=8) :: ifm, niv
        character(len=19) :: modelLigrel
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_82')
        end if
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
        call memame(massOption, model, materialField, materialCoding, &
                    caraElem, time, behaviourField, massElem, &
                    jvBase, modelLigrel)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! asseMass
!
! Assemble elementary matrices for mass
!
! In  lWithDirichlet   : flag to add Dirichlet matrix to mass (for explicit schemes)
! In  listLoad         : name of datastructure for list of loads
! In  numeDof          : name of numbering (NUME_DDL)
! In  numeDofFix       : name of numbering (NUME_DDL)
! In  massElem         : name of elementary matrices for mass
! In  diriElem         : name of elementary matrices for Dirichet
! In  massAsse         : name of assembled matrice for mass
!
! --------------------------------------------------------------------------------------------------
    subroutine asseMass(lWithDirichlet, listLoad, &
                        numeDof, numeDofFix, &
                        diriElemZ, massElemZ, massAsseZ)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        aster_logical, intent(in) :: lWithDirichlet
        character(len=19), intent(in) :: listLoad
        character(len=14), intent(in) :: numeDof, numeDofFix
        character(len=*), intent(in) :: diriElemZ, massElemZ, massAsseZ
! - Local
        character(len=1), parameter :: jvBase = "V"
        integer(kind=8) :: ifm, niv
        character(len=19) :: matrList(2)
        character(len=19) :: massElem
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_72')
        end if
        if (lWithDirichlet) then
            matrList(1) = massElemZ(1:19)
            matrList(2) = diriElemZ(1:19)
            call asmatr(2, matrList, ' ', numeDof, &
                        listLoad, 'ZERO', jvBase, 1, massAsseZ)
        else
            massElem = massElemZ(1:19)
            call asmatr(1, massElem, ' ', numeDofFix, &
                        listLoad, 'ZERO', jvBase, 1, massAsseZ)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! elemDiri
!
! Compute elementary matrices for Dirichet (B matrix for Lagrange multipliers)
!
! In  model            : model
! In  listLoad         : name of datastructure for list of loads
! In  diriElem         : name of elementary matrices for Dirichet
!
! --------------------------------------------------------------------------------------------------
    subroutine elemDiri(model, listLoad, diriElemZ)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=24), intent(in) :: model
        character(len=19), intent(in) :: listLoad
        character(len=*), intent(in) :: diriElemZ
! - Local
        character(len=1), parameter :: jvBase = "V"
        character(len=24) :: diriElem
        integer(kind=8) :: ifm, niv
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_80')
        end if
        diriElem = diriElemZ
        call medime(jvBase, 'ZERO', model, listLoad, diriElem)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! elemElas
!
! Compute elementary matrices for elasticity
!
! In  listFuncActi     : list of active functionnalities
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_system        : datastructure for non-linear system management
! In  sddyna           : name of dynamic parameters datastructure
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
!
! --------------------------------------------------------------------------------------------------
    subroutine elemElas(listFuncActi, &
                        model, caraElem, &
                        ds_material, ds_constitutive, &
                        ds_measure, ds_system, &
                        sddyna, &
                        hval_incr, hval_algo)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: listFuncActi(*)
        character(len=19), intent(in) :: sddyna
        character(len=24), intent(in) :: model, caraElem
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_Measure), intent(inout) :: ds_measure
        type(NL_DS_System), intent(in) :: ds_system
        character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
! - Local
        character(len=16), parameter :: optionElas = "RIGI_MECA"
        integer(kind=8) :: ifm, niv
        integer(kind=8), parameter :: iterNewtPred = 0
        integer(kind=8) :: ldccvg
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_86')
        end if
        call nmrigi(model, caraElem, &
                    ds_material, ds_constitutive, &
                    listFuncActi, iterNewtPred, sddyna, ds_measure, ds_system, &
                    hval_incr, hval_algo, &
                    optionElas, ldccvg)
        if (ldccvg .ne. 0) then
            call utmess('I', 'DAMPING1_1')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! elemDamp
!
! Compute elementary matrices for damp
!
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  materialField    : field for material parameters
! In  materialCoding   : field for coding material parameters
! In  behaviourField   : field for behaviour parameters
! In  vari             : internal state variables
! In  time             : value of time
! In  rigiElem         : name of elementary matrices for rigidity
! In  massElem         : name of elementary matrices for mass
! In  dampElem         : name of elementary matrices for damp
!
! --------------------------------------------------------------------------------------------------
    subroutine elemDamp(model, caraElem, &
                        materialField, materialCoding, &
                        behaviourField, &
                        vari, time, &
                        rigiElem, massElem, &
                        dampElem, sddyna)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=24), intent(in) :: model, caraElem
        character(len=24), intent(in) :: materialField, materialCoding
        character(len=24), intent(in) :: behaviourField
        character(len=24), intent(in) :: vari
        character(len=19), intent(in) :: sddyna
        real(kind=8), intent(in) :: time
        character(len=24), intent(in) :: rigiElem, massElem
        character(len=24), intent(in) :: dampElem
! - Local
        character(len=1), parameter :: jvBase = "V"
        integer(kind=8) :: ifm, niv

!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_83')
        end if
        call meamme(model, &
                    materialField, materialCoding, caraElem, &
                    time, jvBase, &
                    rigiElem, massElem, &
                    dampElem, &
                    vari, behaviourField, sddyna)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! asseDamp
!
! Assemble elementary matrices for damp
!
! In  numeDof          : name of numbering (NUME_DDL)
! In  listLoad         : name of datastructure for list of loads
! In  dampElem         : name of elementary matrices for damp
! In  dampAsse         : name of assembled matrice for damp
!
! --------------------------------------------------------------------------------------------------
    subroutine asseDamp(numeDof, listLoad, dampElemZ, dampAsseZ)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=14), intent(in) :: numeDof
        character(len=19), intent(in) :: listLoad
        character(len=24), intent(in) :: dampElemZ, dampAsseZ
! - Local
        character(len=1), parameter :: jvBase = "V"
        integer(kind=8) :: ifm, niv
        character(len=19) :: dampElem
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_71')
        end if
        dampElem = dampElemZ(1:19)
        call asmatr(1, dampElem, ' ', numeDof, &
                    listLoad, 'ZERO', jvBase, 1, dampAsseZ)
        call mtdscr(dampAsseZ)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! elemNeum
!
! Compute elementary matrices for Neumann loads (undead)
!
! In  listLoad         : name of datastructure for list of loads
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  materialField    : field for material parameters
! In  materialCoding   : field for coding material parameters
! In  behaviourField   : field for behaviour parameters
! In  timePrev         : previous time
! In  timeCurr         : current time
! In  dispPrev         : displacement at beginning of current time
! In  dispCumu         : displacement increment from beginning of current time
! In  neumElem         : name of elementary matrices for Neumann loads
!
! --------------------------------------------------------------------------------------------------
    subroutine elemNeum(listLoad, model, caraElem, &
                        materialField, materialCoding, &
                        behaviourField, &
                        timePrev, timeCurr, &
                        dispPrev, dispCumu, &
                        neumElem)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=19), intent(in) :: listLoad
        character(len=24), intent(in) :: model, caraElem
        character(len=24), intent(in) :: materialField, materialCoding
        character(len=24), intent(in) :: behaviourField
        real(kind=8), intent(in) :: timePrev, timeCurr
        character(len=19), intent(in) :: dispPrev, dispCumu
        character(len=24), intent(in) :: neumElem
! - Local
        integer(kind=8) :: ifm, niv
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_84')
        end if
        call mecgme("S", &
                    model, caraElem, materialField, materialCoding, behaviourField, &
                    listLoad, &
                    timePrev, timeCurr, &
                    dispPrev, dispCumu, &
                    neumElem)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! elemGeom
!
! Compute elementary matrices for geometry (buckling)
!
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  materialCoding   : field for coding material parameters
! In  strx             :
! In  sigm             : stress
! In  geomElem         : name of elementary matrices for geometry (buckling)
!
! --------------------------------------------------------------------------------------------------
    subroutine elemGeom(model, caraElem, &
                        materialCoding, &
                        strx, sigm, &
                        geomElem)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=24), intent(in) :: model, caraElem
        character(len=24), intent(in) :: materialCoding
        character(len=19), intent(in) :: strx, sigm
        character(len=24), intent(in) :: geomElem
! - Local
        character(len=1), parameter :: jvBase = "V"
        integer(kind=8) :: ifm, niv
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_81')
        end if
        call merige(model, caraElem, sigm, strx, geomElem, &
                    jvBase, 0, mateco=materialCoding)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module NonLinearElem_module
