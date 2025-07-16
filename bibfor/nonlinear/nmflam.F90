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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmflam(optionSpec, &
                  model, ds_material, caraElem, listLoad, listFuncActi, &
                  numeDof, ds_system, &
                  ds_constitutive, &
                  sddisc, numeTime, &
                  sddyna, sderro, ds_algopara, &
                  ds_measure, &
                  hval_incr, hval_algo, &
                  hval_meelem, &
                  ds_posttimestep)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/diinst.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmcrel.h"
#include "asterfort/nmflal.h"
#include "asterfort/nmflin.h"
#include "asterfort/nmflma.h"
#include "asterfort/nmop45.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/vpSorensen.h"
#include "asterfort/vpleci.h"
#include "asterfort/nonlinDSPostTimeStepSave.h"
!
    character(len=16), intent(in) :: optionSpec
    character(len=24), intent(in) :: model, caraElem
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=19), intent(in) :: listLoad
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=24), intent(in) :: numeDof
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeTime
    character(len=19), intent(in) :: sddyna
    character(len=24), intent(in) :: sderro
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=19), intent(in) :: hval_meelem(*)
    type(NL_DS_PostTimeStep), intent(inout) :: ds_posttimestep
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Initializations
!
! Spectral analysis (MODE_VIBR/CRIT_STAB)
!
! --------------------------------------------------------------------------------------------------
!
! In  optionSpec       : option to compute (FLAMBSTA/FLAMBDYN/VIBRDYNA)
! In  model            : name of model
! In  ds_material      : datastructure for material parameters
! In  caraElem         : name of elementary characteristics (field)
! In  listLoad         : datastructure for list of loads
! In  listFuncActi     : list of active functionnalities
! In  numeDof          : name of numbering (NUME_DDL)
! In  nume_dof_inva    : name of reference numbering (invariant)
! In  ds_constitutive  : datastructure for constitutive laws management
! In  sddisc           : datastructure for time discretization
! In  numeTime         : index of current time step
! In  sddyna           : datastructure for dynamic
! In  sderro           : datastructure for error management (events)
! In  ds_algopara      : datastructure for algorithm parameters
! In  ds_system        : datastructure for non-linear system management
! IO  ds_measure       : datastructure for measure and statistics management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_meelem      : hat-variable for elementary matrix
! In  hval_measse      : hat-variable for matrix
! IO  ds_posttimestep  : datastructure for post-treatment at each time step
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: linsta, l_hpp, lModiRigi
    integer(kind=8) :: nbFreq, nbFreqCalc
    integer(kind=8) :: nbDofExcl, nbDofStab, coefDimSpace
    integer(kind=8) :: jvPara, iFreq, nfreq_calibr
    real(kind=8) :: freq_calc, freq_mini_abso, freq_abso, freq_mini
    real(kind=8) :: bande(2), r8bid, timeCurr
    character(len=4) :: mod45
    character(len=16) :: optionModal, matrType, calcLevel
    character(len=19) :: matrAsse, matrGeom
    character(len=24) :: k24bid
    character(len=19), parameter :: eigsol = '&&NMFLAM.EIGSOL'
    character(len=8), parameter :: sdmode = '&&NM45BI', sdstab = '&&NM45SI'
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    timeCurr = diinst(sddisc, numeTime)

! - Get parameters
    call nmflal(optionSpec, ds_posttimestep, &
                mod45, l_hpp, &
                nbFreq, coefDimSpace, matrType, optionModal, bande, &
                nbDofExcl, nbDofStab, lModiRigi, calcLevel)

! - Compute matrixes
    call nmflma(matrType, mod45, &
                l_hpp, lModiRigi, &
                listFuncActi, ds_algopara, &
                model, caraElem, &
                ds_material, ds_constitutive, &
                sddyna, listLoad, &
                sddisc, numeTime, &
                ds_posttimestep, nbDofExcl, &
                hval_incr, hval_algo, &
                numeDof, ds_system, &
                ds_measure, hval_meelem, &
                matrAsse, matrGeom)

! - CREATION ET REMPLISSAGE DE LA SD EIGSOL pour GEP SYM REEL RESOLU VIA SORENSEN
    call vpSorensen(mod45, matrAsse, matrGeom, &
                    optionModal, calcLevel, &
                    coefDimSpace, nbFreq, bande, &
                    eigsol)

! - Compute eigen values/vectors
    call nmop45(eigsol, l_hpp, mod45, sdmode, sdstab, ds_posttimestep, nfreq_calibr)
    call vpleci(eigsol, 'I', 1, k24bid, r8bid, nbFreqCalc)
    call detrsd('EIGENSOLVER', eigsol)

! - No eigen value
    if (nbFreqCalc .eq. 0) then
        goto 999
    end if

! - Print info
    do iFreq = 1, nbFreqCalc
        if (mod45 .eq. 'VIBR') then
            call rsadpa(sdmode, 'L', 1, 'FREQ', iFreq, 0, sjv=jvPara)
            call utmess('I', 'MECANONLINE6_10', si=iFreq, sr=zr(jvPara))
        else if (mod45 .eq. 'FLAM') then
            call rsadpa(sdmode, 'L', 1, 'CHAR_CRIT', iFreq, 0, sjv=jvPara)
            call utmess('I', 'MECANONLINE6_11', si=iFreq, sr=zr(jvPara))
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
    iFreq = 1
    if (nbDofStab .ne. 0) then
        call rsadpa(sdstab, 'L', 1, 'CHAR_STAB', iFreq, 0, sjv=jvPara)
        call utmess('I', 'MECANONLINE6_12', si=iFreq, sr=zr(jvPara))
    end if

! - DETECTION INSTABILITE SI DEMANDE
    if (mod45 .eq. 'FLAM') then
        freq_mini_abso = r8maem()
        do iFreq = 1, nbFreqCalc
            call rsadpa(sdmode, 'L', 1, 'CHAR_CRIT', iFreq, 0, sjv=jvPara)
            freq_calc = zr(jvPara)
            freq_abso = abs(freq_calc)
            if (freq_abso .lt. freq_mini_abso) then
                freq_mini_abso = freq_abso
                freq_mini = freq_calc
            end if
        end do
        call nmflin(ds_posttimestep, matrAsse, freq_mini, linsta)
        call nmcrel(sderro, 'CRIT_STAB', linsta)
    end if
!
999 continue

! - Save results
    call nonlinDSPostTimeStepSave(mod45, sdmode, sdstab, &
                                  timeCurr, numeTime, nbFreqCalc, &
                                  nfreq_calibr, ds_posttimestep)

! - Cleaning
    call jedetc('G', sdmode, 1)
    call jedetc('G', sdstab, 1)
!
    call jedema()
end subroutine
