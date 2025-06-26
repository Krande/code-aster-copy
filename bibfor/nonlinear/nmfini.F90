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
subroutine nmfini(sddyna, nlDynaDamping, &
                  valinc, measse, model, ds_material, &
                  caraElem, ds_constitutive, ds_system, &
                  ds_measure, sddisc, numeTime, &
                  solalg, numeDof, listFuncActi)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nonlinNForceCompute.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assvec.h"
!
    character(len=19), intent(in) :: valinc(*), measse(*)
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=24), intent(in) :: model, caraElem
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=24), intent(in) :: numeDof
    character(len=19), intent(in) :: sddisc, solalg(*)
    integer(kind=8), intent(in) :: numeTime
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - CALCUL)
!
! CALCUL DES ENERGIES
! INITIALISATION DES VECTEURS DE FORCE POUR LE CALCUL DES ENERGIES
!
! --------------------------------------------------------------------------------------------------
!
! In  nlDynaDamping    : damping parameters
! In  sddyna           : datastructure for dynamic
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  MODELE : MODELE
! In  listFuncActi           : list of active functionnalities
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_system        : datastructure for non-linear system management
! IO  ds_measure       : datastructure for measure and statistics management
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  NUMINS : NUMERO D'INSTANT
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  NUMEDD : NUME_DDL
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: masse, dampAsse, vitmoi, accmoi
    character(len=19) :: fexmoi, fammoi, flimoi
    integer(kind=8) :: imasse, iamort
    integer(kind=8) :: nbEqua, iEqua
    aster_logical :: lDampMatrix, ldyna
    character(len=19) :: fnomoi
    real(kind=8) :: timePrev, timeCurr
    real(kind=8), pointer :: cv(:) => null()
    real(kind=8), pointer :: ma(:) => null()
    real(kind=8), pointer :: ccmo(:) => null()
    real(kind=8), pointer :: cnfno(:) => null()
    real(kind=8), pointer :: fammo(:) => null()
    real(kind=8), pointer :: fexmo(:) => null()
    real(kind=8), pointer :: flimo(:) => null()
    real(kind=8), pointer :: fnomo(:) => null()
    real(kind=8), pointer :: vitmo(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    ldyna = ndynlo(sddyna, 'DYNAMIQUE')
    lDampMatrix = nlDynaDamping%hasMatrDamp
!
    call nmchex(valinc, 'VALINC', 'FEXMOI', fexmoi)
    call jeveuo(fexmoi//'.VALE', 'E', vr=fexmo)
    call dismoi('NB_EQUA', numeDof, 'NUME_DDL', repi=nbEqua)

! - Get time
    ASSERT(numeTime .eq. 1)
    timePrev = diinst(sddisc, numeTime-1)
    timeCurr = diinst(sddisc, numeTime)

! - AJOUT DE LA FORCE DE LIAISON ET DE LA FORCE D AMORTISSEMENT MODAL
    call nmchex(valinc, 'VALINC', 'FAMMOI', fammoi)
    call jeveuo(fammoi//'.VALE', 'L', vr=fammo)
    call nmchex(valinc, 'VALINC', 'FLIMOI', flimoi)
    call jeveuo(flimoi//'.VALE', 'L', vr=flimo)
    do iEqua = 1, nbEqua
        fexmo(iEqua) = fammo(iEqua)+flimo(iEqua)
    end do

! --- AJOUT DU TERME C.V
    if (lDampMatrix) then
        dampAsse = nlDynaDamping%dampAsse
        call mtdscr(dampAsse)
        call jeveuo(dampAsse//'.&INT', 'L', iamort)
        call nmchex(valinc, 'VALINC', 'VITMOI', vitmoi)
        call jeveuo(vitmoi//'.VALE', 'L', vr=vitmo)
        AS_ALLOCATE(vr=cv, size=nbEqua)
        call mrmult('ZERO', iamort, vitmo, cv, 1, .true._1)
        do iEqua = 1, nbEqua
            fexmo(iEqua) = fexmo(iEqua)+cv(iEqua)
        end do
        AS_DEALLOCATE(vr=cv)
    end if
!
! --- AJOUT DU TERME M.A
!
    if (ldyna) then
        call nmchex(measse, 'MEASSE', 'MEMASS', masse)
        call mtdscr(masse)
        call jeveuo(masse//'.&INT', 'L', imasse)
        call nmchex(valinc, 'VALINC', 'ACCMOI', accmoi)
        call jeveuo(accmoi//'.VALE', 'L', vr=ccmo)
        AS_ALLOCATE(vr=ma, size=nbEqua)
        call mrmult('ZERO', imasse, ccmo, ma, 1, &
                    .true._1)
        do iEqua = 1, nbEqua
            fexmo(iEqua) = fexmo(iEqua)+ma(iEqua)
        end do
        AS_DEALLOCATE(vr=ma)
    end if
!
! - Direct computation (no integration of behaviour)
!
    call nonlinNForceCompute(model, caraElem, listFuncActi, &
                             ds_material, ds_constitutive, &
                             ds_measure, ds_system, &
                             valinc, solalg)
    call assvec('V', ds_system%cnfnod, 1, ds_system%vefnod, [1.d0], &
                ds_system%nume_dof)
    call jeveuo(ds_system%cnfnod//'.VALE', 'L', vr=cnfno)
    do iEqua = 1, nbEqua
        fexmo(iEqua) = fexmo(iEqua)+cnfno(iEqua)
    end do
!
! --- INITIALISATION DES FORCES INTERNES
!
    call nmchex(valinc, 'VALINC', 'FNOMOI', fnomoi)
    call jeveuo(fnomoi//'.VALE', 'E', vr=fnomo)
    do iEqua = 1, nbEqua
        fnomo(iEqua) = cnfno(iEqua)
    end do
!
end subroutine
