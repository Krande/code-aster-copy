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
subroutine nmmatr(phaseType, listFuncActi, listLoad, numeDof, &
                  sddyna, nlDynaDamping, &
                  numeTime, ds_contact, meelem, measse, &
                  matrAsse)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use NonLinearDyna_module, only: shiftMassMatrix
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ascoma.h"
#include "asterfort/cfdisl.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/lccmst.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmasfr.h"
#include "asterfort/nmasun.h"
#include "asterfort/nmchex.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: phaseType, listFuncActi(*)
    character(len=19), intent(in) :: listLoad
    character(len=24), intent(in) :: numeDof
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    integer(kind=8), intent(in) :: numeTime
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: meelem(*), measse(*)
    character(len=19), intent(inout) :: matrAsse
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - CALCUL)
!
! ASSEMBLAGE DE LA MATRICE GLOBALE
!
! --------------------------------------------------------------------------------------------------
!
! In  phaseType        : name of current phase of algorithm
! In  ds_contact       : datastructure for contact management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  numeTime : NUMERO D'INSTANT
! IN  numeDof : NOM DE LA NUMEROTATION MECANIQUE
! IN  listLoad : SD LISTE DES CHARGES
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! IO  MATASS : MATRICE ASSEMBLEE RESULTANTE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    aster_logical :: lDyna, lExpl, lDampMatrix, lNeumUndead
    aster_logical :: lContDiscret, lContLAC
    real(kind=8) :: coefRigi, coefDamp, coefMass
    real(kind=8) :: coefVale(3)
    character(len=24) :: matrName(3)
    character(len=4), parameter :: coefType(3) = (/'R', 'R', 'R'/)
    integer(kind=8) :: nbMatr
    character(len=19) :: rigiAsse, massAsse, dampAsse, matrRefe
    character(len=24) :: matrType
    aster_logical :: lUnil, lUnilPena
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_41')
    end if

! - Get name of matrices
    call nmchex(measse, 'MEASSE', 'MERIGI', rigiAsse)
    call nmchex(measse, 'MEASSE', 'MEMASS', massAsse)
    dampAsse = nlDynaDamping%dampAsse

! - Active functionnalites
    lContDiscret = isfonc(listFuncActi, 'CONT_DISCRET')
    lUnil = isfonc(listFuncActi, 'LIAISON_UNILATER')
    lNeumUndead = isfonc(listFuncActi, 'NEUM_UNDEAD')
    lContLAC = isfonc(listFuncActi, 'CONT_LAC')
    lDyna = isfonc(listFuncActi, 'DYNAMIQUE')
    lDampMatrix = nlDynaDamping%hasMatrDamp
    lExpl = ndynlo(sddyna, 'EXPLICITE')

! - Delete previous matrix
    if (lDyna) then
        call detrsd('MATR_ASSE', matrAsse)
    end if

! - Get coefficients for matrixes
    if (lDyna) then
        coefRigi = ndynre(sddyna, 'COEF_MATR_RIGI')
        coefDamp = ndynre(sddyna, 'COEF_MATR_AMOR')
        coefMass = ndynre(sddyna, 'COEF_MATR_MASS')
    else
        coefRigi = 1.d0
    end if

! - Shift mass matrix
    if (lDyna .and. phaseType .eq. PRED_EULER) then
        call shiftMassMatrix(sddyna, numeTime, measse)
    end if

! - Matrixes and coefficients
    if (lDyna) then
        if (phaseType .eq. ACCEL_INIT) then
            matrName(1) = massAsse
            nbMatr = 1
            coefVale(1) = 1.d0
        else
            if (lExpl) then
                matrName(1) = massAsse
                nbMatr = 1
                coefVale(1) = coefMass
            else
                coefVale(1) = coefRigi
                coefVale(2) = coefMass
                coefVale(3) = coefDamp
                matrName(1) = rigiAsse
                matrName(2) = massAsse
                matrName(3) = dampAsse
                if (lDampMatrix) then
                    nbMatr = 3
                else
                    nbMatr = 2
                end if
            end if
        end if
    end if

! - Define matrix
    if (lDyna) then
        if (phaseType .eq. ACCEL_INIT) then
            call mtdefs(matrAsse, massAsse, 'V', 'R')
        else
            if (lExpl) then
                call mtdefs(matrAsse, massAsse, 'V', 'R')
            else
                matrRefe = rigiAsse
                call dismoi('TYPE_MATRICE', massAsse, 'MATR_ASSE', repk=matrType)
                if (matrType .eq. 'NON_SYM') then
                    matrRefe = massAsse
                end if
                if (lDampMatrix) then
                    call dismoi('TYPE_MATRICE', dampAsse, 'MATR_ASSE', repk=matrType)
                    if (matrType .eq. 'NON_SYM') then
                        matrRefe = dampAsse
                    end if
                end if
                call mtdefs(matrAsse, matrRefe, 'V', 'R')
            end if
        end if
    end if

! - Combination of matrix
    if (lDyna) then
        call mtcmbl(nbMatr, coefType, coefVale, matrName, matrAsse, ' ', ' ', 'ELIM=')
    else
        matrAsse = rigiAsse
    end if
    if (phaseType .eq. ACCEL_INIT) then
        goto 999
    end if

! - Matrix for undead loads
    if (lNeumUndead) then
        call ascoma(meelem, numeDof, listLoad, matrAsse)
    end if

! - Matrix for friction
    if (lContDiscret .and. (phaseType .eq. CORR_NEWTON)) then
        call nmasfr(ds_contact, matrAsse)
    end if

! - Special post-treatment for LAC contact method
    if (lContLAC) then
        call lccmst(ds_contact, matrAsse)
    end if

! - Matrix for unilateral conidition
    if (lUnil .and. (phaseType .eq. CORR_NEWTON)) then
        lUnilPena = cfdisl(ds_contact%sdcont_defi, 'UNIL_PENA')
        if (lUnilPena) then
            call nmasun(ds_contact, matrAsse)
        end if
    end if
!
999 continue
!
end subroutine
