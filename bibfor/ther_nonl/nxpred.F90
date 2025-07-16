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
! aslint: disable=W1504
!
subroutine nxpred(model, mateco, caraElem, listLoad, nume_dof, &
                  solver, l_stat, tpsthe, timeMap, matass, &
                  lonch, maprec, varc_curr, temp_prev, temp_iter, &
                  cn2mbr, hydr_prev, hydr_curr, &
                  comporTher, cndirp, cnchci, vec2nd, vec2ni, &
                  ds_algorom)
!
    use ROM_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asasve.h"
#include "asterfort/ascova.h"
#include "asterfort/copisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nxreso.h"
#include "asterfort/romAlgoNLSystemSolve.h"
#include "asterfort/romAlgoNLCorrEFMatrixModify.h"
#include "asterfort/romAlgoNLCorrEFResiduModify.h"
#include "asterfort/verstp.h"
#include "asterfort/vethbt.h"
#include "asterfort/vethbu.h"
!
    character(len=8), intent(in) :: model, caraElem
    character(len=24), intent(in) :: mateco
    character(len=24), intent(in) :: listLoad
    character(len=24), intent(in) :: nume_dof
    character(len=19), intent(in) :: solver
    real(kind=8), intent(in) :: tpsthe(6)
    character(len=24), intent(in) :: timeMap
    character(len=19), intent(in) :: varc_curr
    integer(kind=8) :: lonch
    character(len=19) :: maprec
    character(len=24) :: matass, cndirp, cnchci, cnresi
    character(len=24) :: temp_iter, temp_prev, vec2nd, vec2ni
    character(len=24) :: hydr_prev, hydr_curr, comporTher
    aster_logical :: l_stat
    character(len=24), intent(in) :: cn2mbr
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
!
! --------------------------------------------------------------------------------------------------
!
! COMMANDE THER_NON_LINE : PHASE DE PREDICTION
!
! --------------------------------------------------------------------------------------------------
!
!     VAR temp_iter : ITERE PRECEDENT DU CHAMP DE TEMPERATURE
!
! In  cn2mbr : name of vector for second member
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iEqua, iMode
    real(kind=8) :: rbid
    real(kind=8) :: time_curr
    character(len=1) :: typres
    character(len=19) :: chsol
    character(len=24) :: bidon, veresi, varesi, vabtla, vebtla
    character(len=24) :: vebuem, vabuem, cnvabt, cnvabu
    character(len=24) :: loadNameJv, loadInfoJv
    real(kind=8), pointer :: v_vec2ni(:) => null()
    real(kind=8), pointer :: v_vec2nd(:) => null()
    real(kind=8), pointer :: v_cn2mbr(:) => null()
    real(kind=8), pointer :: v_cnvabt(:) => null()
    real(kind=8), pointer :: v_cnvabu(:) => null()
    real(kind=8), pointer :: v_cndirp(:) => null()
    real(kind=8), pointer :: v_cnresi(:) => null()
    real(kind=8), pointer :: v_gamma(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    varesi = '&&VARESI'
    cnvabt = ' '
    cnvabu = ' '
    typres = 'R'
    chsol = '&&NXPRED.SOLUTION'
    bidon = '&&FOMULT.BIDON'
    veresi = '&&VERESI           .RELR'
    vebtla = '&&VETBTL           .RELR'
    vabtla = ' '
    vebuem = '&&VEBUEM           .RELR'
    vabuem = ' '
    cnresi = ' '
    time_curr = tpsthe(1)

! - Access to datastructure of list of loads
    loadNameJv = listLoad(1:19)//'.LCHA'
    loadInfoJv = listLoad(1:19)//'.INFC'
!
! --- RECUPERATION D'ADRESSES
!
    call jeveuo(vec2nd(1:19)//'.VALE', 'L', vr=v_vec2nd)
    call jeveuo(cn2mbr(1:19)//'.VALE', 'E', vr=v_cn2mbr)
    call jeveuo(vec2ni(1:19)//'.VALE', 'L', vr=v_vec2ni)
    call jeveuo(cndirp(1:19)//'.VALE', 'L', vr=v_cndirp)
!
    if (l_stat) then
!
!=======================================================================
!  INITIALISATION POUR LE PREMIER PAS DE CALCUL
!=======================================================================
! ----- Compute residual vector (non-linear) - Material and loads
        call verstp(l_stat, &
                    model, caraElem, mateco, &
                    loadNameJv, loadInfoJv, &
                    tpsthe, timeMap, temp_prev, temp_iter, &
                    varc_curr, comporTher, &
                    hydr_prev, hydr_curr, veresi, "V")
        call asasve(veresi, nume_dof, typres, varesi)
        call ascova('D', varesi, bidon, 'INST', rbid, &
                    typres, cnresi)
        call jeveuo(cnresi(1:19)//'.VALE', 'L', vr=v_cnresi)
!
! ----- BT LAMBDA - CALCUL ET ASSEMBLAGE
!
        call vethbt(model, loadNameJv, loadInfoJv, &
                    temp_prev, vebtla, 'V')
        call asasve(vebtla, nume_dof, typres, vabtla)
        call ascova('D', vabtla, bidon, 'INST', rbid, &
                    typres, cnvabt)
        call jeveuo(cnvabt(1:19)//'.VALE', 'L', vr=v_cnvabt)
!
! ----- B . TEMPERATURE - CALCUL ET ASSEMBLAGE
!
        call vethbu(model, matass, loadNameJv, loadInfoJv, &
                    temp_prev, vebuem)
        call asasve(vebuem, nume_dof, typres, vabuem)
        call ascova('D', vabuem, bidon, 'INST', rbid, &
                    typres, cnvabu)
        call jeveuo(cnvabu(1:19)//'.VALE', 'L', vr=v_cnvabu)
!
        do iEqua = 1, lonch
            v_cn2mbr(iEqua) = v_vec2nd(iEqua)-v_cnresi(iEqua)+ &
                              v_cndirp(iEqua)-v_cnvabt(iEqua)-v_cnvabu(iEqua)
        end do
!
! ----- Solve linear system
!
        if (ds_algorom%l_rom .and. ds_algorom%phase .eq. 'HROM') then
            call jeveuo(ds_algorom%gamma, 'E', vr=v_gamma)
            do iMode = 1, ds_algorom%ds_empi%nbMode
                v_gamma(iMode) = 0.d0
            end do
            call copisd('CHAMP_GD', 'V', temp_prev, chsol)
            call romAlgoNLSystemSolve(matass, cn2mbr, cnchci, ds_algorom, chsol)
        else if (ds_algorom%l_rom .and. ds_algorom%phase .eq. 'CORR_EF') then
            call romAlgoNLCorrEFMatrixModify(nume_dof, matass, ds_algorom)
            call romAlgoNLCorrEFResiduModify(cn2mbr, ds_algorom)
            call nxreso(matass, maprec, solver, cnchci, cn2mbr, &
                        chsol)
        else
            call nxreso(matass, maprec, solver, cnchci, cn2mbr, &
                        chsol)
        end if
!
! --- RECOPIE DANS temp_iter DU CHAMP SOLUTION CHSOL
!
        call copisd('CHAMP_GD', 'V', chsol, temp_iter)
!
    else
!
!=======================================================================
!  INITIALISATION POUR LE PREMIER PAS, CALCUL TRANSITOIRE, PAS COURANT
!=======================================================================
!
        do iEqua = 1, lonch
            v_cn2mbr(iEqua) = v_vec2ni(iEqua)+v_cndirp(iEqua)
        end do
!
! ----- Solve linear system
!
        if (ds_algorom%l_rom .and. ds_algorom%phase .eq. 'HROM') then
            call jeveuo(ds_algorom%gamma, 'E', vr=v_gamma)
            do iMode = 1, ds_algorom%ds_empi%nbMode
                v_gamma(iMode) = 0.d0
            end do
            call copisd('CHAMP_GD', 'V', temp_prev, chsol)
            call romAlgoNLSystemSolve(matass, cn2mbr, cnchci, ds_algorom, chsol)
        else if (ds_algorom%l_rom .and. ds_algorom%phase .eq. 'CORR_EF') then
            call romAlgoNLCorrEFMatrixModify(nume_dof, matass, ds_algorom)
            call romAlgoNLCorrEFResiduModify(cn2mbr, ds_algorom)
            call nxreso(matass, maprec, solver, cnchci, cn2mbr, &
                        chsol)
        else
            call nxreso(matass, maprec, solver, cnchci, cn2mbr, &
                        chsol)
        end if

!
! --- RECOPIE DANS temp_iter DU CHAMP SOLUTION CHSOL
!
        call copisd('CHAMP_GD', 'V', chsol, temp_iter)
!
    end if
!
    call jedema()
end subroutine
