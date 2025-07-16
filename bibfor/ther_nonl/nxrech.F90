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
subroutine nxrech(model, mateco, caraElem, listLoad, nume_dof, &
                  tpsthe, timeMap, lonch, comporTher, varc_curr, &
                  temp_iter, vtempp, vtempr, temp_prev, hydr_prev, &
                  hydr_curr, vec2nd, cnvabt, &
                  cnresi, rho, iterho, ds_algopara, l_stat)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/asasve.h"
#include "asterfort/ascova.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/verstp.h"
#include "asterfort/vethbt.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: model, caraElem
    character(len=24), intent(in) :: mateco
    character(len=24), intent(in) :: listLoad
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    real(kind=8), intent(in) :: tpsthe(6)
    character(len=24), intent(in) :: timeMap
    character(len=19), intent(in) :: varc_curr
    integer(kind=8) :: lonch
    real(kind=8) :: rho
    character(len=24) :: temp_prev, vtempr, vtempp, temp_iter, cnvabt, cnresi, vec2nd
    character(len=24) :: hydr_prev, hydr_curr, comporTher
    aster_logical, intent(in) :: l_stat
!
! --------------------------------------------------------------------------------------------------
!
! COMMANDE THER_NON_LINE : RECHERCHE LINEAIRE
! DANS LA DIRECTION DONNEE PAR NEWTON (ON CHERCHE RHO).
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_algopara      : datastructure for algorithm parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: j2nd, jvare, jbtla, rang, i
    real(kind=8) :: rho0, rhot, f0, f1, rhomin, rhomax
    real(kind=8) :: rhof, ffinal
    real(kind=8) :: testm, r8bid, resi_i
    character(len=24) :: vebtla, veresi, varesi, bidon, vabtla
    character(len=1) :: typres
    character(len=8) :: mesh
    character(len=19) :: nume_equa
    integer(kind=8) :: itrmax, iterho
    real(kind=8) :: time_curr
    real(kind=8), pointer :: tempm(:) => null()
    real(kind=8), pointer :: tempp(:) => null()
    real(kind=8), pointer :: tempr(:) => null()
    aster_logical, pointer :: v_pddl(:) => null()
    integer(kind=8), pointer :: v_posdd(:) => null()
    character(len=24) :: loadNameJv, loadInfoJv
    parameter(rhomin=-2.d0, rhomax=2.d0)
    mpi_int :: mrank
    data typres/'R'/
    data bidon/'&&FOMULT.BIDON'/
    data vebtla/'&&VETBTL           .RELR'/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    varesi = '&&VARESI'
    veresi = '&&VERESI'
    time_curr = tpsthe(1)

! - Access to datastructure of list of loads
    loadNameJv = listLoad(1:19)//'.LCHA'
    loadInfoJv = listLoad(1:19)//'.INFC'
!
    call asmpi_info(rank=mrank)
    rang = to_aster_int(mrank)
!
! --- RECUPERATION D'ADRESSES JEVEUX
!
    call jeveuo(temp_iter(1:19)//'.VALE', 'E', vr=tempm)
    call jeveuo(vtempp(1:19)//'.VALE', 'E', vr=tempp)
    call jeveuo(vtempr(1:19)//'.VALE', 'E', vr=tempr)
    call jeveuo(vec2nd(1:19)//'.VALE', 'L', j2nd)
    call jeveuo(cnresi(1:19)//'.VALE', 'L', jvare)
    call jeveuo(cnvabt(1:19)//'.VALE', 'L', jbtla)
!
    call wkvect("&&NXRECH.PDDL", "V V L", lonch, vl=v_pddl)

    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    if (isParallelMesh(mesh)) then
        call dismoi('NUME_EQUA', nume_dof, 'NUME_DDL', repk=nume_equa)
        call jeveuo(nume_equa//'.PDDL', 'L', vi=v_posdd)
!
        do i = 1, lonch
            v_pddl(i) = v_posdd(i) == rang
        end do
    else
        v_pddl = ASTER_TRUE
    end if

!
! --- RECHERCHE LINEAIRE (CALCUL DE RHO) SUR L'INCREMENT VTEMPP
!
    f0 = 0.d0
    do i = 1, lonch
        if (v_pddl(i)) then
            resi_i = zr(j2nd+i-1)-zr(jvare+i-1)-zr(jbtla+i-1)
            f0 = f0+tempp(i)*resi_i
        end if
    end do
    call asmpi_comm_vect('MPI_SUM', 'R', 1, scr=f0)
!
    rho0 = 0.d0
    rho = 1.d0
    itrmax = ds_algopara%line_search%iter_maxi+1
    do iterho = 1, itrmax
        do i = 1, lonch
            tempr(i) = tempm(i)+rho*tempp(i)
        end do
! ----- Compute residual vector (non-linear) - Material and loads
        call verstp(l_stat, &
                    model, caraElem, mateco, &
                    loadNameJv, loadInfoJv, &
                    tpsthe, timeMap, temp_prev, vtempr, &
                    varc_curr, comporTher, &
                    hydr_prev, hydr_curr, veresi, "V")
        call asasve(veresi, nume_dof, typres, varesi)
        call ascova('D', varesi, bidon, 'INST', r8bid, &
                    typres, cnresi)
        call jeveuo(cnresi(1:19)//'.VALE', 'L', jvare)
!
! --- BT LAMBDA - CALCUL ET ASSEMBLAGE
!
        call vethbt(model, loadNameJv, loadInfoJv, &
                    vtempr, vebtla, 'V')
        call asasve(vebtla, nume_dof, typres, vabtla)
        call ascova('D', vabtla, bidon, 'INST', r8bid, &
                    typres, cnvabt)
        call jeveuo(cnvabt(1:19)//'.VALE', 'L', jbtla)
!
        f1 = 0.d0
        testm = 0.d0
!
        do i = 1, lonch
            if (v_pddl(i)) then
                resi_i = zr(j2nd+i-1)-zr(jvare+i-1)-zr(jbtla+i-1)
                f1 = f1+tempp(i)*resi_i
                testm = max(testm, abs(resi_i))
            end if
        end do
!
        call asmpi_comm_vect('MPI_MAX', 'R', 1, scr=testm)
        call asmpi_comm_vect('MPI_SUM', 'R', 1, scr=f1)
!
        if (testm .lt. ds_algopara%line_search%resi_rela) then
            goto 999
        end if
!
        if (iterho .eq. 1) then
            ffinal = f1
            rhof = 1.d0
        end if
        if (abs(f1) .lt. abs(ffinal)) then
            ffinal = f1
            rhof = rho
        end if
        rhot = rho
        if (abs(f1-f0) .gt. r8prem()) then
            rho = -(f0*rhot-f1*rho0)/(f1-f0)
            if (rho .lt. rhomin) rho = rhomin
            if (rho .gt. rhomax) rho = rhomax
            if (abs(rho-rhot) .lt. 1.d-08) then
                goto 40
            end if
        else
            goto 40
        end if
        rho0 = rhot
        f0 = f1
    end do
40  continue
    rho = rhof
    f1 = ffinal
!
999 continue
    iterho = iterho-1
    call jedetr("&&NXRECH.PDDL")
    call jedema()
end subroutine
