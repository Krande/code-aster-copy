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
subroutine chflch(rigthe, vec2nd, listLoad)
!
    use listLoad_module
    use listLoad_type
!
    implicit none
!
#include "asterfort/asasve.h"
#include "asterfort/ascavc.h"
#include "asterfort/ascova.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ntdoch.h"
#include "asterfort/vechth.h"
#include "asterfort/vedith.h"
#include "asterfort/vtcreb.h"
#include "jeveux.h"
!   ------------------------------------------------------------------------------------
!   *** Subroutine arguments
!   ------------------------------------------------------------------------------------
!       --> Input variables
    character(len=8) :: rigthe
!       <-- Output variables
    character(len=24) :: vec2nd
    character(len=24), intent(out) :: listLoad
!
!   ------------------------------------------------------------------------------------
!   *** Definition of local variables
!   ------------------------------------------------------------------------------------
    integer(kind=8) :: neq, jnchtp, j2nd, nbEqua, iEqua, jndirp
    real(kind=8) :: timeCurr
    character(len=8) :: numedd, model
    character(len=24) :: vechtp, cncine
    character(len=24) :: loadNameJv, loadInfoJv, loadFuncJv
    character(len=24) :: vachtp, cnchtp, vediri, vadirp, cndirp, timeMap
    character(len=24), parameter :: mateco = " "
    type(ListLoad_Prep) :: listLoadPrep

!   ------------------------------------------------------------------------------------
    vediri = '&&VETDIR           .RELR'
    vechtp = '&&VETCHA           .RELR'
    listLoad = '&&OP0116_INF_CHARGE'
    cndirp = '    '
    cnchtp = '    '
    timeCurr = 0.d0
    timeMap = ' '
!
    call jemarq()
!
    call dismoi('NOM_MODELE', rigthe, 'MATR_ASSE', repk=model)
    call dismoi('NOM_NUME_DDL', rigthe, 'MATR_ASSE', repk=numedd)
!
    call vtcreb(vec2nd, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)

! - Create list of loads
    listLoadPrep%model = model
    call ntdoch(listLoadPrep, listLoad, "G")
    loadNameJv = listLoad(1:19)//'.LCHA'
    loadInfoJv = listLoad(1:19)//'.INFC'
    loadFuncJv = listLoad(1:19)//'.FCHA'

! - Compute dirichlet RHS
    call vedith(model, loadNameJv, loadInfoJv, timeMap, vediri)
    call asasve(vediri, numedd, 'R', vadirp)
    call ascova('D', vadirp, loadFuncJv, 'INST', timeCurr, &
                'R', cndirp)
    call jeveuo(cndirp(1:19)//'.VALE', 'L', jndirp)

! - Compute dirichlet kinematic elimination
    cncine = ' '
    call ascavc(loadNameJv, loadInfoJv, loadFuncJv, numedd, timeCurr, &
                cncine, l_hho_=ASTER_FALSE)

! - Compute Neumann RHS
    call vechth('STAT', &
                model, mateco, &
                loadNameJv, loadInfoJv, &
                timeCurr, &
                vechtp)
    call asasve(vechtp, numedd, 'R', vachtp)
    call ascova('D', vachtp, loadFuncJv, 'INST', timeCurr, &
                'R', cnchtp)
    call jeveuo(cnchtp(1:19)//'.VALE', 'L', jnchtp)
    call jedetr(vechtp)

! - Final RHS
    call jeveuo(vec2nd(1:19)//'.VALE', 'E', j2nd)
    call jelira(vec2nd(1:19)//'.VALE', 'LONMAX', nbEqua)
    do iEqua = 1, nbEqua
        zr(j2nd+iEqua-1) = zr(jnchtp+iEqua-1)+zr(jndirp+iEqua-1)
    end do
!
    call jedema()
!
end subroutine
