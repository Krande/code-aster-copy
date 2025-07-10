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

subroutine cagrou(load, mesh, vale_type, phenom)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/getnode.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
!  Person in charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: load
    character(len=8), intent(in) :: mesh
    character(len=4), intent(in) :: vale_type
    character(len=4), intent(in) :: phenom
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'LIAISON_UNIF'
!
! --------------------------------------------------------------------------------------------------
!
!
! In mesh      : name of mesh
! In load      : name of load
! In vale_type : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_term
    parameter(nb_term=2)
    real(kind=8) :: coef_real_unit(nb_term)
    complex(kind=8) :: coef_cplx_unit(nb_term)
    character(len=8) :: dof_name(nb_term)
    character(len=8) :: node_name(nb_term)
    integer(kind=8) :: node_nume(nb_term)
    integer(kind=8) :: repe_type(nb_term)
    real(kind=8) :: repe_defi(3, nb_term)
    character(len=4) :: coef_type
    character(len=16) :: keywordfact
    character(len=19) :: list_rela
    integer(kind=8) :: ibid
    integer(kind=8) :: nliai
    real(kind=8) :: vale_real_zero
    character(len=8) :: vale_func_zero
    complex(kind=8) :: vale_cplx_zero
    integer(kind=8) :: iocc, i_no, i_dof
    character(len=24) :: list_node
    integer(kind=8) :: jlino
    integer(kind=8) :: nb_node
    character(len=24) :: list_dof
    integer(kind=8) :: jlidof
    integer(kind=8) :: nb_dof, ier
    aster_logical :: lcolle
!
    data repe_type/0, 0/
    data repe_defi/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    keywordfact = 'LIAISON_UNIF'
    call getfac(keywordfact, nliai)
    if (nliai .eq. 0) goto 999
!
! - Initializations
!
    vale_func_zero = '&FOZERO'
    vale_cplx_zero = (0.d0, 0.d0)
    vale_real_zero = 0.d0
    coef_cplx_unit(1) = (1.0d0, 0.0d0)
    coef_cplx_unit(2) = (-1.0d0, 0.0d0)
    coef_real_unit(1) = 1.d0
    coef_real_unit(2) = -1.d0
    list_rela = '&&CAGROU.RLLISTE'
    list_dof = '&&CAGROU.LIST_DOF'
    lcolle = .false.
    call jeexin(mesh//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
!
! - Initializations of types
!
    if (vale_type .eq. 'COMP') then
        ASSERT(.false.)
    else if (vale_type .eq. 'REEL') then
        coef_type = 'REEL'
    else if (vale_type .eq. 'FONC') then
        coef_type = 'REEL'
    else
        ASSERT(.false.)
    end if
!
    do iocc = 1, nliai
!
! ----- Read mesh affectation
!
        list_node = '&&CAGROU.LIST_NODE'
        call getnode(mesh, keywordfact, iocc, 'F', list_node, &
                     nb_node)
        if (nb_node .lt. 2) then
            call utmess('F', 'CHARGES2_82')
        end if
        call jeveuo(list_node, 'L', jlino)
!
! ----- Get dof
!
        call getvtx(keywordfact, 'DDL', iocc=iocc, nbval=0, nbret=nb_dof)
        ASSERT(nb_dof .ne. 0)
        nb_dof = -nb_dof
        call wkvect(list_dof, 'V V K8', nb_dof, jlidof)
        call getvtx(keywordfact, 'DDL', iocc=iocc, nbval=nb_dof, vect=zk8(jlidof), &
                    nbret=ibid)
!
! ----- First node
!
        node_nume(1) = zi(jlino-1+1)
        node_name(1) = int_to_char8(node_nume(1), lcolle, mesh, 'NOEUD')
!
! ----- Loop on dof
!
        do i_dof = 1, nb_dof
            dof_name(1) = zk8(jlidof-1+i_dof)
            dof_name(2) = zk8(jlidof-1+i_dof)
            do i_no = 2, nb_node
                node_nume(2) = zi(jlino-1+i_no)
                node_name(2) = int_to_char8(node_nume(2), lcolle, mesh, 'NOEUD')
                call afrela(coef_real_unit, coef_cplx_unit, dof_name, node_name, repe_type, &
                            repe_defi, nb_term, vale_real_zero, vale_cplx_zero, vale_func_zero, &
                            coef_type, vale_type, 0.d0, list_rela)
            end do
        end do
!
        call jedetr(list_node)
        call jedetr(list_dof)
    end do
!
! - Final linear relation affectation
!
    if (phenom .eq. 'MECA') then
    end if
    call aflrch(list_rela, load, 'LIN')
!
999 continue
    call jedema()
end subroutine
