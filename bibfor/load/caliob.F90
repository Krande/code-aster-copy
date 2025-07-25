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
subroutine caliob(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8dgrd.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/char_excl_keyw.h"
#include "asterfort/char_read_keyw.h"
#include "asterfort/dismoi.h"
#include "asterfort/getnode.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/matrot.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'LIAISON_OBLIQUE'
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh        : mesh
! In  load        : load
! In  model       : model
! In  valeType    : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: n_max_keyword = 300
    integer(kind=8) :: ddlimp(n_max_keyword)
    real(kind=8) :: valimr(n_max_keyword)
    complex(kind=8) :: valimc(n_max_keyword)
    character(len=8) :: valimf(n_max_keyword)
    character(len=16) :: keywordlist(n_max_keyword)
!
    character(len=24) :: list_node
    integer(kind=8) :: jlino, nb_node
    integer(kind=8) :: ino
    integer(kind=8) :: geomDime, nbec
    integer(kind=8) :: nliai, nume_node
    integer(kind=8) :: i_angle, i_keyword, iocc, i_direct
    real(kind=8) :: coefr, val_r, direct(3)
    character(len=8) :: ddl, coeff, val_f
    complex(kind=8) :: coefc, val_c
    character(len=4) :: typcoe
    character(len=8) :: nomg
    character(len=8) :: name_node
    character(len=16) :: keywordfact, keyword
    integer(kind=8) :: n_keyword
    character(len=19) :: lisrel
    real(kind=8) :: matr_rota(3, 3), rdgd
    real(kind=8) :: zero
    real(kind=8) :: angl_naut(3)
    integer(kind=8) :: n_angle
    character(len=24) :: keywordexcl
    integer(kind=8) :: n_keyexcl
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    keywordfact = 'LIAISON_OBLIQUE'
    call getfac(keywordfact, nliai)
    if (nliai .eq. 0) goto 999
!
! - Initializations
!
    lisrel = '&&CALIOB.RLLISTE'
    zero = 0.d0
    coefc = (1.0d0, 0.0d0)
    coefr = 1.0d0
    coeff = ' '
    rdgd = r8dgrd()
!
    typcoe = 'REEL'
    if (valeType .eq. 'COMP') then
        ASSERT(.false.)
    end if
!
! - Create list of excluded keywords for using in char_read_keyw
!
    keywordexcl = '&&CALIOB.KEYWORDEXCL'
    call char_excl_keyw(keywordfact, keywordexcl, n_keyexcl)
!
! - Information about <GRANDEUR>
!
    nomg = 'DEPL_R'
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 10)

! - Model informations
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)
    if (.not. (geomDime .eq. 2 .or. geomDime .eq. 3)) then
        call utmess('F', 'CHARGES2_6')
    end if
!
! - Loop on factor keyword
!
    do iocc = 1, nliai
!
! ----- Read mesh affectation
!
        list_node = '&&CALIOB.LIST_NODE'
        call getnode(mesh, keywordfact, iocc, 'F', list_node, &
                     nb_node)
        call jeveuo(list_node, 'L', jlino)
!
! ----- Local orientation
!
        angl_naut(1) = zero
        angl_naut(2) = zero
        angl_naut(3) = zero
        call getvr8(keywordfact, 'ANGL_NAUT', iocc=iocc, nbval=3, vect=angl_naut, &
                    nbret=n_angle)
        do i_angle = 1, min(3, abs(n_angle))
            angl_naut(i_angle) = rdgd*angl_naut(i_angle)
        end do
        call matrot(angl_naut, matr_rota)
!
! ----- Read affected components and their values
!
        call char_read_keyw(keywordfact, iocc, valeType, n_keyexcl, keywordexcl, &
                            n_max_keyword, n_keyword, keywordlist, ddlimp, valimr, &
                            valimf, valimc)
!
        do i_keyword = 1, n_keyword
            keyword = keywordlist(i_keyword)
!
! --------- Values
!
            val_r = valimr(i_keyword)
            val_c = valimc(i_keyword)
            val_f = valimf(i_keyword)
            ASSERT(ddlimp(i_keyword) .eq. 1)
!
! --------- Which direction ?
!
            if ((keyword .eq. 'DX') .or. (keyword .eq. 'DRX')) then
                i_direct = 1
            else if ((keyword .eq. 'DY') .or. (keyword .eq. 'DRY')) then
                i_direct = 2
            else if ((keyword .eq. 'DZ') .or. (keyword .eq. 'DRZ')) then
                i_direct = 3
            else
                ASSERT(.false.)
            end if
            direct(1) = matr_rota(i_direct, 1)
            direct(2) = matr_rota(i_direct, 2)
            direct(3) = matr_rota(i_direct, 3)
!
! --------- Which kind of dof ?
!
            if ((keyword .eq. 'DX') .or. (keyword .eq. 'DY') .or. (keyword .eq. 'DZ')) then
                ddl = 'DEPL'
            else if ((keyword .eq. 'DRX') .or. (keyword .eq. 'DRY') .or. (keyword .eq. 'DRZ')) then
                ddl = 'ROTA'
            else
                ASSERT(.false.)
            end if
!
! --------- Affect in direction
!
            do ino = 1, nb_node
                nume_node = zi(jlino+ino-1)
                name_node = int_to_char8(nume_node)
                call afrela([coefr], [coefc], ddl, name_node, [geomDime], &
                            direct, 1, val_r, val_c, val_f, &
                            typcoe, valeType, 0.d0, lisrel)
            end do
        end do
!
        call jedetr(list_node)
    end do
!
! - Final linear relation affectation
!
    call aflrch(lisrel, load, 'LIN')
!
    call jedetc('V', '&&CALIOB.RLLISTE', 1)
    call jedetr(keywordexcl)
!
999 continue
    call jedema()
end subroutine
