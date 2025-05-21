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

subroutine rco3d_apco3d(noma, lismavo, lismaco, nbmavo, nbmaco, epai, &
                        list_pairs, nb_pairs, nt_nodes)
!
    use raco3d_module
    implicit none
!

#include "jeveux.h"
#include "asterfort/jelira.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterc/r8prem.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"

    character(len=8), intent(in) :: noma
    character(len=24), intent(in) :: lismaco, lismavo
    integer, intent(in) :: nbmavo, nbmaco
    real(kind=8), intent(in) :: epai
    integer, intent(out) :: nb_pairs, nt_nodes
    integer, pointer :: list_pairs(:)

! Appariemment COQUE-3D
!
! --------------------------------------------------------------------------------------------------
!
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: desc, vale, connex, typemail
    character(len=8) :: elem_type
    integer ::  i, j, k, l, m, idvale
    integer :: imavo, idconnvo, nbnovo, inovo1, inovo2, jlismavo
    integer :: imaco, idconnco, nbnoco, inoco1, inoco2, jlismaco
    integer :: mesh_typmail, elem_type_nume
    aster_logical :: check
    real(kind=8) :: xv(3), x1, y1, z1, x2, y2, z2, x3, y3, z3
    real(kind=8) :: tab1(3), tab2(3), diam, dist
    real(kind=8) :: coorsegvo(3, 2), coorsegco(3, 2)

    call jemarq()

! INITIALISATIONS
    AS_ALLOCATE(vi=list_pairs, size=3*nbmaco*nbmavo)
    nb_pairs = 0
    nt_nodes = 0

    desc = noma(1:8)//'.COORDO    .DESC'
    vale = noma(1:8)//'.COORDO    .VALE'
    connex = noma(1:8)//'.CONNEX'
    typemail = noma(1:8)//'.TYPMAIL'

    call jeveuo(vale, 'L', idvale)
    call jeveuo(lismavo, 'L', jlismavo)
    call jeveuo(lismaco, 'L', jlismaco)
    call jeveuo(typemail, 'L', mesh_typmail)

    check = ASTER_FALSE
    do i = 1, nbmavo
        imavo = zi(jlismavo-1+i)
        call jeveuo(jexnum(connex, imavo), 'L', idconnvo)
        call jelira(jexnum(connex, imavo), 'LONMAX', nbnovo)
        do l = 1, nbmaco
            imaco = zi(jlismaco-1+l)
            !elem_type_nume = zi(mesh_typmail-1+imaco)
            !call jenuno(jexnum('&CATA.TM.NOMTM', elem_type_nume), elem_type)
            call jeveuo(jexnum(connex, imaco), 'L', idconnco)
            call jelira(jexnum(connex, imaco), 'LONMAX', nbnoco)

            do j = 1, nbnovo-1
                inovo1 = zi(idconnvo+j-1)
                inovo2 = zi(idconnvo+j+1-1)
                do k = 1, 3
                    coorsegvo(k, 1) = zr(idvale-1+3*(inovo1-1)+k)
                    coorsegvo(k, 2) = zr(idvale-1+3*(inovo2-1)+k)
                end do
                do m = 1, nbnoco-1
                    inoco1 = zi(idconnco+m-1)
                    inoco2 = zi(idconnco+m+1-1)
                    do k = 1, 3
                        coorsegco(k, 1) = zr(idvale-1+3*(inoco1-1)+k)
                        coorsegco(k, 2) = zr(idvale-1+3*(inoco2-1)+k)
                    end do
                    !!
                    dist = segseg_distance(coorsegvo, coorsegco)
                    if (dist .le. 0.5d0*epai) then
                        nb_pairs = nb_pairs+1
                        list_pairs(2*(nb_pairs-1)+1) = imaco
                        list_pairs(2*(nb_pairs-1)+2) = imavo
                        nt_nodes = nt_nodes+nbnoco+nbnovo
                        check = ASTER_TRUE
                        exit
                    end if
                end do
                if (check) then
                    ! inutile de comparer les autres segments:
                    ! passer  à l'element suivant
                    exit
                end if
            end do
            check = ASTER_FALSE
        end do
    end do

    call jedema()

end subroutine
