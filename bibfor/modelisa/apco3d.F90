! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1501
!

subroutine apco3d(noma, lismavo, lismaco, nbmavo, nbmaco, &
                list_pairs, nb_pairs, nt_nodes)
    !
    implicit none
    !

#include "jeveux.h"
#include "asterfort/jelira.h"
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
    integer :: imavo, idconnvo, nbnovo, inovo, jlismavo
    integer :: imaco, idconnco, nbnoco, inoco, jlismaco
    integer :: mesh_typmail, elem_type_nume
    aster_logical :: check
    real(kind=8) :: xv(3), x1, y1, z1, x2, y2, z2, x3, y3, z3
    real(kind=8) :: tab1(3), tab2(3), diam, dist, prec


! INITIALISATIONS
    AS_ALLOCATE(vi=list_pairs, size=nbmaco*nbmavo)
    nb_pairs = 0
    nt_nodes = 0

    prec = r8prem()
    !prec = 1.d-4

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
            do j = 1, nbnovo
                inovo = zi(idconnvo+j-1)
                do k = 1, 3
                    xv(k) = zr(idvale-1+3*(inovo-1)+k)
                end do 
                if (.not. check) then
                    imaco = zi(jlismaco-1+l)
                    !elem_type_nume = zi(mesh_typmail-1+imaco)
                    !call jenuno(jexnum('&CATA.TM.NOMTM', elem_type_nume), elem_type)
                    call jeveuo(jexnum(connex, imaco), 'L', idconnco)
                    call jelira(jexnum(connex, imaco), 'LONMAX', nbnoco)
                    do m = 1, nbnoco
                        inoco = zi(idconnco+m-1)
                        
                        if ( m .eq. 1) then
                            x1 = zr(idvale-1+3*(inoco-1)+1)
                            y1 = zr(idvale-1+3*(inoco-1)+2)
                            z1 = zr(idvale-1+3*(inoco-1)+3)
                        end if
                        if ( m .eq. 2) then
                            x2 = zr(idvale-1+3*(inoco-1)+1)
                            y2 = zr(idvale-1+3*(inoco-1)+2)
                            z2 = zr(idvale-1+3*(inoco-1)+3)
                        end if
                        if ( m .eq. 3) then
                            x3 = zr(idvale-1+3*(inoco-1)+1)
                            y3 = zr(idvale-1+3*(inoco-1)+2)
                            z3 = zr(idvale-1+3*(inoco-1)+3)
                        end if

                        !do k = 1, 3
                        !    xc(k) = zr(idvale-1+3*(inoco-1)+k)
                        !end do 
                    end do
                    if (nbnoco .eq. 2) then
                        tab1(1) = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                        diam = sqrt(tab1(1))
                    end if
                    if (nbnoco .eq. 3) then
                        tab1(1) = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                        tab1(2) = (x3-x2)**2+(y3-y2)**2+(z3-z2)**2
                        tab1(3) = (x1-x3)**2+(y1-y3)**2+(z1-z3)**2
                        diam = sqrt(max(tab1(1), tab1(2), tab1(3)))
                    end if
                    check = ASTER_TRUE
                end if
                if (nbnoco .eq. 2) then
                    tab2(1) = (xv(1)-x1)**2+(xv(2)-y1)**2+(xv(3)-z1)**2
                    tab2(2) = (xv(1)-x2)**2+(xv(2)-y2)**2+(xv(3)-z2)**2
                    dist = sqrt(min(tab2(1), tab2(2)))
                end if
                if (dist .le. prec * diam) then
                    nb_pairs = nb_pairs + 1
                    list_pairs(2*(nb_pairs - 1) + 1) = imaco
                    list_pairs(2*(nb_pairs - 1) + 2) = imavo
                    nt_nodes = nt_nodes + nbnoco + nbnovo
                    !write(*,*) "########## ", imaco, "  ###########  ", imavo
                    exit
                end if
            end do
            check = ASTER_FALSE
        end do
    end do

end subroutine