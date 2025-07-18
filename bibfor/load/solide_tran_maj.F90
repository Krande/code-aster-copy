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

subroutine solide_tran_maj(load_name, disp)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/cnocns.h"
#include "asterfort/jedema.h"
#include "asterfort/detrsd.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
!
    character(len=8), intent(in) :: load_name
    character(len=19), intent(in) :: disp
!
! --------------------------------------------------------------------------------------------------
!
! Loads - Computation
!
! Update dualized relations for non-linear Dirichlet boundary conditions (undead) - LIAISON_SOLIDE
!
! --------------------------------------------------------------------------------------------------
!
! In  load_name        : name of load
! In  disp             : displacements
!
!   La programmation de cette routine est directement liee a celle de solide_tran.
!   Les deux programmations doivent evoluer ensemble.
!
! --------------------------------------------------------------------------------------------------

    integer(kind=8) :: numnoe_a, numnoe_b, kdx, kdy, kdz, ncmp1
    character(len=8) :: nomnoe_a, nomnoe_b
    integer(kind=8) :: gd1, gd2, ncmpmx1, ncmpmx2, dim, ka, kb, nb_term
    real(kind=8) :: coer(6)
    real(kind=8) :: xa, ya, za, xb, yb, zb
    real(kind=8) :: ua, va, wa, ub, vb, wb
    real(kind=8) :: c01
    integer(kind=8) :: i_link, nb_link, ideb, nbterm1, k, ico
    character(len=8)  :: mesh
    character(len=13)  :: load_dual
    character(len=19)  :: deplas
    aster_logical :: l3d
    integer(kind=8), pointer :: dual_prdi(:) => null()
    character(len=8), pointer :: dual_prdk(:) => null()
    character(len=8), pointer :: load_type(:) => null()
    character(len=8), pointer :: dual_prdso(:) => null()
    character(len=8), pointer :: deplc(:) => null()
    real(kind=8), pointer :: deplv(:) => null()
    integer(kind=8), pointer :: depld(:) => null()
    aster_logical, pointer :: depll(:) => null()
    real(kind=8), pointer :: mesh_coor(:) => null()
    integer(kind=8), pointer :: cmult_d(:) => null()
    integer(kind=8), pointer :: cimpo_d(:) => null()
    real(kind=8), pointer :: cmult_v(:) => null()
    real(kind=8), pointer :: cimpo_v(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Access to mesh
!
    call dismoi('NOM_MAILLA', load_name, 'CHARGE', repk=mesh)
    call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=mesh_coor)
!
! - Access to load
!
    call jeveuo(load_name//'.TYPE', 'L', vk8=load_type)
    ASSERT(load_type(1) .eq. 'MECA_RE')
    load_dual = load_name//'.DUAL'
    call jeveuo(load_dual//'.PRDK', 'L', vk8=dual_prdk)
    call jelira(load_dual//'.PRDK', 'LONUTI', ival=nb_link)
    call jeveuo(load_dual//'.PRDI', 'L', vi=dual_prdi)
    call jeveuo(load_dual//'.PRDSO', 'L', vk8=dual_prdso)
    ASSERT(size(dual_prdi) .ge. 3*nb_link)
!
!   -- Acces aux cartes .CMULT et .CIMPO :
!   ---------------------------------------
    call jeveuo(load_name//'.CHME.CMULT.DESC', 'L', vi=cmult_d)
    call jeveuo(load_name//'.CHME.CIMPO.DESC', 'L', vi=cimpo_d)
    call jeveuo(load_name//'.CHME.CMULT.VALE', 'E', vr=cmult_v)
    call jeveuo(load_name//'.CHME.CIMPO.VALE', 'E', vr=cimpo_v)
    gd1 = cmult_d(1)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd1), 'LONMAX', ncmpmx1)
    gd2 = cimpo_d(1)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd2), 'LONMAX', ncmpmx2)
!
!   -- On transforme le champ de deplacement en champ simple :
!   -----------------------------------------------------------
    deplas = '&&solide_tran_maj.U'
    call cnocns(disp, 'V', deplas)
    call jeveuo(deplas//'.CNSD', 'L', vi=depld)
    call jeveuo(deplas//'.CNSC', 'L', vk8=deplc)
    call jeveuo(deplas//'.CNSV', 'L', vr=deplv)
    call jeveuo(deplas//'.CNSL', 'L', vl=depll)
    kdx = 1
    kdy = 2
    kdz = 3
    ncmp1 = size(deplc)
    ASSERT(deplc(kdx) .eq. 'DX')
    ASSERT(deplc(kdy) .eq. 'DY')

!   -- on boucle sur les liaisons de load_type '2Dx' ou '3Dx' :
!   -------------------------------------------------------
    do i_link = 1, nb_link
        if (dual_prdk(i_link) (1:2) .ne. '2D' .and. dual_prdk(i_link) (1:2) .ne. '3D') cycle

        l3d = dual_prdk(i_link) (1:1) .eq. '3'
        if (l3d) then
            ASSERT(deplc(kdz) .eq. 'DZ')
        end if
        read (dual_prdk(i_link) (3:3), '(I1)') dim
        ASSERT(dim .ge. 1 .and. dim .le. 3)

        ideb = dual_prdi(3*(i_link-1)+3)
        nbterm1 = dual_prdi(3*(i_link-1)+2)

!       TRES IMPORTANT :
!       ----------------
!       -- Il est indispensable de parcourir les relations non lineaires dans le meme ordre
!       qu'a la creation dans solide_tran.
!       Puis, pour une relation donnee, il faut parcourir les ddls dans le meme ordre.

        ico = ideb-1

        do ka = 1, dim
            nomnoe_a = dual_prdso(4*(i_link-1)+ka)
            numnoe_a = char8_to_int(nomnoe_a)
            xa = mesh_coor(3*(numnoe_a-1)+1)
            ya = mesh_coor(3*(numnoe_a-1)+2)
            if (l3d) za = mesh_coor(3*(numnoe_a-1)+3)

            ua = deplv((numnoe_a-1)*ncmp1+kdx)
            va = deplv((numnoe_a-1)*ncmp1+kdy)
            if (l3d) wa = deplv((numnoe_a-1)*ncmp1+kdz)

            do kb = ka+1, dim+1
                nomnoe_b = dual_prdso(4*(i_link-1)+kb)
                numnoe_b = char8_to_int(nomnoe_b)
                xb = mesh_coor(3*(numnoe_b-1)+1)
                yb = mesh_coor(3*(numnoe_b-1)+2)
                if (l3d) zb = mesh_coor(3*(numnoe_b-1)+3)

                ub = deplv((numnoe_b-1)*ncmp1+kdx)
                vb = deplv((numnoe_b-1)*ncmp1+kdy)
                if (l3d) wb = deplv((numnoe_b-1)*ncmp1+kdz)

                if (l3d) then
                    nb_term = 6
                else
                    nb_term = 4
                end if

!               -- Relation: AB^2 = cste

!               -- Ordre : A,   B,     A,   B     A,   B
!                         'DX','DX',  'DY','DY' ,'DZ','DZ'

                c01 = (yb-ya+vb-va)**2-(yb-ya)**2+(xb-xa+ub-ua)**2-(xb-xa)**2
                if (l3d) c01 = c01+(zb-za+wb-wa)**2-(zb-za)**2

                coer(1) = -2*(xb-xa+ub-ua)
                coer(2) = 2*(xb-xa+ub-ua)
                coer(3) = -2*(yb-ya+vb-va)
                coer(4) = 2*(yb-ya+vb-va)
                if (l3d) then
                    coer(5) = -2*(zb-za+wb-wa)
                    coer(6) = 2*(zb-za+wb-wa)
                end if

                do k = 1, nb_term
                    ico = ico+1
                    cmult_v(ncmpmx1*(ico-1)+1) = coer(k)
                end do
                if (i_link .eq. 1) then
                    ASSERT(ico .le. nbterm1)
                end if
!               -- le coefficient constant est affecte avec le dernier terme de la relation :
                cimpo_v(ncmpmx2*(ico-1)+1) = -c01

            end do
        end do
    end do

    call detrsd('CHAM_NO_S', deplas)
!
    call jedema()
end subroutine
